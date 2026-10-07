"""C3 guest check — web app (Streamlit).

Run locally:   streamlit run app.py
"""
from pathlib import Path

import pandas as pd
import streamlit as st

import c3model_core as core

MAX_MOLECULES = 25
HERE = Path(__file__).parent

st.set_page_config(page_title="C3 guest check", page_icon="⚗️", layout="centered")

if core._host_xyz is None:
    core.load_host((HERE / "host_site.txt").read_text())

MARK = {'✓': '✔️', '✗': '❌', '⚠': '⚠️', '~': '➖', '–': '▫️'}
BOX = {'fits': st.success, 'no': st.error, 'beyond': st.warning, 'border': st.info,
       'outside': st.info, 'unreadable': st.info}


def show(r):
    head = f"**{core.ICON[r['verdict']]} {core.LABEL[r['verdict']]}**"
    if r.get('reasons') and r['verdict'] != 'fits':
        head += " — " + "; ".join(r['reasons'])
    BOX[r['verdict']](f"{head}  \n`{r['smiles']}`")
    if r['status'] == 'done':
        items = [f"Ends: **{r['ends'][0]} … {r['ends'][1]}**, {r['L']:.2f} Å apart "
                 f"(most extended shape within 3 kcal/mol)"]
        items += [f"{MARK.get(m, m)} **{name}** — {text}" for m, name, text in core.checks(r)]
        st.markdown("\n".join(f"- {i}" for i in items))
    for w in r['why']:
        st.markdown(w)


st.title("C3 guest check")
st.markdown("Will a small molecule form the inclusion phase with the PET cyclic trimer (C3)? "
            "Type one or more SMILES and click **Check**.")

with st.form("check"):
    text = st.text_input("SMILES", placeholder="e.g.  Cl/C=C/Cl   ClCCCl   FCCF",
                         help="Several molecules: separate them with spaces or commas.")
    rule = st.checkbox("Carbon ends cannot hold the pocket (recommended)", value=True,
                       help="The pocket is lined by three C–H groups; the end atom needs lone pairs "
                            "(halogen, O, N, S) to accept them.")
    go = st.form_submit_button("Check", type="primary")

if go:
    items = core.split_input(text)
    if not items:
        st.info("Type at least one SMILES.")
    else:
        if len(items) > MAX_MOLECULES:
            st.warning(f"Only the first {MAX_MOLECULES} molecules are checked.")
            items = items[:MAX_MOLECULES]
        results = []
        for smi in items:
            with st.spinner(f"Checking {smi} …"):
                r = core.judge(core.analyse(smi), carbon_end_rule=rule)
            results.append(r)
            st.divider()
            show(r)
        if len(results) > 1:
            st.divider()
            st.subheader("Summary")
            st.dataframe(pd.DataFrame(
                [{"Result": f"{core.ICON[r['verdict']]} {core.LABEL[r['verdict']]}",
                  "SMILES": r['smiles'],
                  "Reason": "; ".join(r['reasons']) if r['verdict'] != 'fits' else ""} for r in results]),
                hide_index=True, width="stretch")

st.divider()
with st.expander("What the results mean"):
    st.markdown("""
* ✅ **FITS** – both ends sit in the pockets of two C3 rings and nothing pushes into the walls.
* ❌ **DOES NOT FIT** – with the reason: an end does not reach its pocket, the molecule hits the walls,
  a carbon end, or the molecule is too large.
* ⚠️ **BEYOND TESTED RANGE** – the ends are larger than bromine; the model cannot say yet.
* ➖ **BORDERLINE** – close to the limits set by the known guests.
* ❔ **OUTSIDE THE MODEL** – salts, mixtures, charged species, unusual elements.
""")
with st.expander("How the model works"):
    st.markdown("""
In the inclusion phase each end of the guest sits in a pocket inside one C3 ring, lined by three aromatic
C–H groups. The two pockets belong to two neighbouring C3 rings; between them is a flat slot.

The host is taken from the DFT-relaxed trans-1,2-dichloroethylene@C3 structure (PBE-D3, VASP) and kept frozen.
For each SMILES the guest is built in 3D (RDKit, MMFF conformers within 3 kcal/mol of the lowest one).
Its two most distant heavy atoms are placed into the two pockets, and the molecule is rotated and shifted
to the best fit. Three things are then measured with van der Waals radii:

* **Pocket contact** – does each end touch the three C–H groups of its pocket?
* **End fit** – how much are the end atoms squeezed? Known guests: Cl 0.12 Å, Br 0.25 Å.
* **Walls** – how much does the rest of the molecule push into the walls of the slot? Known guests: ≤ 0.19 Å.

One chemical rule is applied (it can be switched off): the end atoms must have lone pairs (halogen, O, N, S)
to accept the C–H contacts of the pocket. In DFT, replacing Cl by H weakens binding by about 10 kcal/mol,
so the halogen–pocket contacts carry most of the binding.

**Limits.** The model looks at shape and at this one rule only. It does not compute energies, does not know
whether a compound is stable, and does not consider other crystal phases. Iodine-terminated guests need more
room than bromine; whether the host can provide it has not been checked yet, so they are reported as
*beyond tested range*.
""")
