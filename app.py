"""C3 inclusion-site compatibility screening — web app (Streamlit).

Run locally:   streamlit run app.py
"""
from pathlib import Path

import pandas as pd
import streamlit as st

import c3model_core as core

MAX_MOLECULES = 25
HERE = Path(__file__).parent

st.set_page_config(page_title="C3 Inclusion-Site Compatibility Screening", page_icon="⚗️", layout="centered")

if core._host_xyz is None:
    core.load_host((HERE / "host_site.txt").read_text())

MARK = {'✓': '✔️', '✗': '❌', '⚠': '⚠️', '~': '➖'}
BOX = {'fits': st.success, 'no': st.error, 'beyond': st.warning, 'border': st.info,
       'outside': st.info, 'unreadable': st.info}

LEGEND = """
* ✅ **COMPATIBLE** – both termini occupy their pockets within the reference limits.
* ❌ **NOT COMPATIBLE** – the reason is given for each molecule.
* ⚠️ **OUTSIDE CALIBRATED RANGE** – terminal atoms larger than bromine, the largest characterized terminus.
* ➖ **BORDERLINE** – between the reference limits and the incompatibility limits.
* ❔ **OUTSIDE THE MODEL** – salts, mixtures, charged species, elements other than H, C, N, O, F, S, Cl, Br, I.
"""


def show(r):
    head = f"**{core.ICON[r['verdict']]} {core.LABEL[r['verdict']]}**"
    if r.get('reasons') and r['verdict'] != 'fits':
        head += " — " + "; ".join(r['reasons'])
    BOX[r['verdict']](f"{head}  \n`{r['smiles']}`")
    if r['status'] == 'done':
        items = [f"Termini: **{r['ends'][0]} … {r['ends'][1]}**, {r['L']:.2f} Å apart "
                 f"(most extended conformer within 3 kcal/mol)"]
        items += [f"{MARK.get(m, m)} **{name}** — {text}" for m, name, text in core.checks(r)]
        st.markdown("\n".join(f"- {i}" for i in items))
    for w in r['why']:
        st.markdown(w)


st.title("⚗️ C3 Inclusion-Site Compatibility Screening")
st.markdown("This tool evaluates the geometric and chemical compatibility of small molecules with the "
            "guest-binding site of the PET cyclic trimer (C3) inclusion phase.  \n"
            "Compatibility is assessed based on the structural criteria used in this screening model.")

with st.form("check"):
    text = st.text_input("SMILES", placeholder="e.g.  Cl/C=C/Cl   ClCCCl   FCCF",
                         help="Several molecules: separate them with spaces or commas.")
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
                r = core.judge(core.analyse(smi))
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

with st.expander("Classification categories"):
    st.markdown(LEGEND)

st.divider()
st.markdown((HERE / "about.md").read_text())
