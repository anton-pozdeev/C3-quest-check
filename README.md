# ⚗️ C3 Inclusion-Site Compatibility Screening

This tool evaluates the geometric and chemical compatibility of small molecules with the guest-binding site of the
PET cyclic trimer (C3) inclusion phase. Compatibility is assessed based on the structural criteria used in this
screening model.

**Web app:** https://c3-guest-check.streamlit.app

Open the link, type one or more SMILES (separated by spaces or commas) and click **Check**.
If the page says the app is asleep, click the button to wake it up and wait about a minute.

Examples: `Cl/C=C/Cl` (trans-1,2-dichloroethylene), `Cl/C=C\Cl` (cis), `ClCCCl` (1,2-dichloroethane), `FCCF`.

## Classification categories

* ✅ **COMPATIBLE** – both termini occupy their pockets within the reference limits.
* ❌ **NOT COMPATIBLE** – one or more essential geometric or chemical criteria are not satisfied.
* ➖ **BORDERLINE** – modest deviations from the reference geometric limits.
* ❔ **OUTSIDE THE MODEL** – salts, mixtures, charged species, elements other than H, C, N, O, F, S, Cl, Br, I.

## 🔬 Structural Model

The host geometry is derived from the DFT-optimized trans-1,2-dichloroethylene@C3 structure (PBE-D3, VASP) and remains fixed throughout the analysis.

The binding site consists of two pockets within neighboring C3 rings, connected by a narrow inter-ring region. Each pocket is defined by three aromatic C–H groups oriented toward one end of the guest molecule.

Candidate structures are generated from SMILES using RDKit and MMFF conformational sampling (3 kcal/mol energy window). The two most distant heavy atoms are positioned in opposing pockets, and molecular orientations are systematically sampled to minimize steric overlap.

## 📐 Compatibility Assessment

Geometric compatibility is evaluated using Bondi van der Waals radii and three structural descriptors:

* **Pocket contact** — proximity of each molecular terminus to its binding pocket.
* **Terminal compression** — steric compression experienced by the terminal atoms.
* **Framework overlap** — steric interference with the surrounding host structure.

In addition, each terminal atom must be able to accept C–H···X contacts from its pocket (halogen, O, N or S); carbon termini are classified as incompatible.

The screening criteria are informed by experimentally characterized C3 inclusion complexes and the observed inclusion behavior of guest molecules. These experimental references provide a basis for evaluating geometric compatibility and interpreting the resulting classifications.

## ⚠️ Scope and Limitations

The method assesses molecular compatibility with a specific, experimentally characterized binding geometry. It does not explicitly account for host-framework relaxation, interaction free energies, competing crystalline phases, or inclusion kinetics.

The screening results provide an assessment of structural compatibility, which may help identify promising candidates for inclusion-phase formation.

## Files

* `app.py` – the web app (Streamlit)
* `c3model_core.py` – the screening model
* `host_site.txt` – host atoms around the binding site, from the DFT structure
* `about.md` – method description shown on the web page
* `requirements.txt` – Python packages for the web app
* `C3_guest_check.ipynb` – the same tool as a Google Colab notebook (alternative to the web app):
  [open in Colab](https://colab.research.google.com/github/anton-pozdeev/C3-guest-check/blob/main/C3_guest_check.ipynb),
  then **Runtime → Run all**, type SMILES in the box that appears and click **Check**.

To run the web app on your own computer: `pip install -r requirements.txt`, then `streamlit run app.py`.
