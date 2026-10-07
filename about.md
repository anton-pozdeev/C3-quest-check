### 🔬 Structural Model

The host geometry is derived from the DFT-optimized trans-1,2-dichloroethylene@C3 structure (PBE-D3, VASP) and remains fixed throughout the analysis.

The binding site consists of two pockets within neighboring C3 rings, connected by a narrow inter-ring region. Each pocket is defined by three aromatic C–H groups oriented toward one end of the guest molecule.

Candidate structures are generated from SMILES using RDKit and MMFF conformational sampling (3 kcal/mol energy window). The two most distant heavy atoms are positioned in opposing pockets, and molecular orientations are systematically sampled to minimize steric overlap.

### 📐 Compatibility Assessment

Geometric compatibility is evaluated using Bondi van der Waals radii and three structural descriptors:

* **Pocket contact** — proximity of each molecular terminus to its binding pocket.
* **Terminal compression** — steric compression experienced by the terminal atoms.
* **Framework overlap** — steric interference with the surrounding host structure.

In addition, each terminal atom must be able to accept C–H···X contacts from its pocket (halogen, O, N or S); carbon termini are classified as incompatible.

The screening criteria are informed by the geometric characteristics of experimentally characterized C3 inclusion complexes. Minor deviations from the reference limits are classified as borderline, whereas substantial steric compression or overlap results in an incompatible classification.

### ✅ Evaluation

The screening procedure was additionally evaluated against compounds experimentally reported not to form the corresponding inclusion phase (dichloromethane, methanol, n-hexane, 1,1-dichloroethane, 1,3-dichloropropane and cis-1,2-dichloroethylene). These compounds were classified as incompatible without adjustment of the screening criteria. The classifications remained unchanged when the chemical criterion was not applied.

### ⚠️ Scope and Limitations

The method assesses molecular compatibility with a specific, experimentally characterized binding geometry. It does not explicitly account for host-framework relaxation, interaction free energies, competing crystalline phases, or inclusion kinetics.

The screening results provide an assessment of structural compatibility, which may help identify promising candidates for inclusion-phase formation.
