# ---------------------------------------------------------------------------
# C3 guest check: will a small molecule form the inclusion phase with the PET
# cyclic trimer (C3)?
#
# Model (geometry of the host + one chemical rule):
#   * In the inclusion phase each end of the guest sits in a pocket inside one
#     C3 ring, lined by three aromatic C-H groups. The two pockets belong to two
#     neighbouring C3 rings; between them is a flat slot.
#   * The host is taken from the DFT-relaxed trans-DCE@C3 structure and kept
#     frozen. The guest is built from SMILES (RDKit, MMFF conformers within
#     3 kcal/mol), its two most distant heavy atoms are placed into the two
#     pockets, and it is rotated/shifted to the best fit.
#   * Measured: does each end touch the three C-H of its pocket; how much the
#     ends are squeezed; how much the rest of the molecule pushes into the walls.
#   * Reference values come from the known guests in this host:
#     Cl ends squeezed 0.12 A, Br ends 0.25 A; walls overlap <= 0.19 A.
# ---------------------------------------------------------------------------
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

VDW = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'S': 1.80,
       'Cl': 1.75, 'Br': 1.85, 'I': 1.98}                    # Bondi radii, A
ALLOWED = set(VDW)
LONE_PAIR_ENDS = {'F', 'Cl', 'Br', 'I', 'O', 'N', 'S'}
NAMES = {'F': 'fluorine', 'Cl': 'chlorine', 'Br': 'bromine', 'I': 'iodine',
         'O': 'oxygen', 'N': 'nitrogen', 'S': 'sulfur', 'C': 'carbon'}

# limits (A), set by the guests confirmed in this host
POCKET_TOUCH = 0.05     # end must be within this of contact with all three pocket C-H
SQUEEZE_CL, SQUEEZE_BR = 0.12, 0.25        # largest values among the confirmed guests
SQUEEZE_BORDER, SQUEEZE_MAX = 0.255, 0.30
WALL_OK, WALL_MAX = 0.19, 0.35
MAX_HEAVY = 10

_host_el, _host_xyz, _host_r, _pockets = None, None, None, None


def load_host(text):
    """text: lines 'El x y z' in the site frame (pockets on the z axis at z = +/-2.1406)."""
    global _host_el, _host_xyz, _host_r, _pockets
    el, xyz = [], []
    for line in text.strip().splitlines():
        s, x, y, z = line.split()
        el.append(s); xyz.append((float(x), float(y), float(z)))
    _host_el, _host_xyz = np.array(el), np.array(xyz)
    _host_r = np.array([VDW[s] for s in el])
    H = np.where(_host_el == 'H')[0]
    _pockets = []
    for z in (-2.1406, 2.1406):
        d = np.linalg.norm(_host_xyz[H] - np.array([0, 0, z]), axis=1)
        _pockets.append(_host_xyz[H[np.argsort(d)[:3]]])


def _rot_z_to(v):
    v = v / np.linalg.norm(v); z = np.array([0., 0., 1.])
    c = float(np.dot(v, z)); w = np.cross(v, z)
    if np.linalg.norm(w) < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1., -1., -1.])
    K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
    return np.eye(3) + K + K @ K / (1 + c)


def _rz(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _rx(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _ry(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _screen(smi):
    """Return (mol, problem_message_or_None, verdict_if_rejected)."""
    smi = smi.strip()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, "Could not read this SMILES. Please check the spelling.", "unreadable"
    if len(Chem.GetMolFrags(mol)) > 1:
        return None, "More than one molecule (salt or mixture). The model handles one neutral molecule.", "outside"
    if any(a.GetFormalCharge() != 0 for a in mol.GetAtoms()):
        return None, "Charged species are outside the model (it covers neutral molecules).", "outside"
    bad = sorted({a.GetSymbol() for a in mol.GetAtoms()} - ALLOWED)
    if bad:
        return None, f"Contains {', '.join(bad)}; the model covers H, C, N, O, F, S, Cl, Br, I.", "outside"
    nheavy = mol.GetNumHeavyAtoms()
    if nheavy < 2:
        return None, "A single heavy atom cannot bridge the two pockets.", "small"
    if nheavy > MAX_HEAVY:
        return None, (f"{nheavy} heavy atoms. The cavity holds a rod of about four heavy atoms "
                      f"(X-C-C-X), so this molecule is far too large."), "large"
    return mol, None, None


def analyse(smi, nconf=60, ewin=3.0, seed=7):
    mol, problem, kind = _screen(smi)
    if problem:
        return dict(smiles=smi, status=kind, message=problem)
    m = Chem.AddHs(mol)
    cids = list(AllChem.EmbedMultipleConfs(m, nconf, randomSeed=seed))
    if not cids:
        return dict(smiles=smi, status='outside', message="Could not build a 3D structure for this molecule.")
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    E = np.array([e for _, e in res])
    S = [a.GetSymbol() for a in m.GetAtoms()]
    R = np.array([VDW[s] for s in S])
    heavy = [i for i, s in enumerate(S) if s != 'H']

    # candidate placements: each accessible conformer, both end-to-pocket assignments
    cands = []
    for cid, e in zip(cids, E):
        if e - E.min() > ewin:
            continue
        X = m.GetConformer(cid).GetPositions()
        L, i, j = max((np.linalg.norm(X[a] - X[b]), a, b)
                      for k, a in enumerate(heavy) for b in heavy[k + 1:])
        for a, b in ((i, j), (j, i)):
            Y = X - (X[a] + X[b]) / 2
            cands.append((Y @ _rot_z_to(Y[b] - Y[a]).T, a, b, L, e - E.min()))
    Lmax = max(c[3] for c in cands)
    near = np.linalg.norm(_host_xyz, axis=1) < Lmax / 2 + 7.0
    HP, HR = _host_xyz[near], _host_r[near]
    side_of = lambda a, b: np.array([k for k in range(len(S)) if k not in (a, b)], dtype=int)

    H2 = (HP ** 2).sum(1)

    def score(poses, a, b, chunk=150):
        """poses: (P, n, 3) -> worst overlap per pose, plus per-atom max overlaps.
        Distances via |p|^2 + |h|^2 - 2 p.h in chunks, to keep memory low."""
        P, n, _ = poses.shape
        ov = np.empty((P, n))
        for s in range(0, P, chunk):
            Q = poses[s:s + chunk].reshape(-1, 3)
            d2 = (Q ** 2).sum(1)[:, None] + H2[None, :] - 2.0 * Q @ HP.T
            D = np.sqrt(np.maximum(d2, 0.0)).reshape(-1, n, len(HP))
            ov[s:s + chunk] = (R[None, :, None] + HR[None, None, :] - D).max(2)
        sd = side_of(a, b)
        wall = ov[:, sd].max(1) if len(sd) else np.full(len(poses), -9.0)
        return np.maximum(wall, np.maximum(ov[:, a], ov[:, b])), ov, wall

    def poses_for(Y, phis, shifts, tilts):
        out, par = [], []
        for tx, ty in tilts:
            T = _rx(tx) @ _ry(ty)
            for phi in phis:
                Z = Y @ (T @ _rz(phi)).T
                for sh in shifts:
                    out.append(Z + np.array([0, 0, sh])); par.append((phi, sh, tx, ty))
        return np.array(out), par

    # 1) coarse search
    coarse = []
    for ci, (Y, a, b, L, dE) in enumerate(cands):
        P, par = poses_for(Y, np.radians(np.arange(0, 360, 10)), np.arange(-0.3, 0.31, 0.1), [(0, 0)])
        sc, _, _ = score(P, a, b)
        for k in np.argsort(sc)[:3]:
            coarse.append((sc[k], ci, par[k]))
    coarse.sort(key=lambda t: t[0])

    # 2) refine the best few placements (finer rotation, shift, small tilt of the axis)
    tilts = [(np.radians(x), np.radians(y)) for x in (-6, -3, 0, 3, 6) for y in (-6, -3, 0, 3, 6)]
    best = None
    for s0, ci, (phi0, sh0, _, _) in coarse[:5]:
        Y, a, b, L, dE = cands[ci]
        P, par = poses_for(Y, phi0 + np.radians(np.arange(-6, 6.1, 1.5)),
                           sh0 + np.arange(-0.1, 0.101, 0.025), tilts)
        sc, ov, wall = score(P, a, b)
        k = int(np.argmin(sc))
        if best is None or sc[k] < best['score']:
            W = P[k]
            gapA = (np.linalg.norm(_pockets[0] - W[a], axis=1) - R[a] - 1.20).max()
            gapB = (np.linalg.norm(_pockets[1] - W[b], axis=1) - R[b] - 1.20).max()
            best = dict(score=float(sc[k]), L=L, ends=(S[a], S[b]), squeeze=(ov[k, a], ov[k, b]),
                        gap=(gapA, gapB), wall=float(wall[k]), dE=dE)
    best.update(smiles=smi, status='done')
    return best


def judge(r, carbon_end_rule=True):
    """Adds r['verdict'] in {'fits','no','beyond','border', ...} and r['lines'] (checks) and r['why']."""
    if r['status'] != 'done':
        r['verdict'] = {'large': 'no', 'small': 'no'}.get(r['status'], r['status'])
        r['reasons'] = {'large': ['too large'], 'small': ['too small']}.get(r['status'], [])
        r['lines'], r['why'] = [], [r['message']]
        return r
    e1, e2 = r['ends']; s1, s2 = r['squeeze']; g1, g2 = r['gap']; wall = r['wall']
    fails, warns, borders, lines, why = [], [], [], [], []

    # 1. chemistry of the ends
    carbon = [e for e in (e1, e2) if e not in LONE_PAIR_ENDS]
    if carbon:
        mark = '✗' if carbon_end_rule else '–'
        lines.append(f"{mark} End atoms         {e1} and {e2}: a carbon end has no lone pairs for the pocket C–H"
                     + ("" if carbon_end_rule else "  (rule switched off)"))
        if carbon_end_rule:
            fails.append('carbon end')
            why.append("An end atom is carbon. Each pocket is lined by three aromatic C–H groups that hold the end "
                       "through C–H···X contacts; a carbon end has no lone pairs to accept them "
                       "(in DFT, replacing Cl by H weakens binding by about 10 kcal/mol).")
    else:
        lines.append(f"✓ End atoms         {e1} and {e2}: can accept C–H contacts")

    # 2. pocket contact
    far = [(e, g) for e, g in ((e1, g1), (e2, g2)) if g > POCKET_TOUCH]
    if far:
        if len(far) == 2 and far[0][0] == far[1][0]:
            txt = f"both {far[0][0]} ends: gap {far[0][1]:.2f} / {far[1][1]:.2f} Å"
        else:
            txt = ", ".join(f"{e} end: gap {g:.2f} Å" for e, g in far)
        lines.append(f"✗ Pocket contact    {txt}  (must touch)")
        fails.append('end does not reach pocket')
        why.append("At least one end cannot touch the C–H groups of its pocket while the other end sits in its own "
                   "pocket: the molecule is too short, or its end atom too small.")
    else:
        lines.append(f"✓ Pocket contact    both ends touch the C–H groups of their pockets")

    # 3. squeeze of the ends
    smax = max(s1, s2)
    ref = f"(confirmed: Cl {SQUEEZE_CL:.2f}, Br {SQUEEZE_BR:.2f})"
    if smax > SQUEEZE_MAX:
        lines.append(f"⚠ End fit           squeezed {s1:.2f} / {s2:.2f} Å  {ref}")
        warns.append('ends larger than bromine')
        why.append(f"The ends are larger than bromine (squeezed {smax:.2f} Å vs {SQUEEZE_BR:.2f} Å for "
                   f"1,2-dibromoethane). The host is known to make room for bromine; whether it can open "
                   f"further has not been checked yet.")
    elif smax > SQUEEZE_BORDER:
        lines.append(f"~ End fit           squeezed {s1:.2f} / {s2:.2f} Å  {ref}")
        borders.append('ends at the bromine limit')
    else:
        lines.append(f"✓ End fit           squeezed {max(s1,0):.2f} / {max(s2,0):.2f} Å  {ref}")

    # 4. walls of the slot
    lim = f"(confirmed guests ≤ {WALL_OK:.2f})"
    if wall > WALL_MAX:
        lines.append(f"✗ Walls             overlap {wall:.2f} Å  {lim}")
        fails.append('hits the walls')
        why.append("Part of the molecule pushes into the walls of the flat slot between the two rings — "
                   "usually a side group, a bend, or a molecule that is too long or too bulky.")
    elif wall > WALL_OK:
        lines.append(f"~ Walls             overlap {wall:.2f} Å  {lim}")
        borders.append('close to the walls')
    else:
        lines.append(f"✓ Walls             overlap {max(wall,0):.2f} Å  {lim}")

    if fails:
        r['verdict'] = 'no'; r['reasons'] = fails
    elif warns:
        r['verdict'] = 'beyond'; r['reasons'] = warns
    elif borders:
        r['verdict'] = 'border'; r['reasons'] = borders
        why.append("Close to the limits set by the confirmed guests; the model cannot decide this one reliably.")
    else:
        r['verdict'] = 'fits'; r['reasons'] = []
        why.append("Both ends sit in the C–H pockets of two neighbouring C3 rings, and nothing pushes into "
                   "the walls between them.")
    r['lines'], r['why'] = lines, why
    return r


ICON = {'fits': '✅', 'no': '❌', 'beyond': '⚠️', 'border': '➖', 'outside': '❔', 'unreadable': '❔'}
LABEL = {'fits': 'FITS', 'no': 'DOES NOT FIT', 'beyond': 'BEYOND TESTED RANGE',
         'border': 'BORDERLINE', 'outside': 'OUTSIDE THE MODEL', 'unreadable': 'COULD NOT READ'}


def report(r):
    bar = '═' * 72
    out = [bar, f"{r['smiles']}"]
    head = f"{ICON[r['verdict']]} {LABEL[r['verdict']]}"
    if r.get('reasons') and r['verdict'] != 'fits':
        head += ' — ' + '; '.join(r['reasons'])
    out.append(f"RESULT:  {head}")
    if r['status'] == 'done':
        out.append(f"\n  Ends: {r['ends'][0]} … {r['ends'][1]}, {r['L']:.2f} Å apart (most extended shape within 3 kcal/mol)")
        out += ['  ' + l for l in r['lines']]
    out.append('')
    for w in r['why']:
        out += _wrap(w, 70, '  ')
    return '\n'.join(out)


def _wrap(text, width, indent):
    words, lines, cur = text.split(), [], indent
    for w in words:
        if len(cur) + len(w) + 1 > width + len(indent):
            lines.append(cur.rstrip()); cur = indent
        cur += w + ' '
    lines.append(cur.rstrip())
    return lines


def run(text, carbon_end_rule=True):
    items = [s for s in text.replace(',', ' ').replace(';', ' ').split() if s]
    if not items:
        print("Type one or more SMILES in the box above (separated by spaces or commas), then run again.")
        return
    results = []
    for smi in items:
        r = judge(analyse(smi), carbon_end_rule)
        results.append(r)
        print(report(r), flush=True)
    if len(results) > 1:
        print('═' * 72 + '\nSUMMARY')
        for r in results:
            extra = f"  ({'; '.join(r['reasons'])})" if r.get('reasons') and r['verdict'] != 'fits' else ''
            print(f"  {ICON[r['verdict']]} {LABEL[r['verdict']]:<20s} {r['smiles']}{extra}")


def checks(r):
    """Structured form of r['lines']: list of (mark, check name, text)."""
    return [(l[0], l[2:20].strip(), l[20:].strip()) for l in r.get('lines', [])]


def split_input(text):
    return [s for s in text.replace(',', ' ').replace(';', ' ').split() if s]
