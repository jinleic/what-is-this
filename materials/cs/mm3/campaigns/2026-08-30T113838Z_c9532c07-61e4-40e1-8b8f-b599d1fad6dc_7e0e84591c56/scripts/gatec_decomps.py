#!/usr/bin/env python3
"""Common factor-block loaders for the gate C orientation sweep.

Five decompositions, each normalized to the same interface:
  factors() -> (U, V, Wfac) where
     U, V  : lists of 23 ternary vectors in Z^9  (product r -> A/B coeffs)
     Wfac  : list of 23 ternary vectors in Z^9   (product r -> C-entry coeffs)
             (i.e. Wfac[r][k] = coefficient of product r in output C_k)
  meta    : dict with source paper, ring, anchor count, provenance.
All anchored: 729/729 Brent over Z checked in verify_anchors.py.
"""
import sys
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent          # .../cs/mm3/src
MM3 = HERE.parent                                # .../cs/mm3
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(MM3 / "scratch"))

R = 23


# ---------------------------------------------------------------- paper55 ---
def load_paper55():
    from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED
    U = [tuple(U_BLOCK_PRINTED[i][r] for i in range(9)) for r in range(R)]
    V = [tuple(V_BLOCK_PRINTED[i][r] for i in range(9)) for r in range(R)]
    W = [tuple(W_BLOCK_PRINTED[i][r] for i in range(9)) for r in range(R)]
    return U, V, W


# ------------------------------------------------------------- perminov58 ---
def load_perminov58():
    """Decode cr58_cn122 with Perminov's conventions (gate A's check_perminov),
    then transpose to the (U,V,Wfac) interface. Tensor == paper55 after CPERM."""
    import importlib
    ga = importlib.import_module("gate_a")
    (Up, Vp, Wp), _ = ga.slp_factors()   # 9x23 blocks printed form
    Wp_T = [list(x) for x in zip(*Wp)]   # 23x23? no: Wp is 9 rows x 23 cols
    U = [tuple(Up[i][r] for i in range(9)) for r in range(R)]
    V = [tuple(Vp[i][r] for i in range(9)) for r in range(R)]
    W = [tuple(Wp[i][r] for i in range(9)) for r in range(R)]
    return U, V, W


# ------------------------------------------------------------------ sun56 ---
def load_sun56():
    import sun56_verify as sv
    return [tuple(r) for r in sv.U], [tuple(r) for r in sv.V], [tuple(r) for r in sv.W]


# ------------------------------------------------------------------ mws59 ---
def load_mws59():
    """Table 3 of the MWS paper (layout lines 237-265 of mws59_layout.txt):
    three #-separated blocks. Each block has 23 rows of 9 entries.
    Row r of block 1 = A-coeffs of product r; block 2 = B-coeffs;
    block 3 row r = C-entry coefficients? — resolved by Brent check in
    verify_anchors.py (orientation validated there, not assumed)."""
    lines = (MM3 / "scratch" / "mws59_layout.txt").read_text().splitlines()
    # data lines are layout file lines 237-265 (1-based): 9 rows, '#', 9 rows,
    # '#', 9 rows. Rows = matrix entries, columns = 23 products (SW file format).
    data = lines[236:265]
    blocks, cur = [], []
    for L in data:
        s = L.strip()
        if s.startswith("#"):
            blocks.append(cur); cur = []
        elif s and not s.startswith(("A", "T")):
            try:
                cur.append([int(x) for x in s.split()])
            except ValueError:
                pass
    if cur:
        blocks.append(cur)
    b1, b2, b3 = blocks
    assert all(len(b) == 9 and all(len(row) == 23 for row in b) for b in blocks), \
        [ (len(b), len(b[0]) if b else 0) for b in blocks ]
    # b1[k][r] = coefficient of product r in A-entry k (rows = matrix entries,
    # columns = products, per the SW file format). So the factor target for
    # product r over A entries = column r of b1.
    U = [tuple(b1[k][r] for k in range(9)) for r in range(23)]
    V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]
    W = [tuple(b3[k][r] for k in range(9)) for r in range(23)]
    return U, V, W


# ------------------------------------------------------------- stapleton60 --
def load_stapleton60():
    """sd.build_factors() returns U,V (23 rows x 9) and W (9 rows x 23 cols,
    W[k][r] = coeff of product r in output C_k). Transpose W to Wfac rows."""
    import stapleton60_data as sd
    U, V, W, _ = sd.build_factors()
    Wfac = [tuple(W[k][r] for k in range(9)) for r in range(23)]
    return [tuple(x) for x in U], [tuple(x) for x in V], Wfac


LOADERS = {
    "paper55": load_paper55,
    "perminov58": load_perminov58,
    "sun56": load_sun56,
    "mws59": load_mws59,
    "stapleton60": load_stapleton60,
}

META = {
    "paper55": dict(paper="arXiv:2607.28676 (55)", ring="Z", anchor=55, split=(13, 14, 28)),
    "perminov58": dict(paper="arXiv:2512.21980 (Perminov) cr58_cn122 commit 98ba522",
                       ring="Z", anchor=58, split=None),
    "sun56": dict(paper="arXiv:2604.27645 (Sun)", ring="Z", anchor=56, split=(13, 13, 30)),
    "mws59": dict(paper="arXiv:2601.05272 (MWS) Table 3", ring="Z", anchor=59, split=(15, 15, 29)),
    "stapleton60": dict(paper="arXiv:2508.03857 (Stapleton) Appendix A", ring="Z",
                        anchor=60, split=(16, 16, 28)),
}
