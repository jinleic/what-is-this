"""EXP-045 structure probe: total achievable Delta-bar row span A vs M_bar."""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

from qec_research.codes.bicycle import monomial_matrix  # noqa: E402
from qec_research.codes.pbb_survival import translation_orbit  # noqa: E402
from qec_research.gf2.linalg import matmul as gf2matmul  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

TARGET_FP = "9a7638586033f4c7ae22623a4b53173fa72492f1b1212e94a47e9d8aa2467a4c"
CERT = ROOT / "results" / "partial_runs" / "exp039" / "parent_9a7638586033f4c7.json"

rows = E27.load_catalogue()
row = None
for index, r in enumerate(rows):
    spec, HX, HZ = E27.parent_matrices(r)
    if E27.matrix_fingerprint(HX, HZ) == TARGET_FP:
        row = r
        break
assert row is not None
ell, m = int(row["ell"]), int(row["m"])
dim, n = ell * m, 2 * ell * m
HX, HZ = E27.parent_matrices(row)[1:]
A = HX[:, :dim].copy()
B = HX[:, dim:].copy()
rz = rank_np(HZ)

cert = json.loads(CERT.read_text())
orbit_rows = [translation_orbit(np.asarray(w, dtype=np.uint8).reshape(-1), ell, m)
              for w in cert["witness_vectors"]]

R, pivots = rref_np(HZ)
free = [j for j in range(n) if j not in set(pivots)]
pivrow = {p: i for i, p in enumerate(pivots)}


def pi(mats):
    Z = (np.asarray(mats, dtype=np.uint8) & 1).copy()
    for p, i in pivrow.items():
        sel = Z[:, p].astype(bool)
        Z[sel] ^= R[i]
    return Z[:, free]


def pack(r):
    x = 0
    for j, b in enumerate(r):
        if b & 1:
            x |= 1 << j
    return x


def pivot_basis(ints):
    basis = {}
    for x in ints:
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    return basis


Mbar_bits = [pack(r) for r in pi(np.vstack(orbit_rows))]
mb = pivot_basis(Mbar_bits)
assert len(mb) == cert["T"] == 4

L = nullspace_np(HX.T)
X = [monomial_matrix(ell, m, a, b) for a in range(ell) for b in range(m)]
ii, jj = np.triu_indices(dim, k=1)
cols = np.empty((2 * dim, ii.size), dtype=np.uint8)
for k in range(dim):
    P = gf2matmul(A, X[k].T)
    cols[k] = (P ^ P.T)[ii, jj]
for k in range(dim):
    P = gf2matmul(B, X[k].T)
    cols[dim + k] = (P ^ P.T)[ii, jj]
V = nullspace_np(cols.T)
dimV = int(V.shape[0])
assert dimV == 112

# total achievable Delta-bar row span A_space over ALL of V (span of all rows
# of F(basis_j)); row i of F(v) = pi(L . CD(v))_i
all_rows = []
basis_imgs = []
for j in range(dimV):
    cv = V[j]
    C = np.zeros((dim, dim), np.uint8)
    D = np.zeros((dim, dim), np.uint8)
    for k in range(2 * dim):
        if cv[k]:
            if k < dim:
                C ^= X[k]
            else:
                D ^= X[k - dim]
    proj = pi(gf2matmul(L, np.hstack([C, D])))
    bits = [pack(r) for r in proj]
    basis_imgs.append(bits)
    all_rows.extend(bits)

A_basis = pivot_basis(all_rows)
print("dim A (total achievable Delta-bar row span) =", len(A_basis))

# is M_bar contained in A?
missing = [p for p, v in mb.items() if p not in A_basis]
contained = True
for p, mv in mb.items():
    x = mv
    while x:
        q = x.bit_length() - 1
        if q not in A_basis:
            contained = False
            break
        x ^= A_basis[q]
    if not contained:
        break
print("M_bar contained in A:", contained)

# per-row images: E_i = span of row i over all v in V (same as A for each i?)
for i in range(6):
    Ei = pivot_basis(basis_imgs[j][i] for j in range(dimV))
    print(f"dim E_{i} =", len(Ei))

# distribution over random v in V (sanity: bar values, absorbed rate)
rng = np.random.default_rng(42)
import time
t0 = time.perf_counter()
hist = {}
abs_cnt = 0
for _ in range(20000):
    sel = rng.integers(0, 2, dimV) @ V % 2
    img = [0] * 6
    for k in np.flatnonzero(sel):
        b = basis_imgs[k]
        for t in range(6):
            img[t] ^= b[t]
    bas = pivot_basis(img)
    bar = len(bas)
    hist[bar] = hist.get(bar, 0) + 1
    ab = True
    for mv in mb.values():
        x = mv
        while x:
            q = x.bit_length() - 1
            if q not in bas:
                ab = False
                break
            x ^= bas[q]
        if not ab:
            break
    abs_cnt += ab
print("random-v bar histogram:", dict(sorted(hist.items())),
      "absorbed:", abs_cnt, "of 20000, ", round(time.perf_counter() - t0, 1), "s")
