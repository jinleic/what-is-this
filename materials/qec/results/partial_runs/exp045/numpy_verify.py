"""Independent raw-numpy verification: dim bar_Delta(v) <= 2 for all v in V."""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
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
    _, HX, HZ = E27.parent_matrices(r)
    if E27.matrix_fingerprint(HX, HZ) == TARGET_FP:
        row = r
        break
ell, m = int(row["ell"]), int(row["m"])
dim, n = ell * m, 2 * ell * m
A_ = HX[:, :dim].copy()
B_ = HX[:, dim:].copy()
rz = rank_np(HZ)

cert = json.loads(CERT.read_text())
orbit = np.vstack(
    [translation_orbit(np.asarray(w, dtype=np.uint8).reshape(-1), ell, m)
     for w in cert["witness_vectors"]])

R, pivots = rref_np(HZ)
free = [j for j in range(n) if j not in set(pivots)]
pivrow = {p: i for i, p in enumerate(pivots)}
L = nullspace_np(HX.T)

X = [monomial_matrix(ell, m, a, b) for a in range(ell) for b in range(m)]
ii, jj = np.triu_indices(dim, k=1)
cols = np.empty((2 * dim, ii.size), dtype=np.uint8)
for k in range(dim):
    P = gf2matmul(A_, X[k].T)
    cols[k] = (P ^ P.T)[ii, jj]
for k in range(dim):
    P = gf2matmul(B_, X[k].T)
    cols[dim + k] = (P ^ P.T)[ii, jj]
V = nullspace_np(cols.T)
assert V.shape[0] == 112


def cd(cv):
    C = np.zeros((dim, dim), np.uint8)
    D = np.zeros((dim, dim), np.uint8)
    for k in range(2 * dim):
        if cv[k]:
            if k < dim:
                C ^= X[k]
            else:
                D ^= X[k - dim]
    return np.hstack([C, D])


def pi(mats):
    Z = (np.asarray(mats, dtype=np.uint8) & 1).copy()
    for p, i in pivrow.items():
        sel = Z[:, p].astype(bool)
        Z[sel] ^= R[i]
    return Z[:, free]


# (1) numpy rank of the stacked projected basis-image rows == 2 ?
stacked = []
imgrows = []
for j in range(112):
    proj = pi(gf2matmul(L, cd(V[j])))
    stacked.append(proj)
    imgrows.append(proj)
stacked = np.vstack(stacked)
print("numpy rank of all 672 projected basis rows =", rank_np(stacked))
assert rank_np(stacked) == 2

# per-row numpy ranks: 0,0,1,1,1,1
for i in range(6):
    rows_i = np.vstack([imgrows[j][i] for j in range(112)])
    print(f"numpy rank of row-{i} images =", rank_np(rows_i))

# (2) dense random full-space check: bar_Delta <= 2 and never absorbed
rng = np.random.default_rng(0xE045)
t0 = time.perf_counter()
maxbar = 0
absorbed = 0
N_SAMPLES = 400
for _ in range(N_SAMPLES):
    sel = rng.integers(0, 2, 112).astype(np.uint8)
    cv = (sel @ V) % 2
    CD = cd(cv)
    Delta = gf2matmul(L, CD)
    bar = int(rank_np(np.vstack([HZ, Delta])) - rz)
    maxbar = max(maxbar, bar)
    assert bar <= 2, f"bar_Delta {bar} > 2 observed!"
    ab = int(rank_np(np.vstack([HZ, Delta, orbit]))) == int(rank_np(np.vstack([HZ, Delta])))
    absorbed += ab
    assert not ab
print(f"dense random full-numpy: {N_SAMPLES} draws, max bar_Delta={maxbar}, absorbed={absorbed}, "
      f"{round(time.perf_counter()-t0,1)}s")

# (3) M_bar vs A via raw quotient check on stacked matrix: rank(pi(orbit)) = 4, rank(pi(stacked+orbit)) = 6
proj_orbit = pi(orbit)
rm = rank_np(proj_orbit)
rt = rank_np(np.vstack([stacked, proj_orbit]))
print("rank(pi(orbit)) =", rm, "; rank(pi(stacked + orbit)) =", rt)
assert rm == 4 and rt == 6, "M_bar adds 4 independent dimensions beyond A"
print("EXACT: M_bar NOT contained in A (dim 2); absorption impossible over all of V.")
