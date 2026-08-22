"""EXP-045 prototype: parent 9a7638586033 setup + speed calibration."""
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

t0 = time.perf_counter()

# --- locate the catalogue row(s) for this parent ------------------------------
rows = E27.load_catalogue()
hits = []
for index, row in enumerate(rows):
    spec, HX, HZ = E27.parent_matrices(row)
    if E27.matrix_fingerprint(HX, HZ) == TARGET_FP:
        hits.append((index, E27.catalogue_label(row, index), row, HX, HZ))
print("rows sharing parent:", [(h[0], h[1]) for h in hits])
assert len(hits) >= 1
index, label, row, HX, HZ = hits[0]

cert = json.loads(CERT.read_text())
assert cert["fingerprint"] == TARGET_FP
assert cert["T_is_exact"] and cert["T"] == 4
assert cert["d_z_exact"] and cert["d_z_parent"] == 6
ell, m = int(row["ell"]), int(row["m"])
dim, n = ell * m, 2 * ell * m

A = HX[:, :dim].copy()
B = HX[:, dim:].copy()
assert (HX == np.hstack([A, B])).all()
assert (HZ == np.hstack([B.T, A.T])).all()
rx, rz = rank_np(HX), rank_np(HZ)
k_P = n - rx - rz
print(f"parent {label}: ell={ell} m={m} n={n} rank(HX)={rx} rank(HZ)={rz} k_P={k_P}")
assert k_P == 12 and k_P == cert["k_parent"]

L = nullspace_np(HX.T)
print("dim left-ker(HX) =", L.shape[0], " k_P/2 =", k_P // 2)
assert L.shape[0] * 2 == k_P

# --- module M from the certificate witnesses ---------------------------------
orbit_rows = []
for w in cert["witness_vectors"]:
    v = np.asarray(w, dtype=np.uint8).reshape(-1)
    assert v.shape[0] == n
    orb = translation_orbit(v, ell, m)
    assert (orb.sum(axis=1) == v.sum()).all()
    orbit_rows.append(orb)
M_full = np.vstack([HZ] + orbit_rows)
T_rebuilt = int(rank_np(M_full) - rz)
print("T rebuilt from persisted orbits =", T_rebuilt, "(cert", cert["T"], ")")
assert T_rebuilt == cert["T"] == 4

# --- quotient projection pi: GF(2)^n -> GF(2)^{n - rz}, ker = rowspace(HZ) ----
R, pivots = rref_np(HZ)
free = [j for j in range(n) if j not in set(pivots)]
pivrow = {p: i for i, p in enumerate(pivots)}
assert R.shape[0] == rz
print("quotient dim =", len(free))


def pi(rows: np.ndarray) -> np.ndarray:
    Z = (np.asarray(rows, dtype=np.uint8) & 1).copy()
    for p, i in pivrow.items():
        sel = Z[:, p].astype(bool)
        Z[sel] ^= R[i]
    return Z[:, free]


# check: pi kills rowspace(HZ)
assert not pi(HZ).any()

Mbar = pi(np.vstack(orbit_rows))


def pack(r):
    x = 0
    for j, b in enumerate(r):
        if b & 1:
            x |= 1 << j
    return x


def pivot_basis(rows):
    basis = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    return basis


mb_basis = list(pivot_basis(pack(r) for r in Mbar).values())
mb_basis.sort(reverse=True)
print("dim M_bar =", len(mb_basis), "(must equal T)")
assert len(mb_basis) == T_rebuilt

# --- validity kernel V --------------------------------------------------------
X = [monomial_matrix(ell, m, a, b) for a in range(ell) for b in range(m)]
ii, jj = np.triu_indices(dim, k=1)
cols = np.empty((2 * dim, ii.size), dtype=np.uint8)
for k in range(dim):
    P = gf2matmul(A, X[k].T)
    cols[k] = (P ^ P.T)[ii, jj]
for k in range(dim):
    P = gf2matmul(B, X[k].T)
    cols[dim + k] = (P ^ P.T)[ii, jj]
Nmat = cols.T
V = nullspace_np(Nmat)
dimV = int(V.shape[0])
print("dim V =", dimV, "(exp044 recorded 112)")

# catalogue-own perturbation must lie in V
cvec = np.zeros(2 * dim, dtype=np.uint8)
for a, b in E27.terms(row, "C_terms"):
    cvec[(a % ell) * m + (b % m)] = 1
for a, b in E27.terms(row, "D_terms"):
    cvec[dim + (a % ell) * m + (b % m)] = 1
assert not ((Nmat @ cvec) % 2).any(), "catalogue own perturbation not valid?!"
print("catalogue-own (C,D) in V: OK; weight =", int(cvec.sum()))

# --- per-basis-element packed Delta-bar images --------------------------------
def cd_from_coeff(cv):
    C = np.zeros((dim, dim), np.uint8)
    D = np.zeros((dim, dim), np.uint8)
    for k in range(2 * dim):
        if cv[k]:
            if k < dim:
                C ^= X[k]
            else:
                D ^= X[k - dim]
    return C, D


def classify_full(cv):
    C, D = cd_from_coeff(cv)
    CD = np.hstack([C, D])
    Mskew = gf2matmul(A, C.T) ^ gf2matmul(B, D.T)
    assert not ((Mskew ^ Mskew.T).any())
    Delta = gf2matmul(L, CD)
    bar = int(rank_np(np.vstack([HZ, Delta])) - rz)
    absorbed = int(rank_np(np.vstack([HZ, Delta] + orbit_rows))) == int(
        rank_np(np.vstack([HZ, Delta])))
    return bar, absorbed


bar_own, abs_own = classify_full(cvec)
print(f"catalogue-own: bar_Delta={bar_own} absorbed={abs_own}")

# packed images of V basis vectors
basis_imgs = []
for j in range(dimV):
    C, D = cd_from_coeff(V[j])
    Delta = gf2matmul(L, np.hstack([C, D]))
    proj = pi(Delta)
    basis_imgs.append(tuple(pack(r) for r in proj))

print("setup seconds:", round(time.perf_counter() - t0, 2))


def classify_packed(img):
    """rank of the packed row images + module containment, pivot-basis exact."""
    basis = {}
    for x in img:
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    bar = len(basis)
    absorbed = True
    for mrow in mb_basis:
        x = mrow
        while x:
            p = x.bit_length() - 1
            if p not in basis:
                absorbed = False
                break
            x ^= basis[p]
        if not absorbed:
            break
    return bar, absorbed


# cross-validate packed classifier on catalogue-own + first basis vectors
cv_img = tuple(pack(r) for r in pi(gf2matmul(L, np.hstack(cd_from_coeff(cvec)))))
assert classify_packed(cv_img) == (bar_own, abs_own)
ok = 0
for j in range(min(dimV, 30)):
    cvj = V[j]
    b1, a1 = classify_full(cvj)
    b2, a2 = classify_packed(basis_imgs[j])
    assert (b1, a1) == (b2, a2), (j, b1, a1, b2, a2)
    ok += 1
print(f"cross-validated {ok} basis elements + catalogue-own: OK")

# --- speed calibration ---------------------------------------------------------
import itertools
import math

def bench(wmax, limit=200_000):
    t = time.perf_counter()
    cnt = 0
    for w in range(1, wmax + 1):
        for combo in itertools.combinations(range(dimV), w):
            img = [0] * 6
            for j in combo:
                bj = basis_imgs[j]
                for t2 in range(6):
                    img[t2] ^= bj[t2]
            bar, absorbed = classify_packed(img)
            cnt += 1
            if cnt >= limit:
                return cnt, time.perf_counter() - t
    return cnt, time.perf_counter() - t


cnt, dt = bench(3, limit=200_000)
print(f"bench w<=3: {cnt} instances in {dt:.2f}s -> {dt/cnt*1e6:.2f} us/inst "
      f"(full C(112,<=3)={(sum(1 for _ in itertools.chain.from_iterable(itertools.combinations(range(dimV), w) for w in (1,2,3))))})")

# estimate for w<=4, w<=5
per = dt / cnt
for w in (4, 5):
    total = sum(math.comb(dimV, k) for k in range(1, w + 1))
    print(f"w_max={w}: {total:,} instances ~ {total*per:,.0f}s = {total*per/3600:.2f}h")
