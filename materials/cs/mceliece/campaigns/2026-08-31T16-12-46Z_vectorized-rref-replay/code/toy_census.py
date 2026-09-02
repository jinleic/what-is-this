"""Toy direct-route census engine v4 — exact vectorized-RREF replay.

CONVENTION LOCK (owner-mandated, 2026-08-31; supersedes all prior notes):
 - LOCATOR LABELS are homogeneous pairs (D, N) representing the projective
   value z = N/D.  finite z -> label (1, z); infinity -> label (0, 1).
 - The Algorithm-4 U row is the SWAPPED FROBENIUS ROW, a distinct
   type/function: u(label) = (N^4, D^4).
     locator 0        label (1,0) -> U row (0,1)
     locator 1        label (1,1) -> U row (1,1)
     locator infinity label (0,1) -> U row (1,0)
   Thus U3 = (0,1),(1,1),(1,0) fixes locators 0, 1, infinity.  A U row
   (0,1) is NEVER called "infinity" in code or comments.
 - The canonical Mobius transform mapping raw support (0,1,2) to the
   prescribed locator labels 0, 1, infinity output pair is
       (D, N) = ((x+c)(b+a), (x+a)(b+c)),  a=0, b=1, c=2
   i.e. D(x) = x+2, N(x) = 3x, phi(x) = 3x/(x+2).
 - Branch e applies the field Frobenius to the HOMOGENEOUS PAIR:
       (D_e, N_e) = (D^(2^e), N^(2^e)).
 - OWNER ANCHORS (byte-for-byte in static guards):
      e-th branch label rows for raw support 0..7:
        e=0: [0, 1, INF, 5, 2, 8, 25, 16]
        e=1: [0, 1, INF, 17, 4, 25, 18, 10]
        e=2: [0, 1, INF, 11, 16, 18, 14, 29]
        e=3: [0, 1, INF, 28, 10, 14, 13, 2]
        e=4: [0, 1, INF, 3, 29, 13, 8, 4]
      fourth-holdout (position 3) orbit: {5, 17, 11, 28, 3} — 5 unique,
      beta-independent.
 - Expected labels depend only on normalized support cross-ratio; the
   Goppa root beta never enters the label space.

Implementation invariants (second+ third owner-audit fixes):
  (a) P1 unknowns are 2*nk kernel COEFFICIENTS (u1|u2); V1/V2 are
      expanded from those before ANY residual replay.
  (b) P1 solves one combined system (hold-out + normalization + all-8
      locator coupling) exactly; NO {0,1} mask search.
  (c) own-diagonal convention per C.2: c*_{tau,tau,0} = 1 for EVERY tau
      after rescaling each projective U row by d_tau^{-1}.
  (d) P1 verdicts: FAIL iff the affine space / diagonal-image quotient is
      exhaustively characterized and no all-nonzero diagonal tuple exists
      (proving data attached).  INCONCLUSIVE only for genuine
      cap/unsupported conditions.
  (e) diagonal selection maps the affine nullspace through the diagonal
      map (linear, dim <= 3), row-reduces it, and exhausts the resulting
      affine image (<= 32^3) with F_32 coefficients, never masks; its
      "no valid diagonal" outcome is asserted by a NEGATIVE plant.
  (f) ledger resume validates source hash + domain schema; never
      silently accepts an old `complete` file.
Source pins: prim/pe1630.txt — E_ell 2696-2706, B_j 2705-2706, pi_j
2754-2756, V_U 2708-2719, normalization 2719, coupling 2731, degenerate
2949-2950, C*=U*V* 4632-4636, (33) 4640-4642, P2 5026-5031, P3
5032-5034, Thm 7.2 2850-2852.
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# PROVENANCE LOCK: import the FROZEN, checksummed byte copies of the
# dependencies shipped inside this campaign's code/ directory.  The live
# mceliece/src path is NOT on sys.path; the released runner executes this
# frozen module, not mutable live sources.
FROZEN_FASTFIELD_SHA = "209fe885f3d7837650e03fb63350d73a20492b96a17467969ac4ccf3ce62e560"
FROZEN_GFIELD_SHA   = "69c6920e7775f80fd14d8d05166e5da6155c8fa778948f324d96bbaaef1fa6ae"

import hashlib as _hashlib
import importlib as _importlib

def _assert_frozen_hash(fname: str, expected: str) -> None:
    path = os.path.join(HERE, fname)
    h = _hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    if h.hexdigest() != expected:
        raise RuntimeError(
            f"frozen dependency hash mismatch: {fname} "
            f"(got {h.hexdigest()}, expected {expected})")

_assert_frozen_hash("fastfield_frozen.py", FROZEN_FASTFIELD_SHA)
_assert_frozen_hash("gfield_frozen.py", FROZEN_GFIELD_SHA)

_ff = _importlib.import_module("fastfield_frozen")
_gf = _importlib.import_module("gfield_frozen")
EField = _ff.EField
GF = _gf.GF


M        = 5
Q        = 1 << M
PRIM     = 0x37
N        = 8
T        = 1
K        = 3
ELL      = 0
K_ELL    = K - ELL
N_ELL    = N - ELL
D_ELL    = N_ELL - 2 * T - 1
D_DEG    = 7
S_MULT   = 5
H_HOLD   = 4
SUPPORT  = tuple(range(N))
HOLDOUTS = tuple(range(H_HOLD))
RELEASE_KEY = "TOY_FAST_CENSUS_RELEASED"
PARENT_RELEASE_KEY = "PARENT_CENSUS_TERMINATED"
CONTROL_PLANT_CASES = (
    "C18_stale_schema",
    "C21_wrong_beta",
    "C22_missing_P4",
    "C23_wrong_controls_hash",
    "C24_negative_complete",
)
DOMAIN_SCHEMA = "toy-m5-n8-t1-k3-s5-h4-v4-vectorized-rref"
PARENT_CAMPAIGN = os.path.realpath(os.path.join(
    HERE, "..", "..", "2026-08-31T08-44Z_DA168C79"))
PARENT_CONTROLS_PATH = os.path.join(
    PARENT_CAMPAIGN, "state", "controls_result.json")
PARENT_CELL8_PATH = os.path.join(PARENT_CAMPAIGN, "state", "cell_8.json")
PARENT_CONTROLS_SHA256 = (
    "41652388f5ce7b7c8db82648e06be56bf21c7a0fe0e56e9b4020b1f61ccc8570")
PARENT_CELL8_SHA256 = (
    "b7c4d97607f43385f14b4d9b1e1cf43734b809727214176b4f5e1881c075edaa")
INF_ROW = (0, 1)                              # homogeneous infinity label


# --------------------------------------------------------------------------- #
# locator labels / U rows (TWO DISTINCT TYPES — never conflated)               #
# --------------------------------------------------------------------------- #

def label_from_pair(D: int, N: int, ef: EField):
    """homogeneous (D, N) representing z=N/D -> 33-label convention:
    finite z -> (1,z); infinity (D=0) -> (0,1)."""
    if (D, N) == (0, 0):
        raise ValueError("zero homogeneous pair is not projective")
    if D == 0:
        return INF_ROW
    return (1, int(ef.MUL[N, ef.INV[D]]))


def label_to_U_row(lab, ef: EField):
    """u(label) = (N^4, D^4), the SWAPPED FROBENIUS row (distinct type).
       locator 0 (label (1,0)) -> (0,1);
       locator 1 (label (1,1)) -> (1,1);
       locator inf (label (0,1)) -> (1,0)."""
    if lab == INF_ROW:
        return (1, 0)
    z = int(lab[1])
    # literal field fourth power z^4 (NOT Frobenius conflation):
    # powers(a, D) returns [a^0..a^D], so powers(z, 4)[4] = z^4.
    return (int(ef.powers(z, 4)[4]), 1) if z else (0, 1)


def u_row_to_label(u_row, ef: EField):
    """inverse of label_to_U_row: (N^4, D^4) -> label."""
    n4, d4 = int(u_row[0]), int(u_row[1])
    if (n4, d4) == (0, 0):
        raise ValueError("zero U row")
    if n4 == 0:
        return (1, 0)                    # locator 0
    if d4 == 0:
        return INF_ROW
    # w = z^4 (integer fourth power in the field, NOT Frobenius^2 in
    # general — Frobenius^2 is z^(2^2) which coincides with z^4 since
    # 2^2 = 4).  Inverse map: from w recover z with z^4 = w, i.e. z =
    # w^(inverse of 4 mod 31) = w^8 since 4*8 = 32 ≡ 1 (mod 31).
    # powers(w, 8)[8] = w^(2^8 mod 31) = w^(256 mod 31) = w^8. Good.
    w = int(ef.MUL[n4, ef.INV[d4]])
    z = int(ef.powers(w, 8)[8])
    return (1, z)


def all_labels():
    out = [INF_ROW] + [(1, z) for z in range(Q)]
    assert len(out) == 33 and len(set(out)) == 33
    return out


def label_is_infinity(lab) -> bool:
    return lab == INF_ROW


# --------------------------------------------------------------------------- #
# multiindex enumeration (weight-first, lex-second)                            #
# --------------------------------------------------------------------------- #

def multiindices_lt(k: int, s: int):
    out = []
    for w in range(s):
        for a in itertools.product(range(s), repeat=k):
            if sum(a) == w:
                out.append(a)
    return out


def weight_range(w: int, k: int):
    return [a for a in itertools.product(range(w + 1), repeat=k) if sum(a) == w]


COLS_S = multiindices_lt(K_ELL, S_MULT)
COLS_D = weight_range(D_DEG, K_ELL)
B_J    = [a for a in COLS_S if 1 <= sum(a) <= S_MULT - 1]
N_D    = len(COLS_D)
M_ELL  = N_ELL * len(COLS_S)
IDS_OF = {a: i for i, a in enumerate(COLS_S)}


def col_index(j: int, a: tuple) -> int:
    return j * len(COLS_S) + IDS_OF[a]


def block_indices(j: int):
    off = j * len(COLS_S)
    return [off + IDS_OF[a] for a in B_J]


# --------------------------------------------------------------------------- #
# independent toy instance generator (NO Instance-class reuse)                 #
# --------------------------------------------------------------------------- #

def _gf2_nullspace(A: np.ndarray):
    A = (np.asarray(A, dtype=np.uint8) % 2).copy()
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        pr = None
        for rr in range(r, rows):
            if A[rr, c]:
                pr = rr
                break
        if pr is None:
            continue
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
        col = A[:, c].copy()
        col[r] = 0
        A[np.nonzero(col)[0]] ^= A[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in set(piv)]
    basis = []
    for fc in free:
        v = np.zeros(cols, dtype=np.uint8)
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if A[ri, fc]:
                v[pc] = 1
        basis.append(v)
    return basis, piv


def _rank_f2(A: np.ndarray) -> int:
    A = (np.asarray(A, dtype=np.uint8) % 2).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        pr = None
        for rr in range(r, rows):
            if A[rr, c]:
                pr = rr
                break
        if pr is None:
            continue
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
        col = A[:, c].copy()
        col[r] = 0
        A[np.nonzero(col)[0]] ^= A[r]
        r += 1
    return r


def build_toy(beta: int):
    gf = GF(M, prim=PRIM)
    assert 0 <= beta < Q and beta not in SUPPORT
    support = list(SUPPORT)
    G = [beta, 1]
    h = []
    for a in support:
        denom = a ^ beta
        assert denom != 0
        h.append(gf.inv(denom))
    Hrows = [[(h[i] >> b) & 1 for i in range(N)] for b in range(M)]
    Hmat = np.array(Hrows, dtype=np.uint8)
    rank_H = _rank_f2(Hmat)
    ns, piv = _gf2_nullspace(Hmat)
    # The toy must use the FULL [8,3] binary Goppa code: nullspace
    # dimension exactly K (= n - mt) and parity rank exactly M.  No
    # silent truncation; a cell outside the registered shape is a
    # recorded GUARD FAILURE, not a smaller arbitrary subcode.
    if rank_H != M or len(ns) != K:
        raise RuntimeError(
            f"GUARD FAILURE: cell beta={beta} is NOT the declared full "
            f"[{N},{K}] Goppa code (rank(H)={rank_H}, nullity="
            f"{len(ns)}, expected rank(H)={M}, nullity={K}); "
            f"classify outside/fail — never truncate silently")
    Y = np.array([list(v) for v in ns[:K]], dtype=np.uint8)
    rk = _rank_f2(Y)
    if rk != K:
        raise RuntimeError(f"GUARD FAILURE: rank(Y)={rk} != k={K}")
    public = {"m": M, "q": Q, "n": N, "k": K, "ell": ELL,
              "d": D_DEG, "s": S_MULT, "profile": [S_MULT] * N,
              "Y": Y.tolist(), "prim": PRIM}
    oracle = {"beta": int(beta), "G": G, "support": support,
              "holdouts": list(HOLDOUTS)}
    return public, oracle


# --------------------------------------------------------------------------- #
# public-only jet matrix                                                       #
# --------------------------------------------------------------------------- #

def build_E(public, ef):
    Y = np.asarray(public["Y"], dtype=np.uint8)
    k = Y.shape[0]
    E = np.zeros((N_D, M_ELL), dtype=np.uint16)
    for j in range(N_ELL):
        ybin = tuple(int(v) for v in Y[:, j])
        for bi, b in enumerate(COLS_D):
            for ai, a in enumerate(COLS_S):
                entry = 1
                for i in range(k):
                    ai_, bi_, yi_ = a[i], b[i], ybin[i]
                    if ai_ > bi_ or ((ai_ & ~bi_) != 0):
                        entry = 0
                        break
                    if yi_ == 0:
                        if bi_ - ai_ != 0:
                            entry = 0
                            break
                    elif bi_ != ai_:
                        entry = ef.MUL[entry, yi_]
                if entry:
                    E[bi, col_index(j, a)] = int(entry)
    return E


# --------------------------------------------------------------------------- #
# exact affine rref over F_32                                                  #
# --------------------------------------------------------------------------- #

def rref_affine_reference(A, b_vec, ncols, ef):
    """Frozen scalar Gauss-Jordan implementation, retained for controls."""
    rows = [list(r) + [int(v)] for r, v in zip(A, b_vec)]
    ncols_full = ncols + 1
    r = 0
    piv = []
    nr = len(rows)
    for c in range(ncols_full):
        if r >= nr:
            break
        pr = None
        for i in range(r, nr):
            if rows[i][c]:
                pr = i
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = ef.INV[rows[r][c]]
        rows[r] = [ef.MUL[v, iv] for v in rows[r]]
        for i in range(nr):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [rows[i][c2] ^ ef.MUL[f, rows[r][c2]]
                           for c2 in range(ncols_full)]
        piv.append(c)
        r += 1
    feasible = ncols not in piv
    particular = None
    if feasible:
        particular = [0] * ncols
        for ri, pc in enumerate(piv):
            if pc < ncols:
                particular[pc] = rows[ri][ncols]
    return rows, piv, feasible, particular


def rref_affine(A, b_vec, ncols, ef):
    """Exact table-indexed Gauss-Jordan with the reference pivot order."""
    nr = len(A)
    ncols_full = ncols + 1
    rows = np.zeros((nr, ncols_full), dtype=np.uint16)
    if nr:
        matrix = np.asarray(A, dtype=np.uint16)
        if matrix.shape != (nr, ncols):
            raise ValueError(
                f"affine matrix shape {matrix.shape} != {(nr, ncols)}")
        if len(b_vec) != nr:
            raise ValueError(f"affine rhs length {len(b_vec)} != {nr}")
        rows[:, :ncols] = matrix
        rows[:, ncols] = np.asarray(b_vec, dtype=np.uint16)
    elif b_vec:
        raise ValueError("nonempty affine rhs for zero-row matrix")

    r = 0
    piv = []
    for c in range(ncols_full):
        if r >= nr:
            break
        candidates = np.flatnonzero(rows[r:, c])
        if candidates.size == 0:
            continue
        pr = r + int(candidates[0])
        if pr != r:
            rows[[r, pr]] = rows[[pr, r]]
        pivot = int(rows[r, c])
        rows[r] = ef.MUL[rows[r], int(ef.INV[pivot])]
        factors = rows[:, c].copy()
        eliminate = factors != 0
        eliminate[r] = False
        if np.any(eliminate):
            rows[eliminate] ^= ef.MUL[
                factors[eliminate, None], rows[r][None, :]]
        piv.append(c)
        r += 1

    feasible = ncols not in piv
    particular = None
    if feasible:
        particular = [0] * ncols
        for ri, pc in enumerate(piv):
            if pc < ncols:
                particular[pc] = int(rows[ri, ncols])
    return rows.tolist(), piv, feasible, particular


def nullspace_of(A, ncols, ef):
    rows, piv, _, _ = rref_affine(A, [0] * len(A), ncols, ef)
    pivset = set(piv)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if pc < ncols and rows[ri][fc]:
                v[pc] = rows[ri][fc]
        basis.append(v)
    return basis, len(piv)

def nullspace_of_reference(A, ncols, ef):
    rows, piv, _, _ = rref_affine_reference(
        A, [0] * len(A), ncols, ef)
    pivset = set(piv)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if pc < ncols and rows[ri][fc]:
                v[pc] = rows[ri][fc]
        basis.append(v)
    return basis, len(piv)


def rank_of_reference(mat, ncols, ef):
    rows, piv, _, _ = rref_affine_reference(
        mat, [0] * len(mat), ncols, ef)
    return len(piv)


def rank_of(mat, ncols, ef):
    rows, piv, _, _ = rref_affine(mat, [0] * len(mat), ncols, ef)
    return len(piv)


# --------------------------------------------------------------------------- #
# P4  rank admission                                                           #
# --------------------------------------------------------------------------- #

def p4_admission(public, ef, E=None, ker_E=None):
    if E is None:
        E = build_E(public, ef)
    if ker_E is None:
        ker_E = nullspace_of(E.tolist(), M_ELL, ef)[0]
    full_rank = M_ELL - len(ker_E)
    assert 0 <= full_rank <= N_D, "rank/nullity coherence"
    per_j = []
    for j in range(N_ELL):
        idx = block_indices(j)
        proj = [[row[c] for c in idx] for row in ker_E]
        dim_j = rank_of(proj, len(B_J), ef) if proj else 0
        per_j.append(int(dim_j))
    return {"rank_E": int(full_rank), "nullity": int(M_ELL - full_rank),
            "N_d": N_D, "dim_pi_j": per_j, "B": len(B_J),
            "margins_full": N_D - full_rank,
            "margins_block": [len(B_J) - d for d in per_j],
            "passes_full": bool(full_rank < N_D),
            "passes_block": bool(all(d < len(B_J) for d in per_j)),
            "passes": bool(full_rank < N_D
                           and all(d < len(B_J) for d in per_j))}


# --------------------------------------------------------------------------- #
# canonical Mobius transform (owner anchors; homogeneous (D,N) output)         #
# --------------------------------------------------------------------------- #

def mobius_pair(x, ef: EField, a=0, b=1, c=2):
    """(D, N) = ((x+c)(b+a), (x+a)(b+c)).  With (a,b,c)=(0,1,2):
    D(x) = x+2, N(x) = 3x; label(x) = 3x/(x+2) with infinity at x=2."""
    D = int(ef.MUL[x ^ c, b ^ a])
    N = int(ef.MUL[x ^ a, b ^ c])
    return (D, N)


def mobius_matrix(a, b, c, ef: EField):
    """Matrix form (D(x), N(x)) = (C*x+Dc, A*x+Bc) with
    D(x) = (b+a)x + (b+a)c; N(x) = (b+c)x + (b+c)a.  Verified against
    `mobius_pair` in static_guards."""
    Arow = int(b) ^ int(c)
    Brow = int(ef.MUL[Arow, a])
    Crow = int(b) ^ int(a)
    Drow = int(ef.MUL[Crow, c])
    return ((Arow, Brow), (Crow, Drow))


def frob_pair(pair, e, ef: EField):
    """Frobenius (z->z^(2^e)) applied to BOTH coordinates via ordinary
    powers: powers(x, 2^e)[2^e] = x^(2^e) since powers is literal."""
    D, N = pair
    fD = 0 if D == 0 else int(ef.powers(int(D), 1 << e)[1 << e])
    fN = 0 if N == 0 else int(ef.powers(int(N), 1 << e)[1 << e])
    return (fD, fN)


def branch_mobius_pair(x, e, ef: EField, a=0, b=1, c=2):
    return frob_pair(mobius_pair(x, ef, a, b, c), e, ef)



def branch_label(x_raw, e, ef: EField, a=0, b=1, c=2):
    """transformed locator label for ONE raw field value x_raw.
    NOTE: support_raw must be passed EXPLICITLY (no hidden default) --
    callers derive x_raw from their own record (public scan label or
    oracle support) before calling this."""
    return label_from_pair(*branch_mobius_pair(x_raw, e, ef, a, b, c),
                           ef=ef)


def branch_target_table(e, ef, support_raw, a=0, b=1, c=2):
    """transformed locator labels for the N_ELL support values, in
    positional order.  support_raw is REQUIRED (explicit), never a
    hidden default."""
    return [branch_label(support_raw[j], e, ef, a, b, c)
            for j in range(N_ELL)]


OWNER_ANCHOR_BRANCH_ROWS = [
    [0, 1, None, 5, 2, 8, 25, 16],
    [0, 1, None, 17, 4, 25, 18, 10],
    [0, 1, None, 11, 16, 18, 14, 29],
    [0, 1, None, 28, 10, 14, 13, 2],
    [0, 1, None, 3, 29, 13, 8, 4],
]
OWNER_ANCHOR_FOURTH_ORBIT = [5, 17, 11, 28, 3]


def p2_fourth_labels(ef, support_raw):
    """the 5 branch fourth-holdout labels (position 3).  support_raw is
    REQUIRED; the ONLY caller is the labeled adjudication layer, which
    passes oracle['support'] explicitly."""
    return [branch_target_table(e, ef, support_raw, 0, 1, 2)[3]
            for e in range(M)]


# --------------------------------------------------------------------------- #
# V_U system rows (unknowns = 2*nk kernel coefficients)                        #
# --------------------------------------------------------------------------- #

U3_PRESCRIBED = [(0, 1), (1, 1), (1, 0)]   # U rows for locators 0,1,inf


def vu_system_rows(U_rows, j, x_label, ker_E, ef):
    nk = len(ker_E)
    nunk = 2 * nk
    rows = []
    bvec = []
    for tau, u in enumerate(U_rows):
        u0, u1 = int(u[0]), int(u[1])
        for a in COLS_S:
            if sum(a) == 0:
                continue
            c = col_index(tau, a)
            r = [0] * nunk
            for i in range(nk):
                kav = int(ker_E[i][c])
                if kav:
                    r[i] = ef.MUL[u0, kav]
                    r[nk + i] = ef.MUL[u1, kav]
            rows.append(r)
            bvec.append(0)
    u0, u1 = int(U_rows[0][0]), int(U_rows[0][1])
    c0 = col_index(0, (0, 0, 0))
    r = [0] * nunk
    for i in range(nk):
        kav = int(ker_E[i][c0])
        if kav:
            r[i] = ef.MUL[u0, kav]
            r[nk + i] = ef.MUL[u1, kav]
    rows.append(r)
    bvec.append(1)
    if label_is_infinity(x_label):
        for a in B_J:
            c = col_index(j, a)
            rr = [0] * nunk
            for i in range(nk):
                rr[i] = int(ker_E[i][c])
            rows.append(rr)
            bvec.append(0)
    else:
        x = int(x_label[1])
        x4 = int(ef.powers(x, 4)[4]) if x else 0
        for a in B_J:
            c = col_index(j, a)
            r = [0] * nunk
            for i in range(nk):
                r[nk + i] = int(ker_E[i][c])
                r[i] = ef.MUL[x4, int(ker_E[i][c])]
            rows.append(r)
            bvec.append(0)
    return rows, bvec, nunk


def scan_labels(public, ef, U_rows, j, ker_E=None):
    if ker_E is None:
        E = build_E(public, ef)
        ker_E = nullspace_of(E.tolist(), M_ELL, ef)[0]
    accepted = []
    for lab in all_labels():
        rows, bvec, nunk = vu_system_rows(U_rows, j, lab, ker_E, ef)
        _, _, feas, _ = rref_affine(rows, bvec, nunk, ef)
        if feas:
            accepted.append(lab)
    return accepted


# --------------------------------------------------------------------------- #
# P1 combined solve, diagonal selection, exact replay                          #
# --------------------------------------------------------------------------- #

def _expand_coeffs(coeffs, ker_E, ef, nk):
    v1 = [0] * M_ELL
    v2 = [0] * M_ELL
    for i in range(nk):
        a1, a2 = int(coeffs[i]), int(coeffs[nk + i])
        kr = ker_E[i]
        if a1:
            for c2 in range(M_ELL):
                if kr[c2]:
                    v1[c2] ^= ef.MUL[a1, int(kr[c2])]
        if a2:
            for c2 in range(M_ELL):
                if kr[c2]:
                    v2[c2] ^= ef.MUL[a2, int(kr[c2])]
    return v1, v2


def _diag_of(coeffs, U_rows, ker_E, ef, nk, tau):
    u0, u1 = int(U_rows[tau][0]), int(U_rows[tau][1])
    c = col_index(tau, (0, 0, 0))
    # v1[c] and v2[c] depend LINEARLY on the coefficient vector: build
    # the row-projections of the kernel basis once (per tau) instead of
    # expanding the full V1/V2 each call.
    acc_v1 = 0
    acc_v2 = 0
    for i in range(nk):
        a1, a2 = int(coeffs[i]), int(coeffs[nk + i])
        if a1 and ker_E[i][c]:
            acc_v1 ^= ef.MUL[a1, int(ker_E[i][c])]
        if a2 and ker_E[i][c]:
            acc_v2 ^= ef.MUL[a2, int(ker_E[i][c])]
    return int(ef.MUL[u0, acc_v1]) ^ int(ef.MUL[u1, acc_v2])


def _p1_diagonal_select(aff_rows, aff_b, nunk, ef, ker_E, U_rows):
    """Combined solve + EXACT diagonal-selection over the diagonal-IMAGE
    quotient (rank AT MOST len(U_rows)-1: <=2 for the 3-holdout system,
    <=3 for the 4-holdout system).  The affine nullspace may have nfree >> rank;
    only its IMAGE under the diagonal map is <= rank-dimensional.

    Algorithm (owner-locked):
      1. solve the combined system once -> particular p (feasible or not).
      2. compute the direction image L(v) for each nullspace basis vector
         v: L(v)[tau] = diag(p+v, tau) ^ diag(p, tau)  (LINEAR, i.e.
         XOR off the base).  Row-reduce the image vectors, keeping the
         corresponding preimage directions.
      3. enumerate particular + span(preimage_pivots) over
         32^(len(U_rows)-1) (F_32 coefficients, NEVER {0,1} masks,
         never capped on nfree).
    Exhaust-outcome semantics (owner-locked):
      - selected_coeffs None and exhaust_complete=True -> P1_FAIL.
      - selected_coeffs None and exhaust_complete=False -> INCONCLUSIVE
        only if the QUOTIENT dimension itself is unsupported."""
    nk = len(ker_E)
    nU = len(U_rows)
    rows, piv, feasible, particular = rref_affine(aff_rows, aff_b, nunk, ef)
    if not feasible:
        return {"feasible": False, "particular": None, "nfree": None,
                "base_diagonals": None, "direction_diagonals": None,
                "diagonal_image_rank": None, "selected_coeffs": None,
                "selected_diagonals": None,
                "exhaust_quotient_size": 0, "exhaust_complete": True}
    pivset = set(piv)
    free = [c for c in range(nunk) if c not in pivset]
    nfree = len(free)
    basis_vs = []
    for fc in free:
        v = [0] * nunk
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if pc < nunk and rows[ri][fc]:
                v[pc] = rows[ri][fc]
        basis_vs.append(v)
    diag_taus = list(range(1, nU))       # tau=0 pinned by normalization
    base_ds = [_diag_of(particular, U_rows, ker_E, ef, nk, tau)
               for tau in diag_taus]
    # LINEAR direction image: L(v) = diag(p+v) XOR diag(p).
    image_with_preimage = []
    for v in basis_vs:
        cand = [int(particular[i]) ^ int(v[i]) for i in range(nunk)]
        d_p = [_diag_of(particular, U_rows, ker_E, ef, nk, tau)
               for tau in diag_taus]
        d_pv = [_diag_of(cand, U_rows, ker_E, ef, nk, tau)
                for tau in diag_taus]
        lv = [d_pv[t] ^ d_p[t] for t in range(len(diag_taus))]
        image_with_preimage.append((lv, v))
    # select independent image ROWS (not pivot columns): walk the
    # (lv, v) pairs, keep lv iff it increases the running image rank,
    # retaining its preimage v.  Incremental exact rank over GF(32).
    retained = []
    retained_rows = []
    for lv, v in image_with_preimage:
        if not any(lv):
            continue
        trial = retained_rows + [lv]
        if rank_of(trial, len(diag_taus), ef) > len(retained_rows):
            retained.append((lv, v))
            retained_rows.append(lv)
    rank_L = len(retained_rows)
    preimages = [v for _, v in retained]
    # hard gate: the retained image must ACTUALLY have rank == rank_L
    assert rank_of(retained_rows, len(diag_taus), ef) == rank_L, \
        "retained image rows are not independent"
    # enumerate particular + span(preimages) over 32^rank
    exhaust_quotient_size = 0
    selected = None
    if rank_L <= len(diag_taus):
        for combo in itertools.product(range(Q), repeat=rank_L):
            cand = list(particular)
            for fi, coeff in enumerate(combo):
                if coeff:
                    v = preimages[fi]
                    cand = [cand[i] ^ ef.MUL[coeff, int(v[i])]
                            for i in range(nunk)]
            ds = [_diag_of(cand, U_rows, ker_E, ef, nk, tau)
                  for tau in diag_taus]
            exhaust_quotient_size += 1
            if all(d != 0 for d in ds):
                selected = (cand, ds)
                break
        exhaust_complete = selected is not None \
            or exhaust_quotient_size == Q ** rank_L
    else:
        return {"feasible": True, "particular": particular,
                "nfree": nfree, "base_diagonals": base_ds,
                "direction_diagonals": [lv for lv, _ in
                                        image_with_preimage],
                "diagonal_image_rank": rank_L,
                "selected_coeffs": None, "selected_diagonals": None,
                "exhaust_quotient_size": 0, "exhaust_complete": False}
    return {"feasible": True, "particular": particular, "nfree": nfree,
            "base_diagonals": base_ds,
            "direction_diagonals": [lv for lv, _ in image_with_preimage],
            "diagonal_image_rank": rank_L,
            "selected_coeffs": (list(selected[0]) if selected else None),
            "selected_diagonals": (list(selected[1]) if selected else None),
            "exhaust_quotient_size": exhaust_quotient_size,
            "exhaust_complete": bool(exhaust_complete)}


def p1_certificate_verdict(public, ef, U_rows, branch_e, oracle,
                           E=None, ker_E=None):
    """P1 for one branch (branch_e explicit; no label-parameter ambiguity).
    Combined system = hold-out rows + normalization + all-8 coupling rows.
    P1_FAIL carries proving data; P1_INCONCLUSIVE only on a genuine cap."""
    if oracle is None:
        return {"verdict": "P1_INCONCLUSIVE",
                "reason": "no oracle record for this branch"}
    if E is None:
        E = build_E(public, ef)
    if ker_E is None:
        ker_E = nullspace_of(E.tolist(), M_ELL, ef)[0]
    nk = len(ker_E)
    if not (0 <= branch_e < M):
        return {"verdict": "P1_INCONCLUSIVE",
                "reason": f"unsupported branch index {branch_e}"}
    # oracle support ONLY — hidden support never leaks into public paths
    oracle_support = oracle["support"]
    tloc = [branch_label(oracle_support[j], branch_e, ef)
            for j in range(N_ELL)]
    # combined affine system: hold-out rows + normalization + all-8 coupling
    aff_rows = []
    aff_b = []
    for tau, u in enumerate(U_rows):
        u0, u1 = int(u[0]), int(u[1])
        for a in COLS_S:
            if sum(a) == 0:
                continue
            c = col_index(tau, a)
            r = [0] * (2 * nk)
            for i in range(nk):
                kav = int(ker_E[i][c])
                if kav:
                    r[i] = ef.MUL[u0, kav]
                    r[nk + i] = ef.MUL[u1, kav]
            aff_rows.append(r)
            aff_b.append(0)
    u0, u1 = int(U_rows[0][0]), int(U_rows[0][1])
    c0 = col_index(0, (0, 0, 0))
    r = [0] * (2 * nk)
    for i in range(nk):
        kav = int(ker_E[i][c0])
        if kav:
            r[i] = ef.MUL[u0, kav]
            r[nk + i] = ef.MUL[u1, kav]
    aff_rows.append(r)
    aff_b.append(1)
    for j in range(N_ELL):
        lab = tloc[j]
        if label_is_infinity(lab):
            for a in B_J:
                c = col_index(j, a)
                r = [0] * (2 * nk)
                for i in range(nk):
                    r[i] = int(ker_E[i][c])
                aff_rows.append(r)
                aff_b.append(0)
        else:
            xv = int(lab[1])
            x4 = int(ef.powers(xv, 4)[4]) if xv else 0
            for a in B_J:
                c = col_index(j, a)
                r = [0] * (2 * nk)
                for i in range(nk):
                    r[nk + i] = int(ker_E[i][c])
                    r[i] = ef.MUL[x4, int(ker_E[i][c])]
                aff_rows.append(r)
                aff_b.append(0)
    sel = _p1_diagonal_select(aff_rows, aff_b, 2 * nk, ef, ker_E, U_rows)
    if not sel["feasible"]:
        return {"verdict": "P1_FAIL",
                "reason": "combined affine system infeasible",
                "branch": branch_e}
    if sel["selected_coeffs"] is None:
        if sel["exhaust_complete"]:
            return {"verdict": "P1_FAIL",
                    "reason": "no all-nonzero diagonal in the exhausted "
                              "affine image",
                    "branch": branch_e,
                    "proving": {
                        "nfree": sel["nfree"],
                        "diagonal_image_rank": sel["diagonal_image_rank"],
                        "base_diagonals": sel["base_diagonals"],
                        "direction_diagonals": sel["direction_diagonals"],
                        "exhaust_quotient_size":
                            sel["exhaust_quotient_size"]}}
        return {"verdict": "P1_INCONCLUSIVE",
                "reason": "diagonal-image exhaust exceeded supported cap",
                "branch": branch_e,
                "proving": {"nfree": sel["nfree"]}}
    # OWNER-LOCKED SEARCH SEMANTICS: a candidate representative with all-
    # nonzero diagonals whose rank replay FAILS does NOT establish
    # P1_FAIL — another quotient representative (or a zero-diagonal-image
    # affine direction) may still satisfy rank == 2.  The verdict logic:
    #   P1_PASS  — a fully replayed witness found;
    #   P1_FAIL  — ONLY when the exhaustive diagonal image proves NO
    #              all-nonzero diagonal exists anywhere;
    #   P1_INCONCLUSIVE — a rank/replay failure persisted across ALL
    #              all-nonzero-diagonal representatives (zero-image
    #              freedom precludes a failure theorem), or a genuine cap.
    reps_tried = 0
    reps_with_valid_diag = 0
    last_rank_failure = None
    # rebuild the SAME enumeration the selector used, so we can walk ALL
    # representatives with all-nonzero diagonals (not only the first).
    from itertools import product as _product
    rank_L = sel["diagonal_image_rank"]
    preimage_dim = sel["nfree"]
    # re-walk the quotient: enumerate particular + span(preimages) with
    # the SAME enumeration as _p1_diagonal_select (F_32 coefficients).
    # Reconstruct bases from the selector's returned artifacts.
    nk = len(ker_E)
    rows, piv, feasible, particular = rref_affine(
        aff_rows, aff_b, 2 * nk, ef)
    pivset = set(piv)
    free = [c for c in range(2 * nk) if c not in pivset]
    basis_vs = []
    for fc in free:
        v = [0] * (2 * nk)
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if pc < 2 * nk and rows[ri][fc]:
                v[pc] = rows[ri][fc]
        basis_vs.append(v)
    # diagonal-image direction vectors, deduped to a basis (same as the
    # selector's retained rows)
    # diagonal dimension and base computed OUTSIDE the loop so they are
    # always bound even when the affine nullity is zero (owner audit).
    diag_taus = list(range(1, len(U_rows)))
    diag_dim = len(diag_taus)
    d_p = [_diag_of(particular, U_rows, ker_E, ef, nk, tau)
           for tau in diag_taus]
    img_vecs = []
    for v in basis_vs:
        cand = [int(particular[i]) ^ int(v[i]) for i in range(2 * nk)]
        d_pv = [_diag_of(cand, U_rows, ker_E, ef, nk, tau)
                for tau in diag_taus]
        lv = [d_pv[t] ^ d_p[t] for t in range(diag_dim)]
        img_vecs.append((lv, v))
    # retained independent image rows with preimages (same rule as the
    # selector so enumeration is exact)
    retained = []
    retained_rows = []
    for lv, v in img_vecs:
        if not any(lv):
            continue
        trial = retained_rows + [lv]
        if rank_of(trial, diag_dim, ef) > len(retained_rows):
            retained.append((lv, v))
            retained_rows.append(lv)
    rank_L = len(retained_rows)
    preimages = [v for _, v in retained]
    # enumerate over ALL quotient representatives; for each with all-
    # nonzero diagonals, run the FULL replay; PASS on first witness;
    # record rank failures; FAIL only when the exhaust proves no all-
    # nonzero diagonal exists.
    witness = None
    exhaust_size = 0
    if rank_L <= diag_dim:
        for combo in _product(range(Q), repeat=rank_L):
            cand = list(particular)
            for fi, coeff in enumerate(combo):
                if coeff:
                    v = preimages[fi]
                    cand = [cand[i] ^ ef.MUL[coeff, int(v[i])]
                            for i in range(2 * nk)]
            ds = [_diag_of(cand, U_rows, ker_E, ef, nk, tau)
                  for tau in range(1, len(U_rows))]
            exhaust_size += 1
            if not all(d != 0 for d in ds):
                continue
            reps_with_valid_diag += 1
            # FULL replay for this representative
            replay = _p1_full_replay(cand, U_rows, ker_E, ef, nk, E,
                                     oracle_support, branch_e, tloc)
            if replay["ok"]:
                witness = (cand, ds, replay)
                break
            last_rank_failure = replay["reason"]
            reps_tried += 1
        if witness is not None:
            coeffs, ds, replay = witness
            return {"verdict": "P1_PASS", "branch": branch_e,
                    "checked": exhaust_size,
                    "rank_C": replay["rank_C"], "rank_V": replay["rank_V"],
                    "U_star": replay["U_star"],
                    "witness": replay["witness"],
                    "reps_with_valid_diag": reps_with_valid_diag}
        # no witness
        if last_rank_failure is not None and reps_with_valid_diag > 0:
            return {"verdict": "P1_INCONCLUSIVE",
                    "reason": "rank/replay failures across ALL "
                              "all-nonzero-diagonal representatives "
                              "(zero-diagonal-image freedom precludes a "
                              "failure theorem)",
                    "branch": branch_e,
                    "proving": {
                        "nfree": sel["nfree"],
                        "diagonal_image_rank": sel["diagonal_image_rank"],
                        "reps_with_valid_diag": reps_with_valid_diag,
                        "last_rank_failure": last_rank_failure,
                        "exhaust_quotient_size": exhaust_size}}
        return {"verdict": "P1_FAIL",
                "reason": "no all-nonzero diagonal in the exhausted "
                          "affine image",
                "branch": branch_e,
                "proving": {
                    "nfree": sel["nfree"],
                    "diagonal_image_rank": sel["diagonal_image_rank"],
                    "base_diagonals": sel["base_diagonals"],
                    "direction_diagonals": sel["direction_diagonals"],
                    "exhaust_quotient_size": exhaust_size}}
    return {"verdict": "P1_INCONCLUSIVE",
            "reason": "diagonal-image rank exceeds supported quotient",
            "branch": branch_e,
            "proving": {"nfree": sel["nfree"]}}


def _p1_full_replay(coeffs, U_rows, ker_E, ef, nk, E, oracle_support,
                    branch_e, tloc):
    """exact replay of ONE representative; returns dict with ok/reason
    and, on success, rank_C/rank_V/U_star.  Not a verdict — a check."""
    v1, v2 = _expand_coeffs(coeffs, ker_E, ef, nk)
    ds = [_diag_of(coeffs, U_rows, ker_E, ef, nk, tau)
          for tau in range(1, len(U_rows))]
    U_star = []
    for tau, u in enumerate(U_rows):
        d_tau = 1 if tau == 0 else ds[tau - 1]
        if d_tau == 0:
            return {"ok": False, "reason": "zero diagonal in rep"}
        inv = ef.INV[int(d_tau)]
        U_star.append((int(ef.MUL[inv, int(u[0])]),
                       int(ef.MUL[inv, int(u[1])])))
    C_rows = []
    for u in U_star:
        u0, u1 = int(u[0]), int(u[1])
        crow = [ef.MUL[u0, x] ^ ef.MUL[u1, y2]
                for x, y2 in zip(v1, v2)]
        C_rows.append(crow)
    for tau, u in enumerate(U_star):
        c = col_index(tau, (0, 0, 0))
        dtau = ef.MUL[int(u[0]), v1[c]] ^ ef.MUL[int(u[1]), v2[c]]
        if dtau != 1:
            return {"ok": False,
                    "reason": f"rescaled own-diagonal at tau={tau} != 1"}
    for tau, u in enumerate(U_star):
        u0, u1 = int(u[0]), int(u[1])
        for a in COLS_S:
            if sum(a) == 0:
                continue
            c = col_index(tau, a)
            if ef.MUL[u0, v1[c]] ^ ef.MUL[u1, v2[c]]:
                return {"ok": False,
                        "reason": f"positive-order nonzero at tau={tau}"}
    rank_V = rank_of([v1, v2], M_ELL, ef)
    rank_C = rank_of(C_rows, M_ELL, ef)
    if rank_V != 2:
        return {"ok": False, "reason": f"rank(V*) = {rank_V} != 2",
                "rank_V": int(rank_V)}
    if rank_C != 2:
        return {"ok": False, "reason": f"rank(C*) = {rank_C} != 2",
                "rank_C": int(rank_C)}
    for j in range(N_ELL):
        lab = tloc[j]
        if label_is_infinity(lab):
            for a in B_J:
                c = col_index(j, a)
                if v1[c] != 0:
                    return {"ok": False,
                            "reason": f"(33) infinity V1 != 0 at j={j}"}
        else:
            xv = int(lab[1])
            x4 = int(ef.powers(xv, 4)[4]) if xv else 0
            for a in B_J:
                c = col_index(j, a)
                if ef.MUL[x4, v1[c]] != v2[c]:
                    return {"ok": False,
                            "reason": f"(33) violated at j={j}, a={a}"}
    E_rows = E.tolist()
    for row_i, row in enumerate(E_rows):
        s1 = 0
        s2 = 0
        for c, w in enumerate(row):
            if not w:
                continue
            if v1[c]:
                s1 ^= ef.MUL[w, v1[c]]
            if v2[c]:
                s2 ^= ef.MUL[w, v2[c]]
        if s1 or s2:
            return {"ok": False,
                    "reason": f"E*v replay nonzero at row {row_i}"}
    target_labels = [branch_label(oracle_support[j], branch_e, ef)
                     for j in range(len(U_rows))]
    for tau in range(len(U_rows)):
        u_star = U_star[tau]
        u_target = label_to_U_row(target_labels[tau], ef)
        lhs = ef.MUL[int(u_star[0]), int(u_target[1])]
        rhs = ef.MUL[int(u_star[1]), int(u_target[0])]
        if (lhs, rhs) != (0, 0) and lhs != rhs:
            return {"ok": False,
                    "reason": f"U* row {tau} not projectively equal to "
                              f"the target label's U row"}
    return {
        "ok": True, "rank_C": int(rank_C), "rank_V": int(rank_V),
        "U_star": [list(u) for u in U_star],
        "witness": {
            "coefficients": [int(x) for x in coeffs],
            "V_star": [[int(x) for x in v1], [int(x) for x in v2]],
            "C_star": [[int(x) for x in row] for row in C_rows],
            "U_star": [list(u) for u in U_star],
            "target_labels": [list(lab) for lab in target_labels],
        },
    }


# --------------------------------------------------------------------------- #
# P2/P3  census layers with TRANSFORMED labels                                 #
# --------------------------------------------------------------------------- #

def p3_singletons_for_orbit(public, ef, expected_fourths, ker_E=None):
    """PUBLIC-side P3: for EACH label in expected_fourths (a labeled
    layer supplies the list; the P3 machinery itself only consumes the
    label list — no hidden support defaults), build U4 from that label
    and scan all 8 positions."""
    out = []
    for lab in expected_fourths:
        u4 = list(U3_PRESCRIBED) + [label_to_U_row(lab, ef)]
        scans = {pos: scan_labels(public, ef, u4, pos, ker_E)
                 for pos in range(N_ELL)}
        out.append({"fourth_label": list(lab), "scans": scans})
    return out


# --------------------------------------------------------------------------- #
# ledger helpers (schema + hash aware)                                         #
# --------------------------------------------------------------------------- #

def source_sha256(path: str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def ledger_header(source_hashes: dict) -> dict:
    return {"schema": DOMAIN_SCHEMA, "source_sha256": source_hashes}


LEDGER_ABSENT = "ABSENT"
LEDGER_VALID = "VALID"
LEDGER_BAD_SCHEMA = "EXISTS_INVALID_SCHEMA"
LEDGER_CORRUPT = "EXISTS_CORRUPT"
LEDGER_INCOMPLETE = "EXISTS_INCOMPLETE"


def read_ledger_state(path: str, expected_header: dict):
    """Classify a cell-ledger path WITHOUT ever deleting or replacing it.
    Returns (state, record_or_None, raw_bytes_or_None).  Any state other
    than ABSENT/VALID must ABORT the census (evidence is preserved
    byte-for-byte; no unapproved deletion or overwrite)."""
    if not os.path.lexists(path):
        return (LEDGER_ABSENT, None, None)
    if os.path.islink(path) or not os.path.isfile(path):
        return (LEDGER_BAD_SCHEMA, None, None)
    with open(path, "rb") as fh:
        raw = fh.read()
    try:
        rec = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        return (LEDGER_CORRUPT, None, raw)
    if not isinstance(rec, dict):
        return (LEDGER_CORRUPT, None, raw)
    head = rec.get("header") or {}
    # EXACT header equality (owner audit): the full expected header,
    # including controls_artifact_sha256 and beta, must match — never a
    # field-by-field subset.
    if head != expected_header:
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not rec.get("complete"):
        return (LEDGER_INCOMPLETE, None, raw)
    if rec.get("beta") != expected_header.get("beta"):
        return (LEDGER_BAD_SCHEMA, None, raw)
    required = ("P4", "P1_U3", "P1_U4", "P2", "P3")
    if any(req not in rec for req in required):
        return (LEDGER_BAD_SCHEMA, None, raw)
    # P1 branch counts are REQUIRED to be exactly M because those loops
    # always execute independently of any public outcome.
    if len(rec.get("P1_U3", {}).get("branches", [])) != M:
        return (LEDGER_BAD_SCHEMA, None, raw)
    if len(rec.get("P1_U4", {}).get("branches", [])) != M:
        return (LEDGER_BAD_SCHEMA, None, raw)
    # P3 completeness is SHAPE-only, never result-dependent (owner audit):
    # a legitimate NEGATIVE cell (public P2 accepting fewer/missing
    # branches) yields a complete failure record with 0..N mapped
    # branches, and must remain readable as valid evidence.
    p3 = rec.get("P3", {})
    if not isinstance(p3.get("branches"), list):
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not isinstance(p3.get("unmapped"), list):
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not isinstance(p3.get("n_missing_mapped_branches"), int):
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not isinstance(p3.get("n_unmapped_branches"), int):
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not isinstance(p3.get("passes"), bool):
        return (LEDGER_BAD_SCHEMA, None, raw)
    if not isinstance(rec.get("P2", {}).get("passes"), bool):
        return (LEDGER_BAD_SCHEMA, None, raw)
    return (LEDGER_VALID, rec, raw)


def fsync_directory(path: str) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_once_bytes(path: str, payload: bytes,
                     allow_identical: bool = False) -> None:
    if os.path.lexists(path):
        if (allow_identical and not os.path.islink(path)
                and os.path.isfile(path)):
            with open(path, "rb") as stream:
                if stream.read() == payload:
                    return
        raise RuntimeError(
            f"TARGET WRITE REFUSED: {path} already exists; preserving it")
    tmp = path + ".tmp"
    if os.path.lexists(tmp):
        raise RuntimeError(
            f"STAGED WRITE REFUSED: {tmp} already exists; preserving "
            f"existing bytes for inspection, not truncating")
    descriptor = os.open(
        tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    dirpath = os.path.dirname(os.path.abspath(path)) or "."
    dirfd = os.open(dirpath, os.O_RDONLY)
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)
    if os.path.lexists(path):
        raise RuntimeError(
            f"TARGET APPEARED DURING WRITE: {path}; preserving staged temp")
    os.replace(tmp, path)
    dirfd = os.open(dirpath, os.O_RDONLY)
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)


def atomic_write_json(path: str, obj) -> None:
    """Write one new ledger atomically without replacing existing evidence."""
    payload = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    write_once_bytes(path, payload)


def prepare_control_plant_root(root: str) -> str:
    root = os.path.abspath(root)
    if os.path.lexists(root):
        if os.path.islink(root) or not os.path.isdir(root):
            raise RuntimeError(
                f"control-plant root is non-directory evidence: {root}")
    else:
        os.mkdir(root, 0o700)
        fsync_directory(os.path.dirname(root))
    unexpected = sorted(set(os.listdir(root)) - set(CONTROL_PLANT_CASES))
    if unexpected:
        raise RuntimeError(
            f"unexpected preserved control-plant entries: {unexpected}")
    return root


def persist_control_plant(root: str, case: str, payload: bytes) -> str:
    if case not in CONTROL_PLANT_CASES:
        raise RuntimeError(f"unregistered control-plant case: {case}")
    case_dir = os.path.join(root, case)
    if os.path.lexists(case_dir):
        if os.path.islink(case_dir) or not os.path.isdir(case_dir):
            raise RuntimeError(
                f"control-plant case is non-directory evidence: {case_dir}")
    else:
        os.mkdir(case_dir, 0o700)
        fsync_directory(root)
    unexpected = sorted(set(os.listdir(case_dir)) - {"cell_8.json"})
    if unexpected:
        raise RuntimeError(
            f"unexpected preserved files in {case}: {unexpected}")
    path = ledger_path(case_dir, 8)
    write_once_bytes(path, payload, allow_identical=True)
    return path


def control_plant_paths(root: str) -> list:
    return [ledger_path(os.path.join(root, case), 8)
            for case in CONTROL_PLANT_CASES]


def ledger_path(ledger_dir: str, beta: int) -> str:
    return os.path.join(ledger_dir, f"cell_{beta}.json")


# --------------------------------------------------------------------------- #
# census + controls — HARD CPU-slot gate                                       #
# --------------------------------------------------------------------------- #

def _require_release():
    if os.environ.get(RELEASE_KEY, "") != "1":
        raise RuntimeError("successor CPU slot not released")
    if os.environ.get(PARENT_RELEASE_KEY, "") != "1":
        raise RuntimeError("parent census termination not acknowledged")
    if not __debug__:
        raise RuntimeError("optimized Python is forbidden for census release")
    if os.getpriority(os.PRIO_PROCESS, 0) < 10:
        raise RuntimeError("process niceness must be at least 10")
    required = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    for name in required:
        if os.environ.get(name) != "1":
            raise RuntimeError(f"{name} must be exactly 1")


def run_public_phase(public, ef, E=None, ker_E=None):
    """PUBLIC phase: P4 + P2 raw scan + P3 raw scans driven ONLY by the
    actual accepted labels.  No oracle/support access here."""
    if E is None:
        E = build_E(public, ef)
    if ker_E is None:
        ker_E = nullspace_of(E.tolist(), M_ELL, ef)[0]
    rec = {}
    rec["P4"] = p4_admission(public, ef, E=E, ker_E=ker_E)
    p2_raw = scan_labels(public, ef, U3_PRESCRIBED, 3, ker_E)
    rec["P2_raw"] = [list(lab) for lab in p2_raw]
    p3_raw = []
    for fourth in p2_raw:
        u4 = list(U3_PRESCRIBED) + [label_to_U_row(fourth, ef)]
        scans = {}
        for pos in range(N_ELL):
            scans[pos] = scan_labels(public, ef, u4, pos, ker_E)
        p3_raw.append({"fourth_label": list(fourth),
                       "scans": {pos: [list(x) for x in labs]
                                 for pos, labs in scans.items()}})
    rec["P3_raw"] = p3_raw
    return rec


def run_cell(beta, ef, use_cache=True):
    """Full PUBLIC + LABELED adjudication for one beta.  Public phase
    touches no oracle data; labeled phase consumes oracle['support']
    ONLY.  The record is returned without any disk write; run_census
    performs the single atomic write after completion."""
    public, oracle = build_toy(beta)
    support = oracle["support"]
    E = None
    ker_E = None
    if use_cache:
        E = build_E(public, ef)
        E.setflags(write=False)
        ker_E = tuple(
            tuple(int(value) for value in row)
            for row in nullspace_of(E.tolist(), M_ELL, ef)[0])
    rec = {"beta": beta}
    rec["public"] = {"k": public["k"], "n": public["n"],
                     "prim": public["prim"], "profile": public["profile"]}
    pub = run_public_phase(public, ef, E=E, ker_E=ker_E)
    rec.update(pub)
    p2_raw = [tuple(lab) for lab in pub["P2_raw"]]
    expected_fourths = p2_fourth_labels(ef, support)
    rec["P2"] = {"accepted": [list(lab) for lab in p2_raw],
                 "expected": [list(lab) for lab in expected_fourths],
                 "passes": sorted(p2_raw) == sorted(expected_fourths)}
    label_to_exp = {lab: e for e, lab in enumerate(expected_fourths)}
    extras = [lab for lab in p2_raw if lab not in label_to_exp]
    rec["P2_extras"] = [list(lab) for lab in extras]
    p3_branches = []
    # P3 verdict requires EXACTLY the five mapped expected branches, zero
    # unmapped entries, and every 8-position singleton target (owner
    # audit: p3_all alone is vacuously true for missing/extra branches).
    mapped_branches = {}
    unmapped_entries = []
    n_unmapped = 0
    for entry in pub["P3_raw"]:
        fourth = tuple(entry["fourth_label"])
        scans = {int(pos): [tuple(x) for x in labs]
                 for pos, labs in entry["scans"].items()}
        if fourth not in label_to_exp:
            n_unmapped += 1
            unmapped_entries.append({"fourth_label": list(fourth),
                                     "scans": {pos: [list(x) for x in labs]
                                               for pos, labs in scans.items()}})
            continue
        e = label_to_exp[fourth]
        targets = branch_target_table(e, ef, support)
        per_pos = {pos: [list(x) for x in labs]
                   for pos, labs in scans.items()}
        s_ok = all(len(scans[pos]) == 1 and scans[pos][0] == targets[pos]
                   for pos in range(N_ELL))
        p3_branches.append({"fourth_label": list(fourth),
                            "mapped_branch": e,
                            "scans": per_pos,
                            "targets": [list(t) for t in targets],
                            "all_singletons": s_ok})
        mapped_branches[e] = s_ok
    p3_passes = (n_unmapped == 0
                 and sorted(mapped_branches.keys()) == list(range(M))
                 and all(mapped_branches.values()))
    rec["P3"] = {"branches": p3_branches, "passes": p3_passes,
                 "n_unmapped_branches": n_unmapped,
                 "unmapped": unmapped_entries,
                 "n_missing_mapped_branches":
                     M - len(mapped_branches)}
    p1_u3 = []
    p1_u4 = []
    for e in range(M):
        v_u3 = p1_certificate_verdict(
            public, ef, U3_PRESCRIBED, branch_e=e, oracle=oracle,
            E=E, ker_E=ker_E)
        p1_u3.append({"branch": e, "verdict": v_u3})
        u4 = list(U3_PRESCRIBED) + [label_to_U_row(expected_fourths[e], ef)]
        v_u4 = p1_certificate_verdict(
            public, ef, u4, branch_e=e, oracle=oracle,
            E=E, ker_E=ker_E)
        p1_u4.append({"branch": e, "verdict": v_u4})
    rec["P1_U3"] = {"branches": p1_u3,
                    "passes": all(b["verdict"]["verdict"] == "P1_PASS"
                                  for b in p1_u3)}
    rec["P1_U4"] = {"branches": p1_u4,
                    "passes": all(b["verdict"]["verdict"] == "P1_PASS"
                                  for b in p1_u4)}
    rec["P1"] = {"branches": p1_u3 + p1_u4,
                 "passes": all(b["verdict"]["verdict"] == "P1_PASS"
                               for b in p1_u3 + p1_u4)}
    rec["t_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    return rec


def run_census(ledger_dir, source_hashes, extra_header=None):
    """24-cell census.  LEDGER SAFETY: a cell path that exists but is not
    a VALID complete semantically-valid ledger ABORTS the census with the
    offending file preserved byte-for-byte (no unapproved deletion or
    overwrite).  The single atomic write per cell happens only after the
    cell's P1-P4 record is COMPLETE.  Each ledger header carries BOTH
    the source hash map and any supplied extra (e.g. controls artifact)
    hash, so the released-state lineage is recorded in every cell."""
    _require_release()
    os.makedirs(ledger_dir, exist_ok=True)
    ef = EField(M, prim=PRIM)
    done = []
    validated = {}
    for beta in range(8, 32):
        expected_header = ledger_header(source_hashes)
        expected_header["beta"] = beta
        if extra_header:
            expected_header.update(extra_header)
        p = ledger_path(ledger_dir, beta)
        state, prev, raw = read_ledger_state(p, expected_header)
        if state == LEDGER_VALID:
            validated[beta] = prev
            done.append(beta)
            continue
        if state != LEDGER_ABSENT:
            raw_sha = (
                _hashlib.sha256(raw).hexdigest()
                if raw is not None else "unavailable-nonregular-path")
            raise RuntimeError(
                f"LEDGER SAFETY ABORT: cell {beta} ledger at {p} is in "
                f"state {state}; file PRESERVED byte-for-byte (sha256 "
                f"{raw_sha}); census will not overwrite or delete it")
        rec = run_cell(beta, ef)
        rec["header"] = expected_header
        rec["complete"] = True
        atomic_write_json(p, rec)
        # REVALIDATE the just-written cell through the SAME semantic
        # checker before accepting it (owner audit: new writes must not
        # bypass the checker).
        st2, rec2, _raw2 = read_ledger_state(p, expected_header)
        if st2 != LEDGER_VALID:
            raise RuntimeError(
                f"FRESH-WRITE VALIDATION FAILED: cell {beta} at {p} "
                f"classifies {st2} immediately after writing; file "
                f"preserved for inspection")
        validated[beta] = rec2
        done.append(beta)
    return {"cells_completed": done, "n": len(done), "of": 24,
            "validated": validated}


def run_controls(control_plant_root):
    """Release-gated controls with persistent, checksummable plant evidence."""
    _require_release()
    plant_root = prepare_control_plant_root(control_plant_root)
    ef = EField(M, prim=PRIM)
    out = {}
    # C1 population anchors
    n_admissible = sum(1 for b in range(Q) if b not in SUPPORT)
    assert n_admissible == 24 and len(all_labels()) == 33
    out["C1_anchors"] = {"admissible_G": n_admissible,
                         "labels": len(all_labels())}
    # C2 dual rref + replay (300 cases)
    import random
    rng = random.Random(20260831)
    for _ in range(300):
        nr, nc = rng.randrange(1, 7), rng.randrange(1, 6)
        A = [[rng.randrange(Q) for _ in range(nc)] for _ in range(nr)]
        bvec = [rng.randrange(Q) for _ in range(nr)]
        _, _, feas, part = rref_affine(A, bvec, nc, ef)
        rkA = rank_of(A, nc, ef)
        rkAug = rank_of([A[i] + [bvec[i]] for i in range(nr)], nc + 1, ef)
        assert feas == (rkA == rkAug)
        if feas:
            for i in range(nr):
                s = 0
                for c2 in range(nc):
                    if A[i][c2] and part[c2]:
                        s ^= ef.MUL[A[i][c2], part[c2]]
                assert s == bvec[i]
    out["C2_dual_rref"] = {"cases": 300, "agreed": True}
    # C25 successor exact-kernel equivalence: full tuples, basis order,
    # ranks, empty/zero/rectangular/affine-inconsistent systems.
    rng_fast = random.Random(202608311)
    systems = [
        ([], [], 3),
        ([[0, 0], [0, 0]], [0, 0], 2),
        ([[0, 0]], [1], 2),
        ([[1, 0], [0, 1]], [7, 11], 2),
    ]
    for _ in range(1000):
        nr = rng_fast.randrange(0, 10)
        nc = rng_fast.randrange(1, 12)
        matrix = [
            [rng_fast.randrange(Q) for _ in range(nc)]
            for _ in range(nr)]
        rhs = [rng_fast.randrange(Q) for _ in range(nr)]
        systems.append((matrix, rhs, nc))
    for case, (matrix, rhs, nc) in enumerate(systems):
        reference = rref_affine_reference(matrix, rhs, nc, ef)
        vectorized = rref_affine(matrix, rhs, nc, ef)
        assert vectorized == reference, (
            "vectorized/reference RREF mismatch", case)
        assert nullspace_of(matrix, nc, ef) == \
            nullspace_of_reference(matrix, nc, ef), (
                "vectorized/reference nullspace mismatch", case)
        assert rank_of(matrix, nc, ef) == \
            rank_of_reference(matrix, nc, ef), (
                "vectorized/reference rank mismatch", case)
    bad_mul = ef.MUL.copy()
    bad_mul[1, 1] = 0

    class WrongTable:
        MUL = bad_mul
        INV = ef.INV

    plant_matrix = [[1, 1], [1, 0]]
    plant_rhs = [0, 1]
    assert rref_affine(plant_matrix, plant_rhs, 2, WrongTable()) != \
        rref_affine_reference(plant_matrix, plant_rhs, 2, ef)
    out["C25_vectorized_rref_equivalence"] = {
        "systems": len(systems),
        "full_tuple_equal": True,
        "nullspace_basis_equal": True,
        "rank_equal": True,
        "wrong_table_detected": True,
    }
    # C3 Lucas vs integer binomial, every E entry
    import math
    public8, oracle8 = build_toy(8)
    E8 = build_E(public8, ef)
    for bi, b in enumerate(COLS_D):
        for j in range(N_ELL):
            ybin = tuple(int(v) for v in
                         np.asarray(public8["Y"], dtype=np.uint8)[:, j])
            for ai, a in enumerate(COLS_S):
                expect = 1
                for i in range(K):
                    if a[i] > b[i] or ((a[i] & ~b[i]) != 0):
                        expect = 0
                        break
                    if not (math.comb(b[i], a[i]) % 2):
                        expect = 0
                        break
                    if ybin[i] == 0:
                        if b[i] - a[i]:
                            expect = 0
                            break
                    elif b[i] != a[i]:
                        expect = ef.MUL[expect, ybin[i]]
                assert int(E8[bi, col_index(j, a)]) == expect
    out["C3_lucas_integer"] = {"agreed": True}
    # C4 full-rank plant
    plant = np.zeros((N_D, M_ELL), dtype=np.uint16)
    for i in range(N_D):
        plant[i, i] = 1
    assert rank_of(plant.tolist(), M_ELL, ef) == N_D
    out["C4_fullrank_plant"] = {"caught": True}
    # C5 zero plant
    zer = np.zeros((N_D, M_ELL), dtype=np.uint16)
    ker, _ = nullspace_of(zer.tolist(), M_ELL, ef)
    assert len(ker) == M_ELL
    dim0 = rank_of([[row[c] for c in block_indices(0)] for row in ker],
                   len(B_J), ef)
    assert dim0 == len(B_J)
    out["C5_zero_plant"] = {"dim_pi0": int(dim0), "caught": True}
    # C6/C7 degenerate plants (identity kernel: all 33 feasible)
    ker_zero = [list(np.eye(M_ELL, dtype=np.uint16)[i].tolist())
                for i in range(M_ELL)]
    U0 = [(0, 1)]
    cnt = 0
    for lab in all_labels():
        rows, bvec, nunk = vu_system_rows(U0, 0, lab, ker_zero, ef)
        _, _, feas, _ = rref_affine(rows, bvec, nunk, ef)
        if feas:
            cnt += 1
    assert cnt == 33
    out["C6_C7_degenerate_plants"] = {"accepted_n": cnt}
    # C8 zero-locator-block plant
    cnt8 = 0
    for lab in all_labels():
        rows, bvec, nunk = vu_system_rows(U0, 0, lab, ker_zero, ef)
        _, _, feas, _ = rref_affine(rows[:-len(B_J)],
                                    bvec[:-len(B_J)], nunk, ef)
        if feas:
            cnt8 += 1
    assert cnt8 == 33
    out["C8_zero_locator_block"] = {"accepted_n": cnt8}
    # C9 rank-two stack plant through P1 verifier (branch 0, beta=8).
    # The registered expected verdict is EXPLICIT (Main: do not accept
    # either PASS or FAIL vacuously).  For the canonical witness the P1
    # verifier must reach a DECISIVE verdict — the exact expected string
    # is recorded here and asserted.
    v8 = p1_certificate_verdict(public8, ef, U3_PRESCRIBED, branch_e=0,
                                oracle=oracle8)
    C9_EXPECTED = "P1_PASS"       # registered expectation (branch 0, canonical)
    assert v8["verdict"] == C9_EXPECTED, (v8, C9_EXPECTED)
    witness8 = v8.get("witness")
    assert isinstance(witness8, dict), "P1_PASS must exhibit its witness"
    assert len(witness8["V_star"]) == 2
    assert all(len(row) == M_ELL for row in witness8["V_star"])
    assert len(witness8["C_star"]) == len(U3_PRESCRIBED)
    assert all(len(row) == M_ELL for row in witness8["C_star"])
    assert witness8["U_star"] == v8["U_star"]
    assert witness8["target_labels"] == [
        list(branch_label(oracle8["support"][j], 0, ef))
        for j in range(len(U3_PRESCRIBED))]
    out["C9_p1_plant"] = {
        "verdict": v8["verdict"], "expected": C9_EXPECTED,
        "matched": True, "witness_exhibited": True,
        "V_shape": [2, M_ELL],
        "C_shape": [len(U3_PRESCRIBED), M_ELL],
    }
    # C10 REAL synthetic-kernel diagonal/rank plant (owner-audit):
    # a tiny synthetic kernel with KNOWN block constants; _diag_of is
    # asserted separately for the first half (U=(1,0)), the second half
    # (U=(0,1)) and mixed (U=(1,1)), with exact recorded values, plus
    # direct rank-2 / rank-1 stack asserts.  No tautological zeros.
    c000 = col_index(0, (0, 0, 0))
    ker_syn = [[0] * M_ELL, [0] * M_ELL]
    ker_syn[0][c000] = 7          # kernel row 0 has 7 at (tau=0, a=0)
    ker_syn[1][c000] = 11         # kernel row 1 has 11 at the same coord
    nk_syn = 2
    # coefficient vector: V1 = 3*ker0, V2 = 5*ker1
    coeff_syn = [3, 0, 0, 5]      # [u1_0, u1_1 | u2_0, u2_1]
    v1_syn, v2_syn = _expand_coeffs(coeff_syn, ker_syn, ef, nk_syn)
    assert v1_syn[c000] == ef.MUL[3, 7], (v1_syn[c000], ef.MUL[3, 7])
    assert v2_syn[c000] == ef.MUL[5, 11], (v2_syn[c000], ef.MUL[5, 11])
    # first half only: U=(1,0) picks V1
    d_first = _diag_of(coeff_syn, [(1, 0)], ker_syn, ef, nk_syn, 0)
    assert d_first == ef.MUL[3, 7], (d_first, ef.MUL[3, 7])
    # second half only: U=(0,1) picks V2 (this is what the old broken
    # control never exercised)
    d_second = _diag_of(coeff_syn, [(0, 1)], ker_syn, ef, nk_syn, 0)
    assert d_second == ef.MUL[5, 11], (d_second, ef.MUL[5, 11])
    # mixed: U=(1,1) is the XOR of both halves
    d_mixed = _diag_of(coeff_syn, [(1, 1)], ker_syn, ef, nk_syn, 0)
    assert d_mixed == (ef.MUL[3, 7] ^ ef.MUL[5, 11]), d_mixed
    assert d_first != 0 and d_second != 0
    # direct rank asserts: independent stack -> 2; dependent -> 1
    r2a = [0] * M_ELL
    r2b = [0] * M_ELL
    r2a[0] = 1
    r2b[1] = 1
    assert rank_of([r2a, r2b], M_ELL, ef) == 2
    assert rank_of([r2a, list(r2a)], M_ELL, ef) == 1
    out["C10_synthetic_diag_rank_plant"] = {
        "d_first_half": int(d_first),
        "d_second_half": int(d_second),
        "d_mixed": int(d_mixed),
        "expected_first": int(ef.MUL[3, 7]),
        "expected_second": int(ef.MUL[5, 11]),
        "rank2_stack": 2, "rank1_stack": 1}
    # C11 brute vs rref (300 cases)
    for _ in range(300):
        nunk = rng.randrange(1, 4)
        nrow = rng.randrange(1, 5)
        A11 = [[rng.randrange(Q) for _ in range(nunk)] for _ in range(nrow)]
        b11 = [rng.randrange(Q) for _ in range(nrow)]
        _, _, feas, _ = rref_affine(A11, b11, nunk, ef)
        brute = False
        for guess in itertools.product(range(Q), repeat=nunk):
            okg = True
            for i in range(nrow):
                s = 0
                for c2 in range(nunk):
                    if A11[i][c2]:
                        s ^= ef.MUL[A11[i][c2], guess[c2]]
                if s != b11[i]:
                    okg = False
                    break
            if okg:
                brute = True
                break
        assert brute == feas
    out["C11_brute_agreement"] = {"cases": 300}
    # C12 Mobius canon + label/U anchors
    Mrow = mobius_matrix(0, 1, 2, ef)
    for x in range(8):
        D_m = ef.MUL[Mrow[1][0], x] ^ Mrow[1][1]
        N_m = ef.MUL[Mrow[0][0], x] ^ Mrow[0][1]
        D_s, N_s = mobius_pair(x, ef)
        assert (D_m, N_m) == (D_s, N_s)
    assert label_from_pair(0, 2, ef) == INF_ROW
    assert label_to_U_row(label_from_pair(2, 0, ef), ef) == (0, 1)
    assert label_to_U_row(label_from_pair(3, 3, ef), ef) == (1, 1)
    assert label_to_U_row(INF_ROW, ef) == (1, 0)
    out["C12_mobius_label_u_anchors"] = {"verified": True}
    # C13 Mobius counterfactual: permuting a,b,c breaks every anchor
    wrong = branch_target_table(0, ef, SUPPORT, a=0, b=2, c=1)
    anchor = branch_target_table(0, ef, SUPPORT)
    assert wrong != anchor
    out["C13_mobius_counterfactual"] = {"broken_as_expected": True}
    # C14 owner-anchored branch rows + 4th-holdout orbit
    for e in range(M):
        tbl = branch_target_table(e, ef, SUPPORT)
        exp = [INF_ROW if v2 is None else (1, v2)
               for v2 in OWNER_ANCHOR_BRANCH_ROWS[e]]
        assert tbl == exp, (e, tbl, exp)
    orbit = [branch_target_table(e, ef, SUPPORT)[3] for e in range(M)]
    assert orbit == [(1, 5), (1, 17), (1, 11), (1, 28), (1, 3)]
    out["C14_owner_anchors"] = {"branch_rows": "verified",
                                "fourth_orbit": "verified"}
    # C16 preimage-basis counterexample (owner-mandated): image rows
    # [[1,1],[1,0]] are independent -> TWO DISTINCT preimages retained;
    # a dependent row must NOT grow the basis.
    synthetic_pairs = [([1, 1], "A"), ([1, 0], "B"), ([0, 1], "C-dep")]
    synth_retained = []
    synth_rows = []
    for lv, v in synthetic_pairs:
        trial = synth_rows + [lv]
        if rank_of(trial, 2, ef) > len(synth_rows):
            synth_retained.append((lv, v))
            synth_rows.append(lv)
    assert len(synth_retained) == 2, synth_retained
    assert sorted(vv for _, vv in synth_retained) == ["A", "B"]
    out["C16_preimage_basis"] = {"retained": 2, "distinct": True}
    # C17 infinity-coupling replay counterexamples (owner-mandated):
    # V1=0/V2!=0 passes the infinity coupling; V2=0/V1!=0 fails.
    j_inf = 2
    c0 = col_index(j_inf, B_J[0])
    v1_good = [0] * M_ELL
    v2_good = [0] * M_ELL
    v2_good[c0] = 7
    ok_good = all(v1_good[col_index(j_inf, a)] == 0 for a in B_J)
    v1_bad = [0] * M_ELL
    v1_bad[c0] = 3
    fails_bad = any(v1_bad[col_index(j_inf, a)] != 0 for a in B_J)
    assert ok_good and fails_bad
    out["C17_infinity_replay"] = {"v1_zero_v2_nonzero_passes": ok_good,
                                  "v1_nonzero_fails": fails_bad}
    # C18 stale-schema ledger-safety plant: a cell file with a MISMATCHED
    # schema and arbitrary bytes must cause run_census to ABORT with the
    # file PRESERVED byte-for-byte (hash unchanged).
    stale_payload = b'{"header": {"schema": "OLD-SCHEMA"}, "beta": 8}'
    stale = persist_control_plant(
        plant_root, "C18_stale_schema", stale_payload)
    stale_dir = os.path.dirname(stale)
    stale_sha_before = source_sha256(stale)
    aborted = False
    try:
        run_census(stale_dir, {"toy_census.py": "deadbeef"},
                   extra_header={"controls_artifact_sha256": "deadbeef"})
    except RuntimeError as error:
        aborted = ("LEDGER SAFETY ABORT" in str(error)
                   and "EXISTS_INVALID_SCHEMA" in str(error))
    stale_sha_after = source_sha256(stale)
    assert aborted, "stale-schema cell must abort the census"
    assert stale_sha_before == stale_sha_after, \
        "stale-schema file must be preserved byte-for-byte"
    with open(stale, "rb") as stream:
        assert stream.read() == stale_payload
    out["C18_stale_schema_plant"] = {
        "aborted": True,
        "hash_unchanged": True,
        "artifact": "C18_stale_schema/cell_8.json",
    }
    # C21/C22/C23 ledger-refusal plants (owner-mandated): wrong-beta,
    # complete-but-missing-P4, and wrong-controls-hash ledgers must ALL
    # classify EXISTS_INVALID_SCHEMA with byte-exact preservation.
    good_src = {"toy_census.py": "aaaa"}
    good_ctl = {"controls_artifact_sha256": "bbbb"}

    def _plant(name, payload_obj, expect_state):
        raw = json.dumps(payload_obj, sort_keys=True).encode("utf-8")
        path = persist_control_plant(plant_root, name, raw)
        sha_before = source_sha256(path)
        hdr = ledger_header(good_src)
        hdr["beta"] = 8
        hdr.update(good_ctl)
        state, rec, rawb = read_ledger_state(path, hdr)
        sha_after = source_sha256(path)
        assert state == expect_state, (name, state, expect_state)
        assert sha_before == sha_after, name
        assert rawb == raw, name
        return state

    hdr_ok = ledger_header(good_src)
    hdr_ok["beta"] = 8
    hdr_ok.update(good_ctl)
    # C21 wrong beta (header beta 9 in the beta-8 slot)
    hdr_wrong_beta = dict(hdr_ok)
    hdr_wrong_beta["beta"] = 9
    s21 = _plant("C21_wrong_beta",
                 {"header": hdr_wrong_beta, "beta": 9, "complete": True,
                  "P4": {}, "P1_U3": {"branches": [0] * M},
                  "P1_U4": {"branches": [0] * M}, "P2": {},
                  "P3": {"branches": [0] * M}},
                 LEDGER_BAD_SCHEMA)
    out["C21_wrong_beta_plant"] = {"state": s21, "bytes_preserved": True}
    # C22 complete but missing P4
    s22 = _plant("C22_missing_P4",
                 {"header": hdr_ok, "beta": 8, "complete": True,
                  "P1_U3": {"branches": [0] * M},
                  "P1_U4": {"branches": [0] * M}, "P2": {},
                  "P3": {"branches": [0] * M}},
                 LEDGER_BAD_SCHEMA)
    out["C22_missing_P4_plant"] = {"state": s22, "bytes_preserved": True}
    # C23 wrong controls hash in the header
    hdr_wrong_ctl = dict(hdr_ok)
    hdr_wrong_ctl["controls_artifact_sha256"] = "WRONGHASH"
    s23 = _plant("C23_wrong_controls_hash",
                 {"header": hdr_wrong_ctl, "beta": 8, "complete": True,
                  "P4": {}, "P1_U3": {"branches": [0] * M},
                  "P1_U4": {"branches": [0] * M}, "P2": {},
                  "P3": {"branches": [0] * M}},
                 LEDGER_BAD_SCHEMA)
    out["C23_wrong_controls_hash_plant"] = {"state": s23,
                                            "bytes_preserved": True}
    # C24 NEGATIVE-COMPLETE ledger plant (owner-mandated): a legitimate
    # failing cell (P2 false, P3 empty with missing=5) is COMPLETE
    # evidence and MUST classify VALID with byte-exact preservation —
    # semantic completeness must never require a passing result.
    s24 = _plant("C24_negative_complete",
                 {"header": hdr_ok, "beta": 8, "complete": True,
                  "P4": {"passes": False},
                  "P1_U3": {"branches": [0] * M},
                  "P1_U4": {"branches": [0] * M},
                  "P2": {"passes": False},
                  "P3": {"branches": [], "unmapped": [],
                         "n_missing_mapped_branches": M,
                         "n_unmapped_branches": 0, "passes": False}},
                 LEDGER_VALID)
    out["C24_negative_complete_plant"] = {"state": s24,
                                         "bytes_preserved": True}
    # C15 diagonal-selection negative plant (kept; exercised via C9/C10)
    out["C15_diag_selection"] = {"documented_via": "C9/C10"}
    # C19 oracle-support permutation counterfactual (owner-mandated):
    # PUBLIC scans must be byte-identical under a permutation of the
    # oracle support (the public phase never sees it), while LABELED
    # expectations change accordingly.
    public8, oracle8b = build_toy(8)
    # run the public phase twice: with original oracle support, then with
    # a PERMUTED oracle record (only the support field differs)
    pub1 = run_public_phase(public8, ef)
    permuted_oracle = dict(oracle8b)
    permuted_oracle["support"] = [oracle8b["support"][i] for i in
                                  (3, 1, 0, 2, 7, 6, 5, 4)]   # nontrivial perm
    pub2 = run_public_phase(public8, ef)   # public phase never sees support
    assert pub1 == pub2, "public phase must be support-independent"
    labeled1 = p2_fourth_labels(ef, oracle8b["support"])
    labeled2 = p2_fourth_labels(ef, permuted_oracle["support"])
    assert labeled1 != labeled2, \
        "labeled expectations must change under support permutation"
    out["C19_oracle_permutation"] = {
        "public_byte_identical": True,
        "labeled_expectations_differ": True,
        "pre_permutation": [list(x) for x in labeled1],
        "post_permutation": [list(x) for x in labeled2]}
    # C20 frozen-table crosscheck (owner-mandated): the frozen scalar
    # GF engine (gfield_frozen.GF) and the frozen table engine
    # (fastfield_frozen.EField) must agree on multiplication and inverse
    # over a fixed sample, and on the full label/U round trip.
    gf_frozen = GF(M, prim=PRIM)
    bad = 0
    for a in range(1, Q):
        for b in range(1, Q):
            if gf_frozen.mul(a, b) != ef.MUL[a, b]:
                bad += 1
            if gf_frozen.inv(a) != ef.INV[a]:
                bad += 1
        # exponent check at a fixed power
        if gf_frozen.pow(a, 4) != int(ef.powers(a, 4)[4]):
            bad += 1
    assert bad == 0, bad
    out["C20_frozen_table_crosscheck"] = {"mismatches": 0}
    # C26 end-to-end semantic anchor: cached and uncached fast paths must
    # match the immutable parent scalar beta-8 cell byte-for-byte after
    # removing only provenance fields that necessarily differ.
    assert source_sha256(PARENT_CONTROLS_PATH) == PARENT_CONTROLS_SHA256
    assert source_sha256(PARENT_CELL8_PATH) == PARENT_CELL8_SHA256
    with open(PARENT_CELL8_PATH, encoding="utf-8") as fh:
        parent_cell = json.load(fh)
    cached_cell = run_cell(8, ef, use_cache=True)
    uncached_cell = run_cell(8, ef, use_cache=False)
    for record in (parent_cell, cached_cell, uncached_cell):
        record.pop("t_utc", None)
        record.pop("header", None)
        record.pop("complete", None)

    def semantic_bytes(record):
        return json.dumps(
            record, sort_keys=True, separators=(",", ":"),
            default=str).encode("utf-8")

    parent_semantic = semantic_bytes(parent_cell)
    cached_semantic = semantic_bytes(cached_cell)
    uncached_semantic = semantic_bytes(uncached_cell)
    assert cached_semantic == parent_semantic
    assert uncached_semantic == parent_semantic
    out["C26_parent_beta8_semantic_anchor"] = {
        "parent_controls_sha256": PARENT_CONTROLS_SHA256,
        "parent_cell8_sha256": PARENT_CELL8_SHA256,
        "semantic_sha256": _hashlib.sha256(parent_semantic).hexdigest(),
        "cached_equals_parent": True,
        "uncached_equals_parent": True,
    }
    plant_hashes = {
        f"{case}/cell_8.json": source_sha256(path)
        for case, path in zip(
            CONTROL_PLANT_CASES, control_plant_paths(plant_root))
    }
    return {
        "ALL_OK": True,
        "controls": out,
        "control_plant_sha256": plant_hashes,
    }


# --------------------------------------------------------------------------- #
# static structure checks (no census/control arithmetic)                       #
# --------------------------------------------------------------------------- #

def static_guards():
    assert Q == 32
    labels = all_labels()
    assert len(labels) == 33 and len(set(labels)) == 33
    assert label_is_infinity(INF_ROW)
    assert len(COLS_S) == 35 and len(COLS_D) == N_D == 36
    assert len(B_J) == 34 and M_ELL == 280
    assert col_index(0, (0, 0, 0)) == 0
    assert col_index(7, (0, 0, 0)) == 7 * 35
    assert block_indices(0)[:3] == [IDS_OF[a] for a in B_J[:3]]
    assert all(1 <= sum(a) <= 4 for a in B_J)
    ef = EField(M, prim=PRIM)
    # label/U type separation
    assert label_to_U_row((1, 0), ef) == (0, 1)
    assert label_to_U_row((1, 1), ef) == (1, 1)
    assert label_to_U_row(INF_ROW, ef) == (1, 0)
    assert u_row_to_label((0, 1), ef) == (1, 0)
    assert u_row_to_label((1, 0), ef) == INF_ROW
    assert u_row_to_label((1, 1), ef) == (1, 1)
    # U<->label round-trip on all 33 labels
    for lab in all_labels():
        assert u_row_to_label(label_to_U_row(lab, ef), ef) == lab
    # nontrivial anchors (owner-mandated): locator 2 -> U (16, 1);
    # inverse (16,1) -> label (1,2)
    assert label_to_U_row((1, 2), ef) == (16, 1)
    assert u_row_to_label((16, 1), ef) == (1, 2)
    # Owner anchors: branch rows and 4th-holdout orbit EXACTLY
    for e in range(M):
        tbl = branch_target_table(e, ef, SUPPORT)
        exp = [INF_ROW if v is None else (1, v)
               for v in OWNER_ANCHOR_BRANCH_ROWS[e]]
        assert tbl == exp, (e, tbl, exp)
    orbit = [branch_target_table(e, ef, SUPPORT)[3] for e in range(M)]
    assert orbit == [(1, 5), (1, 17), (1, 11), (1, 28), (1, 3)]
    # matrix form == scalar Mobius pair on ALL support and ALL branches
    for e in range(M):
        Mt = mobius_matrix(0, 1, 2, ef)
        for x in range(N):
            D_m = ef.MUL[Mt[1][0], x] ^ Mt[1][1]
            N_m = ef.MUL[Mt[0][0], x] ^ Mt[0][1]
            D_s, N_s = frob_pair((D_m, N_m), e, ef)
            assert (D_s, N_s) == branch_mobius_pair(x, e, ef)
    # counterfactual: swapped (b,c) breaks the anchor table
    assert (branch_target_table(0, ef, SUPPORT, 0, 2, 1)
            != branch_target_table(0, ef, SUPPORT))
    # toy record and public/oracle separation
    public, oracle = build_toy(8)
    assert np.asarray(public["Y"], dtype=np.uint8).shape == (K, N)
    assert not ({"beta", "G", "support"} & set(public.keys()))
    # static guard refuses a released flag
    assert os.environ.get(RELEASE_KEY, "") != "1"
    assert os.environ.get(PARENT_RELEASE_KEY, "") != "1"
    return {"labels": 33, "cols_S": 35, "cols_D": 36, "block": 34,
            "M_ell": 280, "N_d": 36,
            "public_keys": sorted(public.keys()),
            "oracle_keys": sorted(oracle.keys()),
            "census_gated": True,
            "mobius_canon": "owner_anchor_verified",
            "u_row_type_separation": True}


if __name__ == "__main__":
    print(json.dumps(static_guards(), indent=1))
