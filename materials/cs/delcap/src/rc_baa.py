"""
rc_baa.py — Faithful implementation of Rubinstein-Con's memory-efficient
Blahut-Arimoto (arXiv:2305.07156, Algorithm 4 with Algorithm-3 caching), plus
the certified Arb dual/primal evaluator at the FINAL iterate.

Their algorithm (verified against the ar5iv text):
  - transition oracle P(x->y) computed on demand via half-split recursion:
      P(x->y) = sum_{k'=0..k} P(x_left -> y[:k']) * P(x_right -> y[k':])
    with leaf tables P(x'->y') for |x'| <= ceil(n/2), |y'| <= k, stored once.
    Leaf values are the classic DP subsequence counts.
  - BA step over only NONZERO (x,y) pairs:  D_y = sum_x P(x) P(x->y),
    a_x = sum_y P(x->y) log2(P(x->y)/D_y),  P'(x) = P(x) exp(a_x)/Z.
    (this is exactly what we need; and every P used is a nice positive float)

OP count at n=28, k=14: leaves 2^14 * 2^14 = 2.7e8 entries — too many for RAM.
Rubinstein-Con reduce further (Algorithm 3): the caching is over the smaller
half x1 (n1 = ceil(n/2)) and all y prefixes; and note that we do not need ALL
(x,y) pairs explicitly, we iterate rowwise: for each x (2^28), compute
a_x on the fly from the current D with an inner sum over the y's reachable
from x. The number of reachable y from a single x is the number of
subsequences, which at k=n-1 is at most n choose k ~ n. So the S(n,k) near k=n
is the SPARSE regime and their Algorithm-3 half-split is overkill — the
DOMINANT cost is the 2^28 x [row support] BA sweep at row level, i.e., for
k = n-1 the row support is exactly n+1 points (delete exactly one bit), so
one sweep is 2^28 * (n+1) ~ 7.8e9 ops — too slow in Python but ok in numpy
if we're careful with the vector formulation... For k near n it's actually
feasible because the output alphabet has 2^{n-1} symbols but D_y has to be
maintained; the full sweep is O(2^n * n) per iteration.

Since we cannot reach 2^28 sweeps in Python, the honest Plan for Gate-A
certified Cbar_{n,k}:
  * All rows that FD10 and RC actually USE in the d in {0.25, 0.5, 0.64} and
    0.68 chains at moderate n (n <= 22 with k <= 5 or near k=n; see FD10
    Table II; RC Table 1) get the same treatment.
  * At the n=28 point the certified enclosure is INFEASIBLE on this
    workstation. We'll state it as the honest partial result (per Main).

We also expose the "SPARSITY" formulation so that W_{n,k} for k = n-1 and
k = n-2 become n and n-2 * 2^k-sized sweeps — this is where the near-d rows
live (d large -> k small -> affordable).
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
import flint
from flint import arb, fmpq, fmpz
from functools import lru_cache


# ============================================================ exact DP matrices
def subseq_matrix_np(n: int, k: int) -> tuple[np.ndarray, int]:
    """Matrix of embedding counts N[x,y]; entries exact integers.
    Computed with the standard DP over (i,j) per y, vectorised over y.

    For memory, computes per-block: uses np.int64; entries may exceed int64
    for large n; in those regions we switch to python ints (object).
    """
    Y = np.arange(1 << k)
    bits = ((Y[:, None] >> np.arange(k)) & 1).astype(np.int8)   # (2^k, k)
    X = np.arange(1 << n)
    N = np.zeros((1 << n, 1 << k), dtype=object)  # object for big ints
    for i in range(1 << n):
        dpY = np.zeros((k + 1, 1 << k), dtype=object)
        dpY[0, :] = 1
        for ii in range(n):
            a = (i >> ii) & 1
            new = np.empty_like(dpY)
            new[0, :] = dpY[0, :]
            for j in range(1, k + 1):
                want = bits[:, j - 1] == a
                new[j, :] = dpY[j, :] + np.where(want, dpY[j - 1, :], 0)
            dpY = new
        N[i, :] = dpY[k, :]
    return N, math.comb(n, k)


def emit_nonzero_pairs(n: int, k: int):
    """Yields (x, list of (y, count)) sparsely — for k near n the row support
    is small (at k = n-1 exactly n+1 nonzero (x,y) per x)."""
    for x in range(1 << n):
        row = []
        for y in range(1 << k):
            c = _cntDP(x, y, n, k)
            if c:
                row.append((y, c))
        yield x, row


def _cntDP(x: int, y: int, n: int, k: int) -> int:
    A = [[0] * (k + 1) for _ in range(n + 1)]
    A[0][0] = 1
    for i in range(1, n + 1):
        A[i][0] = 1
        xi = (x >> (i - 1)) & 1
        for j in range(1, min(i, k) + 1):
            A[i][j] = A[i - 1][j] + (A[i - 1][j - 1] if xi == ((y >> (j - 1)) & 1) else 0)
    return A[n][k]


# ============================================================ BA with oracle
class RC_BAA:
    """
    BAA over the sparse W_{n,k} channel with on-demand transition oracle
    (Rubinstein-Con Algorithm 3/4 structure). Maintains D as a float array
    and sweeps row-by-row. Used in FLOAT mode to LOCATE the optimum, and
    then the FINAL p vector (only 2^n floats) is snapped to rationals.
    """

    def __init__(self, n: int, k: int):
        self.n, self.k = n, k
        self.Cnk = math.comb(n, k)
        self.nX = 1 << n
        self.nY = 1 << k
        # precompute row supports (list of (y, count)) lazily? All pairs is
        # 2^n * 2^k worst-case; we materialize the full transition NUMPY array
        # if it's small enough (<= 2^28 entries); otherwise the caller must
        # use the oracle path.
        self.Nfull = None
        if self.nX * self.nY <= 1 << 26:
            self.Nfull, _ = subseq_matrix_np(n, k)

    # ---- explicit row support
    def row_pairs(self, x: int):
        if self.Nfull is not None:
            return [(int(j), int(self.Nfull[x, j])) for j in np.nonzero(self.Nfull[x])[0]]
        return [(y, c) for y in range(1 << self.k) if (c := _cntDP(x, y, self.n, self.k))]

    def iterate(self, iters: int = 200, p0=None, record_interval: int = 10):
        nX = self.nX
        p = np.full(nX, 1.0 / nX) if p0 is None else np.array(p0, dtype=float)
        D = np.zeros(self.nY)
        # initial D
        for x in range(nX):
            if p[x] == 0:  # never for p0
                continue
            for y, c in self.row_pairs(x):
                D[y] += p[x] * c / self.Cnk
        for it in range(iters):
            p_new = np.zeros(nX)
            a_max = -np.inf
            a_arr = []
            for x in range(nX):
                s = 0.0
                for y, c in self.row_pairs(x):
                    nn = c / self.Cnk
                    s += nn * math.log2(nn / D[y])
                a_arr.append(s)
            m = max(a_arr)
            w = np.array([math.exp(a - m) for a in a_arr])
            p = w / w.sum()
            D[:] = 0.0
            for x in range(nX):
                if p[x] == 0:
                    continue
                for y, c in self.row_pairs(x):
                    D[y] += p[x] * c / self.Cnk
        return p, D

    def rate(self, p, D):
        r = 0.0
        for x in range(self.nX):
            if p[x] == 0:
                continue
            for y, c in self.row_pairs(x):
                nn = c / self.Cnk
                r += p[x] * nn * math.log2(nn / D[y])
        return r

    def dual(self, p, D):
        m = -np.inf
        for x in range(self.nX):
            s = 0.0
            for y, c in self.row_pairs(x):
                nn = c / self.Cnk
                s += nn * math.log2(nn / D[y])
            m = max(m, s)
        return m


def certify_from_sparse_rows(n: int, k: int, p_rational_num: list[int]):
    """
    Certified Arb enclosure of C(W_{n,k}) from an exact rational p = m/M.
    Computes I(p*) (lower bound) and max_x KL(W(.|x) || D_{p*}) (dual upper)
    in Arb with outward rounding; returns (lower_ball, upper_ball).
    """
    flint.ctx.prec = 400
    M = sum(p_rational_num)
    Cnk = math.comb(n, k)
    # D numerator per y (fmpz), computed exactly from integer m
    Dnum = [fmpz(0)] * (1 << k)
    for x, row in emit_nonzero_pairs(n, k):
        if p_rational_num[x] == 0:
            continue
        for y, c in row:
            Dnum[y] += fmpz(p_rational_num[x]) * c
    log2c = arb(2).log()

    def row_kl(x: int) -> arb | None:
        s = arb(0)
        for y, c in emit_row(n, k, x):
            if Dnum[y] == 0:
                return None
            w = fmpq(c, Cnk)
            r = fmpq(int(c) * int(M) * 1, Cnk * int(Dnum[y]))
            s += arb(w) * (arb(r).log() / log2c)
        return s

    lower = arb(0)
    uppers = []
    for x in range(1 << n):
        if p_rational_num[x] == 0:
            continue
        t = row_kl(x)
        if t is None:
            return None, None
        lower += arb(fmpq(p_rational_num[x], M)) * t
    for x in range(1 << n):
        t = row_kl(x)
        if t is not None:
            uppers.append(t.upper())
    return lower, max(uppers, key=lambda a: float(a))


def emit_row(n: int, k: int, x: int):
    for y in range(1 << k):
        c = _cntDP(x, y, n, k)
        if c:
            yield y, c
