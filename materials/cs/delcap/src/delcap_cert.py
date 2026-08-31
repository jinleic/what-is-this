"""
delcap_cert.py — Certified enclosures for the uniform-subsequence channel
capacity Cbar_{n,k} (Fertonani-Duman's f(L,R)) and the derived deletion-channel
capacity bounds. Self-contained; python-flint Arb (outward-rounded balls) for
every certified number, mpmath only to LOCATE the BA optimum.

Certified scheme (fixed in cs/delcap/pre_statement.md before any run):
  CS = capacity of DMC = max_p I(p).
  For ANY full-support output distribution D on the output alphabet:
      I(p)                       <= CS  (primal: increasing in p)
      dual(D) := max_x KL(W(.|x) || D) >= CS
                   (Csizar-Tusnady / Blahut 1972: dual upper bound)
  Pipeline: high-precision BA (mpmath, 150 dps) -> snap p to an EXACT rational
  certificate m_x/M on the simplex (integer arithmetic) -> evaluate, in Arb:
      lower = I(p*)          (interval containing the true I(p*))
      upper = dual(D_{p*})   (interval containing the true dual >= CS)
  => CS in [ lower, upper ], a rigorous enclosure. No interval iteration:
  the ONLY interval arithmetic is a one-shot outward-rounded evaluation of
  explicit rational formulas, so nothing can diverge.
"""
from __future__ import annotations
import itertools, math
import numpy as np
import mpmath as mp
import flint
from flint import arb, fmpq, fmpz

MP_DPS = 150


def subseq_matrix(n: int, k: int) -> tuple[np.ndarray, int]:
    """Exact integer embedding-count matrix N (2^n x 2^k): N[x,y] = number of
    k-subsequence position-sets S of x with x|S = y. Row-normalised by C(n,k)
    gives the DMC W_{n,k}."""
    X = np.arange(1 << n, dtype=np.uint32)
    Y = np.arange(1 << k, dtype=np.uint32)
    N = np.zeros((1 << n, 1 << k), dtype=np.int64)
    for i in range(1 << n):
        x = i
        xb = ((x >> np.arange(n)) & 1).astype(np.int8)
        # DP over positions with running pattern match to every y at once
        # classic count-of-subsequences DP, vectorised over y.
        A = np.zeros((1 << k,), dtype=np.int64)
        A[0] = 1
        # full DP table A_prev[j] over y-prefixes: do (i x j) DP over floats of y
        # state: c[j] = #{ways y[0..j) is a subsequence of x[0..i)}
        c = np.zeros((1 << k,), dtype=np.int64)
        c[0] = 1
        # need per-y DP; vectorise by processing y bits column-wise
        # c is indexed by full y (k bits). Transition on consuming x bit a=x_i:
        # c'[y] = c[y] (skip) and for each y with bit j = a:
        #         c'[y] += c[y without bit j set... ] careful: standard DP is over
        # (i,j) per y. Vectorised alternative: iterate over y explicitly when
        # n,k small; for large n use the vectorised two-array trick below.
        yi = np.arange(1 << k)
        bits = ((yi[:, None] >> np.arange(k)) & 1).astype(np.int8)  # (2^k, k)
        # dp[i][j] for each y simultaneously: dpY (j+1 rows)
        dpY = np.zeros((k + 1, 1 << k), dtype=np.int64)
        dpY[0, :] = 1
        for ii in range(n):
            a = (i >> ii) & 1
            new = np.empty_like(dpY)
            new[0, :] = dpY[0, :]
            for j in range(1, k + 1):
                # skip x_ii: dpY[j] accumulated; match: need bit j-1 of y == a
                want = bits[:, j - 1] == a
                new[j, :] = dpY[j, :] + np.where(want, dpY[j - 1, :], 0)
            dpY = new
        N[i, :] = dpY[k, :]
    Cnk = math.comb(n, k)
    return N, Cnk


def baa_vectorised(N: np.ndarray, Cnk: int, iters: int = 4000,
                   dps: int = MP_DPS):
    """High-precision BA on channel N/Cnk. Returns (p, rate, dual, iters_run)."""
    mp.mp.dps = dps
    Cm = mp.mpf(Cnk)
    nX, nY = N.shape
    # sparse row lists of (y, prob)
    rows = []
    for i in range(nX):
        nz = np.nonzero(N[i])[0]
        rows.append([(int(j), mp.mpf(int(N[i, j])) / Cm) for j in nz])
    p = [mp.mpf(1) / nX] * nX

    def marg(p):
        D = [mp.mpf(0)] * nY
        for i, row in enumerate(rows):
            if p[i] == 0:
                continue
            for j, v in row:
                D[j] += p[i] * v
        return D

    def klrow(row, D):
        s = mp.mpf(0)
        for j, v in row:
            s += v * mp.log(v / D[j], 2)
        return s

    def rate(p, D):
        return sum(p[i] * klrow(rows[i], D) for i in range(nX) if p[i] != 0)

    def dual(p, D):
        return max(klrow(rows[i], D) for i in range(nX))

    D = marg(p)
    rprev = rate(p, D)
    it = 0
    for it in range(1, iters + 1):
        a = [klrow(rows[i], D) for i in range(nX)]
        ls = [mp.log(p[i]) + a[i] for i in range(nX)]
        M = max(ls)
        w = [mp.exp(ls[i] - M) for i in range(nX)]
        Z = sum(w)
        p = [wi / Z for wi in w]
        D = marg(p)
        r = rate(p, D)
        if abs(r - rprev) < mp.mpf(10) ** (-dps + 30):
            break
        rprev = r
    return p, r, dual(p, D), it


def snap_to_integers(p, denom: int) -> list[int]:
    """Snap distribution p to m/M with M=denom exactly (round + largest
    remainders). Returns integer list summing to denom."""
    import mpmath as mp
    raw = [mp.mpf(denom) * pi for pi in p]
    m = [int(mp.floor(r)) for r in raw]
    rem = denom - sum(m)
    fr = sorted(range(len(p)), key=lambda i: -(raw[i] - m[i]))
    for i in range(rem):
        m[fr[i]] += 1
    assert sum(m) == denom and all(v >= 0 for v in m)
    return m


def certify_ab(N: np.ndarray, Cnk: int, m: list[int]) -> tuple[arb, arb]:
    """Outward-rounded Arb evaluation of I(p*) and dual(D_{p*}) for the exact
    integer certificate m. Returns (lower_ball, upper_ball) with CS in between."""
    flint.ctx.prec = 400
    M = sum(m)
    nX, nY = N.shape
    log2c = arb(2).log()
    # D_y numerator: sum_x m[x] N[x,y]; exact in fmpz
    Dnum = (np.asarray(m, dtype=object).reshape(-1, 1) * N.astype(object)).sum(axis=0)
    Dden = M * Cnk

    def row_kl(i) -> arb | None:
        s = arb(0)
        for j in np.nonzero(N[i])[0]:
            j = int(j)
            Dn = int(Dnum[j])
            if Dn == 0:
                return None
            # P(x->y)=N/Cnk ; D_y = Dn/Dden ; ratio = N*Dden/(Cnk*Dn)
            r = fmpq(int(N[i, j]) * Dden, Cnk * Dn)
            s += arb(fmpq(int(N[i, j]), Cnk)) * (arb(r).log() / log2c)
        return s

    # primal I(p*)
    I = arb(0)
    for i in range(nX):
        if m[i] == 0:
            continue
        t = row_kl(i)
        assert t is not None
        I += arb(fmpq(m[i], M)) * t
    lower = I
    # dual: CS <= max_x KL(W_x || D) ; the true max lies in the union of the
    # per-row interval evaluations, so CS <= max_i (upper endpoint of row i).
    uppers = []
    for i in range(nX):
        t = row_kl(i)
        if t is not None:
            uppers.append(t.upper())
    upper = max(uppers, key=lambda a: float(a))
    return lower, upper


def certify_capacity(n: int, k: int, iters: int = 3000, denom: int = 1 << 24) -> dict:
    N, Cnk = subseq_matrix(n, k)
    p, rate, dual, it = baa_vectorised(N, Cnk, iters=iters)
    m = snap_to_integers(p, denom)
    lo, hi = certify_ab(N, Cnk, m)
    return dict(n=n, k=k, lo=str(lo), hi=str(hi), lo_f=float(lo), hi_f=float(hi),
                ba_rate=float(rate), ba_dual=float(dual), iters=it, M=denom)


# ----------------------------------------------------------------------
# Combination: U(n,d) = (1/n) sum_k binom(n,k) d^{n-k}(1-d)^k Cbar_{n,k}
# evaluated in outward-rounded Arb from per-k certified intervals.
def bino_arb(nn: int, kk: int) -> arb:
    return arb(math.comb(nn, kk))


def d_pow_arb(dqq: fmpq, e: int) -> arb:
    if e == 0:
        return arb(1)
    return arb(dqq) ** e


def U_bound(d: fmpq, caps: dict[int, tuple[arb, arb]], n: int) -> tuple[arb, arb]:
    """caps[k] = (lower, upper) arb balls for Cbar_{n,k}. Returns Arb interval
    for U(n,d) = (1/n) sum_k binom(n,k) d^{n-k}(1-d)^k Cbar_{n,k}."""
    lo_total = arb(0)
    hi_total = arb(0)
    one = arb(1)
    dlan = arb(d)
    pm = arb(1) - dlan
    for k, (clo, chi) in caps.items():
        w = bino_arb(n, k) * d_pow_arb(d, n - k) * d_pow_arb(1 - d, k) / n
        lo_total += w * clo
        hi_total += w * chi
    return lo_total, hi_total


if __name__ == "__main__":
    r = certify_capacity(3, 2, iters=200, denom=1 << 20)
    print(r)
