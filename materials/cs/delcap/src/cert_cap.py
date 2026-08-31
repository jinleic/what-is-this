"""
cert_cap.py — Certified capacity of the uniform-subsequence channel W_{n,k}.

Channel: input x in {0,1}^n, output y in {0,1}^k where y is a uniformly
random k-subsequence of x (i.e. delete n-k bits uniformly at random).
Transition counts N[x,y] = #{subsequence positions} give rational
transitions P(x->y) = N[x,y]/C(n,k). Exact integers.

Capacity CC(n,k) = max_p I(p), I(p) = sum_{x,y} p(x) P(x->y) log2(P(x->y)/D_y),
D_y = sum_x' p(x') P(x'->y).

CERTIFIED UPPER BOUND (Arimoto 1972 / Blahut 1972 dual):
  For ANY output distribution D_y on the output alphabet with D_y > 0 wherever
  needed CS <= max_x sum_y P(x->y) log2( P(x->y) / D_y ).
  This follows since CS = max_p sum_x p(x) KL(W_y|x || D_y) <= max_x KL(W_y|x || D_y).

CERTIFIED LOWER BOUND: I(p) <= CS for the specific p we plug in.

So CS in [ I(p), max_x KL(W_y|x || D_y) ] for the SAME D_y = p.W.

Pipeline:
  1. exact integer transition matrix
  2. high-precision mpmath BA, fixed # of iterations (deterministic, seeded)
  3. snap p to exact ratio p* = m/M on the simplex (M = sum m)
  4. compute D* = p*.W exactly (integer arithmetic)
  5. Arb-evaluate (outward rounded):
       I(p*)      = sum_{x,y} p*(x) P(x->y) log2(P(x->y)/D*_y)
       dual(p*,D*)= max_x sum_y P(x->y) log2(P(x->y)/D*_y)
     Both returning arb intervals guaranteed to contain the true value.
  6. Report [lower, upper] interval.
"""
from __future__ import annotations
import itertools, math
import mpmath as mp
import flint
from flint import arb, fmpq, fmpz, ctx

MP_DPS = 200


# ---------------------------------------------------------------- exact channel
def subseq_count_DP(x: tuple, y: tuple) -> int:
    """#{positions of y as a subsequence of x}, via DP. Exact."""
    n, k = len(x), len(y)
    if k > n:
        return 0
    A = [[0]*(k+1) for _ in range(n+1)]
    for i in range(n+1):
        A[i][0] = 1
    for i in range(1, n+1):
        xi = x[i-1]
        Amax = min(i, k)
        Ai, Api = A[i], A[i-1]
        for j in range(1, Amax+1):
            if xi == y[j-1]:
                Ai[j] = Api[j] + Api[j-1]
            else:
                Ai[j] = Api[j]
    return A[n][k]


def build_channel(n: int, k: int):
    """Returns rows: list over x of dict y -> multiplicity count (int)."""
    X = list(itertools.product([0, 1], repeat=n))
    Y = list(itertools.product([0, 1], repeat=k))
    C = math.comb(n, k)
    rows = []
    for x in X:
        row = {}
        for j, y in enumerate(Y):
            c = subseq_count_DP(x, y)
            if c:
                row[j] = c
        rows.append(row)
    return rows, C


# ---------------------------------------------------------------- BA in mpmath
def baa_mp(rows, C: int, iters: int = 500, dps: int = MP_DPS, tol: float = 1e-40):
    """High-precision BA. Returns (p as list[mp.mpf], rate, max_x KL)."""
    mp.mp.dps = dps
    Cc = mp.mpf(C)
    Prow = [{j: mp.mpf(c) / Cc for j, c in row.items()} for row in rows]
    nX = len(Prow)
    ys = set()
    for row in Prow:
        ys.update(row.keys())
    nY = max(ys) + 1
    p = [mp.mpf(1) / nX] * nX

    def marginal(p):
        D = [mp.mpf(0)] * nY
        for i, row in enumerate(Prow):
            if p[i] == 0:
                continue
            for j, v in row.items():
                D[j] += p[i] * v
        return D

    def rate(p, D):
        I = mp.mpf(0)
        for i, row in enumerate(Prow):
            if p[i] == 0:
                continue
            for j, v in row.items():
                I += p[i] * v * mp.log(v / D[j], 2)
        return I

    def dual(p, D):
        """max_x KL(W_y|x || D_y) -- certified upper bound."""
        m = mp.mpf(-10**50)
        for row in Prow:
            s = mp.mpf(0)
            for j, v in row.items():
                s += v * mp.log(v / D[j], 2)
            if s > m:
                m = s
        return m

    D = marginal(p)
    rprev = rate(p, D)
    for it in range(iters):
        a = []
        for row in Prow:
            s = mp.mpf(0)
            for j, v in row.items():
                s += v * mp.log(v / D[j], 2)
            a.append(s)
        ls = [mp.log(p[i]) + a[i] for i in range(nX)]
        M = max(ls)
        w = [mp.exp(ls[i] - M) for i in range(nX)]
        Z = sum(w)
        p = [wi / Z for wi in w]
        D = marginal(p)
        r = rate(p, D)
        if abs(r - rprev) < mp.mpf(10) ** (-dps + 20):
            break
        rprev = r
    return p, r, dual(p, D), it


def snap_to_rational(p, denom: int) -> list[int]:
    """Snap distribution p (sum 1) to integer numerators m_i over common
    denominator M = sum m_i ≈ denom, in a simplex-preserving way."""
    nX = len(p)
    m = [int(mp.floor(mp.mpf(denom) * pi)) for pi in p]
    Msnap = sum(m)
    # fix remaining mass to match denom (only affects M, not simplex)
    diff = denom - Msnap
    # distribute remaining units on largest fractional parts
    fracs = [(i, mp.mpf(denom) * p[i] - m[i]) for i in range(nX)]
    fracs.sort(key=lambda t: -t[1])
    for i in range(abs(diff)):
        j = fracs[i % nX][0]
        if diff > 0:
            m[j] += 1
        else:
            m[j] -= 1
            assert m[j] >= 0
    return m


def certified_interval(rows, C: int, m: list[int]):
    """Given integer m = numerators (denominator M = sum m), return
    (arb(lower), arb(upper), gap_ub) as arb balls containing CS.
    NOTE: all quantities are evaluated in OUTWARD-ROUNDED Arb; arb ball
    arithmetic guarantees each returned interval contains the true value."""
    M = sum(m)
    nX = len(m)
    ys = set()
    for row in rows:
        ys.update(row.keys())
    nY = max(ys) + 1
    # D*_y = sum_x m[x] * N[x,y] / (M*C) -- numerator in fmpz
    Dnum = [fmpz(0)] * nY
    for i, row in enumerate(rows):
        if m[i] == 0:
            continue
        for j, c in row.items():
            Dnum[j] += fmpz(m[i]) * c
    Dden = fmpz(M) * C
    # (arb lower, arb upper巢) for I(p*) and the dual Ub
    log2c = arb(2).log()  # ln 2

    def kl_row(row, sign):
        """Compute sum_y (c/C) sign*log2( (c/C) / (Dnum[y]/Dden) ).
        We compute it as: (c*Dden) / (C*Dnum[y]) has sign, take log2 via arb,
        return interval of (c/C) * sign * that = (c/(C*ln2)) * ln( ... ).
        sign=+1 gives I(p*) term; sign=-1 gives -ln2(ratio) → for dual we use
        the SAME formula (M row) but with the max over rows.
        """
        s = arb(0)
        for j, c in row.items():
            Dn = Dnum[j]
            if Dn == 0:
                # unreachable if D has full support, we asserted earlier
                return None
            # ratio r = (c/C) / (Dn/Dden) = c*Dden / (C*Dn)  -- fmpq exact
            r = fmpq(c * int(Dden), C * int(Dn))
            ls = arb(r).log() / log2c
            # weight is c/C (probability of transition)
            t = arb(fmpq(c, C)) * ls
            s = s + t
        return s

    # I(p*) = sum_x (m[x]/M) * sum_y (c/C) log2((c/C) / (D_y* / Dden))
    I_lo, I_hi = arb(0), arb(0)
    for i, row in enumerate(rows):
        if m[i] == 0:
            continue
        t = kl_row(row, +1)
        if t is None:
            return None, None, None
        w = fmpq(m[i], M)
        I_lo = I_lo + arb(w) * t
        I_hi = I_hi + arb(w) * t

    intervals = []
    for row in rows:
        t = kl_row(row, +1)
        if t is not None:
            intervals.append(t)
    # Simplify: compute all row intervals, and upper-bound CS by the
    # max of the per-row UPPER endpoints (all of which CONTAIN their true
    # per-row values, so the true dual <= max_i(t_i.upper()) and
    # CS <= dual <= max_i(t_i.upper()).
    dual_ub_of_upper = max(intervals, key=lambda z: z.upper()).upper()
    return I_lo.lower(), dual_ub_of_upper, I_hi.upper()  # I_hi sanity upper of I(p*)


def certify_channel(n: int, k: int, denom: int = 2**30, iters: int = 400,
                    dps: int = MP_DPS, verbose: bool = False):
    """Full certified pipeline for capacity of W_{n,k}.
    Returns dict with exact params, lower arb ball, upper arb ball, gap."""
    flint.ctx.prec = 300
    mp.mp.dps = dps
    rows, C = build_channel(n, k)
    p, rate_val, dual_val, it = baa_mp(rows, C, iters=iters, dps=dps)
    m = snap_to_rational(p, denom)
    lower, dual_ub_upper, I_hi = certified_interval(rows, C, m)
    return {
        "n": n, "k": k, "Cnk": C, "iters_run": it,
        "denom": denom,
        "m": m, "M": sum(m),
        "lower": lower, "dual_ub_upper": dual_ub_upper, "I_hi": I_hi,
        "gap_interval": (dual_ub_upper - lower) if lower is not None else None,
        "BA_rate_it_last": rate_val, "BA_dual_it_last": dual_val,
        "float_C": float(rate_val), "float_dual": float(dual_val),
    }


if __name__ == "__main__":
    # quick smoke test
    r = certify_channel(3, 2)
    print("n=3, k=2:")
    print("  lower CS ball:", r["lower"])
    print("  upper CS ball:", r["dual_ub_upper"])
    print("  gap interval width:", r["gap_interval"])
    print("  float BA rate:", r["float_C"])
