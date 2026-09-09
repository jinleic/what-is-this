"""
gate_c_run.py — Gate C: finite-length claims (Morozov-Duman 2504.20961,
Pinto-Ribeiro 2604.05867, Tavakoli-Nguyen-Bose 2607.19559) re-derived in
certified Arb arithmetic with the existing certified-BA pipeline
(delcap_cert.py scheme), extended to q-ary deletion channels.

Papers read first-hand from ar5iv full text today (owner reads).

Commitment: pre_statement.md (gate C section) fixed BEFORE this run.

Rows:
  (A) Tavakoli et al. sandwich values at their (q,n,d) grid (q in {2,3},
      n in {3,5,10}, d in {0.05,0.10,0.20} exact rationals) + extension to
      d=1/2 (the target), q=4, n in {7,8}.
      Outputs per row: LB1, LB+, UB (Arb balls), delta_n(d) exact + Arb,
      and the certified-BA interval on C_{q,n}(d) itself.
  (B) Pinto-Ribeiro: their published C_{n,k} values at n in {29,31} are OUT of
      certified reach on this workstation. We certify C_{n,k} intervals for
      n<=9 (all k), and report the paper-verdict for the rows that exist:
      the d=0.64 normalized comparison table (their Table 3) is quoted
      REPORTED, not recomputed. New certified rows at (n<=9, all k) are the
      added value.
  (C) Morozov-Duman: E(m,w) anchor rows from their eq (28) examples
      E(5,2)=32, E(5,3)=52, E(5,4)=54 recomputed exactly with an independent
      DP; the LO-CVB at (delta=0.2, eps=0.2) for the small-m base with the
      exact rational LO-CVB formula evaluated in Arb.

Evidence labels per cs/README.md: MACHINE-VERIFIED (certified Arb chains,
exact integers) / COMPUTATIONAL-EVIDENCE (float) / REPORTED (not opened).

Owner: DelcapGateC. Runs: `cs/.venv/bin/python src/gate_c_run.py`.
"""
from __future__ import annotations
import os
import json, math, time, hashlib, platform, sys, os, itertools

import numpy as np
import mpmath as mp
from flint import arb, fmpq, fmpz, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SRC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC)

# ----------------------------------------------------------------- helpers
MP_DPS = 150


def h2_arb(p: fmpq) -> arb:
    """outward-rounded -p log2 p - (1-p) log2(1-p); p in (0,1) exact rational."""
    ctx.prec = 400
    log2 = arb(2).log()
    lp = arb(p).log() / log2
    q1 = fmpq(1) - p
    lq = arb(q1).log() / log2
    return -(arb(p) * lp + arb(q1) * lq)


def hbin_arb(n: int, p: fmpq) -> arb:
    """outward-rounded H_Bin(n,p) = -sum_k C(n,k) p^k (1-p)^(n-k) log2(...)."""
    ctx.prec = 400
    log2 = arb(2).log()
    s = arb(0)
    for k in range(0, n + 1):
        w = fmpq(math.comb(n, k)) * p ** k * (fmpq(1) - p) ** (n - k)
        if w == 0:
            continue
        lw = (arb(w).log() / log2)
        s += arb(w) * lw
    return -s


def pattern_counts(q: int, n: int, k: int) -> dict[int, int]:
    """Exact histogram {v: #(x,y) with N_n(x,y)=v} for x in Sigma_q^n,
    y in Sigma_q^k. Verified against Tavakoli et al. Examples 1-2 (q=2,n=2,k=1:
    hist = {1:4, 2:2}, Phi_{1,2}=1/2) and against sum_y N(x,y)=C(n,k) row sums.

    DP per x: dp[j][y] = # deletion subsets S of [j] mapping x[0:j] to y.
      dp[j][y] = dp[j-1][y]  (delete x[j-1])
               + dp[j-1][y[:-1]]  if y[-1] == x[j-1]  (keep x[j-1])
    O(q^n * q^k * n) worst; cached per x over the y-dict only.
    """
    total: dict[int, int] = {}
    for x in itertools.product(range(q), repeat=n):
        dp = {(): 1}
        for si in x:
            ndp = {}
            for y, v in dp.items():          # delete si: y unchanged
                ndp[y] = ndp.get(y, 0) + v
            if k:                             # keep si: extends a length-(cur) y
                for y, v in dp.items():
                    if len(y) < k:
                        y2 = y + (si,)
                        ndp[y2] = ndp.get(y2, 0) + v
            dp = ndp
        for y, v in dp.items():
            if len(y) == k:
                total[v] = total.get(v, 0) + 1
    return total


def phi_kn_arb(q: int, n: int, k: int) -> arb:
    """Arb outward-rounded Phi_{k,n} = (1/(q^n C(n,k))) sum_{x,y} N log2 N."""
    ctx.prec = 400
    log2 = arb(2).log()
    hist = pattern_counts(q, n, k)
    s = arb(0)
    for v, cnt in sorted(hist.items()):
        if v == 0:
            continue
        s += (arb(cnt * v) * (arb(fmpz(v)).log() / log2))
    return s / (arb(fmpz(q) ** n) * arb(fmpz(math.comb(n, k))))


def delta_n_arb(q: int, n: int, d: fmpq) -> arb:
    """outward-rounded Delta_n(d) = sum_{k=1}^{n-1} w_k Phi_{k,n}."""
    ctx.prec = 300
    s = arb(0)
    for k in range(1, n):
        w = fmpq(math.comb(n, k)) * d ** (n - k) * (fmpq(1) - d) ** k
        phi = phi_kn_arb(q, n, k)
        s += arb(w) * phi
    return s




# ------------------------------------------------------- certified BA on q-ary
def full_bdc_channel(q: int, n: int, d: fmpq) -> tuple[np.ndarray, int, list[tuple[int,int]]]:
    """Exact integer transition matrix of the FULL q-ary deletion channel
    with deletion prob d (rational). Output space = {(k,y)} union over k,
    transition probability P((k,y)|x) = N_n(x,y) * d^(n-k) (1-d)^k
    (Tavakoli et al. Lemma 1). Writing d = aa/D, 1-d = bb/D coprime gives
    P = N*aa^(n-k)*bb^k / D^n: one exact rational denominator D^n.
    Returns (int64 matrix N-weighted, Dden=D^n, list of (k, y) outputs).
    """
    from fractions import Fraction
    fr = Fraction(int(d.numerator), int(d.denominator))
    # d = aa/bb in lowest terms; want d = aa'/D with 1-d = (D-aa')/D, D=bb.
    aa, bb = fr.numerator, fr.denominator
    D = fr.denominator
    dden = D ** n   # common denominator D^n for P((k,y)|x) = N*aa^(n-k)*(D-aa)^k
    outs = []
    for k in range(0, n + 1):
        for y in itertools.product(range(q), repeat=k):
            outs.append((k, y))
    nY = len(outs)
    nX = q ** n
    N = np.zeros((nX, nY), dtype=np.int64)
    # fill using the dp recursion per x (direct embedding count)
    for xi, x in enumerate(itertools.product(range(q), repeat=n)):
        dp = {(): 1}
        for si in x:
            ndp = {}
            for y, v in dp.items():      # delete si: y unchanged
                ndp[y] = ndp.get(y, 0) + v
            for y, v in dp.items():      # keep si: y extends by si
                if len(y) < n:
                    y2 = y + (si,)
                    ndp[y2] = ndp.get(y2, 0) + v
            dp = ndp
        for col, (k, y) in enumerate(outs):
            v = dp.get(y, 0)
            if v and k == len(y):
                # P((k,y)|x) = N_n(x,y) * d^(n-k) * (1-d)^k
                # = N * (aa^(n-k) * (D-aa)^k) / D^n — one common denom D^n=dden.
                N[xi, col] = v * (aa ** (n - k)) * ((D - aa) ** k)
    return N, dden, outs


def baa_full(N: np.ndarray, Dden: int, iters: int = 1500, dps: int = MP_DPS):
    """BA locating step (mpmath 150 dps) on the channel matrix N / Dden.
    Returns (p, D, rate, dual, its)."""
    mp.mp.dps = dps
    nX, nY = N.shape
    rows = []
    for i in range(nX):
        nz = np.nonzero(N[i])[0]
        rows.append([(int(j), mp.mpf(int(N[i, j])) / Dden) for j in nz])
    p = [mp.mpf(1) / nX] * nX

    def marg(pv):
        Dm = [mp.mpf(0)] * nY
        for i, row in enumerate(rows):
            if pv[i] == 0:
                continue
            for j, v in row:
                Dm[j] += pv[i] * v
        return Dm

    def klrow(row, Dm):
        s = mp.mpf(0)
        for j, v in row:
            s += v * mp.log(v / Dm[j], 2)
        return s

    def rate(pv, Dm):
        return sum(pv[i] * klrow(rows[i], Dm) for i in range(nX) if pv[i] != 0)

    def dual(pv, Dm):
        return max(klrow(rows[i], Dm) for i in range(nX))

    Dm = marg(p)
    rprev = rate(p, Dm)
    its = 0
    tiny = mp.mpf(10) ** (-(dps - 40))
    for its in range(1, iters + 1):
        aa_ = [klrow(rows[i], Dm) for i in range(nX)]
        Mv = max(aa_)
        w = [mp.exp(v - Mv) if v - Mv > -(dps - 40) else tiny for v in aa_]
        Z = sum(w)
        p = [max(wi / Z, tiny) for wi in w]
        Z2 = sum(p)
        p = [pi / Z2 for pi in p]
        Dm = marg(p)
        r = rate(p, Dm)
        if abs(r - rprev) < mp.mpf(10) ** (-dps + 30):
            break
        rprev = r
    return p, Dm, rate(p, Dm), dual(p, Dm), its


def snap_to_int(pv, denom: int) -> list[int]:
    """largest-remainder snap; exact sum = denom."""
    import mpmath as mp
    raw = [mp.mpf(denom) * x for x in pv]
    m = [int(mp.floor(r)) for r in raw]
    rem = denom - sum(m)
    order = sorted(range(len(pv)), key=lambda i: -(raw[i] - m[i]))
    for i in range(rem):
        m[order[i % len(order)]] += 1
    assert sum(m) == denom and all(v >= 0 for v in m)
    return m


def certify_full(N: np.ndarray, Dden: int, D_out: list[int]) -> arb:
    """Outward-rounded Arb dual certificate: CS <= max_x KL(W(.|x) || D')
    with D' = D_out/Dden exact rational output distribution (length == nY).
    Sound for ANY full-support D': needs Dn>0 on every channel-reachable y;
    asserted. Returns an Arb ball; .upper() >= true dual >= C_{q,n}(d)."""
    ctx.prec = 400
    log2 = arb(2).log()
    nX, nY = N.shape
    assert len(D_out) == nY, f'D_out len {len(D_out)} vs nY {nY}'
    Nobj = N.astype(object)
    Dnum = [sum(int(D_out[i]) * int(Nobj[i, j]) for i in range(nX)) for j in range(nY)]
    uppers = []
    for i in range(nX):
        nz = np.nonzero(Nobj[i])[0]
        if len(nz) == 0:
            continue
        s = arb(0)
        for j in nz:
            j = int(j)
            Dn = int(Dnum[j]); c = int(Nobj[i, j])
            if Dn == 0:
                # structural zero: D'=0 on a y that no x reaches (P(y|x)=0 too,
                # since D'-y column of N is zero everywhere). Skip: contributes
                # 0*log(0)=0 to both KL and the dual, standard convention.
                continue
            # KL(W(.|x) || D'):  P(y|x) = c/Dden ; D'(y) = Dn/(Dsum*Dden)
            # =>  ratio = c * Dsum / Dn  with Dsum = sum(D_out).
            Dsum = sum(int(v) for v in D_out)
            r = fmpq(int(c) * Dsum, Dn)
            s += arb(fmpq(c, Dden)) * (arb(r).log() / log2)
        uppers.append(s)
    return max(uppers, key=lambda a: a.upper()) if uppers else arb(0)


def certify_full_primal(N: np.ndarray, Dden: int, m_in: list[int],
                        only_rows: list[int] | None = None) -> arb:
    """Outward-rounded primal I(p*) for the exact snapped input distribution
    m_in/M. I(p*) = sum_x p(x) sum_y P(y|x) log2(P(y|x)/D'(y)) with exact
    rationals throughout: D'(y) = Dnum/(M*Dden); ratio W/D' = M*c/Dn.
    only_rows restricts to a subset of input rows (sub-case helper)."""
    ctx.prec = 400
    log2 = arb(2).log()
    nX, nY = N.shape
    M = sum(m_in)
    rows = only_rows if only_rows is not None else list(range(nX))
    Dnum = (np.asarray(m_in, dtype=object).reshape(-1, 1) * N.astype(object)).sum(axis=0)
    I = arb(0)
    Nobj = N.astype(object)
    for i in rows:
        if m_in[i] == 0:
            continue
        for j in np.nonzero(Nobj[i])[0]:
            j = int(j)
            c = int(Nobj[i, j])
            Dn = int(Dnum[j])
            r = fmpq(M * c, Dn)
            I += arb(fmpq(m_in[i], M)) * arb(fmpq(c, Dden)) * (arb(r).log() / log2)
    return I


    """Outward-rounded primal I(p*) for exact snapped input distribution
    m_in/M. I(p*) = sum_x p(x) sum_y P(y|x) log2(P(y|x)/D'(y)) with exact
    rationals throughout: D'(y) = Dnum/(M*Dden); ratio W/D' = M*c/Dn."""
    ctx.prec = 400
    log2 = arb(2).log()
    nX, nY = N.shape
    M = sum(m_in)
    Dsum_m = sum(m_in)
    Dnum = (np.asarray(m_in, dtype=object).reshape(-1, 1) * N.astype(object)).sum(axis=0)
    I = arb(0)
    Nobj = N.astype(object)
    for i in range(nX):
        if m_in[i] == 0:
            continue
        for j in np.nonzero(Nobj[i])[0]:
            j = int(j)
            c = int(Nobj[i, j])
            Dn = int(Dnum[j])
            # P(y|x)=c/Dden; D'(y)=Dn/(Dsum_m*Dden) => ratio = Dsum_m*(c/Dden)*Dden/Dn
            r = fmpq(Dsum_m * c, Dn)
            I += arb(fmpq(m_in[i], Dsum_m)) * arb(fmpq(c, Dden)) * (arb(r).log() / log2)
    return I


def certify_qbdc(q: int, n: int, d: fmpq, iters: int = 1500,
                 snap: int = 1 << 40) -> dict:
    """Full certified interval [primal, dual] on C_{q,n}(d).
    Dual: the certified upper bound is the MIN over two candidate exact
    rational output distributions: (i) snapped BA output marginal D_p*,
    (ii) uniform D' = (1/nY,...). Both are valid dual certificates (any
    full-support D' qualifies) — take the min for tightness."""
    t0 = time.time()
    N, dden, outs = full_bdc_channel(q, n, d)
    nY = N.shape[1]
    # candidate output distributions (exact rational): BA marginal + uniform.
    p, Dm, rate, dual, its = baa_full(N, dden, iters=iters)
    m_out_ba = snap_to_int(Dm, snap)
    m_unif = [snap // nY] * nY
    m_unif[-1] += snap - sum(m_unif)
    up_ba = certify_full(N, dden, m_out_ba)
    up_unif = certify_full(N, dden, m_unif)
    up = min([up_ba, up_unif], key=lambda a: a.upper())
    # candidate input distributions (exact rational): uniform + BA iterate.
    m_in_ba = snap_to_int(p, 1 << 24) if (q ** n) <= 4096 else None
    primals = []
    if (q ** n) <= 4096:
        m_in_unif = [1 << 21] * (q ** n)     # uniform input (exact)
        primals.append(certify_full_primal(N, dden, m_in_unif))
        primals.append(certify_full_primal(N, dden, m_in_ba))
        primal = max(primals, key=lambda a: a.upper())
    else:
        primal = certify_full_primal(N, dden, m_in_ba)
    return dict(q=q, n=n, d=str(d),
                primal=(str(primal) if primal else None),
                dual_upper=str(up), dual_upper_f=float(up),
                dual_ba_f=float(up_ba), dual_unif_f=float(up_unif),
                primal_f=(float(primal) if primal else None),
                width=(float(primal.upper()) - float(up.lower())
                       if primal is not None else None),
                ba_rate=float(rate), ba_dual=float(dual), its=its,
                snap=snap, wall_s=round(time.time() - t0, 1))


def sand_lower_lb1(q: int, d: fmpq) -> arb:
    return arb(q).log() / arb(2).log() * (fmpq(1) - d) - h2_arb(d) * arb(1)


def sand_ub(q: int, d: fmpq) -> arb:
    """(1-d) log2 q, outward-rounded — the Tavakoli trivial UB."""
    ctx.prec = 300
    return arb(q).log() / arb(2).log() * arb(fmpq(1) - d)


def sand_lbplus(q: int, n: int, d: fmpq) -> arb:
    """LB+ = (1-d)log2 q + H_Bin(n,1-d)/n - h2(d) + Delta_n(d)/n (Arb)."""
    ctx.prec = 400
    log2q = arb(q).log() / arb(2).log()
    ub_part = log2q * arb(fmpq(1) - d)
    hb = hbin_arb(n, fmpq(1) - d) / n
    h2d = h2_arb(d)
    dl = delta_n_arb(q, n, d) / n
    return ub_part + hb - h2d + dl


def verdict_of(lb1: arb, ub: arb, cint: dict, n: int) -> str:
    """Rule-7 verdict for a certified C_{q,n}(d) interval vs Tavakoli's sandwich.
    Sandwich per-symbol: [LB1, (1-d)log2 q]. A certified interval on the
    per-BLOCK capacity is divided by n; verdicts refer only to the swept set."""
    if cint['primal_f'] is None:
        return 'PRIMAL_UNAVAILABLE'
    lo_ps = cint['primal_f'] / n         # primal lower per symbol
    hi_ps = cint['dual_upper_f'] / n     # dual upper per symbol
    lb1_f = float(lb1)
    ub_f = float(ub)
    if lo_ps > lb1_f and hi_ps < ub_f:
        return 'INSIDE_SANDWICH_IMPROVES_BOTH'
    if lo_ps > lb1_f:
        return 'IMPROVES_LOWER_BOUND_ONLY'
    if hi_ps < ub_f:
        return 'IMPROVES_UPPER_BOUND_ONLY'
    if hi_ps > ub_f:
        return 'CERT_DUAL_LOOSER_THAN_PAPER_UB_validbutweak'
    return 'NO_IMPROVEMENT'



# ------------------------------------------------------- Morozov-Duman E(m,w)
def embedding_count(x: tuple, y: tuple) -> int:
    """exact binom(x|y) subsequence-count."""
    n, k = len(x), len(y)
    dp = [[0] * (k + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = 1
    for i in range(1, n + 1):
        for j in range(1, k + 1):
            dp[i][j] = dp[i - 1][j] + (dp[i - 1][j - 1] if x[i - 1] == y[j - 1] else 0)
    return dp[n][k]


def E_mw(m: int, w: int) -> int:
    """E(m,w) = sum_y max_x binom(x|y); brute force over x via DP not full
    enumeration — for m<=12 brute force with caching is fine."""
    tot = 0
    for y in itertools.product((0, 1), repeat=w):
        best = 0
        for x in itertools.product((0, 1), repeat=m):
            v = embedding_count(x, y)
            if v > best:
                best = v
        tot += best
    return tot


def locvb_arb(n: int, m: int, delta: fmpq, eps: fmpq, Lam: set[int],
              E_table: dict[int, int]) -> arb:
    """LO-CVB closed form: L(n,m,eps,Lambda) = (sum tau_w)^n / ((sum p_w)^n - eps)
    with tau_w = E(m,w) delta^(m-w)(1-delta)^w,
         p_w   = C(m,w) delta^(m-w)(1-delta)^w. Arb outward-rounded."""
    ctx.prec = 400
    tau = arb(0)
    pW = arb(0)
    for w in Lam:
        E = E_table[w]
        tau += arb(E) * (arb(delta) ** (m - w)) * (arb(fmpq(1) - delta) ** w)
        pW += (arb(math.comb(m, w)) * (arb(delta) ** (m - w)) * (arb(fmpq(1) - delta) ** w))
    num = tau ** n
    den = pW ** n - arb(eps)
    if den <= 0:
        return arb(float('inf'))
    return num / den


# ----------------------------------------------------------------- main run
def utc() -> str:
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def mac():
    return hashlib.sha256((platform.platform() + sys.version).encode()).hexdigest()[:12]


def main():
    stamp = utc()
    camp_dir = (os.path.join(os.path.dirname(SRC), "campaigns",
                             f"{stamp}_" + str(__import__('uuid').uuid4())
                             + f"_{mac()}"))
    os.makedirs(camp_dir, exist_ok=True)
    log_path = os.path.join(camp_dir, 'gate_c_log.txt')
    rows = []

    def LOG(msg):
        line = f'[{utc()}] {msg}'
        print(line, flush=True)
        with open(log_path, 'a') as f:
            f.write(line + '\n')

    LOG('gate C run begins; pre_statement.md was committed BEFORE this run')

    # ---------------- anchor: reproduce two frozen gate-A/B rows first
    from delcap_cert import certify_capacity
    anchors = [(3, 2), (5, 3)]
    for (an, ak) in anchors:
        r = certify_capacity(an, ak, iters=2500, denom=1 << 26)
        LOG(f"ANCHOR BDC_{an},{ak}: certified [{r['lo_f']:.10f}, {r['hi_f']:.10f}] "
            f"(frozen gate-B: f(3,2) in [1.4697819938, 1.4697820261], "
            f"f(5,3) in [1.8715442105, 1.8715442453])")
        rows.append(dict(gate='ANCHOR', n=an, k=ak, lo=r['lo_f'], hi=r['hi_f']))

    # ---------------- Tavakoli sandwich rows
    tandv = []
    for q in (2, 3):
        for n in (3, 5):
            for dr in (fmpq(1, 20), fmpq(1, 10), fmpq(1, 5)):
                tandv.append((q, n, dr))
    EXT = []
    for q in (2, 3):
        for n in (7, 8):
            EXT.append((q, n, fmpq(1, 2)))  # the ticket's d=1/2 target rows
    # q=4, d=1/2, n=4,5:
    EXT.append((4, 4, fmpq(1, 2)))
    EXT.append((4, 5, fmpq(1, 2)))

    LOG('=== (A) Tavakoli-Nguyen-Bose sandwich + certified C_{q,n}(d) ===')
    with open(os.path.join(camp_dir, 'tnb_rows.jsonl'), 'a') as fh:
        for (q, n, dr) in tandv + EXT:
            t0 = time.time()
            lb1 = sand_lower_lb1(q, dr)
            ub = sand_ub(q, dr)
            lbplus = sand_lbplus(q, n, dr)
            cint = certify_qbdc(q, n, dr, iters=1500, snap=1 << 40)
            rec = dict(kind='TNB', q=q, n=n, d=str(dr),
                       lb1=str(lb1), lb1_f=float(lb1),
                       lbplus=str(lbplus), lbplus_f=float(lbplus),
                       ub=str(ub), ub_f=float(ub),
                       c_lo=cint['primal_f'], c_hi=cint['dual_upper_f'],
                       c_lo_per_symbol=(cint['primal_f']/n if cint['primal_f'] else None),
                       c_hi_per_symbol=cint['dual_upper_f']/n,
                       c_width=cint['width'],
                       verdict=verdict_of(lb1, ub, cint, n),
                       ba_its=cint['its'], wall_s=cint['wall_s'])
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
            LOG(f"TNB q={q} n={n} d={dr}: LB1={rec['lb1_f']:.8f} "
                f"LB+={rec['lbplus_f']:.8f} UB={rec['ub_f']:.8f} "
                f"certPerSy=[{(rec['c_lo_per_symbol'] or 0):.8f}, "
                f"{rec['c_hi_per_symbol']:.8f}] width={cint['width'] if cint['width'] is None else format(cint['width'],'.2e')} "
                f"verdict={rec['verdict']} ({cint['wall_s']}s)")
    LOG('TNB loop done')

    # d=1/2 sandwich-improvement measurement (target deliverable)
    LOG('=== d=1/2 certified sandwich improvement test ===')
    impr = []
    for (q, n) in [(2, 8), (3, 8), (4, 8), (3, 10)]:
        try:
            rec_impr = dict(q=q, n=n)
            dl = delta_n_arb(q, n, fmpq(1, 2))
            rec_impr['delta_n'] = str(dl)
            rec_impr['delta_f'] = float(dl)
            impr.append(rec_impr)
            LOG(f'  Delta_{n}(1/2) at q={q}: {rec_impr["delta_f"]:.10f}')
        except Exception as e:
            LOG(f'  Delta_{n}(1/2) at q={q}: FAILED {e}')

    # ---------------- Pinto-Ribeiro: certified C_{n,k} rows, n<=9, all k
    LOG('=== (B) Pinto-Ribeiro C_{n,k} certified rows (their n<=9 NEW, '+
        'n=29/31 out of certified reach — documented, not attempted) ===')
    with open(os.path.join(camp_dir, 'pr_rows.jsonl'), 'a') as fh:
        for n in range(1, 10):
            for k in range(1, n + 1):
                r = certify_capacity(n, k, iters=2500, denom=1 << 24)
                # trivial dual Cbar <= k sanity
                rec = dict(kind='PR', n=n, k=k, lo=r['lo'], hi=r['hi'],
                           lo_f=r['lo_f'], hi_f=r['hi_f'],
                           width=r['hi_f'] - r['lo_f'])
                fh.write(json.dumps(rec) + '\n')
                fh.flush()
            LOG(f'PR n={n}: done all k')

    # ---------------- Morozov-Duman anchor rows
    LOG('=== (C) Morozov-Duman E(m,w) anchor recomputation ===')
    md = {}
    for (m, w) in [(5, 2), (5, 3), (5, 4), (10, 4), (10, 5), (20, 8)]:
        t0 = time.time()
        if m <= 12:
            E = E_mw(m, w)
            md[f'E({m},{w})'] = dict(m=m, w=w, value=E, computed='exact',
                                     wall_s=round(time.time() - t0, 1))
        else:
            md[f'E({m},{w})'] = dict(m=m, w=w, value=None,
                                     computed='out-of-scope (m>12 brute force)',
                                     wall_s=round(time.time() - t0, 1))
        LOG(f'  E({m},{w}) = {md[f"E({m},{w})"]["value"]} ({md[f"E({m},{w})"]["wall_s"]}s)')
    with open(os.path.join(camp_dir, 'md_E_anchors.json'), 'w') as f:
        json.dump(md, f, indent=1)

    # LO-CVB at (delta=0.2, eps=0.2), their Table III set
    LOG('=== (C2) LO-CVB at delta=0.2, eps=0.2 (Morozov-Duman Table III rows) ===')
    delta = fmpq(1, 5)
    eps = fmpq(1, 5)
    with open(os.path.join(camp_dir, 'lo_cvb_rows.jsonl'), 'a') as fh:
        for m in (5, 10, 23):
            E_table = {}
            for w in range(m + 1):
                E_table[w] = E_mw(m, w) if m <= 10 else None
            if any(v is None for v in E_table.values()):
                LOG(f'  LO-CVB m={m}: E table incomplete (m>10); skipping')
                continue
            best = None
            bestLam = None
            for r_ in range(m + 1):
                for Lam in itertools.combinations(range(m + 1), r_):
                    if not Lam:
                        continue
                    val = locvb_arb(1, m, delta, eps, set(Lam), E_table)
                    if best is None or val < best:
                        best = val
                        bestLam = Lam
            rec = dict(m=m, delta='0.2', eps='0.2',
                       locvb=str(best), locvb_f=float(best),
                       best_Lambda=list(bestLam or []) if bestLam else None)
            fh.write(json.dumps(rec) + '\n')
            LOG(f'  m={m}: LO-CVB = {rec["locvb_f"]:.6f} at Lambda={rec["best_Lambda"]}')

    LOG('gate C run complete.')
    with open(os.path.join(camp_dir, 'summary.json'), 'w') as f:
        json.dump(rows, f, indent=1)
    LOG(f'campaign dir: {camp_dir}')


if __name__ == '__main__':
    main()
