"""
gate_c_md.py — Gate C, Morozov-Duman (arXiv:2504.20961) rows + Pinto-Ribeiro
(arXiv:2604.05867) certified C_{n,k} extension rows.

MORozov-DUMAN, layer-oriented converse bound (LO-CVB), their eqs (24)-(30):

    M(D_{mn}^{(delta)}, eps) <= L(n, m, eps, Lambda)
                             = (sum_{w in Lambda} tau_w)^n
                               / ((sum_{w in Lambda} p_w)^n - eps)
    tau_w = E(m,w) delta^{m-w} (1-delta)^w        (their eq 25)
    p_w   = C(m,w) delta^{m-w} (1-delta)^w        (their eq 26)
    E(m,w) = sum_{y in F_2^w} max_{x in F_2^m} binom(x|y)   (their eq 27)

and their reported quantity (Table III) is the CODE RATE
    r(n, eps) = log2 L(n,m,eps,Lambda*) / (n * m)     [bits per input bit]
minimised over Lambda subset of {0..m}. The n -> infinity row is the
max-oriented limit r(inf) = log2 tau(W) / m  (their eq 17), tau(W) = sum_w tau_w.

Everything here is EXACT: E(m,w) are integers, delta and eps are exact
rationals, so tau_w, p_w, L are exact rationals and the only inexact step is
one outward-rounded Arb evaluation of log2 of an exact rational.

E(m,w) provenance (labels are load-bearing):
  * m = 5: recomputed here from scratch by brute force over all (x,y)
    (their eq (28) prints E(5,2)=32, E(5,3)=52, E(5,4)=54)  -> MACHINE-VERIFIED
  * m = 20..23: typed from their Table I  -> CITED-DEPENDENCY, except the
    columns w in {0, 1, 2, m} which are re-derived here in closed form and
    checked against the typed values -> MACHINE-VERIFIED for those cells.
    (w=0: E=1; w=1: E=2m; w=m: E=2^m -- their eq (28).  w=2: derived here,
     E(m,2) = 2*C(m,2) + 2*floor(m^2/4), proof in `E_m2_closed`, brute-force
     validated at m=4..12 inside the self-test.)

Lambda search (rule 7 -- stated exactly):
  * m = 5: EXHAUSTIVE over all 2^6 = 64 subsets, then Arb-certified.
  * m = 22, 23: ALL 2^{m+1} subsets are enumerated in float64 (subset-sum
    doubling), ranked per n, and the top-K candidates are then certified in
    outward-rounded Arb; the reported number is the min over CERTIFIED
    candidates. Since EVERY Lambda yields a valid converse bound, a
    non-exhaustive certification can only make our reported bound LARGER
    (weaker), never invalid.

PINTO-RIBEIRO rows: their published C_{n,k} live at n in {29, 31} (GPU) and are
out of certified reach on this workstation; we instead certify C_{n,k} for
k <= n <= 9 with the frozen gate-A/B pipeline (delcap_cert.certify_capacity),
which extends the frozen gate-B table (n <= 7) to n = 8, 9.

Owner: DelcapGateC.
"""
from __future__ import annotations
import itertools, json, math, os, sys, time, hashlib, uuid

import numpy as np
from flint import arb, fmpq, fmpz, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PREC = 400
OUT_ROOT = '/Users/jinleic/jinleic-workspace/cs/delcap/campaigns'

# ---------------------------------------------------------------- E(m,w) table
# Typed from Morozov-Duman Table I (complete sets), columns m = 20, 21, 22, 23.
E_TABLE_PAPER = {
    20: [1, 40, 580, 5052, 30932, 142184, 514682, 1481532, 3671204, 7501642,
         12986826, 18226482, 24024636, 27877130, 27614704, 23235832, 17051216,
         11135474, 6263626, 2928320, 1048576],
    21: [1, 42, 640, 5894, 37994, 184954, 708084, 2216868, 5690200, 12575196,
         23446602, 35815610, 49223870, 62086746, 67882662, 63810994, 51413026,
         36780178, 23474060, 12917734, 5933988, 2097152],
    22: [1, 44, 704, 6804, 46148, 237180, 956052, 3191242, 8624830, 20507658,
         40757978, 69815062, 98366644, 130971724, 156701316, 162279170,
         145310300, 112905556, 79003896, 49317072, 26587726, 12015100, 4194304],
    23: [1, 46, 770, 7798, 55508, 299834, 1276366, 4504570, 12874990, 32366540,
         68846108, 126499426, 192944942, 266328578, 341932794, 387072452,
         381739786, 326856142, 246511952, 169050912, 103291092, 54626852,
         24310736, 8388608],
}

# Morozov-Duman Table III: code-rate LO-CVB at delta = 0.2, eps = 0.2.
# (printed to 5 decimals; 'inf' row is their n -> infinity limit)
TABLE_III = {
    5:  {1: 0.71688, 2: 0.69929, 4: 0.81882, 8: 0.81077, 16: 0.80675,
         32: 0.80473, 64: 0.80373, 128: 0.80323, 256: 0.80297, 512: 0.80285,
         1024: 0.80279, 'inf': 0.80272},
    22: {1: 0.55239, 2: 0.58012, 4: 0.61719, 8: 0.62095, 16: 0.65946,
         32: 0.66391, 64: 0.70137, 128: 0.70135, 256: 0.73575, 512: 0.73572,
         1024: 0.73571, 'inf': 0.73569},
    23: {1: 0.54775, 2: 0.57346, 4: 0.59406, 8: 0.62193, 16: 0.66192,
         32: 0.66262, 64: 0.70186, 128: 0.70179, 256: 0.70211, 512: 0.73417,
         1024: 0.73416, 'inf': 0.73414},
}


def embedding_count(x: tuple, y: tuple) -> int:
    n, k = len(x), len(y)
    dp = [1] + [0] * k
    for i in range(1, n + 1):
        for j in range(min(i, k), 0, -1):
            if x[i - 1] == y[j - 1]:
                dp[j] += dp[j - 1]
    return dp[k]


def E_brute(m: int, w: int) -> int:
    """E(m,w) by exhaustive maximisation over x for every y (m <= 12)."""
    tot = 0
    for y in itertools.product((0, 1), repeat=w):
        best = 0
        for x in itertools.product((0, 1), repeat=m):
            v = embedding_count(x, y)
            if v > best:
                best = v
        tot += best
    return tot


def E_m2_closed(m: int) -> int:
    """E(m,2) = 2*C(m,2) + 2*floor(m^2/4).

    Proof. For y = bb (two equal symbols) the count binom(x|y) is the number of
    ordered pairs i<j with x_i = x_j = b, maximised by the constant string:
    C(m,2). For y = bc with b != c the count is #{i<j : x_i=b, x_j=c}; writing
    the positions of b-symbols and c-symbols, the count is maximised by the
    block string b^t c^{m-t}, giving max_t t(m-t) = floor(m^2/4). Summing over
    the four y in F_2^2 (two equal, two unequal) gives the formula."""
    return 2 * math.comb(m, 2) + 2 * (m * m // 4)


def E_table(m: int) -> list[int]:
    """Complete E(m,·) list of length m+1 with provenance recorded by caller."""
    if m <= 12:
        return [E_brute(m, w) for w in range(m + 1)]
    if m in E_TABLE_PAPER:
        return list(E_TABLE_PAPER[m])
    raise KeyError(f'no E table for m={m}')


# -------------------------------------------------------------------- LO-CVB
def tau_p_exact(m: int, delta: fmpq, E: list[int]):
    """Exact rational tau_w and p_w for w = 0..m."""
    one = fmpq(1)
    tau = []
    pw = []
    for w in range(m + 1):
        wt = delta ** (m - w) * (one - delta) ** w
        tau.append(fmpq(E[w]) * wt)
        pw.append(fmpq(math.comb(m, w)) * wt)
    return tau, pw


def locvb_rate_arb(n, m: int, eps: fmpq, tau_sum: fmpq, p_sum: fmpq) -> arb | None:
    """Outward-rounded code rate log2 L /(n m) for one Lambda (given its exact
    tau_sum, p_sum).  n = 'inf' gives the max-oriented limit log2(tau_sum)/m.
    Returns None if the bound is invalid (denominator <= 0)."""
    ctx.prec = PREC
    log2 = arb(2).log()
    if n == 'inf':
        # exact-arithmetic guard: only the full-support Lambda is admissible
        if p_sum != fmpq(1):
            return None
        return (arb(tau_sum).log() / log2) / m
    den = p_sum ** n - eps
    if den <= 0:
        return None
    num_log = n * (arb(tau_sum).log() / log2)
    den_log = arb(den).log() / log2
    return (num_log - den_log) / (n * m)


def all_subset_sums(vals: list[float]) -> np.ndarray:
    """All 2^len(vals) subset sums, by doubling (float64)."""
    out = np.zeros(1, dtype=np.float64)
    for v in vals:
        out = np.concatenate([out, out + v])
    return out


def best_lambda_rows(m: int, delta: fmpq, eps: fmpq, E: list[int],
                     ns: list, topk: int = 24, exhaustive_max_items: int = 12):
    """For each n in ns return the certified min-rate over searched Lambdas."""
    tau, pw = tau_p_exact(m, delta, E)
    items = m + 1
    tau_f = [float(t) for t in tau]
    p_f = [float(t) for t in pw]
    if items <= exhaustive_max_items:
        masks = range(1, 1 << items)
        tau_sums = np.array([sum(tau_f[w] for w in range(items) if msk >> w & 1)
                             for msk in masks])
        p_sums = np.array([sum(p_f[w] for w in range(items) if msk >> w & 1)
                           for msk in masks])
        mask_list = np.array(list(masks))
        exhaustive = True
    else:
        tau_sums = all_subset_sums(tau_f)[1:]
        p_sums = all_subset_sums(p_f)[1:]
        mask_list = np.arange(1, 1 << items)
        exhaustive = False
    eps_f = float(eps)
    out = {}
    for n in ns:
        if n == 'inf':
            # ADMISSIBILITY: L(n,m,eps,Lambda) requires (sum_{Lambda} p_w)^n > eps.
            # As n -> infinity only Lambda with sum p_w = 1 (the full set) stays
            # admissible, and the bound tends to log2(tau(W))/m (their eq (17)).
            # Without this restriction the float ranking would happily pick
            # Lambda={0}, which is NOT a valid bound at any large n.
            score = np.where(p_sums >= 1.0 - 1e-12,
                             np.log2(np.maximum(tau_sums, 1e-300)) / m, np.inf)
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                den = np.power(p_sums, n) - eps_f
                score = np.where(den > 0,
                                 (n * np.log2(np.maximum(tau_sums, 1e-300))
                                  - np.log2(np.maximum(den, 1e-300))) / (n * m),
                                 np.inf)
        order = np.argsort(score)[:topk]
        certified = []
        for idx in order:
            if not np.isfinite(score[idx]):
                continue
            msk = int(mask_list[idx])
            Lam = [w for w in range(items) if msk >> w & 1]
            ts = sum((tau[w] for w in Lam), fmpq(0))
            ps = sum((pw[w] for w in Lam), fmpq(0))
            val = locvb_rate_arb(n, m, eps, ts, ps)
            if val is not None:
                certified.append((val, Lam))
        if not certified:
            out[n] = None
            continue
        val, Lam = min(certified, key=lambda kv: float(kv[0].upper()))
        out[n] = dict(rate_ball=str(val), rate_lo=float(val.lower()),
                      rate_hi=float(val.upper()), rad=float(val.rad()),
                      Lambda=Lam, lambda_search='exhaustive' if exhaustive
                      else f'all-{1 << items}-subsets-float-ranked+top{topk}-certified',
                      is_full_set=(len(Lam) == items))
    return out


# ------------------------------------------------------ greedy achievability
def greedy_avb(m: int, delta: fmpq, max_size: int = 32):
    """Morozov-Duman Alg. 1 (GAVB) for the m-bit deletion channel, EXACT
    rational arithmetic. Returns [(M, eps_M)] with eps_M = 1 - pi(C_M)/M."""
    a = int(delta.numerator)
    D = int(delta.denominator)
    b = D - a
    dden = D ** m
    outs = [(k, y) for k in range(m + 1) for y in itertools.product((0, 1), repeat=k)]
    col = {ky: i for i, ky in enumerate(outs)}
    W = [[0] * len(outs) for _ in range(1 << m)]
    for xi, x in enumerate(itertools.product((0, 1), repeat=m)):
        dp = {(): 1}
        for si in x:
            ndp = dict(dp)
            for y, v in dp.items():
                if len(y) < m:
                    y2 = y + (si,)
                    ndp[y2] = ndp.get(y2, 0) + v
            dp = ndp
        for y, v in dp.items():
            k = len(y)
            W[xi][col[(k, y)]] = v * (a ** (m - k)) * (b ** k)
    YP = [0] * len(outs)
    chosen = []
    frontier = []
    cand = set(range(1 << m))
    cur = 0                                   # pi(C) numerator over dden
    x0 = 0                                    # their all-zero start
    for step in range(1, min(max_size, 1 << m) + 1):
        if step == 1:
            pick = x0
        else:
            pick, gain = None, -1
            for x in cand:
                g = sum(max(W[x][j] - YP[j], 0) for j in range(len(outs)))
                if g > gain:
                    pick, gain = x, g
        cand.discard(pick)
        chosen.append(pick)
        for j in range(len(outs)):
            if W[pick][j] > YP[j]:
                cur += W[pick][j] - YP[j]
                YP[j] = W[pick][j]
        M = len(chosen)
        pc = fmpq(cur, dden * M)               # exact P(correct)
        frontier.append(dict(M=M, log2M=math.log2(M),
                             p_correct=str(pc), fer=str(fmpq(1) - pc),
                             fer_f=float(fmpq(1) - pc)))
    return frontier


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def main():
    stamp = utc()
    code_hash = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:12]
    camp = os.path.join(OUT_ROOT, f'{stamp}_{uuid.uuid4()}_{code_hash}')
    os.makedirs(camp, exist_ok=True)
    logp = os.path.join(camp, 'log.txt')

    def LOG(msg):
        line = f'[{utc()}] {msg}'
        print(line, flush=True)
        with open(logp, 'a') as f:
            f.write(line + '\n')

    LOG(f'gate C Morozov-Duman + Pinto-Ribeiro runner; campaign {camp}')
    LOG(f'code sha256[:12]={code_hash}; Arb prec={PREC} bits')

    # ---------- self-tests / anchors on E(m,w)
    checks = {}
    for (m, w, want, src) in [(5, 2, 32, 'MD eq (28)'), (5, 3, 52, 'MD eq (28)'),
                              (5, 4, 54, 'MD eq (28)')]:
        got = E_brute(m, w)
        checks[f'E({m},{w})'] = dict(brute=got, paper=want, ok=(got == want), src=src)
        LOG(f'  ANCHOR E({m},{w}) brute={got} paper={want} '
            f'{"MATCH" if got == want else "MISMATCH"} [{src}]')
    for m in range(4, 13):
        assert E_brute(m, 2) == E_m2_closed(m), (m, E_brute(m, 2), E_m2_closed(m))
    LOG('  E(m,2) closed form validated by brute force for m = 4..12')
    for m in (20, 21, 22, 23):
        row = E_TABLE_PAPER[m]
        ok0 = row[0] == 1
        ok1 = row[1] == 2 * m
        okm = row[m] == 2 ** m
        ok2 = row[2] == E_m2_closed(m)
        checks[f'paperTableI_m{m}'] = dict(w0=ok0, w1=ok1, w2=ok2, wm=okm,
                                           w2_closed=E_m2_closed(m), w2_typed=row[2])
        LOG(f'  Table-I m={m}: w=0 {ok0}, w=1 (=2m) {ok1}, '
            f'w=2 (closed form {E_m2_closed(m)}) {ok2}, w=m (=2^m) {okm}')
    with open(os.path.join(camp, 'E_checks.json'), 'w') as f:
        json.dump(checks, f, indent=1)

    # ---------- LO-CVB rows at their Table III points
    delta = fmpq(1, 5)
    eps = fmpq(1, 5)
    ns = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 'inf']
    rowsf = os.path.join(camp, 'locvb_rows.jsonl')
    escalations = []
    for m in (5, 22, 23):
        t0 = time.time()
        E = E_table(m)
        prov = ('brute-force MACHINE-VERIFIED' if m <= 12
                else 'typed from MD Table I (CITED-DEPENDENCY; w in {0,1,2,m} '
                     'independently re-derived)')
        res = best_lambda_rows(m, delta, eps, E, ns)
        LOG(f'  m={m}: Lambda search done in {time.time() - t0:.1f}s '
            f'(E provenance: {prov})')
        for n in ns:
            r = res[n]
            printed = TABLE_III[m][n if n != 'inf' else 'inf']
            if r is None:
                LOG(f'    n={n}: no valid Lambda (all denominators <= eps)')
                continue
            win_lo, win_hi = printed - 5e-6, printed + 5e-6
            if r['rate_hi'] < win_lo:
                verdict = 'CERT_STRICTLY_TIGHTER_THAN_PRINTED'
                escalations.append((m, n, printed, r['rate_lo'], r['rate_hi']))
            elif r['rate_lo'] > win_hi:
                verdict = 'CERT_WEAKER_THAN_PRINTED_(search/E-dependency)'
            else:
                verdict = 'REPRODUCES_PRINTED_5dp'
            rec = dict(kind='MD_LOCVB', m=m, n=(n if n != 'inf' else 'inf'),
                       delta='1/5', eps='1/5', printed=printed,
                       cert_lo=r['rate_lo'], cert_hi=r['rate_hi'],
                       width=r['rate_hi'] - r['rate_lo'], rad=r['rad'],
                       Lambda=r['Lambda'], is_full_set=r['is_full_set'],
                       lambda_search=r['lambda_search'],
                       E_provenance=prov, verdict=verdict)
            with open(rowsf, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
            LOG(f'    n={n}: printed={printed:.5f} cert=[{r["rate_lo"]:.8f}, '
                f'{r["rate_hi"]:.8f}] w={rec["width"]:.2e} '
                f'|Lambda|={len(r["Lambda"])} full={r["is_full_set"]} -> {verdict}')

    # ---------- extension: delta the paper only plots (Fig. 2), not tabulates
    ext = os.path.join(camp, 'locvb_extension_rows.jsonl')
    for m in (22, 23):
        E = E_table(m)
        for dd in (fmpq(1, 20), fmpq(1, 2), fmpq(4, 5)):
            res = best_lambda_rows(m, dd, eps, E, [1, 4, 64, 1024, 'inf'])
            for n, r in res.items():
                if r is None:
                    continue
                rec = dict(kind='MD_LOCVB_EXT', m=m, n=n, delta=str(dd),
                           eps='1/5', cert_lo=r['rate_lo'], cert_hi=r['rate_hi'],
                           width=r['rate_hi'] - r['rate_lo'],
                           Lambda=r['Lambda'], is_full_set=r['is_full_set'],
                           lambda_search=r['lambda_search'],
                           note='no printed analogue: MD show delta != 0.2 only '
                                'in Fig. 2 curves')
                with open(ext, 'a') as fh:
                    fh.write(json.dumps(rec) + '\n')
            LOG(f'  EXT m={m} delta={dd}: rows written')

    # ---------- greedy achievability (their Alg. 1) at m = 5
    t0 = time.time()
    fr = greedy_avb(5, fmpq(1, 5), max_size=32)
    with open(os.path.join(camp, 'gavb_m5.json'), 'w') as f:
        json.dump(fr, f, indent=1)
    LOG(f'  GAVB m=5 delta=1/5: {len(fr)} code sizes, exact rational FER; '
        f'M=4 -> FER={fr[3]["fer_f"]:.6f}, M=8 -> FER={fr[7]["fer_f"]:.6f} '
        f'({time.time() - t0:.1f}s)')

    # ---------- Pinto-Ribeiro: certified C_{n,k} rows, n <= 9
    from delcap_cert import certify_capacity
    prf = os.path.join(camp, 'pr_cnk_rows.jsonl')
    for n in range(6, 10):
        for k in range(1, n + 1):
            t0 = time.time()
            try:
                r = certify_capacity(n, k, iters=2000, denom=1 << 24)
            except Exception as e:                    # noqa: BLE001
                LOG(f'  PR C_{n},{k}: FAILED {type(e).__name__}: {e}')
                continue
            rec = dict(kind='PR_CNK', n=n, k=k, lo=r['lo_f'], hi=r['hi_f'],
                       width=r['hi_f'] - r['lo_f'], lo_ball=r['lo'],
                       hi_ball=r['hi'], iters=r['iters'], M=r['M'],
                       wall_s=round(time.time() - t0, 1))
            with open(prf, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
            LOG(f'  PR C_{n},{k} certified [{r["lo_f"]:.10f}, {r["hi_f"]:.10f}] '
                f'w={rec["width"]:.2e} ({rec["wall_s"]}s)')

    if escalations:
        with open(os.path.join(camp, 'ESCALATE.json'), 'w') as f:
            json.dump(escalations, f, indent=1)
        LOG(f'ESCALATION CANDIDATES (certified strictly tighter than printed): '
            f'{escalations}')
    LOG('gate C MD/PR runner complete.')
    print(camp)


if __name__ == '__main__':
    main()
