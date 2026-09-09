"""
gate_c_dhalf.py — Gate C headline rows: certified finite-blocklength enclosure of
the q-ary deletion-channel capacity C_{q,n}(d) at the ticket's target d = 1/2
(plus the paper's own d grid), compared against the Tavakoli-Nguyen-Bose
(arXiv:2607.19559) sandwich

    LB1  = (1-d) log2 q - h2(d)                                  (their Thm 1)
    LB+  = (1-d) log2 q + H_Bin(n,1-d)/n - h2(d) + Delta_n(d)/n   (their Cor 1)
    UB   = (1-d) log2 q                                          (their Thm 1)

Certified route (identical soundness argument to src/delcap_cert.py, which is the
frozen gate-A/B pipeline; the ONLY difference is that the BA *locating* step here
runs in float64 instead of mpmath-150dps — a locating step's rounding cannot
affect validity of either side, only tightness):

  primal:  I(p*) for an EXACT rational input distribution p* = m/M
           (candidates: uniform, and the float-BA iterate snapped to m/M).
           I(p) <= C for every p, so max over candidates is a valid LOWER bound.
  dual:    max_x KL(W(.|x) || D') for an EXACT rational output distribution D'
           (candidates: uniform, and the float-BA output marginal snapped).
           For ANY full-support D', C <= max_x KL(W(.|x)||D')
           (Csiszar-Tusnady / Blahut 1972), so min over candidates is a valid
           UPPER bound.
  Both evaluated ONCE in outward-rounded Arb from exact rationals: every log is
  a log of an exact rational, so there is no interval iteration anywhere.

Channel (exact, Tavakoli Lemma 1): for x in Sigma_q^n and y in Sigma_q^k,
    P((k,y)|x) = N_n(x,y) d^{n-k} (1-d)^k,
with N_n(x,y) the number of deletion subsets carrying x to y. Writing d = a/D in
lowest terms gives one common denominator D^n, so the whole transition matrix is
integer-valued: W[x,(k,y)] = N_n(x,y) a^{n-k} (D-a)^k over D^n.

Delta_n(d) = sum_{k=1}^{n-1} w_k Phi_{k,n}, Phi_{k,n} = (1/(q^n C(n,k))) sum_{x,y}
N log2 N — computed from the EXACT integer histogram of N_n values, in Arb.

Owner: DelcapGateC. Committed grid: see the gate-C section of pre_statement.md.
"""
from __future__ import annotations
import os
import itertools, json, math, os, sys, time, hashlib, platform, uuid

import numpy as np
from flint import arb, fmpq, fmpz, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

PREC = 400
OUT_ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'campaigns')


# ---------------------------------------------------------------- exact channel
def full_bdc_channel(q: int, n: int, d: fmpq):
    """Integer transition matrix of the length-n q-ary deletion channel.
    Returns (N int64 [q^n x sum_k q^k], Dden, outs) with
    P((k,y)|x) = N[x,(k,y)] / Dden exactly."""
    a = int(d.numerator)
    D = int(d.denominator)
    b = D - a                      # numerator of (1-d) over the same D
    dden = D ** n
    outs = [(k, y) for k in range(n + 1) for y in itertools.product(range(q), repeat=k)]
    col_of = {ky: i for i, ky in enumerate(outs)}
    nX, nY = q ** n, len(outs)
    N = np.zeros((nX, nY), dtype=object)
    for xi, x in enumerate(itertools.product(range(q), repeat=n)):
        dp = {(): 1}
        for si in x:
            ndp = dict(dp)                       # delete si
            for y, v in dp.items():              # keep si
                if len(y) < n:
                    y2 = y + (si,)
                    ndp[y2] = ndp.get(y2, 0) + v
            dp = ndp
        for y, v in dp.items():
            k = len(y)
            N[xi, col_of[(k, y)]] = v * (a ** (n - k)) * (b ** k)
    return N, dden, outs


def ba_float(Nf: np.ndarray, iters: int = 4000):
    """float64 BA locating step. Nf = row-stochastic float channel."""
    nX, nY = Nf.shape
    p = np.full(nX, 1.0 / nX)
    mask = Nf > 0
    logW = np.zeros_like(Nf)
    logW[mask] = np.log2(Nf[mask])
    for _ in range(iters):
        Dm = p @ Nf
        Dm = np.maximum(Dm, 1e-300)
        T = np.zeros_like(Nf)
        T[mask] = Nf[mask] * (logW[mask] - np.log2(np.broadcast_to(Dm, Nf.shape))[mask])
        a = T.sum(axis=1)
        w = np.exp2(a - a.max()) * p
        s = w.sum()
        if not np.isfinite(s) or s <= 0:
            break
        p = w / s
        p = np.maximum(p, 1e-24)
        p /= p.sum()
    Dm = np.maximum(p @ Nf, 1e-300)
    T = np.zeros_like(Nf)
    T[mask] = Nf[mask] * (logW[mask] - np.log2(np.broadcast_to(Dm, Nf.shape))[mask])
    rate = float(p @ T.sum(axis=1))
    dual = float(T.sum(axis=1).max())
    return p, Dm, rate, dual


def snap(vec, denom: int) -> list[int]:
    """largest-remainder snap of a float vector to integers summing to denom."""
    raw = [denom * float(v) for v in vec]
    m = [int(math.floor(r)) for r in raw]
    rem = denom - sum(m)
    order = sorted(range(len(raw)), key=lambda i: -(raw[i] - m[i]))
    i = 0
    while rem > 0:
        m[order[i % len(order)]] += 1
        rem -= 1
        i += 1
    assert sum(m) == denom and all(v >= 0 for v in m)
    return m


# ---------------------------------------------------------- certified evaluation
def cert_primal(N, Dden: int, m_in: list[int]) -> arb:
    """Outward-rounded I(p) with p = m_in/M exact. Valid LOWER bound on C."""
    ctx.prec = PREC
    log2 = arb(2).log()
    nX, nY = N.shape
    M = sum(m_in)
    Dnum = [0] * nY
    for i in range(nX):
        if m_in[i] == 0:
            continue
        row = N[i]
        for j in range(nY):
            c = row[j]
            if c:
                Dnum[j] += m_in[i] * int(c)
    I = arb(0)
    for i in range(nX):
        if m_in[i] == 0:
            continue
        pi = arb(fmpq(m_in[i], M))
        row = N[i]
        for j in range(nY):
            c = int(row[j])
            if not c:
                continue
            # W = c/Dden ; D'(y) = Dnum[j]/(M*Dden) ; ratio = M*c/Dnum[j]
            I += pi * arb(fmpq(c, Dden)) * (arb(fmpq(M * c, Dnum[j])).log() / log2)
    return I


def cert_dual(N, Dden: int, m_out: list[int]) -> arb:
    """Outward-rounded max_x KL(W(.|x)||D') with D' = m_out/S exact.
    Valid UPPER bound on C for ANY full-support D' (Csiszar-Tusnady)."""
    ctx.prec = PREC
    log2 = arb(2).log()
    nX, nY = N.shape
    S = sum(m_out)
    best = None
    for i in range(nX):
        row = N[i]
        s = arb(0)
        ok = True
        for j in range(nY):
            c = int(row[j])
            if not c:
                continue
            if m_out[j] == 0:
                ok = False       # D'(y)=0 where W(y|x)>0 -> KL = +infinity
                break
            # ratio = (c/Dden) / (m_out[j]/S) = c*S/(Dden*m_out[j])
            s += arb(fmpq(c, Dden)) * (arb(fmpq(c * S, Dden * m_out[j])).log() / log2)
        if not ok:
            return arb(float('inf'))
        if best is None or s.upper() > best.upper():
            best = s
    return best if best is not None else arb(0)


# ------------------------------------------------------- Tavakoli sandwich exact
def h2_arb(d: fmpq) -> arb:
    ctx.prec = PREC
    log2 = arb(2).log()
    one = fmpq(1)
    return -(arb(d) * (arb(d).log() / log2) + arb(one - d) * (arb(one - d).log() / log2))


def hbin_arb(n: int, p: fmpq) -> arb:
    ctx.prec = PREC
    log2 = arb(2).log()
    s = arb(0)
    for k in range(n + 1):
        w = fmpq(math.comb(n, k)) * p ** k * (fmpq(1) - p) ** (n - k)
        if w == 0:
            continue
        s += arb(w) * (arb(w).log() / log2)
    return -s


def pattern_hist(q: int, n: int, k: int) -> dict[int, int]:
    """Exact histogram {v: #(x,y) : N_n(x,y) = v} over x in Sigma_q^n, y in Sigma_q^k.
    Verified against Tavakoli Example 1 (q=2,n=2,k=1 -> {1:4, 2:2}) and Example 2
    (Phi_{1,2} = 1/q)."""
    total: dict[int, int] = {}
    for x in itertools.product(range(q), repeat=n):
        dp = {(): 1}
        for si in x:
            ndp = dict(dp)
            for y, v in dp.items():
                if len(y) < k:
                    y2 = y + (si,)
                    ndp[y2] = ndp.get(y2, 0) + v
            dp = ndp
        for y, v in dp.items():
            if len(y) == k:
                total[v] = total.get(v, 0) + 1
    return total


def phi_arb(q: int, n: int, k: int) -> arb:
    ctx.prec = PREC
    log2 = arb(2).log()
    s = arb(0)
    for v, cnt in pattern_hist(q, n, k).items():
        if v > 1:
            s += arb(cnt * v) * (arb(fmpz(v)).log() / log2)
    return s / (arb(fmpz(q) ** n) * arb(fmpz(math.comb(n, k))))


def delta_arb(q: int, n: int, d: fmpq) -> arb:
    ctx.prec = PREC
    s = arb(0)
    for k in range(1, n):
        w = fmpq(math.comb(n, k)) * d ** (n - k) * (fmpq(1) - d) ** k
        s += arb(w) * phi_arb(q, n, k)
    return s


def sandwich(q: int, n: int, d: fmpq) -> dict:
    ctx.prec = PREC
    log2q = arb(q).log() / arb(2).log()
    ub = log2q * arb(fmpq(1) - d)
    h2 = h2_arb(d)
    lb1 = ub - h2
    dl = delta_arb(q, n, d)
    lbp = ub + hbin_arb(n, fmpq(1) - d) / n - h2 + dl / n
    return dict(lb1=lb1, lbplus=lbp, ub=ub, delta=dl)


# ------------------------------------------------------------------- one row
def row(q: int, n: int, d: fmpq, iters: int = 4000, snap_in: int = 1 << 30,
        snap_out: int = 1 << 30) -> dict:
    t0 = time.time()
    N, Dden, outs = full_bdc_channel(q, n, d)
    nX, nY = N.shape
    Nf = np.array([[float(N[i, j]) / Dden for j in range(nY)] for i in range(nX)])
    p_ba, D_ba, rate_f, dual_f = ba_float(Nf, iters=iters)

    # --- primal candidates (exact rational input distributions)
    m_unif = [1] * nX
    prim_unif = cert_primal(N, Dden, m_unif)
    m_ba = snap(p_ba, snap_in)
    prim_ba = cert_primal(N, Dden, m_ba) if min(m_ba) >= 0 else None
    prim_cands = [('uniform', prim_unif)]
    if prim_ba is not None:
        prim_cands.append(('ba', prim_ba))
    prim_name, primal = max(prim_cands, key=lambda kv: kv[1].upper())

    # --- dual candidates (exact rational output distributions)
    m_out_unif = [1] * nY
    dual_unif = cert_dual(N, Dden, m_out_unif)
    m_out_ba = snap(D_ba, snap_out)
    # keep support: bump any zero entry that a row needs
    if min(m_out_ba) == 0:
        m_out_ba = [v + 1 for v in m_out_ba]
    dual_ba = cert_dual(N, Dden, m_out_ba)
    dual_name, dualv = min([('uniform', dual_unif), ('ba', dual_ba)],
                           key=lambda kv: kv[1].upper())

    sw = sandwich(q, n, d)
    lo_ps = float(primal.lower()) / n
    hi_ps = float(dualv.upper()) / n
    lb1_f, lbp_f, ub_f = float(sw['lb1']), float(sw['lbplus']), float(sw['ub'])
    # rigorous strictness margins (certified endpoints vs exact-Arb sandwich)
    beats_lb1 = float((primal / n).lower()) > float(sw['lb1'].upper())
    beats_lbp = float((primal / n).lower()) > float(sw['lbplus'].upper())
    beats_ub = float((dualv / n).upper()) < float(sw['ub'].lower())
    verdicts = []
    if beats_lbp:
        verdicts.append('CERT_LOWER_BEATS_LBplus')
    elif beats_lb1:
        verdicts.append('CERT_LOWER_BEATS_LB1_ONLY')
    if beats_ub:
        verdicts.append('CERT_UPPER_BEATS_UB')
    if not verdicts:
        verdicts.append('NO_STRICT_IMPROVEMENT')
    return dict(
        q=q, n=n, d=str(d),
        lb1=lb1_f, lbplus=lbp_f, ub=ub_f, delta_n=float(sw['delta']),
        cert_lo_block=float(primal.lower()), cert_hi_block=float(dualv.upper()),
        cert_lo_per_symbol=lo_ps, cert_hi_per_symbol=hi_ps,
        cert_width_per_symbol=hi_ps - lo_ps,
        primal_ball=str(primal), dual_ball=str(dualv),
        primal_rad=float(primal.rad()), dual_rad=float(dualv.rad()),
        primal_from=prim_name, dual_from=dual_name,
        ba_rate_float=rate_f, ba_dual_float=dual_f,
        snap_in=snap_in, snap_out=snap_out,
        nX=nX, nY=nY, verdict='+'.join(verdicts),
        wall_s=round(time.time() - t0, 1))


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

    LOG(f'gate C d=1/2 + paper-grid runner; campaign {camp}')
    LOG(f'code sha256[:12] = {code_hash}; python-flint prec = {PREC} bits')

    # self-tests against the papers' own printed examples (anchors for THIS code)
    h = pattern_hist(2, 2, 1)
    assert h == {1: 4, 2: 2}, h
    p12 = phi_arb(2, 2, 1)
    assert abs(float(p12) - 0.5) < 1e-30, p12
    d2 = delta_arb(2, 2, fmpq(1, 2))
    assert abs(float(d2) - 0.25) < 1e-30, d2   # their Example 3: Delta_2(d)=d(1-d)
    p12q3 = phi_arb(3, 2, 1)
    assert abs(float(p12q3) - 1.0 / 3) < 1e-30, p12q3   # their Thm 5: Phi_{1,2}=1/q
    LOG('self-test OK: pattern hist {1:4,2:2}, Phi_{1,2}(q=2)=1/2, Phi_{1,2}(q=3)=1/3, '
        'Delta_2(1/2)=1/4 (Tavakoli Examples 1-3, Thm 5)')

    grid = []
    # (i) headline: the ticket's d = 1/2 target
    for n in (2, 3, 4, 5, 6, 7, 8):
        grid.append((2, n, fmpq(1, 2)))
    for n in (2, 3, 4, 5):
        grid.append((3, n, fmpq(1, 2)))
    grid.append((4, 3, fmpq(1, 2)))
    grid.append((4, 4, fmpq(1, 2)))
    # (ii) the paper's own (q,n,d) points, for the recomputation table
    for q in (2, 3):
        for n in (3, 5, 10):
            for dd in (fmpq(1, 20), fmpq(1, 10), fmpq(1, 5)):
                if q == 3 and n == 10:
                    continue          # 3^10 x 88573 — out of this runner's budget
                grid.append((q, n, dd))

    only = os.environ.get('GATEC_ONLY')     # e.g. "2:10,3:3,3:5" -> q:n pairs
    if only:
        keep = {tuple(int(v) for v in tok.split(':')) for tok in only.split(',')}
        grid = [g for g in grid if (g[0], g[1]) in keep]
        LOG(f'GATEC_ONLY={only}: {len(grid)} rows selected')

    # freeze the scripts as run + tool versions BEFORE computing (write-first)
    import shutil, subprocess
    shutil.copy(__file__, os.path.join(camp, os.path.basename(__file__)))
    shutil.copy(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'pre_statement.md'),
                os.path.join(camp, 'pre_statement.md'))
    with open(os.path.join(camp, 'tool_versions.txt'), 'w') as f:
        import flint, numpy, mpmath, platform
        f.write(f'python {sys.version}\nplatform {platform.platform()}\n'
                f'python-flint {flint.__version__}\nnumpy {numpy.__version__}\n'
                f'mpmath {mpmath.__version__}\narb prec bits {PREC}\n')
    subprocess.run(f'cd {camp} && shasum -a 256 * > checksums.sha256',
                   shell=True, check=False)

    rowsf = os.path.join(camp, 'rows.jsonl')
    done = []
    for (q, n, d) in grid:
        try:
            r = row(q, n, d)
        except Exception as e:                      # noqa: BLE001
            LOG(f'ROW FAILED q={q} n={n} d={d}: {type(e).__name__}: {e}')
            continue
        done.append(r)
        with open(rowsf, 'a') as fh:
            fh.write(json.dumps(r) + '\n')
        LOG(f"q={q} n={n} d={d}: LB1={r['lb1']:.8f} LB+={r['lbplus']:.8f} "
            f"UB={r['ub']:.8f} | cert/sym=[{r['cert_lo_per_symbol']:.8f}, "
            f"{r['cert_hi_per_symbol']:.8f}] w={r['cert_width_per_symbol']:.3e} "
            f"(primal:{r['primal_from']}, dual:{r['dual_from']}) "
            f"{r['verdict']} [{r['wall_s']}s]")

    with open(os.path.join(camp, 'summary.json'), 'w') as f:
        json.dump(dict(stamp=stamp, code_sha256_12=code_hash, prec_bits=PREC,
                       rows=done), f, indent=1)
    LOG(f'done: {len(done)} rows -> {rowsf}')
    print(camp)


if __name__ == '__main__':
    main()
