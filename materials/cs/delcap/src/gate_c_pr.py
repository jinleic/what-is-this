"""
gate_c_pr.py — Pinto-Ribeiro (arXiv:2604.05867) C_{n,k} extension rows, k<=n<=9,
certified with the ANCHORED Gate C certificate path, plus a cross-check against
the frozen mpmath-locator rows from delcap_cert.

Why two routes. Pinto-Ribeiro publish C_{n,k} only at n in {29,31}, so these
rows have no published analogue and cannot be checked against a paper. An
unchecked single implementation is exactly what this repo does not accept, so
each row is produced twice by different machinery and the two are compared:

  route 1 (this script): float64 numpy BA locator -> exact rational snap ->
      gate_c_dhalf.cert_primal / cert_dual (Arb, 400 bits)
  route 2 (frozen):      mpmath 150-dps BA locator -> exact rational snap ->
      delcap_cert.certify_ab (Arb, 300 bits)

The two share only the exact integer channel matrix (subseq_matrix), which is
independently anchored by the frozen gate-B table. Different locator, different
snap denominator, different certificate code, different working precision.
Both intervals must contain the same true C_{n,k}, so the check is that they
OVERLAP; a disjoint pair would mean at least one certificate is wrong.
"""
from __future__ import annotations
import hashlib, json, os, sys, time, uuid

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT_ROOT = '/Users/jinleic/jinleic-workspace/cs/delcap/campaigns'
FROZEN_MPMATH = ('/Users/jinleic/jinleic-workspace/cs/delcap/campaigns/'
                 '2026-08-30T12:34:14Z_4cd10de4-efae-43f4-b9a6-39cab1db46b5_'
                 '8ecf82513989/pr_cnk_rows.jsonl')


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def main():
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        'gh', os.path.join(HERE, 'gate_c_dhalf.py'))
    gh = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(gh)
    from delcap_cert import subseq_matrix

    stamp = utc()
    code_hash = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:12]
    camp = os.path.join(OUT_ROOT, f'{stamp}_{uuid.uuid4()}_{code_hash}')
    os.makedirs(camp, exist_ok=True)
    logp = os.path.join(camp, 'log.txt')
    rowsf = os.path.join(camp, 'pr_rows_fastpath.jsonl')

    def LOG(msg):
        line = f'[{utc()}] {msg}'
        print(line, flush=True)
        with open(logp, 'a') as f:
            f.write(line + '\n')

    LOG(f'gate C Pinto-Ribeiro extension rows; campaign {camp}')
    LOG(f'code sha256[:12]={code_hash}')

    frozen = {}
    if os.path.isfile(FROZEN_MPMATH):
        for L in open(FROZEN_MPMATH):
            r = json.loads(L)
            frozen[(r['n'], r['k'])] = (r['lo'], r['hi'])
        LOG(f'loaded {len(frozen)} frozen mpmath-route rows for cross-check')

    done = []
    for n in range(6, 10):
        for k in range(1, n + 1):
            t0 = time.time()
            Ni = np.array(subseq_matrix(n, k)[0], dtype=object)
            Cnk = subseq_matrix(n, k)[1]
            nX, nY = Ni.shape
            Nf = np.array([[float(Ni[i, j]) / Cnk for j in range(nY)]
                           for i in range(nX)])
            p, D, _, _ = gh.ba_float(Nf, iters=4000)
            m_in = gh.snap(p, 1 << 30)
            m_out = gh.snap(D, 1 << 40)
            if min(m_out) == 0:
                m_out = [v + 1 for v in m_out]
            lo = gh.cert_primal(Ni, Cnk, m_in)
            hi_ba = gh.cert_dual(Ni, Cnk, m_out)
            hi_un = gh.cert_dual(Ni, Cnk, [1] * nY)
            hi = min([hi_ba, hi_un], key=lambda a: a.upper())
            L, H = float(lo.lower()), float(hi.upper())
            rec = dict(kind='PR_CNK_FASTPATH', n=n, k=k, cert_lo=L, cert_hi=H,
                       width=H - L, primal_ball=str(lo), dual_ball=str(hi),
                       primal_rad=float(lo.rad()), dual_rad=float(hi.rad()),
                       nX=nX, nY=nY, snap_in=1 << 30, snap_out=1 << 40,
                       wall_s=round(time.time() - t0, 1))
            fz = frozen.get((n, k))
            if fz is not None:
                flo, fhi = float(fz[0]), float(fz[1])
                overlap = not (H < flo or L > fhi)
                rec.update(mpmath_route_lo=flo, mpmath_route_hi=fhi,
                           routes_overlap=overlap,
                           tighter_by=(fhi - flo) - (H - L))
                verdict = ('TWO_ROUTES_OVERLAP' if overlap
                           else 'ROUTES_DISJOINT_INVESTIGATE')
            else:
                verdict = 'FASTPATH_ONLY_no_mpmath_row'
            rec['verdict'] = verdict
            done.append(rec)
            with open(rowsf, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
                fh.flush()
            extra = ''
            if fz is not None:
                extra = (f' | mpmath route [{rec["mpmath_route_lo"]:.10f}, '
                         f'{rec["mpmath_route_hi"]:.10f}] '
                         f'overlap={rec["routes_overlap"]}')
            LOG(f'  C_{n},{k} = [{L:.12f}, {H:.12f}] w={H - L:.2e} '
                f'({rec["wall_s"]}s){extra}')

    disj = [r for r in done if r['verdict'] == 'ROUTES_DISJOINT_INVESTIGATE']
    with open(os.path.join(camp, 'summary.json'), 'w') as f:
        json.dump(dict(stamp=stamp, code_sha256_12=code_hash, rows=done,
                       n_rows=len(done), n_crosschecked=sum(
                           1 for r in done if 'routes_overlap' in r),
                       n_disjoint=len(disj)), f, indent=1)
    LOG(f'done: {len(done)} rows, '
        f'{sum(1 for r in done if r.get("routes_overlap"))} cross-checked and '
        f'overlapping, {len(disj)} disjoint')
    print(camp)


if __name__ == '__main__':
    main()
