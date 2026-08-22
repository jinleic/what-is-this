#!/usr/bin/env python3
"""Wave-3 H-scan: ALTERNATE aligned classes per cell + deep-k rescan.

H needs only ONE aligned class with ONE emergent-free member per cell.
Wave 1-2 scanned the frozen classes (_L11_CLASSES/_L12_ESCAPE) to k<=200.
This driver attacks the open cells two ways:

L11 branch: alternate a values.  For the class to be aligned at the wall
prime w (w | numer(z)), the tower needs w | A = 1+4a^2 (L11h); candidates
are a == r (mod w) with w | 1+4r^2, lifted a = r + w*t or (w-r) + w*t.
For each (a, eps, q1) with an ok _l10_class_cert, scan k <= per_a_k.

ESC branch: alternate f | numer(z) primes and eps flips at a = 1,
certificate mirrors l12_class.l12_class_cert's guard list.  Cheap
self-validation: smooth_emergent recomputes ALL alignment symbols per
member, so a wrongly-built class can only yield 'align-fail' members
(never a false zero); every zero is rechecked by ramified() at freeze.

DEEPK branch: rescan frozen classes of open cells with max-k deep.

Usage:
  python3 l13h_alt.py L11 LO HI OPENJSONLS OUT.jsonl [opts]
  python3 l13h_alt.py ESC LO HI OPENJSONLS OUT.jsonl [opts]
  python3 l13h_alt.py DEEPK L11|ESC OPENJSONLS OUT.jsonl [--max-k 600]
OPENJSONLS: comma-separated wave-2 jsonl files to read open classes from.
"""
import json
import math
import sys
import time
from fractions import Fraction as F

from h10q import (_L11_CLASSES, _L12_ESCAPE, _is_prime, _l10_class_cert,
                  _l10_supp, _l7_tied_status, primerange, ramified,
                  PrimalityBound)
from l13_filter import smooth_emergent, cofactor_decide
from l13h_scan import excluded_primes_L11, class_data

V_NONE = ('align-fail', 'none', 'bad', 'smooth')


def open_indices(jsonls, fam):
    """Open class indices of `fam` across scan-output jsonls: every class
    index seen as a class record of `fam`, minus every CELL a closure
    record exists for (wave 1/2 'status': 'closed' rows and wave-3 alt
    records, which carry no status).  Cell keys are (w, (u1, u2))."""
    keys_l = sorted(_L11_CLASSES) if fam == 'L11' else sorted(_L12_ESCAPE)
    closed_cells = set()
    seen = []
    for path in jsonls.split(','):
        try:
            fh = open(path)
        except FileNotFoundError:
            continue
        for line in fh:
            r = json.loads(line)
            if r.get('type') != 'class':
                continue
            cell = (r['cell'][0], tuple(r['cell'][1]))
            if r.get('status') == 'closed' or 'alt' in str(r.get('family')):
                closed_cells.add(cell)
            elif r.get('family') == fam and r.get('status') == 'open':
                seen.append(r['class_idx'])
    return sorted(idx for idx in set(seen)
                  if (keys_l[idx][0], keys_l[idx][1]) not in closed_cells)


def try_member(fam, key, a, eps, f, N, tau, z, excl, k, Q, deadline):
    if Q in excl:
        return ('excluded', None)
    try:
        if not _is_prime(Q):
            return ('nonprime', None)
    except PrimalityBound:
        return ('nonprime', None)
    b = F(eps*f*Q)
    st, dt = smooth_emergent(a, z, tau, b)
    if st == 'zero':
        return ('zero', (k, Q, 'square', b))
    if st == 'cofactor-big':
        v, info = cofactor_decide(a, z, tau, b, detail=dt, deadline=deadline)
        if v == 'zero':
            return ('zero', (k, Q, info.get('rung'), b))
        return (v, None)
    return (st, None)


def close_cell(fam, key, k, Q, rung, b, a, z, tau):
    tied = None
    try:
        tied = str(_l7_tied_status(a, b, z, tau))
    except Exception as e:
        tied = f"refused:{type(e).__name__}"
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    alpha = -delta*A
    c0 = None
    from h10q import _sun_h
    c0 = _sun_h(a, b, z**3)
    M0 = 16 - delta*c0*c0 - 32*A*b*((a-1)/2)**2
    x0, d0 = alpha*M0, alpha*2*b
    ram_checked = None
    try:
        ram_checked = (ramified(x0, d0) == [])
    except Exception as e:
        ram_checked = f"refused:{type(e).__name__}"
    return {'k_zero': k, 'Q': str(Q), 'rung': rung, 'tied': tied,
            'ramified_empty': ram_checked}


def scan_L11_alt(key, out, o, T_MAX=40, Q1_MAX=4000, PER_A_K=200,
                 ALT_CAP=200, deadline=None):
    w, ut = key
    z = F(w)*F(*ut)
    roots = [r for r in range(1, w) if (1 + 4*r*r) % w == 0]
    a_pool = sorted({abs(r + s*w*t) for r in roots for s in (1, -1)
                     for t in range(T_MAX + 1)} )
    a_pool = [a for a in a_pool if a % 2 == 1]
    if not a_pool:
        a_pool = [1]      # w == 3 (mod 4): no tower a; the a=1 class path
    tried = 0
    for a_ in a_pool:
        a = F(a_)
        A = 1 + 4*a*a
        tau = (1 + 2*a*a)/A
        excl = excluded_primes_L11(a, z, tau)
        for eps in (1, -1):
            for q1 in primerange(41, Q1_MAX + 1):
                if tried >= ALT_CAP:
                    break
                cert = _l10_class_cert(a, z, tau, eps, q1)
                if not cert or not cert.get('ok'):
                    continue
                tried += 1
                N = cert['N']
                o.write(json.dumps({'type': 'alt-class', 'cell': list(key),
                                    'a': a_, 'eps': eps, 'q1': q1,
                                    'N': str(N)}) + '\n')
                for k in range(PER_A_K + 1):
                    if deadline and time.time() > deadline:
                        return ('deadline', tried)
                    Q = q1 + k*N
                    v, hit = try_member('L11', key, a, eps, 1, N, tau, z,
                                        excl, k, Q, deadline)
                    if v == 'zero':
                        rec = close_cell('L11', key, *hit[:3], hit[3],
                                         a, z, tau)
                        rec.update({'type': 'class', 'family': 'L11-alt',
                                    'cell': list(key), 'a': str(a),
                                    'eps': eps, 'q1': q1, 'N': str(N)})
                        o.write(json.dumps(rec) + '\n')
                        o.flush()
                        print(f"*** ALT-CLOSED L11 {key} a={a_} eps={eps}"
                              f" q1={q1} k={hit[0]} rung={hit[2]}",
                              flush=True)
                        return ('closed', tried)
                if tried >= ALT_CAP:
                    break
            if tried >= ALT_CAP:
                break
    return ('open', tried)


def scan_ESC_alt(key, out, o, Q1_MAX=8000, PER_A_K=200, ALT_CAP=400,
                 deadline=None):
    w, ut = key
    z = F(w)*F(*ut)
    a = F(1)
    A = 5
    tau = F(3, 5)
    delta = F(-4, 5)
    alpha = F(4)
    Z = z**3
    D = 1 - Z - Z*Z
    from h10q import factorint, vp
    f_frozen, eps_frozen, _ = _L12_ESCAPE[key]
    f_cands = sorted(p for p in factorint(z.numerator) if p != 2)
    variants = []
    for f_ in f_cands:
        for eps_ in (1, -1):
            if f_ == f_frozen and eps_ == eps_frozen:
                continue            # the wave-1/2 class itself
            variants.append((f_, eps_))
    excl = set()
    for p in factorint(D.numerator):
        excl.add(p)
    for p in factorint(z.numerator):
        excl.add(p)
    tried = 0
    for f_, eps_ in variants:
        if any(vp(t, f_) != 0 for t in (delta, A, a, D)):
            continue                 # mirror the l12 guard: f shares ONLY z
        S = sorted({2, 3, 5, 7, f_})
        from h10q import _l10_P, _l10_exponent
        P = _l10_P(a, Z, D, A, delta, (a - 1)/2)
        ks = {p: _l10_exponent(P, F(eps_*f_), p) for p in S}
        N = 8
        for p, k_ in ks.items():
            N = N * p**k_ // math.gcd(N, p**k_)
        N = N * 4*A // math.gcd(N, 4*A)
        for q1 in primerange(41, Q1_MAX + 1):
            if tried >= ALT_CAP:
                break
            if q1 in S or any(vp(t, q1) != 0 for t in (delta, A, a)):
                continue
            b0 = F(eps_*f_*q1)
            from h10q import hilbert, OO, _sun_h
            M00 = 16 - delta*_sun_h(a, b0, Z)**2 - 32*A*b0*((a-1)/2)**2
            x00, d00 = alpha*M00, alpha*2*b0
            if any(hilbert(x00, d00, p) != 1 for p in S):
                continue
            if hilbert(x00, d00, OO) != 1:
                continue
            tried += 1
            o.write(json.dumps({'type': 'alt-class', 'cell': list(key),
                                'a': 1, 'f': f_, 'eps': eps_, 'q1': q1,
                                'N': str(N)}) + '\n')
            for k in range(PER_A_K + 1):
                if deadline and time.time() > deadline:
                    return ('deadline', tried)
                Q = q1 + k*N
                v, hit = try_member('ESC', key, a, eps_, f_, N, tau, z,
                                    excl, k, Q, deadline)
                if v == 'zero':
                    rec = close_cell('ESC', key, *hit[:3], hit[3], a, z,
                                     tau)
                    rec.update({'type': 'class', 'family': 'ESC-alt',
                                'cell': list(key), 'a': '1', 'f': f_,
                                'eps': eps_, 'q1': q1, 'N': str(N)})
                    o.write(json.dumps(rec) + '\n')
                    o.flush()
                    print(f"*** ALT-CLOSED ESC {key} f={f_} eps={eps_}"
                          f" q1={q1} k={hit[0]} rung={hit[2]}", flush=True)
                    return ('closed', tried)
        if tried >= ALT_CAP:
            break
    return ('open', tried)


def scan_DEEPK(fam, key, out, o, max_k, deadline):
    a, eps, f, q1, N, tau, z, excl = class_data(fam, key)
    for k in range(201, max_k + 1):
        if deadline and time.time() > deadline:
            return ('deadline', k)
        Q = q1 + k*N
        v, hit = try_member(fam, key, a, eps, f, N, tau, z, excl, k, Q,
                            deadline)
        if v == 'zero':
            rec = close_cell(fam, key, *hit[:3], hit[3], a, z, tau)
            rec.update({'type': 'class', 'family': fam + '-deepk',
                        'cell': list(key), 'a': str(a), 'eps': eps,
                        'f': f, 'q1': q1, 'N': str(N)})
            o.write(json.dumps(rec) + '\n')
            o.flush()
            print(f"*** DEEPK-CLOSED {fam} {key} k={hit[0]}"
                  f" rung={hit[2]}", flush=True)
            return ('closed', k)
    return ('open', max_k)


def main():
    mode = sys.argv[1]
    hours = 3.0
    if '--hours' in sys.argv:
        hours = float(sys.argv[sys.argv.index('--hours') + 1])
    T0 = time.time()
    deadline = T0 + hours*3600
    if mode == 'DEEPK':
        fam, jsonls, path = sys.argv[2], sys.argv[3], sys.argv[4]
        max_k = 600
        if '--max-k' in sys.argv:
            max_k = int(sys.argv[sys.argv.index('--max-k') + 1])
        idxs = open_indices(jsonls, fam)
        if '--only' in sys.argv:
            want = {int(x) for x in
                    sys.argv[sys.argv.index('--only') + 1].split(',')}
            idxs = [i for i in idxs if i in want]
        keys_l = sorted(_L11_CLASSES) if fam == 'L11' else sorted(_L12_ESCAPE)
        o = open(path, 'w')
        o.write(json.dumps({'type': 'meta', 'mode': mode, 'family': fam,
                            'open': len(idxs), 'max_k': max_k}) + '\n')
        closed = 0
        for idx in idxs:
            if time.time() > deadline:
                break
            status, k = scan_DEEPK(fam, keys_l[idx], o, o, max_k, deadline)
            closed += status == 'closed'
            print(f"[{idx}] {keys_l[idx]} {status} k_to={k}", flush=True)
        o.write(json.dumps({'type': 'summary', 'closed': closed,
                            'of': len(idxs)}) + '\n')
        o.close()
        print(f"SUMMARY DEEPK {fam} closed {closed}/{len(idxs)}", flush=True)
        return
    lo, hi, jsonls, path = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
    idxs = open_indices(jsonls, mode)[lo:hi]
    keys_l = sorted(_L11_CLASSES) if mode == 'L11' else sorted(_L12_ESCAPE)
    o = open(path, 'w')
    o.write(json.dumps({'type': 'meta', 'mode': mode, 'lo': lo, 'hi': hi,
                        'open': len(idxs)}) + '\n')
    closed = 0
    for i, idx in enumerate(idxs):
        if time.time() > deadline:
            break
        key = keys_l[idx]
        fn = scan_L11_alt if mode == 'L11' else scan_ESC_alt
        status, tried = fn(key, o, o, deadline=deadline)
        closed += status == 'closed'
        o.write(json.dumps({'type': 'cell-done', 'idx': idx,
                            'cell': list(key), 'status': status,
                            'alt_classes_tried': tried}) + '\n')
        o.flush()
        print(f"[{idx}] {key} {status} (alt tried {tried})", flush=True)
    o.write(json.dumps({'type': 'summary', 'closed': closed,
                        'of': len(idxs)}) + '\n')
    o.close()
    print(f"SUMMARY {mode} [{lo}:{hi}] closed {closed}/{len(idxs)}",
          flush=True)


if __name__ == '__main__':
    main()
