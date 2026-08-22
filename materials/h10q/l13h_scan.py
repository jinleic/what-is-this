#!/usr/bin/env python3
"""Per-class emergent-free (zero-bad) member certification for hypothesis H.

For each aligned class of the 353-cell grid (190 L11 prime-b classes + 103
L12 escape classes), scan prime members Q = q1 + k*N and decide each with
the L13f ladder (l13_filter.smooth_emergent + cofactor_decide).  First
PROVED 'zero' member closes the class (zero-bad => all Hilbert symbols +1
=> tied conic globally soluble, Hasse-Minkowski).  'zero-jacobi' rows are
evidence tier and do NOT close a class.

Usage:
  python3 l13h_scan.py L11 LO HI OUT.jsonl [--max-k 200] [--hours H]
  python3 l13h_scan.py ESC LO HI OUT.jsonl [...]

Line schema (one JSON per line):
  {"type":"meta", ...}
  {"type":"member", "class_idx":i, "cell":[w,[u1,u2]], "k":k, "Q":str,
   "verdict":v, "rung":str|None, "digits":int|None}
  {"type":"class", "family":fam, "class_idx":i, "cell":[w,[u1,u2]],
   "status":"closed"|"open", "k_zero":k|None, "Q":str|None,
   "rung":str|None, "verdicts":{...}, "excluded_skips":n,
   "tied_crosscheck":"True"|"refused:..."|None, "elapsed_s":float}
  {"type":"summary", ...}

Authorities: h10q.py::_L11_CLASSES rows ((a,end),eps,q1,N,w);
h10q.py::_L12_ESCAPE (f,eps,q1) + l12_class.l12_class_cert (N, excluded).
No edits to h10q.py or l13_filter.py.
"""
import json
import sys
import time
from fractions import Fraction as F

from h10q import (_L11_CLASSES, _L12_ESCAPE, _is_prime, _l7_tied_status,
                  PrimalityBound)
import l12_class
from l13_filter import smooth_emergent, cofactor_decide


def excluded_primes_L11(a, z, tau):
    from h10q import factorint, _l10_P
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    Z = z**3
    D = 1 - Z - a*a*Z*Z
    parts = [D.numerator, D.denominator, z.numerator, z.denominator,
             a.numerator, a.denominator, delta.numerator, delta.denominator]
    out = set()
    for p in parts:
        for q in factorint(abs(p)):
            out.add(q)
    return out


def class_data(fam, key):
    w, ut = key
    z = F(w)*F(*ut)
    if fam == 'L11':
        row = _L11_CLASSES[key]
        a = F(*row[0])
        eps, q1, N = row[1], row[2], row[3]
        f = 1
        A = 1 + 4*a*a
        tau = (1 + 2*a*a)/A
        excl = excluded_primes_L11(a, z, tau)
    else:
        f_t, eps, q1 = _L12_ESCAPE[key]
        cert = l12_class.l12_class_cert(w, *ut)
        N = cert['N']
        f = cert['f']
        assert f == f_t, (key, f, f_t)
        a = F(1)
        tau = F(3, 5)
        excl = set(cert['excluded'])
    return a, eps, f, q1, N, tau, z, excl


def scan_class(fam, key, idx, out, max_k, class_deadline):
    t0 = time.time()
    a, eps, f, q1, N, tau, z, excl = class_data(fam, key)
    verdicts = {}
    excl_skips = nonprime = 0
    closed = None
    k = -1
    while time.time() < class_deadline and k < max_k:
        k += 1
        Q = q1 + k*N
        if Q in excl:
            excl_skips += 1
            continue
        try:
            if not _is_prime(Q):
                nonprime += 1
                continue
        except PrimalityBound:
            nonprime += 1
            continue
        b = F(eps*f*Q)
        st, dt = smooth_emergent(a, z, tau, b)
        if st == 'zero':
            v, rung, digits = 'zero', 'square', None
        elif st == 'cofactor-big':
            v2, info = cofactor_decide(a, z, tau, b, detail=dt,
                                       deadline=class_deadline)
            v = v2
            rung = info.get('rung')
            digits = max((len(info['rn']), len(info['rd'])))
        else:
            v, rung, digits = st, None, None
        verdicts[v] = verdicts.get(v, 0) + 1
        if v in ('zero', 'bad-1', 'INCONSISTENT', 'align-fail'):
            out.write(json.dumps({
                'type': 'member', 'family': fam, 'class_idx': idx,
                'cell': list(key), 'k': k, 'Q': str(Q), 'verdict': v,
                'rung': rung, 'digits': digits}) + '\n')
            out.flush()
        if v == 'align-fail':
            print(f"*** CLASS-CONTRADICTION {fam} {key} k={k} Q={Q}",
                  flush=True)
        if v == 'zero':
            tied = None
            try:
                tied = str(_l7_tied_status(a, b, z, tau))
            except Exception as e:
                tied = f"refused:{type(e).__name__}"
            closed = (k, Q, rung, tied)
            break
    rec = {'type': 'class', 'family': fam, 'class_idx': idx,
           'cell': list(key),
           'status': 'closed' if closed else 'open',
           'k_zero': closed[0] if closed else None,
           'Q': str(closed[1]) if closed else None,
           'rung': closed[2] if closed else None,
           'tied_crosscheck': closed[3] if closed else None,
           'verdicts': verdicts, 'excluded_skips': excl_skips,
           'nonprime_skips': nonprime,
           'a': str(a), 'eps': eps, 'f': f, 'q1': q1, 'N': str(N),
           'elapsed_s': round(time.time() - t0, 2)}
    out.write(json.dumps(rec) + '\n')
    out.flush()
    return rec


def main():
    fam, lo, hi, path = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    max_k, hours = 200, 2.8
    args = sys.argv[5:]
    if '--max-k' in args:
        max_k = int(args[args.index('--max-k') + 1])
    if '--hours' in args:
        hours = float(args[args.index('--hours') + 1])
    assert fam in ('L11', 'ESC')
    keys = sorted(_L11_CLASSES) if fam == 'L11' else sorted(_L12_ESCAPE)
    keys = keys[lo:hi]
    out = open(path, 'w')
    out.write(json.dumps({'type': 'meta', 'family': fam, 'lo': lo,
                          'hi': hi, 'n_classes': len(keys), 'max_k': max_k,
                          'hours': hours, 'ts': time.time()}) + '\n')
    out.flush()
    T0 = time.time()
    deadline = T0 + hours*3600
    n_closed = 0
    for i, key in enumerate(keys):
        per_class = min(deadline, time.time() + deadline - T0)
        class_deadline = min(deadline, time.time() + 600)
        if time.time() > deadline:
            break
        rec = scan_class(fam, key, lo + i, out, max_k, class_deadline)
        n_closed += rec['status'] == 'closed'
        print(f"[{lo+i}] {rec['cell']} {rec['status']}"
              f" k_zero={rec['k_zero']} verdicts={rec['verdicts']}"
              f" ({rec['elapsed_s']}s)", flush=True)
    out.write(json.dumps({
        'type': 'summary', 'family': fam, 'classes': len(keys),
        'closed': n_closed, 'wall_s': round(time.time() - T0, 1)}) + '\n')
    out.close()
    print(f"SUMMARY {fam}[{lo}:{hi}] closed {n_closed}/{len(keys)}",
          flush=True)


if __name__ == '__main__':
    main()
