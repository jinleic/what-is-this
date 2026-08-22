#!/usr/bin/env python3
"""Regenerate the class table under the CORRECTED controlled set (2026-08-16).

The adversarial audit (L11Audit) found the frozen sets of `_L10_CLASSES` and
`_L11_CLASSES` under-inclusive: L10-0 frees a place only where v_p(x) is EVEN,
and that hypothesis was never checked off S.  Two families of permanent
odd-valuation places were skipped:

  (1) p | A, because delta carries A in its denominator on every tau != 0
      branch (1/A canonical, -4a^4/A square).  130/164 L10 rows and 293/353
      L11 rows had symbol -1 at such a p.
  (2) prime-restricted fixed divisors: P has degree 8, so a p dividing P(u)
      for every unit u without dividing the content needs p - 1 <= 8, i.e.
      p in {2, 3, 5, 7}.  p = 3 really occurs.

`_l10_class_cert` now freezes S = {2,3,5,7} u supp(alpha) u supp(delta).
Under that set the wall at p | A is branch-INDEPENDENT, and a cell is
reachable exactly when z has a numerator prime = 1 mod 4 (THEOREMS L11h):
such a prime admits an admissible rational a with v_p(1+4a^2) > 0, which is
the only way to escape the pairing.  This script builds a by construction
from that prime rather than sampling an s-pool.

Search-side tool: rows are frozen into h10q.py only after in-suite
re-verification.
"""
import json
import sys
import time
from fractions import Fraction as F

from h10q import (_l10_class_cert, _l9_grid, factorint, primerange, vp,
                  FactorBudget, PrimalityBound)


def escape_primes(z):
    """Numerator primes = 1 mod 4 of z: the only escape from the L11g wall."""
    return [p for p in sorted(factorint(abs(z.numerator))) if p % 4 == 1]


def a_candidates(p, limit=600, want=4):
    """Admissible rational a = m/n with v_p(1 + 4a^2) > 0.

    Exists for every p = 1 mod 4: solve n^2 = -4m^2 mod p and pick odd
    representatives, so a is an odd 2-adic unit and A = 1 + 4a^2 = 5 mod 8.
    """
    out = []
    for n in range(1, limit, 2):
        for m in range(1, limit // 2, 2):
            if (n * n + 4 * m * m) % p == 0:
                a = F(m, n)
                if vp(a, 2) == 0 and vp(1 + 4 * a * a, p) > 0:
                    out.append(a)
                    if len(out) >= want:
                        return out
    return out


def scan(qmax=1200, nq=70, out_path='data/l11_classes.jsonl'):
    QS = [q for q in primerange(41, qmax)][:nq]
    grid = sorted(_l9_grid())
    rows, walled, t0 = [], [], time.time()
    for i, (w, ut) in enumerate(grid):
        z = F(w) * F(*ut)
        ps = escape_primes(z)
        if not ps:
            walled.append((w, ut))
            continue
        found = None
        for p in ps:
            for a in a_candidates(p):
                A = 1 + 4 * a * a
                for tau_nm, tau in (('sq', (1 + 2 * a * a) / A),
                                    ('can', 2 * a / A)):
                    for eps in (1, -1):
                        for q1 in QS:
                            try:
                                r = _l10_class_cert(a, z, tau, eps, q1)
                            except (FactorBudget, PrimalityBound):
                                continue
                            if r and r['ok']:
                                found = (a, tau_nm, eps, q1, r['N'], r['S'])
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if found:
            a, tau_nm, eps, q1, N, S = found
            rows.append({'w': w, 'u': list(ut),
                         'a': [a.numerator, a.denominator], 'branch': tau_nm,
                         'eps': eps, 'q1': q1, 'N': N, 'S': [str(x) for x in S],
                         'p': ps[0]})
        else:
            walled.append((w, ut))
        if (i + 1) % 40 == 0:
            print(f"  [{i+1}/{len(grid)}] certified={len(rows)}"
                  f" unreached={len(walled)} {time.time()-t0:.0f}s", flush=True)
    with open(out_path, 'w') as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"certified {len(rows)}/{len(grid)}; {len(walled)} cells have NO"
          f" numerator prime = 1 mod 4 and are walled by L11g"
          f" ({time.time()-t0:.0f}s) -> {out_path}", flush=True)
    return rows, walled


if __name__ == '__main__':
    scan(out_path=sys.argv[1] if len(sys.argv) > 1
         else 'data/l11_classes.jsonl')
