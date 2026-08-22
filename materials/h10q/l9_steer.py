#!/usr/bin/env python3
"""L9 reciprocity-steered assembly probe (2026-08-16).

Search-side tool (peer of l6_search.py): NOT part of the verified suite.
Results it finds are frozen into h10q.py (_L9_STEERED) only after in-suite
re-verification.

Target: the OPEN L6 assembly lemma.  For a cell (w, z), z in m_w, find an
admissible witness (s, b) in Phi_1^{S={2}} (v2(s) >= 0 or s = 0, v2(b) = 0)
and a branch tau in {0, 2a/A} whose TIED conic is soluble.  All certificate
logic lives in h10q._l9_steered_solvable (single source; see THEOREMS.md
L9/L9a and the docstrings there): d resolves completely, x resolves into
proven prime powers times even-exponent composite blobs coprime to
supp(d), every place of T carries exact symbol +1, and with exactly one
wild odd-multiplicity prime its symbol is forced by reciprocity (L9) —
so the search only needs smooth-place alignment, never wild-prime luck.

Each hit is additionally replayed against the independent in-suite path
_l7_tied_status(a, b, z, tau): True is a second proof; False is a hard
contradiction (abort); 'budget' means that path refused and the hit rests
on the steered certificate alone (counted separately).  Refusals/skips
never count as evidence.

Output: one JSON line per cell to data/l9_steer_run.jsonl (run artifact,
search-side provenance only, written incrementally) + a printed summary.
"""
import json
import sys
import time
from fractions import Fraction

from h10q import (primerange, vp, _l7_tied_status, _l9_steered_solvable,
                  _L9_U_POOL)

S_POOL = [Fraction(v) for v in (0, 1, -1, 2, -2, 3, -3)] + \
         [Fraction(1, 3), Fraction(-1, 3), Fraction(4), Fraction(-4),
          Fraction(1, 5), Fraction(2, 5), Fraction(5), Fraction(-5)]
S_POOL_DEEP = S_POOL + \
    [Fraction(v) for v in (6, -6, 7, -7, 8, -8)] + \
    [Fraction(3, 5), Fraction(-3, 5), Fraction(2, 7), Fraction(-2, 7),
     Fraction(1, 7), Fraction(4, 3), Fraction(-4, 3), Fraction(2, 3),
     Fraction(-2, 3), Fraction(1, 9), Fraction(5, 3), Fraction(-5, 3)]


def b_pool(w, deep=False):
    """Odd/odd rationals ordered by height, plus w-dependent entries."""
    pool = []
    nmax, mmax = (99, 20) if deep else (35, 12)
    for n in range(-nmax, nmax + 1, 2):
        for m in range(1, mmax, 2):
            if n in (0, m) or Fraction(n, m) in (0, 1):
                continue
            q = Fraction(n, m)
            if q not in pool:
                pool.append(q)
    extras = [Fraction(w), Fraction(-w), Fraction(3 * w), Fraction(-3 * w),
              Fraction(1, w), Fraction(-1, w)]
    if deep:
        extras += [Fraction(w + 2), Fraction(-w - 2), Fraction(w - 2),
                   Fraction(2 - w), Fraction(5 * w), Fraction(-5 * w),
                   Fraction(3, w), Fraction(-3, w), Fraction(w, 3),
                   Fraction(-w, 3)]
    for extra in extras:
        if extra not in pool and extra not in (0, 1):
            pool.append(extra)
    pool.sort(key=lambda q: (abs(q.numerator) * q.denominator, q < 0))
    return pool


def steer_cell(w, u, max_candidates=20000, time_budget=None, deep=False):
    """Search one cell (w, z = w*u).  Returns a result dict."""
    if time_budget is None:
        time_budget = 25.0 + 0.6 * w
    z = Fraction(w) * u
    assert vp(z, w) >= 1, (w, u)
    t0 = time.time()
    tried = skipped = certified_fail = 0
    pool_b = b_pool(w, deep=deep)
    for s in (S_POOL_DEEP if deep else S_POOL):
        if not (s == 0 or vp(s, 2) >= 0):
            continue
        a = 1 + 2 * s
        A = 1 + 4 * a * a
        for tau in (Fraction(0), 2 * a / A):
            for b in pool_b:
                if vp(b, 2) != 0:
                    continue
                if tried >= max_candidates or time.time() - t0 > time_budget:
                    return {'w': w, 'u': [u.numerator, u.denominator],
                            'ok': False, 'tried': tried, 'skipped': skipped,
                            'certified_fail': certified_fail,
                            'reason': 'exhausted'}
                tried += 1
                ok, mech, wilds = _l9_steered_solvable(a, b, z, tau)
                if ok is None:
                    skipped += 1
                    continue
                if ok is False:
                    certified_fail += 1
                    continue
                # ground-truth replay (independent Hasse-Minkowski path);
                # False would contradict the steered certificate: abort.
                gt = _l7_tied_status(a, b, z, tau)
                assert gt is not False, (w, u, s, b, tau)
                return {'w': w, 'u': [u.numerator, u.denominator],
                        'ok': True,
                        's': [s.numerator, s.denominator],
                        'b': [b.numerator, b.denominator],
                        'tau': [tau.numerator, tau.denominator],
                        'mech': mech,
                        'wilds': [str(q) for q in wilds],
                        'gt_replayed': (gt is True),
                        'tried': tried, 'skipped': skipped,
                        'certified_fail': certified_fail,
                        'secs': round(time.time() - t0, 2)}
    return {'w': w, 'u': [u.numerator, u.denominator], 'ok': False,
            'tried': tried, 'skipped': skipped,
            'certified_fail': certified_fail, 'reason': 'pool_exhausted'}


# Authority moved to h10q.py (_L9_U_POOL) so the verifier and this driver
# cannot drift.  The original literal is PRESERVED below as a drift guard.
_LEGACY_U_POOL = [Fraction(v) for v in (1, -1, 2, -2, 3, -3, 5, -5)] + \
                 [Fraction(n, m) for (n, m) in
                  [(1, 3), (-1, 3), (2, 3), (-1, 5), (7, 3), (1, 7), (-2, 7)]]
U_POOL = _L9_U_POOL
assert U_POOL == _LEGACY_U_POOL, "canonical u-pool drifted from the frozen literal"


def rescue(chunk, nchunks):
    """Retry the failures of data/l9_steer_run.jsonl with deep pools and
    large budgets; process failures[chunk::nchunks]."""
    fails = []
    with open('data/l9_steer_run.jsonl') as fh:
        for line in fh:
            r = json.loads(line)
            if not r['ok']:
                fails.append((r['w'], Fraction(*r['u'])))
    todo = fails[chunk::nchunks]
    print(f"L9 deep rescue chunk {chunk}/{nchunks}: {len(todo)} of"
          f" {len(fails)} failed cells", flush=True)
    ncov = 0
    with open(f'data/l9_rescue_{chunk}.jsonl', 'w') as fh:
        for w, u in todo:
            r = steer_cell(w, u, max_candidates=120000, time_budget=420.0,
                           deep=True)
            fh.write(json.dumps(r) + '\n')
            fh.flush()
            ncov += 1 if r['ok'] else 0
            print(f"  w={w} u={u}: "
                  + (f"COVERED mech={r['mech']} secs={r['secs']}" if r['ok']
                     else f"still open ({r['reason']}, tried={r['tried']},"
                          f" certfail={r['certified_fail']})"), flush=True)
    print(f"rescue chunk {chunk}: {ncov}/{len(todo)} newly covered")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'rescue':
        rescue(int(sys.argv[2]), int(sys.argv[3]))
        return
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    cells = []
    for w in primerange(3, wmax):
        for u in U_POOL:
            z = Fraction(w) * u
            if vp(z, w) >= 1:          # genuine target cell: z in m_w
                cells.append((w, u))
    print(f"L9 steered probe: {len(cells)} cells (odd w < {wmax}; fixed"
          f" {len(U_POOL)}-value u-pool filtered per target by v_w(u) >= 0,"
          " i.e. v_w(w*u) >= 1)", flush=True)
    rows, t0 = [], time.time()
    mechs, failures, gt_budget = {}, [], 0
    with open('data/l9_steer_run.jsonl', 'w') as fh:
        for i, (w, u) in enumerate(cells):
            r = steer_cell(w, u)
            rows.append(r)
            fh.write(json.dumps(r) + '\n')
            fh.flush()
            if r['ok']:
                mechs[r['mech']] = mechs.get(r['mech'], 0) + 1
                gt_budget += 0 if r['gt_replayed'] else 1
            else:
                failures.append(r)
            if (i + 1) % 25 == 0:
                print(f"  [{i+1}/{len(cells)}]"
                      f" covered={sum(1 for x in rows if x['ok'])}"
                      f" mechs={mechs} gt_budget={gt_budget}"
                      f" fails={len(failures)} ({time.time()-t0:.0f}s)",
                      flush=True)
    ncov = sum(1 for r in rows if r['ok'])
    print(f"\ncoverage: {ncov}/{len(cells)} cells")
    print(f"mechanisms: {mechs}; ground-truth replay refused (budget) on"
          f" {gt_budget} hits (steered single-path certificates)")
    if failures:
        print(f"UNCOVERED ({len(failures)}):")
        for r in failures:
            print(f"  w={r['w']} u={Fraction(*r['u'])} tried={r['tried']}"
                  f" skipped={r['skipped']} certfail={r['certified_fail']}"
                  f" reason={r['reason']}")
    else:
        print("ALL CELLS COVERED — every certificate is a complete-"
              "resolution W2 symbol check; ground truth replayed where"
              " the independent path did not refuse")
    print(f"total {time.time()-t0:.0f}s; artifact data/l9_steer_run.jsonl"
          " (search-side provenance, not suite evidence)")


if __name__ == '__main__':
    main()
