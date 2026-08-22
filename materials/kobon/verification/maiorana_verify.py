#!/usr/bin/env python3
"""Fail-closed independent certification of the Maiorana k=14 / 54-triangle
arrangement (sol1). One command, structured JSONL output, nonzero exit on
ANY check failure.

Checks (all on the SAME pinned file object):
  1. SHA-256 of the pinned source file equals the recorded manifest hash.
  2. The 14 lines are pairwise non-proportional (no duplicate lines).
  3. Computed concurrency census equals the file's declared_triple_points.
  4. Exact-rational empty-interior census (record convention) returns 54.
  5. Every selected triangle stays empty under the stricter cevian-aware
     open-interior test (lines through a vertex can still cut the interior).
  6. The constructed 54-selection passes the campaign-audited
     math/kobon/engine.py::verify_selection(n=14, minimum=54).

Usage:
  scratch/kobon-audit/venv/bin/python \\
      scratch/kobon/n14/maiorana_verify.py \\
      scratch/kobon/n14/sol1_lines_rational.json
"""
import json
import sys
import hashlib
import itertools
from collections import defaultdict
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXPECT_SHA = ('83fc26666d73cd39791018a72aba71caeac4a3f7efb27109'
              'ed2955b4540b2523')


def fail(msg, **kw):
    print(json.dumps({'stage': 'FAIL', 'reason': msg, **kw}), flush=True)
    sys.exit(1)


def crossings(L):
    k = len(L)
    Pn = {}
    for i, j in itertools.combinations(range(k), 2):
        a1, b1, c1 = L[i]
        a2, b2, c2 = L[j]
        det = a1 * b2 - a2 * b1
        if det == 0:
            continue
        Pn[(i, j)] = ((c1 * b2 - c2 * b1) / det,
                      (a1 * c2 - a2 * c1) / det)
    return Pn


def census_census(L, Pn):
    pt_lines = defaultdict(set)
    for (i, j), pt in Pn.items():
        pt_lines[pt].add(i)
        pt_lines[pt].add(j)
    return sorted(tuple(sorted(v)) for v in pt_lines.values() if len(v) >= 3)


def select_empty(L, Pn):
    """Record convention: triple selected iff no other line strictly
    separates two of its three vertices."""
    k = len(L)
    sel = []
    for T in itertools.combinations(range(k), 3):
        i, j, l = T
        if not all(x in Pn for x in [(i, j), (i, l), (j, l)]):
            continue
        v = [Pn[(i, j)], Pn[(i, l)], Pn[(j, l)]]
        (x1, y1), (x2, y2), (x3, y3) = v
        if (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) == 0:
            continue
        empty = True
        for m in range(k):
            if m in (i, j, l):
                continue
            a, b, c = L[m]
            s = [a * x + b * y - c for x, y in v]
            if any((s[p] < 0 < s[q]) for p in range(3) for q in range(3)):
                empty = False
                break
        if empty:
            sel.append(T)
    return sel


def cevian_check(L, Pn, sel):
    """Stricter test on the selection: no other line enters the open
    interior, including the subtle case where it only passes through a
    vertex but still separates an open side."""
    k = len(L)
    bad = 0
    for (i, j, l) in sel:
        v = [Pn[(i, j)], Pn[(i, l)], Pn[(j, l)]]
        for m in range(k):
            if m in (i, j, l):
                continue
            a, b, c = L[m]
            s = [a * x + b * y - c for x, y in v]
            pos = any(x > 0 for x in s)
            neg = any(x < 0 for x in s)
            zer = sum(1 for x in s if x == 0)
            if pos and neg:
                bad += 1
                break
            if zer:
                for (Pa, Qa) in [(v[0], v[1]), (v[0], v[2]), (v[1], v[2])]:
                    tn = -(a * Pa[0] + b * Pa[1] - c)
                    td = a * (Qa[0] - Pa[0]) + b * (Qa[1] - Pa[1])
                    if td == 0:
                        continue
                    if 0 < tn / td < 1:
                        bad += 1
                        break
                if bad:
                    break
    return bad


def main():
    if len(sys.argv) != 2:
        fail('usage: maiorana_verify.py solN/lines_rational.json')
    path = Path(sys.argv[1])
    blob = path.read_bytes()
    sha = hashlib.sha256(blob).hexdigest()
    if sha != EXPECT_SHA:
        fail('SHA mismatch', got=sha, want=EXPECT_SHA)
    print(json.dumps({'stage': 'sha_pinned_source', 'ok': True,
                      'sha256': sha}), flush=True)

    data = json.loads(blob)
    L = [tuple(Fr(s) for s in t) for t in data['lines_frac']]
    k = len(L)
    if k != 14:
        fail(f'line count {k} != 14')
    if len(set(L)) != k:
        fail('duplicate (proportional) lines')
    print(json.dumps({'stage': 'lines_distinct', 'ok': True, 'k': k}),
          flush=True)

    Pn = crossings(L)
    trips = census_census(L, Pn)
    declared = sorted(tuple(sorted(t)) for t in
                      data.get('declared_triple_points', []))
    if trips != declared:
        fail('multipoint census != declared_triple_points',
             computed=trips, declared=declared)
    print(json.dumps({'stage': 'multipoint_census', 'ok': True,
                      'triples': [list(t) for t in trips]}), flush=True)

    sel = select_empty(L, Pn)
    if len(sel) != 54 or len(sel) != data['count']:
        fail('record census mismatch', found=len(sel),
             claimed=data['count'])
    print(json.dumps({'stage': 'record_census', 'ok': True,
                      'count': len(sel)}), flush=True)

    bad = cevian_check(L, Pn, sel)
    if bad:
        fail('cevian-aware open-interior check', violating_triangles=bad)
    print(json.dumps({'stage': 'cevian_check', 'ok': True}), flush=True)

    sys.path.insert(0, str(ROOT / 'math' / 'kobon'))
    import engine
    ms = [-a / b for (a, b, c) in L]
    bs = [c / b for (a, b, c) in L]
    ok, why = engine.verify_selection(k, ms, bs, sel, minimum=54)
    if not ok:
        fail('engine.verify_selection', detail=why)
    print(json.dumps({'stage': 'engine_verify_selection', 'ok': True,
                      'detail': why, 'minimum': 54}), flush=True)

    print(json.dumps({'stage': 'CERTIFIED', 'K_gen_14_lower_bound': 54,
                      'encoding_convention': 'campaign engine K_gen'}),
          flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
