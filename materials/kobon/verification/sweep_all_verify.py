#!/usr/bin/env python3
"""Replayable all-15 sweep: verify every Maiorana solN arrangement from the
persisted local files in scratch/kobon/n14/maiorana/ against the persisted
SHA256SUMS.txt manifest. Exact-rational empty-interior census; each run
prints one JSONL line per solution plus a SUMMARY line; exits 1 on any
failure. Convention: record/verifier separation test (documented in
maiorana_verify.py); sol1 additionally gets the cevian-aware open-interior
cross-check and the audited engine verify_selection there.

Usage:
  scratch/kobon-audit/venv/bin/python scratch/kobon/n14/sweep_all_verify.py
"""
import json
import sys
import hashlib
import itertools
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / 'maiorana'


def main():
    man = {}
    for row in (DATA / 'SHA256SUMS.txt').read_text().strip().split('\n'):
        h, f = row.split()
        man[f] = h
    ok_all = True
    for N in range(1, 16):
        fp = DATA / f'sol{N}_lines_rational.json'
        want = man.get(f'sol{N}/lines_rational.json')
        blob = fp.read_bytes()
        sha = hashlib.sha256(blob).hexdigest()
        sha_ok = (sha == want)
        d = json.loads(blob)
        L = [tuple(Fr(s) for s in t) for t in d['lines_frac']]
        k = len(L)
        distinct = len(set(L)) == k
        Pn = {}
        for i, j in itertools.combinations(range(k), 2):
            a1, b1, c1 = L[i]
            a2, b2, c2 = L[j]
            det = a1 * b2 - a2 * b1
            if det == 0:
                continue
            Pn[(i, j)] = ((c1 * b2 - c2 * b1) / det,
                          (a1 * c2 - a2 * c1) / det)
        cnt = 0
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
            cnt += empty
        ok = sha_ok and distinct and k == 14 and cnt == d['count'] == 54
        ok_all = ok_all and ok
        print(json.dumps({'sol': N, 'sha_ok': sha_ok, 'distinct_lines': distinct,
                          'k': k, 'count': cnt, 'claimed': d['count'],
                          'pass': ok}), flush=True)
    print(json.dumps({'stage': 'SUMMARY', 'solutions_verified':
                      sum(1 for N in range(1, 16)), 'all_pass': ok_all}),
          flush=True)
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
