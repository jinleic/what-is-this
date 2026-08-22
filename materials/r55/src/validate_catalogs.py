#!/usr/bin/env python3
"""Gate 1 of the R(5,5) campaign: independently validate McKay's Ramsey
catalogs against their published counts and against provable invariants.

Published counts source: https://users.cecs.anu.edu.au/~bdm/data/ramsey.html
(fetched 2026-08-13). Degree bounds are theorems, not data: in a
Ramsey(s,t,n)-graph, every neighborhood induces a Ramsey(s-1,t)-graph and
every non-neighborhood induces a Ramsey(s,t-1)-graph, so
    n - R(s,t-1) <= deg(v) <= R(s-1,t) - 1
using R(2,5)=5, R(3,4)=9, R(3,5)=14, R(4,4)=18, R(4,5)=25.

Writes data/VALIDATION.json and prints a markdown summary.
"""

import hashlib
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import check_file  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

R35_COUNTS = {1: 1, 2: 2, 3: 3, 4: 7, 5: 13, 6: 32, 7: 71, 8: 179, 9: 290,
              10: 313, 11: 105, 12: 12, 13: 1}
R44_COUNTS = {1: 1, 2: 2, 3: 4, 4: 9, 5: 24, 6: 84, 7: 362, 8: 2079,
              9: 14701, 10: 103706, 11: 546356, 12: 1449166, 13: 1184231,
              14: 130816, 15: 640, 16: 2, 17: 1}

# (file, s, t, expected_count, deg_lo(n), deg_hi(n))
JOBS = []
for n, c in R35_COUNTS.items():
    JOBS.append((f"r35_{n}.g6", 3, 5, c, max(0, n - 9), 4))
for n, c in R44_COUNTS.items():
    JOBS.append((f"r44_{n}.g6", 4, 4, c, max(0, n - 9), 8))
JOBS.append(("r45_24.g6", 4, 5, 352366, 24 - 18, 13))
JOBS.append(("r55_42some.g6", 5, 5, 328, 42 - 25, 24))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    workers = int(os.environ.get("WORKERS", "8"))
    report, failures = [], []
    for fname, s, t, expect, dlo, dhi in JOBS:
        path = os.path.join(DATA, fname)
        if not os.path.exists(path):
            failures.append(f"{fname}: MISSING")
            continue
        t0 = time.time()
        agg = check_file(path, s, t, workers)
        agg["seconds"] = round(time.time() - t0, 2)
        agg["sha256"] = sha256(path)
        agg["expected"] = expect
        agg["count_matches"] = agg["graphs"] == expect
        agg["deg_bounds_theory"] = [dlo, dhi]
        agg["deg_bounds_ok"] = (agg["deg_min"] >= dlo and agg["deg_max"] <= dhi)
        ok = agg["all_ramsey"] and agg["count_matches"] and agg["deg_bounds_ok"]
        agg["ok"] = ok
        if not ok:
            failures.append(f"{fname}: {json.dumps(agg)}")
        report.append(agg)
        print(f"{'PASS' if ok else 'FAIL'} {fname}: {agg['graphs']} graphs "
              f"(expect {expect}), edges [{agg['edge_min']},{agg['edge_max']}], "
              f"deg [{agg['deg_min']},{agg['deg_max']}] within {[dlo, dhi]}, "
              f"{agg['seconds']}s", flush=True)
    out = os.path.join(DATA, "VALIDATION.json")
    with open(out, "w") as f:
        json.dump({"date": "2026-08-13", "workers": workers,
                   "results": report, "failures": failures}, f, indent=1)
    print(f"\nwrote {out}")
    if failures:
        print("FAILURES:\n" + "\n".join(failures))
        sys.exit(1)
    print(f"ALL {len(report)} CATALOGS VALIDATED")


if __name__ == "__main__":
    main()
