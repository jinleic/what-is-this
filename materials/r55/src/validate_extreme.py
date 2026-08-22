#!/usr/bin/env python3
"""Validate the r45extreme catalog: complete sets of Ramsey(4,5,n)-graphs
with extreme edge counts, orders 4..23 (~8M graphs, 247MB).

Each file r45extreme/r45{n}.{e}.g6 claims: all Ramsey(4,5,n)-graphs with
exactly e edges. Checked per file: every graph has n vertices, e edges,
no K4, no I5, and degrees within the provable window
    max(0, n - R(4,4)) <= deg <= R(3,5) - 1 = 13,  R(4,4) = 18.
Counts per file are not published on the source page; the filename encoding
is the cross-check. Writes data/VALIDATION_extreme.json.
"""

import glob
import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import check_file  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def main():
    workers = int(os.environ.get("WORKERS", "8"))
    files = sorted(glob.glob(os.path.join(DATA, "r45extreme", "r45*.g6")))
    if not files:
        sys.exit("no r45extreme files found; extract r45extreme.tar.gz first")
    report, failures, total = [], [], 0
    for path in files:
        # r45{n}.{e}.g6: all Ramsey(4,5,n) with exactly e edges (n=10..23);
        # r45{n}.g6: the complete Ramsey(4,5,n) set (n=4..9), no edge claim.
        m = re.match(r"r45(\d+)(?:\.(\d+))?\.g6$", os.path.basename(path))
        n = int(m.group(1))
        e = int(m.group(2)) if m.group(2) else None
        t0 = time.time()
        agg = check_file(path, 4, 5, workers)
        agg["seconds"] = round(time.time() - t0, 2)
        ok = (agg["all_ramsey"]
              and set(agg["n_values"]) == {n}
              and (e is None or agg["edge_min"] == agg["edge_max"] == e)
              and agg["deg_min"] >= max(0, n - 18) and agg["deg_max"] <= 13)
        agg["ok"] = ok
        total += agg["graphs"]
        if not ok:
            failures.append(os.path.basename(path))
        report.append(agg)
        print(f"{'PASS' if ok else 'FAIL'} r45_{n} e={e if e is not None else 'all'}: "
              f"{agg['graphs']} graphs, {agg['seconds']}s", flush=True)
    out = os.path.join(DATA, "VALIDATION_extreme.json")
    with open(out, "w") as f:
        json.dump({"date": "2026-08-13", "total_graphs": total,
                   "results": report, "failures": failures}, f, indent=1)
    print(f"\n{total} graphs in {len(files)} files; wrote {out}")
    if failures:
        print("FAILURES:", failures)
        sys.exit(1)
    print("ALL EXTREME FILES VALIDATED")


if __name__ == "__main__":
    main()
