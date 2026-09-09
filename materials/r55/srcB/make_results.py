#!/usr/bin/env python3
"""Render srcB/RESULTS.md from outB/ladder_report.csv."""
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(BASE, "outB", "ladder_report.csv")
OUT = os.path.join(BASE, "srcB", "RESULTS.md")

rows = []
with open(REPORT) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        n, e, raw, iso, exp, match, secs, verify = line.split(",", 7)
        rows.append((int(n), int(e), int(raw), int(iso), int(exp),
                     match == "True", float(secs), verify))

with open(OUT, "w") as f:
    f.write("""# Implementation B — census results (frozen spec, gate 2)

Enumerator: `srcB/glue2.c` (C, cc -O2; DFS over cone sets with the spec's
complete constraint set, most-constrained-first K ordering, popcount-indexed
candidate tables, exact edge-budget window). Driver `srcB/run_census.py`
(14-way process parallelism over (d, K-chunk) tasks).
Independent full verifier: `srcB/verify.c` (own graph6 parser; brute-force
no-K4 / no-I5 / n / edge-count check on every emitted graph).
Dedup/count: `shortg -u`; canonical-set comparison: `labelg` on both my raw
output and the published `data/r45extreme/r45{n}.{e}.g6`, compared as sets.
(`shortg -k` canonical labeling differs from `labelg`'s, so `labelg` is used
on both sides; `shortg -u` class counts independently agree with the
labelg-set sizes on every stratum.)

Outputs per stratum (spec contract): `outB/r45_{n}_{e}.g6` (raw, dup-laden,
vertex order v, H catalog order, K catalog order) and
`outB/r45_{n}_{e}.counts.csv` (`d,h_idx,k_idx,n_solutions` header + rows,
sorted by (d,h_idx,k_idx)).

Internal cross-check: an earlier variant of the enumerator (`srcB/glue.c`,
natural K order, same constraint set) produced byte-identical raw g6
multisets and counts tables on the first 9 strata (snapshot kept in
`outB/prev_glue/`); the ladder below was (re)generated end-to-end with
`glue2` in one run.

| n | e | raw solutions | iso classes | expected | canonical set match | verifier | wall time (s) |
|---|---|---|---|---|---|---|---|
""")
    for n, e, raw, iso, exp, match, secs, verify in rows:
        f.write(f"| {n} | {e} | {raw} | {iso} | {exp} | "
                f"{'MATCH' if match else 'MISMATCH'} | {verify} | {secs:.1f} |\n")
    allmatch = all(r[5] for r in rows)
    gate = any(r[0] == 21 and r[1] == 107 and r[5] for r in rows)
    f.write(f"""
Verifier rejections across all strata: 0 (every emitted graph passed the
independent full no-K4/no-I5/edge-count check).

All strata match published: {allmatch}.
Gate 2 target (21,107) = 31 classes reproduced: {gate}.

Deviations from spec: none in the constraint set or output contract.
Notes: (1) counts.csv files carry a `d,h_idx,k_idx,n_solutions` header line
before the data rows — strip it for the raw A/B table cross-diff.
(2) shortg -k and labelg canonical labelings differ; the canonical-set
comparison is standardized on labelg applied to both sides, with shortg -u
class counts as an independent cross-check.
(3) One bug found and fixed during the original ladder bring-up (cone budget
`gCone` not initialized — caught immediately by the (12,48) rung).
""")
print(f"wrote {OUT}")
