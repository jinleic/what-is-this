"""d-range corner-case audit:
 - for each ladder stratum print the spec's Delta_range and implied q=n-1-d,
   flag q>17 (catalog file would not exist) and q<0
 - observed Delta multiset per stratum, confirm within range
 - confirm every r35_d (d<=13) / r44_q (q<=17) catalog is nonempty
"""
import os
from math import ceil
from g6lib import load_g6_file, ecount

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
STRATA = [(12, 48), (13, 53), (13, 52), (14, 60), (15, 66),
          (16, 72), (16, 71), (17, 79), (17, 78), (21, 107)]

for n, e in STRATA:
    lo, hi = ceil(2 * e / n), min(13, n - 1)
    qs = [(d, n - 1 - d) for d in range(lo, hi + 1)]
    bad = [(d, q) for d, q in qs if q > 17 or q < 0]
    deltas = {}
    for gn, adj in load_g6_file(f"{DATA}/r45extreme/r45{n}.{e}.g6"):
        D = max(bin(a).count("1") for a in adj)
        deltas[D] = deltas.get(D, 0) + 1
    inrange = all(lo <= D <= hi for D in deltas)
    print(f"n={n} e={e}: range=[{lo},{hi}] q from {n-1-hi} to {n-1-lo} "
          f"badq={bad or 'none'} observed Delta={sorted(deltas.items())} "
          f"in-range={inrange}")

print()
for d in range(1, 14):
    c = len(load_g6_file(f"{DATA}/r35_{d}.g6"))
    print(f"r35_{d}: {c} graphs{' EMPTY!' if c == 0 else ''}", end="  ")
print()
for q in range(1, 18):
    c = len(load_g6_file(f"{DATA}/r44_{q}.g6"))
    print(f"r44_{q}: {c}{' EMPTY!' if c == 0 else ''}", end="  ")
print()
