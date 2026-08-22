"""Timing probe: how expensive is exact CSS distance certification at n = 162?

Scratch tool for planning EXP-055's certification budget.  Finds the weight-3
pairs on (9,9) with 8 <= k <= 16, ranks them by the EXP-055 certified ceiling,
and times exact_distance_css on the strongest few.
"""
import importlib.util
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

spec = importlib.util.spec_from_file_location(
    "e55", ROOT / "experiments" / "exp055_odd_lattice_sweep.py")
E55 = importlib.util.module_from_spec(spec)
sys.modules["e55"] = E55
spec.loader.exec_module(E55)

from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402

ELL, M = 9, 9
bar = E55.bar_permutation(ELL, M)
sup = [s for s in E55.normalised_supports(ELL, M, 3) if len(s) == 3]
anns = {}
for s in sup:
    A = E55.ann_basis(s, ELL, M)
    if A.shape[0] >= 4:
        anns[tuple(s)] = A
print(f"weight-3 supports: {len(sup)}  with dim Ann >= 4: {len(anns)}", flush=True)

keys = list(anns)
best = []
for i, ka in enumerate(keys):
    for kb in keys[i:]:
        I = E55.intersect(anns[ka], anns[kb])
        k = 2 * I.shape[0]
        if not 8 <= k <= 16:
            continue
        cel = E55.certified_ceiling(I, ELL, M, bar)
        if cel.get("ceiling"):
            best.append((cel["ceiling"], k, list(ka), list(kb), cel["dim_I0"]))
best.sort(key=lambda t: (-t[0], t[1]))
print(f"candidates: {len(best)}  top: {[(b[0], b[1]) for b in best[:10]]}", flush=True)

for ceil, k, A, B, d0 in best[:3]:
    HX, HZ = E55.E53.bb_from_terms(ELL, M, A, B)
    kd = 2 * ELL * M - rank_np(HX) - rank_np(HZ)
    t0 = time.time()
    res = exact_distance_css(HX, HZ, time_limit_s=240, workers=8)
    print(f"n={2*ELL*M} k={k}({kd}) ceiling={ceil} dimI0={d0} -> "
          f"d={res['d']} d_X={res['d_X']} d_Z={res['d_Z']} exact={res['d_exact']} "
          f"[{time.time()-t0:.0f}s] A={A} B={B}", flush=True)
