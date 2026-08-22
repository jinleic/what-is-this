"""Probe: which catalogue parents sit at n=180/360, and do they carry exp039
d_z certificates / exp037 pool bounds?"""
from __future__ import annotations

import importlib.util
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

_spec39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(_spec39)
sys.modules[_spec39.name] = E39
_spec39.loader.exec_module(E39)

rows = E27.load_catalogue()
parents = E39.distinct_parents(rows)
certs = E39.load_certificates()
pool = E39.load_pool_lower_bounds()

by_n = {}
for fp, e in parents.items():
    by_n.setdefault(e["n"], []).append(fp)
for n in sorted(by_n):
    print("n", n, "count", len(by_n[n]))

print()
for n in (180, 360):
    print(f"--- n={n} ---")
    for fp in by_n.get(n, []):
        e = parents[fp]
        c = certs.get(fp)
        dz = c["d_z_parent"] if c and c.get("d_z_exact") else None
        lb = pool.get(fp)
        lab = [m["label"] for m in e["members"][:2]]
        print(fp[:16], e["ell"], e["m"], "k?", lab,
              "d_z_exact=", dz, "pool=", lb)
