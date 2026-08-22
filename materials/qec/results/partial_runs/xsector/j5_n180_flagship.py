"""Generate a fully-verified flagship J.5-collapse witness at 15_6_0218 (n=180,
d_Z(P)=14 exp039 exact) from the syzygy family, dual-path. Writes
j5_witness_15_6_0218.json.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import sys
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE.parent))

_spec = importlib.util.spec_from_file_location("exp027", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = E27; _spec.loader.exec_module(E27)
_spec39 = importlib.util.spec_from_file_location("exp039", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(_spec39); sys.modules[_spec39.name] = E39; _spec39.loader.exec_module(E39)

import j5core

rows = E27.load_catalogue()
parents = E39.distinct_parents(rows)
certs = E39.load_certificates()
import os as _os
fp_prefix = _os.environ.get("J5_FP", "5a61c4ae52436b35")  # 15_6_0219: d_Z exact 14
fp = next(f for f in parents if f.startswith(fp_prefix))
entry = parents[fp]
pool = E39.load_pool_lower_bounds()
cert = certs.get(fp)
d_lb = int(cert["d_z_parent"]) if cert and cert.get("d_z_exact") else int(pool[fp])
ell_p, m_p = entry["ell"], entry["m"]
row0 = next(r for r in rows
            if int(r["ell"]) == ell_p and int(r["m"]) == m_p
            and E27.matrix_fingerprint(*E27.parent_matrices(r)[1:3]) == fp)
lat = j5core.Lattice(ell_p, m_p)
p = j5core.ParentData(lat, row0["A_terms"], row0["B_terms"])
S = p.syzygy
for i in range(S.shape[0]):
    cd = S[i]
    if not cd.any():
        continue
    C, D = j5core.split_CD(p, cd)
    rw = int(p.HX[0].sum())
    demoted, sp = j5core.batch_demotion_test(
        p, C, D, [{"w": rw, "lam": (j,)} for j in range(p.dim)])
    if not demoted:
        continue
    ver = j5core.verify_full(p, cd, demoted[0]["lam"], d_lb)
    if not ver["verified"]:
        continue
    payload = {"schema": "j5-xcollapse-witness-v1",
               "parent_fingerprint": fp, "labels": [mm["label"] for mm in entry["members"][:4]],
               "cert": {"d_z_lower": d_lb,
                        "d_z_source": ("exp039-cert(exact)" if cert and cert.get("d_z_exact")
                                      else "exp037-pool(bound)"),
                        "T": cert.get("T") if cert else None,
                        "path": (f"results/partial_runs/exp039/parent_{fp[:16]}.json"
                                 if cert else f"results/partial_runs/exp037/css_{fp[:16]}.json")},
               "demoted_lam": demoted[0]["lam"], **ver}
    outname = f"j5_witness_{payload['labels'][0]}.json"
    (HERE.parent / outname).write_text(json.dumps(payload, indent=1, default=str))
    print(json.dumps({k: payload[k] for k in
                      ("labels", "demoted_lam")} | {k: ver[k] for k in
                      ("w_x", "w_z", "M_is_zero", "delta_bar", "k_Q", "rho_X", "verified")}))
    break
