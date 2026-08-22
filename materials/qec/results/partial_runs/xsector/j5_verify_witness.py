"""Independent full verification of the J.5 refutation candidate at 12_6_0193.

Two code paths for the demotion condition:
  path A (j5core): z = lam0 [C D] not in S_Z + Delta   (J.2/J.4(ii) rank test)
  path B (exp044): build S_X(Q) basis via Wperp trick, test x0 not in rowspace(S_X(Q))
Plus raw centralizer cross-check of x0 in Xcen(Q), and weight/bound assertions.
"""
from __future__ import annotations

import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
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

from qec_research.gf2.linalg import rank_np, nullspace_np, matmul as gf2matmul  # noqa: E402
import j5core  # noqa: E402

rows = E27.load_catalogue()
parents = E39.distinct_parents(rows)
certs = E39.load_certificates()

TARGET = "12_6_0193"
lat = j5core.Lattice(12, 6)
row = next(r for i, r in enumerate(rows) if E27.catalogue_label(r, i) == TARGET)
p = j5core.ParentData(lat, row["A_terms"], row["B_terms"])
fp = E27.matrix_fingerprint(p.HX, p.HZ)
cert = certs.get(fp)
print("fingerprint", fp[:16], "cert d_z_parent:", cert.get("d_z_parent") if cert else None,
      "exact:", cert.get("d_z_exact") if cert else None, "T:", cert.get("T") if cert else None)

rng = np.random.default_rng(1)
lights = j5core.light_stabilizers(p, cert["d_z_parent"])
probes, cands = j5core.sharp_probe(p, lights, cert["d_z_parent"], rng, max_lams=8, samples_per=4)
print("num cands:", len(cands))

lam0 = np.zeros(72, np.uint8); lam0[cands[0]["lam"]] = 1
cd = np.array(cands[0]["cd"], np.uint8)
C, D = j5core.split_CD(p, cd)
CD = np.hstack([C, D])
x0 = (lam0 @ p.HX) % 2
z0 = (lam0 @ CD) % 2
M = (gf2matmul(p.A, C.T) ^ gf2matmul(p.B, D.T)).astype(np.uint8)

# path A (j5core)
sp = j5core.demotion_spaces(p, C, D)
inS_A = j5core.z_in_SZplusDelta(z0, sp["probe_SZD"])
# path B (exp044 style)
Wperp = nullspace_np(sp["SZD"])
Lam = nullspace_np(((CD @ Wperp.T) % 2).T)
SXQ = (Lam @ p.HX) % 2
ns = nullspace_np(SXQ)
x0_in_SXQ = not bool(((x0 @ ns.T) % 2).sum())
print("pathA z in S_Z+Delta:", inS_A, " pathB x0 in S_X(Q):", x0_in_SXQ,
      " AGREE:", inS_A == x0_in_SXQ)

checks = {
    "M == 0": not bool(M.any()),
    "M symmetric": not bool((M ^ M.T).any()),
    "x0 weight": int(x0.sum()),
    "z0 weight": int(z0.sum()),
    "d_Z(P) certified": int(cert["d_z_parent"]),
    "x0 in ker H_Z": not bool((x0 @ p.HZ.T % 2).any()),
    "x0 in ker [C D]": not bool((x0 @ CD.T % 2).any()),
    "x0 in rowspace H_X": not bool((x0 @ nullspace_np(p.HX).T % 2).sum()),
    "z0 in ker H_X": not bool((z0 @ p.HX.T % 2).any()),
    "z0 NOT in S_Z+Delta (A)": not inS_A,
    "x0 NOT in S_X(Q) (B)": not x0_in_SXQ,
    "delta_bar": sp["delta_bar"],
    "k_Q": int(p.kP - sp["delta_bar"]),
    "A_terms": p.A_terms, "B_terms": p.B_terms,
    "C_terms": lat.terms_of(C), "D_terms": lat.terms_of(D),
    "lam_support": sorted(int(i) for i in np.flatnonzero(lam0)),
    "x0_support": sorted(int(i) for i in np.flatnonzero(x0)),
    "z0_support": sorted(int(i) for i in np.flatnonzero(z0)),
}
verified = (checks["M == 0"] and checks["x0 in ker H_Z"] and checks["x0 in ker [C D]"]
            and checks["x0 in rowspace H_X"] and checks["z0 NOT in S_Z+Delta (A)"]
            and checks["x0 NOT in S_X(Q) (B)"] and checks["x0 weight"] < checks["d_Z(P) certified"])
checks["VERIFIED"] = bool(verified)
print(json.dumps({k: v for k, v in checks.items() if not isinstance(v, list) or len(str(v)) < 200},
                 indent=1, default=str))
if verified:
    out = HERE.parent / "j5_witness_12_6_0193.json"
    out.write_text(json.dumps({"schema": "j5-xcollapse-witness-v1", **checks,
                               "parent_fingerprint": fp, "cert": {"d_z_parent": cert["d_z_parent"],
                                                                  "d_z_exact": cert["d_z_exact"],
                                                                  "T": cert.get("T"),
                                                                  "path": f"results/partial_runs/exp039/parent_{fp[:16]}.json"}},
                              indent=1, default=str))
    print("WITNESS SAVED", out)
