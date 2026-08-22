"""Smoke test of j5core on a catalogue parent: 12_6_0193 (ell=12,m=6) and phase2_58."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
from pathlib import Path
import importlib.util

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE.parent))

_spec = importlib.util.spec_from_file_location("exp027", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

import j5core

rows = E27.load_catalogue()
lat = j5core.Lattice(12, 6)
row = next(r for i, r in enumerate(rows) if E27.catalogue_label(r, i) == "12_6_0193")
p = j5core.ParentData(lat, row["A_terms"], row["B_terms"])
print("k_P", p.kP, "rHX", p.rHX, "rHZ", p.rHZ, "dimV", p.validity.shape[0],
      "dimSyz", p.syzygy.shape[0], "sigma", p.sigma)
# EXP-044 expected: Xcen(P)=78 -> dim ker HZ = 78? (kP=12, rHX=?? ) check J-consistent counts
print("Xcen(P) dim =", 2 * 72 - p.rHZ, " expect 78")

# catalogue's own perturbation: demotion analysis
C = lat.mat(row["C_terms"]); D = lat.mat(row["D_terms"])
sp = j5core.demotion_spaces(p, C, D)
print("M symmetric:", not (sp["M"] ^ sp["M"].T).any(), " delta_bar:", sp["delta_bar"],
      "(expect 0 for this row)", " dim kerM:", sp["kerM"].shape[0])
print("rho_X =", int(__import__("qec_research.gf2.linalg", fromlist=["rank_np"]).rank_np(
    np.vstack([p.HZ, sp["CD"]]))) - p.rHZ, "(expect 46)")

# trivial family (B^T, A^T): M=0 and no demotion at all expected
Ct, Dt = p.B.T.copy(), p.A.T.copy()
spt = j5core.demotion_spaces(p, Ct, Dt)
print("trivial M zero:", not spt["M"].any(), "delta_bar(triv):", spt["delta_bar"], "(expect 0)")
# W = rowspace HX check
kerM = spt["kerM"]
Wv = (kerM @ p.HX) % 2
from qec_research.gf2.linalg import rank_np
print("rank W(triv):", rank_np(Wv), "rank HX:", p.rHX, " V'== full ker?")

# light stabilizers vs d_Z(P) = 12 (certified)
lights = j5core.light_stabilizers(p, 12)
print("lights<12:", [ (r['w'], len(r['lam'])) for r in lights][:20])
rng = np.random.default_rng(1)
probes, cands = j5core.sharp_probe(p, lights, 12, rng, max_lams=8, samples_per=4)
print("probes", probes[:8])
print("cands", cands[:3])
