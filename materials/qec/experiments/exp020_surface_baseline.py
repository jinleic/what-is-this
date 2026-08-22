
import sys; sys.path.insert(0,"src")
import json, platform
from pathlib import Path
from qec_research.codes.bicycle import BRAVYI_BB, BBSpec, build_bb, bb_stabilizer
from qec_research.gf2.linalg import rank_bitset, rows_to_bitsets

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "processed"

def sc_phys_per_logical(d):
    return 2*d*d - 1

# every qLDPC code: n_data, n_checks (= rows of H = one ancilla per check),
# and whether a verified FT circuit was actually built for it
codes = []
# codes with verified circuits (built in exp004/exp007)
codes.append({"name":"[[72,12,6]]", "spec":BRAVYI_BB["[[72,12,6]]"],
              "circuit_built": True, "source": "exp004/exp007 verified circuit"})
codes.append({"name":"[[108,8,10]]", "spec":BRAVYI_BB["[[108,8,10]]"],
              "circuit_built": True, "source": "exp004 verified circuit"})
codes.append({"name":"[[144,12,12]] Gross", "spec":BRAVYI_BB["[[144,12,12]]"],
              "circuit_built": True, "source": "exp004/exp007 verified, Gate 7 passed"})

# codes WITHOUT circuits: compute ancilla count from stabilizer structure
codes.append({"name":"[[288,12,18]]",
              "spec":BBSpec(ell=12,m=12,A=[(3,0),(0,2),(0,7)],B=[(0,3),(1,0),(2,0)]),
              "circuit_built": False, "source": "stabilizer-row count, no circuit built"})
codes.append({"name":"[[288,16,12]]",
              "spec":BBSpec(ell=12,m=12,A=[(0,0),(0,1),(0,2)],B=[(0,0),(2,0),(4,0)]),
              "circuit_built": False, "source": "stabilizer-row count, no circuit built"})

comparisons = []
for c in codes:
    sp = c["spec"]
    HX, HZ = build_bb(sp)
    n_data = sp.n
    n_checks = HX.shape[0] + HZ.shape[0]   # one ancilla per stabilizer row
    total = n_data + n_checks
    code = bb_stabilizer(sp)
    d = None
    # use published distance where available
    pub_d = {"[[72,12,6]]":6, "[[108,8,10]]":10, "[[144,12,12]] Gross":12,
             "[[288,12,18]]":18, "[[288,16,12]]":12}
    d = pub_d[c["name"]]
    k = code.k
    sc = k * sc_phys_per_logical(d)
    comparisons.append({
        "qldpc_code": c["name"], "k": k, "d": d,
        "n_data": n_data, "n_ancilla": n_checks,
        "qldpc_total_physical": total,
        "sc_total_physical": sc,
        "physical_overhead_ratio": round(sc / total, 2),
        "circuit_built": c["circuit_built"],
        "provenance": c["source"],
    })

MATCHING_RULES = {
    "1_same_k": "PERFORMED",
    "2_same_target_pL": "UNPERFORMED — no simulation to a target logical error probability",
    "3_comparable_d": "PERFORMED — same code distance",
    "4_comparable_phys": "RESULT — the overhead ratio is the output",
    "5_cycle_duration": "UNPERFORMED at matched p_L",
}
rec = {
    "experiment": "exp020_surface_baseline",
    "description": ("First-order same-(k,d) physical-qubit comparison. qLDPC totals "
                    "include data + syndrome ancillas (one per stabilizer row)."),
    "method": ("Rotated SC: 2d^2-1 per logical. qLDPC: n_data + n_checks where "
               "n_checks = H_X.rows + H_Z.rows (one ancilla per stabilizer generator). "
               "Codes marked circuit_built=False have no verified FT circuit; their "
               "ancilla count is the stabilizer-row count, not a measured circuit property."),
    "matching_rules": MATCHING_RULES,
    "comparisons": comparisons,
    "explicitly_NOT_claimed": [
        "No target-p_L space-time comparison (rule 2 unperformed).",
        "No cycle-duration comparison at matched p_L (rule 5 unperformed).",
        "No circuit-level simulation of the surface code.",
        "For circuit_built=False codes, ancilla count is structural, not from a built circuit.",
    ],
}
(OUT / "exp020_surface_baseline.json").write_text(json.dumps(rec, indent=2))

print(f"{'code':24s} {'k':>3s} {'d':>3s} {'data':>5s} {'anc':>5s} {'total':>5s} "
      f"{'SC':>6s} {'overhead':>8s} {'circuit':>7s}")
print("-"*72)
for c in comparisons:
    print(f"{c['qldpc_code']:24s} {c['k']:3d} {c['d']:3d} {c['n_data']:5d} "
          f"{c['n_ancilla']:5d} {c['qldpc_total_physical']:5d} {c['sc_total_physical']:6d} "
          f"{c['physical_overhead_ratio']:7.1f}x {'yes' if c['circuit_built'] else 'NO':>7s}")
print(f"\nwrote results/processed/exp020_surface_baseline.json")
