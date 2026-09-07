# -*- coding: utf-8 -*-
"""Gate B: sensitivity sweeps over the QLOPS model.

Axis 1 -- decoder latency: multiply every published reaction time t_r by
m in {0.5, 0.75, 1.0, 1.5, 2.0}; recompute every Eq.(2) output at fixed
k, d, t_SEC.  Falsifier (per pre_statement.md / README gate B): any output
shifting > 2x vs the baseline (m=1) falsifies cross-paper comparability
along this axis.

Axis 2 -- magic-state protocol swap: 7-T/Litinski 15-to-1 baseline vs
zero-level CCZ constants (arXiv:2605.21867, [REPORTED] until reproduced).
Compared on
  (a) feasibility under the paper's own rule p_out <= p0 for each Table-6
      target at physical error p in {1e-3, 1e-4};
  (b) the frozen legacy arithmetic: unit_qubits*plain_cycles (Litinski) vs
      22*24 (zero-level). The latter mixes syndrome-cycle and circuit-layer
      time units; it is preserved only to reproduce the frozen Gate-B artifact.
      zero_level_provenance.py Revision 4 supplies the common-unit comparison.

Run: python3 src/gate_b_sweep.py <out_json>
"""

import json
import math
import re
import sys
from decimal import Decimal as D
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import paper_data as P  # noqa: E402
import formulas as F    # noqa: E402

MULTIPLIERS = P.GATE_B["latency_multipliers"]
FALSIFIER = P.GATE_B["comparability_falsifier"]


def ceil_ratio(a, b):
    return math.ceil(a / b)


def q_at(k, d, t_sec, t_r):
    r = ceil_ratio(D(t_r), D(t_sec))
    return int(k) / ((r + int(d)) * float(D(t_sec)))


SC_DIST = {  # Table-5 published surface-code distance per (code, variant)
    ("72,12,6-Z", "cur"): 13, ("72,12,6-Z", "fut"): 5,
    ("90,8,10-Z", "cur"): 17, ("90,8,10-Z", "fut"): 7,
    ("108,8,10-Z", "cur"): 19, ("108,8,10-Z", "fut"): 7,
    ("144,12,12-Z", "cur"): 19, ("144,12,12-Z", "fut"): 7,
    ("288,12,18-Z", "cur"): 27, ("288,12,18-Z", "fut"): 11,
    ("72,12,6-ALL", "cur"): 15, ("72,12,6-ALL", "fut"): 7,
}


# --------------------------------------------------------------------------
# Axis 1: latency sweep over all Eq.(2) inputs in the paper.
# --------------------------------------------------------------------------
def latency_axis():
    rows = []
    worst = 0.0
    for code, t5 in P.T5.items():
        n, k, d = t5["code"]
        t4 = P.GB_TABLE4[code]
        base = q_at(k, d, t4["t_sec"], t4["tr"])
        for m in MULTIPLIERS:
            q = q_at(k, d, t4["t_sec"], D(str(m)) * D(t4["tr"]))
            ratio = q / base
            worst = max(worst, abs(ratio - 1.0))
            rows.append(dict(target=f"NA {code}", m=m,
                             baseline=base, value=q, ratio=ratio))
    for code, t5 in P.T5.items():
        n, k, _ = t5["code"]
        for variant, tsec, trtab in (
                ("cur", P.T_SEC_SC["current"], P.T_R_SURFACE_CURRENT),
                ("fut", P.T_SEC_SC["future"], P.T_R_SURFACE_FUTURE)):
            d = SC_DIST[(code, variant)]
            tr = trtab.get(d)
            if tr is None:
                continue
            base = q_at(k, d, tsec, tr)
            for m in MULTIPLIERS:
                q = q_at(k, d, tsec, D(str(m)) * D(tr))
                ratio = q / base
                worst = max(worst, abs(ratio - 1.0))
                rows.append(dict(target=f"SC-{variant} {code}", m=m,
                                 baseline=base, value=q, ratio=ratio))
    # RSA systems
    k_sc = P.RSA_SC["n_cold"] + P.RSA_SC["n_active"]
    base = q_at(k_sc, P.RSA_SC["d"], P.RSA_SC["t_sec"], P.RSA_SC["t_r"])
    for m in MULTIPLIERS:
        q = q_at(k_sc, P.RSA_SC["d"], P.RSA_SC["t_sec"],
                 D(str(m)) * D(P.RSA_SC["t_r"]))
        ratio = q / base
        worst = max(worst, abs(ratio - 1.0))
        rows.append(dict(target="RSA-SC", m=m, baseline=base, value=q,
                         ratio=ratio))
    base_na = F.qlops_ignoring_latency(P.RSA_NA["n_logical"],
                                       P.RSA_NA["t_sec"])
    for m in MULTIPLIERS:
        rows.append(dict(target="RSA-NA", m=m, baseline=base_na,
                         value=base_na, ratio=1.0))  # latency ignored here
    return rows, worst


# --------------------------------------------------------------------------
# Axis 2: magic-state protocol swap.
# --------------------------------------------------------------------------
def zero_level_feasibility(p0, p_phys):
    """p_L = c*p^2 <= p0 ?  Returns (p_L_zero, feasible)."""
    z = P.ZERO_LEVEL_CCZ
    pL = z["c"] * p_phys ** z["p_order"]
    return pL, pL <= p0


def litinski_unit_spacetime(dX, dZ, dm):
    return (F.litinski_15to1_unit_qubits(dX, dZ, dm)
            * F.litinski_15to1_plain_cycles(dm))


def msd_axis():
    """Sweep Table-6 rows: 15-to-1 (baseline) vs zero-level CCZ [REPORTED]."""
    rows = []
    z = P.ZERO_LEVEL_CCZ
    for row in P.T6:
        m = re.match(r"15to1_(\d+),(\d+),(\d+)", row["proto"])
        dX, dZ, dm = map(int, m.groups())
        spacetime_base = litinski_unit_spacetime(dX, dZ, dm)
        spacetime_zero = z["legacy_gate_b_spacetime_mixed_units"]
        for p_phys, ptag in ((1e-3, "p=1e-3"), (1e-4, "p=1e-4")):
            pL_zero, feas = zero_level_feasibility(row["p0"], p_phys)
            rows.append(dict(
                target=f"{row['code']}/{row['proto']}", p_phys=ptag,
                p0=row["p0"], p_dist_15to1=row["p_dist"],
                pL_zero_level=pL_zero, zero_level_feasible=feas,
                spacetime_15to1_per_T=spacetime_base,
                spacetime_zero_per_CCZ=spacetime_zero,
                spacetime_15to1_per_CCZ_equiv=7 * spacetime_base,
                spacetime_ratio_zero_vs_litinski_ccz=(
                    spacetime_zero / (7 * spacetime_base)),
                note="zero-level constants [REPORTED] arXiv:2605.21867; "
                     "T->CCZ scale factor 7 arbitrary-but-stated"))
    infeasible_by_p = {}
    for ptag in ("p=1e-3", "p=1e-4"):
        infeasible_by_p[ptag] = sorted(
            {r["target"] for r in rows
             if r["p_phys"] == ptag and not r["zero_level_feasible"]})
    lo = min(r["spacetime_ratio_zero_vs_litinski_ccz"] for r in rows)
    hi = max(r["spacetime_ratio_zero_vs_litinski_ccz"] for r in rows)
    return rows, dict(spacetime_ratio_range=[lo, hi],
                      zero_level_infeasible_by_p=infeasible_by_p)


def main():
    lat_rows, worst_latency_shift = latency_axis()
    msd_rows, msd_summary = msd_axis()

    falsified_latency = worst_latency_shift > FALSIFIER
    infeas_1e4 = msd_summary["zero_level_infeasible_by_p"]["p=1e-4"]
    falsified_msd = bool(infeas_1e4)

    out = dict(
        axis_latency=dict(rows=lat_rows, max_shift=worst_latency_shift,
                          falsified=falsified_latency),
        axis_msd_swap=dict(rows=msd_rows, summary=msd_summary,
                           falsified=falsified_msd,
                           constants_tag=P.ZERO_LEVEL_CCZ["tag"]),
        falsifier=FALSIFIER,
    )
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path("gate_b_results.json")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")

    print(f"gate B latency axis: {len(lat_rows)} evals; "
          f"max |Q(m*tr)/Q(tr)-1| = {worst_latency_shift:.4f}; "
          f"falsified(>2x): {falsified_latency}")
    print(f"gate B msd axis: {len(msd_rows)} evals; "
          f"space-time ratio range {msd_summary['spacetime_ratio_range']}; "
          f"falsified: {falsified_msd}")
    if infeas_1e4:
        print(f"  zero-level infeasible (p_L>p0) at p=1e-4 for "
              f"{len(infeas_1e4)} of {len(P.T6)} targets; see artifact")
    print(f"artifact: {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
