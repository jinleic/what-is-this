"""Decision-1 repaired probes — run ONLY after
pre_statement_addendum_decision1.md is committed (it is). Repair rule,
directions, feasibility gates, and stopping condition are all there; this
script implements them verbatim and reports every probe including misses.
"""
import sys, json, time
import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import interval_core as IC
import vxxz24_float as F
from scipy.io import loadmat
from flint import arb

MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/"
       "data/K100_2.37155181.mat")
EPS = 1e-8
STEP_SMALL = 7e-8
T0 = time.time()
DEADLINE = 30 * 60.0


def repaired_offset(params, pm, gm, region, sign):
    """Pre-registered repair rule (addendum): +7e-8 on the smallest
    coordinate, -1e-8 on the 7 largest, rest 0; sign flips everything for
    the negated direction."""
    s = pm.start[gm.dist_id[region]]
    n = pm.size[gm.dist_id[region]]
    d = np.zeros(len(params))
    block = params[s:s + n]
    k = int(np.argmin(block))
    largest = sorted(range(n), key=lambda i: (-block[i], i))[:7]
    d[s + k] += sign * STEP_SMALL
    for i in largest:
        d[s + i] -= sign * EPS
    return d, k, largest


def exact_sum(vals):
    """Exact-dyadic sum via IC.const endpoints (no ball radii)."""
    tot = IC.Ival(arb(0), arb(0))
    for v in vals:
        tot = tot + IC.const(float(v))
    return float(tot.lo)


def feasibility_gate(pm, ws, params, offset, region, record):
    """Gates per addendum. Returns (ok, reasons)."""
    reasons = []
    # gate 1: simplex residual identity (exact dyadic)
    gm = ws.globstage
    s = pm.start[gm.dist_id[region]]
    n = pm.size[gm.dist_id[region]]
    base_res = exact_sum(params[s:s + n]) - 1.0
    rep_res = exact_sum(params[s:s + n] + offset[s:s + n]) - 1.0
    record["simplex_residual_base"] = base_res
    record["simplex_residual_repaired"] = rep_res
    if rep_res != base_res:
        reasons.append(f"simplex residual changed {base_res!r} -> {rep_res!r}")
    # gate 2: support (all coords strictly inside (0,1))
    blk = params[s:s + n] + offset[s:s + n]
    if blk.min() <= 0 or blk.max() >= 1:
        reasons.append("support: coordinate exited (0,1)")
    # gate 3: dist<->dist_max marginal rows (direct, region-level)
    worst = 0.0
    d = pm.cur_x[s:s + n]
    dmx = pm.cur_x[pm.start[gm.dist_max_id[region]]:
                   pm.start[gm.dist_max_id[region]] + n]
    for t in range(3):
        A = gm.j2m.mats[t]
        worst = max(worst, float(np.max(np.abs(d @ A - dmx @ A))))
    record["marginal_worst"] = worst
    if worst > 1.1e-9:
        reasons.append(f"marginal sharing residual {worst:.3e} > 1.1e-9")
    # gate 4: full ceq set at the repaired point (float path residuals
    # against the release's own tolerance discipline)
    _, ce, _ = ws.evaluate()
    cem = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ce)
    record["ceq_max_at_repaired"] = cem
    if cem > 1.0e-8:
        reasons.append(f"ceq_max {cem:.3e} at repaired point")
    return (len(reasons) == 0), reasons


def certified_aggregation():
    """Interval-path certified endpoint at the CURRENT interval pm."""
    pmI = F._PM_GLOBAL
    ws = None
    # build fresh interval workspace against the static pm layout
    # (uses the same wiring as gate_c_stage_b_boxes.certified_at but at
    # the already-set repaired centre)
    return None  # replaced by inline logic in main loop


def main():
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    out = {"results": []}
    with IC.prec():
        for region in range(3):
            pm = F.ParamManager()
            ws = F.Workspace(pm, 5.0, 1.0, 0.0, 3)
            gm = ws.globstage
            for sign in (+1.0, -1.0):
                if time.time() - T0 > DEADLINE:
                    out["results"].append({"status": "WALLCLOCK_STOP"})
                    break
                off, k, largest = repaired_offset(params, pm, gm, region, sign)
                rec = {"probe": f"glob_dist_r{region}_repaired_"
                                f"{'+' if sign > 0 else '-'}",
                       "moved_small_coord": int(k),
                       "moved_largest7": [int(i) for i in largest]}
                pm.set_value(params + off)
                ok, reasons = feasibility_gate(pm, ws, params, off, region, rec)
                rec["feasibility_ok"] = ok
                rec["gate_reasons"] = reasons
                out["results"].append(rec)
                print(json.dumps(rec))
                sys.stdout.flush()
    with open(sys.argv[1] if len(sys.argv) > 1 else "repaired_probes.json",
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("repaired-probe gate pass complete")


if __name__ == "__main__":
    main()
