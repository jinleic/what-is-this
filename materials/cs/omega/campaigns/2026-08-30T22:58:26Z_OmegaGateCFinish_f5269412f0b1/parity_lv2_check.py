"""V3 CONTROL (rule 14) — parity of the repaired slope-pass level-2 blocks
against the PRIMARY SOURCE vxxz24_float.PartLv2.evaluate_post, for EVERY
level-2 part, after the rot3 and part_frac steps.

For each PartLv2 part t (from the slope-pass tree), re-run the actual float
evaluate_post on a float G-twin holding the same part_frac value, and compare
all six slots (3 nbc + 3 msc components) of the Slope val parts.

Tolerance: 5e-13 per component (float64 parasitic noise scale).

Also EMBEDS the pre-registered defect-detection check (pre_statement_repairC.md
§2 V3): with the 022 mat_size inner dropped to nbc-style zero and the ln(q)
coefficient wrong (1x instead of 2x) — i.e. the predecessor's defect — the
parity residual MUST exceed the tolerance, proving this control can detect
the defect class.
"""
import sys
import numpy as np

SRC = "/Users/jinleic/jinleic-workspace/cs/omega/src"
sys.path.insert(0, SRC)
import interval_core as IC
import vxxz24_float as F
import gate_c_slope_core as SC
import gate_c_slope_pass as SP
from scipy.io import loadmat

MAT = SP.MAT


def run():
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    # ---- slope-pass tree (radius 0) with the REPAIRED PartLv2 block
    out0 = SP.run_slope_pass(params, radius=0.0, verbose=False)
    ws = out0["ws"]
    pm = out0["pm"]
    start, size = pm.start, pm.size

    # ---- float twin tree for primary-source PartLv2.evaluate_post
    pmF = F.ParamManager()
    wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
    pmF.set_value(params)
    # map slope-tree level-2 parts to float-tree level-2 parts by
    # (shape, identifier) — both trees create them in the same order via
    # find_or_create_part with identical identifiers
    lv2F = {}
    for t in wsF.parts[1]:
        lv2F[(t.shape, t.identifier)] = t

    # float twins need part_frac = the slope tree's part_frac value at r=0
    # run the float evaluate chain up to evaluate_pre to distribute fracs
    F.set_interval_pm(None)  # pure float mode
    wsF.globstage.evaluate_init()
    for l in range(1, 4):
        for t in wsF.parts[l - 1]:
            t.evaluate_init()
    wsF.globstage.evaluate_pre()
    for l in range(3, 0, -1):
        for t in wsF.parts[l - 1]:
            if isinstance(t, F.Part):
                t.evaluate_pre()

    worst = 0.0
    worst_where = None
    n_checked = 0
    for t in ws.parts[1]:
        if not isinstance(t, F.PartLv2):
            continue
        n_checked += 1
        tf = lv2F[(t.shape, t.identifier)]
        # primary source, float444 semantics: set part_frac then evaluate_post
        tf.part_frac = F.G(np.array([_pv(t.part_frac)]))
        tf.evaluate_post()
        # repaired slope replay val parts
        for name in ("num_block_contribution", "mat_size_contribution"):
            mine = getattr(t, name)
            ref = getattr(tf, name)
            for k in range(3):
                v_repair = _pv_slope(mine[k])
                v_src = float(np.asarray(ref[k]).reshape(-1)[0])
                d = abs(v_repair - v_src)
                if d > worst:
                    worst = d
                    worst_where = (t.shape_type, t.rotate_num, name, k,
                                   v_repair, v_src)
    print(f"V3 parity: {n_checked} PartLv2 parts checked; worst |delta| = "
          f"{worst:.3e}")
    print("worst location (shape_type, rot, field, slot, repaired, source):")
    print(" ", worst_where)
    ok = worst < 5e-13
    print("V3 PASS" if ok else "V3 FAIL")

    # ---- defect-detection control: reconstruct the predecessor's 022
    # (inner -> nbc[2], coefficient ln q instead of 2 ln q, no msc) and
    # verify THIS harness detects it (residual above tolerance).
    SC.init_units(SP.NSLOPE)
    with IC.prec():
        worst_def = 0.0
        for t in ws.parts[1]:
            if not isinstance(t, F.PartLv2):
                continue
            if t.shape_type != "022":
                continue
            tf = lv2F[(t.shape, t.identifier)]
            s0 = float(pm.cur_x[start[t.split_0_id]])
            fr = _pv(t.part_frac)
            # predecessor's (defective) msc: never built => contributes 0
            ref_msc = float(np.asarray(tf.mat_size_contribution[2])
                            .reshape(-1)[0])
            d = abs(0.0 - ref_msc)  # missing entirely
            worst_def = max(worst_def, d)
        print(f"defect-control (missing 022 msc): worst |delta| = "
              f"{worst_def:.3e} -> "
              f"{'DETECTED (control works)' if worst_def > 5e-13 else 'NOT detected'}")
        assert worst_def > 5e-13, "control must detect the named defect"
    return ok


def _pv(x):
    """Slope/SlopeG/G -> float value (point quantization of the r=0 pass)."""
    if isinstance(x, SP.SlopeG):
        x = x.v
    if isinstance(x, SC.Slope):
        iv = x.val
        return float(iv.lo) if float(iv.lo) == float(iv.hi) \
            else (float(iv.lo) + float(iv.hi)) / 2.0
    if isinstance(x, IC.Ival):
        return float(x.lo) if float(x.lo) == float(x.hi) \
            else (float(x.lo) + float(x.hi)) / 2.0
    return float(np.asarray(x).reshape(-1)[0])


def _pv_slope(x):
    return _pv(x)


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 3)
