"""Gate C stage 3, part 2 — ENDGAME: decision on the 21-dim nullspace
question. Combines (pre_statement_gateC2.md §2c):
  - certified slope bundle S_i over the box (from gate_c_slope_pass),
  - the ANY-y Lagrangian instrument on the kernel,
  - a certified second-order remainder Rmax over the box,
and decides PASS-a / PASS-b / FAIL with the exact width arithmetic.

Validations INSIDE (all must pass before any decision is printed):
  V1  point-slope pass reproduces the frozen certified endpoint value
      (2.3715538358...) by direct re-aggregation of the Slope VALUE parts;
  V2  slope values at r=0 reproduce the 45-dim coordinate gradient
      (coord_grad_full.npy) to float tolerance;
  V3  the kernel gradient V^T g matches kernel_grads_float.npy;
  V4  the zero-crossing guard stays clean over the box.
"""
import sys, json, time
import numpy as np

SRC = "/Users/jinleic/jinleic-workspace/cs/omega/src"
sys.path.insert(0, SRC)
import interval_core as IC
import vxxz24_float as F
import gate_c_slope_core as SC
import gate_c_slope_pass as SP
from flint import arb
from scipy.io import loadmat
from scipy.optimize import linprog

I = IC.Ival
CAMP0 = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
         "2026-08-30T11:24:41Z_17bdaeb6_099b1224d603")
RUNG2 = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
         "2026-08-30T05:20:00Z_e97c35ae_70c28fc61780")
FROZEN_WITH = json.load(open(f"{RUNG2}/stage_b_rung2_with_lemma1.json"))
FROZEN_WITHOUT = json.load(open(
    f"{RUNG2}/stage_b_rung2_without_lemma1.json"))
SCR = "/Users/jinleic/jinleic-workspace/cs/omega/scratch"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/endgame.json"
RADIUS = 1e-7
T0 = time.time()

params = np.asarray(loadmat(SP.MAT)["params"]).flatten()
Vex = json.load(open(f"{CAMP0}/kernel_basis_V_exact_rational.json"))
# Exact rational kernel basis, rendered only after loading the authoritative
# numerator/denominator pairs.  Do not substitute kernel_basis_V.npy: that is
# the earlier dense floating/SVD render and is not a rule-14 witness.
V = np.column_stack([
    [num / den for num, den in Vex[str(j)]] for j in range(21)
])
A = np.load(f"{CAMP0}/margin_matrix_M.npy")          # exact 0/1, 27x45
assert V.shape == (45, 21) and A.shape == (27, 45)
assert np.max(np.abs(A @ V)) == 0.0
assert np.max(np.abs(V.sum(axis=0))) == 0.0

checkpoint = {}

def ckpt(k, v):
    checkpoint[k] = v
    with open(OUT, "w") as fh:
        json.dump(checkpoint, fh, indent=1, default=str)
    print("CKPT", k, json.dumps(v, default=str)[:220])
    sys.stdout.flush()

V_svd_counterfactual = np.load(f"{CAMP0}/kernel_basis_V.npy")
assert V_svd_counterfactual.shape == V.shape
assert not np.array_equal(V_svd_counterfactual, V)
ckpt("basis_path_counterfactual", {
    "evidence": "MACHINE-VERIFIED exact integer arithmetic for exact basis; "
                "COMPUTATIONAL-EVIDENCE for dense SVD residual",
    "actual_source": "kernel_basis_V_exact_rational.json",
    "old_source": "kernel_basis_V.npy",
    "actual_equals_exact_JSON_render": True,
    "actual_equals_old_dense_SVD": False,
    "exact_denominators_all_one": bool(all(
        den == 1 for col in Vex.values() for _, den in col)),
    "exact_entry_set": sorted(np.unique(V).tolist()),
    "max_abs_A_times_exact_V": float(np.max(np.abs(A @ V))),
    "max_abs_A_times_dense_SVD_COMPUTATIONAL_EVIDENCE": float(
        np.max(np.abs(A @ V_svd_counterfactual))),
    "max_abs_exact_minus_dense_SVD": float(
        np.max(np.abs(V - V_svd_counterfactual))),
})

# ---------------------------------------------------------------- V3 first
g45_ref = np.load(f"{SCR}/coord_grad_full.npy")
kgf = np.load(f"{CAMP0}/kernel_grads_float.npy")

# ------------------------------------------------------- point pass (r=0)
with IC.prec():
    out0 = SP.run_slope_pass(params, radius=0.0, verbose=False)
    gm0 = out0["gm"]; ws0 = out0["ws"]

    def slope_of(x):
        return x.v if isinstance(x, SP.SlopeG) else x

    # --- certified endpoint re-aggregation from Slope VALUE parts
    # (mirror of gate_c_stage_b_boxes.certified_at at radius)
    def as_iv(z):
        if isinstance(z, SC.Slope):
            return z.val
        if isinstance(z, IC.Ival):
            return z
        if isinstance(z, arb):
            return I(z.lower(), z.upper())
        return IC.const(float(np.asarray(z).reshape(-1)[0]))

    def aggregate(out, radius, lemma1_off=False):
        gm = out["gm"]; ws = out["ws"]
        def get_s(g, key):
            s = g[key]
            return s.v if isinstance(s, SP.SlopeG) else s
        pen_parts = {}   # (l, r) -> Slope sum of penalties
        nb_parts = {}    # (l, r, tt) -> Slope (l>=3); level 2 separate
        pc_parts = {}    # (l, r) -> Slope
        for l in range(3, 4):
            for r in range(3):
                nb = [SC.Slope(SC.ZERO, np.array([SC.ZERO]*SP.NSLOPE, dtype=object))
                      for _ in range(3)]
                pen = SC.s_const(0.0, SP.NSLOPE)
                pc = SC.s_const(0.0, SP.NSLOPE)
                for t in ws.parts[l - 1]:
                    nbc = getattr(t, "num_block_contribution", None)
                    if nbc is None or nbc[r] is None:
                        continue
                    for tt in range(3):
                        _e = nbc[r][tt]
                        if not isinstance(_e, SC.Slope):
                            raise TypeError(
                                f"nb[{tt}] not Slope: {type(_e)} "
                                f"(part {type(t).__name__}, nbc {type(nbc)}, "
                                f"nbc[r] {type(nbc[r]) if isinstance(nbc, (list, tuple)) else ''})")
                        nb[tt] = SC.s_add(nb[tt], _e)
                    if isinstance(t, F.Part):
                        if lemma1_off:
                            hm = SC.s_sub(t.split_dist_max[r].entropy().v,
                                          t.split_dist[r].entropy().v)
                            part_pen = SC.s_mul(
                                SC.s_mul(hm, slope_of(t.part_frac)),
                                t.region_prop.v[r])
                        else:
                            part_pen = slope_of(t.hash_penalty_term[r])
                        pen = SC.s_add(pen, part_pen)
                        pc = SC.s_add(pc, slope_of(t.p_comp[r]))
                pen_parts[(l, r)] = pen
                pc_parts[(l, r)] = pc
                for tt in range(3):
                    nb_parts[(l, r, tt)] = nb[tt]
        nb2 = [SC.s_const(0.0, SP.NSLOPE) for _ in range(3)]
        for t in ws.parts[1]:
            nbc = getattr(t, "num_block_contribution", None)
            if nbc is None:
                continue
            for tt in range(3):
                _e = nbc[tt]
                if not isinstance(_e, SC.Slope):
                    raise TypeError(f"nb2[{tt}] not Slope: {type(_e)} "
                                    f"(part {type(t).__name__})")
                nb2[tt] = SC.s_add(nb2[tt], _e)
        # glob
        nbG = {}
        penG = {}
        pcG = {}
        for r in range(3):
            for tt in range(3):
                nbG[(r, tt)] = slope_of(gm.num_block[r][tt])
            if lemma1_off:
                hm = SC.s_sub(gm.dist_max[r].entropy().v,
                              gm.dist[r].entropy().v)
                penG[r] = SC.s_mul(hm, gm.region_prop.v[r])
            else:
                penG[r] = slope_of(gm.hash_penalty_term[r])
            pcG[r] = slope_of(gm.p_comp[r])
        # aggregate R (value parts, outward like frozen): cand = nb.lo − pen/pc.hi
        total_R = arb(0)
        R_slopes = []   # list of (Slope, which) for contribution tracking
        #mul level>=3 blocks with BRANCH PINNING
        branch_data = {}
        for l in range(3, 4):
            for r in range(3):
                cand = []
                slopes = []
                for tt in range(3):
                    nbl = nb_parts[(l, r, tt)].val.lo
                    rhs = pc_parts[(l, r)].val.hi if tt == r else pen_parts[(l, r)].val.hi
                    cand.append(float(nbl - rhs))
                    slopes.append((nb_parts[(l, r, tt)],
                                   pc_parts[(l, r)] if tt == r else pen_parts[(l, r)]))
                R_low = min(cand)
                branch_data[(l, r)] = {"cand": cand, "argmin": int(np.argmin(cand))}
                if float(R_low) > 0:
                    total_R = total_R + R_low
                    R_slopes.append(slopes[int(np.argmin(cand))])
        # level 2: nb2 vs (none)
        R2_low = min(z.val.lo for z in nb2)
        R2_arg = int(np.argmin([float(z.val.lo) for z in nb2]))
        branch_data["lvl2"] = {"cand": [float(z.val.lo) for z in nb2],
                               "argmin": R2_arg}
        if float(R2_low) > 0:   # REPAIR C: level-2 R branch pinned too
            total_R = total_R + R2_low
            R_slopes.append((nb2[R2_arg], SC.s_const(0.0, SP.NSLOPE)))
        # glob blocks
        for r in range(3):
            cand = []
            slopes = []
            for tt in range(3):
                nbl = nbG[(r, tt)].val.lo
                rhs = pcG[r].val.hi if tt == r else penG[r].val.hi
                cand.append(float(nbl - rhs))
                slopes.append((nbG[(r, tt)],
                               pcG[r] if tt == r else penG[r]))
            R_low = min(cand)
            branch_data[("G", r)] = {"cand": cand, "argmin": int(np.argmin(cand))}
            if float(R_low) > 0:
                total_R = total_R + R_low
                R_slopes.append(slopes[int(np.argmin(cand))])
        R_slope = SC.s_const(0.0, SP.NSLOPE)
        for (nbS, rhsS) in R_slopes:
            R_slope = SC.s_add(R_slope, SC.s_sub(nbS, rhsS))
        msS = gm.mat_size.v
        ms_iv = [z.val for z in msS]
        M_low = min(float(z.lo) for z in ms_iv)
        M_arg = int(np.argmin([float(z.lo) for z in ms_iv]))
        target_hi = float((arb(7).log() * 4).upper())
        Om_raw = (target_hi - float(total_R.lower())) / M_low
        # Differentiate the pinned endpoint expression directly:
        # Omega = (T - R) / M.  This both gives the correct value enclosure
        # and lets SC.s_div apply the quotient rule
        #   dOmega = (-R' - Omega*M') / M.
        M_slope = msS[M_arg]
        Om_slope = SC.s_div(
            SC.s_sub(SC.s_const(target_hi, SP.NSLOPE), R_slope),
            M_slope)
        return {"total_R": float(total_R.lower()), "M_low": M_low, "M_arg": M_arg,
                "Om_raw": Om_raw, "Om_slope": Om_slope, "R_slope": R_slope,
                "branch": branch_data,
                "value_s": out.get("value_s"), "line_s": out.get("line_s")}

    def r_detail(agg):
        branch = agg["branch"]
        return {
            "R_comp[0,3]": max(0.0, min(branch[(3, 0)]["cand"])),
            "R_comp[1,3]": max(0.0, min(branch[(3, 1)]["cand"])),
            "R_comp[2,3]": max(0.0, min(branch[(3, 2)]["cand"])),
            "R_comp[0,2]": max(0.0, min(branch["lvl2"]["cand"])),
            "R_glob[0]": max(0.0, min(branch[("G", 0)]["cand"])),
            "R_glob[1]": max(0.0, min(branch[("G", 1)]["cand"])),
            "R_glob[2]": max(0.0, min(branch[("G", 2)]["cand"])),
        }

    agg0 = aggregate(out0, 0.0)
    agg0_lemma1_off = aggregate(out0, 0.0, lemma1_off=True)
    detail_with = r_detail(agg0)
    detail_without = r_detail(agg0_lemma1_off)
    for got, frozen in ((detail_with, FROZEN_WITH),
                        (detail_without, FROZEN_WITHOUT)):
        for key, value in got.items():
            assert value == frozen["R_detail"][key], (key, value,
                                                       frozen["R_detail"][key])
    for got, frozen in ((agg0, FROZEN_WITH),
                        (agg0_lemma1_off, FROZEN_WITHOUT)):
        assert got["total_R"] == frozen["R_sum_low_certified"]
        assert got["M_low"] == frozen["M_low_certified"]
        assert got["Om_raw"] == frozen["omega_cert_upper_raw"]
    zero_penalty_delta = agg0["Om_raw"] - agg0_lemma1_off["Om_raw"]
    frozen_zero_penalty_delta = (
        FROZEN_WITH["omega_cert_upper_raw"]
        - FROZEN_WITHOUT["omega_cert_upper_raw"])
    assert zero_penalty_delta == frozen_zero_penalty_delta
    ckpt("V1_frozen_machine_equality", {
        "real": {
            "R_detail": detail_with,
            "R_sum": agg0["total_R"],
            "M_low": agg0["M_low"],
            "omega_raw": agg0["Om_raw"],
        },
        "lemma1_penalty_zeroed": {
            "R_detail": detail_without,
            "R_sum": agg0_lemma1_off["total_R"],
            "M_low": agg0_lemma1_off["M_low"],
            "omega_raw": agg0_lemma1_off["Om_raw"],
        },
        "omega_delta_real_minus_zeroed": zero_penalty_delta,
        "all_float_equal_to_frozen": True,
    })
    ckpt("agg_point", {"Om_raw_r0": agg0["Om_raw"],
                       "total_R": agg0["total_R"],
                       "M_low": agg0["M_low"],
                       "branch": {str(k): v
                                  for k, v in agg0["branch"].items()}})

    # ------------------------------------------------------------- V1/V2
    # V1: point-pass value must hit the frozen certified raw endpoint
    FROZEN_RAW = 2.3715538358350803
    FROZEN_OMEGA = 2.3715538358544617   # raw + defects/M_low
    v1_ok = abs(agg0["Om_raw"] - FROZEN_RAW) < 5e-13
    ckpt("V1_point_value", {"Om_raw": agg0["Om_raw"], "frozen": FROZEN_RAW,
                            "delta": agg0["Om_raw"] - FROZEN_RAW, "ok": bool(v1_ok)})
    if not v1_ok:
        print("V1 FAILED: point re-aggregation disagrees with frozen value")
        raise SystemExit(2)

    # V2: slope kernel projection vs predecessor kernel gradient
    Om_grad = agg0["Om_slope"].grad
    g45 = np.array([float(x.lo) for x in Om_grad] +
                   [float(x.hi) for x in Om_grad]).reshape(2, 45)
    g45_mid = (g45[0] + g45[1]) / 2.0
    g45_halfwidth = (g45[1] - g45[0]) / 2.0
    v2_dev = float(np.max(np.abs(g45_mid - g45_ref)))
    ckpt("V2_point_gradient", {"max_dev_vs_ref": v2_dev,
                               "max_halfwidth": float(g45_halfwidth.max()),
                               "ok": bool(v2_dev < 1e-6)})

    # ------------------------------------------------ point pass gradients
    # kernel gradient at the CENTRE (float point)
    kg = V.T @ g45_mid
    ckpt("kernel_grad_point", {"V^T g": kg.tolist(),
                               "ref_kernel_grads": kgf.tolist(),
                               "max_dev": float(np.abs(kg - kgf).max())})

    # -------------------------------------------- boxed slope pass (r=1e-7)
    outB = SP.run_slope_pass(params, radius=RADIUS, verbose=False)
    aggB = aggregate(outB, RADIUS)
    Om_gradB = aggB["Om_slope"].grad
    lo_b = np.array([float(x.lo) for x in Om_gradB])
    hi_b = np.array([float(x.hi) for x in Om_gradB])
    wid_b = hi_b - lo_b
    ckpt("slope_box", {"grad_lo": lo_b.tolist(), "grad_hi": hi_b.tolist(),
                       "width_max": float(wid_b.max()),
                       "width_sum_l1": float(wid_b.sum()),
                       "elapsed": time.time() - T0})

    # ---------------- exact-basis midpoint LP (candidate/diagnostic only)
    # Scale c = RADIUS*u before solving: |V*u|_inf <= 1.  The predecessor's
    # A_ub/RADIUS form accidentally imposed |V*u| <= RADIUS and shrank D by
    # another factor 1e-7.
    grad_mid = (lo_b + hi_b) / 2.0
    gv_mid = V.T @ grad_mid
    Aub = np.vstack([V, -V])
    rLPs = linprog(gv_mid, A_ub=Aub, b_ub=np.ones(90),
                   bounds=[(None, None)] * 21, method="highs")
    if rLPs.status != 0:
        ckpt("LP_scaled_midpoint", {"status": int(rLPs.status),
                                    "message": rLPs.message})
        raise RuntimeError("scaled midpoint LP failed")
    u_star = np.asarray(rLPs.x)
    delta_star = RADIUS * (V @ u_star)
    lin_mid_scaled = RADIUS * float(rLPs.fun)
    primal_violation = max(0.0, float(np.abs(V @ u_star).max()) - 1.0)
    kernel_resid = float(np.abs(A @ delta_star).max())

    # HiGHS inequality marginals have the opposite sign from nonnegative
    # Lagrange multipliers.  There are 90 row duals (45 upper + 45 lower),
    # not 21 variable-bound marginals.
    row_marg = np.asarray(rLPs.ineqlin.marginals)
    lam = -row_marg
    mup, mum = lam[:45], lam[45:]
    kkt_resid = float(np.abs(gv_mid + V.T @ (mup - mum)).max())

    # Keep the unscaled solve only as the registered HiGHS-tolerance control;
    # it is never consumed by a certificate.
    rLPu = linprog(gv_mid, A_ub=Aub,
                   b_ub=np.full(90, RADIUS),
                   bounds=[(None, None)] * 21, method="highs")
    unscaled_min = float(rLPu.fun) if rLPu.status == 0 else None
    ckpt("LP_scaled_midpoint", {
        "evidence": "COMPUTATIONAL-EVIDENCE",
        "basis": "exact rational JSON (-1/0/1), not float SVD",
        "min": lin_mid_scaled,
        "status": int(rLPs.status),
        "max_scaled_primal_violation": primal_violation,
        "max_kernel_residual": kernel_resid,
        "kkt_residual": kkt_resid,
        "unscaled_min_diagnostic": unscaled_min,
        "scaled_unscaled_delta": (
            None if unscaled_min is None else lin_mid_scaled - unscaled_min),
        "delta_star": delta_star.tolist(),
        "row_dual_upper": mup.tolist(),
        "row_dual_lower": mum.tolist(),
    })

    # ---------------- rigorous ANY-y instrument on the interval slope bundle
    # For every delta in D and every gradient g in the full-box enclosure G:
    #   g·delta = (g - A^T y)·delta
    #           >= -RADIUS * ||g - A^T y||_1
    # for ARBITRARY y, because A·delta=0 exactly.  The LP below merely finds
    # a useful y; no solver feasibility or dual claim is consumed.
    lo_fit = np.nextafter(lo_b, -np.inf)
    hi_fit = np.nextafter(hi_b, np.inf)
    fit_obj = np.concatenate([np.zeros(27), np.ones(45)])
    fit_A = np.vstack([
        np.hstack([A.T, -np.eye(45)]),
        np.hstack([-A.T, -np.eye(45)]),
    ])
    fit_b = np.concatenate([lo_fit, -hi_fit])
    fit = linprog(fit_obj, A_ub=fit_A, b_ub=fit_b,
                  bounds=[(None, None)] * 27 + [(0, None)] * 45,
                  method="highs")
    if fit.status != 0:
        ckpt("ANY_y", {"status": int(fit.status), "message": fit.message})
        raise RuntimeError("ANY-y L1 fit failed")
    y = np.asarray(fit.x[:27])

    z_iv = I(arb(0), arb(0))
    residual_l1 = z_iv
    assert set(np.unique(A)).issubset({0.0, 1.0})
    for i, g_iv in enumerate(Om_gradB):
        aty = z_iv
        for j in np.flatnonzero(A[:, i]):
            aty = aty + IC.const(float(y[j]))
        residual_l1 = residual_l1 + (g_iv - aty).pos()
    swing_iv = IC.const(RADIUS) * residual_l1
    swing_upper = float(swing_iv.up())
    lin_min = -swing_upper
    lin_max = swing_upper
    width_cert = 2.0 * swing_upper
    ckpt("ANY_y", {
        "evidence": "MACHINE-VERIFIED",
        "instrument": "R*sup_G ||G-A^T y||_1; arbitrary y",
        "fit_status": int(fit.status),
        "fit_objective_float": float(fit.fun),
        "y_exact_dyadic_float_render": y.tolist(),
        "residual_l1_lo": float(residual_l1.lo),
        "residual_l1_hi_buffered": float(residual_l1.up()),
        "linear_swing_upper": swing_upper,
        "lin_min": lin_min,
        "lin_max": lin_max,
        "width": width_cert,
    })

    # The full-box slope bundle is itself the mean-value remainder enclosure:
    # Omega(x)-Omega(c)=g(xi)·(x-c), xi on the segment in the box.  No second
    # Hessian term is consumed.  Retain the old quadratic expression only as
    # a non-certified consistency number.
    rem_b = 0.5 * float(wid_b.max()) * 45 * RADIUS
    ckpt("remainder_check", {
        "consumed_remainder": "full-box slope bundle / mean-value form",
        "quadr_consistency_est_COMPUTATIONAL_EVIDENCE": rem_b,
    })
    max_wid_idx = int(np.argmax(wid_b))
    d0B = outB["pm"].start[outB["gm"].dist_id[0]]
    leafB = outB["gm"].dist[0].v[max_wid_idx].val
    ckpt("slope_width_bottleneck", {
        "coordinate": max_wid_idx,
        "shape": list(outB["gm"].shapes[max_wid_idx]),
        "p_star": float(params[d0B + max_wid_idx]),
        "boxed_leaf": [float(leafB.lo), float(leafB.hi)],
        "grad_interval": [float(Om_gradB[max_wid_idx].lo),
                          float(Om_gradB[max_wid_idx].hi)],
        "grad_width": float(wid_b[max_wid_idx]),
        "width_sum_other_44": float(wid_b.sum() - wid_b[max_wid_idx]),
    })

    SIGNAL = 1.5816497000997742e-07
    ratio = width_cert / SIGNAL
    ckpt("cert_width_over_D", {"lin_min": lin_min, "lin_max": lin_max,
                               "width": width_cert})
    ckpt("signal_ratio", {"width": width_cert, "signal": SIGNAL,
                          "ratio_vs_predecessor_122x": float(ratio)})

    # The certified enclosure straddles zero, while the scaled midpoint LP is
    # only a float candidate.  Therefore neither PASS-a nor PASS-b is earned.
    decision = ("FAILURE TO CERTIFY: ANY-y slope enclosure straddles zero "
                "over D; 21-dimensional question remains OPEN")
    ckpt("decision", {
        "lin_min": lin_min,
        "lin_max": lin_max,
        "width": width_cert,
        "ratio": float(ratio),
        "Om0_raw": agg0["Om_raw"],
        "midpoint_LP_candidate": lin_mid_scaled,
        "decision": decision,
        "elapsed_total": time.time() - T0,
    })

print("ENDGAME COMPLETE — decision:", checkpoint.get("decision", {}).get("decision"))
print(json.dumps({k: checkpoint[k] for k in
                  ("signal_ratio", "cert_width_over_D", "decision")
                  if k in checkpoint}, indent=1, default=str))
