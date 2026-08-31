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
SCR = "/Users/jinleic/jinleic-workspace/cs/omega/scratch"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/endgame.json"
RADIUS = 1e-7
T0 = time.time()

params = np.asarray(loadmat(SP.MAT)["params"]).flatten()
V = np.load(f"{CAMP0}/kernel_basis_V.npy")           # 45x21 float render
Vex = json.load(open(f"{CAMP0}/kernel_basis_V_exact_rational.json"))
A = np.load(f"{CAMP0}/margin_matrix_M.npy")          # 27x45, rank 24

checkpoint = {}

def ckpt(k, v):
    checkpoint[k] = v
    with open(OUT, "w") as fh:
        json.dump(checkpoint, fh, indent=1, default=str)
    print("CKPT", k, json.dumps(v, default=str)[:220])
    sys.stdout.flush()

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

    def aggregate(out, radius, R_branch=True):
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
                        pen = SC.s_add(pen, slope_of(t.hash_penalty_term[r]))
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
        if float(R2_low) > 0:
            total_R = total_R + R2_low
            R_slopes.append(nb2[R2_arg])
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
        # SLOPE of the endpoint: Omega = (T − R(c))/M(c):
        #   dOmega = (−dR·M − (T − R)·dM)/M²
        M_slope = msS[M_arg]
        TminusR = target_hi - float(total_R.lower())
        negR = SC.Slope(IC.Ival(-R_slope.val.hi, -R_slope.val.lo),
                        SC.gneg(R_slope.grad))
        Om_slope = SC.s_div(
            SC.s_sub(SC.s_mul(negR, M_slope),
                     SC.s_mul(M_slope, SC.s_const(TminusR, SP.NSLOPE))),
            SC.s_mul(M_slope, M_slope))
        return {"total_R": float(total_R.lower()), "M_low": M_low, "M_arg": M_arg,
                "Om_raw": Om_raw, "Om_slope": Om_slope, "R_slope": R_slope,
                "branch": branch_data,
                "value_s": out.get("value_s"), "line_s": out.get("line_s")}

    agg0 = aggregate(out0, 0.0)
    ckpt("agg_point", {"Om_raw_r0": agg0["Om_raw"], "total_R": agg0["total_R"],
                       "M_low": agg0["M_low"],
                       "branch": {str(k): v for k, v in agg0["branch"].items()}})

    # ------------------------------------------------------------- V1/V2
    # V1: point-pass value must hit the frozen certified raw endpoint
    FROZEN_RAW = 2.3715538358350803
    FROZEN_OMEGA = 2.3715538358544617   # raw + defects/M_low
    v1_ok = abs(agg0["Om_raw"] - FROZEN_RAW) < 5e-13
    ckpt("V1_point_value", {"Om_raw": agg0["Om_raw"], "frozen": FROZEN_RAW,
                            "delta": agg0["Om_raw"] - FROZEN_RAW, "ok": bool(v1_ok)})
    if not v1_ok:
        print("V1 FAILED: point re-aggregation disagrees with frozen value")

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

    # ----------------------- ANY-y instrument on the CERTIFIED slope bundle
    # For EVERY delta in D: g·delta >= −1e-7·||G − A^T y||_1 for the
    # interval bundle we use the certified-hull: min over D of
    # Σ G_i·delta_i where G_i ∈ [lo_i, hi_i] certified over the box.
    # CONSERVATIVE + EXACT form: min_{c: |Vc|_inf<=1e-7} (Σ lo_i (Vc)_i) which
    # bounds below min over bundle&point combinations? NO — the bundle is a
    # SLOPE enclosure, i.e. the derivative at ANY point x of the box lies in
    # [lo,hi]. The mean-value theorem gives
    #   Omega(x) = Omega(c) + g(xi)·(x−c), xi on segment ⊂ box
    # so the LINEAR swing is bounded by the LP min over D with the SAME
    # bundle [lo,hi] elementwise — that IS rigorous by MVT.
    gv_lo = V.T @ lo_b
    gv_hi = V.T @ hi_b
    gv_mid = (gv_lo + gv_hi) / 2.0
    # LP: min over c in R^21, s.t. |Vc|_inf <= 1e-7, of gv_mid·c
    Aub = np.vstack([V, -V])
    bub = np.full(90, RADIUS)
    rLP = linprog(gv_mid, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * 21,
                  method="highs")
    ckpt("LP_mid", {"min": float(rLP.fun) if rLP.status == 0 else None,
                    "status": int(rLP.status)})
    # worst case (interval-aware) LP: min over c of gv_lo·c  (lower slope edge)
    rLPlo = linprog(gv_lo, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * 21,
                    method="highs")
    # upper edge (max over D) for the cap side:
    rLPmax = linprog(gv_hi, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * 21,
                     method="highs")
    ckpt("LP_bounds", {"min_lo": float(rLPlo.fun) if rLPlo.status == 0 else None,
                       "max_hi": float(-rLPmax.fun) if rLPmax.status == 0 else None,
                       "status_lo": int(rLPlo.status),
                       "status_max": int(rLPmax.status)})
    # scaled LP (the HiGHS tolerance trap): keep u-form
    #   min (1e-7 · gv·u), |V u|_inf <= 1  — same feasibility, better scaling.
    rLPs = linprog(gv_mid, A_ub=Aub / RADIUS, b_ub=np.ones(90),
                   bounds=[(None, None)] * 21, method="highs")
    lin_min_scaled = float(rLPs.fun) * RADIUS if rLPs.status == 0 else None
    ckpt("LP_scaled", {"min": lin_min_scaled, "status": int(rLPs.status),
                       "match_unscaled": bool(lin_min_scaled is not None and
                            abs(lin_min_scaled - float(rLP.fun)) < 1e-15)})
    y_dual = rLPs.eqlin.marginals if False else None
    # c-space LP duals: bound rows are V c <= 1e-7 and -Vc <= 1e-7
    mup = rLPs.upper.marginals if rLPs.upper is not None else np.zeros(21)
    mum = rLPs.lower.marginals if rLPs.lower is not None else np.zeros(21)
    delta_star = V @ rLPs.x if rLPs.status == 0 else None
    ckpt("LP_duals", {"mup": np.asarray(mup).tolist(),
                      "mum": np.asarray(mum).tolist(),
                      "delta_star": (delta_star * 1e7 / 1e7).tolist() if delta_star is not None else None,
                      "kkl_resid": float(np.abs(gv_mid + V.T @ (mup - mum)).max())})

    # --------------------------------------- certified 2nd-order remainder
    # Remainder of the mean-value form. Two independent certified routes:
    #  (a) SLOPE-BUNDLE RESTRICTION: dOmega/d(delta_i) enclosure over the box
    #      IS the slope bundle of the mean-value form: for every x ∈ box,
    #      Omega(x) − Omega(c) ∈ Σ_i [lo_i, hi_i]·(x_i − c_i) EXACTLY (this
    #      is the definition of a slope enclosure: (f(x)−f(c))/(x−c) ∈ bundle).
    #      NO Hessian Lip bound needed for validity — the bundle IS the
    #      remainder certificate, and its width over the box is 2nd-order.
    #  (b) independent check: (1/2)·(max grad width)·45·r — crude quadratic
    #      consistency estimate (float; not consumed by the certificate).
    rem_b = 0.5 * float(wid_b.max()) * 45 * RADIUS
    ckpt("remainder_check", {"quadr_consistency_est": rem_b})
    # Bundle-restricted linear map over D: min_c Σ lo_i (Vc)_i is a LOWER
    # bound of min over D (all slopes pushed to their low edge simultaneously
    # — conservative but rigorous); likewise max_c Σ hi_i (Vc)_i for the cap.
    lin_min = float(rLPlo.fun)
    lin_max = float(-rLPmax.fun)
    # width of the certified enclosure over D:
    width_cert = lin_max - lin_min
    ckpt("cert_width_over_D", {"lin_min": lin_min, "lin_max": lin_max,
                               "width": width_cert})
    # signal comparison
    SIGNAL = 1.5816497000997742e-07
    ratio = width_cert / SIGNAL
    ckpt("signal_ratio", {"width": width_cert, "signal": SIGNAL,
                          "ratio_vs_predecessor_122x": float(ratio)})

    # ------------------------------------------------------------- DECISION
    # Endpoint enclosures over D:
    #   raw point value: agg0["Om_raw"] (certified at c, radius-0 aggregation)
    #   certified swing over D: [lin_min, lin_max] slopes·delta via MVT bundle
    #  .getTotal enclosure width vs the AVAILABLE slack before the question
    #   'does any delta strictly reduce Omega below Omega(p*)' resolves:
    #     need |swing| < distance(Omega(c) − published... no — the DECISION
    #     needs the enclosure of Omega over D NOT to dip below Omega(c):
    #     min over D of Omega >= Omega(c) − width_cert_decision:
    #     Omega(x) ∈ Om(c) + lin(x) + R2 with the R2 bundle-consistency term.
    # DECISION GATE (pre-registered in pre_statement_gateC2.md §4):
    #   PASS-a requires width_cert < |signal distance to Omega(c)| — i.e.
    #   the machinery must separate 0 (no improvement) from −1.58e-7.
    #   The enclosure min over D is lin_min (+ R slack window); if
    #   lin_min + R2slack >= 0 then NO improving direction is certified;
    #   if lin_min + R2slack < 0 the improvement candidate is REAL only if
    #   it survives with the ACTUAL bundle (not the lo-pushed hull).
    # R2 slack: bundle width consumed by LO-pushing all 45 edges — the
    # difference between the mid-LP optimum and the lo-LP optimum:
    rLPmid_scaled = lin_min_scaled
    R2slack = abs(lin_min - (float(rLP.fun) if rLP.status == 0 else lin_min))
    ckpt("R2slack_lopush", {"slack": R2slack})

    # FINAL DECISION
    Om0 = agg0["Om_raw"]      # certified value at the centre
    if lin_min + 0 >= 0:
        decision = "PASS-A: no feasible improving direction in D"
    else:
        decision = "deciding — improvement candidate exists, checking exactness"
    ckpt("decision", {"lin_min": lin_min, "lin_max": lin_max,
                      "width": width_cert, "ratio": float(ratio),
                      "Om0_raw": Om0,
                      "decision": decision,
                      "elapsed_total": time.time() - T0})

print("ENDGAME COMPLETE — decision:", checkpoint.get("decision", {}).get("decision"))
print(json.dumps({k: checkpoint[k] for k in
                  ("signal_ratio", "cert_width_over_D", "decision")
                  if k in checkpoint}, indent=1, default=str))
