"""Stage (b) rung 1 — slack-absorption certified enclosure (Addendum A3).

Quantities:
 - omega_pub = 2.371339 (published display of rung 1, SODA'25 Alman et al.)
 - Schonhage line: R_total + M*omega >= 2^(L-1) ln(q+2), L=3, q=5
 - slack absorption: build the certified-feasible retained sum R_sum by
   charging P_true = (H^max - H(rho)) enhanced by the Lemma-1 dual slack 2*eps
   (from shipped Lagrange multipliers), then report
       omega_cert_upper = (2^(L-1) ln(q+2) - R_sum_low) / M_low
   together with eps = the explicit relaxation size.

Per Main's directive the feasibility defects folded in are:
     c_max  = 1.137e-10 (inequality violation at the float point)
     ceq_max = 2.390e-11 (equality/Lagrange violation)
     schonhage = -6.2e-15 (half-ulp float noise on the Schonhage line)
"""

import sys, json, hashlib
import numpy as np
from scipy.io import loadmat
from flint import arb

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import interval_core as IC
import alman25_float as F

MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/alman_code/"
       "data/W1.00_2.371339.mat")
OMEGA_PUB = 2.371339
OMEGA_PRIOR_RUNG = 2.371552   # VXXZ24 rung for context comparison


def rung1_stage_b():
    with IC.prec():
        pmF = F.ParamManager()
        wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
        params = np.asarray(loadmat(MAT)['params']).flatten()
        pmF.set_value(params)
        cF, ceF, vF = wsF.evaluate()
        vF_f = float(np.asarray(vF).reshape(-1)[0])
        cF_max = max(float(np.asarray(z).reshape(-1)[0]) for z in cF)
        ceF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ceF)

        tree = IC.IntervalTree(MAT)
        ws = tree.ws
        c, ce, value = ws.evaluate()
        gm = ws.globstage
        st, sz = pmF.start, pmF.size

        def vecf(gid):
            return pmF.cur_x[st[gid]:st[gid] + sz[gid]]

        report = {"blocks": {}}

        # ---- per-block containment: num_block --------------------------------
        jm = gm.j2m.mats
        n_ok = 0
        for r in range(6):
            dist_f = vecf(gm.dist_id[r]).astype(float)
            for d in range(3):
                marg_f = dist_f @ jm[d]
                H_f = float(F.entropy_vec(marg_f))
                marg_iv = []
                M = jm[d]
                for j in range(M.shape[1]):
                    idxs = np.nonzero(M[:, j])[0]
                    s = IC.const(0.0)
                    for i in idxs:
                        s = s + gm.dist[r].v[i]
                    marg_iv.append(s)
                H_iv = IC.ent_vec(marg_iv)
                assert H_iv.contains(H_f), (r, d, H_iv, H_f)
                n_ok += 1
        report["blocks"]["num_block_containment"] = f"{n_ok}/18"

        # ---- penalty: H(dm) - H(dist) ----------------------------------------
        p_ok = 0
        for r in range(6):
            Hdm_f = float(F.entropy_vec(vecf(gm.dist_max_id[r])))
            Hsd_f = float(F.entropy_vec(vecf(gm.dist_id[r])))
            Hdm_iv = IC.ent_vec(gm.dist_max[r].v)
            Hsd_iv = IC.ent_vec(gm.dist[r].v)
            P_f = Hdm_f - Hsd_f
            P_iv = Hdm_iv - Hsd_iv
            tol = 1e-12 * max(1.0, abs(P_f))
            assert float(P_iv.lo) - tol <= P_f <= float(P_iv.hi) + tol
            p_ok += 1
        report["blocks"]["penalty_containment"] = f"{p_ok}/6"

        # ---- Lemma-1 eps per region ------------------------------------------
        eps_list = []
        for r in range(6):
            dm = gm.dist_max[r].v
            eps = IC.const(0.0)
            for i, shp in enumerate(gm.shapes):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = gm.lam_sum[r].v[0].lo - 1
                for d in range(3):
                    g_val = g_val + gm.lam_margin[(r, d)].v[shp[d]].lo
                diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
                e = diff.pos()
                eps = IC.Ival(min(eps.lo, e.lo), max(eps.hi, e.hi))
            eps_list.append(eps)
        report["lemma1_eps_upper_per_region"] = [float(e.hi) for e in eps_list]
        eps_max = max(float(e.hi) for e in eps_list)
        report["lemma1_eps_max"] = eps_max

        # ---- ceq containment --------------------------------------------------
        ceq_f = [float(np.asarray(z).reshape(-1)[0]) for z in ceF]
        okc = 0
        for k, z in enumerate(ce):
            if isinstance(z, IC.Ival):
                assert z.contains(ceq_f[k])
                okc += 1
            elif isinstance(z, (tuple, list)):
                lo, hi = float(min(z)), float(max(z))
                assert lo - 1e-60 <= ceq_f[k] <= hi + 1e-60
                okc += 1
        report["blocks"]["ceq_containment"] = f"{okc}"

        # ---- value containment ------------------------------------------------
        v_lo, v_hi = (float(value.lo), float(value.hi)) if hasattr(value, 'lo') else value
        tol = 1e-12 * max(1.0, abs(vF_f))
        assert v_lo - tol <= vF_f <= v_hi + tol
        report["blocks"]["value_containment"] = "OK"
        report["value_published_point_interval"] = [v_lo, v_hi]
        report["stage_a_value_float64"] = vF_f

        # ---- Schonhage slack absorption ---------------------------------------
        # In interval mode ws.evaluate() already accounts for:
        #   - P_true from the interval entropies and Lemma-1 eps (via
        #     hash_penalty_term), which reduces the certified retained sum;
        #   - the Schonhage line sv = target - value (interval) as c_viol[-1].
        # The certified-feasible Omega at the released point is
        #     Omega_cert = (target - R_sum) / M
        # recomputed from the interval bounds.

        # the interval value already includes M*omega with omega as loaded
        # (2.3713389 dyadic).  To get Omega_cert as the *bound consistent with
        # the SlackA3 statement*: solve for Omega from
        #   R_sum + M*Omega >= target,  with R_sum from the certified lower,
        #   M from the certified-feasible maximum (lower end of M interval).
        # Rebuild: value = R_sum + M*omega_loaded; so R_sum = value - M*omega.
        # We recompute M (mat_size min) interval, and R_sum = value - M*omega.
        Om_loaded = IC.const(2.3713389005434182)  # the released float
        ms_iv = gm.mat_size
        # M is an endpoint-triple (single_mat_size uses min of three); from
        # the float stage, single_mat_size = min(ms)/K.
        # We take M_low = the minimum over the three interval LOWERS.
        def _iv_lo(z):
            if isinstance(z, IC.Ival):
                return float(z.lo)
            if isinstance(z, tuple):
                return float(min(z))
            return float(np.asarray(z).reshape(-1)[0])
        M_low = min(_iv_lo(ms_iv[t]) for t in range(3))
        # value interval:
        if hasattr(value, 'lo'):
            val_lo = float(value.lo)
        elif isinstance(value, tuple):
            val_lo = float(value[0])
        else:
            val_lo = float(value)
        R_sum_low = val_lo - M_low * 2.3713389005434182
        # sanity (Addendum A3 R2): retained sum cannot be negative
        assert R_sum_low > 1.0, f"R_sum_low = {R_sum_low} <= 1: impossible"
        report["M_low"] = M_low
        report["R_sum_low"] = R_sum_low

        # certified Omega:
        target = 4.0 * np.log(7.0)
        Om_cert_upper = (target - R_sum_low) / M_low
        Om_cert_lower = (target - (val_lo + (v_hi - v_lo)) ) / (M_low + 0.0)
        report["omega_cert_upper"] = Om_cert_upper
        report["omega_cert_lower"] = Om_cert_lower
        report["target_4ln7"] = target

        # feasibility of the published point: is lhs >= target at omega_pub?
        lhs_pub_low = R_sum_low + M_low * OMEGA_PUB
        feasible_at_pub = lhs_pub_low >= target
        report["lhs_pub_low"] = lhs_pub_low
        report["feasible_at_published_omega"] = bool(feasible_at_pub)
        report["omega_published"] = OMEGA_PUB
        report["eps_abs"] = max(0.0, Om_cert_upper - OMEGA_PUB)
        report["improves_over_vxxz24"] = bool(Om_cert_upper < OMEGA_PRIOR_RUNG)
        report["prior_rung"] = OMEGA_PRIOR_RUNG

        # feasibility defects from stage (a)
        report["stage_a_c_max"] = cF_max
        report["stage_a_ceq_max"] = ceF_max
        report["stage_a_schonhage_violation"] = vF_f - target  # expected ~ -6e-15
        return report


if __name__ == "__main__":
    rep = rung1_stage_b()
    print(json.dumps(rep, indent=1, default=str))
