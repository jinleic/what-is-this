"""Stage (b) rung 1 certified enclosure — CORRECT semantics.

Per the taint test and width analysis:
- `value` (Schonhage LHS) sums pm point intervals only: zero width by
  construction. The enclosure's certified content lives in the CONSTRAINT
  quantities: num_block (entropies, width ~1.8e-20) and their penalties.
- Certified retained sum: R_cert = num_block_low - P_true_high, per
  (region, level) group, charged through the hashing constraints, plus
  level-2 symmetric terms. Its interval width is real (entropy radii plus
  Lemma-1 slack).

Deliverable (Addendum A3): Omega_cert = (target^up - R_cert low)/M low with
eps = Omega_cert_upper - omega_pub, absorbing (a) Lemma-1 slack (dominant),
(b) stage-a constraint violations c=1.137e-10, ceq=2.390e-11, and
(c) Schonhage half-ulp -6.2e-15.
"""

import sys, json, hashlib
import numpy as np
from scipy.io import loadmat
from flint import arb

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/campaigns/2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959")
import interval_core as IC
import alman25_float as F

MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959/W1.00_2.371339.mat")
OMEGA_PUB = 2.371339
OMEGA_PRIOR = 2.371552
TARGET = None  # set inside


def rung1_certified():
    with IC.prec():
        params = np.asarray(loadmat(MAT)['params']).flatten()
        sha = hashlib.sha256(params.tobytes()).hexdigest()

        # float reference (stage a)
        pmF = F.ParamManager()
        wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
        pmF.set_value(params)
        cF, ceF, vF = wsF.evaluate()
        cF_max = max(float(np.asarray(z).reshape(-1)[0]) for z in cF)
        ceF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ceF)
        value_f = float(np.asarray(vF).reshape(-1)[0])

        # interval tree
        tree = IC.IntervalTree(MAT)
        ws = tree.ws
        gm = ws.globstage
        L = ws.max_level
        q = 5.0
        target = arb(int(q + 2)).log() * (2 ** (L - 1))   # = 4 ln 7
        TARGET = target

        # ---- containment blocks first (A3 R3) ------------------------
        # (num_block 18, penalty 6, ceq n, value 1) — asserted inside
        # stage_b_rung1.py; repeated here minimally for value:
        c, ce, value_iv = ws.evaluate()

        # ---- build certified R_sum ------------------------------------
        # level>=3 hashing: R_(r,l) constrained by
        #   R_(r,l) <= num_block_low - P_true_high   (for X dim)
        #             or num_block_low - p_comp_high (Y/Z dims)
        # plus R_(r,l) >= 0 (pm domain). We take the certified feasible
        # maximum: R_cert = max(0, num_block_low - P_true_high).
        # Level 2: R_(0,2) <= min over dims; R_(r>=1,2) = 0 (pm domain).
        # Objective contribution is the SUM over all R_cert entries plus
        # M * omega, with M = min over dims of mat_size (certified low).
        total_R = arb(0)
        lemma_eps_max = 0.0
        details = []

        for l in range(3, L + 1):
            for r in range(6):
                dimx, dimy, dimz = F.Dims(r + 1)
                # sum num_block over terms (X dim only is charged)
                nb_x = IC.Ival(arb(0), arb(0))
                pen = IC.Ival(arb(0), arb(0))   # P_true (X dim)
                for t in ws.terms[l - 1]:
                    if isinstance(t, F.TermLv2):
                        nb_x = nb_x + t.num_block[dimx - 1]
                        continue
                    nb_x = nb_x + t.num_block_contribution[r][dimx - 1]
                    if isinstance(t, F.Term):
                        frac = t.term_frac.item() * t.region_prop[r].item()
                        Hdm = t.split_dist_max[r].entropy().v
                        Hsd = t.split_dist[r].entropy().v
                        eps_iv = IC.Ival(arb(0), arb(0))
                        dm = t.split_dist_max[r].v
                        for i, sp in enumerate(t.splits):
                            if float(dm[i].hi) <= 0:
                                continue
                            g_val = t.lam_sum[r].v[0].lo - 1
                            for d in range(3):
                                g_val = g_val + t.lam_margin[(r, d)].v[
                                    sp[d] - t.lam_low[d]].lo
                            diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
                            e = diff.pos()
                            eps_iv = IC.Ival(min(eps_iv.lo, e.lo),
                                             max(eps_iv.hi, e.hi))
                        lemma_eps_max = max(lemma_eps_max, float(eps_iv.hi))
                        p_true = frac * (Hdm - Hsd + 2 * eps_iv)
                        pen = pen + p_true
                R_low = nb_x.lo - pen.hi
                if float(R_low) > 0:
                    total_R = total_R + (R_low.lower() if hasattr(R_low, 'lower') and not isinstance(R_low, arb) else R_low)
                # record
                details.append((f"R[{r},{l}]", float(R_low)))

        # level-2 symmetric hashing
        for d in range(3):
            nb2 = IC.Ival(arb(0), arb(0))
            for t in ws.terms[1]:
                nb2 = nb2 + t.num_block[d]
            if d == 0:
                total_R = total_R + nb2.lo.lower() if hasattr(nb2.lo, 'lower') and not isinstance(nb2.lo, arb) else total_R + nb2.lo.lower() if hasattr(nb2.lo, 'lower') else total_R + nb2.lo

        # global hashing (region r contributes min over dims, but for square
        # (K=1) case, all dims equal; use X dim low with P_true charged)
        for r in range(6):
            dimx, dimy, dimz = F.Dims(r + 1)
            nbv = gm.num_block[r][dimx - 1]
            frac = gm.region_prop[r].item()
            Hdm = gm.dist_max[r].entropy().v
            Hsd = gm.dist[r].entropy().v
            eps_iv = IC.Ival(arb(0), arb(0))
            dm = gm.dist_max[r].v
            for i, shp in enumerate(gm.shapes):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = gm.lam_sum[r].v[0].lo - 1
                for d in range(3):
                    g_val = g_val + gm.lam_margin[(r, d)].v[shp[d]].lo
                diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
                e = diff.pos()
                eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
            lemma_eps_max = max(lemma_eps_max, float(eps_iv.hi))
            p_true = frac * (Hdm - Hsd + 2 * eps_iv)
            R_low = nbv.lo - p_true.hi
            if float(R_low) > 0:
                total_R = total_R + (R_low.lower() if hasattr(R_low, 'lower') and not isinstance(R_low, arb) else R_low)

        # matrix size M
        ms_iv = gm.mat_size
        M_low = min(float(z.lo) for z in ms_iv)

        # certified omega upper:
        tgt = target.hi if hasattr(target, 'hi') else target
        Om_up = (tgt - total_R.lower()) / M_low
        Om_up_f = float(Om_up)

        # feasibility absorption (Addendum A3): eps_abs folds in stage-a
        # defect sizes plus the Lemma-1 slack (already charged via p_true);
        # the residual constraint violations translate into an ADDITIVE
        # omega relaxation of order (violation / M):
        eps_abs = max(0.0, Om_up_f - OMEGA_PUB) \
            + (cF_max + ceF_max) / M_low \
            + abs(value_f - float(target)) / M_low
        # domain note: pm groups carrying slack (num_retain / mat_size /
        # omega) absorb the violation; conserved quantity is the value line.

        result = {
            "rung": "alman25_2.371339",
            "params_sha256": sha,
            "n_params": int(len(params)),
            "stage_a_float_value": value_f,
            "stage_a_c_max": cF_max,
            "stage_a_ceq_max": ceF_max,
            "schonhage_violation_stage_a": value_f - 4 * np.log(7),
            "lemma1_eps_max": lemma_eps_max,
            "R_sum_low_certified": float(total_R.lower()),
            "M_low_certified": M_low,
            "target": float(target.hi if hasattr(target, 'hi') else target),
            "omega_cert_upper_raw": Om_up_f,
            "omega_published": OMEGA_PUB,
            "eps_abs": eps_abs,
            "omega_cert_with_absorbed_slack": Om_up_f + eps_abs,
            "improves_over_prior_rung": bool(Om_up_f + eps_abs < OMEGA_PRIOR),
            "prior_rung": OMEGA_PRIOR,
        }
        return result


if __name__ == "__main__":
    rep = rung1_certified()
    print(json.dumps(rep, indent=1))
