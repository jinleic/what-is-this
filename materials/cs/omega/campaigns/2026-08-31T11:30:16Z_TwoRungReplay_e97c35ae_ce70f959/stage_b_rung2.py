"""Stage (b) rung 2 — VXXZ24 K100_2.37155181 certified enclosure.

Follows the certified rung-1 protocol (campaign 2026-08-30T03:20:00Z_c9d4e2a8,
src/stage_b_rung1_certified.py + CORRECTION_coincidence.md), adapted to the
3-region single-p_comp program (src/vxxz24_float.py). Per pre_statement
Addendum A4 the sequence is:
  1. stage-(a) float64 reference recomputed here (must match standalone run);
  2. interval tree evaluation under endpoint-pair semantics (interval_core);
  3. DOUBLE-COUNT TEST: raw certified endpoint with and without the Lemma-1
     charge — if the endpoint moves by ~ the charge, the gap is already
     inside the raw endpoint and must not be added again;
  4. certified-feasible retained sum R_sum (interval lows; per-dim bounds
     num_block - penalty for t != r, num_block - p_comp for t == r), M at
     interval low, Omega_cert = (target_hi - R_sum)/M_low;
  5. explicit absorbed feasibility slack eps (stage-(a) defects / M_low),
     reported separately from the Lemma-1 charge; per-group Lemma-1 residual
     eps upper bounds reported next to every certified quantity.

Per-block containment is asserted inside run() before any aggregate is
trusted (num_block 9 blocks, penalty 3, ceq n, value line 1).

Run: OMP_NUM_THREADS=1 nice -n 10 python src/stage_b_rung2.py [without]
argv[1] == "without" toggles OFF the Lemma-1 charge (double-count test).
"""

import sys, json, hashlib
import numpy as np
from scipy.io import loadmat
from flint import arb

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/campaigns/2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959")
import interval_core as IC
import vxxz24_float as F

MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959/K100_2.37155181.mat")
OMEGA_PUB = 2.37155181
OMEGA_RUNG1_CERT = 2.3713400836689   # corrected rung-1 endpoint
FLOAT_ULPS = 64   # documented containment tolerance for stage-(a) noise


def ulp_tol(x):
    return FLOAT_ULPS * 2.0 ** -52 * max(1.0, abs(x))


def contains_safe(iv, x):
    """Containment with documented float-noise tolerance (R3 fix): the
    interval encloses the TRUE value at the exact dyadic inputs; the
    float64 reference carries summation noise of a few ulps (observed
    max 1.74 ulp on the K100 entropy blocks). Transcription errors would
    appear at >= 1e-6 relative scale, i.e. 9+ orders above this tolerance."""
    return float(iv.lo) - ulp_tol(x) <= x <= float(iv.hi) + ulp_tol(x)

def as_iv(x):
    """Coerce float/np/arb/G scalar into an Ival endpoint pair."""
    if x.__class__.__name__ == "Ival":
        return x
    if isinstance(x, arb):
        return IC.Ival(x.lower(), x.upper())
    return IC.const(float(np.asarray(x).reshape(-1)[0]))


def run(include_lemma1=True):
    with IC.prec():
        # ---- stage (a): float reference -------------------------------
        pmF = F.ParamManager()
        wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
        params = np.asarray(loadmat(MAT)['params']).flatten()
        sha = hashlib.sha256(params.tobytes()).hexdigest()
        pmF.set_value(params)
        cF, ceF, vF = wsF.evaluate()
        cF_max = max(float(np.asarray(z).reshape(-1)[0]) for z in cF)
        ceF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ceF)
        value_f = float(np.asarray(vF).reshape(-1)[0])
        omega_file = wsF.omega_value()
        st, sz = pmF.start, pmF.size

        def vecf(gid):
            return pmF.cur_x[st[gid]:st[gid] + sz[gid]]

        # ---- stage (b): interval tree ---------------------------------
        tree = IC.IntervalTree(MAT, float_module=F)
        ws = tree.ws
        gm = ws.globstage
        L = ws.max_level
        c, ce, value_iv = ws.evaluate()
        target = arb(7).log() * 4   # 2^(L-1) ln(q+2) for L=3, q=5

        report = {
            "rung": "vxxz24_2.37155181",
            "params_sha256": sha,
            "n_params": int(len(params)),
            "n_registered": int(pmF.num_input),
            "stage_a_float_value": value_f,
            "stage_a_omega_in_file": omega_file,
            "stage_a_c_max": cF_max,
            "stage_a_ceq_max": ceF_max,
            "target_4ln7": float(target),
        }

        # ---- containment (A3 R3) ---------------------------------------
        blocks = {}
        n_ok = 0
        for r in range(3):
            dist_f = vecf(gm.dist_id[r]).astype(float)
            marg_ivs = gm.j2m.apply(gm.dist[r].v)
            for d in range(3):
                H_f = float(F.entropy_vec(dist_f @ gm.j2m.mats[d]))
                H_iv = IC.ent_vec(marg_ivs[d])
                # R3 with documented float noise: exact-dyadic entropy vs
                # float64 reference differs by ~1.7 ulp (pairwise sum).
                assert contains_safe(H_iv, H_f), \
                    ("num_block", r, d, str(H_iv), H_f)
                n_ok += 1
        blocks["global_num_block"] = f"{n_ok}/9"
        p_ok = 0
        for r in range(3):
            Hdm_f = float(F.entropy_vec(vecf(gm.dist_max_id[r])))
            Hsd_f = float(F.entropy_vec(vecf(gm.dist_id[r])))
            P_f = Hdm_f - Hsd_f
            P_iv = IC.ent_vec(gm.dist_max[r].v) - IC.ent_vec(gm.dist[r].v)
            assert contains_safe(P_iv, P_f), ("penalty", r, str(P_iv), P_f)
            p_ok += 1
        blocks["global_penalty"] = f"{p_ok}/3"
        # ceq (Lagrange) consistency: the interval path emits LOG-FORM pairs
        # (lo,hi) = ln(dmv)-(tv-1); the float path emits EXP-FORM residuals
        # exp(tv-1)-dmv. Consistency is |exp-form| <= |log-form| * dmv with
        # dmv <= 1 — verified pairwise per constraint below (form-to-form,
        # not raw comparison which would be a category error).
        nceq = 0
        for k, z in enumerate(ce):
            if isinstance(z, (tuple, list)) and len(z) == 2:
                fv = float(np.asarray(ceF[k]).reshape(-1)[0])
                Lw = max(abs(z[0]), abs(z[1]))
                # |exp(tv-1)-dmv| = dmv*|e^{-L}-1| <= |L| (dmv<=1, small |L|)
                assert abs(fv) <= Lw * 1.000001 + 1e-15, \
                    ("ceq", k, z, fv)
                nceq += 1
        blocks["ceq_logform_consistency"] = f"{nceq}"
        v_lo, v_hi = float(value_iv.lo), float(value_iv.hi)
        assert contains_safe(IC.Ival(arb(v_lo), arb(v_hi)), value_f), \
            ("value", v_lo, value_f, v_hi)
        report["blocks"] = blocks
        report["value_interval"] = [v_lo, v_hi]
        report["value_width"] = v_hi - v_lo

        # ---- certified aggregation -------------------------------------
        details = {}
        lemma_eps = {}

        def gm_penalty_true(r, include):
            Hdm = gm.dist_max[r].entropy().v
            Hsd = gm.dist[r].entropy().v
            prop = gm.region_prop[r].item()
            if not include:
                return prop * (Hdm - Hsd), None
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
            lemma_eps[f"glob_r{r}"] = float(eps_iv.hi)
            return prop * (Hdm - Hsd + 2 * eps_iv), eps_iv

        def part_penalty_true(t, r, include):
            Hdm = t.split_dist_max[r].entropy().v
            Hsd = t.split_dist[r].entropy().v
            frac = t.part_frac.item() * t.region_prop[r].item()
            if not (include and isinstance(t, F.Part)):
                return frac * (Hdm - Hsd), None
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
                eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
            key = f"part{t.part_id}_r{r}"
            if (key not in lemma_eps) or float(eps_iv.hi) > lemma_eps[key]:
                lemma_eps[key] = float(eps_iv.hi)
            return frac * (Hdm - Hsd + 2 * eps_iv), eps_iv

        total_R = arb(0)

        # level >= 3 per (r, l)
        for l in range(3, L + 1):
            # accumulate once per region: nb over all parts; pen, pc over Parts
            for r in range(3):
                nb = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
                pen = IC.Ival(arb(0), arb(0))
                pc = IC.Ival(arb(0), arb(0))
                for t in ws.parts[l - 1]:
                    for tt in range(3):
                        val = t.num_block_contribution[r][tt]
                        nb[tt] = nb[tt] + as_iv(val)
                    if isinstance(t, F.Part):
                        pt, _ = part_penalty_true(t, r, include_lemma1)
                        pen = pen + as_iv(pt)
                        pc = pc + as_iv(t.p_comp[r].item())
                cand = [nb[tt].lo - (pc.hi if tt == r else pen.hi)
                        for tt in range(3)]
                R_low = min(cand)
                if float(R_low) > 0:
                    total_R = total_R + R_low
                details[f"R_comp[{r},{l}]"] = float(R_low)

        # level-2 symmetric: only region 0 group
        nb2 = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
        for t in ws.parts[1]:
            for tt in range(3):
                nb2[tt] = nb2[tt] + as_iv(t.num_block_contribution[tt])
        R2_low = min(z.lo for z in nb2)
        if float(R2_low) > 0:
            total_R = total_R + R2_low
        details["R_comp[0,2]"] = float(R2_low)

        # global hashing per region
        for r in range(3):
            nbv = gm.num_block[r]
            pen_iv, _ = gm_penalty_true(r, include_lemma1)
            pcm = gm.p_comp[r].item()
            pc_iv = as_iv(pcm)
            cand = [nbv[tt].lo - (pc_iv.hi if tt == r else pen_iv.hi)
                    for tt in range(3)]
            R_low = min(cand)
            if float(R_low) > 0:
                total_R = total_R + R_low
            details[f"R_glob[{r}]"] = float(R_low)

        # matrix size + Omega
        ms_iv = gm.mat_size
        M_low = min(float(z.lo) for z in ms_iv)
        M_hi_at_min = min(float(z.hi) for z in ms_iv)

        tgt_hi = float(target.upper()) if hasattr(target, "upper") \
            else float(target)
        Om_raw = (tgt_hi - float(total_R.lower())) / M_low

        eps_abs = (cF_max + ceF_max) / M_low \
            + max(0.0, tgt_hi - value_f) / M_low

        report.update({
            "include_lemma1": include_lemma1,
            "R_sum_low_certified": float(total_R.lower()),
            "R_detail": details,
            "M_low_certified": M_low,
            "M_hi_at_min_dim": M_hi_at_min,
            "omega_cert_upper_raw": Om_raw,
            "lemma1_eps_by_group": lemma_eps,
            "lemma1_eps_max": (max(lemma_eps.values()) if lemma_eps else 0.0),
            "eps_abs_defects_only": eps_abs,
            "omega_with_absorbed_defects": Om_raw + eps_abs,
            "omega_published": OMEGA_PUB,
            "rung1_cert": OMEGA_RUNG1_CERT,
            "two_rung_gap": Om_raw + eps_abs - OMEGA_RUNG1_CERT,
        })
        return report


if __name__ == "__main__":
    include_lemma1 = (len(sys.argv) < 2) or (sys.argv[1] != "without")
    rep = run(include_lemma1=include_lemma1)
    print(json.dumps(rep, indent=1, default=str))
