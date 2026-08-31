"""Part 1 — decompose the 2.0258545e-6 rung-2 gap by precision.

Re-runs the SAME interval code path (stage_b_rung2.run) at several working
precisions and reports: certified endpoint, gap to published 2.37155181,
interval widths (R_sum, M, value). If the gap shrinks with precision it is
outward-rounding accumulation; if flat it is structural (Lemma-1 dual
residual of the shipped multipliers is precision-independent).
"""
import sys, json, gc
import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import interval_core as IC
import vxxz24_float as F

PRECISIONS = [300, 200, 128, 96, 64]
OMEGA_PUB = 2.37155181


def run_at_prec(m_bits, include_lemma1=True):
    IC.MID = m_bits
    with IC.prec():
        rep = _run_bare(include_lemma1)
    return rep


def _run_bare(include_lemma1):
    """stage_b_rung2.run body, using IC.MID set above (skips float refeps)."""
    import hashlib
    from scipy.io import loadmat
    from flint import arb
    MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/"
           "data/K100_2.37155181.mat")
    pmF = F.ParamManager()
    wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
    params = np.asarray(loadmat(MAT)['params']).flatten()
    pmF.set_value(params)
    cF, ceF, vF = wsF.evaluate()
    cF_max = max(float(np.asarray(z).reshape(-1)[0]) for z in cF)
    ceF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ceF)
    value_f = float(np.asarray(vF).reshape(-1)[0])

    tree = IC.IntervalTree(MAT, float_module=F)
    ws = tree.ws
    gm = ws.globstage
    L = ws.max_level
    c, ce, value_iv = ws.evaluate()
    target = arb(7).log() * 4

    # containment asserts first (lightweight: only num_block + value + pen)
    for r in range(3):
        dist_f = pmF.cur_x[pmF.start[gm.dist_id[r]]:
                           pmF.start[gm.dist_id[r]] + pmF.size[gm.dist_id[r]]]
        marg_ivs = gm.j2m.apply(gm.dist[r].v)
        for d in range(3):
            H_f = float(F.entropy_vec(dist_f.astype(float) @ gm.j2m.mats[d]))
            H_iv = IC.ent_vec(marg_ivs[d])
            assert float(H_iv.lo) - 1e-11 <= H_f <= float(H_iv.hi) + 1e-11, \
                ("num_block", r, d)
    v_lo, v_hi = float(value_iv.lo), float(value_iv.hi)
    assert v_lo - 1e-11 <= value_f <= v_hi + 1e-11, "value"

    def gm_penalty_true(r, include):
        Hdm = gm.dist_max[r].entropy().v
        Hsd = gm.dist[r].entropy().v
        prop = gm.region_prop[r].item()
        if not include:
            return prop * (Hdm - Hsd)
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
        return prop * (Hdm - Hsd + 2 * eps_iv)

    def part_penalty_true(t, r, include):
        Hdm = t.split_dist_max[r].entropy().v
        Hsd = t.split_dist[r].entropy().v
        frac = t.part_frac.item() * t.region_prop[r].item()
        if not (include and isinstance(t, F.Part)):
            return frac * (Hdm - Hsd)
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
        return frac * (Hdm - Hsd + 2 * eps_iv)

    total_R = arb(0)
    R_width_sum = 0.0
    for l in range(3, L + 1):
        for r in range(3):
            nb = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
            pen = IC.Ival(arb(0), arb(0))
            pc = IC.Ival(arb(0), arb(0))
            for t in ws.parts[l - 1]:
                for tt in range(3):
                    val = t.num_block_contribution[r][tt]
                    nb[tt] = nb[tt] + (val if val.__class__.__name__ == "Ival"
                                       else IC.const(float(np.asarray(val).reshape(-1)[0])))
                if isinstance(t, F.Part):
                    pt = part_penalty_true(t, r, include_lemma1)
                    pen = pen + (pt if pt.__class__.__name__ == "Ival"
                                 else IC.const(float(pt)))
                    pv = t.p_comp[r].item()
                    pc = pc + (pv if pv.__class__.__name__ == "Ival"
                               else IC.const(float(np.asarray(pv).reshape(-1)[0])))
            cand = [nb[tt].lo - (pc.hi if tt == r else pen.hi)
                    for tt in range(3)]
            R_low = min(cand)
            if float(R_low) > 0:
                total_R = total_R + R_low

    nb2 = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
    for t in ws.parts[1]:
        for tt in range(3):
            v = t.num_block_contribution[tt]
            nb2[tt] = nb2[tt] + (v if v.__class__.__name__ == "Ival"
                                 else IC.const(float(np.asarray(v).reshape(-1)[0])))
    R2_low = min(z.lo for z in nb2)
    if float(R2_low) > 0:
        total_R = total_R + R2_low

    for r in range(3):
        nbv = gm.num_block[r]
        pen_iv = gm_penalty_true(r, include_lemma1)
        pcm = gm.p_comp[r].item()
        pc_iv = (pcm if pcm.__class__.__name__ == "Ival"
                 else IC.const(float(np.asarray(pcm).reshape(-1)[0])))
        cand = [nbv[tt].lo - (pc_iv.hi if tt == r else pen_iv.hi)
                for tt in range(3)]
        R_low = min(cand)
        if float(R_low) > 0:
            total_R = total_R + R_low

    ms_iv = gm.mat_size
    M_low = min(float(z.lo) for z in ms_iv)
    M_hi_at_min = min(float(z.hi) for z in ms_iv)

    tgt_hi = float(target.upper())
    Om_raw = (tgt_hi - float(total_R.lower())) / M_low
    eps_abs = (cF_max + ceF_max) / M_low + max(0.0, tgt_hi - value_f) / M_low

    return {
        "prec_bits": IC.MID,
        "omega_cert_upper_raw": Om_raw,
        "omega_with_absorbed_defects": Om_raw + eps_abs,
        "eps_abs": eps_abs,
        "R_sum_low": float(total_R.lower()),
        "R_sum_hi": float(total_R.upper()),
        "R_sum_width": float(total_R.upper() - total_R.lower()),
        "M_low": M_low,
        "M_width": M_hi_at_min - M_low,
        "value_width": float(value_iv.hi - value_iv.lo),
        "include_lemma1": include_lemma1,
        "cF_max": cF_max,
        "ceF_max": ceF_max,
    }


if __name__ == "__main__":
    arm = (len(sys.argv) > 1 and sys.argv[1] == "without")
    out = {"arm": "without_lemma1" if arm else "with_lemma1", "sweep": []}
    for bits in PRECISIONS:
        rep = run_at_prec(bits, include_lemma1=not arm)
        rep["gap_to_published"] = rep["omega_with_absorbed_defects"] - OMEGA_PUB
        out["sweep"].append(rep)
        print(json.dumps(rep))
        sys.stdout.flush()
        gc.collect()
    with open(sys.argv[2] if len(sys.argv) > 2 else
              "precision_sweep.json", "w") as fh:
        json.dump(out, fh, indent=1)
