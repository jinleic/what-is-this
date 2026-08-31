
import sys, json
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import numpy as np
from scipy.io import loadmat
from flint import arb
import interval_core as IC
import alman25_float as F
MAT = "/Users/jinleic/jinleic-workspace/cs/omega/scratch/alman_code/data/W1.00_2.371339.mat"
INCLUDE = (sys.argv[1] == "with")
with IC.prec():
    tree = IC.IntervalTree(MAT)
    ws = tree.ws
    gm = ws.globstage
    ws.evaluate()    # full run populates terms' post results
    L = ws.max_level
    total_R = arb(0)
    for l in range(3, L + 1):
        for r in range(6):
            dimx, dimy, dimz = F.Dims(r + 1)
            nb_x = IC.Ival(arb(0), arb(0))
            pen = IC.Ival(arb(0), arb(0))
            for t in ws.terms[l - 1]:
                if isinstance(t, F.TermLv2):
                    nb_x = nb_x + t.num_block[dimx - 1]
                    continue
                elif isinstance(t, F.TermZero):
                    continue
                nb_x = nb_x + t.num_block_contribution[r][dimx - 1]
                if isinstance(t, F.Term):
                    frac = t.term_frac.item() * t.region_prop[r].item()
                    Hdm = t.split_dist_max[r].entropy().v
                    Hsd = t.split_dist[r].entropy().v
                    if INCLUDE:
                        eps_iv = IC.Ival(arb(0), arb(0))
                        dm = t.split_dist_max[r].v
                        for i, sp in enumerate(t.splits):
                            if float(dm[i].hi) <= 0: continue
                            g_val = t.lam_sum[r].v[0].lo - 1
                            for d in range(3):
                                g_val = g_val + t.lam_margin[(r, d)].v[sp[d] - t.lam_low[d]].lo
                            diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
                            e = diff.pos()
                            eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
                        p_true = frac * (Hdm - Hsd + 2 * eps_iv)
                    else:
                        p_true = frac * (Hdm - Hsd)
                    pen = pen + p_true
            R_low = nb_x.lo - pen.hi
            if float(R_low) > 0:
                total_R = total_R + R_low
    for d in [0]:
        nb2 = IC.Ival(arb(0), arb(0))
        for t in ws.terms[1]:
            nb2 = nb2 + t.num_block[d]
        total_R = total_R + nb2.lo
    for r in range(6):
        dimx, dimy, dimz = F.Dims(r + 1)
        nbv = gm.num_block[r][dimx - 1]
        frac = gm.region_prop[r].item()
        Hdm = gm.dist_max[r].entropy().v
        Hsd = gm.dist[r].entropy().v
        if INCLUDE:
            eps_iv = IC.Ival(arb(0), arb(0))
            dm = gm.dist_max[r].v
            for i, shp in enumerate(gm.shapes):
                if float(dm[i].hi) <= 0: continue
                g_val = gm.lam_sum[r].v[0].lo - 1
                for d2 in range(3):
                    g_val = g_val + gm.lam_margin[(r, d2)].v[shp[d2]].lo
                diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
                e = diff.pos()
                eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
            p_true = frac * (Hdm - Hsd + 2 * eps_iv)
        else:
            p_true = frac * (Hdm - Hsd)
        R_low = nbv.lo - p_true.hi
        if float(R_low) > 0:
            total_R = total_R + R_low
    M_low = min(float(z.lo) for z in gm.mat_size)
    tgt = (arb(7).log() * 4)
    up = (tgt - total_R) / M_low
    print(json.dumps({"include_lemma1": INCLUDE,
                      "R_sum_low": float(total_R.lower()),
                      "M_low": M_low,
                      "omega_up": float(up)}))
