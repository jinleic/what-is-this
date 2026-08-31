"""V6 CONTROL — full-tree value parity of the repaired slope replay against
the frozen interval protocol (vxxz24_float + stage_b_boxes-certified_at
aggregation), at the centre (radius 0). Checks every aggregation BLOCK
(level-3 per r: nb tt candidates with pen/pc; level-2 nb2; glob per r;
mat_size per tt; value/Schoenage line). Tolerance 5e-13 per component.
"""
import sys
sys.path.insert(0, "src")
import numpy as np
import gate_c_slope_pass as SP
import gate_c_slope_core as SC
import interval_core as IC
import vxxz24_float as F
from scipy.io import loadmat
from flint import arb

TOL = 5e-13


def as_lo(z):
    if isinstance(z, SP.SlopeG):
        z = z.v
    if isinstance(z, SC.Slope):
        return float(z.val.lo)
    if isinstance(z, IC.Ival):
        return float(z.lo)
    if hasattr(z, "v"):
        z = z.v
    if isinstance(z, IC.Ival):
        return float(z.lo)
    if isinstance(z, (int, float)):
        return float(z)
    return float(np.asarray(z).reshape(-1)[0])


def as_hi(z):
    if isinstance(z, SP.SlopeG):
        z = z.v
    if isinstance(z, SC.Slope):
        return float(z.val.hi)
    if isinstance(z, IC.Ival):
        return float(z.hi)
    if hasattr(z, "v"):
        z = z.v
    if isinstance(z, IC.Ival):
        return float(z.hi)
    if isinstance(z, (int, float)):
        return float(z)
    return float(np.asarray(z).reshape(-1)[0])


params = np.asarray(loadmat(SP.MAT)["params"]).flatten()

# ---------------- slope replay (radius 0)
out0 = SP.run_slope_pass(params, radius=0.0, verbose=False)
ws = out0["ws"]; pm = out0["pm"]; gm = ws.globstage

# ---------------- frozen protocol tree
pmF = F.ParamManager(); wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
stF, szF = pmF.start, pmF.size
pmF.set_value(params)
pmF.get = lambda gid: [IC.point(IC.const(v).lo)
                       for v in pmF.cur_x[stF[gid]:stF[gid] + szF[gid]]]
pmF.get_scalar = lambda gid: float(pmF.cur_x[stF[gid]])
F.set_interval_pm(pmF)
gmF = wsF.globstage
wsF.evaluate()


def gm_penalty_true_frozen(g, r):
    Hdm = g.dist_max[r].entropy().v
    Hsd = g.dist[r].entropy().v
    prop = g.region_prop[r].item()
    eps_iv = IC.Ival(arb(0), arb(0))
    dm = [x.v if hasattr(x, "v") else x for x in g.dist_max[r].v]
    dm = [z if isinstance(z, IC.Ival) else IC.const(as_lo(z)) for z in dm]
    for i, shp in enumerate(g.shapes):
        if float(dm[i].hi) <= 0:
            continue
        g_val = as_lo(g.lam_sum[r].v[0]) - 1
        for d in range(3):
            g_val = g_val + as_lo(g.lam_margin[(r, d)].v[shp[d]])
        e = (IC.ln_iv(dm[i]) - IC.const(g_val)).pos()
        eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
    Hdm_l = float(dm_iv_lo(g.dist_max[r])) if False else None
    return prop, eps_iv


def dm_iv_lo(gvec):
    return None


# simpler: compute frozen penalty components as floats
worst = 0.0
worst_where = None


def upd(name, a, b):
    global worst, worst_where
    d = abs(a - b)
    if d > worst:
        worst = d
        worst_where = (name, a, b)


# ---- glob blocks
for r in range(3):
    nbS = [as_lo(gm.num_block[r][tt]) for tt in range(3)]
    nbF = [as_lo(gmF.num_block[r][tt]) for tt in range(3)]
    for tt in range(3):
        upd(f"glob nb r{r} tt{tt}", nbS[tt], nbF[tt])
    penS = as_hi(gm.hash_penalty_term[r])
    penF = as_hi(gmF.hash_penalty_term[r])
    upd(f"glob pen r{r}", penS, penF)
    pcS = as_hi(gm.p_comp[r])
    pcF = as_hi(gmF.p_comp[r])
    upd(f"glob pc r{r}", pcS, pcF)

# ---- level-3 blocks
for r in range(3):
    nbS = [0.0, 0.0, 0.0]
    penS = 0.0
    pcS = 0.0
    nbF = [0.0, 0.0, 0.0]
    penF = 0.0
    pcF = 0.0
    for t in ws.parts[2]:
        nbc = getattr(t, "num_block_contribution", None)
        if nbc is None or nbc[r] is None:
            continue
        for tt in range(3):
            nbS[tt] += as_lo(nbc[r][tt])
        if isinstance(t, F.Part):
            penS += as_hi(t.hash_penalty_term[r])
            pcS += as_hi(t.p_comp[r])
    for t in wsF.parts[2]:
        nbc = getattr(t, "num_block_contribution", None)
        if nbc is None or nbc[r] is None:
            continue
        for tt in range(3):
            nbF[tt] += as_lo(nbc[r][tt])
        if isinstance(t, F.Part):
            penF += as_hi(t.hash_penalty_term[r])
            pcF += as_hi(t.p_comp[r])
    for tt in range(3):
        upd(f"L3 nb r{r} tt{tt}", nbS[tt], nbF[tt])
    upd(f"L3 pen r{r}", penS, penF)
    upd(f"L3 pc r{r}", pcS, pcF)

# ---- level-2
nb2S = [0.0, 0.0, 0.0]
nb2F = [0.0, 0.0, 0.0]
for t in ws.parts[1]:
    nbc = getattr(t, "num_block_contribution", None)
    if nbc is None:
        continue
    for tt in range(3):
        nb2S[tt] += as_lo(nbc[tt])
for t in wsF.parts[1]:
    nbc = getattr(t, "num_block_contribution", None)
    if nbc is None:
        continue
    for tt in range(3):
        nb2F[tt] += as_lo(nbc[tt])
for tt in range(3):
    upd(f"L2 nb tt{tt}", nb2S[tt], nb2F[tt])

# ---- mat_size
msS = [as_lo(z) for z in gm.mat_size.v]
msF = [as_lo(z) for z in gmF.mat_size]
for tt in range(3):
    upd(f"ms tt{tt}", msS[tt], msF[tt])

# ---- value + Schoenage line
vs = as_lo(out0["value_s"])
vl = as_lo(out0["line_s"])
valF = as_lo(wsF.value)
tgt_hi = float((arb(7).log() * 4).upper())
upd("value", vs, valF)
upd("line", vl, tgt_hi - valF)

print(f"V6 full-tree parity: worst |delta| = {worst:.3e}")
print("  worst at:", worst_where)
ok = worst < TOL
print("V6 PASS" if ok else "V6 FAIL")
sys.exit(0 if ok else 3)
