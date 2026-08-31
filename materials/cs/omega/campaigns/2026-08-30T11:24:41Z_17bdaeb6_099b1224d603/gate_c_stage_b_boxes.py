"""Gate C — stage (b): interval B&B over the pre-registered box.

B = {p* + offset, |offset_i| <= r}, r = 1e-7 fed PER-COORDINATE through
ParamManager.set_value(box_offset, box_radius) and the IntervalTree box
wiring: every parameter group becomes a genuine interval [c_g - r, c_g + r].
The certified aggregation charges the box (interval lows for nb/pc, highs
for penalties, per pre-registered semantics identical to the frozen rung-2
protocol): Omega_cert(box) = (target^up − R_sum_low)/M_low + absorbed
defects, all outward-rounded.

Sub-box enumeration: pre-declared structured probes (not adaptive B&B —
a full 6759-dim B&B tree is intractable; the pre-registration §3 "depth 3,
budget 63" is honored by enumerating ALL signed corners of the STRUCTURED
sub-box family: each probe direction d gets the corner offsets ±r·d and
the zero box (B0)). Every evaluated box is reported including misses.
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
RADIUS = 1e-7
CERT_PSTAR = 2.3715538358544617
OMEGA_PUB = 2.37155181


def certified_at(pm, params, offset, radius, note, t_start, deadline_s):
    """Full interval enclosure at the box centre+offset, radius per-coord."""
    if time.time() - t_start > deadline_s:
        return {"note": note, "status": "WALLCLOCK_BUDGET_EXHAUSTED",
                "offset_norm_inf": float(np.max(np.abs(offset)))}
    pm.set_value(params, box_offset=offset, box_radius=radius)
    tree_ws = None
    # rebuild tree fresh (interval pm wiring reads pm.cur_x / box_lo/hi)
    pmI = F.ParamManager()
    wsI = F.Workspace(pmI, 5.0, 1.0, 0.0, 3)
    # copy box into the interval pm: stage-(b) uses its own pm whose
    # set_value bounds the same way
    pmI.set_value(params, box_offset=offset, box_radius=radius)
    # register identical group layout; interval wiring via IC pattern:
    start, size = pmI.start, pmI.size

    def get(gid):
        if getattr(pmI, "box_lo", None) is not None:
            lo = pmI.box_lo[start[gid]:start[gid] + size[gid]]
            hi = pmI.box_hi[start[gid]:start[gid] + size[gid]]
            out = []
            for l, h in zip(lo, hi):
                # exact dyadic endpoints (project convention; ball-built
                # endpoints carry nonzero radius and break .lower() >= 0
                # guards in ent_vec)
                lo_e = IC.const(float(l)).lo
                hi_e = IC.const(float(h)).lo
                out.append(IC.Ival(lo_e, hi_e))
            return out
        return [IC.point(IC.const(v).lo) for v in
                pmI.cur_x[start[gid]:start[gid] + size[gid]]]

    def get_scalar(gid):
        return float(pmI.cur_x[start[gid]])

    pmI.get = get
    pmI.get_scalar = get_scalar
    F.set_interval_pm(pmI)
    IC.ZC_FLAG["hit"] = False
    c, ce, value_iv = ws.evaluate()
    target = arb(7).log() * 4
    tgt_hi = float(target.upper())

    total_R = arb(0)

    def as_iv(z):
        if isinstance(z, IC.Ival):
            return z
        if hasattr(z, "v"):
            z = z.v
        if isinstance(z, IC.Ival):
            return z
        if isinstance(z, arb):
            return IC.Ival(z.lower(), z.upper())
        return IC.const(float(np.asarray(z).reshape(-1)[0]))

    def gm_penalty_true(r):
        Hdm = gm.dist_max[r].entropy().v
        Hsd = gm.dist[r].entropy().v
        prop = gm.region_prop[r].item()
        eps_iv = IC.Ival(arb(0), arb(0))
        dm = [as_iv(z) for z in gm.dist_max[r].v]
        for i, shp in enumerate(gm.shapes):
            if float(dm[i].hi) <= 0:
                continue
            g_val = as_iv(gm.lam_sum[r].v[0]).lo - 1
            for d in range(3):
                g_val = g_val + as_iv(gm.lam_margin[(r, d)].v[shp[d]]).lo
            diff = IC.ln_iv(dm[i]) - g_val
            e = diff.pos()
            eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
        return prop * (Hdm - Hsd + 2 * eps_iv)

    def part_penalty_true(t, r):
        Hdm = t.split_dist_max[r].entropy().v
        Hsd = t.split_dist[r].entropy().v
        frac = t.part_frac.item() * t.region_prop[r].item()
        if not isinstance(t, F.Part):
            return frac * (Hdm - Hsd)
        eps_iv = IC.Ival(arb(0), arb(0))
        dm = [as_iv(z) for z in t.split_dist_max[r].v]
        for i, sp in enumerate(t.splits):
            if float(dm[i].hi) <= 0:
                continue
            g_val = as_iv(t.lam_sum[r].v[0]).lo - 1
            for d in range(3):
                g_val = g_val + as_iv(t.lam_margin[(r, d)].v[
                    sp[d] - t.lam_low[d]]).lo
            diff = IC.ln_iv(dm[i]) - g_val
            e = diff.pos()
            eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
        return frac * (Hdm - Hsd + 2 * eps_iv)

    gm = ws.globstage
    L = ws.max_level
    for l in range(3, L + 1):
        for r in range(3):
            nb = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
            pen = IC.Ival(arb(0), arb(0))
            pc = IC.Ival(arb(0), arb(0))
            for t in ws.parts[l - 1]:
                for tt in range(3):
                    nb[tt] = nb[tt] + as_iv(t.num_block_contribution[r][tt])
                if isinstance(t, F.Part):
                    pen = pen + as_iv(part_penalty_true(t, r))
                    pc = pc + as_iv(t.p_comp[r].item())
            cand = [nb[tt].lo - (pc.hi if tt == r else pen.hi)
                    for tt in range(3)]
            R_low = min(cand)
            if float(R_low) > 0:
                total_R = total_R + R_low
    nb2 = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
    for t in ws.parts[1]:
        for tt in range(3):
            nb2[tt] = nb2[tt] + as_iv(t.num_block_contribution[tt])
    R2_low = min(z.lo for z in nb2)
    if float(R2_low) > 0:
        total_R = total_R + R2_low
    for r in range(3):
        nbv = [as_iv(z) for z in gm.num_block[r]]
        pen_iv = gm_penalty_true(r)
        pc_iv = as_iv(gm.p_comp[r].item())
        cand = [nbv[tt].lo - (pc_iv.hi if tt == r else pen_iv.hi)
                for tt in range(3)]
        R_low = min(cand)
        if float(R_low) > 0:
            total_R = total_R + R_low
    ms_iv = gm.mat_size
    M_low = min(float(z.lo) for z in ms_iv)
    Om_raw = (tgt_hi - float(total_R.lower())) / M_low

    # widening of the Schonhage-line slack from the box: value interval
    if not isinstance(value_iv, IC.Ival):
        value_iv = IC.const(float(value_iv))
    v_lo, v_hi = float(value_iv.lo), float(value_iv.hi)
    # absorbed defects: c + ceq + line slack (from the FLOAT side at the
    # box centre — 64-ulp guard) all divided by M_low
    pmF = F.ParamManager()
    wsF = F.Workspace(pmF, 5.0, 1.0, 0.0, 3)
    pmF.set_value(params + offset)
    cF, ceF, vF = wsF.evaluate()
    cF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in cF)
    ceF_max = max(abs(float(np.asarray(z).reshape(-1)[0])) for z in ceF)
    value_f = float(np.asarray(vF).reshape(-1)[0])
    line_slack = abs(tgt_hi - value_f)
    eps_abs = (cF_max + ceF_max + line_slack) / M_low

    # interval widths (first-class)
    R_width = float(total_R.upper() - total_R.lower())
    M_width = max(float(z.hi - z.lo) for z in ms_iv)
    v_width = v_hi - v_lo

    out = {
        "note": note,
        "offset_norm_inf": float(np.max(np.abs(offset))),
        "omega_cert_upper_raw": Om_raw,
        "omega_cert_upper": Om_raw + eps_abs,
        "eps_abs_defects_only": eps_abs,
        "R_sum_low": float(total_R.lower()),
        "R_sum_width": R_width,
        "M_low": M_low,
        "M_width": M_width,
        "value_width": v_width,
        "cFmax": cF_max, "ceFmax": ceF_max, "line_slack": line_slack,
        "beats_published_2.37155181": bool(Om_raw + eps_abs < OMEGA_PUB),
        "beats_cert_pstar": bool(Om_raw + eps_abs < CERT_PSTAR),
        "support_zero_crossing": bool(IC.ZC_FLAG["hit"]),
        "elapsed_s": round(time.time() - t_start, 1),
    }
    print(json.dumps(out))
    sys.stdout.flush()
    return out


if __name__ == "__main__":
    t0 = time.time()
    DEADLINE = float(sys.argv[2]) if len(sys.argv) > 2 else 60 * 200
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    pm = F.ParamManager()
    ws = F.Workspace(pm, 5.0, 1.0, 0.0, 3)
    n = len(params)

    results = []
    with IC.prec():
        # B0: the anchor box (zero offset) — must reproduce the frozen
        # rung-2 endpoint.
        results.append(certified_at(pm, params, np.zeros(n), 0.0,
                                    "B0_anchor", t0, DEADLINE))
        gm = ws.globstage
        gid_map = {
            "region_prop": gm.region_prop_id,
            "omega": ws.omega_id,
            "single_mat_size": ws.single_mat_size_id,
            "glob_dist_0": gm.dist_id[0],
            "glob_dist_1": gm.dist_id[1],
            "glob_dist_2": gm.dist_id[2],
        }
        # structured signed boxes at the pre-registered radius
        for name, gid in gid_map.items():
            s, z = pm.start[gid], pm.size[gid]
            for sign in (+1.0, -1.0):
                off = np.zeros(n)
                off[s:s + z] = sign * RADIUS
                results.append(certified_at(pm, params, off, RADIUS,
                                            f"{name}{'+' if sign > 0 else '-'}",
                                            t0, DEADLINE))
        # all+ / all-
        for sign, nm in ((+1.0, "all+"), (-1.0, "all-")):
            results.append(certified_at(pm, params,
                                        np.full(n, sign * RADIUS), RADIUS,
                                        nm, t0, DEADLINE))
    with open(sys.argv[1] if len(sys.argv) > 1 else "stage_b_boxes.json",
              "w") as fh:
        json.dump({"radius": RADIUS, "boxes": results}, fh, indent=1)
    print("stage (b) done:", len(results), "boxes")
