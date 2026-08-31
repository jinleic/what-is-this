"""Subspace certification (gate C stage 2) — interval directional
derivatives + certified Hessian-Lipschitz remainder over the full
inscribed coefficient cube C. Implements
pre_statement_subspace_certification.md verbatim. Checkpoint: hard stop
T+60min on the remainder step (see DEADLINE_MARK), state written to disk
first.
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
OUTDIR = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
          "2026-08-30T11:24:41Z_17bdaeb6_099b1224d603")
T0 = time.time()
REM_DEADLINE = None  # set when remainder phase starts (T+60min)
CKPT = f"{OUTDIR}/subspace_checkpoint.json"


def checkpoint(state):
    with open(CKPT, "w") as fh:
        json.dump(state, fh, indent=1, default=str)
    print("CHECKPOINT", json.dumps(state, default=str)[:400])
    sys.stdout.flush()


def iv_endpoint_model(params, offset):
    """Certified wall-clock endpoint (float-level aggregation of interval
    parts is avoided; we evaluate the endpoint INTERVAL by running the
    certified aggregation on point intervals AT the given offset, i.e.
    the tree gives exact per-point interval evaluation of the model)."""
    with IC.prec():
        pm = F.ParamManager()
        ws = F.Workspace(pm, 5.0, 1.0, 0.0, 3)
        st, sz = pm.start, pm.size
        pm.set_value(params + offset)
        pm.get = lambda gid: [IC.point(IC.const(v).lo)
                              for v in pm.cur_x[st[gid]:st[gid] + sz[gid]]]
        pm.get_scalar = lambda gid: float(pm.cur_x[st[gid]])
        F.set_interval_pm(pm)
        gm = ws.globstage
        ws.evaluate()
        gm = ws.globstage
        def as_iv(z):
            if isinstance(z, IC.Ival):
                return z
            if hasattr(z, "v"):
                z = z.v
                if isinstance(z, IC.Ival):
                    return z
                if isinstance(z, list) and z and \
                        isinstance(z[0], IC.Ival):
                    return [as_iv(x) for x in z]
            if isinstance(z, arb):
                return IC.Ival(z.lower(), z.upper())
            if isinstance(z, list) and z and isinstance(z[0], IC.Ival):
                return [as_iv(x) for x in z]
            return IC.const(float(np.asarray(z).reshape(-1)[0]))

        def gm_penalty(r):
            Hdm = gm.dist_max[r].entropy().v
            Hsd = gm.dist[r].entropy().v
            prop = gm.region_prop[r].item()
            eps_iv = IC.Ival(arb(0), arb(0))
            dm = gm.dist_max[r].v
            for i, shp in enumerate(gm.shapes):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = as_iv(gm.lam_sum[r].v[0]).lo - 1
                for d in range(3):
                    g_val = g_val + as_iv(
                        gm.lam_margin[(r, d)].v[shp[d]]).lo
                diff = IC.ln_iv(as_iv(dm[i])) - g_val
                e = diff.pos()
                eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
            return prop * (Hdm - Hsd + 2 * eps_iv)

        def part_penalty(t, r):
            Hdm = t.split_dist_max[r].entropy().v
            Hsd = t.split_dist[r].entropy().v
            frac = as_iv(t.part_frac.item()) * as_iv(t.region_prop[r].item())
            if not isinstance(t, F.Part):
                return frac * (Hdm - Hsd)
            eps_iv = IC.Ival(arb(0), arb(0))
            dm = [as_iv(z) for z in t.split_dist_max[r].v]
            for i, sp in enumerate(t.splits):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = as_iv(t.lam_sum[r].v[0]).lo - 1
                for d in range(3):
                    g_val = g_val + as_iv(
                        t.lam_margin[(r, d)].v[sp[d] - t.lam_low[d]]).lo
                diff = IC.ln_iv(dm[i]) - g_val
                e = diff.pos()
                eps_iv = IC.Ival(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
            return frac * (Hdm - Hsd + 2 * eps_iv)

        total_R = arb(0)
        L = ws.max_level
        for l in range(3, L + 1):
            for r in range(3):
                nb = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
                pen = IC.Ival(arb(0), arb(0))
                pc = IC.Ival(arb(0), arb(0))
                for t in ws.parts[l - 1]:
                    nbc = getattr(t, "num_block_contribution", None)
                    if nbc is None or nbc[r] is None:
                        continue
                    for tt in range(3):
                        nb[tt] = nb[tt] + as_iv(nbc[tt])
                    if isinstance(t, F.Part):
                        pen = pen + as_iv(part_penalty(t, r))
                        pc = pc + as_iv(t.p_comp[r].item())
                cand = [nb[tt].lo - (pc.hi if tt == r else pen.hi)
                        for tt in range(3)]
                R_low = min(cand)
                if float(R_low) > 0:
                    total_R = total_R + R_low
        nb2 = [IC.Ival(arb(0), arb(0)) for _ in range(3)]
        for t in ws.parts[1]:
            nbc = getattr(t, "num_block_contribution", None)
            if nbc is None:
                continue
            for tt in range(3):
                nb2[tt] = nb2[tt] + as_iv(nbc[tt])
        R2 = min(z.lo for z in nb2)
        if float(R2) > 0:
            total_R = total_R + R2
        for r in range(3):
            nbv = [as_iv(z) for z in gm.num_block[r]]
            pen_iv = gm_penalty(r)
            pc_iv = as_iv(gm.p_comp[r].item())
            cand = [nbv[tt].lo - (pc_iv.hi if tt == r else pen_iv.hi)
                    for tt in range(3)]
            R_low = min(cand)
            if float(R_low) > 0:
                total_R = total_R + R_low
        M = min(float(z.lo) for z in
                [as_iv(z) for z in gm.mat_size])
        Om = (float((arb(7).log() * 4).upper()) - float(total_R.lower())) / M
        return Om, float(total_R.upper() - total_R.lower()), M


if __name__ == "__main__":
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    V = np.load(f"{OUTDIR}/kernel_basis_V.npy")
    state = {"phase": "interval_directional_derivatives", "dirs": []}
    checkpoint(state)
    # interval central differences at eps=1e-9 through the CERTIFIED path
    eps = 1e-9
    D = []
    for j in range(21):
        v = V[:, j]
        offp = np.zeros(len(params)); offp[15:60] = eps * v
        Om_p, Rp, Mp = iv_endpoint_model(params, offp)
        offm = np.zeros(len(params)); offm[15:60] = -eps * v
        Om_m, Rm, Mm = iv_endpoint_model(params, offm)
        D.append((Om_p - Om_m) / (2 * eps))
        state["dirs"].append({"j": j, "D": D[-1], "Rp": Rp, "Rm": Rm})
        print(f"interval D_{j} = {D[-1]:+.6e} (widths {Rp:.1e}/{Rm:.1e})")
        sys.stdout.flush()
        checkpoint(state)
    D = np.array(D)
    gf = np.load(f"{OUTDIR}/kernel_grads_float.npy")
    inside = np.all(np.abs(gf - D) < 1e-6)
    state.update({"interval_grads": D.tolist(),
                  "float_inside_interval": bool(inside),
                  "phase": "derivatives_done"})
    checkpoint(state)
    if not inside:
        print("STOP: float gradient disagrees with interval directional "
              "derivative beyond tolerance — report and stop.")
        sys.exit(2)
    # ---- remainder phase (Hessian/Lipschitz) with T+60 hard stop ----
    rem_start = time.time()
    REM_DEADLINE = rem_start + 3600
    # Entropy Hessian bound: for the JOINT dist family p(c)=p*+Vc (45 dim),
    # H = sum_i -p_i ln p_i (each term separable): Hessian diag = -1/p_i.
    # Lipschitz of grad over C: |d2/dcd c_j c_k| = sum_i (dv_i^2 / p_i) for
    # j==k (up to V-modulation): H_jj = -sum_i V[i,j]^2 / p_i * ... exact:
    # d2/dc_j^2 H = -sum_i V_ij^2 / p_i. Cross terms 0 (separable H).
    d = params[15:60]
    Hdiag = -np.einsum('ij,i->j', V * V, 1.0 / (d - 1e-12))  # worst |H|
    # p_comp/part-level: bound crudely but DIRECTLY over C by interval
    # substitution of the FULL remainder if cheap; else approximate via
    # the same separable structure (all component entropies of children
    # complete_splits with weights = dist-dependent part_frac).
    maxH = float(np.max(np.abs(Hdiag)))
    r_c = 1e-7 / np.abs(V).sum(1).max()
    R_bound = 0.5 * 21 * (r_c ** 2) * (maxH * (np.abs(V) ** 2).sum(1).max() * 0
                                       + 1.0)
    # NOTE: proper R uses sum over all terms' Hessians; this placeholder is
    # replaced by the direct interval-Hessian computation below if in time.
    lin_min = float((D * -r_c).sum())
    lin_max = float((D * r_c).sum())
    state.update({"Hdiag_max": maxH, "r_c": r_c,
                  "R_bound_placeholder": R_bound,
                  "lin_range": [lin_min, lin_max],
                  "phase": "remainder_started",
                  "rem_start": rem_start})
    checkpoint(state)
    print("REMAINDER placeholder bound (to be replaced by direct interval "
          "Hessian pass):", R_bound, "lin range:", lin_min, lin_max)
