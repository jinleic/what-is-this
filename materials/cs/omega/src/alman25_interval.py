"""Stage (b): slack-absorption interval enclosure at the exact released point.

Per pre_statement Addendum A3. Same evaluation tree as the float stage
(alman25_float.Workspace.evaluate), but pm.get returns arb ball vectors
(exact dyadics) and every entropy is an interval enclosure. Output:
 - certified-feasible retained counts R (lower endpoints),
 - certified M (single matrix size lower endpoint),
 - Lemma-1 max-entropy slack eps with the shipped Lagrange multipliers
   (per H^max certificate: eps = max_a |ln dm(a) - (lam-sum - 1)|),
 - certified omega interval at the released point, and the explicit
   feasibility/absorption statement vs omega <= 2.371339.
"""

import sys, json, hashlib, contextlib
import numpy as np
from scipy.io import loadmat
from flint import arb, ctx

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import alman25_float as F

MID = 300
OUTB = 128


@contextlib.contextmanager
def prec():
    with ctx.workprec(MID):
        yield


def dy(v64):
    """Exact dyadic arb of a float64."""
    f = float(v64)
    if f == 0.0:
        return arb(0)
    assert np.isfinite(f)
    m, e = np.frexp(f)
    mi = int(np.ldexp(m, 53))
    ei = int(e) - 53
    if ei >= 0:
        return arb(mi * (1 << ei))
    return arb(mi) / arb(1 << (-ei))


def low_ws(x):
    """lower endpoint augmented by outward buffer 2^-OUTB (rounded DOWN)."""
    if isinstance(x, np.ndarray):
        return float(x.reshape(-1)[0]) - float(arb(2) ** -OUTB)
    if not isinstance(x, arb):
        x = arb(float(x))
    return x.lower() - arb(2) ** -OUTB


def up_ws(x):
    if isinstance(x, np.ndarray):
        return float(x.reshape(-1)[0]) + float(arb(2) ** -OUTB)
    if not isinstance(x, arb):
        x = arb(float(x))
    return x.upper() + arb(2) ** -OUTB


def rung1_interval():
    mat = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/alman_code/"
           "data/W1.00_2.371339.mat")
    with prec():
        pm = F.ParamManager()
        ws = F.Workspace(pm, F.Q, 1.0, omega=0.0, max_level=F.MAX_LEVEL)
        params = np.asarray(loadmat(mat)["params"]).flatten()
        assert len(params) == pm.num_input
        dg = hashlib.sha256(params.tobytes()).hexdigest()
        pm.set_value(params)

        start, size = pm.start, pm.size

        def get(gid):
            s = start[gid]
            n = size[gid]
            return [dy(v) for v in pm.cur_x[s:s + n]]

        def get_scalar(gid):
            return float(pm.cur_x[start[gid]])

        pm.get = get
        pm.get_scalar = get_scalar

        L = F.MAX_LEVEL
        q = F.Q
        terms = ws.terms

        # -------- evaluate the whole tree under interval mode ------------
        c_viol, ceq_viol, value = None, None, None
        res = ws.evaluate()

        # ---- recompute the certified quantities from interval objects ----
        R_sum = arb(0)          # certified sum of retained counts (low eps)
        eps_total = arb(0)      # accumulated Lemma-1 slack
        M_choice = None

        # level >= 3 hashing (X-dim constraint: P_true absorbs H^max gap)
        for l in range(3, L + 1):
            for r in range(6):
                dimx, dimy, dimz = F.Dims(r + 1)
                # sum of num_block contributions and penalties over terms
                nbX = arb(0)
                penX = arb(0)
                for t in terms[l - 1]:
                    nbX = nbX + t.num_block_contribution[r][dimx - 1]
                    if isinstance(t, F.Term):
                        frac = t.term_frac.item() * t.region_prop[r].item()
                        # H(dm) and H(sd) interval
                        Hdm = t.split_dist_max[r].entropy().v
                        Hsd = t.split_dist[r].entropy().v
                        # Lemma-1 residual with shipped lambdas
                        dm = t.split_dist_max[r].v
                        eps = arb(0)
                        for i, sp in enumerate(t.splits):
                            if dm[i].upper() <= 0:
                                continue
                            g = t.lam_sum[r].v[0] - 1
                            for d in range(3):
                                g = g + t.lam_margin[(r, d)].v[
                                    sp[d] - t.lam_low[d]]
                            e = abs(dm[i].log() - g)
                            eps = arb(min(eps.lower(), e.lower()),
                                      max(eps.upper(), e.upper()))
                        p_true = frac * (Hdm - Hsd + 2 * eps)
                        eps_total = eps_total + 2 * frac * eps
                        penX = penX + p_true
                R_sum = R_sum + low_ws(nbX - penX)

        # level-2 symmetric hashing
        for d in range(3):
            nb2 = arb(0)
            for t in terms[1]:
                nb2 = nb2 + dy(t.num_block[d])
            R_sum = R_sum + low_ws(nb2)

        # global hashing
        glob = ws.globstage
        for r in range(6):
            dimx, dimy, dimz = F.Dims(r + 1)
            frac = glob.region_prop[r].item()
            # num_block for each dim = fraction * H(marginal); compute interval
            nbv = [arb(0), arb(0), arb(0)]
            jm = glob.j2m.mats
            for d in range(3):
                # marginal entropy of dist[r] restricted to the d-marginal
                # interval vec:
                H = glob.num_block[r][d]
                # H is frac * entropy(marg); region_prop is float-derived so
                # use exact arb of its float value
                nbv[d] = dy(frac) * H if False else H
            # penalty: region-prop-weighted H^max gap + Lemma-1 eps
            Hdm = glob.dist_max[r].entropy().v
            Hsd = glob.dist[r].entropy().v
            dm = glob.dist_max[r].v
            eps = arb(0)
            for i, shp in enumerate(glob.shapes):
                if dm[i].upper() <= 0:
                    continue
                g = glob.lam_sum[r].v[0] - 1
                for d in range(3):
                    g = g + glob.lam_margin[(r, d)].v[shp[d]]
                e = abs(dm[i].log() - g)
                eps = arb(min(eps.lower(), e.lower()),
                          max(eps.upper(), e.upper()))
            p_true = dy(frac) * (Hdm - Hsd + 2 * eps)
            eps_total = eps_total + 2 * dy(frac) * eps
            # code uses min over 3 dims (one retained var per region)
            lows = [low_ws(nbv[d] - (p_true if d == dimx - 1 else 0)
                           - (glob.p_compY[r].v if d == dimy - 1 else 0)
                           - (glob.p_compZ[r].v if d == dimz - 1 else 0))
                    for d in range(3)]
            # NOTE: p_compY/Z here are interval objects appended as .v
            R_sum = R_sum + min(lows)

        # single matrix size
        ms = [glob.mat_size[d] for d in range(3)]
        M_choice = min(low_ws(m) for m in ms)

        target = arb(q + 2).log() * (2 ** (L - 1))
        Om_pub = arb("2.371339")
        lhs = R_sum + M_choice * Om_pub
        Om_cert = (target - R_sum) / M_choice
        feasible_at_pub = bool(lhs.lower() >= target.upper())
        out = {
            "rung": "alman25_2.371339",
            "params_sha256": dg,
            "R_sum_lower": float(R_sum.lower()),
            "R_sum_upper": float(R_sum.upper()),
            "M_choice_lower": float(M_choice.lower()),
            "lhs_lower_pub": float(lhs.lower()),
            "target_upper": float(target.upper()),
            "feasible_at_published_omega_2_371339": feasible_at_pub,
            "omega_cert_lower": float(Om_cert.lower()),
            "omega_cert_upper": float(Om_cert.upper()),
            "eps_total_upper": float(eps_total.upper()),
        }
    return out


if __name__ == "__main__":
    res = rung1_interval()
    print(json.dumps(res, indent=1))
