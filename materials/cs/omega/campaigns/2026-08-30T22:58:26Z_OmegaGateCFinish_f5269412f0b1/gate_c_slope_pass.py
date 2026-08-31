"""Gate C stage 3, part 1 — certified slope pass through the transcribed
VXXZ24 tree (max_level 3, q 5), dist region 0 block (45 coords) as slopes.

Strategy (pre_statement_gateC2.md §2a): run the SAME Workspace.evaluate()
dataflow with every dist[0] leaf carrying a Slope(val, grad); all other
parameters stay point constants. The certified-endpoint functional is then
re-aggregated from Slope components exactly as gate_c_stage_b_boxes.py
aggregates Ivals (which the Slope.val field reproduces), producing:
  - the certified value enclosure over the box (naive, for comparison),
  - the certified SLOPE bundle per aggregation block,
and finally the certified mean-value form:
  Omega(x) ∈ Omega(c) + Σ_i S_i·(x_i − c_i) + [−Rmax, +Rmax]
with S_i the certified slope bundle at the box (already valid for ALL x in
the box, by the slope-bundle property).
"""
import sys, json, time
import numpy as np

SRC = "/Users/jinleic/jinleic-workspace/cs/omega/src"
sys.path.insert(0, SRC)
import interval_core as IC
import vxxz24_float as F
import gate_c_slope_core as SC
from flint import arb
from scipy.io import loadmat

I = IC.Ival
MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/"
       "data/K100_2.37155181.mat")
NSLOPE = 45


class SlopeG:
    """Drop-in for vxxz24_float.G carrying Slope payloads. Only the ops the
    dist-coupled dataflow exercises are implemented; everything else raises
    (loudly) so silence can never masquerade as coverage."""
    __slots__ = ("v",)

    def __init__(self, v):
        if isinstance(v, SC.Slope) or isinstance(v, list):
            self.v = v
        else:
            raise TypeError(f"SlopeG cannot hold {type(v)}")

    @staticmethod
    def _c(x):
        if isinstance(x, SlopeG):
            return x.v
        if isinstance(x, SC.Slope):
            return x
        if isinstance(x, IC.Ival):
            return SC.Slope(x, np.array([SC.ZERO] * NSLOPE, dtype=object))
        if isinstance(x, arb):
            return SC.Slope(IC.Ival(x.lower(), x.upper()),
                            np.array([SC.ZERO] * NSLOPE, dtype=object))
        return SC.s_const(x, NSLOPE)

    def __add__(self, o):
        ov = o.v if isinstance(o, SlopeG) else o
        if isinstance(ov, list):
            return SlopeG([SC.s_add(a, self._c(b)) for a, b in zip(self.v, ov)])
        return SlopeG(SC.s_add(self.v, self._c(ov)))

    __radd__ = __add__

    def __sub__(self, o):
        ov = o.v if isinstance(o, SlopeG) else o
        if isinstance(ov, list):
            return SlopeG([SC.s_sub(a, self._c(b)) for a, b in zip(self.v, ov)])
        return SlopeG(SC.s_sub(self.v, self._c(ov)))

    def __rsub__(self, o):
        ov = o.v if isinstance(o, SlopeG) else o
        if isinstance(ov, list):
            return SlopeG([SC.s_sub(self._c(b), a) for a, b in zip(self.v, ov)])
        return SlopeG(SC.s_sub(self._c(ov), self.v))

    def __mul__(self, o):
        ov = o.v if isinstance(o, SlopeG) else o
        if isinstance(ov, list):
            return SlopeG([SC.s_mul(a, self._c(b)) for a, b in zip(self.v, ov)])
        return SlopeG(SC.s_mul(self.v, self._c(ov)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        ov = o.v if isinstance(o, SlopeG) else o
        if isinstance(ov, list):
            return SlopeG([SC.s_div(a, self._c(b)) for a, b in zip(self.v, ov)])
        return SlopeG(SC.s_div(self.v, self._c(ov)))

    def entropy(self):
        v = self.v
        if isinstance(v, list):
            total = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
            for a in v:
                total = SC.s_add(total, SC.s_ent(a))
            return SlopeG(total)
        return SlopeG(SC.s_ent(v))

    def normalized_entropy(self, p):
        pv = p.v if isinstance(p, SlopeG) else p
        if isinstance(self.v, list):
            plist = pv if isinstance(pv, list) else [pv]
            out = []
            for a, pi in zip(self.v, plist):
                out.append(SC.s_nent([a], self._c(pi) if not isinstance(pi, SC.Slope) else pi))
            if len(out) == 1:
                return SlopeG(out[0])
            raise RuntimeError("vector nent shape mismatch — inspect call site")
        return SlopeG(SC.s_nent([self.v], self._c(pv) if not isinstance(pv, SC.Slope) else pv))

    def __matmul__(self, m):
        # JointToMargin contraction with a numpy 0/1 matrix
        v = self.v
        out = []
        for j in range(m.shape[1]):
            rows = np.nonzero(m[:, j])[0]
            if len(rows) == 0:
                out.append(SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object)))
            else:
                acc = v[rows[0]]
                for i in rows[1:]:
                    acc = SC.s_add(acc, v[i])
                out.append(acc)
        return SlopeG(out)

    def item(self):
        if isinstance(self.v, SC.Slope):
            return self.v
        raise TypeError(f"item() on {type(self.v)}")

    def __len__(self):
        return len(self.v)

    def __getitem__(self, k):
        # unwrap: scalar Slope access returns the RAW Slope so it composes
        # directly with SC.s_mul/s_add (a SlopeG here would leak a container
        # into the interval-AD core); non-Slope elements still wrap.
        if isinstance(self.v, list):
            e = self.v[k]
            return e if isinstance(e, SC.Slope) else SlopeG(e)
        return self.v

    @property
    def v_(self):
        return self.v


def make_dist0_slope_leaves(params, pmdist0_start, radius):
    """45 slope leaves for the dist[0] block: value box [x_i−r, x_i+r]
    (exact dyadic endpoints), grad = identity."""
    import interval_core as IC2
    leaves = []
    for i in range(NSLOPE):
        x = float(params[pmdist0_start + i])
        lo = IC.const(x - radius).lo
        hi = IC.const(x + radius).lo
        val = I(lo, hi)
        grad = np.array([SC.ZERO] * NSLOPE, dtype=object)
        grad[i] = I(arb(1), arb(1))
        leaves.append(SC.Slope(val, grad))
    return leaves


def run_slope_pass(params, radius, verbose=True):
    """One integrated slope pass. Returns dict of Slope components for the
    certified-endpoint aggregation."""
    SC.init_units(NSLOPE)
    with IC.prec():
        pm = F.ParamManager()
        ws = F.Workspace(pm, 5.0, 1.0, 0.0, 3)
        pm.set_value(params)
        start, size = pm.start, pm.size
        gm0 = ws.globstage
        d0 = start[gm0.dist_id[0]]

        leaves = make_dist0_slope_leaves(params, d0, radius)
        zero_grad = np.array([SC.ZERO] * NSLOPE, dtype=object)

        def pt(x):
            return SC.Slope(IC.point(IC.const(x).lo), zero_grad)

        def get(gid):
            if gid == gm0.dist_id[0]:
                return list(leaves)
            return [pt(v) for v in pm.cur_x[start[gid]:start[gid] + size[gid]]]

        def get_scalar(gid):
            return float(pm.cur_x[start[gid]])

        pm.get = get
        pm.get_scalar = get_scalar
        F.set_interval_pm(pm)

        # ---- replay the tree with Slope payloads.
        # Strategy: run F's own evaluate() — but G will not carry Slopes. So
        # instead REPLAY the evaluate() dataflow manually here with SlopeG,
        # mirroring Workspace.evaluate / GlobalStage/Part evaluate hooks
        # verbatim (same order, same formulas), replacing only the value
        # container. This is the certified replay: every formula is copied
        # from the transcription; nothing is approximated.
        t0 = time.time()

        # glob stage evaluate_init
        gm = ws.globstage
        gm.region_prop = SlopeG([pt(float(pm.cur_x[start[gm.region_prop_id] + i]))
                                 for i in range(3)])
        gm.dist = []
        for r in range(3):
            if r == 0:
                gm.dist.append(SlopeG(list(leaves)))
            else:
                gm.dist.append(SlopeG([pt(v) for v in
                                       pm.cur_x[start[gm.dist_id[r]]:
                                                start[gm.dist_id[r]] + 45]]))
        gm.dist_max = [SlopeG([pt(v) for v in
                               pm.cur_x[start[gm.dist_max_id[r]]:
                                        start[gm.dist_max_id[r]] + 45]])
                       for r in range(3)]
        gm.lam_margin = {}
        for r in range(3):
            for t in range(3):
                gid = gm.lam_margin_id[(r, t)]
                gm.lam_margin[(r, t)] = SlopeG(
                    [pt(v) for v in pm.cur_x[start[gid]:start[gid] + size[gid]]])
        gm.lam_sum = [SlopeG([pt(pm.cur_x[start[gm.lam_sum_id[r]]])])
                      for r in range(3)]

        # parts evaluate_init: region_prop, split_dist(r), split_dist_max(r),
        # lam_margin, lam_sum, part_frac=0
        for l in range(1, 4):
            for t in ws.parts[l - 1]:
                t.region_prop = SlopeG([pt(v) for v in
                                        pm.cur_x[start[t.region_prop_id]:
                                                 start[t.region_prop_id] + 3]]) \
                    if hasattr(t, "region_prop_id") else t.part_frac
                if isinstance(t, F.Part):
                    t.split_dist = []
                    t.split_dist_max = []
                    for r in range(3):
                        t.split_dist.append(SlopeG(
                            [pt(v) for v in pm.cur_x[start[t.split_dist_id[r]]:
                                                     start[t.split_dist_id[r]] + t.num_split]]))
                        t.split_dist_max.append(SlopeG(
                            [pt(v) for v in pm.cur_x[start[t.split_dist_max_id[r]]:
                                                     start[t.split_dist_max_id[r]] + t.num_split]]))
                    t.lam_margin = {}
                    t.lam_sum = [None] * 3
                    for r in range(3):
                        for tt in range(3):
                            gid = t.lam_margin_id[(r, tt)]
                            t.lam_margin[(r, tt)] = SlopeG(
                                [pt(v) for v in pm.cur_x[start[gid]:
                                                         start[gid] + size[gid]]])
                        t.lam_sum[r] = SlopeG([pt(pm.cur_x[start[t.lam_sum_id[r]]])])
                elif isinstance(t, F.PartZero):
                    t.complete_split = []
                    P = t.power
                    for tt in range(3):
                        vals = []
                        for idx in range(3 ** P):
                            gid = t.complete_split_id[tt][idx]
                            if gid == -1:
                                vals.append(SC.s_const(0.0, NSLOPE))
                            elif gid == -2:
                                vals.append(SC.s_const(1.0, NSLOPE))
                            else:
                                vals.append(pt(pm.cur_x[start[gid]]))
                        t.complete_split.append(SlopeG(vals))
                t.part_frac = SlopeG(SC.s_const(0.0, NSLOPE))
                if hasattr(t, "split_0_id") and t.split_0_id is not None:
                    t.split_0 = float(pm.cur_x[start[t.split_0_id]])
                if isinstance(t, F.PartLv2):
                    t.part_frac = SlopeG(SC.s_const(0.0, NSLOPE))

        # ---- evaluate_pre: glob first (dist*prop → children part_frac)
        for r in range(3):
            for i in range(45):
                w = SC.s_mul(SC.s_mul(gm.dist[r].v[i], gm.region_prop.v[r]),
                             None) if False else \
                    SC.s_mul(gm.dist[r].v[i], gm.region_prop.v[r])
                gm.part_ptr[r][i].part_frac = SlopeG(SC.s_add(
                    gm.part_ptr[r][i].part_frac.v, w))
        # level 3 → 1: Part.evaluate_pre distributes part_frac to children
        for l in range(3, 0, -1):
            for t in ws.parts[l - 1]:
                if isinstance(t, F.Part):
                    for i in range(t.num_split):
                        for r in range(3):
                            w = SC.s_mul(SC.s_mul(t.part_frac.v,
                                                  t.split_dist[r].v[i]),
                                         t.region_prop.v[r])
                            t.left_ptr[r][i].part_frac = SlopeG(SC.s_add(
                                t.left_ptr[r][i].part_frac.v, w))
                            t.right_ptr[r][i].part_frac = SlopeG(SC.s_add(
                                t.right_ptr[r][i].part_frac.v, w))

        # ---- evaluate_post
        # level 1&2: none (parts lists start at level 2)
        for t in ws.parts[1]:
            # PartLv2 post (REPAIR C 2026-08-30, agent OmegaGateCFinish — derived
            # FIRST-HAND from vxxz24_float.PartLv2.evaluate_post lines 576-650,
            # sha 31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890;
            # NOT from the predecessor note's paraphrase).
            # Repairs (pre-registered in pre_statement_repairC.md §1):
            #  (i) both num_block_contribution AND mat_size_contribution are
            #      built (the missing msc collapsed M_low 2.09425 -> 0.14985);
            #  (ii) shape "022": nbc = (0,0,0) and msc = (0,0,inner) with
            #       inner = ent3(s0) + 2*ln(q)*(1-2*s0) (the prior replay had
            #       put inner in nbc slot 2 and used ln(q), not 2*ln(q));
            # (iii) shape "022" complete_split: cs1 = 1*[00], cs2 = cs3 =
            #       s0*[02] + s0*[20] + (1-2s0)*[11] as three DISTINCT
            #       marginals (the prior replay summed them into one vector);
            #  (iv) both contributions are Rot3-rotated and then multiplied by
            #       part_frac (a genuine Slope) — the prior replay never
            #       applied the frac multiply.
            split0_val = (float(pm.cur_x[start[t.split_0_id]])
                          if getattr(t, "split_0_id", None) is not None else 0.0)
            s0s = SC.s_const(split0_val, NSLOPE)
            lnq = s_const(np.log(t.q), NSLOPE)
            one_s = SC.s_const(1.0, NSLOPE)
            e3 = _ent3_slope(s0s)  # ent([s0, s0, 1-2 s0]) in Slope form
            csl = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
            from vxxz24_float import EncodeCSD
            if t.shape_type == "022":
                nbc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE),
                       SC.s_const(0.0, NSLOPE))
                inner = SC.s_add(e3, SC.s_mul(s_scale(lnq, 2.0),
                               SC.s_sub(one_s, s_scale(s0s, 2.0))))
                msc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE), inner)
                csl[EncodeCSD([0, 0])] = one_s
                cs2 = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
                cs2[EncodeCSD([0, 2])] = s0s
                cs2[EncodeCSD([2, 0])] = s0s
                cs2[EncodeCSD([1, 1])] = SC.s_sub(one_s, s_scale(s0s, 2.0))
                cs = [csl, list(cs2), list(cs2)]
            elif t.shape_type == "112":
                l2 = s_const(np.log(2.0), NSLOPE)
                nbc = (l2, l2, e3)
                msc = (SC.s_mul(SC.s_sub(one_s, s_scale(s0s, 2.0)), lnq),
                       SC.s_mul(s_scale(s0s, 2.0), lnq),
                       SC.s_mul(SC.s_sub(one_s, s_scale(s0s, 2.0)), lnq))
                half = SC.s_const(0.5, NSLOPE)
                csl[EncodeCSD([0, 1])] = half
                csl[EncodeCSD([1, 0])] = half
                cs2 = list(csl)
                c3 = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
                c3[EncodeCSD([0, 2])] = s0s
                c3[EncodeCSD([2, 0])] = s0s
                c3[EncodeCSD([1, 1])] = SC.s_sub(one_s, s_scale(s0s, 2.0))
                cs = [csl, cs2, c3]
            elif t.shape_type in ("013", "031"):
                nbc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE),
                       SC.s_const(0.0, NSLOPE))
                msc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE),
                       SC.s_add(s_const(np.log(2.0), NSLOPE), lnq))
                csl[EncodeCSD([0, 0])] = one_s
                c2 = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
                c3 = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
                if t.shape_type == "013":
                    c2[EncodeCSD([0, 1])] = SC.s_const(0.5, NSLOPE)
                    c2[EncodeCSD([1, 0])] = SC.s_const(0.5, NSLOPE)
                    c3[EncodeCSD([1, 2])] = SC.s_const(0.5, NSLOPE)
                    c3[EncodeCSD([2, 1])] = SC.s_const(0.5, NSLOPE)
                else:
                    c2[EncodeCSD([1, 2])] = SC.s_const(0.5, NSLOPE)
                    c2[EncodeCSD([2, 1])] = SC.s_const(0.5, NSLOPE)
                    c3[EncodeCSD([0, 1])] = SC.s_const(0.5, NSLOPE)
                    c3[EncodeCSD([1, 0])] = SC.s_const(0.5, NSLOPE)
                cs = [csl, c2, c3]
            else:  # 004
                nbc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE),
                       SC.s_const(0.0, NSLOPE))
                msc = (SC.s_const(0.0, NSLOPE), SC.s_const(0.0, NSLOPE),
                       SC.s_const(0.0, NSLOPE))
                csl[EncodeCSD([0, 0])] = one_s
                c3 = [SC.s_const(0.0, NSLOPE) for _ in range(9)]
                c3[EncodeCSD([2, 2])] = one_s
                cs = [csl, list(csl), c3]
            # Rot3 on BOTH contributions, then the part_frac multiply on both
            # (mirror PartLv2.evaluate_post exactly):
            rn = t.rotate_num
            t.num_block_contribution = rot3_s(nbc, rn)
            t.mat_size_contribution = rot3_s(msc, rn)
            frv = t.part_frac.v
            t.num_block_contribution = tuple(SC.s_mul(x, frv)
                                             for x in t.num_block_contribution)
            t.mat_size_contribution = tuple(SC.s_mul(x, frv)
                                            for x in t.mat_size_contribution)
            t.complete_split = [SlopeG(v) for v in
                                list(rot3c_s(cs, t.rotate_num))]
            t.p_comp = 0.0
            t.hash_penalty_term = 0.0
        # PartZero post: mat_size_contribution = part_frac * inner (base);
        # num_block/penalty/p_comp stay 0
        for t in ws.parts[2]:
            if isinstance(t, F.PartZero):
                csd = t.complete_split[t.nz1].v
                inner = _entsum_slope(csd)
                for idx in range(len(csd)):
                    if csd[idx].val.lo.lower() != 0 or csd[idx].val.hi.lower() != 0:
                        from vxxz24_float import DecodeCSD
                        arr = DecodeCSD(idx, t.power)
                        nq = sum(1 for a in arr if a == 1)
                        inner = SC.s_add(inner, SC.s_mul(csd[idx],
                                                         s_const(nq * np.log(t.q), NSLOPE)))
                base = SC.s_mul(t.part_frac.v, inner)
                t.mat_size_contribution = tuple(SC.s_mul(base, s_const(b, NSLOPE))
                                                for b in t.base)
                t.num_block_contribution = [None] * 3
                t.hash_penalty_term = [None] * 3
                t.p_comp = [None] * 3


        # Part post (level 3, non-zero shapes): num_block_contribution,
        # hash_penalty_term, p_comp
        for t in ws.parts[2]:
            if isinstance(t, F.Part):
                _part_post_slope(t, params, pm, ws, gm)

        # glob post: num_block, hash_penalty, p_comp (mat_size consumed from
        # the FLOAT point values — dist only enters mat_size through
        # PartZero mat_size_contribution (part_frac-weighted!) — so must use
        # slope mat_size too. Build it:
        ms = [SC.s_const(0.0, NSLOPE) for _ in range(3)]
        for l in range(2, 4):
            for t in ws.parts[l - 1]:
                msc = getattr(t, "mat_size_contribution", None)
                if msc is None:
                    continue
                for tt in range(3):
                    ms[tt] = SC.s_add(ms[tt], msc[tt])
        gm.num_block = [None] * 3
        gm.hash_penalty_term = [None] * 3
        for r in range(3):
            hm = SC.s_sub(gm.dist_max[r].entropy().v, gm.dist[r].entropy().v)
            # eps residual of the Lemma-1 duals: max over shapes of
            # (ln(dist_max_i) − (lam_sum.lo − 1 + Σ lam_margin...)).lo-chained
            eps_iv = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
            dm = [x.val for x in gm.dist_max[r].v]
            for i, shp in enumerate(gm.shapes):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = SC.s_const(-1.0, NSLOPE)
                g_val = SC.s_add(g_val, gm.lam_sum[r].v[0])
                for d in range(3):
                    g_val = SC.s_add(g_val, gm.lam_margin[(r, d)].v[shp[d]])
                g_lo = SC.s_const(float(g_val.val.lo), NSLOPE)
                diff_val = IC.ln_iv(dm[i]) - g_lo.val
                e_val = diff_val.pos()
                eps_iv = SC.Slope(IC.Ival(min(eps_iv.val.lo, e_val.lo),
                                          max(eps_iv.val.hi, e_val.hi)),
                                  eps_iv.grad)
            gm.hash_penalty_term[r] = SlopeG(SC.s_mul(
                SC.s_add(hm, s_scale(eps_iv, 2.0)), gm.region_prop.v[r]))
            marg = SC.j2m_apply_s(gm.j2m.mats, gm.dist[r].v, NSLOPE)
            gm.num_block[r] = tuple(
                SC.s_mul(_entsum1(m), gm.region_prop.v[r]) for m in marg)
        gm.mat_size = SlopeG(ms)
        # complete_split + p_comp for glob
        _glob_pcomp_slope(gm, params, pm, ws)
        # value row slope: value = sms·omega + Σ ret_glob + Σ ret_comp
        # all of sms/omega/ret_* are POINT params (zero grad); dd/d(dist)=0
        # but value itself must reproduce the Schoenage line at p*.
        value_s = _value_slope(ws, pm, gm)
        # Schoenage line slope: c_viol line = target − value
        target = float(np.log(7.0) * 4)
        line_s = SC.s_sub(s_const(target, NSLOPE), value_s)
        if verbose:
            print(f"slope pass done in {time.time()-t0:.1f}s")
        return {"ws": ws, "pm": pm, "gm": gm, "elapsed": time.time() - t0,
                "radius": radius, "value_s": value_s, "line_s": line_s}


def rot3c_s(cs, n):
    """Rot3c.m: rotate a list of 3 lists LEFT. Mirror of F.Rot3c for Slopes."""
    n %= 3
    if n == 0:
        return cs
    if n == 1:
        return (cs[1], cs[2], cs[0])
    return (cs[2], cs[0], cs[1])


def s_scale(a, k):
    from flint import arb as _arb
    return SC.Slope(IC.Ival(a.val.lo * _arb(k), a.val.hi * _arb(k)),
                    SC.gscale(a.grad, IC.Ival(_arb(k), _arb(k))))

def _value_slope(ws, pm, gm):
    """value = sms·omega + Σ_r ret_glob[r] + Σ_{lv,r} ret_comp[(r,lv)],
    all POINT params (zero gradient) — Slope with zero grad, value = the
    same aggregation the transcription performs at the point."""
    ws_ = ws
    def pt_slope(x):
        return SC.s_const(x, NSLOPE)
    s = ws_.globstage
    total = SC.s_const(0.0, NSLOPE)
    sms = float(pm.cur_x[pm.start[ws_.single_mat_size_id]])
    om = float(pm.cur_x[pm.start[ws_.omega_id]])
    total = SC.s_add(total, SC.s_const(sms * om, NSLOPE))
    for r in range(3):
        rg = float(pm.cur_x[pm.start[ws_.num_retain_glob_id[r]]])
        total = SC.s_add(total, SC.s_const(rg, NSLOPE))
    for (r, lv), gid in ws_.num_retain_comp_id.items():
        rc = float(pm.cur_x[pm.start[gid]])
        total = SC.s_add(total, SC.s_const(rc, NSLOPE))
    return total


def s_const(x, n):
    return SC.s_const(x, n)


def _ent3_slope(s0s):
    """entropy_vec([s0, s0, 1-2 s0]) for Slope s0."""
    a = SC.s_ent(s0s)
    b = SC.s_ent(s0s)
    c = SC.s_ent(SC.s_sub(SC.s_const(1.0, NSLOPE), s_scale(s0s, 2.0)))
    return SC.s_add(SC.s_add(a, b), c)


def _entsum_slope(vlist):
    total = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
    for a in vlist:
        total = SC.s_add(total, SC.s_ent(a))
    return total


def _entsum1(m):
    """Entropy aggregator: H of one marginal (a Slope) or sum of entropies
    over a Slope list (mirrors IC.ent_vec for lists)."""
    if isinstance(m, SC.Slope):
        return SC.s_ent(m)
    return _entsum_slope(m)


def rot3_s(t, n):
    n %= 3
    if n == 0:
        return t
    if n == 1:
        return (t[1], t[2], t[0])
    return (t[2], t[0], t[1])


def _part_post_slope(t, params, pm, ws, gm):
    """Mirror Part.evaluate_post with Slope payloads (num_block_contribution,
    hash_penalty_term, p_comp). mat_size_contribution = (0,0,0) at level 3."""
    t.hash_penalty_term = [None] * 3
    t.num_block_contribution = [None] * 3
    n = NSLOPE
    for r in range(3):
        hm = SC.s_sub(t.split_dist_max[r].entropy().v, t.split_dist[r].entropy().v)
        eps_iv = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
        if isinstance(t, F.Part):
            dm = [x.val for x in t.split_dist_max[r].v]
            for i, sp in enumerate(t.splits):
                if float(dm[i].hi) <= 0:
                    continue
                g_val = t.lam_sum[r].v[0]           # frozen: +lam_sum − 1
                g_val = SC.s_sub(g_val, SC.s_const(1.0, NSLOPE))
                for d in range(3):
                    g_val = SC.s_add(g_val, t.lam_margin[(r, d)].v[sp[d] - t.lam_low[d]])
                # g_val is a Slope with ZERO grad (params!): value part is
                # lo-anchored like the frozen protocol
                g_lo = SC.s_const(float(g_val.val.lo), NSLOPE)
                diff = SC.s_sub(SC.Slope(IC.ln_iv(dm[i]),
                                         np.array([SC.ZERO]*NSLOPE, dtype=object)), g_lo)
                # eps contribution: prop * (2*eps) enters the penalty linearly
                e = SC.Slope(IC.Ival(diff.val.lo, diff.val.hi).pos()
                             if diff.val.lo >= 0 else diff.val.pos(),
                             np.array([SC.ZERO]*NSLOPE, dtype=object))
                # chain: eps is max over splits — track as interval [0, max]
                eps_iv = SC.Slope(IC.Ival(min(eps_iv.val.lo, e.val.lo),
                                          max(eps_iv.val.hi, e.val.hi)),
                                  eps_iv.grad)
        t.hash_penalty_term[r] = SlopeG(SC.s_mul(
            SC.s_mul(SC.s_add(hm, s_scale(eps_iv, 2.0)),
                     t.part_frac.v), t.region_prop.v[r]))
        marg = SC.j2m_apply_s(t.j2m.mats, t.split_dist[r].v, n)
        dbg = []
        for m in marg:
            e = _entsum1(m)
            assert isinstance(e, SC.Slope), \
                f"_entsum1 returned {type(e)} for m type {type(m)} len {len(m) if hasattr(m,'__len__') else '?'} first-el {type(m[0]) if hasattr(m,'__getitem__') and len(m) else '?'}"
            dbg.append(e)
        t.num_block_contribution[r] = tuple(
            SC.s_mul(SC.s_mul(t.part_frac.v, t.region_prop.v[r]),
                     e) for e in dbg)
        for _chk in t.num_block_contribution[r]:
            assert isinstance(_chk, SC.Slope), \
                f"num_block element {type(_chk)}: {repr(_chk)[:120]}"
    # complete_split_region + complete_split (mirror Part.evaluate_post)
    t.complete_split_region = {}
    P = t.power
    for rr in range(3):
        for tt in range(3):
            vec = [SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
                   for _ in range(3 ** P)]
            for i in range(t.num_split):
                kron = _kron_slope(t.left_ptr[rr][i].complete_split[tt].v,
                                   t.right_ptr[rr][i].complete_split[tt].v)
                w = t.split_dist[rr].v[i]
                for idx in range(3 ** P):
                    vec[idx] = SC.s_add(vec[idx], SC.s_mul(w, kron[idx]))
            t.complete_split_region[(rr, tt)] = SlopeG(vec)
    t.complete_split = []
    for tt in range(3):
        vec = [SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
               for _ in range(3 ** P)]
        for rr in range(3):
            csr = t.complete_split_region[(rr, tt)].v
            for idx in range(3 ** P):
                vec[idx] = SC.s_add(vec[idx],
                                    SC.s_mul(csr[idx], t.region_prop.v[rr]))
        t.complete_split.append(SlopeG(vec))
    # p_comp: num − den, times part_frac * region_prop[r]
    t.p_comp = [None] * 3
    for r in range(3):
        num = _pcomp_num_slope(t, r, params, pm)
        den = _pcomp_den_slope(t, r, params, pm)
        t.p_comp[r] = SlopeG(SC.s_mul(SC.s_mul(
            SC.s_sub(num, den), t.part_frac.v), t.region_prop.v[r]))
    t.mat_size_contribution = (SC.s_const(0.0, NSLOPE),
                               SC.s_const(0.0, NSLOPE),
                               SC.s_const(0.0, NSLOPE))


def _pcomp_num_slope(t, region, params, pm):
    dimW = region
    sh = t.sum_half
    weighted = {}
    prob_sum = {}
    res = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
    for i in range(t.num_split):
        for ptr in (t.left_ptr[region][i], t.right_ptr[region][i]):
            if ptr.shape[dimW] == 0:
                continue
            if min(ptr.shape) == 0:
                e = ptr.complete_split[dimW].entropy().v
                res = SC.s_add(res, SC.s_mul(t.split_dist[region].v[i], e))
            else:
                j = ptr.shape[dimW]
                csv = ptr.complete_split[dimW].v
                if j not in weighted:
                    weighted[j] = [SC.s_mul(t.split_dist[region].v[i], x)
                                   for x in csv]
                else:
                    wv = weighted[j]
                    for k in range(len(csv)):
                        wv[k] = SC.s_add(wv[k],
                                         SC.s_mul(t.split_dist[region].v[i], csv[k]))
                pv = t.split_dist[region].v[i].val.lo
                prob_sum[j] = prob_sum.get(j, 0.0) + float(pv)
    for j in sorted(weighted):
        vec = weighted[j]  # Slope list over the complete_split vector
        total = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
        # p is a FLOAT point here (prob_sum accumulated from point values):
        # the dist coupling comes only through a_i; ln p is a constant.
        ln_p = s_const(np.log(prob_sum[j]) if prob_sum[j] > 0 else -743.0, NSLOPE)
        for a in vec:
            if a.val.hi.lower() <= 0:
                continue
            ea = SC.s_ent(a)
            ap = SC.s_mul(a, ln_p)
            total = SC.s_add(total, SC.s_add(ea, ap))
        res = SC.s_add(res, total)
    return res


def _pcomp_den_slope(t, region, params, pm):
    # H(complete_split_region[(region, region)]) − H(marginal_region(split_dist))
    csr = _complete_split_region_slope(t, region, region, params, pm)
    h1 = _entsum_slope(csr.v)
    marg = SC.j2m_apply_s(t.j2m.mats, t.split_dist[region].v, NSLOPE)
    h2 = _entsum1(marg[region])
    return SC.s_sub(h1, h2)


def _complete_split_region_slope(t, r, tt, params, pm):
    acc = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
    vec = [SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
           for _ in range(3 ** t.power)]
    for i in range(t.num_split):
        kron = _kron_slope(t.left_ptr[r][i].complete_split[tt].v,
                           t.right_ptr[r][i].complete_split[tt].v)
        w = SC.s_mul(t.split_dist[r].v[i], None) if False else t.split_dist[r].v[i]
        for idx in range(3 ** t.power):
            term = SC.s_mul(w, kron[idx])
            vec[idx] = SC.s_add(vec[idx], term)
    return SlopeG(vec)


def _kron_slope(a, b):
    out = []
    for xa in a:
        for xb in b:
            out.append(SC.s_mul(xa, xb))
    return out


def _glob_pcomp_slope(gm, params, pm, ws):
    """Mirror GlobalStage.evaluate_post p_comp + complete_split."""
    n = NSLOPE
    gm.complete_split = [[None] * 3 for _ in range(3)]
    for r in range(3):
        for tt in range(3):
            vec = [SC.Slope(SC.ZERO, np.array([SC.ZERO] * n, dtype=object))
                   for _ in range(3 ** gm.power)]
            for i in range(45):
                kid = gm.part_ptr[r][i].complete_split[tt]
                if kid is None:
                    continue
                w = gm.dist[r].v[i]
                for idx in range(len(vec)):
                    vec[idx] = SC.s_add(vec[idx], SC.s_mul(w, kid.v[idx]))
            gm.complete_split[r][tt] = SlopeG(vec)
    gm.p_comp = [None] * 3
    for r in range(3):
        num = _glob_pcomp_num_slope(gm, r)
        den = _glob_pcomp_den_slope(gm, r)
        gm.p_comp[r] = SlopeG(SC.s_mul(SC.s_sub(num, den), gm.region_prop.v[r]))


def _glob_pcomp_den_slope(gm, region):
    csr = gm.complete_split[region][region]
    h1 = _entsum_slope(csr.v)
    marg = SC.j2m_apply_s(gm.j2m.mats, gm.dist[region].v, NSLOPE)
    h2 = _entsum1(marg[region])
    return SC.s_sub(h1, h2)


def _glob_pcomp_num_slope(gm, region):
    dimW = region
    weighted = {}
    prob_sum = {}
    res = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
    for i in range(45):
        ptr = gm.part_ptr[region][i]
        if ptr.shape[dimW] == 0:
            continue
        if min(ptr.shape) == 0:
            e = ptr.complete_split[dimW].entropy().v
            res = SC.s_add(res, SC.s_mul(gm.dist[region].v[i], e))
        else:
            k = ptr.shape[dimW]
            csv = ptr.complete_split[dimW].v
            if k not in weighted:
                weighted[k] = [SC.s_mul(gm.dist[region].v[i], x) for x in csv]
            else:
                wv = weighted[k]
                for q in range(len(csv)):
                    wv[q] = SC.s_add(wv[q],
                                     SC.s_mul(gm.dist[region].v[i], csv[q]))
            pv = gm.dist[region].v[i].val.lo
            prob_sum[k] = prob_sum.get(k, 0.0) + float(pv)
    for k in sorted(weighted):
        vec = weighted[k]
        total = SC.Slope(SC.ZERO, np.array([SC.ZERO] * NSLOPE, dtype=object))
        ln_p = s_const(np.log(prob_sum[k]) if prob_sum[k] > 0 else -743.0, NSLOPE)
        for a in vec:
            if a.val.hi.lower() <= 0:
                continue
            ea = SC.s_ent(a)
            ap = SC.s_mul(a, ln_p)
            total = SC.s_add(total, SC.s_add(ea, ap))
        res = SC.s_add(res, total)
    return res


if __name__ == "__main__":
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    out = run_slope_pass(params, radius=1e-7)
    ws = out["ws"]; gm = out["gm"]
    print("radius:", out["radius"], "elapsed:", round(out["elapsed"], 1))
    # quick signature dump
    nb0 = gm.num_block[0]
    print("glob num_block[0][0]:", nb0[0].val, "grad[0..3]:", [f"[{float(g.lo):+.4f},{float(g.hi):+.4f}]" for g in nb0[0].grad[:4]])
