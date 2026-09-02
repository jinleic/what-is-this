"""Stage (b) core v2 — endpoint-pair interval representation (v11).

v11 provenance: 4 hard-hash dependency chain (alman25_float + interval_core
+ vxxz24_float + gate_c_candidate_witness_dependency); exact finite Arb
point extraction supports arbitrary exact finite dyadics (>53-bit);
binary64-only conversion is deprecated in favor of arb_exact_fraction.


ROOT CAUSE FIXED: python-flint's `arb(a, b)` constructor is (mid, radius),
NOT an interval [a, b]. The v1 stores used it as an interval constructor,
which is why the aggregation produced impossible values. This module uses an
explicit (lo, hi) endpoint-pair representation throughout; balls are only
created at the inside of transcendental-function evaluations and immediately
converted back to endpoint pairs via their outward endpoints.

Representation: Ival = (lo, hi) with both arb, lo <= hi (enforced).
All arithmetic on endpoint pairs uses exact outward interval rules.

Boundary & sanity rules from Main directive:
 R1 explicit type-wrapping at boundary; assert types at aggregation entries.
 R2 sanity asserts on aggregates (retained sum >= 0, sizes > 0, eps >= 0).
 R3 per-block containment check vs stage (a) float values, before margins.
 R4 outward endpoints for every consumed inequality: low()/up() add an
    explicit buffer of 2^-OUTB (in addition to exact endpoint arithmetic).
"""

import sys, os, hashlib, contextlib

if not __debug__:  # v11 (Main D3): asserts must never strip in
    raise RuntimeError("production run requires __debug__; "
                       "python -O would strip interval asserts")

import numpy as np
from scipy.io import loadmat
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # v11 campaign-local perimeter
import alman25_float as F

MID = 300
OUTB = 64


@contextlib.contextmanager
def prec():
    with ctx.workprec(MID):
        yield


# ------------------------------------------------------------- Ival type

class Ival:
    """Endpoint-pair interval over the reals (lo, hi), both arb."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        assert isinstance(lo, arb) and isinstance(hi, arb), \
            f"Ival endpoints must be arb, got {type(lo)}, {type(hi)}"
        assert lo <= hi, f"inverted Ival [{lo}, {hi}]"
        self.lo = lo
        self.hi = hi

    def __add__(self, o):
        if isinstance(o, Ival):
            return from_two(self.lo + o.lo, self.hi + o.hi)
        o = const(o)
        return from_two(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, Ival):
            return from_two(self.lo - o.hi, self.hi - o.lo)
        o = const(o)
        return from_two(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        o = const(o)
        return from_two(o.lo - self.hi, o.hi - self.lo)

    def __mul__(self, o):
        if not isinstance(o, Ival):
            o = const(o)
        a = self.lo * o.lo; b = self.lo * o.hi
        c = self.hi * o.lo; d = self.hi * o.hi
        # v11 (B): endpoint-based hull — pick min/max of the four
        # candidate LOWER/UPPER endpoints separately (handles
        # nonzero-radius Arb balls whose own lo/hi may overlap).
        lo_candidates = [x.lower() for x in (a, b, c, d)]
        hi_candidates = [x.upper() for x in (a, b, c, d)]
        lo_pick = lo_candidates[0]
        hi_pick = hi_candidates[0]
        for cand_lo in lo_candidates[1:]:
            if cand_lo < lo_pick:
                lo_pick = cand_lo
        for cand_hi in hi_candidates[1:]:
            if cand_hi > hi_pick:
                hi_pick = cand_hi
        return Ival(lo_pick, hi_pick)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, Ival):
            o = const(o)
        assert o.lo > 0 or o.hi < 0, "division by interval containing 0"
        cands = []
        for x in (self.lo, self.hi):
            for y in (o.lo, o.hi):
                cands.append(x / y)
        # v11 (B): endpoint-based hull for nonzero-radius divisor balls.
        lo_candidates = [x.lower() for x in cands]
        hi_candidates = [x.upper() for x in cands]
        lo_pick = lo_candidates[0]
        hi_pick = hi_candidates[0]
        for cand_lo in lo_candidates[1:]:
            if cand_lo < lo_pick:
                lo_pick = cand_lo
        for cand_hi in hi_candidates[1:]:
            if cand_hi > hi_pick:
                hi_pick = cand_hi
        return Ival(lo_pick, hi_pick)

    def pos(self):
        """|self| as interval.  v11 (B): endpoint-based hull for
        nonzero-radius balls — compute all candidates as Arb endpoints
        and pick exact min/max, never via ball-internal min/max."""
        if self.lo >= 0:
            return Ival(self.lo, self.hi)
        if self.hi <= 0:
            return Ival(-self.hi, -self.lo)
        # straddles zero: lower is exactly 0; upper = max(|lo|, |hi|)
        upper_candidates = [-self.lo, self.hi]
        upper_pick = upper_candidates[0]
        for cand in upper_candidates[1:]:
            if cand > upper_pick:
                upper_pick = cand
        return Ival(arb(0), upper_pick)

    def contains(self, x: float) -> bool:
        eps = 2.0 ** -52
        return float(self.lo) <= x + eps and x <= float(self.hi) + eps

    def low(self):
        """Outward-rounded lower endpoint (buffered down)."""
        return self.lo - arb(2) ** -OUTB

    def up(self):
        """Outward-rounded upper endpoint (buffered up)."""
        return self.hi + arb(2) ** -OUTB

    def __repr__(self):
        return f"I[{float(self.lo):.6e}, {float(self.hi):.6e}]"


def from_two(a, b):
    """From two (possibly ball/float) values: take exact endpoints.
    v11 fix (Main D3): the UPPER endpoint uses b.upper(); the previous
    b.lower() truncated the top of a nonzero-radius ball and could fail
    to enclose (regression plant: exp of a wide interval)."""
    lo = a.lower() if isinstance(a, arb) else arb(float(a))
    hi = b.upper() if isinstance(b, arb) else arb(float(b))
    return Ival(lo, hi)


def const(x):
    if isinstance(x, Ival):
        return x
    if isinstance(x, arb):
        return Ival(x.lower(), x.upper())
    f = float(x)
    assert np.isfinite(f)
    # exact dyadic point
    if f == 0.0:
        z = arb(0)
        return Ival(z, z)
    m, e = np.frexp(f)
    mi = int(np.ldexp(m, 53))
    ei = int(e) - 53
    v = arb(mi * (1 << ei)) if ei >= 0 else arb(mi) / arb(1 << (-ei))
    return Ival(v, v)


def point(dyadic_arb):
    return Ival(dyadic_arb, dyadic_arb)


# ----------------------------------------------------- transcendental ops

def ln_iv(x: Ival) -> Ival:
    """ln over interval; requires x.hi > 0. Encloses via ball eval."""
    assert x.hi.lower() > 0, "ln over nonpositive interval"
    lo_part = x.lo if x.lo.lower() > 0 else arb(0)
    with ctx.workprec(MID):
        lo_b = x.lo.log() if x.lo.lower() > 0 else None
        hi_b = x.hi.log()
    lo_e = (lo_b.lower() - arb(2) ** -OUTB if lo_b is not None
            else -(arb(10) ** 40))  # v11 fix: positive magnitude, negative sign
    hi_e = hi_b.upper() + arb(2) ** -OUTB
    return Ival(lo_e, hi_e)


def ent_iv_point(x: Ival) -> Ival:
    """interval for -t ln t, t >= 0 (concave). Compute at both endpoints via
    outward-rounded ln; add the interior maximum 1/e (max of -t ln t on
    t >= 0) when the interval covers 1/e."""
    zero = arb(0)
    if x.hi.lower() <= 0:
        return Ival(zero, zero)
    lo = x.lo if x.lo.lower() > 0 else zero
    hi = x.hi

    def f_at(tv):
        l = ln_iv(Ival(tv, tv))
        a = -(tv * l.lo)
        b = -(tv * l.hi)
        # v11 (B): endpoint-based lo/hi (exact arb comparisons).
        if a.lower() < b.lower():
            lo_e = a
        else:
            lo_e = b
        if a.upper() > b.upper():
            hi_e = a
        else:
            hi_e = b
        return Ival(lo_e.lower(), hi_e.upper())

    fa = f_at(lo) if lo > 0 else Ival(zero, zero)
    fb = f_at(hi)
    # v11 (B): endpoint-based lo/hi
    if fa.lo < fb.lo:
        lo_e = fa.lo
    else:
        lo_e = fb.lo
    if fa.hi > fb.hi:
        hi_e = fa.hi
    else:
        hi_e = fb.hi
    # v11 fix (Main D3, refine B): 1/e via pure Arb conservative overlap
    # on x.lo/x.hi (no float anywhere in the decision).
    with ctx.workprec(MID):
        one_over_e_ball = (arb(-1).exp())  # = 1/e via exp(-1)
        one_over_e = Ival(one_over_e_ball.lower() - arb(2) ** -OUTB,
                          one_over_e_ball.upper() + arb(2) ** -OUTB)
    # Conservative overlap: t = 1/e; include max iff [x.lo, x.hi] can
    # possibly contain t, using ONLY Arb bounds (no float).
    t = arb(-1).exp()
    if x.lo.lower() <= t.upper() and t.lower() <= x.hi.upper():
        hi_e = max(hi_e, one_over_e.hi)
    return Ival(lo_e, hi_e)


ZC_FLAG = {"hit": False}


def ent_vec(ws) -> Ival:
    """sum of -w_i ln w_i over list of Ivals.

    Gate-C B&B semantics (pre-registered): a box may push a coordinate's
    nominal value (all >= 1e-3 in the released vector) only to
    radius 1e-7, so a zero-crossing here certifies that the box left the
    simplex support; we avoid silently clipping by clamping the interval
    to [0, hi] (valid because -t ln t >= 0 on t >= 0 and the enclosure of
    -t ln t for t in [0, hi] is exactly that of [0, hi] with lo clamped;
    on t < 0 the true -t ln t is NOT defined for the probability program,
    so a crossing marks the box infeasible-by-support) and FLAG the hit so
    the report records the box as support-invalid."""
    total = Ival(arb(0), arb(0))
    for wi in ws:
        assert isinstance(wi, Ival), f"ent_vec not-Ival: {wi!r}"
        if wi.lo.lower() < 0:
            if wi.hi.lower() <= 0:
                continue  # entire weight <= 0: contributes its [0, hi] side
            ZC_FLAG["hit"] = True
            wi = Ival(arb(0), wi.hi)
        total = total + ent_iv_point(wi)
    return total


def nent_vec(ws, p: Ival) -> Ival:
    """sum -a_i ln(a_i/p)."""
    assert p.hi.lower() > 0
    total = Ival(arb(0), arb(0))
    ln_p = ln_iv(p)
    for wi in ws:
        if wi.hi.lower() <= 0:
            continue
        # -a ln a + a ln p
        t = ent_iv_point(wi)
        term1 = t
        term2 = wi * ln_p          # a*ln(p)
        total = total + term1 + term2
    return total


# ------------------------------------------------------------- tree load


class IntervalTree:
    def __init__(self, mat_path, q=5.0, K=1.0, max_level=3, float_module=None):
        self.mat = mat_path
        self.q = q
        self.K = K
        self.L = max_level
        FM = float_module if float_module is not None else F
        self.FM = FM
        pm = FM.ParamManager()
        ws = FM.Workspace(pm, q, K, omega=0.0, max_level=max_level)
        params = np.asarray(loadmat(mat_path)["params"]).flatten()
        assert len(params) == pm.num_input
        self.params_sha = hashlib.sha256(params.tobytes()).hexdigest()
        pm.set_value(params)
        self.pm = pm
        self.ws = ws
        start, size = pm.start, pm.size

        def get(gid):
            if getattr(pm, "box_lo", None) is not None:
                lo = pm.box_lo[start[gid]:start[gid] + size[gid]]
                hi = pm.box_hi[start[gid]:start[gid] + size[gid]]
                out = []
                for l, h in zip(lo, hi):
                    # interval_core project convention: endpoints are exact
                    # dyadics from float64 (IC.const), never float-converted
                    # arb balls — a ball-built endpoint has nonzero radius
                    # and breaks .lower() >= 0 guards in ent_vec.
                    lo_e = const(float(l)).lo
                    hi_e = const(float(h)).lo
                    out.append(Ival(lo_e, hi_e))
                return out
            return [point(const(v).lo) for v in pm.cur_x[start[gid]:start[gid] + size[gid]]]

        def get_scalar(gid):
            return float(pm.cur_x[start[gid]])

        pm.get = get
        pm.get_scalar = get_scalar
        FM.set_interval_pm(pm)

    def evaluate(self):
        self.ws.evaluate()



def iv_vec(lst_of_arb_points):
    return [point(v) for v in lst_of_arb_points]
