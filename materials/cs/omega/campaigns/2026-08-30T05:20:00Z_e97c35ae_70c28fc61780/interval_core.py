"""Stage (b) core v2 — endpoint-pair interval representation.

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

import sys, hashlib, contextlib
import numpy as np
from scipy.io import loadmat
from flint import arb, ctx

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
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
        assert float(lo) <= float(hi), f"inverted Ival [{lo}, {hi}]"
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
        return from_two(min(a, b, c, d), max(a, b, c, d))

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, Ival):
            o = const(o)
        assert o.lo > 0 or o.hi < 0, "division by interval containing 0"
        cands = []
        for x in (self.lo, self.hi):
            for y in (o.lo, o.hi):
                cands.append(x / y)
        return from_two(min(cands), max(cands))

    def pos(self):
        """|self| as interval."""
        if self.lo >= 0:
            return Ival(self.lo, self.hi)
        if self.hi <= 0:
            return Ival(-self.hi, -self.lo)
        return Ival(arb(0), max(-self.lo, self.hi))

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
    """From two (possibly ball/float) values: take exact endpoints."""
    lo = a.lower() if isinstance(a, arb) else arb(float(a))
    hi = b.lower() if isinstance(b, arb) else arb(float(b))
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
    lo_e = lo_b.lower() - arb(2) ** -OUTB if lo_b is not None else arb(-10) ** 40
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
        lo_e = min(a, b)
        hi_e = max(a, b)
        return Ival(lo_e.lower(), hi_e.upper())

    fa = f_at(lo) if float(lo) > 0 else Ival(zero, zero)
    fb = f_at(hi)
    lo_e = min(fa.lo, fb.lo)
    hi_e = max(fa.hi, fb.hi)
    with ctx.workprec(MID):
        one_over_e = (arb(1) / arb(np.e))
    if x.lo.lower() <= 0.36787944117144233 <= x.hi.upper():
        hi_e = max(hi_e, one_over_e.upper() + arb(2) ** -OUTB)
    return Ival(lo_e, hi_e)


def ent_vec(ws) -> Ival:
    """sum of -w_i ln w_i over list of Ivals (each >= 0)."""
    total = Ival(arb(0), arb(0))
    for wi in ws:
        assert isinstance(wi, Ival) and wi.lo.lower() >= 0, \
            f"ent_vec negative/prob not-Ival: {wi!r}"
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
