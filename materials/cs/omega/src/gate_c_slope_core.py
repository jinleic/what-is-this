"""Gate C stage 3 —-certified slope (forward-mode interval AD) through the
transcribed VXXZ24 tree. Implements pre_statement_gateC2.md.

Representation: Slope(val, grad)
  val  : IC.Ival  — enclosure of the quantity over the CURRENT box
  grad : np.ndarray(object, shape (n_slopes,)) of IC.Ival — enclosure of the
         derivative of the quantity w.r.t. each slope coordinate, over the
         SAME box (slope bundle: (f(x) - f(c))/(x - c) style enclosure; for
         leaves val box = centre ± r, grad is the identity — a genuine
         derivative at the centre is the degenerate radius-0 case).

All arithmetic keeps interval endpoints outward-rounded inside IC (arb
endpoint pairs, buffered low()/up()); every log/entropy applies the chain
rule with a certified derivative enclosure over the value interval.
NO dual feasibility, NO float approximations inside the certified pass:
the float layer is used ONLY for fixtures and cross-checks.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import interval_core as IC
from flint import arb

I = IC.Ival
ZERO = None  # set after import: I(arb(0), arb(0))
ONE = None


def init_units(n):
    global ZERO, ONE
    ZERO = I(arb(0), arb(0))
    ONE = I(arb(1), arb(1))


class Slope:
    __slots__ = ("val", "grad")
    def __init__(self, val, grad):
        assert isinstance(val, I), f"Slope.val must be Ival, got {type(val)}"
        self.val = val
        self.grad = grad


def slopify(z, n, wid):
    """Coerce a leaf quantity (IC value modes used by vxxz24_float interval
    path) to Slope with zero gradient."""
    if isinstance(z, Slope):
        return z
    if isinstance(z, I):
        return Slope(z, np.array([ZERO] * n, dtype=object))
    if isinstance(z, arb):
        z = IC.Ival(z.lower(), z.upper())
        return Slope(z, np.array([ZERO] * n, dtype=object))
    if isinstance(z, list):
        return [slopify(x, n, wid) for x in z]
    # float / numpy scalar
    return Slope(IC.const(float(z)), np.array([ZERO] * n, dtype=object))


def s_const(x, n):
    return Slope(IC.const(x), np.array([ZERO] * n, dtype=object))
def init_units(n):
    global ZERO, ONE
    ZERO = I(arb(0), arb(0))
    ONE = I(arb(1), arb(1))


def gzeros(n):
    return np.array([ZERO] * n, dtype=object)


def gscale(g, iv):
    """Ival * grad-array elementwise (grad entries are Ivals)."""
    zz = I(arb(0), arb(0))
    return np.array([(iv * x) if True else None for x in g], dtype=object)


def gadd(g1, g2):
    return np.array([x + y for x, y in zip(g1, g2)], dtype=object)


def gneg(g):
    z = I(arb(0), arb(0))
    return np.array([z - x for x in g], dtype=object)


def gdiv(g, iv):
    return np.array([x / iv for x in g], dtype=object)
def s_add(a, b):
    if isinstance(a, Slope) and isinstance(b, Slope):
        return Slope(a.val + b.val, gadd(a.grad, b.grad))
    if isinstance(a, Slope):
        return Slope(a.val + IC.const(b), a.grad)
    if isinstance(b, Slope):
        return Slope(IC.const(a) + b.val, b.grad)
    return IC.const(a) + IC.const(b)


def s_sub(a, b):
    if isinstance(a, Slope) and isinstance(b, Slope):
        return Slope(a.val - b.val, gadd(a.grad, gneg(b.grad)))
    if isinstance(a, Slope):
        return Slope(a.val - IC.const(b), a.grad)
    if isinstance(b, Slope):
        return Slope(IC.const(a) - b.val, gneg(b.grad))
    return IC.const(a) - IC.const(b)


def s_mul(a, b):
    if not isinstance(a, Slope) and not isinstance(b, Slope):
        return IC.const(a) * IC.const(b)
    if not isinstance(a, Slope):
        return Slope(IC.const(a) * b.val, gscale(b.grad, IC.const(a)))
    if not isinstance(b, Slope):
        return Slope(a.val * IC.const(b), gscale(a.grad, IC.const(b)))
    # product rule: df = x dy + y dx (Leibniz)
    return Slope(a.val * b.val,
                 gadd(gscale(b.grad, a.val), gscale(a.grad, b.val)))


def s_div(a, b):
    if not isinstance(a, Slope) and not isinstance(b, Slope):
        return IC.const(a) / IC.const(b)
    if not isinstance(b, Slope):
        return Slope(a.val / IC.const(b), gdiv(a.grad, IC.const(b)))
    if not isinstance(a, Slope):
        f = IC.const(a) / b.val
        neg_fb = IC.Ival(-(f / b.val).hi, -(f / b.val).lo)
        return Slope(f, gscale(b.grad, neg_fb))
    # f = a / b: df = (da − f·db)/b ; negate f as interval: [−hi, −lo]
    f = a.val / b.val
    neg_f = IC.Ival(-f.hi, -f.lo)
    return Slope(f, gdiv(gadd(a.grad, gscale(b.grad, neg_f)), b.val))


def s_inv(b):
    f = IC.Ival(arb(1), arb(1)) / b.val
    neg_fb = IC.Ival(-(f / b.val).hi, -(f / b.val).lo)
    return Slope(f, gscale(b.grad, neg_fb))


def s_pos(a):
    """|a| enclosures of value AND gradient (chain rule: sign(x) in [-1,1])."""
    v = a.val.pos()
    if a.val.lo >= 0:
        s = I(arb(1), arb(1))
    elif a.val.hi <= 0:
        s = I(arb(-1), arb(-1))
    else:
        s = I(arb(-1), arb(1))
    return Slope(v, gscale(a.grad, s))


def s_ent(w):
    """f = -t ln t on t >= 0: value 0 at t = 0; f' = -(ln t + 1) certified
    over t-box (strictly positive t)."""
    t = w.val
    if t.hi.lower() <= 0:
        return Slope(IC.Ival(arb(0), arb(0)), gzeros(len(w.grad)))
    l = IC.ln_iv(t)
    d = IC.Ival((-l.hi) - 1, (-l.lo) - 1)
    return Slope(IC.ent_iv_point(t), gscale(w.grad, d))


def s_nent(w, p):
    """f = Σ_i -a_i ln(a_i / p) for vector-over-slopes w and Slope p.
    In vxxz24_float.G.normalized_entropy the vector path is
    IC.nent_vec(list_of_Ivals, p) = Σ [ent(a_i) + a_i ln p]  (a_i >= 0).
    Grad: w.r.t. a_i: -(ln(a_i/p) + 1) = -(ln a_i + 1) + ln p
          w.r.t. p:    Σ_i a_i / p.
    Note: in the tree, nent_vec is applied to a python-LIST of slopes of the
    same grad length; p may be a Slope (when p is dist-coupled) or const."""
    ws, pw = w, p
    total = Slope(ZERO, gzeros(len(ws[0].grad)))
    ln_p = pw.val if isinstance(pw, Slope) else IC.ln_iv(pw)
    for a in ws:
        if a.val.hi.lower() <= 0:
            continue
        ea = IC.ent_iv_point(a.val)
        val_term = ea + a.val * ln_p
        la = IC.ln_iv(a.val)
        da = IC.Ival((-la.hi) - 1, (-la.lo) - 1) + ln_p
        total = Slope(total.val + val_term, gadd(total.grad, gscale(a.grad, da)))
    if isinstance(pw, Slope):
        asum = ZERO
        for a in ws:
            if a.val.hi.lower() <= 0:
                continue
            asum = asum + a.val
        gp = (asum / pw.val)  # d/dp Σ a_i ln p = Σ a_i / p
        total = Slope(total.val, gadd(total.grad, gscale(pw.grad, gp)))
    return total


def s_ln_max_form(dm):
    """IC.ln_iv on a vector (list) — vxxz24 uses IC.ln_iv(dm[i]) elementwise."""
    return IC.ln_iv(dm)


# ------------------------------------------------------------- vector ops

def kron_csd_s(a, b, n, wid):
    """vxxz24_float kron_csd on Slope lists."""
    from vxxz24_float import DecodeCSD
    len_a, len_b = len(a), len(b)
    out = [Slope(ZERO, np.array([ZERO] * n, dtype=object)) for _ in range(len_a * len_b)]
    total = 0
    for ia in range(len_a):
        arr = DecodeCSD(ia, total) if False else None
    # fall back: mirror vxxz24_float.kron_csd arithmetic
    import vxxz24_float as F
    return F.kron_csd(a, b)


def j2m_apply_s(mats, dist_slopes, n):
    """JointToMargin.apply for slopes: each output coord j is a fixed 0/1
    sum over input coords (exact dyadic coefficients)."""
    out = []
    for m in mats:
        v = []
        for j in range(m.shape[1]):
            rows = np.nonzero(m[:, j])[0]
            if len(rows) == 0:
                v.append(Slope(ZERO, np.array([ZERO] * n, dtype=object)))
            else:
                acc = dist_slopes[rows[0]]
                for i in rows[1:]:
                    acc = s_add(acc, dist_slopes[i])
                v.append(acc)
        out.append(v)
    return out
