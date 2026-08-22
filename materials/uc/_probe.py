from flint import arb, ctx

ctx.prec = 160
ONE = arb(1)
TWO = arb(2)
HALF = arb(1) / 2
LOG2 = TWO.log()


def hull(*bs):
    o = bs[0]
    for b in bs[1:]:
        o = o.union(b)
    return o


def h_pt(x):
    x = arb(x)
    if x <= 0 or x >= 1:
        return arb(0)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def h_encl(u, v):
    u, v = arb(u), arb(v)
    e = hull(h_pt(u), h_pt(v))
    return e if (v <= HALF or u >= HALF) else hull(e, ONE)


def amax(a, b):
    a, b = arb(a), arb(b)
    if a >= b:
        return a
    if b >= a:
        return b
    return hull(a, b)


def amin(a, b):
    a, b = arb(a), arb(b)
    if a <= b:
        return a
    if b <= a:
        return b
    return hull(a, b)


def rh_encl(u, v, gmax=None):
    u, v = arb(u), arb(v)
    p = hull(u, v)
    hp = h_encl(u, v)
    xl, xh = (ONE - v) * (ONE - v), (ONE - u) * (ONE - u)
    h2 = h_encl(xl, xh)
    e = TWO * (ONE - p) * hp - h2
    lo = amax(e.lower(), arb(0))
    hi = e.upper() if gmax is None else amin(e.upper(), gmax)
    return hull(lo, hi)


def sqrt_rh_enc(u, v, gmax=None):
    hi = rh_encl(u, v, gmax).upper()
    return hull(arb(0), hi.sqrt())


GMAX = arb("0.234318694")
print("=== full trace of rh_encl(0.95, 1.0) ===")
u, v = arb("0.95"), arb(1.0)
p = hull(u, v)
hp = h_encl(u, v)
xl, xh = (ONE - v) * (ONE - v), (ONE - u) * (ONE - u)
h2 = h_encl(xl, xh)
e = TWO * (ONE - p) * hp - h2
print("  p  =", p)
print("  hp =", hp, " lower=%.3e upper=%.3e" % (float(hp.lower()), float(hp.upper())))
print("  h2 =", h2, " lower=%.3e upper=%.3e" % (float(h2.lower()), float(h2.upper())))
print("  e  =", e, " lower=%.3e upper=%.3e" % (float(e.lower()), float(e.upper())))
r = rh_encl("0.95", 1.0, GMAX)
print("  rh_encl =", r, " lower=%.3e upper=%.3e" % (float(r.lower()), float(r.upper())))
print("  e.lower() is finite?", (e.lower() == e.lower()))
print("  r.upper() is finite?", (r.upper() == r.upper()))
hi = r.upper()
print("  hi =", hi, " finite?", (hi == hi))
s = hi.sqrt()
print("  hi.sqrt =", s, " finite?", (s == s))

print("\n=== sqrt_rh_enc(0.95,1.0) ===")
val = sqrt_rh_enc("0.95", 1.0, GMAX)
print("  =", val, " finite?", (val == val),
      " lower=%.4e upper=%.4e" % (float(val.lower()), float(val.upper())))

print("\n=== does NaN come from hull of a straddling ball's sqrt? ===")
straddle = hull(arb(0), arb("0.0286"))   # dips ~1e-10 below 0
print("  straddle =", straddle, " lower=%.3e" % float(straddle.lower()))
print("  straddle.sqrt =", straddle.sqrt(), " finite?",
      (straddle.sqrt() == straddle.sqrt()))
