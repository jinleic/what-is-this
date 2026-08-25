"""CERTIFIED strict mixing gap between the two c* equality mechanisms.

At mean c* the Sawin/Cambie functional has two zero-cost laws:

  bad entropy-active law:  mu* = a* delta_1 + (1-a*) delta_b*,
  deterministic law:      nu  = c* delta_1 + (1-c*) delta_0.

Both have mean c*, and F(mu*) = F(nu) = 0.  `reimer_entropy.py` shows that a
fixed coordinate order can balance the Reimer deficit of mu* with the free
surplus of nu, so scalar Reimer information cannot break c*.  This file proves
that EVERY nontrivial mixture inside one conditional-probability law pays a
uniform quadratic gap:

    F((1-theta)nu + theta mu*) >= 0.1475 theta(1-theta),  0<=theta<=1.

Scope is essential.  Averaging fixed-order entropy proofs produces
E_pi F(mu_{i,pi}), not F(E_pi mu_{i,pi}).  If changing the order merely switches
one local law from nu to mu*, the averaged endpoint cost is still zero; random
order alone does not invoke this theorem.  The gap becomes operative only after
a structural argument forces a nontrivial mixture inside some mu_{i,pi}, or a
separate inequality compares the direct order average with such mixtures.

What remains is therefore structural and quantitative: prove within-order
mixing from join-homomorphic dependencies, or bound the direct permutation
average away from the endpoint equality mechanisms.

Why the proof is finite.  The mixture lives on atoms {0,b*,1}.  Write their
weights as x,y,z.  The coupled cost matrix is

              0       b       1
        0     0      h(b)     0
        b    h(b)     1       0
        1     0       0       0.

An optimal transportation coupling greedily pairs b-mass with atom 1 (cost 0),
then atom 0 (cost h(b)<1), and uses b-b (cost 1) only when non-b capacity is
exhausted.  Exchange arguments make this exact: replacing any b-b + 0-1 pair by
b-1 + 0-b changes cost from 1 to h(b)<1; replacing b-0 + 1-anything by b-1
never increases cost.  Therefore

    C(theta) = 0                           if y <= z,
             = 2(y-z)h(b*)                if z < y <= 1/2,
             = (2y-1) + 2x h(b*)          if y > 1/2.

The switches are theta_1=c*/(1-2a*+c*) and theta_2=1/[2(1-a*)].  IID entropy is
Q=2xy h(b*)+y^2 h(2b*-b*^2), and L=y h(b*), so on each regime
F-0.1475 theta(1-theta) is a quadratic.  Exact rational endpoints bracket the
unique defining root by an Arb sign change and a positive derivative.  The
active endpoint zero follows algebraically from the root equation (Q=C=L).
Every regime is certified concave, so its minimum is controlled by those exact
endpoint identities and positive switch values, not a grid approximation.

Run: ./.venv/bin/python uc/shapley_gap.py
"""

from flint import arb, ctx
import sympy as sp

ctx.prec = 384
ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
LOG2 = TWO.log()
ALPHA = arb("0.0356069")
KAPPA = arb("0.1475")


def h(x):
    if x == ZERO or x == ONE:
        return ZERO
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def hprime(x):
    return ((ONE - x) / x).log() / LOG2


def root_residual(x):
    hx = h(x)
    return hx * (TWO - hx) - h(TWO * x - x * x)


def root_derivative(x):
    hx = h(x)
    y = TWO * x - x * x
    return TWO * (ONE - hx) * hprime(x) \
        - TWO * (ONE - x) * hprime(y)


# Exact rational decimal endpoints.  Arb proves a sign change and a strictly
# positive derivative throughout this interval, hence one unique defining root.
B_LO = arb(
    "0.329454738503036972391705383877134130751839116916698512937367708874881240708840123088556565701276198872"
)
B_HI = arb(
    "0.329454738503036972391705383877134130751839116916698512937367708874881240708840123088556565701276198874"
)
ROOT_LO_RESIDUAL = root_residual(B_LO)
ROOT_HI_RESIDUAL = root_residual(B_HI)
BSTAR = B_LO.union(B_HI)
ROOT_DERIVATIVE = root_derivative(BSTAR)
assert ROOT_LO_RESIDUAL.upper() < 0
assert ROOT_HI_RESIDUAL.lower() > 0
assert ROOT_DERIVATIVE.lower() > 0

HB = h(BSTAR)
H2 = h(TWO * BSTAR - BSTAR * BSTAR)

# At the exact root H2=HB(2-HB).  With w=1/(2-HB), the active endpoint has
# Q=w^2 H2=w HB=L and C=2w-1=w HB=L, so F=(1-alpha)Q+alpha C-L=0 exactly.
_hs, _alpha = sp.symbols("h alpha")
_w = 1 / (2 - _hs)
_q_eq = _w ** 2 * _hs * (2 - _hs)
_c_eq = 2 * _w - 1
_l_eq = _w * _hs
assert sp.cancel(_q_eq - _l_eq) == 0
assert sp.cancel(_c_eq - _l_eq) == 0
assert sp.cancel((1 - _alpha) * _q_eq + _alpha * _c_eq - _l_eq) == 0
ENDPOINT_IDENTITY_PROVED = True
ASTAR = ONE - ONE / (TWO - HB)
CSTAR = ASTAR + (ONE - ASTAR) * BSTAR
Y1 = ONE - ASTAR
X0 = ONE - CSTAR
THETA1 = CSTAR / (ONE - TWO * ASTAR + CSTAR)
THETA2 = ONE / (TWO * Y1)


def coupled(theta, regime):
    x = (ONE - theta) * X0
    y = theta * Y1
    z = (ONE - theta) * CSTAR + theta * ASTAR
    if regime == 1:
        return ZERO
    if regime == 2:
        return TWO * (y - z) * HB
    if regime == 3:
        return TWO * y - ONE + TWO * x * HB
    raise ValueError(regime)


def gap(theta, regime):
    x = (ONE - theta) * X0
    y = theta * Y1
    q = TWO * x * y * HB + y * y * H2
    ell = y * HB
    functional = (ONE - ALPHA) * q + ALPHA * coupled(theta, regime) - ell
    return functional - KAPPA * theta * (ONE - theta)


def quadratic_coefficients(regime):
    """Recover c+b*t+a*t^2 exactly in Arb arithmetic from t=0,1,-1."""
    c = gap(ZERO, regime)
    p1 = gap(ONE, regime)
    pm = gap(-ONE, regime)
    b = (p1 - pm) / TWO
    a = (p1 + pm) / TWO - c
    return a, b, c


def lower_point(value):
    return value.lower()


def certify_regime(regime, lo, hi):
    a, b, c = quadratic_coefficients(regime)
    candidates = [("left", lo, gap(lo, regime)),
                  ("right", hi, gap(hi, regime))]
    # Float guidance is harmless: the candidate value is reevaluated as an Arb
    # ball and endpoint derivative signs certify cases where the vertex is out.
    af, bf = float(a.mid()), float(b.mid())
    if abs(af) > 1e-30:
        tv = -bf / (2.0 * af)
        if float(lo.mid()) < tv < float(hi.mid()):
            tball = arb(str(tv))
            candidates.append(("vertex", tball, gap(tball, regime)))
    minimum = min(candidates, key=lambda item: float(item[2].lower()))
    dlo = TWO * a * lo + b
    dhi = TWO * a * hi + b

    # Every regime is concave (a<0), so its minimum on a closed interval is
    # attained at an endpoint.  This is stronger and correct; an earlier draft
    # incorrectly required monotonicity in regimes 1 and 3.
    assert a.upper() < 0
    if regime == 1:
        assert gap(ZERO, regime) == ZERO
        assert gap(hi, regime).lower() > 0
    elif regime == 2:
        assert gap(lo, regime).lower() > 0
        assert gap(hi, regime).lower() > 0
    elif regime == 3:
        assert gap(lo, regime).lower() > 0
        assert ENDPOINT_IDENTITY_PROVED
    else:
        raise ValueError(regime)
    return a, b, c, candidates, minimum, dlo, dhi

if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    print("PROVED [rational root bracket + monotonicity + exact identities]")
    print("  g(B_LO) =", ROOT_LO_RESIDUAL, "< 0")
    print("  g(B_HI) =", ROOT_HI_RESIDUAL, "> 0")
    print("  g'(B_LO..B_HI) =", ROOT_DERIVATIVE, "> 0")
    print("  unique b* bracket =", BSTAR)
    print("  exact endpoint identity: Q=C=L from g(b*)=0")
    print("  a* =", ASTAR)
    print("  c* =", CSTAR)
    print("  h(b*) < 1:", HB, "; certified:", HB.upper() < 1)
    print("  theta_1 (y=z)   =", THETA1)
    print("  theta_2 (y=1/2) =", THETA2)
    assert ZERO < THETA1 < THETA2 < ONE
    assert HB.upper() < 1
    assert gap(ZERO, 1) == ZERO
    assert ENDPOINT_IDENTITY_PROVED
    # Consistency only; endpoint equality is the exact algebra above.
    assert gap(ONE, 3).contains(0)

    intervals = ((1, ZERO, THETA1),
                 (2, THETA1, THETA2),
                 (3, THETA2, ONE))
    for regime, lo, hi in intervals:
        a, b, c, candidates, minimum, dlo, dhi = certify_regime(regime, lo, hi)
        print()
        print("PROVED [Arb quadratic regime %d] interval %s .. %s" %
              (regime, lo, hi))
        print("  coefficients a,b,c:", a, b, c)
        print("  derivative at endpoints:", dlo, dhi)
        for name, theta, value in candidates:
            print("  candidate %-6s theta=%s gap=%s" % (name, theta, value))
        print("  minimum candidate:", minimum[0], minimum[2])

    # Sharp candidate (NUMERICAL, not needed for the certified 0.1475 theorem).
    sharp = gap(THETA1, 1) + KAPPA * THETA1 * (ONE - THETA1)
    sharp_ratio = sharp / (THETA1 * (ONE - THETA1))
    print()
    print("NUMERICAL orientation (Arb-enclosed evaluation of exact expression)")
    print("  sharp ratio candidate at theta_1 =", sharp_ratio)
    print()
    print("NONCOMMUTATION CHECK")
    half = ONE / TWO
    external_average = ZERO  # Exact: both endpoint costs vanish algebraically.
    endpoint_consistency = gap(ONE, 3)
    internal_mixture = gap(half, 2) + KAPPA * half * (ONE - half)
    assert ENDPOINT_IDENTITY_PROVED
    assert endpoint_consistency.contains(0)
    assert internal_mixture.lower() > KAPPA / arb(4)
    print("  (F(nu)+F(mu*))/2 =", external_average, "[exact identity]")
    print("  F((nu+mu*)/2)    =", internal_mixture)
    print("  These are different operations; order averaging gives the first.")
    print()
    print("FINAL VERDICT")
    print("PROVED: F((1-theta)nu + theta mu*) >= 0.1475 theta(1-theta)")
    print("for every theta in [0,1], with equality only at theta=0 or 1.")
    print("SCOPE: role switching across orders does not itself earn this gap.")
    print("The frontier is a direct order-average or within-order mixing theorem")
    print("for join-homomorphic deterministic dependencies.")
