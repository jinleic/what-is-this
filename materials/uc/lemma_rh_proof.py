"""Machine-checked proof of rh(x) <= 2x on 0 <= x <= 1/2.

Put F(x) = ln(2) f(x).  The near-zero argument is analytic rather than an
interval evaluation at zero: exact symbolic algebra cancels the apparent
x*ln(x) scale, and an explicit logarithmic Taylor remainder gives a positive
lower bound.  Away from zero, Arb interval arithmetic covers the remaining
compact interval adaptively.  ``cert2.rho_upper`` and ``cert2._rho_global``
consume this half-interval theorem.  On ``[1/2,1]`` the companion bound
``rh(x) <= 2(1-x)`` follows directly from ``h <= 1`` and the nonnegative term
subtracted in the definition of ``rh``.
"""

from fractions import Fraction

import sympy as sp
from flint import arb, ctx


ctx.prec = 160

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
LN2 = TWO.log()
EPS_Q = Fraction(1, 16)
HALF_Q = Fraction(1, 2)
MAX_DEPTH = 40
MAX_PROCESSED = 2_000_000


def checked_print(condition, text, label):
    """Print a mathematical claim only after its arithmetic assertion passes."""
    assert bool(condition), label
    print(text, flush=True)


def qarb(value):
    """Convert a Fraction to an Arb enclosure without binary floats."""
    return arb(value.numerator) / arb(value.denominator)


# ---------------------------------------------------------------------------
# Part 1: exact series and an analytic, explicit remainder bound near zero.
# ---------------------------------------------------------------------------

sx = sp.symbols("x", positive=True)
sy = 1 - sx
one_minus_y2_identity = sp.simplify(1 - sy**2 - sx * (2 - sx)) == 0

# This is exactly ln(2) times the f in the theorem, with
# 1 - (1-x)^2 replaced by the identical x(2-x).
F_original = (
    -(sy**2 * sp.log(sy**2) + sx * (2 - sx) * sp.log(sx * (2 - sx)))
    + 2 * sy * (sx * sp.log(sx) + sy * sp.log(sy))
    + 2 * sx * sp.log(2)
)

# SymPy expands the positive factors.  The second replacement is the exact
# factorization 2-x = 2(1-x/2), valid on the proof domain.
F_expanded = sp.expand(sp.expand_log(F_original, force=True))
F_expanded = F_expanded.xreplace(
    {sp.log(2 - sx): sp.log(2) + sp.log(1 - sx / 2)}
)
F_exact = sx**2 * sp.log(2 / sx) - sx * (2 - sx) * sp.log(1 - sx / 2)
exact_identity = (
    sp.simplify(sp.expand_log(F_expanded - F_exact, force=True)) == 0
)

series_actual = sp.series(F_exact, sx, 0, 7).removeO().expand()
series_expected = (
    sx**2 * (-sp.log(sx) + sp.log(2) + 1)
    - sx**3 / 4
    - sx**4 / 24
    - sx**5 / 96
    - sx**6 / 320
)
series_identity = sp.simplify(series_actual - series_expected) == 0
leading_limit = sp.limit(F_exact / (sx**2 * sp.log(1 / sx)), sx, 0, dir="+")
little_o_limit = sp.limit((F_exact / sp.log(2)) / sx, sx, 0, dir="+")
leading_behavior = leading_limit == 1 and little_o_limit == 0

# Let rho(u) = -ln(1-u)-u.  Its integral representation proves, for
# 0 <= u <= U < 1,
#
#   0 <= rho(u) = integral_0^u t/(1-t) dt <= u^2/[2(1-U)].
#
# The following exact identities machine-check every algebraic step.  The sign
# follows because t, U-t, 1-U, and 1-t are all nonnegative on that domain.
su, st, sU = sp.symbols("u t U", nonnegative=True)
rho = -sp.log(1 - su) - su
rho_derivative_identity = sp.simplify(
    sp.diff(rho, su) - su / (1 - su)
) == 0
rho_origin_identity = rho.subs(su, 0) == 0
integrand_gap = st / (1 - sU) - st / (1 - st)
integrand_gap_factored = st * (sU - st) / ((1 - sU) * (1 - st))
integrand_order_identity = sp.simplify(
    integrand_gap - integrand_gap_factored
) == 0
integral_upper_identity = sp.simplify(
    sp.integrate(st / (1 - sU), (st, 0, su))
    - su**2 / (2 * (1 - sU))
) == 0

leading_form = sx**2 * (-sp.log(sx) + sp.log(2) + 1)
rho_at_x = -sp.log(1 - sx / 2) - sx / 2
remainder_form = -sx**3 / 2 + sx * (2 - sx) * rho_at_x
leading_remainder_identity = sp.simplify(
    sp.expand_log(F_exact - leading_form - remainder_form, force=True)
) == 0

# These identities expose the two remaining order steps as products of
# nonnegative factors on 0 < x <= eps: x(2-x) <= 2x, and the lower-bound
# bracket is decreasing, so its minimum is attained at eps.
product_upper_identity = sp.simplify(
    2 * sx - sx * (2 - sx) - sx**2
) == 0
near_lower_form = sx**2 * (sp.log(2 / sx) + 1 - sx / 2)
near_excess_identity = sp.simplify(
    sp.expand_log(
        F_exact - near_lower_form - sx * (2 - sx) * rho_at_x,
        force=True,
    )
) == 0
near_bracket = sp.log(2 / sx) + 1 - sx / 2
near_bracket_derivative_identity = sp.simplify(
    sp.diff(near_bracket, sx) + 1 / sx + sp.Rational(1, 2)
) == 0
eps_symbolic = sp.Rational(EPS_Q.numerator, EPS_Q.denominator)
u_domain_identity = sp.simplify(
    eps_symbolic / 2 - sx / 2 - (eps_symbolic - sx) / 2
) == 0
near_endpoint_identity = sp.simplify(
    near_bracket.subs(sx, eps_symbolic)
    - (sp.log(2 / eps_symbolic) + 1 - eps_symbolic / 2)
) == 0

EPS = qarb(EPS_Q)
U_MAX = EPS / TWO
remainder_lo = -ONE / TWO
remainder_hi = -ONE / TWO + ONE / (arb(4) * (ONE - U_MAX))
remainder_hi_exact = sp.simplify(
    -sp.Rational(1, 2)
    + 1 / (4 * (1 - sp.Rational(EPS_Q.numerator, EPS_Q.denominator) / 2))
) == -sp.Rational(15, 62)
remainder_constants_valid = (
    U_MAX > ZERO
    and U_MAX < ONE
    and remainder_lo <= remainder_hi
    and remainder_hi < ZERO
    and remainder_hi_exact
)

# From the remainder enclosure,
#   F >= x^2[-ln x + ln 2 + 1] - x^3/2
#     >= x^2[ln(2/eps) + 1 - eps/2] = K x^2.
# Every operation defining K remains in Arb.
K_NEAR = (TWO / EPS).log() + ONE - EPS / TWO
K_FOR_F = K_NEAR / LN2
near_chain_valid = (
    EPS > ZERO
    and EPS <= qarb(HALF_Q)
    and LN2 > ZERO
    and K_NEAR.is_finite()
    and K_NEAR > ZERO
    and K_FOR_F > ZERO
)


def h_point(z):
    """Binary entropy at an exact endpoint, with h(0)=h(1)=0."""
    z = arb(z)
    if z <= ZERO or z >= ONE:
        return ZERO
    return -(z * z.log() + (ONE - z) * (ONE - z).log()) / LN2


f_at_zero = h_point(ONE) - TWO * h_point(ZERO) + ZERO
endpoint_valid = f_at_zero == ZERO

near_checks = [
    ("1-(1-x)^2 = x(2-x)", one_minus_y2_identity),
    ("exact simplification of ln(2) f", exact_identity),
    ("SymPy series through order x^6", series_identity),
    ("leading asymptotic and little-o limit", leading_behavior),
    ("rho derivative identity", rho_derivative_identity),
    ("rho value at zero", rho_origin_identity),
    ("integrand comparison factorization", integrand_order_identity),
    ("integrated remainder upper bound", integral_upper_identity),
    ("leading-plus-remainder identity", leading_remainder_identity),
    ("x(2-x) upper-bound factorization", product_upper_identity),
    ("direct near-zero lower-bound identity", near_excess_identity),
    ("near-zero bracket derivative", near_bracket_derivative_identity),
    ("u=x/2 respects the epsilon domain", u_domain_identity),
    ("near-zero bracket endpoint", near_endpoint_identity),
    ("remainder constants on 0 < x <= eps", remainder_constants_valid),
    ("final Arb lower-bound chain", near_chain_valid),
    ("endpoint x=0", endpoint_valid),
]
near_gap = None
for near_label, near_condition in near_checks:
    try:
        assert bool(near_condition), near_label
    except AssertionError:
        near_gap = near_label
        break
near_ok = near_gap is None

if near_ok:
    checked_print(
        exact_identity and one_minus_y2_identity,
        "SYMPY EXACT IDENTITY: ln(2) f(x) = x^2 ln(2/x) "
        "- x(2-x) ln(1-x/2).",
        "exact identity did not pass",
    )
    checked_print(
        series_identity,
        "SYMPY SERIES: ln(2) f(x) = x^2[-ln(x)+ln(2)+1] "
        "- x^3/4 - x^4/24 - x^5/96 - x^6/320 + O(x^7).",
        "series identity did not pass",
    )
    checked_print(
        leading_behavior,
        "LEADING ORDER: f(x) ~ x^2 log_2(1/x); the apparent x ln(x) "
        "terms cancel and f(x) = o(x).",
        "leading behavior did not pass",
    )
    checked_print(
        rho_derivative_identity
        and rho_origin_identity
        and integrand_order_identity
        and integral_upper_identity
        and leading_remainder_identity
        and product_upper_identity
        and u_domain_identity
        and remainder_constants_valid,
        "EXPLICIT REMAINDER ON 0 < x <= 1/16: writing ln(2)f as its "
        "displayed x^2 leading expression plus R(x), "
        "-x^3/2 <= R(x) <= -(15/62)x^3.",
        "explicit remainder enclosure did not pass",
    )
    checked_print(
        near_chain_valid
        and near_excess_identity
        and rho_derivative_identity
        and rho_origin_identity
        and near_bracket_derivative_identity
        and near_endpoint_identity,
        "NEAR ZERO: ln(2)f(x) >= K x^2 on 0 < x <= 1/16, "
        "where the certified Arb K is " + str(K_NEAR.lower()) + " > 0.",
        "near-zero Arb lower bound did not pass",
    )
    checked_print(
        endpoint_valid,
        "ENDPOINT: f(0) = 0 exactly.",
        "endpoint identity did not pass",
    )
else:
    checked_print(
        not near_ok and near_gap is not None,
        "PART 1 FAILED: exact unchecked gap is '" + str(near_gap) + "'.",
        "near-zero failure was not identified",
    )


# ---------------------------------------------------------------------------
# Part 2: adaptive direct Arb interval evaluation on [eps, 1/2].
# ---------------------------------------------------------------------------


def f_direct_enclosure(lo_q, hi_q):
    """Direct natural-log interval enclosure of f on [lo_q, hi_q]."""
    lo = qarb(lo_q)
    hi = qarb(hi_q)
    assert lo > ZERO and hi < ONE and lo <= hi, "input endpoints not in (0,1)"
    x = lo.union(hi)
    y = ONE - x
    y2 = y * y
    x2c = ONE - y2  # equals x(2-x) pointwise, with a sharper positive ball
    assert x > ZERO and x < ONE, "x ball not strictly inside (0,1)"
    assert y > ZERO and y < ONE, "1-x ball not strictly inside (0,1)"
    assert y2 > ZERO and y2 < ONE, "(1-x)^2 ball not inside (0,1)"
    assert x2c > ZERO and x2c < ONE, "x(2-x) ball not inside (0,1)"
    natural = (
        -(y2 * y2.log() + x2c * x2c.log())
        + TWO * y * (x * x.log() + y * y.log())
        + TWO * x * LN2
    )
    enclosure = natural / LN2
    assert enclosure.is_finite(), "direct f enclosure is non-finite"
    return enclosure


stack = [(EPS_Q, HALF_Q, 0)]
accepted = []
accepted_lowers = []
accepted_count = 0
processed_count = 0
min_certified_lower = None
bulk_gap = None

while stack:
    lo_q, hi_q, depth = stack.pop()
    processed_count += 1
    try:
        enclosure = f_direct_enclosure(lo_q, hi_q)
    except AssertionError as exc:
        bulk_gap = "non-finite or out-of-domain interval at [{}, {}]: {}".format(
            lo_q, hi_q, exc
        )
        break

    if enclosure > ZERO:
        lower = enclosure.lower()
        assert lower > ZERO
        accepted.append((lo_q, hi_q))
        accepted_lowers.append(lower)
        accepted_count += 1
        if min_certified_lower is None or lower < min_certified_lower:
            min_certified_lower = lower
        else:
            assert lower >= min_certified_lower
    else:
        if depth >= MAX_DEPTH:
            bulk_gap = "depth cap reached at [{}, {}]".format(lo_q, hi_q)
            break
        if processed_count >= MAX_PROCESSED:
            bulk_gap = "processed-interval cap reached with uncleared work"
            break
        mid_q = (lo_q + hi_q) / 2
        assert lo_q < mid_q < hi_q
        stack.append((mid_q, hi_q, depth + 1))
        stack.append((lo_q, mid_q, depth + 1))

    if processed_count % 100000 == 0:
        checked_print(
            processed_count >= accepted_count >= 0
            and accepted_count == len(accepted_lowers)
            and all(lower_bound > ZERO for lower_bound in accepted_lowers),
            "BULK PROGRESS: processed {} intervals and certified {} leaf "
            "subintervals so far.".format(processed_count, accepted_count),
            "bulk progress counters are inconsistent",
        )

ordered = sorted(accepted)
coverage_valid = bool(ordered)
if coverage_valid:
    coverage_valid = ordered[0][0] == EPS_Q and ordered[-1][1] == HALF_Q
if coverage_valid:
    for left, right in zip(ordered, ordered[1:]):
        if left[1] != right[0]:
            coverage_valid = False
            break

all_cleared_valid = (
    accepted_count == len(accepted_lowers)
    and all(lower_bound > ZERO for lower_bound in accepted_lowers)
)
minimum_record_valid = (
    min_certified_lower is not None
    and all(
        min_certified_lower <= lower_bound
        for lower_bound in accepted_lowers
    )
)

bulk_ok = (
    bulk_gap is None
    and not stack
    and coverage_valid
    and accepted_count == len(accepted)
    and all_cleared_valid
    and minimum_record_valid
    and min_certified_lower is not None
    and min_certified_lower > ZERO
)

if bulk_ok:
    checked_print(
        bulk_ok,
        "BULK: adaptive Arb coverage proved f(x) > 0 on [1/16, 1/2] "
        "with {} certified subintervals.".format(accepted_count),
        "bulk coverage did not pass",
    )
    checked_print(
        bulk_ok and minimum_record_valid and min_certified_lower > ZERO,
        "BULK MINIMUM CERTIFIED LOWER BOUND: "
        + str(min_certified_lower)
        + ".",
        "bulk minimum lower bound was not positive",
    )
else:
    if bulk_gap is None:
        bulk_gap = "accepted leaves did not form an exact contiguous cover"
    checked_print(
        not bulk_ok and bulk_gap is not None,
        "PART 2 FAILED: exact unchecked gap is '" + str(bulk_gap) + "'.",
        "bulk failure was not identified",
    )


# ---------------------------------------------------------------------------
# Formal conclusion, emitted only behind the two successful proof gates.
# ---------------------------------------------------------------------------

proof_complete = near_ok and bulk_ok
if proof_complete:
    checked_print(
        endpoint_valid and near_ok and bulk_ok,
        "THEOREM: f(x) := h((1-x)^2) - 2(1-x)h(x) + 2x >= 0 "
        "for every x in [0, 1/2].",
        "the theorem gates did not all pass",
    )
    checked_print(
        endpoint_valid and near_ok and bulk_ok,
        "PROVED DOMAIN SPLIT: f(0)=0; f(x)>0 on (0,1/16] by the "
        "explicit remainder bound; f(x)>0 on [1/16,1/2] by Arb coverage.",
        "domain split gates did not all pass",
    )
    checked_print(
        proof_complete,
        "PROOF COMPLETE: machine-checked",
        "proof was not complete",
    )
else:
    checked_print(
        not proof_complete,
        "PARTIAL RESULT ONLY: the theorem on [0,1/2] was not claimed.",
        "partial-result gate was inconsistent",
    )
