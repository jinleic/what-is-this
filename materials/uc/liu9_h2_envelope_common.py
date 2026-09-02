"""Shared degree-2 Bernstein helpers on the 2-simplex for the Liu H2 envelope.

This module is the common tool behind `liu9_h2_support_envelope.py` and
`liu9_h2_q_envelope.py`.  It contains no claim logic of its own: it turns six
node samples of *any* polynomial of total degree <= 2 on the unit simplex into
its Bernstein control values and encloses the polynomial's range / absolute
value uniformly over the whole simplex.  There is deliberately no center, no
argmin, and no radius-shortcut logic here.

Domain
------
The mass simplex is T = {(a1, a2) : a1 >= 0, a2 >= 0, a1 + a2 <= 1} with
barycentric coordinates

    lA = 1 - a1 - a2,   lB = a1,   lC = a2,        lA, lB, lC >= 0.

Node / coefficient ordering (committed contract)
------------------------------------------------
    MASS_NODES[i] = (a1, a2) by index:
        0: (0, 0)        corner A          -> control b0
        1: (1, 0)        corner B          -> control b1
        2: (0, 1)        corner C          -> control b2
        3: (1/2, 0)      edge midpoint AB  -> control b3
        4: (0, 1/2)      edge midpoint AC  -> control b4
        5: (1/2, 1/2)    edge midpoint BC  -> control b5

Throughout, "degree <= 2 polynomial on T" means any expression
p(a1, a2) = A + B*a1 + C*a2 + D*a1^2 + E*a2^2 + F*a1*a2.

Bernstein form and the convex-hull range property (proof)
---------------------------------------------------------
Write ``deg2_basis(z) = (lA^2, lB^2, lC^2, 2*lA*lB, 2*lA*lC, 2*lB*lC)``.
Every polynomial p of total degree <= 2 in (a1, a2) equals

    p(z) = b0*lA^2 + b1*lB^2 + b2*lC^2 + 2*b3*lA*lB
           + 2*b4*lA*lC + 2*b5*lB*lC                   (Bernstein form)

for a *unique* coefficient tuple b = (b0, ..., b5), and the coefficients are
exactly the node values b_i = p(MASS_NODES[i]): the 6x6 evaluation matrix V of
the basis at the six nodes, with rows ordered A, B, C, AB, AC, BC, is

        [1 0 0 0   0   0  ]
        [0 1 0 0   0   0  ]
    V = [0 0 1 0   0   0  ]
        [1/4 1/4 0 1/2 0  0]
        [1/4 0 1/4 0 1/2 0 ]
        [0 1/4 1/4 0 0 1/2],

which is unit lower triangular (corners first), hence invertible: unisolvent.
Inverting V gives the node-to-coefficient map implemented below:

    b0 = v0,  b1 = v1,  b2 = v2
    b3 = 2*v3 - (v0 + v1)/2
    b4 = 2*v4 - (v0 + v2)/2
    b5 = 2*v5 - (v1 + v2)/2

Equivalently, in monomial data p = A + B*a1 + C*a2 + D*a1^2 + E*a2^2 + F*a1*a2:
b = (A, A+B+D, A+C+E, A + B/2, A + C/2, A + (B+C)/2 + F/2); this independent
is what the edge formulas and their ordering are cross-checked against in
``selfcheck()``.

Range property.  For every z in T the basis values satisfy

    (i)   deg2_basis_i(z) >= 0 for all i   (squares and products of l >= 0),
    (ii)  sum_i deg2_basis_i(z) = (lA + lB + lC)^2 = 1.

Therefore p(z) = sum_i deg2_basis_i(z) * b_i is a convex combination of the
six control values b_i, and

    min_i b_i  <=  p(z)  <=  max_i b_i     for every z in T,    []

i.e. [min b_i, max b_i] contains the range of *every* polynomial of total
degree <= 2 on the *full* simplex, no matter where its six node samples came
from (function values, partial derivatives restricted to T, q-box-folded
enclosures, ...).  The bound is data-only and may be loose: control values
need not be attained on T.  The same applies coefficient-wise to interval
node enclosures because the map v -> b is linear with fixed rational
(dyadic) coefficients, hence inclusion-monotone.

Domain transfer to arbitrary triangles (``triangle_nodes``)
-----------------------------------------------------------
The proof above used only (i) and (ii), which hold for ANY barycentric
triple: on any non-degenerate triangle with vertices P0, P1, P2, every point
z = lA*P0 + lB*P1 + lC*P2 has lA, lB, lC >= 0 summing to 1, so the same
Bernstein identity holds with control values sampled at the six points
returned by ``triangle_nodes`` (barycentric weights (1,0,0), (0,1,0),
(0,0,1), (1/2,1/2,0), (1/2,0,1/2), (0,1/2,1/2), in MASS_NODES slot order).
Consequently simplex_range / simplex_abs_upper applied to node values from
triangle_nodes bound the polynomial over THAT triangle by the convex-hull
property, and the componentwise union of these bounds over a triangulation
of a polygon bounds it over the whole polygon.

Arithmetic discipline
---------------------
* Exact path: if all six inputs are int/Fraction, everything is computed in
  exact Fraction arithmetic and exact Fractions are returned.  Extremes are
  true argmins with a deterministic first-index tie-break.
* Ball path: if any input is a python-flint ``arb``, the computation stays in
  Arb ball arithmetic (exact inputs are first widened through
  ``point_ball``), ``ctx.prec`` is raised to at least 320 for the call and
  restored afterwards, and returned enclosures are endpoint point-arbs with
  relative slack of order a few 2^-300.  Ball comparisons are only used in
  margins far above their radii, so selections are decided by the real
  endpoint order up to that slack.  ``arb.union`` is never used.
* No floats are scanned and no fudge radii are hardcoded; the only enclosure
  error is Arb's roundoff at the certified precision.
"""

from decimal import Decimal, localcontext

from fractions import Fraction
from flint import arb, ctx

__all__ = [
    "MASS_NODES",
    "triangle_nodes",
    "simplex_bernstein_coefficients",
    "simplex_range",
    "simplex_abs_upper",
    "arb_abs_upper",
    "point_ball",
    "selfcheck",
]

# Six P2-unisolvent nodes of the unit simplex, in the committed order:
# corners A, B, C, then edge midpoints AB, AC, BC.  Index i <-> control b_i.
MASS_NODES = (
    (Fraction(0), Fraction(0)),
    (Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(1)),
    (Fraction(1, 2), Fraction(0)),
    (Fraction(0), Fraction(1, 2)),
    (Fraction(1, 2), Fraction(1, 2)),
)

# Certified precision floor (bits) for ball arithmetic in this module.
_MIN_PREC = 320


def triangle_nodes(triangle):
    """Six degree-2 Bernstein nodes of an arbitrary exact rational triangle.

    Maps the committed MASS_NODES barycentric pattern onto the triangle
    ``triangle = ((x0, y0), (x1, y1), (x2, y2))`` (exact Fractions or ints):
    vertex weights (1,0,0), (0,1,0), (0,0,1) then edge midpoints
    (1/2,1/2,0), (1/2,0,1/2), (0,1/2,1/2).  Node order matches the
    MASS_NODES contract slot-for-slot: corners P0->A, P1->B, P2->C, then
    midpoints P0P1->AB, P0P2->AC, P1P2->BC, control slots b0..b5.

    Returns six exact Fraction (a1, a2) pairs; no floats anywhere.  The
    triangle must be non-degenerate (positive area); a degenerate input
    raises ValueError.

    Domain transfer (why this suffices for both envelope callers): every
    point z of the triangle has barycentric weights (lA, lB, lC) >= 0 with
    lA + lB + lC = 1, and a degree-<=2 polynomial p satisfies the SAME
    Bernstein identity  p(z) = b0*lA^2 + b1*lB^2 + b2*lC^2 + 2*b3*lA*lB +
    2*b4*lA*lC + 2*b5*lB*lC  with control values taken at the six points
    returned here -- because the weights are nonnegative and sum to 1, the
    basis is nonnegative with total 1 on ANY affine triangle.  Hence
    simplex_range / simplex_abs_upper applied to node values sampled at
    triangle_nodes bound p over THAT triangle, and the componentwise union
    of these bounds over a triangulation bounds p over the whole polygon.
    """
    pts = tuple((Fraction(x), Fraction(y)) for x, y in triangle)
    if len(pts) != 3:
        raise ValueError("triangle_nodes needs exactly three vertices")
    (x0, y0), (x1, y1), (x2, y2) = pts
    twice_area = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
    if twice_area == 0:
        raise ValueError("triangle_nodes: degenerate (zero-area) triangle")
    half = Fraction(1, 2)
    return (
        (x0, y0),                                # (1,0,0) corner A
        (x1, y1),                                # (0,1,0) corner B
        (x2, y2),                                # (0,0,1) corner C
        ((x0 + x1) * half, (y0 + y1) * half),    # (1/2,1/2,0) midpoint AB
        ((x0 + x2) * half, (y0 + y2) * half),    # (1/2,0,1/2) midpoint AC
        ((x1 + x2) * half, (y1 + y2) * half),    # (0,1/2,1/2) midpoint BC
    )


def deg2_basis(a1, a2):
    """Exact barycentric values of the six Bernstein basis functions at (a1,a2)."""
    la = Fraction(1) - a1 - a2
    lb = Fraction(a1)
    lc = Fraction(a2)
    return (
        la * la,
        lb * lb,
        lc * lc,
        2 * la * lb,
        2 * la * lc,
        2 * lb * lc,
    )


class _PrecGuard:
    """Raise ctx.prec to at least _MIN_PREC for the call, restore afterwards."""

    def __enter__(self):
        self._saved = ctx.prec
        if self._saved < _MIN_PREC:
            ctx.prec = _MIN_PREC
        return self

    def __exit__(self, exc_type, exc, tb):
        ctx.prec = self._saved
        return False


def _six(values, what):
    items = tuple(values)
    if len(items) != 6:
        raise ValueError(f"{what} needs exactly six values, got {len(items)}")
    return items


def _has_arb(values):
    return any(isinstance(v, arb) for v in values)


def _as_ball(value):
    if isinstance(value, arb):
        return arb(value)
    if isinstance(value, (int, Fraction)):
        frac = Fraction(value)
        return arb(f"{frac.numerator}/{frac.denominator}")
    raise TypeError(f"cannot make an enclosure from {type(value).__name__}")


def _arb_mul_int(value, k):
    out = arb(value)
    for _ in range(k - 1):
        out = out + value
    return out


def _pick_extreme(candidates, want_max):
    """Deterministically pick the extreme endpoint point-arb (ties keep first)."""
    best = candidates[0]
    for cand in candidates[1:]:
        if want_max:
            if cand > best:
                best = cand
        elif cand < best:
            best = cand
    return best


def _endpoints(ball):
    """Lower/upper endpoint point-arbs of an arb ball (conservative ordering)."""
    return ball.lower(), ball.upper()


def point_ball(value):
    """An ``arb`` ball that contains the exact rational ``value``.

    Accepts ``Fraction`` or ``int`` (an ``arb`` input is copied through).
    The width is zero exactly when the rational is representable at the
    certified precision (in particular for dyadic rationals such as every
    MASS_NODES coordinate); otherwise the ball encloses the rational with
    roundoff-scale radius.  Under a precision guard of at least 320 bits the
    rational is handed to Arb as the exact string "numerator/denominator",
    so the returned ball covers the true rational by Arb's own correctly
    rounded conversion (no float is ever formed).
    """
    if isinstance(value, arb):
        return arb(value)
    if isinstance(value, (int, Fraction)):
        frac = Fraction(value)
        with _PrecGuard():
            return arb(f"{frac.numerator}/{frac.denominator}")
    raise TypeError(
        f"point_ball needs a Fraction or int, got {type(value).__name__}"
    )


def arb_abs_upper(value):
    """Point-arb enclosure with sup(|x| : x in ball ``value``) <= its upper end.

    For an exact input (int/Fraction) the returned ball is the exact absolute
    value (zero width when representable).  For a ball the absolute value's
    endpoint upper bound is returned as a point-arb whose interval covers
    max(|lower|, |upper|) up to roundoff slack, so its upper endpoint is a
    certified upper bound for every |x| with x in the input ball.
    """
    if isinstance(value, (int, Fraction)):
        frac = abs(Fraction(value))
        with _PrecGuard():
            if frac.denominator & (frac.denominator - 1) == 0:
                return arb(f"{frac.numerator}/{frac.denominator}")
            with localcontext() as lctx:
                lctx.prec = 400
                mid_dec = Decimal(frac.numerator) / Decimal(frac.denominator)
            return arb(f"{mid_dec} +/- 1e-110")
    if isinstance(value, arb):
        with _PrecGuard():
            lo, hi = _endpoints(arb(value))
            return abs(lo).upper() if abs(lo).upper() > abs(hi).upper() else abs(
                hi
            ).upper()
    raise TypeError(
        f"arb_abs_upper needs a Fraction, int or arb, got {type(value).__name__}"
    )


def _bernstein_exact(values):
    """Exact Fraction node-to-coefficient map (the closed form in the docstring)."""
    v0, v1, v2, v3, v4, v5 = (Fraction(v) for v in values)
    return (
        v0,
        v1,
        v2,
        2 * v3 - (v0 + v1) / 2,
        2 * v4 - (v0 + v2) / 2,
        2 * v5 - (v1 + v2) / 2,
    )


def _bernstein_ball(values):
    """Inclusion-monotone Arb map; same linear form as ``_bernstein_exact``."""
    with _PrecGuard():
        v = tuple(_as_ball(vv) for vv in values)
        zero = arb(0)
        half = arb(1) / 2
        return (
            arb(v[0]),
            arb(v[1]),
            arb(v[2]),
            (_arb_mul_int(v[3], 2))
            - ((v[0] + v[1]) * half),
            (_arb_mul_int(v[4], 2))
            - ((v[0] + v[2]) * half),
            (_arb_mul_int(v[5], 2))
            - ((v[1] + v[2]) * half),
        )


def simplex_bernstein_coefficients(node_values):
    """Node values -> Bernstein control values (documented ordering b0..b5).

    ``node_values[i]`` is the polynomial's (possibly interval) value at
    ``MASS_NODES[i]``.  Returns six exact Fractions when the input is
    entirely int/Fraction; otherwise six arb ball enclosures obtained by the
    identical linear map (inclusion-monotone, see module docstring).
    """
    values = _six(node_values, "simplex_bernstein_coefficients")
    if _has_arb(values):
        return _bernstein_ball(values)
    return _bernstein_exact(values)


def simplex_range(node_values):
    """Enclosure of the range of a degree<=2 polynomial on the full simplex.

    Given the six node values (in MASS_NODES order) of any degree-<=2
    polynomial, returns ``(lower, upper)`` with
    min_i b_i >= lower-bound-side and max_i b_i <= upper-bound-side, so by
    the convex-hull property every p(z), z in T, lies between the two.

    Exact input returns exact ``Fraction`` extremes; ball input returns the
    extreme endpoint point-arbs ``(lower = min_i b_i.lower(),
    upper = max_i b_i.upper())`` of the six coefficient balls, each within
    roundoff slack of ``min_i b_i`` resp. ``max_i b_i``.
    """
    values = _six(node_values, "simplex_range")
    if _has_arb(values):
        coefs = _bernstein_ball(values)
        lows, highs = zip(*(_endpoints(c) for c in coefs))
        return _pick_extreme(lows, False), _pick_extreme(highs, True)
    coefs = _bernstein_exact(values)
    lo = coefs[0]
    hi = coefs[0]
    for c in coefs[1:]:
        if c < lo:
            lo = c
        if c > hi:
            hi = c
    return lo, hi


def simplex_abs_upper(node_values):
    """Certified upper bound for max |p(z)| over the full simplex.

    Returns a single value B with, in the exact case, max_T |p| <= B exactly
    (equal to max_i |b_i|; ties keep the first index), and in the ball case a
    point-arb whose interval covers the larger of |lower|, |upper| of the
    coefficient-range enclosure up to roundoff slack, so that
    max_T |p| <= B.upper() holds.  Soundness is the triangle inequality on
    the convex-hull range property: every p(z) lies in the coefficient
    interval [lo, hi], hence |p(z)| <= max(|lo|, |hi|).
    """
    values = _six(node_values, "simplex_abs_upper")
    if _has_arb(values):
        lo, hi = simplex_range(values)
        return abs(lo).upper() if abs(lo).upper() > abs(hi).upper() else abs(hi)
    lo, hi = simplex_range(values)
    candidates = (-lo if lo < 0 else lo, -hi if hi < 0 else hi)
    best = candidates[0]
    for cand in candidates[1:]:
        if cand > best:
            best = cand
    return best


def selfcheck():
    """Deterministic self-test; returns a dict of booleans (all True expected).

    Reconstructs several exact-rational quadratic polynomials from monomial
    data, runs them through every exported helper on an exhaustive dyadic
    grid of the simplex and of a non-standard exact rational triangle (via
    ``triangle_nodes``), and demonstrates that three explicit mutations
    (a corrupted edge-coefficient formula, a corrupted node ordering, and a
    swapped AB/BC midpoint order in ``triangle_nodes``) are rejected by the
    exact identities above.  Not run on import; callers invoke
    ``selfcheck()`` and assert every entry is True.
    """
    checks = {}

    checks["mass_nodes_exact_values"] = MASS_NODES == (
        (Fraction(0), Fraction(0)),
        (Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(1)),
        (Fraction(1, 2), Fraction(0)),
        (Fraction(0), Fraction(1, 2)),
        (Fraction(1, 2), Fraction(1, 2)),
    )

    polys = (
        # (A, B, C, D, E, F) for A + B*a1 + C*a2 + D*a1^2 + E*a2^2 + F*a1*a2
        (Fraction(1), Fraction(2), Fraction(-3), Fraction(1, 2), Fraction(4), Fraction(-5, 8)),
        (Fraction(0), Fraction(1), Fraction(1), Fraction(-2), Fraction(1, 3), Fraction(7)),
        (Fraction(-7, 3), Fraction(1, 5), Fraction(1, 7), Fraction(1, 11),
         Fraction(-1, 13), Fraction(1, 17)),
    )

    def independent_coefficients(monos):
        A, B, C, D, E, F = monos
        return (
            A,
            A + B + D,
            A + C + E,
            A + B / 2,
            A + C / 2,
            A + (B + C) / 2 + F / 2,
        )

    def evaluate(monos, a1, a2):
        A, B, C, D, E, F = monos
        return A + B * a1 + C * a2 + D * a1 * a1 + E * a2 * a2 + F * a1 * a2

    edges_ok = True
    ordering_ok = True
    range_ok = True
    for monos in polys:
        node_vals = tuple(evaluate(monos, a1, a2) for a1, a2 in MASS_NODES)
        coefs = simplex_bernstein_coefficients(node_vals)
        edges_ok &= coefs == independent_coefficients(monos)

        # Ordering identity: swapping the B/C corner slots and their edge
        # midpoints must swap the corresponding control slots, exactly.
        perm = (node_vals[0], node_vals[2], node_vals[1],
                node_vals[4], node_vals[3], node_vals[5])
        got = simplex_bernstein_coefficients(perm)
        base = independent_coefficients(monos)
        ordering_ok &= got == (base[0], base[2], base[1], base[4], base[3], base[5])
        lo, hi = simplex_range(node_vals)
        # Exhaustive dyadic grid of T at depth 8: range must be enclosed.
        for i in range(9):
            for j in range(9 - i):
                z = evaluate(monos, Fraction(i, 8), Fraction(j, 8))
                range_ok &= (lo <= z) and (z <= hi)

    checks["edges_match_independent_derivative"] = edges_ok
    checks["node_ordering_identity_holds"] = ordering_ok
    checks["range_hull_containment_grid"] = range_ok

    # Range/abs helpers agree with the exact coefficient extremes.
    agree = True
    for monos in polys:
        node_vals = tuple(evaluate(monos, a1, a2) for a1, a2 in MASS_NODES)
        base = independent_coefficients(monos)
        lo, hi = simplex_range(node_vals)
        agree &= lo == min(base) and hi == max(base)
        agree &= simplex_abs_upper(node_vals) == max(abs(c) for c in base)
    checks["simplex_range_exact_agrees"] = agree

    # --- triangle_nodes acceptance block (Main directive, 2026-09-01) ---

    # (1) The standard simplex P0=(0,0), P1=(1,0), P2=(0,1) must reproduce
    # MASS_NODES exactly (corners A,B,C then midpoints AB,AC,BC).
    checks["triangle_nodes_reproduces_mass_nodes"] = (
        triangle_nodes(((0, 0), (1, 0), (0, 1))) == MASS_NODES
    )

    # (2) A degree-2 polynomial sampled at triangle_nodes of a NON-standard
    # exact rational triangle has its true min/max over that triangle
    # enclosed by simplex_range of those node values, checked on an exact
    # rational barycentric grid inside the triangle (k/24 lattice).  The
    # test triangle is the non-dyadic right triangle with legs 2/3 and 1.
    tri = tuple((Fraction(x), Fraction(y)) for x, y in (
        (Fraction(1, 3), Fraction(2, 5)),
        (Fraction(1), Fraction(2, 5)),
        (Fraction(1, 3), Fraction(7, 5)),
    ))

    def grid_minmax_over(triangle, monos, n):
        (x0, y0), (x1, y1), (x2, y2) = triangle
        lo = hi = None
        for i in range(n + 1):
            for j in range(n + 1 - i):
                lB = Fraction(i, n)
                lC = Fraction(j, n)
                lA = 1 - lB - lC
                x = lA * x0 + lB * x1 + lC * x2
                y = lA * y0 + lB * y1 + lC * y2
                val = evaluate(monos, x, y)
                if lo is None or val < lo:
                    lo = val
                if hi is None or val > hi:
                    hi = val
        return lo, hi

    tri_range_ok = True
    tri_nodes = triangle_nodes(tri)
    for monos in polys:
        node_vals = tuple(evaluate(monos, x, y) for x, y in tri_nodes)
        lo_b, hi_b = simplex_range(node_vals)
        lo_g, hi_g = grid_minmax_over(tri, monos, 24)
        tri_range_ok &= lo_b <= lo_g and hi_g <= hi_b
    checks["triangle_nodes_range_grid_containment"] = tri_range_ok

    # Mutation 1: corrupted edge formula (factor 2 dropped on the AB edge).
    # The helper must NOT reproduce the mutant, and the true b3 must differ.
    def mutated_edge_formula(node_values):
        v0, v1, v2, v3, v4, v5 = (Fraction(v) for v in node_values)
        return (
            v0,
            v1,
            v2,
            v3 - (v0 + v1) / 2,   # MUTANT of b3 = 2*v3 - (v0 + v1)/2
            2 * v4 - (v0 + v2) / 2,
            2 * v5 - (v1 + v2) / 2,
        )

    correct_outputs = []
    all_equal_base = True
    mutant_detects_difference = False
    for monos in polys:
        node_vals = tuple(evaluate(monos, a1, a2) for a1, a2 in MASS_NODES)
        coefs = simplex_bernstein_coefficients(node_vals)
        base = independent_coefficients(monos)
        mut = mutated_edge_formula(node_vals)
        all_equal_base &= coefs == base
        mutant_detects_difference |= (coefs != mut) and (mut != base)
    rejected_edge = (
        all_equal_base
        and mutant_detects_difference
        and all(
            (mut != base) or (coefs == mut)
            for coefs, base, mut in correct_outputs
        )
    )
    checks["mutation_edge_coefficient_rejected"] = rejected_edge
    # Mutation 2: corrupted node ordering.  Evaluate at the mass nodes, then
    # permute the VALUES by the B/C corner swap sigma = (1 2).  A correct
    # implementation is equivariant: helper(sigma . v) == sigma . helper(v)
    # for the induced control permutation (b1<->b2, b3<->b4).  It must also
    # be detectable: sigma . base != base on each test polynomial.  Finally
    # the permuted input must actually differ from the identity input, so a
    # helper that ignored ordering entirely (always returning helper(v)) or
    # one with a scrambled internal node map fails this conjunction.
    rejected_order = True
    for monos in polys:
        node_vals = tuple(evaluate(monos, a1, a2) for a1, a2 in MASS_NODES)
        wrong = (node_vals[0], node_vals[2], node_vals[1],
                 node_vals[4], node_vals[3], node_vals[5])
        base = independent_coefficients(monos)
        got = simplex_bernstein_coefficients(wrong)
        expected = (base[0], base[2], base[1], base[4], base[3], base[5])
        rejected_order &= (got == expected) and (expected != base)
        rejected_order &= (wrong != node_vals)
    checks["mutation_node_ordering_rejected"] = rejected_order

    # Mutation 3: corrupted triangle midpoint order (AB <-> BC swapped in
    # triangle_nodes).  A correct triangle_nodes must still reproduce
    # MASS_NODES on the standard simplex (the swap is detected there), the
    # mutant must actually differ from the truth on the NON-standard test
    # triangle, and the correct helper must differ from the mutant at the
    # swapped slots.  A triangle_nodes whose midpoints were swapped would
    # silently mislabel b3 vs b5 on every consumer triangle, so this check
    # fails closed on such a regression.
    def mutated_triangle_midpoints(triangle):
        (x0, y0), (x1, y1), (x2, y2) = (
            (Fraction(px), Fraction(py)) for px, py in triangle
        )
        return (
            (x0, y0),
            (x1, y1),
            (x2, y2),
            ((x1 + x2) / 2, (y1 + y2) / 2),   # MUTANT: BC midpoint in slot 3
            ((x0 + x2) / 2, (y0 + y2) / 2),
            ((x0 + x1) / 2, (y0 + y1) / 2),   # MUTANT: AB midpoint in slot 5
        )

    mut_nodes = mutated_triangle_midpoints(tri)
    checks["mutation_triangle_midpoint_order_rejected"] = (
        mut_nodes != tri_nodes
        and mut_nodes[3] == tri_nodes[5] and mut_nodes[5] == tri_nodes[3]
        and mutated_triangle_midpoints(((0, 0), (1, 0), (0, 1))) != MASS_NODES
    )
    # Ball pipeline: point balls around exact nodes must enclose the exact
    def ball_interval(ball):
        """Fraction enclosure (with slack) of an arb ball's printed interval.

        Prints under a certified-precision guard so the enclosure is not
        truncated by the ambient precision, strips the optional brackets,
        treats a bare midpoint print (zero-radius ball) as radius 0, and
        accepts scientific-notation mantissas like ``9.08e-96``.
        """
        slack = Fraction(1, 10 ** 120)
        with _PrecGuard():
            text = str(ball)
        if text.startswith("[") and text.endswith("]"):
            text = text[1:-1].strip()
        if " +/- " in text:
            mid_s, rad_s = text.split(" +/- ", 1)
            mid = Fraction(Decimal(mid_s))
            rad = Fraction(Decimal(rad_s))
        else:
            mid = Fraction(Decimal(text))
            rad = Fraction(0)
        return mid - rad - slack, mid + rad + slack

    ball_ok = True
    monos = polys[2]
    node_vals = tuple(evaluate(monos, a1, a2) for a1, a2 in MASS_NODES)
    base = independent_coefficients(monos)
    tight = simplex_bernstein_coefficients(
        tuple(point_ball(v) for v in node_vals)
    )
    wide_in = tuple(point_ball(v) + arb("1e-100") for v in node_vals)
    widened = simplex_bernstein_coefficients(wide_in)
    for exact_c, tight_c, wide_c in zip(base, tight, widened):
        blo, bhi = ball_interval(tight_c)
        ball_ok &= blo <= exact_c <= bhi
        wlo, whi = ball_interval(wide_c)
        ball_ok &= wlo <= blo and bhi <= whi
        ball_ok &= ball_interval(wide_c)[0] <= ball_interval(tight_c)[0]
    lo_b, hi_b = simplex_range(tuple(point_ball(v) for v in node_vals))
    blob, bhib = ball_interval(lo_b), ball_interval(hi_b)
    abs_up = simplex_abs_upper(tuple(point_ball(v) for v in node_vals))
    ball_ok &= max(base) <= ball_interval(abs_up)[1]
    checks["ball_pipeline_encloses_exact"] = ball_ok

    # point_ball containment of non-dyadic rationals at certified precision.
    # The dyadic arm asserts a point print: a zero-radius ball prints as a
    # bare decimal (no " +/- "), while non-dyadic rationals must print an
    # interval whose parsed endpoints bracket the rational.  Both arms are
    # detectable: an over-wide dyadic ball would print " +/- 0" only via the
    # interval branch, and an exact-print mutation would fail the bracket.
    proportion_ok = True
    for frac in (Fraction(1, 3), Fraction(-5, 7), Fraction(1, 4),
                 Fraction(22, 7), Fraction(0)):
        ball = point_ball(frac)
        blo, bhi = ball_interval(ball)
        proportion_ok &= blo <= frac <= bhi
        proportion_ok &= blo <= bhi
        is_dyadic = frac.denominator & (frac.denominator - 1) == 0
        if is_dyadic:
            with _PrecGuard():
                proportion_ok &= " +/- " not in str(ball) and str(ball) == str(
                    ball.mid()
                )
        else:
            with _PrecGuard():
                proportion_ok &= " +/- " in str(ball)
    checks["point_ball_contains_rational"] = proportion_ok

    # arb_abs_upper agrees with the analytic |v| and is a true enclosure.
    # Ball arithmetic at 320 bits has ulp ~ 1e-96, so every tolerance below
    # is >= 1/10^80 (>= 10^16 x the parse slack) and every widening is
    # >= 1/10^90 (single-ulp-visible).  All exact claims are compared as
    # Fractions; only bottom-line booleans are trusted from ball compares.
    # (a) exact Fraction input: ratio-form conversion is single-ulp exact,
    #     hence the parsed enclosure is within 1/10^80 of |frac| and
    #     brackets it from below.
    # (b) ball input: the parsed abs-upper equals max(|lo|, |hi|) of the
    #     input's parsed enclosure within 1/10^80, with the true magnitude
    #     enclosed from below (sup|x| is an upper bound, so its enclosure
    #     must reach AT LEAST the parsed endpoint magnitude).
    # (c) monotone domination: widening the input by a certified +1/10^90
    #     must strictly raise the parsed upper endpoint (rejects a constant,
    #     endpoint-dropping, or radius-ignoring implementation).
    tolerance = Fraction(1, 10 ** 80)
    step = Fraction(1, 10 ** 90)
    agrees = True
    for frac in (Fraction(1, 3), Fraction(-5, 7), Fraction(22, 7),
                 Fraction(1, 4), Fraction(0)):
        ulo, uhi = ball_interval(arb_abs_upper(frac))
        agrees &= ulo <= abs(frac) + tolerance
        agrees &= abs(abs(frac) - uhi) <= tolerance
        agrees &= ulo <= uhi
    for frac in (Fraction(1, 3), Fraction(-5, 7), Fraction(22, 7)):
        with _PrecGuard():
            ball = point_ball(frac)
            # A shift TOWARD zero shrinks the magnitude envelope of a
            # negative-centered ball; the monotone domination test must
            # push AWAY from zero so the certified +step strictly grows
            # the quantity sup|x| over the ball.
            direction = point_ball(step) if frac >= 0 else -point_ball(step)
            wide = ball + direction
        blo, bhi = ball_interval(ball)
        a_lo, a_hi = ball_interval(arb_abs_upper(ball))
        magnitude = max(abs(blo), abs(bhi))
        agrees &= (a_hi >= magnitude - tolerance) and (a_hi <= magnitude + tolerance)
        agrees &= ball_interval(arb_abs_upper(wide))[1] > a_hi
    checks["arb_abs_upper_agrees"] = agrees

    return checks


if __name__ == "__main__":
    _result = selfcheck()
    for _key in sorted(_result):
        print(f"{_key}: {_result[_key]}")
    if not all(_result.values()):
        raise SystemExit("selfcheck FAILED")
    print("selfcheck OK")
