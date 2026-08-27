"""Exact local bridge from the certified pair-orbit functional to Cambie's Q2.

Scope and evidence level
------------------------
HUMAN PROOFS WITH EXACT FINITE IMPLEMENTATION CONTROLS:

1. For ``p,r in [0,1]``, the repository's clipped-median formula

       min(max(1/2,p,r), min(p+r,1))

   equals Cambie's Question-2 formula

       max(p,r,min(p+r,1/2)).

   Put ``m=max(p,r)`` and ``s=p+r``.  The following three exhaustive cases
   prove the identity without numerical approximation:

   * ``s <= 1/2``: ``m <= s <= 1/2``, and both formulas equal ``s``;
   * ``s >= 1/2`` and ``m <= 1/2``: both formulas equal ``1/2``;
   * ``m >= 1/2``: ``s >= m`` and ``1 >= m``, and both equal ``m``.

2. If a coupling ``pi`` has equal first and second marginals ``mu``, then
   ``(pi + transpose(pi))/2`` has both marginals ``mu`` and preserves every
   symmetric cost.  Pushing this symmetric coupling through
   ``(p,r) -> (min(p,r),max(p,r))`` produces a probability measure ``nu`` on
   ``Delta={p<=r}``; its induced marginal is ``mu``, and its cost is unchanged.
   This is an identity of integrals, so it applies to arbitrary Borel
   probability measures, not only finite support.  The symbolic finite-matrix
   schema below checks the coefficient identities behind this argument.

3. Consequently the exact functional minimized in ``thmB3_proof.py`` is the
   left side minus the right side of Cambie's Question 2:

       (1-alpha) E h(p+q-pq) + alpha E h(s*(p,r)) - E h(p),

   where ``p,q`` are iid with law ``mu`` and ``(p,r)`` has any equal-marginal
   coupling.  The support reduction, Margin Lemma, and interval certificate are
   separate pinned proof artifacts.

CITED ORIGIN, RECONSTRUCTED LOCALLY: Cambie v2, Question 2 / equation (1)
states the hypothesis for expectations **strictly below** ``c``. Section 4's
wording switches to “at most” while deriving the union-closed consequence.
The closed-domain certificate is stronger than Question 2's domain, and the
sequential Bernoulli-coupling induction in ``AUDIT.md`` and ``paper/main.tex``
resolves the boundary by proving a strict positive first-coordinate margin.
The exact exhaustive control below constructs that coupling for every nonempty
set family on at most three coordinates and checks all prefix and final
marginals.  This guards the finite implementation; the displayed induction is
the proof for arbitrary finite families.

Primary source (fixed version, retrievable from arXiv):
    Stijn Cambie, "Better bounds for the union-closed sets conjecture using the
    entropy approach", arXiv:2212.12500v2, Question 2 and Section 4,
    https://arxiv.org/abs/2212.12500v2

No sampled check below is promoted to a universal claim.  The universal
clipped-median and measure identities are the displayed algebraic proofs; the
exact controls guard their implementation and sign conventions.
"""

from fractions import Fraction

import sympy as sp


SOURCE = "https://arxiv.org/abs/2212.12500v2"
HALF = Fraction(1, 2)
ONE = Fraction(1)


def cambie_sstar(p, r):
    """Cambie Q2's ``max(p,r,min(p+r,1/2))`` on exact rationals."""
    return max(p, r, min(p + r, HALF))


def repository_sstar(p, r):
    """The repository's ``min(max(1/2,p,r),min(p+r,1))``."""
    return min(max(HALF, p, r), min(p + r, ONE))


def proved_case_value(p, r):
    """Apply the three exhaustive cases proved in the module docstring."""
    assert 0 <= p <= 1 and 0 <= r <= 1
    m = max(p, r)
    s = p + r
    if s <= HALF:
        assert m <= s <= HALF
        value = s
        case = "sum<=half"
    elif m <= HALF:
        assert s >= HALF and m <= HALF
        value = HALF
        case = "straddles-half"
    else:
        assert m >= HALF and s >= m and ONE >= m
        value = m
        case = "max>=half"
    assert cambie_sstar(p, r) == value
    assert repository_sstar(p, r) == value
    return case


def exact_sstar_control(denominator=64):
    """Exact rational implementation control; the proof is the case split."""
    counts = {"sum<=half": 0, "straddles-half": 0, "max>=half": 0}
    for i in range(denominator + 1):
        for j in range(denominator + 1):
            p = Fraction(i, denominator)
            r = Fraction(j, denominator)
            counts[proved_case_value(p, r)] += 1
    assert all(counts.values())
    return counts


def bernoulli_union_identities():
    """Exact symbolic OR identities for independent and coupled Bernoullis."""
    p, q, r, z = sp.symbols("p q r z")

    def bernoulli_prob(parameter, bit):
        return parameter if bit else 1 - parameter

    independent_or = sum(
        int(bool(a or b)) * bernoulli_prob(p, a) * bernoulli_prob(q, b)
        for a in (0, 1) for b in (0, 1)
    )
    assert sp.expand(independent_or - (p + q - p * q)) == 0

    # A general Bernoulli(p), Bernoulli(r) coupling is determined by
    # z=P(1,1): P(1,0)=p-z, P(0,1)=r-z, P(0,0)=1-p-r+z.
    coupled_or = (p - z) + (r - z) + z
    assert sp.expand(coupled_or - (p + r - z)) == 0
    assert sp.expand(z + (p - z) + (r - z) + (1 - p - r + z) - 1) == 0
    return True


def _prefix_parameter(family, prefix, coordinate):
    """Exact next-bit probability under the uniform law on ``family``."""
    mask = (1 << coordinate) - 1
    eligible = tuple(member for member in family if member & mask == prefix)
    assert eligible
    bit = 1 << coordinate
    return Fraction(sum(bool(member & bit) for member in eligible), len(eligible))


def _add_mass(mapping, key, mass):
    mapping[key] = mapping.get(key, Fraction(0)) + mass


def sequential_coupling_control(max_coordinates=3):
    """Exhaust the sequential coupling on every small nonempty set family.

    This is an exact finite implementation control, not the universal proof.
    The universal statement follows from the same one-step marginal calculation
    and induction written in the manuscript.
    """
    summaries = {}
    for coordinate_count in range(1, max_coordinates + 1):
        universe_size = 1 << coordinate_count
        family_count = 0
        prefix_state_count = 0
        transition_count = 0
        for family_mask in range(1, 1 << universe_size):
            family = tuple(
                member for member in range(universe_size)
                if family_mask & (1 << member)
            )
            family_count += 1
            uniform_mass = Fraction(1, len(family))
            joint = {(0, 0): ONE}

            for coordinate in range(coordinate_count):
                next_joint = {}
                law_p = {}
                law_r = {}
                frequency = Fraction(
                    sum(bool(member & (1 << coordinate)) for member in family),
                    len(family),
                )
                for (a_prefix, c_prefix), state_mass in joint.items():
                    prefix_state_count += 1
                    p = _prefix_parameter(family, a_prefix, coordinate)
                    r = _prefix_parameter(family, c_prefix, coordinate)
                    _add_mass(law_p, p, state_mass)
                    _add_mass(law_r, r, state_mass)
                    s = cambie_sstar(p, r)
                    probabilities = {
                        (1, 1): p + r - s,
                        (1, 0): s - r,
                        (0, 1): s - p,
                        (0, 0): ONE - s,
                    }
                    assert sum(probabilities.values(), Fraction(0)) == ONE
                    assert all(probability >= 0 for probability in probabilities.values())
                    assert sum(
                        probability for (a_bit, _), probability in probabilities.items()
                        if a_bit
                    ) == p
                    assert sum(
                        probability for (_, c_bit), probability in probabilities.items()
                        if c_bit
                    ) == r
                    assert sum(
                        probability for bits, probability in probabilities.items()
                        if any(bits)
                    ) == s
                    for (a_bit, c_bit), probability in probabilities.items():
                        if not probability:
                            continue
                        transition_count += 1
                        next_a = a_prefix | (a_bit << coordinate)
                        next_c = c_prefix | (c_bit << coordinate)
                        _add_mass(
                            next_joint,
                            (next_a, next_c),
                            state_mass * probability,
                        )

                assert law_p == law_r
                assert sum(value * mass for value, mass in law_p.items()) == frequency
                assert sum(value * mass for value, mass in law_r.items()) == frequency
                assert sum(next_joint.values(), Fraction(0)) == ONE
                joint = next_joint

                prefix_mask = (1 << (coordinate + 1)) - 1
                for prefix in range(1 << (coordinate + 1)):
                    expected = uniform_mass * sum(
                        member & prefix_mask == prefix for member in family
                    )
                    first_mass = sum(
                        mass for (a_prefix, _), mass in joint.items()
                        if a_prefix == prefix
                    )
                    second_mass = sum(
                        mass for (_, c_prefix), mass in joint.items()
                        if c_prefix == prefix
                    )
                    assert first_mass == expected
                    assert second_mass == expected

            for member in range(universe_size):
                expected = uniform_mass if member in family else Fraction(0)
                assert sum(
                    mass for (a_value, _), mass in joint.items()
                    if a_value == member
                ) == expected
                assert sum(
                    mass for (_, c_value), mass in joint.items()
                    if c_value == member
                ) == expected

        summaries[coordinate_count] = {
            "families": family_count,
            "prefix_states": prefix_state_count,
            "positive_transitions": transition_count,
        }
    return summaries


def symbolic_fold_schema(n):
    """Exact coefficient identities for symmetrization and folding."""
    entries = sp.symbols("x0:%d" % (n * n))
    coupling = sp.Matrix(n, n, entries)
    symmetric = (coupling + coupling.T) / 2

    # A generic symmetric cost matrix.
    cost_symbols = {}
    for i in range(n):
        for j in range(i, n):
            cost_symbols[i, j] = sp.Symbol("c%d_%d" % (i, j))

    def cost(i, j):
        return cost_symbols[min(i, j), max(i, j)]

    old_cost = sum(coupling[i, j] * cost(i, j)
                   for i in range(n) for j in range(n))
    sym_cost = sum(symmetric[i, j] * cost(i, j)
                   for i in range(n) for j in range(n))
    assert sp.expand(sym_cost - old_cost) == 0

    for i in range(n):
        row = sum(coupling[i, j] for j in range(n))
        column = sum(coupling[j, i] for j in range(n))
        sym_row = sum(symmetric[i, j] for j in range(n))
        assert sp.expand(sym_row - (row + column) / 2) == 0

        # Folding gives diagonal mass S_ii and off-diagonal orbit mass 2S_ij.
        # Splitting each off-diagonal orbit equally between its two endpoints
        # reproduces the symmetrized first marginal.
        induced = symmetric[i, i] + sum(
            symmetric[min(i, j), max(i, j)] for j in range(n) if j != i)
        assert sp.expand(induced - sym_row) == 0
    return True


def functional_identity():
    """Exact finite-schema control of all three Q2/local term identifications.

    The universal identity follows directly from the product law for the iid
    pair and from symmetrization/folding of the equal-marginal coupling.  This
    deliberately non-tautological control uses distinct orbit weights and all
    three branches of ``s*``; a missing half-weight, transpose, or cost term
    changes at least one asserted expression.
    """
    h_symbol = sp.Function("h")
    orbits = (
        (Fraction(1, 10), Fraction(1, 5), Fraction(1, 4)),
        (Fraction(2, 5), Fraction(4, 5), Fraction(1, 3)),
        (Fraction(3, 5), Fraction(3, 5), Fraction(5, 12)),
    )
    assert sum(weight for _, _, weight in orbits) == 1

    def sympy_q(value):
        return sp.Rational(value.numerator, value.denominator)

    def entropy(value):
        return h_symbol(sympy_q(value))

    # Local pair-orbit representation: each off-diagonal orbit contributes
    # half its mass to each endpoint of the induced marginal.
    local_marginal = {}
    for p, r, weight in orbits:
        local_marginal[p] = local_marginal.get(p, Fraction(0)) + weight / 2
        local_marginal[r] = local_marginal.get(r, Fraction(0)) + weight / 2

    # Q2 representation: unfold each orbit to a symmetric coupling, then derive
    # its first marginal rather than reusing ``local_marginal``.
    coupling = {}
    for p, r, weight in orbits:
        if p == r:
            coupling[p, r] = coupling.get((p, r), Fraction(0)) + weight
        else:
            coupling[p, r] = coupling.get((p, r), Fraction(0)) + weight / 2
            coupling[r, p] = coupling.get((r, p), Fraction(0)) + weight / 2
    q2_marginal = {}
    for (p, _r), weight in coupling.items():
        q2_marginal[p] = q2_marginal.get(p, Fraction(0)) + weight
    assert local_marginal == q2_marginal
    assert sum(local_marginal.values()) == 1
    assert sum(coupling.values()) == 1

    local_q = sum(
        sympy_q(wx * wy) * entropy(x + y - x * y)
        for x, wx in local_marginal.items()
        for y, wy in local_marginal.items()
    )
    q2_q = sum(
        sympy_q(wx * wy) * entropy(x + y - x * y)
        for x, wx in q2_marginal.items()
        for y, wy in q2_marginal.items()
    )
    local_c = sum(
        sympy_q(weight) * entropy(repository_sstar(p, r))
        for p, r, weight in orbits
    )
    q2_c = sum(
        sympy_q(weight) * entropy(cambie_sstar(p, r))
        for (p, r), weight in coupling.items()
    )
    local_l = sum(sympy_q(weight) * entropy(p)
                  for p, weight in local_marginal.items())
    q2_l = sum(sympy_q(weight) * entropy(p)
               for p, weight in q2_marginal.items())
    assert sp.expand(local_q - q2_q) == 0
    assert sp.expand(local_c - q2_c) == 0
    assert sp.expand(local_l - q2_l) == 0

    alpha = sp.Symbol("alpha")
    local_phi_exact = (1 - alpha) * local_q + alpha * local_c - local_l
    cambie_q2_gap = (1 - alpha) * q2_q + alpha * q2_c - q2_l
    assert sp.expand(local_phi_exact - cambie_q2_gap) == 0
    return True


def main():
    counts = exact_sstar_control()
    assert bernoulli_union_identities()
    for n in range(1, 6):
        assert symbolic_fold_schema(n)
    assert functional_identity()
    sequential = sequential_coupling_control()
    print("HUMAN PROOF BRIDGE: Cambie Q2 s* equals repository s* by 3 exhaustive cases")
    print("EXACT IMPLEMENTATION CONTROL: %d rational pairs; cases=%s" %
          (sum(counts.values()), counts))
    print("HUMAN PROOF BRIDGE: symmetrize+fold preserves equal marginal and symmetric cost")
    print("EXACT FINITE SCHEMA: Q, coupling, and marginal terms match Cambie Q2")
    print("EXACT SEQUENTIAL COUPLING CONTROL: %s" % sequential)
    print("HUMAN-AUDITED UC IMPLICATION: reconstructed in uc/AUDIT.md and paper/main.tex")
    print("PRIMARY SOURCE: %s" % SOURCE)
    print("BRIDGE CHECK PASS")


if __name__ == "__main__":
    main()
