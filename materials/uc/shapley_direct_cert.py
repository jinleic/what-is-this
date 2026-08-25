"""CERTIFIED direct-average values for the n=4 and n=5 Shapley extremizers.

`shapley_direct.py` identifies the n=4 family

    B = {0,1,2,4,7,9,10,12},

while `shapley_n5_complete.py` finds two lower normalized n=5 orbits, represented
by a duplicated-column family and a separating family.  This file certifies the
three displayed values with no float transport optimization.

Each displayed family has only six or seven exact laws mu_{i,S}, supported
on at most three rational atoms.  For each law, enumerate every vertex of its
transportation polytope: a vertex with positive marginals has between m and
2m-1 support edges, and its support columns are linearly independent.  Every
candidate plan is solved and checked over the rationals.  A linear transport
cost attains its minimum at one of these exhaustively enumerated vertices.
Binary-entropy costs and the final minimum are enclosed with Arb at 256 bits.

The same exact machinery also refutes a tempting proof shortcut.  Two
conditionally independent one-sided signals of a Bernoulli variable can make
both complementary conditional laws equal to the same negative-cost law, even
when every marginal is at most 0.38261.  Thus complementary predecessor sets
cannot be bounded pairwise; Reimer must enter the aggregate argument.

This certifies the values of the displayed families, not their minimality in the
complete enumerations; those comparisons still use float64 HiGHS in
`shapley_direct.py` and `shapley_n5_complete.py`.

Run: ./.venv/bin/python uc/shapley_direct_cert.py
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations

from flint import arb, ctx
import sympy as sp

from shapley_direct import BASE4, predecessor_law, predecessor_weight

ctx.prec = 256
ZERO = arb(0)
ONE = arb(1)
LOG2 = arb(2).log()
ALPHA = arb("0.0356069")
DUPLICATE5 = (0, 1, 2, 4, 7, 25, 26, 28)
SEPARATING5 = (0, 1, 6, 10, 13, 19, 20, 24)



def aq(value):
    return arb(value.numerator) / value.denominator


def h(value):
    if value == 0 or value == 1:
        return ZERO
    x = aq(value)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def sstar(p, q):
    return max(p, q, min(p + q, Fraction(1, 2)))


def transport_vertices(weights):
    """Enumerate all vertices of the equal-marginal transport polytope."""
    m = len(weights)
    edges = tuple((i, j) for i in range(m) for j in range(m))
    rhs = [sp.Rational(w.numerator, w.denominator) for w in weights]
    # The last column equation is redundant with all row equations.
    rhs += [sp.Rational(weights[j].numerator, weights[j].denominator)
            for j in range(m - 1)]
    vertices = set()

    for support_size in range(m, 2 * m):
        for support in combinations(edges, support_size):
            equations = []
            for i in range(m):
                equations.append([int(edge[0] == i) for edge in support])
            for j in range(m - 1):
                equations.append([int(edge[1] == j) for edge in support])
            matrix = sp.Matrix(equations)
            vector = sp.Matrix(rhs)
            if (matrix.rank() != support_size
                    or matrix.row_join(vector).rank() != support_size):
                continue
            solution = next(iter(sp.linsolve((matrix, vector))))
            if any(value < 0 for value in solution):
                continue

            plan = [Fraction() for _ in edges]
            for edge, value in zip(support, solution):
                plan[edges.index(edge)] = Fraction(int(value.p), int(value.q))
            if not all(sum(plan[i * m + j] for j in range(m)) == weights[i]
                       for i in range(m)):
                continue
            if not all(sum(plan[i * m + j] for i in range(m)) == weights[j]
                       for j in range(m)):
                continue
            vertices.add(tuple(plan))
    return tuple(sorted(vertices))


def amin_endpoint(values, endpoint):
    best = endpoint(values[0])
    for value in values[1:]:
        candidate = endpoint(value)
        if candidate < best:
            best = candidate
    return best


def certified_transport(atoms, weights):
    vertices = transport_vertices(weights)
    m = len(atoms)
    costs = []
    for plan in vertices:
        cost = ZERO
        for i in range(m):
            for j in range(m):
                cost += aq(plan[i * m + j]) * h(sstar(atoms[i], atoms[j]))
        costs.append(cost)
    lower = amin_endpoint(costs, lambda value: value.lower())
    upper = amin_endpoint(costs, lambda value: value.upper())
    return lower.union(upper), len(vertices)


def certified_functional(law):
    atoms = tuple(atom for atom, _ in law)
    weights = tuple(weight for _, weight in law)
    entropy = sum((aq(weight) * h(atom) for atom, weight in law), ZERO)
    iid = sum(
        (aq(weights[i] * weights[j])
         * h(atoms[i] + atoms[j] - atoms[i] * atoms[j])
         for i in range(len(law)) for j in range(len(law))),
        ZERO,
    )
    coupled, vertex_count = certified_transport(atoms, weights)
    value = (ONE - ALPHA) * iid + ALPHA * coupled - entropy
    return value, coupled, vertex_count


def base_law_coefficients():
    coefficients = defaultdict(Fraction)
    for coordinate in range(4):
        others = tuple(j for j in range(4) if j != coordinate)
        for size in range(4):
            for predecessors in combinations(others, size):
                law = predecessor_law(BASE4, coordinate, predecessors)
                coefficients[law] += predecessor_weight(4, size)
    assert sum(coefficients.values(), Fraction()) == 4
    assert len(coefficients) == 6
    return coefficients

def law_coefficients(mem, n):
    coefficients = defaultdict(Fraction)
    for coordinate in range(n):
        others = tuple(j for j in range(n) if j != coordinate)
        for size in range(n):
            for predecessors in combinations(others, size):
                law = predecessor_law(mem, coordinate, predecessors)
                coefficients[law] += predecessor_weight(n, size)
    assert sum(coefficients.values(), Fraction()) == n
    return coefficients


def certified_direct_value(mem, n):
    coefficients = law_coefficients(mem, n)
    total = ZERO
    vertex_count = 0
    for law, coefficient in coefficients.items():
        value, _coupled, vertices = certified_functional(law)
        total += aq(coefficient) * value
        vertex_count += vertices
    return total, len(coefficients), vertex_count

def signal_pair_refutation():
    """Return an exact negative complementary-pair law and Reimer deficit."""
    c = Fraction(1913, 5000)
    signal_rate = Fraction(1, 5)
    signal_mass = c * signal_rate
    posterior_zero = c * (1 - signal_rate) / (1 - signal_mass)
    law = (
        (posterior_zero, 1 - signal_mass),
        (Fraction(1), signal_mass),
    )
    assert sum(atom * weight for atom, weight in law) == c
    value, coupled, vertex_count = certified_functional(law)
    assert value.upper() < 0

    # X~Bernoulli(c); U,V are conditionally iid Bernoulli(signal_rate) when
    # X=1 and are zero when X=0.  Then law(X|U)=law(X|V)=law above.
    joint_entropy = h(c) + 2 * aq(c) * h(signal_rate)
    twice_mean_sum = 2 * aq(c * (1 + 2 * signal_rate))
    reimer_slack = twice_mean_sum - joint_entropy
    assert reimer_slack.upper() < 0
    return law, value, coupled, vertex_count, reimer_slack


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    coefficients = base_law_coefficients()
    total = ZERO
    total_vertices = 0
    for law in sorted(coefficients, key=lambda item: (len(item), item)):
        value, coupled, vertex_count = certified_functional(law)
        total += aq(coefficients[law]) * value
        total_vertices += vertex_count
        print("law coefficient=%s atoms=%s" % (coefficients[law], law))
        print("  transport vertices:", vertex_count)
        print("  C =", coupled)
        print("  F =", value)

    gamma = total / 3  # H(uniform on BASE4) = log_2(8) = 3 exactly.
    assert total.lower() > arb("0.12987354972888714")
    assert total.upper() < arb("0.12987354972888715")
    assert gamma.lower() > arb("0.04329118324296238")
    assert gamma.upper() < arb("0.04329118324296239")

    print()
    print("PROVED [exact rational transport vertices + Arb entropy]")
    print("  distinct local laws:", len(coefficients))
    print("  total enumerated transport vertices:", total_vertices)
    print("  A(BASE4) =", total)
    print("  A(BASE4)/H(BASE4) =", gamma)
    print("  certified decimal enclosure:")
    print("    0.04329118324296238 < gamma < 0.04329118324296239")
    duplicate_value, duplicate_laws, duplicate_vertices = \
        certified_direct_value(DUPLICATE5, 5)
    duplicate_ratio = duplicate_value / 3
    separating_value, separating_laws, separating_vertices = \
        certified_direct_value(SEPARATING5, 5)
    separating_ratio = separating_value / 3
    assert duplicate_ratio.lower() > arb("0.04078419694312722")
    assert duplicate_ratio.upper() < arb("0.04078419694312723")
    assert separating_ratio.lower() > arb("0.04245552114301732")
    assert separating_ratio.upper() < arb("0.04245552114301734")
    print()
    print("PROVED [n=5 representative values]")
    print("  duplicated-column family =", DUPLICATE5)
    print("    local laws / vertex sum =", duplicate_laws, "/", duplicate_vertices)
    print("    A =", duplicate_value)
    print("    A/H =", duplicate_ratio)
    print("  separating family =", SEPARATING5)
    print("    local laws / vertex sum =", separating_laws, "/", separating_vertices)
    print("    A =", separating_value)
    print("    A/H =", separating_ratio)
    print("  Both certified ratios are strictly below the n=4 gamma.")
    law, pair_value, pair_coupled, pair_vertices, pair_slack = \
        signal_pair_refutation()
    print()
    print("PROVED BOUNDARY [complement-pair shortcut is false]")
    print("  X marginal = 1913/5000; each signal marginal = 1913/25000")
    print("  law(X|U) = law(X|V) =", law)
    print("  transport vertices:", pair_vertices)
    print("  C =", pair_coupled)
    print("  F(law) =", pair_value, "< 0")
    print("  F(law(X|U))+F(law(X|V)) =", 2 * pair_value, "< 0")
    print("  Reimer slack of (X,U,V) =", pair_slack, "< 0")
    print("  Therefore any surviving global nonnegativity theorem must use")
    print("  Reimer in aggregate, not nonnegative complementary pairs.")
    print()
    print("SCOPE: displayed-family values are certified; enumerated minimality")
    print("remains the complete float64 comparison in the companion scripts.")
