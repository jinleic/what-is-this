"""GLOBAL join-loss identity on uniform union-closed families.

Let X,Y be independent uniform rows of a finite family F and Z=X union Y.  For
a coordinate order pi, write W_j=(X_j,Y_j).  The iid Shapley term is

    Q_pi = sum_i H(Z_i | W_pred(i)).

Since Z_pred is a function of W_pred,

    Q_pi = H(Z) - J_pi,
    J_pi = sum_i I(Z_i ; W_pred(i) | Z_pred(i)) >= 0.

If F is union-closed, Z is supported on F.  Diagonal pairs show every row of F
occurs, and therefore

    D_join = D_2(P_Z || Uniform(F)) = log2|F| - H(Z) >= 0.

Averaging over all orders gives the exact full-law identity

    log2|F| - Q_Shapley = D_join + J_Shapley.                 (1)

Thus union closure moves the iid Q term in the wrong direction: Q<=H.  Any
beyond-c* proof must make the globally consistent dependent strategy pay both
the join-concentration loss D_join and the hidden-prefix leakage J_Shapley.
This identity is not projection-local and is untouched by balanced labels.

The script exhausts all 4,959 union-closed families on [4].  It checks (1), the
nonnegativity of both losses, and the discovery inequality

    J_Shapley <= D_join/2.

That inequality holds throughout [4] with maximum ratio about 0.443589, but an
explicit union-closed family on [5] refutes it with ratio 0.551652.  The weaker
candidate J_Shapley<=D_join survives this probe but remains conjectural.
Entropy is evaluated in float64; pair counts and conditional laws are exact.

Run:
    ./.venv/bin/python uc/shapley_join_loss.py
"""

from fractions import Fraction
from itertools import combinations
from math import comb, log2

from shapley_direct import predecessor_law
from shapley_entropy import union_closed

N = 4
SUBSETS = 1 << N


def binary_entropy(value):
    value = float(value)
    if value <= 0.0 or value >= 1.0:
        return 0.0
    return -value * log2(value) - (1.0 - value) * log2(1.0 - value)


def family_members(mask):
    return tuple(row for row in range(SUBSETS) if (mask >> row) & 1)


def iid_local_cost(law):
    return sum(
        float(left_mass * right_mass)
        * binary_entropy(left + right - left * right)
        for left, left_mass in law
        for right, right_mass in law
    )


def shapley_iid(family, dimension):
    total = 0.0
    for coordinate in range(dimension):
        others = tuple(i for i in range(dimension) if i != coordinate)
        for size in range(dimension):
            weight = Fraction(
                1, dimension * comb(dimension - 1, size))
            for predecessors in combinations(others, size):
                total += float(weight) * iid_local_cost(
                    predecessor_law(family, coordinate, predecessors)
                )
    return total


def join_entropy(family):
    pair_count = len(family) ** 2
    counts = {}
    for left in family:
        for right in family:
            joined = left | right
            counts[joined] = counts.get(joined, 0) + 1
    assert set(counts) == set(family)
    return -sum(
        Fraction(count, pair_count) * log2(count / pair_count)
        for count in counts.values()
    )


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    count = 0
    maximum_residual = 0.0
    maximum_ratio = (-1.0, None, None, None)
    maximum_hidden = (-1.0, None)
    maximum_divergence = (-1.0, None)

    for mask in range(1, 1 << SUBSETS):
        family = family_members(mask)
        if not union_closed(family):
            continue
        count += 1
        entropy = log2(len(family))
        output_entropy = join_entropy(family)
        iid = shapley_iid(family, N)
        divergence = entropy - output_entropy
        hidden = output_entropy - iid
        residual = abs((entropy - iid) - (divergence + hidden))
        maximum_residual = max(maximum_residual, residual)
        assert divergence >= -2e-14
        assert hidden >= -2e-14
        if divergence > 1e-13:
            ratio = hidden / divergence
            if ratio > maximum_ratio[0]:
                maximum_ratio = (ratio, family, hidden, divergence)
        if hidden > maximum_hidden[0]:
            maximum_hidden = (hidden, family)
        if divergence > maximum_divergence[0]:
            maximum_divergence = (divergence, family)

    assert count == 4959
    assert maximum_residual < 3e-15
    assert maximum_ratio[0] < 0.5
    assert union_closed(maximum_ratio[1])
    assert abs(maximum_ratio[0] - 0.443588993565) < 1e-12

    counterexample = (0, 13, 17, 26, 27, 29, 31)
    assert union_closed(counterexample)
    counter_entropy = log2(len(counterexample))
    counter_output_entropy = join_entropy(counterexample)
    counter_iid = shapley_iid(counterexample, 5)
    counter_divergence = counter_entropy - counter_output_entropy
    counter_hidden = counter_output_entropy - counter_iid
    counter_ratio = counter_hidden / counter_divergence
    assert counter_ratio > 0.55
    assert abs(counter_ratio - 0.551652224828) < 1e-12

    print("PROVED [identity] + COMPLETE NUMERICAL [4] CENSUS")
    print("  union-closed families:", count)
    print("  maximum identity residual = %.3e" % maximum_residual)
    print("  max J/D = %.12f" % maximum_ratio[0])
    print("    family =", maximum_ratio[1])
    print("    J / D = %.12f / %.12f" %
          (maximum_ratio[2], maximum_ratio[3]))
    print("  max J = %.12f at %s" % maximum_hidden)
    print("  max D = %.12f at %s" % maximum_divergence)
    print()
    print("FINITE VERDICT")
    print("  J<=D/2 holds on every union-closed family on [4], but is")
    print("  REFUTED on [5] by", counterexample)
    print("  counterexample J / D / ratio = %.12f / %.12f / %.12f" %
          (counter_hidden, counter_divergence, counter_ratio))
    print("  The weaker J<=D remains a live conjecture; no theorem is claimed.")
