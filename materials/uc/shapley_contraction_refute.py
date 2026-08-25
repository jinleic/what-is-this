"""CERTIFIED uniform-family refutation of 2/3 projection contraction.

The complete [5] projection lattice suggested

    A(proj_S X) >= (2/3) A(proj_{S\\{j}} X).

This file gives an explicit admissible uniform family on 13 coordinates whose
3-coordinate projection violates that inequality.

Take all 120 weight-3 labels in {0,1}^10, once each. Attach projected state 000
to 75 labels, 101 to 22 labels, and 110 to 23 labels. The resulting 120 full
rows are distinct. Every label coordinate has frequency 36/120; the projected
coordinates have frequencies 22/120, 23/120, and 45/120. Thus the maximum is
3/8, and the average row weight is 3.75, giving positive Reimer slack.

The 3-coordinate projection has law

    000: 75/120, 101: 22/120, 110: 23/120.

Deleting its OR coordinate gives 00:75/120, 01:22/120, 10:23/120.  All local
laws have at most three rational atoms. Transportation vertices are exhausted
over the rationals and entropy is enclosed with Arb.

This edge refutes quantitative 2/3 contraction but has two positive direct
values; projection nonnegativity for unrestricted admissible parents is refuted
separately by `shapley_uniform_lift.py`. The family is not union-closed. The
displayed KL correction repairs only this edge and is not a surviving universal
method.

Run: ./.venv/bin/python uc/shapley_contraction_refute.py
"""

from fractions import Fraction
from itertools import combinations
from math import comb

from flint import arb

from shapley_direct_cert import LOG2, ZERO, aq, certified_functional


def weighted_law(states, probabilities, coordinate, predecessors):
    groups = {}
    for state, probability in zip(states, probabilities):
        pattern = tuple((state >> j) & 1 for j in predecessors)
        total, ones = groups.get(pattern, (Fraction(), Fraction()))
        groups[pattern] = (
            total + probability,
            ones + probability * ((state >> coordinate) & 1),
        )
    mass = {}
    for total, ones in groups.values():
        atom = ones / total
        mass[atom] = mass.get(atom, Fraction()) + total
    return tuple(sorted(mass.items()))


def certified_direct(states, probabilities, dimension):
    total = ZERO
    vertex_sum = 0
    for coordinate in range(dimension):
        others = tuple(j for j in range(dimension) if j != coordinate)
        for size in range(dimension):
            coefficient = Fraction(1, dimension * comb(dimension - 1, size))
            for predecessors in combinations(others, size):
                law = weighted_law(
                    states, probabilities, coordinate, predecessors)
                value, _coupled, vertices = certified_functional(law)
                total += aq(coefficient) * value
                vertex_sum += vertices
    return total, vertex_sum


def build_family():
    labels = []
    for support in combinations(range(10), 3):
        labels.append(sum(1 << coordinate for coordinate in support))
    assert len(labels) == 120
    rows = []
    for index, label in enumerate(labels):
        projected = 0 if index < 75 else (5 if index < 97 else 6)
        rows.append(label | (projected << 10))
    return tuple(rows)


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    family = build_family()
    assert len(family) == len(set(family)) == 120
    counts = tuple(sum((row >> coordinate) & 1 for row in family)
                   for coordinate in range(13))
    assert counts[:10] == (36,) * 10
    assert counts[10:] == (22, 23, 45)
    assert max(counts) * 8 == 3 * len(family)

    twice_average_size = arb(2 * sum(counts)) / len(family)
    support_entropy = arb(len(family)).log() / arb(2).log()
    reimer_slack = twice_average_size - support_entropy
    assert reimer_slack.lower() > arb("0.59")

    parent_states = (0, 5, 6)
    child_states = (0, 1, 2)
    probabilities = (Fraction(75, 120), Fraction(22, 120), Fraction(23, 120))
    parent, parent_vertices = certified_direct(
        parent_states, probabilities, 3)
    child, child_vertices = certified_direct(
        child_states, probabilities, 2)
    margin = parent - aq(Fraction(2, 3)) * child
    ratio = parent / child
    entropy = -sum((aq(p) * aq(p).log() / LOG2 for p in probabilities), ZERO)
    uniformity_defect = arb(3).log() / LOG2 - entropy
    corrected_margin = (
        parent + aq(Fraction(1, 10)) * uniformity_defect
        - aq(Fraction(2, 3))
        * (child + aq(Fraction(1, 10)) * uniformity_defect)
    )
    assert parent.lower() > 0
    assert child.lower() > 0
    assert margin.upper() < arb("-0.0077")
    assert ratio.upper() < arb("0.651")
    assert corrected_margin.lower() > arb("0.00079")

    print("PROVED [explicit simple family + rational LP vertices + Arb]")
    print("  family size / coordinates = 120 / 13")
    print("  coordinate counts =", counts)
    print("  Reimer slack =", reimer_slack)
    print("  projected parent law =", (parent_states, probabilities))
    print("  projected child law  =", (child_states, probabilities))
    print("  transport-vertex sums parent / child =",
          parent_vertices, "/", child_vertices)
    print("  A(parent) =", parent)
    print("  A(child)  =", child)
    print("  ratio     =", ratio)
    print("  A(parent)-(2/3)A(child) =", margin, "< 0")
    print("  uniformity defect log|supp|-H =", uniformity_defect)
    print("  KL-corrected 2/3 margin =", corrected_margin, "> 0")
    print()
    print("FINAL VERDICT")
    print("  REFUTED: universal 2/3 projection contraction.")
    print("  This edge has positive values; unrestricted admissible-parent")
    print("  projection nonnegativity is refuted by uniform_lift.py.")
    print("  The KL correction repairs this edge only; no universal claim survives.")
