"""CERTIFIED impossibility of any constant KL-corrected 2/3 contraction.

Consider

    B_lambda(P) = A(P) + lambda (log_2 |supp P| - H(P)).

Two explicit projection edges of uniform admissible families impose incompatible
requirements on lambda.

LOWER-BOUND EDGE. Take 63 distinct weight-2 labels in {0,1}^12 and attach four
projected states with multiplicities 14,25,14,10. The projected child rows are
100,000,010,011 and the parent adds a deterministic feature, giving rows
1100,0000,0010,1011. Every full-family coordinate frequency is at most 24/63,
and Reimer holds. Its corrected contraction margin is nonnegative only for
lambda >= 0.774...

UPPER-BOUND EDGE. A complete-[5] projection edge from the admissible parent mask
0x40031443 has parent law

    states 0,3,5,6,8,15 with weights 1/4,1/8,1/8,1/8,1/4,1/8

and child law

    states 0,3,5,6,7 with weights 1/2,1/8,1/8,1/8,1/8.

Its corrected margin is nonnegative only for lambda <= 0.418...

All local laws have at most four rational atoms. Transportation vertices are
exhausted over the rationals and entropy is enclosed with Arb. Therefore no
constant lambda can make KL-corrected 2/3 projection contraction universal.
These edges stay positive; projection nonnegativity for unrestricted admissible
parents is refuted separately by `shapley_uniform_lift.py`.

Run: ./.venv/bin/python uc/shapley_kl_potential.py
"""

from fractions import Fraction
from itertools import combinations

from flint import arb

from shapley_contraction_refute import certified_direct
from shapley_direct_cert import LOG2, ZERO, aq
from shapley_n5_complete import exact_admissible, members
from shapley_n5_projections import project


def shannon(probabilities):
    return -sum((aq(p) * aq(p).log() / LOG2 for p in probabilities), ZERO)


def defect(probabilities):
    return arb(len(probabilities)).log() / LOG2 - shannon(probabilities)


def build_lower_family():
    labels = []
    for support in combinations(range(12), 2):
        labels.append(sum(1 << coordinate for coordinate in support))
    labels = labels[:63]
    assert len(labels) == len(set(labels)) == 63

    projected_states = (12, 0, 2, 11)
    multiplicities = (14, 25, 14, 10)
    rows = []
    offset = 0
    for projected, multiplicity in zip(projected_states, multiplicities):
        for label in labels[offset:offset + multiplicity]:
            rows.append(label | (projected << 12))
        offset += multiplicity
    assert offset == 63
    return tuple(rows), projected_states, multiplicities


def threshold(parent, child, parent_defect, child_defect):
    bare = parent - aq(Fraction(2, 3)) * child
    slope = parent_defect - aq(Fraction(2, 3)) * child_defect
    return bare, slope, -bare / slope


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    lower_family, lower_states, lower_multiplicities = build_lower_family()
    assert len(lower_family) == len(set(lower_family)) == 63
    lower_counts = tuple(sum((row >> coordinate) & 1 for row in lower_family)
                         for coordinate in range(16))
    assert max(lower_counts) == 24
    lower_reimer = (
        arb(2 * sum(lower_counts)) / len(lower_family)
        - arb(len(lower_family)).log() / LOG2
    )
    assert lower_reimer.lower() > arb("0.30")

    lower_probabilities = tuple(
        Fraction(multiplicity, 63) for multiplicity in lower_multiplicities)
    lower_child_states = (4, 0, 2, 3)
    lower_parent, _ = certified_direct(
        lower_states, lower_probabilities, 4)
    lower_child, _ = certified_direct(
        lower_child_states, lower_probabilities, 3)
    assert lower_parent.lower() > 0
    assert lower_child.lower() > 0
    lower_defect = defect(lower_probabilities)
    lower_bare, lower_slope, lambda_lower = threshold(
        lower_parent, lower_child, lower_defect, lower_defect)
    assert lower_bare.upper() < arb("-0.021")
    assert lower_slope.lower() > 0
    assert lambda_lower.lower() > arb("0.77")

    upper_mask = 0x40031443
    upper_family = members(upper_mask)
    assert exact_admissible(upper_family)
    upper_parent_states, upper_parent_probabilities = project(
        upper_family, (1, 2, 3, 4))
    upper_child_states, upper_child_probabilities = project(
        upper_family, (1, 2, 3))
    assert upper_parent_states == (0, 3, 5, 6, 8, 15)
    assert upper_parent_probabilities == (
        Fraction(1, 4), Fraction(1, 8), Fraction(1, 8),
        Fraction(1, 8), Fraction(1, 4), Fraction(1, 8))
    assert upper_child_states == (0, 3, 5, 6, 7)
    assert upper_child_probabilities == (
        Fraction(1, 2), Fraction(1, 8), Fraction(1, 8),
        Fraction(1, 8), Fraction(1, 8))

    upper_parent, _ = certified_direct(
        upper_parent_states, upper_parent_probabilities, 4)
    upper_child, _ = certified_direct(
        upper_child_states, upper_child_probabilities, 3)
    assert upper_parent.lower() > 0
    assert upper_child.lower() > 0
    upper_parent_defect = defect(upper_parent_probabilities)
    upper_child_defect = defect(upper_child_probabilities)
    upper_bare, upper_slope, lambda_upper = threshold(
        upper_parent, upper_child,
        upper_parent_defect, upper_child_defect)
    assert upper_bare.lower() > arb("0.05")
    assert upper_slope.upper() < 0
    assert lambda_upper.upper() < arb("0.42")
    assert lambda_upper.upper() < lambda_lower.lower()

    truncation = arb("0.1")
    assert lower_defect.upper() < truncation
    assert upper_parent_defect.upper() < truncation
    assert upper_child_defect.lower() > truncation
    lower_truncated = lower_bare + lower_defect / 3
    upper_truncated = (
        upper_bare + upper_parent_defect
        - aq(Fraction(2, 3)) * truncation
    )
    assert lower_truncated.lower() > arb("0.006")
    assert upper_truncated.lower() > arb("0.07")

    print("PROVED [two simple-family projection edges + Arb]")
    print("  lower-bound family size / coordinates = 63 / 16")
    print("  coordinate count maximum =", max(lower_counts), "/ 63")
    print("  Reimer slack =", lower_reimer)
    print("  lower edge bare / KL slope =", lower_bare, "/", lower_slope)
    print("  requires lambda >=", lambda_lower)
    print()
    print("  upper edge parent mask = 0x%08x" % upper_mask)
    print("  upper edge bare / KL slope =", upper_bare, "/", upper_slope)
    print("  requires lambda <=", lambda_upper)
    print()
    print("  local two-edge control phi(u)=min(u,0.1):")
    print("    lower-edge corrected margin =", lower_truncated)
    print("    upper-edge corrected margin =", upper_truncated)
    print()
    print("FINAL VERDICT")
    print("  REFUTED: every constant KL weight lambda.")
    print("  The admissible intervals are disjoint: lambda>=0.77 and <=0.42.")
    print("  Truncation repairs these edges only; the unrestricted lift refutes it.")
