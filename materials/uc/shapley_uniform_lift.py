"""PROVED lift: any rational cap-respecting law is an admissible projection.

This kills projection invariants based only on uniformity, the cap, and Reimer.

GENERAL LIFT. Let P be a rational distribution on distinct binary states with
all marginals at most c, and let M be a common denominator of its atom
probabilities, enlarged if necessary.
Create M copies according to those multiplicities. Choose 0<w<M and append to
the copies the M cyclic shifts of the binary word 1^w 0^(M-w). The shifts are
distinct: the cyclic word has exactly two transitions, so any proper period
would replicate the transitions. Each added label coordinate occurs in exactly
w rows. Hence the full rows are distinct, their projection is exactly P, and
all label frequencies equal w/M.

Taking w=floor(cM) preserves a marginal cap c. For sufficiently large M,
2w>=log_2 M, so the labels alone make the uniform lifted family satisfy Reimer.
Thus uniformity, the cap, and Reimer put no further restriction on rational
projected laws. The lifted parent is not union-closed in general.

CERTIFIED INSTANCE. Lift the k=900,d=479 shared-signal distribution from
`shapley_shared_signal.py`. Its projection has A<-0.15. The exact common
denominator M=25000*23087^900 is divisible by 5000, so
w=(1913/5000)M is integral, all full-family marginals are <=0.3826, and
2w>>log_2 M. Even the truncated candidate
A+min(log|supp|-H,0.1) is negative on this projection. Its one-coordinate
projection onto a U copy has positive truncated value. Therefore a universal
2/3 contraction along coordinate-addition chains would preserve positivity and
contradict the negative endpoint.

This does not refute A>=0 for the full uniform family: only its projection is
negative. Nor does it refute projection nonnegativity for uniform
**union-closed** parents. Distinct constant-weight cyclic labels have a union
of larger weight that is absent from the lift, so this parent is not
union-closed. The theorem rules out unrestricted admissible-parent invariants;
a union-closed projection argument may still use the cross-fiber constraint
L_g join L_h subseteq L_{g join h}. The original full-law and UC-restricted
projection questions remain open.

Run: ./.venv/bin/python uc/shapley_uniform_lift.py
"""

from fractions import Fraction
from flint import arb

from shapley_shared_signal import (
    A, B, C, LOG2, certificate, compact, functional, h, law)
# The imported module initializes the shared 4096-bit Arb context.

K = 900
DUPLICATES = 479


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    denominator = A.denominator * B.denominator ** K
    assert denominator == 25000 * 23087 ** K
    assert denominator % A.denominator == 0
    assert denominator % 5000 == 0

    # Every U=0 pattern probability has denominator dividing this product;
    # the binomial theorem checks that all integer multiplicities sum to M.
    total_probability = A + (1 - A) * (B + (1 - B)) ** K
    assert total_probability == 1
    assert (denominator * A).denominator == 1
    for ones in (0, 1, K // 2, K - 1, K):
        probability = (1 - A) * B ** ones * (1 - B) ** (K - ones)
        assert (denominator * probability).denominator == 1

    label_weight = C * denominator
    assert label_weight.denominator == 1
    label_weight = label_weight.numerator
    assert 0 < label_weight < denominator
    assert 2 * label_weight >= denominator.bit_length()
    direct, _reimer, _support_slack, _bad_value, entropy = certificate(
        K, DUPLICATES)
    assert direct.upper() < -0.15
    projected_defect = (arb(2) ** K + 1).log() / LOG2 - entropy
    assert projected_defect.lower() > arb("0.1")
    truncated = direct + arb("0.1")
    assert truncated.upper() < -0.05

    base_direct = functional(law(((A, Fraction(1)),)))
    base_defect = arb(1) - h(A)
    assert base_direct.lower() > 0
    assert base_defect.lower() > arb("0.1")
    base_truncated = base_direct + arb("0.1")
    assert base_truncated.lower() > 0

    projected_coordinates = K + DUPLICATES + 1
    print("PROVED [exact denominator + cyclic balanced-label lift]")
    print("  projected coordinates =", projected_coordinates)
    print("  log10 family size is between", len(str(denominator)) - 1,
          "and", len(str(denominator)))
    print("  full coordinate count = projected coordinates + family size")
    print("  label frequency =", C, "= 0.3826")
    print("  2*label row weight >= bit_length(M):",
          2 * label_weight >= denominator.bit_length())
    print("  projected uniformity defect =", compact(projected_defect), "> 0.1")
    print("  projected direct A =", compact(direct), "< 0")
    print("  projected truncated potential =", compact(truncated), "< 0")
    print("  one-U-coordinate direct A =", compact(base_direct), "> 0")
    print("  one-U-coordinate truncated potential =",
          compact(base_truncated), "> 0")
    print()
    print("FINAL VERDICT")
    print("  REFUTED: projection nonnegativity from cap/Reimer/uniformity alone.")
    print("  REFUTED: unrestricted truncated 2/3 projection contraction.")
    print("  NOT REFUTED: A>=0 for the original full uniform law.")
    print("  NOT REFUTED: projection nonnegativity for union-closed parents.")
