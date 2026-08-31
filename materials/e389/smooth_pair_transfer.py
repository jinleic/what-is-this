#!/usr/bin/env python3
"""Transfer Pascadi's 5/8 theorem to dense consecutive smooth pairs.

Yang (arXiv:2607.16032, Theorem 1.5) uses the 66/107 exponent of distribution
for smooth numbers to prove a positive density of consecutive integers whose
largest prime factors are below x^(41/107+epsilon).  Pascadi's stronger
Theorem 1.5 (arXiv:2505.00653v2) supplies the same absolute-value
Bombieri--Vinogradov estimate through x^(5/8-epsilon).  The identical two-prime
construction therefore gives the sharper exponent

    1 - 5/8 = 3/8.

This producer certifies the parameter algebra exactly.  For eta=5/16 and a
small delta, take

    p1 in (x^(eta-3 delta), x^(eta-2 delta)],
    p2 in (x^(eta-2 delta), x^(eta-delta)].

Then p1*p2 <= x^(5/8-3 delta), strictly inside Pascadi's level
x^(5/8-2 delta), while the remaining cofactor of n+1 is at most
x^(3/8+5 delta).  Choosing delta < epsilon/6 puts every prime factor below
x^(3/8+epsilon).  Mertens gives a positive product of two logarithms and
Pascadi makes the total progression error o(Psi(x,y)); hence the family has
positive density.  This is a theorem derived from a cited external input, not a
finite experiment.

The remaining reports expand the selected-prime mask into zero, one-sided and
primitive bilinear frequencies; certify the $M^2$ lift and the exact paired
q-adic Poisson/reciprocity identities; audit Pascadi,
Drappeau--Shparlinski, square-moduli large sieves, and Chen switching at the
actual variables; price Parseval, tails, and automatic-compensation fallbacks;
and run two exact finite censuses.  The experiments are evidence only and are
never inputs to the smooth-pair theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import SPFFactorizer
from good_density_decomposition import classify_term, compensation_count

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "smooth_pair_transfer.json"

YANG_LEVEL = Fraction(66, 107)
YANG_SMOOTH_EXPONENT = Fraction(41, 107)
PASCADI_LEVEL = Fraction(5, 8)
PAIR_PRIME_EXPONENT = Fraction(5, 16)
TRANSFER_SMOOTH_EXPONENT = Fraction(3, 8)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def rational_power_floor(value: int, numerator: int, denominator: int) -> int:
    """Largest integer y with ``y**denominator <= value**numerator``."""
    _require_int("value", value, 1)
    _require_int("numerator", numerator, 1)
    _require_int("denominator", denominator, 1)
    target = value**numerator
    low, high = 0, 1
    while high**denominator <= target:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**denominator <= target:
            low = middle
        else:
            high = middle
    return low



def drappeau_shparlinski_linear_saving(
    modulus_exponent: Fraction,
) -> Fraction | None:
    """Power saving from Theorem 1.1 for y=x^(1/C), C large.

    The three relevant savings are 1/5, theta/2 and (1-theta)/2.
    A modulus exponent theta >= 1 is outside the nontrivial range.
    """
    if not isinstance(modulus_exponent, Fraction) or modulus_exponent <= 0:
        raise ValueError("modulus_exponent must be a positive Fraction")
    if modulus_exponent >= 1:
        return None
    return min(
        Fraction(1, 5),
        modulus_exponent / 2,
        (1 - modulus_exponent) / 2,
    )


def transfer_parameters(epsilon: Fraction = Fraction(1, 200)) -> dict[str, Any]:
    """Exact exponent certificate for the 5/8 -> 3/8 transfer."""
    if not isinstance(epsilon, Fraction) or not 0 < epsilon < Fraction(1, 8):
        raise ValueError("epsilon must be a Fraction in (0, 1/8)")
    delta = min(epsilon / 6, Fraction(1, 100))
    eta = PAIR_PRIME_EXPONENT
    p1_lower = eta - 3 * delta
    p1_upper = eta - 2 * delta
    p2_lower = p1_upper
    p2_upper = eta - delta
    modulus_lower = p1_lower + p2_lower
    modulus_upper = p1_upper + p2_upper
    distribution_limit = PASCADI_LEVEL - 2 * delta
    cofactor_upper = 1 - modulus_lower
    target_upper = TRANSFER_SMOOTH_EXPONENT + epsilon
    improves_yang = target_upper < YANG_SMOOTH_EXPONENT
    max_one_sided_lift = max(
        2 * p1_upper + p2_upper,
        p1_upper + 2 * p2_upper,
    )
    one_sided_individual_saving = drappeau_shparlinski_linear_saving(
        max_one_sided_lift
    )
    if one_sided_individual_saving is None:
        raise AssertionError("one-sided lift left the nontrivial theorem range")
    one_sided_naive_sum_exponent = (
        1 - one_sided_individual_saving + modulus_upper
    )
    one_sided_naive_sum_excess = one_sided_naive_sum_exponent - 1
    max_bilinear_lift = 2 * modulus_upper
    bilinear_lift_excess = max_bilinear_lift - 1
    if drappeau_shparlinski_linear_saving(max_bilinear_lift) is not None:
        raise AssertionError("bilinear lift unexpectedly entered theorem range")

    checks = {
        "disjoint_prime_ranges": p1_upper == p2_lower,
        "at_most_three_selected_primes": p1_lower > Fraction(1, 4),
        "modulus_inside_pascadi_level": modulus_upper < distribution_limit,
        "cofactor_below_target": cofactor_upper < target_upper,
    }
    if not all(checks.values()):
        raise AssertionError(f"invalid transfer parameters: {checks}")

    # Mertens: sum 1/p over (x^a, x^b] -> log(b/a).  The product is the
    # positive coefficient multiplying Psi(x, x^(1/C)); C itself is supplied
    # existentially by Pascadi's theorem and is deliberately not invented here.
    harmonic_p1 = math.log(float(p1_upper / p1_lower))
    harmonic_p2 = math.log(float(p2_upper / p2_lower))
    main_coefficient_without_dickman = 0.5 * harmonic_p1 * harmonic_p2
    if main_coefficient_without_dickman <= 0.0:
        raise AssertionError("Mertens main coefficient is not positive")

    return {
        "epsilon": str(epsilon),
        "delta": str(delta),
        "eta": str(eta),
        "p1_exponents": [str(p1_lower), str(p1_upper)],
        "p2_exponents": [str(p2_lower), str(p2_upper)],
        "modulus_exponents": [str(modulus_lower), str(modulus_upper)],
        "pascadi_distribution_limit": str(distribution_limit),
        "cofactor_upper_exponent": str(cofactor_upper),
        "target_smooth_exponent": str(target_upper),
        "geometric_C_lower_bound": str(1 / p1_lower),
        "main_coefficient_without_dickman": main_coefficient_without_dickman,
        "improves_yang_exponent_at_this_epsilon": improves_yang,
        "max_one_sided_lift_exponent": str(max_one_sided_lift),
        "one_sided_individual_sum_saving": str(one_sided_individual_saving),
        "one_sided_naive_sum_exponent": str(one_sided_naive_sum_exponent),
        "one_sided_naive_sum_excess_over_main": str(one_sided_naive_sum_excess),
        "max_bilinear_lift_exponent": str(max_bilinear_lift),
        "bilinear_lift_excess_over_n_range": str(bilinear_lift_excess),
        "checks": checks,
    }


def _first_obstruction(
    value: int, factorization: dict[int, int], target_m: int = 1
) -> tuple[int, int] | None:
    """Largest obstructing prime and its exponent, or ``None``."""
    for prime, exponent in sorted(factorization.items(), reverse=True):
        if prime <= target_m:
            continue
        cofactor = value // prime**exponent
        if compensation_count(cofactor, prime) < exponent:
            return prime, exponent
    return None


def _band(value: int, prime: int) -> int:
    """Number of usable levels for an exponent-one prime at this scale."""
    count, power = 0, prime * prime
    while power <= 2 * value:
        count += 1
        power *= prime
    return count


def measure_block(
    *, start: int, count: int, target_m: int = 1
) -> dict[str, Any]:
    """Exact census of consecutive 3/8-smooth pairs and level survival."""
    _require_int("start", start, 2)
    _require_int("count", count, 1)
    _require_int("target_m", target_m, 1)
    top = start + count + 1
    smooth_bound = rational_power_floor(top, 3, 8)
    factorizer = SPFFactorizer(top + 1)
    factorizations: list[dict[int, int]] = []
    goodness: list[bool] = []
    for value in range(start, top):
        factorization = factorizer.factor(value)
        factorizations.append(factorization)
        short, level = classify_term(value, target_m, factorizer)
        goodness.append(not short and not level)

    smooth_pairs = 0
    smooth_good_pairs = 0
    all_good_pairs = 0
    failure_patterns: dict[str, int] = {}
    blocker_tallies: dict[tuple[int, int, int], int] = {}
    for offset in range(count):
        left_factorization = factorizations[offset]
        right_factorization = factorizations[offset + 1]
        left_smooth = max(left_factorization, default=1) <= smooth_bound
        right_smooth = max(right_factorization, default=1) <= smooth_bound
        all_good_pairs += goodness[offset] and goodness[offset + 1]
        if not (left_smooth and right_smooth):
            continue
        smooth_pairs += 1
        smooth_good_pairs += goodness[offset] and goodness[offset + 1]
        local = []
        for side, (value, factorization) in enumerate(
            ((start + offset, left_factorization),
             (start + offset + 1, right_factorization))
        ):
            obstruction = _first_obstruction(value, factorization, target_m)
            if obstruction is None:
                local.append("good")
                continue
            prime, exponent = obstruction
            band = _band(value, prime)
            local.append(f"b{band}e{exponent}")
            key = (side, band, exponent)
            blocker_tallies[key] = blocker_tallies.get(key, 0) + 1
        pattern = "/".join(local)
        failure_patterns[pattern] = failure_patterns.get(pattern, 0) + 1

    if smooth_pairs == 0:
        raise AssertionError("screened block contains no 3/8-smooth pairs")
    return {
        "start": start,
        "count": count,
        "theta": "3/8",
        "smooth_bound": smooth_bound,
        "smooth_pairs": smooth_pairs,
        "smooth_pair_density": smooth_pairs / count,
        "smooth_good_pairs": smooth_good_pairs,
        "smooth_and_good_pair_density": smooth_good_pairs / count,
        "conditional_good_given_smooth_pair": smooth_good_pairs / smooth_pairs,
        "all_good_pair_density": all_good_pairs / count,
        "failure_pattern_count": len(failure_patterns),
        "failure_patterns": [
            {"pattern": pattern, "count": tally}
            for pattern, tally in sorted(
                failure_patterns.items(), key=lambda item: (-item[1], item[0])
            )[:20]
        ],
        "first_blocker_class_count": len(blocker_tallies),
        "first_blockers": [
            {"side": key[0], "band": key[1], "exponent": key[2], "count": tally}
            for key, tally in sorted(
                blocker_tallies.items(), key=lambda item: (-item[1], item[0])
            )[:20]
        ],
    }


def half_lift_fourier_l1_bound(prime: int) -> Fraction:
    """Bound the normalized Fourier l1 norm of a centered half-lift."""
    _require_int("prime", prime, 3)
    if prime % 2 == 0:
        raise ValueError("prime must be odd")
    return sum(
        (Fraction(1, index) for index in range(1, (prime - 1) // 2 + 1)),
        start=Fraction(0),
    )



def double_mask_frequency_classes(first_prime: int, second_prime: int) -> dict[str, Any]:
    """Exact zero/one-sided/bilinear frequency partition of the double mask."""
    _require_int("first_prime", first_prime, 3)
    _require_int("second_prime", second_prime, 3)
    if (
        first_prime % 2 == 0
        or second_prime % 2 == 0
        or math.gcd(first_prime, second_prime) != 1
    ):
        raise ValueError("inputs must be odd and coprime")
    modulus = first_prime * second_prime
    primitive_bilinear = 0
    for first_frequency in range(1, first_prime):
        for second_frequency in range(1, second_prime):
            numerator = (
                first_frequency * second_prime**2
                + second_frequency * first_prime**2
            )
            if math.gcd(numerator, modulus) != 1:
                raise AssertionError("bilinear frequency is not primitive")
            primitive_bilinear += 1
    classes = {
        "zero": 1,
        "one_sided_first": first_prime - 1,
        "one_sided_second": second_prime - 1,
        "genuinely_bilinear": primitive_bilinear,
    }
    if sum(classes.values()) != modulus:
        raise AssertionError("frequency partition does not cover the full grid")
    return {
        "first_prime": first_prime,
        "second_prime": second_prime,
        "quotient_modulus": modulus,
        "classes": classes,
        "bilinear_frequencies_primitive_mod_product": True,
        "one_sided_quotient_periods": [first_prime, second_prime],
        "bilinear_quotient_period": modulus,
        "one_sided_n_lift_moduli": [
            first_prime**2 * second_prime,
            first_prime * second_prime**2,
        ],
        "bilinear_n_lift_modulus": modulus**2,
        "phase": (
            "e_{p1*p2}((h1*p2^2+h2*p1^2)*ell); h2=0 reduces to p1, "
            "h1=0 reduces to p2, h1*h2!=0 is primitive mod p1*p2"
        ),
    }




def q_adic_lift_report(modulus: int, frequency: int) -> dict[str, Any]:
    """Exact support of a quotient character after lifting from M to M^2.

    For gcd(frequency, modulus)=1,

      1_{s=-1 mod M} e_M(frequency*(s+1)/M)
        = M^{-1} sum_{j mod M} e_{M^2}((frequency+jM)*(s+1)).

    Every lifted frequency remains primitive modulo M^2.
    """
    _require_int("modulus", modulus, 2)
    _require_int("frequency", frequency, 1)
    if math.gcd(modulus, frequency) != 1:
        raise ValueError("frequency must be coprime to modulus")
    lifted = [frequency + index * modulus for index in range(modulus)]
    primitive = all(
        math.gcd(value, modulus * modulus) == 1 for value in lifted
    )
    if not primitive:
        raise AssertionError("a lifted primitive frequency became imprimitive")
    return {
        "base_modulus": modulus,
        "base_frequency": frequency,
        "lift_modulus": modulus * modulus,
        "lift_frequency_count": modulus,
        "all_lift_frequencies_primitive": primitive,
        "identity": (
            "1{s=-1 mod M}*e_M(a*(s+1)/M)="
            "(1/M)*sum_{j mod M} e_{M^2}((a+jM)*(s+1))"
        ),
    }


def q_adic_poisson_reindex_report(
    modulus: int, multiplier: int, cutoff: int
) -> dict[str, Any]:
    """Certify the exact reindexing in the q-adic Poisson formula.

    If K*multiplier=1 mod modulus^2, each dual integer j corresponds to
    a=-j*K mod modulus and h=(j+a*multiplier)/modulus, so
    j=modulus*h-a*multiplier.  This is the moving-center term missed by the
    unweighted completion.
    """
    _require_int("modulus", modulus, 2)
    _require_int("multiplier", multiplier, 2)
    _require_int("cutoff", cutoff, 0)
    if math.gcd(modulus, multiplier) != 1:
        raise ValueError("multiplier must be coprime to modulus")
    lift_modulus = modulus * modulus
    inverse = pow(multiplier, -1, lift_modulus)
    inverse_modulus = inverse % modulus
    inverse_lift_mod_multiplier = pow(lift_modulus, -1, multiplier)
    residue_counts = [0] * modulus
    reciprocity_holds = True
    for dual in range(-cutoff, cutoff + 1):
        frequency = (-dual * inverse_modulus) % modulus
        numerator = dual + frequency * multiplier
        if numerator % modulus:
            raise AssertionError("dual frequency did not reindex integrally")
        poisson_frequency = numerator // modulus
        if dual != modulus * poisson_frequency - frequency * multiplier:
            raise AssertionError("moving-center identity failed")
        residue_counts[frequency] += 1

        # -K*j/r^2 = j*bar(r^2)/k - j/(k*r^2) (mod 1).
        reciprocity_numerator = (
            -inverse * dual * multiplier
            - dual * inverse_lift_mod_multiplier * lift_modulus
            + dual
        )
        reciprocity_holds &= (
            reciprocity_numerator % (multiplier * lift_modulus) == 0
        )
    second_multiplier = multiplier + modulus
    second_quotient_base = (
        1 - inverse * second_multiplier
    ) // modulus
    pair_phase_holds = True
    for first_frequency in range(modulus):
        for second_frequency in range(modulus):
            for poisson_frequency in range(-2, 3):
                dual = (
                    modulus * poisson_frequency
                    - first_frequency * multiplier
                    + second_frequency * second_multiplier
                )
                original_phase = (
                    -inverse * dual
                    - second_frequency * modulus * second_quotient_base
                ) % lift_modulus
                reduced_phase = (
                    first_frequency
                    - second_frequency
                    - inverse * poisson_frequency * modulus
                ) % lift_modulus
                pair_phase_holds &= original_phase == reduced_phase
    return {
        "modulus": modulus,
        "multiplier": multiplier,
        "lift_modulus": lift_modulus,
        "inverse_mod_lift": inverse,
        "cutoff": cutoff,
        "dual_count": 2 * cutoff + 1,
        "max_repetitions_per_mask_frequency": max(residue_counts),
        "min_repetitions_per_mask_frequency": min(residue_counts),
        "all_dual_frequencies_reindexed": True,
        "additive_reciprocity_verified": reciprocity_holds,
        "pair_poisson_phase_verified": pair_phase_holds,
        "pair_poisson_identity": (
            "center=(a*k1-b*k2)/r and phase="
            "e_{r^2}(a-b)*e_r(-h*inverse(k1))"
        ),
        "poisson_identity": (
            "sum_{mk=-1 mod r} f(m)W((mk+1)/r)="
            "(1/r)sum_j fhat(j/r^2)e_{r^2}(-Kj)What(-jK mod r)"
        ),
        "reciprocity_identity": (
            "e_{r^2}(-Kj)=e_k(j*inverse(r^2))*e(-j/(k*r^2))"
        ),
    }


def mixed_dispersion_exponents(
    level: Fraction = PASCADI_LEVEL,
) -> dict[str, Any]:
    """Exact exponent map after inserting the mask before Poisson."""
    if (
        not isinstance(level, Fraction)
        or not Fraction(1, 2) < level < Fraction(2, 3)
    ):
        raise ValueError("level must be a Fraction in (1/2, 2/3)")
    factor = 1 - level
    native_dual = 2 * level - 1
    selected_prime = level / 2
    spectral_level = 2 * native_dual
    mixed_support = native_dual + level
    support_inflation = mixed_support - spectral_level
    pre_poisson_gap = selected_prime - native_dual
    additive_saving_benchmark = native_dual
    residual_after_additive_benchmark = (
        support_inflation - additive_saving_benchmark
    )
    return {
        "outer_level": str(level),
        "factor_M_exponent": str(factor),
        "factor_N_exponent": str(native_dual),
        "factor_L_exponent": str(factor),
        "complement_K_exponent": str(level),
        "native_poisson_length_exponent": str(native_dual),
        "balanced_mask_period_exponent": str(selected_prime),
        "pre_poisson_frequency_gap": str(pre_poisson_gap),
        "spectral_level_Q_exponent": str(spectral_level),
        "unmasked_off_diagonal_support_exponent": str(spectral_level),
        "masked_moving_center_support_exponent": str(mixed_support),
        "poisson_j_support_exponent": str(mixed_support),
        "post_cauchy_support_inflation": str(support_inflation),
        "native_additive_saving_benchmark": str(additive_saving_benchmark),
        "support_gap_after_native_additive_benchmark": str(
            residual_after_additive_benchmark
        ),
        "single_mask_dual_l1": "O(H*log(p1)*log(p2))",
        "paired_mask_dual_l1": "O(H*log(p1)^2*log(p2)^2)",
        "paired_mask_dual_l2_squared": "O(H)",
        "existing_additive_large_sieve_applies": False,
        "existing_watt_large_sieve_applies": False,
        "reason": (
            "the h-interval keeps length x^(2*level-1) but is centered at "
            "(a*k1-b*k2)/r and therefore depends on the same complementary "
            "variables that later enter the Kloosterman modulus"
        ),
    }


def square_moduli_large_sieve_ledger() -> dict[str, Any]:
    """Price the best generic square-moduli large sieve at Input VII."""
    prime = PAIR_PRIME_EXPONENT
    pair = PASCADI_LEVEL

    def delta(root_exponent: Fraction) -> Fraction:
        return max(
            3 * root_exponent,
            Fraction(1),
            min(
                2 * root_exponent + Fraction(1, 2),
                root_exponent / 2 + 1,
            ),
        )

    one_sided_delta = delta(prime)
    one_sided_error = (one_sided_delta + 1) / 2
    one_sided_conjectural_delta = max(3 * prime, Fraction(1))
    one_sided_conjectural_error = (one_sided_conjectural_delta + 1) / 2
    bilinear_delta = delta(pair)
    bilinear_error = (bilinear_delta + 1) / 2
    bilinear_conjectural_delta = max(3 * pair, Fraction(1))
    bilinear_conjectural_error = (bilinear_conjectural_delta + 1) / 2
    return {
        "one_sided_root_modulus_exponent": str(prime),
        "one_sided_best_known_delta_exponent": str(one_sided_delta),
        "one_sided_cauchy_error_exponent": str(one_sided_error),
        "one_sided_excess": str(one_sided_error - 1),
        "one_sided_zhao_conjectural_delta_exponent": str(
            one_sided_conjectural_delta
        ),
        "one_sided_zhao_conjectural_error_exponent": str(
            one_sided_conjectural_error
        ),
        "one_sided_projection_assumption": (
            "optimistic benchmark: treats the other selected-prime weight as "
            "a common arbitrary sequence, although it actually depends on p"
        ),
        "bilinear_root_modulus_exponent": str(pair),
        "bilinear_best_known_delta_exponent": str(bilinear_delta),
        "bilinear_cauchy_error_exponent": str(bilinear_error),
        "bilinear_excess": str(bilinear_error - 1),
        "bilinear_zhao_conjectural_delta_exponent": str(
            bilinear_conjectural_delta
        ),
        "bilinear_zhao_conjectural_error_exponent": str(
            bilinear_conjectural_error
        ),
        "lifted_fourier_l2_exponent": "0",
        "verdict": (
            "even the optimistic generic square-moduli Cauchy benchmark is "
            "dominated by the order-x trivial bound; Zhao's conjectural delta "
            "reaches only the one-sided boundary and leaves bilinear excess "
            "7/16"
        ),
    }


def unbalanced_p2_report() -> dict[str, Any]:
    """Exact dichotomy for replacing a selected prime by P2=r*s."""
    return {
        "prime_prime_route": (
            "if w+1=b*r*s with r>=s prime, the second coordinate has a "
            "short-cofactor prime exactly when r>2*b*s; standard Chen P2 "
            "does not impose this factor imbalance"
        ),
        "input_vii_route": (
            "if r>2*b*s the candidate already has a short-cofactor failure; "
            "otherwise r and s add two q-adic digit masks instead of one"
        ),
        "fully_nonzero_lift_exponent_after_replacement": "5/4",
        "mask_count_change": "two selected-prime masks become three",
        "nominal_mask_mass_change": "1/4 -> 1/8",
        "substitutes_for_input_vii": False,
        "verdict": (
            "Chen switching may address the already-closed prime-prime margin "
            "only with an additional unbalanced-P2 lower bound; it does not "
            "reduce the mixed q-adic estimate"
        ),
    }


def parseval_pair_bound() -> dict[str, Any]:
    """Exact limiting-exponent accounting after Fourier l2/Parseval."""
    prime = PAIR_PRIME_EXPONENT
    quotient = TRANSFER_SMOOTH_EXPONENT
    pair_count = 2 * prime

    # One-sided Parseval:
    # sum_h |S(h)|^2 <= p*N + N^2.  After the square root and summing
    # p1,p2, the N term has exponent pair_count+quotient = 1, while the
    # p*N term has exponent pair_count+(prime+quotient)/2 = 31/32.
    one_sided_diagonal = pair_count + quotient
    one_sided_off_diagonal = pair_count + (prime + quotient) / 2
    one_sided_exponent = max(one_sided_diagonal, one_sided_off_diagonal)

    # Bilinear frequencies are primitive mod M=p1*p2.  Since N=x/M<M,
    # orthogonality gives sum_a |S_M(a)|^2=M*N<=x.  Each pair therefore costs
    # at most x^(1/2), and there are x^(5/8+o(1)) pairs.
    bilinear_per_pair = (pair_count + quotient) / 2
    bilinear_exponent = pair_count + bilinear_per_pair

    trivial_absolute_exponent = pair_count + quotient
    effective_bilinear_exponent = min(
        bilinear_exponent, trivial_absolute_exponent
    )
    return {
        "prime_exponent": str(prime),
        "quotient_exponent": str(quotient),
        "pair_count_exponent": str(pair_count),
        "one_sided_diagonal_exponent": str(one_sided_diagonal),
        "one_sided_off_diagonal_exponent": str(one_sided_off_diagonal),
        "one_sided_parseval_exponent": str(one_sided_exponent),
        "bilinear_per_pair_exponent": str(bilinear_per_pair),
        "bilinear_parseval_exponent": str(bilinear_exponent),
        "trivial_absolute_exponent": str(trivial_absolute_exponent),
        "effective_bilinear_exponent": str(effective_bilinear_exponent),
        "main_term_exponent": "1",
        "one_sided_excess": str(one_sided_exponent - 1),
        "bilinear_parseval_excess": str(bilinear_exponent - 1),
        "effective_bilinear_excess": str(effective_bilinear_exponent - 1),
        "naive_individual_ds_excess": "19/32",
        "parseval_beats_trivial": bilinear_exponent < trivial_absolute_exponent,
        "status": (
            "FAILED APPROACH: Parseval reaches x^(9/8) on the bilinear "
            "frequencies and is dominated by the trivial order-x bound; "
            "neither gives the o(x) cancellation required by (67b)"
        ),
    }


def analytic_fit_report() -> dict[str, Any]:
    """Classify which frequency families fit the cited analytic theorems."""
    return {
        "pascadi_smooth_proposition": {
            "outer_modulus_weight": (
                "arbitrary scalar lambda_r is allowed; triply-well-factorable "
                "weights belong to Pascadi's prime theorem, not the smooth theorem"
            ),
            "zero_frequency": "APPLIES: exactly the base class n=-1 mod p1*p2",
            "nonzero_frequency": (
                "DOES_NOT_APPLY AS STATED: the mask depends jointly on the "
                "factorization r=p1*p2 and the quotient (n+1)/r, so it is "
                "neither a scalar outer weight nor an inner coefficient "
                "alpha_m beta_n gamma_l independent of r"
            ),
            "direct_one_sided_lift": (
                "FAILS proposition range: limiting modulus p1^2*p2 has "
                "exponent 15/16 > 5/8"
            ),
            "direct_bilinear_lift": (
                "FAILS even as a generic n-modulus: (p1*p2)^2 has limiting "
                "exponent 5/4 > 1"
            ),
        },
        "drappeau_shparlinski": {
            "one_sided": (
                "PARTIAL: individual additive sums at limiting modulus 15/16 "
                "save 1/32, but naive prime-pair summation loses 19/32"
            ),
            "bilinear": (
                "DOES NOT MATCH: in ell the phase is primitive mod p1*p2, "
                "but the weight is 1{P+(p1*p2*ell-1)<=y}, not a smooth-number "
                "weight in ell; in n the modulus is beyond the summation range"
            ),
        },
        "pascadi_completion": {
            "native_poisson_length_exponent": "1/4",
            "mask_period_exponent": "5/16",
            "pre_poisson_frequency_gap": "1/16",
            "masked_dual_support_exponent": "7/8",
            "unmasked_off_diagonal_support_exponent": "1/2",
            "post_cauchy_support_inflation": "3/8",
            "support_gap_after_native_additive_benchmark": "1/8",
            "after_cauchy_sequence": (
                "moving h-centers (a*k1-b*k2)/r depend on the same "
                "complementary variables that later enter the Kloosterman "
                "fraction, producing the mixed h*k*ell-h'*k'*ell' coefficient"
            ),
            "primary_source_verdict": (
                "OPEN: Pascadi explicitly states that no corresponding "
                "exceptional-spectrum large-sieve inequality is known for "
                "this mixed sequence"
            ),
            "post_2025_search": (
                "no unconditional resolution located; Baier arXiv:2503.18009 "
                "gives only conditional improvements for square moduli"
            ),
        },
        "irreducible_estimate": (
            "average the primitive bilinear quotient phases over p1,p2 before "
            "absolute values, with the shifted-smooth weight retained"
        ),
    }


def level_interface() -> dict[str, Any]:
    """Exact exponent gap for the first unresolved q-adic half-lift."""
    modulus_range_prime_exponent = PASCADI_LEVEL / 2
    largest_first_digit_modulus = 2 * TRANSFER_SMOOTH_EXPONENT
    return {
        "smooth_pair_exponent": str(TRANSFER_SMOOTH_EXPONENT),
        "pascadi_distribution_exponent": str(PASCADI_LEVEL),
        "modulus_range_reaches_q_squared_through": str(
            modulus_range_prime_exponent
        ),
        "largest_first_digit_modulus_exponent": str(
            largest_first_digit_modulus
        ),
        "distribution_gap_at_top": str(
            largest_first_digit_modulus - PASCADI_LEVEL
        ),
        "uncontrolled_prime_strip": [
            str(modulus_range_prime_exponent),
            str(TRANSFER_SMOOTH_EXPONENT),
        ],
        "uncontrolled_strip_width": str(
            TRANSFER_SMOOTH_EXPONENT - modulus_range_prime_exponent
        ),
        "half_lift_classes": (
            "H_q={-1+q*r mod q^2 : q/2<r<q}; count smooth n in H_q "
            "relative to smooth n congruent to -1 mod q"
        ),
        "centered_half_lift_fourier_l1": (
            "at most H_((q-1)/2) <= 1+log(q); the selected-prime double "
            "mask costs O(log(p1)*log(p2))"
        ),
        "yang_selected_prime_exponent": str(PAIR_PRIME_EXPONENT),
        "yang_quotient_exponent": str(TRANSFER_SMOOTH_EXPONENT),
        "quotient_over_prime_period_margin": str(
            TRANSFER_SMOOTH_EXPONENT - PAIR_PRIME_EXPONENT
        ),
        "limiting_one_sided_n_lift_modulus": str(
            3 * PAIR_PRIME_EXPONENT
        ),
        "limiting_one_sided_lift_beyond_pascadi": str(
            3 * PAIR_PRIME_EXPONENT - PASCADI_LEVEL
        ),
        "limiting_one_sided_ds_saving": "1/32",
        "limiting_one_sided_sum_exponent": "31/32",
        "limiting_one_sided_naive_pair_sum_exponent": "51/32",
        "limiting_one_sided_naive_excess": "19/32",
        "limiting_bilinear_quotient_modulus": str(
            2 * PAIR_PRIME_EXPONENT
        ),
        "limiting_bilinear_n_lift_modulus": str(
            4 * PAIR_PRIME_EXPONENT
        ),
        "bilinear_n_lift_excess_over_length": "1/4",
        "double_half_lift_mask": (
            "D(p1,p2,l)=1{2*(p2*l mod p1)>p1}*"
            "1{2*(p1*l mod p2)>p2}; expected unweighted density 1/4"
        ),
        "verdict": (
            "OPEN structured half-lift dispersion. Zero frequency is Pascadi's "
            "base progression. One-sided frequencies lift n to limiting "
            "modulus 15/16. Frequencies nonzero on both sides are primitive "
            "mod p1*p2 in ell and lift n to limiting modulus 5/4, beyond the "
            "n-range; a proof must preserve the factorized ell variable and "
            "average over both selected primes before absolute values"
        ),
    }





def k_factor_construction_report(selected_factor_count: int) -> dict[str, Any]:
    """Price splitting the 5/8 selected modulus among k equal prime factors."""
    _require_int("selected_factor_count", selected_factor_count, 2)
    count = selected_factor_count
    base = PASCADI_LEVEL
    prime = base / count
    quotient = 1 - base
    one_sided_lift = base + prime
    one_sided_saving = (1 - one_sided_lift) / 2
    one_sided_naive_excess = base - one_sided_saving
    return {
        "selected_factor_count": count,
        "base_modulus_exponent": str(base),
        "each_prime_exponent": str(prime),
        "quotient_exponent": str(quotient),
        "one_sided_lift_exponent": str(one_sided_lift),
        "one_sided_individual_saving": str(one_sided_saving),
        "one_sided_naive_excess": str(one_sided_naive_excess),
        "limiting_excess_as_k_grows": "7/16",
        "nominal_independent_success_mass": str(Fraction(1, 2**count)),
        "fully_bilinear_quotient_modulus": "5/8",
        "fully_bilinear_n_lift_modulus": "5/4",
        "closes": False,
    }


def automatic_compensation_report() -> dict[str, Any]:
    """Proved obstructions to the two lower-complexity construction changes."""
    return {
        "minimum_congruence_modulus": {
            "prime_adic_exponent": 2,
            "proof": (
                "inside w=-1 mod p, the first digit is (w+1)/p mod p; "
                "any progression whose modulus has v_p<=1 cycles through all "
                "p lifts and contains both successes and failures"
            ),
        },
        "size_forced": {
            "condition": "w+1=p*c with p/2<c<p",
            "consequence": "p^2/2<w+1<p^2, so p is confined to sqrt(w)",
            "candidate_order": "x/log(x)",
            "positive_density": False,
        },
        "factor_splitting": [
            k_factor_construction_report(count) for count in (2, 3, 4, 8)
        ],
        "near_diagonal_cut": {
            "gap_cutoff": "p2-p1<=x^(5/16-kappa)",
            "discarded_representation_exponent": "1+3*delta-kappa",
            "negligible_when": "kappa>3*delta",
            "purpose": (
                "removes bounded-gap and near-resonant prime pairs without "
                "changing the positive main density"
            ),
        },
        "status": (
            "FAILED APPROACH: no lower-p-adic congruence or fixed-k factor "
            "split makes compensation automatic at positive density; this "
            "does not exclude a different algebraic construction"
        ),
    }


def remaining_tail_report(
    epsilon: Fraction = Fraction(1, 200),
) -> dict[str, Any]:
    """Exact exponent ledger after the selected-prime mask is certified."""
    parameters = transfer_parameters(epsilon)
    p1_lower = Fraction(parameters["p1_exponents"][0])
    p2_lower = Fraction(parameters["p2_exponents"][0])
    pair_upper = Fraction(parameters["modulus_exponents"][1])
    quotient_upper = Fraction(parameters["cofactor_upper_exponent"])
    selected_repeat_exponent = max(
        1 - p1_lower,
        1 - p2_lower,
        pair_upper,
    )
    large_repeat_floor_error = pair_upper + quotient_upper / 2
    return {
        "selected_prime_repetition_exponent": str(selected_repeat_exponent),
        "limiting_selected_prime_repetition_exponent": "11/16",
        "selected_prime_repetition_is_negligible": selected_repeat_exponent < 1,
        "large_repeated_factor_main_exponent": "1-zeta",
        "large_repeated_factor_floor_error_exponent": str(
            large_repeat_floor_error
        ),
        "large_repeated_factor_floor_error_is_negligible": (
            large_repeat_floor_error < 1
        ),
        "left_coordinate": (
            "P+(n)<=x^(1/C); all remaining primes are ultra-small but their "
            "q-adic prefix automata still require a uniform smooth-weighted tail"
        ),
        "quotient_coordinate": (
            "for q^e||ell, compensation tests "
            "P=(p1*p2)*(ell/q^e) modulo q^r; exponent-one q is another "
            "mixed half-lift, while repeated q>x^zeta are o(x)"
        ),
        "status": (
            "PROVED selected-prime repetitions and all repeated factors above "
            "x^zeta are negligible; OPEN small repeated factors and "
            "exponent-one quotient factors"
        ),
    }


def measure_constructed_probe(
    *, start: int = 10**6, count: int = 4_000_000
) -> dict[str, Any]:
    """Finite friendly-width probe of the structured mask in (67b).

    The theorem lets delta be arbitrarily small.  A finite block cannot resolve
    those asymptotically narrow prime ranges, so this probe uses delta=1/50 and
    P+(n)<=top^(1/4).  It measures representations, not distinct n, and is
    explicitly experimental.
    """
    _require_int("start", start, 2)
    _require_int("count", count, 1)
    top = start + count
    delta = Fraction(1, 50)
    eta = PAIR_PRIME_EXPONENT
    exponents = (
        eta - 3 * delta,
        eta - 2 * delta,
        eta - delta,
    )
    bounds = [
        rational_power_floor(top, exponent.numerator, exponent.denominator)
        for exponent in exponents
    ]
    left_smooth_bound = rational_power_floor(top, 1, 4)
    factorizer = SPFFactorizer(top + 1)
    representations = 0
    selected_prime_passes = 0
    fully_good_representations = 0
    selected_and_fully_good = 0
    selected_and_left_good = 0
    selected_and_right_good = 0
    for value in range(start, top):
        left_factorization = factorizer.factor(value)
        if max(left_factorization, default=1) > left_smooth_bound:
            continue
        right_factorization = factorizer.factor(value + 1)
        first_primes = [
            prime
            for prime in right_factorization
            if bounds[0] < prime <= bounds[1]
        ]
        second_primes = [
            prime
            for prime in right_factorization
            if bounds[1] < prime <= bounds[2]
        ]
        if not first_primes or not second_primes:
            continue
        left_short, left_level = classify_term(value, 1, factorizer)
        right_short, right_level = classify_term(value + 1, 1, factorizer)
        left_good = not (left_short or left_level)
        right_good = not (right_short or right_level)
        fully_good = left_good and right_good
        for first_prime in first_primes:
            for second_prime in second_primes:
                quotient = (value + 1) // (first_prime * second_prime)
                selected_pass = (
                    2 * ((second_prime * quotient) % first_prime) > first_prime
                    and 2 * ((first_prime * quotient) % second_prime) > second_prime
                )
                representations += 1
                selected_prime_passes += selected_pass
                fully_good_representations += fully_good
                selected_and_fully_good += selected_pass and fully_good
                selected_and_left_good += selected_pass and left_good
                selected_and_right_good += selected_pass and right_good
    if representations == 0:
        raise AssertionError("constructed probe has no representations")
    return {
        "status": "FINITE_EXPERIMENT_NOT_THEOREM_PARAMETERS",
        "start": start,
        "count": count,
        "delta": str(delta),
        "left_smooth_exponent": "1/4",
        "prime_exponents": [str(exponent) for exponent in exponents],
        "prime_bounds": bounds,
        "representations": representations,
        "selected_prime_passes": selected_prime_passes,
        "selected_prime_pass_rate": selected_prime_passes / representations,
        "expected_unweighted_pass_rate": 0.25,
        "fully_good_representations": fully_good_representations,
        "selected_and_fully_good": selected_and_fully_good,
        "selected_and_left_good": selected_and_left_good,
        "selected_and_right_good": selected_and_right_good,
    }




def analyze(*, start: int, count: int, target_m: int) -> dict[str, Any]:
    started = time.perf_counter()
    parameters = transfer_parameters()
    block = measure_block(start=start, count=count, target_m=target_m)
    return {
        "theorem": {
            "status": "PROVED_FROM_CITED_PASCADI_THEOREM_1_5",
            "statement": (
                "for every epsilon>0, a positive density of n satisfies "
                "P+(n)<P+(n+1)<x^(3/8+epsilon)"
            ),
            "yang_arxiv": "2607.16032",
            "yang_theorem": "Theorem 1.5 (41/107 + epsilon)",
            "pascadi_arxiv": "2505.00653v2",
            "pascadi_theorem": "Theorem 1.5 (absolute smooth BV to 5/8 - epsilon)",
            "exponent_identity": "1 - 5/8 = 3/8",
        },
        "parameters": parameters,
        "level_interface": level_interface(),
        "frequency_decomposition": double_mask_frequency_classes(101, 127),
        "q_adic_lift": q_adic_lift_report(101 * 127, 1),
        "analytic_fit": analytic_fit_report(),
        "parseval_bound": parseval_pair_bound(),
        "q_adic_poisson": q_adic_poisson_reindex_report(15, 2, 60),
        "mixed_dispersion": mixed_dispersion_exponents(),
        "square_moduli_large_sieve": square_moduli_large_sieve_ledger(),
        "unbalanced_p2": unbalanced_p2_report(),
        "remaining_tail": remaining_tail_report(),
        "automatic_compensation": automatic_compensation_report(),
        "finite_evidence": block,
        "constructed_probe": measure_constructed_probe(),
        "status": "PASS_DENSE_SMOOTH_PAIR_TRANSFER",
        "implementation_sha256": {
            "smooth_pair_transfer.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "good_density_decomposition.py": sha256_file(
                ROOT / "good_density_decomposition.py"
            ),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "resources": {"wall_seconds": time.perf_counter() - started},
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=10**6)
    parser.add_argument("--count", type=int, default=500_000)
    parser.add_argument("--target-m", type=int, default=1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = analyze(start=args.start, count=args.count, target_m=args.target_m)
    atomic_write_json(args.output, report)
    block = report["finite_evidence"]
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(args.output),
                "smooth_pairs": block["smooth_pairs"],
                "smooth_good_pairs": block["smooth_good_pairs"],
                "conditional_level_survival": block[
                    "conditional_good_given_smooth_pair"
                ],
                "wall_seconds": report["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
