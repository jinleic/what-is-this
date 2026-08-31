#!/usr/bin/env python3
"""Locate exactly where the elementary route to Conjecture R(L) dies.

Section 17 proved R(1) by buying back the union-bound overshoot ``u - 1`` with
a lower bound on the Bonferroni overlap, where ``u = log 2 + Lambda_B``. The
obvious next question is R(2). This producer prices the next rungs and records
the boundaries of the methods tried here.

Four exact facts are certified here.

1.  THE SMOOTHNESS COLLAPSE.  When ``m >= sqrt(2w)`` every prime ``p > m``
    dividing ``w`` has ``p^2 > 2w``, so by (16a) it is a short-cofactor failure
    outright and no level test survives.  Hence in that regime

        w is m-compensation-good   <==>   P^+(w) <= m,

    an equivalence with no digit content at all.  It is checked exhaustively,
    not sampled.

2.  THE UNION-BOUND THRESHOLD.  In the collapsed regime the failure density of
    one term tends to ``log(1/alpha)`` by Mertens, where
    ``m = (2X)^alpha`` and ``alpha > 1/2``.  A run of ``L`` terms therefore has
    density at least ``1 - L log(1/alpha)``, whose positivity threshold is
    ``alpha > e^(-1/L)``.  For ``L = 2`` this is
    ``e^(-1/2) = 0.606531``.  The theorem is about growing ``m``; it cannot
    reach Erdos #389, where ``m`` is fixed and ``alpha -> 0``.

3.  THE DEFICIT JUMP.  At ``L = 1`` the first union bound misses by ``u - 1``.
    At ``L = 2`` it misses by ``2u - 1``, larger by a factor of about 66.
    Section 17 supplied ``0.045`` of diagonal overlap against a demand of
    ``0.0157``; the cross-coordinate terms are indispensable at ``L = 2``.

4.  THE MOMENT COUNTERMODEL.  Exact rational joint laws have marginal
    ``31/100`` at every coordinate but zero mass on a full run, for
    ``L = 2, 3, 4``.  At ``L = 2`` the second moment is algebraically identical
    to the missing pair correlation, so it is not an independent route.

The producer also records the two routes that were tried and failed, each with
the number that kills it:

  parity          the density route needs d(A_0 and A_1) > 2 log 2 - 1; that
                  count is #{pa + 1 = p'b} with p, p' > sqrt(2x), i.e. moduli
                  pp' > 2x strictly beyond the range - a binary problem.
  thin family     the shifted squares (x^2 - 1, x^2) make BOTH coordinates
                  short-cofactor-free unconditionally (Lemma 18.5, elementary).
                  For a fixed band r the number of parameters carrying q is
                  smaller than the required q^r residue period by at least
                  X^(r/(r+2)).  The local q-test is deterministic for every
                  prime q in (X, sqrt2 X]: x = q + 1 fails it and x = q - 1
                  passes it, with no exceptions.
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
from good_density_decomposition import compensation_count

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "run_threshold.json"

#: sum_{r>=1} 2^-r log((r+2)/(r+1)), the level-union constant of Lemma 17.3.
LAMBDA_B = sum(2.0**-r * math.log((r + 2) / (r + 1)) for r in range(1, 400))
#: u = log 2 + Lambda_B, the first-order union-bound total per coordinate.
U_TOTAL = math.log(2.0) + LAMBDA_B


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def is_good(value: int, factorisation: dict[int, int], target_m: int) -> bool:
    """(18): ``C_p(w) >= v_p(w)`` for every prime ``p > m`` dividing ``w``."""
    for prime, exponent in factorisation.items():
        if prime <= target_m:
            continue
        if compensation_count(value // prime**exponent, prime) < exponent:
            return False
    return True


def largest_prime_factor(factorisation: dict[int, int]) -> int:
    return max(factorisation) if factorisation else 1


def smoothness_collapse(start: int, span: int, alpha: float) -> dict[str, Any]:
    """Fact 1: for ``m >= sqrt(2w)``, m-goodness is exactly m-smoothness.

    Exhaustive over the declared block; any mismatch is an assertion failure,
    never a reported statistic.
    """
    _require_int("start", start, 2)
    _require_int("span", span, 1)
    if not isinstance(alpha, float) or not 0.5 <= alpha < 1.0:
        raise ValueError("alpha must be a float in [0.5, 1)")
    target_m = int((2 * start) ** alpha)
    top = start + span
    if target_m * target_m < 2 * top:
        raise ValueError("alpha too small: the collapse hypothesis m >= sqrt(2w) fails")
    factorizer = SPFFactorizer(top + 1)
    checked = 0
    for value in range(start, top):
        factorisation = factorizer.factor(value)
        smooth = largest_prime_factor(factorisation) <= target_m
        if smooth != is_good(value, factorisation, target_m):
            raise AssertionError(f"collapse fails at w={value}, m={target_m}")
        checked += 1
    return {
        "start": start,
        "span": span,
        "alpha": alpha,
        "m": target_m,
        "hypothesis_m_sq_ge_2top": True,
        "checked": checked,
        "mismatches": 0,
    }


def run_density(start: int, span: int, alpha: float, length: int) -> dict[str, Any]:
    """Fact 2: measured density of ``L``-runs against ``1 - L log(1/alpha)``."""
    _require_int("length", length, 1)
    target_m = int((2 * start) ** alpha)
    top = start + span + length
    if target_m * target_m < 2 * top:
        raise ValueError("alpha too small: the collapse hypothesis m >= sqrt(2w) fails")
    factorizer = SPFFactorizer(top + 1)
    smooth = [
        largest_prime_factor(factorizer.factor(value)) <= target_m
        for value in range(start, top)
    ]
    runs = sum(all(smooth[i : i + length]) for i in range(span))
    bound = 1.0 - length * math.log(1.0 / alpha)
    measured = runs / span
    if measured < bound:
        raise AssertionError(f"measured run density {measured} below proved bound {bound}")
    return {
        "alpha": alpha,
        "length": length,
        "m": target_m,
        "measured_run_density": measured,
        "proved_lower_bound": bound,
        "threshold_alpha": math.exp(-1.0 / length),
        "bound_is_positive": bound > 0.0,
    }


def deficit_ladder(max_length: int) -> dict[str, Any]:
    """Fact 3: what the union bound misses at each run length, and the cap."""
    _require_int("max_length", max_length, 1)
    rungs = []
    for length in range(1, max_length + 1):
        rungs.append(
            {
                "length": length,
                "union_total": length * U_TOTAL,
                "deficit": length * U_TOTAL - 1.0,
                "bonferroni_mass_cap": 2.0 * U_TOTAL if length >= 2 else U_TOTAL,
            }
        )
    return {
        "u": U_TOTAL,
        "rungs": rungs,
        "jump_factor_L1_to_L2": (2 * U_TOTAL - 1.0) / (U_TOTAL - 1.0),
    }


def moment_obstruction(max_length: int) -> dict[str, Any]:
    """Exact countermodels: marginals below 1-1/L do not force an L-run.

    The proved density ceiling is ``1-log(2) < 0.31``.  For each L >= 2 put
    mass ``g/(L-1)`` on each of the L subsets missing one coordinate and the
    remaining mass on the empty subset, with ``g = 31/100``.  Every coordinate
    then has marginal exactly g but the all-ones event has mass zero.  Fractions
    are returned as strings so the certificate is exact, not floating point.
    """
    _require_int("max_length", max_length, 2)
    marginal = Fraction(31, 100)
    if not 1.0 - math.log(2.0) < float(marginal):
        raise AssertionError("31/100 no longer upper-bounds the proved ceiling")
    rows = []
    for length in range(2, max_length + 1):
        subset_weight = marginal / (length - 1)
        empty_weight = 1 - length * subset_weight
        reconstructed = (length - 1) * subset_weight
        if empty_weight < 0 or reconstructed != marginal:
            raise AssertionError(f"countermodel failed at L={length}")
        rows.append(
            {
                "length": length,
                "marginal": str(marginal),
                "weight_each_subset_missing_one": str(subset_weight),
                "empty_weight": str(empty_weight),
                "reconstructed_marginal": str(reconstructed),
                "all_ones_weight": "0",
                "second_moment_identity_L2": (
                    "E[S^2] = E[S] + 2 E[G0 G1]" if length == 2 else None
                ),
            }
        )
    return {
        "proved_density_ceiling": 1.0 - math.log(2.0),
        "rational_upper_bound": str(marginal),
        "rows": rows,
        "conclusion": (
            "marginals alone cannot force a full run; at L=2 the second "
            "moment is algebraically identical to the missing pair correlation"
        ),
    }



def _log_interval_integer(value: int, terms: int = 24) -> tuple[Fraction, Fraction]:
    """Exact rational enclosure of ``log(value)``.

    Reduce to ``value / 2^a in [1, 2)`` and use
    ``log y = 2 atanh((y-1)/(y+1))``.  Both resulting atanh arguments are at
    most ``1/3``, and the omitted positive tail is bounded geometrically.
    """
    _require_int("value", value, 1)
    _require_int("terms", terms, 1)

    def atanh_interval(argument: Fraction) -> tuple[Fraction, Fraction]:
        total = Fraction(0)
        power = argument
        for index in range(terms):
            total += power / (2 * index + 1)
            power *= argument * argument
        lower = 2 * total
        tail = 2 * power / ((2 * terms + 1) * (1 - argument * argument))
        return lower, lower + tail

    exponent = value.bit_length() - 1
    power_of_two = 1 << exponent
    log_two_lower, log_two_upper = atanh_interval(Fraction(1, 3))
    reduced_lower, reduced_upper = atanh_interval(
        Fraction(value - power_of_two, value + power_of_two)
    )
    return (
        exponent * log_two_lower + reduced_lower,
        exponent * log_two_upper + reduced_upper,
    )


def _exp_reciprocal_interval(
    denominator: int, terms: int = 20
) -> tuple[Fraction, Fraction]:
    """Exact rational enclosure of ``exp(1 / denominator)``."""
    _require_int("denominator", denominator, 1)
    _require_int("terms", terms, 1)
    argument = Fraction(1, denominator)
    term = Fraction(1)
    total = term
    for index in range(1, terms + 1):
        term *= argument / index
        total += term
    next_term = term * argument / (terms + 1)
    tail = next_term / (1 - argument / (terms + 2))
    return total, total + tail


def _minimum_size_admissible_k(target_m: int) -> int:
    """Smallest k not excluded by ``binomial ratio < 2``, with exact bounds."""
    log_two_lower, log_two_upper = _log_interval_integer(2)
    candidate = 1
    while Fraction(candidate * candidate, target_m + 1) <= log_two_upper:
        candidate += 1
    previous = Fraction((candidate - 1) ** 2, target_m + 1)
    if previous > log_two_lower:
        raise AssertionError("log(2) enclosure too wide to decide the size floor")
    return candidate


def large_m_at_once_certificate(finite_max_m: int = 2047) -> dict[str, Any]:
    """The positive first-order run regime contains no possible witness.

    A witness has binomial ratio at least two, forcing
    ``k^2 / (m+1) > log 2``.  The positive regime of Theorem 18.2, applied at
    the first bad-window term ``X = k + floor(m/2) + 1``, instead requires
    ``log(2X) < exp(1/L) log m``, ``L = ceil(m/2)``.  Exact rational intervals
    prove these incompatible for every ``3 <= m <= finite_max_m``.

    Section 18.8 gives the elementary continuation for all ``m >= 2048``.
    """
    _require_int("finite_max_m", finite_max_m, 3)
    minimum_margin: Fraction | None = None
    minimum_margin_m = 0
    maximum_size_floor = 0
    for target_m in range(3, finite_max_m + 1):
        length = (target_m + 1) // 2
        candidate_k = _minimum_size_admissible_k(target_m)
        first_term = candidate_k + target_m // 2 + 1
        log_scale_lower, _ = _log_interval_integer(2 * first_term)
        _, log_m_upper = _log_interval_integer(target_m)
        _, exp_upper = _exp_reciprocal_interval(length)
        margin = log_scale_lower - exp_upper * log_m_upper
        if margin <= 0:
            raise AssertionError(
                f"large-m route not excluded at m={target_m}, k={candidate_k}"
            )
        if minimum_margin is None or margin < minimum_margin:
            minimum_margin = margin
            minimum_margin_m = target_m
        maximum_size_floor = max(maximum_size_floor, candidate_k)

    assert minimum_margin is not None
    log_two_lower, _ = _log_interval_integer(2)
    _, log_2048_upper = _log_interval_integer(2048)
    analytic_base = 2048 * log_two_lower > 16 * log_2048_upper * log_2048_upper
    if not analytic_base:
        raise AssertionError("analytic continuation base inequality failed")
    return {
        "finite_range": [3, finite_max_m],
        "checked_m": finite_max_m - 2,
        "minimum_log_margin": str(minimum_margin),
        "minimum_log_margin_decimal": float(minimum_margin),
        "minimum_margin_at_m": minimum_margin_m,
        "maximum_size_floor": maximum_size_floor,
        "analytic_start_m": 2048,
        "analytic_base_inequality": analytic_base,
        "conclusion": (
            "the parameter range where the first-order run bound is positive "
            "contains no k large enough for binomial divisibility"
        ),
    }

def parity_requirement() -> dict[str, Any]:
    """The density route: what a lower bound on d(A_0 and A_1) must beat."""
    trivial = 2 * math.log(2.0) - 1.0
    independent = math.log(2.0) ** 2
    return {
        "bonferroni_floor": trivial,
        "independence_prediction": independent,
        "required_gap_over_floor": independent - trivial,
        "moduli_exceed_range": True,
        "reason": (
            "p > sqrt(2x) and p' > sqrt(2(x+1)) force pp' > 2x, so the "
            "congruence defining the cross term has modulus strictly larger "
            "than the range; the count is #{pa + 1 = p'b} with "
            "a < sqrt(x/2), b < sqrt((x+1)/2), a binary problem of Chen type"
        ),
    }


def shifted_square_family(start: int, span: int, target_m: int) -> dict[str, Any]:
    """Measure a short-cofactor-free family lacking residue-density control.

    Lemma 18.5 is verified exhaustively (never sampled): for every x >= 3 the
    largest prime factor of x^2 - 1 is at most x + 1 <= sqrt(2(x^2 - 1)), and
    x^2 is a square, so neither coordinate can carry a short cofactor.
    """
    _require_int("start", start, 3)
    _require_int("span", span, 1)
    _require_int("target_m", target_m, 1)
    top = start + span
    factorizer = SPFFactorizer(top + 2)
    low_good = square_good = both = 0
    union_terms = 0
    for x in range(start, top):
        value = x * x - 1
        factorisation = dict(factorizer.factor(x - 1))
        for prime, exponent in factorizer.factor(x + 1).items():
            factorisation[prime] = factorisation.get(prime, 0) + exponent
        # Lemma 18.5, checked rather than assumed.
        if largest_prime_factor(factorisation) ** 2 > 2 * value:
            raise AssertionError(f"x^2-1 carries a short cofactor at x={x}")
        ok_low = True
        for prime, exponent in factorisation.items():
            if prime <= target_m:
                continue
            if compensation_count(value // prime**exponent, prime) < exponent:
                ok_low = False
                union_terms += 1
        square_factorisation = {p: 2 * e for p, e in factorizer.factor(x).items()}
        if largest_prime_factor(square_factorisation) ** 2 > 2 * x * x:
            raise AssertionError(f"x^2 carries a short cofactor at x={x}")
        ok_square = is_good(x * x, square_factorisation, target_m)
        low_good += ok_low
        square_good += ok_square
        both += ok_low and ok_square
    return {
        "start": start,
        "span": span,
        "target_m": target_m,
        "short_cofactor_free_both_coordinates": True,
        "density_x2_minus_1_good": low_good / span,
        "density_x2_good": square_good / span,
        "density_both_good": both / span,
        "union_terms_per_x": union_terms / span,
        "free_interval_level_constant_2_lambda_b": 2 * LAMBDA_B,
        "measured_union_exceeds_free_interval_constant": union_terms / span > 2 * LAMBDA_B,
    }


def deterministic_bias(scale: int) -> dict[str, Any]:
    """Why Section 17 cannot count the thin family by residue density.

    For x in [X, 2X) and a prime q in (X, sqrt2 X], the only parameters with
    q dividing x^2 - 1 are x = q +- 1.  Then (x^2 - 1)/q is q + 2 or q - 2,
    so the local q-test is decided outright: x = q + 1 fails it and
    x = q - 1 passes it.  A residue heuristic would predict one half each.
    """
    _require_int("scale", scale, 1000)
    top = int(math.isqrt(2 * scale * scale)) + 1
    factorizer = SPFFactorizer(top + 2)
    primes = [
        q
        for q in range(scale + 1, top)
        if factorizer.factor(q) == {q: 1}
    ]
    plus_fail = sum(
        1 for q in primes if compensation_count(((q + 1) ** 2 - 1) // q, q) < 1
    )
    minus_fail = sum(
        1 for q in primes if compensation_count(((q - 1) ** 2 - 1) // q, q) < 1
    )
    return {
        "scale": scale,
        "primes_tested": len(primes),
        "x_eq_q_plus_1_failures": plus_fail,
        "x_eq_q_minus_1_failures": minus_fail,
        "heuristic_failure_rate": 0.5,
        "observed_plus_rate": plus_fail / len(primes) if primes else None,
        "observed_minus_rate": minus_fail / len(primes) if primes else None,
        "deterministic": plus_fail == len(primes) and minus_fail == 0,
    }


def analyze(
    *,
    collapse_start: int,
    collapse_span: int,
    run_start: int,
    run_span: int,
    family_start: int,
    family_span: int,
    bias_scale: int,
    max_length: int,
) -> dict[str, Any]:
    started = time.perf_counter()
    collapse = [
        smoothness_collapse(collapse_start, collapse_span, alpha)
        for alpha in (0.55, 0.70, 0.85)
    ]
    runs = [
        run_density(run_start, run_span, alpha, length)
        for length, alpha in ((1, 0.55), (2, 0.75), (2, 0.62), (3, 0.80))
    ]
    family = shifted_square_family(family_start, family_span, 1)
    return {
        "smoothness_collapse": collapse,
        "run_threshold": runs,
        "deficit_ladder": deficit_ladder(max_length),
        "parity_requirement": parity_requirement(),
        "moment_obstruction": moment_obstruction(max_length),
        "large_m_at_once": large_m_at_once_certificate(),
        "shifted_square_family": family,
        "deterministic_bias": deterministic_bias(bias_scale),
        "constants": {"lambda_b": LAMBDA_B, "u": U_TOTAL},
        "status": "PASS_EXACT_RUN_THRESHOLD",
        "implementation_sha256": {
            "run_threshold.py": sha256_file(Path(__file__)),
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
    parser.add_argument("--collapse-start", type=int, default=10**6)
    parser.add_argument("--collapse-span", type=int, default=40_000)
    parser.add_argument("--run-start", type=int, default=10**6)
    parser.add_argument("--run-span", type=int, default=60_000)
    parser.add_argument("--family-start", type=int, default=10**6)
    parser.add_argument("--family-span", type=int, default=30_000)
    parser.add_argument("--bias-scale", type=int, default=10**6)
    parser.add_argument("--max-length", type=int, default=4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = analyze(
        collapse_start=args.collapse_start,
        collapse_span=args.collapse_span,
        run_start=args.run_start,
        run_span=args.run_span,
        family_start=args.family_start,
        family_span=args.family_span,
        bias_scale=args.bias_scale,
        max_length=args.max_length,
    )
    atomic_write_json(args.output, report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(args.output),
                "threshold_L2": math.exp(-0.5),
                "jump_factor": report["deficit_ladder"]["jump_factor_L1_to_L2"],
                "family_both_good": report["shifted_square_family"]["density_both_good"],
                "bias_deterministic": report["deterministic_bias"]["deterministic"],
                "wall_seconds": report["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
