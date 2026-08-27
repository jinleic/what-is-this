#!/usr/bin/env python3
"""Certify the local large-prime structure behind Erdős #389 repairs.

For ``p > m``, every negative slack belongs to the unique bad-window integer
that ``p`` divides.  This module turns that fact into an exact term predicate,
checks its separation from the small-prime tier, and quantifies why retaining
all currently obstructing primes forces enormous CRT displacements.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from math import gcd
from pathlib import Path
from typing import Any, Iterable

from artifact_io import atomic_write_json
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    SPFFactorizer,
    TrialFactorizer,
    carry_count,
    integer_valuation,
    primes_up_to,
    slack,
    witness_certificate,
)
from residue_control_analysis import (
    crt_coprime,
    least_power_above,
    natural_shift_spike_formula,
    upper_half_prefix_count,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_MOVING = ROOT / "data" / "moving_bad_window_m27_h100000000.json"
DEFAULT_EXTENDED_MOVING = (
    ROOT / "data" / "moving_bad_window_m27_h100000000_199999999.json"
)
DEFAULT_NEAR_MISS = (
    ROOT / "data" / "bad_window_near_miss_m27_h100000000.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "data"
    / "compensation_structure_m1_20_k5000_m27_h100000000.json"
)

Factorizer = TrialFactorizer | SPFFactorizer


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _term_formula_from_exponent(term: int, p: int, exponent: int) -> dict[str, int]:
    power = p**exponent
    cofactor = term // power
    positive_levels = upper_half_prefix_count(cofactor, p)
    return {
        "term": term,
        "prime": p,
        "exponent": exponent,
        "cofactor": cofactor,
        "positive_prefix_levels": positive_levels,
        "predicted_slack": positive_levels - exponent,
    }


def term_compensation_formula(term: int, p: int) -> dict[str, int]:
    """Return the exact ``p > m`` slack attached to one bad-window term.

    Primality is supplied by callers' exact factorizations.  The returned
    formula itself is independent of ``m`` and of the term's bad-window
    position; those values only determine whether the formula applies.
    """
    _require_int("term", term, 1)
    _require_int("p", p, 2)
    if term % p:
        raise ValueError("p must divide term")
    exponent = integer_valuation(term, p)
    return _term_formula_from_exponent(term, p, exponent)

def short_cofactor_compensation_formula(term: int, p: int) -> dict[str, int]:
    """Specialize the exact term formula when ``term / p**e < p``."""
    record = term_compensation_formula(term, p)
    cofactor = record["cofactor"]
    if cofactor >= p:
        raise ValueError("require the p-free cofactor to be smaller than p")
    indicator = int(2 * cofactor > p)
    if record["positive_prefix_levels"] != indicator:
        raise AssertionError("short-cofactor prefix count formula failed")
    return {
        **record,
        "upper_half_indicator": indicator,
    }


def large_prime_window_records(
    m: int,
    k: int,
    factorizer: Factorizer,
) -> tuple[dict[str, int], ...]:
    """Return every ``p > m`` term formula in the bad window of ``(m, k)``."""
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    first_position = m // 2 + 1
    records: list[dict[str, int]] = []
    seen: set[int] = set()
    for position in range(first_position, m + 1):
        term = k + position
        for p, exponent in factorizer.factor(term).items():
            if p <= m:
                continue
            if p in seen:
                raise AssertionError("a prime larger than m hit two bad-window terms")
            seen.add(p)
            records.append(
                {
                    "position": position,
                    **_term_formula_from_exponent(term, p, exponent),
                }
            )
    return tuple(sorted(records, key=lambda record: record["prime"]))


def large_prime_window_obstructions(
    m: int,
    k: int,
    factorizer: Factorizer,
) -> tuple[dict[str, int], ...]:
    """Return exactly the negative-slack primes larger than ``m``."""
    return tuple(
        record
        for record in large_prime_window_records(m, k, factorizer)
        if record["predicted_slack"] < 0
    )


def flexible_prefix_residue(p: int, exponent: int, residue: int) -> bool:
    """Whether an odd-prime cofactor residue has ``exponent`` safe prefixes.

    If this holds and ``w = p**exponent * c`` with the given residue for
    ``c modulo p**exponent``, the first ``exponent`` positive levels cancel
    the ``exponent`` negative levels at ``p``.  This is the full prefix cone
    containing the single mirror residue ``c = -1``.
    """
    _require_int("p", p, 3)
    _require_int("exponent", exponent, 1)
    _require_int("residue", residue, 0)
    if p % 2 == 0:
        raise ValueError("the strict first-prefix cone requires odd p")
    modulus = p**exponent
    residue %= modulus
    power = p
    for _ in range(exponent):
        if 2 * (residue % power) <= power:
            return False
        power *= p
    return True


def flexible_prefix_residue_count(p: int, exponent: int) -> int:
    """Count the odd-prime safe prefix cone modulo ``p**exponent`` exactly."""
    _require_int("p", p, 3)
    _require_int("exponent", exponent, 1)
    if p % 2 == 0:
        raise ValueError("the strict first-prefix cone requires odd p")
    return ((p - 1) // 2) * ((p + 1) // 2) ** (exponent - 1)


def retained_prime_displacement(
    k: int,
    prime_positions: Iterable[tuple[int, int]],
) -> dict[str, int]:
    """Return the exact forward spacing forced by retaining attached primes.

    Each pair is ``(p, i)`` and must satisfy ``p | k+i``.  Pairwise-coprime
    primes force every candidate retaining all attachments to be congruent to
    ``k`` modulo their product.
    """
    _require_int("k", k, 1)
    controls = tuple(prime_positions)
    if not controls:
        raise ValueError("at least one retained prime is required")
    product = 1
    seen: set[int] = set()
    for p, position in controls:
        _require_int("p", p, 2)
        _require_int("position", position, 1)
        if p in seen or gcd(product, p) != 1:
            raise ValueError("retained primes must be distinct and coprime")
        if (k + position) % p:
            raise ValueError("each retained prime must divide its attached term")
        seen.add(p)
        product *= p
    return {
        "prime_count": len(controls),
        "prime_product": product,
        "least_forward_delta": product,
        "least_forward_candidate": k + product,
    }

def blocker_eviction_crt(
    m: int,
    k: int,
    prime_positions: Iterable[tuple[int, int]],
) -> dict[str, int]:
    """Combine zero-carry small-prime control with eviction of old blockers.

    The chosen large-prime residue translates the old bad window by exactly
    its length.  Every old attached prime is therefore absent from the new
    window.  The returned class count includes every eviction residue while
    fixing the small-prime tier to the zero-carry class.
    """
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    controls = tuple(prime_positions)
    if not controls:
        raise ValueError("at least one blocker is required")

    first_position = m // 2 + 1
    bad_window_length = m - first_position + 1
    small_modulus = 1
    for small_p in primes_up_to(m):
        _, prime_power = least_power_above(m, small_p)
        small_modulus *= prime_power

    congruences: list[tuple[int, int]] = [(0, small_modulus)]
    blocker_product = 1
    eviction_class_count = 1
    seen: set[int] = set()
    for p, position in controls:
        _require_int("p", p, 2)
        _require_int("position", position, 1)
        if p <= m:
            raise ValueError("blocker primes must be larger than m")
        if not first_position <= position <= m:
            raise ValueError("each position must lie in the bad window")
        if p in seen or gcd(blocker_product, p) != 1:
            raise ValueError("blocker primes must be distinct and coprime")
        if (k + position) % p:
            raise ValueError("each blocker must divide its attached term")
        seen.add(p)
        blocker_product *= p
        eviction_class_count *= p - bad_window_length
        congruences.append(((k + bad_window_length) % p, p))

    canonical_residue, combined_modulus = crt_coprime(congruences)
    if combined_modulus != small_modulus * blocker_product:
        raise AssertionError("eviction CRT modulus product mismatch")
    candidate = canonical_residue
    if candidate <= k:
        candidate += ((k - candidate) // combined_modulus + 1) * combined_modulus

    for small_p in primes_up_to(m):
        if carry_count(candidate, m, small_p):
            raise AssertionError("zero-carry CRT class lost small-prime control")
    for p, _ in controls:
        if any(
            (candidate + position) % p == 0
            for position in range(first_position, m + 1)
        ):
            raise AssertionError("canonical CRT class retained an old blocker")

    return {
        "source_k": k,
        "bad_window_first_position": first_position,
        "bad_window_length": bad_window_length,
        "small_prime_zero_carry_modulus": small_modulus,
        "blocker_prime_count": len(controls),
        "blocker_prime_product": blocker_product,
        "combined_modulus": combined_modulus,
        "combined_eviction_class_count": eviction_class_count,
        "canonical_residue": canonical_residue,
        "least_forward_candidate": candidate,
        "least_forward_delta": candidate - k,
    }


def shift_first_bad_term_formula(m: int, k: int, p: int) -> dict[str, int]:
    """Identify an even-to-odd shift spike with its first bad-window term.

    Here ``m`` is the even source index.  For odd ``p > m+1`` dividing
    ``x=m+2k``, the first bad-window term after shifting to ``(m+1,k-1)`` is
    ``x/2`` and has exactly the shifted ``p``-slack.
    """
    _require_int("m", m, 2)
    _require_int("k", k, 2)
    _require_int("p", p, 3)
    if m % 2 or p % 2 == 0 or p <= m + 1:
        raise ValueError("require even m and odd p > m+1")
    x = m + 2 * k
    if x % p:
        raise ValueError("p must divide m+2k")
    first_bad_term = x // 2
    term_record = term_compensation_formula(first_bad_term, p)
    spike_record = natural_shift_spike_formula(m, k, p)
    shifted = spike_record["predicted_shifted_slack"]
    if term_record["predicted_slack"] != shifted:
        raise AssertionError("shift spike and first-term formulas disagree")
    return {
        "source_m": m,
        "source_k": k,
        "target_m": m + 1,
        "target_k": k - 1,
        "prime": p,
        "first_bad_term": first_bad_term,
        "exponent": term_record["exponent"],
        "cofactor": term_record["cofactor"],
        "positive_prefix_levels": term_record["positive_prefix_levels"],
        "first_bad_term_slack": term_record["predicted_slack"],
        "predicted_shifted_slack": shifted,
    }


def _verify_prefix_cone() -> dict[str, Any]:
    counts = Counter()
    rows: list[dict[str, int]] = []
    for p in (3, 5, 7, 11):
        for exponent in range(1, 5):
            modulus = p**exponent
            actual_count = 0
            for residue in range(modulus):
                safe = flexible_prefix_residue(p, exponent, residue)
                counts["enumerated_residues"] += 1
                if not safe:
                    continue
                actual_count += 1
                counts["safe_residue_formula_checks"] += 1
                term = p**exponent * residue
                record = term_compensation_formula(term, p)
                if record["exponent"] != exponent:
                    raise AssertionError("safe residue changed the requested valuation")
                if record["predicted_slack"] < 0:
                    raise AssertionError("safe prefix cone did not repair its prime")
                if slack(1, term - 1, p) != record["predicted_slack"]:
                    raise AssertionError("safe prefix cone disagrees with exact slack")
            expected_count = flexible_prefix_residue_count(p, exponent)
            if actual_count != expected_count:
                raise AssertionError("safe prefix residue count formula failed")
            rows.append(
                {
                    "prime": p,
                    "exponent": exponent,
                    "modulus": modulus,
                    "safe_residue_count": actual_count,
                }
            )
    return {"counts": dict(sorted(counts.items())), "rows": rows}


def _verify_local_separation(max_m: int, max_k: int) -> dict[str, int]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    counts = Counter()
    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            certificate = witness_certificate(m, k, factorizer)
            records = large_prime_window_records(m, k, factorizer)
            local_obstructions = {
                record["prime"]: record["predicted_slack"]
                for record in records
                if record["predicted_slack"] < 0
            }
            exact_obstructions = {
                p: value
                for p, value in certificate.obstructions.items()
                if p > m
            }
            if local_obstructions != exact_obstructions:
                raise AssertionError("local term criterion missed a large-prime slack")
            for record in records:
                if slack(m, k, record["prime"]) != record["predicted_slack"]:
                    raise AssertionError("local term formula disagrees with exact slack")
                if record["cofactor"] < record["prime"]:
                    short = short_cofactor_compensation_formula(
                        record["term"], record["prime"]
                    )
                    if short["predicted_slack"] != record["predicted_slack"]:
                        raise AssertionError("short-cofactor formula disagrees")
                    counts["short_cofactor_checks"] += 1
                    if record["predicted_slack"] < 0:
                        counts["short_cofactor_obstructions"] += 1
                    if record["exponent"] >= 2:
                        if record["predicted_slack"] >= 0:
                            raise AssertionError(
                                "short cofactor cannot repair a repeated prime"
                            )
                        counts["short_cofactor_prime_power_spikes"] += 1
            small_failure = any(p <= m for p in certificate.obstructions)
            separated_witness = not local_obstructions and not small_failure
            if certificate.is_witness != separated_witness:
                raise AssertionError("small/local separation changed witness truth")
            counts["pairs"] += 1
            counts["large_prime_factor_checks"] += len(records)
            counts["large_prime_obstructions"] += len(local_obstructions)
            if not local_obstructions:
                counts["large_prime_good_windows"] += 1
                if small_failure:
                    counts["large_good_but_small_prime_failed"] += 1
            if certificate.is_witness:
                counts["witnesses"] += 1
    return dict(sorted(counts.items()))


def _verify_shift_identity(max_m: int, max_k: int) -> dict[str, int]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    counts = Counter()
    for m in range(2, max_m + 1, 2):
        for k in range(2, max_k + 1):
            x = m + 2 * k
            for p in factorizer.factor(x):
                if p <= m + 1:
                    continue
                record = shift_first_bad_term_formula(m, k, p)
                exact = slack(m + 1, k - 1, p)
                if record["first_bad_term_slack"] != exact:
                    raise AssertionError("first bad term disagrees with shifted slack")
                if record["cofactor"] < p:
                    short = short_cofactor_compensation_formula(
                        record["first_bad_term"], p
                    )
                    if short["predicted_slack"] != exact:
                        raise AssertionError("short-cofactor shift formula disagrees")
                    counts["short_cofactor_checks"] += 1
                    if exact < 0:
                        counts["short_cofactor_spikes"] += 1
                    if record["exponent"] >= 2:
                        if exact >= 0:
                            raise AssertionError(
                                "short cofactor cannot repair a repeated shift prime"
                            )
                        counts["short_cofactor_prime_power_spikes"] += 1
                counts["prime_checks"] += 1
                if exact < 0:
                    counts["spikes"] += 1
    return dict(sorted(counts.items()))

def _verify_blocker_evictions(
    max_m: int,
    max_k: int,
    candidate_max_x: int,
) -> dict[str, Any]:
    if candidate_max_x > MAX_TRIAL_FACTORIZER_VALUE:
        raise ValueError("eviction candidate bound exceeds the exact factorizer limit")
    source_factorizer = SPFFactorizer(max_m + 2 * max_k)
    candidate_factorizer = TrialFactorizer(candidate_max_x)
    counts = Counter()
    representatives: list[dict[str, Any]] = []
    maximum_candidate_digits = 0

    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            old = large_prime_window_obstructions(m, k, source_factorizer)
            if not old:
                continue
            controls = tuple((record["prime"], record["position"]) for record in old)
            eviction = blocker_eviction_crt(m, k, controls)
            candidate = eviction["least_forward_candidate"]
            candidate_x = m + 2 * candidate
            maximum_candidate_digits = max(
                maximum_candidate_digits, len(str(candidate))
            )
            counts["source_systems"] += 1
            counts["old_blocker_prime_checks"] += len(old)
            if candidate_x > candidate_max_x:
                counts["canonical_candidates_beyond_exact_check_bound"] += 1
                continue

            certificate = witness_certificate(m, candidate, candidate_factorizer)
            small_obstructions = {
                p: value
                for p, value in certificate.obstructions.items()
                if p <= m
            }
            if small_obstructions:
                raise AssertionError("zero-carry eviction candidate failed small tier")
            new = large_prime_window_obstructions(m, candidate, candidate_factorizer)
            new_map = {
                record["prime"]: record["predicted_slack"] for record in new
            }
            exact_large = {
                p: value
                for p, value in certificate.obstructions.items()
                if p > m
            }
            if new_map != exact_large:
                raise AssertionError("eviction candidate large-prime mismatch")
            old_primes = {record["prime"] for record in old}
            if old_primes.intersection(new_map):
                raise AssertionError("an evicted blocker re-entered the new window")

            counts["canonical_candidates_exactly_checked"] += 1
            counts["new_large_obstruction_prime_checks"] += len(new)
            if certificate.is_witness:
                counts["canonical_candidate_witnesses"] += 1
            else:
                if not new:
                    raise AssertionError("nonwitness eviction candidate lacked blocker")
                counts["canonical_candidates_with_new_prime_interference"] += 1
            for record in new:
                if record["cofactor"] >= record["prime"]:
                    continue
                short = short_cofactor_compensation_formula(
                    record["term"], record["prime"]
                )
                if short["predicted_slack"] != record["predicted_slack"]:
                    raise AssertionError("new short-cofactor obstruction mismatch")
                counts["new_short_cofactor_obstructions"] += 1
                if record["exponent"] >= 2:
                    counts["new_short_cofactor_prime_power_spikes"] += 1

            if len(representatives) < 20:
                representatives.append(
                    {
                        "m": m,
                        "source_k": k,
                        "old_blockers": [
                            {
                                "prime": record["prime"],
                                "position": record["position"],
                                "slack": record["predicted_slack"],
                            }
                            for record in old
                        ],
                        "eviction": eviction,
                        "canonical_candidate_x": candidate_x,
                        "canonical_candidate_is_witness": certificate.is_witness,
                        "new_large_obstructions": list(new),
                    }
                )

    return {
        "parameters": {
            "min_m": 1,
            "max_m": max_m,
            "min_k": 1,
            "max_k": max_k,
            "candidate_max_x": candidate_max_x,
        },
        "counts": dict(sorted(counts.items())),
        "maximum_canonical_candidate_decimal_digits": maximum_candidate_digits,
        "reported_representatives": representatives,
    }


def _verify_moving_survivors(moving: dict[str, Any]) -> dict[str, Any]:
    target_m = moving["parameters"]["target_m"]
    moving_records = moving["sieve_survivors"]["reported_records"]
    moving_by_k = {record["k"]: record for record in moving_records}
    survivor_count = moving["counts"].get("single_level_sieve_survivors", 0)
    if survivor_count != len(moving_by_k):
        raise AssertionError("moving survivor count does not match its records")

    rows: list[dict[str, Any]] = []
    counts = Counter()
    factorizer: TrialFactorizer | None = None
    if moving_by_k:
        max_x = max(target_m + 2 * k for k in moving_by_k)
        factorizer = TrialFactorizer(max_x)
    for k in sorted(moving_by_k):
        assert factorizer is not None
        certificate = witness_certificate(target_m, k, factorizer)
        moving_obstructions = {
            record["prime"]: record["slack"]
            for record in moving_by_k[k]["obstructions"]
        }
        if moving_obstructions != certificate.obstructions:
            raise AssertionError("moving survivor certificate disagrees")
        local = large_prime_window_obstructions(target_m, k, factorizer)
        local_map = {
            record["prime"]: record["predicted_slack"] for record in local
        }
        exact_large = {
            p: value
            for p, value in certificate.obstructions.items()
            if p > target_m
        }
        if local_map != exact_large:
            raise AssertionError("survivor local obstructions are incomplete")
        small = [
            {"prime": p, "slack": value}
            for p, value in certificate.obstructions.items()
            if p <= target_m
        ]

        row: dict[str, Any] = {
            "k": k,
            "offset": moving_by_k[k]["offset"],
            "is_witness": certificate.is_witness,
            "exact_obstruction_count": len(certificate.obstructions),
            "large_prime_obstruction_count": len(local),
            "small_prime_obstructions": small,
            "obstructions": list(local),
        }
        if certificate.is_witness:
            counts["witnesses"] += 1
        if small:
            counts["survivors_with_small_prime_obstructions"] += 1
        if local:
            controls = tuple(
                (record["prime"], record["position"]) for record in local
            )
            displacement = retained_prime_displacement(k, controls)
            eviction = blocker_eviction_crt(target_m, k, controls)
            counts["zero_carry_eviction_crt_systems"] += 1
            prime_product = displacement["prime_product"]
            flexible_modulus = 1
            flexible_choices = 1
            for record in local:
                p = record["prime"]
                exponent = record["exponent"]
                flexible_modulus *= p ** (2 * exponent)
                flexible_choices *= flexible_prefix_residue_count(p, exponent)
                counts["obstruction_prime_checks"] += 1
                if record["cofactor"] < p:
                    short = short_cofactor_compensation_formula(
                        record["term"], p
                    )
                    if short["predicted_slack"] != record["predicted_slack"]:
                        raise AssertionError(
                            "survivor short-cofactor formula disagrees"
                        )
                    counts["short_cofactor_obstruction_checks"] += 1
                    if record["exponent"] >= 2:
                        counts["short_cofactor_prime_power_spikes"] += 1
            next_k = displacement["least_forward_candidate"]
            for p, position in controls:
                if (next_k + position) % p:
                    raise AssertionError(
                        "product displacement lost a retained prime"
                    )
            row.update(
                {
                    "retained_prime_product": prime_product,
                    "retained_prime_product_decimal_digits": len(
                        str(prime_product)
                    ),
                    "least_forward_retaining_candidate": next_k,
                    "least_forward_retaining_candidate_decimal_digits": len(
                        str(next_k)
                    ),
                    "exceeds_factorizer_bound": (
                        target_m + 2 * next_k > MAX_TRIAL_FACTORIZER_VALUE
                    ),
                    "flexible_prefix_modulus": flexible_modulus,
                    "flexible_prefix_choice_count": flexible_choices,
                    "flexible_prefix_density": {
                        "numerator": flexible_choices,
                        "denominator": flexible_modulus,
                    },
                    "zero_carry_blocker_eviction": eviction,
                    "zero_carry_eviction_candidate_decimal_digits": len(
                        str(eviction["least_forward_candidate"])
                    ),
                    "zero_carry_eviction_candidate_exceeds_factorizer_bound": (
                        target_m + 2 * eviction["least_forward_candidate"]
                        > MAX_TRIAL_FACTORIZER_VALUE
                    ),
                }
            )
            counts["survivors_with_large_prime_obstructions"] += 1
        rows.append(row)
        counts["survivors"] += 1

    rejected = moving["counts"]["single_level_rejected_offsets"]
    if rejected + counts["survivors"] != moving["parameters"]["offset_count"]:
        raise AssertionError("moving artifact does not partition its interval")
    counts["single_level_rejected_offsets"] = rejected
    large_obstructed = (
        rejected + counts["survivors_with_large_prime_obstructions"]
    )
    counts["offsets_with_certified_large_prime_obstruction"] = large_obstructed
    if large_obstructed == moving["parameters"]["offset_count"]:
        counts["all_offsets_with_certified_large_prime_obstruction"] = (
            large_obstructed
        )
    return {"counts": dict(sorted(counts.items())), "rows": rows}


def _verify_survivor_bounds(
    moving: dict[str, Any],
    near_miss: dict[str, Any],
) -> dict[str, Any]:
    bounds = _verify_moving_survivors(moving)
    near_by_k = {
        record["k"]: record
        for record in near_miss["minimum_window_representatives"]
    }
    rows_by_k = {record["k"]: record for record in bounds["rows"]}
    if set(rows_by_k) != set(near_by_k):
        raise AssertionError("moving and near-miss survivor sets disagree")
    combined_by_k = {
        system["source_k"]: system
        for system in near_miss["combined_crt_systems"]
    }
    for k, row in rows_by_k.items():
        near_obstructions = {
            record["prime"]: record["slack"]
            for record in near_by_k[k]["exact_obstructions"]
        }
        exact_obstructions = {
            record["prime"]: record["predicted_slack"]
            for record in row["obstructions"]
        }
        exact_obstructions.update(
            {
                record["prime"]: record["slack"]
                for record in row["small_prime_obstructions"]
            }
        )
        if near_obstructions != exact_obstructions:
            raise AssertionError("near-miss survivor certificate disagrees")
        system = combined_by_k[k]
        if system["combined_solution"] == k:
            raise AssertionError("the obstructed source cannot be a repair")
        row["canonical_mirror_crt_solution"] = system["combined_solution"]
    return bounds


def analyze(
    *,
    max_m: int,
    max_k: int,
    shift_max_m: int,
    shift_max_k: int,
    eviction_max_m: int,
    eviction_max_k: int,
    eviction_candidate_max_x: int,
    moving_path: Path,
    near_miss_path: Path,
    extended_moving_path: Path | None = None,
) -> dict[str, Any]:
    if min(
        max_m,
        max_k,
        shift_max_m,
        shift_max_k,
        eviction_max_m,
        eviction_max_k,
        eviction_candidate_max_x,
    ) < 1:
        raise ValueError("all finite bounds must be positive")
    moving_path = moving_path.resolve()
    near_miss_path = near_miss_path.resolve()
    moving_bytes = moving_path.read_bytes()
    near_miss_bytes = near_miss_path.read_bytes()
    moving = json.loads(moving_bytes)
    near_miss = json.loads(near_miss_bytes)
    if moving["status"] != "EXACT_BOUNDED_MOVING_BAD_WINDOW_SEARCH":
        raise ValueError("unexpected moving-window artifact status")
    if near_miss["status"] != "EXACT_BOUNDED_BAD_WINDOW_NEAR_MISS_ANALYSIS":
        raise ValueError("unexpected near-miss artifact status")

    result = {
        "schema_version": 3,
        "status": "PASS_EXACT_COMPENSATION_STRUCTURE_ANALYSIS",
        "claim_boundary": (
            "The local decomposition, short-cofactor formula, prefix-cone count, "
            "retained-prime product bound, blocker-eviction CRT, and shift identity "
            "are proved algebraically in THEOREMS.md. The reported counts are exact "
            "only over their declared finite domains. No witness or existence claim "
            "is made outside those domains."
        ),
        "parameters": {
            "local_rectangle": {
                "min_m": 1,
                "max_m": max_m,
                "min_k": 1,
                "max_k": max_k,
            },
            "shift_rectangle": {
                "even_m_min": 2,
                "even_m_max": shift_max_m,
                "min_k": 2,
                "max_k": shift_max_k,
            },
            "eviction_rectangle": {
                "min_m": 1,
                "max_m": eviction_max_m,
                "min_k": 1,
                "max_k": eviction_max_k,
                "candidate_max_x": eviction_candidate_max_x,
            },
        },
        "flexible_prefix_cone": _verify_prefix_cone(),
        "small_and_large_prime_separation": {
            "counts": _verify_local_separation(max_m, max_k)
        },
        "even_shift_first_bad_term_identity": {
            "counts": _verify_shift_identity(shift_max_m, shift_max_k)
        },
        "zero_carry_blocker_eviction": _verify_blocker_evictions(
            eviction_max_m,
            eviction_max_k,
            eviction_candidate_max_x,
        ),
        "m27_survivor_retained_prime_bounds": _verify_survivor_bounds(
            moving, near_miss
        ),
        "input": {
            "moving_window": {
                "logical_name": moving_path.name,
                "sha256": hashlib.sha256(moving_bytes).hexdigest(),
            },
            "near_miss": {
                "logical_name": near_miss_path.name,
                "sha256": hashlib.sha256(near_miss_bytes).hexdigest(),
            },
        },
        "implementation_sha256": {
            "compensation_structure_analysis.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
        },
    }
    if extended_moving_path is not None:
        extended_moving_path = extended_moving_path.resolve()
        extended_bytes = extended_moving_path.read_bytes()
        extended = json.loads(extended_bytes)
        if extended["status"] != "EXACT_BOUNDED_MOVING_BAD_WINDOW_SEARCH":
            raise ValueError("unexpected extended moving-window artifact status")
        if extended["parameters"]["target_m"] != moving["parameters"]["target_m"]:
            raise ValueError("moving-window target indices disagree")
        expected_start = (
            moving["parameters"]["start_k"]
            + moving["parameters"]["offset_count"]
        )
        if extended["parameters"]["start_k"] != expected_start:
            raise ValueError("moving-window artifacts are not contiguous")
        result["parameters"]["contiguous_moving_offset_count"] = (
            moving["parameters"]["offset_count"]
            + extended["parameters"]["offset_count"]
        )
        result["m27_extended_survivor_retained_prime_bounds"] = (
            _verify_moving_survivors(extended)
        )
        result["input"]["extended_moving_window"] = {
            "logical_name": extended_moving_path.name,
            "sha256": hashlib.sha256(extended_bytes).hexdigest(),
        }
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=20)
    parser.add_argument("--max-k", type=int, default=5_000)
    parser.add_argument("--shift-max-m", type=int, default=20)
    parser.add_argument("--shift-max-k", type=int, default=2_000)
    parser.add_argument("--eviction-max-m", type=int, default=8)
    parser.add_argument("--eviction-max-k", type=int, default=500)
    parser.add_argument(
        "--eviction-candidate-max-x",
        type=int,
        default=100_000_000,
    )
    parser.add_argument("--moving", type=Path, default=DEFAULT_MOVING)
    parser.add_argument("--near-miss", type=Path, default=DEFAULT_NEAR_MISS)
    parser.add_argument(
        "--extended-moving",
        type=Path,
        help=f"Optional contiguous artifact; canonical: {DEFAULT_EXTENDED_MOVING.name}",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        max_m=args.max_m,
        max_k=args.max_k,
        shift_max_m=args.shift_max_m,
        shift_max_k=args.shift_max_k,
        eviction_max_m=args.eviction_max_m,
        eviction_max_k=args.eviction_max_k,
        eviction_candidate_max_x=args.eviction_candidate_max_x,
        moving_path=args.moving,
        near_miss_path=args.near_miss,
        extended_moving_path=args.extended_moving,
    )
    result["resources"] = {
        "wall_seconds": time.perf_counter() - started_wall,
        "process_cpu_seconds": time.process_time() - started_cpu,
        "processes": 1,
    }
    atomic_write_json(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "local_counts": result["small_and_large_prime_separation"][
                    "counts"
                ],
                "survivor_counts": result[
                    "m27_survivor_retained_prime_bounds"
                ]["counts"],
                "extended_survivor_counts": result.get(
                    "m27_extended_survivor_retained_prime_bounds", {}
                ).get("counts"),
                "wall_seconds": result["resources"]["wall_seconds"],
                "process_cpu_seconds": result["resources"][
                    "process_cpu_seconds"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
