#!/usr/bin/env python3
"""Exact finite analysis of local CRT controls and prime-power compensation.

The proved statements emitted by this script are elementary consequences of
Kummer's and Legendre's formulas.  Every aggregate count is bounded by the
parameters recorded in the output artifact.
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
    SPFFactorizer,
    carry_count,
    integer_valuation,
    primes_up_to,
    slack,
    zone_contribution,
)
from shift_repair_analysis import scan_witness_bits

ROOT = Path(__file__).resolve().parent
DEFAULT_ATLAS = ROOT / "data" / "atlas_m1_50_k50000.json"
DEFAULT_OUTPUT = ROOT / "data" / "residue_control_analysis_m1_50_k50000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def least_power_above(m: int, p: int, extra_digits: int = 0) -> tuple[int, int]:
    """Return ``(a, p**a)`` with minimal ``p**a > m``, plus extra digits."""
    if m < 0 or p < 2 or extra_digits < 0:
        raise ValueError("require m>=0, p>=2, and extra_digits>=0")
    exponent = 1
    power = p
    while power <= m:
        exponent += 1
        power *= p
    return exponent + extra_digits, power * p**extra_digits


def local_safe_residues(
    m: int, p: int, extra_digits: int = 0
) -> tuple[int, tuple[int, ...]]:
    """Residues that certify ``s_p(m, k) >= 0`` for every congruent ``k``.

    The modulus is a power ``q=p**a>m``.  Requiring ``r+m<q`` prevents a carry
    from the low ``a`` digits into the unrestricted high digits when adding
    ``m``.  The low-digit right carry count can then be compared with the
    complete left carry count.
    """
    _, modulus = least_power_above(m, p, extra_digits)
    residues = tuple(
        r
        for r in range(modulus - m)
        if carry_count(r, r + m, p) >= carry_count(r, m, p)
    )
    return modulus, residues


def carry_free_residue_count(m: int, p: int) -> int:
    """Count residues modulo the least power above ``m`` with no left carry."""
    _, modulus = least_power_above(m, p)
    count = 1
    value = m
    while modulus > 1:
        count *= p - value % p
        value //= p
        modulus //= p
    return count


def local_safe_profile(m: int, extra_digits: int = 0) -> dict[str, Any]:
    modulus_product = 1
    class_product = 1
    rows: list[dict[str, Any]] = []
    for p in primes_up_to(m):
        exponent, modulus = least_power_above(m, p, extra_digits)
        actual_modulus, residues = local_safe_residues(m, p, extra_digits)
        assert actual_modulus == modulus and residues
        carry_free_count = carry_free_residue_count(m, p) if extra_digits == 0 else None
        rows.append(
            {
                "prime": p,
                "exponent": exponent,
                "modulus": modulus,
                "safe_residue_count": len(residues),
                "carry_free_residue_count": carry_free_count,
                "safe_residue_representatives": list(residues[:12]),
            }
        )
        modulus_product *= modulus
        class_product *= len(residues)
    divisor = gcd(class_product, modulus_product)
    return {
        "m": m,
        "extra_digits": extra_digits,
        "crt_modulus": modulus_product,
        "crt_safe_class_count": class_product,
        "density": {
            "numerator": class_product // divisor,
            "denominator": modulus_product // divisor,
        },
        "prime_rows": rows,
    }

def prime_power_band(m: int, p: int, exponent: int) -> tuple[int, int]:
    """Inclusive k-band where ``p**exponent <= m+2k < p**(exponent+1)``."""
    if m < 0 or p < 2 or exponent < 1:
        raise ValueError("require m>=0, p>=2, and exponent>=1")
    power = p**exponent
    lower = max(1, (power - m + 1) // 2)
    upper = (power * p - 1 - m) // 2
    return lower, upper


def band_residue_slack(m: int, p: int, exponent: int, residue: int) -> int:
    """Slack on a prime-power band as a function modulo ``p**exponent``."""
    modulus = p**exponent
    residue %= modulus
    representative = residue if residue else modulus
    power = p
    total = 0
    for _ in range(exponent):
        total += zone_contribution(m, representative, power)
        power *= p
    return total


def band_good_residues(m: int, p: int, exponent: int) -> tuple[int, ...]:
    """Exact residues with nonnegative p-slack on one prime-power band."""
    modulus = p**exponent
    return tuple(
        residue
        for residue in range(modulus)
        if band_residue_slack(m, p, exponent, residue) >= 0
    )


def upper_half_prefix_count(cofactor: int, p: int) -> int:
    """Count base-``p`` prefixes ``u_r`` satisfying ``2*u_r > p**r``."""
    if cofactor < 1 or p < 2:
        raise ValueError("require cofactor>=1 and p>=2")
    count = 0
    power = p
    while power <= 2 * cofactor:
        if 2 * (cofactor % power) > power:
            count += 1
        power *= p
    return count


def large_prime_bad_window_formula(m: int, k: int, p: int) -> dict[str, int]:
    """Return the exact slack formula when ``p>m`` divides the bad window."""
    if m < 1 or k < 1 or p <= m:
        raise ValueError("require m>=1, k>=1, and p>m")
    positions = [
        i
        for i in range(m // 2 + 1, m + 1)
        if (k + i) % p == 0
    ]
    if len(positions) != 1:
        raise ValueError("p must divide exactly one bad-window term")
    position = positions[0]
    term = k + position
    exponent = integer_valuation(term, p)
    cofactor = term // p**exponent
    positive_levels = upper_half_prefix_count(cofactor, p)
    return {
        "position": position,
        "term": term,
        "exponent": exponent,
        "cofactor": cofactor,
        "positive_prefix_levels": positive_levels,
        "predicted_slack": positive_levels - exponent,
    }


def natural_shift_spike_formula(m: int, k: int, p: int) -> dict[str, int]:
    """Exact source/shifted slack at an odd ``p>m`` dividing ``x=m+2k``."""
    if m < 1 or k < 2 or p <= m or p % 2 == 0:
        raise ValueError("require m>=1, k>=2, and odd p>m")
    x = m + 2 * k
    if x % p:
        raise ValueError("p must divide m+2k")
    exponent = integer_valuation(x, p)
    cofactor = x // p**exponent
    positive_levels = sum(
        (cofactor // p**r) % 2
        for r in range(1, cofactor.bit_length() + 1)
        if p**r <= cofactor
    )
    source_slack = positive_levels + (exponent if m % 2 else 0)
    added = integer_valuation(m + 1, p)
    return {
        "exponent": exponent,
        "cofactor": cofactor,
        "positive_tail_levels": positive_levels,
        "m_plus_one_valuation": added,
        "predicted_source_slack": source_slack,
        "predicted_shifted_slack": source_slack + added - exponent,
    }


def crt_coprime(congruences: Iterable[tuple[int, int]]) -> tuple[int, int]:
    """Combine pairwise-coprime ``(residue, modulus)`` congruences."""
    value = 0
    modulus = 1
    for residue, next_modulus in congruences:
        if next_modulus < 1 or gcd(modulus, next_modulus) != 1:
            raise ValueError("CRT moduli must be positive and pairwise coprime")
        step = ((residue - value) * pow(modulus, -1, next_modulus)) % next_modulus
        value += modulus * step
        modulus *= next_modulus
    return value % modulus, modulus

def fixed_class_single_level_obstruction(
    m: int,
    modulus: int,
    residue: int,
    position: int,
    prime: int,
) -> dict[str, int]:
    """Construct one Dirichlet-progression obstruction in a fixed CRT class.

    ``prime`` is assumed prime and must already lie in the required reduced
    progression.  Dirichlet's theorem supplies infinitely many such primes.
    """
    if m < 1 or modulus < 1 or prime < 2:
        raise ValueError("require m>=1, modulus>=1, and prime>=2")
    if not m // 2 < position <= m:
        raise ValueError("position must lie in the bad window")
    residue %= modulus
    divisor = gcd(residue + position, modulus)
    reduced_modulus = modulus // divisor
    reduced_residue = (residue + position) // divisor
    if prime % reduced_modulus != reduced_residue % reduced_modulus:
        raise ValueError("prime is not in the required reduced progression")
    if prime <= max(m, 2 * divisor):
        raise ValueError("prime is too small for the single-level certificate")
    k = divisor * prime - position
    x = m + 2 * k
    if k < 1 or k % modulus != residue or prime * prime <= x:
        raise AssertionError("fixed-class construction invariant failed")
    return {
        "m": m,
        "modulus": modulus,
        "residue": residue,
        "position": position,
        "divisor": divisor,
        "prime": prime,
        "k": k,
        "x": x,
        "predicted_slack": -1,
    }


def _check_fixed_class_obstructions() -> list[dict[str, int]]:
    examples = (
        (2, 6, 1, 2, 101),
        (4, 30, 7, 3, 103),
        (9, 210, 13, 7, 107),
        (5, 12, 5, 3, 101),
    )
    records: list[dict[str, int]] = []
    for arguments in examples:
        record = fixed_class_single_level_obstruction(*arguments)
        if slack(record["m"], record["k"], record["prime"]) != -1:
            raise AssertionError("fixed-class single-level obstruction failed")
        records.append(record)
    return records

def _check_band_periodicity() -> tuple[int, list[dict[str, int]]]:
    cases = (
        (3, 5, 3),
        (2, 7, 2),
        (4, 5, 2),
        (5, 7, 2),
        (3, 3, 3),
        (3, 2, 4),
        (9, 5, 2),
        (9, 2, 6),
    )
    checks = 0
    records: list[dict[str, int]] = []
    for m, p, exponent in cases:
        lower, upper = prime_power_band(m, p, exponent)
        modulus = p**exponent
        good = band_good_residues(m, p, exponent)
        for k in range(lower, upper + 1):
            expected = band_residue_slack(m, p, exponent, k)
            if slack(m, k, p) != expected:
                raise AssertionError("prime-power band periodicity failed")
            checks += 1
        records.append(
            {
                "m": m,
                "prime": p,
                "exponent": exponent,
                "band_lower_k": lower,
                "band_upper_k": upper,
                "modulus": modulus,
                "good_residue_count": len(good),
            }
        )
    return checks, records


def _check_local_safe_profiles(max_m: int) -> tuple[list[dict[str, Any]], int]:
    profiles: list[dict[str, Any]] = []
    extension_checks = 0
    for m in range(1, max_m + 1):
        profile = local_safe_profile(m)
        profiles.append(profile)
        for row in profile["prime_rows"]:
            p = row["prime"]
            modulus, residues = local_safe_residues(m, p)
            carry_free = sum(carry_count(r, m, p) == 0 for r in range(modulus))
            if carry_free != row["carry_free_residue_count"]:
                raise AssertionError("carry-free digit count mismatch")
            for r in residues:
                left_low = carry_count(r, m, p)
                right_low = carry_count(r, r + m, p)
                if r + m >= modulus or right_low < left_low:
                    raise AssertionError("invalid local-safe residue")
                for high in range(4):
                    k = r + high * modulus
                    if carry_count(k, k + m, p) < carry_count(k, m, p):
                        raise AssertionError("local-safe residue did not extend")
                    extension_checks += 1
    return profiles, extension_checks


def _check_general_large_prime_formula(
    max_m: int, max_k: int
) -> tuple[dict[str, int], list[dict[str, int]]]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    counts = Counter()
    compensated_examples: list[dict[str, int]] = []
    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            seen: set[int] = set()
            for position in range(m // 2 + 1, m + 1):
                for p in factorizer.factor(k + position):
                    if p <= m or p in seen:
                        continue
                    seen.add(p)
                    record = large_prime_bad_window_formula(m, k, p)
                    actual = slack(m, k, p)
                    if record["predicted_slack"] != actual:
                        raise AssertionError(
                            f"large-prime formula mismatch at {(m, k, p)}"
                        )
                    counts["formula_checks"] += 1
                    category = (
                        "negative"
                        if actual < 0
                        else "zero" if actual == 0 else "positive"
                    )
                    counts[f"slack_{category}"] += 1
                    if actual >= 0 and len(compensated_examples) < 40:
                        compensated_examples.append(
                            {"m": m, "k": k, "prime": p, **record}
                        )
    return dict(sorted(counts.items())), compensated_examples


def _check_mirror_lifts(max_m: int) -> tuple[int, list[dict[str, int]]]:
    primes = primes_up_to(2 * max_m + 100)
    checks = 0
    examples: list[dict[str, int]] = []
    for m in range(1, max_m + 1):
        large_primes = [p for p in primes if p > m][:2]
        for position in range(m // 2 + 1, m + 1):
            for p in large_primes:
                for exponent in (1, 2):
                    power = p**exponent
                    lift_digits = exponent if p % 2 else exponent + 1
                    k = power * (p**lift_digits - 1) - position
                    record = large_prime_bad_window_formula(m, k, p)
                    if record["exponent"] != exponent:
                        raise AssertionError("mirror lift changed the prescribed valuation")
                    if record["predicted_slack"] != 0 or slack(m, k, p) != 0:
                        raise AssertionError("mirror lift did not exactly cancel")
                    checks += 1
                    if len(examples) < 24:
                        examples.append({"m": m, "k": k, "prime": p, **record})
    return checks, examples


def _check_combined_crt(max_m: int) -> tuple[int, list[dict[str, Any]]]:
    prime_pool = primes_up_to(2 * max_m + 100)
    checks = 0
    examples: list[dict[str, Any]] = []
    for m in range(2, max_m + 1):
        congruences: list[tuple[int, int]] = []
        controlled_small: list[int] = []
        for p in primes_up_to(m):
            modulus, residues = local_safe_residues(m, p)
            residue = residues[-1]
            congruences.append((residue, modulus))
            controlled_small.append(p)
        positions = list(range(m // 2 + 1, m + 1))
        large = [p for p in prime_pool if p > m][: len(positions)]
        controlled_large: list[dict[str, int]] = []
        for position, p in zip(positions, large, strict=True):
            exponent = 1
            power = p**exponent
            modulus = p ** (2 * exponent)
            residue = (-position - power) % modulus
            congruences.append((residue, modulus))
            controlled_large.append(
                {"prime": p, "position": position, "exponent": exponent}
            )
        k, modulus = crt_coprime(congruences)
        if k == 0:
            k = modulus
        for p in controlled_small:
            if slack(m, k, p) < 0:
                raise AssertionError("combined CRT lost a small-prime control")
            checks += 1
        for item in controlled_large:
            record = large_prime_bad_window_formula(m, k, item["prime"])
            if record["predicted_slack"] < 0 or slack(m, k, item["prime"]) < 0:
                raise AssertionError("combined CRT lost a mirror-lift control")
            checks += 1
        if len(examples) < 12:
            examples.append(
                {
                    "m": m,
                    "least_positive_solution": k,
                    "modulus": modulus,
                    "controlled_small_primes": controlled_small,
                    "controlled_large_prime_lifts": controlled_large,
                }
            )
    return checks, examples


def _analyze_natural_shift_spikes(
    min_m: int,
    max_m: int,
    max_k: int,
    witness_bits: list[bytearray | None],
) -> dict[str, Any]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    rows: list[dict[str, int]] = []
    totals = Counter()
    examples: list[dict[str, Any]] = []
    for m in range(min_m, max_m):
        source = witness_bits[m]
        target = witness_bits[m + 1]
        assert source is not None and target is not None
        counts = Counter()
        for k in range(2, max_k + 1):
            if not source[k]:
                continue
            success = bool(target[k - 1])
            counts["eligible"] += 1
            counts["successful" if success else "failed"] += 1
            x = m + 2 * k
            exact_large: list[dict[str, int | bool]] = []
            old_fatal: list[int] = []
            for p, exponent in sorted(factorizer.factor(x).items()):
                if p <= m or p % 2 == 0:
                    continue
                record = natural_shift_spike_formula(m, k, p)
                actual_source = slack(m, k, p)
                actual_shifted = slack(m + 1, k - 1, p)
                if record["exponent"] != exponent:
                    raise AssertionError("factorization valuation mismatch")
                if record["predicted_source_slack"] != actual_source:
                    raise AssertionError("source spike formula mismatch")
                if record["predicted_shifted_slack"] != actual_shifted:
                    raise AssertionError("shifted spike formula mismatch")
                counts["formula_prime_checks"] += 1
                is_old_fatal = (
                    m % 2 == 0
                    and p > m + 1
                    and p ** (exponent + 1) > x + m
                )
                if is_old_fatal:
                    old_fatal.append(p)
                if actual_shifted < 0:
                    exact_large.append(
                        {"prime": p, "old_fatal_criterion": is_old_fatal, **record}
                    )
            if old_fatal:
                counts["sources_with_old_fatal_spike"] += 1
            if exact_large:
                if success:
                    raise AssertionError("exact large spike survived the shift")
                counts["failed_with_exact_large_spike"] += 1
                if not old_fatal:
                    counts["additional_failures_from_multilevel_formula"] += 1
                counts["exact_large_obstruction_primes"] += len(exact_large)
                if len(examples) < 50:
                    examples.append({"m": m, "k": k, "x": x, "spikes": exact_large})
            elif not success:
                counts["failed_without_exact_large_spike"] += 1
        row = {"m": m, **dict(sorted(counts.items()))}
        rows.append(row)
        totals.update(counts)
    return {
        "totals": dict(sorted(totals.items())),
        "rows": rows,
        "examples": examples,
    }


def analyze(
    *,
    min_m: int,
    max_m: int,
    max_k: int,
    atlas_path: Path,
    formula_max_m: int = 8,
    formula_max_k: int = 2_000,
) -> dict[str, Any]:
    if not 1 <= min_m < max_m or max_k < 2:
        raise ValueError("require 1<=min_m<max_m and max_k>=2")
    if formula_max_m < 1 or formula_max_k < 1:
        raise ValueError("formula bounds must be positive")
    atlas_path = atlas_path.resolve()
    atlas_bytes = atlas_path.read_bytes()
    atlas = json.loads(atlas_bytes)
    for name, expected in (("min_m", min_m), ("max_m", max_m), ("max_k", max_k)):
        if atlas["parameters"][name] != expected:
            raise AssertionError(f"atlas parameter mismatch: {name}")
    if atlas["schema_version"] != 2:
        raise AssertionError("control analysis requires repaired atlas schema 2")

    witness_bits = scan_witness_bits(min_m, max_m, max_k)
    profiles, extension_checks = _check_local_safe_profiles(max_m)
    deep_profiles = [local_safe_profile(27, extra) for extra in range(3)] if max_m >= 27 else []
    formula_counts, compensated_examples = _check_general_large_prime_formula(
        formula_max_m, formula_max_k
    )
    mirror_checks, mirror_examples = _check_mirror_lifts(min(max_m, 20))
    crt_checks, crt_examples = _check_combined_crt(min(max_m, 16))
    spike_analysis = _analyze_natural_shift_spikes(
        min_m, max_m, max_k, witness_bits
    )
    fixed_class_examples = _check_fixed_class_obstructions()
    band_checks, band_records = _check_band_periodicity()

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_BOUNDED_CONTROL_ANALYSIS",
        "claim_boundary": (
            "The band, local-safe, bad-window compensation, mirror-lift, and "
            "shift-spike formulas are proved elementary implications. Fixed-class "
            "insufficiency additionally invokes Dirichlet's theorem. Aggregate counts "
            "and searches are exact only in their declared finite domains; no global "
            "existence claim is made."
        ),
        "parameters": {
            "min_m": min_m,
            "max_m": max_m,
            "max_k": max_k,
            "general_formula_max_m": formula_max_m,
            "general_formula_max_k": formula_max_k,
            "local_safe_extension_high_digits": [0, 1, 2, 3],
        },
        "proved_lemmas": [
            {
                "name": "prime_power_band_periodicity",
                "statement": (
                    "On p^L<=m+2k<p^(L+1), s_p(m,k) depends only on "
                    "k modulo p^L, so its nonnegative set is an exact union "
                    "of prime-power residue classes."
                ),
            },
            {
                "name": "local_safe_residue",
                "statement": (
                    "For q=p^a>m, if r+m<q and carries_p(r,r+m)>=carries_p(r,m), "
                    "then every k congruent to r modulo q has s_p(m,k)>=0."
                ),
            },
            {
                "name": "fixed_crt_class_insufficiency",
                "statement": (
                    "Dirichlet's theorem implies that every residue class modulo "
                    "every fixed modulus contains arbitrarily large k with a "
                    "single-level bad-window obstruction."
                ),
            },
            {
                "name": "large_prime_bad_window_compensation",
                "statement": (
                    "If p>m and p^e exactly divides k+i for a unique bad-window "
                    "position i, then s_p=-e plus the number of r>=1 for which "
                    "2*((k+i)/p^e mod p^r)>p^r."
                ),
            },
            {
                "name": "mirror_lift",
                "statement": (
                    "For odd p>m, the lift k+i congruent to -p^e modulo "
                    "p^(2e) contributes e strict-good levels after e bad levels. "
                    "For p=2 (only possible here when m=1), modulus 2^(2e+1) "
                    "is needed instead."
                ),
            },
            {
                "name": "natural_shift_spike_formula",
                "statement": (
                    "For odd p>m, p^a exactly dividing x=m+2k, and y=x/p^a, "
                    "s_p(m,k)=a*(m mod 2)+sum_{r>=1}(floor(y/p^r) mod 2); "
                    "the shifted slack adds v_p(m+1)-a."
                ),
            },
        ],
        "prime_power_band_crt": {
            "periodicity_checks": band_checks,
            "profiles": band_records,
        },
        "small_prime_local_safe_crt": {
            "profiles": profiles,
            "m27_deepening_profiles": deep_profiles,
            "extension_checks": extension_checks,
        },
        "fixed_class_obstruction_checks": {
            "exact_checks": len(fixed_class_examples),
            "examples": fixed_class_examples,
        },
        "general_large_prime_compensation": {
            "finite_domain": {"max_m": formula_max_m, "max_k": formula_max_k},
            "counts": formula_counts,
            "compensated_examples": compensated_examples,
        },
        "mirror_lift_checks": {
            "exact_checks": mirror_checks,
            "examples": mirror_examples,
        },
        "combined_crt_checks": {
            "controlled_prime_checks": crt_checks,
            "examples": crt_examples,
        },
        "natural_shift_spike_analysis": spike_analysis,
        "input": {
            "atlas_logical_name": atlas_path.name,
            "atlas_sha256": hashlib.sha256(atlas_bytes).hexdigest(),
            "atlas_schema_version": atlas["schema_version"],
        },
        "implementation_sha256": {
            "residue_control_analysis.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "shift_repair_analysis.py": sha256_file(ROOT / "shift_repair_analysis.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-m", type=int, default=1)
    parser.add_argument("--max-m", type=int, default=50)
    parser.add_argument("--max-k", type=int, default=50_000)
    parser.add_argument("--formula-max-m", type=int, default=8)
    parser.add_argument("--formula-max-k", type=int, default=2_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        min_m=args.min_m,
        max_m=args.max_m,
        max_k=args.max_k,
        atlas_path=args.atlas,
        formula_max_m=args.formula_max_m,
        formula_max_k=args.formula_max_k,
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
                "finite_checks": {
                    "local_extensions": result["small_prime_local_safe_crt"][
                        "extension_checks"
                    ],
                    "large_prime_formula": result[
                        "general_large_prime_compensation"
                    ]["counts"]["formula_checks"],
                    "mirror_lifts": result["mirror_lift_checks"]["exact_checks"],
                },
                "spike_totals": result["natural_shift_spike_analysis"]["totals"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
