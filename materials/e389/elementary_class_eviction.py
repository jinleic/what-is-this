#!/usr/bin/env python3
"""Evict every residue class from Erdos #389 by an elementary large-prime construction.

Necessity (25) makes a bad-window term fatal as soon as one prime ``p > m`` has
``p**v_p(w) > sqrt(2w)``.  Section 15 turns that necessity into a construction.
A class ``c mod M`` pins

    d_i = prod_{p | M} p**min(v_p(c+i), a_p),   q_i = M / d_i,
    a_i = ((c+i)/d_i) mod q_i     with gcd(a_i, q_i) = 1,

on the window term ``k+i`` for every ``k = c (mod M)``; the cofactor
``u_i = (k+i)/d_i`` is free inside the reduced class ``a_i mod q_i``.  Any
factorisation ``u_i = p*t`` with ``p`` prime and ``p > 2 d_i t`` is fatal, because
then ``p**2 > 2 d_i p t = 2(k+i)``.

Choosing ``p`` first and solving ``t = a_i p^{-1} (mod q_i)`` needs no prime
localisation: one prime in an interval of ratio two, and one residue in an
interval of length at least ``q_i``.  That is why the construction beats both
Linnik (which places its prime below the target block) and GRH (which reaches
only ``q <= x^{1/2-eps}``).

Checked here exactly, with no sampling inside a declared rectangle:

  (a) the fatality criterion against the repository's Legendre/Kummer verifier;
  (b) Theorem 15.3 - every arithmetic progression of length at least the least
      prime above sqrt(2*top) contains a non-witness, whatever its difference;
  (c) its corollaries: a non-witness below 2*M*max(2M, m) in every class, and
      one inside [N, 2N) whenever M*P0(N,m) <= N;
  (d) the measured barrier past the proved threshold, by scanning primes;
  (e) the arithmetic of the modulus family {q_i} that any analytic route faces.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import sys
import time
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
# ``_is_prime_64`` is the repository's deterministic 64-bit Miller--Rabin.
# It is imported rather than reimplemented: a second primality routine would
# be a second authority for the same fact.
from erdos389 import TrialFactorizer, slack
from erdos389 import _is_prime_64 as is_prime_64
from exact_eviction_certificate import bad_window_positions

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "elementary_class_eviction.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def least_prime_above(n: int) -> int:
    """Least prime strictly greater than ``n`` (exact, deterministic)."""
    _require_int("n", n, 0)
    candidate = n + 1
    if candidate <= 2:
        return 2
    if candidate % 2 == 0:
        candidate += 1
    while not is_prime_64(candidate):
        candidate += 2
    return candidate


def factor_modulus(modulus: int) -> dict[int, int]:
    """Exact factorisation of a modulus by trial division."""
    _require_int("modulus", modulus, 1)
    factors: dict[int, int] = {}
    remaining, prime = modulus, 2
    while prime * prime <= remaining:
        while remaining % prime == 0:
            factors[prime] = factors.get(prime, 0) + 1
            remaining //= prime
        prime += 1 if prime == 2 else 2
    if remaining > 1:
        factors[remaining] = factors.get(remaining, 0) + 1
    return factors


def pinned_split(
    residue: int, position: int, modulus: int, factors: dict[int, int]
) -> tuple[int, int, int]:
    """Return ``(d_i, q_i, a_i)`` for the class ``residue mod modulus`` at ``position``.

    ``d_i`` is the part of ``k+i`` that the class pins for every ``k`` in it,
    ``q_i = modulus/d_i`` is the modulus the free cofactor still runs over, and
    ``a_i`` is its residue, which is always coprime to ``q_i``.
    """
    pinned = 1
    for prime, exponent in factors.items():
        value, level = residue + position, 0
        while level < exponent and value % prime == 0:
            value //= prime
            level += 1
        pinned *= prime**level
    cofactor_modulus = modulus // pinned
    cofactor_residue = ((residue + position) // pinned) % cofactor_modulus
    return pinned, cofactor_modulus, cofactor_residue


def fatal_prime_of(m: int, k: int, factorizer: TrialFactorizer) -> dict[str, Any] | None:
    """First violation of (25) in the bad window, found by real factorisation."""
    for position in bad_window_positions(m):
        term = k + position
        for prime, exponent in factorizer.factor(term).items():
            if prime > m and (prime**exponent) ** 2 > 2 * term:
                return {"i": position, "p": prime, "exponent": exponent, "term": term}
    return None


def progression_prime(m: int, top: int) -> int:
    """Least prime with ``p**2 > 2*top`` and ``p > m``.

    Both conditions are needed: (25) only speaks about primes above ``m``, and
    ``p**2 > 2*top`` does not imply ``p > m`` at small scales (``m = 5``,
    ``top = 10`` gives ``p = 5``, which (25) does not reach).
    """
    _require_int("m", m, 1)
    _require_int("top", top, 1)
    return least_prime_above(max(math.isqrt(2 * top), m))


def minimal_progression_length(m: int, start: int, step: int) -> int:
    """Least ``length`` for which (36) holds at this ``start`` and ``step``.

    ``top`` grows with ``length``, so the condition is a fixed point; it is
    reached from below in a few steps because ``P0`` grows like a square root.
    """
    _require_int("m", m, 1)
    _require_int("start", start, 1)
    _require_int("step", step, 1)
    length = progression_prime(m, start + m)
    for _ in range(64):
        needed = progression_prime(m, start + (length - 1) * step + m)
        if needed <= length:
            return length
        length = needed
    raise AssertionError(f"no fixed point for (36) at {(m, start, step)}")


def evict_progression(
    m: int, start: int, step: int, length: int
) -> dict[str, Any] | None:
    """Theorem 15.3: a non-witness inside ``{start + j*step : 0 <= j < length}``.

    Succeeds whenever ``length >= progression_prime(m, top)``; the returned
    record names the position ``i``, the index ``j``, the fatal prime, and the
    factorisation ``k+i = d*p*t`` that Lemma 15.2 consumes.
    """
    _require_int("m", m, 1)
    _require_int("start", start, 1)
    _require_int("step", step, 1)
    _require_int("length", length, 1)
    top = start + (length - 1) * step + m
    prime = progression_prime(m, top)
    if length < prime:
        return None
    factors = factor_modulus(step)
    best: dict[str, Any] | None = None
    for position in bad_window_positions(m):
        pinned, cofactor_modulus, _ = pinned_split(start, position, step, factors)
        head = (start + position) // pinned
        if cofactor_modulus % prime == 0:
            raise AssertionError(f"p divides q, excluded by the proof at {(m, start, step)}")
        index = (-head * pow(cofactor_modulus, -1, prime)) % prime
        if index >= length:
            continue
        k = start + index * step
        term = k + position
        if term % (pinned * prime):
            raise AssertionError(f"construction missed its factorisation at {(m, k)}")
        record = {
            "i": position,
            "j": index,
            "d": pinned,
            "q": cofactor_modulus,
            "p": prime,
            "t": term // (pinned * prime),
            "k": k,
        }
        if best is None or record["k"] < best["k"]:
            best = record
    return best


def evict_initial(m: int, modulus: int, residue: int, factors: dict[int, int]) -> dict[str, Any]:
    """Corollary 15.4(a): a fatal ``k = residue (mod modulus)`` below ``2*M*max(2M, m)``."""
    prime = least_prime_above(max(2 * modulus, m))
    best: dict[str, Any] | None = None
    for position in bad_window_positions(m):
        pinned, cofactor_modulus, cofactor_residue = pinned_split(
            residue, position, modulus, factors
        )
        if cofactor_modulus == 1:
            cofactor = 1
        else:
            cofactor = (cofactor_residue * pow(prime, -1, cofactor_modulus)) % cofactor_modulus
            if cofactor == 0:
                cofactor = cofactor_modulus
        k = pinned * prime * cofactor - position
        if k >= 1 and (best is None or k < best["k"]):
            best = {
                "i": position,
                "d": pinned,
                "q": cofactor_modulus,
                "a": cofactor_residue,
                "p": prime,
                "t": cofactor,
                "k": k,
            }
    if best is None:  # unreachable: p > m forces k = d*p*t - i >= p - m >= 1
        raise AssertionError(f"no eviction constructed at {(m, modulus, residue)}")
    return best


def block_prime(m: int, block_start: int) -> int:
    """Least prime with ``p**2 > 2*(2N+m)`` and ``p > m``."""
    return progression_prime(m, 2 * block_start + m)


def covering_spans(m: int, block_start: int, modulus: int, prime: int) -> list[tuple[int, int]]:
    """Intervals of ``x = c * p^{-1} mod M`` whose class the prime ``p`` evicts.

    A class ``c`` is evicted by ``p`` exactly when some ``z`` has
    ``p*z - i`` in ``[N,2N)`` and ``p*z = c+i (mod M)`` for a window position
    ``i``; the pinned parts of (34) cancel, since only ``p | k+i`` is used.
    In the coordinate ``x = c*p^{-1}`` the condition is
    ``x in Z_i - i*p^{-1}`` with ``Z_i`` an interval of integers, so eviction by
    one prime is a union of ``L`` intervals modulo ``M``.
    """
    _require_int("modulus", modulus, 1)
    inverse = pow(prime, -1, modulus) if modulus > 1 else 0
    spans: list[tuple[int, int]] = []
    for position in bad_window_positions(m):
        low = -(-(block_start + position) // prime)
        high = (2 * block_start - 1 + position) // prime
        if high < low:
            continue
        spans.append((((low - position * inverse) % modulus), high - low + 1))
    return spans


def union_measure(spans: list[tuple[int, int]], modulus: int) -> int:
    """Exact size of a union of cyclic intervals modulo ``modulus``."""
    if any(length >= modulus for _, length in spans):
        return modulus
    pieces: list[tuple[int, int]] = []
    for start, length in spans:
        if start + length <= modulus:
            pieces.append((start, start + length))
        else:
            pieces.append((start, modulus))
            pieces.append((0, start + length - modulus))
    pieces.sort()
    total = 0
    current: tuple[int, int] | None = None
    for start, end in pieces:
        if current is None:
            current = (start, end)
        elif start <= current[1]:
            current = (current[0], max(current[1], end))
        else:
            total += current[1] - current[0]
            current = (start, end)
    return total + (current[1] - current[0] if current is not None else 0)


def covering_prime(
    m: int, block_start: int, modulus: int, *, prime_budget: int = 64
) -> dict[str, Any] | None:
    """First admissible prime whose window union covers every class mod ``M``.

    When one exists the eviction of *every* class from ``[N,2N)`` is a theorem
    for that ``(m, N, M)``, proved by an ``O(L log L)`` computation. The reach is
    ``M <= L * N / P0``, about ``ceil(m/2) * sqrt(N)/2``, since the union of
    ``L`` intervals of length ``N/P0`` cannot exceed that.
    """
    prime = block_prime(m, block_start)
    for tried in range(1, prime_budget + 1):
        if union_measure(covering_spans(m, block_start, modulus, prime), modulus) == modulus:
            return {"p": prime, "primes_tried": tried}
        prime = least_prime_above(prime)
    return None


def evict_block(
    m: int,
    modulus: int,
    residue: int,
    block_start: int,
    factors: dict[int, int],
    *,
    prime_budget: int = 1,
    first_prime: int | None = None,
) -> dict[str, Any] | None:
    """Corollary 15.4(b): a fatal ``k = residue (mod modulus)`` inside ``[N, 2N)``.

    With ``prime_budget == 1`` this is exactly the proved construction, which
    succeeds whenever ``modulus * block_prime(m, N) <= N``.  A larger budget
    scans further primes and measures how far past that threshold the same
    construction keeps working.
    """
    prime = block_prime(m, block_start) if first_prime is None else first_prime
    best: dict[str, Any] | None = None
    for scanned in range(prime_budget):
        for position in bad_window_positions(m):
            pinned, cofactor_modulus, cofactor_residue = pinned_split(
                residue, position, modulus, factors
            )
            step = pinned * prime
            low = -(-(block_start + position) // step)
            if cofactor_modulus == 1:
                cofactor = max(low, 1)
            else:
                target = (cofactor_residue * pow(prime, -1, cofactor_modulus)) % cofactor_modulus
                cofactor = low + ((target - low) % cofactor_modulus)
                if cofactor < 1:
                    cofactor += cofactor_modulus
            k = step * cofactor - position
            if block_start <= k < 2 * block_start and (best is None or k < best["k"]):
                best = {
                    "i": position,
                    "d": pinned,
                    "q": cofactor_modulus,
                    "a": cofactor_residue,
                    "p": prime,
                    "t": cofactor,
                    "k": k,
                    "primes_scanned": scanned + 1,
                }
        if best is not None:
            return best
        prime = least_prime_above(prime)
    return None


def verify_eviction(
    m: int, modulus: int, residue: int, record: dict[str, Any], block_start: int | None = None
) -> None:
    """Exact independent check of one eviction; raises on any failure."""
    k, position, prime = record["k"], record["i"], record["p"]
    term = k + position
    exponent = 0
    value = term
    while value % prime == 0:
        value //= prime
        exponent += 1
    if k % modulus != residue % modulus:
        raise AssertionError(f"eviction left the class at {(m, modulus, residue)}")
    if position not in bad_window_positions(m):
        raise AssertionError(f"eviction used a position outside the window at {(m, k)}")
    if exponent == 0:
        raise AssertionError(f"eviction prime does not divide its term at {(m, k)}")
    if not (prime > m and (prime**exponent) ** 2 > 2 * term):
        raise AssertionError(f"eviction prime is not fatal by (25) at {(m, k)}")
    if slack(m, k, prime) >= 0:
        raise AssertionError(f"fatal prime left nonnegative slack at {(m, k, prime)}")
    if block_start is not None and not block_start <= k < 2 * block_start:
        raise AssertionError(f"eviction left the block at {(m, modulus, residue)}")


def criterion_agreement(max_m: int, max_k: int) -> dict[str, Any]:
    """(25)-fatality by factorisation versus the repository witness verifier."""
    _require_int("max_m", max_m, 1)
    _require_int("max_k", max_k, 1)
    factorizer = TrialFactorizer(max_m + max_k)
    checked = fatal = witnesses = 0
    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            checked += 1
            violation = fatal_prime_of(m, k, factorizer)
            if violation is None:
                continue
            fatal += 1
            if slack(m, k, violation["p"]) >= 0:
                raise AssertionError(f"(25) violation with nonnegative slack at {(m, k)}")
    from erdos389 import is_witness

    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            if is_witness(m, k, factorizer):
                witnesses += 1
                if fatal_prime_of(m, k, factorizer) is not None:
                    raise AssertionError(f"witness carries a fatal prime at {(m, k)}")
    return {
        "checked_pairs": checked,
        "pairs_with_fatal_prime": fatal,
        "witnesses": witnesses,
        "max_m": max_m,
        "max_k": max_k,
        "claim": (
            "every (25) violation found by factorisation has negative Legendre "
            "slack at that prime, and no witness in the rectangle carries one"
        ),
    }


def progression_screen(
    max_m: int, max_step: int, starts: tuple[int, ...]
) -> dict[str, Any]:
    """Theorem 15.3 at its exact minimal length, plus one step below it."""
    _require_int("max_m", max_m, 1)
    _require_int("max_step", max_step, 1)
    checked = short_failures = short_successes = 0
    for start in starts:
        for step in range(1, max_step + 1):
            for m in range(1, max_m + 1):
                length = minimal_progression_length(m, start, step)
                record = evict_progression(m, start, step, length)
                if record is None:
                    raise AssertionError(
                        f"15.3 failed at its own minimal length {(m, start, step, length)}"
                    )
                k = record["k"]
                if not start <= k <= start + (length - 1) * step:
                    raise AssertionError(f"eviction left the progression at {(m, k)}")
                if (k - start) % step:
                    raise AssertionError(f"eviction left the class at {(m, k)}")
                if slack(m, k, record["p"]) >= 0:
                    raise AssertionError(f"progression eviction not fatal at {(m, k)}")
                if record["p"] ** 2 <= 2 * (k + record["i"]):
                    raise AssertionError(f"prime too small for (25) at {(m, k)}")
                checked += 1
                if length > 1:
                    if evict_progression(m, start, step, length - 1) is None:
                        short_failures += 1
                    else:
                        short_successes += 1
    return {
        "exhaustive": True,
        "max_m": max_m,
        "max_step": max_step,
        "starts": list(starts),
        "checked_progressions": checked,
        "one_shorter_no_construction": short_failures,
        "one_shorter_still_worked": short_successes,
        "claim": (
            "at the minimal length of (36) every progression was evicted with an "
            "exact certificate; the one-shorter counts show the hypothesis is the "
            "binding constraint on the construction, not on the conclusion"
        ),
    }


def covering_screen(
    ms: tuple[int, ...],
    block_start: int,
    prime_budget: int,
    classes: int,
    seed: int,
) -> dict[str, Any]:
    """Corollary 15.4(e): bisect the covering reach and verify the certificate."""
    rng = random.Random(seed)
    rows = []
    for m in ms:
        window = len(list(bad_window_positions(m)))
        low, high = 0.50, 0.75
        for _ in range(24):
            mid = (low + high) / 2
            modulus = max(2, int(round(block_start**mid)))
            if covering_prime(m, block_start, modulus, prime_budget=prime_budget):
                low = mid
            else:
                high = mid
        modulus = max(2, int(round(block_start**low)))
        certificate = covering_prime(m, block_start, modulus, prime_budget=prime_budget)
        if certificate is None:
            raise AssertionError(f"bisection kept an uncertified modulus at {(m, low)}")
        factors = factor_modulus(modulus)
        verified = 0
        for _ in range(classes):
            residue = rng.randrange(modulus)
            record = evict_block(
                m, modulus, residue, block_start, factors, first_prime=certificate["p"]
            )
            if record is None:
                raise AssertionError(
                    f"covering certificate missed a class at {(m, modulus, residue)}"
                )
            verify_eviction(m, modulus, residue, record, block_start=block_start)
            verified += 1
        ceiling = math.log(window * (block_start / block_prime(m, block_start))) / math.log(
            block_start
        )
        rows.append(
            {
                "m": m,
                "window_length": window,
                "largest_exponent_certified": low,
                "measure_ceiling_exponent": ceiling,
                "fraction_of_ceiling_used": low / ceiling,
                "M": modulus,
                "certificate": certificate,
                "classes_verified": verified,
            }
        )
    return {
        "block_start": block_start,
        "prime_budget": prime_budget,
        "classes_per_m": classes,
        "seed": seed,
        "claim": (
            "each row's exponent is certified for that (m, N, M) by an exact "
            "O(L log L) covering check, then spot-verified end to end; the "
            "ceiling is the measure bound L*floor(N/P0), which no covering "
            "argument can pass"
        ),
        "rows": rows,
    }


def initial_segment_screen(max_m: int, max_modulus: int) -> dict[str, Any]:
    """Corollary 15.4(a) on every class of every modulus up to ``max_modulus``."""
    _require_int("max_m", max_m, 1)
    _require_int("max_modulus", max_modulus, 1)
    worst_ratio = 0.0
    worst: dict[str, Any] = {}
    checked = 0
    for modulus in range(1, max_modulus + 1):
        factors = factor_modulus(modulus)
        for m in range(1, max_m + 1):
            bound = 2 * modulus * max(2 * modulus, m)
            for residue in range(modulus):
                record = evict_initial(m, modulus, residue, factors)
                verify_eviction(m, modulus, residue, record)
                if record["k"] >= bound:
                    raise AssertionError(f"15a bound violated at {(m, modulus, residue)}")
                checked += 1
                ratio = record["k"] / bound
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    worst = {
                        "m": m,
                        "M": modulus,
                        "c": residue,
                        "k": record["k"],
                        "bound": bound,
                    }
    return {
        "exhaustive": True,
        "max_m": max_m,
        "max_modulus": max_modulus,
        "checked_classes": checked,
        "worst_ratio_to_bound": worst_ratio,
        "worst_case": worst,
    }


def block_screen(max_m: int, max_modulus: int, block_starts: tuple[int, ...]) -> dict[str, Any]:
    """Corollary 15.4(b) on every class of every modulus meeting its condition."""
    _require_int("max_m", max_m, 1)
    _require_int("max_modulus", max_modulus, 1)
    checked = skipped = 0
    for block_start in block_starts:
        _require_int("block_start", block_start, 2)
        for modulus in range(1, max_modulus + 1):
            factors = factor_modulus(modulus)
            for m in range(1, max_m + 1):
                if modulus * block_prime(m, block_start) > block_start:
                    skipped += 1
                    continue
                for residue in range(modulus):
                    record = evict_block(m, modulus, residue, block_start, factors)
                    if record is None:
                        raise AssertionError(
                            f"15b failed inside its hypothesis at "
                            f"{(m, modulus, residue, block_start)}"
                        )
                    verify_eviction(m, modulus, residue, record, block_start=block_start)
                    checked += 1
    return {
        "exhaustive": True,
        "max_m": max_m,
        "max_modulus": max_modulus,
        "block_starts": list(block_starts),
        "checked_classes": checked,
        "skipped_outside_hypothesis": skipped,
        "claim": "every class satisfying M * P0(N,m) <= N was evicted inside [N,2N)",
    }


def threshold_sharpness(
    m: int, block_start: int, classes: int, ratios: tuple[float, ...], seed: int
) -> dict[str, Any]:
    """Where the one-prime construction stops, against the proved threshold."""
    rng = random.Random(seed)
    rows = []
    for ratio in ratios:
        modulus = max(1, int(ratio * math.isqrt(block_start)))
        factors = factor_modulus(modulus)
        inside = modulus * block_prime(m, block_start) <= block_start
        successes = 0
        for _ in range(classes):
            residue = rng.randrange(modulus)
            record = evict_block(m, modulus, residue, block_start, factors)
            if record is None:
                continue
            verify_eviction(m, modulus, residue, record, block_start=block_start)
            successes += 1
        rows.append(
            {
                "ratio_M_over_sqrt_N": ratio,
                "M": modulus,
                "inside_proved_hypothesis": inside,
                "classes_tested": classes,
                "successes": successes,
            }
        )
    inside_rows = [row for row in rows if row["inside_proved_hypothesis"]]
    if any(row["successes"] != row["classes_tested"] for row in inside_rows):
        raise AssertionError("one-prime construction failed inside its own hypothesis")
    return {
        "m": m,
        "block_start": block_start,
        "classes_per_ratio": classes,
        "sampling": "uniform residues, seeded; the rows are a measurement, not a claim",
        "seed": seed,
        "rows": rows,
        "largest_ratio_all_success": max(
            (row["ratio_M_over_sqrt_N"] for row in rows if row["successes"] == classes),
            default=None,
        ),
    }


def barrier_scan(
    m: int,
    block_start: int,
    exponents: tuple[float, ...],
    classes: int,
    prime_budget: int,
    seed: int,
) -> dict[str, Any]:
    """How far past sqrt(N) prime scanning still evicts, and at what cost."""
    rng = random.Random(seed)
    rows = []
    for exponent in exponents:
        modulus = max(2, int(round(block_start**exponent)))
        factors = factor_modulus(modulus)
        successes = 0
        scans = []
        for _ in range(classes):
            residue = rng.randrange(modulus)
            record = evict_block(
                m, modulus, residue, block_start, factors, prime_budget=prime_budget
            )
            if record is None:
                continue
            verify_eviction(m, modulus, residue, record, block_start=block_start)
            successes += 1
            scans.append(record["primes_scanned"])
        window = len(list(bad_window_positions(m)))
        predicted = modulus / (window * (block_start / block_prime(m, block_start)))
        rows.append(
            {
                "log_M_over_log_N": exponent,
                "M": modulus,
                "classes_tested": classes,
                "successes": successes,
                "max_primes_scanned": max(scans) if scans else None,
                "mean_primes_scanned": (sum(scans) / len(scans)) if scans else None,
                "equidistribution_predicted_primes": max(1.0, predicted),
            }
        )
    return {
        "m": m,
        "block_start": block_start,
        "prime_budget": prime_budget,
        "classes_per_exponent": classes,
        "seed": seed,
        "sampling": "uniform residues, seeded; a measurement of the barrier, not a claim",
        "proved_exponent_ceiling": 0.5,
        "cost_model": (
            "If the residues p^{-1} mod q behaved like independent uniform draws, "
            "one prime would place k in the block with probability about "
            "L*(N/p)/q, so the expected number of primes to scan is "
            "q/(L*N/p) = M*P0/(L*N). The measured means are compared against that "
            "prediction: agreement says the obstruction is purely the missing "
            "equidistribution statement, not an arithmetic conspiracy."
        ),
        "rows": rows,
    }


def modulus_family_stats(
    max_m: int, moduli: tuple[int, ...], classes: int, seed: int
) -> dict[str, Any]:
    """Arithmetic of the family {q_i = M/d_i} any analytic route must handle."""
    rng = random.Random(seed)
    rows = []
    for modulus in moduli:
        factors = factor_modulus(modulus)
        trivial = 0
        total = 0
        distinct_counts = []
        smallest_ratio = 1.0
        for m in range(1, max_m + 1):
            positions = list(bad_window_positions(m))
            for _ in range(classes):
                residue = rng.randrange(modulus)
                pinned_parts = [
                    pinned_split(residue, position, modulus, factors)[0]
                    for position in positions
                ]
                total += 1
                if all(part == 1 for part in pinned_parts):
                    trivial += 1
                distinct_counts.append(len({modulus // part for part in pinned_parts}))
                smallest_ratio = min(smallest_ratio, min(1 / part for part in pinned_parts))
        rows.append(
            {
                "M": modulus,
                "sampled_classes": total,
                "fraction_all_pinned_parts_trivial": trivial / total,
                "mean_distinct_moduli_per_class": sum(distinct_counts) / len(distinct_counts),
                "max_distinct_moduli_per_class": max(distinct_counts),
                "smallest_q_over_M": smallest_ratio,
            }
        )
    return {
        "max_m": max_m,
        "classes_per_modulus_and_m": classes,
        "seed": seed,
        "sampling": "uniform residues, seeded",
        "interpretation": (
            "d_i * q_i = M always, so the i-freedom never lowers the product; it "
            "only redistributes it. A route that needs a small modulus cannot get "
            "one by choosing i."
        ),
        "rows": rows,
    }


def analyze(
    *,
    criterion_max_m: int,
    criterion_max_k: int,
    progression_max_m: int,
    progression_max_step: int,
    progression_starts: tuple[int, ...],
    initial_max_m: int,
    initial_max_modulus: int,
    block_max_m: int,
    block_max_modulus: int,
    block_starts: tuple[int, ...],
    covering_ms: tuple[int, ...],
    covering_block: int,
    covering_classes: int,
    sharpness_m: int,
    sharpness_block: int,
    sharpness_classes: int,
    barrier_block: int,
    barrier_classes: int,
    prime_budget: int,
    family_moduli: tuple[int, ...],
    family_classes: int,
    seed: int,
) -> dict[str, Any]:
    return {
        "status": "PASS_EXACT_ELEMENTARY_CLASS_EVICTION",
        "criterion_agreement": criterion_agreement(criterion_max_m, criterion_max_k),
        "progression": progression_screen(
            progression_max_m, progression_max_step, progression_starts
        ),
        "initial_segment": initial_segment_screen(initial_max_m, initial_max_modulus),
        "block": block_screen(block_max_m, block_max_modulus, block_starts),
        "covering": covering_screen(
            covering_ms, covering_block, prime_budget, covering_classes, seed
        ),
        "threshold_sharpness": threshold_sharpness(
            sharpness_m,
            sharpness_block,
            sharpness_classes,
            (0.10, 0.25, 0.40, 0.49, 0.51, 0.60, 0.80, 1.00, 1.50,
             2.00, 3.00, 5.00, 8.00, 14.00, 20.00, 40.00),
            seed,
        ),
        "barrier": barrier_scan(
            sharpness_m,
            barrier_block,
            (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95),
            barrier_classes,
            prime_budget,
            seed,
        ),
        "modulus_family": modulus_family_stats(
            block_max_m, family_moduli, family_classes, seed
        ),
        "claim_boundary": (
            "Theorem 15a and Theorem 15b are proved for every class; the screens "
            "exhaust declared rectangles and are regression evidence, not the proof. "
            "The sharpness and barrier blocks are seeded measurements of where the "
            "one-prime construction stops, and are labelled as measurements."
        ),
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "elementary_class_eviction.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "exact_eviction_certificate.py": sha256_file(ROOT / "exact_eviction_certificate.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--criterion-max-m", type=int, default=12)
    parser.add_argument("--criterion-max-k", type=int, default=300)
    parser.add_argument("--progression-max-m", type=int, default=12)
    parser.add_argument("--progression-max-step", type=int, default=24)
    parser.add_argument(
        "--progression-starts", type=int, nargs="+", default=[1, 10**4, 10**6, 10**9]
    )
    parser.add_argument("--initial-max-m", type=int, default=16)
    parser.add_argument("--initial-max-modulus", type=int, default=90)
    parser.add_argument("--block-max-m", type=int, default=20)
    parser.add_argument("--block-max-modulus", type=int, default=48)
    parser.add_argument("--block-starts", type=int, nargs="+", default=[10**6, 10**12])
    parser.add_argument("--covering-ms", type=int, nargs="+", default=[5, 13, 27, 51])
    parser.add_argument("--covering-block", type=int, default=10**12)
    parser.add_argument("--covering-classes", type=int, default=40)
    parser.add_argument("--sharpness-m", type=int, default=27)
    parser.add_argument("--sharpness-block", type=int, default=10**12)
    parser.add_argument("--sharpness-classes", type=int, default=200)
    parser.add_argument("--barrier-block", type=int, default=10**12)
    parser.add_argument("--barrier-classes", type=int, default=100)
    parser.add_argument("--prime-budget", type=int, default=2000)
    parser.add_argument("--family-moduli", type=int, nargs="+", default=[720, 5040, 999983])
    parser.add_argument("--family-classes", type=int, default=40)
    parser.add_argument("--seed", type=int, default=389)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        criterion_max_m=args.criterion_max_m,
        criterion_max_k=args.criterion_max_k,
        progression_max_m=args.progression_max_m,
        progression_max_step=args.progression_max_step,
        progression_starts=tuple(args.progression_starts),
        initial_max_m=args.initial_max_m,
        initial_max_modulus=args.initial_max_modulus,
        block_max_m=args.block_max_m,
        block_max_modulus=args.block_max_modulus,
        block_starts=tuple(args.block_starts),
        covering_ms=tuple(args.covering_ms),
        covering_block=args.covering_block,
        covering_classes=args.covering_classes,
        sharpness_m=args.sharpness_m,
        sharpness_block=args.sharpness_block,
        sharpness_classes=args.sharpness_classes,
        barrier_block=args.barrier_block,
        barrier_classes=args.barrier_classes,
        prime_budget=args.prime_budget,
        family_moduli=tuple(args.family_moduli),
        family_classes=args.family_classes,
        seed=args.seed,
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
                "progressions": result["progression"]["checked_progressions"],
                "covering": [
                    {
                        "m": row["m"],
                        "theta": round(row["largest_exponent_certified"], 4),
                        "of_ceiling": round(row["fraction_of_ceiling_used"], 4),
                    }
                    for row in result["covering"]["rows"]
                ],
                "initial_classes": result["initial_segment"]["checked_classes"],
                "block_classes": result["block"]["checked_classes"],
                "largest_ratio_all_success": result["threshold_sharpness"][
                    "largest_ratio_all_success"
                ],
                "barrier_rows": [
                    {
                        "exponent": row["log_M_over_log_N"],
                        "successes": row["successes"],
                        "max_primes": row["max_primes_scanned"],
                    }
                    for row in result["barrier"]["rows"]
                ],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
