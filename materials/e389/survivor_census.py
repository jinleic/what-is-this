#!/usr/bin/env python3
"""Census the survivors that a Brun--Titchmarsh constant would have to bound.

Section 15 evicts every class of modulus ``M <= sqrt(N)/2`` from ``[N,2N)`` by
construction. Above that the literature offers a second, entirely different
route (Section 16): count the *survivors*

    S(N, q, a) = #{ n in [N,2N) : n = a (mod q), no prime p | n has p^2 > 2n },

and prove ``S(N,q,a) < #{n in [N,2N): n = a (mod q)}``. Then some member of the
class carries a fatal prime and the class is evicted, with no prime
localisation at all. Balog--Pomerance (1992) and Shiu (1980) both bound
``S(N,q,a)`` for an individual modulus in exactly this range, but with an
unspecified absolute constant; the Dickman heuristic says the truth is
``rho(2) = 1 - log 2 = 0.3069`` per class, so the needed inequality has a
factor ``1/rho(2) = 3.26`` of room against the class size.

This producer measures that margin exactly. It builds the fatality mask of
``[N,2N)`` by sieve — ``n`` is fatal iff some prime ``p | n`` has ``p^2 > 2n`` —
and then, for each tested modulus, reports every class's survivor count. Three
things are extracted:

  * the overall survivor density, against ``rho(2)``;
  * the worst class: ``max_a S(N,q,a) / (class size)``, the quantity an
    effective Balog--Pomerance or Shiu constant must push below ``1``;
  * the number of classes with no fatal member at all. Zero such classes is a
    finite certificate that every class mod ``q`` is evicted from ``[N,2N)``,
    unconditionally and far above ``sqrt(N)``.

Moduli are tested in three shapes — prime, smooth, and highly composite — since
the only known structural obstruction to per-class upper bounds (concentration
into a subgroup of index bounded in terms of the least non-residue) is
shape-sensitive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "survivor_census.json"
DICKMAN_RHO_2 = 1.0 - math.log(2.0)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def primes_below(limit: int) -> list[int]:
    """Exact sieve of Eratosthenes with byte-slice striking."""
    _require_int("limit", limit, 2)
    flags = bytearray(b"\x01") * limit
    flags[0:2] = b"\x00\x00"
    for candidate in range(2, math.isqrt(limit - 1) + 1):
        if flags[candidate]:
            start = candidate * candidate
            flags[start:limit:candidate] = bytes(
                len(range(start, limit, candidate))
            )
    return [index for index in range(limit) if flags[index]]


def fatality_mask(block_start: int, max_m: int) -> bytearray:
    """``mask[n - N] = 1`` iff ``n`` violates (25) for every ``m <= max_m``.

    (25) fails at ``n`` when some prime power ``p**v`` exactly dividing ``n``
    with ``p > m`` satisfies ``(p**v)**2 > 2n``. Prime *powers* must be included:
    at ``q`` prime with ``N < q**2 < 2N`` the member ``n = q**2`` of the class
    ``0 mod q`` is fatal through ``q**2`` while no prime factor of any member is.
    Only prime powers above ``sqrt(2N)`` can qualify, and for each the
    qualifying multiples are exactly those below ``(p**v)**2/2``; ``p**v | n``
    already implies ``p**{v_p(n)} >= p**v``, so taking a maximum over powers is
    unnecessary and one byte-slice per prime power suffices.
    """
    _require_int("block_start", block_start, 4)
    _require_int("max_m", max_m, 1)
    top = 2 * block_start
    mask = bytearray(block_start)
    for prime in primes_below(top):
        if prime <= max_m:
            continue
        power = prime
        while power < top:
            if power * power > 2 * block_start:
                limit = min(top, (power * power + 1) // 2)
                first = ((block_start + power - 1) // power) * power
                if first < limit:
                    lo = first - block_start
                    hi = limit - block_start
                    mask[lo:hi:power] = b"\x01" * len(range(lo, hi, power))
            power *= prime
    return mask


def window_survivors(has_fatal: list[bool], modulus: int, m: int) -> int:
    """Number of k-classes mod ``modulus`` that no window position can evict.

    Eviction needs only one position, so ``c`` survives exactly when every
    ``c+i``, ``i`` in the bad window, is a fatal-free n-class. The window is an
    interval of ``L = ceil(m/2)`` consecutive integers, so survivors are the
    starts of cyclic runs of ``L`` fatal-free residues.
    """
    _require_int("m", m, 1)
    positions = range(m // 2 + 1, m + 1)
    length = len(positions)
    if not any(has_fatal):
        return modulus
    if length > modulus:
        # The window covers every residue at least once, so one fatal class
        # evicts every k-class.
        return 0
    # Survivors are the starts of cyclic runs of ``length`` fatal-free residues.
    # Walking from a fatal residue keeps every run intact across the wrap.
    start = has_fatal.index(True)
    survivors = 0
    run = 0
    for step in range(1, modulus + 1):
        if has_fatal[(start + step) % modulus]:
            if run >= length:
                survivors += run - length + 1
            run = 0
        else:
            run += 1
    return survivors


def class_census(mask: bytearray, block_start: int, modulus: int) -> dict[str, Any]:
    """Survivor and fatal counts for every class modulo ``modulus``."""
    _require_int("modulus", modulus, 1)
    worst_ratio = -1.0
    worst_class = -1
    empty_classes = 0
    has_fatal = [False] * modulus
    fatal_total = 0
    smallest_fatal = None
    for residue in range(modulus):
        offset = (residue - block_start) % modulus
        members = len(range(offset, block_start, modulus))
        if members == 0:
            continue
        fatal = sum(mask[offset:block_start:modulus])
        fatal_total += fatal
        if smallest_fatal is None or fatal < smallest_fatal:
            smallest_fatal = fatal
        has_fatal[residue] = fatal > 0
        if fatal == 0:
            empty_classes += 1
        ratio = (members - fatal) / members
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_class = residue
    return {
        "M": modulus,
        "classes": modulus,
        "min_fatal_per_class": smallest_fatal,
        "classes_with_no_fatal_member": empty_classes,
        "every_class_evicted": empty_classes == 0,
        "worst_class": worst_class,
        "worst_class_survivor_fraction": worst_ratio,
        "overall_survivor_fraction": (block_start - fatal_total) / block_start,
        "has_fatal": has_fatal,
    }


def modulus_menu(block_start: int, exponents: tuple[float, ...]) -> list[dict[str, Any]]:
    """Prime, smooth, and highly composite moduli near each target exponent."""
    smooth_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    menu: list[dict[str, Any]] = []
    for exponent in exponents:
        target = block_start**exponent
        prime = 2
        for candidate in range(int(target), int(target) + 10_000):
            if candidate < 2:
                continue
            if all(candidate % p for p in range(2, math.isqrt(candidate) + 1)):
                prime = candidate
                break
        smooth = 1
        index = 0
        while index < 400:
            base = smooth_primes[index % len(smooth_primes)]
            if smooth * base > target:
                index += 1
                if index % len(smooth_primes) == 0 and smooth * 2 > target:
                    break
                continue
            smooth *= base
            index += 1
        composite = 1
        for base in smooth_primes:
            if composite * base > target:
                break
            composite *= base
        while composite * 2 <= target:
            composite *= 2
        menu.append(
            {
                "exponent": exponent,
                "prime": prime,
                "smooth": max(2, smooth),
                "primorial": max(2, composite),
            }
        )
    return menu


def analyze(
    *,
    block_start: int,
    exponents: tuple[float, ...],
    max_m: int,
    window_ms: tuple[int, ...],
) -> dict[str, Any]:
    mask = fatality_mask(block_start, max_m)
    overall = sum(mask) / block_start
    menu = modulus_menu(block_start, exponents)
    rows = []
    for entry in menu:
        for shape in ("prime", "smooth", "primorial"):
            modulus = entry[shape]
            if modulus >= block_start:
                continue
            census = class_census(mask, block_start, modulus)
            has_fatal = census.pop("has_fatal")
            census["window_survivors_by_m"] = {
                str(m): window_survivors(has_fatal, modulus, m) for m in window_ms
            }
            census["every_k_class_evicted_by_m"] = {
                str(m): count == 0
                for m, count in zip(
                    window_ms, census["window_survivors_by_m"].values()
                )
            }
            census["shape"] = shape
            census["target_exponent"] = entry["exponent"]
            census["log_M_over_log_N"] = math.log(modulus) / math.log(block_start)
            rows.append(census)
    surviving = [row for row in rows if not row["every_class_evicted"]]
    certified = [row for row in rows if row["every_class_evicted"]]
    per_m = {}
    for m in window_ms:
        good = [
            row["log_M_over_log_N"]
            for row in rows
            if row["every_k_class_evicted_by_m"][str(m)]
        ]
        bad = [
            row["log_M_over_log_N"]
            for row in rows
            if not row["every_k_class_evicted_by_m"][str(m)]
        ]
        per_m[str(m)] = {
            "largest_exponent_certified": max(good, default=None),
            "smallest_exponent_failed": min(bad, default=None),
        }
    return {
        "status": "PASS_EXACT_SURVIVOR_CENSUS",
        "block_start": block_start,
        "block": [block_start, 2 * block_start],
        "max_m": max_m,
        "parameters": {
            "block_start": block_start,
            "target_exponents": list(exponents),
            "max_m": max_m,
            "window_ms": list(window_ms),
        },
        "fatal_density": overall,
        "survivor_density": 1.0 - overall,
        "dickman_rho_2": DICKMAN_RHO_2,
        "survivor_density_over_rho2": (1.0 - overall) / DICKMAN_RHO_2,
        "window_ms": list(window_ms),
        "multi_position_by_m": per_m,
        "rows": rows,
        "single_position_largest_exponent_all_classes_evicted": max(
            (row["log_M_over_log_N"] for row in certified), default=None
        ),
        "single_position_smallest_exponent_with_a_surviving_class": min(
            (row["log_M_over_log_N"] for row in surviving), default=None
        ),
        "worst_class_survivor_fraction_overall": max(
            row["worst_class_survivor_fraction"] for row in rows
        ),
        "claim_boundary": (
            "Every row is an exact count over the whole block, not a sample. A "
            "row with every_class_evicted = true is a finite certificate for "
            "that (N, M) and every m <= max_m; it is not a theorem for "
            "other N. The survivor fractions are the quantity that an effective "
            "Balog-Pomerance or Shiu constant must bound below 1."
        ),
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "survivor_census.py": sha256_file(Path(__file__)),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-start", type=int, default=10**7)
    parser.add_argument("--max-m", type=int, default=60)
    parser.add_argument(
        "--window-ms", type=int, nargs="+", default=[1, 3, 5, 13, 27, 51]
    )
    parser.add_argument(
        "--exponents",
        type=float,
        nargs="+",
        default=[
            0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99
        ],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        block_start=args.block_start,
        exponents=tuple(args.exponents),
        max_m=args.max_m,
        window_ms=tuple(args.window_ms),
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
                "survivor_density": result["survivor_density"],
                "survivor_density_over_rho2": result["survivor_density_over_rho2"],
                "multi_position_by_m": result["multi_position_by_m"],
                "worst_class_survivor_fraction_overall": result[
                    "worst_class_survivor_fraction_overall"
                ],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
