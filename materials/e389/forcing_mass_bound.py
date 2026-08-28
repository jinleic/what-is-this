#!/usr/bin/env python3
"""Price the modulus that congruence-forced powersmooth mass costs.

Call a pair ``(M, k)`` *size-forcing* when the part of every bad-window term
that ``k mod M`` pins already meets the powersmooth threshold of (25), so no
unforced cofactor can carry a fatal prime whatever it turns out to be. Section
14 proves that for such a pair

    log M >= (L/4) log(k/2) - (L/2)(log L + 1),      L = ceil(m/2).

This prices one pair; it is not a claim that some class certifies unbounded
``k``. That reading is already settled by the Dirichlet argument of Section 8,
which empties the hypothesis. The reading with content is a dyadic range,
where a size-forcing modulus must exceed the range itself. Small moduli are
excluded inside a range only in the Siegel-Walfisz regime, or under GRH;
Linnik bounds the least prime of the cofactor progression, so it places a fatal
k before the range, not inside it. The interval between is open; see the scope
note in Section 14.

Three things are checked exactly here: the window counting lemma that supplies
the mass bound, the Mertens constant used for the small-prime correction, and
the resulting inequality against explicit size-forcing moduli built from real
factorizations of published witness windows. Nothing is sampled.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import time
from math import isqrt
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import TrialFactorizer, primes_up_to
from exact_eviction_certificate import bad_window_positions

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "forcing_mass_bound.json"
DEFAULT_WITNESS_CSV = ROOT / "data" / "oeis_a375071.csv"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def window_length(m: int) -> int:
    _require_int("m", m, 1)
    return (m + 1) // 2


def window_multiplicity_sum(m: int, p: int, exponent: int, residue: int) -> int:
    """Return ``sum_{j<=a} #{i in I_m : p^j | residue + i}`` exactly."""
    _require_int("m", m, 1)
    _require_int("p", p, 2)
    _require_int("exponent", exponent, 1)
    _require_int("residue", residue, 0)
    total = 0
    power = 1
    for _ in range(exponent):
        power *= p
        total += sum(1 for i in bad_window_positions(m) if (residue + i) % power == 0)
    return total


def counting_lemma_bound(m: int, p: int, exponent: int) -> float:
    """Return the Section 14 bound ``a + L/(p-1)`` on the multiplicity sum."""
    return exponent + window_length(m) / (p - 1)


def verify_counting_lemma(max_m: int, primes: tuple[int, ...], max_exponent: int) -> dict[str, Any]:
    """Check the multiplicity bound for every residue on a declared rectangle."""
    _require_int("max_m", max_m, 1)
    _require_int("max_exponent", max_exponent, 1)
    checks = 0
    worst = None
    for m in range(1, max_m + 1):
        for p in primes:
            for exponent in range(1, max_exponent + 1):
                bound = counting_lemma_bound(m, p, exponent)
                for residue in range(p**exponent):
                    observed = window_multiplicity_sum(m, p, exponent, residue)
                    if observed > bound:
                        raise AssertionError(
                            f"counting lemma failed at {(m, p, exponent, residue)}: "
                            f"{observed} > {bound}"
                        )
                    slack = bound - observed
                    if worst is None or slack < worst["slack"]:
                        worst = {
                            "m": m,
                            "p": p,
                            "exponent": exponent,
                            "residue": residue,
                            "observed": observed,
                            "bound": bound,
                            "slack": slack,
                        }
                    checks += 1
    return {"checked_residues": checks, "tightest_case": worst}


def verify_mertens_constant(limit: int) -> dict[str, Any]:
    """Check ``sum_{p<=P} log p/(p-1) <= log P + 1`` for every prime ``P<=limit``."""
    _require_int("limit", limit, 2)
    running = 0.0
    worst = None
    for p in primes_up_to(limit):
        running += math.log(p) / (p - 1)
        margin = math.log(p) + 1.0 - running
        if worst is None or margin < worst["margin"]:
            worst = {"prime": p, "sum": running, "margin": margin}
    if worst["margin"] < 0:
        raise AssertionError("Mertens correction bound failed")
    return {"limit": limit, "tightest": worst, "final_sum": running}


def modulus_log_lower_bound(m: int, k: int) -> float:
    """Return the Section 14 lower bound on ``log M``."""
    _require_int("m", m, 1)
    _require_int("k", k, 2)
    length = window_length(m)
    return (length / 4.0) * math.log(k / 2.0) - (length / 2.0) * (
        math.log(length) + 1.0
    )


def exceeds_target_threshold(m: int) -> int | None:
    """Smallest ``k`` for which the bound already forces ``M > k``, if any."""
    length = window_length(m)
    if length <= 4:
        return None
    low, high = 3, 10**60
    while low < high:
        mid = (low + high) // 2
        if modulus_log_lower_bound(m, mid) > math.log(mid):
            high = mid
        else:
            low = mid + 1
    return low


def greedy_size_forcing_modulus(
    m: int, k: int, factorizer: TrialFactorizer
) -> dict[str, Any]:
    """Build a modulus making ``(M, k)`` size-forcing.

    Greedy over the largest prime powers of each term, sharing exponents across
    the window. The result is an upper bound for the least size-forcing
    modulus, so comparing it with :func:`modulus_log_lower_bound` brackets the
    truth and shows the hypothesis of (33) is satisfiable.
    """
    _require_int("m", m, 1)
    _require_int("k", k, 2)
    needed: dict[int, int] = {}
    per_term = []
    for i in bad_window_positions(m):
        term = k + i
        target = isqrt(term // 2) + 1
        factors = sorted(
            factorizer.factor(term).items(), key=lambda item: -(item[0] ** item[1])
        )
        forced = 1
        used: dict[int, int] = {}
        for prime, exponent in factors:
            if forced > target:
                break
            for _ in range(exponent):
                forced *= prime
                used[prime] = used.get(prime, 0) + 1
                if forced > target:
                    break
        if forced <= target:
            raise AssertionError(
                f"term {term} has no divisor above sqrt(w/2) inside its factorization"
            )
        for prime, exponent in used.items():
            needed[prime] = max(needed.get(prime, 0), exponent)
        per_term.append(
            {"position": i, "term": term, "forced_divisor_decimal_digits": len(str(forced))}
        )
    modulus = 1
    for prime, exponent in sorted(needed.items()):
        modulus *= prime**exponent
    return {
        "modulus_decimal_digits": len(str(modulus)),
        "log_modulus": math.log(modulus),
        "distinct_primes": len(needed),
        "terms": per_term,
    }


def read_witness_rows(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [(int(row["m"]), int(row["k"])) for row in csv.DictReader(stream)]


def analyze(
    *,
    witness_csv: Path,
    counting_max_m: int,
    counting_primes: tuple[int, ...],
    counting_max_exponent: int,
    mertens_limit: int,
    greedy_max_k: int,
) -> dict[str, Any]:
    witness_csv = witness_csv.resolve()
    rows = read_witness_rows(witness_csv)
    factorizer = TrialFactorizer(10**14)

    published: list[dict[str, Any]] = []
    bracketed = 0
    for m, k in rows:
        if m < 1 or k < 2:
            continue
        bound = modulus_log_lower_bound(m, k)
        record = {
            "m": m,
            "k": k,
            "k_decimal_digits": len(str(k)),
            "window_length": window_length(m),
            "log_modulus_lower_bound": bound,
            "modulus_decimal_digits_lower_bound": (
                math.floor(bound / math.log(10)) if bound > 0 else 0
            ),
            "bound_exceeds_k": bound > math.log(k),
        }
        if k <= greedy_max_k:
            greedy = greedy_size_forcing_modulus(m, k, factorizer)
            if greedy["log_modulus"] < bound:
                raise AssertionError(
                    f"greedy size-forcing modulus below the proved bound at {(m, k)}"
                )
            record["greedy_size_forcing_modulus"] = greedy
            record["bracket_ratio"] = (
                greedy["log_modulus"] / bound if bound > 0 else None
            )
            bracketed += 1
        published.append(record)

    thresholds = [
        {
            "m": m,
            "window_length": window_length(m),
            "least_k_forcing_modulus_above_k": exceeds_target_threshold(m),
        }
        for m in range(1, 31)
    ]

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_FORCING_MASS_BOUND",
        "claim_boundary": (
            "The counting lemma is verified for every residue class on its "
            "declared rectangle, and the Mertens correction is verified at "
            "every prime below its declared limit; both are exact finite "
            "checks, not samples. The greedy moduli are upper bounds for the "
            "least size-forcing modulus, computed from real factorizations, "
            "and are only required to dominate the proved lower bound. The "
            "bound prices one pair (M, k) of congruence-forced powersmooth "
            "mass. It is not a statement about a class certifying unbounded "
            "k: Section 8 settles that reading by Dirichlet and empties this "
            "hypothesis. It says nothing about small moduli, which need a "
            "separate argument, and it proves nothing about the existence or "
            "the rarity of witnesses."
        ),
        "parameters": {
            "counting_max_m": counting_max_m,
            "counting_primes": list(counting_primes),
            "counting_max_exponent": counting_max_exponent,
            "mertens_limit": mertens_limit,
            "greedy_max_k": greedy_max_k,
        },
        "counting_lemma": verify_counting_lemma(
            counting_max_m, counting_primes, counting_max_exponent
        ),
        "mertens_correction": verify_mertens_constant(mertens_limit),
        "published_witnesses": published,
        "bracketed_witness_count": bracketed,
        "threshold_by_m": thresholds,
        "input": {
            "witness_csv_logical_name": witness_csv.name,
            "witness_csv_sha256": sha256_file(witness_csv),
            "published_witness_rows": len(rows),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "forcing_mass_bound.py": sha256_file(Path(__file__)),
            "exact_eviction_certificate.py": sha256_file(
                ROOT / "exact_eviction_certificate.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counting-max-m", type=int, default=30)
    parser.add_argument("--counting-primes", type=int, nargs="+", default=[2, 3, 5, 7, 11])
    parser.add_argument("--counting-max-exponent", type=int, default=4)
    parser.add_argument("--mertens-limit", type=int, default=2_000_000)
    parser.add_argument("--greedy-max-k", type=int, default=10**9)
    parser.add_argument("--witness-csv", type=Path, default=DEFAULT_WITNESS_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        witness_csv=args.witness_csv,
        counting_max_m=args.counting_max_m,
        counting_primes=tuple(args.counting_primes),
        counting_max_exponent=args.counting_max_exponent,
        mertens_limit=args.mertens_limit,
        greedy_max_k=args.greedy_max_k,
    )
    result["resources"] = {
        "wall_seconds": time.perf_counter() - started_wall,
        "process_cpu_seconds": time.process_time() - started_cpu,
        "processes": 1,
    }
    atomic_write_json(args.output, result)
    largest = max(result["published_witnesses"], key=lambda row: row["k"])
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "checked_residues": result["counting_lemma"]["checked_residues"],
                "mertens_margin": result["mertens_correction"]["tightest"]["margin"],
                "bracketed_witness_count": result["bracketed_witness_count"],
                "largest_witness": {
                    "m": largest["m"],
                    "k_decimal_digits": largest["k_decimal_digits"],
                    "modulus_decimal_digits_lower_bound": largest[
                        "modulus_decimal_digits_lower_bound"
                    ],
                },
                "witnesses_with_bound_above_k": sum(
                    1 for row in result["published_witnesses"] if row["bound_exceeds_k"]
                ),
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
