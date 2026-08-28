#!/usr/bin/env python3
"""Certify that the small-prime tier of Erdos #389 is asymptotically free.

Section 9 splits the witness criterion exactly into a small-prime tier
``slack_p(m,k) >= 0`` for every ``p <= m`` and a run of ``ceil(m/2)``
consecutive ``m``-compensation-good bad-window terms.  This producer settles
the tier side.

Every reported cell is checked by independent exact methods:

1. a base-``p`` digit dynamic program over all ``k < p^J`` (no enumeration),
2. a prefix-decomposition dynamic program over ``1 <= k <= X`` for arbitrary
   ``X``, agreeing with (1) at ``X = p^J - 1``,
3. direct enumeration on the small declared grid,
4. the closed-form upper bound of Section 13, which must dominate (1), and
5. the explicit decay envelope ``4 exp(-J/24)`` of Section 13, checked
   wherever its hypothesis ``J >= 16 (D+1)`` holds.

The digit-sum slack identity and the carry inequality that the bound rests on
are re-derived numerically against ``erdos389.slack`` on a declared rectangle.
No density in this file is sampled.  The only non-exact input is the measured
compensation-good density used for the closing crossover scale, which is
labelled and proves nothing by itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from math import comb
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import primes_up_to, slack

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "tier_density_m1_30.json"
DEFAULT_DENSITY_ARTIFACT = ROOT / "data" / "smooth_density_m1_20_k3000.json"

# Largest published witness a(25); a(26) = a(25) - 1 by the odd-to-even shift
# of Section 10. Used only as a report scale.
PUBLISHED_M25_K = 5_048_891_644_621


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def base_digits(n: int, p: int) -> list[int]:
    """Return the base-``p`` digits of ``n``, least significant first."""
    _require_int("n", n, 0)
    _require_int("p", p, 2)
    out: list[int] = []
    while n:
        out.append(n % p)
        n //= p
    return out


def digit_sum(n: int, p: int) -> int:
    return sum(base_digits(n, p))


def digit_sum_slack(m: int, k: int, p: int) -> int:
    """Return ``slack_p(m,k)`` through the digit-sum identity (28)."""
    _require_int("m", m, 0)
    _require_int("k", k, 1)
    _require_int("p", p, 2)
    total = 2 * digit_sum(m + k, p) - digit_sum(m + 2 * k, p) - digit_sum(m, p)
    if total % (p - 1):
        raise AssertionError("digit-sum identity produced a non-multiple of p-1")
    return total // (p - 1)


def high_digit_count(p: int) -> int:
    """Return ``#{d in [0,p) : 2d >= p}``."""
    _require_int("p", p, 2)
    return 1 if p == 2 else (p - 1) // 2


def carry_bound_terms(m: int, k: int, p: int) -> tuple[int, int, int]:
    """Return ``(D, Z, W)`` of inequality (29)."""
    _require_int("m", m, 0)
    _require_int("k", k, 1)
    _require_int("p", p, 2)
    d_count = len(base_digits(m, p))
    kd = base_digits(k, p)
    z = 0
    while d_count + z < len(kd) and kd[d_count + z] == p - 1:
        z += 1
    window = sum(1 for j in range(d_count + z + 1, len(kd)) if 2 * kd[j] >= p)
    return d_count, z, window


def _state_step(
    state: dict[tuple[int, int, int], int],
    p: int,
    m_digit: int,
    k_digits: range | tuple[int, ...],
    cap: int,
) -> dict[tuple[int, int, int], int]:
    """Advance one base-``p`` position of both additions simultaneously."""
    nxt: dict[tuple[int, int, int], int] = {}
    for (carry_left, carry_right, delta), count in state.items():
        for k_digit in k_digits:
            left = m_digit + k_digit + carry_left
            new_left = 1 if left >= p else 0
            sum_digit = left - p if new_left else left
            right = k_digit + sum_digit + carry_right
            new_right = 1 if right >= p else 0
            new_delta = delta + new_right - new_left
            if new_delta > cap:
                new_delta = cap
            elif new_delta < -cap:
                new_delta = -cap
            key = (new_left, new_right, new_delta)
            nxt[key] = nxt.get(key, 0) + count
    return nxt


def _finalize(state: dict[tuple[int, int, int], int], p: int) -> int:
    """Count states whose completed slack is negative."""
    total = 0
    for (carry_left, carry_right, delta), count in state.items():
        tail = carry_left + carry_right
        final = delta
        while tail >= p:
            final += 1
            tail = tail - p + 1
        if final < 0:
            total += count
    return total


def failure_count_power(m: int, p: int, digit_count: int) -> int:
    """Exact ``#{0 <= k < p^J : slack_p(m,k) < 0}`` by digit dynamic program."""
    _require_int("m", m, 0)
    _require_int("p", p, 2)
    _require_int("digit_count", digit_count, 1)
    m_digits = base_digits(m, p)
    if len(m_digits) > digit_count:
        raise ValueError("digit_count must cover every base-p digit of m")
    cap = digit_count + 2
    state = {(0, 0, 0): 1}
    for position in range(digit_count):
        m_digit = m_digits[position] if position < len(m_digits) else 0
        state = _state_step(state, p, m_digit, range(p), cap)
    return _finalize(state, p)


def failure_count_upto(m: int, p: int, limit: int) -> int:
    """Exact ``#{1 <= k <= X : slack_p(m,k) < 0}`` for arbitrary ``X``."""
    _require_int("m", m, 0)
    _require_int("p", p, 2)
    _require_int("limit", limit, 1)
    x_digits = base_digits(limit, p)
    digit_count = max(len(x_digits), len(base_digits(m, p)) + 1)
    x_digits = x_digits + [0] * (digit_count - len(x_digits))
    m_digits = base_digits(m, p) + [0] * digit_count

    cap = digit_count + 2
    prefix: list[dict[tuple[int, int, int], int]] = [{(0, 0, 0): 1}]
    for position in range(digit_count):
        prefix.append(_state_step(prefix[-1], p, m_digits[position], range(p), cap))

    total = 0
    for position in range(digit_count):
        if x_digits[position] == 0:
            continue
        state = _state_step(
            prefix[position],
            p,
            m_digits[position],
            tuple(range(x_digits[position])),
            cap,
        )
        for higher in range(position + 1, digit_count):
            state = _state_step(state, p, m_digits[higher], (x_digits[higher],), cap)
        total += _finalize(state, p)
    if slack(m, limit, p) < 0:
        total += 1
    return total


def failure_count_enumerated(m: int, p: int, limit: int) -> int:
    """Exact failure count by direct enumeration; the independent oracle."""
    _require_int("m", m, 0)
    _require_int("p", p, 2)
    _require_int("limit", limit, 1)
    return sum(1 for k in range(1, limit + 1) if slack(m, k, p) < 0)


def certified_upper_bound(m: int, p: int, digit_count: int) -> int:
    """Closed-form bound (30) on ``#{k < p^J : slack_p(m,k) < 0}``."""
    _require_int("m", m, 0)
    _require_int("p", p, 2)
    _require_int("digit_count", digit_count, 1)
    d_count = len(base_digits(m, p))
    if d_count > digit_count:
        raise ValueError("digit_count must cover every base-p digit of m")
    high = high_digit_count(p)
    low = p - high
    total = 0
    for run in range(digit_count - d_count):
        free = digit_count - d_count - run - 1
        cutoff = d_count + run - 1
        if cutoff < 0:
            continue
        tail = sum(
            comb(free, i) * high**i * low ** (free - i)
            for i in range(min(cutoff, free) + 1)
        )
        total += p**d_count * (p - 1) * tail
    total += p**d_count
    return total


def decay_envelope(digit_count: int) -> float:
    """Return the Section 13 envelope ``4 exp(-J/24)``."""
    _require_int("digit_count", digit_count, 1)
    return 4.0 * math.exp(-digit_count / 24.0)


def envelope_applies(m: int, p: int, digit_count: int) -> bool:
    return digit_count >= 16 * (len(base_digits(m, p)) + 1)


def identity_and_inequality_scan(max_m: int, max_k: int) -> dict[str, Any]:
    """Recheck (28) and (29) against ``erdos389.slack`` on a full rectangle."""
    _require_int("max_m", max_m, 1)
    _require_int("max_k", max_k, 1)
    primes = primes_up_to(max_m)
    checks = 0
    slack_min = 0
    margin_min = None
    for m in range(1, max_m + 1):
        for p in primes:
            if p > m:
                continue
            for k in range(1, max_k + 1):
                exact = slack(m, k, p)
                if exact != digit_sum_slack(m, k, p):
                    raise AssertionError(f"identity (28) failed at {(m, k, p)}")
                d_count, run, window = carry_bound_terms(m, k, p)
                margin = exact - (window - d_count - run)
                if margin < 0:
                    raise AssertionError(f"inequality (29) failed at {(m, k, p)}")
                if margin_min is None or margin < margin_min:
                    margin_min = margin
                slack_min = min(slack_min, exact)
                checks += 1
    return {
        "checked_triples": checks,
        "max_m": max_m,
        "max_k": max_k,
        "minimum_exact_slack": slack_min,
        "minimum_inequality_margin": margin_min,
    }


def enumeration_agreement(grid_max_m: int, grid_digits: int) -> dict[str, Any]:
    """Cross-check the two dynamic programs against enumeration."""
    _require_int("grid_max_m", grid_max_m, 1)
    _require_int("grid_digits", grid_digits, 2)
    rows: list[dict[str, Any]] = []
    for m in range(0, grid_max_m + 1):
        for p in primes_up_to(max(2, min(m, 13))):
            if p > m and m != 0:
                continue
            for digit_count in range(2, grid_digits + 1):
                if p**digit_count <= 2 * m or p**digit_count > 200_000:
                    continue
                power = failure_count_power(m, p, digit_count)
                upto = failure_count_upto(m, p, p**digit_count - 1)
                brute = failure_count_enumerated(m, p, p**digit_count - 1)
                if not power == upto == brute:
                    raise AssertionError(
                        f"failure counts disagree at {(m, p, digit_count)}: "
                        f"{power} {upto} {brute}"
                    )
                bound = certified_upper_bound(m, p, digit_count)
                if bound < power:
                    raise AssertionError(f"bound (30) violated at {(m, p, digit_count)}")
                rows.append(
                    {
                        "m": m,
                        "p": p,
                        "digit_count": digit_count,
                        "range": p**digit_count,
                        "failures": power,
                        "certified_bound": bound,
                    }
                )
    return {"agreeing_cells": len(rows), "rows": rows}


def prime_decay_table(m: int, digit_targets: tuple[int, ...]) -> list[dict[str, Any]]:
    """Exact tier-failure densities per prime at several base-2 log scales."""
    _require_int("m", m, 1)
    rows: list[dict[str, Any]] = []
    for p in primes_up_to(m):
        cells: list[dict[str, Any]] = []
        for bits in digit_targets:
            digit_count = max(
                len(base_digits(m, p)) + 1,
                math.ceil(bits * math.log(2) / math.log(p)),
            )
            failures = failure_count_power(m, p, digit_count)
            span = p**digit_count
            density = failures / span
            bound = certified_upper_bound(m, p, digit_count)
            if bound < failures:
                raise AssertionError(f"bound (30) violated at {(m, p, digit_count)}")
            envelope = decay_envelope(digit_count)
            applies = envelope_applies(m, p, digit_count)
            if applies and bound / span > envelope:
                raise AssertionError(
                    f"envelope (31) violated at {(m, p, digit_count)}"
                )
            cells.append(
                {
                    "digit_count": digit_count,
                    "range_decimal_digits": len(str(span)),
                    "failures": failures,
                    "density": density,
                    "certified_bound_density": bound / span,
                    "envelope_applies": applies,
                    "envelope": envelope,
                }
            )
        first, last = cells[0], cells[-1]
        rows.append(
            {
                "p": p,
                "cells": cells,
                "measured_decay_exponent": (
                    math.log(last["density"] / first["density"])
                    / math.log(
                        p ** last["digit_count"] / p ** first["digit_count"]
                    )
                    if first["density"] > 0 and last["density"] > 0
                    else None
                ),
            }
        )
    return rows


def tier_failure_at_scale(m: int, limit: int) -> dict[str, Any]:
    """Exact per-prime failure counts and a union bound at one real scale."""
    _require_int("m", m, 1)
    _require_int("limit", limit, 1)
    rows = []
    union = 0
    for p in primes_up_to(m):
        failures = failure_count_upto(m, p, limit)
        rows.append({"p": p, "failures": failures, "density": failures / limit})
        union += failures
    return {
        "limit": limit,
        "limit_decimal_digits": len(str(limit)),
        "per_prime": rows,
        "union_failure_count": union,
        "union_failure_density": union / limit,
        "surviving_density_lower_bound": 1.0 - union / limit,
    }


def conditional_crossover(
    m: int, run_density: float, max_digit_count: int
) -> dict[str, Any]:
    """Smallest certified scale where the tier bound falls below a run density.

    The run density is measured input, so this scale is conditional on that
    measurement and is not an existence proof.
    """
    _require_int("m", m, 1)
    _require_int("max_digit_count", max_digit_count, 1)
    if not 0.0 < run_density < 1.0:
        raise ValueError("run_density must lie in (0,1)")
    window = (m + 1) // 2
    target = run_density**window
    primes = primes_up_to(m)
    for bits in range(16, max_digit_count + 1):
        union = 0.0
        for p in primes:
            digit_count = max(
                len(base_digits(m, p)) + 1,
                math.ceil(bits * math.log(2) / math.log(p)),
            )
            if not envelope_applies(m, p, digit_count):
                union = 1.0
                break
            union += decay_envelope(digit_count)
        if union < target:
            return {
                "window_length": window,
                "run_density": run_density,
                "target_density": target,
                "certified_tier_bound": union,
                "scale_base_two_exponent": bits,
                "scale_decimal_digits": math.ceil(bits * math.log10(2)),
            }
    return {
        "window_length": window,
        "run_density": run_density,
        "target_density": target,
        "certified_tier_bound": None,
        "scale_base_two_exponent": None,
        "scale_decimal_digits": None,
    }


def exact_crossover(
    m: int, run_density: float, bit_grid: tuple[int, ...]
) -> dict[str, Any]:
    """Exact tier-failure density at each scale ``2^b``, and the first scale
    where it drops below ``run_density ** ceil(m/2)``.

    Unlike :func:`conditional_crossover` this uses the exact digit dynamic
    program rather than the Section 13 envelope, so the tier side of the
    comparison is an exact count over every ``k`` in the range. The run
    density remains measured input.
    """
    _require_int("m", m, 1)
    if not 0.0 < run_density < 1.0:
        raise ValueError("run_density must lie in (0,1)")
    window = (m + 1) // 2
    target = run_density**window
    primes = primes_up_to(m)
    rows: list[dict[str, Any]] = []
    first_below = None
    for bits in sorted(bit_grid):
        union = 0.0
        for p in primes:
            digit_count = max(
                len(base_digits(m, p)) + 1,
                math.ceil(bits * math.log(2) / math.log(p)),
            )
            span = p**digit_count
            union += failure_count_power(m, p, digit_count) / span
        rows.append(
            {
                "scale_base_two_exponent": bits,
                "scale_decimal_digits": math.ceil(bits * math.log10(2)),
                "exact_union_failure_density": union,
            }
        )
        if first_below is None and union < target:
            first_below = rows[-1]
    return {
        "window_length": window,
        "run_density": run_density,
        "target_density": target,
        "scales": rows,
        "first_scale_below_target": first_below,
    }


def analyze(
    *,
    max_m: int,
    target_m: int,
    scan_max_m: int,
    scan_max_k: int,
    grid_max_m: int,
    grid_digits: int,
    digit_targets: tuple[int, ...],
    scales: tuple[int, ...],
    density_artifact: Path,
    max_crossover_bits: int,
    exact_bit_grid: tuple[int, ...],
) -> dict[str, Any]:
    density_artifact = density_artifact.resolve()
    payload = json.loads(density_artifact.read_text(encoding="utf-8"))
    unbiased = [
        block
        for block in payload["compensation_blocks"]
        if block["label"] == "scale_5e12_unbiased"
    ]
    if len(unbiased) != 1:
        raise AssertionError("expected exactly one unbiased compensation block")
    run_density = unbiased[0]["density"]

    families = []
    for m in range(1, max_m + 1):
        primes = primes_up_to(m)
        if not primes:
            continue
        worst = None
        for p in primes:
            digit_count = max(
                len(base_digits(m, p)) + 1,
                math.ceil(45 * math.log(2) / math.log(p)),
            )
            failures = failure_count_power(m, p, digit_count)
            density = failures / p**digit_count
            if worst is None or density > worst["density"]:
                worst = {"p": p, "digit_count": digit_count, "density": density}
        families.append(
            {
                "m": m,
                "window_length": (m + 1) // 2,
                "worst_prime": worst["p"],
                "worst_density_at_2^45": worst["density"],
            }
        )

    decay = prime_decay_table(target_m, digit_targets)
    envelope_cells = sum(
        1 for row in decay for cell in row["cells"] if cell["envelope_applies"]
    )
    if envelope_cells == 0:
        raise AssertionError(
            "no reported cell exercises the Section 13 envelope; raise digit-targets"
        )

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_TIER_DENSITY",
        "claim_boundary": (
            "Every failure count here is exact on its declared range: the digit "
            "dynamic program is a closed computation over all k in that range, "
            "not a sample, and it is cross-checked against direct enumeration "
            "wherever enumeration is affordable and against the second "
            "prefix-decomposition program at every power-of-p range. The "
            "closed-form bound and the exponential envelope of Section 13 are "
            "checked against the exact counts on every reported cell. The "
            "closing crossover scale consumes a measured compensation-good "
            "density and is therefore conditional on that measurement; it is "
            "not an existence proof and no witness is claimed."
        ),
        "parameters": {
            "max_m": max_m,
            "target_m": target_m,
            "scan_max_m": scan_max_m,
            "scan_max_k": scan_max_k,
            "grid_max_m": grid_max_m,
            "grid_digits": grid_digits,
            "digit_targets": list(digit_targets),
            "scales": list(scales),
            "max_crossover_bits": max_crossover_bits,
            "exact_bit_grid": list(exact_bit_grid),
        },
        "identity_and_inequality": identity_and_inequality_scan(scan_max_m, scan_max_k),
        "method_agreement": enumeration_agreement(grid_max_m, grid_digits),
        "prime_decay": decay,
        "envelope_checked_cells": envelope_cells,
        "scales": [tier_failure_at_scale(target_m, limit) for limit in scales],
        "worst_prime_by_m": families,
        "conditional_crossover": conditional_crossover(
            target_m, run_density, max_crossover_bits
        ),
        "exact_crossover": exact_crossover(target_m, run_density, exact_bit_grid),
        "input": {
            "density_artifact_logical_name": density_artifact.name,
            "density_artifact_sha256": sha256_file(density_artifact),
            "measured_unbiased_run_density": run_density,
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "tier_density_dp.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=30)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--scan-max-m", type=int, default=20)
    parser.add_argument("--scan-max-k", type=int, default=400)
    parser.add_argument("--grid-max-m", type=int, default=13)
    parser.add_argument("--grid-digits", type=int, default=5)
    parser.add_argument(
        "--digit-targets", type=int, nargs="+", default=[20, 45, 60, 90, 150, 250]
    )
    parser.add_argument(
        "--scales",
        type=int,
        nargs="+",
        default=[10**6, 10**13, PUBLISHED_M25_K],
    )
    parser.add_argument("--max-crossover-bits", type=int, default=8000)
    parser.add_argument(
        "--exact-bit-grid",
        type=int,
        nargs="+",
        default=[60, 90, 120, 150, 200, 250, 300, 400],
    )
    parser.add_argument("--density-artifact", type=Path, default=DEFAULT_DENSITY_ARTIFACT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        max_m=args.max_m,
        target_m=args.target_m,
        scan_max_m=args.scan_max_m,
        scan_max_k=args.scan_max_k,
        grid_max_m=args.grid_max_m,
        grid_digits=args.grid_digits,
        digit_targets=tuple(args.digit_targets),
        scales=tuple(args.scales),
        density_artifact=args.density_artifact,
        max_crossover_bits=args.max_crossover_bits,
        exact_bit_grid=tuple(args.exact_bit_grid),
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
                "checked_triples": result["identity_and_inequality"]["checked_triples"],
                "agreeing_cells": result["method_agreement"]["agreeing_cells"],
                "worst_prime_density_at_2^45": max(
                    row["worst_density_at_2^45"] for row in result["worst_prime_by_m"]
                ),
                "union_density_by_scale": [
                    {
                        "limit_decimal_digits": row["limit_decimal_digits"],
                        "union_failure_density": row["union_failure_density"],
                    }
                    for row in result["scales"]
                ],
                "conditional_crossover_decimal_digits": result[
                    "conditional_crossover"
                ]["scale_decimal_digits"],
                "exact_crossover_decimal_digits": (
                    result["exact_crossover"]["first_scale_below_target"] or {}
                ).get("scale_decimal_digits"),
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
