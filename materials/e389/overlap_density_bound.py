#!/usr/bin/env python3
"""Certify that the compensation-good integers have positive lower density.

Section 17 proves Conjecture R(1). The first-order union bound over failure
events tends to ``log 2 + 0.322521 = 1.015668 > 1`` and therefore proves
nothing. The second Bonferroni term repairs it: the two failure families

  A   some prime p > m with p | w and p^2 > 2w      (short cofactor)
  B   some prime q > m with q^e || w, q^2 <= 2w and C_q(w/q^e) < e   (level)

satisfy the exact identity ``#good = x - #A - #B + #(A and B)`` on any block,
so a *lower* bound on the overlap buys back exactly what the union bound
overshoots. Only ``0.015668`` of overlap is needed.

This producer certifies three numbers, each as a one-sided bound with an
explicit tail, never as a measurement:

  level_union_bound       upper bound for B, sum_r 2^-r log((r+2)/(r+1))
  exponent_tail           upper bound for the e >= 2 part of B, which -> 0
  overlap_lower_bound     lower bound for #(A and B), from the disjoint
                          families D_1 and D_2 of Section 17

The overlap families are counted on the *prime* side: the free variable is the
large prime p, and the level condition is a congruence on p modulo q or q^2.
Vinogradov's exponential-sum bound applies to each individual modulus in the
admissible range, so no averaging over moduli and no Bombieri-Vinogradov input
is required. Bounding the same overlap from the smooth side instead would need
smooth numbers in progressions to moduli near sqrt(x), which is why the earlier
attempt in Section 13 stalled.

Bands 1 and 2 are used because there the family member is unique: a band-r
prime q > (2w)^(1/(r+2)) together with p > (2w)^(1/2) forces
(2w)^(1/2 + 2/(r+2)) < w, impossible for r <= 2, so no integer is counted twice
and the two families are disjoint.
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
from itertools import combinations
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import SPFFactorizer, primes_up_to
from good_density_decomposition import classify_term, compensation_count

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "overlap_density_bound.json"

#: Bands whose overlap family is automatically multiplicity-free.
CERTIFIED_BANDS = (1, 2)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _require_unit(name: str, value: float) -> None:
    if not isinstance(value, float) or not 0.0 <= value < 0.5:
        raise ValueError(f"{name} must be a float in [0, 0.5)")


def level_union_bound(terms: int = 200) -> dict[str, Any]:
    """Upper bound for the density of level failures with exponent one.

    A prime with exactly ``r`` usable levels satisfies
    ``(2w)^(1/(r+2)) < q <= (2w)^(1/(r+1))``; Mertens gives ``sum 1/q =
    log((r+2)/(r+1))`` over that range and the ``r`` level tests cost ``2^-r``.
    The discarded tail is bounded by ``log((r+2)/(r+1)) <= 1/(r+1)``, so
    ``sum_{r>terms} 2^-r/(r+1) <= 2^-terms/(terms+2)``.
    """
    _require_int("terms", terms, 1)
    partial = sum(2.0 ** (-r) * math.log((r + 2) / (r + 1)) for r in range(1, terms + 1))
    tail = 2.0 ** (-terms) / (terms + 2)
    return {"terms": terms, "partial": partial, "tail_bound": tail, "bound": partial + tail}


def exponent_tail_bound(
    scale: float, target_m: int, prime_cutoff: int = 10**6
) -> dict[str, Any]:
    """Upper bound for level failures carrying exponent ``e >= 2`` at scale ``w ~ scale``.

    For ``q^e || w`` the density is at most ``q^-e`` and the failure needs fewer
    than ``e`` of the ``R = floor(log_q 2w) - e`` levels to land in the upper
    half, so it costs at most ``sum_{i<e} C(R,i) 2^-R``. Primes above
    ``prime_cutoff`` are discarded by ``sum_{n>Q} n^-2 (1-1/n)^-1 < 2/Q``.

    The bound tends to zero with the scale: every fixed ``q`` gains levels, and
    the large-``q`` end is killed by ``q^-e``.
    """
    if not isinstance(scale, float) or scale <= 4.0:
        raise ValueError("scale must be a float > 4")
    _require_int("target_m", target_m, 0)
    _require_int("prime_cutoff", prime_cutoff, 100)
    log_scale = math.log(2.0 * scale)
    total = 0.0
    for q in primes_up_to(prime_cutoff):
        if q <= target_m:
            continue
        levels = int(log_scale / math.log(q))
        for exponent in range(2, levels + 1):
            remaining = max(0, levels - exponent)
            cost = sum(
                math.comb(remaining, i) for i in range(min(exponent, remaining + 1))
            ) * 2.0 ** (-remaining)
            total += q ** (-exponent) * min(1.0, cost)
    tail = 2.0 / prime_cutoff
    return {
        "scale": scale,
        "prime_cutoff": prime_cutoff,
        "enumerated": total,
        "tail_bound": tail,
        "bound": total + tail,
    }


def band_integrand(beta: float, band: int, eta: float) -> float:
    """Inner ``alpha`` measure of the band-``band`` overlap family at ``beta``.

    ``w = p q t`` contributes when ``alpha = log p / log 2w`` satisfies
    ``alpha >= 1/2 + eta`` (so ``p^2 > 2w``), ``alpha >= band*beta*(1+eta)`` (so
    Vinogradov's modulus ``q^band`` stays below the length of the ``p`` range)
    and ``alpha <= 1 - beta`` (so ``t >= 1``). The measure is ``d alpha / alpha``.
    """
    floor = max(0.5 + eta, band * beta * (1.0 + eta))
    ceiling = 1.0 - beta
    if ceiling <= floor:
        return 0.0
    return math.log(ceiling / floor) / beta


def band_overlap_lower_bound(band: int, eta: float, steps: int = 20_000) -> dict[str, Any]:
    """Certified lower bound for the density of the band-``band`` overlap family.

    The integrand is a positive decreasing function of ``beta`` on the band
    (both ``1/beta`` and ``log((1-beta)/max(...))`` decrease), so the
    right-endpoint Riemann sum is a rigorous lower bound for the integral.
    Monotonicity is asserted, not assumed.
    """
    _require_int("band", band, 1)
    _require_unit("eta", eta)
    _require_int("steps", steps, 16)
    lo, hi = 1.0 / (band + 2), 1.0 / (band + 1)
    width = (hi - lo) / steps
    total = 0.0
    previous = None
    for i in range(1, steps + 1):
        beta = lo + i * width
        value = band_integrand(beta, band, eta)
        if previous is not None and value > previous + 1e-12:
            raise AssertionError(
                f"integrand not decreasing on band {band} at beta={beta}"
            )
        previous = value
        total += value * width
    return {
        "band": band,
        "eta": eta,
        "steps": steps,
        "beta_range": [lo, hi],
        "level_cost": 2.0 ** (-band),
        "bound": total * 2.0 ** (-band),
    }


def overlap_lower_bound(eta: float, steps: int = 20_000) -> dict[str, Any]:
    """Sum the certified band bounds over the multiplicity-free bands."""
    bands = [band_overlap_lower_bound(band, eta, steps) for band in CERTIFIED_BANDS]
    return {"bands": bands, "bound": sum(row["bound"] for row in bands)}


def band_of(value: int, prime: int) -> int | None:
    """Return ``r`` with ``prime^(r+1) <= 2*value < prime^(r+2)``, else ``None``."""
    _require_int("value", value, 2)
    _require_int("prime", prime, 2)
    if prime * prime > 2 * value:
        return None
    band, power = 0, prime
    while power * prime <= 2 * value:
        band += 1
        power *= prime
    return band


def overlap_family_member(value: int, target_m: int, factors: dict[int, int]) -> int | None:
    """Return the band of ``value`` in the Section 17 overlap family, else ``None``.

    Membership requires a short-cofactor prime and a band-1 or band-2 prime
    whose every level lands in the lower half.
    """
    if not any(
        prime > target_m and exponent == 1 and prime * prime > 2 * value
        for prime, exponent in factors.items()
    ):
        return None
    found = None
    for prime, exponent in factors.items():
        if prime <= target_m or exponent != 1:
            continue
        band = band_of(value, prime)
        if band in CERTIFIED_BANDS and compensation_count(value // prime, prime) == 0:
            if found is not None:
                raise AssertionError(f"overlap family member {value} counted twice")
            found = band
    return found


def census_block(start: int, count: int, target_m: int) -> dict[str, Any]:
    """Exhaustively verify the Section 17 inclusion on ``[start, start+count)``."""
    _require_int("start", start, 8)
    _require_int("count", count, 1)
    _require_int("target_m", target_m, 0)
    factorizer = SPFFactorizer(start + count + 1)
    counts = {"good": 0, "short": 0, "level": 0, "both": 0}
    family = {band: 0 for band in CERTIFIED_BANDS}
    for value in range(start, start + count):
        short, level = classify_term(value, target_m, factorizer)
        if short:
            counts["short"] += 1
        if level:
            counts["level"] += 1
        if short and level:
            counts["both"] += 1
        if not short and not level:
            counts["good"] += 1
        band = overlap_family_member(value, target_m, factorizer.factor(value))
        if band is not None:
            if not (short and level):
                raise AssertionError(
                    f"overlap family member {value} is outside A and B"
                )
            family[band] += 1
    identity = count - counts["short"] - counts["level"] + counts["both"]
    if identity != counts["good"]:
        raise AssertionError("inclusion-exclusion identity failed on the block")
    return {
        "start": start,
        "count": count,
        "counts": counts,
        "family_counts": {str(band): value for band, value in family.items()},
        "densities": {name: value / count for name, value in counts.items()},
        "family_density": sum(family.values()) / count,
        "identity_holds": True,
    }


def pair_accounting(start: int, count: int, target_m: int) -> dict[str, Any]:
    """Exact four-event inclusion-exclusion for ``w`` and ``w+1`` both good.

    The four events are ``A``, ``B`` at ``w`` and at ``w+1``. Truncating
    inclusion-exclusion after three sums is an upper bound for the union
    (Bonferroni), so ``S1 - S2 + S3 < 1`` would give R(2) the way ``(47)`` gives
    R(1). Under independence across the two terms that total collapses to

        S1 - S2 + S3 = 1 + (u - 1) * (2*pi - (u - 1)),   u = log2 + Lambda_B,

    so it drops below one exactly when ``pi < (u-1)/2``. The overlap that proves
    R(1) is ``pi >= 0.045``, six times too large: the same term that rescues one
    coordinate defeats two, because it enters with the opposite sign.
    """
    _require_int("start", start, 8)
    _require_int("count", count, 1)
    _require_int("target_m", target_m, 0)
    factorizer = SPFFactorizer(start + count + 2)
    flags: list[tuple[bool, bool, bool, bool]] = []
    for value in range(start, start + count):
        short0, level0 = classify_term(value, target_m, factorizer)
        short1, level1 = classify_term(value + 1, target_m, factorizer)
        flags.append((short0, level0, short1, level1))
    sums = {}
    for depth in range(1, 5):
        sums[depth] = sum(
            sum(1 for row in flags if all(row[i] for i in idx))
            for idx in combinations(range(4), depth)
        ) / count
    union = sum(1 for row in flags if any(row)) / count
    depth3 = sums[1] - sums[2] + sums[3]
    if depth3 < union:
        raise AssertionError("Bonferroni depth 3 must dominate the union")
    if abs(sums[1] - sums[2] + sums[3] - sums[4] - union) > 1e-12:
        raise AssertionError("four-event inclusion-exclusion is not exact")
    return {
        "start": start,
        "count": count,
        "S1": sums[1],
        "S2": sums[2],
        "S3": sums[3],
        "S4": sums[4],
        "exact_union": union,
        "both_good": 1.0 - union,
        "bonferroni_depth3": depth3,
        "single_term_total": sums[1] / 2.0,
    }


def pair_obstruction(overlap: float, level_bound: float) -> dict[str, Any]:
    """Where the depth-3 route stands asymptotically, in closed form."""
    u = math.log(2) + level_bound
    return {
        "u": u,
        "overlap": overlap,
        "depth3_limit": 1.0 + (u - 1.0) * (2.0 * overlap - (u - 1.0)),
        "overlap_threshold_for_closure": (u - 1.0) / 2.0,
        "shortfall_factor": overlap / ((u - 1.0) / 2.0),
        "closes": 1.0 + (u - 1.0) * (2.0 * overlap - (u - 1.0)) < 1.0,
    }


def finitary_threshold(target_m: int, run_density: float) -> dict[str, Any]:
    """Where the finitary form of the Section 13 reduction starts to bite.

    (32) bounds tier failures by ``4*pi(m)*m*X^(1-c)``, ``c = 1/(24 log m)``, so
    a single scale with more good runs than that yields a witness. The criterion
    is vacuous until that bound falls below ``X``, and useful only once it falls
    below the run density.
    """
    _require_int("target_m", target_m, 2)
    if not isinstance(run_density, float) or not 0.0 < run_density < 1.0:
        raise ValueError("run_density must be a float in (0, 1)")
    decay = 1.0 / (24.0 * math.log(target_m))
    constant = 4 * len(primes_up_to(target_m)) * target_m
    return {
        "target_m": target_m,
        "decay_exponent": decay,
        "constant": constant,
        "run_density": run_density,
        "log10_scale_nonvacuous": math.log10(constant) / decay,
        "log10_scale_beats_run_density": (
            math.log10(constant) - math.log10(run_density)
        ) / decay,
    }


def analyze(
    *,
    blocks: tuple[tuple[int, int], ...],
    target_m: int,
    eta: float,
    scales: tuple[float, ...],
    steps: int,
) -> dict[str, Any]:
    level = level_union_bound()
    overlap = overlap_lower_bound(eta, steps)
    tails = [exponent_tail_bound(scale, target_m) for scale in scales]
    if any(a["bound"] <= b["bound"] for a, b in zip(tails, tails[1:])):
        raise AssertionError("the exponent tail bound must decrease with scale")
    limiting_tail = tails[-1]["bound"]
    first_order = math.log(2) + level["bound"]
    margin = 1.0 - first_order - limiting_tail + overlap["bound"]
    if first_order <= 1.0:
        raise AssertionError("first-order union bound no longer exceeds one")
    if margin <= 0.0:
        raise AssertionError("certified margin is not positive")
    rows = [census_block(start, count, target_m) for start, count in blocks]
    pairs = [pair_accounting(start, count, target_m) for start, count in blocks]
    obstruction = pair_obstruction(overlap["bound"], level["bound"])
    finitary = finitary_threshold(27, 0.123664**14)
    if obstruction["closes"]:
        raise AssertionError(
            "the depth-3 route now closes R(2); the recorded obstruction is stale"
        )
    for row in rows:
        if row["family_density"] > row["densities"]["both"]:
            raise AssertionError("overlap family exceeds the measured overlap")
    return {
        "schema_version": 1,
        "status": "PASS_CERTIFIED_POSITIVE_GOOD_DENSITY",
        "claim_boundary": (
            "The three assembled numbers are one-sided bounds with explicit "
            "tails, not measurements: level_union_bound and exponent_tail are "
            "upper bounds for the two level-failure families, and "
            "overlap_lower_bound is a right-endpoint Riemann sum, hence a lower "
            "bound, for the asymptotic density of the Section 17 overlap "
            "families. The asymptotic evaluation of those families uses "
            "Mertens for the prime sums and Vinogradov's exponential-sum bound "
            "for the level congruence; both are cited theorems, not verified "
            "here. The block census is exhaustive and verifies the "
            "inclusion-exclusion identity and the family inclusion exactly, "
            "but it is finite evidence and proves nothing asymptotic."
        ),
        "parameters": {
            "target_m": target_m,
            "eta": eta,
            "steps": steps,
            "scales": list(scales),
            "blocks": [list(block) for block in blocks],
            "certified_bands": list(CERTIFIED_BANDS),
        },
        "level_union_bound": level,
        "exponent_tails": tails,
        "overlap_lower_bound": overlap,
        "certificate": {
            "log_two": math.log(2),
            "first_order_union_total": first_order,
            "first_order_deficit": first_order - 1.0,
            "limiting_exponent_tail": limiting_tail,
            "overlap_bound": overlap["bound"],
            "certified_lower_density": margin,
        },
        "blocks": rows,
        "pair_accounting": pairs,
        "pair_obstruction": obstruction,
        "finitary_reduction": finitary,
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "overlap_density_bound.py": sha256_file(Path(__file__)),
            "good_density_decomposition.py": sha256_file(
                ROOT / "good_density_decomposition.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=1)
    parser.add_argument("--eta", type=float, default=0.005)
    parser.add_argument("--steps", type=int, default=20_000)
    parser.add_argument("--scales", type=float, nargs="+", default=[1e8, 1e20, 1e40])
    parser.add_argument("--block-starts", type=int, nargs="+", default=[10**6, 10**7])
    parser.add_argument("--block-count", type=int, default=200_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        blocks=tuple((start, args.block_count) for start in args.block_starts),
        target_m=args.target_m,
        eta=args.eta,
        scales=tuple(args.scales),
        steps=args.steps,
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
                "certified_lower_density": result["certificate"][
                    "certified_lower_density"
                ],
                "first_order_deficit": result["certificate"]["first_order_deficit"],
                "overlap_bound": result["certificate"]["overlap_bound"],
                "block_family_densities": [
                    row["family_density"] for row in result["blocks"]
                ],
                "r2_depth3_limit": result["pair_obstruction"]["depth3_limit"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
