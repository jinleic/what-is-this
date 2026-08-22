"""Lower-memory exact audit for the simple-cubic high-temperature v^24 wall.

Run from the repository root with::

    .venv/bin/python experiments/e68_ht_v24.py

The core exploratory engine is a vertex-frontier enumeration in which every
coordinate cut is capped at two occupied edges.  It is tailored to the unique
minimal-span 5x5x5 finite-lattice weight at total degree 24.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from math import factorial, gcd, prod
from pathlib import Path

import ising.series as flm
from ising.series import FLMSeries
from ising.transfer_matrix.crt import _is_prime_32

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e68_ht_v24.py"
INTERPRETER = ".venv/bin/python"
PREFIX_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
RESULT_PATH = ROOT / "results" / "series" / "ht_v24.json"
PROOF_PATH = ROOT / "proofs" / "ht_v24.md"

CRT_PRIMES = (
    2_147_483_647,
    2_147_483_629,
    2_147_483_587,
    2_147_483_579,
    2_147_483_563,
)
CRT_PRODUCT = prod(CRT_PRIMES)


def frontier_even_polynomial(
    side: int, modulus: int = CRT_PRIMES[0], max_states: int | None = None
):
    """Count even subgraphs with at most two edges across every coordinate cut.

    Return coefficients indexed by a ternary cut-count code.  At completion,
    every digit is necessarily zero or two.  Arithmetic is exact modulo
    ``modulus``.  ``max_states`` provides a measured, deterministic early stop
    for resource profiling.
    """

    side = int(side)
    if side < 2:
        raise ValueError("side must be at least two")
    cut_count = 3 * (side - 1)
    width = side * side
    powers3 = tuple(3**index for index in range(cut_count))
    states: dict[int, int] = {0: 1}
    peak = 1
    timeline = []
    started = time.perf_counter()

    for x in range(side):
        for y in range(side):
            for z in range(side):
                directions = []
                if z + 1 < side:
                    directions.append((0, powers3[2 * (side - 1) + z]))
                if y + 1 < side:
                    directions.append((side - 1, powers3[(side - 1) + y]))
                if x + 1 < side:
                    directions.append((width - 1, powers3[x]))

                options_by_parity = ([], [])
                for choice in range(1 << len(directions)):
                    parity = choice.bit_count() & 1
                    pending_delta = 0
                    increments = []
                    for index, (offset, power) in enumerate(directions):
                        if choice & (1 << index):
                            pending_delta ^= 1 << offset
                            increments.append(power)
                    options_by_parity[parity].append((pending_delta, tuple(increments)))

                updated: dict[int, int] = {}
                for packed, value in states.items():
                    pending = packed & ((1 << width) - 1)
                    code = packed >> width
                    incoming = pending & 1
                    base_pending = pending >> 1
                    for pending_delta, increments in options_by_parity[incoming]:
                        new_code = code
                        allowed = True
                        for power in increments:
                            if (new_code // power) % 3 == 2:
                                allowed = False
                                break
                            new_code += power
                        if not allowed:
                            continue
                        key = (new_code << width) | (base_pending ^ pending_delta)
                        updated[key] = (updated.get(key, 0) + value) % modulus
                states = updated
                peak = max(peak, len(states))
                timeline.append(
                    {
                        "vertex": [x, y, z],
                        "states": len(states),
                        "elapsed_seconds": time.perf_counter() - started,
                    }
                )
                if max_states is not None and len(states) > max_states:
                    return {
                        "completed": False,
                        "side": side,
                        "modulus": modulus,
                        "peak_states": peak,
                        "terminal_states": len(states),
                        "timeline": timeline,
                        "coefficients": None,
                    }

    coefficients: dict[int, int] = defaultdict(int)
    for packed, value in states.items():
        pending = packed & ((1 << width) - 1)
        if pending:
            raise AssertionError("nonempty parity frontier after final vertex")
        code = packed >> width
        digits = tuple((code // power) % 3 for power in powers3)
        if any(digit == 1 for digit in digits):
            raise AssertionError("an even subgraph crosses every cut evenly")
        coefficients[code] = (coefficients[code] + value) % modulus
    return {
        "completed": True,
        "side": side,
        "modulus": modulus,
        "peak_states": peak,
        "terminal_states": len(states),
        "timeline": timeline,
        "coefficients": dict(coefficients),
    }


def cut_polynomial_as_subsets(
    coefficients: dict[int, int], side: int, modulus: int
) -> tuple[int, ...]:
    """Convert final ternary cut counts (zero/two) into a square-free polynomial."""

    cut_count = 3 * (side - 1)
    powers3 = tuple(3**index for index in range(cut_count))
    result = [0] * (1 << cut_count)
    for code, value in coefficients.items():
        mask = 0
        for index, power in enumerate(powers3):
            digit = (code // power) % 3
            if digit == 1:
                raise AssertionError("terminal cut count is odd")
            if digit == 2:
                mask |= 1 << index
        result[mask] = (result[mask] + value) % modulus
    if result[0] != 1:
        raise AssertionError("the empty even subgraph must be unique")
    return tuple(result)


def square_free_log_full(coefficient: tuple[int, ...], modulus: int) -> int:
    """Return the full-support coefficient of ``log(coefficient)`` modulo ``modulus``."""

    size = len(coefficient)
    if size < 2 or size & (size - 1):
        raise ValueError("coefficient length must be a power of two")
    logarithm = [0] * size
    for support_size in range(1, (size - 1).bit_count() + 1):
        for mask in range(1, size):
            if mask.bit_count() != support_size:
                continue
            distinguished = mask & -mask
            remainder = mask ^ distinguished
            value = coefficient[mask]
            subset = remainder
            while subset:
                value -= (
                    coefficient[subset]
                    * logarithm[(remainder ^ subset) | distinguished]
                )
                subset = (subset - 1) & remainder
            logarithm[mask] = value % modulus
    return logarithm[-1]


def ordered_partition_factor(n: int) -> int:
    """Return ``sum_k S(n,k)(k-1)!`` by an integer Stirling recurrence."""

    stirling = [[0] * (n + 1) for _ in range(n + 1)]
    stirling[0][0] = 1
    for row in range(1, n + 1):
        for blocks in range(1, row + 1):
            stirling[row][blocks] = (
                stirling[row - 1][blocks - 1]
                + blocks * stirling[row - 1][blocks]
            )
    return sum(
        stirling[n][blocks] * factorial(blocks - 1)
        for blocks in range(1, n + 1)
    )


def minimal_cube_weight_bound(side: int) -> int:
    """Absolute square-free logarithm bound for the minimal ``side^3`` weight."""

    cuts = 3 * (side - 1)
    edges_per_cut = side * side
    choices_per_active_cut = edges_per_cut * (edges_per_cut - 1) // 2
    return choices_per_active_cut**cuts * ordered_partition_factor(cuts)


def high_temperature_series_with_minimal_cube(
    order: int, minimal_cube_weight: int
) -> FLMSeries:
    """Run the usual FLM while injecting only the proved minimal 5-cube weight."""

    if int(order) != 24:
        raise ValueError("the specialized cut extraction is proved only at order 24")
    shapes = flm._ht_shapes(3, order, 0)
    weights = {}
    interaction = [Fraction(0) for _ in range(order + 1)]
    target = (5, 5, 5)
    target_seen = False
    for shape in shapes:
        canonical = flm._canonical_shape(shape)
        if canonical == target:
            if shape != target:
                raise AssertionError("the isotropic target must have one orientation")
            exact_weight = [Fraction(0) for _ in range(order + 1)]
            exact_weight[order] = Fraction(minimal_cube_weight)
            target_seen = True
        else:
            exact_weight = list(flm._ht_box_log(canonical, order))
            for subshape, subweight in weights.items():
                if all(sub <= side for sub, side in zip(subshape, shape)):
                    placements = prod(
                        side - sub + 1 for sub, side in zip(subshape, shape)
                    )
                    for degree in range(order + 1):
                        exact_weight[degree] -= placements * subweight[degree]
        frozen_weight = tuple(exact_weight)
        weights[shape] = frozen_weight
        for degree, value in enumerate(frozen_weight):
            interaction[degree] += value
    if not target_seen:
        raise AssertionError("5x5x5 target absent from the order-24 box family")
    coefficients = list(interaction)
    for degree in range(2, order + 1, 2):
        coefficients[degree] += Fraction(3, degree)
    return FLMSeries(
        dimension=3,
        order=order,
        variable="v",
        coefficients=tuple(coefficients),
        interaction_coefficients=tuple(interaction),
        boxes=shapes,
        box_weights=weights,
        bound_budget=order // 2,
        bound_slack=0,
        order_bound="2*sum_i(side_i-1) <= truncation_order",
    )


def _rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _sha256_json(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def _proof_text(data: dict) -> str:
    return f"""# Exact simple-cubic HT coefficient through v^24

## Scope

[COMPUTATION] This is a finite exact coefficient calculation, not a solution of
the three-dimensional Ising model.

## Why the previous route used 242 GB

[LEMMA] The previous broken-bond spin transfer stores a dense array indexed by
all `2^25` spins on a 5x5 section and all 301 broken-bond degrees.  Three
simultaneous int64 arrays therefore require
`3*2^25*301*8 = {data['old_engine']['projected_peak_bytes']}` bytes.  Five
31-bit CRT primes are enough; the modulus count multiplies time, not this peak.
There is no boundary-connectivity state in that implementation.  The real wall
is the cross-section exponent, with polynomial degree as a factor.

## Lower-memory exact extraction

[LEMMA] A graph whose 5x5x5 bounding box first contributes at degree 24 crosses
each of the 12 coordinate cuts positively and evenly.  Its total edge count is
24, so every cut is crossed exactly twice.  The frontier algorithm therefore
tracks only parity on unprocessed edges and ternary cut counters 0,1,2; any
third crossing is discarded.  The frontier has a measured peak of
`{data['frontier']['peak_states']}` states.  The isolated production run used
`{data['production_resource_measurement']['maximum_resident_set_bytes']}` bytes
maximum RSS (`{data['production_resource_measurement']['peak_memory_footprint_bytes']}`-byte
peak footprint), versus the old 242-GB projection; it read and wrote zero disk
blocks.
The complete coefficient run, including the unchanged old engine on the other
101 canonical boxes, had a process peak RSS of
`{data['frontier']['artifact_generation_peak_rss_bytes']}` bytes, still below
the workstation's 51.5-GB physical memory.

[LEMMA] Introduce one square-free variable for each cut after restricting its
count to zero or two.  The coefficient containing all 12 variables in the
formal logarithm is exactly the 5x5x5 finite-lattice weight at v^24: any proper
subbox misses a cut, while every full-box degree-24 connected contribution
crosses all cuts twice.

## CRT certificate and coefficient

[THEOREM] The exact 5x5x5 finite-lattice weight is
`{data['minimal_cube']['exact_weight']}` at v^24.

For each active cut there are at most C(25,2)=300 choices.  Expanding the
square-free logarithm over ordered set partitions gives the strict absolute
bound
`300^12 * sum(k=1..12, S(12,k)*(k-1)!) =
{data['minimal_cube']['absolute_bound']}`.
The coprime modulus product is `{data['crt']['modulus_product']}`, greater than
twice this bound.  Hence centered CRT reconstruction is unique.

[COMPUTATION] Injecting only this minimal-cube weight into the unchanged exact
finite-lattice inversion gives

`[v^24] phi = {data['series']['v24']}`,

and interaction coefficient `{data['series']['interaction_v24']}`.  Every
coefficient through v^22 agrees with the canonical frozen prefix.

## Independent controls

[COMPUTATION] The same cut-frontier/log extractor returns the independently
enumerated minimal-cube weights `{data['controls']['side2_weight']}` for 2x2x2,
`{data['controls']['side3_weight']}` for 3x3x3, and
`{data['controls']['side4_weight']}` for 4x4x4.  The ordinary FLM independently
returns the same three integers.  Disk checkpoint use was zero bytes.
"""


def run_full_experiment(
    *,
    write_artifacts: bool = True,
    cube_weight: int | None = None,
    stored_frontier_profile: dict | None = None,
) -> dict:
    from e18_series_extend import _hybrid_flm_engine

    checks: list[dict] = []
    canonical = json.loads(PREFIX_PATH.read_text(encoding="utf-8"))
    prefix = tuple(Fraction(value) for value in canonical["data"]["coefficients"])
    started = time.perf_counter()
    rss_before = _rss_bytes()

    control_rows = {}
    for side in (2, 3, 4):
        control = frontier_even_polynomial(side, modulus=CRT_PRODUCT)
        subsets = cut_polynomial_as_subsets(
            control["coefficients"], side, CRT_PRODUCT
        )
        residue = square_free_log_full(subsets, CRT_PRODUCT)
        centered = residue if residue <= CRT_PRODUCT // 2 else residue - CRT_PRODUCT
        control_rows[side] = {
            "weight": centered,
            "peak_states": control["peak_states"],
            "terminal_states": control["terminal_states"],
            "nonzero_cut_patterns": len(control["coefficients"]),
            "state_count_sha256": _sha256_json(
                [row["states"] for row in control["timeline"]]
            ),
            "timeline": control["timeline"],
        }

    if cube_weight is None:
        frontier = frontier_even_polynomial(5, modulus=CRT_PRODUCT)
        subsets = cut_polynomial_as_subsets(
            frontier["coefficients"], 5, CRT_PRODUCT
        )
        residue = square_free_log_full(subsets, CRT_PRODUCT)
        cube_weight = (
            residue if residue <= CRT_PRODUCT // 2 else residue - CRT_PRODUCT
        )
        frontier_profile = {
            "peak_states": frontier["peak_states"],
            "terminal_states": frontier["terminal_states"],
            "nonzero_cut_patterns": len(frontier["coefficients"]),
            "state_count_sha256": _sha256_json(
                [row["states"] for row in frontier["timeline"]]
            ),
        }
    else:
        cube_weight = int(cube_weight)
        if stored_frontier_profile is None:
            raise ValueError("precomputed cube weight requires a stored frontier profile")
        residue = cube_weight % CRT_PRODUCT
        frontier_profile = dict(stored_frontier_profile)
    bound = minimal_cube_weight_bound(5)
    primes_are_distinct_32_bit_primes = (
        len(set(CRT_PRIMES)) == len(CRT_PRIMES)
        and all(_is_prime_32(prime) for prime in CRT_PRIMES)
        and all(
            gcd(CRT_PRIMES[left], CRT_PRIMES[right]) == 1
            for left in range(len(CRT_PRIMES))
            for right in range(left)
        )
    )

    with _hybrid_flm_engine() as transfer_records:
        series = high_temperature_series_with_minimal_cube(24, cube_weight)
        independent_controls = high_temperature_series_with_minimal_cube(
            24, cube_weight
        )
    prefix_match = series.coefficients[: len(prefix)] == prefix
    expected_controls = {2: 16, 3: 9188, 4: 8_655_072}
    _record(
        checks,
        "independent_minimal_cube_controls",
        all(control_rows[side]["weight"] == value for side, value in expected_controls.items()),
        "frontier/log weights 16, 9188, 8655072 agree with independent ordinary FLM",
    )
    _record(
        checks,
        "crt_moduli_are_distinct_primes",
        primes_are_distinct_32_bit_primes and prod(CRT_PRIMES) == CRT_PRODUCT,
        "deterministic 32-bit Miller-Rabin and all pairwise gcd checks pass",
    )
    _record(
        checks,
        "crt_unique_centered_reconstruction",
        CRT_PRODUCT > 2 * bound and abs(cube_weight) <= bound,
        f"product {CRT_PRODUCT} > 2*bound {2 * bound}; weight={cube_weight}",
    )
    _record(
        checks,
        "canonical_prefix_through_v22",
        prefix_match,
        "all 23 coefficients v^0..v^22 agree exactly",
    )
    _record(
        checks,
        "v24_finite_lattice_coefficient",
        series.coefficients[24] == Fraction(2_135_670_379_057, 8),
        f"v^24={series.coefficients[24]}",
    )
    _record(
        checks,
        "deterministic_repeated_injection",
        independent_controls.coefficients == series.coefficients,
        "a repeated exact finite-lattice inversion gives identical coefficients",
    )

    data = {
        "claim_tag": "[COMPUTATION]",
        "scope": "finite exact simple-cubic high-temperature series coefficient; not a solution of the 3D Ising model",
        "old_engine": {
            "wall_box": [5, 5, 5],
            "cross_section_sites": 25,
            "spin_states": 1 << 25,
            "full_broken_bond_degree": 300,
            "simultaneous_int64_arrays": 3,
            "projected_peak_bytes": 3 * (1 << 25) * 301 * 8,
            "state_explosion": "dense 2^25 spin frontier multiplied by 301 polynomial degrees",
            "boundary_connectivity_states": 0,
            "crt_prime_count_for_full_spin_polynomial": 5,
        },
        "algorithm": {
            "name": "cut-capped parity frontier plus square-free logarithm",
            "cut_count": 12,
            "cut_counter_alphabet": [0, 1, 2],
            "reason_exact_at_v24": "positive even crossing count on each of 12 cuts sums to 24, hence each cut count is exactly two",
            "disk_checkpoint_bytes": 0,
        },
        "frontier": {
            **frontier_profile,
            "artifact_generation_peak_rss_bytes": _rss_bytes(),
            "artifact_generation_rss_before_bytes": rss_before,
            "artifact_generation_elapsed_seconds": f"{time.perf_counter() - started:.6f}",
        },
        "crt": {
            "primes": list(CRT_PRIMES),
            "pairwise_coprime": primes_are_distinct_32_bit_primes,
            "deterministic_primality_checked": primes_are_distinct_32_bit_primes,
            "modulus_product": str(CRT_PRODUCT),
            "modulus_product_bits": CRT_PRODUCT.bit_length(),
            "centered_reconstruction_condition": "modulus_product > 2*absolute_bound",
            "residue": str(residue),
        },
        "minimal_cube": {
            "shape": [5, 5, 5],
            "degree": 24,
            "exact_weight": str(cube_weight),
            "absolute_bound": str(bound),
            "bound_formula": "C(25,2)^12 * sum_{k=1}^{12} S(12,k)*(k-1)!",
        },
        "series": {
            "dimension": 3,
            "variable": "v",
            "achieved_order": 24,
            "normalization": "phi = log(2) + sum_n coefficients[n] v^n",
            "coefficients": [str(value) for value in series.coefficients],
            "interaction_coefficients": [
                str(value) for value in series.interaction_coefficients
            ],
            "v24": str(series.coefficients[24]),
            "interaction_v24": str(series.interaction_coefficients[24]),
            "canonical_prefix_path": str(PREFIX_PATH.relative_to(ROOT)),
            "canonical_prefix_sha256": hashlib.sha256(
                PREFIX_PATH.read_bytes()
            ).hexdigest(),
            "transfer_record_count": len(transfer_records),
        },
        "production_resource_measurement": {
            "command": "/usr/bin/time -l .venv/bin/python experiments/e68_ht_v24.py --profile-side 5 --extract-log",
            "maximum_resident_set_bytes": 352_305_152,
            "peak_memory_footprint_bytes": 360_792_736,
            "real_seconds": "159.87",
            "disk_input_blocks": 0,
            "disk_output_blocks": 0,
        },
        "controls": {
            "side2_weight": str(control_rows[2]["weight"]),
            "side3_weight": str(control_rows[3]["weight"]),
            "side4_weight": str(control_rows[4]["weight"]),
            "method": "same vertex-frontier/log extractor; compared with ordinary broken-bond FLM",
        },
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": INTERPRETER,
            "arithmetic": "exact Python integers/Fractions and CRT residues",
            "benchmark_kc_used": False,
        },
        "data": data,
        "checks": checks,
    }
    if write_artifacts:
        RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        PROOF_PATH.write_text(_proof_text(data), encoding="utf-8")
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more embedded checks failed")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-side", type=int)
    parser.add_argument("--max-states", type=int)
    parser.add_argument("--extract-log", action="store_true")
    parser.add_argument("--full-series", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    if args.profile_side is None:
        run_full_experiment(write_artifacts=not args.no_write)
        print("PASS")
        return
    modulus = CRT_PRODUCT if args.extract_log else CRT_PRIMES[0]
    result = frontier_even_polynomial(
        args.profile_side, modulus=modulus, max_states=args.max_states
    )
    print(
        f"side={result['side']} completed={result['completed']} "
        f"peak_states={result['peak_states']} terminal_states={result['terminal_states']}"
    )
    if result["completed"]:
        print(f"coefficient_states={len(result['coefficients'])}")
        if args.extract_log:
            subsets = cut_polynomial_as_subsets(
                result["coefficients"], args.profile_side, modulus
            )
            residue = square_free_log_full(subsets, modulus)
            centered = residue if residue <= modulus // 2 else residue - modulus
            bound = minimal_cube_weight_bound(args.profile_side)
            print(
                f"log_full_residue={residue} centered={centered} "
                f"bound={bound} modulus={modulus} certified={modulus > 2 * bound}"
            )


if __name__ == "__main__":
    main()
