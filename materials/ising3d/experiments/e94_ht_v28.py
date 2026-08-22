"""Exact simple-cubic high-temperature series through v^28.

The producer uses a cut-capped vertex frontier.  Unlike the v^26 producer, the
cap is a vector: a requested half-crossing profile ``alpha`` uses cut caps
``2*alpha``.  This retains the ordinary mixed-radix logarithm terms needed for
both a six-crossing cut and two distinct four-crossing cuts at order 28.

Run from the repository root with::

    .venv/bin/python experiments/e94_ht_v28.py

The worker mode exists so every memory-intensive frontier runs in an isolated
process and reports an exact modular residue plus resource telemetry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, prod
from pathlib import Path
from typing import Sequence

import ising.series as flm
from ising.series import FLMSeries
from ising.transfer_matrix.crt import _descending_primes, _is_prime_32

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e94_ht_v28.py"
INTERPRETER = ".venv/bin/python"
ORDER = 28
V24_PATH = ROOT / "results" / "series" / "ht_v24.json"
V26_PATH = ROOT / "results" / "series" / "ht_v26.json"
RESULT_PATH = ROOT / "results" / "series" / "ht_v28.json"
PROOF_PATH = ROOT / "proofs" / "ht_v28.md"
MEMORY_BUDGET_BYTES = 40_000_000_000
# A worker trips at 90% of this process-local budget while both source and
# destination dictionaries coexist.  Sampling inside a vertex transition avoids
# an allocator kill before a post-vertex measurement can emit an unresolved
# certificate.
MEMORY_GUARD_TRIP_BYTES = MEMORY_BUDGET_BYTES * 9 // 10
MEMORY_RSS_SAMPLE_INSERTS = 65_536

# Six deterministically verified 31-bit primes already exceed every proved
# order-28 profile bound below.  Direct arithmetic modulo their product is
# equivalent to retaining the six CRT residues and avoids six frontier passes.
CRT_PRIMES = tuple(int(value) for value in _descending_primes(6))
CRT_PRODUCT = prod(CRT_PRIMES)

# [EXTERNAL] Arisue--Fujiwara Table I, arXiv:hep-lat/0209002.  This is checked
# only after the finite-lattice calculation has been completed; it is never an
# input to profile extraction or FLM assembly.
EXTERNAL_A28 = Fraction(525_549_581_866_326, 7)
EXTERNAL_TABLE = {
    "arxiv": "hep-lat/0209002",
    "table_i_sha256": "bf79ec7918694768b57e6793ed4b0c8f581aaea6a705b4273ecab42641723afc",
    "interaction_v28": str(EXTERNAL_A28),
}

WALLED_SHAPES = (
    (4, 6, 6),
    (4, 6, 7),
    (5, 5, 5),
    (5, 5, 6),
    (5, 5, 7),
    (5, 6, 6),
)

# The frontier advances in its first direction.  These permutations put a
# largest direction first, minimizing the pending parity register b*c.  A
# profile is generated in canonical sorted coordinates and then reoriented.
SWEEP_AXES = {
    (4, 6, 6): (1, 2, 0),  # (6, 6, 4), width 24
    (4, 6, 7): (2, 1, 0),  # (7, 6, 4), width 24
    (5, 5, 5): (0, 1, 2),
    (5, 5, 6): (2, 0, 1),  # (6, 5, 5), width 25
    (5, 5, 7): (2, 0, 1),  # (7, 5, 5), width 25
    (5, 6, 6): (1, 2, 0),  # (6, 6, 5), width 30
}


class FrontierResourceWall(RuntimeError):
    """A worker reached an explicit state or RSS guard before completion."""

    def __init__(self, record: dict):
        self.record = record
        reason = record.get("wall_reason", "unknown_resource_wall")
        super().__init__(
            "frontier did not complete profile "
            f"{record['shape']} degree {record['degree']} at "
            f"{record['peak_states']} states ({reason})"
        )


def _shape_sides(shape: Sequence[int]) -> tuple[int, int, int]:
    sides = tuple(int(side) for side in shape)
    if len(sides) != 3 or any(side < 1 for side in sides):
        raise ValueError("shape must have exactly three positive side lengths")
    return sides


def _cut_inventory(shape: Sequence[int]) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]]:
    """Return coordinate cuts and their crossing-edge counts."""

    a, b, c = _shape_sides(shape)
    sides = (a, b, c)
    cuts = tuple(
        (axis, position)
        for axis, side in enumerate(sides)
        for position in range(side - 1)
    )
    edge_counts = tuple(
        prod(sides[other] for other in range(3) if other != axis)
        for axis, _ in cuts
    )
    return cuts, edge_counts


def _radix_powers(radices: Sequence[int]) -> tuple[int, ...]:
    powers = []
    current = 1
    for radix in radices:
        radix = int(radix)
        if radix < 2:
            raise ValueError("mixed-radix bases must be at least two")
        powers.append(current)
        current *= radix
    return tuple(powers)


def _digits_of(code: int, radices: Sequence[int]) -> tuple[int, ...]:
    digits = []
    rest = int(code)
    for radix in radices:
        digits.append(rest % radix)
        rest //= radix
    if rest:
        raise AssertionError("mixed-radix code exceeds its declared range")
    return tuple(digits)


def _mixed_code(digits: Sequence[int], powers: Sequence[int]) -> int:
    if len(digits) != len(powers):
        raise ValueError("digits and powers have incompatible lengths")
    return sum(int(digit) * int(power) for digit, power in zip(digits, powers))


def _validate_cut_caps(shape: Sequence[int], cut_caps: Sequence[int] | None) -> tuple[int, ...]:
    cut_count = len(_cut_inventory(shape)[0])
    if cut_caps is None:
        return (2,) * cut_count
    caps = tuple(int(cap) for cap in cut_caps)
    if len(caps) != cut_count:
        raise ValueError("cut cap vector must have one entry per coordinate cut")
    if any(cap < 2 or cap % 2 for cap in caps):
        raise ValueError("every cut cap must be a positive even integer at least two")
    return caps


def frontier_cut_polynomial(
    shape: Sequence[int],
    *,
    modulus: int = CRT_PRODUCT,
    cut_caps: Sequence[int] | None = None,
    max_states: int | None = None,
    memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES,
) -> dict:
    """Enumerate even subgraphs with a separately capped counter per cut.

    The pending register has width ``b*c`` for a sweep in the first coordinate.
    Bit ``j`` is the parity already owed by an occupied edge ending ``j+1``
    vertices in the future.  Thus exactly one incoming parity bit is consumed at
    each vertex.  The cut counters are ordinary mixed-radix digits and are
    capped before a transition is admitted.  During a transition both the old
    and new dictionaries are live, so an RSS sample is taken periodically
    while the destination dictionary is built.
    """

    a, b, c = _shape_sides(shape)
    modulus = int(modulus)
    if modulus <= 1:
        raise ValueError("modulus must exceed one")
    if memory_budget_bytes is not None and int(memory_budget_bytes) < 1:
        raise ValueError("memory budget must be positive or None")
    memory_budget = (
        None if memory_budget_bytes is None else int(memory_budget_bytes)
    )
    memory_guard = (
        None
        if memory_budget is None
        else min(MEMORY_GUARD_TRIP_BYTES, memory_budget * 9 // 10)
    )
    caps = _validate_cut_caps((a, b, c), cut_caps)
    radices = tuple(cap + 1 for cap in caps)
    powers = _radix_powers(radices)
    width = b * c
    pending_mask = (1 << width) - 1
    states: dict[int, int] = {0: 1}
    peak = 1
    timeline = []
    started = time.perf_counter()

    def unresolved(
        reason: str,
        vertex: tuple[int, int, int],
        source_states: int,
        partial_updated_states: int,
    ) -> dict:
        return {
            "completed": False,
            "wall_reason": reason,
            "shape": [a, b, c],
            "cut_caps": list(caps),
            "modulus": str(modulus),
            "memory_budget_bytes": memory_budget,
            "memory_guard_trip_bytes": memory_guard,
            "rss_bytes": _rss_bytes(),
            "vertex": list(vertex),
            "source_states": source_states,
            "partial_updated_states": partial_updated_states,
            "peak_states": peak,
            "terminal_states": partial_updated_states,
            "timeline": timeline,
            "coefficients": None,
        }

    for x in range(a):
        for y in range(b):
            for z in range(c):
                directions: list[tuple[int, int]] = []
                if z + 1 < c:
                    directions.append((0, (a - 1) + (b - 1) + z))
                if y + 1 < b:
                    directions.append((c - 1, (a - 1) + y))
                if x + 1 < a:
                    directions.append((width - 1, x))

                options_by_parity: tuple[list[tuple[int, tuple[tuple[int, int, int], ...]]], list[tuple[int, tuple[tuple[int, int, int], ...]]]] = ([], [])
                for choice in range(1 << len(directions)):
                    parity = choice.bit_count() & 1
                    pending_delta = 0
                    increments = []
                    for index, (offset, cut) in enumerate(directions):
                        if choice & (1 << index):
                            pending_delta ^= 1 << offset
                            increments.append((powers[cut], radices[cut], caps[cut]))
                    options_by_parity[parity].append(
                        (pending_delta, tuple(increments))
                    )

                source_states = len(states)
                inserted_states = 0
                updated: dict[int, int] = {}
                for packed, value in states.items():
                    pending = packed & pending_mask
                    code = packed >> width
                    incoming = pending & 1
                    base_pending = pending >> 1
                    for pending_delta, increments in options_by_parity[incoming]:
                        new_code = code
                        allowed = True
                        for power, radix, cap in increments:
                            if (new_code // power) % radix == cap:
                                allowed = False
                                break
                            new_code += power
                        if not allowed:
                            continue
                        key = (new_code << width) | (base_pending ^ pending_delta)
                        previous = updated.get(key)
                        if previous is None:
                            updated[key] = value
                            inserted_states += 1
                            if (
                                memory_guard is not None
                                and inserted_states % MEMORY_RSS_SAMPLE_INSERTS == 0
                                and _rss_bytes() >= memory_guard
                            ):
                                peak = max(peak, len(updated))
                                return unresolved(
                                    "rss_budget",
                                    (x, y, z),
                                    source_states,
                                    len(updated),
                                )
                        else:
                            total = previous + value
                            updated[key] = total if total < modulus else total - modulus
                peak = max(peak, len(updated))
                if max_states is not None and len(updated) > int(max_states):
                    return unresolved(
                        "state_cap",
                        (x, y, z),
                        source_states,
                        len(updated),
                    )
                if memory_guard is not None and _rss_bytes() >= memory_guard:
                    return unresolved(
                        "rss_budget",
                        (x, y, z),
                        source_states,
                        len(updated),
                    )
                states = updated
                timeline.append(
                    {
                        "vertex": [x, y, z],
                        "states": len(states),
                        "elapsed_seconds": round(time.perf_counter() - started, 6),
                    }
                )

    coefficients: dict[int, int] = defaultdict(int)
    for packed, value in states.items():
        pending = packed & pending_mask
        if pending:
            raise AssertionError("nonempty parity register after the final vertex")
        code = packed >> width
        digits = _digits_of(code, radices)
        if any(digit % 2 for digit in digits):
            raise AssertionError("an even subgraph crossed a coordinate cut oddly")
        prior = coefficients[code]
        total = prior + value
        coefficients[code] = total if total < modulus else total - modulus
    return {
        "completed": True,
        "shape": [a, b, c],
        "cut_caps": list(caps),
        "modulus": str(modulus),
        "memory_budget_bytes": memory_budget,
        "memory_guard_trip_bytes": memory_guard,
        "peak_states": peak,
        "terminal_states": len(states),
        "timeline": timeline,
        "coefficients": dict(coefficients),
    }


def cut_polynomial_as_half_crossings(
    coefficients: dict[int, int],
    shape: Sequence[int],
    cut_caps: Sequence[int],
    modulus: int,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Convert terminal even crossing counters to an ordinary exponent quotient."""

    caps = _validate_cut_caps(shape, cut_caps)
    terminal_radices = tuple(cap + 1 for cap in caps)
    exponent_radices = tuple(cap // 2 + 1 for cap in caps)
    powers = _radix_powers(exponent_radices)
    result = [0] * prod(exponent_radices)
    for code, value in coefficients.items():
        counts = _digits_of(code, terminal_radices)
        if any(count % 2 for count in counts):
            raise AssertionError("terminal cut count must be even")
        exponent_code = _mixed_code(tuple(count // 2 for count in counts), powers)
        total = result[exponent_code] + value
        result[exponent_code] = total if total < modulus else total - modulus
    if result[0] != 1:
        raise AssertionError("the empty even subgraph must be unique")
    return tuple(result), exponent_radices


def truncated_mixed_log(
    coefficient: Sequence[int], radices: Sequence[int], modulus: int
) -> list[int]:
    """Compute formal log in a coordinatewise-complete mixed-radix quotient.

    The recurrence comes from ``partial_j P = (partial_j log(P))*P``.  Unlike a
    square-free quotient it retains, for example, a 2+2 split of one cut across
    two logarithm factors.  A cap vector supplies every exponent coordinatewise
    below the requested profile, so no omitted monomial can enter the recurrence.
    """

    radices = tuple(int(radix) for radix in radices)
    expected_size = prod(radices)
    if len(coefficient) != expected_size or coefficient[0] % modulus != 1:
        raise ValueError("mixed logarithm requires a complete unit-constant quotient")
    powers = _radix_powers(radices)
    digits_by_code = tuple(_digits_of(code, radices) for code in range(expected_size))
    logarithm = [0] * expected_size
    ordered_codes = sorted(
        range(1, expected_size),
        key=lambda code: (sum(digits_by_code[code]), code),
    )
    for alpha_code in ordered_codes:
        alpha = digits_by_code[alpha_code]
        axis = next(index for index, digit in enumerate(alpha) if digit)
        if gcd(alpha[axis], modulus) != 1:
            raise AssertionError("a log denominator is not invertible modulo the CRT product")
        subtotal = 0
        gamma = [0] * len(radices)

        def visit(index: int, gamma_code: int) -> None:
            nonlocal subtotal
            if index == len(radices):
                if gamma_code != alpha_code:
                    delta_code = alpha_code - gamma_code
                    subtotal += (
                        gamma[axis]
                        * logarithm[gamma_code]
                        * coefficient[delta_code]
                    )
                return
            lower = 1 if index == axis else 0
            for digit in range(lower, alpha[index] + 1):
                gamma[index] = digit
                visit(index + 1, gamma_code + digit * powers[index])

        visit(0, 0)
        logarithm[alpha_code] = (
            coefficient[alpha_code]
            - subtotal * pow(alpha[axis], -1, modulus)
        ) % modulus
    return logarithm


def ordered_partition_factor(slot_count: int) -> int:
    """Return ``sum_k S(slot_count,k)*(k-1)!`` exactly."""

    if slot_count < 1:
        raise ValueError("a profile needs at least one half-crossing slot")
    stirling = [[0] * (slot_count + 1) for _ in range(slot_count + 1)]
    stirling[0][0] = 1
    for row in range(1, slot_count + 1):
        for blocks in range(1, row + 1):
            stirling[row][blocks] = (
                stirling[row - 1][blocks - 1]
                + blocks * stirling[row - 1][blocks]
            )
    return sum(
        stirling[slot_count][blocks] * factorial(blocks - 1)
        for blocks in range(1, slot_count + 1)
    )


def frontier_profile_bound(shape: Sequence[int], profile: Sequence[int]) -> int:
    """A strict labelled-slot absolute bound for one formal-log profile."""

    _, edge_counts = _cut_inventory(shape)
    profile = tuple(int(exponent) for exponent in profile)
    if len(profile) != len(edge_counts) or any(exponent < 1 for exponent in profile):
        raise ValueError("profile must contain one positive exponent per cut")
    # For all production cuts E>=24 and exponents <=3, C(E,2r) increases with r.
    # Label each of the r half-crossing slots and bound every slot by C(E,2r).
    # Ordered set partitions of all slots contribute the O_q factor.
    if any(2 * exponent > edge_count for edge_count, exponent in zip(edge_counts, profile)):
        raise ValueError("profile asks for more crossing edges than a cut has")
    return ordered_partition_factor(sum(profile)) * prod(
        comb(edge_count, 2 * exponent) ** exponent
        for edge_count, exponent in zip(edge_counts, profile)
    )


def _centered(residue: int, modulus: int) -> int:
    return residue if residue <= modulus // 2 else residue - modulus


def extract_frontier_profile(
    shape: Sequence[int],
    profile: Sequence[int],
    *,
    modulus: int = CRT_PRODUCT,
    max_states: int | None = None,
    memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES,
) -> dict:
    """Extract one exact log coefficient indexed by a primitive half profile."""

    shape = _shape_sides(shape)
    profile = tuple(int(exponent) for exponent in profile)
    cut_count = len(_cut_inventory(shape)[0])
    if len(profile) != cut_count or any(exponent < 1 for exponent in profile):
        raise ValueError("profile must contain one positive half-crossing exponent per cut")
    if gcd(*profile) != 1:
        raise ValueError("log-profile exponents must be primitive")
    caps = tuple(2 * exponent for exponent in profile)
    started = time.perf_counter()
    frontier = frontier_cut_polynomial(
        shape,
        modulus=modulus,
        cut_caps=caps,
        max_states=max_states,
        memory_budget_bytes=memory_budget_bytes,
    )
    wall = time.perf_counter() - started
    common = {
        "shape": frontier["shape"],
        "profile": list(profile),
        "cut_caps": list(caps),
        "modulus": str(modulus),
        "peak_states": frontier["peak_states"],
        "terminal_states": frontier["terminal_states"],
        "timeline": frontier["timeline"],
        "wall_seconds": f"{wall:.6f}",
        "bound": str(frontier_profile_bound(shape, profile)),
        "primitive_profile": True,
        "memory_budget_bytes": frontier["memory_budget_bytes"],
        "memory_guard_trip_bytes": frontier["memory_guard_trip_bytes"],
    }
    if not frontier["completed"]:
        return {
            "completed": False,
            **common,
            **{
                key: frontier[key]
                for key in (
                    "wall_reason",
                    "memory_budget_bytes",
                    "memory_guard_trip_bytes",
                    "rss_bytes",
                    "vertex",
                    "source_states",
                    "partial_updated_states",
                )
            },
        }
    half_polynomial, exponent_radices = cut_polynomial_as_half_crossings(
        frontier["coefficients"], shape, caps, modulus
    )
    logarithm = truncated_mixed_log(half_polynomial, exponent_radices, modulus)
    target = _mixed_code(profile, _radix_powers(exponent_radices))
    residue = logarithm[target]
    return {
        "completed": True,
        **common,
        "residue": str(residue),
        "centered": str(_centered(residue, modulus)),
    }


def _cut_actions(shape: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """All box isometries: independent reflections and equal-side permutations."""

    sides = _shape_sides(shape)
    cuts, _ = _cut_inventory(sides)
    lookup = {cut: index for index, cut in enumerate(cuts)}
    actions = []
    for axis_permutation in permutations(range(3)):
        if any(sides[axis] != sides[axis_permutation[axis]] for axis in range(3)):
            continue
        for reflected in product((False, True), repeat=3):
            actions.append(
                tuple(
                    lookup[
                        (
                            axis_permutation[axis],
                            sides[axis] - 2 - position if reflected[axis] else position,
                        )
                    ]
                    for axis, position in cuts
                )
            )
    return tuple(actions)


def _act_profile(profile: Sequence[int], action: Sequence[int]) -> tuple[int, ...]:
    output = [0] * len(profile)
    for source, target in enumerate(action):
        output[target] = int(profile[source])
    return tuple(output)


def _positive_profiles(cut_count: int, extras: int) -> set[tuple[int, ...]]:
    """All positive profiles needed here (only 0, 1, or 2 excess slots)."""

    if extras < 0 or extras > 2:
        raise ValueError("order-28 producer only supports zero, one, or two excess slots")
    profiles: set[tuple[int, ...]] = set()
    if extras == 0:
        profiles.add((1,) * cut_count)
    elif extras == 1:
        for index in range(cut_count):
            profile = [1] * cut_count
            profile[index] = 2
            profiles.add(tuple(profile))
    else:
        for index in range(cut_count):
            profile = [1] * cut_count
            profile[index] = 3
            profiles.add(tuple(profile))
        for left, right in combinations(range(cut_count), 2):
            profile = [1] * cut_count
            profile[left] = profile[right] = 2
            profiles.add(tuple(profile))
    return profiles


def profile_orbits(shape: Sequence[int], degree: int) -> list[dict]:
    """Group every positive degree profile under exact box symmetries."""

    shape = _shape_sides(shape)
    cut_count = len(_cut_inventory(shape)[0])
    if degree % 2:
        raise ValueError("even-subgraph degree must be even")
    extras = degree // 2 - cut_count
    remaining = _positive_profiles(cut_count, extras)
    actions = _cut_actions(shape)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = {_act_profile(seed, action) for action in actions}
        representative = min(orbit)
        orbits.append(
            {
                "profile": representative,
                "multiplicity": len(orbit),
                "orbit": tuple(sorted(orbit)),
            }
        )
        remaining -= orbit
    return sorted(orbits, key=lambda row: (row["profile"], row["multiplicity"]))


def _reorient_profile(
    profile: Sequence[int], source_shape: Sequence[int], destination_axes: Sequence[int]
) -> tuple[int, ...]:
    """Reorder profile segments when the sweep directions are permuted."""

    source = _shape_sides(source_shape)
    axes = tuple(int(axis) for axis in destination_axes)
    if sorted(axes) != [0, 1, 2]:
        raise ValueError("destination axes must be a permutation of 0,1,2")
    if len(profile) != sum(side - 1 for side in source):
        raise ValueError("profile length does not match source shape")
    segments = []
    position = 0
    for side in source:
        segments.append(tuple(int(value) for value in profile[position : position + side - 1]))
        position += side - 1
    return tuple(value for axis in axes for value in segments[axis])


def _sweep_shape(canonical_shape: Sequence[int]) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    canonical = _shape_sides(canonical_shape)
    axes = SWEEP_AXES.get(canonical, (0, 1, 2))
    return tuple(canonical[axis] for axis in axes), axes


def _rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _sha256_json(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def _profile_record(
    canonical_shape: tuple[int, int, int],
    degree: int,
    canonical_profile: tuple[int, ...],
    multiplicity: int,
    sweep_shape: tuple[int, int, int],
    sweep_profile: tuple[int, ...],
    payload: dict,
) -> dict:
    row = {
        "shape": list(canonical_shape),
        "degree": int(degree),
        "profile": list(canonical_profile),
        "multiplicity": int(multiplicity),
        "sweep_shape": list(sweep_shape),
        "sweep_profile": list(sweep_profile),
        "cut_caps": payload["cut_caps"],
        "completed": payload["completed"],
        "peak_states": payload["peak_states"],
        "terminal_states": payload["terminal_states"],
        "timeline": payload["timeline"],
        "wall_seconds": payload["wall_seconds"],
        "max_rss_bytes": payload["max_rss_bytes"],
        "state_count_sha256": _sha256_json(
            [entry["states"] for entry in payload["timeline"]]
        ),
        "bound": payload["bound"],
        "modulus": payload["modulus"],
    }
    for key in (
        "primitive_profile",
        "memory_budget_bytes",
        "memory_guard_trip_bytes",
        "wall_reason",
        "rss_bytes",
        "vertex",
        "source_states",
        "partial_updated_states",
    ):
        if key in payload:
            row[key] = payload[key]
    if "parent_wall_seconds" in payload:
        row["parent_wall_seconds"] = payload["parent_wall_seconds"]
    if payload["completed"]:
        row["residue"] = payload["residue"]
        row["centered"] = payload["centered"]
    return row


def run_frontier_worker(
    shape: Sequence[int],
    profile: Sequence[int],
    *,
    timeout_seconds: int,
    max_states: int | None = None,
    memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES,
    one_prime: bool = False,
) -> dict:
    """Run one profile extraction outside the producer's accumulated RSS."""

    command = [
        sys.executable,
        str(ROOT / SCRIPT),
        "--worker",
        "--shape",
        ",".join(str(side) for side in shape),
        "--half-profile",
        ",".join(str(exponent) for exponent in profile),
    ]
    if max_states is not None:
        command += ["--max-states", str(max_states)]
    if memory_budget_bytes is not None:
        command += ["--memory-budget-bytes", str(memory_budget_bytes)]
    if one_prime:
        command.append("--one-prime")
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            cwd=ROOT,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"frontier worker exceeded its {timeout_seconds}-second timeout for {shape}"
        ) from error
    parent_wall = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"frontier worker failed for {shape}: {completed.stderr.strip()[-2000:]}"
        )
    try:
        payload = json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"frontier worker emitted no JSON for {shape}: {completed.stdout[-1000:]}"
        ) from error
    payload["parent_wall_seconds"] = f"{parent_wall:.6f}"
    payload["worker_command"] = command
    return payload


def _run_profile_family(
    canonical_shape: Sequence[int],
    degree: int,
    *,
    worker: bool,
    timeout_seconds: int = 14_400,
    max_states: int | None = None,
    memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES,
) -> tuple[int, list[dict]]:
    """Compute and symmetry-sum every profile contributing at one box degree."""

    canonical = tuple(sorted(_shape_sides(canonical_shape)))
    sweep, axes = _sweep_shape(canonical)
    records = []
    total = 0
    for orbit in profile_orbits(canonical, degree):
        profile = orbit["profile"]
        sweep_profile = _reorient_profile(profile, canonical, axes)
        if worker:
            payload = run_frontier_worker(
                sweep,
                sweep_profile,
                timeout_seconds=timeout_seconds,
                max_states=max_states,
                memory_budget_bytes=memory_budget_bytes,
            )
        else:
            payload = extract_frontier_profile(
                sweep,
                sweep_profile,
                max_states=max_states,
                memory_budget_bytes=memory_budget_bytes,
            )
            payload["max_rss_bytes"] = _rss_bytes()
        row = _profile_record(
            canonical,
            degree,
            profile,
            orbit["multiplicity"],
            sweep,
            sweep_profile,
            payload,
        )
        records.append(row)
        if not payload["completed"]:
            raise FrontierResourceWall(row)
        total += orbit["multiplicity"] * int(payload["centered"])
    return total, records


def high_temperature_series_with_injections(
    order: int, injections: dict[tuple[int, int, int], dict[int, int]]
) -> FLMSeries:
    """Exact FLM inversion with only cross-section-walled weights injected."""

    order = int(order)
    shapes = flm._ht_shapes(3, order, 0)
    normalized = {
        tuple(sorted(tuple(int(side) for side in shape))): {
            int(degree): int(value) for degree, value in values.items()
        }
        for shape, values in injections.items()
    }
    weights = {}
    interaction = [Fraction(0) for _ in range(order + 1)]
    injected_seen = set()
    for shape in shapes:
        canonical = flm._canonical_shape(shape)
        minimum_degree = 2 * sum(side - 1 for side in canonical)
        if canonical in normalized:
            exact_weight = [Fraction(0) for _ in range(order + 1)]
            for degree, value in normalized[canonical].items():
                if degree % 2 or not minimum_degree <= degree <= order:
                    raise AssertionError("injected weight contradicts the cut-parity degree lemma")
                exact_weight[degree] = Fraction(value)
            injected_seen.add(canonical)
        else:
            exact_weight = list(flm._ht_box_log(canonical, order))
            for subshape, subweight in weights.items():
                if all(sub <= side for sub, side in zip(subshape, shape)):
                    placements = prod(
                        side - sub + 1 for sub, side in zip(subshape, shape)
                    )
                    for degree in range(order + 1):
                        exact_weight[degree] -= placements * subweight[degree]
        frozen = tuple(exact_weight)
        weights[shape] = frozen
        for degree, value in enumerate(frozen):
            interaction[degree] += value
    if injected_seen != set(normalized):
        raise AssertionError("an injected walled box was absent from the FLM family")
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


def _assemble(order: int, injections: dict[tuple[int, int, int], dict[int, int]]) -> tuple[FLMSeries, dict]:
    """Use the existing CRT spin-transfer route for all non-walled boxes."""

    from e18_series_extend import _hybrid_flm_engine

    with _hybrid_flm_engine() as transfer_records:
        series = high_temperature_series_with_injections(order, injections)
        records = dict(transfer_records)
    return series, records


def _engine_projection(shape: Sequence[int]) -> dict:
    sides = tuple(sorted(_shape_sides(shape)))
    cross = sides[0] * sides[1]
    columns = flm._box_bond_count(sides) + 1
    return {
        "shape": list(sides),
        "cross_section_sites": cross,
        "projected_dense_bytes": 3 * (1 << cross) * columns * 8,
    }


def _profile_summary(records: Sequence[dict]) -> dict:
    return {
        "profile_count": len(records),
        "max_peak_states": max(record["peak_states"] for record in records),
        "max_worker_rss_bytes": max(record["max_rss_bytes"] for record in records),
        "sum_worker_wall_seconds": f"{sum(float(record['wall_seconds']) for record in records):.6f}",
    }


def _proof_text(data: dict) -> str:
    box = data["box_class_lemma"]
    weights = data.get("weights", {})
    series = data.get("series", {})
    crt = data["crt"]
    controls = data.get("controls", {})
    resource = data["production_resource_measurement"]
    if data.get("status") == "[UNRESOLVED]":
        return f"""# Exact simple-cubic HT v^28 frontier status

## Scope

[UNRESOLVED] The order-28 finite calculation did not complete within the
recorded frontier resource guard.  This note makes no claim for `[v^28]`.

## Validated lower-order generalized-register work

[COMPUTATION] The generic cut-cap vector reproduced the stored exact v^24 and
v^26 series before the order-28 worker reached its wall.  The partial profile
telemetry and exact state cap are recorded in
`results/series/ht_v28.json`.

## Box-class lemma

[LEMMA] Every edge crosses exactly one coordinate cut and an even subgraph
crosses every cut evenly.  Therefore a connected exact-bounding-box graph has
at least `2*sum_i(side_i-1)` edges.  At order 28 only boxes with
`a+b+c <= 17` can contribute.

## Resource certificate

[COMPUTATION] {data['resource_wall']}
"""
    return f"""# Exact simple-cubic HT coefficient through v^28

## Scope

[COMPUTATION] This is a finite exact coefficient calculation, not a solution of
the three-dimensional Ising model.

## Exact order-28 box-class bound

[LEMMA] Enting's finite-lattice identity is
`L(A)=sum_(R<=A) product_i(A_i-R_i+1) W(R)`.  Every edge of an
`a*b*c` box crosses exactly one of its `a+b+c-3` coordinate cuts.  An even
subgraph crosses each cut evenly; an exact-bounding-box connected contribution
crosses each cut positively.  Hence it has at least
`2*(a+b+c-3)` edges.  At truncation order 28 this proves the exact class bound
`a+b+c <= 17`, not merely a heuristic cutoff.  The artifact enumerates all
{box['contributing_canonical_class_count']} canonical classes, of which
{box['minimal_span_class_count']} have sum 17 and first occur at v^28.

[COMPUTATION] The six and only six classes whose canonical open spin-transfer
cross-sections exceed 22 are `{box['walled']}` with cross-sections
`{box['walled_cross_sections']}`.  Their dense broken-bond projections are
`{box['walled_projected_bytes']}` bytes, so their weights are instead computed
by the sparse parity frontier.

## Generalized cut-profile frontier

[LEMMA] For a requested positive half-crossing profile `alpha`, cap cut `i` at
`2*alpha_i` crossings and retain the complete coordinatewise quotient with
radix `alpha_i+1`.  The pending parity register has width `b*c` in a sweep of
an `a*b*c` box; it enforces even vertex degree exactly.  The formal-log
recurrence in this ordinary mixed-radix quotient includes all decompositions of
`alpha`, including a two-plus-two split of a four-crossing cut.  Counts above a
cap cannot enter the target coefficient, so the cap is exact.

[LEMMA] At degree 28, a sum-17 box has the all-ones profile; each sum-16 box
has exactly one exponent two; and the 5x5x5 box has either one exponent three
(one six-crossing cut) or two distinct exponents two (two four-crossing cuts).
The latter requires exponent radix four.  Box reflections and permutations of
equal sides give two one-six-crossing cube orbits (edge/middle, multiplicity
six each) and seven two-four-crossing orbits: same-axis opposite edge-edge
(multiplicity 3), same-axis middle-middle (3), same-axis edge-middle adjacent
(6), same-axis edge-middle nonadjacent (6), and different-axis edge-edge (12),
middle-middle (12), edge-middle (24).  The cut-coordinate distance distinction
is retained; no inequivalent profiles were merged.

## Primitive log profiles

[LEMMA] Let `P` be an integer-coefficient formal power series with constant
term one.  If a multi-index `alpha` is primitive
(`gcd(alpha_i)=1`), then `[x^alpha] log P` is an integer.  Indeed, in the
`k`-fold term of the logarithm expansion, cyclic rotation acts on the ordered
`k` factors.  A nontrivial stabilizer would repeat a shorter word and hence
divide every coordinate of `alpha`; primitivity makes every orbit have size
`k`, cancelling the displayed denominator `1/k` orbit by orbit.

[COMPUTATION] Every listed v^24, v^26, and v^28 cut profile has gcd one.  Thus
the reconstructed centered residues are integral formal-log coefficients before
their exact symmetry multiplicities are applied.

## CRT certificate

[LEMMA] Put `O_q=sum_(k=1..q) S(q,k)(k-1)!`.  For a profile `alpha`, label its
`q=sum_i alpha_i` half-crossing slots.  Ignoring vertex parity, a slot on a cut
with `E_i` crossing edges is bounded by `C(E_i,2*alpha_i)`.  Expanding the
formal logarithm over ordered set partitions gives the strict bound
`O_q product_i C(E_i,2*alpha_i)^alpha_i`.  The profile-wise values are stored
verbatim in the JSON artifact.  The product of the six independently primality
checked 31-bit primes is `{crt['modulus_product']}`, and it exceeds twice every
stored profile bound; centered reconstruction is therefore unique.

## Lower-order reproduction and independent controls

[COMPUTATION] Before any degree-28 extraction, the new cap-vector pipeline
recomputed the v^24/v^26 walled weights and fresh FLM assemblies.  It reproduced
`[v^24] phi = {data['reproductions']['v24']['v']}` and
`[v^26] phi = {data['reproductions']['v26']['v']}` exactly.  Independent
spin-transfer FLM controls for cap-2, cap-4, and the new cap-6/pair profile
families are `{controls['independent_spin_transfer']}`.

## Exact coefficient

[COMPUTATION] The injected walled weights are `{weights['injections']}`.
The exact finite-lattice assembly gives

`[v^28] phi = {series['v28']}`,

with interaction coefficient `a_28 = {series['interaction_v28']}`.  The odd
coefficient v^27 vanishes and every coefficient through v^26 reproduces the
previous exact artifact.

[EXTERNAL] Arisue and Fujiwara, arXiv:hep-lat/0209002 Table I tabulates
`a_28 = {EXTERNAL_TABLE['interaction_v28']}`.  This value was used only after
computation as a falsification witness; it agrees exactly.  Recorded PDF sha256:
`{EXTERNAL_TABLE['table_i_sha256']}`.

## Resource telemetry

[COMPUTATION] Frontier workers were isolated.  Their summary is
`{data['frontier']['summary']}`.  The producer peak RSS was
`{resource['maximum_resident_set_bytes']}` bytes and elapsed wall time was
`{resource['real_seconds']}` seconds.  No disk checkpoint was used.

[COMPUTATION] The finalized worker samples RSS while both source and
destination frontier dictionaries coexist and stops with an `[UNRESOLVED]`
certificate at `{MEMORY_GUARD_TRIP_BYTES}` bytes rather than relying on an
allocator failure.  This post-production guard is not used as evidence for the
coefficient above.
"""


def _base_box_data() -> dict:
    canonical = tuple(
        sorted({flm._canonical_shape(shape) for shape in flm._ht_shapes(3, ORDER, 0)})
    )
    minimal = tuple(shape for shape in canonical if sum(shape) == 17)
    walled = tuple(shape for shape in canonical if shape[0] * shape[1] > 22)
    if walled != WALLED_SHAPES:
        raise AssertionError("walled class inventory disagrees with the transfer guard")
    return {
        "order": ORDER,
        "contributing_classes": "a+b+c <= 17",
        "contributing_canonical_class_count": len(canonical),
        "contributing_canonical_classes": [list(shape) for shape in canonical],
        "minimal_span_class_count": len(minimal),
        "minimal_span_classes": [list(shape) for shape in minimal],
        "walled": [list(shape) for shape in walled],
        "walled_cross_sections": [shape[0] * shape[1] for shape in walled],
        "walled_projected_bytes": [_engine_projection(shape)["projected_dense_bytes"] for shape in walled],
    }


def _read_reference_series(path: Path) -> tuple[Fraction, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return tuple(Fraction(value) for value in payload["data"]["series"]["coefficients"])


def _lower_controls(
    *, memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES
) -> tuple[dict, list[dict]]:
    """Run generic caps before any walled order-28 frontier is launched."""

    cube333, rows333 = _run_profile_family(
        (3, 3, 3), 16, worker=False, memory_budget_bytes=memory_budget_bytes
    )
    cube444_minimal, rows444_minimal = _run_profile_family(
        (4, 4, 4), 18, worker=False, memory_budget_bytes=memory_budget_bytes
    )
    cube444_cap4, rows444_cap4 = _run_profile_family(
        (4, 4, 4), 20, worker=False, memory_budget_bytes=memory_budget_bytes
    )
    # All these boxes fit the independent established spin-transfer FLM route.
    spin333 = flm.high_temperature_free_energy(3, 16).box_weights[(3, 3, 3)][16]
    from e18_series_extend import _hybrid_flm_engine

    with _hybrid_flm_engine():
        spin444_minimal = flm.high_temperature_free_energy(3, 18).box_weights[(4, 4, 4)][18]
        spin444_cap4 = flm.high_temperature_free_energy(3, 20).box_weights[(4, 4, 4)][20]
    controls = {
        "cube333_degree16": {
            "aggregate": str(cube333),
            "spin_transfer": str(spin333),
            "profile_orbit_count": len(rows333),
        },
        "cube444_degree18": {
            "aggregate": str(cube444_minimal),
            "spin_transfer": str(spin444_minimal),
            "profile_orbit_count": len(rows444_minimal),
        },
        "cube444_degree20": {
            "aggregate": str(cube444_cap4),
            "spin_transfer": str(spin444_cap4),
            "profile_orbit_count": len(rows444_cap4),
        },
    }
    return controls, rows333 + rows444_minimal + rows444_cap4


def _make_payload(data: dict, checks: list[dict]) -> dict:
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": INTERPRETER,
            "arithmetic": "exact Python integers/Fractions and direct CRT-product residues",
            "benchmark_kc_used": False,
        },
        "data": data,
        "checks": checks,
    }


def run_full_experiment(
    *,
    write_artifacts: bool = True,
    max_states: int | None = None,
    memory_budget_bytes: int | None = MEMORY_BUDGET_BYTES,
) -> dict:
    """Reproduce v24/v26 first, then compute and assemble the v28 coefficient."""

    checks: list[dict] = []
    started = time.perf_counter()
    rss_before = _rss_bytes()
    box_data = _base_box_data()
    v24_reference = _read_reference_series(V24_PATH)
    v26_reference = _read_reference_series(V26_PATH)

    # ---- generic-register controls before production -----------------------
    controls, control_rows = _lower_controls(
        memory_budget_bytes=memory_budget_bytes
    )
    controls_ok = all(
        Fraction(row["aggregate"]) == Fraction(row["spin_transfer"])
        for row in controls.values()
    )
    _record(
        checks,
        "independent_generic_profile_controls",
        controls_ok,
        "generic cap-6/pair, cap-2, and cap-4 profile aggregates agree with spin-transfer FLM",
    )
    if not controls_ok:
        raise AssertionError("generic-register control failed before walled production")

    # ---- recompute every lower-order walled input from the new code --------
    profile_records: list[dict] = []
    w555_24, rows = _run_profile_family(
        (5, 5, 5),
        24,
        worker=True,
        max_states=max_states,
        memory_budget_bytes=memory_budget_bytes,
    )
    profile_records += rows
    w555_26, rows = _run_profile_family(
        (5, 5, 5),
        26,
        worker=True,
        max_states=max_states,
        memory_budget_bytes=memory_budget_bytes,
    )
    profile_records += rows
    w556_26, rows = _run_profile_family(
        (5, 5, 6),
        26,
        worker=True,
        max_states=max_states,
        memory_budget_bytes=memory_budget_bytes,
    )
    profile_records += rows
    w466_26, rows = _run_profile_family(
        (4, 6, 6),
        26,
        worker=True,
        max_states=max_states,
        memory_budget_bytes=memory_budget_bytes,
    )
    profile_records += rows

    v26_injections = {
        (5, 5, 5): {24: w555_24, 26: w555_26},
        (5, 5, 6): {26: w556_26},
        (4, 6, 6): {26: w466_26},
    }
    # One fresh degree-26 FLM inversion contains the complete degree-24 prefix,
    # so this independently reproduces both stored lower-order artifacts without
    # paying a second transfer-matrix sweep.
    series26, _ = _assemble(26, v26_injections)
    lower_match = (
        series26.coefficients[: len(v24_reference)] == v24_reference
        and series26.coefficients == v26_reference
    )
    _record(
        checks,
        "fresh_v24_v26_pipeline_reproduction",
        lower_match,
        f"fresh cap-vector FLM gives v24={series26.coefficients[24]} and v26={series26.coefficients[26]}",
    )
    if not lower_match:
        raise AssertionError("new generic pipeline failed lower-order reproduction")

    # Existing certified values are comparison targets only, after extraction.
    old_v24 = json.loads(V24_PATH.read_text(encoding="utf-8"))["data"]["minimal_cube"]
    old_v26 = json.loads(V26_PATH.read_text(encoding="utf-8"))["data"]["weights"]
    lower_weight_match = (
        w555_24 == int(old_v24["exact_weight"])
        and w555_26 == int(old_v26["w_555_26"])
        and w556_26 == int(old_v26["w_655"])
        and w466_26 == int(old_v26["w_664"])
    )
    _record(
        checks,
        "fresh_walled_weight_reproduction",
        lower_weight_match,
        "new generalized frontier reproduces all v24/v26 walled weights before v28",
    )
    if not lower_weight_match:
        raise AssertionError("new generic frontier failed lower-order walled-weight reproduction")

    # ---- order-28 profiles: all needed walled classes ----------------------
    try:
        w467_28, rows = _run_profile_family(
            (4, 6, 7),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
        w557_28, rows = _run_profile_family(
            (5, 5, 7),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
        w566_28, rows = _run_profile_family(
            (5, 6, 6),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
        w466_28, rows = _run_profile_family(
            (4, 6, 6),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
        w556_28, rows = _run_profile_family(
            (5, 5, 6),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
        w555_28, rows = _run_profile_family(
            (5, 5, 5),
            28,
            worker=True,
            max_states=max_states,
            memory_budget_bytes=memory_budget_bytes,
        )
        profile_records += rows
    except FrontierResourceWall as wall:
        profile_records.append(wall.record)
        data = {
            "status": "[UNRESOLVED]",
            "claim_tag": "[UNRESOLVED]",
            "scope": "generic exact cut-profile frontier reached an explicit resource cap; no v28 coefficient is claimed",
            "box_class_lemma": box_data,
            "profiles": profile_records,
            "crt": {
                "primes": list(CRT_PRIMES),
                "modulus_product": str(CRT_PRODUCT),
                "centered_reconstruction_condition": "modulus_product > 2*absolute_profile_bound",
            },
            "controls": {"independent_spin_transfer": controls},
            "reproductions": {
                "v24": {
                    "v": str(series26.coefficients[24]),
                    "coefficients": [
                        str(value) for value in series26.coefficients[: len(v24_reference)]
                    ],
                },
                "v26": {"v": str(series26.coefficients[26]), "coefficients": [str(value) for value in series26.coefficients]},
            },
            "resource_wall": {
                "shape": wall.record["shape"],
                "degree": wall.record["degree"],
                "profile": wall.record["profile"],
                "peak_states": wall.record["peak_states"],
                "terminal_states": wall.record["terminal_states"],
                "max_rss_bytes": wall.record["max_rss_bytes"],
                "wall_seconds": wall.record["wall_seconds"],
                "max_states": max_states,
                "state_density_ceiling_bytes": (
                    (wall.record["max_rss_bytes"] + wall.record["peak_states"] - 1)
                    // wall.record["peak_states"]
                ),
                "projected_bytes_at_peak_states": (
                    ((wall.record["max_rss_bytes"] + wall.record["peak_states"] - 1)
                    // wall.record["peak_states"])
                    * wall.record["peak_states"]
                ),
            },
            "production_resource_measurement": {
                "command": ".venv/bin/python experiments/e94_ht_v28.py",
                "maximum_resident_set_bytes": _rss_bytes(),
                "real_seconds": f"{time.perf_counter() - started:.6f}",
                "disk_checkpoint_bytes": 0,
            },
        }
        payload = _make_payload(data, checks)
        if write_artifacts:
            RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            PROOF_PATH.write_text(_proof_text(data), encoding="utf-8")
        return payload

    final_injections = {
        (5, 5, 5): {24: w555_24, 26: w555_26, 28: w555_28},
        (5, 5, 6): {26: w556_26, 28: w556_28},
        (4, 6, 6): {26: w466_26, 28: w466_28},
        (4, 6, 7): {28: w467_28},
        (5, 5, 7): {28: w557_28},
        (5, 6, 6): {28: w566_28},
    }
    series28, transfer_records = _assemble(28, final_injections)

    primes_ok = (
        len(set(CRT_PRIMES)) == len(CRT_PRIMES)
        and all(_is_prime_32(prime) for prime in CRT_PRIMES)
        and all(
            gcd(CRT_PRIMES[left], CRT_PRIMES[right]) == 1
            for left in range(len(CRT_PRIMES))
            for right in range(left)
        )
    )
    crt_ok = all(
        CRT_PRODUCT > 2 * int(row["bound"])
        and abs(int(row["centered"])) <= int(row["bound"])
        and int(row["residue"]) == int(row["centered"]) % CRT_PRODUCT
        for row in profile_records
    )
    _record(
        checks,
        "profile_wise_crt_centered_reconstruction",
        primes_ok and crt_ok,
        f"M={CRT_PRODUCT} exceeds twice every exact labelled-slot profile bound",
    )
    primitive_profiles_ok = all(
        row.get("primitive_profile") is True
        and gcd(*tuple(row["profile"])) == 1
        for row in profile_records
    )
    _record(
        checks,
        "primitive_log_profile_integrality",
        primitive_profiles_ok,
        "every extracted profile is primitive, so its formal logarithm coefficient is integral",
    )
    _record(
        checks,
        "order28_cube_orbit_inventory",
        len(profile_orbits((5, 5, 5), 28)) == 9
        and sum(item["multiplicity"] for item in profile_orbits((5, 5, 5), 28)) == 78,
        "two six-crossing and seven distance-sensitive pair-cut cube orbits cover 78 profiles",
    )
    _record(
        checks,
        "order28_rectangular_orbit_inventory",
        len(profile_orbits((4, 6, 6), 28)) == 5
        and len(profile_orbits((5, 5, 6), 28)) == 5,
        "each span-13 rectangle has five reflection/equal-axis one-cap-four orbits",
    )
    _record(
        checks,
        "v28_prefix_and_odd_parity",
        series28.coefficients[: len(v26_reference)] == v26_reference
        and series28.coefficients[27] == 0
        and series28.interaction_coefficients[27] == 0,
        "fresh order-28 assembly reproduces v0..v26 and has zero odd v27 coefficient",
    )
    external_match = series28.interaction_coefficients[28] == EXTERNAL_A28
    _record(
        checks,
        "external_arisue_fujiwara_v28_witness",
        external_match,
        f"computed a28={series28.interaction_coefficients[28]}; external witness={EXTERNAL_A28}",
    )
    _record(
        checks,
        "v28_normalization",
        series28.coefficients[28]
        == series28.interaction_coefficients[28] + Fraction(3, 28),
        f"[v^28] phi={series28.coefficients[28]} = a28 + 3/28",
    )

    elapsed = time.perf_counter() - started
    run_usage = resource.getrusage(resource.RUSAGE_SELF)
    data = {
        "status": "[COMPUTATION]",
        "claim_tag": "[COMPUTATION]",
        "scope": "finite exact simple-cubic high-temperature series coefficient; not a solution of the 3D Ising model",
        "box_class_lemma": box_data,
        "algorithm": {
            "name": "box-general parity frontier with per-cut cap vector and ordinary mixed-radix logarithm",
            "pending_register": "b*c bit parity register for an a-by-b-by-c sweep",
            "profile_caps": "profile alpha uses terminal cut cap 2*alpha and exponent radix alpha+1",
            "reason_exact": "every degree edge crosses one cut, positive even cut totals determine every requested profile",
            "primitive_log_profiles": "every requested alpha has gcd(alpha)=1; the formal-log coefficient is therefore integral",
            "rss_guard": (
                "workers sample RSS while source and destination frontiers coexist; "
                f"they stop with [UNRESOLVED] at {MEMORY_GUARD_TRIP_BYTES} bytes"
            ),
            "disk_checkpoint_bytes": 0,
        },
        "profiles": profile_records,
        "frontier": {
            "summary": _profile_summary(profile_records),
            "control_profile_count": len(control_rows),
            "memory_budget_bytes": MEMORY_BUDGET_BYTES,
        },
        "crt": {
            "primes": list(CRT_PRIMES),
            "pairwise_coprime": primes_ok,
            "deterministic_primality_checked": primes_ok,
            "modulus_product": str(CRT_PRODUCT),
            "modulus_product_bits": CRT_PRODUCT.bit_length(),
            "centered_reconstruction_condition": "modulus_product > 2*absolute_profile_bound",
            "ordered_partition_factors": {
                "O_12": str(ordered_partition_factor(12)),
                "O_13": str(ordered_partition_factor(13)),
                "O_14": str(ordered_partition_factor(14)),
            },
            "bound_formula": "O_sum(alpha_i) * product_i C(E_i,2*alpha_i)^alpha_i",
        },
        "weights": {
            "w555_24": str(w555_24),
            "w555_26": str(w555_26),
            "w555_28": str(w555_28),
            "w556_26": str(w556_26),
            "w556_28": str(w556_28),
            "w466_26": str(w466_26),
            "w466_28": str(w466_28),
            "w467_28": str(w467_28),
            "w557_28": str(w557_28),
            "w566_28": str(w566_28),
            "injections": {
                str(list(shape)): {str(degree): str(value) for degree, value in values.items()}
                for shape, values in final_injections.items()
            },
        },
        "reproductions": {
            "v24": {
                "v": str(series26.coefficients[24]),
                "interaction": str(series26.interaction_coefficients[24]),
                "coefficients": [
                    str(value) for value in series26.coefficients[: len(v24_reference)]
                ],
            },
            "v26": {
                "v": str(series26.coefficients[26]),
                "interaction": str(series26.interaction_coefficients[26]),
                "coefficients": [str(value) for value in series26.coefficients],
            },
        },
        "series": {
            "dimension": 3,
            "variable": "v",
            "achieved_order": ORDER,
            "normalization": "phi = log(2) + sum_n coefficients[n] v^n",
            "coefficients": [str(value) for value in series28.coefficients],
            "interaction_coefficients": [str(value) for value in series28.interaction_coefficients],
            "v28": str(series28.coefficients[28]),
            "interaction_v28": str(series28.interaction_coefficients[28]),
            "v24_artifact_sha256": hashlib.sha256(V24_PATH.read_bytes()).hexdigest(),
            "v26_artifact_sha256": hashlib.sha256(V26_PATH.read_bytes()).hexdigest(),
            "transfer_record_count": len(transfer_records),
            "transfer_records_large_cross_section": [
                record
                for record in transfer_records.values()
                if record["cross_section_sites"] >= 16
            ],
        },
        "controls": {
            "independent_spin_transfer": controls,
            "external_arisue_fujiwara": {
                "tag": "[EXTERNAL]",
                **EXTERNAL_TABLE,
                "agreement": external_match,
                "relation": "checked only after exact finite-lattice assembly",
            },
        },
        "production_resource_measurement": {
            "command": ".venv/bin/python experiments/e94_ht_v28.py",
            "measurement_method": "resource.getrusage(RUSAGE_SELF) in producer; isolated worker RSS in every profile record",
            "maximum_resident_set_bytes": _rss_bytes(),
            "rss_before_bytes": rss_before,
            "real_seconds": f"{elapsed:.6f}",
            "disk_input_blocks_before_artifact_write": int(run_usage.ru_inblock),
            "disk_output_blocks_before_artifact_write": int(run_usage.ru_oublock),
            "disk_checkpoint_bytes": 0,
        },
    }
    payload = _make_payload(data, checks)
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more embedded checks failed")
    if write_artifacts:
        RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        PROOF_PATH.write_text(_proof_text(data), encoding="utf-8")
    return payload


def _worker_main(args: argparse.Namespace) -> None:
    shape = tuple(int(side) for side in args.shape.split(","))
    profile = tuple(int(value) for value in args.half_profile.split(","))
    modulus = CRT_PRIMES[0] if args.one_prime else CRT_PRODUCT
    payload = extract_frontier_profile(
        shape,
        profile,
        modulus=modulus,
        max_states=args.max_states,
        memory_budget_bytes=args.memory_budget_bytes,
    )
    payload["max_rss_bytes"] = _rss_bytes()
    print(json.dumps(payload, separators=(",", ":")))


def _profile_main(args: argparse.Namespace) -> None:
    shape = tuple(int(side) for side in args.shape.split(","))
    profile = tuple(int(value) for value in args.half_profile.split(","))
    modulus = CRT_PRIMES[0] if args.one_prime else CRT_PRODUCT
    payload = extract_frontier_profile(
        shape,
        profile,
        modulus=modulus,
        max_states=args.max_states,
        memory_budget_bytes=args.memory_budget_bytes,
    )
    print(
        f"shape={payload['shape']} profile={payload['profile']} completed={payload['completed']} "
        f"peak_states={payload['peak_states']} terminal_states={payload['terminal_states']} "
        f"wall_seconds={payload['wall_seconds']}"
    )
    if payload["completed"]:
        print(
            f"centered={payload['centered']} bound={payload['bound']} "
            f"certified={modulus > 2 * int(payload['bound'])}"
        )
    else:
        print(
            f"wall_reason={payload['wall_reason']} rss_bytes={payload['rss_bytes']} "
            f"memory_budget_bytes={payload['memory_budget_bytes']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--profile-run", action="store_true")
    parser.add_argument("--shape", type=str)
    parser.add_argument("--half-profile", type=str)
    parser.add_argument("--one-prime", action="store_true")
    parser.add_argument("--max-states", type=int)
    parser.add_argument(
        "--memory-budget-bytes", type=int, default=MEMORY_BUDGET_BYTES
    )
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    if args.worker:
        if args.shape is None or args.half_profile is None:
            parser.error("--worker requires --shape and --half-profile")
        _worker_main(args)
        return
    if args.profile_run:
        if args.shape is None or args.half_profile is None:
            parser.error("--profile-run requires --shape and --half-profile")
        _profile_main(args)
        return
    run_full_experiment(
        write_artifacts=not args.no_write,
        max_states=args.max_states,
        memory_budget_bytes=args.memory_budget_bytes,
    )
    print("PASS")


if __name__ == "__main__":
    main()
