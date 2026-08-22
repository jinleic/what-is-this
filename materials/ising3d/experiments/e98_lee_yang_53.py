"""Exact streamed Lee--Yang circle evaluations with cross-section symmetries.

The graded 5x5 layer transfer from e88 needs a 2**25 by 63 residue buffer.
This experiment instead evaluates the reduced circle polynomial Q(t) at one
rational t at a time.  It works in F_p[z]/(z**2 - 2*t*z + 1), stores one
algebra element per spatial-and-spin orbit, and streams the factorized
interlayer butterflies through bounded state-pair strips.

The default run does not claim a 5x5x5 value.  It validates the new evaluator
at 4x4x4 and at the axis-permuted 5x4x4 box, measures the 4x4x4 strip wall,
and records a transparent 5x5x5 extrapolation.

Run from the repository root:
    timeout 7200 .venv/bin/python experiments/e98_lee_yang_53.py

All decisive arithmetic is integer modular arithmetic plus CRT and Fraction.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

from ising.lee_yang import _CRT_PRIMES, rational_field_polynomial_transfer
from ising.transfer_matrix import layer_bonds

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e98_lee_yang_53.py"
OUTPUT = ROOT / "results" / "lee_yang" / "zeros_53.json"
PRIOR_ARTIFACT = ROOT / "results" / "lee_yang" / "zeros_l5.json"

COUPLING = (2, 3)
VALIDATION_T_VALUES = (Fraction(0), Fraction(1, 2))
DEFAULT_STRIP_PAIRS = 1 << 20
VALIDATED_STATE_LIMIT = 1 << 20


def fraction_text(value: Fraction) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_box_bond_count(shape: tuple[int, int, int]) -> int:
    a, b, c = shape
    return (a - 1) * b * c + a * (b - 1) * c + a * b * (c - 1)


# ---------------------------------------------------------------------------
# Cross-section symmetry bookkeeping.
#
# A group element is a bit permutation: permutation[source_bit] is its target
# bit.  Rectangles have row/column reflections; squares additionally have the
# two diagonal reflections and quarter turns (D8).
# ---------------------------------------------------------------------------


def cross_permutations(cross: tuple[int, int]) -> tuple[tuple[int, ...], ...]:
    rows, cols = cross

    def bit(row: int, col: int) -> int:
        return row * cols + col

    if rows == cols:
        side = rows
        transforms = (
            lambda r, c: (r, c),
            lambda r, c: (c, side - 1 - r),
            lambda r, c: (side - 1 - r, side - 1 - c),
            lambda r, c: (side - 1 - c, r),
            lambda r, c: (side - 1 - r, c),
            lambda r, c: (r, side - 1 - c),
            lambda r, c: (c, r),
            lambda r, c: (side - 1 - c, side - 1 - r),
        )
    else:
        transforms = (
            lambda r, c: (r, c),
            lambda r, c: (rows - 1 - r, c),
            lambda r, c: (r, cols - 1 - c),
            lambda r, c: (rows - 1 - r, cols - 1 - c),
        )

    permutations: list[tuple[int, ...]] = []
    for transform in transforms:
        permutation = tuple(
            bit(*transform(row, col))
            for row in range(rows)
            for col in range(cols)
        )
        if permutation not in permutations:
            permutations.append(permutation)
    return tuple(permutations)


def permutation_cycle_count(permutation: tuple[int, ...]) -> int:
    seen = [False] * len(permutation)
    cycles = 0
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cycles += 1
        point = start
        while not seen[point]:
            seen[point] = True
            point = permutation[point]
    return cycles


def permutation_has_only_even_cycles(permutation: tuple[int, ...]) -> bool:
    seen = [False] * len(permutation)
    for start in range(len(permutation)):
        if seen[start]:
            continue
        point = start
        length = 0
        while not seen[point]:
            seen[point] = True
            point = permutation[point]
            length += 1
        if length % 2:
            return False
    return True


def burnside_orbit_counts(cross: tuple[int, int]) -> dict[str, int]:
    """Exact orbit counts for spatial and spatial-plus-spin actions."""

    permutations = cross_permutations(cross)
    spatial_fixed_sum = sum(1 << permutation_cycle_count(item) for item in permutations)
    complement_fixed_sum = sum(
        (1 << permutation_cycle_count(item))
        for item in permutations
        if permutation_has_only_even_cycles(item)
    )
    group_order = len(permutations)
    if spatial_fixed_sum % group_order or (spatial_fixed_sum + complement_fixed_sum) % (2 * group_order):
        raise ArithmeticError("Burnside orbit count was not integral")
    return {
        "state_count": 1 << (cross[0] * cross[1]),
        "spatial_group_order": group_order,
        "spatial_orbits": spatial_fixed_sum // group_order,
        "spatial_plus_spin_orbits": (spatial_fixed_sum + complement_fixed_sum) // (2 * group_order),
        "spatial_fixed_sum": spatial_fixed_sum,
        "spin_complement_fixed_sum": complement_fixed_sum,
    }


def permute_state_vector(states: np.ndarray, permutation: tuple[int, ...]) -> np.ndarray:
    """Apply a bit permutation to a uint32 vector of cross-section states."""

    out = np.zeros_like(states)
    one = np.uint32(1)
    for source, target in enumerate(permutation):
        out |= ((states >> np.uint32(source)) & one) << np.uint32(target)
    return out


def state_statistics_for_states(
    states: np.ndarray, cross: tuple[int, int]
) -> tuple[np.ndarray, np.ndarray]:
    """Broken in-layer bonds and down-spin counts for selected state integers."""

    ns = cross[0] * cross[1]
    one = np.uint32(1)
    occupied = np.zeros(states.shape, dtype=np.uint8)
    for bit in range(ns):
        occupied += ((states >> np.uint32(bit)) & one).astype(np.uint8)
    down = np.uint8(ns) - occupied
    broken = np.zeros(states.shape, dtype=np.uint8)
    for left, right in layer_bonds(cross, (False, False)):
        broken += (
            ((states >> np.uint32(left)) ^ (states >> np.uint32(right))) & one
        ).astype(np.uint8)
    return broken, down


def build_symmetry_data(
    cross: tuple[int, int], *, allow_large: bool = False
) -> dict:
    """Build exact spatial-plus-spin orbit maps for one cross-section.

    The stored representative is the least spatial representative of a
    spatial orbit chosen from each spin-complement pair.  ``flipped`` records
    whether a raw state is in the complementary spatial orbit.  At layer j,
    its value is z**(j*ns) times the z -> z**-1 conjugate of the stored value.
    """

    ns = cross[0] * cross[1]
    n_states = 1 << ns
    if n_states > VALIDATED_STATE_LIMIT and not allow_large:
        raise ValueError(
            f"{cross} has {n_states:,} states; pass allow_large=True only after a budget check"
        )

    permutations = cross_permutations(cross)
    states = np.arange(n_states, dtype=np.uint32)
    canonical = states.copy()
    for permutation in permutations[1:]:
        image = permute_state_vector(states, permutation)
        np.minimum(canonical, image, out=canonical)

    spatial_representatives = np.unique(canonical)
    sentinel = np.uint32(np.iinfo(np.uint32).max)
    representative_to_spatial = np.full(n_states, sentinel, dtype=np.uint32)
    representative_to_spatial[spatial_representatives] = np.arange(
        len(spatial_representatives), dtype=np.uint32
    )
    spatial_for_state = representative_to_spatial[canonical]
    if np.any(spatial_for_state == sentinel):
        raise ArithmeticError("spatial canonical map was incomplete")

    mask = np.uint32(n_states - 1)
    complement_states = mask ^ spatial_representatives
    complement_spatial = representative_to_spatial[canonical[complement_states]]
    spatial_indices = np.arange(len(spatial_representatives), dtype=np.uint32)
    leaders = spatial_indices <= complement_spatial
    selected_spatial = spatial_indices[leaders]

    spatial_to_orbit = np.full(len(spatial_representatives), sentinel, dtype=np.uint32)
    spatial_to_orbit[selected_spatial] = np.arange(len(selected_spatial), dtype=np.uint32)
    smaller_partner = np.minimum(spatial_indices, complement_spatial)
    orbit_for_spatial = spatial_to_orbit[smaller_partner]
    if np.any(orbit_for_spatial == sentinel):
        raise ArithmeticError("spin-paired canonical map was incomplete")
    flipped_for_spatial = spatial_indices > complement_spatial
    state_to_orbit = orbit_for_spatial[spatial_for_state]
    flipped = flipped_for_spatial[spatial_for_state]
    representatives = spatial_representatives[selected_spatial]

    count0 = np.bincount(
        state_to_orbit[~flipped], minlength=len(representatives)
    ).astype(np.uint64)
    count1 = np.bincount(
        state_to_orbit[flipped], minlength=len(representatives)
    ).astype(np.uint64)
    broken, down = state_statistics_for_states(representatives, cross)

    formula = burnside_orbit_counts(cross)
    if len(spatial_representatives) != formula["spatial_orbits"]:
        raise ArithmeticError("enumerated spatial orbit count disagrees with Burnside")
    if len(representatives) != formula["spatial_plus_spin_orbits"]:
        raise ArithmeticError("enumerated spin-paired orbit count disagrees with Burnside")
    if int(count0.sum() + count1.sum()) != n_states:
        raise ArithmeticError("orbit multiplicities did not cover every state")

    # Arrays that existed only to construct the maps can be discarded before
    # the streamed transfer allocates its working vectors.
    del states, canonical, representative_to_spatial, spatial_for_state
    del complement_states, complement_spatial, spatial_indices, spatial_to_orbit
    del smaller_partner, orbit_for_spatial, flipped_for_spatial
    gc.collect()

    return {
        "cross": cross,
        "ns": ns,
        "n_states": n_states,
        "permutations": permutations,
        "representatives": representatives,
        "state_to_orbit": state_to_orbit,
        "flipped": flipped,
        "flip_indices": np.flatnonzero(flipped).astype(np.uint32),
        "count0": count0,
        "count1": count1,
        "broken": broken,
        "down": down,
        "formula": formula,
    }


# ---------------------------------------------------------------------------
# Exact quadratic-algebra arithmetic and striped butterflies.
# ---------------------------------------------------------------------------


def pair_multiply(
    left_a: int, left_b: int, right_a: int, right_b: int, two_t: int, modulus: int
) -> tuple[int, int]:
    """Multiply (a+bz) pairs modulo z^2-2tz+1 and ``modulus``."""

    return (
        (left_a * right_a - left_b * right_b) % modulus,
        (left_a * right_b + left_b * right_a + two_t * left_b * right_b) % modulus,
    )


def pair_conjugate(value_a: int, value_b: int, two_t: int, modulus: int) -> tuple[int, int]:
    """Apply z -> z**-1 = 2t-z to one algebra element."""

    return ((value_a + two_t * value_b) % modulus, (-value_b) % modulus)


def z_power_pair(power: int, two_t: int, modulus: int) -> tuple[int, int]:
    """Return z**power in the quadratic algebra by an exact recurrence."""

    value_a, value_b = 1, 0
    for _ in range(power):
        value_a, value_b = (-value_b) % modulus, (value_a + two_t * value_b) % modulus
    return value_a, value_b


def butterfly_blocks(ns: int, bit: int, strip_pairs: int):
    """Yield bounded (hi, lo) rectangles for one tensor-product butterfly."""

    if strip_pairs < 1:
        raise ValueError("strip_pairs must be positive")
    lo = 1 << bit
    hi = 1 << (ns - bit - 1)
    if lo <= strip_pairs:
        rows = max(1, strip_pairs // lo)
        for first_hi in range(0, hi, rows):
            yield first_hi, min(hi, first_hi + rows), 0, lo
    else:
        for first_hi in range(hi):
            for first_lo in range(0, lo, strip_pairs):
                yield first_hi, first_hi + 1, first_lo, min(lo, first_lo + strip_pairs)


def butterfly_strip_count(ns: int, strip_pairs: int) -> int:
    return sum(1 for bit in range(ns) for _ in butterfly_blocks(ns, bit, strip_pairs))


def apply_interlayer_component(
    source: np.ndarray,
    scratch: np.ndarray,
    ns: int,
    agree: np.uint64,
    disagree: np.uint64,
    modulus: np.uint64,
    strip_pairs: int,
    strip_order: str,
) -> tuple[np.ndarray, int]:
    """Apply all interlayer butterflies to one algebra component exactly."""

    if strip_order not in ("forward", "reverse"):
        raise ValueError("strip_order must be 'forward' or 'reverse'")
    current, updated = source, scratch
    strips = 0
    for bit in range(ns):
        lo = 1 << bit
        hi = 1 << (ns - bit - 1)
        current_view = current.reshape(hi, 2, lo)
        updated_view = updated.reshape(hi, 2, lo)
        blocks = butterfly_blocks(ns, bit, strip_pairs)
        if strip_order == "reverse":
            blocks = reversed(tuple(blocks))
        for first_hi, last_hi, first_lo, last_lo in blocks:
            lower = current_view[first_hi:last_hi, 0, first_lo:last_lo]
            upper = current_view[first_hi:last_hi, 1, first_lo:last_lo]
            updated_view[first_hi:last_hi, 0, first_lo:last_lo] = (
                agree * lower + disagree * upper
            ) % modulus
            updated_view[first_hi:last_hi, 1, first_lo:last_lo] = (
                disagree * lower + agree * upper
            ) % modulus
            strips += 1
        current, updated = updated, current
    return current, strips


def expand_spin_paired_values(
    current_a: np.ndarray,
    current_b: np.ndarray,
    symmetry: dict,
    layer_count: int,
    two_t: int,
    modulus: int,
    raw_a: np.ndarray,
    raw_b: np.ndarray,
) -> None:
    """Expand one stored orbit vector to raw states without materializing grades."""

    np.take(current_a, symmetry["state_to_orbit"], out=raw_a)
    np.take(current_b, symmetry["state_to_orbit"], out=raw_b)
    power_a, power_b = z_power_pair(layer_count * symmetry["ns"], two_t, modulus)
    conjugate_b_coefficient = (two_t * power_a + power_b) % modulus
    negative_power_a = (-power_a) % modulus
    flip_indices = symmetry["flip_indices"]
    chunk = 1 << 20
    for start in range(0, len(flip_indices), chunk):
        indices = flip_indices[start : start + chunk]
        old_a = raw_a[indices]
        old_b = raw_b[indices]
        raw_a[indices] = (power_a * old_a + conjugate_b_coefficient * old_b) % modulus
        raw_b[indices] = (power_b * old_a + negative_power_a * old_b) % modulus


def field_factors_for_representatives(
    symmetry: dict, numerator: int, denominator: int, two_t: int, modulus: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return w_s z**down(s) for every stored representative."""

    ns = symmetry["ns"]
    powers = [z_power_pair(power, two_t, modulus) for power in range(ns + 1)]
    bonds_per_layer = len(layer_bonds(symmetry["cross"], (False, False)))
    first = np.empty(len(symmetry["representatives"]), dtype=np.uint64)
    second = np.empty_like(first)
    for index, (broken, down) in enumerate(zip(symmetry["broken"], symmetry["down"])):
        weight = (
            pow(numerator, int(broken), modulus)
            * pow(denominator, bonds_per_layer - int(broken), modulus)
        ) % modulus
        z_a, z_b = powers[int(down)]
        first[index] = (weight * z_a) % modulus
        second[index] = (weight * z_b) % modulus
    return first, second


def allocate_vector_pair(
    length: int, disk_dir: Path | None, label: str
) -> tuple[np.ndarray, np.ndarray, tuple[Path, Path] | None]:
    if disk_dir is None:
        return np.empty(length, dtype=np.uint64), np.empty(length, dtype=np.uint64), None
    disk_dir.mkdir(parents=True, exist_ok=True)
    first_path = disk_dir / f"{label}_a.u64"
    second_path = disk_dir / f"{label}_b.u64"
    first = np.memmap(first_path, mode="w+", dtype=np.uint64, shape=(length,))
    second = np.memmap(second_path, mode="w+", dtype=np.uint64, shape=(length,))
    return first, second, (first_path, second_path)


def release_vector_pair(
    first: np.ndarray, second: np.ndarray, paths: tuple[Path, Path] | None
) -> None:
    if isinstance(first, np.memmap):
        first.flush()
    if isinstance(second, np.memmap):
        second.flush()
    del first, second
    if paths is not None:
        for path in paths:
            path.unlink(missing_ok=True)


def modular_reduced_circle_residue(
    shape: tuple[int, int, int],
    numerator: int,
    denominator: int,
    t: Fraction,
    modulus: int,
    symmetry: dict,
    *,
    strip_pairs: int = DEFAULT_STRIP_PAIRS,
    strip_order: str = "forward",
    disk_dir: Path | None = None,
    label: str = "residue",
) -> dict:
    """One exact modular residue of Q(t), with streamed state-pair butterflies."""

    cross = (shape[0], shape[1])
    length = shape[2]
    ns = cross[0] * cross[1]
    n_sites = ns * length
    if symmetry["cross"] != cross:
        raise ValueError("symmetry data does not match transfer cross-section")
    if modulus <= 2 or t.denominator % modulus == 0:
        raise ValueError("modulus is not admissible for t")

    t_mod = (t.numerator % modulus) * pow(t.denominator % modulus, -1, modulus) % modulus
    two_t = (2 * t_mod) % modulus
    if (two_t + 2) % modulus == 0:
        raise ValueError("t == -1 modulo modulus makes division by z+1 singular")

    factor_a, factor_b = field_factors_for_representatives(
        symmetry, numerator, denominator, two_t, modulus
    )
    field_b_coefficient = (factor_a + two_t * factor_b) % modulus
    negative_factor_b = (np.uint64(modulus) - factor_b) % np.uint64(modulus)
    n_orbits = len(symmetry["representatives"])

    current_a, current_b, current_paths = allocate_vector_pair(
        n_orbits, disk_dir, f"{label}_layer1"
    )
    current_a[:] = factor_a
    current_b[:] = factor_b
    raw = np.empty(symmetry["n_states"], dtype=np.uint64)
    scratch = np.empty_like(raw)
    transformed_a = np.empty(n_orbits, dtype=np.uint64)
    transformed_b = np.empty_like(transformed_a)
    agree = np.uint64(denominator % modulus)
    disagree = np.uint64(numerator % modulus)
    modulus_u = np.uint64(modulus)
    strip_total = 0

    try:
        for layer in range(1, length):
            expand_spin_paired_values(
                current_a,
                current_b,
                symmetry,
                layer,
                two_t,
                modulus,
                raw,
                scratch,
            )
            transformed, strips = apply_interlayer_component(
                raw,
                scratch,
                ns,
                agree,
                disagree,
                modulus_u,
                strip_pairs,
                strip_order,
            )
            transformed_a[:] = transformed[symmetry["representatives"]]
            strip_total += strips

            expand_spin_paired_values(
                current_a,
                current_b,
                symmetry,
                layer,
                two_t,
                modulus,
                raw,
                scratch,
            )
            # ``scratch`` is the expanded b component.  Copy it into the
            # source buffer; the butterfly immediately overwrites scratch as
            # its output buffer, so only two raw vectors are live.
            raw[:] = scratch
            transformed, strips = apply_interlayer_component(
                raw,
                scratch,
                ns,
                agree,
                disagree,
                modulus_u,
                strip_pairs,
                strip_order,
            )
            transformed_b[:] = transformed[symmetry["representatives"]]
            strip_total += strips

            next_a, next_b, next_paths = allocate_vector_pair(
                n_orbits, disk_dir, f"{label}_layer{layer + 1}"
            )
            next_a[:] = (factor_a * transformed_a + negative_factor_b * transformed_b) % modulus
            next_b[:] = (factor_b * transformed_a + field_b_coefficient * transformed_b) % modulus
            if isinstance(next_a, np.memmap):
                next_a.flush()
                next_b.flush()
            release_vector_pair(current_a, current_b, current_paths)
            current_a, current_b, current_paths = next_a, next_b, next_paths

        final_power_a, final_power_b = z_power_pair(n_sites, two_t, modulus)
        conjugate_b_coefficient = (two_t * final_power_a + final_power_b) % modulus
        negative_final_power_a = (-final_power_a) % modulus
        flipped_a = (
            final_power_a * current_a + conjugate_b_coefficient * current_b
        ) % modulus
        flipped_b = (final_power_b * current_a + negative_final_power_a * current_b) % modulus
        total_a = np.sum(
            (symmetry["count0"] * current_a + symmetry["count1"] * flipped_a) % modulus,
            dtype=np.uint64,
        ) % modulus_u
        total_b = np.sum(
            (symmetry["count0"] * current_b + symmetry["count1"] * flipped_b) % modulus,
            dtype=np.uint64,
        ) % modulus_u
        polynomial_a, polynomial_b = int(total_a), int(total_b)

        if n_sites % 2:
            inverse = pow((two_t + 2) % modulus, -1, modulus)
            divisor_a = ((two_t + 1) * inverse) % modulus
            divisor_b = (-inverse) % modulus
            polynomial_a, polynomial_b = pair_multiply(
                polynomial_a,
                polynomial_b,
                divisor_a,
                divisor_b,
                two_t,
                modulus,
            )
        center = n_sites // 2
        center_a, center_b = z_power_pair(center, two_t, modulus)
        inverse_center_a, inverse_center_b = pair_conjugate(
            center_a, center_b, two_t, modulus
        )
        reduced_a, reduced_b = pair_multiply(
            polynomial_a,
            polynomial_b,
            inverse_center_a,
            inverse_center_b,
            two_t,
            modulus,
        )
        if reduced_b % modulus:
            raise ArithmeticError(
                "spin-paired evaluation did not reduce to the base field: "
                f"P=({polynomial_a},{polynomial_b}), "
                f"z^-m=({inverse_center_a},{inverse_center_b}), "
                f"Q=({reduced_a},{reduced_b}) mod {modulus}"
            )
        return {
            "residue": reduced_a,
            "butterfly_strips": strip_total,
            "n_orbits": n_orbits,
        }
    finally:
        release_vector_pair(current_a, current_b, current_paths)


# ---------------------------------------------------------------------------
# CRT reconstruction and independent reference evaluation.
# ---------------------------------------------------------------------------


def circle_value_bound(
    shape: tuple[int, int, int], numerator: int, denominator: int, t: Fraction
) -> int:
    """A safe signed bound for den(t)**floor(N/2) Q(t)."""

    n_sites = math.prod(shape)
    n_bonds = open_box_bond_count(shape)
    coefficient_sum_bound = (1 << n_sites) * max(numerator, denominator) ** n_bonds
    # Dividing an odd palindrome by z+1 gives alternating partial sums, each
    # bounded by the coefficient sum.  The Chebyshev expression has at most
    # N coefficients and |T_j(t)| <= 1 for rational t in [-1, 1].
    return (n_sites + 1) * coefficient_sum_bound * t.denominator ** (n_sites // 2)


def exact_reduced_circle_value(
    shape: tuple[int, int, int],
    numerator: int,
    denominator: int,
    t: Fraction,
    *,
    strip_pairs: int = DEFAULT_STRIP_PAIRS,
    strip_order: str = "forward",
    allow_large: bool = False,
    disk_dir: Path | None = None,
) -> dict:
    """CRT-reconstruct the exact rational reduced-circle value Q(t)."""

    shape = tuple(map(int, shape))
    if len(shape) != 3 or any(side < 1 for side in shape):
        raise ValueError("shape must be a positive three-dimensional box")
    t = Fraction(t)
    if not -1 < t < 1:
        raise ValueError("the circle evaluator is restricted to rational -1 < t < 1")
    symmetry = build_symmetry_data((shape[0], shape[1]), allow_large=allow_large)
    bound = circle_value_bound(shape, numerator, denominator, t)
    scaled_value = 0
    combined = 1
    modular_records: list[dict] = []
    started = time.perf_counter()
    for prime_index, modulus in enumerate(_CRT_PRIMES, start=1):
        prime_started = time.perf_counter()
        record = modular_reduced_circle_residue(
            shape,
            numerator,
            denominator,
            t,
            modulus,
            symmetry,
            strip_pairs=strip_pairs,
            strip_order=strip_order,
            disk_dir=disk_dir,
            label=f"p{prime_index}",
        )
        scaled_residue = record["residue"] * pow(t.denominator, shape[0] * shape[1] * shape[2] // 2, modulus) % modulus
        inverse = pow(combined % modulus, -1, modulus)
        scaled_value += combined * ((scaled_residue - scaled_value % modulus) * inverse % modulus)
        combined *= modulus
        record["prime"] = modulus
        record["wall_seconds"] = time.perf_counter() - prime_started
        modular_records.append(record)
        if combined > 2 * bound:
            break
    else:
        raise ArithmeticError("checked-in CRT prime list did not exceed the signed value bound")

    signed = scaled_value if scaled_value <= combined // 2 else scaled_value - combined
    if abs(signed) > bound:
        raise ArithmeticError("centered CRT value exceeded its proved bound")
    value = Fraction(signed, t.denominator ** (math.prod(shape) // 2))
    return {
        "shape": list(shape),
        "x": fraction_text(Fraction(numerator, denominator)),
        "t": fraction_text(t),
        "value": value,
        "scaled_numerator": signed,
        "scaled_denominator": t.denominator ** (math.prod(shape) // 2),
        "value_bound_bits": bound.bit_length(),
        "combined_modulus_bits": combined.bit_length(),
        "crt_primes_used": len(modular_records),
        "primes": [item["prime"] for item in modular_records],
        "butterfly_strips": sum(item["butterfly_strips"] for item in modular_records),
        "wall_seconds": time.perf_counter() - started,
        "per_prime_wall_seconds": [item["wall_seconds"] for item in modular_records],
        "symmetry": symmetry["formula"],
        "strip_pairs": strip_pairs,
    }


def reduced_circle_value_from_coefficients(coefficients: list[int], t: Fraction) -> Fraction:
    """Independent exact Q(t) evaluation from an ascending field polynomial."""

    values = [int(value) for value in coefficients]
    degree = len(values) - 1
    if values != values[::-1]:
        raise ValueError("field polynomial is not palindromic")
    if degree % 2:
        quotient = [0] * degree
        quotient[0] = values[0]
        for index in range(1, degree):
            quotient[index] = values[index] - quotient[index - 1]
        if values[-1] != quotient[-1]:
            raise ArithmeticError("z+1 did not divide an odd palindrome")
        values = quotient
        degree -= 1
    center = degree // 2
    result = Fraction(values[center])
    chebyshev_previous = Fraction(1)
    chebyshev_current = t
    for power in range(1, center + 1):
        if power == 1:
            chebyshev = chebyshev_current
        else:
            chebyshev = 2 * t * chebyshev_current - chebyshev_previous
            chebyshev_previous, chebyshev_current = chebyshev_current, chebyshev
        result += 2 * values[center + power] * chebyshev
    return result


def prior_5x4x4_coefficients() -> list[int]:
    artifact = json.loads(PRIOR_ARTIFACT.read_text())
    return [
        int(value)
        for value in artifact["data"]["exact_zero_data"]["x_2_over_3"]["ladder16"]["4x4x5"]["A_k"]
    ]


def compact_value_record(record: dict) -> dict:
    return {
        "shape": record["shape"],
        "x": record["x"],
        "t": record["t"],
        "Q_t": fraction_text(record["value"]),
        "scaled_numerator": str(record["scaled_numerator"]),
        "scaled_denominator": str(record["scaled_denominator"]),
        "value_bound_bits": record["value_bound_bits"],
        "combined_modulus_bits": record["combined_modulus_bits"],
        "crt_primes_used": record["crt_primes_used"],
        "crt_primes": record["primes"],
        "butterfly_strips": record["butterfly_strips"],
        "wall_seconds": record["wall_seconds"],
        "symmetry": record["symmetry"],
        "strip_pairs": record["strip_pairs"],
    }


def validate_small_boxes(strip_pairs: int) -> tuple[dict, list[dict]]:
    """Reproduce independent exact Q(t) values for 4^3 and 5x4x4."""

    numerator, denominator = COUPLING
    reference_4cube = rational_field_polynomial_transfer((4, 4, 4), numerator, denominator)
    reference_5x4x4 = prior_5x4x4_coefficients()
    cases = (
        ("cube_4x4x4", (4, 4, 4), reference_4cube),
        ("box_5x4x4", (5, 4, 4), reference_5x4x4),
    )
    checks: list[dict] = []
    output: dict[str, dict] = {}
    for label, shape, coefficients in cases:
        output[label] = {}
        for t in VALIDATION_T_VALUES:
            record = exact_reduced_circle_value(
                shape, numerator, denominator, t, strip_pairs=strip_pairs
            )
            expected = reduced_circle_value_from_coefficients(coefficients, t)
            passed = record["value"] == expected
            checks.append(
                {
                    "name": f"{label}_Q_at_{fraction_text(t)}",
                    "passed": passed,
                    "detail": "streamed orbit evaluator equals independent full-coefficient evaluation",
                }
            )
            if not passed:
                raise AssertionError(f"{label} mismatch at t={t}: {record['value']} != {expected}")
            output[label][f"t_{t.numerator}_{t.denominator}"] = compact_value_record(record)
    return output, checks


def projected_53_budget(validation_4cube_t0: dict, strip_pairs: int) -> dict:
    """Scale only measured fixed-evaluation butterfly work; do not run 5^3."""

    source_shape = (4, 4, 4)
    target_shape = (5, 5, 5)
    source_work = 2 * (source_shape[2] - 1) * (source_shape[0] * source_shape[1]) * (1 << 16)
    target_work = 2 * (target_shape[2] - 1) * (target_shape[0] * target_shape[1]) * (1 << 25)
    work_ratio = Fraction(target_work, source_work)
    target_counts = burnside_orbit_counts((5, 5))
    source_wall = validation_4cube_t0["wall_seconds"]
    source_primes = validation_4cube_t0["crt_primes_used"]
    source_per_prime = source_wall / source_primes
    target_per_prime = float(work_ratio) * source_per_prime
    target_t0_bound = circle_value_bound(target_shape, *COUPLING, Fraction(0))
    target_half_bound = circle_value_bound(target_shape, *COUPLING, Fraction(1, 2))
    target_t0_primes = next(
        index
        for index in range(1, len(_CRT_PRIMES) + 1)
        if math.prod(_CRT_PRIMES[:index]) > 2 * target_t0_bound
    )
    target_half_primes = next(
        index
        for index in range(1, len(_CRT_PRIMES) + 1)
        if math.prod(_CRT_PRIMES[:index]) > 2 * target_half_bound
    )
    interpolation_t = Fraction(1, 32)
    interpolation_bound = circle_value_bound(target_shape, *COUPLING, interpolation_t)
    interpolation_primes = next(
        index
        for index in range(1, len(_CRT_PRIMES) + 1)
        if math.prod(_CRT_PRIMES[:index]) > 2 * interpolation_bound
    )
    per_component_strips = butterfly_strip_count(25, strip_pairs)
    target_strips_per_prime = 2 * (target_shape[2] - 1) * per_component_strips
    source_strips_per_prime = validation_4cube_t0["butterfly_strips"] // source_primes
    raw_bytes = 2 * (1 << 25) * 8
    map_bytes = (1 << 25) * (np.dtype(np.uint32).itemsize + np.dtype(np.bool_).itemsize)
    orbit_count = target_counts["spatial_plus_spin_orbits"]
    orbit_component_bytes = orbit_count * 10 * np.dtype(np.uint64).itemsize
    orbit_metadata_bytes = (
        (1 << 24) * np.dtype(np.uint32).itemsize
        + orbit_count
        * (
            np.dtype(np.uint32).itemsize
            + 2 * np.dtype(np.uint64).itemsize
            + 2 * np.dtype(np.uint8).itemsize
        )
    )
    return {
        "claim_tag": "COMPUTATION",
        "statement": "Measured 4^3 streamed fixed-circle passes are linearly extrapolated by exact butterfly-slot count only; no 5^3 pass was run.",
        "coupling": "2/3",
        "circle_points": {
            "t_0": {
                "signed_bound_bits": target_t0_bound.bit_length(),
                "required_existing_30bit_primes": target_t0_primes,
                "linear_estimated_seconds": target_per_prime * target_t0_primes,
            },
            "t_1_over_2": {
                "signed_bound_bits": target_half_bound.bit_length(),
                "required_existing_30bit_primes": target_half_primes,
                "linear_estimated_seconds": target_per_prime * target_half_primes,
            },
        },
        "coefficient_reconstruction_via_evaluations": {
            "reduced_Q_degree": 62,
            "distinct_rational_values_required": 63,
            "example_grid": "j/32 for -31 <= j <= 31",
            "worst_grid_signed_bound_bits": interpolation_bound.bit_length(),
            "worst_grid_required_existing_30bit_primes": interpolation_primes,
            "linear_estimated_seconds_for_63_values": target_per_prime
            * interpolation_primes
            * 63,
            "scope": "Exact rational interpolation would reconstruct Q, but this is only a linear work extrapolation and excludes 5x5 symmetry-map construction, disk traffic, and subsequent Sturm work.",
        },
        "validated_source": {
            "shape": list(source_shape),
            "wall_seconds_all_primes": source_wall,
            "primes": source_primes,
            "wall_seconds_per_prime": source_per_prime,
            "butterfly_strips_per_prime": source_strips_per_prime,
        },
        "exact_work": {
            "source_scalar_butterfly_slots": source_work,
            "target_scalar_butterfly_slots": target_work,
            "target_over_source_ratio": f"{work_ratio.numerator}/{work_ratio.denominator}",
            "target_butterfly_strips_per_prime": target_strips_per_prime,
            "strip_pairs": strip_pairs,
        },
        "target_symmetry": target_counts,
        "streamed_memory_layout_bytes": {
            "two_raw_component_buffers": raw_bytes,
            "persistent_state_to_orbit_and_flip_maps": map_bytes,
            "ten_long_lived_orbit_component_vectors": orbit_component_bytes,
            "persistent_orbit_metadata": orbit_metadata_bytes,
            "explicit_subtotal_before_numpy_expression_temporaries_and_allocator_overhead": (
                raw_bytes + map_bytes + orbit_component_bytes + orbit_metadata_bytes
            ),
            "scope": "The ten component vectors are factor_a/factor_b, field_b_coefficient/negative_factor_b, current_a/current_b, transformed_a/transformed_b, and next_a/next_b (or final flipped_a/flipped_b).",
        },
        "scope": "This is a fixed-Q(t) evaluator, not a reconstructed coefficient vector or an exact Sturm root count. A first-zero claim needs additional exact zero-counting certification.",
    }


def correction_route_status() -> dict:
    return {
        "claim_tag": "UNRESOLVED",
        "statement": "No correction-amplitude majorant from 5x4x4/5x5x4 slabs to the isotropic 5^3 first zero is claimed.",
        "certificate": [
            "The stored slab fourth points remain compatible with both p=2 and p=4 after the edge is allowed to move (proofs/lee_yang_zeros.md, four-size analysis).",
            "The L=2,3,4 isotropic cube data have certified globally decreasing p=2 and p=4 continuations with a common positive edge and exponent separation 3/4 (proofs/lee_yang_prg.md).",
            "No repository theorem controls complex first Lee--Yang zeros under transverse layer addition, bond addition, or box inclusion; real-field correlation monotonicity is not such a bridge.",
        ],
        "conclusion": "The correction-bound alternative is not established; the artifact records only the streamed-transfer fallback.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strip-pairs", type=int, default=DEFAULT_STRIP_PAIRS)
    parser.add_argument("--evaluate-53-t", type=str)
    parser.add_argument("--disk-dir", type=Path)
    args = parser.parse_args()

    if args.evaluate_53_t:
        numerator_text, denominator_text = args.evaluate_53_t.split("/", 1)
        target_t = Fraction(int(numerator_text), int(denominator_text))
        record = exact_reduced_circle_value(
            (5, 5, 5),
            *COUPLING,
            target_t,
            strip_pairs=args.strip_pairs,
            allow_large=True,
            disk_dir=args.disk_dir,
        )
        print(json.dumps(compact_value_record(record), indent=2, sort_keys=True))
        print("PASS")
        return

    started = time.time()
    validation, checks = validate_small_boxes(args.strip_pairs)
    source_record = validation["cube_4x4x4"]["t_0_1"]
    budget = projected_53_budget(source_record, args.strip_pairs)
    artifact = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "prior_zeros_l5_sha256": sha256_file(PRIOR_ARTIFACT),
            "exact_arithmetic": "uint64 residues below 2**30; exact CRT; Fraction reference Q(t)",
            "benchmarks_used": "none",
        },
        "definitions": {
            "x": "e^{-2K}, evaluated at the exact rational 2/3",
            "z": "field fugacity satisfying z+z^{-1}=2t",
            "Q_t": "z^{-N/2}P(z) for even N; z^{-(N-1)/2}P(z)/(z+1) for odd N",
            "state_reduction": "D8 spatial symmetry for square cross-sections (row/column/diagonal reflections and rotations), plus spin-flip conjugation in F_p[z]/(z^2-2tz+1)",
        },
        "data": {
            "validation": validation,
            "symmetry_counts": {
                "cross_4x4": burnside_orbit_counts((4, 4)),
                "cross_5x4": burnside_orbit_counts((5, 4)),
                "cross_5x5": burnside_orbit_counts((5, 5)),
            },
            "budget_check_5x5x5": budget,
            "correction_amplitude_route": correction_route_status(),
        },
        "checks": checks,
        "wall_seconds": time.time() - started,
    }
    if not all(item["passed"] for item in checks):
        raise AssertionError("validation check failed")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    for item in checks:
        print(("PASS" if item["passed"] else "FAIL"), item["name"], "--", item["detail"], flush=True)
    print("PASS")


if __name__ == "__main__":
    main()
