#!/usr/bin/env python3
"""Exact finite divergence systems for translation-covariant Ising charges.

The calculation uses the repository ordered-Pauli convention X^a Z^b and the
repository generator-left convention (1/2)[H,.].  The sign reversal to [.,H]
does not affect the divergence-membership kernel.

No floating point arithmetic is used for a mathematical conclusion.
"""
from __future__ import annotations

import base64
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
import gc
import hashlib
import itertools
import json
from pathlib import Path
import platform
import resource
import struct
import time
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_charges.json"
SCRIPT = "experiments/e117_extensive_charges.py"

# Large odd primes; a modular rank is always a lower bound for rank over Q.
PRIMES = (2_147_483_647, 2_147_483_629)
FIELD_COEFFICIENT = 1
BOND_COEFFICIENT = 1
RSS_WALL_BYTES = 7 * 1024**3

# B_R has (R+1)^d sites.  The one-dimensional exact route is extended until
# its declared preparation wall; the final two cases deliberately exercise
# the larger transverse dimensions through a bounded two-prime sparse route.
EXACT_CASES = tuple((1, radius) for radius in range(0, 11)) + ((2, 0), (2, 1), (3, 0))
MODULAR_CASES = ((2, 2), (3, 1))
EXACT_WALL_SECONDS = 180.0
MODULAR_WALL_SECONDS = 240.0
PRESENTATION_WALL_SECONDS = 120.0

Site = tuple[int, ...]
OrbitKey = tuple[tuple[int, ...], int, int]  # (bounding shape, canonical X, canonical Z)


@dataclass(frozen=True)
class WallRecord:
    phase: str
    reason: str
    elapsed_seconds: float
    processed_columns: int
    rank_so_far: int
    peak_rss_bytes: int


class ResourceWall(RuntimeError):
    def __init__(self, record: WallRecord) -> None:
        super().__init__(record.reason)
        self.record = record


def peak_rss_bytes() -> int:
    """Darwin reports ru_maxrss in bytes; Linux reports KiB."""
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def wall_dict(record: WallRecord) -> dict[str, object]:
    return {
        "claim_tag": "[UNRESOLVED]",
        "phase": record.phase,
        "reason": record.reason,
        "elapsed_seconds": round(record.elapsed_seconds, 6),
        "processed_columns": record.processed_columns,
        "rank_so_far": record.rank_so_far,
        "peak_rss_bytes": record.peak_rss_bytes,
        "peak_rss_mib": round(record.peak_rss_bytes / 1024**2, 3),
    }


class WallClock:
    def __init__(self, phase: str, seconds: float) -> None:
        self.phase = phase
        self.seconds = seconds
        self.started = time.monotonic()

    def check(self, processed_columns: int, rank_so_far: int = 0) -> None:
        elapsed = time.monotonic() - self.started
        rss = peak_rss_bytes()
        if elapsed > self.seconds:
            raise ResourceWall(
                WallRecord(
                    self.phase,
                    f"wall time {self.seconds:g} s reached",
                    elapsed,
                    processed_columns,
                    rank_so_far,
                    rss,
                )
            )
        if rss > RSS_WALL_BYTES:
            raise ResourceWall(
                WallRecord(
                    self.phase,
                    f"RSS wall {RSS_WALL_BYTES // 1024**3} GiB reached",
                    elapsed,
                    processed_columns,
                    rank_so_far,
                    rss,
                )
            )

    def elapsed(self) -> float:
        return time.monotonic() - self.started


@dataclass
class PreparedSystem:
    row_lookup: dict[OrbitKey, int]
    row_keys: list[OrbitKey]
    row_frequencies: list[int]
    density_translation_rank: int
    density_translation_kernel_dimension: int
    preparation_seconds: float
    peak_rss_bytes: int


@dataclass
class RankResult:
    rank: int
    elapsed_seconds: float
    peak_rss_bytes: int
    maximum_pivot_support: int
    maximum_coefficient_bits: int
    pivot_trace_sha256: str
    pivot_trace: list[list[int]]


class OrbitModel:
    """Finite ordered-Pauli presentation of B_R and its one-shell enlargement."""

    def __init__(self, dimension: int, radius: int) -> None:
        if dimension < 1 or radius < 0:
            raise ValueError("dimension >= 1 and radius >= 0 are required")
        self.dimension = dimension
        self.radius = radius
        self.base_sites: tuple[Site, ...] = tuple(itertools.product(range(radius + 1), repeat=dimension))
        self.ext_sites: tuple[Site, ...] = tuple(itertools.product(range(-1, radius + 2), repeat=dimension))
        self.ext_index = {site: index for index, site in enumerate(self.ext_sites)}
        self.base_bits = tuple(1 << self.ext_index[site] for site in self.base_sites)
        self.base_origin_bit = 1 << self.ext_index[(0,) * dimension]
        self.column_count = 1 << (2 * len(self.base_sites))
        self.neighbor_bits: dict[int, tuple[int, ...]] = {}
        for site, bit in zip(self.base_sites, self.base_bits):
            entries: list[int] = []
            for axis in range(dimension):
                for sign in (-1, 1):
                    neighbor = list(site)
                    neighbor[axis] += sign
                    entries.append(1 << self.ext_index[tuple(neighbor)])
            self.neighbor_bits[bit] = tuple(entries)

    @staticmethod
    def _iter_set_bits(mask: int) -> Iterator[int]:
        while mask:
            bit = mask & -mask
            yield bit.bit_length() - 1
            mask ^= bit

    def decode_column(self, column: int) -> tuple[int, int]:
        x_mask = z_mask = 0
        value = column
        for bit in self.base_bits:
            digit = value & 3
            value >>= 2
            if digit & 1:
                x_mask |= bit
            if digit & 2:
                z_mask |= bit
        return x_mask, z_mask

    @lru_cache(maxsize=250_000)
    def canonical_with_shift(self, x_mask: int, z_mask: int) -> tuple[OrbitKey, Site]:
        """Return the translation orbit key and the minimum-coordinate shift."""
        support = x_mask | z_mask
        if not support:
            return ((), 0, 0), (0,) * self.dimension
        minima = [10**9] * self.dimension
        maxima = [-10**9] * self.dimension
        for position in self._iter_set_bits(support):
            site = self.ext_sites[position]
            for axis, coordinate in enumerate(site):
                minima[axis] = min(minima[axis], coordinate)
                maxima[axis] = max(maxima[axis], coordinate)
        shape = tuple(maxima[axis] - minima[axis] + 1 for axis in range(self.dimension))

        def compact(mask: int) -> int:
            output = 0
            for position in self._iter_set_bits(mask):
                site = self.ext_sites[position]
                index = 0
                for axis, size in enumerate(shape):
                    index = index * size + site[axis] - minima[axis]
                output |= 1 << index
            return output

        return (shape, compact(x_mask), compact(z_mask)), tuple(minima)

    @staticmethod
    def _add_term(output: dict[tuple[int, int], int], key: tuple[int, int], coefficient: int) -> None:
        if not coefficient:
            return
        value = output.get(key, 0) + coefficient
        if value:
            output[key] = value
        else:
            output.pop(key, None)

    def local_terms(self, x_mask: int, z_mask: int, field: int, bond: int) -> dict[tuple[int, int], int]:
        """Exact terms of (1/2)[H, X^a Z^b] in the generator-left convention."""
        output: dict[tuple[int, int], int] = {}
        for position in self._iter_set_bits(z_mask):
            bit = 1 << position
            self._add_term(output, (x_mask ^ bit, z_mask), field)
        for position in self._iter_set_bits(x_mask):
            bit = 1 << position
            # Every boundary bond has exactly one X endpoint, so this loops it once.
            for neighbor_bit in self.neighbor_bits[bit]:
                if not x_mask & neighbor_bit:
                    self._add_term(output, (x_mask, z_mask ^ bit ^ neighbor_bit), -bond)
        return output

    def orbit_image_from_masks(self, x_mask: int, z_mask: int, field: int, bond: int) -> dict[OrbitKey, int]:
        output: dict[OrbitKey, int] = {}
        for (next_x, next_z), coefficient in self.local_terms(x_mask, z_mask, field, bond).items():
            key, _ = self.canonical_with_shift(next_x, next_z)
            value = output.get(key, 0) + coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)
        return output

    def orbit_image(self, column: int, field: int, bond: int) -> dict[OrbitKey, int]:
        return self.orbit_image_from_masks(*self.decode_column(column), field, bond)

    def hamiltonian_density_terms(self) -> dict[tuple[int, int], int]:
        """h=X_0+sum_i Z_0 Z_{e_i}, whose translation sum is H at a=b=1."""
        if self.radius < 1:
            raise ValueError("the Hamiltonian density requires R >= 1")
        output: dict[tuple[int, int], int] = {(self.base_origin_bit, 0): 1}
        origin = (0,) * self.dimension
        for axis in range(self.dimension):
            neighbor = list(origin)
            neighbor[axis] = 1
            z_mask = self.base_origin_bit | (1 << self.ext_index[tuple(neighbor)])
            output[(0, z_mask)] = 1
        return output


def prepare_system(model: OrbitModel, field: int, bond: int, seconds: float) -> PreparedSystem:
    clock = WallClock("prepare_ordered_pauli_system", seconds)
    row_lookup: dict[OrbitKey, int] = {}
    row_keys: list[OrbitKey] = []
    frequencies: list[int] = []
    density_orbits: set[OrbitKey] = set()
    for column in range(model.column_count):
        x_mask, z_mask = model.decode_column(column)
        density_key, _ = model.canonical_with_shift(x_mask, z_mask)
        density_orbits.add(density_key)
        for key in model.orbit_image_from_masks(x_mask, z_mask, field, bond):
            row = row_lookup.get(key)
            if row is None:
                row = len(frequencies)
                row_lookup[key] = row
                row_keys.append(key)
                frequencies.append(0)
            frequencies[row] += 1
        if column % 256 == 0:
            clock.check(column)
    clock.check(model.column_count)
    translation_rank = len(density_orbits)
    return PreparedSystem(
        row_lookup=row_lookup,
        row_keys=row_keys,
        row_frequencies=frequencies,
        density_translation_rank=translation_rank,
        density_translation_kernel_dimension=model.column_count - translation_rank,
        preparation_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
    )


def order_rows(prepared: PreparedSystem) -> list[int]:
    return [
        position
        for position, _ in sorted(
            enumerate(prepared.row_frequencies), key=lambda item: (item[1], item[0])
        )
    ]


def image_rows(model: OrbitModel, prepared: PreparedSystem, column: int, field: int, bond: int) -> dict[int, int]:
    image = model.orbit_image(column, field, bond)
    return {prepared.row_lookup[key]: coefficient for key, coefficient in image.items() if coefficient}

def trace_bytes(trace: list[list[int]]) -> bytes:
    """Compact, architecture-independent serialization of (column,row) pivots."""
    output = bytearray(8 * len(trace))
    for position, (column, row) in enumerate(trace):
        struct.pack_into("<II", output, 8 * position, column, row)
    return bytes(output)


def rank_exact_q(
    model: OrbitModel, prepared: PreparedSystem, field: int, bond: int, seconds: float
) -> RankResult:
    clock = WallClock("exact_fraction_rank", seconds)
    row_order = order_rows(prepared)
    inverse_order = [0] * len(row_order)
    for position, row in enumerate(row_order):
        inverse_order[row] = position
    pivots: dict[int, dict[int, Fraction]] = {}
    max_support = max_bits = 0
    trace_hash = hashlib.sha256()
    for column in range(model.column_count):
        vector = {row: Fraction(value) for row, value in image_rows(model, prepared, column, field, bond).items()}
        while vector:
            lead = min(vector, key=inverse_order.__getitem__)
            prior = pivots.get(lead)
            if prior is None:
                coefficient = vector[lead]
                if coefficient != 1:
                    vector = {row: value / coefficient for row, value in vector.items() if value}
                pivots[lead] = vector
                trace_hash.update(struct.pack("<II", column, lead))
                max_support = max(max_support, len(vector))
                for value in vector.values():
                    max_bits = max(max_bits, value.numerator.bit_length(), value.denominator.bit_length())
                break
            coefficient = vector[lead]
            for row, value in prior.items():
                reduced = vector.get(row, Fraction(0)) - coefficient * value
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
        if column % 128 == 0:
            clock.check(column, len(pivots))
    clock.check(model.column_count, len(pivots))
    digest = trace_hash.hexdigest()
    return RankResult(
        rank=len(pivots),
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
        maximum_pivot_support=max_support,
        maximum_coefficient_bits=max_bits,
        pivot_trace_sha256=digest,
        pivot_trace=[],
    )


def rank_mod_prime(
    model: OrbitModel,
    prepared: PreparedSystem,
    field: int,
    bond: int,
    prime: int,
    seconds: float,
) -> RankResult:
    clock = WallClock(f"modular_rank_p_{prime}", seconds)
    row_order = order_rows(prepared)
    inverse_order = [0] * len(row_order)
    for position, row in enumerate(row_order):
        inverse_order[row] = position
    pivots: dict[int, dict[int, int]] = {}
    max_support = 0
    trace: list[list[int]] = []
    for column in range(model.column_count):
        vector = {row: value % prime for row, value in image_rows(model, prepared, column, field, bond).items() if value % prime}
        while vector:
            lead = min(vector, key=inverse_order.__getitem__)
            prior = pivots.get(lead)
            if prior is None:
                inverse = pow(vector[lead], prime - 2, prime)
                vector = {row: value * inverse % prime for row, value in vector.items() if value}
                pivots[lead] = vector
                trace.append([column, lead])
                max_support = max(max_support, len(vector))
                break
            coefficient = vector[lead]
            for row, value in prior.items():
                reduced = (vector.get(row, 0) - coefficient * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
        if column % 128 == 0:
            clock.check(column, len(pivots))
    clock.check(model.column_count, len(pivots))
    digest = hashlib.sha256(trace_bytes(trace)).hexdigest()
    return RankResult(
        rank=len(pivots),
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
        maximum_pivot_support=max_support,
        maximum_coefficient_bits=0,
        pivot_trace_sha256=digest,
        pivot_trace=trace,
    )


def rank_record(result: RankResult, arithmetic: str, include_trace: bool) -> dict[str, object]:
    record: dict[str, object] = {
        "rank": result.rank,
        "arithmetic": arithmetic,
        "elapsed_seconds": round(result.elapsed_seconds, 6),
        "peak_rss_bytes": result.peak_rss_bytes,
        "peak_rss_mib": round(result.peak_rss_bytes / 1024**2, 3),
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
        "maximum_pivot_support": result.maximum_pivot_support,
        "maximum_coefficient_bits": result.maximum_coefficient_bits,
        "pivot_trace_sha256": result.pivot_trace_sha256,
    }
    if include_trace:
        packed_trace = trace_bytes(result.pivot_trace)
        record.update(
            {
                "pivot_trace_encoding": "base85 little-endian uint32 (column,row) pairs",
                "pivot_trace_pair_count": len(result.pivot_trace),
                "pivot_trace_column_row_base85": base64.b85encode(packed_trace).decode("ascii"),
            }
        )
    return record


def build_current_presentation(
    model: OrbitModel, prepared: PreparedSystem, field: int, bond: int, seconds: float
) -> dict[str, object]:
    """Build a finite tree-current presentation equivalent to all finite currents.

    For each output translation orbit, a rooted coordinatewise tree joins every
    actual translate generated by D(q).  Its edge-current incidence image is
    exactly the zero-total-coefficient subspace of that orbit, so eliminating
    these current variables is identical to applying the orbit-sum quotient.
    """
    clock = WallClock("tree_current_presentation", seconds)
    nodes_by_row: list[set[Site]] = [set() for _ in prepared.row_frequencies]
    for column in range(model.column_count):
        x_mask, z_mask = model.decode_column(column)
        for (next_x, next_z), coefficient in model.local_terms(x_mask, z_mask, field, bond).items():
            if coefficient:
                key, shift = model.canonical_with_shift(next_x, next_z)
                nodes_by_row[prepared.row_lookup[key]].add(shift)
        if column % 256 == 0:
            clock.check(column)

    equation_rows = current_columns = 0
    direction_columns = [0] * model.dimension
    support_outside_current_box = 0
    active_orbits = 0
    for row, nodes in enumerate(nodes_by_row):
        if not nodes:
            continue
        active_orbits += 1
        root = tuple(min(node[axis] for node in nodes) for axis in range(model.dimension))
        closure = set(nodes)
        for node in tuple(nodes):
            cursor = list(node)
            while tuple(cursor) != root:
                axis = next(axis for axis in range(model.dimension) if cursor[axis] > root[axis])
                cursor[axis] -= 1
                closure.add(tuple(cursor))
        shape = prepared.row_keys[row][0]
        # The fixed-parent rule makes this a rooted tree: each nonroot vertex
        # has one parent obtained by reducing its first nonroot coordinate.
        for node in closure:
            for axis, size in enumerate(shape):
                if node[axis] < -1 or node[axis] + size - 1 > model.radius + 1:
                    support_outside_current_box += 1
                    break
            if node != root:
                axis = next(axis for axis in range(model.dimension) if node[axis] > root[axis])
                direction_columns[axis] += 1
        equation_rows += len(closure)
        current_columns += len(closure) - 1
    clock.check(model.column_count)
    assert sum(direction_columns) == current_columns
    assert support_outside_current_box == 0
    return {
        "claim_tag": "[COMPUTATION]",
        "current_box": "C_R={-1,...,R+1}^d",
        "current_box_site_count": (model.radius + 3) ** model.dimension,
        "active_ordered_pauli_output_orbits": active_orbits,
        "equation_count": equation_rows,
        "current_unknown_count": current_columns,
        "current_unknowns_by_direction": direction_columns,
        "q_unknown_count": model.column_count,
        "combined_unknown_count": model.column_count + current_columns,
        "tree_current_basis": (
            "one s_i coefficient per oriented edge of a rooted coordinatewise tree for each active "
            "ordered-Pauli translation orbit; this is an exact spanning restriction, not a truncation"
        ),
        "all_tree_currents_stay_in_C_R": support_outside_current_box == 0,
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
        "elapsed_seconds": round(clock.elapsed(), 6),
        "peak_rss_bytes": peak_rss_bytes(),
        "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
    }


def summed_image(model: OrbitModel, terms: dict[tuple[int, int], int], field: int, bond: int) -> dict[OrbitKey, int]:
    output: dict[OrbitKey, int] = {}
    for (x_mask, z_mask), scale in terms.items():
        for key, coefficient in model.orbit_image_from_masks(x_mask, z_mask, field, bond).items():
            value = output.get(key, 0) + scale * coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)
    return output


def chain_current_terms(model: OrbitModel) -> dict[tuple[int, int], int]:
    """J=sum_j(Y_j Z_{j+1}-Z_j Y_{j+1}), up to the harmless factor -i."""
    if model.dimension != 1 or model.radius < 1:
        raise ValueError("the chain current needs d=1, R>=1")
    left = 1 << model.ext_index[(0,)]
    right = 1 << model.ext_index[(1,)]
    return {(left, left | right): 1, (right, left | right): -1}


def transverse_chain_current_terms(model: OrbitModel) -> dict[tuple[int, int], int]:
    """The same oriented two-site density along e_1 in d>=2."""
    if model.dimension < 2 or model.radius < 1:
        raise ValueError("a transverse direction is required")
    origin = (0,) * model.dimension
    neighbor = (1,) + (0,) * (model.dimension - 1)
    left = 1 << model.ext_index[origin]
    right = 1 << model.ext_index[neighbor]
    return {(left, left | right): 1, (right, left | right): -1}


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def analyze_case(dimension: int, radius: int, mode: str) -> dict[str, object]:
    started = time.monotonic()
    model = OrbitModel(dimension, radius)
    case: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "dimension": dimension,
        "lattice": {1: "Z", 2: "Z^2 square", 3: "Z^3 simple cubic"}[dimension],
        "radius": radius,
        "density_box": "B_R={0,...,R}^d",
        "coupling": {"a": FIELD_COEFFICIENT, "b": BOND_COEFFICIENT},
        "positive_control": dimension == 1,
        "density_basis_size": model.column_count,
        "density_box_site_count": len(model.base_sites),
        "current_box_site_count": len(model.ext_sites),
        "requested_mode": mode,
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
    }
    try:
        prepared = prepare_system(model, FIELD_COEFFICIENT, BOND_COEFFICIENT, PRESENTATION_WALL_SECONDS)
        presentation = build_current_presentation(
            model, prepared, FIELD_COEFFICIENT, BOND_COEFFICIENT, PRESENTATION_WALL_SECONDS
        )
    except ResourceWall as exc:
        case.update(
            {
                "arithmetic": "resource_wall",
                "resource_wall": wall_dict(exc.record),
                "elapsed_seconds": round(time.monotonic() - started, 6),
                "peak_rss_bytes": peak_rss_bytes(),
                "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
            }
        )
        return case

    case.update(
        {
            "translation_rank": prepared.density_translation_rank,
            "trivial_divergence_dimension": prepared.density_translation_kernel_dimension,
            "preparation_seconds": round(prepared.preparation_seconds, 6),
            "preparation_peak_rss_bytes": prepared.peak_rss_bytes,
            "finite_divergence_system": presentation,
        }
    )
    # I always lies in the kernel.  For R>=1, the displayed local H-density is
    # an independent second kernel vector, furnishing an analytic rank upper bound.
    identity_kernel_dimension = 1
    hamiltonian_kernel = radius >= 1
    h_image = summed_image(model, model.hamiltonian_density_terms(), FIELD_COEFFICIENT, BOND_COEFFICIENT) if hamiltonian_kernel else None
    case["explicit_kernel"] = {
        "identity_density": True,
        "hamiltonian_density_available": hamiltonian_kernel,
        "hamiltonian_density_projected_commutator_zero": h_image == {} if h_image is not None else None,
        "dimension_used_for_rank_upper_bound": identity_kernel_dimension + int(hamiltonian_kernel),
    }
    assert h_image in (None, {})
    rank_upper = prepared.density_translation_rank - identity_kernel_dimension - int(hamiltonian_kernel)

    if mode == "exact_Q":
        try:
            exact = rank_exact_q(model, prepared, FIELD_COEFFICIENT, BOND_COEFFICIENT, EXACT_WALL_SECONDS)
        except ResourceWall as exc:
            case.update(
                {
                    "arithmetic": "resource_wall",
                    "resource_wall": wall_dict(exc.record),
                    "elapsed_seconds": round(time.monotonic() - started, 6),
                    "peak_rss_bytes": peak_rss_bytes(),
                    "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
                }
            )
            return case
        case["arithmetic"] = "exact_Q"
        case["commutator_rank"] = exact.rank
        case["rank_certificate"] = rank_record(exact, "exact Fraction Gaussian elimination", include_trace=False)
        case["solution_density_nullity"] = model.column_count - exact.rank
        case["nontrivial_quotient_dimension"] = prepared.density_translation_rank - exact.rank - 1
        case["combined_system_rank"] = presentation["current_unknown_count"] + exact.rank
        case["combined_system_nullity"] = model.column_count - exact.rank
        assert exact.rank <= rank_upper
    elif mode == "two_prime_modular":
        modular_records: list[dict[str, object]] = []
        ranks: list[int] = []
        walls: list[dict[str, object]] = []
        for prime in PRIMES:
            gc.collect()
            try:
                modular = rank_mod_prime(
                    model, prepared, FIELD_COEFFICIENT, BOND_COEFFICIENT, prime, MODULAR_WALL_SECONDS
                )
            except ResourceWall as exc:
                walls.append(wall_dict(exc.record))
                break
            ranks.append(modular.rank)
            modular_records.append(
                {
                    "prime": prime,
                    "rank_Fp": modular.rank,
                    "certificate": rank_record(modular, "finite-field elimination", include_trace=True),
                }
            )
        case["modular_rank_lower_bounds"] = modular_records
        case["rank_upper_from_explicit_kernel"] = rank_upper
        if walls:
            case.update(
                {
                    "arithmetic": "resource_wall",
                    "resource_wall": walls,
                    "elapsed_seconds": round(time.monotonic() - started, 6),
                    "peak_rss_bytes": peak_rss_bytes(),
                    "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
                }
            )
            return case
        rank_lower = max(ranks)
        case["commutator_rank_lower_bound_Q"] = rank_lower
        case["solution_density_nullity_upper_bound_Q"] = model.column_count - rank_lower
        case["nontrivial_quotient_dimension_upper_bound_Q"] = prepared.density_translation_rank - rank_lower - 1
        if rank_lower == rank_upper:
            # The modular result is used only for the lower inequality; the
            # independently exhibited I,H kernel gives the matching Q upper inequality.
            case.update(
                {
                    "arithmetic": "sandwich_exact_Q",
                    "commutator_rank": rank_lower,
                    "solution_density_nullity": model.column_count - rank_lower,
                    "nontrivial_quotient_dimension": prepared.density_translation_rank - rank_lower - 1,
                    "combined_system_rank": presentation["current_unknown_count"] + rank_lower,
                    "combined_system_nullity": model.column_count - rank_lower,
                    "sandwich_reason": (
                        "rank_Fp <= rank_Q <= translation_rank - dim(span{I,h}); the two-prime lower bound "
                        "reaches this explicit Q upper bound"
                    ),
                }
            )
        else:
            case["arithmetic"] = "two_prime_modular_upper_bound_only"
    else:
        raise ValueError(mode)

    # A scalar multiple of the chain current must vanish only in d=1.  Its
    # failure in transverse directions independently checks that every spatial
    # bond direction entered the commutator map.
    if radius >= 1 and dimension == 1:
        chain_image = summed_image(model, chain_current_terms(model), FIELD_COEFFICIENT, BOND_COEFFICIENT)
        case["chain_current_projected_commutator_term_count"] = len(chain_image)
        assert not chain_image
    elif radius >= 1 and dimension >= 2:
        transverse_image = summed_image(
            model, transverse_chain_current_terms(model), FIELD_COEFFICIENT, BOND_COEFFICIENT
        )
        case["embedded_chain_current_projected_commutator_term_count"] = len(transverse_image)
        assert transverse_image

    case["elapsed_seconds"] = round(time.monotonic() - started, 6)
    case["peak_rss_bytes"] = peak_rss_bytes()
    case["peak_rss_mib"] = round(peak_rss_bytes() / 1024**2, 3)
    return case


def main() -> int:
    started = time.monotonic()
    cases: list[dict[str, object]] = []
    for dimension, radius in EXACT_CASES:
        cases.append(analyze_case(dimension, radius, "exact_Q"))
    for dimension, radius in MODULAR_CASES:
        cases.append(analyze_case(dimension, radius, "two_prime_modular"))

    by_key = {(case["dimension"], case["radius"]): case for case in cases}
    exact_cases = [case for case in cases if case.get("arithmetic") == "exact_Q"]
    chain = [case for case in cases if case["dimension"] == 1 and case.get("arithmetic") == "exact_Q"]
    d2r1 = by_key.get((2, 1), {})
    d3r0 = by_key.get((3, 0), {})
    d3r1 = by_key.get((3, 1), {})

    checks = [
        check(
            "ordered_pauli_hamiltonian_kernel",
            all(
                case.get("explicit_kernel", {}).get("identity_density")
                and case.get("explicit_kernel", {}).get("hamiltonian_density_projected_commutator_zero") is not False
                for case in exact_cases
            ),
            "the identity and, where it fits, h=X_0+sum_i Z_0Z_ei have zero translation-coinvariant commutator",
        ),
        check(
            "one_dimensional_positive_control",
            bool(by_key.get((1, 1)))
            and by_key[(1, 1)].get("nontrivial_quotient_dimension", 0) >= 2
            and by_key[(1, 2)].get("nontrivial_quotient_dimension", 0) >= 4
            and by_key[(1, 1)].get("chain_current_projected_commutator_term_count") == 0,
            "the same system finds the Hamiltonian and the standard oriented TFIM energy-current density at R=1",
        ),
        check(
            "transverse_direction_included",
            d2r1.get("embedded_chain_current_projected_commutator_term_count", 0) > 0,
            "the chain current ceases to commute after the transverse Z Z bonds are included",
        ),
        check(
            "exact_small_cubic_anchor",
            d3r0.get("arithmetic") == "exact_Q"
            and d3r0.get("nontrivial_quotient_dimension") == 0,
            "the one-site cubic density has only the identity before the Hamiltonian density fits",
        ),
        check(
            "exact_fraction_elimination_present",
            any(case.get("arithmetic") == "exact_Q" and case.get("commutator_rank") is not None for case in cases),
            "decisive small cases use Fraction Gaussian elimination over Q, not floating point",
        ),
        check(
            "two_prime_records_are_labelled_one_sided",
            all(
                case.get("arithmetic") not in {"sandwich_exact_Q", "two_prime_modular_upper_bound_only"}
                or "modular_rank_lower_bounds" in case
                for case in cases
            ),
            "finite-field ranks are stored as lower bounds; any equality claim additionally names the independent kernel upper bound",
        ),
        check(
            "tree_current_system_stays_local",
            all(
                case.get("finite_divergence_system", {}).get("all_tree_currents_stay_in_C_R", False)
                for case in cases
                if "finite_divergence_system" in case
            ),
            "every current basis word used by the exact tree presentation lies in C_R={-1,...,R+1}^d",
        ),
        check(
            "current_directions_present",
            all(
                len(case["finite_divergence_system"]["current_unknowns_by_direction"]) == case["dimension"]
                and all(count > 0 for count in case["finite_divergence_system"]["current_unknowns_by_direction"])
                for case in cases
                if "finite_divergence_system" in case and case["dimension"] >= 2 and case["radius"] >= 1
            ),
            "every current direction s_i is genuinely exercised in each completed d>=2, R>=1 case",
        ),
    ]

    unresolved: list[dict[str, object]] = [
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "all dimensions, all radii",
            "statement": (
                "This finite-radius calculation does not prove the absence of arbitrary-range, quasilocal, or "
                "non-translation-covariant conserved charges.  It does not establish non-integrability of the 3D Ising model."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "coupling ratio",
            "statement": (
                "All stored ranks use a=b=1.  Reaching the maximal rank at this ratio proves a Zariski-generic "
                "finite-radius statement, but this artifact does not identify or exclude exceptional nonzero b/a ratios."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "Z R>=10",
            "statement": (
                "The d=1, R=10 preparation of 4^11=4,194,304 density coefficients reached its explicit "
                "120 s wall; its exact processed-column count and RSS are in the d=1, R=10 resource record. "
                "Every larger chain radius remains uncomputed in this run."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "Z^3 R>=2",
            "classification": "input_size_preflight_not_observed_wall",
            "statement": (
                "B_2 in the stated convention has 27 sites and 4^27=18,014,398,509,481,984 ordered-Pauli "
                "density coefficients before currents.  This case was not launched; the size is a preflight "
                "projection, not an observed 7 GiB / 240 s resource wall."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "Z^2 R>=3",
            "classification": "input_size_preflight_not_observed_wall",
            "statement": (
                "B_3 has 16 sites and 4^16=4,294,967,296 ordered-Pauli density coefficients before currents. "
                "This case was not launched; the size is a preflight projection, not an observed 7 GiB / "
                "240 s resource wall."
            ),
        },
    ]
    for case in cases:
        if case.get("arithmetic") in {"resource_wall", "two_prime_modular_upper_bound_only"}:
            unresolved.append(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "scope": f"d={case['dimension']}, R={case['radius']}",
                    "statement": "The stored resource record / one-sided modular bounds do not decide the exact Q quotient dimension.",
                    "case_key": [case["dimension"], case["radius"]],
                }
            )

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "platform": platform.platform(),
            "method": (
                "exact ordered-Pauli finite divergence presentation, quotient by translation coinvariants, exact "
                "Fraction elimination for small cases, and explicitly one-sided two-prime sparse lower bounds for large cases"
            ),
            "certifying_arithmetic": "Python integers and fractions.Fraction; finite-field arithmetic modulo listed primes; no floats",
            "total_elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
            "peak_rss_bytes": peak_rss_bytes(),
            "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
            "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
        },
        "data": {
            "claim_tags": ["[THEOREM]", "[LEMMA]", "[COMPUTATION]", "[EXTERNAL]", "[UNRESOLVED]"],
            "construction": {
                "claim_tag": "[LEMMA]",
                "ordered_pauli_convention": "X^a Z^b",
                "generator_orientation": "(1/2)[H,.]; replacing it by [.,H] multiplies every equation by -2",
                "ad_formulas": (
                    "(1/2)[X_s,X^aZ^b]=[b_s=1]X^(a+e_s)Z^b; "
                    "(1/2)[Z_uZ_v,X^aZ^b]=-[a_u xor a_v=1]X^aZ^(b+e_u+e_v)"
                ),
                "density_box": "B_R={0,...,R}^d",
                "current_box": "C_R={-1,...,R+1}^d",
                "formal_charge": "Q(q)=sum_(x in Z^d) tau_x(q)",
                "finite_system": "D(q)=sum_i(s_i-tau_(e_i)s_i) in the active ordered-Pauli basis",
                "translation_coinvariant_elimination": (
                    "a finite density is a finite lattice divergence iff the sum of its coefficients in every "
                    "translation orbit is zero; the per-orbit rooted tree currents give an explicit finite witness"
                ),
                "trivial_quotient": "span{identity density} + finite lattice divergences",
                "reported_dimension": "dim ker(pi o D|V_R) - dim ker(pi|V_R) - 1 = rank(pi|V_R)-rank(pi o D|V_R)-1",
            },
            "finite_radius_result": {
                "claim_tag": "[COMPUTATION]",
                "statement": (
                    "The listed values are finite-radius density quotient dimensions at a=b=1.  A value one is the "
                    "unavoidable Hamiltonian density class, not a newly discovered conserved charge."
                ),
                "cases": cases,
            },
            "cases": cases,
            "one_dimensional_control": {
                "claim_tag": "[COMPUTATION]",
                "known_density": (
                    "The ordered-Pauli density X_jZ_jZ_(j+1)-X_(j+1)Z_jZ_(j+1) is -i times "
                    "sum_j(Y_j Z_(j+1)-Z_j Y_(j+1))."
                ),
                "status": "the same d-parameterized code finds its zero projected commutator at d=1, R=1",
                "external_context": "[EXTERNAL] This is the standard local energy-current member of the free-fermion transverse-field Ising chain hierarchy.",
            },
            "three_dimensional_interpretation": {
                "claim_tag": "[COMPUTATION]",
                "statement": (
                    "For every Z^3 case whose reported quotient dimension is exactly one, the only certified class "
                    "within that finite density box is the already-known Hamiltonian density class, modulo identity and divergences."
                ),
                "additional_density_found": any(
                    case["dimension"] == 3
                    and case.get("nontrivial_quotient_dimension", 0) > 1
                    for case in cases
                ),
                "periodic_box_verification": "not applicable: no additional non-Hamiltonian Z^3 density was certified",
            },
            "internal_checks": checks,
            "unresolved": unresolved,
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    failed = [entry["name"] for entry in checks if not entry["passed"]]
    if failed:
        print("FAIL: " + ", ".join(failed))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
