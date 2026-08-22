#!/usr/bin/env python3
"""Support-size-stratified extensive-charge census past the wave-12 frontier.

Wave 12 (`experiments/e117_extensive_charges.py`, `proofs/extensive_charges.md`)
certified the translation-coinvariant density quotient

    dim ker(pi o D|_V) / (ker(pi|_V) + Q I) = rank S_R - rank M_R - 1

for full boxes up to Z^3 R=1 and Z^2 R=2 and stopped: Z^3 R=2 has 4^27
density columns and Z^2 R=3 has 4^16, both beyond any streaming budget.

This script pushes the 3D census into the Z^3 R=2 (3x3x3) box and the Z^2
R=3 (4x4) box by STRATIFYING BY SUPPORT SIZE instead of enumerating the full
ordered-Pauli box.  For s >= 0 let C(d,R,s) be the rational span of all
ordered-Pauli words supported in B_R={0..R}^d whose word support (number of
non-identity sites) is at most s.  Lemma A of `proofs/extensive_3d_r2.md`
proves the stratified system exact: D maps a size-<=s word to size-<=s+1
words inside C_R={-1..R+1}^d, so pi o D|_C is presented exactly by the
translation-orbit rows that the class columns actually hit.  The class
quotient is again rank S_C - rank M_C - 1 (Lemma B: I lies in every class
and {q in C: pi(q) in Q e_empty} = Q I + (C cap ker pi)).

The rank over Q is certified by the wave-12 sandwich, not by a bare modular
equality: rank_Fp <= rank_Q always, and the analytically exhibited kernel
{I, h} (h = X_0 + sum_i Z_0 Z_{e_i}, of word sizes 1 and 2) gives
rank_Q M_C <= rank S_C - 2.  When the two-prime lower bound meets that upper
bound, rank_Q is certified and the class quotient follows.

Cases (a=b=1, exact integer arithmetic only):
  * wave-12 full-box anchors Z^3 R=1 and Z^2 R=2 (must reproduce the frozen
    values 64813/64811 and 254209/254207 exactly);
  * Z^3 R=2 classes s<=2 (exact Fraction rank), s<=3, s<=4;
  * Z^2 R=3 classes s<=4, s<=5;
  * a time-walled OBSERVED probe of Z^3 R=2 s<=5 (21,121,156 columns) that
    records measured throughput and RSS rather than guessing;
  * preflight arithmetic for the unreachable full boxes, including the exact
    inclusion-exclusion anchored-orbit count of the full Z^3 R=2 box and the
    Burnside lower bound showing point-group symmetry reduction (factor <=48)
    cannot make that box finite work.

Run: PYTHONPATH=src .venv/bin/python experiments/e125_extensive_3d_r2.py
"""
from __future__ import annotations

import base64
import gc
import hashlib
import itertools
import json
import math
import platform
import resource
import struct
import time
from array import array
from dataclasses import dataclass, field
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_3d_r2.json"
FROZEN_WAVE12 = ROOT / "results" / "integrability" / "extensive_charges.json"
SCRIPT = "experiments/e125_extensive_3d_r2.py"

PRIMES = (2_147_483_647, 2_147_483_629)
FIELD_COEFFICIENT = 1
BOND_COEFFICIENT = 1
RSS_WALL_BYTES = 5 * 1024**3
RSS_CHECK_COLUMN_INTERVAL = 8_192
PROBE_WALL_SECONDS = 480.0
# Store full pivot traces only below this rank; bigger cases store the SHA-256.
TRACE_RANK_LIMIT = 300_000

Site = tuple[int, ...]
OrbitKey = tuple[tuple[int, ...], int, int]
MASKS = tuple[int, int]


@dataclass(frozen=True)
class WallRecord:
    phase: str
    reason: str
    processed_columns: int
    total_columns: int
    elapsed_seconds: float
    peak_rss_bytes: int

    def as_dict(self) -> dict[str, object]:
        return {
            "phase": self.phase,
            "reason": self.reason,
            "processed_columns": self.processed_columns,
            "total_columns": self.total_columns,
            "elapsed_seconds": round(self.elapsed_seconds, 6),
            "peak_rss_bytes": self.peak_rss_bytes,
            "peak_rss_mib": round(self.peak_rss_bytes / 1024**2, 3),
            "measured_columns_per_second": (
                round(self.processed_columns / self.elapsed_seconds, 3)
                if self.elapsed_seconds > 0
                else None
            ),
        }


class ResourceWall(RuntimeError):
    def __init__(self, record: WallRecord) -> None:
        super().__init__(record.reason)
        self.record = record


def peak_rss_bytes() -> int:
    """Darwin reports ru_maxrss in bytes; Linux reports KiB."""
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class WallClock:
    def __init__(self, phase: str, seconds: float) -> None:
        self.phase = phase
        self.limit = seconds
        self.started = time.monotonic()

    def elapsed(self) -> float:
        return time.monotonic() - self.started

    def check(self, processed: int, total: int) -> None:
        spent = self.elapsed()
        if spent > self.limit:
            raise ResourceWall(
                WallRecord(
                    phase=self.phase,
                    reason=f"explicit {self.limit:.1f}s wall exceeded after {spent:.1f}s",
                    processed_columns=processed,
                    total_columns=total,
                    elapsed_seconds=spent,
                    peak_rss_bytes=peak_rss_bytes(),
                )
            )


def check_rss(phase: str, processed: int, total: int, elapsed_seconds: float) -> None:
    rss = peak_rss_bytes()
    if rss > RSS_WALL_BYTES:
        raise ResourceWall(
            WallRecord(
                phase=phase,
                reason=f"peak RSS {rss / 1024**3:.3f} GiB exceeded the {RSS_WALL_BYTES / 1024**3:.1f} GiB wall",
                processed_columns=processed,
                total_columns=total,
                elapsed_seconds=elapsed_seconds,
                peak_rss_bytes=rss,
            )
        )


def iter_set_bits(mask: int) -> Iterator[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


class ClassModel:
    """Ordered-Pauli presentation of the size-<=smax word class in B_R."""

    def __init__(self, dimension: int, radius: int, smax: int) -> None:
        if dimension < 1 or radius < 0 or smax < 0:
            raise ValueError("dimension >= 1, radius >= 0, smax >= 0 required")
        if smax > (radius + 1) ** dimension:
            raise ValueError("smax beyond the box site count is the full box; use full-box mode")
        self.dimension = dimension
        self.radius = radius
        self.smax = smax
        self.base_sites: tuple[Site, ...] = tuple(
            itertools.product(range(radius + 1), repeat=dimension)
        )
        self.ext_sites: tuple[Site, ...] = tuple(
            itertools.product(range(-1, radius + 2), repeat=dimension)
        )
        self.ext_index = {site: index for index, site in enumerate(self.ext_sites)}
        self.base_bits = tuple(1 << self.ext_index[site] for site in self.base_sites)
        self.neighbor_bits: dict[int, tuple[int, ...]] = {}
        for site, bit in zip(self.base_sites, self.base_bits):
            entries: list[int] = []
            for axis in range(dimension):
                for sign in (-1, 1):
                    neighbor = list(site)
                    neighbor[axis] += sign
                    entries.append(1 << self.ext_index[tuple(neighbor)])
            self.neighbor_bits[bit] = tuple(entries)

    def column_count(self) -> int:
        n = len(self.base_sites)
        return sum(math.comb(n, k) * 3**k for k in range(self.smax + 1))

    def iter_columns(self) -> Iterator[MASKS]:
        """Every (x_mask, z_mask) word on B_R with support size <= smax."""
        n = len(self.base_sites)
        yield (0, 0)
        for size in range(1, self.smax + 1):
            for positions in itertools.combinations(range(n), size):
                bits = [self.base_bits[p] for p in positions]
                for assignment in itertools.product((1, 2, 3), repeat=size):
                    x_mask = 0
                    z_mask = 0
                    for bit, pauli in zip(bits, assignment):
                        if pauli & 1:
                            x_mask |= bit
                        if pauli & 2:
                            z_mask |= bit
                    yield (x_mask, z_mask)

    def local_terms(self, x_mask: int, z_mask: int) -> dict[MASKS, int]:
        """Exact terms of (1/2)[H, X^a Z^b] at a=b=1, generator-left convention.

        (1/2)[X_s, X^a Z^b] = [b_s=1] X^(a+e_s) Z^b
        (1/2)[Z_u Z_v, X^a Z^b] = -[a_u xor a_v] X^a Z^(b+e_u+e_v)
        """
        output: dict[MASKS, int] = {}
        mask = z_mask
        while mask:
            bit = mask & -mask
            mask ^= bit
            key = (x_mask ^ bit, z_mask)
            output[key] = output.get(key, 0) + FIELD_COEFFICIENT
        mask = x_mask
        while mask:
            bit = mask & -mask
            mask ^= bit
            for neighbor_bit in self.neighbor_bits[bit]:
                if not x_mask & neighbor_bit:
                    key = (x_mask, z_mask ^ bit ^ neighbor_bit)
                    output[key] = output.get(key, 0) - BOND_COEFFICIENT
        return {key: value for key, value in output.items() if value}

    def anchor(self, x_mask: int, z_mask: int) -> OrbitKey:
        """Translation-orbit key: anchored at min coordinates, compacted."""
        support = x_mask | z_mask
        if not support:
            return ((), 0, 0)
        minima = [10**9] * self.dimension
        maxima = [-10**9] * self.dimension
        for position in iter_set_bits(support):
            site = self.ext_sites[position]
            for axis, coordinate in enumerate(site):
                if coordinate < minima[axis]:
                    minima[axis] = coordinate
                if coordinate > maxima[axis]:
                    maxima[axis] = coordinate
        shape = tuple(maxima[axis] - minima[axis] + 1 for axis in range(self.dimension))
        compact_x = 0
        compact_z = 0
        for position in iter_set_bits(support):
            site = self.ext_sites[position]
            index = 0
            for axis, size in enumerate(shape):
                index = index * size + site[axis] - minima[axis]
            bit = 1 << index
            if x_mask & (1 << position):
                compact_x |= bit
            if z_mask & (1 << position):
                compact_z |= bit
        return (shape, compact_x, compact_z)

    def support_size(self, masks: MASKS) -> int:
        return bin(masks[0] | masks[1]).count("1")

    def support_inside_current_box(self, masks: MASKS) -> bool:
        for position in iter_set_bits(masks[0] | masks[1]):
            site = self.ext_sites[position]
            for coordinate in site:
                if coordinate < -1 or coordinate > self.radius + 1:
                    return False
        return True

    def orbit_image(self, x_mask: int, z_mask: int) -> dict[OrbitKey, int]:
        output: dict[OrbitKey, int] = {}
        for (next_x, next_z), coefficient in self.local_terms(x_mask, z_mask).items():
            key = self.anchor(next_x, next_z)
            value = output.get(key, 0) + coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)
        return output

    def hamiltonian_density_terms(self) -> dict[MASKS, int]:
        """h = X_0 + sum_i Z_0 Z_{e_i}; word support sizes 1 and 2."""
        if self.radius < 1:
            raise ValueError("the Hamiltonian density requires R >= 1")
        origin_bit = self.base_bits[self.base_sites.index((0,) * self.dimension)]
        output: dict[MASKS, int] = {(origin_bit, 0): 1}
        for axis in range(self.dimension):
            neighbor = [0] * self.dimension
            neighbor[axis] = 1
            bit = self.base_bits[self.base_sites.index(tuple(neighbor))]
            output[(0, origin_bit | bit)] = 1
        return output

    def chain_current_terms(self) -> dict[MASKS, int]:
        """j = X_0 Z_0 Z_1 - X_1 Z_0 Z_1 along e_1 (word support size 3)."""
        if self.radius < 1:
            raise ValueError("the chain current needs R >= 1")
        origin = (0,) * self.dimension
        right = (1,) + (0,) * (self.dimension - 1)
        left_bit = self.base_bits[self.base_sites.index(origin)]
        right_bit = self.base_bits[self.base_sites.index(right)]
        return {(left_bit, left_bit | right_bit): 1, (right_bit, left_bit | right_bit): -1}

    def summed_orbit_image(self, terms: dict[MASKS, int]) -> dict[OrbitKey, int]:
        output: dict[OrbitKey, int] = {}
        for (x_mask, z_mask), scale in terms.items():
            for key, coefficient in self.orbit_image(x_mask, z_mask).items():
                value = output.get(key, 0) + scale * coefficient
                if value:
                    output[key] = value
                else:
                    output.pop(key, None)
        return output


@dataclass
class ClassSystem:
    """CSR presentation of pi o D restricted to the class columns."""

    column_count: int
    orbit_count: int  # rank of S|_C = number of translation orbits met by the class
    row_count: int  # number of output orbits hit (rows of the exact presentation)
    row_frequencies: array
    csr_rows: array
    csr_coeffs: array
    csr_offsets: array
    max_output_support_size: int
    outputs_outside_current_box: int
    elapsed_seconds: float
    peak_rss_bytes: int


def build_class_system(model: ClassModel, wall_seconds: float) -> ClassSystem:
    """Enumerate class columns once; cache coinvariant images in CSR arrays.

    Lemma A's two structural facts are asserted per column: every output word
    has support size <= smax+1 and support inside C_R.
    """
    total = model.column_count()
    clock = WallClock("enumerate_class_columns", wall_seconds)
    row_lookup: dict[OrbitKey, int] = {}
    frequencies = array("q")
    column_orbits: set[OrbitKey] = set()
    csr_rows = array("i")
    csr_coeffs = array("b")
    csr_offsets = array("q", [0])
    ext_mask = (1 << len(model.ext_sites)) - 1
    max_output_support_size = 0
    outside = 0
    processed = 0
    for x_mask, z_mask in model.iter_columns():
        column_orbits.add(model.anchor(x_mask, z_mask))
        per_row: dict[int, int] = {}
        for (next_x, next_z), coefficient in model.local_terms(x_mask, z_mask).items():
            size = (next_x | next_z).bit_count()
            if size > max_output_support_size:
                max_output_support_size = size
            if size > model.smax + 1:
                raise AssertionError(
                    f"Lemma A violated: output support size {size} > smax+1 = {model.smax + 1}"
                )
            if (next_x | next_z) & ~ext_mask:
                outside += 1
                raise AssertionError("Lemma A violated: output support outside C_R")
            key = model.anchor(next_x, next_z)
            row = row_lookup.get(key)
            if row is None:
                row = len(frequencies)
                row_lookup[key] = row
                frequencies.append(0)
            frequencies[row] += 1
            per_row[row] = per_row.get(row, 0) + coefficient
        for row, coefficient in per_row.items():
            if coefficient:
                csr_rows.append(row)
                csr_coeffs.append(coefficient)
        csr_offsets.append(len(csr_rows))
        processed += 1
        if processed % 1024 == 0:
            clock.check(processed, total)
            if processed % RSS_CHECK_COLUMN_INTERVAL == 0:
                check_rss("enumerate_class_columns", processed, total, clock.elapsed())
    clock.check(processed, total)
    if processed != total:
        raise AssertionError(f"enumerated {processed} columns, expected closed form {total}")
    system = ClassSystem(
        column_count=processed,
        orbit_count=len(column_orbits),
        row_count=len(frequencies),
        row_frequencies=frequencies,
        csr_rows=csr_rows,
        csr_coeffs=csr_coeffs,
        csr_offsets=csr_offsets,
        max_output_support_size=max_output_support_size,
        outputs_outside_current_box=outside,
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
    )
    del row_lookup, column_orbits
    gc.collect()
    return system


def row_order_arrays(system: ClassSystem) -> list[int]:
    """Rare rows first (wave-12 ordering) to keep elimination fill small."""
    order = [
        row
        for row, _ in sorted(
            enumerate(system.row_frequencies), key=lambda item: (item[1], item[0])
        )
    ]
    return order


def column_vector(system: ClassSystem, column: int) -> dict[int, int]:
    start = system.csr_offsets[column]
    end = system.csr_offsets[column + 1]
    return {
        system.csr_rows[pos]: system.csr_coeffs[pos]
        for pos in range(start, end)
    }


@dataclass
class RankResult:
    rank: int
    elapsed_seconds: float
    peak_rss_bytes: int
    maximum_pivot_support: int
    pivot_trace_sha256: str
    pivot_trace_pair_count: int
    pivot_trace: list[list[int]] = field(default_factory=list)


def trace_digest_update(digest: "hashlib._Hash", column: int, row: int) -> None:
    digest.update(struct.pack("<II", column, row))


def rank_mod_prime(system: ClassSystem, prime: int, wall_seconds: float) -> RankResult:
    clock = WallClock(f"modular_rank_p_{prime}", wall_seconds)
    order = row_order_arrays(system)
    inverse_order = [0] * len(order)
    for position, row in enumerate(order):
        inverse_order[row] = position
    get_position = inverse_order.__getitem__
    pivots: dict[int, dict[int, int]] = {}
    max_support = 0
    trace: list[list[int]] = []
    keep_trace = True
    digest = hashlib.sha256()
    for column in range(system.column_count):
        vector = {
            row: value % prime
            for row, value in column_vector(system, column).items()
            if value % prime
        }
        while vector:
            lead = min(vector, key=get_position)
            prior = pivots.get(lead)
            if prior is None:
                inverse = pow(vector[lead], prime - 2, prime)
                vector = {
                    row: value * inverse % prime for row, value in vector.items() if value
                }
                pivots[lead] = vector
                trace_digest_update(digest, column, lead)
                if keep_trace:
                    trace.append([column, lead])
                    if len(trace) > TRACE_RANK_LIMIT:
                        keep_trace = False
                        trace = []
                max_support = max(max_support, len(vector))
                break
            scale = vector[lead]
            for row, value in prior.items():
                reduced = (vector.get(row, 0) - scale * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
        if column % 4096 == 0:
            clock.check(column, system.column_count)
            if column % RSS_CHECK_COLUMN_INTERVAL == 0:
                check_rss(
                    f"modular_rank_p_{prime}", column, system.column_count, clock.elapsed()
                )
    clock.check(system.column_count, system.column_count)
    return RankResult(
        rank=len(pivots),
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
        maximum_pivot_support=max_support,
        pivot_trace_sha256=digest.hexdigest(),
        pivot_trace_pair_count=len(pivots),
        pivot_trace=trace,
    )


def rank_exact_q(system: ClassSystem, wall_seconds: float) -> RankResult:
    clock = WallClock("exact_fraction_rank", wall_seconds)
    order = row_order_arrays(system)
    inverse_order = [0] * len(order)
    for position, row in enumerate(order):
        inverse_order[row] = position
    get_position = inverse_order.__getitem__
    pivots: dict[int, dict[int, Fraction]] = {}
    max_support = 0
    max_bits = 0
    digest = hashlib.sha256()
    trace: list[list[int]] = []
    for column in range(system.column_count):
        vector = {
            row: Fraction(value)
            for row, value in column_vector(system, column).items()
        }
        while vector:
            lead = min(vector, key=get_position)
            prior = pivots.get(lead)
            if prior is None:
                scale = vector[lead]
                if scale != 1:
                    vector = {
                        row: value / scale for row, value in vector.items() if value
                    }
                pivots[lead] = vector
                trace.append([column, lead])
                trace_digest_update(digest, column, lead)
                max_support = max(max_support, len(vector))
                for value in vector.values():
                    max_bits = max(
                        max_bits, value.numerator.bit_length(), value.denominator.bit_length()
                    )
                break
            scale = vector[lead]
            for row, value in prior.items():
                reduced = vector.get(row, Fraction(0)) - scale * value
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
        if column % 512 == 0:
            clock.check(column, system.column_count)
            if column % RSS_CHECK_COLUMN_INTERVAL == 0:
                check_rss("exact_fraction_rank", column, system.column_count, clock.elapsed())
    clock.check(system.column_count, system.column_count)
    return RankResult(
        rank=len(pivots),
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
        maximum_pivot_support=max_support,
        pivot_trace_sha256=digest.hexdigest(),
        pivot_trace_pair_count=len(pivots),
        pivot_trace=trace,
    )


def trace_bytes(trace: list[list[int]]) -> bytes:
    output = bytearray(8 * len(trace))
    for position, (column, row) in enumerate(trace):
        struct.pack_into("<II", output, 8 * position, column, row)
    return bytes(output)


def rank_certificate(result: RankResult, arithmetic: str) -> dict[str, object]:
    record: dict[str, object] = {
        "rank": result.rank,
        "arithmetic": arithmetic,
        "elapsed_seconds": round(result.elapsed_seconds, 6),
        "peak_rss_bytes": result.peak_rss_bytes,
        "peak_rss_mib": round(result.peak_rss_bytes / 1024**2, 3),
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
        "maximum_pivot_support": result.maximum_pivot_support,
        "pivot_trace_sha256": result.pivot_trace_sha256,
        "pivot_trace_pair_count": result.pivot_trace_pair_count,
        "pivot_trace_encoding": "base85 little-endian uint32 (column,row) pairs",
        "pivot_trace_stored_in_full": bool(result.pivot_trace),
    }
    if result.pivot_trace:
        record["pivot_trace_column_row_base85"] = base64.b85encode(
            trace_bytes(result.pivot_trace)
        ).decode("ascii")
    return record


def fullbox_anchored_orbit_count(dimension: int, radius: int) -> int:
    """Exact inclusion-exclusion count of translation orbits met by the full box.

    A word anchored at the origin (min coordinate 0 in every axis, among its
    support sites) is the unique representative of its orbit inside the box;
    the identity is one extra orbit.  Not touching the plane x_i = 0 confines
    the support to a box with one fewer layer per excluded axis.
    """
    side = radius + 1
    count = 0
    for excluded in range(dimension + 1):
        remaining_sites = (side - 1) ** excluded * side ** (dimension - excluded)
        count += (-1) ** excluded * math.comb(dimension, excluded) * 4**remaining_sites
    # The identity word has empty support and forms its own orbit.
    return count + 1


def analyze_class_case(dimension: int, radius: int, smax: int, mode: str) -> dict[str, object]:
    started = time.monotonic()
    model = ClassModel(dimension, radius, smax)
    n_sites = len(model.base_sites)
    case: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "case_kind": "support_size_class",
        "dimension": dimension,
        "lattice": {1: "Z", 2: "Z^2 square", 3: "Z^3 simple cubic"}[dimension],
        "radius": radius,
        "density_box": "B_R={0,...,R}^d",
        "class_definition": (
            f"rational span of ordered-Pauli words supported in B_{radius} with word "
            f"support size (non-identity sites) <= {smax}"
        ),
        "smax": smax,
        "coupling": {"a": FIELD_COEFFICIENT, "b": BOND_COEFFICIENT},
        "density_box_site_count": n_sites,
        "current_box": "C_R={-1,...,R+1}^d",
        "closed_form_column_count": model.column_count(),
        "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
    }
    wall = max(120.0, model.column_count() / 4000.0)
    try:
        system = build_class_system(model, wall)
    except ResourceWall as exc:
        case.update(
            {
                "arithmetic": "resource_wall",
                "resource_wall": exc.record.as_dict(),
                "wall_classification": (
                    "observed_wall" if exc.record.elapsed_seconds > 0 else "observed_rss_wall"
                ),
                "elapsed_seconds": round(time.monotonic() - started, 6),
                "peak_rss_bytes": peak_rss_bytes(),
                "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
            }
        )
        return case

    case.update(
        {
            "density_basis_size": system.column_count,
            "translation_rank": system.orbit_count,
            "trivial_divergence_dimension": system.column_count - system.orbit_count,
            "active_output_orbit_rows": system.row_count,
            "maximum_output_support_size": system.max_output_support_size,
            "outputs_outside_current_box": system.outputs_outside_current_box,
            "enumeration_seconds": round(system.elapsed_seconds, 6),
            "enumeration_peak_rss_bytes": system.peak_rss_bytes,
        }
    )
    assert system.column_count == model.column_count()
    assert system.outputs_outside_current_box == 0
    assert system.max_output_support_size <= smax + 1

    # Explicit kernel: the identity and the Hamiltonian density h (word sizes
    # 1 and 2) lie in every class with smax >= 2 and R >= 1; both must have
    # exactly zero coinvariant commutator.
    hamiltonian_available = radius >= 1 and smax >= 2
    h_image = (
        model.summed_orbit_image(model.hamiltonian_density_terms())
        if hamiltonian_available
        else None
    )
    assert h_image in (None, {})
    case["explicit_kernel"] = {
        "identity_density": True,
        "hamiltonian_density_available": hamiltonian_available,
        "hamiltonian_density_word_support_sizes": [1] + [2] * dimension,
        "hamiltonian_density_projected_commutator_zero": h_image == {},
        "dimension_used_for_rank_upper_bound": 1 + int(hamiltonian_available),
    }
    rank_upper = system.orbit_count - 1 - int(hamiltonian_available)

    if mode == "exact_Q":
        result = rank_exact_q(system, wall)
        case["arithmetic"] = "exact_Q"
        case["commutator_rank"] = result.rank
        case["rank_certificate"] = rank_certificate(result, "exact Fraction Gaussian elimination")
        case["solution_density_nullity"] = system.column_count - result.rank
        case["nontrivial_quotient_dimension"] = system.orbit_count - result.rank - 1
        assert result.rank <= rank_upper
    elif mode == "two_prime_modular":
        records: list[dict[str, object]] = []
        ranks: list[int] = []
        walls: list[dict[str, object]] = []
        for prime in PRIMES:
            gc.collect()
            try:
                result = rank_mod_prime(system, prime, wall)
            except ResourceWall as exc:
                walls.append(exc.record.as_dict())
                break
            ranks.append(result.rank)
            records.append({"prime": prime, "rank_Fp": result.rank, "certificate": rank_certificate(result, "finite-field elimination")})
        case["modular_rank_lower_bounds"] = records
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
        case["solution_density_nullity_upper_bound_Q"] = system.column_count - rank_lower
        case["nontrivial_quotient_dimension_upper_bound_Q"] = system.orbit_count - rank_lower - 1
        if rank_lower == rank_upper:
            case.update(
                {
                    "arithmetic": "sandwich_exact_Q",
                    "commutator_rank": rank_lower,
                    "solution_density_nullity": system.column_count - rank_lower,
                    "nontrivial_quotient_dimension": system.orbit_count - rank_lower - 1,
                    "sandwich_reason": (
                        "rank_Fp <= rank_Q <= orbit_count - dim span{I,h}; the two-prime "
                        "lower bound meets the analytic kernel upper bound"
                    ),
                }
            )
        else:
            case["arithmetic"] = "two_prime_modular_upper_bound_only"
    else:
        raise ValueError(mode)

    # Directional control: the oriented chain current (word size 3) must be a
    # null density along e_1 only when no transverse bonds act, and a nonzero
    # coinvariant image once transverse directions exist.
    if radius >= 1 and smax >= 3:
        j_image = model.summed_orbit_image(model.chain_current_terms())
        if dimension == 1:
            case["chain_current_projected_commutator_term_count"] = len(j_image)
            assert not j_image
        else:
            case["embedded_chain_current_projected_commutator_term_count"] = len(j_image)
            assert j_image

    case["elapsed_seconds"] = round(time.monotonic() - started, 6)
    case["peak_rss_bytes"] = peak_rss_bytes()
    case["peak_rss_mib"] = round(peak_rss_bytes() / 1024**2, 3)
    return case


def analyze_fullbox_anchor(dimension: int, radius: int) -> dict[str, object]:
    """Reproduce a frozen wave-12 full-box value inside this script's pipeline."""
    started = time.monotonic()
    model = ClassModel(dimension, radius, (radius + 1) ** dimension)
    total = 1 << (2 * len(model.base_sites))
    case: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "case_kind": "full_box_wave12_anchor",
        "dimension": dimension,
        "radius": radius,
        "smax": (radius + 1) ** dimension,
        "density_box_site_count": len(model.base_sites),
        "closed_form_column_count": total,
        "coupling": {"a": FIELD_COEFFICIENT, "b": BOND_COEFFICIENT},
    }
    # Full-box enumeration by integer decode, matching the wave-12 basis.
    clock = WallClock("enumerate_full_box", 240.0)
    row_lookup: dict[OrbitKey, int] = {}
    frequencies = array("q")
    column_orbits: set[OrbitKey] = set()
    csr_rows = array("i")
    csr_coeffs = array("b")
    csr_offsets = array("q", [0])
    base_bits = model.base_bits
    for column in range(total):
        x_mask = 0
        z_mask = 0
        value = column
        for bit in base_bits:
            digit = value & 3
            value >>= 2
            if digit & 1:
                x_mask |= bit
            if digit & 2:
                z_mask |= bit
        column_orbits.add(model.anchor(x_mask, z_mask))
        per_row: dict[int, int] = {}
        for (next_x, next_z), coefficient in model.local_terms(x_mask, z_mask).items():
            key = model.anchor(next_x, next_z)
            row = row_lookup.get(key)
            if row is None:
                row = len(frequencies)
                row_lookup[key] = row
                frequencies.append(0)
            frequencies[row] += 1
            per_row[row] = per_row.get(row, 0) + coefficient
        for row, coefficient in per_row.items():
            if coefficient:
                csr_rows.append(row)
                csr_coeffs.append(coefficient)
        csr_offsets.append(len(csr_rows))
        if column % 1024 == 0:
            clock.check(column, total)
    system = ClassSystem(
        column_count=total,
        orbit_count=len(column_orbits),
        row_count=len(frequencies),
        row_frequencies=frequencies,
        csr_rows=csr_rows,
        csr_coeffs=csr_coeffs,
        csr_offsets=csr_offsets,
        max_output_support_size=0,
        outputs_outside_current_box=0,
        elapsed_seconds=clock.elapsed(),
        peak_rss_bytes=peak_rss_bytes(),
    )
    del row_lookup, column_orbits
    gc.collect()
    records = []
    ranks = []
    for prime in PRIMES:
        result = rank_mod_prime(system, prime, 240.0)
        ranks.append(result.rank)
        records.append({"prime": prime, "rank_Fp": result.rank, "certificate": rank_certificate(result, "finite-field elimination")})
    case.update(
        {
            "density_basis_size": total,
            "translation_rank": system.orbit_count,
            "active_output_orbit_rows": system.row_count,
            "inclusion_exclusion_orbit_count": fullbox_anchored_orbit_count(dimension, radius),
            "modular_rank_lower_bounds": records,
            "arithmetic": "two_prime_modular",
            "commutator_rank_lower_bound_Q": max(ranks),
        }
    )
    assert system.orbit_count == fullbox_anchored_orbit_count(dimension, radius)
    case["elapsed_seconds"] = round(time.monotonic() - started, 6)
    case["peak_rss_mib"] = round(peak_rss_bytes() / 1024**2, 3)
    return case


def probe_observed_wall(dimension: int, radius: int, smax: int, wall_seconds: float) -> dict[str, object]:
    """Launch the next class deliberately under an explicit wall; record facts."""
    started = time.monotonic()
    model = ClassModel(dimension, radius, smax)
    total = model.column_count()
    record: dict[str, object] = {
        "claim_tag": "[UNRESOLVED]",
        "case_kind": "observed_wall_probe",
        "dimension": dimension,
        "radius": radius,
        "smax": smax,
        "closed_form_column_count": total,
        "launched": True,
        "wall_seconds": wall_seconds,
    }
    try:
        system = build_class_system(model, wall_seconds)
    except ResourceWall as exc:
        wall = exc.record.as_dict()
        record.update(
            {
                "resource_wall": wall,
                "wall_classification": "observed_wall",
                "statement": (
                    f"The d={dimension}, R={radius}, smax<={smax} class ({total:,} columns by "
                    f"closed form) was launched and stopped by its explicit {wall_seconds:.0f}s "
                    f"wall after {wall['processed_columns']:,} columns "
                    f"({wall['measured_columns_per_second']} columns/s measured, "
                    f"{wall['peak_rss_mib']} MiB process-lifetime RSS at the wall). "
                    "The projected completion time at the measured rate is recorded below; "
                    "no rank or quotient claim is made for this case."
                ),
                "projected_completion_seconds_at_measured_rate": (
                    round(total / wall["measured_columns_per_second"], 1)
                    if wall["measured_columns_per_second"]
                    else None
                ),
                "elapsed_seconds": round(time.monotonic() - started, 6),
            }
        )
        return record
    # Unexpected completion within the wall: keep the honest result.
    record.update(
        {
            "wall_classification": "completed_within_wall",
            "column_count": system.column_count,
            "translation_rank": system.orbit_count,
            "active_output_orbit_rows": system.row_count,
            "enumeration_seconds": round(system.elapsed_seconds, 6),
            "statement": "The probe completed enumeration inside its wall; no rank claim is attached to the probe itself.",
        }
    )
    return record


def symmetry_and_fullbox_preflight() -> dict[str, object]:
    """Arithmetic (no launch) for the boxes this route deliberately does not run."""
    full_columns_z3r2 = 4**27
    full_columns_z2r3 = 4**16
    point_group_order = 48  # octahedral group O_h acting on the 27 sites
    burnside_orbits = -(-full_columns_z3r2 // point_group_order)  # ceiling
    return {
        "claim_tag": "[COMPUTATION]",
        "kind": "input_size_preflight_not_observed_wall",
        "statement": (
            "Counts and projections computed by exact integer arithmetic; nothing at this "
            "size was launched. These are input-size preflights, not observed resource walls."
        ),
        "full_box_z3_r2": {
            "density_coefficients": full_columns_z3r2,
            "exact_translation_orbit_count_inclusion_exclusion": fullbox_anchored_orbit_count(3, 2),
            "symmetry_reduction_note": (
                "The bond-preserving site-permutation group of the 3x3x3 box is the octahedral "
                "group of order 48. By Burnside, the invariant-density class still has at least "
                "ceil(4^27/48) = "
                f"{burnside_orbits:,} independent coefficients, so point-group symmetry "
                "reduction alone cannot make the full Z^3 R=2 box finite work at any "
                "achievable rate; symmetry sectors of the class computation remain finite only "
                "because of the support-size stratification."
            ),
            "invariant_class_lower_bound_burnside": burnside_orbits,
        },
        "full_box_z2_r3": {
            "density_coefficients": full_columns_z2r3,
            "exact_translation_orbit_count_inclusion_exclusion": fullbox_anchored_orbit_count(2, 3),
        },
        "class_z3_r2_s5_columns": sum(
            math.comb(27, k) * 3**k for k in range(6)
        ),
    }


def load_frozen_wave12_values() -> dict[tuple[int, int], tuple[int, int]]:
    if not FROZEN_WAVE12.exists():
        return {}
    payload = json.loads(FROZEN_WAVE12.read_text(encoding="utf-8"))
    output: dict[tuple[int, int], tuple[int, int]] = {}
    for case in payload.get("data", {}).get("cases", []):
        key = (case.get("dimension"), case.get("radius"))
        translation = case.get("translation_rank")
        rank = case.get("commutator_rank")
        if translation is not None and rank is not None:
            output[key] = (translation, rank)
    return output


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    started = time.monotonic()
    frozen = load_frozen_wave12_values()

    anchor_31 = analyze_fullbox_anchor(3, 1)
    anchor_22 = analyze_fullbox_anchor(2, 2)
    cases = [
        anchor_31,
        anchor_22,
        analyze_class_case(3, 2, 2, "exact_Q"),
        analyze_class_case(3, 2, 3, "two_prime_modular"),
        analyze_class_case(3, 2, 4, "two_prime_modular"),
        analyze_class_case(2, 3, 4, "two_prime_modular"),
        analyze_class_case(2, 3, 5, "two_prime_modular"),
    ]
    probe = probe_observed_wall(3, 2, 5, PROBE_WALL_SECONDS)
    preflight = symmetry_and_fullbox_preflight()

    by_class = {
        (case["dimension"], case["radius"], case["smax"]): case
        for case in cases
        if case.get("case_kind") == "support_size_class"
    }
    z3_s4 = by_class.get((3, 2, 4), {})
    z2_s5 = by_class.get((2, 3, 5), {})

    checks = [
        check(
            "wave12_anchor_z3_r1_reproduced",
            anchor_31.get("translation_rank") == 64813
            and anchor_31.get("commutator_rank_lower_bound_Q") == 64811
            and frozen.get((3, 1), (64813, 64811)) == (64813, 64811),
            "this script's pipeline independently reproduces the frozen Z^3 R=1 full-box values 64813/64811",
        ),
        check(
            "wave12_anchor_z2_r2_reproduced",
            anchor_22.get("translation_rank") == 254209
            and anchor_22.get("commutator_rank_lower_bound_Q") == 254207
            and frozen.get((2, 2), (254209, 254207)) == (254209, 254207),
            "this script's pipeline independently reproduces the frozen Z^2 R=2 full-box values 254209/254207",
        ),
        check(
            "anchored_orbit_count_matches_inclusion_exclusion",
            anchor_31.get("translation_rank") == anchor_31.get("inclusion_exclusion_orbit_count")
            and anchor_22.get("translation_rank") == anchor_22.get("inclusion_exclusion_orbit_count"),
            "the orbit count equals the closed-form anchored inclusion-exclusion count for both anchors",
        ),
        check(
            "hamiltonian_kernel_present_in_every_class",
            all(
                case.get("explicit_kernel", {}).get("hamiltonian_density_projected_commutator_zero")
                is True
                for case in by_class.values()
            ),
            "h = X_0 + sum_i Z_0 Z_ei has exactly zero coinvariant commutator in every class",
        ),
        check(
            "stratified_rows_exact",
            all(
                case.get("outputs_outside_current_box") == 0
                and case.get("maximum_output_support_size", 0) <= case["smax"] + 1
                for case in by_class.values()
            ),
            "Lemma A machine check: every output of every class column stays in C_R with support size <= smax+1",
        ),
        check(
            "headline_z3_r2_class_s4_certified",
            z3_s4.get("arithmetic") == "sandwich_exact_Q"
            and z3_s4.get("nontrivial_quotient_dimension") == 1,
            "Z^3 R=2 class smax<=4: two-prime rank meets the analytic kernel upper bound; quotient = 1",
        ),
        check(
            "z2_r3_class_s5_certified",
            z2_s5.get("arithmetic") == "sandwich_exact_Q"
            and z2_s5.get("nontrivial_quotient_dimension") == 1,
            "Z^2 R=3 class smax<=5: two-prime rank meets the analytic kernel upper bound; quotient = 1",
        ),
        check(
            "transverse_direction_control",
            any(
                case.get("embedded_chain_current_projected_commutator_term_count", 0) > 0
                for case in by_class.values()
            ),
            "the oriented chain current fails to be conserved once transverse bonds act (bond-direction trap)",
        ),
        check(
            "probe_is_an_observed_wall",
            probe.get("wall_classification") in {"observed_wall", "completed_within_wall"},
            "the smax<=5 Z^3 R=2 probe was launched and records measured time/RSS, not a projection",
        ),
    ]

    class_cases = [case for case in cases if case.get("case_kind") == "support_size_class"]
    unresolved: list[dict[str, object]] = [
        probe,
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "full boxes Z^3 R>=2 and Z^2 R>=3",
            "classification": "input_size_preflight_not_observed_wall",
            "statement": (
                "The full ordered-Pauli boxes (4^27 resp. 4^16 density coefficients) were not "
                "launched; the preflight counts in this artifact are arithmetic projections."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "classes with smax >= 5 in Z^3 R=2",
            "classification": "observed_wall_probe" if probe.get("wall_classification") == "observed_wall" else "input_size_preflight_not_observed_wall",
            "statement": (
                "Only the launched probe data constrain Z^3 R=2 classes with smax >= 5; no rank "
                "or quotient is claimed for them. Densities whose words have support size >= 5 "
                "remain unchecked in the 3x3x3 box."
            ),
        },
        {
            "claim_tag": "[UNRESOLVED]",
            "scope": "all-range and coupling-ratio statements",
            "statement": (
                "All certified values are finite-box, finite-class, a=b=1 statements. Nothing "
                "here rules out quasilocal, non-translation-covariant, larger-support, or "
                "exceptional-ratio conserved quantities, and nothing decides integrability of "
                "the 3D Ising model."
            ),
        },
    ]
    for case in class_cases:
        if case.get("arithmetic") not in {"sandwich_exact_Q", "exact_Q"}:
            unresolved.append(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "scope": f"class d={case['dimension']}, R={case['radius']}, smax<={case['smax']}",
                    "statement": "The stored one-sided bounds do not certify this class quotient.",
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
                "support-size-stratified exact ordered-Pauli class systems; translation "
                "coinvariant (orbit-sum) rows; exact Fraction elimination for the smallest "
                "class; streamed CSR column storage with two-prime finite-field lower bounds "
                "sandwiched against the analytic {I,h} kernel upper bound"
            ),
            "certifying_arithmetic": "Python integers and fractions.Fraction; finite-field arithmetic modulo listed primes; no floats",
            "rss_wall_bytes": RSS_WALL_BYTES,
            "total_elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
            "peak_rss_bytes": peak_rss_bytes(),
            "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
            "rss_measurement": "process-lifetime ru_maxrss at measurement end; conservative",
        },
        "data": {
            "claim_tags": ["[THEOREM]", "[LEMMA]", "[COMPUTATION]", "[UNRESOLVED]"],
            "construction": {
                "claim_tag": "[LEMMA]",
                "ordered_pauli_convention": "X^a Z^b",
                "generator_orientation": "(1/2)[H,.]; the divergence criterion of proofs/extensive_charges.md is invariant under this orientation",
                "ad_formulas": (
                    "(1/2)[X_s,X^aZ^b]=[b_s=1]X^(a+e_s)Z^b; "
                    "(1/2)[Z_uZ_v,X^aZ^b]=-[a_u xor a_v=1]X^aZ^(b+e_u+e_v)"
                ),
                "class": "C(d,R,s) = span of words supported in B_R={0..R}^d with word support size <= s",
                "class_quotient": (
                    "dim ker(pi o D|_C)/(C cap ker pi + Q I) = rank S_C - rank M_C - 1, exact by "
                    "Lemmas A and B of proofs/extensive_3d_r2.md"
                ),
                "sandwich": "rank_Fp(M_C) <= rank_Q(M_C) <= rank S_C - 2 (kernel I and h)",
                "reported_dimension": "rank S_C - rank M_C - 1",
            },
            "class_cases": class_cases,
            "wave12_anchor_reproductions": [anchor_31, anchor_22],
            "cases": cases,
            "observed_wall_probe": probe,
            "preflight_arithmetic": preflight,
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
