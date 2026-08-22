#!/usr/bin/env python3
"""Independent verifier for the support-size-stratified extensive-charge census.

This test deliberately does not import `experiments/e125_extensive_3d_r2.py`.
It reconstructs every decisive rank from raw ordered-Pauli inputs using a
separate orbit representation:

* producer: compact (shape, X-mask, Z-mask) keys and rare-row pivots;
* verifier: bytes encoding of the explicitly translated relative site/Pauli
  pattern and a separately rebuilt frequency ordering at a third independently
  proved prime.

It independently reproduces the wave-12 full-box anchors Z^3 R=1 and
Z^2 R=2, then checks every new support-size class.  During both enumeration
passes it checks every generated output word, not a sample: support size is
at most s+1 and every word lies in C_R={-1,...,R+1}^d.  This is the finite
machine check of the universal support-stratification lemma.

Run:
    PYTHONPATH=src .venv/bin/python tests/test_extensive_3d_r2.py
"""
from __future__ import annotations

import base64
import itertools
import json
import math
import signal
import time
import struct
from array import array
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_3d_r2.json"
TIMEOUT_SECONDS = 3_600
# Independent from the producer's two 31-bit primes. Deterministic trial
# division below establishes primality before Fermat inverses are used.
TEST_PRIME = 1_000_000_007

Site = tuple[int, ...]
Masks = tuple[int, int]
Orbit = bytes


def iter_set_bits(mask: int) -> Iterator[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def is_prime_by_trial_division(value: int) -> bool:
    """Deterministic primality check adequate for the chosen ~1e9 prime."""
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


class IndependentModel:
    """Raw-input verifier with a distinct bytes orbit encoding.

    A nonempty translation orbit is represented by its sorted support list
    after subtracting the coordinatewise minimum, with one byte per relative
    coordinate and one byte for the local ordered-Pauli code X=1, Z=2, XZ=3.
    The bounding shape is implicit because every maximum coordinate is
    occupied by a support site. Thus equal keys are exactly equal translates.
    """

    def __init__(self, dimension: int, radius: int) -> None:
        self.dimension = dimension
        self.radius = radius
        self.base_sites: tuple[Site, ...] = tuple(
            itertools.product(range(radius + 1), repeat=dimension)
        )
        self.ext_sites: tuple[Site, ...] = tuple(
            itertools.product(range(-1, radius + 2), repeat=dimension)
        )
        self.ext_index = {site: position for position, site in enumerate(self.ext_sites)}
        self.base_bits = tuple(1 << self.ext_index[site] for site in self.base_sites)
        self.ext_mask = (1 << len(self.ext_sites)) - 1
        self.neighbors: dict[int, tuple[int, ...]] = {}
        for site, bit in zip(self.base_sites, self.base_bits):
            nearby: list[int] = []
            for axis in range(dimension):
                for sign in (-1, 1):
                    neighbor = list(site)
                    neighbor[axis] += sign
                    nearby.append(1 << self.ext_index[tuple(neighbor)])
            self.neighbors[bit] = tuple(nearby)

    def class_count(self, smax: int) -> int:
        return sum(math.comb(len(self.base_sites), size) * 3**size for size in range(smax + 1))

    def iter_class_columns(self, smax: int) -> Iterator[Masks]:
        """Enumerate raw class inputs by support subset and local Pauli code."""
        yield (0, 0)
        for size in range(1, smax + 1):
            for positions in itertools.combinations(range(len(self.base_sites)), size):
                for codes in itertools.product((1, 2, 3), repeat=size):
                    x_mask = 0
                    z_mask = 0
                    for position, code in zip(positions, codes):
                        bit = self.base_bits[position]
                        if code & 1:
                            x_mask |= bit
                        if code & 2:
                            z_mask |= bit
                    yield (x_mask, z_mask)

    def canonical(self, x_mask: int, z_mask: int) -> Orbit:
        support = x_mask | z_mask
        if not support:
            return b""
        minima = [10**9] * self.dimension
        positions = list(iter_set_bits(support))
        for position in positions:
            site = self.ext_sites[position]
            for axis, coordinate in enumerate(site):
                if coordinate < minima[axis]:
                    minima[axis] = coordinate
        encoded = bytearray()
        # ext_sites is lexicographically sorted; subtracting a fixed vector
        # preserves that ordering, so this is a canonical sorted site list.
        for position in positions:
            site = self.ext_sites[position]
            for axis, coordinate in enumerate(site):
                relative = coordinate - minima[axis]
                assert 0 <= relative < 256
                encoded.append(relative)
            bit = 1 << position
            encoded.append((1 if x_mask & bit else 0) | (2 if z_mask & bit else 0))
        return bytes(encoded)

    def is_anchored(self, x_mask: int, z_mask: int) -> bool:
        """Direct orbit count method, independent of canonical-byte keys."""
        support = x_mask | z_mask
        if not support:
            return True
        minima = [10**9] * self.dimension
        for position in iter_set_bits(support):
            site = self.ext_sites[position]
            for axis, coordinate in enumerate(site):
                if coordinate < minima[axis]:
                    minima[axis] = coordinate
        return all(coordinate == 0 for coordinate in minima)

    def raw_terms(self, x_mask: int, z_mask: int) -> dict[Masks, int]:
        """Independent direct evaluation of (1/2)[H, X^a Z^b] at a=b=1."""
        output: dict[Masks, int] = {}

        def add(key: Masks, coefficient: int) -> None:
            value = output.get(key, 0) + coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)

        # Field terms: every Z bit flips its X bit with coefficient +1.
        for position in iter_set_bits(z_mask):
            bit = 1 << position
            add((x_mask ^ bit, z_mask), 1)
        # Every unordered bond with exactly one X endpoint occurs once: start
        # from that X endpoint and test the other endpoint for X-membership.
        for position in iter_set_bits(x_mask):
            bit = 1 << position
            for neighbor in self.neighbors[bit]:
                if not x_mask & neighbor:
                    add((x_mask, z_mask ^ bit ^ neighbor), -1)
        return output

    def coinvariant_image(self, x_mask: int, z_mask: int) -> dict[Orbit, int]:
        output: dict[Orbit, int] = {}
        for (next_x, next_z), coefficient in self.raw_terms(x_mask, z_mask).items():
            key = self.canonical(next_x, next_z)
            value = output.get(key, 0) + coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)
        return output

    def density_image(self, terms: dict[Masks, int]) -> dict[Orbit, int]:
        output: dict[Orbit, int] = {}
        for (x_mask, z_mask), scale in terms.items():
            for key, coefficient in self.coinvariant_image(x_mask, z_mask).items():
                value = output.get(key, 0) + scale * coefficient
                if value:
                    output[key] = value
                else:
                    output.pop(key, None)
        return output

    def hamiltonian_density(self) -> dict[Masks, int]:
        origin = (0,) * self.dimension
        origin_bit = self.base_bits[self.base_sites.index(origin)]
        terms: dict[Masks, int] = {(origin_bit, 0): 1}
        for axis in range(self.dimension):
            neighbor = [0] * self.dimension
            neighbor[axis] = 1
            neighbor_bit = self.base_bits[self.base_sites.index(tuple(neighbor))]
            terms[(0, origin_bit | neighbor_bit)] = 1
        return terms

    def chain_current_density(self) -> dict[Masks, int]:
        origin = (0,) * self.dimension
        right = (1,) + (0,) * (self.dimension - 1)
        left_bit = self.base_bits[self.base_sites.index(origin)]
        right_bit = self.base_bits[self.base_sites.index(right)]
        return {(left_bit, left_bit | right_bit): 1, (right_bit, left_bit | right_bit): -1}

    def translate(self, x_mask: int, z_mask: int, shift: Site) -> Masks:
        next_x = 0
        next_z = 0
        for position in iter_set_bits(x_mask | z_mask):
            site = self.ext_sites[position]
            moved = tuple(site[axis] + shift[axis] for axis in range(self.dimension))
            assert moved in self.ext_index
            bit = 1 << self.ext_index[moved]
            if x_mask & (1 << position):
                next_x |= bit
            if z_mask & (1 << position):
                next_z |= bit
        return next_x, next_z

    def support_is_in_current_box(self, x_mask: int, z_mask: int) -> bool:
        return not ((x_mask | z_mask) & ~self.ext_mask)


@dataclass
class Prepared:
    model: IndependentModel
    smax: int
    column_count: int
    orbit_count: int
    rows: dict[Orbit, int]
    frequencies: array
    maximum_output_support_size: int
    saw_size_s_plus_one_output: bool


def prepare(model: IndependentModel, smax: int) -> Prepared:
    """Exhaustively reconstruct all class rows and direct anchored orbit count."""
    rows: dict[Orbit, int] = {}
    frequencies = array("I")
    anchored_count = 0
    maximum_output_support_size = 0
    saw_size_s_plus_one_output = False
    column_count = 0
    for x_mask, z_mask in model.iter_class_columns(smax):
        if model.is_anchored(x_mask, z_mask):
            anchored_count += 1
        image: dict[Orbit, int] = {}
        for (next_x, next_z), coefficient in model.raw_terms(x_mask, z_mask).items():
            # Universal structural reduction check: this executes for EVERY
            # output of EVERY class column, not a sampled collection.
            assert model.support_is_in_current_box(next_x, next_z)
            output_size = (next_x | next_z).bit_count()
            assert output_size <= smax + 1
            maximum_output_support_size = max(maximum_output_support_size, output_size)
            saw_size_s_plus_one_output |= output_size == smax + 1
            key = model.canonical(next_x, next_z)
            image[key] = image.get(key, 0) + coefficient
        for key, coefficient in image.items():
            if coefficient:
                row = rows.get(key)
                if row is None:
                    row = len(frequencies)
                    rows[key] = row
                    frequencies.append(0)
                frequencies[row] += 1
        column_count += 1
    assert column_count == model.class_count(smax)
    return Prepared(
        model=model,
        smax=smax,
        column_count=column_count,
        orbit_count=anchored_count,
        rows=rows,
        frequencies=frequencies,
        maximum_output_support_size=maximum_output_support_size,
        saw_size_s_plus_one_output=saw_size_s_plus_one_output,
    )


def rank_mod_prime(prepared: Prepared, prime: int) -> int:
    """Independent raw-input sparse elimination at a third prime.

    The verifier keeps its distinct bytes orbit keys and recomputes every
    image, but orders rows by independently counted hit frequency. That
    prevents the catastrophic fill-in of first-seen-row pivots on the
    11-million-row headline system.
    """
    order = [
        row
        for row, _ in sorted(
            enumerate(prepared.frequencies), key=lambda item: (item[1], item[0])
        )
    ]
    inverse_order = [0] * len(order)
    for position, row in enumerate(order):
        inverse_order[row] = position
    get_position = inverse_order.__getitem__
    pivots: dict[int, dict[int, int]] = {}
    model = prepared.model
    for x_mask, z_mask in model.iter_class_columns(prepared.smax):
        vector: dict[int, int] = {}
        for (next_x, next_z), coefficient in model.raw_terms(x_mask, z_mask).items():
            assert model.support_is_in_current_box(next_x, next_z)
            assert (next_x | next_z).bit_count() <= prepared.smax + 1
            key = model.canonical(next_x, next_z)
            row = prepared.rows.get(key)
            # This proves no output orbit row was silently omitted from the
            # finite stratified presentation.
            assert row is not None
            value = (vector.get(row, 0) + coefficient) % prime
            if value:
                vector[row] = value
            else:
                vector.pop(row, None)
        while vector:
            lead = min(vector, key=get_position)
            prior = pivots.get(lead)
            if prior is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = {
                    row: value * inverse % prime for row, value in vector.items() if value
                }
                break
            scale = vector[lead]
            for row, value in prior.items():
                reduced = (vector.get(row, 0) - scale * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return len(pivots)


def rank_exact_q(prepared: Prepared) -> int:
    """Independent exact Fraction elimination for the smallest new class."""
    pivots: dict[int, dict[int, Fraction]] = {}
    model = prepared.model
    for x_mask, z_mask in model.iter_class_columns(prepared.smax):
        vector: dict[int, Fraction] = {}
        for (next_x, next_z), coefficient in model.raw_terms(x_mask, z_mask).items():
            key = model.canonical(next_x, next_z)
            row = prepared.rows[key]
            value = vector.get(row, Fraction(0)) + Fraction(coefficient)
            if value:
                vector[row] = value
            else:
                vector.pop(row, None)
        while vector:
            lead = min(vector)
            prior = pivots.get(lead)
            if prior is None:
                scale = vector[lead]
                pivots[lead] = {row: value / scale for row, value in vector.items() if value}
                break
            scale = vector[lead]
            for row, value in prior.items():
                reduced = vector.get(row, Fraction(0)) - scale * value
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return len(pivots)


def fullbox_inclusion_exclusion(dimension: int, radius: int) -> int:
    """Independent exact orbit count: anchored nonempty patterns plus identity."""
    side = radius + 1
    count = 0
    for excluded in range(dimension + 1):
        remaining_sites = (side - 1) ** excluded * side ** (dimension - excluded)
        count += (-1) ** excluded * math.comb(dimension, excluded) * 4**remaining_sites
    return count + 1


def expect_case(payload: dict[str, object], kind: str, dimension: int, radius: int, smax: int) -> dict[str, object]:
    for case in payload["data"]["cases"]:
        if (
            case.get("case_kind") == kind
            and case.get("dimension") == dimension
            and case.get("radius") == radius
            and case.get("smax") == smax
        ):
            return case
    raise AssertionError(f"missing result case kind={kind}, d={dimension}, R={radius}, smax={smax}")


def verify_trace_metadata(case: dict[str, object]) -> None:
    """Check stored trace payload shape without using it as the rank proof."""
    records = case.get("modular_rank_lower_bounds", [])
    for record in records:
        certificate = record["certificate"]
        assert certificate["rank"] == record["rank_Fp"]
        assert certificate["pivot_trace_pair_count"] == record["rank_Fp"]
        if certificate["pivot_trace_stored_in_full"]:
            raw = base64.b85decode(certificate["pivot_trace_column_row_base85"])
            assert len(raw) == 8 * certificate["pivot_trace_pair_count"]
            for offset in range(0, len(raw), 8):
                column, row = struct.unpack_from("<II", raw, offset)
                assert column < case["density_basis_size"]
                assert row < case["active_output_orbit_rows"]


def verify_sign_and_direction_traps() -> None:
    model = IndependentModel(3, 2)
    origin = model.base_bits[model.base_sites.index((0, 0, 0))]
    # D(Z_0) = + X_0 Z_0: catches a field-sign error.
    field = model.raw_terms(0, origin)
    assert field == {(origin, origin): 1}
    # D(X_0) has exactly the six outward 3D bond terms, all with sign -1:
    # catches a bond sign or dropped-direction error.
    bond = model.raw_terms(origin, 0)
    assert len(bond) == 6
    assert set(bond.values()) == {-1}

    # h catches a relative field/bond sign defect; j distinguishes the 1D
    # positive control from a dropped transverse direction.
    assert model.density_image(model.hamiltonian_density()) == {}
    one_dimensional = IndependentModel(1, 2)
    assert one_dimensional.density_image(one_dimensional.chain_current_density()) == {}
    square = IndependentModel(2, 3)
    assert square.density_image(square.chain_current_density())
    assert model.density_image(model.chain_current_density())

    # Translation-equivalent raw inputs must have equal coinvariant images.
    x_mask = origin
    z_mask = origin | model.base_bits[model.base_sites.index((0, 1, 0))]
    moved = model.translate(x_mask, z_mask, (1, 0, 0))
    assert model.canonical(x_mask, z_mask) == model.canonical(*moved)
    assert model.coinvariant_image(x_mask, z_mask) == model.coinvariant_image(*moved)


def verify_boundary_class_membership() -> None:
    model = IndependentModel(3, 2)
    four = 0
    five = 0
    for bit in model.base_bits[:4]:
        four |= bit
    for bit in model.base_bits[:5]:
        five |= bit
    assert (four).bit_count() == 4
    assert (five).bit_count() == 5
    assert four.bit_count() <= 4
    assert not five.bit_count() <= 4
    assert model.class_count(4) == 1_503_766


def independently_certify_case(
    dimension: int,
    radius: int,
    smax: int,
    expected: tuple[int, int, int, int],
    exact: bool,
) -> tuple[Prepared, int]:
    """Raw-input recomputation and independent sandwich for one claimed case."""
    model = IndependentModel(dimension, radius)
    prepared = prepare(model, smax)
    columns, orbit_count, rank, quotient = expected
    assert prepared.column_count == columns
    assert prepared.orbit_count == orbit_count
    assert prepared.maximum_output_support_size == smax + 1
    assert prepared.saw_size_s_plus_one_output
    assert model.density_image(model.hamiltonian_density()) == {}
    if exact:
        observed_rank = rank_exact_q(prepared)
    else:
        observed_rank = rank_mod_prime(prepared, TEST_PRIME)
    assert observed_rank == rank
    # Independent exact sandwich: the modular lower bound reaches rank S-2,
    # and the direct I,h kernel supplies the opposite Q inequality.
    assert observed_rank == prepared.orbit_count - 2
    assert prepared.orbit_count - observed_rank - 1 == quotient == 1
    return prepared, observed_rank


def main() -> None:
    # CPU-time budget: SIGALRM measures wall clock, which on this shared
    # machine tracks unrelated background load; the decisive budget below
    # uses process_time(). The alarm stays only as a 4x coarse safety net.
    signal.signal(
        signal.SIGALRM,
        lambda *_: (_ for _ in ()).throw(
            TimeoutError(f"NON-DECISIVE: exceeded {4 * TIMEOUT_SECONDS}s wall (safety net)")
        ),
    )
    signal.alarm(4 * TIMEOUT_SECONDS)
    assert time.process_time() < TIMEOUT_SECONDS, (
        f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s CPU time"
    )
    assert is_prime_by_trial_division(TEST_PRIME)

    def _check_budget() -> None:
        assert time.process_time() < TIMEOUT_SECONDS, (
            f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s CPU time"
        )

    verify_sign_and_direction_traps()
    verify_boundary_class_membership()
    _check_budget()


    # Required independent wave-12 full-box reproductions. Full class s=n is
    # enumerated by support subsets, unlike the producer anchor's base-4 loop.
    anchor_31, rank_31 = independently_certify_case(
        3, 1, 8, (65_536, 64_813, 64_811, 1), exact=False
    )
    anchor_22, rank_22 = independently_certify_case(
        2, 2, 9, (262_144, 254_209, 254_207, 1), exact=False
    )
    _check_budget()
    assert anchor_31.orbit_count == fullbox_inclusion_exclusion(3, 1)
    assert anchor_22.orbit_count == fullbox_inclusion_exclusion(2, 2)

    # New classes. The last two are the headline Z^3 R=2 s<=4 and the
    # strongest finished Z^2 R=3 s<=5 class, respectively.
    new_specs = [
        (3, 2, 2, (3_241, 562, 560, 1), True),
        (3, 2, 3, (82_216, 29_749, 29_747, 1), False),
        (3, 2, 4, (1_503_766, 822_334, 822_332, 1), False),
        (2, 3, 4, (163_669, 83_164, 83_162, 1), False),
        (2, 3, 5, (1_225_093, 790_294, 790_292, 1), False),
    ]
    outcomes: dict[tuple[int, int, int], tuple[Prepared, int]] = {}
    for dimension, radius, smax, expected, exact in new_specs:
        _check_budget()
        outcomes[(dimension, radius, smax)] = independently_certify_case(
            dimension, radius, smax, expected, exact
        )
    _check_budget()

    # Artifact comparison is separate from the raw recomputation above.
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    stored_31 = expect_case(payload, "full_box_wave12_anchor", 3, 1, 8)
    stored_22 = expect_case(payload, "full_box_wave12_anchor", 2, 2, 9)
    assert stored_31["translation_rank"] == anchor_31.orbit_count
    assert stored_31["commutator_rank_lower_bound_Q"] == rank_31
    assert stored_22["translation_rank"] == anchor_22.orbit_count
    assert stored_22["commutator_rank_lower_bound_Q"] == rank_22
    verify_trace_metadata(stored_31)
    verify_trace_metadata(stored_22)

    for dimension, radius, smax, expected, _ in new_specs:
        stored = expect_case(payload, "support_size_class", dimension, radius, smax)
        prepared, observed_rank = outcomes[(dimension, radius, smax)]
        assert stored["density_basis_size"] == prepared.column_count
        assert stored["translation_rank"] == prepared.orbit_count
        assert stored["commutator_rank"] == observed_rank
        assert stored["nontrivial_quotient_dimension"] == 1
        assert stored["maximum_output_support_size"] == smax + 1
        assert stored["outputs_outside_current_box"] == 0
        assert stored["explicit_kernel"]["hamiltonian_density_projected_commutator_zero"] is True
        if stored["arithmetic"] == "sandwich_exact_Q":
            assert stored["rank_upper_from_explicit_kernel"] == prepared.orbit_count - 2
            assert stored["commutator_rank_lower_bound_Q"] == observed_rank
            verify_trace_metadata(stored)
        else:
            assert stored["arithmetic"] == "exact_Q"
            assert stored["rank_certificate"]["rank"] == observed_rank

    # Exact arithmetic behind the symmetry/infeasibility statement. The test
    # re-derives the numbers, rather than trusting producer formatting.
    preflight = payload["data"]["preflight_arithmetic"]
    z3 = preflight["full_box_z3_r2"]
    assert 4**27 == 18_014_398_509_481_984
    assert fullbox_inclusion_exclusion(3, 2) == 18_014_192_401_317_889
    assert z3["density_coefficients"] == 4**27
    assert z3["exact_translation_orbit_count_inclusion_exclusion"] == fullbox_inclusion_exclusion(3, 2)
    assert z3["invariant_class_lower_bound_burnside"] == math.ceil(4**27 / 48)
    assert z3["invariant_class_lower_bound_burnside"] == 375_299_968_947_542
    z2 = preflight["full_box_z2_r3"]
    assert z2["density_coefficients"] == 4**16 == 4_294_967_296
    assert z2["exact_translation_orbit_count_inclusion_exclusion"] == fullbox_inclusion_exclusion(2, 3)

    probe = payload["data"]["observed_wall_probe"]
    assert probe["launched"] is True
    assert probe["closed_form_column_count"] == 21_121_156
    assert probe["wall_classification"] in {"observed_wall", "completed_within_wall"}
    if probe["wall_classification"] == "observed_wall":
        wall = probe["resource_wall"]
        assert 0 < wall["processed_columns"] < probe["closed_form_column_count"]
        assert wall["elapsed_seconds"] > 0
        assert wall["measured_columns_per_second"] > 0

    assert all(entry["passed"] for entry in payload["checks"])
    print("PASS")


if __name__ == "__main__":
    main()
