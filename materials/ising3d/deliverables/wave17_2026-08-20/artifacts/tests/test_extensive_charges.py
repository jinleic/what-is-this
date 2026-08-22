#!/usr/bin/env python3
"""Independent standalone verifier for the extensive-charge certificate.

This verifier deliberately does not import the producer.  It rebuilds the
ordered-Pauli commutator in its own tuple representation, takes the finite
translation-coinvariant quotient directly, and compares the decisive exact
small-box ranks recorded by the artifact.
"""
from __future__ import annotations

import base64
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import signal
import struct

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "extensive_charges.json"
TIMEOUT_SECONDS = 900

Site = tuple[int, ...]
Word = tuple[tuple[Site, int, int], ...]


def box(dim: int, radius: int) -> tuple[Site, ...]:
    """The producer's explicit convention: B_R={0,...,R}^dim."""
    return tuple(itertools.product(range(radius + 1), repeat=dim))


def decode(column: int, sites: tuple[Site, ...]) -> tuple[frozenset[Site], frozenset[Site]]:
    x_sites: set[Site] = set()
    z_sites: set[Site] = set()
    value = column
    for site in sites:
        digit = value & 3
        value >>= 2
        if digit & 1:
            x_sites.add(site)
        if digit & 2:
            z_sites.add(site)
    return frozenset(x_sites), frozenset(z_sites)


def translate_word(x_sites: frozenset[Site], z_sites: frozenset[Site], shift: Site) -> Word:
    entries: list[tuple[Site, int, int]] = []
    for site in sorted(x_sites | z_sites):
        moved = tuple(site[i] + shift[i] for i in range(len(site)))
        entries.append((moved, int(site in x_sites), int(site in z_sites)))
    return tuple(entries)


def canonical(x_sites: frozenset[Site], z_sites: frozenset[Site], dim: int) -> Word:
    support = x_sites | z_sites
    if not support:
        return ()
    minima = tuple(min(site[i] for site in support) for i in range(dim))
    return translate_word(x_sites, z_sites, tuple(-coordinate for coordinate in minima))

def canonical_with_shift(
    x_sites: frozenset[Site], z_sites: frozenset[Site], dim: int
) -> tuple[Word, Site]:
    support = x_sites | z_sites
    if not support:
        return (), (0,) * dim
    minima = tuple(min(site[i] for site in support) for i in range(dim))
    return canonical(x_sites, z_sites, dim), minima


def local_commutator(
    x_sites: frozenset[Site], z_sites: frozenset[Site], dim: int
) -> dict[tuple[frozenset[Site], frozenset[Site]], int]:
    """Unquotiented (1/2)[H, X^a Z^b] at a=b=1, in actual coordinates."""
    output: dict[tuple[frozenset[Site], frozenset[Site]], int] = {}

    def add(x_part: frozenset[Site], z_part: frozenset[Site], coefficient: int) -> None:
        key = (x_part, z_part)
        value = output.get(key, 0) + coefficient
        if value:
            output[key] = value
        else:
            output.pop(key, None)

    for site in z_sites:
        add(x_sites ^ {site}, z_sites, 1)
    edges: set[tuple[Site, Site]] = set()
    for site in x_sites:
        for axis in range(dim):
            for sign in (-1, 1):
                neighbor = list(site)
                neighbor[axis] += sign
                edges.add(tuple(sorted((site, tuple(neighbor)))))
    for left, right in edges:
        if (left in x_sites) ^ (right in x_sites):
            add(x_sites, z_sites ^ {left, right}, -1)
    return output


def combine_coinvariant(
    terms: dict[tuple[frozenset[Site], frozenset[Site]], int], dim: int
) -> dict[Word, int]:
    """Coinvariant commutator of a signed density, term by term."""
    output: dict[Word, int] = {}
    for (x_sites, z_sites), scale in terms.items():
        for key, coefficient in commutator_coinvariant(x_sites, z_sites, dim).items():
            value = output.get(key, 0) + scale * coefficient
            if value:
                output[key] = value
            else:
                output.pop(key, None)
    return output


def orbit_sum_vector(
    terms: dict[tuple[frozenset[Site], frozenset[Site]], int], dim: int
) -> dict[Word, int]:
    """The class representative in V/ker(sum-per-translation-orbit)."""
    output: dict[Word, int] = {}
    for (x_sites, z_sites), coefficient in terms.items():
        key = canonical(x_sites, z_sites, dim)
        value = output.get(key, 0) + coefficient
        if value:
            output[key] = value
        else:
            output.pop(key, None)
    return output


def commutator_coinvariant(
    x_sites: frozenset[Site], z_sites: frozenset[Site], dim: int
) -> dict[Word, int]:
    """Independent calculation of (1/2)[H, X^a Z^b] at a=b=1."""
    output: dict[Word, int] = {}

    def add(key: Word, coefficient: int) -> None:
        value = output.get(key, 0) + coefficient
        if value:
            output[key] = value
        else:
            output.pop(key, None)

    for site in z_sites:
        add(canonical(x_sites ^ {site}, z_sites, dim), 1)

    edges: set[tuple[Site, Site]] = set()
    for site in x_sites:
        for axis in range(dim):
            for sign in (-1, 1):
                neighbor = list(site)
                neighbor[axis] += sign
                edge = tuple(sorted((site, tuple(neighbor))))
                edges.add(edge)
    for left, right in edges:
        if (left in x_sites) ^ (right in x_sites):
            add(canonical(x_sites, z_sites ^ {left, right}, dim), -1)
    return output


def rank_q(columns: list[dict[Word, int]]) -> int:
    """Exact Fraction elimination, intentionally unlike the producer's sparse order."""
    pivots: dict[Word, dict[Word, Fraction]] = {}
    for column in columns:
        vector = {row: Fraction(value) for row, value in column.items() if value}
        while vector:
            pivot = min(vector)
            prior = pivots.get(pivot)
            if prior is None:
                scale = vector[pivot]
                pivots[pivot] = {row: value / scale for row, value in vector.items() if value}
                break
            scale = vector[pivot]
            for row, value in prior.items():
                updated = vector.get(row, Fraction(0)) - scale * value
                if updated:
                    vector[row] = updated
                else:
                    vector.pop(row, None)
    return len(pivots)


def independent_ranks(dim: int, radius: int) -> dict[str, int]:
    sites = box(dim, radius)
    canonical_density = set()
    images: list[dict[Word, int]] = []
    for column in range(1 << (2 * len(sites))):
        x_sites, z_sites = decode(column, sites)
        canonical_density.add(canonical(x_sites, z_sites, dim))
        images.append(commutator_coinvariant(x_sites, z_sites, dim))
    rank_commutator = rank_q(images)
    rank_translation = len(canonical_density)
    return {
        "density_basis_size": 1 << (2 * len(sites)),
        "translation_rank": rank_translation,
        "commutator_rank": rank_commutator,
        "nontrivial_quotient_dimension": rank_translation - rank_commutator - 1,
    }

def independent_modular_ranks(dim: int, radius: int, prime: int) -> dict[str, int]:
    """A second, reverse-column-order finite-field rank computation for the large cubic box."""
    sites = box(dim, radius)
    column_count = 1 << (2 * len(sites))
    row_lookup: dict[Word, int] = {}
    frequencies: list[int] = []
    canonical_density: set[Word] = set()

    for column in range(column_count):
        x_sites, z_sites = decode(column, sites)
        canonical_density.add(canonical(x_sites, z_sites, dim))
        for row in commutator_coinvariant(x_sites, z_sites, dim):
            row_id = row_lookup.get(row)
            if row_id is None:
                row_id = len(frequencies)
                row_lookup[row] = row_id
                frequencies.append(0)
            frequencies[row_id] += 1

    order = [
        row_id
        for row_id, _ in sorted(
            enumerate(frequencies), key=lambda item: (item[1], item[0])
        )
    ]
    inverse_order = [0] * len(order)
    for position, row_id in enumerate(order):
        inverse_order[row_id] = position

    pivots: dict[int, dict[int, int]] = {}
    for column in reversed(range(column_count)):
        x_sites, z_sites = decode(column, sites)
        vector = {
            row_lookup[row]: coefficient % prime
            for row, coefficient in commutator_coinvariant(x_sites, z_sites, dim).items()
            if coefficient % prime
        }
        while vector:
            pivot = min(vector, key=inverse_order.__getitem__)
            prior = pivots.get(pivot)
            if prior is None:
                inverse = pow(vector[pivot], prime - 2, prime)
                pivots[pivot] = {
                    row: value * inverse % prime for row, value in vector.items() if value
                }
                break
            scale = vector[pivot]
            for row, value in prior.items():
                updated = (vector.get(row, 0) - scale * value) % prime
                if updated:
                    vector[row] = updated
                else:
                    vector.pop(row, None)

    rank_commutator = len(pivots)
    rank_translation = len(canonical_density)
    return {
        "density_basis_size": column_count,
        "translation_rank": rank_translation,
        "commutator_rank": rank_commutator,
        "nontrivial_quotient_dimension": rank_translation - rank_commutator - 1,
    }


def locate_case(payload: dict[str, object], dimension: int, radius: int) -> dict[str, object]:
    for case in payload["data"]["cases"]:
        if case["dimension"] == dimension and case["radius"] == radius:
            return case
    raise AssertionError(f"missing d={dimension}, R={radius}")


def verify_exact_small_cases(payload: dict[str, object]) -> None:
    # In d=1, R=1 is the first box exposing a class beyond the Hamiltonian
    # density (the oriented energy current); R=2 confirms the growth pattern.
    chain_one = independent_ranks(1, 1)
    assert chain_one["nontrivial_quotient_dimension"] == 2, chain_one
    stored_chain_one = locate_case(payload, 1, 1)
    assert stored_chain_one["arithmetic"] == "exact_Q"
    for key, value in chain_one.items():
        assert stored_chain_one[key] == value, (key, stored_chain_one[key], value)

    chain = independent_ranks(1, 2)
    assert chain["nontrivial_quotient_dimension"] >= 2, chain
    stored_chain = locate_case(payload, 1, 2)
    assert stored_chain["arithmetic"] == "exact_Q"
    for key, value in chain.items():
        assert stored_chain[key] == value, (key, stored_chain[key], value)

    # The first square box independently tests that the transverse bond
    # direction is retained; only the Hamiltonian class survives.
    square = independent_ranks(2, 1)
    assert square["nontrivial_quotient_dimension"] == 1, square
    stored_square = locate_case(payload, 2, 1)
    assert stored_square["arithmetic"] == "exact_Q"
    for key, value in square.items():
        assert stored_square[key] == value, (key, stored_square[key], value)

    # The other decisive sandwich is independently rebuilt over a third prime.
    # This deliberately uses reverse column order and never imports producer
    # code or its stored pivot data.
    square_large = independent_modular_ranks(2, 2, 1_000_000_007)
    assert square_large["nontrivial_quotient_dimension"] == 1, square_large
    stored_square_large = locate_case(payload, 2, 2)
    assert stored_square_large["arithmetic"] == "sandwich_exact_Q"
    for key, value in square_large.items():
        assert stored_square_large[key] == value, (key, stored_square_large[key], value)

    # This independently catches a dimension/boundary-direction mistake in the
    # cubic construction without importing any producer implementation.
    cubic_one_site = independent_ranks(3, 0)
    stored_cubic = locate_case(payload, 3, 0)
    assert stored_cubic["arithmetic"] == "exact_Q"
    for key, value in cubic_one_site.items():
        assert stored_cubic[key] == value, (key, stored_cubic[key], value)

    # This is the decisive cubic R=1 rank, recomputed over a third prime with
    # reverse column order.  It imports neither the experiment nor its pivot
    # certificate.
    cubic_large = independent_modular_ranks(3, 1, 1_000_000_007)
    assert cubic_large["nontrivial_quotient_dimension"] == 1, cubic_large
    stored_large = locate_case(payload, 3, 1)
    assert stored_large["arithmetic"] == "sandwich_exact_Q"
    for key, value in cubic_large.items():
        assert stored_large[key] == value, (key, stored_large[key], value)



def verify_compact_pivot_traces(payload: dict[str, object]) -> None:
    for case in payload["data"]["cases"]:
        active_rows = case.get("finite_divergence_system", {}).get(
            "active_ordered_pauli_output_orbits", 0
        )
        for modular in case.get("modular_rank_lower_bounds", []):
            certificate = modular["certificate"]
            assert "pivot_trace_column_row" not in certificate
            assert certificate["pivot_trace_encoding"] == "base85 little-endian uint32 (column,row) pairs"
            raw = base64.b85decode(certificate["pivot_trace_column_row_base85"].encode("ascii"))
            count = certificate["pivot_trace_pair_count"]
            assert len(raw) == 8 * count == 8 * certificate["rank"]
            assert hashlib.sha256(raw).hexdigest() == certificate["pivot_trace_sha256"]
            for column, row in struct.iter_unpack("<II", raw):
                assert column < case["density_basis_size"]
                assert row < active_rows


def verify_energy_current_density(payload: dict[str, object]) -> None:
    """Evaluate the explicit d=1 energy-current density, not just its label."""
    # j = X_0 Z_0 Z_1 - X_1 Z_0 Z_1 (the scalar -i relative to the Hermitian
    # Y-form is immaterial for commutation and independence).
    j_terms = {
        (frozenset({(0,)}), frozenset({(0,), (1,)})): 1,
        (frozenset({(1,)}), frozenset({(0,), (1,)})): -1,
    }
    h_terms = {
        (frozenset({(0,)}), frozenset()): 1,
        (frozenset(), frozenset({(0,), (1,)})): 1,
    }
    identity_terms = {(frozenset(), frozenset()): 1}

    # (i) The unquotiented local commutator of j is genuinely nonzero: the
    # cancellation is a telescoping orbit effect, not a termwise one.
    raw: dict[tuple[frozenset[Site], frozenset[Site]], int] = {}
    for (x_sites, z_sites), scale in j_terms.items():
        for key, coefficient in local_commutator(x_sites, z_sites, 1).items():
            value = raw.get(key, 0) + scale * coefficient
            if value:
                raw[key] = value
            else:
                raw.pop(key, None)
    assert raw, "energy current must have nonzero unquotiented commutator"

    # (ii) Its coinvariant (orbit-summed) commutator vanishes in d=1.
    assert not combine_coinvariant(j_terms, 1)
    assert not combine_coinvariant(h_terms, 1)
    assert not combine_coinvariant(identity_terms, 1)

    # (iii) Its class is independent of the identity and of the Hamiltonian
    # density in V/ker(orbit sum): the three class vectors have exact rank 3.
    class_vectors = [
        orbit_sum_vector(identity_terms, 1),
        orbit_sum_vector(h_terms, 1),
        orbit_sum_vector(j_terms, 1),
    ]
    assert rank_q(class_vectors) == 3, "I, h, j must span independent classes"

    # (iv) The same density embedded along e_1 in d=2 no longer commutes once
    # the transverse bonds are included.
    j_terms_d2 = {
        (frozenset({(0, 0)}), frozenset({(0, 0), (1, 0)})): 1,
        (frozenset({(1, 0)}), frozenset({(0, 0), (1, 0)})): -1,
    }
    assert combine_coinvariant(j_terms_d2, 2)

    stored_chain = locate_case(payload, 1, 1)
    assert stored_chain["chain_current_projected_commutator_term_count"] == 0
    stored_square = locate_case(payload, 2, 1)
    assert stored_square["embedded_chain_current_projected_commutator_term_count"] > 0


def reconstruct_current_system(dim: int, radius: int) -> dict[str, object]:
    """Independently rebuild the per-orbit rooted-tree current presentation."""
    sites = box(dim, radius)
    nodes_by_orbit: dict[Word, set[Site]] = {}
    for column in range(1 << (2 * len(sites))):
        x_sites, z_sites = decode(column, sites)
        for next_x, next_z in local_commutator(x_sites, z_sites, dim):
            key, shift = canonical_with_shift(next_x, next_z, dim)
            nodes_by_orbit.setdefault(key, set()).add(shift)

    equation_rows = 0
    current_columns = 0
    direction_columns = [0] * dim
    incidences: list[tuple[int, list[dict[int, int]]]] = []
    for nodes in nodes_by_orbit.values():
        root = tuple(min(node[axis] for node in nodes) for axis in range(dim))
        closure = set(nodes)
        for node in tuple(nodes):
            cursor = list(node)
            while tuple(cursor) != root:
                axis = next(axis for axis in range(dim) if cursor[axis] > root[axis])
                cursor[axis] -= 1
                closure.add(tuple(cursor))
        index = {node: position for position, node in enumerate(sorted(closure))}
        columns: list[dict[int, int]] = []
        for node in sorted(closure):
            if node == root:
                continue
            axis = next(axis for axis in range(dim) if node[axis] > root[axis])
            parent = list(node)
            parent[axis] -= 1
            direction_columns[axis] += 1
            columns.append({index[node]: 1, index[tuple(parent)]: -1})
        equation_rows += len(closure)
        current_columns += len(columns)
        incidences.append((len(closure), columns))
    return {
        "equation_count": equation_rows,
        "current_unknown_count": current_columns,
        "current_unknowns_by_direction": direction_columns,
        "incidences": incidences,
    }


def verify_current_incidence(payload: dict[str, object]) -> None:
    """Tree incidence spans exactly the zero-orbit-sum subspace, per orbit."""
    for dim, radius in ((1, 1), (2, 1)):
        rebuilt = reconstruct_current_system(dim, radius)
        stored = locate_case(payload, dim, radius)["finite_divergence_system"]
        assert rebuilt["equation_count"] == stored["equation_count"]
        assert rebuilt["current_unknown_count"] == stored["current_unknown_count"]
        assert rebuilt["current_unknowns_by_direction"] == stored["current_unknowns_by_direction"]
        assert len(stored["current_unknowns_by_direction"]) == dim

        for size, columns in rebuilt["incidences"]:
            if size == 1:
                assert not columns
                continue
            # Every current column is a pure difference (t - tau t).
            assert all(sum(column.values()) == 0 for column in columns)
            # The rooted tree has incidence rank |vertices| - 1, so its column
            # space is exactly the zero-sum subspace of that orbit.
            assert rank_q(columns) == size - 1
            # A nonzero orbit sum is NOT solvable with currents on this orbit.
            augmented = list(columns) + [{0: 1}]
            assert rank_q(augmented) == size
        if dim == 2:
            # Both transverse current directions are genuinely exercised.
            assert all(count > 0 for count in rebuilt["current_unknowns_by_direction"])

def verify_metadata(payload: dict[str, object]) -> None:
    data = payload["data"]
    construction = data["construction"]
    assert construction["ordered_pauli_convention"] == "X^a Z^b"
    assert construction["density_box"] == "B_R={0,...,R}^d"
    assert construction["trivial_quotient"] == "span{identity density} + finite lattice divergences"
    checks = data["internal_checks"]
    assert checks and all(item["passed"] for item in checks), checks
    assert any(case["dimension"] == 1 and case["positive_control"] for case in data["cases"])
    assert any(case["dimension"] == 3 for case in data["cases"])
    assert "-i" in data["one_dimensional_control"]["known_density"]
    assert payload["provenance"]["rss_measurement"] == "process-lifetime ru_maxrss at measurement end; conservative"
    assert all(
        case["rss_measurement"] == "process-lifetime ru_maxrss at measurement end; conservative"
        for case in data["cases"]
    )
    unresolved_scopes = {entry["scope"] for entry in data["unresolved"]}
    assert {"Z R>=10", "Z^2 R>=3", "Z^3 R>=2"}.issubset(unresolved_scopes)
    unresolved_by_scope = {entry["scope"]: entry for entry in data["unresolved"]}
    for scope in ("Z^2 R>=3", "Z^3 R>=2"):
        assert unresolved_by_scope[scope]["classification"] == "input_size_preflight_not_observed_wall"
    chain_wall = locate_case(payload, 1, 10)
    assert chain_wall["arithmetic"] == "resource_wall"
    assert chain_wall["resource_wall"]["phase"] == "prepare_ordered_pauli_system"
    assert chain_wall["resource_wall"]["processed_columns"] < chain_wall["density_basis_size"]


def main() -> None:
    signal.alarm(TIMEOUT_SECONDS)
    assert RESULT.exists(), "producer artifact is missing"
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    verify_metadata(payload)
    verify_compact_pivot_traces(payload)
    verify_exact_small_cases(payload)
    verify_energy_current_density(payload)
    verify_current_incidence(payload)
    print("PASS")


if __name__ == "__main__":
    main()
