#!/usr/bin/env python3
"""Exact support-resolved commutants for finite Ising layer generator pairs.

For A=sum_v X_v and B=sum_{uv in E} Z_u Z_v, enumerate connected site sets
S and solve exactly for operators O=O_S tensor I_(V\\S) with [O,A]=[O,B]=0.
The main scan uses two prime fields.  A prime nullity of one is a rigorous
certificate over Q because I is an explicit rational kernel vector.  Every
nontrivial kernel is rebuilt by sparse Fraction elimination and substituted
in the integer equations.
"""
from __future__ import annotations

import itertools
import json
import platform
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "flux_support.json"
SCRIPT = "experiments/e57_flux_support.py"
PRIMES = (2_147_483_647, 2_147_483_629)

SparseRow = dict[int, int]
RationalVector = dict[int, Fraction]


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        [(r * columns + c, (r + 1) * columns + c) for r in range(rows - 1) for c in range(columns)]
        + [(r * columns + c, r * columns + c + 1) for r in range(rows) for c in range(columns - 1)]
    )


def chain_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, site + 1) for site in range(n - 1))


def cycle_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, (site + 1) % n) for site in range(n))


def graph_automorphisms(rows: int, columns: int) -> tuple[tuple[int, ...], ...]:
    maps: list[tuple[int, ...]] = []
    if rows == columns:
        for transpose in (False, True):
            for flip_rows in (False, True):
                for flip_columns in (False, True):
                    image = []
                    for row in range(rows):
                        for column in range(columns):
                            r, c = (column, row) if transpose else (row, column)
                            if flip_rows:
                                r = rows - 1 - r
                            if flip_columns:
                                c = columns - 1 - c
                            image.append(r * columns + c)
                    maps.append(tuple(image))
    else:
        for flip_rows in (False, True):
            for flip_columns in (False, True):
                maps.append(
                    tuple(
                        (rows - 1 - r if flip_rows else r) * columns
                        + (columns - 1 - c if flip_columns else c)
                        for r in range(rows)
                        for c in range(columns)
                    )
                )
    return tuple(dict.fromkeys(maps))


def transform_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[site] for site in range(len(permutation)) if (mask >> site) & 1)


def is_connected(mask: int, n: int, edges: tuple[tuple[int, int], ...]) -> bool:
    if not mask:
        return False
    seen = mask & -mask
    while True:
        old = seen
        for left, right in edges:
            if ((seen >> left) & 1) and ((mask >> right) & 1):
                seen |= 1 << right
            if ((seen >> right) & 1) and ((mask >> left) & 1):
                seen |= 1 << left
        if seen == old:
            return seen == mask


def connected_support_representatives(
    rows: int, columns: int, max_support: int
) -> tuple[list[tuple[int, int]], tuple[tuple[int, ...], ...]]:
    n = rows * columns
    edges = grid_edges(rows, columns)
    automorphisms = graph_automorphisms(rows, columns)
    representatives: list[tuple[int, int]] = []
    for size in range(2, max_support + 1):
        for sites in itertools.combinations(range(n), size):
            mask = sum(1 << site for site in sites)
            if not is_connected(mask, n, edges):
                continue
            orbit = {transform_mask(mask, permutation) for permutation in automorphisms}
            if mask == min(orbit):
                representatives.append((mask, len(orbit)))
    return representatives, automorphisms


def local_commutator_rows(
    n: int, edges: tuple[tuple[int, int], ...], support_mask: int
) -> tuple[list[SparseRow], dict[str, int | list[int]]]:
    """Rows of (ad_A/2, ad_B/2) restricted to End(H_S).

    Local columns are ordered monomials Q_(a,b)=X^a Z^b on S.  Output rows
    retain global Pauli masks, so boundary-bond outputs carrying an outside Z
    can never be cancelled accidentally by an internal output.
    """
    sites = [site for site in range(n) if (support_mask >> site) & 1]
    dimension = 1 << (2 * len(sites))
    output: defaultdict[tuple[int, int], SparseRow] = defaultdict(dict)
    incident_edges = [
        edge for edge in edges if ((support_mask >> edge[0]) & 1) or ((support_mask >> edge[1]) & 1)
    ]
    for column in range(dimension):
        digits = column
        a_mask = 0
        b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a_mask |= 1 << site
            if label & 2:
                b_mask |= 1 << site
        pauli = a_mask | (b_mask << n)
        for site in sites:
            if (b_mask >> site) & 1:
                target = pauli ^ (1 << site)
                output[(0, target)][column] = 1
        for left, right in incident_edges:
            if ((a_mask >> left) ^ (a_mask >> right)) & 1:
                target = pauli ^ (1 << (n + left)) ^ (1 << (n + right))
                output[(1, target)][column] = -1
    rows = [row for _, row in sorted(output.items()) if row]
    return rows, {
        "sites": sites,
        "support_mask": support_mask,
        "unknowns": dimension,
        "nonzero_equations": len(rows),
        "nonzero_entries": sum(len(row) for row in rows),
        "incident_bonds": len(incident_edges),
        "boundary_bonds": sum(((support_mask >> u) & 1) ^ ((support_mask >> v) & 1) for u, v in edges),
    }


def local_commutator_columns(
    n: int, edges: tuple[tuple[int, int], ...], support_mask: int
) -> tuple[list[SparseRow], dict[str, int | list[int]]]:
    """Sparse columns of the same stacked map, for fast modular elimination."""
    sites = [site for site in range(n) if (support_mask >> site) & 1]
    dimension = 1 << (2 * len(sites))
    incident_edges = [
        edge for edge in edges if ((support_mask >> edge[0]) & 1) or ((support_mask >> edge[1]) & 1)
    ]
    columns: list[SparseRow] = []
    nonzero_entries = 0
    target_rows: set[tuple[int, int]] = set()
    for column in range(dimension):
        digits = column
        a_mask = 0
        b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a_mask |= 1 << site
            if label & 2:
                b_mask |= 1 << site
        pauli = a_mask | (b_mask << n)
        vector: SparseRow = {}
        for site in sites:
            if (b_mask >> site) & 1:
                target = (0, pauli ^ (1 << site))
                vector[target] = 1
                target_rows.add(target)
        for left, right in incident_edges:
            if ((a_mask >> left) ^ ((a_mask >> right) & 1)) & 1:
                target = (1, pauli ^ (1 << (n + left)) ^ (1 << (n + right)))
                vector[target] = -1
                target_rows.add(target)
        nonzero_entries += len(vector)
        columns.append(vector)
    return columns, {
        "sites": sites,
        "support_mask": support_mask,
        "unknowns": dimension,
        "nonzero_equations": len(target_rows),
        "nonzero_entries": nonzero_entries,
        "incident_bonds": len(incident_edges),
        "boundary_bonds": sum(((support_mask >> u) & 1) ^ ((support_mask >> v) & 1) for u, v in edges),
    }


def column_rank_mod(columns: Iterable[SparseRow], prime: int) -> int:
    """Sparse rare-output-first column rank over F_p."""
    materialized = list(columns)
    frequencies = Counter(row for column in materialized for row in column)
    row_order = {
        row: index for index, row in enumerate(sorted(frequencies, key=lambda row: (frequencies[row], row)))
    }
    pivots: dict[tuple[int, int], dict[tuple[int, int], int]] = {}
    rank = 0
    for source in sorted(materialized, key=len):
        vector = {row: value % prime for row, value in source.items() if value % prime}
        while vector:
            lead = min(vector, key=row_order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = {row: value * inverse % prime for row, value in vector.items()}
                rank += 1
                break
            coefficient = vector[lead]
            for row, value in old.items():
                reduced = (vector.get(row, 0) - coefficient * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return rank


def rank_mod(rows: Iterable[SparseRow], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for source in rows:
        row = {column: value % prime for column, value in source.items() if value % prime}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(row[lead], prime - 2, prime)
                pivots[lead] = {column: value * inverse % prime for column, value in row.items()}
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = (row.get(column, 0) - coefficient * value) % prime
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return len(pivots)


def echelon_q(rows: Iterable[SparseRow]) -> dict[int, dict[int, Fraction]]:
    pivots: dict[int, dict[int, Fraction]] = {}
    for source in rows:
        row = {column: Fraction(value) for column, value in source.items() if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                coefficient = row[lead]
                pivots[lead] = {column: value / coefficient for column, value in row.items()}
                break
            coefficient = row[lead]
            for column, value in old.items():
                reduced = row.get(column, Fraction(0)) - coefficient * value
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return pivots


def nullspace_q(pivots: dict[int, dict[int, Fraction]], dimension: int) -> list[RationalVector]:
    basis: list[RationalVector] = []
    for free in range(dimension):
        if free in pivots:
            continue
        vector: RationalVector = {free: Fraction(1)}
        for pivot in sorted(pivots, reverse=True):
            value = sum(
                coefficient * vector.get(column, Fraction(0))
                for column, coefficient in pivots[pivot].items()
                if column != pivot
            )
            if value:
                vector[pivot] = -value
        basis.append(vector)
    return basis


def verify_kernel(rows: list[SparseRow], basis: list[RationalVector]) -> bool:
    return all(
        sum(Fraction(coefficient) * vector.get(column, Fraction(0)) for column, coefficient in row.items()) == 0
        for vector in basis
        for row in rows
    )


def vector_rank_q(vectors: list[RationalVector]) -> int:
    return len(echelon_q(vectors))


def encode_vector(vector: RationalVector, sites: list[int], n: int) -> list[list[int]]:
    terms: list[list[int]] = []
    for column, coefficient in sorted(vector.items()):
        digits = column
        a_mask = 0
        b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a_mask |= 1 << site
            if label & 2:
                b_mask |= 1 << site
        terms.append([a_mask, b_mask, coefficient.numerator, coefficient.denominator])
    return terms


def scan_layer(rows: int, columns: int, max_support: int) -> dict[str, object]:
    started = time.monotonic()
    n = rows * columns
    edges = grid_edges(rows, columns)
    representatives, automorphisms = connected_support_representatives(rows, columns, max_support)
    table: list[dict[str, object]] = []
    exceptional: list[dict[str, object]] = []
    for support_mask, orbit_size in representatives:
        matrix_columns, metadata = local_commutator_columns(n, edges, support_mask)
        dimension = int(metadata["unknowns"])
        modular = []
        for prime in PRIMES:
            rank = column_rank_mod(matrix_columns, prime)
            modular.append({"prime": prime, "rank": rank, "nullity": dimension - rank})
        nullities = [entry["nullity"] for entry in modular]
        record: dict[str, object] = {
            "support_mask": support_mask,
            "sites": metadata["sites"],
            "size": support_mask.bit_count(),
            "symmetry_orbit_size": orbit_size,
            "unknowns": dimension,
            "nonzero_equations": metadata["nonzero_equations"],
            "nonzero_entries": metadata["nonzero_entries"],
            "incident_bonds": metadata["incident_bonds"],
            "boundary_bonds": metadata["boundary_bonds"],
            "modular": modular,
            "kernel_dimension_Q": 1 if nullities[0] == 1 else None,
            "rational_certificate": (
                "nullity_Fp=1 and explicit I imply nullity_Q=1"
                if nullities[0] == 1
                else "rebuilt by exact Fraction elimination"
            ),
        }
        if len(set(nullities)) != 1:
            raise ArithmeticError(f"prime nullities disagree for {rows}x{columns}, mask {support_mask}: {nullities}")
        matrix_rows: list[SparseRow] | None = None
        if nullities[0] > 1:
            matrix_rows, row_metadata = local_commutator_rows(n, edges, support_mask)
            if any(row_metadata[key] != metadata[key] for key in ("unknowns", "nonzero_equations", "nonzero_entries")):
                raise ArithmeticError(f"row/column map metadata disagree for support mask {support_mask}")
            pivots = echelon_q(matrix_rows)
            basis = nullspace_q(pivots, dimension)
            if not verify_kernel(matrix_rows, basis) or vector_rank_q(basis) != len(basis):
                raise ArithmeticError(f"invalid exact kernel for support mask {support_mask}")
            record["kernel_dimension_Q"] = len(basis)
            record["exact_basis_support_sizes"] = [len(vector) for vector in basis]
            record["exact_basis_coefficient_set"] = sorted({str(value) for vector in basis for value in vector.values()})
            record["exact_basis"] = [
                {"name": "identity" if index == 0 else f"local_charge_{index}", "terms_a_b_num_den": encode_vector(vector, list(metadata["sites"]), n)}
                for index, vector in enumerate(basis)
            ]
            record["exact_substitution_passed"] = True
            record["basis_exact_rank_Q"] = len(basis)
            record["parent_review_flag"] = "LOCAL CONSERVED CHARGE: non-scalar exact kernel"
            exceptional.append(record)
        table.append(record)
    by_size = []
    for size in range(2, max_support + 1):
        rows_at_size = [entry for entry in table if entry["size"] == size]
        expanded = sum(int(entry["symmetry_orbit_size"]) for entry in rows_at_size)
        dimensions = Counter(int(entry["kernel_dimension_Q"]) for entry in rows_at_size)
        by_size.append(
            {
                "size": size,
                "orbit_representatives": len(rows_at_size),
                "all_connected_supports": expanded,
                "kernel_dimension_histogram_on_orbit_representatives": {str(key): dimensions[key] for key in sorted(dimensions)},
            }
        )
    scalar_through = max(
        size
        for size in range(1, max_support + 1)
        if all(int(entry["kernel_dimension_Q"]) == 1 for entry in table if int(entry["size"]) <= size)
    )
    return {
        "claim_tag": "[THEOREM]",
        "graph": f"{rows}x{columns}_open_grid",
        "shape": [rows, columns],
        "n": n,
        "edges": [list(edge) for edge in edges],
        "support_budget": max_support,
        "automorphism_group_size": len(automorphisms),
        "connected_supports_are_quotiented_by_full_graph_automorphisms": True,
        "scalar_only_for_every_connected_support_through": scalar_through,
        "theorem": (
            f"[THEOREM] Every operator supported on a connected set of at most {scalar_through} sites and commuting with A and B is scalar."
        ),
        "extension_outcome": (
            "[COMPUTATION] all tested supports remain scalar-only"
            if not exceptional
            else f"[COMPUTATION] non-scalar kernels occur at {len(exceptional)} support-orbit representatives beyond the theorem cutoff"
        ),
        "by_size": by_size,
        "table": table,
        "exceptional_support_masks": [int(entry["support_mask"]) for entry in exceptional],
        "elapsed_seconds": time.monotonic() - started,
    }


def control_case(name: str, n: int, edges: tuple[tuple[int, int], ...], support_mask: int) -> dict[str, object]:
    rows, metadata = local_commutator_rows(n, edges, support_mask)
    pivots = echelon_q(rows)
    basis = nullspace_q(pivots, int(metadata["unknowns"]))
    modular = [
        {"prime": prime, "nullity": int(metadata["unknowns"]) - rank_mod(rows, prime)} for prime in PRIMES
    ]
    return {
        "claim_tag": "[COMPUTATION]",
        "name": name,
        "n": n,
        "edges": [list(edge) for edge in edges],
        "support_mask": support_mask,
        "sites": metadata["sites"],
        "kernel_dimension_Q": len(basis),
        "modular": modular,
        "exact_basis_support_sizes": [len(vector) for vector in basis],
        "exact_substitution_passed": verify_kernel(rows, basis),
        "basis_exact_rank_Q": vector_rank_q(basis),
        "exact_basis": [
            {"name": "identity" if index == 0 else f"control_charge_{index}", "terms_a_b_num_den": encode_vector(vector, list(metadata["sites"]), n)}
            for index, vector in enumerate(basis)
        ],
    }


def main() -> None:
    layers = [scan_layer(2, 4, 7), scan_layer(3, 3, 7), scan_layer(3, 4, 7)]
    controls = [
        control_case("open_chain_P4_full_support", 4, chain_edges(4), (1 << 4) - 1),
        control_case("cycle_C4_three_site_support", 4, cycle_edges(4), 0b0111),
        control_case("cycle_C4_full_support", 4, cycle_edges(4), (1 << 4) - 1),
    ]
    anomalies = [entry for layer in layers for entry in layer["table"] if int(entry["kernel_dimension_Q"]) > 1]
    checks = [
        {
            "name": "two_prime_agreement",
            "passed": all(
                len({int(row["nullity"]) for row in entry["modular"]}) == 1
                for layer in layers
                for entry in layer["table"]
            ),
            "detail": f"all support rows agree at primes {PRIMES}",
        },
        {
            "name": "layer_theorem_cutoffs",
            "passed": [layer["scalar_only_for_every_connected_support_through"] for layer in layers] == [7, 6, 7],
            "detail": "scalar-only cutoffs are 7 on 2x4, 6 on 3x3, and 7 on 3x4",
        },
        {
            "name": "all_nontrivial_kernels_verified_over_Q",
            "passed": all(entry.get("exact_substitution_passed") and entry.get("basis_exact_rank_Q") == entry["kernel_dimension_Q"] for entry in anomalies),
            "detail": f"{len(anomalies)} exceptional support orbits have exact Fraction bases and zero integer residuals",
        },
        {
            "name": "three_by_three_local_charges_flagged",
            "passed": sorted((entry["support_mask"], entry["kernel_dimension_Q"]) for entry in anomalies) == [(239, 4), (254, 4), (367, 2), (381, 2)],
            "detail": "four seven-site support orbits are non-scalar and require parent review",
        },
        {
            "name": "one_dimensional_free_fermion_control",
            "passed": controls[0]["kernel_dimension_Q"] == 5 and controls[0]["exact_substitution_passed"],
            "detail": "the full P4 joint commutant has exact dimension 5, reproducing the finite free-fermion charge structure",
        },
        {
            "name": "C4_extra_symmetry_control",
            "passed": [controls[1]["kernel_dimension_Q"], controls[2]["kernel_dimension_Q"]] == [3, 27]
            and all(control["exact_substitution_passed"] for control in controls[1:]),
            "detail": "C4 has dimension 3 already on a three-site support and full dimension 27",
        },
    ]
    artifact = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "method": (
                "ordered-Pauli sparse commutator maps; exact ranks over two primes; prime-nullity-one plus explicit identity gives exact Q scalar kernels; every excess kernel rebuilt and substituted with Fraction arithmetic"
            ),
        },
        "data": {
            "claim_tags": ["[LEMMA]", "[THEOREM]", "[COMPUTATION]", "[UNRESOLVED]"],
            "local_condition": {
                "claim_tag": "[LEMMA]",
                "statement": (
                    "[LEMMA] For O=O_S tensor I_(V\\S), [O,A]=[O,A_S] and [O,B] is the sum over internal and boundary bonds incident to S; bonds outside S vanish. Hence [O,A]=[O,B]=0 is exactly the stacked local linear system used here."
                ),
                "ordered_pauli_formula": (
                    "[LEMMA] ad_A(Q_(a,b))/2=sum_{v in S:b_v=1} Q_(a xor e_v,b); ad_B(Q_(a,b))/2=-sum_{uv incident to S:a_u xor a_v=1} Q_(a,b xor e_u xor e_v), with a,b zero outside S"
                ),
                "unknowns": "[COMPUTATION] 4^|S| coefficients over Q; global output masks retain outside Z factors on boundary bonds",
            },
            "transfer_family_scope": (
                "[LEMMA] A coupling-independent operator commutes with the analytic one-parameter family generated by A and B iff it commutes with both generators; the artifact certifies the joint commutant on each support."
            ),
            "layers": layers,
            "controls": controls,
            "predecessor": {
                "claim_tag": "[THEOREM]",
                "reference": "[THEOREM] proofs/algebraic_obstruction.md, Theorem 11 and Lemma 12",
                "relationship": (
                    "[THEOREM] the predecessor excludes individual nontrivial Pauli strings; this finite support calculation allows arbitrary rational linear combinations of all 4^|S| Pauli operators"
                ),
            },
            "scope_and_wall": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "[UNRESOLVED] The scan stops at seven sites. At eight sites each support has 65536 unknowns; 3x4 has 66 symmetry-inequivalent connected eight-site supports. The available sparse column elimination showed severe fill-in (a single full 2x4 modular kernel took about 156 seconds), so no size-eight 3x4 theorem is claimed."
                ),
                "not_claimed": (
                    "[UNRESOLVED] No thermodynamic non-integrability theorem, no disconnected-support theorem, no size-eight-or-larger theorem, and no assertion that the exact seven-site 3x3 charges extend with layer size."
                ),
            },
        },
        "checks": checks,
    }
    if not all(check["passed"] for check in checks):
        raise AssertionError([check for check in checks if not check["passed"]])
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"wrote {RESULT.relative_to(ROOT)}")
    print("PASS")


if __name__ == "__main__":
    main()
