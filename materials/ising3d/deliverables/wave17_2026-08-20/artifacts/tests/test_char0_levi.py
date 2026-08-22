#!/usr/bin/env python3
"""Standalone exact re-derivation for the 2x3 characteristic-zero Levi result."""
from __future__ import annotations

import importlib.util
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("e45", ROOT / "experiments" / "e45_char0_levi.py")
e45 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(e45)


def check(name, condition, detail=""):
    print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not condition:
        raise AssertionError(name)


def independent_full_orbit_actions():
    """Cross-check the compressed action against an independently assembled map.

    This routine intentionally does not call exact_pauli_closure or its nested
    action: it expands each orbit sum in all 1056 Pauli strings, brackets every
    source with every term of A/B, and checks representative independence.
    """
    xs = tuple(1 << i for i in range(e45.N))
    zz = tuple(((1 << u) | (1 << v)) << e45.N for u, v in e45.BONDS)
    singles = xs + zz
    columns = e45.support_columns(singles)
    permutations = e45.site_permutations()
    seen, orbits = set(), []
    for value in columns:
        if value in seen:
            continue
        orbit = tuple(sorted({e45.permute_pauli(value, p) for p in permutations}))
        seen.update(orbit)
        orbits.append(orbit)
    index = {value: i for i, orbit in enumerate(orbits) for value in orbit}
    maps = []
    for group in (xs, zz):
        action = []
        for source_orbit in orbits:
            coefficients = [dict() for _ in orbits]
            for source in source_orbit:
                for generator in group:
                    sign = e45.bracket_sign(generator, source)
                    if sign:
                        target = source ^ generator
                        coefficients[index[target]][target] = coefficients[index[target]].get(target, 0) + sign
            row = []
            for target_orbit, values in zip(orbits, coefficients):
                coefficient_set = {values.get(target, 0) for target in target_orbit}
                if len(coefficient_set) != 1:
                    raise AssertionError("orbit coefficient is representative-dependent")
                row.append(coefficient_set.pop())
            action.append(row)
        maps.append(action)
    return orbits, maps


def independent_pauli_dimension(orbits, maps):
    xs = tuple(1 << i for i in range(e45.N))
    zz = tuple(((1 << u) | (1 << v)) << e45.N for u, v in e45.BONDS)
    index = {value: i for i, orbit in enumerate(orbits) for value in orbit}

    def seed(group):
        vector = [0] * len(orbits)
        for value in group:
            vector[index[value]] += 1
        return vector

    def apply(mapping, vector):
        return [sum(vector[source] * mapping[source][target] for source in range(len(orbits))) for target in range(len(orbits))]

    rows, pivots, frontier = [], {}, []
    for vector in (seed(xs), seed(zz)):
        assert e45.add_fraction_row(rows, pivots, vector)
        frontier.append(vector)
    head = 0
    while head < len(frontier):
        vector = frontier[head]
        head += 1
        for mapping in maps:
            image = apply(mapping, vector)
            if e45.add_fraction_row(rows, pivots, image):
                frontier.append(image)
    return len(rows)


def main():
    orbits, independent_maps = independent_full_orbit_actions()
    production_groups, production_orbits, production_maps = e45.orbit_pauli_actions()
    check(
        "full orbit action maps agree entrywise",
        tuple(tuple(sorted(orbit)) for orbit in production_orbits) == tuple(orbits)
        and all(independent == production for independent, production in zip(independent_maps, production_maps)),
    )
    independent_dimension = independent_pauli_dimension(orbits, independent_maps)
    check("independent full-orbit Pauli closure dimension", independent_dimension == 263)

    # Rebuild everything rather than reading the result JSON.  This independently
    # exercises exact QQ ranks, the full Killing matrix, radical, and factors.
    artifact = e45.build_structure()
    data = artifact["data"]
    check("exact dim_Q g", data["dimension_Q"] == 263)
    check("faithful matrix dimension", data["faithful_matrix_dimension_Q"] == independent_dimension)
    check("Killing form is 263x263 symmetric integer", data["Killing"]["size"] == [263, 263] and data["Killing"]["symmetric"] and data["Killing"]["entry_domain"].startswith("Z"))
    check("exact Killing rank", data["Killing"]["rank_Q"] == 262)
    check("exact radical dimension", data["solvable_radical"]["dimension_Q"] == 1 and data["solvable_radical"]["equals_center"])
    check("exact derived dimension", data["derived_dimension_Q"] == 262)
    factors = data["semisimple_quotient"]["factors"]
    check("Levi factor dimensions", [x["dimension"] for x in factors] == [105, 21, 21, 80, 35] and data["semisimple_quotient"]["selected_joint_projection_rank_Q"] == 262)
    largest = factors[0]
    check("largest factor dimension", largest["dimension"] == 105)
    check("largest factor simplicity certificate", largest["type"] == "C7" and largest["certificate"]["image_dimension"] == largest["certificate"]["ambient_classical_dimension"] == 105 and largest["certificate"]["invariant_form_determinant"] == 16384)
    check("B7 versus C7", data["largest_factor"]["rank"] == 7 and "14-dimensional" in data["largest_factor"]["B7_vs_C7_invariant"])
    check("all generated checks", all(item["passed"] for item in artifact["checks"]))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
