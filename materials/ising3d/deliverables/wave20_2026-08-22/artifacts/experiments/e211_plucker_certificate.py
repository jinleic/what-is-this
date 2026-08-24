#!/usr/bin/env python3
"""Produce the exact Pluecker/Hirota boundary-tensor certificate.

[COMPUTATION] A 2x3 planar rectangle selects one coefficient vector from a
predeclared 125-element class; a fresh 3x3 planar rectangle checks it without
participating in selection.  The same vector is then evaluated on an elementary
cube and a 2x2x3 slab using exact repository lattice tensors.

[THEOREM] The emitted nonzero residuals are finite counterexamples to this one
Hadamard-basis four-channel matchgate identity on the declared terminal sets.

[UNRESOLVED] The result makes no claim against arbitrary determinant/Pfaffian
representations, other basis changes or terminal orders, higher identities, or
an all-size solution of the simple-cubic Ising model.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from ising.lattices import Lattice, hyperrect  # noqa: E402

from e209_plucker_tensor import (  # noqa: E402
    Budget,
    CPU_BUDGET_SECONDS,
    PLUCKER_COEFFICIENTS,
    RSS_CAP_BYTES,
    boundary_partition_tensor,
    edge_subset_currents,
    lower_priority,
    max_rss_bytes,
    monomial,
    plucker_residual,
    poly_at_fraction,
    poly_mul,
    poly_pow,
    walsh_currents,
)
from e210_plucker_search import (  # noqa: E402
    COEFFICIENT_BOUND,
    bounded_search,
    coefficient_candidates,
    first_witness,
    mpoly_digest,
    mpoly_linear_combination,
    multivariate_channels,
    multivariate_currents,
    uniform_specialization,
)

ARTIFACT = ROOT / "results" / "integrability" / "plucker_hirota.json"
V_HALF = Fraction(1, 2)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def first_nonzero(poly: Sequence[int]) -> dict[str, int] | None:
    for degree, coefficient in enumerate(poly):
        if coefficient:
            return {"degree": degree, "coefficient": int(coefficient)}
    return None


def fraction_record(value: Fraction) -> dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
        "nonzero": value != 0,
    }


def coordinates_of_edges(lat: Lattice) -> list[dict[str, object]]:
    return [
        {
            "edge_index": edge_index,
            "left_site": left,
            "right_site": right,
            "left_coordinate": list(lat.coord(left)),
            "right_coordinate": list(lat.coord(right)),
        }
        for edge_index, (left, right) in enumerate(lat.bonds)
    ]


def build_tensor_case(
    name: str,
    shape: Sequence[int],
    terminal_coordinates: Sequence[Sequence[int]],
    budget: Budget,
) -> tuple[dict[str, object], Lattice, tuple[int, ...]]:
    lat = hyperrect(tuple(shape), periodic=False)
    coordinates = tuple(tuple(int(value) for value in coord) for coord in terminal_coordinates)
    terminals = tuple(lat.index(coord) for coord in coordinates)
    tensor = boundary_partition_tensor(lat, terminals, budget)
    walsh = walsh_currents(tensor, lat.n_sites)
    direct = edge_subset_currents(lat, terminals, budget)
    residual = plucker_residual(walsh)
    odd_currents_zero = all(walsh[mask] == [0] for mask in range(16) if mask.bit_count() & 1)
    global_flip_symmetric = all(tensor[mask] == tensor[mask ^ 0b1111] for mask in range(16))
    expected_v0 = 1 << (lat.n_sites - len(terminals))
    v0_normalized = all(poly[0] == expected_v0 for poly in tensor)
    return (
        {
            "tag": "[COMPUTATION]",
            "name": name,
            "shape": list(lat.shape),
            "periodic": list(lat.periodic),
            "lattice_description": lat.describe(),
            "site_count": lat.n_sites,
            "bond_count": lat.n_bonds,
            "terminal_order_coordinates": [list(coord) for coord in coordinates],
            "terminal_order_site_indices": list(terminals),
            "boundary_bit_convention": "bit 1 is spin +1, matching ising.lattices repository convention",
            "tensor_definition": (
                "T_b(v)=sum over interior spins product_edges(1+v*sigma_u*sigma_w); "
                "physical fixed-boundary Z_b=(cosh K)^|E| T_b(tanh K)"
            ),
            "boundary_partition_tensor_coefficients_ascending": tensor,
            "normalized_walsh_current_coefficients_ascending": walsh,
            "walsh_equals_independent_edge_parity_dp": walsh == direct,
            "global_spin_flip_symmetric": global_flip_symmetric,
            "odd_current_entries_zero": odd_currents_zero,
            "all_tensor_entries_at_v0": [poly[0] for poly in tensor],
            "expected_tensor_entry_at_v0": expected_v0,
            "v0_normalization_passed": v0_normalized,
            "plucker_coefficients": list(PLUCKER_COEFFICIENTS),
            "plucker_residual_coefficients_ascending": residual,
            "plucker_residual_first_nonzero": first_nonzero(residual),
            "plucker_residual_at_v_half": fraction_record(poly_at_fraction(residual, V_HALF)),
        },
        lat,
        terminals,
    )


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def coefficient_histogram(poly: dict[int, int]) -> dict[str, int]:
    counts = Counter(poly.values())
    return {str(coefficient): counts[coefficient] for coefficient in sorted(counts)}


def run_all() -> tuple[dict[str, object], list[dict[str, object]], Budget, int | None]:
    niceness = lower_priority()
    budget = Budget(CPU_BUDGET_SECONDS, RSS_CAP_BYTES)
    checks: list[dict[str, object]] = []

    training, training_lat, training_terminals = build_tensor_case(
        "planar_training_2x3",
        (2, 3),
        ((0, 0), (0, 2), (1, 2), (1, 0)),
        budget,
    )
    planar_holdout, holdout_lat, holdout_terminals = build_tensor_case(
        "planar_holdout_3x3",
        (3, 3),
        ((0, 0), (0, 2), (2, 2), (2, 0)),
        budget,
    )
    cube_cofacial, _, _ = build_tensor_case(
        "cube_cofacial_control_2x2x2",
        (2, 2, 2),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
        budget,
    )
    cube_skew, cube_lat, cube_terminals = build_tensor_case(
        "cube_skew_certificate_2x2x2",
        (2, 2, 2),
        ((0, 0, 0), (0, 0, 1), (1, 1, 1), (1, 1, 0)),
        budget,
    )
    slab, _, _ = build_tensor_case(
        "slab_surface_holdout_2x2x3",
        (2, 2, 3),
        ((0, 0, 0), (0, 0, 2), (0, 1, 2), (0, 1, 0)),
        budget,
    )
    cases = (training, planar_holdout, cube_cofacial, cube_skew, slab)

    for case in cases:
        add_check(
            checks,
            f"[COMPUTATION] {case['name']}: Walsh/current agreement",
            bool(case["walsh_equals_independent_edge_parity_dp"]),
            "[LEMMA] fixed-spin enumeration and edge-parity DP give identical integer polynomials",
        )
        add_check(
            checks,
            f"[COMPUTATION] {case['name']}: tensor invariants",
            bool(case["global_spin_flip_symmetric"])
            and bool(case["odd_current_entries_zero"])
            and bool(case["v0_normalization_passed"]),
            (
                f"[LEMMA] global flip holds, odd currents vanish, and each v=0 entry is "
                f"{case['expected_tensor_entry_at_v0']}"
            ),
        )

    training_multi = multivariate_currents(training_lat, training_terminals, budget)
    training_channels = multivariate_channels(training_multi)
    search = bounded_search([training_channels], COEFFICIENT_BOUND)
    selected = [tuple(row) for row in search["survivors"]]

    holdout_multi = multivariate_currents(holdout_lat, holdout_terminals, budget)
    holdout_multi_residual = mpoly_linear_combination(
        multivariate_channels(holdout_multi), PLUCKER_COEFFICIENTS
    )
    cube_multi = multivariate_currents(cube_lat, cube_terminals, budget)
    cube_channels = multivariate_channels(cube_multi)
    cube_multi_residual = mpoly_linear_combination(cube_channels, PLUCKER_COEFFICIENTS)
    joint_survivors = [
        coefficients
        for coefficients in coefficient_candidates(COEFFICIENT_BOUND)
        if not mpoly_linear_combination(training_channels, coefficients)
        and not mpoly_linear_combination(cube_channels, coefficients)
    ]

    cube_uniform_residual = cube_skew["plucker_residual_coefficients_ascending"]
    cube_factor = poly_mul(monomial(6, 4), poly_pow([1, 0, -1], 6))
    slab_uniform_residual = slab["plucker_residual_coefficients_ascending"]
    slab_factor = poly_mul(
        poly_mul(
            poly_mul(monomial(8, 4), poly_pow([1, 0, -1], 8)),
            [1, 0, 1],
        ),
        [3, 0, 6, 0, -1],
    )
    cube_multivariate_uniform = uniform_specialization(cube_multi_residual, cube_lat.n_bonds)
    witness = first_witness(cube_multi_residual, cube_lat.n_bonds)
    if witness is None:
        witness_edges: list[dict[str, object]] = []
    else:
        edge_records = coordinates_of_edges(cube_lat)
        witness_edges = [
            {**edge_records[index], "exponent": exponent}
            for index, exponent in enumerate(witness["edge_exponents"])
            if exponent
        ]

    expected_candidate_count = (2 * COEFFICIENT_BOUND + 1) ** 3
    add_check(
        checks,
        "[COMPUTATION] bounded planar coefficient search",
        search["candidate_count"] == expected_candidate_count
        and selected == [PLUCKER_COEFFICIENTS],
        (
            f"[COMPUTATION] {search['candidate_count']} normalized vectors searched; "
            f"survivors={search['survivors']}"
        ),
    )
    add_check(
        checks,
        "[COMPUTATION] fresh planar holdout",
        planar_holdout["plucker_residual_coefficients_ascending"] == [0]
        and not holdout_multi_residual,
        "[COMPUTATION] both uniform and independent-edge 3x3 residuals vanish exactly",
    )
    add_check(
        checks,
        "[COMPUTATION] cofacial cube scope control",
        cube_cofacial["plucker_residual_coefficients_ascending"] == [0],
        "[COMPUTATION] the cyclic terminals on one elementary-cube face still pass",
    )
    add_check(
        checks,
        "[THEOREM] skew elementary-cube residual certificate",
        cube_uniform_residual == cube_factor
        and poly_at_fraction(cube_uniform_residual, V_HALF) != 0,
        (
            f"[THEOREM] R_cube(v)=4*v^6*(1-v^2)^6; "
            f"R_cube(1/2)={poly_at_fraction(cube_uniform_residual, V_HALF)}"
        ),
    )
    add_check(
        checks,
        "[COMPUTATION] multivariate cube certificate",
        bool(cube_multi_residual)
        and cube_multivariate_uniform == cube_uniform_residual
        and witness is not None
        and int(witness["coefficient"]) != 0,
        (
            f"[COMPUTATION] {len(cube_multi_residual)} nonzero monomials, "
            f"digest={mpoly_digest(cube_multi_residual)}, witness={witness}"
        ),
    )
    add_check(
        checks,
        "[THEOREM] 2x2x3 surface-slab holdout certificate",
        slab_uniform_residual == slab_factor
        and poly_at_fraction(slab_uniform_residual, V_HALF) != 0,
        (
            "[THEOREM] R_slab(v)=4*v^8*(1-v^2)^8*(1+v^2)*(3+6*v^2-v^4); "
            f"R_slab(1/2)={poly_at_fraction(slab_uniform_residual, V_HALF)}"
        ),
    )
    add_check(
        checks,
        "[THEOREM] bounded-class planar/cube no-go",
        not joint_survivors,
        (
            f"[THEOREM] no vector in the declared {expected_candidate_count}-element class "
            "passes both the 2x3 planar training tensor and the skew cube"
        ),
    )

    budget.check("final resource audit")
    add_check(
        checks,
        "[COMPUTATION] process-time and RSS budgets",
        budget.used() < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        (
            f"[COMPUTATION] CPU={budget.used():.3f}s<{CPU_BUDGET_SECONDS}; "
            f"RSS={max_rss_bytes()}<{RSS_CAP_BYTES}"
        ),
    )

    search_record = {
        **search,
        "training_only": "planar_training_2x3 with independent edge variables",
        "not_used_for_selection": [
            "planar_holdout_3x3",
            "cube_cofacial_control_2x2x2",
            "cube_skew_certificate_2x2x2",
            "slab_surface_holdout_2x2x3",
        ],
        "channel_order": [
            "C_empty*C_1234",
            "C_12*C_34",
            "C_13*C_24",
            "C_14*C_23",
        ],
        "joint_planar_training_and_cube_survivors": [list(row) for row in joint_survivors],
        "scope": (
            "[UNRESOLVED] exhaustive only for the 125 normalized integer coefficient vectors; "
            "not a search over arbitrary determinant or Pfaffian representations"
        ),
    }
    cube_skew["exact_factorization"] = "4*v^6*(1-v^2)^6"
    slab["exact_factorization"] = (
        "4*v^8*(1-v^2)^8*(1+v^2)*(3+6*v^2-v^4)"
    )
    cube_multivariate_record = {
        "tag": "[COMPUTATION]",
        "edge_variable_order": coordinates_of_edges(cube_lat),
        "nonzero_term_count": len(cube_multi_residual),
        "coefficient_histogram": coefficient_histogram(cube_multi_residual),
        "canonical_sha256": mpoly_digest(cube_multi_residual),
        "first_minimum_degree_witness": witness,
        "witness_nonzero_edges": witness_edges,
        "uniform_specialization_coefficients_ascending": cube_multivariate_uniform,
    }

    data = {
        "identity": {
            "tag": "[LEMMA]",
            "name": "four-terminal Pfaffian Pluecker relation in the terminal Walsh basis",
            "formula": (
                "C_empty*C_1234-C_12*C_34+C_13*C_24-C_14*C_23=0"
            ),
            "coefficient_origin": (
                "[LEMMA] Pf([[0,a12,a13,a14],...])="
                "a12*a34-a13*a24+a14*a23; no 3D value fixes a sign"
            ),
            "basis": (
                "C_S=2^(-|V|) sum_b chi_S(b) T_b, with bit 1 denoting spin +1"
            ),
        },
        "bounded_identity_search": search_record,
        "planar_training": training,
        "planar_holdout": planar_holdout,
        "cofacial_cube_scope_control": cube_cofacial,
        "skew_cube_certificate": cube_skew,
        "skew_cube_multivariate_certificate": cube_multivariate_record,
        "surface_slab_holdout": slab,
        "conclusion": {
            "exact_finite_statement": (
                "[THEOREM] The standard four-channel Pfaffian Pluecker relation, with signs "
                "fixed before all 3D evaluations, has residual 4*v^6*(1-v^2)^6 on the "
                "declared skew-terminal elementary cube and the displayed nonzero residual on "
                "the declared 2x2x3 surface slab."
            ),
            "physical_interval": (
                "[THEOREM] Both displayed factorizations are positive for 0<v<1; for the slab, "
                "3+6*v^2-v^4>0 on that interval."
            ),
            "bounded_no_go": (
                "[THEOREM] Within the explicitly searched 125 normalized coefficient vectors, "
                "none passes both the multivariate 2x3 planar control and the skew cube."
            ),
            "counterexample_to_dimension_only_overstatement": (
                "[COMPUTATION] A cofacial cyclic terminal order on the same 2x2x2 cube has zero "
                "residual, so merely calling a graph three-dimensional does not force failure."
            ),
            "scope": (
                "[UNRESOLVED] These finite certificates do not exclude arbitrary determinant/"
                "Pfaffian representations, other local transforms, other terminal orders, "
                "higher Pluecker/Dodgson/octahedron systems, or any all-size identity."
            ),
        },
    }
    return data, checks, budget, niceness


def make_artifact(
    data: dict[str, object],
    checks: list[dict[str, object]],
    budget: Budget,
    niceness: int | None,
) -> dict[str, object]:
    sources = [
        ROOT / "experiments" / "e209_plucker_tensor.py",
        ROOT / "experiments" / "e210_plucker_search.py",
        Path(__file__),
    ]
    return {
        "meta": {
            "experiment": "e211_plucker_certificate",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "provenance": "experiments/e211_plucker_certificate.py",
            "command": "nice -n 10 .venv/bin/python experiments/e211_plucker_certificate.py",
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "process_id": os.getpid(),
            "arithmetic": (
                "exact Python integers and fractions.Fraction; fixed-spin enumeration, exact "
                "Walsh division, edge-parity dynamic programming, and bounded multivariate search"
            ),
            "parallel_processes": 1,
            "observed_niceness": niceness,
            "benchmark_use": "[COMPUTATION] K_c=0.221654626 was not read, selected, or fitted",
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "observed_process_time_seconds": round(budget.used(), 6),
                "observed_peak_rss_bytes": max_rss_bytes(),
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in sources
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks, budget, niceness = run_all()
    artifact = make_artifact(data, checks, budget, niceness)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"[COMPUTATION] wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if all(bool(item["passed"]) for item in checks):
        print("PASS e211 Pluecker/Hirota certificate", flush=True)
        return 0
    print("FAIL e211 Pluecker/Hirota certificate", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
