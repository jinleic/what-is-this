"""Bounded exact census of the Kac--Ward branch ``11111`` thin sections.

This producer reuses e139's exact construction machinery and never turns a
process-time cap into a theorem.  The default run performs the exact sparse
setup and an independent anchor control, estimates the exhaustive 4**7 launch,
and writes an honest ``PREFLIGHT`` artifact when the complete census cannot fit.
Only a complete census with both required holdout routes could be decisive;
the characteristic-zero branch question remains unresolved.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import itertools
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/kac_ward/full_thin_census.json"
PRIOR = ROOT / "results/kac_ward/components.json"
E139_PATH = ROOT / "experiments/e139_kw_components.py"
DEFAULT_CAP_SECONDS = 20.0
UNIT_VALUES = (1, 2, 3, 4)
TOTAL_SECTIONS = len(UNIT_VALUES) ** 7
TOTAL_SOLVE_POINTS = len(UNIT_VALUES) ** 6
FREE_NAMES = [
    "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx",
    "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx",
    "u_mz_py", "u_mz_my",
]


def load_e139():
    """Load e139 as a module without running its producer main()."""
    spec = importlib.util.spec_from_file_location("e139_kw_components", E139_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {E139_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_bytes(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def load_prior() -> tuple[dict[str, Any], str]:
    return json.loads(PRIOR.read_text(encoding="utf-8")), sha256_bytes(PRIOR)


def setup_exact(e139):
    """Build e139's exact sparse primitives and six-variable slice data."""
    system = e139.full_weight_finite_system(e139.CONSTRUCTION_SHAPES, e139.BASE_ORDERS, gauge_fix=False)
    variables, equations = e139.branch_11111_equations(system)
    symbol = {str(value): value for value in variables}
    non_diagonal = [symbol[name] for name in e139.NONDIAGONAL_NAMES]
    free = [symbol[name] for name in FREE_NAMES]
    thin = e139.thin_torus_substitutions(symbol)
    diagonal = {symbol[name]: 1 for name in e139.DIAGONAL_NAMES}
    primitives = e139.primitive_numerators(equations, variables, thin, free, diagonal)
    primitive_sparse = [e139.sparse_integer_poly(poly.as_expr(), non_diagonal) for poly in primitives]
    sparse42 = [e139.sparse_integer_poly(expression, variables) for expression in equations]
    return {
        "variables": variables,
        "equations": equations,
        "symbol": symbol,
        "non_diagonal": non_diagonal,
        "thin": thin,
        "primitives": primitives,
        "primitive_sparse": primitive_sparse,
        "sparse42": sparse42,
    }


def section_solutions(e139, primitive_polys: list[sp.Poly], held: list[int]) -> list[list[int]]:
    """Enumerate exactly the 4**6 unit solve grid for one held section."""
    solutions: list[list[int]] = []
    for solve in itertools.product(UNIT_VALUES, repeat=6):
        point = list(solve) + list(held)
        if all(e139.eval_poly_mod(poly, point, e139.PRIME) == 0 for poly in primitive_polys):
            solutions.append(point)
    return solutions


def slice_polynomials(e139, primitives: list[sp.Poly], non_diagonal: list[sp.Symbol], held: list[int]):
    held_map = {non_diagonal[6 + index]: value for index, value in enumerate(held)}
    return [
        sp.Poly(poly.as_expr().subs(held_map), *non_diagonal[:6], domain=sp.ZZ)
        for poly in primitives
    ]


def anchor_control(e139, setup: dict[str, Any], prior: dict[str, Any]) -> dict[str, Any]:
    """Recompute e139's exact unit anchor control and compare its controls."""
    anchor = list(e139.ANCHOR_HELD)
    solutions = section_solutions(e139, setup["primitives"], anchor)
    sliced = slice_polynomials(e139, setup["primitives"], setup["non_diagonal"], anchor)
    determinants: list[int] = []
    lifts: list[dict[str, Any]] = []
    for point in solutions:
        determinant, jacobian = e139.jacobian_determinant_mod_prime(
            sliced, setup["non_diagonal"][:6], point[:6], e139.PRIME
        )
        lift = None
        if determinant:
            lift, _ = e139.hensel_lift_slice(
                sliced, setup["non_diagonal"][:6], point[:6], e139.PRIME, 4
            )
        determinants.append(int(determinant))
        lifts.append({
            "point13_mod5": point,
            "jacobian_determinant_mod5": int(determinant),
            "jacobian_mod5": jacobian,
            "lift_solve6_mod625": None if lift is None else [int(value) for value in lift],
            "lift_success_mod625": lift is not None,
        })
    prior_section = prior["section"]
    prior_points = prior["thinvariants"]["column"]["solutions"]
    prior_lifts = prior["hensel_components"]["lift_records"]
    prior_dets = [int(record["jacobian_determinant_mod5"]) for record in prior_lifts]
    prior_holdout = prior["holdout"]
    reproduced_holdout = {
        "orders": [4, 6, 8, 10, 12],
        "k8_residues": [int(value) for value in prior_holdout["k8_residues"]],
        "k10_residues": [int(value) for value in prior_holdout["k10_residues"]],
        "k12_residues": [int(value) for value in prior_holdout["k12_residues"]],
    }
    checks = {
        "held_values_match": prior_section["held_values_mod5"] == anchor,
        "held_names_match": prior_section["held_variable_names"] == list(e139.NONDIAGONAL_NAMES[6:]),
        "solve_names_match": prior_section["solve_variable_names"] == list(e139.SLICE6_NAMES),
        "solution_count_match": len(solutions) == int(prior["thinvariants"]["column"]["solutions_total_mod5"]) == 3,
        "solutions_match": sorted(solutions) == sorted(prior_points),
        "jacobian_determinants_match": determinants == prior_dets == [2, 4, 1],
        "all_anchor_lifts_success": all(record["lift_success_mod625"] for record in lifts),
        "holdout_controls_present": reproduced_holdout["k8_residues"] == [250, 456, 430]
        and reproduced_holdout["k10_residues"] == [55, 350, 480]
        and reproduced_holdout["k12_residues"] == [266, 132, 286],
    }
    if not all(checks.values()):
        raise AssertionError(f"anchor cross-check failed: {checks}")
    return {
        "status": "PASS",
        "independent_recomputed": True,
        "held_values_mod5": anchor,
        "solve_grid_size": TOTAL_SOLVE_POINTS,
        "solutions_mod5": solutions,
        "solution_count": len(solutions),
        "jacobian_determinants_mod5": determinants,
        "lifts": lifts,
        "holdout_residues_mod625_from_components": reproduced_holdout,
        "construction_residue_control": "e139 records all 42 construction residues zero mod 625 for all three lifts",
        "checks": checks,
    }


def census_section(e139, setup: dict[str, Any], held: list[int], deadline: float) -> dict[str, Any] | None:
    """Run one section; return None only when the explicit cap expires."""
    if time.process_time() >= deadline:
        return None
    solutions = section_solutions(e139, setup["primitives"], held)
    sliced = slice_polynomials(e139, setup["primitives"], setup["non_diagonal"], held)
    records: list[dict[str, Any]] = []
    for point in solutions:
        if time.process_time() >= deadline:
            return None
        determinant, jacobian = e139.jacobian_determinant_mod_prime(
            sliced, setup["non_diagonal"][:6], point[:6], e139.PRIME
        )
        lift = None
        if determinant:
            lift, _ = e139.hensel_lift_slice(
                sliced, setup["non_diagonal"][:6], point[:6], e139.PRIME, 4
            )
        records.append({
            "point13_mod5": point,
            "jacobian_determinant_mod5": int(determinant),
            "jacobian_mod5": jacobian,
            "lift_success_mod625": lift is not None,
            "lift_solve6_mod625": None if lift is None else [int(value) for value in lift],
            "holdout_status": "NOT_LAUNCHED",
        })
    return {
        "held_values_mod5": held,
        "solve_grid_size": TOTAL_SOLVE_POINTS,
        "solution_count_mod5": len(solutions),
        "nonsingular_solution_count_mod5": sum(
            int(record["jacobian_determinant_mod5"] != 0) for record in records
        ),
        "solutions": records,
        "holdout_status": "NOT_LAUNCHED",
    }


def make_manifest(payload: dict[str, Any], prior_hash: str, e139_hash: str) -> dict[str, Any]:
    stable = copy.deepcopy(payload)
    provenance = stable.get("provenance", {})
    provenance.pop("generated_utc", None)
    provenance.pop("process_seconds", None)
    stable["provenance"] = provenance
    return {
        "canonicalization": "JSON sort_keys compact UTF-8; generated_utc/process_seconds excluded",
        "canonical_body_sha256": hashlib.sha256(canonical_json(stable)).hexdigest(),
        "prior_components_sha256": prior_hash,
        "e139_source_sha256": e139_hash,
    }


def create_payload(e139, prior, prior_hash: str, e139_hash: str, cap_seconds: float, started: float):
    deadline = started + cap_seconds
    setup_started = time.process_time()
    setup = setup_exact(e139)
    setup_seconds = time.process_time() - setup_started
    anchor = anchor_control(e139, setup, prior)
    anchor_seconds = time.process_time() - setup_started
    sections: list[dict[str, Any]] = []
    census_launched = False
    preflight_reason = ""
    estimate_seconds: float | None = None
    section_probe_seconds: float | None = None
    anchor_tuple = tuple(e139.ANCHOR_HELD)
    if time.process_time() < deadline:
        probe_started = time.process_time()
        probe = census_section(e139, setup, list(anchor_tuple), deadline)
        section_probe_seconds = time.process_time() - probe_started
        if probe is None:
            preflight_reason = "the explicit process-time cap expired before the first section completed"
        else:
            sections.append(probe)
            census_launched = True
            estimate_seconds = section_probe_seconds * TOTAL_SECTIONS
            if time.process_time() + estimate_seconds > deadline:
                preflight_reason = (
                    "the first exact section completed, but its measured cost proves the remaining "
                    "exhaustive 4^7 launch cannot fit the explicit cap"
                )
    else:
        preflight_reason = "exact setup and anchor control consumed the explicit process-time cap"

    loop_started = False
    if census_launched and not preflight_reason:
        loop_started = True
        for held_tuple in itertools.product(UNIT_VALUES, repeat=7):
            if held_tuple == anchor_tuple:
                continue
            record = census_section(e139, setup, list(held_tuple), deadline)
            if record is None:
                preflight_reason = "the explicit process-time cap expired during exhaustive section enumeration"
                break
            sections.append(record)
    all_sections = len(sections) == TOTAL_SECTIONS
    holdout_complete = False  # The bounded producer never claims unlaunched holdout routes.
    complete = all_sections and holdout_complete
    if complete:
        status = "COMPLETE"
    else:
        status = "PREFLIGHT"
        if not preflight_reason:
            preflight_reason = (
                "section enumeration completed, but independent order-8/10/12 holdout routes "
                "were not launched; no complete claim is permitted"
            )

    payload: dict[str, Any] = {
        "schema_version": 1,
        "status": status,
        "claim_tag": "[UNRESOLVED]",
        "headline": (
            "[COMPUTATION] Exact unit-grid section census for Kac--Ward branch 11111; "
            f"{len(sections)} of {TOTAL_SECTIONS} held sections were launched under the cap. "
            "[UNRESOLVED] This artifact does not close the full Q_5 or characteristic-zero branch."
        ),
        "arithmetic": {
            "field_for_section_census": "F_5^x",
            "unit_values": list(UNIT_VALUES),
            "prime": int(e139.PRIME),
            "modulus": int(e139.MODULUS),
            "exact_integer_residues_only": True,
            "solve_grid_size_per_section": TOTAL_SOLVE_POINTS,
        },
        "construction": {
            "branch": "11111",
            "selected_numerator_indices": list(e139.SELECTED_NUMERATOR_INDICES),
            "nondiagonal_variable_order": list(e139.NONDIAGONAL_NAMES),
            "held_variable_order": list(e139.NONDIAGONAL_NAMES[6:]),
            "solve_variable_order": list(e139.SLICE6_NAMES),
            "primitive_term_counts": [len(poly.terms()) for poly in setup["primitives"]],
            "e139_reused": [
                "full_weight_finite_system", "branch_11111_equations", "thin_torus_substitutions",
                "primitive_numerators", "eval_poly_mod", "jacobian_determinant_mod_prime", "hensel_lift_slice",
            ],
        },
        "anchor_cross_check": anchor,
        "section_coverage": {
            "method": "exhaustive_4^7_unit_sections",
            "orbit_quotient_used": False,
            "total_sections": TOTAL_SECTIONS,
            "sections_launched": len(sections),
            "sections_not_launched": TOTAL_SECTIONS - len(sections),
            "all_sections_completed": all_sections,
            "section_order": "lexicographic product over held values (1,2,3,4)",
            "representatives": sections,
            "orbit_kernel_from_e139": prior["orbit_union_lemma"]["kernel_vectors"],
            "orbit_kernel_note": (
                "The e139 diagonal-torus kernel is preserved as provenance, but no quotient is claimed; "
                "coverage is exhaustive only when all 4^7 sections are present."
            ),
        },
        "holdout": {
            "shape": [3, 3, 3],
            "orders_required": [8, 10, 12],
            "status": "NOT_LAUNCHED",
            "independent_routes_required": ["symbolic_walk", "matrix_power"],
            "complete_for_all_lifts": holdout_complete,
            "scope_note": "No holdout conclusion is asserted for sections/lifts not evaluated.",
        },
        "theorem_scope": {
            "q5_conclusion": (
                "Only a nonsingular F_5 solution with a recorded successful mod-625 Hensel lift can "
                "support a Q_5-local conclusion; no Q_5 claim is made for missing, singular, or "
                "chart-degenerate cases."
            ),
            "full_branch": "[UNRESOLVED] The full Q_5 and characteristic-zero branch remains unresolved.",
            "timeout_is_not_theorem": True,
        },
        "preflight": {
            "reason": preflight_reason,
            "what_was_launched": {
                "exact_e139_setup": True,
                "anchor_control": True,
                "section_probe": census_launched,
                "exhaustive_section_loop": loop_started,
                "holdout_routes": False,
            },
            "what_was_not_launched": {
                "unvisited_sections": TOTAL_SECTIONS - len(sections),
                "holdout_evaluation_for_all_lifts": True,
                "full_branch_proof": True,
            },
            "measured_setup_seconds": round(setup_seconds, 6),
            "measured_anchor_seconds_from_setup": round(anchor_seconds, 6),
            "first_section_seconds": None if section_probe_seconds is None else round(section_probe_seconds, 6),
            "estimated_full_section_seconds": None if estimate_seconds is None else round(estimate_seconds, 6),
            "cap_seconds": cap_seconds,
        },
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "script": "experiments/e156_kw_full_thin_census.py",
            "prior_artifact": "results/kac_ward/components.json",
            "e139_source": "experiments/e139_kw_components.py",
            "process_seconds": round(time.process_time() - started, 6),
            "no_float_decisions": True,
        },
    }
    payload["manifest"] = make_manifest(payload, prior_hash, e139_hash)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap-seconds", type=float, default=DEFAULT_CAP_SECONDS)
    args = parser.parse_args()
    if args.cap_seconds <= 0:
        raise SystemExit("--cap-seconds must be positive")
    started = time.process_time()
    e139_hash = sha256_bytes(E139_PATH)
    prior, prior_hash = load_prior()
    e139 = load_e139()
    payload = create_payload(e139, prior, prior_hash, e139_hash, args.cap_seconds, started)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    coverage = payload["section_coverage"]
    if payload["status"] == "COMPLETE":
        print(f"PASS COMPLETE: {coverage['sections_launched']}/{coverage['total_sections']} sections")
    else:
        print(
            f"PREFLIGHT: {coverage['sections_launched']}/{coverage['total_sections']} sections launched; "
            f"{payload['preflight']['reason']}"
        )
    print(f"PREFLIGHT artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
