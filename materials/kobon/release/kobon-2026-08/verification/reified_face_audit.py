#!/usr/bin/env python3
"""Audit exact face reification against embedded rational certificates.

The audit fixes every crossing/parallel, concurrency, strict line-order and face
literal from exact rational geometry.  SAT then checks that the reified CNF,
including endpoint closure, multipoint incidence, and K(4) cuts, accepts the
known n=10 optimum and n=12 38-face record.  The independently evaluated
direct-gap face set must also equal the embedded certificate set exactly.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
GAP_PATH = ROOT / "scratch" / "kobon" / "gap_faces.py"
if not (ENGINE_DIR / "engine.py").is_file():
    bundle = Path(__file__).resolve().parents[1]
    ENGINE_DIR = bundle / "scripts"
    GAP_PATH = bundle / "scripts" / "gap_faces.py"
sys.path.insert(0, str(ENGINE_DIR))
import engine  # noqa: E402

SPEC = importlib.util.spec_from_file_location("gap_faces", GAP_PATH)
gap_faces = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(gap_faces)


def slope_sorted(lines, selected):
    order = sorted(range(len(lines)), key=lambda index: lines[index][0])
    inverse = {old: new for new, old in enumerate(order)}
    sorted_lines = [lines[old] for old in order]
    sorted_selected = {
        tuple(sorted(inverse[index] for index in triple))
        for triple in selected
    }
    return sorted_lines, sorted_selected


def direct_faces(lines):
    n = len(lines)
    ms = [engine.Fr(slope) for slope, _ in lines]
    bs = [engine.Fr(intercept) for _, intercept in lines]
    points = engine.crossings(n, ms, bs)
    faces = set()
    for triple in combinations(range(n), 3):
        if not engine.tri_ok(points, triple):
            continue
        is_face = True
        for side in triple:
            endpoints = [
                engine.cross(points, side, line)
                for line in triple if line != side
            ]
            lower, upper = sorted((endpoints[0][0], endpoints[1][0]))
            for outside in range(n):
                if outside in triple:
                    continue
                point = engine.cross(points, side, outside)
                if point is not None and lower < point[0] < upper:
                    is_face = False
                    break
            if not is_face:
                break
        if is_face:
            faces.add(triple)
    return points, faces


def fixed_geometry_assumptions(pool, points, faces, n):
    assumptions = []
    for first, second in combinations(range(n), 2):
        literal = pool.id(("P", first, second))
        assumptions.append(
            literal if engine.cross(points, first, second) is not None
            else -literal)

    for triple in combinations(range(n), 3):
        crossings = [
            engine.cross(points, first, second)
            for first, second in combinations(triple, 2)
        ]
        is_concurrent = (
            crossings[0] is not None
            and crossings[0] == crossings[1] == crossings[2]
        )
        literal = pool.id(("C",) + triple)
        assumptions.append(literal if is_concurrent else -literal)

    for line in range(n):
        others = [other for other in range(n) if other != line]
        for first in others:
            for second in others:
                if first == second:
                    continue
                first_point = engine.cross(points, line, first)
                second_point = engine.cross(points, line, second)
                precedes = (
                    first_point is not None
                    and second_point is not None
                    and first_point != second_point
                    and first_point[0] < second_point[0]
                )
                literal = pool.id(("X", line, first, second))
                assumptions.append(literal if precedes else -literal)

    for triple in combinations(range(n), 3):
        literal = pool.id(("S",) + triple)
        assumptions.append(literal if triple in faces else -literal)
    return assumptions


def audit_case(
        name, lines, selected, subarrangement_bounds, simple_upper):
    sorted_lines, expected_faces = slope_sorted(lines, selected)
    points, faces = direct_faces(sorted_lines)
    n = len(sorted_lines)
    target = len(expected_faces)
    point_lines = {}
    for first, second in combinations(range(n), 2):
        point = engine.cross(points, first, second)
        if point is not None:
            point_lines.setdefault(point, set()).update((first, second))
    multipoint_incident_faces = sum(
        any(
            len(point_lines[engine.cross(points, first, second)]) >= 3
            for first, second in combinations(face, 2)
        )
        for face in faces
    )
    required_incident_faces = max(0, target - simple_upper)
    cnf, pool = gap_faces.build(
        n,
        target,
        subarrangement_bounds=subarrangement_bounds,
        endpoint_closure=True,
        reify_faces=True,
        simple_upper=simple_upper,
        k4_bound=True,
    )
    assumptions = fixed_geometry_assumptions(pool, points, faces, n)
    with engine.Solver(
            name="cadical195", bootstrap_with=cnf.clauses) as solver:
        satisfiable = solver.solve(assumptions=assumptions)
        stats = solver.accum_stats()
    return {
        "case": name,
        "n": n,
        "target": target,
        "direct_face_count": len(faces),
        "face_set_matches_certificate": faces == expected_faces,
        "fixed_geometry_cnf_sat": satisfiable,
        "multipoint_incident_faces": multipoint_incident_faces,
        "required_multipoint_incident_faces": required_incident_faces,
        "simple_perturbation_bound_holds": (
            multipoint_incident_faces >= required_incident_faces),
        "variables": cnf.nv,
        "clauses": len(cnf.clauses),
        "solver_stats": stats,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    cases = [
        audit_case(
            "K10=25",
            engine.N10_LOWER_BOUND_LINES,
            engine.N10_LOWER_BOUND_TRIANGLES,
            {9: 21, 8: 15},
            25,
        ),
        audit_case(
            "K12>=38",
            engine.N12_LOWER_BOUND_LINES,
            engine.N12_LOWER_BOUND_TRIANGLES,
            {11: 32, 10: 25},
            37,
        ),
    ]
    passed = all(
        case["face_set_matches_certificate"]
        and case["fixed_geometry_cnf_sat"]
        and case["simple_perturbation_bound_holds"]
        for case in cases
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "predicate": "exact direct-gap face reification",
        "cases": cases,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
