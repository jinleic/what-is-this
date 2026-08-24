#!/usr/bin/env python3
"""Probe a straight-line incidence axiom absent from the order relaxation.

Dual Pappus: two distinct pencils A,B,C and a,b,c, together with the three
cross-connector lines x,y,z, force x,y,z to be concurrent projectively.  In an
affine chart the conclusion is finite concurrency or one parallel class.

This script exhibits a SAT model of the current abstract order axioms with the
nondegenerate Pappus antecedent and a false conclusion, then checks that the
three projective-conclusion clauses reject exactly that model class.
"""
from __future__ import annotations
import argparse

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
if not (ENGINE_DIR / "engine.py").is_file():
    ENGINE_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(ENGINE_DIR))
import engine  # noqa: E402

ROLES = tuple(range(9))


def variables(pool):
    A, B, C, a, b, c, x, y, z = ROLES
    concurrent = lambda *triple: pool.id(
        ("C",) + tuple(sorted(triple)))
    crossing = lambda left, right: pool.id(
        ("P",) + tuple(sorted((left, right))))
    antecedent = [
        concurrent(A, B, C),
        concurrent(a, b, c),
        concurrent(A, b, x),
        concurrent(a, B, x),
        concurrent(A, c, y),
        concurrent(a, C, y),
        concurrent(B, c, z),
        concurrent(b, C, z),
    ]
    distinct_centres = -concurrent(A, B, a)
    conclusion = concurrent(x, y, z)
    conclusion_pairs = (
        crossing(x, y), crossing(x, z), crossing(y, z))
    return antecedent, distinct_centres, conclusion, conclusion_pairs


def add_dual_pappus_cut(cnf, pool):
    antecedent, distinct_centres, conclusion, conclusion_pairs = variables(pool)
    guard_failure = [-literal for literal in antecedent]
    # distinct_centres is the literal -C(A,B,a); its failure is +C(A,B,a).
    guard_failure.append(-distinct_centres)
    for pair_crossing in conclusion_pairs:
        cnf.append(guard_failure + [conclusion, -pair_crossing])


def solve(cnf, assumptions):
    with engine.Solver(name="cadical195", bootstrap_with=cnf.clauses) as solver:
        return solver.solve(assumptions=assumptions), solver.accum_stats()


def main() -> int:
    cnf, pool = engine.build_model(9, 1, rules=())
    antecedent, distinct_centres, conclusion, conclusion_pairs = variables(pool)
    violating = antecedent + [distinct_centres, -conclusion, conclusion_pairs[0]]
    base_sat, base_stats = solve(cnf, violating)

    add_dual_pappus_cut(cnf, pool)
    cut_sat, cut_stats = solve(cnf, violating)
    result = {
        "status": "PASS" if base_sat and not cut_sat else "FAIL",
        "roles": ["A", "B", "C", "a", "b", "c", "x", "y", "z"],
        "labels": list(ROLES),
        "nondegenerate_distinct_pencils_guard": "not C(A,B,a)",
        "base_relaxation_admits_false_pappus_conclusion": base_sat,
        "pappus_cut_rejects_false_conclusion": not cut_sat,
        "base_stats": base_stats,
        "cut_stats": cut_stats,
        "base_dimensions": {
            "variables": cnf.nv,
            "clauses_after_cut": len(cnf.clauses),
        },
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    if result["status"] != "PASS":
        raise AssertionError("Pappus relaxation probe failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
