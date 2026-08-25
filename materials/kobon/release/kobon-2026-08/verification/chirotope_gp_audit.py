#!/usr/bin/env python3
"""Exact audit and relaxation-gap probe for projective chirotope cuts."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
if not (ENGINE_DIR / "engine.py").is_file():
    ENGINE_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(ENGINE_DIR))
import engine  # noqa: E402

MODES = (
    ("generic", 0, 0),
    ("parallel_pair", 1, 0),
    ("two_parallel_pairs", 2, 0),
    ("triple_point", 0, 1),
    ("two_triple_points", 0, 2),
    ("parallel_and_triple", 1, 1),
)


def sign(value):
    return (value > 0) - (value < 0)


def dual_determinant(lines, triple):
    first, second, third = triple
    m0, b0 = lines[first]
    m1, b1 = lines[second]
    m2, b2 = lines[third]
    return (m1 - m0) * (b2 - b0) - (m2 - m0) * (b1 - b0)


def encoded_sign(lines, triple):
    first, middle, last = triple
    m0, b0 = lines[first]
    m1, b1 = lines[middle]
    m2, b2 = lines[last]
    lower_parallel = m0 == m1
    upper_parallel = m1 == m2
    if lower_parallel and upper_parallel:
        return 0
    if lower_parallel:
        return 1 if b0 > b1 else -1
    if upper_parallel:
        return 1 if b2 > b1 else -1

    x01 = (b1 - b0) / (m0 - m1)
    x02 = (b2 - b0) / (m0 - m2)
    if x01 == x02:
        return 0
    return 1 if x01 > x02 else -1


def ordered_sign(signs, ordered):
    inversions = sum(
        ordered[left] > ordered[right]
        for left in range(3) for right in range(left + 1, 3)
    )
    value = signs[tuple(sorted(ordered))]
    return value if inversions % 2 == 0 else -value


def exact_audit(samples_per_mode):
    counts = {
        "arrangements": 0,
        "triples": 0,
        "zero_triples": 0,
        "gp_relations": 0,
        "sign_mismatches": 0,
        "gp_violations": 0,
    }
    examples = []
    for n in range(5, 11):
        for mode_index, (mode, parallel_pairs, triple_points) in enumerate(MODES):
            for sample in range(samples_per_mode):
                seed = 20_000_000 * n + 100_000 * mode_index + sample
                slopes, intercepts = engine.gen_config(
                    n,
                    n_par_pairs=parallel_pairs,
                    n_conc=triple_points,
                    coord=10**5,
                    seed=seed,
                )
                lines = tuple(zip(slopes, intercepts))
                signs = {}
                counts["arrangements"] += 1
                for triple in combinations(range(n), 3):
                    exact = sign(dual_determinant(lines, triple))
                    encoded = encoded_sign(lines, triple)
                    signs[triple] = exact
                    counts["triples"] += 1
                    counts["zero_triples"] += exact == 0
                    if exact != encoded:
                        counts["sign_mismatches"] += 1
                        examples.append((
                            mode, n, seed, "sign", triple, exact, encoded))

                for five in combinations(range(n), 5):
                    for anchor in five:
                        first, second, third, fourth = (
                            line for line in five if line != anchor)
                        products = (
                            ordered_sign(signs, (anchor, first, second))
                            * ordered_sign(signs, (anchor, third, fourth)),
                            -ordered_sign(signs, (anchor, first, third))
                            * ordered_sign(signs, (anchor, second, fourth)),
                            ordered_sign(signs, (anchor, first, fourth))
                            * ordered_sign(signs, (anchor, second, third)),
                        )
                        counts["gp_relations"] += 1
                        has_positive = any(value > 0 for value in products)
                        has_negative = any(value < 0 for value in products)
                        if has_positive != has_negative:
                            counts["gp_violations"] += 1
                            examples.append((
                                mode, n, seed, "gp", five, anchor, products))
    return counts, examples


def relaxation_gap_probe():
    for n in range(5, 8):
        cnf, pool = engine.build_model(n, 1, rules=())
        _, relations = engine.add_projective_chirotope_constraints(
            cnf, pool, n, enforce=False)
        with engine.Solver(
                name="cadical195", bootstrap_with=cnf.clauses) as solver:
            witness = None
            for relation_index, products in enumerate(relations):
                for product_index, (positive, _) in enumerate(products):
                    assumptions = [positive] + [
                        -negative for _, negative in products]
                    if solver.solve(assumptions=assumptions):
                        witness = (
                            relation_index, product_index, "positive", assumptions)
                        break
                if witness is not None:
                    break
                for product_index, (_, negative) in enumerate(products):
                    assumptions = [negative] + [
                        -positive for positive, _ in products]
                    if solver.solve(assumptions=assumptions):
                        witness = (
                            relation_index, product_index, "negative", assumptions)
                        break
                if witness is not None:
                    break
            base_stats = solver.accum_stats()
        if witness is None:
            continue

        engine.add_projective_chirotope_constraints(cnf, pool, n, enforce=True)
        with engine.Solver(
                name="cadical195", bootstrap_with=cnf.clauses) as solver:
            rejected = not solver.solve(assumptions=witness[3])
            cut_stats = solver.accum_stats()
        return {
            "base_relaxation_admits_gp_violation": True,
            "gp_cut_rejects_witness": rejected,
            "n": n,
            "relation_index": witness[0],
            "product_index": witness[1],
            "one_sided_sign": witness[2],
            "base_stats": base_stats,
            "cut_stats": cut_stats,
        }
    return {
        "base_relaxation_admits_gp_violation": False,
        "gp_cut_rejects_witness": False,
        "n_range": [5, 7],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-mode", type=int, default=100)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    counts, examples = exact_audit(args.samples_per_mode)
    probe = relaxation_gap_probe()
    passed = (
        not examples
        and (
            not probe["base_relaxation_admits_gp_violation"]
            or probe["gp_cut_rejects_witness"]
        )
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "samples_per_mode": args.samples_per_mode,
        "n_range": [5, 10],
        "modes": [mode for mode, _, _ in MODES],
        "counts": counts,
        "relaxation_gap_probe": probe,
        "semantic_strength": (
            "GAP_FOUND_AND_CLOSED"
            if probe["base_relaxation_admits_gp_violation"]
            else "NO_GAP_FOUND_THROUGH_N7_PROPAGATION_ONLY"),
        "violation_examples": examples[:20],
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
