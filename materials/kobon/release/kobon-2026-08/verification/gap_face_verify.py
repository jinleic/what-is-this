#!/usr/bin/env python3
"""Exact audit of the direct-gap triangular-face criterion.

Compares two independent predicates on every nondegenerate line triple:
(1) no outside line strictly straddles the triangle's vertices; and
(2) no outside-line intersection lies strictly inside any of its three sides.
All arithmetic is Fraction.  Random suites include generic arrangements,
parallel pairs, planted triple points, mixed degeneracies and planted 4-fold
points.  Any mismatch exits nonzero and prints the first counterexample.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402


def straddle_face(n, ms, bs, X, triple):
    vertices = engine.tri_verts(X, triple)
    for outside in range(n):
        if outside in triple:
            continue
        signs = [y - ms[outside] * x - bs[outside] for x, y in vertices]
        if min(signs) < 0 < max(signs):
            return False
    return True


def direct_gap_face(n, X, triple):
    for side in triple:
        endpoints = [engine.cross(X, side, x) for x in triple if x != side]
        x0, x1 = sorted((endpoints[0][0], endpoints[1][0]))
        for outside in range(n):
            if outside in triple:
                continue
            point = engine.cross(X, side, outside)
            if point is not None and x0 < point[0] < x1:
                return False
    return True


def planted_fourfold(n, seed):
    rng = random.Random(seed)
    ms = sorted(F(x) for x in rng.sample(range(-100, 101), n))
    bs = [F(rng.randrange(-100, 101)) for _ in range(n)]
    x0, y0 = F(rng.randrange(-7, 8)), F(rng.randrange(-7, 8))
    for line in rng.sample(range(n), 4):
        bs[line] = y0 - ms[line] * x0
    return ms, bs


def configurations(n, trials):
    for seed in range(trials):
        yield "generic", seed, engine.gen_config(n, seed=10_000 * n + seed)
        yield "parallel", seed, engine.gen_config(
            n, n_par_pairs=1, seed=20_000 * n + seed)
        yield "triple", seed, engine.gen_config(
            n, n_conc=1, seed=30_000 * n + seed)
        yield "parallel+triple", seed, engine.gen_config(
            n, n_par_pairs=1, n_conc=1, seed=40_000 * n + seed)
        yield "fourfold", seed, planted_fourfold(n, 50_000 * n + seed)


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    checked_configurations = 0
    checked_triples = 0
    mode_counts = {}
    square_violations = 0
    minimum_square_slack = None
    for n in range(5, 9):
        for mode, seed, (ms, bs) in configurations(n, trials):
            if len(set(zip(ms, bs))) != n:
                continue
            X = engine.crossings(n, ms, bs)
            checked_configurations += 1
            mode_counts[mode] = mode_counts.get(mode, 0) + 1
            face_count = 0
            for triple in combinations(range(n), 3):
                if not engine.tri_ok(X, triple):
                    continue
                strict = straddle_face(n, ms, bs, X, triple)
                direct = direct_gap_face(n, X, triple)
                checked_triples += 1
                if direct:
                    face_count += 1
                if strict != direct:
                    print(json.dumps({
                        "status": "MISMATCH", "n": n, "mode": mode,
                        "seed": seed, "triple": triple,
                        "straddle_face": strict, "direct_gap_face": direct,
                        "lines": [[str(m), str(b)] for m, b in zip(ms, bs)],
                    }, indent=2))
                    return 1

            parallel_pairs = 0
            point_lines = {}
            for i, j in combinations(range(n), 2):
                point = engine.cross(X, i, j)
                if point is None:
                    parallel_pairs += 1
                else:
                    point_lines.setdefault(point, set()).update((i, j))
            multiplicities = sorted(
                len(lines) for lines in point_lines.values() if len(lines) >= 3
            )
            square_rhs = (
                n * (n - 2)
                - 2 * parallel_pairs
                - sum((k - 2) ** 2 for k in multiplicities)
            )
            square_slack = square_rhs - 3 * face_count
            record = {
                "slack": square_slack, "n": n, "mode": mode, "seed": seed,
                "faces": face_count, "parallel_pairs": parallel_pairs,
                "multiplicities": multiplicities,
            }
            if (minimum_square_slack is None
                    or square_slack < minimum_square_slack["slack"]):
                minimum_square_slack = record
            if square_slack < 0:
                square_violations += 1
                print(json.dumps({
                    "status": "SQUARE_PENALTY_VIOLATION", **record,
                    "lines": [[str(m), str(b)] for m, b in zip(ms, bs)],
                }, indent=2))
                return 1
    print(json.dumps({
        "status": "PASS",
        "trials_per_n_per_mode": trials,
        "configurations": checked_configurations,
        "nondegenerate_triples": checked_triples,
        "modes": mode_counts,
        "mismatches": 0,
        "square_penalty_violations": square_violations,
        "minimum_square_penalty_slack": minimum_square_slack,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
