#!/usr/bin/env python3
"""Search for a combinatorial counterexample to a weighted face penalty.

For an essential n-line arrangement, with Q parallel pairs, triangular-face
count F, and finite multipoints p of multiplicity k_p >= 3, test

    3 F + 2 Q + c * sum_p (k_p - 2)^2 <= n (n - 2),

where c is a nonnegative rational coefficient (default 1). The direct-gap
clauses make every selected triple a triangular arrangement face. Selecting
any subset is sufficient: a violating arrangement can select all of its
faces, while a violating selected subset already proves that the arrangement
violates the inequality.

The order-table model is a necessary relaxation for straight-line
arrangements. Therefore SAT is a combinatorial candidate requiring
straightening; UNSAT becomes a theorem only after its emitted proof is checked.

Usage: square_penalty_sat.py N [--coefficient NUM/DEN]
                             [--out FILE] [--solve]
                             [--solver cadical195]
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path
import resource
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "math" / "kobon"
if not (ENGINE_DIR / "engine.py").is_file():
    ENGINE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ENGINE_DIR))
import engine  # noqa: E402

from gap_faces import add_direct_face_clauses  # noqa: E402


def rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def add_square_violation(
        cnf, pool, n: int, coefficient: Fraction) -> list[int]:
    """Require a strict violation and return its scaled score literals.

    Multiplying by the coefficient denominator makes every weight integral.
    A k-line point has one canonical first triple and u = k - 3 later lines.
    Since (k - 2)^2 = (u + 1)^2 = 1 + 3u + 2*C(u, 2), duplicating literals
    gives an exact unit-weight cardinality encoding. The canonical point
    variables use the same concurrency closure as ``engine.add_face_bound``.
    """
    sp = lambda *xs: tuple(sorted(xs))
    crosses = lambda i, j: pool.id(("P",) + sp(i, j))
    concurrent = lambda *xs: pool.id(("C",) + sp(*xs))
    selected = lambda triple: pool.id(("S",) + tuple(sorted(triple)))

    score: list[int] = []
    numerator, denominator = coefficient.numerator, coefficient.denominator
    for triple in combinations(range(n), 3):
        score.extend([selected(triple)] * (3 * denominator))
    for i, j in combinations(range(n), 2):
        score.extend([-crosses(i, j)] * (2 * denominator))

    for a, b, c in combinations(range(n), 3):
        is_concurrent = concurrent(a, b, c)
        earlier = [
            concurrent(a, b, d)
            for d in range(c)
            if d not in (a, b)
        ]
        first = pool.id(("SQP", a, b, c))
        cnf.append([-first, is_concurrent])
        for old in earlier:
            cnf.append([-first, -old])
        cnf.append([-is_concurrent] + earlier + [first])
        score.extend([first] * numerator)

        later = []
        for d in range(c + 1, n):
            through = concurrent(a, b, d)
            extra = pool.id(("SQU", a, b, c, d))
            cnf.append([-extra, first])
            cnf.append([-extra, through])
            cnf.append([-first, -through, extra])
            score.extend([extra] * (3 * numerator))
            later.append(extra)
        for d, e in combinations(later, 2):
            pair = pool.id(("SQV", a, b, c, d, e))
            cnf.append([-pair, d])
            cnf.append([-pair, e])
            cnf.append([-d, -e, pair])
            score.extend([pair] * (2 * numerator))

    threshold = denominator * n * (n - 2) + 1
    if threshold > len(score):
        cnf.append([])
    else:
        cnf.extend(CardEnc.atleast(
            lits=score,
            bound=threshold,
            vpool=pool,
            encoding=EncType.seqcounter,
        ))
    return score


def build(n: int, coefficient: Fraction = Fraction(1)):
    if n < 3:
        raise ValueError("n must be at least 3")
    if coefficient < 0:
        raise ValueError("coefficient must be nonnegative")
    cnf, pool = engine.build_model(n, 0, rules=())
    cnf.append([
        pool.id(("P", i, j)) for i, j in combinations(range(n), 2)
    ])
    add_direct_face_clauses(cnf, pool, n)
    score = add_square_violation(cnf, pool, n, coefficient)
    return cnf, pool, score


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    n = int(sys.argv[1])
    output = None
    solver = "cadical195"
    coefficient = Fraction(1)
    solve = "--solve" in sys.argv
    for i, arg in enumerate(sys.argv):
        if arg == "--out":
            output = Path(sys.argv[i + 1])
        elif arg == "--solver":
            solver = sys.argv[i + 1]
        elif arg == "--coefficient":
            coefficient = Fraction(sys.argv[i + 1])

    started = time.time()
    cnf, pool, score = build(n, coefficient)
    print(json.dumps({
        "stage": "built",
        "formulation": "direct-gap-weighted-square-penalty-violation",
        "n": n,
        "coefficient": str(coefficient),
        "threshold": coefficient.denominator * n * (n - 2) + 1,
        "score_occurrences": len(score),
        "vars": cnf.nv,
        "clauses": len(cnf.clauses),
        "seconds": round(time.time() - started, 3),
        "rss_gb": round(rss_gb(), 3),
    }), flush=True)

    if output is not None:
        cnf.to_file(str(output))
        print(json.dumps({
            "stage": "written",
            "path": str(output),
            "bytes": output.stat().st_size,
        }), flush=True)
    if not solve:
        return 0

    solve_started = time.time()
    with Solver(name=solver, bootstrap_with=cnf.clauses) as sat:
        satisfiable = sat.solve()
        stats = sat.accum_stats()
        model = sat.get_model() if satisfiable else None
    result = {
        "stage": "VERDICT",
        "n": n,
        "coefficient": str(coefficient),
        "result": "SAT-candidate" if satisfiable else "UNSAT-discovery",
        "seconds": round(time.time() - solve_started, 3),
        "stats": stats,
    }
    if model is not None:
        truth = {literal for literal in model if literal > 0}
        result["encoded_score"] = sum(abs(lit) in truth if lit > 0 else abs(lit) not in truth
                                      for lit in score)
        result.update(engine.extract(model, pool, n))
    print(json.dumps(result), flush=True)
    return 10 if satisfiable else 20


if __name__ == "__main__":
    sys.exit(main())
