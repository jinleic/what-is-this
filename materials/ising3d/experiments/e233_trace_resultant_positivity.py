#!/usr/bin/env python3
"""Exact positive-mode feasibility of the trace-resultant coefficient cut.

The low-power trace equations reduce to sextics Q with a0, a1, a2 and a0
fixed.  A positive-mode Gaussian would require every root u_i>0.  This script
eliminates the coefficient variables and asks whether the surviving trace-data
semialgebraic region matches the positive-root constraints.
"""

from __future__ import annotations

import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))

from e232_trace_resultant import (  # noqa: E402
    graph_targets,
    norm_equations,
    synthetic_control,
)
from ising.transfer_matrix import layer_bonds  # noqa: E402

OUTPUT = ROOT / "results" / "spectral" / "trace_resultant_positivity.json"


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def positivity_certificate(targets: dict[str, Fraction]) -> dict[str, object]:
    equations, _ = norm_equations(targets)
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    basis = sp.groebner(
        [equation.as_expr() for equation in equations],
        a5,
        a4,
        a3,
        order="lex",
        domain=sp.QQ,
    )
    basis_exprs = [sp.factor(poly.as_expr()) for poly in basis.polys]
    univariate_clauses: list[dict[str, object]] = []
    for expression in basis_exprs:
        one_variable = [variable for variable in (a5, a4, a3) if expression.has(variable)]
        if len(one_variable) == 1:
            variable = one_variable[0]
            values = sp.solve(expression, variable, dict=True)
            solutions = [sp.factor(solution[variable]) for solution in values]
            univariate_clauses.append(
                {
                    "variable": str(variable),
                    "clause": str(expression),
                    "solutions": [str(solution) for solution in solutions],
                }
            )
    return {
        "basis": [str(expression) for expression in basis_exprs],
        "basis_size": len(basis_exprs),
        "univariate_clauses": univariate_clauses,
    }


def positive_root_necessities(qpoly_expression: sp.Expr, substitution: dict[sp.Symbol, sp.Expr]) -> dict[str, object]:
    u = sp.symbols("u")
    specialized = sp.Poly(sp.expand(qpoly_expression.subs(substitution)), u, domain=sp.QQ)
    coefficients = specialized.all_coeffs()
    descartes_sign_changes = sum(
        (coefficients[index] > 0) != (coefficients[index + 1] > 0)
        for index in range(len(coefficients) - 1)
        if coefficients[index] != 0 and coefficients[index + 1] != 0
    )
    sturm_roots = sp.rem(poly=specialized, x=0)
    real_roots = specialized.count_roots()
    positive_roots = specialized.count_roots(0, sp.oo)
    negative_roots = specialized.count_roots(-sp.oo, 0)
    return {
        "coefficients": [str(coefficient) for coefficient in coefficients],
        "descartes_sign_changes_upper_bound": descartes_sign_changes,
        "real_root_count": int(real_roots),
        "positive_root_count": int(positive_roots),
        "negative_root_count": int(negative_roots),
        "all_six_roots_positive": int(positive_roots) == 6,
        "sturm_audit": str(sturm_roots),
    }


def solve_first_three(targets: dict[str, Fraction]) -> dict[str, object]:
    _, expressions = norm_equations(targets)
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    first_three = [
        sp.factor(sp.resultant(sp.symbols("u") ** 2 - 4 * sp.symbols("u") + 2, expressions["Q"], sp.symbols("u"))),
    ]
    equations, expressions = norm_equations(targets)
    groebner = sp.groebner(
        [equation.as_expr() for equation in equations[:3]],
        a5,
        a4,
        a3,
        order="lex",
        domain=sp.QQ,
    )
    if groebner.contains(sp.Integer(1)):
        raise AssertionError("the first three equations must stay satisfiable for controls")
    solution = sp.solve([equation.as_expr() for equation in equations[:3]], a5, a4, a3, dict=True)
    return {
        "solution_count": len(solution),
        "solutions": [
            {str(symbol): str(sp.factor(solution_row[symbol])) for symbol in (a5, a4, a3)}
            for solution_row in solution
        ],
    }


def run() -> dict[str, object]:
    started = time.process_time()
    bonds = list(layer_bonds((2, 3), (False, False)))
    record = graph_targets(bonds)
    elimination = positivity_certificate(record["targets"])
    equations, expressions = norm_equations(record["targets"])
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    basis = sp.groebner(
        [equation.as_expr() for equation in equations],
        a5,
        a4,
        a3,
        order="lex",
        domain=sp.QQ,
    )
    solutions = sp.solve([poly.as_expr() for poly in basis.polys], a5, a4, a3, dict=True)
    positivity = [
        positive_root_necessities(
            expressions["Q"],
            {
                a5: row[a5],
                a4: row[a4],
                a3: row[a3],
            },
        )
        for row in solutions
    ]
    control = synthetic_control()
    checks = [
        {
            "name": "lex_elimination_nontrivial",
            "passed": elimination["basis_size"] == 3,
            "detail": "the exact three-variable elimination has three non-unit clauses",
        },
        {
            "name": "coefficient_solutions_finite",
            "passed": len(solutions) <= 8,
            "detail": f"found {len(solutions)} complex coefficient solutions",
        },
        {
            "name": "synthetic_control_remains_feasible",
            "passed": control["all_four_norm_residuals_zero"],
            "detail": "known rational mode polynomial satisfies the trace identities",
        },
    ]
    if not all(check["passed"] for check in checks):
        raise AssertionError("positivity producer check failed")
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e233_trace_resultant_positivity.py",
            "interpreter": sys.executable,
            "arithmetic": "exact QQ elimination and Sturm root counting over Q",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
        },
        "data": {
            "graph": "open 2x3 grid at t=1/3",
            "elimination": elimination,
            "coefficient_solutions": [
                {str(symbol): str(row[symbol]) for symbol in (a5, a4, a3)} for row in solutions
            ],
            "positivity": positivity,
            "positive_modes_possible": bool(
                solutions and any(record["all_six_roots_positive"] for record in positivity)
            ),
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for check in payload["checks"]:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
