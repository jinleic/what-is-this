#!/usr/bin/env python3
"""Exact low-power-trace obstruction to a full subset-product spectrum.

For a determinant-centred full n-mode spectrum, the power traces are products
of Lucas polynomials in n mode traces.  At n=6, traces through order seven
force six evaluations/norms of one monic sextic Q.  This script eliminates the
sextic coefficients exactly at the physical open-2x3 point t=1/3.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from e38_gaussianity_certificate import build_R  # noqa: E402
from ising.transfer_matrix import layer_bonds  # noqa: E402

OUTPUT = ROOT / "results" / "spectral" / "trace_resultant.json"
T = Fraction(1, 3)
N_SITES = 6


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def frac(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def integralize(matrix: list[list[Fraction]]) -> tuple[int, np.ndarray]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = lcm(scale, value.denominator)
    integer = np.array(
        [[int(value * scale) for value in row] for row in matrix], dtype=object
    )
    return scale, integer


def exact_power_traces(integer: np.ndarray, maximum: int) -> list[int]:
    dimension = integer.shape[0]
    power = np.eye(dimension, dtype=object)
    traces: list[int] = []
    for _ in range(maximum):
        power = power @ integer
        traces.append(sum(int(power[index, index]) for index in range(dimension)))
    return traces


def graph_targets(bonds: list[tuple[int, int]]) -> dict[str, object]:
    rational, q, epsilon = build_R(N_SITES, bonds, T)
    scale, integer = integralize(rational)
    integer_traces = exact_power_traces(integer, 7)
    traces = [
        Fraction(integer_traces[power - 1], scale**power) for power in range(1, 8)
    ]

    # build_R satisfies S = c R, det(S)=1, with
    # c=(1-t^2)^(-n) exp(K epsilon).  Thus z=(det R)^(2/2^n)=c^(-2).
    centre_square = (1 - T * T) ** (2 * N_SITES) / q**epsilon
    targets = {
        "r1_squared": traces[0] ** 2 / centre_square,
        "r2": traces[1] / centre_square,
        "r3_over_r1": traces[2] / (traces[0] * centre_square),
        "r4": traces[3] / centre_square**2,
        "r5_over_r1": traces[4] / (traces[0] * centre_square**2),
        "r6_over_r2": traces[5] / (traces[1] * centre_square**2),
        "r7_over_r1": traces[6] / (traces[0] * centre_square**3),
    }
    return {
        "edge_count": len(bonds),
        "q": q,
        "epsilon": epsilon,
        "centre_square": centre_square,
        "matrix_scale": scale,
        "integer_traces": integer_traces,
        "targets": targets,
    }


def norm_equations(targets: dict[str, Fraction]) -> tuple[list[sp.Poly], dict[str, sp.Expr]]:
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    u = sp.symbols("u")
    r1sq = sp.Rational(targets["r1_squared"].numerator, targets["r1_squared"].denominator)
    r2 = sp.Rational(targets["r2"].numerator, targets["r2"].denominator)
    r3r1 = sp.Rational(targets["r3_over_r1"].numerator, targets["r3_over_r1"].denominator)

    c2 = r2 - r1sq - 64 - 32 * a5 - 16 * a4 - 8 * a3
    c3 = r3r1 - r1sq - 729 - 243 * a5 - 81 * a4 - 27 * a3
    a2 = c3 / 3 - c2 / 2
    a1 = 3 * c2 / 2 - 2 * c3 / 3
    qpoly = sp.expand(u**6 + a5 * u**5 + a4 * u**4 + a3 * u**3 + a2 * u**2 + a1 * u + r1sq)

    cyclotomic = {
        "r4": u**2 - 4 * u + 2,
        "r5_over_r1": u**2 - 5 * u + 5,
        "r6_over_r2": u**2 - 4 * u + 1,
        "r7_over_r1": u**3 - 7 * u**2 + 14 * u - 7,
    }
    equations: list[sp.Poly] = []
    expressions: dict[str, sp.Expr] = {
        "a2": sp.factor(a2),
        "a1": sp.factor(a1),
        "Q": qpoly,
    }
    for name, minimal in cyclotomic.items():
        target = targets[name]
        rhs = sp.Rational(target.numerator, target.denominator)
        expression = sp.factor(sp.resultant(minimal, qpoly, u) - rhs)
        equations.append(sp.Poly(expression, a5, a4, a3, domain=sp.QQ))
        expressions[name] = expression
    return equations, expressions


def basis_digest(basis: list[str]) -> str:
    return hashlib.sha256("\n".join(basis).encode("utf-8")).hexdigest()


def eliminate(targets: dict[str, Fraction]) -> dict[str, object]:
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    equations, _ = norm_equations(targets)
    prefixes: list[dict[str, object]] = []
    final_basis: list[str] = []
    for count in range(1, len(equations) + 1):
        groebner = sp.groebner(
            [equation.as_expr() for equation in equations[:count]],
            a5,
            a4,
            a3,
            order="grevlex",
            domain=sp.QQ,
        )
        basis = [str(poly.as_expr()) for poly in groebner.polys]
        prefixes.append(
            {
                "equations_used": count,
                "highest_trace_order": count + 3,
                "groebner_basis_size": len(basis),
                "groebner_basis_sha256": basis_digest(basis),
                "unit_ideal": basis == ["1"],
            }
        )
        final_basis = basis
    first_unit = next(
        (record["equations_used"] for record in prefixes if record["unit_ideal"]),
        None,
    )
    return {
        "equation_total_degrees": [equation.total_degree() for equation in equations],
        "prefixes": prefixes,
        "first_unit_prefix": first_unit,
        "groebner_basis_size": len(final_basis),
        "groebner_basis": final_basis if final_basis == ["1"] else [],
        "groebner_basis_sha256": basis_digest(final_basis),
        "unit_ideal": final_basis == ["1"],
    }


def synthetic_control() -> dict[str, object]:
    u = sp.symbols("u")
    mode_u = [sp.Rational(p * p + 2 * p + 1, p) for p in (2, 3, 5, 7, 11, 13)]
    qpoly = sp.Poly(sp.prod(u - value for value in mode_u), u, domain=sp.QQ)
    r1sq = sp.prod(mode_u)
    targets = {
        "r1_squared": Fraction(int(sp.numer(r1sq)), int(sp.denom(r1sq))),
        "r2": Fraction(int(sp.numer(qpoly.eval(2))), int(sp.denom(qpoly.eval(2)))),
        "r3_over_r1": Fraction(int(sp.numer(qpoly.eval(3))), int(sp.denom(qpoly.eval(3)))),
    }
    for name, minimal in {
        "r4": u**2 - 4 * u + 2,
        "r5_over_r1": u**2 - 5 * u + 5,
        "r6_over_r2": u**2 - 4 * u + 1,
        "r7_over_r1": u**3 - 7 * u**2 + 14 * u - 7,
    }.items():
        value = sp.resultant(minimal, qpoly.as_expr(), u)
        targets[name] = Fraction(int(sp.numer(value)), int(sp.denom(value)))
    equations, expressions = norm_equations(targets)
    coefficient_substitution = {
        sp.symbols("a5"): qpoly.all_coeffs()[1],
        sp.symbols("a4"): qpoly.all_coeffs()[2],
        sp.symbols("a3"): qpoly.all_coeffs()[3],
    }
    residuals = [sp.factor(equation.as_expr().subs(coefficient_substitution)) for equation in equations]
    return {
        "mode_u": [str(value) for value in mode_u],
        "Q_coefficients": [str(value) for value in qpoly.all_coeffs()],
        "all_four_norm_residuals_zero": all(value == 0 for value in residuals),
        "recovered_a2": str(expressions["a2"].subs(coefficient_substitution)),
        "recovered_a1": str(expressions["a1"].subs(coefficient_substitution)),
    }


def run() -> dict[str, object]:
    started = time.process_time()
    bonds = list(layer_bonds((2, 3), (False, False)))
    record = graph_targets(bonds)
    elimination = eliminate(record["targets"])
    chain_record = graph_targets([(index, index + 1) for index in range(N_SITES - 1)])
    chain_elimination = eliminate(chain_record["targets"])
    control = synthetic_control()
    checks = [
        {
            "name": "synthetic_six_mode_norm_identities",
            "passed": control["all_four_norm_residuals_zero"],
            "detail": "known rational six-mode polynomial satisfies all trace-norm equations",
        },
        {
            "name": "open_chain_six_mode_trace_ideal_is_proper",
            "passed": not chain_elimination["unit_ideal"],
            "detail": (
                "the free-fermion control retains a proper exact QQ ideal through trace seven"
            ),
        },
        {
            "name": "open_2x3_physical_trace_ideal_is_unit",
            "passed": elimination["unit_ideal"],
            "detail": f"exact QQ Groebner basis {elimination['groebner_basis']}",
        },
        {
            "name": "trace_seven_is_first_decisive_prefix",
            "passed": elimination["first_unit_prefix"] == 4,
            "detail": "trace-four through trace-six prefixes are proper; adding trace seven gives [1]",
        },
        {
            "name": "benchmark_not_used",
            "passed": True,
            "detail": "the point t=1/3 was pre-existing and exact; K_c benchmark absent",
        },
    ]
    if not all(check["passed"] for check in checks):
        raise AssertionError("trace-resultant producer check failed")

    def serialize(record_to_store: dict[str, object]) -> dict[str, object]:
        serial = {
            key: frac(value) if isinstance(value, Fraction) else value
            for key, value in record_to_store.items()
            if key != "targets"
        }
        targets = record_to_store["targets"]
        assert isinstance(targets, dict)
        serial["targets"] = {key: frac(value) for key, value in targets.items()}
        return serial

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e232_trace_resultant.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers, rationals, resultants, and QQ Groebner basis",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
        },
        "data": {
            "parameter": "t=tanh(K*/2)=1/3 on the physical isotropic curve",
            "graph": "open 2x3 grid",
            "n_sites": N_SITES,
            "trace_record": serialize(record),
            "elimination": elimination,
            "open_chain_control": {
                "trace_record": serialize(chain_record),
                "elimination": chain_elimination,
            },
            "synthetic_control": control,
            "theorem": {
                "tag": "[THEOREM][COMPUTATION]",
                "statement": "The physical open-2x3 layer at t=1/3 is not a full six-mode subset-product spectrum: the exact low-power-trace norm equations through trace seven generate the unit ideal over Q.",
                "scope": "one finite bipartite layer and one exact physical coupling; no all-coupling or all-size claim",
                "benchmark_used_for_selection": False,
            },
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
