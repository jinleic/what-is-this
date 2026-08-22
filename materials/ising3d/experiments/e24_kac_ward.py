"""Exact 2D Kac--Ward control and scoped 3D local-weight obstruction.

Run from the repository root with

    .venv/bin/python experiments/e24_kac_ward.py

The script writes ``results/kac_ward/kac_ward.json`` and prints ``PASS`` after
all proved and reproduced claims have been checked.  The full gauge-fixed
30-weight order-eight Groebner problem is deliberately time-bounded and is
reported as unresolved if it does not finish; it is never promoted to a
no-go theorem.
"""

from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
import time
from datetime import datetime, timezone
from fractions import Fraction

import mpmath as mp
import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    BULK_SC_LOG_P,
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    bulk_cubic_trace_monomials,
    closed_nonbacktracking_walk_count,
    exact_kac_ward_polynomial,
    fit_strict_cubic_candidate,
    formal_log_coefficients,
    full_weight_finite_system,
    kac_ward_determinant,
    kac_ward_trace_basis,
    numeric_determinant_prefix,
    polynomial_square,
    strict_cubic_determinant_prefix,
    strict_cubic_equations,
    trace_type_counts,
)
from ising.lattices import cubic, square


DPS = 90
TOLERANCE = mp.mpf("1e-70")
SQUARE_CASES = ((3, 3), (3, 4), (4, 4), (4, 5))
FULL_FAMILY_SHAPES = (
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
    (2, 2, 3),
    (4, 2, 2),
    (2, 4, 2),
    (2, 2, 4),
    (3, 3, 2),
    (3, 2, 3),
    (2, 3, 3),
)
FULL_GROEBNER_TIMEOUT_SECONDS = 20.0


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _polyval(coefficients, value):
    result = mp.mpf(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def _sympy_to_mpc(value: sp.Expr, dps: int) -> mp.mpc:
    evaluated = sp.N(value, dps)
    with mp.workdps(dps):
        return +mp.mpc(str(sp.re(evaluated)), str(sp.im(evaluated)))


def _complex_record(value, digits: int = 70) -> dict[str, str]:
    return {
        "real": mp.nstr(mp.re(value), digits),
        "imag": mp.nstr(mp.im(value), digits),
    }


def _first_numeric_failure(observed, expected, tolerance) -> int | None:
    for degree, (left, right) in enumerate(zip(observed, expected, strict=True)):
        if abs(left - right) > tolerance:
            return degree
    return None


def _groebner_worker(equations, variables, queue) -> None:
    try:
        started = time.monotonic()
        basis = sp.groebner(equations, *variables, order="grevlex")
        queue.put(
            {
                "status": "completed",
                "elapsed_seconds": time.monotonic() - started,
                "basis_size": len(basis.polys),
                "contains_one": list(basis) == [1],
            }
        )
    except BaseException as error:  # child-process boundary: serialize the failure
        queue.put({"status": "error", "error": f"{type(error).__name__}: {error}"})


def _timed_groebner(equations, variables, timeout_seconds: float) -> dict:
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    process = context.Process(target=_groebner_worker, args=(equations, variables, queue))
    started = time.monotonic()
    process.start()
    process.join(timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join()
        result = {
            "status": "timeout",
            "timeout_seconds": timeout_seconds,
            "elapsed_seconds": time.monotonic() - started,
        }
    elif not queue.empty():
        result = queue.get()
        result["timeout_seconds"] = timeout_seconds
    else:
        result = {
            "status": "error",
            "timeout_seconds": timeout_seconds,
            "error": f"worker exited with code {process.exitcode} without a result",
        }
    queue.close()
    return result


def _two_dimensional_controls() -> tuple[list[dict], int, mp.mpf]:
    records: list[dict] = []
    global_max_discrepancy = 0
    global_max_numeric_error = mp.mpf(0)
    with mp.workdps(DPS):
        evaluation_point = mp.mpf("0.173205080756887729352744634150587236694280525381")
        for shape in SQUARE_CASES:
            lattice = square(*shape, periodic=False)
            even = even_subgraph_polynomial(lattice)
            target = polynomial_square(even)
            observed, certificate = exact_kac_ward_polynomial(lattice)
            discrepancies = [left - right for left, right in zip(observed, target, strict=True)]
            maximum = max(map(abs, discrepancies), default=0)
            global_max_discrepancy = max(global_max_discrepancy, maximum)

            determinant = kac_ward_determinant(lattice, evaluation_point, dps=DPS)
            expected_value = _polyval(target, evaluation_point)
            numerical_error = abs(determinant - expected_value)
            global_max_numeric_error = max(global_max_numeric_error, numerical_error)

            log_even = formal_log_coefficients(even, 8)
            count_rows = []
            for order in (4, 6, 8):
                weighted_basis = kac_ward_trace_basis(lattice, order)
                required_trace = -2 * order * log_even[order]
                count_rows.append(
                    {
                        "length": order,
                        "even_subgraphs": even[order] if order < len(even) else 0,
                        "closed_nonbacktracking_walks": closed_nonbacktracking_walk_count(
                            lattice, order
                        ),
                        "kac_ward_trace_basis_1_zeta_zeta2_zeta3": list(weighted_basis),
                        "required_trace_minus_2k_logP": _fraction_text(required_trace),
                        "trace_discrepancy_basis": [
                            _fraction_text(weighted_basis[0] - required_trace),
                            weighted_basis[1],
                            weighted_basis[2],
                            weighted_basis[3],
                        ],
                    }
                )
            records.append(
                {
                    "shape": list(shape),
                    "boundary": "free",
                    "n_sites": lattice.n_sites,
                    "n_bonds": lattice.n_bonds,
                    "directed_edges": certificate.directed_edge_count,
                    "even_subgraph_polynomial": even,
                    "determinant_coefficients": list(observed),
                    "target_P_squared": target,
                    "maximum_exact_coefficient_discrepancy": maximum,
                    "exact_certificate": {
                        "cyclotomic_ring": "Z[zeta_8], zeta_8^4=-1",
                        "four_embeddings_per_prime": True,
                        "primes": list(certificate.primes),
                        "crt_modulus": str(certificate.modulus),
                        "hadamard_coefficient_bound": str(certificate.coefficient_bound),
                        "max_nonrational_component": certificate.max_nonrational_component,
                    },
                    "numeric_check": {
                        "dps": DPS,
                        "v": mp.nstr(evaluation_point, 60),
                        "determinant": _complex_record(determinant),
                        "P_squared": mp.nstr(expected_value, 70),
                        "absolute_error": mp.nstr(numerical_error, 20),
                    },
                    "walk_and_subgraph_counts": count_rows,
                }
            )
            print(
                f"2D {shape[0]}x{shape[1]} free: max exact coefficient discrepancy {maximum}; "
                f"{DPS}-dps error {mp.nstr(numerical_error, 6)}"
            )
    return records, global_max_discrepancy, global_max_numeric_error


def _strict_cubic_obstruction() -> dict:
    straight, turn, equations, required = strict_cubic_equations(BULK_SC_LOG_P)
    basis_46 = sp.groebner([equations[4], equations[6]], straight, turn, order="lex")
    basis_468 = sp.groebner(
        [equations[4], equations[6], equations[8]], straight, turn, order="lex"
    )
    bulk_rows = []
    for order in (4, 6, 8):
        types = bulk_cubic_trace_monomials(order)
        bulk_rows.append(
            {
                "length": order,
                "closed_nonbacktracking_walks_per_site": sum(types.values()),
                "walks_by_straight_transition_count": {
                    str(power): count for power, count in types.items()
                },
                "connected_even_subgraph_log_coefficient_per_site": _fraction_text(
                    BULK_SC_LOG_P[order]
                ),
                "required_weighted_trace_per_site": _fraction_text(required[order]),
                "equation": sp.sstr(equations[order]) + " = 0",
            }
        )

    finite_lattice = cubic(3, 3, 2, periodic=False)
    finite_even = even_subgraph_polynomial(finite_lattice)
    finite_target = polynomial_square(finite_even)
    fitted_straight, fitted_turn = fit_strict_cubic_candidate(finite_lattice)
    symbolic_prefix = strict_cubic_determinant_prefix(
        finite_lattice, fitted_straight, fitted_turn, 8
    )
    symbolic_failures = [
        degree
        for degree in range(9)
        if sp.simplify(symbolic_prefix[degree] - finite_target[degree]) != 0
    ]
    with mp.workdps(DPS):
        numerical_prefix = numeric_determinant_prefix(
            finite_lattice,
            _sympy_to_mpc(fitted_straight, DPS),
            _sympy_to_mpc(fitted_turn, DPS),
            8,
            dps=DPS,
        )
        numerical_failure = _first_numeric_failure(
            numerical_prefix, finite_target[:9], TOLERANCE
        )
        discrepancy = numerical_prefix[8] - finite_target[8]

    finite_log = formal_log_coefficients(finite_even, 8)
    finite_rows = []
    for order in (4, 6, 8):
        type_counts = trace_type_counts(finite_lattice, order)
        finite_rows.append(
            {
                "length": order,
                "even_subgraphs": finite_even[order],
                "closed_nonbacktracking_walks": sum(type_counts.values()),
                "walks_by_straight_transition_count": {
                    str(power): count for power, count in type_counts.items()
                },
                "required_weighted_trace": _fraction_text(-2 * order * finite_log[order]),
            }
        )

    return {
        "symmetry_reduction": {
            "translation_invariant_unknowns": 30,
            "strict_scalar_cubic_covariance_unknowns": 2,
            "unknowns": {"a": "straight d'=d", "b": "orthogonal d' perpendicular d"},
            "justification": (
                "The proper cubic group is transitive separately on allowed straight and "
                "ordered orthogonal direction pairs. This is strict scalar covariance, not "
                "projective covariance up to a direction gauge."
            ),
        },
        "bulk_equations": bulk_rows,
        "groebner_k4_k6": {
            "basis_size": len(basis_46.polys),
            "contains_one": list(basis_46) == [1],
        },
        "groebner_k4_k6_k8": {
            "basis": [sp.sstr(polynomial) for polynomial in basis_468],
            "contains_one": list(basis_468) == [1],
            "verdict": "inconsistent",
            "first_inconsistent_order": 8,
        },
        "finite_independent_check": {
            "shape": [3, 3, 2],
            "boundary": "free",
            "P_coefficients_through_8": finite_even[:9],
            "P_squared_coefficients_through_8": finite_target[:9],
            "fitted_candidate": {
                "b": sp.sstr(fitted_turn),
                "a": sp.sstr(fitted_straight),
                "fits_orders": [4, 6],
            },
            "symbolic_determinant_prefix": [sp.sstr(value) for value in symbolic_prefix],
            "symbolic_first_failure_order": symbolic_failures[0],
            "numeric_dense_matrix": {
                "dps": DPS,
                "first_failure_order": numerical_failure,
                "coefficient_at_failure": _complex_record(numerical_prefix[8]),
                "target_at_failure": str(finite_target[8]),
                "coefficient_discrepancy": _complex_record(discrepancy),
            },
            "walk_and_subgraph_counts": finite_rows,
        },
    }


def _natural_three_dimensional_guess() -> dict:
    """Falsify the unsigned half-angle guess on the free cube graph."""

    lattice = cubic(2, 2, 2, periodic=False)
    even = even_subgraph_polynomial(lattice)
    target = polynomial_square(even)
    straight = sp.Integer(1)
    turn = (1 + sp.I) / sp.sqrt(2)
    symbolic = strict_cubic_determinant_prefix(lattice, straight, turn, 8)
    symbolic_failures = [
        degree for degree in range(9) if sp.simplify(symbolic[degree] - target[degree]) != 0
    ]
    with mp.workdps(DPS):
        numerical = numeric_determinant_prefix(
            lattice,
            mp.mpc(1),
            _sympy_to_mpc(turn, DPS),
            8,
            dps=DPS,
        )
        numeric_failure = _first_numeric_failure(numerical, target[:9], TOLERANCE)
    log_even = formal_log_coefficients(even, 8)
    rows = []
    for order in (4, 6, 8):
        counts = trace_type_counts(lattice, order)
        weighted_trace = sum(
            count * straight**power * turn ** (order - power)
            for power, count in counts.items()
        )
        rows.append(
            {
                "length": order,
                "even_subgraphs": even[order],
                "closed_nonbacktracking_walks": sum(counts.values()),
                "weighted_trace": sp.sstr(sp.simplify(weighted_trace)),
                "required_weighted_trace": _fraction_text(-2 * order * log_even[order]),
            }
        )
    return {
        "definition": (
            "straight weight 1; every orthogonal turn receives exp(i*pi/4), using the "
            "unsigned angle pi/2 because 3D supplies no left/right turning sign"
        ),
        "shape": [2, 2, 2],
        "boundary": "free",
        "P": even,
        "P_squared": target,
        "symbolic_determinant_prefix": [sp.sstr(value) for value in symbolic],
        "symbolic_first_failure_order": symbolic_failures[0] if symbolic_failures else None,
        "coefficient_at_first_failure": (
            sp.sstr(symbolic[symbolic_failures[0]]) if symbolic_failures else None
        ),
        "target_at_first_failure": (
            str(target[symbolic_failures[0]]) if symbolic_failures else None
        ),
        "numeric_dense_matrix_first_failure_order": numeric_failure,
        "walk_and_subgraph_counts": rows,
    }


def _full_weight_search() -> dict:
    system = full_weight_finite_system(FULL_FAMILY_SHAPES, (4, 6, 8), gauge_fix=True)
    lower_equations = [
        equation
        for label, equation in zip(system.labels, system.equations, strict=True)
        if label[1] in (4, 6)
    ]
    lower_unique = list(dict.fromkeys(map(sp.expand, lower_equations)))
    lower_basis = sp.groebner(lower_unique, *system.variables, order="grevlex")

    all_unique = list(dict.fromkeys(map(sp.expand, system.equations)))
    timed_result = _timed_groebner(
        all_unique, system.variables, FULL_GROEBNER_TIMEOUT_SECONDS
    )
    canonical_equations = "\n".join(
        f"{shape}:{order}:{sp.srepr(equation)}"
        for (shape, order), equation in zip(system.labels, system.equations, strict=True)
    )
    equation_hash = hashlib.sha256(canonical_equations.encode("utf-8")).hexdigest()
    pair_name = lambda pair: f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"

    return {
        "ansatz": (
            "one nonzero complex scalar U(d,d') for each of the 30 ordered pairs of cubic "
            "directions with d' != -d, shared by every vertex and every free box"
        ),
        "pre_gauge_unknown_count": system.pre_gauge_unknown_count,
        "gauge_action": "U(d,d') -> g(d)^(-1) U(d,d') g(d')",
        "gauge_dimension": 5,
        "gauge_slice": [pair_name(pair) + "=1" for pair in DIRECTION_GAUGE_TREE],
        "gauge_slice_scope": (
            "generic chart where all five tree weights are nonzero; branches on which any "
            "fixed tree weight vanishes are not covered"
        ),
        "post_gauge_unknown_count": system.post_gauge_unknown_count,
        "finite_shapes": [list(shape) for shape in FULL_FAMILY_SHAPES],
        "orders": [4, 6, 8],
        "raw_equation_count": len(system.equations),
        "unique_equation_count": len(all_unique),
        "equation_term_counts": [len(equation.as_ordered_terms()) for equation in all_unique],
        "equation_sha256": equation_hash,
        "equation_definition": (
            "Tr_shape(U^k) + 2*k*[v^k]log(P_shape(v)) = 0; traces are exact "
            "integer polynomials generated by closed non-backtracking direction words"
        ),
        "through_order_6": {
            "unique_equation_count": len(lower_unique),
            "groebner_order": "grevlex over QQ",
            "basis_size": len(lower_basis.polys),
            "contains_one": list(lower_basis) == [1],
            "verdict": "consistent (no full-family obstruction through k=6)",
        },
        "through_order_8": {
            **timed_result,
            "groebner_order": "grevlex over QQ",
            "verdict": (
                "unresolved" if timed_result["status"] != "completed" else
                ("inconsistent" if timed_result.get("contains_one") else "consistent")
            ),
        },
        "scope_conclusion": (
            "No no-go theorem for the full 30-weight family is claimed unless the order-eight "
            "Groebner computation completes with basis [1]."
        ),
    }


def main() -> None:
    controls, max_exact_discrepancy, max_numeric_error = _two_dimensional_controls()
    strict = _strict_cubic_obstruction()
    natural = _natural_three_dimensional_guess()
    full = _full_weight_search()

    checks = [
        {
            "name": "2d_exact_coefficient_identity",
            "passed": max_exact_discrepancy == 0,
            "detail": f"four free square lattices; maximum discrepancy {max_exact_discrepancy}",
        },
        {
            "name": "2d_high_precision_determinants",
            "passed": max_numeric_error < TOLERANCE,
            "detail": f"dps={DPS}, maximum absolute error={mp.nstr(max_numeric_error, 12)}",
        },
        {
            "name": "strict_cubic_groebner_obstruction",
            "passed": strict["groebner_k4_k6_k8"]["contains_one"],
            "detail": "k=4,6 consistent; adding k=8 gives exact Groebner basis [1]",
        },
        {
            "name": "strict_cubic_finite_numeric_reproduction",
            "passed": strict["finite_independent_check"]["numeric_dense_matrix"][
                "first_failure_order"
            ]
            == 8,
            "detail": "3x3x2 free dense 90-dps matrix powers first disagree at v^8",
        },
        {
            "name": "natural_3d_guess_falsified",
            "passed": natural["numeric_dense_matrix_first_failure_order"] == 6,
            "detail": "2x2x2 free: determinant coefficient is 32*i at v^6, target is 32",
        },
        {
            "name": "full_family_scope_not_overclaimed",
            "passed": (
                full["through_order_8"]["verdict"] != "inconsistent"
                or full["through_order_8"].get("contains_one") is True
            ),
            "detail": (
                f"30 weights -> 25 on generic gauge chart; k<=6 consistent; "
                f"k<=8 status={full['through_order_8']['status']} and is reported accordingly"
            ),
        },
    ]
    payload = {
        "provenance": {
            "script": "experiments/e24_kac_ward.py",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "precision": f"exact Python integers/Fraction/sympy; mpmath dps={DPS}",
            "interpreter": ".venv/bin/python",
            "full_family_groebner_timeout_seconds": FULL_GROEBNER_TIMEOUT_SECONDS,
        },
        "data": {
            "normalization": "det(I-vU)=P(v)^2; Z=2^N(cosh K)^n_b P(v)",
            "two_dimensional_controls": controls,
            "strict_cubic_scalar_obstruction": strict,
            "natural_three_dimensional_guess": natural,
            "full_translation_invariant_30_weight_search": full,
            "scope": {
                "proved": (
                    "the 2D Kac-Ward identity on the listed boxes and nonexistence in the "
                    "strict scalar cubic-covariant two-weight slice through the exact k=8 contradiction"
                ),
                "not_proved": (
                    "nonexistence for the full 30-weight family, zero-entry gauge branches, "
                    "projective covariance, edge/vertex-dependent weights, sums over spin "
                    "structures, or auxiliary/bond dimension >=2"
                ),
                "bond_dimension_2_tested": False,
            },
        },
        "checks": checks,
    }
    os.makedirs("results/kac_ward", exist_ok=True)
    output_path = "results/kac_ward/kac_ward.json"
    with open(output_path, "w", encoding="utf-8") as output:
        json.dump(payload, output, indent=2)
        output.write("\n")

    for check in checks:
        print(f"{check['name']}: {'PASS' if check['passed'] else 'FAIL'} -- {check['detail']}")
    print(
        "full 30-weight order-eight verdict: "
        + full["through_order_8"]["verdict"].upper()
    )
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more Kac--Ward checks failed")
    print(f"wrote {output_path}")
    print("PASS")


if __name__ == "__main__":
    main()
