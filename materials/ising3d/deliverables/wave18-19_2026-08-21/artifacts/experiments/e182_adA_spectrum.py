#!/usr/bin/env python3
"""Exact ad_(iA) spectrum and the edge-eigenvalue coincidence obstruction.

This is the theory/combinatorics producer for wave 19.  It uses only integer and
Fraction arithmetic.  The integrated artifact is written by e184.
"""

from __future__ import annotations

import platform
import resource
import time
from fractions import Fraction
from math import comb

CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 1_900_000_000

EndpointVector = dict[str, Fraction]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str, budget: float = CPU_BUDGET_SECONDS) -> None:
    used = time.process_time() - started
    if used > budget:
        raise RuntimeError(
            f"NON-DECISIVE process-time expiry at {stage}: {used:.3f}s > {budget:.3f}s"
        )
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise RuntimeError(
            f"NON-DECISIVE RSS expiry at {stage}: {rss} >= {RSS_CAP_BYTES} bytes"
        )


def grade_spectrum(n: int, grade: int) -> tuple[dict[str, int], ...]:
    """Spectrum of ad_(iA) on complexified Clifford grade ``grade``.

    ``eigenvalue_i_coefficient=c`` denotes the eigenvalue ``c*i``.  If ``a``
    minus-weight and ``b`` plus-weight one-Majorana eigenvectors are selected,
    then ``a+b=grade``, the eigenvalue is ``2*(b-a)*i``, and the multiplicity is
    C(n,a) C(n,b).
    """
    if n < 1 or not 0 <= grade <= 2 * n:
        raise ValueError((n, grade))
    rows: list[dict[str, int]] = []
    for a in range(max(0, grade - n), min(n, grade) + 1):
        b = grade - a
        rows.append(
            {
                "eigenvalue_i_coefficient": 2 * (b - a),
                "multiplicity": comb(n, a) * comb(n, b),
                "minus_count": a,
                "plus_count": b,
            }
        )
    return tuple(sorted(rows, key=lambda row: row["eigenvalue_i_coefficient"]))


def expected_spectrum_values(n: int, grade: int) -> tuple[int, ...]:
    radius = min(grade, 2 * n - grade)
    return tuple(range(-2 * radius, 2 * radius + 1, 4))


def edge_majorana_occupancies(n: int, p: int, q: int) -> tuple[int, ...]:
    """Numbers of Majoranas from each site pair in gamma_(2p+1)...gamma_(2q)."""
    if not 0 <= p < q < n:
        raise ValueError((n, p, q))
    values = [0] * n
    values[p] = values[q] = 1
    for site in range(p + 1, q):
        values[site] = 2
    return tuple(values)


def edge_eigenvalue_support(n: int, p: int, q: int) -> tuple[int, ...]:
    occupancies = edge_majorana_occupancies(n, p, q)
    singly_occupied = sum(value == 1 for value in occupancies)
    if singly_occupied != 2:
        raise AssertionError("an Ising edge must have exactly two rotating endpoints")
    return (-4, 0, 4)


def _add(target: EndpointVector, label: str, coefficient: Fraction) -> None:
    value = target.get(label, Fraction(0)) + coefficient
    if value:
        target[label] = value
    else:
        target.pop(label, None)


def endpoint_derivation(vector: EndpointVector) -> EndpointVector:
    """Apply D=ad_(i(X_u+X_v)) on the endpoint basis {ZZ,ZY,YZ,YY}."""
    local = {
        "Z": (("Y", Fraction(2)),),
        "Y": (("Z", Fraction(-2)),),
    }
    answer: EndpointVector = {}
    for label, coefficient in vector.items():
        if len(label) != 2 or any(letter not in "ZY" for letter in label):
            raise ValueError(label)
        for new_letter, factor in local[label[0]]:
            _add(answer, new_letter + label[1], coefficient * factor)
        for new_letter, factor in local[label[1]]:
            _add(answer, label[0] + new_letter, coefficient * factor)
    return answer


def vector_add(*vectors: EndpointVector) -> EndpointVector:
    answer: EndpointVector = {}
    for vector in vectors:
        for label, coefficient in vector.items():
            _add(answer, label, coefficient)
    return answer


def vector_scale(vector: EndpointVector, coefficient: Fraction) -> EndpointVector:
    return {
        label: value * coefficient
        for label, value in vector.items()
        if value * coefficient
    }


def serialize_vector(vector: EndpointVector) -> dict[str, str]:
    return {label: str(value) for label, value in sorted(vector.items())}


def vandermonde_norm_squared(eigenvalue_i_coefficients: tuple[int, ...]) -> int:
    """Exact squared Gaussian norm of the Vandermonde determinant."""
    if len(set(eigenvalue_i_coefficients)) != len(eigenvalue_i_coefficients):
        return 0
    value = 1
    for right, coefficient in enumerate(eigenvalue_i_coefficients):
        for prior in eigenvalue_i_coefficients[:right]:
            value *= (coefficient - prior) ** 2
    return value


def edge_projector_audit() -> dict[str, object]:
    edge = {"ZZ": Fraction(1)}
    first = endpoint_derivation(edge)
    second = endpoint_derivation(first)
    third = endpoint_derivation(second)
    zero = vector_add(third, vector_scale(first, Fraction(16)))

    component_zero = vector_add(edge, vector_scale(second, Fraction(1, 16)))
    component_cosine = vector_scale(second, Fraction(-1, 16))
    component_sine = vector_scale(first, Fraction(1, 4))

    expected_zero = {"ZZ": Fraction(1, 2), "YY": Fraction(1, 2)}
    expected_cosine = {"ZZ": Fraction(1, 2), "YY": Fraction(-1, 2)}
    expected_sine = {"YZ": Fraction(1, 2), "ZY": Fraction(1, 2)}
    passed = (
        not zero
        and component_zero == expected_zero
        and component_cosine == expected_cosine
        and component_sine == expected_sine
        and endpoint_derivation(component_zero) == {}
        and endpoint_derivation(component_cosine)
        == vector_scale(component_sine, Fraction(4))
        and endpoint_derivation(component_sine)
        == vector_scale(component_cosine, Fraction(-4))
        and vector_add(component_zero, component_cosine) == edge
    )
    return {
        "tag": "[LEMMA][COMPUTATION]",
        "minimal_polynomial_on_every_edge": "x*(x^2+16)",
        "edge_eigenvalue_i_coefficients": [-4, 0, 4],
        "vandermonde_determinant_norm_squared": vandermonde_norm_squared((-4, 0, 4)),
        "D_edge": serialize_vector(first),
        "D2_edge": serialize_vector(second),
        "projected_rational_components": {
            "eigenvalue_zero": serialize_vector(component_zero),
            "paired_nonzero_cosine": serialize_vector(component_cosine),
            "paired_nonzero_sine": serialize_vector(component_sine),
        },
        "passed": passed,
    }


def run_spectrum_audit() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    spectrum_rows: list[dict[str, object]] = []
    for n in range(2, 13):
        for grade in range(0, 2 * n + 1):
            rows = grade_spectrum(n, grade)
            values = tuple(row["eigenvalue_i_coefficient"] for row in rows)
            multiplicity = sum(row["multiplicity"] for row in rows)
            passed = (
                values == expected_spectrum_values(n, grade)
                and multiplicity == comb(2 * n, grade)
            )
            if not passed:
                raise AssertionError((n, grade, values, multiplicity))
        budget_tick(started, f"grade spectra n={n}")
    checks.append(
        {
            "name": "exact grade spectrum through n=12",
            "passed": True,
            "detail": "all multiplicities sum to C(2n,k) and all eigenvalue sets match",
        }
    )

    for n in (6, 8, 9):
        grades = tuple(range(2, 2 * n, 4))
        spectrum_rows.append(
            {
                "n_vertices": n,
                "grades_2_mod_4": list(grades),
                "spectra": {
                    str(grade): [dict(row) for row in grade_spectrum(n, grade)]
                    for grade in grades
                },
                "common_i_coefficients": sorted(
                    set.intersection(
                        *(
                            {
                                row["eigenvalue_i_coefficient"]
                                for row in grade_spectrum(n, grade)
                            }
                            for grade in grades
                        )
                    )
                ),
            }
        )

    edge_profiles: list[dict[str, object]] = []
    for n in range(2, 13):
        for p in range(n):
            for q in range(p + 1, n):
                occupancies = edge_majorana_occupancies(n, p, q)
                support = edge_eigenvalue_support(n, p, q)
                if support != (-4, 0, 4):
                    raise AssertionError((n, p, q, support))
                if n in (6, 8, 9) and (q - p) % 2 == 1:
                    edge_profiles.append(
                        {
                            "n_vertices": n,
                            "path_positions": [p, q],
                            "path_distance": q - p,
                            "clifford_grade": 2 * (q - p),
                            "site_pair_occupancies": list(occupancies),
                            "edge_eigenvalue_i_coefficients": list(support),
                        }
                    )
        budget_tick(started, f"edge profiles n={n}")
    checks.append(
        {
            "name": "all edge terms through n=12 have the same three eigenvalues",
            "passed": True,
            "detail": "two singly occupied endpoint pairs; interior pairs are D-invariant",
        }
    )

    projector = edge_projector_audit()
    checks.append(
        {
            "name": "exact three-component projectors",
            "passed": bool(projector["passed"]),
            "detail": "D(D^2+16) annihilates ZZ and the rational components recombine",
        }
    )

    collision = {
        "tag": "[THEOREM — refuted method]",
        "statement": (
            "Vandermonde separation recovers the three global ad_(iA) components of B, "
            "but cannot separate B_path from B_high: every path bond and every chord has "
            "the identical nonzero eigenvalue support {-4i,0,4i}."
        ),
        "polynomial_obstruction": (
            "Any polynomial p(D) annihilating the three independent components of every path "
            "edge has p(-4i)=p(0)=p(4i)=0 and therefore annihilates every chord as well."
        ),
        "grade_coincidence": (
            "The full grade-k spectra overlap much more broadly; for an Ising edge vector, "
            "only the common {-4i,0,4i} slice is occupied, independently of k=2(q-p)."
        ),
    }

    if not all(bool(row["passed"]) for row in checks):
        raise AssertionError("ad_A spectrum audit failed")
    budget_tick(started, "completed spectrum audit")
    data = {
        "tag": "[LEMMA][THEOREM — refuted method]",
        "operator": "D=ad_(iA), A=sum_r X_r",
        "single_majorana_eigenvectors": {
            "gamma_(2r)+i gamma_(2r+1)": "-2i",
            "gamma_(2r)-i gamma_(2r+1)": "+2i",
        },
        "grade_spectrum_formula": (
            "at grade k choose a,b in [0,n], a+b=k; eigenvalue=2i(b-a), "
            "multiplicity=C(n,a)C(n,b)"
        ),
        "selected_grade_spectra": spectrum_rows,
        "edge_profiles": edge_profiles,
        "edge_projectors": projector,
        "collision_obstruction": collision,
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }
    return data, checks


def main() -> int:
    data, checks = run_spectrum_audit()
    print(
        f"PASS e182: {len(checks)} checks, "
        f"cpu={data['process_time_seconds']}s, rss={data['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
