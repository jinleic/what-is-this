"""A determinant-center-free spectral invariant on the isotropic curve.

For N=2^n and M=N/2, every invertible n-mode subset-product multiset obeys

    e_N^(M-1) e_1^M - e_(N-1)^M = 0.

The physical determinant-one transfer matrix reduces this invariant to a power
difference of the alignment polynomial Q_G and its reversal.  Exact finite graph
rows below audit the algebra; the non-bipartite all-size conclusion is proved
without finite extrapolation in proofs/isotropic_invariant.md.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from math import prod

from e200_isotropic_puiseux import (
    SAMPLE_GRAPHS,
    alignment_counts,
    fraction_record,
    linear_power,
    physical_clear,
    poly_add,
    poly_eval,
    poly_mul,
    poly_pow,
    poly_scale,
    poly_sub,
    reverse_alignment,
    trace_skew_from_spins,
    trim,
)

Polynomial = list[int]
Edge = tuple[int, int]


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def polynomial_digest(poly: Polynomial) -> str:
    payload = ",".join(str(coefficient) for coefficient in trim(poly)).encode()
    return hashlib.sha256(payload).hexdigest()


def subset_product_slots(modes: int) -> list[Fraction]:
    base = Fraction(2, 3)
    factors = [Fraction(index + 3, index + 2) for index in range(modes)]
    slots: list[Fraction] = []
    for mask in range(1 << modes):
        value = base
        for index, factor in enumerate(factors):
            if (mask >> index) & 1:
                value *= factor
        slots.append(value)
    return slots


def elementary_edge_data(slots: list[Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    determinant = prod(slots, start=Fraction(1))
    first = sum(slots, start=Fraction(0))
    penultimate = determinant * sum((1 / value for value in slots), start=Fraction(0))
    return first, penultimate, determinant


def spectral_invariant(slots: list[Fraction]) -> Fraction:
    if len(slots) % 2:
        raise ValueError("the complement invariant requires an even number of slots")
    half_dimension = len(slots) // 2
    first, penultimate, determinant = elementary_edge_data(slots)
    return determinant ** (half_dimension - 1) * first**half_dimension - penultimate**half_dimension


def geometric_power_factor(left: Polynomial, right: Polynomial, exponent: int) -> Polynomial:
    """Return sum_{j=0}^{M-1} left^(M-1-j) right^j."""
    result = [0]
    for right_power in range(exponent):
        term = poly_mul(poly_pow(left, exponent - 1 - right_power), poly_pow(right, right_power))
        result = poly_add(result, term)
    return trim(result)


def lollipop_graph(tail_edges: int) -> tuple[int, tuple[Edge, ...]]:
    if tail_edges < 0:
        raise ValueError("tail length must be nonnegative")
    edges: list[Edge] = [(0, 1), (1, 2), (2, 0)]
    previous = 0
    for offset in range(tail_edges):
        vertex = 3 + offset
        edges.append((previous, vertex))
        previous = vertex
    return 3 + tail_edges, tuple(edges)


def graph_invariant_row(label: str, sites: int, edges: tuple[Edge, ...]) -> tuple[dict[str, object], list[bool]]:
    edge_count = len(edges)
    dimension = 1 << sites
    half_dimension = dimension // 2
    alignments = alignment_counts(sites, edges)
    reversed_alignments = reverse_alignment(alignments)
    trace_skew = trace_skew_from_spins(alignments)
    psi = poly_sub(poly_pow(alignments, half_dimension), poly_pow(reversed_alignments, half_dimension))
    factor = geometric_power_factor(alignments, reversed_alignments, half_dimension)
    factorization_holds = poly_mul(trace_skew, factor) == psi
    is_zero = psi == [0]
    degree = None if is_zero else len(psi) - 1
    expected_degree = edge_count * half_dimension
    y0 = Fraction(5, 3)
    psi_at_y0 = poly_eval(psi, y0)
    trace_skew_at_y0 = poly_eval(trace_skew, y0)
    cleared = [0] if is_zero else physical_clear(psi, expected_degree)
    t0 = Fraction(1, 3)
    cleared_relation = is_zero or poly_eval(cleared, t0) == (2 * t0) ** expected_degree * psi_at_y0
    row = {
        "tag": "[COMPUTATION]",
        "label": label,
        "sites": sites,
        "edge_count": edge_count,
        "dimension": dimension,
        "half_dimension_M": half_dimension,
        "spectral_polynomial": "Psi_G(y)=Q_G(y)^M-Q_G^*(y)^M",
        "psi_zero": is_zero,
        "psi_degree": degree,
        "psi_sha256": polynomial_digest(psi),
        "psi_nonzero_coefficient_count": sum(coefficient != 0 for coefficient in psi),
        "psi_constant_coefficient": psi[0],
        "psi_leading_coefficient": psi[-1],
        "psi_at_y_5_over_3": fraction_record(psi_at_y0),
        "trace_skew_at_y_5_over_3": fraction_record(trace_skew_at_y0),
        "physical_cleared_degree": None if is_zero else len(cleared) - 1,
        "physical_cleared_sha256": polynomial_digest(cleared),
        "complex_exception_root_bound": None if is_zero else 2 * expected_degree,
        "nonzero_modulo_physical_curve": not is_zero,
    }
    audits = [
        factorization_holds,
        is_zero == (trace_skew == [0]),
        is_zero or degree == expected_degree,
        is_zero or len(cleared) - 1 == 2 * expected_degree,
        cleared_relation,
        (psi_at_y0 == 0) == is_zero,
        (psi_at_y0 > 0) == (not is_zero),
    ]
    return row, audits


def run_invariant() -> dict[str, object]:
    checks: list[dict[str, object]] = []

    subset_rows: list[dict[str, object]] = []
    for modes in range(1, 5):
        slots = subset_product_slots(modes)
        first, penultimate, determinant = elementary_edge_data(slots)
        invariant = spectral_invariant(slots)
        subset_rows.append(
            {
                "modes": modes,
                "slots": len(slots),
                "half_dimension_M": len(slots) // 2,
                "invariant_zero": invariant == 0,
                "e1_nonzero": first != 0,
                "e_N_minus_1_nonzero": penultimate != 0,
                "determinant_nonzero": determinant != 0,
            }
        )
    _check(
        checks,
        "subset-product complement invariant vanishes exactly in rational audits",
        all(bool(row["invariant_zero"]) for row in subset_rows),
        "checked from all slots for n=1..4, without symbolic parameter cancellation",
    )

    nonsubset_slots = [Fraction(2), Fraction(3), Fraction(1, 5), Fraction(5, 6)]
    nonsubset_first, nonsubset_penultimate, nonsubset_det = elementary_edge_data(nonsubset_slots)
    nonsubset_invariant = spectral_invariant(nonsubset_slots)
    _check(
        checks,
        "the invariant is not an identically zero symmetric polynomial",
        nonsubset_det == 1 and nonsubset_invariant != 0,
        f"det={nonsubset_det}, e1={nonsubset_first}, e3={nonsubset_penultimate}, Phi={nonsubset_invariant}",
    )

    invariant_rows: list[dict[str, object]] = []
    for label, (sites, edges) in SAMPLE_GRAPHS.items():
        row, audits = graph_invariant_row(label, sites, edges)
        invariant_rows.append(row)
        _check(
            checks,
            f"{label}: determinant-center-free graph invariant",
            all(audits),
            (
                f"Psi zero={row['psi_zero']}, degree={row['psi_degree']}, "
                f"complex-root bound={row['complex_exception_root_bound']}"
            ),
        )

    leaf_pairs = (("triangle", "paw"), ("paw", "lollipop_2"))
    for old_label, new_label in leaf_pairs:
        old_sites, old_edges = SAMPLE_GRAPHS[old_label]
        new_sites, new_edges = SAMPLE_GRAPHS[new_label]
        old_q = alignment_counts(old_sites, old_edges)
        old_q_star = reverse_alignment(old_q)
        old_m = 1 << (old_sites - 1)
        new_m = 1 << (new_sites - 1)
        new_q = alignment_counts(new_sites, new_edges)
        new_q_star = reverse_alignment(new_q)
        new_psi = poly_sub(poly_pow(new_q, new_m), poly_pow(new_q_star, new_m))
        old_psi = poly_sub(poly_pow(old_q, old_m), poly_pow(old_q_star, old_m))
        expected_new = poly_mul(
            poly_pow([1, 1], new_m),
            poly_mul(old_psi, poly_add(poly_pow(old_q, old_m), poly_pow(old_q_star, old_m))),
        )
        _check(
            checks,
            f"leaf propagation {old_label} -> {new_label}",
            new_m == 2 * old_m
            and new_q == poly_mul([1, 1], old_q)
            and trim(new_q_star) == poly_mul([1, 1], old_q_star)
            and new_psi == expected_new,
            f"M doubles from {old_m} to {new_m}; Psi degree={len(new_psi) - 1}",
        )

    lollipop_rows: list[dict[str, object]] = []
    triangle_q = [0, 6, 0, 2]
    triangle_q_star = list(reversed(triangle_q))
    triangle_skew = poly_scale(poly_pow([-1, 1], 3), 2)
    for tail_edges in range(0, 6):
        sites, edges = lollipop_graph(tail_edges)
        alignments = alignment_counts(sites, edges)
        expected_q = poly_mul(poly_pow([1, 1], tail_edges), triangle_q)
        expected_q_star = poly_mul(poly_pow([1, 1], tail_edges), triangle_q_star)
        trace_skew = trace_skew_from_spins(alignments)
        expected_skew = poly_mul(poly_pow([1, 1], tail_edges), triangle_skew)
        physical_skew = physical_clear(trace_skew, len(edges))
        expected_physical = poly_scale(
            poly_mul(linear_power(-1, 6), linear_power(1, 2 * tail_edges)), 2
        )
        row_passed = (
            alignments == expected_q
            and trim(reverse_alignment(alignments)) == expected_q_star
            and trace_skew == expected_skew
            and physical_skew == expected_physical
        )
        lollipop_rows.append(
            {
                "tail_edges_r": tail_edges,
                "sites": sites,
                "edge_count": len(edges),
                "checks_passed": row_passed,
                "trace_skew_sha256": polynomial_digest(trace_skew),
                "physical_trace_skew_sha256": polynomial_digest(physical_skew),
            }
        )
    _check(
        checks,
        "lollipop leaf recurrence finite lemma tests",
        all(bool(row["checks_passed"]) for row in lollipop_rows),
        "r=0..5 agree with I_r=2(y-1)^3(y+1)^r and J_r=2(t-1)^6(t+1)^(2r)",
    )

    nonbipartite_rows = [row for row in invariant_rows if not bool(row["psi_zero"])]
    bipartite_rows = [row for row in invariant_rows if bool(row["psi_zero"])]
    _check(
        checks,
        "finite graph rows exercise both exact sides of the route boundary",
        len(nonbipartite_rows) > 0
        and len(bipartite_rows) > 0
        and all(bool(row["nonzero_modulo_physical_curve"]) for row in nonbipartite_rows)
        and all(not bool(row["nonzero_modulo_physical_curve"]) for row in bipartite_rows),
        f"nonbipartite rows={len(nonbipartite_rows)}, bipartite rows={len(bipartite_rows)}",
    )

    data = {
        "tag": "[THEOREM]",
        "characteristic_coefficient_invariant": {
            "tag": "[LEMMA]",
            "dimension": "N=2^n, M=N/2",
            "definition": "Phi_N(W)=e_N(W)^(M-1)e_1(W)^M-e_(N-1)(W)^M",
            "necessary_condition": (
                "every invertible full n-mode subset-product spectrum has Phi_N=0; complement "
                "pairing gives e_N=c^M and e_1=c sum(lambda^-1)"
            ),
            "determinant_one_form": "Phi_N(W)=tr(W)^M-tr(W^-1)^M when det(W)=1",
            "exact_slot_audits": subset_rows,
        },
        "physical_graph_reduction": {
            "tag": "[LEMMA]",
            "operator": "V_G(a,b)=exp(a sum_v X_v) exp(b sum_uv Z_u Z_v)",
            "multiplicative_coordinates": "x=exp(a), y=exp(2b)",
            "trace_formulas": (
                "tr(V)=cosh(a)^n y^(-m/2)Q_G(y), "
                "tr(V^-1)=cosh(a)^n y^(-m/2)Q_G^*(y)"
            ),
            "polynomial_representative": "Psi_G(y)=Q_G(y)^M-Q_G^*(y)^M",
            "graph_rows": invariant_rows,
        },
        "curve_nonvanishing": {
            "tag": "[THEOREM]",
            "curve": "F_phys(x,y)=(y-1)x^2-(y+1)=0",
            "exact_witness": "t=1/3 gives (x,y)=(2,5/3) and Psi_G(5/3)>0 for every non-bipartite G",
            "modulus_conclusion": "Psi_G is nonzero modulo F_phys for every non-bipartite G",
            "finite_exception_polynomial": (
                "K_G(t)=(2t)^(mM) Psi_G((1+t^2)/(2t)) is a nonzero integer polynomial "
                "of degree 2mM"
            ),
            "complex_exception_bound": "at most 2mM roots, on the nonsingular multiplicative physical curve",
        },
        "all_size_result": {
            "tag": "[THEOREM]",
            "statement": (
                "for every finite simple non-bipartite graph G and every real physical "
                "0<t<1, the isotropic layer spectrum is not a full n-mode subset-product multiset"
            ),
            "reason": (
                "q(t)>1; the odd-Eulerian expansion gives Q_G(q)>Q_G^*(q)>0, hence "
                "Psi_G(q)>0 while every subset-product spectrum requires Phi_N=0"
            ),
            "bipartite_boundary": (
                "for every bipartite G, Q_G=Q_G^* and this entire trace-reciprocity route "
                "vanishes identically at every coupling"
            ),
        },
        "leaf_extension": {
            "tag": "[LEMMA]",
            "recurrence": (
                "adding one leaf gives Q_new=(y+1)Q, Q_new^*=(y+1)Q^*, M_new=2M, and "
                "Psi_new=(y+1)^(2M) Psi_old (Q^M+(Q^*)^M)"
            ),
            "infinite_family": (
                "the triangle with a tail of r>=1 edges is connected, branching, and has "
                "I_r=2(y-1)^3(y+1)^r and cleared J_r=2(t-1)^6(t+1)^(2r)"
            ),
            "finite_lemma_tests": lollipop_rows,
        },
        "scope": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "the theorem covers all non-bipartite layers and an infinite connected branching "
                "family, but it gives no obstruction for bipartite branching layers such as open "
                "rectangular grids; finite rows are lemma tests only"
            ),
        },
    }
    return {"data": data, "checks": checks}


def main() -> int:
    result = run_invariant()
    for item in result["checks"]:
        print(f"[{'PASS' if item['passed'] else 'FAIL'}] {item['name']}: {item['detail']}")
    if all(bool(item["passed"]) for item in result["checks"]):
        print("PASS e201 isotropic trace invariant")
        return 0
    print("FAIL e201 isotropic trace invariant")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
