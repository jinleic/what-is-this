"""Ordering-free parity-projected Gaussian obstruction for the 2x3 layer.

At t = tanh(K*/2) = 1/3, construct the exact rational symmetrized transfer
representative R, restrict it to the physical P = product_v X_v sectors, and
apply the monic integral finite-field pair-product certificate separately to
P=+, P=-, and to their cross product.  No eigenvalue ordering or inertia is
used anywhere in this certificate.
"""

from __future__ import annotations

import json
import os
import platform
import signal
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from e38_gaussianity_certificate import build_R  # noqa: E402
from e107_pair_product_obstruction import (  # noqa: E402
    INT64_MAX,
    coeff_digest,
    extend_power_sums,
    factorize,
    horner_matrix_zero,
    integralize,
    is_prime,
    newton_from_power_sums,
    poly_deriv_desc,
    poly_gcd,
    poly_is_zero,
    poly_rem,
    run_pipeline,
    sequential_power_traces,
    trace_power_binary,
)
from ising.transfer_matrix import layer_bonds  # noqa: E402

N_MODES = 6
SECTOR_DIM = 1 << (N_MODES - 1)
WITHIN_PAIR_DEGREE = SECTOR_DIM * (SECTOR_DIM - 1) // 2
CROSS_PAIR_DEGREE = SECTOR_DIM * SECTOR_DIM
PRIMES = (1_000_003, 2_000_003)
WALL_SECONDS = 300
CHECKS: list[dict[str, object]] = []
FAILURES: list[str] = []


def hard_timeout(_signum, _frame):
    raise TimeoutError(f"NON-DECISIVE: computation exceeded the {WALL_SECONDS}-second wall")


def check(name: str, passed: bool, detail: str) -> bool:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)
    CHECKS.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        FAILURES.append(name)
    return bool(passed)


def within_exponent_class_count(modes: int) -> int:
    """Classes epsilon+delta from two distinct slots of one parity half."""
    return sum(comb(modes, ones) * 2 ** (modes - ones) for ones in range(2, modes + 1, 2))


def cross_exponent_class_count(modes: int) -> int:
    """Classes epsilon+delta from an even slot and an odd slot."""
    return sum(comb(modes, ones) * 2 ** (modes - ones) for ones in range(1, modes + 1, 2))


def threshold_data(modes: int = N_MODES) -> dict[str, int]:
    slots = 1 << (modes - 1)
    within_pairs = slots * (slots - 1) // 2
    cross_pairs = slots * slots
    within_classes = within_exponent_class_count(modes)
    cross_classes = cross_exponent_class_count(modes)
    return {
        "modes": modes,
        "slots_per_parity_half": slots,
        "within_pair_slots": within_pairs,
        "cross_pair_slots": cross_pairs,
        "even_even_max_exponent_classes": within_classes,
        "odd_odd_max_exponent_classes": within_classes,
        "even_odd_max_exponent_classes": cross_classes,
        "even_even_gcd_degree_floor": within_pairs - within_classes,
        "odd_odd_gcd_degree_floor": within_pairs - within_classes,
        "even_odd_gcd_degree_floor": cross_pairs - cross_classes,
    }


def physical_flip_map(n: int) -> list[int]:
    """Computational-basis permutation of P=product_v X_v."""
    mask = (1 << n) - 1
    return [mask ^ state for state in range(1 << n)]


def split_physical_parity(R: list[list[Fraction]]) -> dict[int, list[list[Fraction]]]:
    """Exact coordinate restrictions in b_a^s=e_a+s e_{P a}, a<2^(n-1)."""
    dim = len(R)
    assert dim and dim & (dim - 1) == 0
    half = dim // 2
    mask = dim - 1
    sectors: dict[int, list[list[Fraction]]] = {}
    for sign in (1, -1):
        block = [[Fraction(0)] * half for _ in range(half)]
        for a in range(half):
            pa = mask ^ a
            for b in range(half):
                pb = mask ^ b
                block[a][b] = (R[a][b] + sign * R[a][pb] + sign * R[pa][b] + R[pa][pb]) / 2
        sectors[sign] = block
    return sectors

SYMBOLIC_CLEARING_T_POWER = 3
SYMBOLIC_CLEARING_ONE_PLUS_T2_POWER = 4
SYMBOLIC_ENTRY_DEGREE_BOUND = 26
WEDGE_ENTRY_DEGREE_BOUND = 2 * SYMBOLIC_ENTRY_DEGREE_BOUND
GENERIC_RESULTANT_DEGREE_BOUND = 2 * WEDGE_ENTRY_DEGREE_BOUND * WITHIN_PAIR_DEGREE * (WITHIN_PAIR_DEGREE - 1)
RATIONAL_GRID = (
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
    Fraction(2, 3),
)


def trim_polynomial(coefficients: list[int | Fraction]) -> tuple[int | Fraction, ...]:
    """Canonical ascending coefficient tuple over Q[t]."""
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients or [0])


def polynomial_degree(coefficients: tuple[int | Fraction, ...]) -> int:
    return -1 if coefficients == (0,) else len(coefficients) - 1


def one_plus_t_squared_power(power: int) -> tuple[int, ...]:
    coefficients = [0] * (2 * power + 1)
    for index in range(power + 1):
        coefficients[2 * index] = comb(power, index)
    return tuple(coefficients)


def shifted_scaled_polynomial(coefficients: tuple[int, ...], scale: int, shift: int) -> tuple[int, ...]:
    return tuple([0] * shift + [scale * coefficient for coefficient in coefficients])


def clearing_factor_polynomial() -> tuple[int, ...]:
    """C(t)=8 t^3 (1+t^2)^4, clearing every q(t)^m denominator for m in [-4,3]."""
    return shifted_scaled_polynomial(
        one_plus_t_squared_power(SYMBOLIC_CLEARING_ONE_PLUS_T2_POWER),
        8,
        SYMBOLIC_CLEARING_T_POWER,
    )


def scaled_diagonal_polynomial(exponent: int) -> tuple[int, ...]:
    """C(t) q(t)^exponent as an integer polynomial for q=(1+t^2)/(2t)."""
    if exponent >= 0:
        assert exponent <= 3
        return shifted_scaled_polynomial(
            one_plus_t_squared_power(4 + exponent),
            2 ** (3 - exponent),
            3 - exponent,
        )
    magnitude = -exponent
    assert magnitude <= 4
    return shifted_scaled_polynomial(
        one_plus_t_squared_power(4 - magnitude),
        2 ** (3 + magnitude),
        3 + magnitude,
    )


def evaluate_polynomial(coefficients: tuple[int | Fraction, ...], t: Fraction) -> Fraction:
    return sum((Fraction(coefficient) * t**power for power, coefficient in enumerate(coefficients)), Fraction(0))


def symbolic_plus_sector(bonds: list[tuple[int, int]]) -> tuple[list[list[tuple[int | Fraction, ...]]], dict]:
    """Build C(t)R_+(t) over Q[t] without a CAS or floating-point arithmetic."""
    dim = 1 << N_MODES
    bond_sums = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (N_MODES - 1 - index)) & 1) for index in range(N_MODES)]
        bond_sums.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = bond_sums[0] % 2
    exponents = [(value - parity) // 2 for value in bond_sums]
    assert all((value - parity) % 2 == 0 for value in bond_sums)
    assert min(exponents) == -4 and max(exponents) == 3
    diagonal = {exponent: scaled_diagonal_polynomial(exponent) for exponent in set(exponents)}
    full: list[list[tuple[int, ...]]] = [[(0,)] * dim for _ in range(dim)]
    for left in range(dim):
        for right in range(left, dim):
            coefficients = [0] * (SYMBOLIC_ENTRY_DEGREE_BOUND + 1)
            for state, exponent in enumerate(exponents):
                shift = (left ^ state).bit_count() + (right ^ state).bit_count()
                for degree, coefficient in enumerate(diagonal[exponent]):
                    coefficients[degree + shift] += coefficient
            entry = trim_polynomial(coefficients)
            assert polynomial_degree(entry) <= SYMBOLIC_ENTRY_DEGREE_BOUND
            full[left][right] = full[right][left] = entry
    plus: list[list[tuple[int | Fraction, ...]]] = [[(0,)] * SECTOR_DIM for _ in range(SECTOR_DIM)]
    mask = dim - 1
    for left in range(SECTOR_DIM):
        paired_left = mask ^ left
        for right in range(SECTOR_DIM):
            paired_right = mask ^ right
            entries = (
                full[left][right],
                full[left][paired_right],
                full[paired_left][right],
                full[paired_left][paired_right],
            )
            coefficients = [0] * max(len(entry) for entry in entries)
            for entry in entries:
                for degree, coefficient in enumerate(entry):
                    coefficients[degree] += coefficient
            plus[left][right] = trim_polynomial([Fraction(coefficient, 2) for coefficient in coefficients])
    denominators = sorted({coefficient.denominator for row in plus for entry in row for coefficient in entry if isinstance(coefficient, Fraction)})
    actual_degree = max(polynomial_degree(entry) for row in plus for entry in row)
    return plus, {
        "clearing_factor": "8*t^3*(1+t^2)^4",
        "clearing_factor_coefficients_ascending": list(clearing_factor_polynomial()),
        "bond_parity": parity,
        "diagonal_q_exponent_range": [min(exponents), max(exponents)],
        "sector_dimension": SECTOR_DIM,
        "entry_degree_bound": SYMBOLIC_ENTRY_DEGREE_BOUND,
        "actual_max_entry_degree": actual_degree,
        "sector_coefficient_denominators": denominators,
    }


def genericity_data(bonds: list[tuple[int, int]], thresholds: dict[str, int]) -> dict:
    """Exact symbolic family metadata plus independent rational point certificates."""
    symbolic_plus, symbolic = symbolic_plus_sector(bonds)
    specialization = Fraction(1, 3)
    R, _, _ = build_R(6, bonds, specialization)
    rational_plus = split_physical_parity(R)[1]
    clearing_value = evaluate_polynomial(clearing_factor_polynomial(), specialization)
    specialization_matches = all(
        evaluate_polynomial(symbolic_plus[left][right], specialization)
        == clearing_value * rational_plus[left][right]
        for left in range(SECTOR_DIM)
        for right in range(SECTOR_DIM)
    )
    grid = []
    for t in RATIONAL_GRID:
        point_R, q, parity = build_R(6, bonds, t)
        point_plus = split_physical_parity(point_R)[1]
        scale, point_integer = integralize(point_plus)
        point_row = public_pipeline_row(run_pipeline(point_integer, PRIMES[0]))
        grid.append(
            {
                "t": str(t),
                "exp_2K": str(q),
                "bond_parity": parity,
                "integral_scale": str(scale),
                "plus_within_per_prime": {str(PRIMES[0]): point_row},
                "valid_point_certificate": bool(
                    point_row["pair_poly_degree"] == WITHIN_PAIR_DEGREE
                    and point_row["gcd_divides_both"]
                    and point_row["gcd_monic"]
                ),
            }
        )
    specialization_row = next(sample["plus_within_per_prime"][str(PRIMES[0])] for sample in grid if sample["t"] == "1/3")
    point_upper_bound = specialization_row["gcd_degree"]
    assert point_upper_bound == 177
    return {
        "family": {
            **symbolic,
            "A_description": "A(z,t)=charpoly(wedge^2(C(t) R_+(t))) is monic of z-degree 496 in Q[t][z]",
            "B_description": "B(z,t)=partial_z A(z,t)",
            "A_z_degree": WITHIN_PAIR_DEGREE,
            "B_z_degree": WITHIN_PAIR_DEGREE - 1,
            "B_leading_z_coefficient": WITHIN_PAIR_DEGREE,
            "wedge_entry_degree_bound": WEDGE_ENTRY_DEGREE_BOUND,
            "coefficient_t_degree_slope_bound": WEDGE_ENTRY_DEGREE_BOUND,
            "symbolic_specialization_t": "1/3",
            "symbolic_specialization_matches_rational_sector_exactly": specialization_matches,
        },
        "L1_specialization": {
            "point_t": "1/3",
            "modular_gcd_degree": point_upper_bound,
            "conclusion": "deg gcd over Q(t) is at most 177, because monic generic gcd specializes at every t.",
        },
        "L2_genericity": {
            "rho_description": "Res_z(A/G, B/G); H=A/G is monic, J=B/G lies in Q[t][z] with leading z coefficient 496, and H,J are coprime over Q(t).",
            "H_is_monic": True,
            "J_is_in_Q_t_z": True,
            "J_leading_z_coefficient": WITHIN_PAIR_DEGREE,
            "H_J_coprime_over_Q_t": True,
            "generic_gcd_degree_upper_bound": point_upper_bound,
            "generic_distinct_pair_product_lower_bound": WITHIN_PAIR_DEGREE - point_upper_bound,
            "gaussian_within_half_max_exponent_classes": thresholds["even_even_max_exponent_classes"],
            "resultant_degree_bound": GENERIC_RESULTANT_DEGREE_BOUND,
            "denominator_roots_added_to_exceptional_set": ["0", "I", "-I"],
        },
        "exceptional_set_cardinality_bound_including_denominator_roots": GENERIC_RESULTANT_DEGREE_BOUND + 3,
        "rational_grid": grid,
        "rational_grid_all_gcd_degrees_equal_point_bound": all(
            sample["plus_within_per_prime"][str(PRIMES[0])]["gcd_degree"] == point_upper_bound
            for sample in grid
        ),
    }


def public_pipeline_row(row: dict) -> dict:
    """Remove producer-only coefficient vectors while retaining independently checkable digests."""
    return {key: value for key, value in row.items() if not key.startswith("_")}


def sector_trace_data(M_rows: list[list[int]], upto: int, p: int) -> dict:
    """Charpoly/trace recurrence for one integral sector matrix over F_p."""
    dim = len(M_rows)
    assert dim == len(M_rows[0])
    assert p > dim
    assert dim * (p - 1) ** 2 < INT64_MAX
    matrix = np.array([[int(entry) % p for entry in row] for row in M_rows], dtype=np.int64)
    traces = sequential_power_traces(matrix, dim, p)
    inverse = [0] + [pow(k, -1, p) for k in range(1, upto + 1)]
    characteristic = newton_from_power_sums(np.array(traces, dtype=np.int64), dim, p, inverse)
    horner_ok = horner_matrix_zero(characteristic, matrix, p)
    traces = extend_power_sums(traces, characteristic, upto, p)
    spots = []
    for power in (33, 65, 128, 256, 512, 992, 1024):
        if dim < power <= upto:
            direct = trace_power_binary(matrix, power, p)
            spots.append({"k": power, "direct": direct, "recurrence": int(traces[power]), "match": direct == int(traces[power])})
    return {
        "matrix": matrix,
        "traces": traces,
        "characteristic": characteristic,
        "charpoly_degree": dim,
        "charpoly_monic": int(characteristic[0]) == 1,
        "charpoly_sha256": coeff_digest([int(c) for c in characteristic[::-1]]),
        "horner_charpoly_zero_matrix": bool(horner_ok),
        "spot_trace_checks": spots,
        "spot_trace_checks_all_match": all(item["match"] for item in spots),
    }


def cross_pipeline(left_rows: list[list[int]], right_rows: list[list[int]], p: int) -> dict:
    """Monic gcd computation for prod_{i in left,j in right}(z-rho_i sigma_j)."""
    started = time.time()
    left_dim = len(left_rows)
    right_dim = len(right_rows)
    degree = left_dim * right_dim
    assert p > degree and p % 2 == 1
    assert degree * (p - 1) ** 2 < INT64_MAX
    left = sector_trace_data(left_rows, degree, p)
    right = sector_trace_data(right_rows, degree, p)
    pair_powers = np.zeros(degree + 1, dtype=np.int64)
    for power in range(1, degree + 1):
        pair_powers[power] = int(left["traces"][power]) * int(right["traces"][power]) % p
    inverses = [0] + [pow(k, -1, p) for k in range(1, degree + 1)]
    polynomial = newton_from_power_sums(pair_powers, degree, p, inverses)
    derivative = poly_deriv_desc(polynomial, p)
    gcd = poly_gcd(polynomial, derivative, p)
    gcd_degree = 0 if poly_is_zero(gcd) else int(gcd.size) - 1
    divides = poly_is_zero(poly_rem(polynomial, gcd, p)) and poly_is_zero(poly_rem(derivative, gcd, p))
    return {
        "prime": p,
        "left_sector_dimension": left_dim,
        "right_sector_dimension": right_dim,
        "pair_poly_degree": degree,
        "pair_poly_monic": int(polynomial[0]) == 1,
        "pair_poly_sha256": coeff_digest([int(c) for c in polynomial[::-1]]),
        "gcd_degree": gcd_degree,
        "gcd_sha256": coeff_digest([int(c) for c in gcd[::-1]]),
        "gcd_monic": bool(poly_is_zero(gcd) or int(gcd[0]) == 1),
        "gcd_divides_both": bool(divides),
        "left_charpoly": {
            key: value
            for key, value in left.items()
            if key not in {"matrix", "traces", "characteristic"}
        },
        "right_charpoly": {
            key: value
            for key, value in right.items()
            if key not in {"matrix", "traces", "characteristic"}
        },
        "elapsed_seconds": round(time.time() - started, 3),
    }


def run_three_pair_polynomials(sectors: dict[int, list[list[int]]], thresholds: dict[str, int]) -> dict:
    """Run within P=+, within P=-, and P+ x P- at every certificate prime."""
    within: dict[str, dict] = {"+1": {"per_prime": {}}, "-1": {"per_prime": {}}}
    cross = {"per_prime": {}}
    for p in PRIMES:
        assert is_prime(p), f"chosen modulus {p} is not prime"
        plus = public_pipeline_row(run_pipeline(sectors[1], p))
        minus = public_pipeline_row(run_pipeline(sectors[-1], p))
        pm = cross_pipeline(sectors[1], sectors[-1], p)
        within["+1"]["per_prime"][str(p)] = plus
        within["-1"]["per_prime"][str(p)] = minus
        cross["per_prime"][str(p)] = pm
    for sign in ("+1", "-1"):
        values = [within[sign]["per_prime"][str(p)]["gcd_degree"] for p in PRIMES]
        within[sign]["gcd_degree_by_prime"] = values
        within[sign]["required_gcd_degree_floor"] = thresholds["even_even_gcd_degree_floor"]
        within[sign]["meets_floor_at_every_prime"] = all(value >= thresholds["even_even_gcd_degree_floor"] for value in values)
        within[sign]["violates_floor_at_some_prime"] = any(value < thresholds["even_even_gcd_degree_floor"] for value in values)
    cross_values = [cross["per_prime"][str(p)]["gcd_degree"] for p in PRIMES]
    cross["gcd_degree_by_prime"] = cross_values
    cross["required_gcd_degree_floor"] = thresholds["even_odd_gcd_degree_floor"]
    cross["meets_floor_at_every_prime"] = all(value >= thresholds["even_odd_gcd_degree_floor"] for value in cross_values)
    cross["violates_floor_at_some_prime"] = any(value < thresholds["even_odd_gcd_degree_floor"] for value in cross_values)
    return {"within_sector": within, "cross_sector": cross}


def assignment_rows(pair_data: dict, thresholds: dict[str, int]) -> list[dict]:
    """State both P-to-fermionic parity conventions explicitly."""
    rows = []
    for even_sign, odd_sign in ((1, -1), (-1, 1)):
        even_name = "+1" if even_sign == 1 else "-1"
        odd_name = "+1" if odd_sign == 1 else "-1"
        even_within = pair_data["within_sector"][even_name]
        odd_within = pair_data["within_sector"][odd_name]
        cross = pair_data["cross_sector"]
        violations = []
        if even_within["violates_floor_at_some_prime"]:
            violations.append(f"assigned even P={even_name} within-pair gcd degree is below {thresholds['even_even_gcd_degree_floor']}")
        if odd_within["violates_floor_at_some_prime"]:
            violations.append(f"assigned odd P={odd_name} within-pair gcd degree is below {thresholds['odd_odd_gcd_degree_floor']}")
        if cross["violates_floor_at_some_prime"]:
            violations.append(f"P+ x P- cross-pair gcd degree is below {thresholds['even_odd_gcd_degree_floor']}")
        rows.append(
            {
                "physical_to_fermionic_assignment": {
                    "+1": "even" if even_sign == 1 else "odd",
                    "-1": "odd" if even_sign == 1 else "even",
                },
                "assigned_fermionic_even_physical_sign": even_sign,
                "assigned_fermionic_odd_physical_sign": odd_sign,
                "required_gcd_floors": {
                    "even_even": thresholds["even_even_gcd_degree_floor"],
                    "odd_odd": thresholds["odd_odd_gcd_degree_floor"],
                    "even_odd": thresholds["even_odd_gcd_degree_floor"],
                },
                "violations": violations,
                "excluded": bool(violations),
            }
        )
    return rows


def rational_case(name: str, bonds: list[tuple[int, int]], thresholds: dict[str, int]) -> dict:
    """Build and split one exact rational transfer operator, then certify its three polynomials."""
    started = time.time()
    R, q, eps = build_R(6, bonds, Fraction(1, 3))
    flip = physical_flip_map(6)
    exact_flip = (
        sorted(flip) == list(range(64))
        and all(flip[flip[state]] == state for state in range(64))
        and all(flip[state] != state for state in range(64))
    )
    commutes = all(R[flip[i]][flip[j]] == R[i][j] for i in range(64) for j in range(64))
    rational_sectors = split_physical_parity(R)
    scales: dict[str, int] = {}
    integral_sectors: dict[int, list[list[int]]] = {}
    for sign in (1, -1):
        scale, matrix = integralize(rational_sectors[sign])
        scales[str(sign)] = scale
        integral_sectors[sign] = matrix
    pair_data = run_three_pair_polynomials(integral_sectors, thresholds)
    assignments = assignment_rows(pair_data, thresholds)
    return {
        "name": name,
        "kind": "exact_rational_transfer_operator",
        "n_sites": 6,
        "n_bonds": len(bonds),
        "bonds": [list(edge) for edge in bonds],
        "t_tanh_half_Kstar": "1/3",
        "exp_2K": str(q),
        "bond_parity": eps,
        "physical_flip_is_fixed_point_free_involution": exact_flip,
        "R_commutes_with_physical_flip_exactly": commutes,
        "sector_matrix_convention": "(1/2) B_s^T R B_s in b_a^s=e_a+s e_{Pa}; B_s^T B_s=2I",
        "integral_scales": {sign: str(scale) for sign, scale in scales.items()},
        "integral_scale_factorizations": {sign: {str(base): exponent for base, exponent in factorize(scale).items()} for sign, scale in scales.items()},
        **pair_data,
        "assignments": assignments,
        "both_assignments_excluded": all(row["excluded"] for row in assignments),
        "all_thresholds_met_exactly": all(
            entry["gcd_degree_by_prime"] == [entry["required_gcd_degree_floor"]] * len(PRIMES)
            for entry in (*pair_data["within_sector"].values(), pair_data["cross_sector"])
        ),
        "elapsed_seconds": round(time.time() - started, 3),
    }


def synthetic_case(thresholds: dict[str, int]) -> dict:
    """A true six-mode Gaussian, separated into its actual subset-cardinality halves."""
    factors = (2, 3, 5, 7, 11, 13)
    even_values = []
    odd_values = []
    for slot in range(1 << N_MODES):
        value = 1
        for mode, factor in enumerate(factors):
            if (slot >> mode) & 1:
                value *= factor
        (even_values if slot.bit_count() % 2 == 0 else odd_values).append(value)
    sectors = {
        1: [[value if row == column else 0 for column, value in enumerate(even_values)] for row in range(SECTOR_DIM)],
        -1: [[value if row == column else 0 for column, value in enumerate(odd_values)] for row in range(SECTOR_DIM)],
    }
    pair_data = run_three_pair_polynomials(sectors, thresholds)
    return {
        "name": "synthetic six-mode Gaussian split into true even/odd halves",
        "kind": "synthetic_true_gaussian_control",
        "mode_factors": list(factors),
        "true_physical_labeling": {"+1": "even", "-1": "odd"},
        "integral_scales": {"+1": "1", "-1": "1"},
        **pair_data,
        "all_thresholds_met_exactly": all(
            entry["gcd_degree_by_prime"] == [entry["required_gcd_degree_floor"]] * len(PRIMES)
            for entry in (*pair_data["within_sector"].values(), pair_data["cross_sector"])
        ),
    }


def main() -> int:
    signal.signal(signal.SIGALRM, hard_timeout)
    signal.alarm(WALL_SECONDS)
    started = time.time()
    thresholds = threshold_data()
    check("two certificate primes are prime", all(is_prime(p) for p in PRIMES), f"primes={list(PRIMES)}")
    check(
        "all Newton and int64 product bounds hold",
        all(p > CROSS_PAIR_DEGREE and CROSS_PAIR_DEGREE * (p - 1) ** 2 < INT64_MAX for p in PRIMES),
        f"degree={CROSS_PAIR_DEGREE}, max guard={CROSS_PAIR_DEGREE * (max(PRIMES) - 1) ** 2} < {INT64_MAX}",
    )
    check(
        "m=6 parity exponent thresholds are EE=OO=301, EO=364 and gcd floors 195,195,660",
        thresholds["even_even_max_exponent_classes"] == 301
        and thresholds["odd_odd_max_exponent_classes"] == 301
        and thresholds["even_odd_max_exponent_classes"] == 364
        and thresholds["even_even_gcd_degree_floor"] == 195
        and thresholds["odd_odd_gcd_degree_floor"] == 195
        and thresholds["even_odd_gcd_degree_floor"] == 660,
        str(thresholds),
    )

    layer_bond_list = list(layer_bonds((2, 3), (False, False)))
    print("constructing symbolic P=+ family and rational genericity grid ...", flush=True)
    genericity = genericity_data(layer_bond_list, thresholds)
    check(
        "symbolic C(t)R_+(t) has integral coefficients, degree at most 26, and specializes exactly at t=1/3",
        genericity["family"]["symbolic_specialization_matches_rational_sector_exactly"]
        and genericity["family"]["actual_max_entry_degree"] <= SYMBOLIC_ENTRY_DEGREE_BOUND
        and genericity["family"]["sector_coefficient_denominators"] == [1],
        str({
            "actual_max_entry_degree": genericity["family"]["actual_max_entry_degree"],
            "coefficient_denominators": genericity["family"]["sector_coefficient_denominators"],
        }),
    )
    check(
        "L1/L2 genericity input gives >=319 distinct P=+ pair products outside a finite exceptional set",
        genericity["L1_specialization"]["modular_gcd_degree"] == 177
        and genericity["L2_genericity"]["generic_distinct_pair_product_lower_bound"] == 319
        and genericity["L2_genericity"]["generic_distinct_pair_product_lower_bound"] > thresholds["even_even_max_exponent_classes"],
        str(genericity["L2_genericity"]),
    )
    quotient_metadata = genericity["L2_genericity"]
    check(
        "generic resultant metadata has monic H and nonmonic J with leading coefficient 496",
        genericity["family"]["A_z_degree"] == WITHIN_PAIR_DEGREE
        and genericity["family"]["B_leading_z_coefficient"] == WITHIN_PAIR_DEGREE
        and quotient_metadata["H_is_monic"]
        and quotient_metadata["J_is_in_Q_t_z"]
        and quotient_metadata["J_leading_z_coefficient"] == WITHIN_PAIR_DEGREE
        and quotient_metadata["H_J_coprime_over_Q_t"],
        str({
            "A_z_degree": genericity["family"]["A_z_degree"],
            "B_leading_z_coefficient": genericity["family"]["B_leading_z_coefficient"],
            "H_is_monic": quotient_metadata["H_is_monic"],
            "J_leading_z_coefficient": quotient_metadata["J_leading_z_coefficient"],
        }),
    )
    check(
        "all deterministic rational grid samples have the same P=+ modular gcd degree 177",
        genericity["rational_grid_all_gcd_degrees_equal_point_bound"],
        str([(sample["t"], sample["plus_within_per_prime"][str(PRIMES[0])]["gcd_degree"]) for sample in genericity["rational_grid"]]),
    )

    print("certifying true six-mode Gaussian control ...", flush=True)
    synthetic = synthetic_case(thresholds)
    check("synthetic true-parity Gaussian control meets all three thresholds exactly", synthetic["all_thresholds_met_exactly"], str({key: value["gcd_degree_by_prime"] for key, value in synthetic["within_sector"].items()} | {"cross": synthetic["cross_sector"]["gcd_degree_by_prime"]}))

    print("certifying n=6 open-chain sector control ...", flush=True)
    chain = rational_case("n=6 open chain sectors", [(i, i + 1) for i in range(5)], thresholds)
    check("open-chain sectors meet all three thresholds exactly", chain["all_thresholds_met_exactly"], str({key: value["gcd_degree_by_prime"] for key, value in chain["within_sector"].items()} | {"cross": chain["cross_sector"]["gcd_degree_by_prime"]}))

    print("certifying 2x3 physical-parity sectors ...", flush=True)
    layer = rational_case("2x3 open layer sectors", layer_bond_list, thresholds)
    check(
        "2x3 rational operator has exact physical P split",
        layer["physical_flip_is_fixed_point_free_involution"] and layer["R_commutes_with_physical_flip_exactly"],
        f"integral scales={layer['integral_scales']}",
    )
    check(
        "both physical-to-fermionic parity assignments are excluded without ordering",
        layer["both_assignments_excluded"],
        str([{ "assignment": row["physical_to_fermionic_assignment"], "violations": row["violations"] } for row in layer["assignments"]]),
    )

    artifact = {
        "meta": {
            "provenance": "experiments/e119_parity_pair_product.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "arithmetic": "Exact Fraction construction and integral monic models; finite-field Newton identities and Euclidean gcd at listed primes. Decimal arithmetic is absent.",
            "method": "Ordering-free sectorwise and cross-sector pair-product collision obstruction; no eigenvalue ordering or inertia.",
            "resource_budget": {"wall_seconds": WALL_SECONDS, "sympy_or_groebner_calls": 0, "resource_wall_hit": False},
            "scope": "Exact point certificate at t=1/3 plus an algebraic all-but-finitely-many-t physical-P genericity theorem for the open 2x3 layer only.",
        },
        "data": {
            "thresholds": thresholds,
            "primes": list(PRIMES),
            "genericity": genericity,
            "cases": {
                "synthetic_six_mode": synthetic,
                "chain_n6": chain,
                "layer_2x3": layer,
            },
        },
        "checks": CHECKS,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    os.makedirs(ROOT / "results" / "spectral", exist_ok=True)
    output = ROOT / "results" / "spectral" / "parity_pair_product.json"
    output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    signal.alarm(0)
    if FAILURES:
        print(f"FAIL: {FAILURES}")
        return 1
    print("PASS: sectorwise/cross-sector pair-product thresholds reject both physical parity assignments.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
