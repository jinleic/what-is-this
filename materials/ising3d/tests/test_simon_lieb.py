"""Standalone acceptance checks for the Simon--Lieb finite-box bound."""

from __future__ import annotations

import math

import mpmath as mp

from ising.rigorous_bounds.simon_lieb import (
    boundary_multiplicities,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
    kappa_float_from_v,
    kappa_mpf_from_v,
    mpmath_relative_roundoff_bound,
    solve_box_root_float,
)


def _scaled_polynomial(coefficients: tuple[int, ...], p: int, q: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * p**power * q ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def main() -> None:
    # B={0}: P=1, Q=6, hence kappa(v)=6v and the root is exactly v=1/6.
    singleton_p, singleton_q = exact_pq_at_rational((1, 1, 1), 1, 6)
    assert (singleton_p, singleton_q) == (1, 6)
    assert exact_criterion_residual((1, 1, 1), 1, 6) == 0
    assert exact_criterion_residual((1, 1, 1), 1, 7) < 0
    assert exact_criterion_residual((1, 1, 1), 1, 5) > 0
    assert kappa_float_from_v((1, 1, 1), 1.0 / 6.0) == 1.0
    singleton_v, singleton_k = solve_box_root_float((1, 1, 1))
    assert abs(singleton_v - 1.0 / 6.0) < 2e-16
    assert abs(singleton_k - math.atanh(1.0 / 6.0)) < 2e-16
    print("singleton criterion 6*tanh(K)<1 reproduced exactly: PASS")

    # The side-three cube is enumerated over all 2^26 configurations with
    # sigma_origin fixed; global spin flip supplies the other half exactly.
    p_poly, q_poly = enumerate_pq_polynomials((3, 3, 3), chunk_size=1 << 21)
    assert len(p_poly) == len(q_poly) == 55  # 54 internal bonds.
    assert p_poly[0] == 1 and all(coefficient >= 0 for coefficient in p_poly)
    assert q_poly[0] == 0 and all(coefficient >= 0 for coefficient in q_poly)
    cycle_count = 1 << (54 - 27 + 1)
    assert sum(p_poly) == cycle_count
    assert sum(q_poly) == 54 * cycle_count

    p, q = 1, 5
    transfer_p, transfer_q = exact_pq_at_rational((3, 3, 3), p, q)
    assert transfer_p == _scaled_polynomial(p_poly, p, q)
    assert transfer_q == _scaled_polynomial(q_poly, p, q)
    assert p * transfer_q - q * transfer_p == exact_criterion_residual(
        (3, 3, 3), p, q
    )
    for shape in ((1, 1, 1), (2, 2, 2), (2, 3, 2)):
        rational_p, rational_q = exact_pq_at_rational(shape, 2, 9)
        assert exact_criterion_residual(shape, 2, 9) == (
            2 * rational_q - 9 * rational_p
        )
    print("3x3x3 exact enumeration agrees with independent parity transfer: PASS")

    shapes = (
        (1, 1, 1),
        (2, 2, 2),
        (3, 3, 3),
        (4, 4, 4),
        (4, 4, 8),
        (4, 4, 16),
        (4, 4, 32),
    )
    expected_k = (
        0.16823611831060642,
        0.1834637215521519,
        0.19740388284055016,
        0.20193444723608656,
        0.20465062730735273,
        0.20515984606810345,
        0.20517823158101176,
    )
    roots = [solve_box_root_float(shape) for shape in shapes]
    k_bounds = [root[1] for root in roots]
    assert all(
        right >= left for left, right in zip(k_bounds, k_bounds[1:])
    )
    for shape, actual, expected in zip(shapes, k_bounds, expected_k, strict=True):
        assert abs(actual - expected) < 2e-13, (shape, actual, expected)
        assert sum(boundary_multiplicities(shape)) == 2 * (
            shape[0] * shape[1]
            + shape[0] * shape[2]
            + shape[1] * shape[2]
        )
    assert k_bounds[-1] < mp.mpf("0.221654626")
    print("finite-box lower-bound sequence is non-decreasing: PASS")
    print("largest test bound remains below 0.221654626: PASS")

    with mp.workdps(70):
        value = kappa_mpf_from_v((1, 1, 1), mp.mpf(1) / 7, dps=60)
        assert abs(value - mp.mpf(6) / 7) < mp.mpf("1e-55")
        error_bound = mpmath_relative_roundoff_bound((4, 4, 32), 50)
        assert mp.mpf(0) < error_bound < mp.mpf("1e-35")
    print("multiprecision evaluator and conservative roundoff bound: PASS")
    print("ALL Simon--Lieb checks: PASS")


if __name__ == "__main__":
    main()
