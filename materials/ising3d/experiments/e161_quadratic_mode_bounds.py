"""Exact dimension-to-mode bounds for the all-L local-term no-go.

No floating point is used.  ``minimum_modes`` returns the least integer m with
m(2m-1) at least the supplied proved Lie-algebra dimension lower bound.
"""
from __future__ import annotations

from math import comb, isqrt


def gaussian_lie_dimension(m: int) -> int:
    if m < 0:
        raise ValueError(m)
    return m * (2 * m - 1)


def minimum_modes(dimension_lower_bound: int) -> int:
    """Least integer m satisfying m(2m-1) >= dimension_lower_bound."""

    if dimension_lower_bound < 0:
        raise ValueError(dimension_lower_bound)
    discriminant = 1 + 8 * dimension_lower_bound
    candidate = max(0, (1 + isqrt(discriminant)) // 4)
    while gaussian_lie_dimension(candidate) < dimension_lower_bound:
        candidate += 1
    while candidate and gaussian_lie_dimension(candidate - 1) >= dimension_lower_bound:
        candidate -= 1
    return candidate


def physical_qubits(L: int) -> int:
    if L < 2:
        raise ValueError("the theorem starts at L=2")
    return 2 * L


def physical_gaussian_ceiling(L: int) -> int:
    n = physical_qubits(L)
    return gaussian_lie_dimension(n)


def polynomial_family_bound(L: int) -> int:
    """All path bilinears plus the unused closing rung."""

    return physical_gaussian_ceiling(L) + 1


def ring_family_bound(L: int) -> int:
    """Grades 2 and 4L-2 from the perimeter cycle."""

    return 2 * physical_gaussian_ceiling(L)


def materialised_hypercube_count(L: int) -> int:
    n = physical_qubits(L)
    if L == 2:
        return 0
    return (1 << n) if L % 2 else (1 << (n - 1))


def materialised_family_bound(L: int) -> int:
    return ring_family_bound(L) + materialised_hypercube_count(L)


def compact_central_grade(L: int) -> int:
    n = physical_qubits(L)
    return n if L % 2 else n + 2


def compact_central_orbit_bound(L: int) -> int:
    """Half the middle binomial for odd L; the full non-middle binomial for even L."""

    n = physical_qubits(L)
    value = comb(2 * n, compact_central_grade(L))
    return value // 2 if L % 2 else value


def strongest_dimension_bound(L: int) -> int:
    """Disjoint-grade union, except at L=2 where the central grade is already in the ring."""

    if L == 2:
        return ring_family_bound(L)
    return ring_family_bound(L) + compact_central_orbit_bound(L)


def strict_sqrt_rational_lower(numerator: int, denominator: int) -> int:
    """Least integer q with denominator*q^2 > numerator."""

    if numerator < 0 or denominator <= 0:
        raise ValueError((numerator, denominator))
    q = isqrt(numerator // denominator)
    while denominator * q * q <= numerator:
        q += 1
    while q and denominator * (q - 1) * (q - 1) > numerator:
        q -= 1
    return q


def elementary_exponential_mode_bound(L: int) -> int:
    """Integer form of m > 2^(2L)/sqrt(5(4L+1))."""

    physical_qubits(L)
    return strict_sqrt_rational_lower(1 << (4 * L), 5 * (4 * L + 1))


def size_row(L: int) -> dict:
    strongest = strongest_dimension_bound(L)
    modes = minimum_modes(strongest)
    return {
        "L": L,
        "physical_qubits_n": physical_qubits(L),
        "physical_gaussian_ceiling": physical_gaussian_ceiling(L),
        "polynomial_family_bound": polynomial_family_bound(L),
        "ring_family_bound": ring_family_bound(L),
        "materialised_hypercube_count": materialised_hypercube_count(L),
        "materialised_family_bound": materialised_family_bound(L),
        "compact_central_grade": compact_central_grade(L),
        "compact_central_orbit_bound": compact_central_orbit_bound(L),
        "strongest_dimension_bound": strongest,
        "minimum_modes_from_strongest_bound": modes,
        "elementary_exponential_mode_bound": elementary_exponential_mode_bound(L),
        "predecessor_gaussian_dimension": gaussian_lie_dimension(modes - 1),
        "attained_gaussian_dimension": gaussian_lie_dimension(modes),
    }


def main() -> None:
    rows = [size_row(L) for L in range(2, 129)]
    for row in rows:
        L = row["L"]
        assert row["polynomial_family_bound"] > row["physical_gaussian_ceiling"]
        assert row["predecessor_gaussian_dimension"] < row["strongest_dimension_bound"]
        assert row["attained_gaussian_dimension"] >= row["strongest_dimension_bound"]
        orbit = row["compact_central_orbit_bound"]
        assert 5 * (4 * L + 1) * orbit >= (1 << (4 * L + 1))
        assert row["minimum_modes_from_strongest_bound"] >= row["elementary_exponential_mode_bound"]
    for row in rows[:7]:
        print(
            f"L={row['L']}: dim>={row['strongest_dimension_bound']} "
            f"m>={row['minimum_modes_from_strongest_bound']}"
        )
    print("PASS")


if __name__ == "__main__":
    main()
