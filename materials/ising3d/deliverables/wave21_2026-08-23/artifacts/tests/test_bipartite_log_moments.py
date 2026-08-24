#!/usr/bin/env python3
"""Independent verifier for the e228 artifact, not cross-wave universal coverage.

It owns only the finite log-moment certificate.  It imports nothing from the
producer; the lead's broad/project verifier execution is intentionally left
outside this file's scope.

Clean-room verifier for the bipartite log-moment certificate.

This file imports nothing from experiments/e228_bipartite_log_moments.py.  It
independently

* proves the stated cumulant normalization on exact finite Rademacher sums and
  checks the Hankel minors against explicit sums of squares;
* rebuilds every rational Ising layer matrix directly from its entry formula;
* replays the stored approximate-basis certificate using exact integer matrix
  products and the generalized min-max bounds, including degenerate spectra;
* recomputes centered logarithms and cumulants with a separate 2^-192 outward
  interval implementation; and
* reconstructs the uniform rational-t perturbation bounds for both claimed
  bipartite grids.

Float fields in the JSON are ignored.  Every sign verdict below uses integers
or Fractions only.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "bipartite_log_moments.json"
POINT_T = Fraction(1, 3)
RADIUS = Fraction(1, 1_000_000_000)
VERIFY_BITS = 192
VERIFY_SCALE = 1 << VERIFY_BITS
VERIFY_LOG_TERMS = 64
ALLOWED_TAGS = {
    "[THEOREM]",
    "[LEMMA]",
    "[COMPUTATION]",
    "[CONJECTURE]",
    "[EXTERNAL]",
    "[UNRESOLVED]",
}
FAILS: list[str] = []

Edge = tuple[int, int]
Interval = tuple[int, int]  # endpoints over VERIFY_SCALE


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILS.append(name)


def chain_edges(sites: int) -> tuple[Edge, ...]:
    return tuple((index, index + 1) for index in range(sites - 1))


def rectangle_edges(rows: int, columns: int) -> tuple[Edge, ...]:
    horizontal = tuple(
        (row * columns + column, row * columns + column + 1)
        for row in range(rows)
        for column in range(columns - 1)
    )
    vertical = tuple(
        (row * columns + column, (row + 1) * columns + column)
        for row in range(rows - 1)
        for column in range(columns)
    )
    return tuple(
        sorted(horizontal + vertical, key=lambda edge: (edge[0], edge[1]))
    )


def producer_order_rectangle_edges(rows: int, columns: int) -> tuple[Edge, ...]:
    # Used only to check the artifact's declared ordering.  Matrix construction
    # itself is insensitive to edge order and uses rectangle_edges above.
    edges: list[Edge] = []
    for row in range(rows):
        for column in range(columns):
            site = row * columns + column
            if row + 1 < rows:
                edges.append((site, (row + 1) * columns + column))
            if column + 1 < columns:
                edges.append((site, site + 1))
    return tuple(edges)


def bond_exponents(sites: int, edges: tuple[Edge, ...]) -> tuple[int, tuple[int, ...]]:
    energies: list[int] = []
    for state in range(1 << sites):
        spins = tuple(
            1 - 2 * ((state >> (sites - 1 - vertex)) & 1)
            for vertex in range(sites)
        )
        energies.append(sum(spins[left] * spins[right] for left, right in edges))
    parity = energies[0] % 2
    assert all((energy - parity) % 2 == 0 for energy in energies)
    return parity, tuple((energy - parity) // 2 for energy in energies)


def rebuild_layer(
    sites: int, unordered_edges: tuple[Edge, ...], t: Fraction
) -> tuple[list[list[Fraction]], int, list[list[int]], Fraction, int]:
    """Direct entry-sum construction, without materializing P_t."""

    edges = tuple(sorted(tuple(sorted(edge)) for edge in unordered_edges))
    parity, exponents = bond_exponents(sites, edges)
    dimension = 1 << sites
    q = (1 + t * t) / (2 * t)
    matrix: list[list[Fraction]] = []
    for row in range(dimension):
        output_row: list[Fraction] = []
        for column in range(dimension):
            value = Fraction(0)
            for state, exponent in enumerate(exponents):
                hamming = (row ^ state).bit_count() + (column ^ state).bit_count()
                value += t**hamming * q**exponent
            output_row.append(value)
        matrix.append(output_row)
    assert all(
        matrix[row][column] == matrix[column][row]
        for row in range(dimension)
        for column in range(dimension)
    )
    denominator = 1
    for row in matrix:
        for value in row:
            denominator = math.lcm(denominator, value.denominator)
    integral = [[int(denominator * value) for value in row] for row in matrix]
    return matrix, denominator, integral, q, parity


def digest_integer_matrix(matrix: list[list[int]]) -> str:
    encoded = ";".join(",".join(map(str, row)) for row in matrix).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def replay_basis_certificate(
    row: dict[str, object], integral: list[list[int]], denominator: int
) -> list[tuple[Fraction, Fraction]]:
    certificate = row["spectral_certificate"]
    dimension = len(integral)
    scale = int(certificate["vector_scale"])
    assert int(certificate["diagonal_scale"]) == scale
    vectors = [[int(value) for value in source] for source in certificate["eigenbasis_rows"]]
    diagonal = [int(value) for value in certificate["diagonal_numerators"]]
    assert len(vectors) == dimension
    assert all(len(source) == dimension for source in vectors)
    assert len(diagonal) == dimension and diagonal == sorted(diagonal)

    # Compute G=V^T V and B=V^T M V.  Store M V once: this reads the same exact
    # integers as the producer while using a different accumulation layout.
    matrix_times_vectors = [
        [
            sum(integral[row_index][target] * vectors[target][column] for target in range(dimension))
            for column in range(dimension)
        ]
        for row_index in range(dimension)
    ]
    gram = [[0 for _ in range(dimension)] for _ in range(dimension)]
    congruence = [[0 for _ in range(dimension)] for _ in range(dimension)]
    for left in range(dimension):
        for right in range(dimension):
            total = 0
            weighted_total = 0
            for source in range(dimension):
                total += vectors[source][left] * vectors[source][right]
                weighted_total += vectors[source][left] * matrix_times_vectors[source][right]
            gram[left][right] = total
            congruence[left][right] = weighted_total
    assert all(
        gram[left][right] == gram[right][left]
        and congruence[left][right] == congruence[right][left]
        for left in range(dimension)
        for right in range(dimension)
    )

    scale_squared = scale * scale
    f_numerator = max(
        sum(
            abs(gram[row_index][column] - (scale_squared if row_index == column else 0))
            for column in range(dimension)
        )
        for row_index in range(dimension)
    )
    e_numerator = max(
        sum(
            abs(
                scale * congruence[row_index][column]
                - (
                    diagonal[row_index] * scale_squared
                    if row_index == column
                    else 0
                )
            )
            for column in range(dimension)
        )
        for row_index in range(dimension)
    )
    assert f_numerator == int(certificate["gram_error_numerator"])
    assert e_numerator == int(certificate["congruence_error_numerator"])
    f = Fraction(f_numerator, scale_squared)
    e = Fraction(e_numerator, scale_squared * scale)
    assert str(f) == certificate["gram_error_exact"]
    assert str(e) == certificate["congruence_error_exact"]
    assert f < 1

    intervals = []
    for index, diagonal_numerator in enumerate(diagonal):
        center = Fraction(diagonal_numerator, scale)
        lower = (center - e) / (1 + f) / denominator
        upper = (center + e) / (1 - f) / denominator
        stored = certificate["ordered_eigenvalue_intervals"][index]
        assert lower == Fraction(stored["lower"])
        assert upper == Fraction(stored["upper"])
        assert 0 < lower <= upper
        intervals.append((lower, upper))
    return intervals


# ---------------------------------------------------------------------------
# Independent outward dyadic interval engine (2^-192, 64 atanh terms).
# ---------------------------------------------------------------------------
def floor_div(numerator: int, denominator: int) -> int:
    return numerator // denominator


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def point_interval(value: Fraction) -> Interval:
    scaled = value.numerator * VERIFY_SCALE
    return floor_div(scaled, value.denominator), ceil_div(scaled, value.denominator)


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def iv_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def iv_sub(left: Interval, right: Interval) -> Interval:
    return iv_add(left, iv_neg(right))


def iv_scale(value: Interval, scalar: Fraction | int) -> Interval:
    factor = Fraction(scalar)
    if factor < 0:
        return iv_neg(iv_scale(value, -factor))
    return (
        floor_div(value[0] * factor.numerator, factor.denominator),
        ceil_div(value[1] * factor.numerator, factor.denominator),
    )


def iv_multiply(left: Interval, right: Interval) -> Interval:
    candidates = [a * b for a in left for b in right]
    return floor_div(min(candidates), VERIFY_SCALE), ceil_div(max(candidates), VERIFY_SCALE)


def iv_power(value: Interval, exponent: int) -> Interval:
    assert exponent >= 0
    if exponent == 0:
        return VERIFY_SCALE, VERIFY_SCALE
    denominator = VERIFY_SCALE ** (exponent - 1)
    if exponent % 2:
        return (
            floor_div(value[0] ** exponent, denominator),
            ceil_div(value[1] ** exponent, denominator),
        )
    low = 0 if value[0] <= 0 <= value[1] else min(value[0] ** exponent, value[1] ** exponent)
    high = max(value[0] ** exponent, value[1] ** exponent)
    return floor_div(low, denominator), ceil_div(high, denominator)


def log_near_one(z: Fraction) -> Interval:
    """Enclose 2*atanh(z), 0<=z<=1/3, by a positive-tail series."""

    assert 0 <= z <= Fraction(1, 3)
    z_iv = point_interval(z)
    square = iv_multiply(z_iv, z_iv)
    power = z_iv
    answer = (0, 0)
    for term in range(VERIFY_LOG_TERMS):
        answer = iv_add(answer, iv_scale(power, Fraction(2, 2 * term + 1)))
        power = iv_multiply(power, square)
    tail = iv_scale((0, power[1]), Fraction(9, 4 * (2 * VERIFY_LOG_TERMS + 1)))
    return answer[0], answer[1] + tail[1]


LOG_TWO = log_near_one(Fraction(1, 3))


def rational_log(value: Fraction) -> Interval:
    assert value > 0
    normalized = value
    binary_exponent = 0
    while normalized >= 2:
        normalized /= 2
        binary_exponent += 1
    while normalized < 1:
        normalized *= 2
        binary_exponent -= 1
    z = (normalized - 1) / (normalized + 1)
    return iv_add(log_near_one(z), iv_scale(LOG_TWO, binary_exponent))


def recompute_log_data(
    eigenvalue_intervals: list[tuple[Fraction, Fraction]]
) -> dict[str, dict[int | str, Interval]]:
    logs = []
    for lower, upper in eigenvalue_intervals:
        lower_log = rational_log(lower)
        upper_log = rational_log(upper)
        logs.append((lower_log[0], upper_log[1]))
    dimension = len(logs)
    mean = iv_scale(
        (sum(value[0] for value in logs), sum(value[1] for value in logs)),
        Fraction(1, dimension),
    )
    centered = [(value[0] - mean[1], value[1] - mean[0]) for value in logs]
    moments = {
        order: iv_scale(
            (
                sum(iv_power(value, order)[0] for value in centered),
                sum(iv_power(value, order)[1] for value in centered),
            ),
            Fraction(1, dimension),
        )
        for order in (2, 4, 6, 8)
    }
    k2 = moments[2]
    k4 = iv_sub(moments[4], iv_scale(iv_power(moments[2], 2), 3))
    k6 = iv_add(
        iv_sub(moments[6], iv_scale(iv_multiply(moments[4], moments[2]), 15)),
        iv_scale(iv_power(moments[2], 3), 30),
    )
    k8 = iv_sub(
        iv_add(
            iv_sub(
                iv_sub(moments[8], iv_scale(iv_multiply(moments[6], moments[2]), 28)),
                iv_scale(iv_power(moments[4], 2), 35),
            ),
            iv_scale(iv_multiply(moments[4], iv_power(moments[2], 2)), 420),
        ),
        iv_scale(iv_power(moments[2], 4), 630),
    )
    cumulants = {2: k2, 4: k4, 6: k6, 8: k8}
    p = {
        1: k2,
        2: iv_scale(k4, Fraction(-1, 2)),
        3: iv_scale(k6, Fraction(1, 16)),
        4: iv_scale(k8, Fraction(-1, 272)),
    }
    p1, p2, p3, p4 = (p[index] for index in range(1, 5))
    inequalities = {
        "hankel_n_p2_minus_p1_squared": iv_sub(iv_scale(p2, dimension), iv_power(p1, 2)),
        "hankel_p1_p3_minus_p2_squared": iv_sub(iv_multiply(p1, p3), iv_power(p2, 2)),
        "hankel_p2_p4_minus_p3_squared": iv_sub(iv_multiply(p2, p4), iv_power(p3, 2)),
        "newton_2e2": iv_sub(iv_power(p1, 2), p2),
        "newton_6e3": iv_add(
            iv_sub(iv_power(p1, 3), iv_scale(iv_multiply(p1, p2), 3)),
            iv_scale(p3, 2),
        ),
        "newton_24e4": iv_sub(
            iv_add(
                iv_add(
                    iv_sub(iv_power(p1, 4), iv_scale(iv_multiply(iv_power(p1, 2), p2), 6)),
                    iv_scale(iv_power(p2, 2), 3),
                ),
                iv_scale(iv_multiply(p1, p3), 8),
            ),
            iv_scale(p4, 6),
        ),
    }
    return {"cumulants": cumulants, "powers": p, "inequalities": inequalities}


def exact_sign(value: Interval) -> str:
    if value[0] > 0:
        return "strictly_positive"
    if value[1] < 0:
        return "strictly_negative"
    if value == (0, 0):
        return "zero"
    return "contains_zero"


def rebuild_uniform_eta(
    point_matrix: list[list[Fraction]],
    sites: int,
    edges: tuple[Edge, ...],
    lower_t: Fraction,
    upper_t: Fraction,
) -> Fraction:
    dimension = 1 << sites
    _, exponents = bond_exponents(sites, tuple(sorted(edges)))
    q_lower = (1 + upper_t * upper_t) / (2 * upper_t)
    q_upper = (1 + lower_t * lower_t) / (2 * lower_t)
    rows = [Fraction(0) for _ in range(dimension)]
    for left in range(dimension):
        for right in range(left, dimension):
            lower_entry = Fraction(0)
            upper_entry = Fraction(0)
            for state, exponent in enumerate(exponents):
                hamming = (left ^ state).bit_count() + (right ^ state).bit_count()
                if exponent >= 0:
                    q_lo, q_hi = q_lower**exponent, q_upper**exponent
                else:
                    q_lo, q_hi = q_upper**exponent, q_lower**exponent
                lower_entry += lower_t**hamming * q_lo
                upper_entry += upper_t**hamming * q_hi
            center = point_matrix[left][right]
            assert lower_entry <= center <= upper_entry
            error = max(center - lower_entry, upper_entry - center)
            rows[left] += error
            if right != left:
                rows[right] += error
    return max(rows)


def elementary(values: tuple[Fraction, ...], degree: int) -> Fraction:
    return sum(
        (math.prod(choice, start=Fraction(1)) for choice in itertools.combinations(values, degree)),
        start=Fraction(0),
    )


def exact_rademacher_audit() -> bool:
    for modes in range(1, 6):
        energies = tuple(Fraction(index + 1, index + 3) for index in range(modes))
        samples = tuple(
            sum((sign * energy for sign, energy in zip(signs, energies, strict=True)), Fraction(0))
            for signs in itertools.product((-1, 1), repeat=modes)
        )
        moments = {
            order: sum((sample**order for sample in samples), Fraction(0)) / len(samples)
            for order in (2, 4, 6, 8)
        }
        kappa = {
            2: moments[2],
            4: moments[4] - 3 * moments[2] ** 2,
            6: moments[6] - 15 * moments[4] * moments[2] + 30 * moments[2] ** 3,
            8: (
                moments[8]
                - 28 * moments[6] * moments[2]
                - 35 * moments[4] ** 2
                + 420 * moments[4] * moments[2] ** 2
                - 630 * moments[2] ** 4
            ),
        }
        y = tuple(energy * energy for energy in energies)
        p = {
            power: sum((value**power for value in y), Fraction(0))
            for power in range(1, 5)
        }
        if not (
            kappa[2] == p[1]
            and kappa[4] == -2 * p[2]
            and kappa[6] == 16 * p[3]
            and kappa[8] == -272 * p[4]
        ):
            return False
        h0 = modes * p[2] - p[1] ** 2
        h1 = p[1] * p[3] - p[2] ** 2
        h2 = p[2] * p[4] - p[3] ** 2
        if h0 != sum((a - b) ** 2 for a, b in itertools.combinations(y, 2)):
            return False
        if h1 != sum(a * b * (a - b) ** 2 for a, b in itertools.combinations(y, 2)):
            return False
        if h2 != sum(a * a * b * b * (a - b) ** 2 for a, b in itertools.combinations(y, 2)):
            return False
        if p[1] ** 2 - p[2] != 2 * elementary(y, 2):
            return False
        if p[1] ** 3 - 3 * p[1] * p[2] + 2 * p[3] != 6 * elementary(y, 3):
            return False
        if p[1] ** 4 - 6 * p[1] ** 2 * p[2] + 3 * p[2] ** 2 + 8 * p[1] * p[3] - 6 * p[4] != 24 * elementary(y, 4):
            return False
    return True


def collect_tags(value: object) -> list[str]:
    tags: list[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "tag":
                tags.append(str(item))
            else:
                tags.extend(collect_tags(item))
    elif isinstance(value, list):
        for item in value:
            tags.extend(collect_tags(item))
    return tags


def main() -> int:
    payload = json.loads(ARTIFACT.read_text())
    check("artifact top-level schema", set(payload) == {"meta", "data", "checks"})
    check(
        "artifact checks are unique and true",
        len({item["name"] for item in payload["checks"]}) == len(payload["checks"])
        and all(item["passed"] is True for item in payload["checks"]),
        f"{len(payload['checks'])} stored checks",
    )
    tags = collect_tags(payload["data"])
    check("claim tags are allowed", bool(tags) and set(tags) <= ALLOWED_TAGS, str(sorted(set(tags))))
    check(
        "exact universal Rademacher and sum-of-squares audit",
        exact_rademacher_audit(),
        "n=1..5: cumulant constants, Newton numerators, and three Hankel SOS identities",
    )

    expected = {
        "open_chain_P4": (4, chain_edges(4), "strictly_positive"),
        "open_chain_P6": (6, chain_edges(6), "strictly_positive"),
        "open_grid_2x2": (4, producer_order_rectangle_edges(2, 2), "strictly_negative"),
        "open_grid_2x3": (6, producer_order_rectangle_edges(2, 3), "strictly_negative"),
    }
    cases = {row["label"]: row for row in payload["data"]["physical_cases"]}
    check(
        "finite physical universal-claim coverage is exact",
        set(cases) == set(expected),
        "P4, P6, open 2x2, open 2x3; no implicit 2x4 or all-size row",
    )

    rebuilt: dict[str, tuple[list[list[Fraction]], list[tuple[Fraction, Fraction]]]] = {}
    for label, (sites, declared_edges, expected_h1_sign) in expected.items():
        row = cases[label]
        stored_edges = tuple(tuple(edge) for edge in row["edges"])
        matrix, denominator, integral, q, parity = rebuild_layer(
            sites, tuple(sorted(stored_edges)), POINT_T
        )
        matrix_record = row["exact_matrix"]
        matrix_ok = (
            row["sites"] == sites
            and row["dimension"] == 1 << sites
            and stored_edges == declared_edges
            and denominator == int(matrix_record["common_denominator"])
            and digest_integer_matrix(integral) == matrix_record["integral_matrix_sha256"]
            and q == Fraction(5, 3)
            and parity == row["point"]["bond_energy_parity"]
        )
        check(f"{label} independent exact matrix", matrix_ok, f"D={denominator}")
        intervals = replay_basis_certificate(row, integral, denominator)
        check(
            f"{label} exact generalized-minmax replay",
            len(intervals) == 1 << sites and all(lower > 0 for lower, _ in intervals),
            "full ordered spectrum enclosed without a simplicity assumption",
        )
        independent = recompute_log_data(intervals)
        h1 = independent["inequalities"]["hankel_p1_p3_minus_p2_squared"]
        h2 = independent["inequalities"]["hankel_p2_p4_minus_p3_squared"]
        stored_h1 = row["log_moments"]["inequalities"]["hankel_p1_p3_minus_p2_squared"]
        stored_h2 = row["log_moments"]["inequalities"]["hankel_p2_p4_minus_p3_squared"]
        sign_ok = (
            exact_sign(h1) == expected_h1_sign
            and stored_h1["sign"] == expected_h1_sign
            and (
                int(stored_h1["lower_numerator"]) > 0
                if expected_h1_sign == "strictly_positive"
                else int(stored_h1["upper_numerator"]) < 0
            )
        )
        secondary_ok = (
            exact_sign(h2) == expected_h1_sign
            and stored_h2["sign"] == expected_h1_sign
            and (
                int(stored_h2["lower_numerator"]) > 0
                if expected_h1_sign == "strictly_positive"
                else int(stored_h2["upper_numerator"]) < 0
            )
        )
        check(
            f"{label} independent 2^-192 log-cumulant sign",
            sign_ok and secondary_ok,
            f"p1*p3-p2^2 sign={exact_sign(h1)}, p2*p4-p3^2 sign={exact_sign(h2)}",
        )
        necessities_ok = all(
            independent["powers"][order][0] > 0 for order in range(1, 5)
        ) and all(
            independent["inequalities"][name][0] > 0
            for name in ("newton_2e2", "newton_6e3", "newton_24e4")
        )
        check(
            f"{label} independent power-sum and Newton signs",
            necessities_ok,
            "p_j=sum_i(E_i^2)^j is kept distinct from kappa_(2j)",
        )
        rebuilt[label] = (matrix, intervals)

    uniform_rows = {
        row["label"]: row for row in payload["data"]["uniform_interval_certificates"]
    }
    check(
        "uniform theorem coverage is exact",
        set(uniform_rows) == {"open_grid_2x2", "open_grid_2x3"},
        "both and only the grids named by the finite theorem are present",
    )
    lower_t, upper_t = POINT_T - RADIUS, POINT_T + RADIUS
    for label in ("open_grid_2x2", "open_grid_2x3"):
        sites, edges, _ = expected[label]
        row = uniform_rows[label]
        assert Fraction(row["t_interval"]["lower"]) == lower_t
        assert Fraction(row["t_interval"]["upper"]) == upper_t
        point_matrix, point_intervals = rebuilt[label]
        eta = rebuild_uniform_eta(point_matrix, sites, tuple(sorted(edges)), lower_t, upper_t)
        stored_eta = Fraction(row["entrywise_to_operator_norm_bound"]["exact"])
        enlarged = [(lower - eta, upper + eta) for lower, upper in point_intervals]
        positivity = all(lower > 0 for lower, _ in enlarged)
        independent = recompute_log_data(enlarged)
        h1 = independent["inequalities"]["hankel_p1_p3_minus_p2_squared"]
        stored_h1 = row["log_moments"]["inequalities"]["hankel_p1_p3_minus_p2_squared"]
        check(
            f"{label} independent uniform perturbation",
            eta == stored_eta and positivity,
            f"eta={eta}",
        )
        check(
            f"{label} independent uniform strict sign",
            h1[1] < 0
            and row["strict_obstruction"] is True
            and stored_h1["sign"] == "strictly_negative"
            and int(stored_h1["upper_numerator"]) < 0,
            f"independent 2^-192 upper numerator={h1[1]}",
        )

    theorem = payload["data"]["finite_result"]
    scope = payload["data"]["route_scope"]
    check(
        "claim scope remains finite and benchmark-independent",
        theorem["tag"] == "[THEOREM]"
        and scope["tag"] == "[UNRESOLVED]"
        and payload["data"]["physical_representative"]["benchmark_used_for_selection"] is False
        and payload["data"]["physical_representative"]["K_c_used_for_selection_or_tuning"] is False
        and "No all-coupling" in scope["statement"]
        and "2x4" in scope["statement"],
        "two finite grids on one explicit rational interval; P4/P6 controls; no extrapolation",
    )

    if FAILS:
        print(f"FAIL: {len(FAILS)} checks failed: {', '.join(FAILS)}")
        return 1
    print("PASS: independent bipartite log-moment verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
