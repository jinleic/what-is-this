#!/usr/bin/env python3
"""Clean-room standalone verifier for ``e85_kc_interval2``.

This verifier never imports the producer.  It independently rebuilds the
strongest full-transfer Simon--Lieb signs through the stable library API,
checks every displayed finite-torus LP certificate by exact primal--dual
feasibility, and recomputes the small FKG-free SAW controls.
"""
from __future__ import annotations

from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.rigorous_bounds.simon_lieb import (  # noqa: E402
    atanh_rational_interval,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
)


RESULT = ROOT / "results" / "bounds" / "kc_interval2.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
RATIONAL_COSINES = {
    4: (Fraction(1), Fraction(0), Fraction(-1), Fraction(0)),
    6: (
        Fraction(1),
        Fraction(1, 2),
        Fraction(-1, 2),
        Fraction(-1),
        Fraction(-1, 2),
        Fraction(1, 2),
    ),
}
EXTERNAL_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: C36,
}
EXPECTED_ROUTE1_SHAPES = (
    (3, 6, 4),
    (3, 6, 12),
    (4, 5, 4),
    (4, 5, 16),
)
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    payload = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(payload).hexdigest(),
    }


def frac_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def floor_decimal(value: str, places: int = 40) -> str:
    with localcontext() as context:
        context.prec = max(110, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


# ---------------------------------------------------------------------------
# Independent finite-torus LP reconstruction
# ---------------------------------------------------------------------------


def orbit(point: tuple[int, int, int], side: int) -> tuple[tuple[int, int, int], ...]:
    images: set[tuple[int, int, int]] = set()
    for coordinates in set(itertools.permutations(point)):
        for signs in itertools.product((-1, 1), repeat=3):
            images.add(tuple((sign * coordinate) % side for sign, coordinate in zip(signs, coordinates, strict=True)))
    return tuple(sorted(images))


def all_orbits(side: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    unseen = set(itertools.product(range(side), repeat=3))
    result = []
    while unseen:
        current = min(unseen)
        current_orbit = orbit(current, side)
        unseen.difference_update(current_orbit)
        result.append(current_orbit)
    return tuple(result)


def green_values(side: int) -> dict[tuple[int, int, int], Fraction]:
    cosines = RATIONAL_COSINES[side]
    sites = list(itertools.product(range(side), repeat=3))
    modes = [mode for mode in sites if mode != (0, 0, 0)]
    output: dict[tuple[int, int, int], Fraction] = {}
    for site in sites:
        value = Fraction(0)
        for mode in modes:
            lam = 3 - sum(cosines[index] for index in mode)
            phase = sum(site[axis] * mode[axis] for axis in range(3)) % side
            value += cosines[phase] / lam
        output[site] = value / side**3
    return output


def lp_data(side: int, coupling: Fraction) -> tuple[list[list[Fraction]], list[Fraction], list[Fraction], list[str], dict[tuple[int, int, int], Fraction], list[int], list[Fraction]]:
    cosines = RATIONAL_COSINES[side]
    volume = side**3
    nonzero_orbits = tuple(item for item in all_orbits(side) if item != ((0, 0, 0),))
    multiplicities = [len(item) for item in nonzero_orbits]
    dispersions = [3 - sum(cosines[index] for index in item[0]) for item in nonzero_orbits]
    phase = []
    for site_orbit in nonzero_orbits:
        z = site_orbit[0]
        phase.append([
            sum(cosines[sum(mode[axis] * z[axis] for axis in range(3)) % side] for mode in mode_orbit)
            for mode_orbit in nonzero_orbits
        ])

    a: list[list[Fraction]] = []
    b: list[Fraction] = []
    labels: list[str] = []
    for j, lam in enumerate(dispersions):
        row = [Fraction(0)] * len(nonzero_orbits)
        row[j] = Fraction(1)
        a.append(row)
        b.append(Fraction(1, 2) / (coupling * lam))
        labels.append(f"ceiling_mode_{j}")
    for site_orbit, phase_row in zip(nonzero_orbits, phase, strict=True):
        z = site_orbit[0]
        delta = [Fraction(m) - coeff for m, coeff in zip(multiplicities, phase_row, strict=True)]
        a.append(delta)
        b.append(Fraction(volume))
        labels.append("gks_nonnegative_" + "_".join(map(str, z)))
        a.append([-entry for entry in delta])
        b.append(Fraction(0))
        labels.append("spin_upper_" + "_".join(map(str, z)))
    a.append([Fraction(m) for m in multiplicities])
    b.append(Fraction(volume))
    labels.append("zero_mode_nonnegative")
    e_index = next(index for index, item in enumerate(nonzero_orbits) if (1, 0, 0) in item)
    delta_e = [Fraction(m) - coeff for m, coeff in zip(multiplicities, phase[e_index], strict=True)]
    a.append(delta_e)
    b.append(Fraction(volume - 1, 1) / (6 * coupling))
    labels.append("energy_nearest_neighbour")
    return a, b, [Fraction(m) for m in multiplicities], labels, green_values(side), multiplicities, dispersions


def verify_lp_certificate(side: int, certificate: dict[str, object]) -> bool:
    coupling = Fraction(str(certificate["coupling"]))
    a, b, c, labels, _, _, _ = lp_data(side, coupling)
    x = [Fraction(str(item)) for item in certificate["primal_x"]]
    y = [Fraction(str(item)) for item in certificate["dual_y"]]
    if len(x) != len(c) or len(y) != len(b) or labels != certificate["constraint_labels"]:
        return False
    primal_ok = all(value >= 0 for value in x) and all(
        sum(coefficient * value for coefficient, value in zip(row, x, strict=True)) <= bound
        for row, bound in zip(a, b, strict=True)
    )
    dual_ok = all(value >= 0 for value in y) and all(
        sum(y[i] * a[i][j] for i in range(len(a))) >= c[j]
        for j in range(len(c))
    )
    primal_obj = sum(value * coefficient for value, coefficient in zip(x, c, strict=True))
    dual_obj = sum(value * bound for value, bound in zip(y, b, strict=True))
    return (
        primal_ok
        and dual_ok
        and primal_obj == dual_obj
        and frac_text(primal_obj) == certificate["objective_sum_nonzero_modes"]
        and frac_text(primal_obj / side**3) == certificate["S_star"]
        and frac_text(Fraction(1) - primal_obj / side**3) == certificate["magnetization_square_floor"]
    )


# ---------------------------------------------------------------------------
# Independent small SAW enumeration
# ---------------------------------------------------------------------------


def saw_counts(max_steps: int, forbidden: tuple[int, int, int] | None = None) -> list[int]:
    counts = [0] * (max_steps + 1)
    moves = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    seen = {(0, 0, 0)}

    def walk(position: tuple[int, int, int], depth: int) -> None:
        counts[depth] += 1
        if depth == max_steps:
            return
        for dx, dy, dz in moves:
            following = (position[0] + dx, position[1] + dy, position[2] + dz)
            if following != forbidden and following not in seen:
                seen.add(following)
                walk(following, depth + 1)
                seen.remove(following)

    walk((0, 0, 0), 0)
    return counts


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == "experiments/e85_kc_interval2.py"
        and bool(artifact["checks"])
        and all(item["passed"] for item in artifact["checks"]),
    )
    data = artifact["data"]

    # Route 1: rebuild the strongest certificate and its challenge sign using
    # the full stable transfer, never the producer.
    route1 = data["route1_wider_simon_lieb"]
    records = {tuple(item["shape"]): item for item in route1["certificates"]}
    check("route1 fixed family", tuple(records) == EXPECTED_ROUTE1_SHAPES and all(item["cross_section_sites"] > 16 for item in records.values()))
    best = route1["best_exact_safe_certificate"]
    best_shape = tuple(best["shape"])
    check("route1 best belongs to family", best_shape in records and best == records[best_shape])
    for shape in (best_shape, (3, 6, 12)):
        record = records[shape]
        p, q = map(int, record["safe_exact_t"].split("/"))
        safe = exact_criterion_residual(shape, p, q)
        cp, cq = map(int, record["challenge_exact_t"].split("/"))
        challenge = exact_criterion_residual(shape, cp, cq)
        interval = atanh_rational_interval(p, q, dps=90)
        check(
            "route1 exact residual " + "x".join(map(str, shape)),
            safe < 0
            and fingerprint(safe) == record["safe_residual_tQ_minus_P"]
            and floor_decimal(interval[0]) == record["safe_K_decimal_lower"],
            "library full transfer independently reproduces safe sign and directed endpoint",
        )
        check(
            "route1 challenge residual " + "x".join(map(str, shape)),
            (challenge < 0) == record["challenge_passes"]
            and fingerprint(challenge) == record["challenge_residual_tQ_minus_P"],
            "library full transfer independently reproduces the challenge sign",
        )
    bp, bq = map(int, best["safe_exact_t"].split("/"))
    route1_improved = bool(route1["global_interval_improved"])
    check(
        "route1 exact incumbent comparison",
        (bp**36 * C36 > bq**36) == route1_improved,
        "36th-power comparison is independent of decimal output",
    )
    p_poly, q_poly = enumerate_pq_polynomials((2, 2, 2))
    p_value, q_value = exact_pq_at_rational((2, 2, 2), 1, 5)
    check(
        "route1 enumeration control",
        p_value == sum(value * 5 ** (len(p_poly) - 1 - index) for index, value in enumerate(p_poly))
        and q_value == sum(value * 5 ** (len(q_poly) - 1 - index) for index, value in enumerate(q_poly)),
    )

    # Route 2: direct mode sums and exact LP primal--dual certificates.
    route2 = data["route2_finite_volume_floor"]
    for record in route2["records"]:
        side = int(record["side"])
        green = green_values(side)
        check(
            f"route2 exact C_L table L{side}",
            frac_text(green[(0, 0, 0)]) == record["C_L_at_origin"]
            and frac_text(min(green.values())) == record["minimum_C_L"]
            and sum(green.values(), Fraction(0)) == 0,
        )
        for certificate in record["finite_floor_certificates"]:
            check(
                f"route2 LP primal-dual L{side} K{certificate['coupling']}",
                verify_lp_certificate(side, certificate),
                "independently rebuilt Fraction constraints have matching feasible primal and dual",
            )
        coupling = Fraction(6, 25)
        a, b, _, _, green_again, multiplicities, dispersions = lp_data(side, coupling)
        all_ceiling = [Fraction(1, 2) / (coupling * lam) for lam in dispersions]
        threshold = (green_again[(0, 0, 0)] - min(green_again.values())) / 2
        check(
            f"route2 all-ceiling witness L{side}",
            coupling >= threshold
            and all(sum(coefficient * value for coefficient, value in zip(row, all_ceiling, strict=True)) <= bound for row, bound in zip(a, b, strict=True))
            and sum(Fraction(m) * value for m, value in zip(multiplicities, all_ceiling, strict=True)) / side**3
            == green_again[(0, 0, 0)] / (2 * coupling),
        )

    # Route 3: recompute all small inputs and every decisive integer comparison.
    counts = saw_counts(9)
    avoid = saw_counts(8, forbidden=(1, 0, 0))
    check("route3 small SAW rebuild", tuple(counts) == (1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886))
    check("route3 source hash", hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest() == C36_SHA256)
    all_counts = {index: value for index, value in enumerate(counts)} | EXTERNAL_COUNTS
    route3 = data["route3_kesten_attempt"]
    check("route3 stored small sequences", route3["small_exact_counts"] == counts and route3["avoid_fixed_neighbour_counts_a_n"] == avoid)
    for total in range(2, 10):
        for n in range(1, total):
            m = total - n
            check(
                f"route3 injection m{m} n{n}",
                all_counts[total] <= all_counts[m] * (all_counts[n] - avoid[n - 1]),
            )
    for m, n in ((30, 6), (31, 5), (32, 4), (33, 3), (34, 2), (35, 1)):
        check(
            f"route3 mixed injection m{m} n{n}",
            all_counts[m + n] <= all_counts[m] * (all_counts[n] - avoid[n - 1]),
        )
    rows = route3["fekete_certified_finite_ratio_instances"]["rows_n_0_through_7"]
    for row in rows:
        n = int(row["n"])
        passed = C36 * all_counts[n] ** 18 <= all_counts[n + 2] ** 18
        check(f"route3 Fekete n{n}", passed == row["proved_by_c36_fekete"])
    check(
        "route3 conditional n30 improvement", all_counts[32] ** 18 < all_counts[30] ** 18 * C36
    )

    final = data["final_certified_interval"]
    check(
        "honest final conclusion",
        bool(final["improved"]) == route1_improved
        and str(final["lower"]) == (best["safe_K_decimal_lower"] if route1_improved else "0.2122119011661678393310862783954278184914")
        and ("[THEOREM" in data["classification"] if route1_improved else "[UNRESOLVED]" in data["classification"]),
    )
    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {FAILURES}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
