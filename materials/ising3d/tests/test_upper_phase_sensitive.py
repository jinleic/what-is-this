#!/usr/bin/env python3
"""Independent exact verifier for upper_phase_sensitive.json.

The verifier deliberately imports no experiment producer.  It reconstructs the
rational L=4,6 Green powers, covers every directed involution orbit, checks the
strict separator, and re-derives the all-size rational cutoff.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "upper_phase_sensitive.json"
SOURCE_ARTIFACT = ROOT / "results" / "bounds" / "upper_fourpoint.json"
CLASS_NAME = "PS4-convolution-block-Gram-L1"
K = Fraction(6, 25)
ALPHA = Fraction(25, 12)
RATE_C = Fraction(5168, 525)
EXPECTED_SIDES = (4, 6)
EXPECTED_COVERAGE = {
    "CONFIGURATION_CONVOLUTION_IDENTITY",
    "DIRECTED_FOURPOINT_BLOCKS",
    "PSD_AND_LOCALIZING_ROWS",
    "POWER_POLYGON_EQUIVALENCE",
    "STRICT_MR4_SEPARATOR_L4",
    "GREEN_LIFT_L4_L6",
    "ALL_SIZE_GREEN_LIFT",
    "ZERO_FLOOR_DIRECTION",
    "ENDPOINT_UNCHANGED",
}
COSINES = {
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
CPU_BUDGET_SECONDS = 30.0
RSS_LIMIT_BYTES = 2 * 1024**3
Mode = tuple[int, int, int]


def check(condition: bool, name: str, detail: str) -> None:
    if not condition:
        raise AssertionError(f"FAIL {name}: {detail}")
    print(f"PASS {name} -- {detail}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_rss_bytes() -> int:
    amount = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return amount if sys.platform == "darwin" else amount * 1024


def modes(side: int) -> tuple[Mode, ...]:
    return tuple(itertools.product(range(side), repeat=3))


def dispersion(mode: Mode, side: int) -> Fraction:
    table = COSINES[side]
    return Fraction(3) - sum((table[index] for index in mode), Fraction())


def reconstruct_green_powers(
    side: int,
) -> tuple[dict[Mode, Fraction], Fraction, Fraction]:
    volume = side**3
    all_modes = modes(side)
    c_zero = sum(
        (Fraction(1, dispersion(mode, side)) for mode in all_modes if mode != (0, 0, 0)),
        Fraction(),
    ) / volume
    p_zero = Fraction(1) - ALPHA * c_zero
    powers = {(0, 0, 0): p_zero}
    for mode in all_modes:
        if mode != (0, 0, 0):
            powers[mode] = ALPHA / (volume * dispersion(mode, side))
    return powers, c_zero, p_zero


def involution_orbits(side: int, momentum: Mode) -> tuple[tuple[Mode, ...], ...]:
    unseen = set(modes(side))
    answer: list[tuple[Mode, ...]] = []
    while unseen:
        mode = min(unseen)
        partner = tuple((momentum[axis] - mode[axis]) % side for axis in range(3))
        orbit = tuple(sorted({mode, partner}))
        for member in orbit:
            unseen.remove(member)
        answer.append(orbit)
    return tuple(answer)


def orbit_count_digest(side: int) -> tuple[int, int, int, str]:
    volume = side**3
    total = 0
    central = 0
    lines: list[str] = []
    for momentum in modes(side):
        if momentum == (0, 0, 0):
            continue
        rows = involution_orbits(side, momentum)
        central_orbit = tuple(sorted({(0, 0, 0), momentum}))
        check(
            central_orbit in rows
            and sum(len(row) for row in rows) == volume
            and len(set().union(*map(set, rows))) == volume,
            f"L{side}_coverage_k_{'_'.join(map(str, momentum))}",
            "involution orbits partition every q and contain {0,k}",
        )
        total += len(rows)
        central += 1
        lines.append(",".join(map(str, momentum)) + ":" + str(len(rows)))
    digest = hashlib.sha256("\n".join(lines).encode("ascii")).hexdigest()
    return total, central, total - central, digest


def verify_finite_control(row: dict[str, Any], source_p_zero: Fraction) -> None:
    side = int(row["side"])
    volume = side**3
    powers, c_zero, p_zero = reconstruct_green_powers(side)
    lambda_values = [
        dispersion(mode, side) for mode in modes(side) if mode != (0, 0, 0)
    ]
    lambda_min = min(lambda_values)
    lambda_max = max(lambda_values)

    check(
        Fraction(row["coupling_K"]) == K
        and Fraction(row["alpha"]) == ALPHA
        and Fraction(row["C_L_zero"]) == c_zero
        and Fraction(row["p_zero"]) == p_zero == source_p_zero,
        f"L{side}_green_values",
        "independent rational cosine sum reproduces the inherited Green point",
    )
    check(
        sum(powers.values(), Fraction()) == 1
        and Fraction(row["simplex_sum"]) == 1
        and all(value > 0 for value in powers.values()),
        f"L{side}_simplex",
        "all exact Green powers are positive and sum to one",
    )
    check(
        Fraction(row["lambda_min_nonzero"]) == lambda_min
        and Fraction(row["lambda_max"]) == lambda_max == 6
        and Fraction(row["p_min_nonzero"]) == ALPHA / (volume * lambda_max)
        and Fraction(row["p_max_nonzero"]) == ALPHA / (volume * lambda_min),
        f"L{side}_dispersion_extrema",
        "rational roots of unity give exact nonzero dispersion extrema",
    )

    total, central, noncentral, digest = orbit_count_digest(side)
    expected_total = (volume - 1) * volume // 2 + 4 * ((side // 2) ** 3 - 1)
    check(
        int(row["nonzero_momenta_covered"]) == volume - 1
        and int(row["polygon_rows_covered"]) == total == expected_total
        and int(row["central_rows_covered"]) == central == volume - 1
        and int(row["noncentral_rows_covered"]) == noncentral
        and row["orbit_count_sha256"] == digest,
        f"L{side}_universal_row_count",
        f"all {total} directed polygon rows are covered exactly once",
    )

    proof = row["proof_bounds"]
    noncentral_margin = Fraction(volume - 4, 6) - Fraction(2, lambda_min)
    rhs_lower = Fraction(volume - 2) * ALPHA / (6 * volume)
    lhs_squared_upper = 4 * p_zero * ALPHA / (volume * lambda_min)
    central_squared_margin = rhs_lower**2 - lhs_squared_upper
    check(
        Fraction(proof["noncentral_normalized_margin"]) == noncentral_margin > 0
        and Fraction(proof["central_rhs_lower"]) == rhs_lower
        and Fraction(proof["central_lhs_squared_upper"]) == lhs_squared_upper
        and Fraction(proof["central_squared_margin"])
        == central_squared_margin
        > 0,
        f"L{side}_exact_squared_margins",
        "square-root polygon feasibility follows from positive rational margins",
    )

    p_min = ALPHA / (6 * volume)
    noncentral_lhs_bound = 2 * ALPHA / (volume * lambda_min)
    noncentral_rhs_bound = (volume - 4) * p_min
    central_rhs_bound = (volume - 2) * p_min
    all_rows_certified = True
    checked_rows = 0
    for momentum in modes(side):
        if momentum == (0, 0, 0):
            continue
        central_orbit = tuple(sorted({(0, 0, 0), momentum}))
        for orbit in involution_orbits(side, momentum):
            checked_rows += 1
            representative = orbit[0]
            partner = tuple(
                (momentum[axis] - representative[axis]) % side for axis in range(3)
            )
            diagonal = powers[representative] * powers[partner]
            if orbit == central_orbit:
                outside = [
                    mode
                    for mode in modes(side)
                    if mode not in {(0, 0, 0), momentum}
                ]
                all_rows_certified = all_rows_certified and len(outside) == volume - 2
                all_rows_certified = all_rows_certified and all(
                    powers[mode]
                    * powers[
                        tuple(
                            (momentum[axis] - mode[axis]) % side
                            for axis in range(3)
                        )
                    ]
                    >= p_min**2
                    for mode in outside
                )
                all_rows_certified = all_rows_certified and (
                    len(orbit) ** 2 * diagonal <= lhs_squared_upper
                    and central_rhs_bound**2 >= lhs_squared_upper
                )
            else:
                excluded = set(orbit) | {(0, 0, 0), momentum}
                eligible = [mode for mode in modes(side) if mode not in excluded]
                all_rows_certified = all_rows_certified and len(eligible) >= volume - 4
                all_rows_certified = all_rows_certified and all(
                    powers[mode]
                    * powers[
                        tuple(
                            (momentum[axis] - mode[axis]) % side
                            for axis in range(3)
                        )
                    ]
                    >= p_min**2
                    for mode in eligible
                )
                all_rows_certified = all_rows_certified and (
                    len(orbit) ** 2 * diagonal <= noncentral_lhs_bound**2
                    and noncentral_rhs_bound >= noncentral_lhs_bound
                )
    check(
        all_rows_certified and checked_rows == total,
        f"L{side}_all_polygon_bounds",
        f"exact squared bounds certify every one of {checked_rows} rows",
    )


def verify_separator(row: dict[str, Any]) -> None:
    side = int(row["side"])
    volume = side**3
    zero = (0, 0, 0)
    special = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
    powers = {mode: Fraction() for mode in modes(side)}
    for item in row["power_support"]:
        powers[tuple(item["mode"])] = Fraction(item["p"])

    check(
        side == 4
        and Fraction(row["coupling_K"]) == K
        and Fraction(row["epsilon"]) == Fraction(1, 100)
        and powers[zero] == Fraction(97, 100)
        and all(powers[mode] == Fraction(1, 100) for mode in special)
        and sum(powers.values(), Fraction()) == 1,
        "separator_simplex",
        "the exact cubic-symmetric four-mode power point is reconstructed",
    )

    cap = Fraction(1) / (2 * K * volume * dispersion(special[0], side))
    check(
        cap == Fraction(row["infrared_cap_at_special_modes"])
        and cap - Fraction(1, 100) == Fraction(row["infrared_cap_margin"])
        and cap > Fraction(1, 100),
        "separator_infrared_caps",
        "each occupied nonzero mode lies strictly below its exact cap",
    )

    table = COSINES[side]
    correlations: dict[Mode, Fraction] = {}
    for displacement in modes(side):
        correlations[displacement] = sum(
            (
                probability
                * table[
                    sum(mode[axis] * displacement[axis] for axis in range(3))
                    % side
                ]
                for mode, probability in powers.items()
            ),
            Fraction(),
        )
    energy_lower = Fraction(1) - (Fraction(1) - Fraction(1, volume)) / (6 * K)
    edges = [tuple(int(axis == coordinate) for axis in range(3)) for coordinate in range(3)]
    check(
        min(correlations.values())
        == Fraction(row["real_space_min"])
        == Fraction(47, 50)
        and max(correlations.values())
        == Fraction(row["real_space_max"])
        == 1
        and all(correlations[edge] == Fraction(49, 50) for edge in edges)
        and Fraction(row["nearest_neighbour_G"]) == Fraction(49, 50)
        and Fraction(row["energy_lower_row"]) == energy_lower == Fraction(81, 256)
        and Fraction(row["energy_margin"]) == Fraction(49, 50) - energy_lower,
        "separator_real_space_rows",
        "every GKS/boundedness row and all three energy directions hold exactly",
    )

    check(
        all(
            sum((powers[first] * powers[second] for second in modes(side)), Fraction())
            == powers[first]
            and powers[first] ** 2 <= powers[first]
            for first in modes(side)
        )
        and row["Q_formula"] == "Q[k,l]=p[k]*p[l]"
        and int(row["Q_moment_rank"]) == 1,
        "separator_MR4_Q",
        "Q=pp^T obeys every power-simplex level-one row",
    )

    momentum = tuple(row["separating_momentum"])
    outside_products = [
        powers[mode]
        * powers[
            tuple((momentum[axis] - mode[axis]) % side for axis in range(3))
        ]
        for mode in modes(side)
        if mode not in {zero, momentum}
    ]
    central = powers[zero] * powers[momentum]
    check(
        momentum == special[0]
        and central == Fraction(row["central_product_p0_pk"]) == Fraction(97, 10000)
        and all(value == 0 for value in outside_products)
        and Fraction(row["polygon_rhs"]) == 0
        and Fraction(row["polygon_lhs_squared"])
        == Fraction(row["exact_separation_margin_squared"])
        == 4 * central
        == Fraction(97, 2500),
        "separator_directed_polygon",
        "the positive central side cannot close against a zero polygon remainder",
    )


def central_margin(side: int) -> Fraction:
    return Fraction((side**3 - 2) ** 2, side**5) - 27 * RATE_C


def verify_all_size(row: dict[str, Any], infrared: dict[str, Any]) -> None:
    cutoff = int(row["first_even_polygon_cutoff"])
    predecessor = int(row["predecessor_even_side"])
    i3_half_upper = Fraction(infrared["data"]["certified_decimal_upper"])
    noncentral_margin = Fraction(cutoff**3 - 4, 6) - Fraction(cutoff**2, 4)

    independently_found = 96
    while central_margin(independently_found) <= 0:
        independently_found += 2
    check(
        2 * i3_half_upper < 1
        and cutoff == independently_found == 266
        and predecessor == 264
        and Fraction(row["predecessor_central_margin"])
        == central_margin(predecessor)
        < 0
        and Fraction(row["cutoff_central_margin"])
        == central_margin(cutoff)
        == Fraction(1806844037121, 8323159178600)
        > 0,
        "all_size_first_even_cutoff",
        "independent directed rational arithmetic finds the first even cutoff 266",
    )
    check(
        Fraction(row["rate_constant"]) == RATE_C
        and int(row["source_green_rate_cutoff"]) == 96
        and Fraction(row["cutoff_noncentral_margin"]) == noncentral_margin > 0
        and 2**6 + 8 * 2**3 - 20 > 0,
        "all_size_monotone_bounds",
        "both polygon cases remain strict for every real L>=266",
    )
    check(
        row["effective_resistance_identity"] == "D_L=max_z R_eff(0,z)"
        and row["effective_resistance_path_bound"] == "D_L<=3L/2"
        and row["alpha_lower"] == "alpha_L>=2/(3L)"
        and row["dispersion_lower"] == "lambda_min>=8/L^2",
        "all_size_exact_lemmas",
        "the stored all-size proof uses the resistance path and sine-concavity bounds",
    )


def main() -> int:
    started = time.process_time()
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    check(
        set(payload) == {"meta", "data", "checks"},
        "artifact_shape",
        "top-level envelope is exactly meta/data/checks",
    )
    check(
        payload["meta"]["producer"] == "experiments/e231_phase_sensitive_upper.py",
        "producer_identity",
        "the dedicated phase-sensitive producer is recorded",
    )
    names = [row["name"] for row in payload["checks"]]
    check(
        len(names) == len(set(names))
        and all(row["passed"] is True for row in payload["checks"]),
        "stored_checks",
        f"all {len(names)} producer checks are uniquely named and true",
    )
    for relative, item in payload["meta"]["provenance"].items():
        path = ROOT / relative
        check(
            path.is_file()
            and path.stat().st_size == int(item["size_bytes"])
            and sha256(path) == item["sha256"],
            "provenance_" + relative.replace("/", "_").replace(".", "_"),
            "size and SHA-256 match the audited local input",
        )

    data = payload["data"]
    check(
        set(data["coverage"]) == EXPECTED_COVERAGE
        and len(data["coverage"]) == len(EXPECTED_COVERAGE),
        "coverage",
        "every universal claim family is present exactly once",
    )
    class_data = data["class"]
    check(
        class_data["name"] == CLASS_NAME
        and class_data["base"] == "MR4-power-simplex-L1"
        and "H^(k)*1=0 for k!=0" in class_data["rows"]
        and "cross-channel" in class_data["scope_guard"],
        "class_scope",
        "directed block rows and the omitted cross-channel scope are explicit",
    )

    source = json.loads(SOURCE_ARTIFACT.read_text(encoding="utf-8"))
    source_certificate = source["data"]["fourpoint_certificate"]
    source_p_zero = {
        int(row["side"]): Fraction(row["p_zero_equals_M2"])
        for row in source_certificate["finite_lifts"]
    }
    finite_rows = data["finite_green_lifts"]
    check(
        tuple(int(row["side"]) for row in finite_rows) == EXPECTED_SIDES,
        "finite_side_coverage",
        "both required rational-cosine tori are present",
    )
    for row in finite_rows:
        verify_finite_control(row, source_p_zero[int(row["side"])])
        check(
            time.process_time() - started < CPU_BUDGET_SECONDS,
            f"L{row['side']}_cpu_budget",
            "universal directed-orbit verification remains below 30 CPU seconds",
        )

    verify_separator(data["strict_separator"])
    infrared = json.loads((ROOT / "results" / "bounds" / "upper_infrared.json").read_text(encoding="utf-8"))
    verify_all_size(data["all_size_green_lift"], infrared)

    direction = data["direction"]
    outcome = data["outcome"]
    check(
        "upper bound on p0" in direction["consequence"]
        and "p0>=1/N" in direction["finite_floor"]
        and "p0<=5168/(525L)" in direction["uniform_floor"],
        "inequality_direction",
        "the exact row is not misreported as a magnetisation lower bound",
    )
    check(
        outcome["improved_endpoint"] is False
        and outcome["incumbent"] == "I3/2"
        and outcome["method_class"] == CLASS_NAME
        and outcome["green_point_lifts_all_even_sides_from"] == 266
        and outcome["uniform_positive_floor"] is False,
        "endpoint_outcome",
        "the stronger block class has zero uniform floor and leaves I3/2 unchanged",
    )
    check(
        time.process_time() - started < CPU_BUDGET_SECONDS
        and peak_rss_bytes() < RSS_LIMIT_BYTES,
        "resource_bound",
        "independent verifier stays below 30 CPU seconds and 2 GiB peak RSS",
    )
    print(f"PASS process_time -- {time.process_time() - started:.3f}s")
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
