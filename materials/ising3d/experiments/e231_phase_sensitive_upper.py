#!/usr/bin/env python3
"""Produce exact certificates for directed Fourier-convolution Gram rows."""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e231_phase_sensitive_upper.py"
OUTPUT = ROOT / "results" / "bounds" / "upper_phase_sensitive.json"
SOURCE_ARTIFACT = ROOT / "results" / "bounds" / "upper_fourpoint.json"
INFRARED_ARTIFACT = ROOT / "results" / "bounds" / "upper_infrared.json"
INPUTS = (
    ROOT / SCRIPT,
    ROOT / "proofs" / "upper_phase_sensitive.md",
    ROOT / "tests" / "test_upper_phase_sensitive.py",
    ROOT / "proofs" / "upper_fourpoint.md",
    SOURCE_ARTIFACT,
    ROOT / "proofs" / "upper_infrared.md",
    INFRARED_ARTIFACT,
)

CLASS_NAME = "PS4-convolution-block-Gram-L1"
CHALLENGE_K = Fraction(6, 25)
FINITE_SIDES = (4, 6)
FINITE_ALPHA = Fraction(25, 12)
RATE_C = Fraction(5168, 525)
BASE_ALL_SIZE_CUTOFF = 96
CPU_BUDGET_SECONDS = 30.0
RSS_LIMIT_BYTES = 2 * 1024**3
COVERAGE = (
    "CONFIGURATION_CONVOLUTION_IDENTITY",
    "DIRECTED_FOURPOINT_BLOCKS",
    "PSD_AND_LOCALIZING_ROWS",
    "POWER_POLYGON_EQUIVALENCE",
    "STRICT_MR4_SEPARATOR_L4",
    "GREEN_LIFT_L4_L6",
    "ALL_SIZE_GREEN_LIFT",
    "ZERO_FLOOR_DIRECTION",
    "ENDPOINT_UNCHANGED",
)
COSINE_TABLES = {
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
Mode = tuple[int, int, int]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _peak_rss_bytes() -> int:
    amount = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return amount if sys.platform == "darwin" else amount * 1024


def _record(
    checks: list[dict[str, Any]], name: str, condition: bool, detail: str
) -> None:
    passed = bool(condition)
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    checks.append({"name": name, "passed": passed, "detail": detail})


def _modes(side: int) -> tuple[Mode, ...]:
    return tuple(itertools.product(range(side), repeat=3))


def _dispersion(mode: Mode, side: int) -> Fraction:
    table = COSINE_TABLES[side]
    return Fraction(3) - sum((table[value] for value in mode), Fraction())


def _green_powers(side: int) -> tuple[dict[Mode, Fraction], Fraction, Fraction]:
    modes = _modes(side)
    volume = side**3
    inverse_sum = sum(
        (Fraction(1, _dispersion(mode, side)) for mode in modes if mode != (0, 0, 0)),
        Fraction(),
    )
    c_zero = inverse_sum / volume
    p_zero = Fraction(1) - FINITE_ALPHA * c_zero
    powers = {(0, 0, 0): p_zero}
    for mode in modes:
        if mode != (0, 0, 0):
            powers[mode] = FINITE_ALPHA / (volume * _dispersion(mode, side))
    return powers, c_zero, p_zero


def _involution_orbits(side: int, momentum: Mode) -> tuple[tuple[Mode, ...], ...]:
    unseen = set(_modes(side))
    answer: list[tuple[Mode, ...]] = []
    while unseen:
        mode = min(unseen)
        partner = tuple((momentum[i] - mode[i]) % side for i in range(3))
        orbit = tuple(sorted({mode, partner}))
        for item in orbit:
            unseen.remove(item)
        answer.append(orbit)
    return tuple(answer)


def _finite_control(
    side: int, source_p_zero: Fraction
) -> tuple[dict[str, Any], dict[str, bool]]:
    modes = _modes(side)
    volume = side**3
    powers, c_zero, p_zero = _green_powers(side)
    nonzero_lambdas = [
        _dispersion(mode, side) for mode in modes if mode != (0, 0, 0)
    ]
    lambda_min = min(nonzero_lambdas)
    lambda_max = max(nonzero_lambdas)

    total_rows = 0
    central_rows = 0
    orbit_counts: list[str] = []
    coverage_ok = True
    for momentum in modes:
        if momentum == (0, 0, 0):
            continue
        orbits = _involution_orbits(side, momentum)
        central = tuple(sorted({(0, 0, 0), momentum}))
        coverage_ok = coverage_ok and central in orbits
        coverage_ok = coverage_ok and sum(len(orbit) for orbit in orbits) == volume
        total_rows += len(orbits)
        central_rows += int(central in orbits)
        orbit_counts.append(
            ",".join(map(str, momentum)) + ":" + str(len(orbits))
        )

    noncentral_rows = total_rows - central_rows
    noncentral_margin = Fraction(volume - 4, 6) - Fraction(2, lambda_min)
    central_rhs_lower = Fraction(volume - 2) * FINITE_ALPHA / (6 * volume)
    central_squared_margin = central_rhs_lower**2 - (
        4 * p_zero * FINITE_ALPHA / (volume * lambda_min)
    )
    p_min_nonzero = FINITE_ALPHA / (volume * lambda_max)
    p_max_nonzero = FINITE_ALPHA / (volume * lambda_min)

    record = {
        "side": side,
        "volume": volume,
        "coupling_K": str(CHALLENGE_K),
        "alpha": str(FINITE_ALPHA),
        "C_L_zero": str(c_zero),
        "p_zero": str(p_zero),
        "lambda_min_nonzero": str(lambda_min),
        "lambda_max": str(lambda_max),
        "p_min_nonzero": str(p_min_nonzero),
        "p_max_nonzero": str(p_max_nonzero),
        "simplex_sum": str(sum(powers.values(), Fraction())),
        "nonzero_momenta_covered": volume - 1,
        "polygon_rows_covered": total_rows,
        "central_rows_covered": central_rows,
        "noncentral_rows_covered": noncentral_rows,
        "orbit_count_sha256": hashlib.sha256(
            "\n".join(orbit_counts).encode("ascii")
        ).hexdigest(),
        "proof_bounds": {
            "noncentral_normalized_margin": str(noncentral_margin),
            "central_rhs_lower": str(central_rhs_lower),
            "central_squared_margin": str(central_squared_margin),
            "central_lhs_squared_upper": str(
                4 * p_zero * FINITE_ALPHA / (volume * lambda_min)
            ),
        },
    }
    conditions = {
        "source_p_zero": p_zero == source_p_zero,
        "probability_simplex": p_zero > 0
        and all(value > 0 for value in powers.values())
        and sum(powers.values(), Fraction()) == 1,
        "rational_cosine_algebra": lambda_min > 0
        and lambda_max == 6
        and all(isinstance(value, Fraction) for value in powers.values()),
        "universal_polygon_coverage": coverage_ok
        and central_rows == volume - 1
        and total_rows == (volume - 1) * volume // 2
        + 4 * ((side // 2) ** 3 - 1),
        "noncentral_polygon_bound": noncentral_margin > 0,
        "central_polygon_squared_bound": central_squared_margin > 0,
    }
    return record, conditions


def _cubic_images(mode: Mode, side: int) -> set[Mode]:
    images: set[Mode] = set()
    for permutation in itertools.permutations(range(3)):
        permuted = tuple(mode[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=3):
            images.add(
                tuple((signs[index] * permuted[index]) % side for index in range(3))
            )
    return images


def _strict_separator() -> tuple[dict[str, Any], dict[str, bool]]:
    side = 4
    volume = side**3
    modes = _modes(side)
    zero = (0, 0, 0)
    special = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
    epsilon = Fraction(1, 100)
    powers = {mode: Fraction() for mode in modes}
    powers[zero] = Fraction(97, 100)
    for mode in special:
        powers[mode] = epsilon

    cap = Fraction(1, 1) / (
        2 * CHALLENGE_K * volume * _dispersion(special[0], side)
    )
    real_space: dict[Mode, Fraction] = {}
    table = COSINE_TABLES[side]
    for displacement in modes:
        value = Fraction()
        for mode, probability in powers.items():
            phase = sum(mode[i] * displacement[i] for i in range(3)) % side
            value += probability * table[phase]
        real_space[displacement] = value

    energy_lower = Fraction(1) - (Fraction(1) - Fraction(1, volume)) / (
        6 * CHALLENGE_K
    )
    edge_values = [real_space[tuple(int(i == axis) for i in range(3))] for axis in range(3)]
    momentum = special[0]
    central_product = powers[zero] * powers[momentum]
    outside_products = [
        powers[mode]
        * powers[tuple((momentum[i] - mode[i]) % side for i in range(3))]
        for mode in modes
        if mode not in {zero, momentum}
    ]

    support = {mode for mode, value in powers.items() if value}
    cubic_support = set().union(*(_cubic_images(mode, side) for mode in special))
    q_row_sums_ok = all(
        sum((powers[first] * powers[second] for second in modes), Fraction())
        == powers[first]
        for first in modes
    )
    q_diag_ok = all(value * value <= value for value in powers.values())

    record = {
        "side": side,
        "volume": volume,
        "coupling_K": str(CHALLENGE_K),
        "epsilon": str(epsilon),
        "power_support": [
            {"mode": list(mode), "p": str(powers[mode])}
            for mode in (zero,) + special
        ],
        "infrared_cap_at_special_modes": str(cap),
        "infrared_cap_margin": str(cap - epsilon),
        "real_space_min": str(min(real_space.values())),
        "real_space_max": str(max(real_space.values())),
        "nearest_neighbour_G": str(edge_values[0]),
        "energy_lower_row": str(energy_lower),
        "energy_margin": str(edge_values[0] - energy_lower),
        "Q_formula": "Q[k,l]=p[k]*p[l]",
        "Q_moment_rank": 1,
        "separating_momentum": list(momentum),
        "central_product_p0_pk": str(central_product),
        "polygon_lhs_squared": str(4 * central_product),
        "polygon_rhs": "0",
        "exact_separation_margin_squared": str(4 * central_product),
    }
    conditions = {
        "base_simplex_and_symmetry": sum(powers.values(), Fraction()) == 1
        and support == {zero, *special}
        and cubic_support == set(special),
        "base_infrared_caps": epsilon < cap
        and all(
            mode == zero
            or probability
            <= Fraction(1, 1)
            / (2 * CHALLENGE_K * volume * _dispersion(mode, side))
            for mode, probability in powers.items()
        ),
        "base_real_space_rows": min(real_space.values()) == Fraction(47, 50)
        and max(real_space.values()) == 1,
        "base_energy_row": edge_values == [Fraction(49, 50)] * 3
        and energy_lower == Fraction(81, 256)
        and edge_values[0] >= energy_lower,
        "base_rank_one_Q_rows": q_row_sums_ok and q_diag_ok,
        "phase_polygon_strict_separation": central_product > 0
        and all(value == 0 for value in outside_products)
        and 4 * central_product == Fraction(97, 2500),
    }
    return record, conditions


def _central_cutoff_margin(side: int) -> Fraction:
    return Fraction((side**3 - 2) ** 2, side**5) - 27 * RATE_C


def _all_size_control(i3_half_upper: Fraction) -> tuple[dict[str, Any], dict[str, bool]]:
    cutoff = BASE_ALL_SIZE_CUTOFF
    if cutoff % 2:
        cutoff += 1
    while _central_cutoff_margin(cutoff) <= 0:
        cutoff += 2

    noncentral_margin = Fraction(cutoff**3 - 4, 6) - Fraction(cutoff**2, 4)
    central_margin = _central_cutoff_margin(cutoff)
    previous_margin = _central_cutoff_margin(cutoff - 2)
    derivative_numerator_at_two = 2**6 + 8 * 2**3 - 20

    record = {
        "profile": "alpha_L=min(I3^-1,D_L^-1)",
        "applicable_couplings": "every 0<K<=I3/2",
        "source_green_rate_cutoff": BASE_ALL_SIZE_CUTOFF,
        "source_zero_mode_rate": f"p0<={RATE_C}/L",
        "rate_constant": str(RATE_C),
        "effective_resistance_identity": "D_L=max_z R_eff(0,z)",
        "effective_resistance_path_bound": "D_L<=3L/2",
        "alpha_lower": "alpha_L>=2/(3L)",
        "dispersion_lower": "lambda_min>=8/L^2",
        "lambda_upper": "lambda(k)<=6",
        "noncentral_sufficient_inequality": "(L^3-4)/6>=L^2/4",
        "central_sufficient_inequality": "L*(1-2/L^3)^2>27*(5168/525)",
        "first_even_polygon_cutoff": cutoff,
        "predecessor_even_side": cutoff - 2,
        "predecessor_central_margin": str(previous_margin),
        "cutoff_central_margin": str(central_margin),
        "cutoff_noncentral_margin": str(noncentral_margin),
        "monotonicity_derivative": "1+8/L^3-20/L^6>0 for L>=2",
        "finite_GKS_floor": "p0=(1/N)sum_z G(z)>=1/N",
        "thermodynamic_floor": "inf_class p0<=5168/(525L)->0",
    }
    conditions = {
        "watson_upper_below_one": 2 * i3_half_upper < 1,
        "cutoff_is_first_even": cutoff == 266
        and previous_margin < 0
        and central_margin > 0,
        "noncentral_all_size_bound": noncentral_margin > 0,
        "central_monotonicity": derivative_numerator_at_two > 0,
        "fixed_eta_obstruction": RATE_C / 266 > 0,
    }
    return record, conditions


def main() -> int:
    started = time.process_time()
    checks: list[dict[str, Any]] = []

    source = json.loads(SOURCE_ARTIFACT.read_text(encoding="utf-8"))
    source_certificate = source["data"]["fourpoint_certificate"]
    source_p_zero = {
        int(row["side"]): Fraction(row["p_zero_equals_M2"])
        for row in source_certificate["finite_lifts"]
    }
    _record(
        checks,
        "source::challenge",
        Fraction(source_certificate["challenge"]["K"]) == CHALLENGE_K,
        "the inherited exact Green finite controls use K=6/25",
    )
    _record(
        checks,
        "source::rate",
        Fraction(source_certificate["all_size_certificate"]["rate_constant"])
        == RATE_C
        and int(
            source_certificate["all_size_certificate"]["challenge_even_cutoff"]
        )
        == BASE_ALL_SIZE_CUTOFF,
        "the exact Green theorem supplies p0<=5168/(525L) from even L>=96",
    )

    finite_controls: list[dict[str, Any]] = []
    for side in FINITE_SIDES:
        control, conditions = _finite_control(side, source_p_zero[side])
        finite_controls.append(control)
        for name, condition in conditions.items():
            _record(
                checks,
                f"finite::L{side}_{name}",
                condition,
                f"exact directed polygon control on all nonzero momenta of T_{side}",
            )

    separator, separator_conditions = _strict_separator()
    for name, condition in separator_conditions.items():
        _record(
            checks,
            f"separator::{name}",
            condition,
            "exact L=4 MR4 point and directed k=(2,0,0) polygon row",
        )

    infrared = json.loads(INFRARED_ARTIFACT.read_text(encoding="utf-8"))
    i3_half_upper = Fraction(infrared["data"]["certified_decimal_upper"])
    all_size, all_size_conditions = _all_size_control(i3_half_upper)
    for name, condition in all_size_conditions.items():
        _record(
            checks,
            f"all_size::{name}",
            condition,
            "exact rational sufficient bounds for every even L>=266",
        )

    _record(
        checks,
        "cross::coverage_unique",
        len(COVERAGE) == len(set(COVERAGE)),
        ",".join(COVERAGE),
    )
    _record(
        checks,
        "cross::resource_bound",
        time.process_time() - started < CPU_BUDGET_SECONDS
        and _peak_rss_bytes() < RSS_LIMIT_BYTES,
        "producer stays below 30 CPU seconds and 2 GiB peak RSS",
    )

    provenance: dict[str, dict[str, Any]] = {}
    for path in INPUTS:
        stat = path.stat()
        provenance[str(path.relative_to(ROOT))] = {
            "sha256": _sha256(path),
            "size_bytes": stat.st_size,
        }

    meta = {
        "producer": SCRIPT,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": provenance,
        "environment": {
            "python": sys.version,
            "executable": sys.executable,
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "hash_seed": os.environ.get("PYTHONHASHSEED", "unset"),
        },
        "arithmetic": (
            "finite controls use only integers and Fraction on L=4,6; the "
            "all-size theorem uses directed exact inequalities and the exact I3 symbol"
        ),
        "benchmark_policy": (
            "K_c=0.221654626 is comparison-only and has no selection, tuning, "
            "fitting, or validation role"
        ),
    }
    data = {
        "headline": (
            "[THEOREM] PS4-convolution-block-Gram-L1 is strictly stronger than "
            "MR4-power-simplex-L1, but the exact Green profile lifts for every "
            "even L>=266 with p0<=5168/(525L); its uniform floor is zero."
        ),
        "classification": {
            "strictness": (
                "[THEOREM] an exact L=4,K=6/25 MR4 point violates a directed "
                "convolution polygon row"
            ),
            "limitation": (
                "[THEOREM] exact zero-uniform-floor limitation for the named "
                "phase-sensitive block relaxation"
            ),
            "finite_controls": (
                "[COMPUTATION] exhaustive directed orbit coverage on rational "
                "cosine tori L=4,6"
            ),
            "endpoint": "[UNRESOLVED] the incumbent certified endpoint remains I3/2",
        },
        "coverage": list(COVERAGE),
        "class": {
            "name": CLASS_NAME,
            "base": "MR4-power-simplex-L1",
            "amplitude": "f_k=N^-1 sum_x sigma_x exp(i k.x)",
            "configuration_identity": "sum_q f_q*f_(k-q)=delta_(k,0)",
            "directed_monomial": "x^(k)_q=f_q*f_(k-q)",
            "phase_block": (
                "H^(k)_[q,r]=E[x^(k)_q*conj(x^(k)_r)]"
            ),
            "rows": [
                "H^(k) is Hermitian positive semidefinite",
                "diag H^(k)_q=Q_[q,k-q]",
                "rows q and k-q agree",
                "H^(0)=Q and H^(0)*1=p",
                "H^(k)*1=0 for k!=0",
                "H^(-k)_[-q,-r]=conj(H^(k)_[q,r])",
                "|H_[q,r]|^2<=Q_[q,k-q]*Q_[r,k-r]",
            ],
            "polygon": (
                "for every k!=0 and involution orbit O under q->k-q, "
                "m_O*sqrt(Q_[q,k-q])<=sum_(O'!=O) "
                "m_O'*sqrt(Q_[r,k-r])"
            ),
            "polygon_equivalence": (
                "the complete directed polygon family is necessary and sufficient "
                "for a PSD block with the prescribed diagonal, equal involution "
                "rows, and H*1=0"
            ),
            "scope_guard": (
                "different total-momentum blocks are not identified by full "
                "cross-channel permutation constraints; no DLR or random-current "
                "row is included"
            ),
        },
        "strict_separator": separator,
        "finite_green_lifts": finite_controls,
        "all_size_green_lift": all_size,
        "direction": {
            "central_rank_one_row": (
                "4*p0*p_k<=(sum_(q notin {0,k}) sqrt(p_q*p_(k-q)))^2"
            ),
            "consequence": (
                "for p_k>0 this is an upper bound on p0, not a positive lower bound"
            ),
            "finite_floor": (
                "the inherited GKS rows alone give p0>=1/N, so an exact finite "
                "p0=0 point is impossible in this class"
            ),
            "uniform_floor": (
                "the lifted Green sequence has p0<=5168/(525L), hence no "
                "L-independent positive zero-mode floor"
            ),
        },
        "outcome": {
            "improved_endpoint": False,
            "incumbent": "I3/2",
            "method_class": CLASS_NAME,
            "green_point_lifts_on_tested_tori": True,
            "green_point_lifts_all_even_sides_from": 266,
            "uniform_positive_floor": False,
            "physicality": (
                "the Green lift is a block pseudo-moment construction, not an "
                "asserted spin or Gibbs measure"
            ),
            "remaining_scope": (
                "[UNRESOLVED] full cross-channel fourth-moment consistency, "
                "higher localizers, DLR identities, and sourced or multi-edge "
                "current identities may still separate the Green sequence"
            ),
        },
    }
    payload = {"meta": meta, "data": data, "checks": checks}
    _record(
        checks,
        "cross::artifact_shape",
        set(payload) == {"meta", "data", "checks"},
        "top-level artifact envelope is exactly meta/data/checks",
    )
    names = [row["name"] for row in checks]
    _record(
        checks,
        "cross::stored_check_names_unique",
        len(names) == len(set(names)),
        f"{len(names)} producer check names are unique before the final row",
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in checks:
        print(f"PASS {row['name']} -- {row['detail']}")
    print(f"wrote {OUTPUT}")
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
