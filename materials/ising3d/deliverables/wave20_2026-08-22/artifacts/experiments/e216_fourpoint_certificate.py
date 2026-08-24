#!/usr/bin/env python3
"""Exact MR4 power-simplex moment certificates at a rational challenge.

The class studied here keeps every Fourier-power mode and its full four-point
second-moment matrix.  Its limitation is proved by an exact deterministic
moment lift, not by a floating-point SDP optimum.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import e215_fourpoint_kernel as kernel_component  # noqa: E402
from e93_mag_floor import (  # noqa: E402
    I3_LOW,
    Quad,
    cos_table,
    field_text,
    torus_orbits,
)

SCRIPT = "experiments/e216_fourpoint_certificate.py"
CLASS_NAME = "MR4-power-simplex-L1"
CHALLENGE_K = Fraction(6, 25)
CHALLENGE_ETA = Fraction(1, 32)
RIEMANN_ERROR_A = Fraction(992015, 408608)
FLOOR_RATE_C = Fraction(5168, 525)
CPU_BUDGET_SECONDS = 120.0

Field = Fraction | Quad


def _record(
    checks: list[dict[str, Any]], name: str, condition: bool, detail: str
) -> None:
    passed = bool(condition)
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    checks.append({"name": name, "passed": passed, "detail": detail})


def _ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _least_even_at_least(value: int) -> int:
    return value if value % 2 == 0 else value + 1


def _least_even_strictly_above(value: Fraction) -> int:
    integer = value.numerator // value.denominator + 1
    return _least_even_at_least(integer)


def _fraction_kernel_map(row: dict[str, Any]) -> dict[tuple[int, int, int], Fraction]:
    side = int(row["side"])
    expected_orbits = torus_orbits(side)
    stored = {tuple(item["representative"]): item for item in row["orbits"]}
    green: dict[tuple[int, int, int], Fraction] = {}
    for orbit in expected_orbits:
        representative = orbit[0]
        item = stored[representative]
        if int(item["multiplicity"]) != len(orbit):
            raise AssertionError(f"L={side}: stored orbit multiplicity changed")
        value = Fraction(item["value"])
        for site in orbit:
            green[site] = value
    return green


def _mode_lift(
    kernel_row: dict[str, Any], coupling: Fraction
) -> tuple[dict[str, Any], dict[str, bool]]:
    side = int(kernel_row["side"])
    volume = side**3
    c0 = Fraction(kernel_row["C0"])
    minimum = Fraction(kernel_row["minimum"])
    diameter = Fraction(kernel_row["diameter_D"])
    green = _fraction_kernel_map(kernel_row)

    # alpha is the largest coefficient on the Green ray allowed jointly by
    # the infrared ceiling alpha<=1/(2K) and GKS alpha<=1/D.
    ceiling_scale = Fraction(1, 1) / (2 * coupling)
    gks_scale = Fraction(1, 1) / diameter
    alpha = min(ceiling_scale, gks_scale)
    regime = "infrared-all-ceiling" if ceiling_scale <= gks_scale else "GKS-diameter"
    p0 = Fraction(1) - alpha * c0

    table, _, _ = cos_table(side)
    mode_orbits = torus_orbits(side)
    full_p: dict[tuple[int, int, int], Field] = {(0, 0, 0): p0}
    mode_rows: list[dict[str, Any]] = []
    for orbit in mode_orbits:
        if orbit == ((0, 0, 0),):
            continue
        representative = orbit[0]
        lam: Field = table[0] * 0 + 3
        for coordinate in representative:
            lam = lam - table[coordinate]
        probability = alpha / (volume * lam)
        for mode in orbit:
            full_p[mode] = probability
        mode_rows.append(
            {
                "representative": list(representative),
                "multiplicity": len(orbit),
                "lambda": field_text(lam),
                "p": field_text(probability),
            }
        )

    sum_p: Field = sum(full_p.values(), table[0] * 0)
    nonnegative = all(value >= 0 for value in full_p.values())
    at_most_one = all(value <= 1 for value in full_p.values())
    cap_ok = True
    conjugacy_ok = True
    for mode, probability in full_p.items():
        if mode == (0, 0, 0):
            continue
        lam: Field = table[0] * 0 + 3
        for coordinate in mode:
            lam = lam - table[coordinate]
        cap_ok = cap_ok and probability <= Fraction(1, 1) / (
            2 * coupling * volume * lam
        )
        inverse = tuple((-coordinate) % side for coordinate in mode)
        conjugacy_ok = conjugacy_ok and probability == full_p[inverse]

    # The direct inverse transform checks every stored real-space orbit.  It
    # does not rely on the formula p0+alpha*C after forming the mode vector.
    inverse_transform_ok = True
    gks_ok = True
    spectral_min: Field | None = None
    spectral_max: Field | None = None
    for item in kernel_row["orbits"]:
        site = tuple(item["representative"])
        direct: Field = table[0] * 0 + p0
        for mode, probability in full_p.items():
            if mode == (0, 0, 0):
                continue
            phase = sum(
                mode[axis] * site[axis] for axis in range(3)
            ) % side
            direct = direct + probability * table[phase]
        expected: Field = p0 + alpha * green[site]
        inverse_transform_ok = inverse_transform_ok and direct == expected
        gks_ok = gks_ok and 0 <= direct <= 1
        if spectral_min is None or direct < spectral_min:
            spectral_min = direct
        if spectral_max is None or direct > spectral_max:
            spectral_max = direct

    # q_{kl}=p_k p_l is represented by its exact factor p.  These checks are
    # the complete simplex-localizing rows: Q1=p, q>=0, and diag(Q)<=p.
    q_row_sums_ok = sum_p == 1 and all(
        probability * sum_p == probability for probability in full_p.values()
    )
    q_nonnegative = nonnegative
    q_diagonal_ok = all(
        probability - probability * probability >= 0
        for probability in full_p.values()
    )

    edge = (1 % side, 0, 0)
    energy_lower = Fraction(1) - (Fraction(1) - Fraction(1, volume)) / (
        6 * coupling
    )
    energy_ok = p0 + alpha * green[edge] >= energy_lower

    checks = {
        "probability_simplex": sum_p == 1 and nonnegative and at_most_one,
        "mode_conjugacy": conjugacy_ok,
        "infrared_caps": cap_ok,
        "inverse_transform": inverse_transform_ok,
        "gks_range": gks_ok and spectral_min is not None and spectral_max == 1,
        "energy_row": energy_ok,
        "fourpoint_row_sums": q_row_sums_ok,
        "fourpoint_nonnegative": q_nonnegative,
        "fourpoint_diagonal_localizers": q_diagonal_ok,
        "moment_psd_factorization": q_row_sums_ok and len(full_p) == volume,
    }

    record = {
        "side": side,
        "volume": volume,
        "coupling_K": str(coupling),
        "regime": regime,
        "alpha": str(alpha),
        "p_zero_equals_M2": str(p0),
        "nonzero_mass": str(Fraction(1) - p0),
        "C0": str(c0),
        "minimum_C": str(minimum),
        "diameter_D": str(diameter),
        "mode_orbits": mode_rows,
        "q_representation": {
            "formula": "q[k,l]=p[k]*p[l]",
            "moment_matrix_factor": "[1;p] [1;p]^T",
            "covariance": "Q-p*p^T=0",
            "rank": 1,
        },
        "real_space_G_min": field_text(spectral_min),
        "real_space_G_max": field_text(spectral_max),
        "energy_lower_row": str(energy_lower),
    }
    return record, checks


def build_certificate_data(
    kernel_data: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    started = time.process_time()
    if kernel_data is None:
        kernel_data, _ = kernel_component.build_kernel_data()

    checks: list[dict[str, Any]] = []
    lifts: list[dict[str, Any]] = []
    for kernel_row in kernel_data["sides"]:
        lift, lift_checks = _mode_lift(kernel_row, CHALLENGE_K)
        lifts.append(lift)
        for name, condition in lift_checks.items():
            _record(
                checks,
                f"L{lift['side']}_{name}",
                condition,
                f"exact {CLASS_NAME} row at K={CHALLENGE_K}",
            )
        if time.process_time() - started > CPU_BUDGET_SECONDS:
            raise RuntimeError(
                f"certificate CPU budget exceeded after L={lift['side']}"
            )

    twice_k = 2 * CHALLENGE_K
    endpoint_gap_lower = I3_LOW - twice_k
    cutoff_raw = _ceil_fraction(RIEMANN_ERROR_A / endpoint_gap_lower)
    challenge_cutoff = _least_even_at_least(cutoff_raw)
    eta_cutoff = max(
        challenge_cutoff,
        _least_even_strictly_above(FLOOR_RATE_C / CHALLENGE_ETA),
    )

    _record(
        checks,
        "challenge_strictly_below_Watson",
        endpoint_gap_lower > 0,
        f"I3-2K >= {endpoint_gap_lower} > 0",
    )
    _record(
        checks,
        "all_size_challenge_cutoff",
        I3_LOW - RIEMANN_ERROR_A / challenge_cutoff >= twice_k,
        f"C_L(0)>=I3-A/L>=2K for every even L>={challenge_cutoff}",
    )
    _record(
        checks,
        "all_size_eta_obstruction",
        eta_cutoff >= challenge_cutoff
        and FLOOR_RATE_C / eta_cutoff < CHALLENGE_ETA,
        f"p0<=({FLOOR_RATE_C})/L<{CHALLENGE_ETA} for every even L>={eta_cutoff}",
    )

    side12 = next(row for row in lifts if row["side"] == 12)
    side12_p0 = Fraction(side12["p_zero_equals_M2"])
    side12_gap = CHALLENGE_ETA - side12_p0
    _record(
        checks,
        "L12_fixed_eta_primal_obstruction",
        side12_gap > 0,
        f"eta-p0={side12_gap} > 0",
    )
    _record(
        checks,
        "finite_rank_one_lifts_complete",
        all(
            row["q_representation"]["rank"] == 1
            and Fraction(row["p_zero_equals_M2"]) >= 0
            for row in lifts
        ),
        f"{len(lifts)} exact mode-resolved outer-product lifts",
    )

    data = {
        "class": {
            "name": CLASS_NAME,
            "physical_intensity": "P_k=|sigma_hat(k)|^2/N^2",
            "two_point_variables": "p_k=E[P_k]=G_hat(k)/N",
            "four_point_variables": "q_kl=E[P_k P_l]",
            "constraints": [
                "p>=0, sum_k p_k=1, p_k=p_-k",
                "p_k<=1/(2*K*N*lambda(k)) for k!=0",
                "0<=G_p(z)<=1 and the inherited nearest-neighbour energy row",
                "q>=0, Q*1=p, diag(Q)<=p",
                "[[1,p^T],[p,Q]] is positive semidefinite",
            ],
            "exact_projection_theorem": (
                "for every feasible two-point p, q=p*p^T satisfies every listed "
                "four-point row; hence the projection and the optimum p0 are "
                "exactly those of the inherited two-point relaxation"
            ),
            "psd_identity": (
                "[[1,p^T],[p,p*p^T]]=[1;p][1;p]^T and Q-p*p^T=0"
            ),
        },
        "challenge": {
            "K": str(CHALLENGE_K),
            "twice_K": str(twice_k),
            "eta": str(CHALLENGE_ETA),
            "benchmark_role": "none",
            "I3_lower_input": str(I3_LOW),
            "I3_minus_twice_K_lower": str(endpoint_gap_lower),
            "strictly_below_Watson_endpoint": True,
        },
        "finite_lifts": lifts,
        "fixed_eta_L12": {
            "eta": str(CHALLENGE_ETA),
            "p0": str(side12_p0),
            "exact_gap_eta_minus_p0": str(side12_gap),
        },
        "all_size_certificate": {
            "C_L0_lower": "C_L(0)>=I3-A/L",
            "A": str(RIEMANN_ERROR_A),
            "negative_kernel_bound": "-min_z C_L(z)<=sqrt(Sigma4)/(4L)",
            "floor_rate": f"p0<=({FLOOR_RATE_C})/L",
            "rate_constant": str(FLOOR_RATE_C),
            "challenge_even_cutoff": challenge_cutoff,
            "fixed_eta_even_cutoff": eta_cutoff,
            "witness_scale": "alpha_L=min(1/(2K),1/D_L)",
            "witness_zero_mode": "p0_L=1-alpha_L*C_L(0)",
            "thermodynamic_conclusion": (
                "infimum p0 of MR4-power-simplex-L1 tends to zero; no uniform "
                "positive magnetisation floor is certifiable in this class"
            ),
        },
        "scope": {
            "negative_theorem": (
                "exact limitation of MR4-power-simplex-L1, including all modes "
                "and the level-one power-intensity moment matrix"
            ),
            "not_covered": [
                "phase-sensitive Fourier convolution identities from sigma_x^2=1",
                "higher localizing levels or full spin moment matrices",
                "DLR conditional identities",
                "sourced or multi-current switching inequalities",
                "the true physical value of M_L^2",
            ],
        },
        "process_time_seconds": f"{time.process_time() - started:.6f}",
    }
    return data, checks


def main() -> int:
    data, checks = build_certificate_data()
    for check in checks:
        print(f"PASS {check['name']} -- {check['detail']}")
    print(
        json.dumps(
            {
                "script": SCRIPT,
                "class": data["class"]["name"],
                "K": data["challenge"]["K"],
                "check_count": len(checks),
            },
            sort_keys=True,
        )
    )
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
