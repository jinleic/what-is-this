"""Exact interlayer expansion and independent anisotropic-series validation.

Run from the repository root with
    .venv/bin/python experiments/e31_interlayer.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np

from ising.interlayer import (
    anisotropic_flm_c2_series,
    anisotropic_box_even_subgraph,
    critical_coupling_2d,
    decoupled_chain_phi,
    exact_2d_overlap_series,
    log_cosh_w_series,
    torus_overlap_sum_float,
    torus_overlap_sum_mp,
)

SCRIPT = "experiments/e31_interlayer.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "interlayer" / "interlayer_expansion.json"
V_ORDER = 6
W_CONTROL_ORDER = 12
DPS = 70


def _fraction_strings(values):
    return [str(Fraction(value)) for value in values]


def _record(checks, name, passed, detail):
    row = {"name": name, "passed": bool(passed), "detail": detail}
    checks.append(row)
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _aitken(values):
    """Three-size geometric-tail extrapolation with an explicitly heuristic error."""

    first, second, third = (float(value) for value in values)
    delta_one = second - first
    delta_two = third - second
    if delta_one == 0:
        return third, 0.0, 0.0
    ratio = delta_two / delta_one
    if abs(ratio) >= 1:
        return third, abs(delta_two), ratio
    estimate = third + delta_two * ratio / (1 - ratio)
    # Twice the extrapolation correction is deliberately conservative, but is
    # an empirical finite-size error estimate, not a rigorous interval.
    uncertainty = 2 * abs(estimate - third) + 32 * np.finfo(float).eps * abs(estimate)
    return estimate, uncertainty, ratio


def main():
    checks = []
    timestamp = datetime.now(timezone.utc).isoformat()

    overlap = exact_2d_overlap_series(V_ORDER)
    correlation_total = tuple(value / 2 for value in overlap)
    correlation_residual = list(correlation_total)
    correlation_residual[0] -= Fraction(1, 2)
    correlation_residual = tuple(correlation_residual)

    flm = anisotropic_flm_c2_series(V_ORDER)
    enlarged = anisotropic_flm_c2_series(V_ORDER, bound_slack=1)
    c1_vanishes = all(
        all(row[1] == 0 for row in anisotropic_box_even_subgraph(shape, V_ORDER, 2))
        for shape in flm.boxes
    )
    _record(
        checks,
        "spin-reversal c1 control",
        c1_vanishes,
        "every exact anisotropic box polynomial has identically zero w^1 coefficient",
    )
    _record(
        checks,
        "2D correlations equal anisotropic 3D FLM coefficient",
        flm.residual == correlation_residual,
        "exact Fraction identity through v^6; nonzero residual orders v^2, v^4, v^6",
    )
    _record(
        checks,
        "anisotropic FLM box-bound stability",
        enlarged.residual == flm.residual,
        f"{len(flm.boxes)} minimal boxes and {len(enlarged.boxes)} enlarged boxes give the same series",
    )
    _record(
        checks,
        "prefactor reconciliation",
        all(flm.total[k] == correlation_total[k] for k in range(V_ORDER + 1)),
        "total [w^2](phi3D-phi2D) differs from residual c2 only by G(0)^2/2=1/2",
    )

    expected_prefactor = tuple(
        Fraction(1, degree) if degree >= 2 and degree % 2 == 0 else Fraction(0)
        for degree in range(W_CONTROL_ORDER + 1)
    )
    prefactor = log_cosh_w_series(W_CONTROL_ORDER)
    _record(
        checks,
        "K=0 chain w-series control",
        prefactor == expected_prefactor and flm.residual[0] == 0,
        "through w^12 the residual vanishes and log cosh(atanh w) has coefficient 1/(2m) at w^(2m)",
    )
    with mp.workdps(DPS):
        kz_control = mp.mpf("0.317")
        chain_phi = decoupled_chain_phi(kz_control, DPS)
        independent_chain_phi = mp.log(mp.e**kz_control + mp.e ** (-kz_control))
        chain_error = abs(chain_phi - independent_chain_phi)
    _record(
        checks,
        "K=0 exact thermodynamic chain control",
        chain_error < mp.mpf("1e-65"),
        f"phi3D(0,Kz)=log(2 cosh Kz); error at Kz=0.317 is {mp.nstr(chain_error, 8)}",
    )

    kc = critical_coupling_2d(90)
    critical_rows = [torus_overlap_sum_float(side, float(kc)) for side in range(3, 11)]
    sides = np.array([row["side"] for row in critical_rows], dtype=float)
    sums = np.array([row["overlap_sum"] for row in critical_rows], dtype=float)
    fit_count = 4
    raw_exponent = float(
        np.polyfit(np.log(sides[-fit_count:]), np.log(sums[-fit_count:]), 1)[0]
    )
    correction_design = np.column_stack(
        [np.ones(fit_count), np.log(sides[-fit_count:]), 1 / sides[-fit_count:]]
    )
    corrected_exponent = float(
        np.linalg.lstsq(correction_design, np.log(sums[-fit_count:]), rcond=None)[0][1]
    )
    exponent_uncertainty = max(
        abs(raw_exponent - corrected_exponent),
        abs(
            raw_exponent
            - float(np.polyfit(np.log(sides[-5:]), np.log(sums[-5:]), 1)[0])
        ),
    )
    _record(
        checks,
        "critical finite-size exponent agrees with eta=1/4 prediction",
        abs(raw_exponent - 1.5) <= exponent_uncertainty,
        f"raw L=7..10 exponent {raw_exponent:.12g}, 1/L-corrected {corrected_exponent:.12g}, systematic {exponent_uncertainty:.3g}; exact prediction 3/2",
    )

    mp_reference = torus_overlap_sum_mp(5, kc, dps=DPS)
    float_reference = next(row for row in critical_rows if row["side"] == 5)
    precision_difference = abs(
        float(mp_reference["overlap_sum"]) - float_reference["overlap_sum"]
    )
    _record(
        checks,
        "binary64 finite-size transfer agrees with multiprecision transfer",
        precision_difference < 1e-12,
        f"L=5 at exact Kc2D; absolute difference {precision_difference:.3g} against {DPS}-dps evaluation",
    )
    maximum_numeric_residual = max(
        max(
            float(row["identity_residual"]),
            float(row["maximum_correlation_asymmetry"]),
        )
        for row in critical_rows
    )
    _record(
        checks,
        "finite-torus correlation identities",
        maximum_numeric_residual < 1e-12,
        f"max of G(0)-1 and inversion-asymmetry residual is {maximum_numeric_residual:.3g}",
    )

    subcritical_rows = []
    for coupling_text in ("0.1", "0.15", "0.2", "0.25", "0.3"):
        coupling = float(coupling_text)
        finite = [torus_overlap_sum_float(side, coupling) for side in (6, 8, 10)]
        total_values = [row["total_c2"] for row in finite]
        estimate, uncertainty, ratio = _aitken(total_values)
        subcritical_rows.append(
            {
                "K": coupling_text,
                "torus_sides": [6, 8, 10],
                "finite_total_c2": [format(value, ".17g") for value in total_values],
                "extrapolated_total_c2": format(estimate, ".17g"),
                "extrapolated_residual_c2": format(estimate - 0.5, ".17g"),
                "empirical_absolute_error": format(uncertainty, ".6g"),
                "successive_correction_ratio": format(ratio, ".10g"),
                "error_scope": "heuristic finite-size estimate, not a rigorous interval",
            }
        )

    exact_data = {
        "conventions": {
            "residual_c2": "[w^2](phi_3D-phi_2D-log(cosh K_z)) = (1/2) sum_{r != 0} G(r)^2",
            "total_c2": "[w^2](phi_3D-phi_2D) = (1/2) sum_r G(r)^2",
            "reconciliation": "total_c2-residual_c2=1/2=[w^2]log(cosh(atanh w))=G(0)^2/2",
        },
        "c1_residual_and_total": "0",
        "v_order": V_ORDER,
        "two_dimensional_overlap_sum": _fraction_strings(overlap),
        "correlation_total_c2": _fraction_strings(correlation_total),
        "correlation_residual_c2": _fraction_strings(correlation_residual),
        "anisotropic_flm_total_c2": _fraction_strings(flm.total),
        "anisotropic_flm_residual_c2": _fraction_strings(flm.residual),
        "nonzero_residual_coefficients": {
            str(degree): str(value)
            for degree, value in enumerate(flm.residual)
            if value
        },
        "minimal_box_count": len(flm.boxes),
        "minimal_boxes": [list(shape) for shape in flm.boxes],
        "enlarged_box_count": len(enlarged.boxes),
        "log_cosh_w_through_12": _fraction_strings(prefactor),
        "independent_methods": {
            "2D": "exact integer polynomial row transfer on a 7x13 cylinder, source six rows from either open boundary; no finite-lattice inversion",
            "3D": "exact anisotropic open-box density-of-states transform plus rectangular finite-lattice Moebius inversion",
        },
    }
    numerical_data = {
        "method": (
            "algebraically exact symmetric row transfer matrix on periodic LxL tori; "
            "binary64 for L=3..10, independently checked at 70 dps for L=5"
        ),
        "two_dimensional_critical_coupling": mp.nstr(kc, 80),
        "subcritical_extrapolations": subcritical_rows,
        "critical_finite_size": [
            {
                "L": row["side"],
                "overlap_sum": format(float(row["overlap_sum"]), ".17g"),
                "total_c2": format(float(row["total_c2"]), ".17g"),
                "overlap_over_L_to_3_over_2": format(
                    float(row["overlap_sum"]) / row["side"] ** 1.5, ".17g"
                ),
            }
            for row in critical_rows
        ],
        "critical_exponent_fit": {
            "fit_window": [7, 8, 9, 10],
            "model": "log S_L = intercept + p log L",
            "raw_p": format(raw_exponent, ".17g"),
            "one_over_L_corrected_p": format(corrected_exponent, ".17g"),
            "systematic_error_from_fit_variation": format(exponent_uncertainty, ".6g"),
            "exact_eta_prediction": "2-2*eta=3/2 for eta=1/4",
        },
        "multiprecision_L5_overlap_sum": mp.nstr(mp_reference["overlap_sum"], 60),
        "binary64_vs_multiprecision_absolute_difference": format(
            precision_difference, ".17g"
        ),
        "infinite_lattice_scope": (
            "For K<Kc the quoted extrapolations have empirical, not certified, finite-size errors. "
            "At K=Kc the infinite sum is divergent; finite L demonstrates L^(3/2) growth and is not an estimate of a finite limit."
        ),
    }

    payload = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": timestamp,
            "precision": {
                "formal_series": "exact Python int/Fraction",
                "multiprecision_transfer_dps": DPS,
                "finite_size_scaling": "IEEE-754 binary64, independently checked against mpmath",
            },
        },
        "data": {
            "exact_series": exact_data,
            "numerics": numerical_data,
            "scope": {
                "expansion_point": "K_z=0",
                "radius_of_convergence": "unknown",
                "isotropic_point": "not reached or solved by this calculation",
                "D_finiteness": "undetermined from three nonzero v coefficients; no ODE claim",
                "crossover": "no rigorous anisotropic-to-isotropic continuation follows",
            },
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT.relative_to(ROOT)}")
    print("PASS e31_interlayer")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e31_interlayer: {error}")
        raise
