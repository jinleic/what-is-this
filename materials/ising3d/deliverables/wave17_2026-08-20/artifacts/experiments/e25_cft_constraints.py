"""Universal-data and singularity constraints for candidate 3-D Ising free energies.

Run from the repository root with

    .venv/bin/python experiments/e25_cft_constraints.py

The conformal-bootstrap inputs are quoted from the locally archived papers.  Derived
critical exponents use rectangular propagation of the quoted marginal bounds (no
unpublished covariance is assumed).  Candidate singularities are found from the callable
alone by an adaptive positive-real-axis curvature scan; claimed critical points are used
only after the scan as independent comparisons.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from statistics import median
from typing import Callable

import mpmath as mp

from ising.onsager import onsager_free_energy

PRECISION = 50
mp.mp.dps = PRECISION
SCRIPT = "experiments/e25_cft_constraints.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "cft" / "cft_constraints.json"
CLAIMS_PATH = ROOT / "results" / "reference_data" / "claimed_solutions.json"

# A deliberately conservative operational constraint requested for the candidate checker.
# The higher-precision bootstrap-derived value is recorded separately below.
ALPHA_CONSTRAINT = mp.mpf("0.110")
ALPHA_CONSTRAINT_UNCERTAINTY = mp.mpf("0.001")

# The scan starts above the removable high-temperature endpoint K=0 and covers all
# positive-real singularities claimed in the repository plus the exact 2-D control.
SCAN_LOW = mp.mpf("0.03")
SCAN_HIGH = mp.mpf("0.58")
INITIAL_INTERVALS = 28
REFINEMENT_INTERVALS = 12
REFINEMENT_LEVELS = 7
OFFSET_LEVELS = 8


CITATIONS: dict[str, dict[str, str]] = {
    "kos2016": {
        "authors": "F. Kos, D. Poland, D. Simmons-Duffin, and A. Vichi",
        "title": "Precision Islands in the Ising and O(N) Models",
        "venue": "JHEP 08 (2016) 036",
        "doi": "10.1007/JHEP08(2016)036",
        "arxiv": "1603.04436",
        "local_path": "sources/fulltext/num_kos2016.pdf",
        "location": "abstract, introduction, and Eqs. (3.1)-(3.2)",
    },
    "campostrini2002": {
        "authors": "M. Campostrini, A. Pelissetto, P. Rossi, and E. Vicari",
        "title": "25th-Order High-Temperature Expansion Results for Three-Dimensional Ising-Like Systems on the Simple-Cubic Lattice",
        "venue": "Phys. Rev. E 65 (2002) 066127",
        "doi": "10.1103/PhysRevE.65.066127",
        "arxiv": "cond-mat/0201180",
        "local_path": "sources/fulltext/num_campostrini2002.pdf",
        "location": "abstract; Table 1; Sec. 5 and Tables 7-8",
    },
    "ferrenberg2018": {
        "authors": "A. M. Ferrenberg, J. Xu, and D. P. Landau",
        "title": "Pushing the Limits of Monte Carlo Simulations for the Three-Dimensional Ising Model",
        "venue": "Phys. Rev. E 97 (2018) 043301",
        "doi": "10.1103/PhysRevE.97.043301",
        "arxiv": "1806.03558",
        "local_path": "sources/fulltext/num_ferrenberg2018.pdf",
        "location": "abstract; Eqs. (39), (46), and (47)",
    },
    "onsager1944": {
        "authors": "L. Onsager",
        "title": "Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition",
        "venue": "Phys. Rev. 65 (1944) 117-149",
        "doi": "10.1103/PhysRev.65.117",
        "arxiv": "",
        "local_path": "metadata in sources/manifest.yaml; independently reproduced in src/ising/onsager/",
        "location": "thermodynamic-limit square-lattice free energy",
    },
    "zhang2007": {
        "authors": "Z.-D. Zhang",
        "title": "Conjectures on Exact Solution of Three-Dimensional Simple Orthorhombic Ising Lattices",
        "venue": "Philosophical Magazine 87 (2007) 5309-5419",
        "doi": "10.1080/14786430701646325",
        "arxiv": "0705.1045",
        "local_path": "sources/fulltext/cla_zhang2007.pdf",
        "location": "Eqs. (3.37), (3.60), and (3.62)",
    },
    "dzhang2021": {
        "authors": "D. Zhang",
        "title": "Exact Solution for Three-Dimensional Ising Model",
        "venue": "Symmetry 13 (2021) 1837",
        "doi": "10.3390/sym13101837",
        "arxiv": "2110.11233",
        "local_path": "sources/fulltext/cla_zhang2021.pdf",
        "location": "Eqs. (32), (35)-(36), and (44)-(45)",
    },
}


def _mpstr(value: mp.mpf | mp.mpc, digits: int = 32) -> str:
    """Serialize an mpmath real without conversion through binary float."""

    real = mp.mpf(mp.re(value))
    if real == 0:
        return "0"
    return mp.nstr(real, digits)


def _rectangular_bound(
    function: Callable[[mp.mpf, mp.mpf], mp.mpf],
    delta_sigma: mp.mpf,
    delta_epsilon: mp.mpf,
    sigma_uncertainty: mp.mpf,
    epsilon_uncertainty: mp.mpf,
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Return central value, conservative half-width, minimum, and maximum over a box."""

    central = function(delta_sigma, delta_epsilon)
    corners = [
        function(delta_sigma + sign_sigma * sigma_uncertainty, delta_epsilon + sign_epsilon * epsilon_uncertainty)
        for sign_sigma, sign_epsilon in product((-1, 1), repeat=2)
    ]
    lower, upper = min(corners), max(corners)
    uncertainty = max(abs(lower - central), abs(upper - central))
    return central, uncertainty, lower, upper


def universal_data() -> dict[str, dict[str, object]]:
    """Assemble quoted CFT/RG data and scaling-derived exponents."""

    delta_sigma = mp.mpf("0.5181489")
    delta_sigma_uncertainty = mp.mpf("0.0000010")
    delta_epsilon = mp.mpf("1.412625")
    delta_epsilon_uncertainty = mp.mpf("0.000010")

    definitions: dict[str, tuple[str, Callable[[mp.mpf, mp.mpf], mp.mpf]]] = {
        "nu": ("nu = 1/(3 - Delta_epsilon)", lambda ds, de: 1 / (3 - de)),
        "eta": ("eta = 2*Delta_sigma - 1", lambda ds, de: 2 * ds - 1),
        "alpha": ("alpha = 2 - 3*nu", lambda ds, de: 2 - 3 / (3 - de)),
        "beta": ("beta = nu*(1 + eta)/2 = Delta_sigma/(3 - Delta_epsilon)", lambda ds, de: ds / (3 - de)),
        "gamma": ("gamma = nu*(2 - eta) = (3 - 2*Delta_sigma)/(3 - Delta_epsilon)", lambda ds, de: (3 - 2 * ds) / (3 - de)),
        "delta": ("delta = (5 - eta)/(1 + eta) = (3 - Delta_sigma)/Delta_sigma", lambda ds, de: (3 - ds) / ds),
    }

    data: dict[str, dict[str, object]] = {
        "Delta_sigma": {
            "value": _mpstr(delta_sigma),
            "uncertainty": _mpstr(delta_sigma_uncertainty),
            "quoted": "0.5181489(10)",
            "uncertainty_kind": "quoted numerical-bootstrap island half-width; not interpreted as a Gaussian standard deviation",
            "citation": "kos2016",
        },
        "Delta_epsilon": {
            "value": _mpstr(delta_epsilon),
            "uncertainty": _mpstr(delta_epsilon_uncertainty),
            "quoted": "1.412625(10)",
            "uncertainty_kind": "quoted numerical-bootstrap island half-width; not interpreted as a Gaussian standard deviation",
            "citation": "kos2016",
        },
    }

    displayed = {
        "nu": "0.629971(4)",
        "eta": "0.0362978(20)",
        "alpha": "0.110087(12)",
        "beta": "0.326419(3)",
        "gamma": "1.237075(9)",
        "delta": "4.789841(12)",
    }
    for name, (relation, function) in definitions.items():
        central, uncertainty, lower, upper = _rectangular_bound(
            function,
            delta_sigma,
            delta_epsilon,
            delta_sigma_uncertainty,
            delta_epsilon_uncertainty,
        )
        data[name] = {
            "value": _mpstr(central),
            "uncertainty": _mpstr(uncertainty),
            "interval": [_mpstr(lower), _mpstr(upper)],
            "quoted_or_rounded": displayed[name],
            "derivation": relation,
            "uncertainty_kind": "conservative rectangular propagation of the two quoted bootstrap marginal bounds; covariance unavailable",
            "citation": "kos2016",
        }

    # The downloaded precision-island paper does not determine the leading irrelevant
    # scalar.  This is therefore the best sourced omega estimate in the local full-text
    # corpus, and its visibly larger uncertainty is retained rather than improved by fiat.
    data["omega"] = {
        "value": "0.83",
        "uncertainty": "0.05",
        "quoted": "0.83(5)",
        "related_confluent_exponent": "Delta = 0.52(3) = omega*nu",
        "uncertainty_kind": "series-analysis uncertainty quoted by the source",
        "citation": "campostrini2002",
    }
    return data


def scaling_relation_checks(data: dict[str, dict[str, object]]) -> dict[str, object]:
    """Check internal identities and independent published observables."""

    values = {name: mp.mpf(entry["value"]) for name, entry in data.items() if "value" in entry}
    internal_relations = {
        "alpha = 2 - 3 nu": values["alpha"] - (2 - 3 * values["nu"]),
        "gamma = nu (2 - eta)": values["gamma"] - values["nu"] * (2 - values["eta"]),
        "beta = nu (1 + eta) / 2": values["beta"] - values["nu"] * (1 + values["eta"]) / 2,
        "delta = (5 - eta) / (1 + eta)": values["delta"] - (5 - values["eta"]) / (1 + values["eta"]),
    }

    # These are genuinely cross-source comparisons rather than a tautological check of
    # exponents derived from the same two dimensions.  Covariances are unavailable, so
    # quadrature pulls are diagnostics only and are labelled accordingly.
    observations = [
        ("alpha = 2 - 3 nu", mp.mpf("0.110"), mp.mpf("0.002"), values["alpha"], mp.mpf(data["alpha"]["uncertainty"]), "campostrini2002", "direct HT alpha 0.110(2), Table 1"),
        ("gamma = nu (2 - eta)", mp.mpf("1.23708"), mp.mpf("0.00033"), values["gamma"], mp.mpf(data["gamma"]["uncertainty"]), "ferrenberg2018", "Monte Carlo gamma, Eq. (46)"),
        ("beta = nu (1 + eta) / 2", mp.mpf("0.32630"), mp.mpf("0.00022"), values["beta"], mp.mpf(data["beta"]["uncertainty"]), "ferrenberg2018", "Monte Carlo beta, Eq. (47)"),
        ("delta = (5 - eta) / (1 + eta)", mp.mpf("4.7893"), mp.mpf("0.0008"), values["delta"], mp.mpf(data["delta"]["uncertainty"]), "campostrini2002", "reported delta 4.7893(8), Table 1; source marks it scaling-derived"),
    ]
    independent: list[dict[str, object]] = []
    for relation, observed, observed_uncertainty, predicted, predicted_uncertainty, citation, detail in observations:
        residual = observed - predicted
        combined = mp.sqrt(observed_uncertainty**2 + predicted_uncertainty**2)
        independent.append(
            {
                "relation": relation,
                "published_lhs": _mpstr(observed),
                "published_lhs_uncertainty": _mpstr(observed_uncertainty),
                "rhs_from_bootstrap_dimensions": _mpstr(predicted),
                "rhs_uncertainty": _mpstr(predicted_uncertainty),
                "residual_lhs_minus_rhs": _mpstr(residual),
                "combined_uncertainty_if_independent": _mpstr(combined),
                "absolute_pull_if_independent": _mpstr(abs(residual) / combined),
                "citation_lhs": citation,
                "citation_rhs": "kos2016",
                "detail": detail,
            }
        )

    return {
        "internal_bootstrap_derived": [
            {
                "relation": relation,
                "residual": _mpstr(residual),
                "note": "algebraic consistency check; zero is expected because the recommended exponents are derived from the two dimensions",
                "citation": "kos2016",
            }
            for relation, residual in internal_relations.items()
        ],
        "cross_source": independent,
    }


def amplitude_ratios() -> dict[str, dict[str, str]]:
    """Return three standard universal amplitude ratios with source definitions."""

    return {
        "A_plus_over_A_minus": {
            "value": "0.532",
            "uncertainty": "0.003",
            "quoted": "0.532(3)",
            "definition": "U0 = A_+/A_- for C_H = A_+ t^(-alpha) above and A_- |t|^(-alpha) below Tc",
            "citation": "campostrini2002",
            "location": "Table 1; definition in Table 7 and Eq. (5.1)",
        },
        "C_plus_over_C_minus": {
            "value": "4.76",
            "uncertainty": "0.02",
            "quoted": "4.76(2)",
            "definition": "U2 = C_+/C_- for the zero-field susceptibility amplitudes",
            "citation": "campostrini2002",
            "location": "Table 1; definition in Table 7",
        },
        "R_chi": {
            "value": "1.660",
            "uncertainty": "0.004",
            "quoted": "1.660(4)",
            "definition": "R_chi = C_+ B^(delta-1) / B_c^delta",
            "citation": "campostrini2002",
            "location": "Table 1; definition in Table 7; k=3 estimate in Table 8",
        },
    }


def _real_phi(phi: Callable[[mp.mpf], mp.mpf | mp.mpc], coupling: mp.mpf) -> mp.mpf:
    value = phi(coupling)
    real, imaginary = mp.re(value), abs(mp.im(value))
    tolerance = mp.mpf(10) ** (-(PRECISION - 12)) * max(1, abs(real))
    if imaginary > tolerance:
        raise ValueError(f"candidate free energy is non-real at K={coupling}: Im(phi)={imaginary}")
    if not mp.isfinite(real):
        raise ValueError(f"candidate free energy is non-finite at K={coupling}")
    return mp.mpf(real)


def _scan_curvature_peak(
    phi: Callable[[mp.mpf], mp.mpf], lower: mp.mpf, upper: mp.mpf, intervals: int
) -> tuple[mp.mpf, mp.mpf, list[mp.mpf], list[mp.mpf]]:
    step = (upper - lower) / intervals
    points = [lower + index * step for index in range(intervals + 1)]
    values = [phi(point) for point in points]
    curvatures = [
        abs((values[index + 1] - 2 * values[index] + values[index - 1]) / step**2)
        for index in range(1, intervals)
    ]
    peak_offset = max(range(len(curvatures)), key=curvatures.__getitem__)
    return points[peak_offset + 1], step, points, curvatures


def _locate_nearest_positive_singularity(
    phi: Callable[[mp.mpf], mp.mpf]
) -> dict[str, object]:
    center, step, points, curvatures = _scan_curvature_peak(
        phi, SCAN_LOW, SCAN_HIGH, INITIAL_INTERVALS
    )
    background = mp.mpf(median(curvatures))
    local_peaks: list[int] = []
    for index, curvature in enumerate(curvatures):
        left_ok = index == 0 or curvature >= curvatures[index - 1]
        right_ok = index == len(curvatures) - 1 or curvature >= curvatures[index + 1]
        if left_ok and right_ok and curvature >= mp.mpf("2.5") * background:
            local_peaks.append(index)
    if local_peaks:
        nearest_index = min(local_peaks, key=lambda index: points[index + 1])
        center = points[nearest_index + 1]
        initial_significance = curvatures[nearest_index] / background
    else:
        initial_significance = max(curvatures) / background

    history: list[dict[str, str]] = [
        {
            "center": _mpstr(center),
            "grid_step": _mpstr(step),
            "peak_over_median_curvature": _mpstr(initial_significance),
        }
    ]
    previous_center = center
    for _ in range(REFINEMENT_LEVELS):
        lower = center - mp.mpf("1.5") * step
        upper = center + mp.mpf("1.5") * step
        previous_center = center
        center, step, _, _ = _scan_curvature_peak(phi, lower, upper, REFINEMENT_INTERVALS)
        history.append({"center": _mpstr(center), "grid_step": _mpstr(step)})

    uncertainty = max(2 * step, abs(center - previous_center) + step)
    return {
        "K": center,
        "uncertainty": uncertainty,
        "interval": (center - uncertainty, center + uncertainty),
        "initial_significant_peak_count": len(local_peaks),
        "initial_peak_over_median_curvature": initial_significance,
        "history": history,
    }


def _five_point_curvature(
    phi: Callable[[mp.mpf], mp.mpf], coupling: mp.mpf, step: mp.mpf
) -> mp.mpf:
    values = [phi(coupling + offset * step) for offset in (-2, -1, 0, 1, 2)]
    return (-values[4] + 16 * values[3] - 30 * values[2] + 16 * values[1] - values[0]) / (12 * step**2)


def _linear_relative_rms(xs: list[mp.mpf], ys: list[mp.mpf]) -> mp.mpf:
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denominator = sum((value - mean_x) ** 2 for value in xs)
    if denominator == 0:
        return mp.inf
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denominator
    intercept = mean_y - slope * mean_x
    rms = mp.sqrt(sum((y - intercept - slope * x) ** 2 for x, y in zip(xs, ys)) / len(xs))
    scale = max(ys) - min(ys)
    return rms / scale if scale else mp.inf


def _extract_specific_heat_exponent(
    phi: Callable[[mp.mpf], mp.mpf], singularity: dict[str, object]
) -> dict[str, object]:
    critical = mp.mpf(singularity["K"])
    location_uncertainty = mp.mpf(singularity["uncertainty"])
    offsets = [critical * mp.mpf("0.08") / (2**level) for level in range(OFFSET_LEVELS)]
    side_data: dict[str, dict[str, object]] = {}
    final_effective: list[mp.mpf] = []
    previous_effective: list[mp.mpf] = []
    log_residuals: list[mp.mpf] = []
    power_residuals: list[mp.mpf] = []

    for side_name, sign in (("low_K", -1), ("high_K", 1)):
        curvatures = [
            _five_point_curvature(phi, critical + sign * offset, offset / 8)
            for offset in offsets
        ]
        differences = [
            curvatures[index + 1] - curvatures[index]
            for index in range(len(curvatures) - 1)
        ]
        effective = [
            mp.log(abs(differences[index + 1] / differences[index])) / mp.log(2)
            for index in range(len(differences) - 1)
        ]
        final_effective.append(effective[-1])
        previous_effective.append(effective[-2])

        tail_offsets = offsets[2:]
        tail_curvatures = curvatures[2:]
        log_residual = _linear_relative_rms([mp.log(offset) for offset in tail_offsets], tail_curvatures)
        power_residual = _linear_relative_rms(
            [offset ** (-ALPHA_CONSTRAINT) for offset in tail_offsets], tail_curvatures
        )
        log_residuals.append(log_residual)
        power_residuals.append(power_residual)
        side_data[side_name] = {
            "reduced_offsets_abs_K_minus_Kc": [_mpstr(offset / critical) for offset in offsets],
            "phi_second_derivative": [_mpstr(value) for value in curvatures],
            "successive_curvature_differences": [_mpstr(value) for value in differences],
            "effective_alpha_from_difference_ratios": [_mpstr(value) for value in effective],
            "log_fit_relative_rms": _mpstr(log_residual),
            "fixed_alpha_0.110_fit_relative_rms": _mpstr(power_residual),
        }

    estimate = sum(final_effective) / 2
    previous_estimate = sum(previous_effective) / 2
    side_asymmetry = abs(final_effective[0] - final_effective[1]) / 2
    convergence_drift = abs(estimate - previous_estimate)
    location_component = location_uncertainty / offsets[-1]
    numerical_floor = mp.mpf("0.002")
    uncertainty = max(side_asymmetry, convergence_drift, location_component, numerical_floor)

    mean_log_residual = sum(log_residuals) / len(log_residuals)
    mean_power_residual = sum(power_residuals) / len(power_residuals)
    logarithmic = (
        abs(estimate) <= uncertainty
        and abs(estimate) + uncertainty < mp.mpf("0.06")
        and mean_log_residual * mp.mpf("1.5") < mean_power_residual
    )
    classification = "logarithmic_alpha_0" if logarithmic else "power_law"
    lower, upper = estimate - uncertainty, estimate + uncertainty
    benchmark_lower = ALPHA_CONSTRAINT - ALPHA_CONSTRAINT_UNCERTAINTY
    benchmark_upper = ALPHA_CONSTRAINT + ALPHA_CONSTRAINT_UNCERTAINTY
    interval_overlap = max(lower, benchmark_lower) <= min(upper, benchmark_upper)
    compatible = classification == "power_law" and interval_overlap

    return {
        "classification": classification,
        "alpha_estimate": _mpstr(mp.mpf(0) if logarithmic else estimate),
        "raw_two_sided_alpha_estimate": _mpstr(estimate),
        "alpha_uncertainty": _mpstr(uncertainty),
        "alpha_interval": [_mpstr(lower), _mpstr(upper)],
        "error_budget": {
            "low_high_side_asymmetry": _mpstr(side_asymmetry),
            "last_level_convergence_drift": _mpstr(convergence_drift),
            "singularity_location_uncertainty_over_smallest_offset": _mpstr(location_component),
            "numerical_floor": _mpstr(numerical_floor),
            "rule": "maximum of the four displayed components",
        },
        "model_comparison": {
            "mean_log_fit_relative_rms": _mpstr(mean_log_residual),
            "mean_fixed_alpha_0.110_fit_relative_rms": _mpstr(mean_power_residual),
        },
        "three_dimensional_constraint": {
            "alpha": _mpstr(ALPHA_CONSTRAINT),
            "uncertainty": _mpstr(ALPHA_CONSTRAINT_UNCERTAINTY),
            "quoted": "0.110(1)",
            "citation": "campostrini2002",
            "compatible": compatible,
            "reason": (
                "power-law exponent interval overlaps the conservative 3-D interval"
                if compatible
                else "logarithmic alpha=0 or exponent interval does not overlap the conservative 3-D interval"
            ),
        },
        "side_data": side_data,
    }


def check_candidate(phi_callable: Callable[[mp.mpf], mp.mpf | mp.mpc], name: str) -> dict[str, object]:
    """Numerically test a candidate reduced free energy ``phi(K)``.

    The function scans the positive real interval recorded in the report, locates the
    nearest significant curvature peak without using a supplied ``K_c``, and extracts
    the specific-heat exponent from background-cancelling curvature differences on both
    sides.  A logarithm has asymptotically constant differences and hence effective
    alpha zero; a power law has difference ratios approaching ``2**alpha``.
    """

    mp.mp.dps = PRECISION
    cache: dict[mp.mpf, mp.mpf] = {}

    def cached(coupling: mp.mpf) -> mp.mpf:
        key = mp.mpf(coupling)
        if key not in cache:
            cache[key] = _real_phi(phi_callable, key)
        return cache[key]

    singularity = _locate_nearest_positive_singularity(cached)
    extraction = _extract_specific_heat_exponent(cached, singularity)
    return {
        "name": name,
        "status": "checked",
        "nearest_detected_positive_real_singularity": {
            "K": _mpstr(mp.mpf(singularity["K"])),
            "uncertainty": _mpstr(mp.mpf(singularity["uncertainty"])),
            "interval": [_mpstr(value) for value in singularity["interval"]],
            "search_interval": [_mpstr(SCAN_LOW), _mpstr(SCAN_HIGH)],
            "initial_grid_intervals": INITIAL_INTERVALS,
            "refinement_grid_intervals": REFINEMENT_INTERVALS,
            "refinement_levels": REFINEMENT_LEVELS,
            "initial_significant_peak_count": singularity["initial_significant_peak_count"],
            "initial_peak_over_median_curvature": _mpstr(mp.mpf(singularity["initial_peak_over_median_curvature"])),
            "history": singularity["history"],
            "method": "adaptive maximum of absolute centered second difference; earliest coarse peak above 2.5 times median curvature",
        },
        "specific_heat_singularity": extraction,
        "function_evaluations": len(cache),
    }


def _onsager_control(coupling: mp.mpf) -> mp.mpf:
    return onsager_free_energy(coupling, coupling)


def _zhang_2007_reduced(coupling: mp.mpf) -> mp.mpf:
    """Zhang's recorded double integral after exact evaluation of one angle."""

    coupling = mp.mpf(coupling)
    sinh_two, sinh_six = mp.sinh(2 * coupling), mp.sinh(6 * coupling)
    cosh_product = mp.cosh(2 * coupling) * mp.cosh(6 * coupling)
    # (A-b-c)(A+b+c) = (b*c-1)^2: cancellation-free at the claimed singularity.
    determinant_gap = (sinh_two * sinh_six - 1) ** 2 / (
        cosh_product + sinh_two + sinh_six
    )

    def integrated_angle(theta: mp.mpf) -> mp.mpf:
        a_minus_c = determinant_gap + 2 * sinh_two * mp.sin(theta / 2) ** 2
        a = a_minus_c + sinh_six
        radical = mp.sqrt(a_minus_c * (a_minus_c + 2 * sinh_six))
        return mp.log((a + radical) / 2)

    integral = mp.quad(integrated_angle, [0, mp.pi / 2, mp.pi])
    return mp.log(2) + integral / (2 * mp.pi)


def _dzhang_2021(coupling: mp.mpf) -> mp.mpf:
    """Direct stable transcription of the one-angle expression recorded in the JSON."""

    coupling = mp.mpf(coupling)
    dual = mp.atanh(mp.exp(-2 * coupling))
    difference = dual - coupling
    sinh_product = mp.sinh(2 * coupling) * mp.sinh(2 * difference)

    def integrated_angle(theta: mp.mpf) -> mp.mpf:
        if sinh_product >= 0:
            argument_minus_one = (
                mp.cosh(2 * (coupling - difference))
                - 1
                + 2 * sinh_product * mp.sin(theta / 2) ** 2
            )
        else:
            argument_minus_one = (
                mp.cosh(2 * (coupling + difference))
                - 1
                - 2 * sinh_product * mp.cos(theta / 2) ** 2
            )
        return mp.acosh(1 + argument_minus_one)

    integral = mp.quad(integrated_angle, [0, mp.pi / 2, mp.pi])
    return mp.log(2 * mp.sinh(2 * coupling)) / 2 + integral / (2 * mp.pi)


def _series_assessment() -> dict[str, object]:
    high_path = ROOT / "results" / "series" / "sc_ht_free_energy.json"
    low_path = ROOT / "results" / "series" / "sc_lt_free_energy.json"
    high = json.loads(high_path.read_text())
    low = json.loads(low_path.read_text())
    return {
        "high_temperature_free_energy": {
            "path": str(high_path.relative_to(ROOT)),
            "variable": high["data"]["variable"],
            "achieved_order": high["data"]["achieved_order"],
            "citation": "local exact artifact generated by experiments/e03_flm_series.py",
        },
        "low_temperature_free_energy": {
            "path": str(low_path.relative_to(ROOT)),
            "variable": low["data"]["variable"],
            "achieved_order": low["data"]["achieved_order"],
            "citation": "local exact artifact generated by experiments/e03_flm_series.py",
        },
        "testable_amplitude_ratios": [],
        "conclusion": "none",
        "reason": (
            "A+/A- would require controlled two-sided critical continuation and amplitude extraction; "
            "the short zero-field free-energy series do not support it. C+/C- and R_chi require "
            "field derivatives (susceptibility and magnetization amplitudes), which these artifacts do not contain."
        ),
    }


def _synthetic_power_law(coupling: mp.mpf) -> mp.mpf:
    critical = mp.mpf("0.35")
    return coupling**2 + abs(coupling - critical) ** (2 - ALPHA_CONSTRAINT)


def run_analysis(write_result: bool = True) -> dict[str, object]:
    mp.mp.dps = PRECISION
    data = universal_data()
    scaling = scaling_relation_checks(data)
    ratios = amplitude_ratios()
    series = _series_assessment()

    onsager = check_candidate(_onsager_control, "2D Onsager isotropic control")
    onsager_exact_kc = mp.log(1 + mp.sqrt(2)) / 2
    onsager["citation"] = "onsager1944"
    onsager["exact_Kc"] = _mpstr(onsager_exact_kc)

    claims = json.loads(CLAIMS_PATH.read_text())
    implementations: dict[str, Callable[[mp.mpf], mp.mpf]] = {
        "zhang2007": _zhang_2007_reduced,
        "dzhang2021": _dzhang_2021,
    }
    claimed_reports: list[dict[str, object]] = []
    for record in claims:
        key = record["key"]
        expression = record.get("free_energy_python")
        if expression is None:
            claimed_reports.append(
                {
                    "key": key,
                    "status": "not_testable_no_free_energy_callable",
                    "reason": "the reference record supplies a critical-point conjecture only",
                    "source_record": str(CLAIMS_PATH.relative_to(ROOT)),
                    "claimed_Kc": _mpstr(mp.mpf(str(record["claimed_kc_value"]))),
                    "citation": "zhang2007" if key == "rosengren1985_kc" else key,
                }
            )
            continue
        if key not in implementations:
            raise KeyError(f"no stable callable transcription for claimed solution {key}")
        report = check_candidate(implementations[key], key)
        located = mp.mpf(report["nearest_detected_positive_real_singularity"]["K"])
        claimed_kc = mp.mpf(str(record["claimed_kc_value"]))
        report.update(
            {
                "key": key,
                "citation": key,
                "source_record": str(CLAIMS_PATH.relative_to(ROOT)),
                "source_expression_sha256": hashlib.sha256(expression.encode()).hexdigest(),
                "callable_transcription": (
                    "exact analytic evaluation of one angular integral from the recorded expression"
                    if key == "zhang2007"
                    else "stable direct transcription of the recorded one-angular-integral expression"
                ),
                "claimed_Kc_not_used_by_locator": _mpstr(claimed_kc),
                "located_minus_claimed_Kc": _mpstr(located - claimed_kc),
            }
        )
        claimed_reports.append(report)

    calibration = check_candidate(_synthetic_power_law, "synthetic alpha=0.110 power-law calibration")
    calibration["provenance"] = "constructed by this script to verify that the extractor distinguishes a small power from a logarithm"
    calibration["benchmark_citation"] = "campostrini2002"

    onsager_interval = [
        mp.mpf(value)
        for value in onsager["nearest_detected_positive_real_singularity"]["interval"]
    ]
    evaluable_claims = [report for report in claimed_reports if report["status"] == "checked"]
    internal_residuals = [
        abs(mp.mpf(item["residual"])) for item in scaling["internal_bootstrap_derived"]
    ]
    cross_source_pulls = [
        mp.mpf(item["absolute_pull_if_independent"]) for item in scaling["cross_source"]
    ]
    calibration_extraction = calibration["specific_heat_singularity"]
    calibration_estimate = mp.mpf(calibration_extraction["alpha_estimate"])
    calibration_uncertainty = mp.mpf(calibration_extraction["alpha_uncertainty"])

    checks = [
        {
            "name": "bootstrap-derived exponents satisfy all four scaling identities",
            "passed": max(internal_residuals) < mp.mpf("1e-28"),
            "detail": f"maximum absolute residual {_mpstr(max(internal_residuals))}",
        },
        {
            "name": "cross-source exponent observations are consistent with scaling predictions",
            "passed": max(cross_source_pulls) < 1,
            "detail": f"maximum diagnostic pull {_mpstr(max(cross_source_pulls))}",
        },
        {
            "name": "Onsager control singularity location contains exact square-lattice Kc",
            "passed": onsager_interval[0] <= onsager_exact_kc <= onsager_interval[1],
            "detail": (
                f"located interval [{_mpstr(onsager_interval[0])}, {_mpstr(onsager_interval[1])}], "
                f"exact {_mpstr(onsager_exact_kc)}"
            ),
        },
        {
            "name": "Onsager control is identified as logarithmic alpha=0",
            "passed": (
                onsager["specific_heat_singularity"]["classification"] == "logarithmic_alpha_0"
                and mp.mpf(onsager["specific_heat_singularity"]["alpha_interval"][0]) <= 0
                <= mp.mpf(onsager["specific_heat_singularity"]["alpha_interval"][1])
            ),
            "detail": (
                f"classification {onsager['specific_heat_singularity']['classification']}, "
                f"raw alpha {onsager['specific_heat_singularity']['raw_two_sided_alpha_estimate']} "
                f"+/- {onsager['specific_heat_singularity']['alpha_uncertainty']}"
            ),
        },
        {
            "name": "extractor distinguishes a synthetic alpha=0.110 power from a logarithm",
            "passed": (
                calibration_extraction["classification"] == "power_law"
                and abs(calibration_estimate - ALPHA_CONSTRAINT)
                <= calibration_uncertainty + mp.mpf("0.002")
            ),
            "detail": (
                f"classification {calibration_extraction['classification']}, alpha "
                f"{calibration_extraction['alpha_estimate']} +/- {calibration_extraction['alpha_uncertainty']}"
            ),
        },
        {
            "name": "every claimed-solution reference record is accounted for",
            "passed": len(claimed_reports) == len(claims),
            "detail": f"{len(claimed_reports)} of {len(claims)} records reported",
        },
        {
            "name": "all claimed free-energy formulas fail the 3D alpha constraint",
            "passed": (
                len(evaluable_claims) == 2
                and all(
                    report["specific_heat_singularity"]["classification"] == "logarithmic_alpha_0"
                    and not report["specific_heat_singularity"]["three_dimensional_constraint"]["compatible"]
                    for report in evaluable_claims
                )
            ),
            "detail": ", ".join(
                f"{report['key']}: {report['specific_heat_singularity']['classification']}"
                for report in evaluable_claims
            ),
        },
        {
            "name": "current exact zero-field series support no amplitude-ratio extraction",
            "passed": series["testable_amplitude_ratios"] == [],
            "detail": "none; the available observables and orders are insufficient",
        },
    ]

    result: dict[str, object] = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision_decimal_digits": PRECISION,
            "arithmetic": "mpmath for transcendental numerics; exact decimal strings for literature inputs",
            "source_manifest": "sources/manifest.yaml",
            "literature_matrix": "notes/literature_matrix.md",
        },
        "data": {
            "citations": CITATIONS,
            "universal_data": data,
            "scaling_relation_checks": scaling,
            "specific_heat_constraint": {
                "alpha": _mpstr(ALPHA_CONSTRAINT),
                "uncertainty": _mpstr(ALPHA_CONSTRAINT_UNCERTAINTY),
                "quoted": "0.110(1)",
                "citation": "campostrini2002",
                "note": "conservative checker interval requested by the task; wider than the bootstrap-derived interval",
            },
            "universal_amplitude_ratios": ratios,
            "series_testability": series,
            "candidate_checks": {
                "onsager_2d_control": onsager,
                "claimed_solutions": claimed_reports,
                "extractor_calibration": calibration,
            },
        },
        "checks": checks,
    }

    if write_result:
        RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
        RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    result = run_analysis(write_result=True)
    for check in result["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"{status}: {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in result["checks"])
    print(("PASS" if passed else "FAIL") + f": wrote {RESULT_PATH.relative_to(ROOT)}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
