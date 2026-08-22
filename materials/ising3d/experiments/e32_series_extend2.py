"""Profile, extend, independently witness, and analyze the cubic Ising series.

Run from the repository root:
    .venv/bin/python experiments/e32_series_extend2.py

The resource profile is printed before either finite-lattice computation starts.
"""

from __future__ import annotations

import errno
import json
import math
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import prod
from pathlib import Path

import mpmath as mp
from sympy import Matrix, Rational

import ising.series as flm
from e18_series_extend import _hybrid_flm_engine, _series_payload
from e19_series_analysis import _prefix_estimate
from ising.series import high_temperature_free_energy, low_temperature_free_energy
from ising.series.analysis import differential_approximant_scan, even_to_squared_variable

SCRIPT = "experiments/e32_series_extend2.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "series"
NOTE_PATH = ROOT / "notes" / "series_extension2.md"
FROZEN_HT_PATH = RESULT_DIR / "extended_sc_ht_free_energy.json"
FROZEN_LT_PATH = RESULT_DIR / "extended_sc_lt_free_energy.json"
HT_PATH = RESULT_DIR / "extended2_sc_ht_free_energy.json"
LT_PATH = RESULT_DIR / "extended2_sc_lt_free_energy.json"
ANALYSIS_PATH = RESULT_DIR / "dfinite_holdout.json"
HT_TARGET = 22
LT_TARGET = 32
PRECISION = 80
CRT_CROSS_SECTION_GUARD = 22
SOURCE_KEY = "guttmann_enting1993"
SOURCE_PATH = "sources/fulltext/num_guttmann_enting1993.pdf"
SOURCE_SHA256 = "b866f359ca4bc0f75d534f8561cc2a286edd0680a0b37e3f7d8a741830609bff"


def _fraction_strings(values) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _load_payload(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _physical_memory_bytes() -> int | None:
    if sys.platform != "darwin":
        return None
    try:
        return int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip())
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


def _open_bonds(shape: tuple[int, int, int]) -> int:
    return sum((shape[axis] - 1) * prod(shape[:axis] + shape[axis + 1 :]) for axis in range(3))


def _canonical_shapes(side: str, order: int) -> tuple[tuple[int, int, int], ...]:
    if side == "HT":
        budget = order // 2
        return tuple(
            shape
            for shape in combinations_with_replacement(range(1, budget + 2), 3)
            if sum(value - 1 for value in shape) <= budget
        )
    if side == "LT":
        budget = (order + 6) // 4
        return tuple(
            shape
            for shape in combinations_with_replacement(range(1, budget + 1), 3)
            if sum(shape) <= budget
        )
    raise ValueError(f"unknown series side {side}")


def _ordered_shape_count(side: str, order: int) -> int:
    if side == "HT":
        budget = order // 2
        return sum(
            1
            for shape in product(range(1, budget + 2), repeat=3)
            if sum(value - 1 for value in shape) <= budget
        )
    budget = (order + 6) // 4
    return sum(1 for shape in product(range(1, budget + 1), repeat=3) if sum(shape) <= budget)


def _profile(side: str, order: int) -> dict:
    canonical = _canonical_shapes(side, order)
    budget = order // 2 if side == "HT" else (order + 6) // 4
    records = []
    for shape in canonical:
        sites = prod(shape)
        cross_section = prod(shape[:-1])
        open_bonds = _open_bonds(shape)
        full_degree = open_bonds if side == "HT" else 6 * sites - open_bonds
        peak_bytes = 3 * (1 << cross_section) * (full_degree + 1) * 8
        records.append(
            {
                "shape": list(shape),
                "sites": sites,
                "cross_section_sites": cross_section,
                "full_polynomial_degree": full_degree,
                "estimated_peak_bytes": peak_bytes,
                "within_crt_cross_section_guard": cross_section <= CRT_CROSS_SECTION_GUARD,
            }
        )
    memory_wall = max(records, key=lambda row: row["estimated_peak_bytes"])
    maximum_sites = max(records, key=lambda row: row["sites"])
    maximum_cross = max(records, key=lambda row: row["cross_section_sites"])
    if side == "HT":
        frontier = [list(shape) for shape in canonical if sum(value - 1 for value in shape) == budget]
        bound = "sum_i(side_i-1) <= order/2"
    else:
        frontier = [list(shape) for shape in canonical if sum(shape) == budget]
        bound = "sum_i(side_i) <= floor((order+6)/4)"
    return {
        "side": side,
        "order": order,
        "box_budget": budget,
        "order_bound": bound,
        "ordered_box_count": _ordered_shape_count(side, order),
        "canonical_box_count": len(canonical),
        "canonical_boxes": [list(shape) for shape in canonical],
        "new_budget_frontier_boxes": frontier,
        "maximum_sites_box": maximum_sites,
        "maximum_cross_section_box": maximum_cross,
        "memory_wall_box": memory_wall,
        "crt_cross_section_guard": CRT_CROSS_SECTION_GUARD,
        "reachable_by_current_guard": all(row["within_crt_cross_section_guard"] for row in records),
    }


def _print_profiles(profiles: list[dict], physical_memory: int | None) -> None:
    print("PROFILE FIRST (no coefficient computation has started)", flush=True)
    for item in profiles:
        wall = item["memory_wall_box"]
        print(
            "PROFILE "
            f"{item['side']} order={item['order']}: canonical={item['canonical_box_count']}, "
            f"frontier={item['new_budget_frontier_boxes']}, wall={wall['shape']} "
            f"N={wall['sites']} cross={wall['cross_section_sites']} "
            f"peak={wall['estimated_peak_bytes']} bytes, "
            f"guard_ok={item['reachable_by_current_guard']}",
            flush=True,
        )
    print(
        "PROFILE DECISION: extend LT through x^32 first (two new nonzero even "
        "orders at a 1,339,392-byte wall); then attempt bonus HT v^22. "
        "HT v^24 is blocked by cross-section 25 > 22 and its hypothetical "
        "242,397,216,768-byte propagation peak.",
        flush=True,
    )
    if physical_memory is not None:
        print(f"PROFILE physical_memory_bytes={physical_memory}", flush=True)


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def _clear_series_caches() -> None:
    flm._free_broken_polynomial.cache_clear()
    flm._plus_broken_polynomial.cache_clear()
    flm._ht_box_log.cache_clear()
    flm._lt_box_log.cache_clear()


def _formal_log(values: list[int], order: int) -> tuple[Fraction, ...]:
    if len(values) < order + 1 or values[0] != 1:
        raise ValueError("formal logarithm requires a complete unit-constant prefix")
    result = [Fraction(0) for _ in range(order + 1)]
    for degree in range(1, order + 1):
        convolution = sum(
            (index * result[index] * values[degree - index] for index in range(1, degree)),
            Fraction(0),
        )
        result[degree] = Fraction(values[degree]) - convolution / degree
    return tuple(result)


def _external_primary_witness() -> dict:
    # Instantiated only after local series derivations finish: validation, never input.
    published_lambda = {
        0: 1, 3: 1, 4: 0, 5: 3, 6: -3, 7: 15, 8: -30, 9: 101,
        10: -261, 11: 807, 12: -2308, 13: 7065, 14: -21171,
        15: 65337, 16: -200934,
    }
    lambda_series = [published_lambda.get(index, 0) for index in range(17)]
    published_a = {
        0: 1, 4: 3, 6: 22, 8: 192, 10: 2046, 12: 24853,
        14: 329334, 16: 4649601, 18: 68884356, 20: 1059830112,
        22: 16809862992,
    }
    a_series = [published_a.get(index, 0) for index in range(23)]
    log_lambda = _formal_log(lambda_series, 16)
    log_a = _formal_log(a_series, 22)
    return {
        "provenance_class": "external_primary_source",
        "manifest_key": SOURCE_KEY,
        "title": "Series studies of the Potts model. I. The simple cubic Ising model",
        "authors": ["A. J. Guttmann", "I. G. Enting"],
        "journal": "Journal of Physics A: Mathematical and General 26 (1993) 807-821",
        "doi": "10.1088/0305-4470/26/4/010",
        "arxiv": "hep-lat/9212032",
        "local_path": SOURCE_PATH,
        "sha256": SOURCE_SHA256,
        "exact_location": "Table 2, arXiv PDF page 21 of 24; normalizations in equations (9) and (17)",
        "raw_LT_partition_coefficients_lambda_n": {
            str(index): str(value) for index, value in published_lambda.items()
        },
        "raw_HT_interaction_partition_coefficients_a_n": {
            str(index): str(value) for index, value in published_a.items()
        },
        "exact_formal_log_witnesses": {
            "x^30_equals_log_Lambda0_u^15": str(log_lambda[15]),
            "x^32_equals_log_Lambda0_u^16": str(log_lambda[16]),
            "v^22_interaction_equals_log_Phi_v^22": str(log_a[22]),
            "v^22_total_adds_3_over_22": str(log_a[22] + Fraction(3, 22)),
        },
        "coefficient_source_role": (
            "validation only; locally derived coefficients are computed before this "
            "table is instantiated and are never seeded from these values"
        ),
    }


def _sympy_row(series: tuple[Fraction, ...], n: int, order: int, degree: int):
    row = []
    for derivative_order in range(order + 1):
        for shift in range(degree + 1):
            source = n - shift
            value = series[source] * (source**derivative_order) if source >= 0 else Fraction(0)
            row.append(Rational(value.numerator, value.denominator))
    return row


def _primitive_fraction_vector(vector) -> tuple[Fraction, ...]:
    fractions = [Fraction(int(Rational(value).p), int(Rational(value).q)) for value in vector]
    common_denominator = math.lcm(*(value.denominator for value in fractions))
    integers = [value.numerator * (common_denominator // value.denominator) for value in fractions]
    common_divisor = 0
    for value in integers:
        common_divisor = math.gcd(common_divisor, abs(value))
    if common_divisor:
        integers = [value // common_divisor for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(Fraction(value) for value in integers)


def _full_dfinite_scan(series: tuple[Fraction, ...]) -> dict:
    coefficient_count = len(series)
    pairs = []
    for order in range(coefficient_count):
        for degree in range(coefficient_count):
            unknown_count = (order + 1) * (degree + 1)
            if unknown_count > coefficient_count:
                continue
            matrix = Matrix([_sympy_row(series, n, order, degree) for n in range(coefficient_count)])
            rank = int(matrix.rank())
            pairs.append(
                {
                    "order": order,
                    "degree": degree,
                    "parameter_budget": unknown_count,
                    "data_count": coefficient_count,
                    "rank": rank,
                    "nullity": unknown_count - rank,
                    "consistent_nonzero_relation": rank < unknown_count,
                }
            )
    smallest_vacuous_budget = coefficient_count + 1
    vacuous_pairs = [
        [order, degree]
        for order in range(smallest_vacuous_budget)
        for degree in range(smallest_vacuous_budget)
        if (order + 1) * (degree + 1) == smallest_vacuous_budget
    ]
    found = [item for item in pairs if item["consistent_nonzero_relation"]]
    return {
        "analyzed_series": "G(z)=(phi_HT-log(2))/z with z=v^2",
        "equivalence_note": (
            "Writing F=z*G conjugates theta on F to theta+1 on G, so Euler-ODE "
            "order and polynomial degree, hence the parameter budget, are preserved "
            "in both directions"
        ),
        "operator_basis": "Euler theta=z*d/dz",
        "ansatz": "sum_j Q_j(z) theta^j G(z)=0, deg Q_j<=degree",
        "known_coefficient_count": coefficient_count,
        "tested_pair_count": len(pairs),
        "tested_pairs": pairs,
        "consistent_relation_count": len(found),
        "largest_fully_nonvacuous_parameter_budget": coefficient_count,
        "smallest_vacuous_parameter_budget": smallest_vacuous_budget,
        "first_vacuous_order_degree_pairs": vacuous_pairs,
        "conclusion": (
            f"no D-finite recurrence of budget <= {coefficient_count} is consistent with the known coefficients"
            if not found
            else "at least one searched finite-budget relation is consistent"
        ),
        "scope": "finite exact scan only; it does not disprove D-finiteness at larger order or degree",
    }


def _predict_ode(
    training: tuple[Fraction, ...],
    total_count: int,
    order: int,
    degree: int,
    vector: tuple[Fraction, ...],
) -> tuple[Fraction, ...] | None:
    predicted = list(training)
    for n in range(len(training), total_count):
        leading = sum(vector[j * (degree + 1)] * (n**j) for j in range(order + 1))
        if leading == 0:
            return None
        remainder = Fraction(0)
        for j in range(order + 1):
            for shift in range(1, degree + 1):
                source = n - shift
                if source >= 0:
                    remainder += (
                        vector[j * (degree + 1) + shift]
                        * (source**j)
                        * predicted[source]
                    )
        predicted.append(-remainder / leading)
    return tuple(predicted)


def _ode_holdouts(series: tuple[Fraction, ...]) -> list[dict]:
    results = []
    total_count = len(series)
    first_m = max(5, total_count - 6)
    for m in range(first_m, total_count):
        candidates = []
        interpolation_budget = m + 1
        for order in range(interpolation_budget):
            for degree in range(interpolation_budget):
                unknown_count = (order + 1) * (degree + 1)
                if unknown_count != interpolation_budget:
                    continue
                matrix = Matrix([_sympy_row(series, n, order, degree) for n in range(m)])
                nullspace = matrix.nullspace()
                if len(nullspace) != 1:
                    continue
                vector = _primitive_fraction_vector(nullspace[0])
                prediction = _predict_ode(series[:m], total_count, order, degree, vector)
                if prediction is None:
                    candidates.append(
                        {
                            "order": order,
                            "degree": degree,
                            "parameter_budget": unknown_count,
                            "predictive": False,
                            "reason": "coefficient multiplying the first held-out term vanishes",
                        }
                    )
                    continue
                errors = tuple(prediction[index] - series[index] for index in range(m, total_count))
                candidates.append(
                    {
                        "order": order,
                        "degree": degree,
                        "parameter_budget": unknown_count,
                        "predictive": True,
                        "operator_vector": _fraction_strings(vector),
                        "predicted_remaining": _fraction_strings(prediction[m:]),
                        "actual_remaining": _fraction_strings(series[m:]),
                        "signed_errors": _fraction_strings(errors),
                        "absolute_errors": _fraction_strings(abs(value) for value in errors),
                        "relative_errors": _fraction_strings(
                            abs(error / actual) if actual else Fraction(0)
                            for error, actual in zip(errors, series[m:])
                        ),
                        "all_holdouts_exact": all(value == 0 for value in errors),
                    }
                )
        predictive = [item for item in candidates if item.get("predictive")]
        selected = (
            min(
                predictive,
                key=lambda item: (
                    abs(item["order"] - item["degree"]),
                    item["order"],
                    item["degree"],
                ),
            )
            if predictive
            else None
        )
        results.append(
            {
                "m_training_coefficients": m,
                "held_out_coefficients": total_count - m,
                "fit_budget": interpolation_budget,
                "selection_rule": (
                    "among prefix-only predictive interpolants, minimize |order-degree|, "
                    "then lexicographic (order,degree); holdouts are not consulted"
                ),
                "selected_pair": [selected["order"], selected["degree"]] if selected else None,
                "selected_absolute_errors": selected["absolute_errors"] if selected else [],
                "selected_relative_errors": selected["relative_errors"] if selected else [],
                "candidates": candidates,
            }
        )
    return results


def _predict_differential(fit, training: tuple[Fraction, ...], total_count: int):
    q1 = tuple(Fraction(value) for value in fit.derivative_polynomial)
    q0 = tuple(Fraction(value) for value in fit.function_polynomial)
    forcing = tuple(Fraction(value) for value in fit.inhomogeneous_polynomial)
    predicted = list(training)
    for n in range(len(training), total_count):
        equation_degree = n - 1
        leading = n * q1[0]
        if leading == 0:
            return None
        right = forcing[equation_degree] if equation_degree < len(forcing) else Fraction(0)
        for shift in range(1, len(q1)):
            source = n - shift
            if source >= 0:
                right -= q1[shift] * source * predicted[source]
        for shift in range(len(q0)):
            source = equation_degree - shift
            if source >= 0:
                right -= q0[shift] * predicted[source]
        predicted.append(right / leading)
    return tuple(predicted)


def _differential_holdouts(series: tuple[Fraction, ...]) -> list[dict]:
    results = []
    total_count = len(series)
    for m in range(max(7, total_count - 5), total_count):
        fits = [
            fit
            for fit in differential_approximant_scan(series[:m], precision=PRECISION, balanced=True)
            if 0 < fit.singularity < 1 and 0 < fit.singular_exponent < 4
        ]
        models = []
        for fit in fits:
            prediction = _predict_differential(fit, series[:m], total_count)
            if prediction is None:
                continue
            errors = tuple(prediction[index] - series[index] for index in range(m, total_count))
            models.append(
                {
                    "degrees_Q1_Q0_P": [
                        fit.derivative_degree,
                        fit.function_degree,
                        fit.inhomogeneous_degree,
                    ],
                    "z_singularity": mp.nstr(fit.singularity, 35),
                    "implied_K_singularity": mp.nstr(mp.atanh(mp.sqrt(fit.singularity)), 35),
                    "predicted_remaining": _fraction_strings(prediction[m:]),
                    "actual_remaining": _fraction_strings(series[m:]),
                    "signed_errors": _fraction_strings(errors),
                    "absolute_errors": _fraction_strings(abs(value) for value in errors),
                    "relative_errors": _fraction_strings(
                        abs(error / actual) if actual else Fraction(0)
                        for error, actual in zip(errors, series[m:])
                    ),
                    "all_holdouts_exact": all(value == 0 for value in errors),
                }
            )
        selected = (
            min(
                models,
                key=lambda item: (
                    max(item["degrees_Q1_Q0_P"]) - min(item["degrees_Q1_Q0_P"]),
                    item["degrees_Q1_Q0_P"],
                ),
            )
            if models
            else None
        )
        results.append(
            {
                "m_training_coefficients": m,
                "held_out_coefficients": total_count - m,
                "selection_rule": (
                    "balanced first-order inhomogeneous approximants with a training-only "
                    "positive physical pole and exponent in (0,4); minimize degree spread "
                    "then lexicographic degrees; holdouts are not consulted"
                ),
                "admissible_model_count": len(models),
                "selected_degrees": selected["degrees_Q1_Q0_P"] if selected else None,
                "selected_absolute_errors": selected["absolute_errors"] if selected else [],
                "selected_relative_errors": selected["relative_errors"] if selected else [],
                "models": models,
            }
        )
    return results


def _critical_estimate(coefficients: tuple[Fraction, ...]) -> dict:
    achieved_order = len(coefficients) - 1
    orders = tuple(order for order in (16, 18, 20, 22) if order <= achieved_order)
    prefixes = [_prefix_estimate(coefficients, order) for order in orders]
    current = prefixes[-1]
    all_k = [value for result in prefixes for value in result["k_estimates"]]
    all_alpha = [value for result in prefixes for value in result["alpha_estimates"]]
    central_k = current["central_k"]
    central_alpha = current["central_alpha"]
    k_uncertainty = max(abs(value - central_k) for value in all_k)
    alpha_uncertainty = max(abs(value - central_alpha) for value in all_alpha)
    return {
        "benchmark_used_in_fit_or_selection": False,
        "selection_rule": (
            "at the longest series, median of highest-complexity nondefective "
            "near-diagonal Dlog-Pade poles and balanced first-order inhomogeneous "
            "differential-approximant poles; only series coefficients enter"
        ),
        "uncertainty_rule": (
            "maximum displacement from the longest-series central value among every "
            "retained series-only approximant at truncations "
            + ", ".join(f"v^{order}" for order in orders)
        ),
        "K_c": mp.nstr(central_k, 35),
        "K_c_uncertainty": mp.nstr(k_uncertainty, 35),
        "v_c": mp.nstr(current["central_v"], 35),
        "alpha": mp.nstr(central_alpha, 35),
        "alpha_uncertainty": mp.nstr(alpha_uncertainty, 35),
        "truncation_stability": [
            {
                "v_order": result["order"],
                "central_K_c": mp.nstr(result["central_k"], 35),
                "central_alpha": mp.nstr(result["central_alpha"], 35),
                "dlog_count": len(result["dlog"]),
                "differential_count": len(result["differential"]),
                "all_K_estimates": [mp.nstr(value, 35) for value in result["k_estimates"]],
            }
            for result in prefixes
        ],
    }


def _series_data_from_frozen(payload: dict) -> dict:
    data = dict(payload["data"])
    data["extension_status"] = "not attempted successfully; frozen prefix retained"
    return data


def _write_json(path: Path, provenance: dict, data: dict, checks: list[dict]) -> None:
    payload = {"provenance": provenance, "data": data, "checks": checks}
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {path.relative_to(ROOT)}", flush=True)


def _human_bytes(value: int) -> str:
    return f"{value / (1 << 30):.3f} GiB" if value >= 1 << 30 else f"{value / (1 << 20):.3f} MiB"


def _write_note(profiles: list[dict], lt_data: dict, ht_data: dict, analysis: dict, witness: dict) -> None:
    by_key = {(item["side"], item["order"]): item for item in profiles}
    profile_rows = []
    for side, order in (("HT", 22), ("HT", 24), ("LT", 30), ("LT", 32)):
        item = by_key[(side, order)]
        wall = item["memory_wall_box"]
        profile_rows.append(
            f"| {side} `{('v' if side == 'HT' else 'x')}^{order}` | "
            f"{item['canonical_box_count']} | `{item['new_budget_frontier_boxes']}` | "
            f"`{wall['shape']}` | {wall['sites']} | {wall['cross_section_sites']} | "
            f"{wall['estimated_peak_bytes']} ({_human_bytes(wall['estimated_peak_bytes'])}) | "
            f"{'yes' if item['reachable_by_current_guard'] else 'no'} |"
        )
    lt_coefficients = lt_data["coefficients"]
    ht_coefficients = ht_data["coefficients"]
    ode_rows = []
    for row in analysis["strict_holdout"]["d_finite_interpolants"]:
        ode_rows.append(
            f"| {row['m_training_coefficients']} | {row['held_out_coefficients']} | "
            f"`{row['selected_pair']}` | `{row['selected_absolute_errors']}` | "
            f"`{row['selected_relative_errors']}` |"
        )
    da_rows = []
    for row in analysis["strict_holdout"]["differential_approximants"]:
        da_rows.append(
            f"| {row['m_training_coefficients']} | {row['held_out_coefficients']} | "
            f"`{row['selected_degrees']}` | `{row['selected_absolute_errors']}` | "
            f"`{row['selected_relative_errors']}` |"
        )
    kc = analysis["critical_estimate"]
    benchmark = analysis["final_comparison_only"]
    ht_reached = ht_data["achieved_order"] >= 22
    note = f"""# Series extension, resource wall, and strict holdout analysis

[COMPUTATION] The resource profile below was printed before either coefficient computation. The finite-lattice computation then extended the low-temperature reduced free energy through `x^32`; the same run {'reached `v^22`' if ht_reached else 'did not reach `v^22`'} on the high-temperature side.

## Profile first

[LEMMA] The HT box bound is `sum(side_i-1) <= order/2`; the LT bound is `sum(side_i) <= floor((order+6)/4)`. Sorting a box puts its largest side in the transfer direction, so the transfer cross-section is the product of its two smallest sides. The memory column is the CRT engine's conservative three-array estimate `3*2^cross_section*(full_degree+1)*8` bytes.

| target | canonical boxes | new budget-frontier canonical boxes | wall box | N | cross-section | peak memory | guard <=22 |
|---|---:|---|---|---:|---:|---:|---|
{chr(10).join(profile_rows)}

[COMPUTATION] `HT v^22` fits the current guard at cross-section 20. `HT v^24` does not: the required canonical `5x5x5` box has cross-section 25, beyond the hard guard 22, and its hypothetical three-array peak is 242,397,216,768 bytes (225.750 GiB). Thus `v^24` is the quantified resource wall.

[COMPUTATION] The frozen `x^28` FLM artifact actually has canonical wall `2x3x3`, cross-section 6. The 20-spin item stored in the frozen HT artifact is explicitly a separate `4x5x2` engine probe, not an LT-required box. This corrects the resource-wall inventory without changing any coefficient.

## Exact extensions and witnesses

[COMPUTATION] The locally derived LT coefficients are `c_30={lt_coefficients[30]}` and `c_32={lt_coefficients[32]}` in `phi=3K+sum c_n x^n`. They were obtained from the local plus-boundary finite-lattice calculation before the external table was instantiated. Every coefficient through frozen order `x^28` agrees exactly.

[EXTERNAL] Guttmann and Enting, *Series studies of the Potts model. I. The simple cubic Ising model*, J. Phys. A 26 (1993) 807-821, DOI `10.1088/0305-4470/26/4/010`, Table 2 (arXiv PDF page 21/24), print `lambda_15=65337` and `lambda_16=-200934` for `Lambda_0(u)=sum lambda_n u^n`, `u=x^2`. An exact local formal logarithm of their printed table gives `[{witness['exact_formal_log_witnesses']['x^30_equals_log_Lambda0_u^15']}, {witness['exact_formal_log_witnesses']['x^32_equals_log_Lambda0_u^16']}]`, agreeing coefficient-by-coefficient with the local `x^30,x^32` result. This is an external primary-source witness, not the source of the local coefficients.

[LEMMA] The internal `4x4xc` periodic-slab route in `experiments/e26_lt_series_independent.py` is rigorous only below its first transverse wrapping term `x^(4L)`. To witness through `x^32` without wrapping requires `4L>32`, hence `L>=9` in both periodic transverse directions and an 81-spin cross-section. That exceeds the current guard 22. Therefore no independent internal slab witness is available for `x^30,x^32`; the honest witness status is external-primary-source only.

[COMPUTATION] The locally derived HT `v^22` coefficient is `{ht_coefficients[22] if len(ht_coefficients) > 22 else 'not reached'}`. When reached, the same Table 2 value `a_22=16809862992` for `Phi(v)` gives the exact formal-log total-free-energy witness `{witness['exact_formal_log_witnesses']['v^22_total_adds_3_over_22']}` after adding the known `3/22` from `3 log(cosh K)`. The recorded process peak RSS was {ht_data.get('resource_usage', {}).get('process_peak_rss_bytes', 'not available')} bytes against the propagation-array estimate {by_key[('HT', 22)]['memory_wall_box']['estimated_peak_bytes']} bytes.

## Strict holdout prediction

[COMPUTATION] For the homogeneous Euler-ODE holdout, each model is fitted only to the first `m` coefficients of `G(z)=(phi-log 2)/z`, `z=v^2`; all remaining coefficients are recursively predicted. The fixed selector minimizes `|order-degree|` and then uses lexicographic degrees, without seeing any holdout.

| m | held out | selected (order,degree) | absolute errors | relative errors |
|---:|---:|---|---|---|
{chr(10).join(ode_rows)}

[COMPUTATION] Balanced first-order inhomogeneous differential approximants were likewise fitted to only the first `m` coefficients of `F(z)=phi-log 2`; the table reports genuine recursive predictions of every remaining coefficient.

| m | held out | selected `(deg Q1,deg Q0,deg P)` | absolute errors | relative errors |
|---:|---:|---|---|---|
{chr(10).join(da_rows)}

## Finite D-finiteness budget

[COMPUTATION] Every `(order,degree)` pair with parameter budget `(order+1)(degree+1) <= {analysis['d_finiteness']['largest_fully_nonvacuous_parameter_budget']}` was tested by exact rational rank on all known coefficients of `G`. The result is: **{analysis['d_finiteness']['conclusion']}**. Budget {analysis['d_finiteness']['smallest_vacuous_parameter_budget']} is the first vacuous budget because it has more operator coefficients than the {analysis['d_finiteness']['known_coefficient_count']} data. This finite negative result is not a proof that the full series is non-D-finite.

## Untuned series-only critical estimate

[COMPUTATION] With the model family and pole filters fixed without external input, the extended HT series gives `K_c={kc['K_c']} +/- {kc['K_c_uncertainty']}`. The uncertainty is the largest displacement from the longest-series median among every retained approximant at the recorded truncations; it is a method-spread envelope, not a statistical confidence interval.

[EXTERNAL] Final comparison only: the benchmark `K_c={benchmark['benchmark_K_c']}` differs from the untuned series-only central estimate by `{benchmark['absolute_difference']}` and was not used in fitting, filtering, or model selection.
"""
    NOTE_PATH.write_text(note, encoding="utf-8")
    print(f"WROTE {NOTE_PATH.relative_to(ROOT)}", flush=True)


def main() -> bool:
    mp.mp.dps = PRECISION
    timestamp = datetime.now(timezone.utc).isoformat()
    physical_memory = _physical_memory_bytes()
    profiles = [_profile("HT", 22), _profile("HT", 24), _profile("LT", 30), _profile("LT", 32)]
    _print_profiles(profiles, physical_memory)

    frozen_ht = _load_payload(FROZEN_HT_PATH)
    frozen_lt = _load_payload(FROZEN_LT_PATH)
    old_ht = tuple(Fraction(value) for value in frozen_ht["data"]["coefficients"])
    old_lt = tuple(Fraction(value) for value in frozen_lt["data"]["coefficients"])

    lt_started = time.perf_counter()
    lt_rss_before = _rss_bytes()
    with _hybrid_flm_engine() as lt_transfer_records:
        lt_series = low_temperature_free_energy(3, LT_TARGET)
        lt_enlarged = low_temperature_free_energy(3, LT_TARGET, bound_slack=1)
    lt_elapsed = time.perf_counter() - lt_started
    lt_peak_rss = _rss_bytes()

    ht_started = time.perf_counter()
    ht_rss_before = _rss_bytes()
    ht_series = None
    ht_transfer_records = {}
    ht_failure = None
    try:
        with _hybrid_flm_engine() as records:
            ht_series = high_temperature_free_energy(3, HT_TARGET)
        ht_transfer_records = records
    except MemoryError as error:
        ht_failure = f"MemoryError: {error}"
        _clear_series_caches()
    except OSError as error:
        if error.errno != errno.ENOMEM:
            raise
        ht_failure = f"OSError ENOMEM: {error}"
        _clear_series_caches()
    ht_elapsed = time.perf_counter() - ht_started
    ht_peak_rss = _rss_bytes()

    witness = _external_primary_witness()
    lt_log_witness = witness["exact_formal_log_witnesses"]

    lt_checks: list[dict] = []
    _record(
        lt_checks,
        "frozen_LT_prefix_exact",
        lt_series.coefficients[: len(old_lt)] == old_lt,
        f"all {len(old_lt)} coefficients through x^{len(old_lt)-1} agree exactly",
    )
    _record(
        lt_checks,
        "two_new_nonzero_LT_coefficients",
        bool(lt_series.coefficients[30]) and bool(lt_series.coefficients[32]),
        f"x^30={lt_series.coefficients[30]}, x^32={lt_series.coefficients[32]}",
    )
    lt_extra = set(lt_enlarged.boxes) - set(lt_series.boxes)
    _record(
        lt_checks,
        "LT_enlarged_box_self_consistency",
        bool(lt_extra)
        and lt_enlarged.coefficients == lt_series.coefficients
        and all(not any(lt_enlarged.box_weights[shape]) for shape in lt_extra),
        f"{len(lt_extra)} extra ordered boxes at side-sum budget 10 have zero retained weight",
    )
    _record(
        lt_checks,
        "LT_x30_external_primary_witness",
        str(lt_series.coefficients[30]) == lt_log_witness["x^30_equals_log_Lambda0_u^15"],
        "local FLM x^30 compared after derivation with exact log of published Table 2",
    )
    _record(
        lt_checks,
        "LT_x32_external_primary_witness",
        str(lt_series.coefficients[32]) == lt_log_witness["x^32_equals_log_Lambda0_u^16"],
        "local FLM x^32 compared after derivation with exact log of published Table 2",
    )

    lt_data = _series_payload(lt_series)
    lt_data.update(
        {
            "lattice": "simple cubic",
            "normalization": "phi = 3*K + sum_n coefficients[n] * x^n",
            "coefficient_provenance": "derived_locally",
            "extension": {
                "previous_order": len(old_lt) - 1,
                "new_order": LT_TARGET,
                "new_nonzero_coefficients": {
                    "30": str(lt_series.coefficients[30]),
                    "32": str(lt_series.coefficients[32]),
                },
            },
            "profile_completed_before_computation": True,
            "resource_profile": next(item for item in profiles if item["side"] == "LT" and item["order"] == 32),
            "resource_usage": {
                "elapsed_seconds": f"{lt_elapsed:.6f}",
                "rss_before_bytes": lt_rss_before,
                "process_peak_rss_bytes": lt_peak_rss,
                "transfer_record_count": len(lt_transfer_records),
            },
            "external_witness": witness,
            "internal_witness_status": {
                "available": False,
                "method_considered": "periodic transverse slab difference from experiments/e26_lt_series_independent.py",
                "reason": (
                    "validity through x^32 requires transverse circumference L>=9 in both "
                    "directions (4L>32), hence cross-section 81 > CRT guard 22"
                ),
            },
            "self_consistency": {
                "extra_ordered_boxes": len(lt_extra),
                "extra_budget": 10,
                "all_extra_weights_zero": all(not any(lt_enlarged.box_weights[shape]) for shape in lt_extra),
            },
        }
    )

    ht_checks: list[dict] = []
    if ht_series is not None:
        ht_coefficients = ht_series.coefficients
        _record(
            ht_checks,
            "frozen_HT_prefix_exact",
            ht_coefficients[: len(old_ht)] == old_ht,
            f"all {len(old_ht)} coefficients through v^{len(old_ht)-1} agree exactly",
        )
        _record(ht_checks, "HT_v22_nonzero", bool(ht_coefficients[22]), f"v^22={ht_coefficients[22]}")
        _record(
            ht_checks,
            "HT_v22_external_primary_witness",
            str(ht_coefficients[22]) == witness["exact_formal_log_witnesses"]["v^22_total_adds_3_over_22"],
            "local FLM v^22 compared after derivation with exact log of published Table 2",
        )
        ht_data = _series_payload(ht_series)
        ht_data.update(
            {
                "lattice": "simple cubic",
                "normalization": "phi = log(2) + sum_n coefficients[n] * v^n",
                "interaction_normalization": (
                    "phi = log(2) + 3*log(cosh K) + sum_n interaction_coefficients[n] * v^n"
                ),
                "coefficient_provenance": "derived_locally",
                "extension": {
                    "previous_order": len(old_ht) - 1,
                    "new_order": HT_TARGET,
                    "new_nonzero_coefficients": {"22": str(ht_coefficients[22])},
                },
                "external_witness": witness,
            }
        )
    else:
        ht_coefficients = old_ht
        ht_data = _series_data_from_frozen(frozen_ht)
        _record(
            ht_checks,
            "frozen_HT_prefix_retained_after_resource_failure",
            True,
            ht_failure or "HT extension unavailable",
        )
    ht_data.update(
        {
            "profile_completed_before_computation": True,
            "resource_profile_v22": next(item for item in profiles if item["side"] == "HT" and item["order"] == 22),
            "resource_wall_v24": next(item for item in profiles if item["side"] == "HT" and item["order"] == 24),
            "resource_usage": {
                "attempted": True,
                "succeeded": ht_series is not None,
                "failure": ht_failure,
                "elapsed_seconds": f"{ht_elapsed:.6f}",
                "rss_before_bytes": ht_rss_before,
                "process_peak_rss_bytes": ht_peak_rss,
                "profile_estimated_peak_bytes": next(
                    item for item in profiles if item["side"] == "HT" and item["order"] == 22
                )["memory_wall_box"]["estimated_peak_bytes"],
                "physical_memory_bytes": physical_memory,
                "transfer_record_count": len(ht_transfer_records),
                "crt_box_records": sorted(
                    (
                        record
                        for record in ht_transfer_records.values()
                        if record["engine"] == "31-bit-prime CRT"
                    ),
                    key=lambda record: (sum(record["shape"]), record["shape"]),
                ),
            },
        }
    )

    z_series = even_to_squared_variable(ht_coefficients)
    g_series = tuple(z_series[1:])
    dfinite = _full_dfinite_scan(g_series)
    ode_holdout = _ode_holdouts(g_series)
    da_holdout = _differential_holdouts(z_series)
    critical = _critical_estimate(ht_coefficients)

    analysis_checks: list[dict] = []
    _record(
        analysis_checks,
        "all_supported_D_finite_budgets_rejected",
        dfinite["consistent_relation_count"] == 0,
        dfinite["conclusion"],
    )
    _record(
        analysis_checks,
        "strict_ODE_holdouts_predict_remaining",
        bool(ode_holdout) and all(row["held_out_coefficients"] >= 1 for row in ode_holdout),
        f"reported prefix-only recursive predictions for {len(ode_holdout)} values of m",
    )
    _record(
        analysis_checks,
        "strict_differential_holdouts_predict_remaining",
        bool(da_holdout) and all(row["held_out_coefficients"] >= 1 for row in da_holdout),
        f"reported prefix-only recursive predictions for {len(da_holdout)} values of m",
    )
    _record(
        analysis_checks,
        "benchmark_not_used_in_fit_or_selection",
        not critical["benchmark_used_in_fit_or_selection"],
        "critical model family, filters, central value, and uncertainty were fixed from series coefficients only",
    )

    benchmark_kc = mp.mpf("0.221654626")
    central_kc = mp.mpf(critical["K_c"])
    final_comparison = {
        "benchmark_K_c": mp.nstr(benchmark_kc, 12),
        "absolute_difference": mp.nstr(abs(central_kc - benchmark_kc), 35),
        "used_in_fit_or_model_selection": False,
    }
    analysis_data = {
        "input": {
            "artifact": str(HT_PATH.relative_to(ROOT)),
            "known_v_order": len(ht_coefficients) - 1,
            "transformation": "z=v^2; G(z)=(phi_HT-log(2))/z",
            "z_coefficients": _fraction_strings(z_series),
            "G_coefficients": _fraction_strings(g_series),
            "mpmath_decimal_precision": PRECISION,
        },
        "strict_holdout": {
            "rule": (
                "fit only the first m coefficients, recursively predict every remaining "
                "coefficient, and never use holdout error or an external benchmark to select a model"
            ),
            "d_finite_interpolants": ode_holdout,
            "differential_approximants": da_holdout,
        },
        "d_finiteness": dfinite,
        "critical_estimate": critical,
        "final_comparison_only": final_comparison,
    }

    exact_provenance = {
        "script": SCRIPT,
        "timestamp": timestamp,
        "precision": (
            "exact Python int/Fraction; numpy.int64 propagation for N<=62; "
            "deterministic 31-bit-prime CRT with product >2^N beyond 62 sites"
        ),
    }
    analysis_provenance = {
        "script": SCRIPT,
        "timestamp": timestamp,
        "precision": (
            f"exact Fraction/SymPy rational rank and prediction; mpmath {PRECISION} dps for polynomial roots"
        ),
    }
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    _write_json(LT_PATH, exact_provenance, lt_data, lt_checks)
    _write_json(HT_PATH, exact_provenance, ht_data, ht_checks)
    _write_json(ANALYSIS_PATH, analysis_provenance, analysis_data, analysis_checks)
    _write_note(profiles, lt_data, ht_data, analysis_data, witness)

    passed = all(
        check["passed"]
        for checks in (lt_checks, ht_checks, analysis_checks)
        for check in checks
    )
    if passed:
        print(
            f"PASS e32_series_extend2: LT x^{LT_TARGET}; HT v^{len(ht_coefficients)-1}; "
            f"exact external witnesses; strict holdouts; "
            f"D-finite budget <= {dfinite['largest_fully_nonvacuous_parameter_budget']} rejected",
            flush=True,
        )
    else:
        print("FAIL e32_series_extend2: one or more recorded checks failed", flush=True)
    return passed


if __name__ == "__main__":
    try:
        ok = main()
    except Exception as error:
        print(f"FAIL e32_series_extend2: {type(error).__name__}: {error}", flush=True)
        raise
    raise SystemExit(0 if ok else 1)
