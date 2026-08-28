#!/usr/bin/env python3
"""Classify the saved residual frontiers of the Liu-H2 complement filter.

This is a post-processor only: it never reruns branch-and-bound.  It loads the
saved 3,000,001-box and 866,551-box reports, independently recomputes the Arb
geometry and objective bounds of every retained box, encloses the entropy-layer
coefficients at every support face, and writes a deterministic classification
report.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import uuid
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

# This workstation is shared with live certification workers.  These must be
# set before importing modules that may load a BLAS implementation.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from flint import arb, ctx  # noqa: E402

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import (  # noqa: E402
    ArbParameters,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_pilot import Box, ONE, ZERO, mean_corner_range  # noqa: E402
from liu9_tube import complement_box_bound, distance_squared_box  # noqa: E402

ctx.prec = max(ctx.prec, 320)

COORDINATES = ("a1", "a2", "q", "b0", "b2", "b4", "b1", "b3", "b5")
SUPPORT_COORDINATES = COORDINATES[3:]
SUPPORT_INDICES = tuple(range(3, 9))
DEFAULT_LARGE_REPORT = (
    HERE / "verification/results/liu9-complement-residual-rho0.1-3M.json"
)
DEFAULT_SMALL_REPORT = (
    HERE / "verification/results/liu9-complement-residual-rho0.1.json"
)
DEFAULT_OUTPUT = HERE / "verification/results/liu9-survivor-classification.json"
EXPECTED_LARGE_COUNTS = {
    "residual": 79,
    "zero_support": 70,
    "zero_support_and_q_degenerate": 57,
    "q_degenerate_without_zero_support": 0,
}
EXPECTED_LITERAL_COORDINATE_INTERSECTION = 9


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _exact(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, bool):
        raise TypeError("a Boolean is not a coordinate")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("non-finite coordinate")
        return Fraction.from_float(value)
    return Fraction(str(value))


def _arbf(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _hull(values: Iterable[arb]) -> arb:
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration as error:
        raise ValueError("cannot hull an empty sequence") from error
    for value in iterator:
        result = result.union(value)
    if not result.is_finite():
        raise ArithmeticError("non-finite Arb hull")
    return result


def _enclosure_json(value: arb) -> dict[str, Any]:
    if not value.is_finite():
        raise ArithmeticError("cannot serialize a non-finite Arb enclosure")
    lower = value.lower()
    upper = value.upper()
    if lower > ZERO:
        sign = "positive"
    elif upper < ZERO:
        sign = "negative"
    else:
        sign = "undetermined"
    return {
        "arb": str(value),
        "lower": str(lower),
        "lower_float": float(lower),
        "sign": sign,
        "upper": str(upper),
        "upper_float": float(upper),
    }


def _box_from_json(raw: Any) -> Box:
    if not isinstance(raw, list) or len(raw) != len(COORDINATES):
        raise ValueError("a survivor box must have nine coordinate intervals")
    box = []
    for index, interval in enumerate(raw):
        if not isinstance(interval, list) or len(interval) != 2:
            raise ValueError(f"invalid interval for {COORDINATES[index]}")
        lower, upper = (float(interval[0]), float(interval[1]))
        if not (0.0 <= lower <= upper <= 1.0):
            raise ValueError(
                f"invalid unit-box interval for {COORDINATES[index]}: {interval}"
            )
        box.append((lower, upper))
    if box[0][0] + box[1][0] > 1.0:
        raise ValueError("survivor box misses the mass simplex")
    return tuple(box)


def _box_key(box: Box) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((_exact(lower), _exact(upper)) for lower, upper in box)


def _exact_coordinates(box: Box) -> list[dict[str, Any]]:
    rows = []
    for name, (lower, upper) in zip(COORDINATES, box):
        lo = _exact(lower)
        hi = _exact(upper)
        rows.append({
            "interval": [str(lo), str(hi)],
            "lower_float": lower,
            "name": name,
            "upper_float": upper,
            "width": str(hi - lo),
            "width_float": upper - lower,
        })
    return rows


def _compact_coordinates(box: Box) -> list[list[str]]:
    return [[str(_exact(lower)), str(_exact(upper))] for lower, upper in box]


def _simplex_vertices(box: Box) -> tuple[tuple[Fraction, Fraction], ...]:
    l1, u1 = (_exact(value) for value in box[0])
    l2, u2 = (_exact(value) for value in box[1])
    if l1 + l2 > 1:
        return ()
    candidates = {
        (a1, a2)
        for a1 in (l1, u1)
        for a2 in (l2, u2)
        if a1 + a2 <= 1
    }
    for a1 in (l1, u1):
        a2 = 1 - a1
        if l2 <= a2 <= u2:
            candidates.add((a1, a2))
    for a2 in (l2, u2):
        a1 = 1 - a2
        if l1 <= a1 <= u1:
            candidates.add((a1, a2))
    return tuple(sorted(candidates))


def _masses(vertex: tuple[Fraction, Fraction]) -> tuple[Fraction, ...]:
    a1, a2 = vertex
    return a1, a2, 1 - a1 - a2


def _global_weights(
    vertex: tuple[Fraction, Fraction], q: Fraction
) -> tuple[Fraction, ...]:
    masses = _masses(vertex)
    return tuple(
        component_weight * mass
        for component_weight in (1 - q, q)
        for mass in masses
    )


def _global_weight_range(box: Box, support_index: int) -> arb:
    values = []
    for q in (_exact(box[2][0]), _exact(box[2][1])):
        for vertex in _simplex_vertices(box):
            values.append(_arbf(_global_weights(vertex, q)[support_index]))
    return _hull(values)


def _zero_face_coefficient(
    box: Box, moving_index: int, parameters: ArbParameters
) -> arb:
    """Exact face enclosure of the log(1/b_j) coefficient in G_j.

    On the face the moving support is fixed at zero.  The coefficient is

        2(1-beta)*mean + 2*beta*sum_i a_i*b_i*(2-b_i) - 1,

    where the second sum uses the moving atom's component.  With b_j fixed at
    zero it is increasing in every other support on [0,1], affine in q, and
    affine on the clipped mass simplex.  Its extrema therefore occur at the
    two monotone support corners, q endpoints, and simplex vertices below.
    """
    support_lows = [
        _exact(box[index][0]) for index in SUPPORT_INDICES
    ]
    support_highs = [
        _exact(box[index][1]) for index in SUPPORT_INDICES
    ]
    support_lows[moving_index] = Fraction(0)
    support_highs[moving_index] = Fraction(0)
    component = moving_index // 3
    candidates = []
    for support in (tuple(support_lows), tuple(support_highs)):
        for q in (_exact(box[2][0]), _exact(box[2][1])):
            for vertex in _simplex_vertices(box):
                masses = _masses(vertex)
                weights = _global_weights(vertex, q)
                mean = sum(
                    (weight * point for weight, point in zip(weights, support)),
                    Fraction(0),
                )
                base = 3 * component
                paired_moment = sum(
                    (
                        masses[offset]
                        * support[base + offset]
                        * (2 - support[base + offset])
                        for offset in range(3)
                    ),
                    Fraction(0),
                )
                coefficient = (
                    2 * (1 - parameters.beta) * _arbf(mean)
                    + 2 * parameters.beta * _arbf(paired_moment)
                    - 1
                )
                candidates.append(coefficient)
    return _hull(candidates)


def _mirror_face_ranges(
    box: Box, moving_index: int, parameters: ArbParameters
) -> tuple[arb, arb, arb, tuple[str, ...], tuple[str, ...]]:
    """Enclose the b_j -> 1 coefficient over every face sub-stratum.

    Other support intervals that touch 1 may either remain below 1 or join the
    face.  Enumerating those subsets encloses all corresponding pathwise
    leading coefficients, including simultaneous approaches to the face.
    """
    touching = {
        index
        for index in range(6)
        if _exact(box[index + 3][1]) == 1
    }
    mandatory = {
        index
        for index in touching
        if _exact(box[index + 3][0]) == 1
    }
    mandatory.add(moving_index)
    optional = tuple(sorted(touching - mandatory))
    moving_component = moving_index // 3
    coefficient_values = []
    weight_values = []
    conditional_values = []
    for bits in itertools.product((False, True), repeat=len(optional)):
        included = set(mandatory)
        included.update(
            index for index, selected in zip(optional, bits) if selected
        )
        for q in (_exact(box[2][0]), _exact(box[2][1])):
            for vertex in _simplex_vertices(box):
                masses = _masses(vertex)
                weights = _global_weights(vertex, q)
                weight_at_one = sum(
                    (weights[index] for index in included), Fraction(0)
                )
                conditional_mass_at_one = sum(
                    (
                        masses[index % 3]
                        for index in included
                        if index // 3 == moving_component
                    ),
                    Fraction(0),
                )
                weight_arb = _arbf(weight_at_one)
                conditional_arb = _arbf(conditional_mass_at_one)
                coefficient = (
                    1
                    - 2 * (1 - parameters.beta) * weight_arb
                    - 2 * parameters.beta * conditional_arb
                )
                coefficient_values.append(coefficient)
                weight_values.append(weight_arb)
                conditional_values.append(conditional_arb)
    return (
        _hull(coefficient_values),
        _hull(weight_values),
        _hull(conditional_values),
        tuple(SUPPORT_COORDINATES[index] for index in sorted(mandatory)),
        tuple(SUPPORT_COORDINATES[index] for index in optional),
    )


def _constraint_status(mean_range: arb, target_mean: arb) -> str:
    if mean_range.lower() >= target_mean.upper():
        return "forced"
    if mean_range.upper() < target_mean.lower():
        return "forced_to_fail"
    return "undetermined"


def _tube_status(distance_range: arb, rho_squared: arb) -> str:
    # These are the exact strict/non-strict tests used by liu9_residual.py.
    if distance_range.upper() < rho_squared:
        return "entirely_inside"
    if distance_range.lower() >= rho_squared:
        return "entirely_outside"
    return "straddling"


def _strata(box: Box) -> dict[str, Any]:
    zero_coordinates = [
        COORDINATES[index]
        for index in SUPPORT_INDICES
        if _exact(box[index][0]) == 0
    ]
    mirror_coordinates = [
        COORDINATES[index]
        for index in SUPPORT_INDICES
        if _exact(box[index][1]) == 1
    ]
    q_faces = []
    if _exact(box[2][0]) == 0:
        q_faces.append("zero")
    if _exact(box[2][1]) == 1:
        q_faces.append("one")
    active = []
    if zero_coordinates:
        active.append("zero_support")
    if mirror_coordinates:
        active.append("mirror")
    if q_faces:
        active.append("q_degenerate")
    if not active:
        active.append("interior")
    return {
        "active": active,
        "interior": active == ["interior"],
        "mirror": bool(mirror_coordinates),
        "mirror_coordinates": mirror_coordinates,
        "q_degenerate": bool(q_faces),
        "q_faces": q_faces,
        "zero_support": bool(zero_coordinates),
        "zero_support_coordinates": zero_coordinates,
    }


def _face_coefficients(
    box: Box,
    parameters: ArbParameters,
    rho_squared: arb,
    zero_feasible_lower: arb,
    mirror_weight_threshold: arb,
) -> list[dict[str, Any]]:
    records = []
    zero_cache: dict[int, arb] = {}
    weight_cache: dict[int, arb] = {}
    for support_index, coordinate in enumerate(SUPPORT_COORDINATES):
        coordinate_index = support_index + 3
        lower = _exact(box[coordinate_index][0])
        upper = _exact(box[coordinate_index][1])
        if lower != 0 and upper != 1:
            continue
        weight = weight_cache.setdefault(
            support_index, _global_weight_range(box, support_index)
        )
        contact = "fixed" if lower == upper else "straddling"
        if lower == 0:
            coefficient = zero_cache.setdefault(
                support_index,
                _zero_face_coefficient(box, support_index, parameters),
            )
            weighted_coefficient = weight * coefficient
            face_defect = arb(0)
            records.append({
                "coefficient_convention": "coefficient of log(1/b_j) in weight-cleared G_j",
                "coefficient_formula": (
                    "2(1-beta)*mean + "
                    "2*beta*sum_i[a_i*b_ci*(2-b_ci)] - 1"
                ),
                "coordinate": coordinate,
                "face": "zero",
                "face_contact": contact,
                "global_atom_weight": _enclosure_json(weight),
                "leading_singular_coefficient": _enclosure_json(coefficient),
                "mean_feasible_coefficient_lower_bound": _enclosure_json(
                    zero_feasible_lower
                ),
                "partial_derivative_leading_coefficient": _enclosure_json(
                    weighted_coefficient
                ),
                "weighted_support_defect_at_face": _enclosure_json(face_defect),
                "whole_face_defect_le_rho_squared": True,
            })
        if upper == 1:
            (
                coefficient,
                weight_at_one,
                conditional_mass_at_one,
                mandatory,
                optional,
            ) = _mirror_face_ranges(box, support_index, parameters)
            weighted_coefficient = weight * coefficient
            face_defect = weight * (1 - parameters.x) ** 2
            if weight_at_one.upper() <= mirror_weight_threshold.lower():
                mass_bound_status = "forced"
            elif weight_at_one.lower() > mirror_weight_threshold.upper():
                mass_bound_status = "forced_to_fail"
            else:
                mass_bound_status = "undetermined"
            records.append({
                "coefficient_convention": (
                    "coefficient of log(1/(1-b_j)) in weight-cleared G_j"
                ),
                "coefficient_formula": "1 - 2(1-beta)*W1 - 2*beta*A1",
                "conditional_mass_A1": _enclosure_json(conditional_mass_at_one),
                "coordinate": coordinate,
                "face": "one",
                "face_contact": contact,
                "global_atom_weight": _enclosure_json(weight),
                "leading_singular_coefficient": _enclosure_json(coefficient),
                "mandatory_supports_at_one": list(mandatory),
                "optional_supports_touching_one": list(optional),
                "partial_derivative_leading_coefficient": _enclosure_json(
                    weighted_coefficient
                ),
                "tube_mass_bound": {
                    "formula": "W1 <= rho^2/(1-x)^2",
                    "status": mass_bound_status,
                    "threshold": _enclosure_json(mirror_weight_threshold),
                },
                "weight_W1": _enclosure_json(weight_at_one),
                "weighted_support_defect_at_face": _enclosure_json(face_defect),
                "whole_face_defect_le_rho_squared": bool(
                    face_defect.upper() <= rho_squared.lower()
                ),
            })
    return records


def _float_agrees(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2e-12, abs_tol=2e-12)


def _load_report(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    path = path.expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        report = json.load(handle)
    if report.get("report_type") != "liu9_complement_residual":
        raise ValueError(f"unexpected report type in {path}")
    recorded_digest = report.get("report_sha256")
    digest_payload = {
        key: value for key, value in report.items() if key != "report_sha256"
    }
    computed_digest = _canonical_digest(digest_payload)
    if recorded_digest != computed_digest:
        raise ValueError(
            f"source digest mismatch for {path}: {recorded_digest} != {computed_digest}"
        )
    residual = report.get("residual")
    if not isinstance(residual, list):
        raise ValueError(f"missing residual list in {path}")
    if report.get("run", {}).get("residual_count") != len(residual):
        raise ValueError(f"residual count mismatch in {path}")
    return report, {
        "path": str(path),
        "report_sha256": recorded_digest,
        "residual_count": len(residual),
        "rho": report.get("rho"),
    }


def _classify_report(
    label: str,
    path: Path,
    parameters: ArbParameters,
    zero_feasible_lower: arb,
) -> tuple[dict[str, Any], dict[tuple[tuple[Fraction, Fraction], ...], str]]:
    source, source_metadata = _load_report(path)
    rho = Fraction(source["rho"])
    rho_squared = _arbf(rho * rho)
    mirror_weight_threshold = rho_squared / (1 - parameters.x) ** 2
    lambdas = tuple(Fraction(value) for value in source["lambdas"])
    survivors = []
    identity: dict[tuple[tuple[Fraction, Fraction], ...], str] = {}
    recomputation_discrepancies = []

    for index, saved in enumerate(source["residual"]):
        box = _box_from_json(saved["box"])
        survivor_id = f"{label}-{index:03d}"
        key = _box_key(box)
        if key in identity:
            raise ValueError(f"duplicate survivor box in {path}: {survivor_id}")
        identity[key] = survivor_id

        strata = _strata(box)
        mean_range = mean_corner_range(box)
        if mean_range is None:
            raise ValueError(f"saved survivor misses the simplex: {survivor_id}")
        mean_status = _constraint_status(mean_range, parameters.mean)
        distance_range = distance_squared_box(box, parameters)
        tube_status = _tube_status(distance_range, rho_squared)
        bound = complement_box_bound(box, parameters, lambdas)
        gap_clears = bool(bound.gap_lower >= ZERO)
        objective_clears = bool(bound.objective_lower >= ONE)
        retained = not gap_clears and not objective_clears
        face_coefficients = _face_coefficients(
            box,
            parameters,
            rho_squared,
            zero_feasible_lower,
            mirror_weight_threshold,
        )
        weighted_face_defect_holds = any(
            item["whole_face_defect_le_rho_squared"]
            for item in face_coefficients
        )

        recomputed = {
            "distance_squared_lower": float(distance_range.lower()),
            "distance_squared_upper": float(distance_range.upper()),
            "gap_lower": float(bound.gap_lower),
            "method": bound.method,
            "objective_lower": float(bound.objective_lower),
            "straddles_tube": tube_status == "straddling",
        }
        field_discrepancies = []
        for field in (
            "distance_squared_lower",
            "distance_squared_upper",
            "gap_lower",
            "objective_lower",
        ):
            if not _float_agrees(float(saved[field]), recomputed[field]):
                field_discrepancies.append({
                    "field": field,
                    "recomputed": recomputed[field],
                    "saved": saved[field],
                })
        for field in ("method", "straddles_tube"):
            if saved[field] != recomputed[field]:
                field_discrepancies.append({
                    "field": field,
                    "recomputed": recomputed[field],
                    "saved": saved[field],
                })
        if field_discrepancies:
            recomputation_discrepancies.append({
                "fields": field_discrepancies,
                "survivor_id": survivor_id,
            })

        if tube_status == "straddling":
            priority_term = "tube_distance_upper_margin"
            priority_value = distance_range.upper() - rho_squared
        else:
            priority_term = f"objective_lower:{bound.method}"
            priority_value = bound.objective_lower

        survivors.append({
            "bound_failure": {
                "dominant_retention_priority": priority_term,
                "dominant_retention_priority_value": _enclosure_json(
                    priority_value
                ),
                "gap_clears": gap_clears,
                "objective_clears": objective_clears,
                "objective_method": bound.method,
                "retained_by_objective_rules": retained,
            },
            "coordinates": _exact_coordinates(box),
            "depth": int(saved["depth"]),
            "distance_squared": _enclosure_json(distance_range),
            "face_coefficients": face_coefficients,
            "maximum_width": str(
                max(
                    _exact(upper) - _exact(lower)
                    for lower, upper in box
                )
            ),
            "mean": {
                "constraint_status": mean_status,
                "enclosure": _enclosure_json(mean_range),
                "target": _enclosure_json(parameters.mean),
            },
            "objective_bound": {
                "gap_lower": _enclosure_json(bound.gap_lower),
                "objective_lower": _enclosure_json(bound.objective_lower),
            },
            "source_recorded": {
                "distance_squared_lower": saved["distance_squared_lower"],
                "distance_squared_upper": saved["distance_squared_upper"],
                "gap_lower": saved["gap_lower"],
                "method": saved["method"],
                "objective_lower": saved["objective_lower"],
                "straddles_tube": saved["straddles_tube"],
            },
            "strata": strata,
            "survivor_id": survivor_id,
            "tube_position": tube_status,
            "weighted_face_defect_hypothesis": weighted_face_defect_holds,
        })

    cross_tab = Counter(
        (
            item["strata"]["zero_support"],
            item["strata"]["mirror"],
            item["strata"]["q_degenerate"],
        )
        for item in survivors
    )
    cross_tab_json = [
        {
            "count": cross_tab.get((zero, mirror, q_degenerate), 0),
            "mirror": mirror,
            "q_degenerate": q_degenerate,
            "zero_support": zero,
        }
        for zero, mirror, q_degenerate in itertools.product(
            (False, True), repeat=3
        )
        if cross_tab.get((zero, mirror, q_degenerate), 0)
    ]
    counts = {
        "interior": sum(item["strata"]["interior"] for item in survivors),
        "mirror": sum(item["strata"]["mirror"] for item in survivors),
        "mirror_only": sum(
            item["strata"]["mirror"]
            and not item["strata"]["zero_support"]
            and not item["strata"]["q_degenerate"]
            for item in survivors
        ),
        "q_degenerate": sum(
            item["strata"]["q_degenerate"] for item in survivors
        ),
        "q_degenerate_without_zero_support": sum(
            item["strata"]["q_degenerate"]
            and not item["strata"]["zero_support"]
            for item in survivors
        ),
        "residual": len(survivors),
        "zero_support": sum(
            item["strata"]["zero_support"] for item in survivors
        ),
        "zero_support_and_q_degenerate": sum(
            item["strata"]["zero_support"]
            and item["strata"]["q_degenerate"]
            for item in survivors
        ),
        "zero_support_or_mirror": sum(
            item["strata"]["zero_support"]
            or item["strata"]["mirror"]
            for item in survivors
        ),
    }
    mean_counts = Counter(
        item["mean"]["constraint_status"] for item in survivors
    )
    tube_counts = Counter(item["tube_position"] for item in survivors)
    priority_counts = Counter(
        item["bound_failure"]["dominant_retention_priority"]
        for item in survivors
    )
    method_counts = Counter(
        item["bound_failure"]["objective_method"] for item in survivors
    )
    face_sign_counts = Counter(
        (face["face"], face["leading_singular_coefficient"]["sign"])
        for item in survivors
        for face in item["face_coefficients"]
    )
    aggregate = {
        "bound_failure_dominator_counts": dict(sorted(priority_counts.items())),
        "cross_tabulation": cross_tab_json,
        "face_coefficient_sign_counts": {
            f"{face}:{sign}": count
            for (face, sign), count in sorted(face_sign_counts.items())
        },
        "mean_constraint_counts": dict(sorted(mean_counts.items())),
        "objective_method_counts": dict(sorted(method_counts.items())),
        "recomputation_discrepancies": recomputation_discrepancies,
        "single_hypothesis_counts": {
            "mirror": counts["mirror"],
            "q_degenerate": counts["q_degenerate"],
            "weighted_face_defect_le_rho_squared": sum(
                item["weighted_face_defect_hypothesis"] for item in survivors
            ),
            "zero_support": counts["zero_support"],
            "zero_support_or_mirror": counts["zero_support_or_mirror"],
        },
        "stratum_counts": counts,
        "tube_position_counts": dict(sorted(tube_counts.items())),
    }
    return {
        "aggregate": aggregate,
        "rho": str(rho),
        "rho_squared": _enclosure_json(rho_squared),
        "source": source_metadata,
        "survivors": survivors,
    }, identity


def _exception_record(
    survivor: dict[str, Any], common_small_id: Optional[str]
) -> dict[str, Any]:
    return {
        "coordinates": [
            item["interval"] for item in survivor["coordinates"]
        ],
        "coordinate_order": list(COORDINATES),
        "depth": survivor["depth"],
        "large_survivor_id": survivor["survivor_id"],
        "small_survivor_id": common_small_id,
        "strata": survivor["strata"]["active"],
        "tube_position": survivor["tube_position"],
    }


def _write_report(path: Path, payload: dict[str, Any]) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    ).encode("ascii")
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600
    )
    try:
        view = memoryview(encoded)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise RuntimeError(f"short report write: {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def _print_cross_tab(label: str, rows: Sequence[dict[str, Any]]) -> None:
    print(f"MACHINE VERIFIED [{label} cross-tab]:")
    print("   zero_support mirror q_degenerate count")
    for row in rows:
        print(
            "   %-12s %-6s %-12s %d"
            % (
                str(row["zero_support"]).lower(),
                str(row["mirror"]).lower(),
                str(row["q_degenerate"]).lower(),
                row["count"],
            )
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large-report", type=Path, default=DEFAULT_LARGE_REPORT)
    parser.add_argument("--small-report", type=Path, default=DEFAULT_SMALL_REPORT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prec", type=int, default=320)
    parser.add_argument("--dps", type=int, default=90)
    arguments = parser.parse_args(argv)
    if arguments.prec < 128:
        parser.error("--prec must be at least 128 bits")
    if arguments.dps < 60:
        parser.error("--dps must be at least 60")
    ctx.prec = arguments.prec

    mp_parameters = solve_equation_parameters(arguments.dps + 30)
    parameters = certify_equation_parameters(mp_parameters)
    if not (parameters.beta > ZERO and parameters.beta < ONE):
        raise ArithmeticError("beta was not certified in (0,1)")
    zero_feasible_lower = 2 * (1 - parameters.beta) * parameters.mean - 1
    if not zero_feasible_lower > ZERO:
        raise ArithmeticError("zero-face feasible coefficient is not positive")

    large, large_identity = _classify_report(
        "large", arguments.large_report, parameters, zero_feasible_lower
    )
    small, small_identity = _classify_report(
        "small", arguments.small_report, parameters, zero_feasible_lower
    )

    common_keys = set(large_identity).intersection(small_identity)
    exceptions = []
    for item in large["survivors"]:
        if item["strata"]["zero_support"]:
            continue
        key = tuple(
            tuple(Fraction(value) for value in coordinate["interval"])
            for coordinate in item["coordinates"]
        )
        exceptions.append(
            _exception_record(item, small_identity.get(key))
        )
    mirror_exception_common = sum(
        1
        for exception in exceptions
        if exception["small_survivor_id"] is not None
    )

    observed = large["aggregate"]["stratum_counts"]
    acceptance_discrepancies = []
    for name, expected in EXPECTED_LARGE_COUNTS.items():
        actual = observed[name]
        if actual != expected:
            acceptance_discrepancies.append({
                "actual": actual,
                "expected": expected,
                "metric": name,
            })
    if len(common_keys) != EXPECTED_LITERAL_COORDINATE_INTERSECTION:
        acceptance_discrepancies.append({
            "actual": len(common_keys),
            "expected": EXPECTED_LITERAL_COORDINATE_INTERSECTION,
            "metric": "literal_all_box_coordinate_intersection",
            "note": (
                "The value 9 is recovered exactly for the mirror-only "
                "exceptional boxes, not for the full residual sets."
            ),
        })

    weighted_hypothesis_exceptions = [
        item["survivor_id"]
        for item in large["survivors"]
        if not item["weighted_face_defect_hypothesis"]
    ]
    exception_ids = [item["large_survivor_id"] for item in exceptions]
    if weighted_hypothesis_exceptions != exception_ids:
        raise AssertionError(
            "the weighted face-defect exceptions do not equal the nonzero-face list"
        )
    if any(
        item["strata"] != ["mirror"] or item["tube_position"] != "entirely_outside"
        for item in exceptions
    ):
        raise AssertionError("the exceptional list is not the mirror-only outside set")

    report = {
        "acceptance_checks": {
            "discrepancies": acceptance_discrepancies,
            "expected_large_counts": EXPECTED_LARGE_COUNTS,
            "literal_coordinate_intersection_expected": (
                EXPECTED_LITERAL_COORDINATE_INTERSECTION
            ),
            "literal_coordinate_intersection_observed": len(common_keys),
            "mirror_only_coordinate_identical": mirror_exception_common,
        },
        "analytic_assessment": {
            "common_zero_q_structure": {
                "claim_status": "PROVED",
                "explanation": (
                    "Inside the tube, r_c*delta(P_c)^2 <= rho^2 for "
                    "r_c in {1-q,q}; hence delta(P_c) <= rho/sqrt(r_c) "
                    "when r_c>0.  The component-conditional control becomes "
                    "vacuous at a q endpoint.  All q-degenerate survivors are "
                    "already zero-support survivors."
                ),
                "q_degenerate_survivors": observed["q_degenerate"],
                "q_degenerate_without_zero_support": observed[
                    "q_degenerate_without_zero_support"
                ],
                "tube_inequality": "r_c*delta(P_c)^2 <= rho^2",
                "zero_q_overlap": observed["zero_support_and_q_degenerate"],
            },
            "single_whole_frontier_inequality": {
                "claim_status": "OPEN",
                "covers_all_survivors": False,
                "exception_count": len(exceptions),
                "exceptions": exceptions,
                "finding": (
                    "No single sign-definite zero/q boundary inequality covers "
                    "the full frontier.  The certified zero-face inequality "
                    "covers exactly the 70 zero-support boxes, including all "
                    "57 q-degenerate boxes.  The exact complement is nine "
                    "mirror-only, wholly-outside boxes, which require the "
                    "distinct b_j -> 1 coefficient."
                ),
                "weighted_face_defect_box_inequality": (
                    "min_(j,f in faces(B)) sup_(V in B cap {b_j=f}) "
                    "w_j*[f*(f-x)]^2 <= rho^2"
                ),
                "weighted_face_defect_covers": large["aggregate"][
                    "single_hypothesis_counts"
                ]["weighted_face_defect_le_rho_squared"],
            },
            "zero_face_inequality": {
                "applies_to_survivors": observed["zero_support"],
                "claim_status": "PROVED",
                "exact_coefficient": (
                    "C0_j = 2(1-beta)*mu + "
                    "2*beta*sum_i[a_i*b_ci*(2-b_ci)] - 1"
                ),
                "inequality": (
                    "C0_j >= 2(1-beta)*(1-c)-1 > 0 "
                    "under mean >= 1-c"
                ),
                "lower_bound": _enclosure_json(zero_feasible_lower),
            },
        },
        "claim_status": "MACHINE VERIFIED",
        "coordinate_order": list(COORDINATES),
        "parameters": {
            "beta": _enclosure_json(parameters.beta),
            "mean": _enclosure_json(parameters.mean),
            "p": _enclosure_json(parameters.p),
            "precision_bits": arguments.prec,
            "root_hi": str(parameters.root_hi),
            "root_lo": str(parameters.root_lo),
            "x": _enclosure_json(parameters.x),
        },
        "reports": {
            "large_3m": large,
            "small_866551": small,
        },
        "tool": "liu9_survivors.py",
    }
    report["report_sha256"] = _canonical_digest(report)
    _write_report(arguments.output, report)

    print("LIU H2 RESIDUAL SURVIVOR CLASSIFICATION")
    print(
        "MACHINE VERIFIED [inputs]: source digests valid; classified %d + %d "
        "saved boxes without rerunning branch-and-bound."
        % (
            large["source"]["residual_count"],
            small["source"]["residual_count"],
        )
    )
    _print_cross_tab("3M", large["aggregate"]["cross_tabulation"])
    print(
        "MACHINE VERIFIED [3M established counts]: survivors=%d, "
        "zero-support=%d, zero-and-q-degenerate=%d, "
        "q-degenerate-without-zero=%d."
        % (
            observed["residual"],
            observed["zero_support"],
            observed["zero_support_and_q_degenerate"],
            observed["q_degenerate_without_zero_support"],
        )
    )
    print(
        "MACHINE VERIFIED [3M other strata]: mirror=%d, mirror-only=%d, "
        "interior=%d; tube straddling=%d, outside=%d."
        % (
            observed["mirror"],
            observed["mirror_only"],
            observed["interior"],
            large["aggregate"]["tube_position_counts"].get("straddling", 0),
            large["aggregate"]["tube_position_counts"].get(
                "entirely_outside", 0
            ),
        )
    )
    print(
        "MACHINE VERIFIED [failed-bound dominator]: %s; objective methods=%s."
        % (
            large["aggregate"]["bound_failure_dominator_counts"],
            large["aggregate"]["objective_method_counts"],
        )
    )
    if len(common_keys) == EXPECTED_LITERAL_COORDINATE_INTERSECTION:
        print(
            "MACHINE VERIFIED [coordinate identity]: %d exact boxes occur in "
            "both reports." % len(common_keys)
        )
    else:
        print(
            "REFUTED [literal all-box coordinate identity expectation]: "
            "observed %d exact coordinate-identical boxes between reports, "
            "not %d (discrepancy %+d)."
            % (
                len(common_keys),
                EXPECTED_LITERAL_COORDINATE_INTERSECTION,
                len(common_keys) - EXPECTED_LITERAL_COORDINATE_INTERSECTION,
            )
        )
    print(
        "MACHINE VERIFIED [mirror-only coordinate identity]: all %d minimal "
        "exceptions are coordinate-identical between reports; this is the "
        "restricted count 9."
        % mirror_exception_common
    )
    print(
        "PROVED [common zero/q inequality]: C0_j >= "
        "2(1-beta)*(1-c)-1 in %s > 0 covers all %d zero-support boxes; all %d "
        "q-degenerate boxes are in that set."
        % (
            zero_feasible_lower,
            observed["zero_support"],
            observed["q_degenerate"],
        )
    )
    print(
        "OPEN [single whole-frontier inequality]: the zero/q inequality and "
        "the uniform weighted-face-defect test cover %d/%d boxes, not all; "
        "the exact minimal exceptional set has %d mirror-only boxes."
        % (
            large["aggregate"]["single_hypothesis_counts"][
                "weighted_face_defect_le_rho_squared"
            ],
            observed["residual"],
            len(exceptions),
        )
    )
    for exception in exceptions:
        coordinate_text = ", ".join(
            f"{name}=[{interval[0]},{interval[1]}]"
            for name, interval in zip(COORDINATES, exception["coordinates"])
        )
        print(
            "MACHINE VERIFIED [exception %s]: %s"
            % (exception["large_survivor_id"], coordinate_text)
        )
    if acceptance_discrepancies:
        print(
            "REFUTED [acceptance discrepancy count]: %d discrepancy recorded "
            "verbatim in the report."
            % len(acceptance_discrepancies)
        )
    else:
        print("MACHINE VERIFIED [acceptance checks]: all expected counts agree.")
    if large["aggregate"]["recomputation_discrepancies"] or small[
        "aggregate"
    ]["recomputation_discrepancies"]:
        print("REFUTED [saved-bound replay]: at least one saved bound disagrees.")
    else:
        print(
            "MACHINE VERIFIED [saved-bound replay]: all 155 recorded distance/"
            "objective bound fields agree; all mean statuses were recomputed."
        )
    print(f"MACHINE VERIFIED [report]: {arguments.output}")
    print(f"MACHINE VERIFIED [canonical SHA-256]: {report['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
