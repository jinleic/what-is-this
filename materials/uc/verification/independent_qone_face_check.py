#!/usr/bin/env python3
"""Independent q=1 active-face sampling check for Liu Hypothesis 2.

The source report is parsed and authenticated directly.  Its q-at-one boxes are
projected to (a1, a2, b1, b3, b5), the only coordinates with nonzero weight at
q=1.  A deterministic, feasibility-gated coordinate descent then starts from
an actual feasible corner of every projected box.  Raw gap evaluations use the
shared high-precision ``gap_mp`` transcription; the objective pencil is never
used as a substitute for the raw gap.

This is a numerical cross-check, not an exhaustive proof of the q=1 face.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Sequence

# Keep this checker harmless alongside the live certification campaigns.  These
# must be set before importing modules that can transitively load SciPy/BLAS.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

import mpmath  # noqa: E402

from liu9_binding import MPParameters, solve_equation_parameters  # noqa: E402
from liu9_boundary_layer import gap_mp  # noqa: E402


SOURCE_REPORT = (
    HERE / "results" / "liu9-complement-residual-rho1-1728-100k.json"
)
QONE_REPORT = HERE / "results" / "liu9-qone-face-rho1-1728-corrected-100k.json"
ACTIVE_INDICES = (0, 1, 6, 7, 8)
ACTIVE_NAMES = ("a1", "a2", "b1", "b3", "b5")
EXPECTED_Q_AT_ONE = 48_996
EXPECTED_UNIQUE_ACTIVE = 230
DPS = 90
DESCENT_LEVELS = 4
ZERO_TOLERANCE_TEXT = "1e-70"
NEGATIVE_TOLERANCE_TEXT = "1e-60"

ActiveBox = tuple[tuple[float, float], ...]
MPBox = tuple[tuple[mpmath.mpf, mpmath.mpf], ...]
Point = tuple[mpmath.mpf, ...]


class VerificationError(RuntimeError):
    """A required independent-check invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def load_digest_checked(path: Path) -> tuple[dict[str, Any], str]:
    """Load a report and replay its canonical object digest."""
    require(path.is_file(), f"missing report: {path}")
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    require(isinstance(payload, dict), f"report root is not an object: {path}")
    recorded = payload.get("report_sha256")
    require(
        isinstance(recorded, str) and len(recorded) == 64,
        f"report has no canonical SHA-256: {path}",
    )
    unsigned = dict(payload)
    unsigned.pop("report_sha256", None)
    computed = canonical_digest(unsigned)
    require(
        recorded == computed,
        f"canonical digest mismatch for {path}: {recorded} != {computed}",
    )
    return payload, computed


def parse_box(row: Any, row_number: int) -> tuple[tuple[float, float], ...]:
    require(isinstance(row, dict), f"residual row {row_number} is not an object")
    raw_box = row.get("box")
    require(
        isinstance(raw_box, list) and len(raw_box) == 9,
        f"residual row {row_number} does not contain a nine-coordinate box",
    )
    parsed: list[tuple[float, float]] = []
    for coordinate, pair in enumerate(raw_box):
        require(
            isinstance(pair, list) and len(pair) == 2,
            f"row {row_number}, coordinate {coordinate} is not an interval",
        )
        require(
            all(
                isinstance(value, (int, float)) and not isinstance(value, bool)
                for value in pair
            ),
            f"row {row_number}, coordinate {coordinate} is not numeric",
        )
        lo, hi = float(pair[0]), float(pair[1])
        require(
            math.isfinite(lo) and math.isfinite(hi) and 0.0 <= lo <= hi <= 1.0,
            f"row {row_number}, coordinate {coordinate} has invalid bounds",
        )
        parsed.append((lo, hi))
    return tuple(parsed)


@dataclass(frozen=True)
class Projection:
    boxes: tuple[ActiveBox, ...]
    q_at_one_count: int
    residual_count: int
    q_intervals: dict[str, int]


def project_q_at_one(source: dict[str, Any]) -> Projection:
    """Quotient q and the three zero-weight P0 support coordinates at q=1."""
    require(
        source.get("report_type") == "liu9_complement_residual",
        "source report type is not liu9_complement_residual",
    )
    residual = source.get("residual")
    require(isinstance(residual, list), "source residual is not a list")

    active_keys: set[ActiveBox] = set()
    q_at_one_count = 0
    q_intervals: dict[str, int] = {}
    for row_number, row in enumerate(residual):
        box = parse_box(row, row_number)
        qlo, qhi = box[2]
        if "q_low" in row:
            require(float(row["q_low"]) == qlo, f"row {row_number} q_low disagrees")
        if "q_high" in row:
            require(float(row["q_high"]) == qhi, f"row {row_number} q_high disagrees")
        if qhi != 1.0 or qlo <= 0.0:
            continue
        q_at_one_count += 1
        interval_label = str([qlo, qhi])
        q_intervals[interval_label] = q_intervals.get(interval_label, 0) + 1
        active_keys.add(tuple(box[index] for index in ACTIVE_INDICES))

    run = source.get("run")
    summary = source.get("summary")
    require(isinstance(run, dict), "source run metadata is missing")
    require(isinstance(summary, dict), "source summary metadata is missing")
    require(
        run.get("residual_count") == len(residual)
        and summary.get("residual_count") == len(residual),
        "source residual counts disagree with the residual array",
    )
    require(
        summary.get("q_at_one_only") == q_at_one_count,
        "source q-at-one summary disagrees with direct box classification",
    )
    return Projection(
        boxes=tuple(sorted(active_keys)),
        q_at_one_count=q_at_one_count,
        residual_count=len(residual),
        q_intervals=dict(sorted(q_intervals.items())),
    )


def cross_check_qone_report(
    path: Path,
    source: dict[str, Any],
    source_digest: str,
    projection: Projection,
) -> str | None:
    """Replay the optional q-one report digest and compare its source envelope."""
    if not path.is_file():
        return None
    report, report_digest = load_digest_checked(path)
    metadata = report.get("source")
    require(isinstance(metadata, dict), "q-one report source metadata is missing")
    require(
        metadata.get("q_at_one_source_boxes") == projection.q_at_one_count,
        "q-one report q-at-one source count disagrees",
    )
    require(
        metadata.get("unique_active_boxes") == len(projection.boxes),
        "q-one report unique active-box count disagrees",
    )
    require(
        metadata.get("q_intervals") == projection.q_intervals,
        "q-one report q-interval counts disagree",
    )
    expected_duplication = projection.q_at_one_count / len(projection.boxes)
    require(
        metadata.get("duplication_factor") == expected_duplication,
        "q-one report duplication factor disagrees",
    )

    nested = metadata.get("source")
    require(isinstance(nested, dict), "q-one report nested source envelope is missing")
    require(
        nested.get("report_sha256") == source_digest,
        "q-one report names a different source digest",
    )
    require(
        nested.get("residual_count") == projection.residual_count,
        "q-one report names a different source residual count",
    )
    require(
        nested.get("rho") == source.get("rho"),
        "q-one report names a different source rho",
    )
    return report_digest


def to_mp_box(box: ActiveBox) -> MPBox:
    return tuple(
        (mpmath.mpf(repr(lo)), mpmath.mpf(repr(hi))) for lo, hi in box
    )


def qone_mean(point: Sequence[mpmath.mpf]) -> mpmath.mpf:
    a1, a2, b1, b3, b5 = point
    return a1 * b1 + a2 * b3 + (1 - a1 - a2) * b5


def classify_qone(
    point: Sequence[mpmath.mpf],
    box: MPBox,
    target_mean: mpmath.mpf,
) -> tuple[bool, str, mpmath.mpf]:
    """Classify before any raw-gap call; no feasibility tolerance is used."""
    if len(point) != 5:
        return False, "shape", mpmath.nan
    for value, (lo, hi) in zip(point, box):
        if not mpmath.isfinite(value) or value < lo or value > hi:
            return False, "bounds", mpmath.nan
    a1, a2 = point[0], point[1]
    if a1 < 0 or a2 < 0 or a1 + a2 > 1:
        return False, "simplex", mpmath.nan
    mean = qone_mean(point)
    if mean < target_mean:
        return False, "mean", mean
    return True, "feasible", mean


def first_feasible_corner(box: MPBox, target_mean: mpmath.mpf) -> tuple[Point, mpmath.mpf]:
    """Return the lexicographically first exact corner satisfying every constraint."""
    for bits in product((0, 1), repeat=5):
        point = tuple(box[index][bit] for index, bit in enumerate(bits))
        feasible, _, mean = classify_qone(point, box, target_mean)
        if feasible:
            return point, mean - target_mean
    raise VerificationError("an active box has no genuinely feasible q=1 corner")


def qone_values(
    point: Sequence[mpmath.mpf],
    inactive: Sequence[mpmath.mpf] | None = None,
) -> tuple[mpmath.mpf, ...]:
    """Lift five active coordinates to q=1 with canonical inactive gauges."""
    a1, a2, b1, b3, b5 = point
    if inactive is None:
        inactive = (mpmath.mpf("0.5"),) * 3
    return (a1, a2, mpmath.mpf(1), *inactive, b1, b3, b5)


@dataclass
class SearchResult:
    seed_count: int
    feasible_evaluations: int
    accepted_moves: int
    accepted_results: int
    accepted_infeasible: int
    invalid_rejects: dict[str, int]
    minimum_gap: mpmath.mpf
    minimum_box: int
    minimum_point: Point
    minimum_mean: mpmath.mpf
    minimum_seed_slack: mpmath.mpf
    strict_negative_evaluations: int
    material_negative_found: bool


def deterministic_descent(
    boxes: Sequence[ActiveBox], parameters: MPParameters
) -> SearchResult:
    """Run a small coordinate descent whose gap gate accepts only feasible points."""
    mp_boxes = tuple(to_mp_box(box) for box in boxes)
    zero = mpmath.mpf(0)
    negative_tolerance = mpmath.mpf(NEGATIVE_TOLERANCE_TEXT)
    invalid_rejects = {"shape": 0, "bounds": 0, "simplex": 0, "mean": 0}
    feasible_evaluations = 0
    accepted_moves = 0
    accepted_results = 0
    accepted_infeasible = 0
    strict_negative_evaluations = 0
    minimum_gap = mpmath.inf
    minimum_box = -1
    minimum_point: Point = (zero,) * 5
    minimum_mean = mpmath.nan
    minimum_seed_slack = mpmath.inf

    def inspect_if_feasible(
        point: Point, box: MPBox, box_number: int
    ) -> tuple[mpmath.mpf, mpmath.mpf] | None:
        nonlocal feasible_evaluations
        nonlocal strict_negative_evaluations
        nonlocal minimum_gap, minimum_box, minimum_point, minimum_mean

        feasible, reason, mean = classify_qone(point, box, parameters.mean)
        if not feasible:
            invalid_rejects[reason] += 1
            return None
        # This is deliberately below the feasibility gate.  Infeasible search
        # proposals are never passed to Liu's raw-gap transcription.
        raw_gap = gap_mp(qone_values(point), parameters.beta)
        require(mpmath.isfinite(raw_gap), "a feasible raw-gap evaluation was non-finite")
        feasible_evaluations += 1
        if raw_gap < zero:
            strict_negative_evaluations += 1
        if raw_gap < minimum_gap:
            minimum_gap = +raw_gap
            minimum_box = box_number
            minimum_point = tuple(point)
            minimum_mean = +mean
        return raw_gap, mean

    for box_number, box in enumerate(mp_boxes):
        seed, seed_slack = first_feasible_corner(box, parameters.mean)
        minimum_seed_slack = min(minimum_seed_slack, seed_slack)
        inspected = inspect_if_feasible(seed, box, box_number)
        require(inspected is not None, "a selected feasible seed failed its own gate")
        current = list(seed)
        current_gap = inspected[0]
        widths = tuple(hi - lo for lo, hi in box)

        for level in range(DESCENT_LEVELS):
            divisor = mpmath.mpf(2 ** (level + 1))
            for coordinate in range(5):
                step = widths[coordinate] / divisor
                if step == zero:
                    continue
                for direction in (-1, 1):
                    proposal_list = list(current)
                    proposal_list[coordinate] += direction * step
                    proposal = tuple(proposal_list)
                    candidate = inspect_if_feasible(proposal, box, box_number)
                    if candidate is None:
                        continue
                    candidate_gap = candidate[0]
                    if candidate_gap < current_gap:
                        current = proposal_list
                        current_gap = candidate_gap
                        accepted_moves += 1

        final_point = tuple(current)
        final_feasible, _, _ = classify_qone(final_point, box, parameters.mean)
        if not final_feasible:
            accepted_infeasible += 1
        accepted_results += 1

    require(minimum_box >= 0, "the descent performed no feasible raw-gap evaluation")
    return SearchResult(
        seed_count=len(mp_boxes),
        feasible_evaluations=feasible_evaluations,
        accepted_moves=accepted_moves,
        accepted_results=accepted_results,
        accepted_infeasible=accepted_infeasible,
        invalid_rejects=invalid_rejects,
        minimum_gap=minimum_gap,
        minimum_box=minimum_box,
        minimum_point=minimum_point,
        minimum_mean=minimum_mean,
        minimum_seed_slack=minimum_seed_slack,
        strict_negative_evaluations=strict_negative_evaluations,
        material_negative_found=minimum_gap < -negative_tolerance,
    )


@dataclass(frozen=True)
class Calibration:
    candidate_gap: mpmath.mpf
    candidate_gauge_spread: mpmath.mpf
    binary_effective_gap: mpmath.mpf
    binary_all_one_gap: mpmath.mpf


def calibrate_raw_gap(parameters: MPParameters) -> Calibration:
    """Check Liu's q=1 structural zero and two feasible binary zeros."""
    zero = mpmath.mpf(0)
    half = mpmath.mpf("0.5")
    one = mpmath.mpf(1)
    tolerance = mpmath.mpf(ZERO_TOLERANCE_TEXT)
    r = (one - parameters.p) / 2

    candidate = (parameters.p, r, parameters.x, zero, zero)
    gauge_choices = (
        (zero, zero, zero),
        (half, half, half),
        (one, one, one),
    )
    candidate_gaps = tuple(
        gap_mp(qone_values(candidate, gauges), parameters.beta)
        for gauges in gauge_choices
    )
    candidate_gap = candidate_gaps[1]
    candidate_spread = max(candidate_gaps) - min(candidate_gaps)
    require(
        abs(qone_mean(candidate) - parameters.mean) <= tolerance,
        "Liu's q=1 candidate does not calibrate to the mean boundary",
    )
    require(
        max(abs(value) for value in candidate_gaps) <= tolerance,
        "Liu's q=1 candidate raw gap did not calibrate to zero",
    )
    require(
        abs(candidate_spread) <= tolerance,
        "q=1 raw gap depends on an inactive P0 support gauge",
    )

    # At the known sharp endpoint, the x-support has zero mass and the
    # positive-mass support is literally binary: mean mass at 1, the rest at 0.
    binary_effective = (zero, parameters.mean, parameters.x, one, zero)
    binary_effective_gap = gap_mp(qone_values(binary_effective), parameters.beta)
    require(
        qone_mean(binary_effective) >= parameters.mean,
        "effective-binary endpoint is not q=1 mean-feasible",
    )
    require(
        abs(binary_effective_gap) <= tolerance,
        "effective-binary boundary raw gap did not calibrate to zero",
    )

    binary_all_one = (zero, zero, one, one, one)
    binary_all_one_gap = gap_mp(qone_values(binary_all_one), parameters.beta)
    require(
        qone_mean(binary_all_one) >= parameters.mean,
        "all-one binary endpoint is not q=1 mean-feasible",
    )
    require(
        abs(binary_all_one_gap) <= tolerance,
        "all-one binary boundary raw gap did not calibrate to zero",
    )
    return Calibration(
        candidate_gap=candidate_gap,
        candidate_gauge_spread=candidate_spread,
        binary_effective_gap=binary_effective_gap,
        binary_all_one_gap=binary_all_one_gap,
    )


def aggregate_mean(values: Sequence[mpmath.mpf]) -> mpmath.mpf:
    """Direct nine-coordinate mean, used only to classify the negative control."""
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    a3 = 1 - a1 - a2
    p0_mean = a1 * b0 + a2 * b2 + a3 * b4
    p1_mean = a1 * b1 + a2 * b3 + a3 * b5
    return (1 - q) * p0_mean + q * p1_mean


@dataclass(frozen=True)
class NegativeControl:
    raw_gap: mpmath.mpf
    mean: mpmath.mpf
    mean_deficit: mpmath.mpf
    q: mpmath.mpf


def reproduce_negative_control(parameters: MPParameters) -> NegativeControl:
    """Replay the known negative raw gap only after classifying it infeasible."""
    zero = mpmath.mpf(0)
    q = mpmath.mpf("0.5")
    epsilon = mpmath.mpf(1) / 1024
    y = mpmath.mpf(1) / 32
    values = (
        parameters.p - epsilon,
        epsilon,
        q,
        parameters.x,
        y,
        zero,
        parameters.x,
        y,
        zero,
    )
    mean = aggregate_mean(values)
    # Classification is intentionally completed before the diagnostic gap call.
    require(mean < parameters.mean, "negative control unexpectedly meets the mean constraint")
    require(q != 1, "negative control unexpectedly lies on the q=1 face")
    raw_gap = gap_mp(values, parameters.beta)
    require(raw_gap < mpmath.mpf("-5e-4"), "known ambient negative control was not reproduced")
    return NegativeControl(
        raw_gap=raw_gap,
        mean=mean,
        mean_deficit=mean - parameters.mean,
        q=q,
    )


def fmt(value: mpmath.mpf, digits: int = 32) -> str:
    return mpmath.nstr(value, digits, strip_zeros=False)


def main() -> int:
    started = time.monotonic()
    mpmath.mp.dps = DPS
    require(
        "liu9_qone_face" not in sys.modules,
        "forbidden q-one implementation was imported",
    )

    source, source_digest = load_digest_checked(SOURCE_REPORT)
    projection = project_q_at_one(source)
    require(
        projection.q_at_one_count == EXPECTED_Q_AT_ONE,
        f"expected {EXPECTED_Q_AT_ONE} q-at-one boxes, got {projection.q_at_one_count}",
    )
    require(
        len(projection.boxes) == EXPECTED_UNIQUE_ACTIVE,
        f"expected {EXPECTED_UNIQUE_ACTIVE} unique active boxes, got {len(projection.boxes)}",
    )
    qone_digest = cross_check_qone_report(
        QONE_REPORT, source, source_digest, projection
    )

    parameters = solve_equation_parameters(DPS)
    calibration = calibrate_raw_gap(parameters)
    search = deterministic_descent(projection.boxes, parameters)
    require(
        search.seed_count == EXPECTED_UNIQUE_ACTIVE,
        "not every active box supplied a feasible corner seed",
    )
    require(
        search.accepted_results == EXPECTED_UNIQUE_ACTIVE,
        "not every active box supplied a local-search result",
    )
    require(
        search.accepted_infeasible == 0,
        "an infeasible local-search result was accepted",
    )
    control = reproduce_negative_control(parameters)

    invalid_total = sum(search.invalid_rejects.values())
    print("LIU H2 INDEPENDENT Q=1 ACTIVE-FACE CHECK")
    print("single-core thread caps             = 1")
    print("source canonical digest             =", source_digest)
    print("source residual boxes               =", projection.residual_count)
    print("q-at-one source boxes               =", projection.q_at_one_count)
    print("unique active boxes                 =", len(projection.boxes))
    print("active quotient                     =", ",".join(ACTIVE_NAMES))
    print("q interval counts                   =", json.dumps(projection.q_intervals, sort_keys=True))
    if qone_digest is None:
        print("stored q-one report cross-check     = NOT PRESENT (optional)")
    else:
        print("stored q-one report digest          =", qone_digest)
        print("stored q-one source cross-check     = PASS")

    print("\n[raw-gap calibrations via gap_mp]")
    print("Liu q=1 candidate raw gap           =", fmt(calibration.candidate_gap))
    print("inactive-gauge raw-gap spread       =", fmt(calibration.candidate_gauge_spread))
    print("effective-binary boundary raw gap   =", fmt(calibration.binary_effective_gap))
    print("all-one binary boundary raw gap     =", fmt(calibration.binary_all_one_gap))

    print("\n[feasibility-gated deterministic descent]")
    print("feasible corner seeds               =", search.seed_count)
    print("minimum feasible-corner mean slack  =", fmt(search.minimum_seed_slack))
    print("feasible raw-gap evaluations        =", search.feasible_evaluations)
    print("accepted descent moves              =", search.accepted_moves)
    print("accepted final results              =", search.accepted_results)
    print("accepted infeasible results         =", search.accepted_infeasible)
    print("invalid-result rejects before gap   =", invalid_total)
    print("invalid reject breakdown            =", json.dumps(search.invalid_rejects, sort_keys=True))
    print("minimum feasible raw gap            =", fmt(search.minimum_gap, 50))
    print("minimum active box (one-based)      =", search.minimum_box + 1)
    print("minimum active point                =", tuple(fmt(v, 18) for v in search.minimum_point))
    print("minimum point mean                  =", fmt(search.minimum_mean))
    print("strict negative numerical values    =", search.strict_negative_evaluations)
    print(
        "negative feasible raw gap found      =",
        "YES" if search.material_negative_found else "NO",
        f"(threshold -{NEGATIVE_TOLERANCE_TEXT})",
    )

    print("\n[negative control: rejected before interpretation]")
    print("ambient control q                   =", fmt(control.q))
    print("ambient control mean deficit        =", fmt(control.mean_deficit))
    print("ambient control raw gap             =", fmt(control.raw_gap))
    print("q=1 feasible result                 = NO")
    print("interpretation                      = mean-infeasible ambient negative control only")

    elapsed = time.monotonic() - started
    print("\nelapsed seconds                     = %.3f" % elapsed)
    require(elapsed < 240.0, "standalone checker exceeded four minutes")
    if search.material_negative_found:
        print("FINAL VERDICT: FEASIBLE NEGATIVE FOUND")
        return 2
    print("FINAL VERDICT: PASS; all accepted results are q=1 feasible and none is negative.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except VerificationError as error:
        print(f"FINAL VERDICT: FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
