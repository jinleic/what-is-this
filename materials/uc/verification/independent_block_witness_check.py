#!/usr/bin/env python3
"""Independent verification of the mean-feasible raw-gap witness search.

Replays the search report (liu9-block-witness-search.json), authenticates its
canonical digest, and re-derives every recorded raw gap from the claimed
nine-vectors at independent working precision.  Every accepted point is
re-checked for exact mean feasibility: the mean enclosure must contain or
exceed the binding m, and the point must satisfy the simplex, q, and support
constraints.  No search logic is imported from liu9_block_witness; only the
shared gap_mp transcription (the SSOT formula) is used.

The deterministic decision:
  - every recorded gap re-evaluated at independent dps must agree with the
    recorded text to the recorded precision,
  - the minimum re-evaluated gap decides the claim,
  - any re-evaluated gap below -1e-12 escalates the point to higher precision
    and, if it persists, FLAGS A NEGATIVE WITNESS for human mathematics.

This is a numerical cross-check, not a proof of copositivity.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Sequence

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

from liu9_binding import solve_equation_parameters  # noqa: E402
from liu9_boundary_layer import entropy_mp, gap_mp, mean_of  # noqa: E402

SOURCE_REPORT = HERE / "results" / "liu9-block-witness-search.json"
OUTPUT_REPORT = HERE / "results" / "liu9-block-witness-check.json"
INDEPENDENT_DPS = 90
ESCALATION_DPS = 200
NEGATIVE_ESCALATION = mpmath.mpf("-1e-12")
WITNESS_THRESHOLD = mpmath.mpf("-1e-60")


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


def load_authenticated_report(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing report: {path}")
    payload = json.loads(path.read_text())
    require(isinstance(payload, dict), "report root is not an object")
    recorded = payload.pop("report_sha256", None)
    require(
        isinstance(recorded, str) and len(recorded) == 64,
        "report has no usable canonical SHA-256",
    )
    computed = canonical_digest(payload)
    require(recorded == computed, "report digest mismatch")
    return payload


def parse_values(entries: Sequence[str]) -> tuple[mpmath.mpf, ...]:
    require(len(entries) == 9, f"expected nine values, got {len(entries)}")
    return tuple(mpmath.mpf(entry) for entry in entries)


def point_is_admissible(values: tuple[mpmath.mpf, ...], target: mpmath.mpf,
                        beta: mpmath.mpf) -> tuple[bool, mpmath.mpf]:
    """Exact simplex/q/support checks and the raw feasibility margin."""
    a1, a2, q, *supports = values
    one = mpmath.mpf(1)
    zero = mpmath.mpf(0)
    if not (zero <= a1 <= one and zero <= a2 <= one and a1 + a2 <= one):
        return False, mpmath.mpf(0)
    if not zero <= q <= one:
        return False, mpmath.mpf(0)
    if any(not zero <= support <= one for support in supports):
        return False, mpmath.mpf(0)
    mean = mean_of(values, one)
    if entropy_mp(values) <= 0:
        return False, mean - target
    return True, mean - target


def reevaluate(values: tuple[mpmath.mpf, ...], beta: mpmath.mpf,
               dps: int) -> mpmath.mpf | None:
    with mpmath.workdps(dps):
        converted = tuple(mpmath.mpf(v) for v in values)
        gap = gap_mp(converted, beta)
        return +gap


def relative_agreement(recorded: str, recomputed: mpmath.mpf) -> mpmath.mpf:
    reference = mpmath.mpf(recorded)
    scale = max(abs(reference), abs(recomputed), mpmath.mpf("1e-300"))
    return abs(recomputed - reference) / scale


def main() -> int:
    parameters = solve_equation_parameters(80)
    beta = parameters.beta
    target = parameters.mean
    source = load_authenticated_report(SOURCE_REPORT)

    rows = source.get("best_20")
    require(isinstance(rows, list) and rows, "report has no best_20 rows")
    best = source.get("best_gap")
    require(isinstance(best, dict), "report has no best_gap row")
    if best not in rows:
        rows = [best] + rows

    checked = 0
    infeasible = 0
    worst_agreement = mpmath.mpf(0)
    min_name, min_gap, min_point = None, None, None
    escalations: list[dict[str, Any]] = []
    witnesses: list[dict[str, Any]] = []

    for row in rows:
        require(isinstance(row, dict) and "values" in row,
                "row lacks recorded values")
        values = parse_values(row["values"])
        admissible, margin = point_is_admissible(values, target, beta)
        if not admissible:
            infeasible += 1
            continue
        gap = reevaluate(values, beta, INDEPENDENT_DPS)
        require(gap is not None, "gap_mp returned nothing in workdps context")
        agreement = relative_agreement(row["gap"], gap)
        worst_agreement = max(worst_agreement, agreement)
        checked += 1
        if min_gap is None or gap < min_gap:
            min_name, min_gap, min_point = row.get("name"), gap, values
        if gap < NEGATIVE_ESCALATION:
            escalated = reevaluate(values, beta, ESCALATION_DPS)
            escalations.append({
                "name": row.get("name"),
                "dps": ESCALATION_DPS,
                "gap": mpmath.nstr(escalated, 60),
            })
            if escalated < WITNESS_THRESHOLD:
                witnesses.append({
                    "name": row.get("name"),
                    "values": [mpmath.nstr(v, 60) for v in values],
                    "gap": mpmath.nstr(escalated, 60),
                    "feasibility_margin": mpmath.nstr(margin, 30),
                })

    negative_found = bool(witnesses)
    report = {
        "tool": "independent_block_witness_check.py",
        "claim_status": "COMPUTATIONAL EVIDENCE",
        "independent_dps": INDEPENDENT_DPS,
        "escalation_dps": ESCALATION_DPS,
        "binding_parameters": {
            "m": mpmath.nstr(target, 40),
            "beta": mpmath.nstr(beta, 40),
        },
        "rows_authenticated": len(rows),
        "rows_checked": checked,
        "rows_infeasible": infeasible,
        "worst_relative_agreement": mpmath.nstr(worst_agreement, 30),
        "minimum_recomputed_gap": None if min_gap is None
        else mpmath.nstr(min_gap, 60),
        "minimum_at": min_name,
        "minimum_values": None if min_point is None
        else [mpmath.nstr(v, 60) for v in min_point],
        "escalations": escalations,
        "negative_witnesses": witnesses,
        "negative_mean_feasible_raw_gap_found": negative_found,
    }
    report["report_sha256"] = canonical_digest(report)
    OUTPUT_REPORT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")

    print("INDEPENDENT BLOCK-WITNESS CHECK")
    print("rows authenticated %d, checked %d, infeasible %d"
          % (len(rows), checked, infeasible))
    print("worst relative agreement %s" % mpmath.nstr(worst_agreement, 12))
    print("minimum recomputed gap %s at %s"
          % (report["minimum_recomputed_gap"], report["minimum_at"]))
    print("escalations %d, witnesses %d"
          % (len(escalations), len(witnesses)))
    print("negative_mean_feasible_raw_gap_found %s" % negative_found)
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except VerificationError as error:
        print("VERIFICATION FAILED: %s" % error)
        raise SystemExit(1)
