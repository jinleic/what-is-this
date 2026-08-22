"""EXP-038: zero-search structural domination via minimum-weight survival.

Tests the criterion in ``qec_research.codes.pbb_survival`` against every row
EXP-036 has already decided, then applies it to the rows EXP-036 left open.

The falsification requirement is the point of the experiment: the criterion is
*only* trustworthy if it certifies domination on rows independently proved
dominated and **fails on all seven certified reversals**.  A criterion that
"proves" domination for a known reversal is refuted, not useful.

Inputs are rebuilt from the catalogue every run; the parent's minimum-weight
Z-logical comes from a certified exact CSS distance computation, never from a
catalogue field.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.pbb_survival import (  # noqa: E402
    dressing_space,
    survival_test,
)
from qec_research.distance.sat_decide import (  # noqa: E402
    css_side_instance,
    decide_weight_bounded,
    verify_witness_two_paths,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E27
_SPEC.loader.exec_module(E27)

OUT = ROOT / "results" / "partial_runs" / "exp038_survival_criterion.json"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp038"
E36_DIR = ROOT / "results" / "partial_runs" / "exp036"
SOLVER_NAME = "cadical195"
SCHEMA = "exp038-survival-v1"
MAX_CLIMB = 64


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def cached_parent_distance(fingerprint: str) -> dict[str, Any] | None:
    """A certified-exact CSS distance for this parent from EXP-037's pool.

    ``certified_css_distance`` proves UNSAT at ``distance - 1`` on *both* sides
    before setting ``exact``, so an exact record certifies ``d_Z >= distance``
    regardless of which side supplied the witness.  That lower bound is the
    expensive half; the matching Z-side witness is then a single SAT call at cap
    ``distance``, which must succeed at weight exactly ``distance`` when the Z
    side is minimising and is UNSAT otherwise.
    """

    path = ROOT / "results" / "partial_runs" / "exp037" / f"css_{fingerprint[:16]}.json"
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    if not record.get("exact") or record.get("encoding_version") != "sat-decide-v2":
        return None
    return record


def parent_min_weight_logical(
    HX: np.ndarray, HZ: np.ndarray, hint: int, block_length: int
) -> tuple[np.ndarray, int, list[dict[str, Any]]]:
    """A verified minimum-weight Z-logical of the parent, with its exact weight.

    Climb to the first SAT weight, then prove the weight below is UNSAT, so the
    returned weight is exactly ``d_Z(parent)`` -- never an upper bound.  The
    sound BB translation symmetry break is enabled: it is worth ~6x on exactly
    these UNSAT calls, which dominate the runtime.
    """

    instance = css_side_instance(HX, HZ, "z", block_length=block_length)
    calls: list[dict[str, Any]] = []
    found: dict[str, Any] | None = None
    cap = max(1, int(hint))
    while cap <= MAX_CLIMB:
        record = decide_weight_bounded(
            instance, cap, solver_name=SOLVER_NAME, conflict_budget=0
        )
        calls.append(
            {
                "weight_cap": cap,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
            }
        )
        if record["status"] == "SAT":
            found = record
            break
        if record["status"] != "UNSAT":
            raise RuntimeError(f"parent d_Z undecided at cap {cap}")
        cap += 1
    if found is None:
        raise RuntimeError(f"no parent Z-logical up to weight {MAX_CLIMB}")

    weight = int(found["weight"])
    if weight > 1:
        below = decide_weight_bounded(
            instance, weight - 1, solver_name=SOLVER_NAME, conflict_budget=0
        )
        calls.append(
            {
                "weight_cap": weight - 1,
                "status": below["status"],
                "cnf_sha256": below["cnf_sha256"],
                "wall_time_s": below["solver"]["wall_time_s"],
            }
        )
        if below["status"] != "UNSAT":
            raise RuntimeError(
                f"weight {weight} not minimal: cap {weight - 1} returned "
                f"{below['status']}"
            )

    vector = np.asarray(found["vector"], dtype=np.uint8)
    checked = verify_witness_two_paths(instance, vector)
    if not checked["valid"]:
        raise RuntimeError(f"parent witness failed verification: {checked}")
    # css_side_instance works in the n-bit Z sector already.
    return vector, weight, calls


def row_state_path(index: int) -> Path:
    return STATE_DIR / f"row_{index:04d}.json"


def known_verdicts() -> dict[int, str]:
    """EXP-036's independently established verdicts, keyed by catalogue index."""

    out: dict[int, str] = {}
    if not E36_DIR.exists():
        return out
    for path in sorted(E36_DIR.glob("row_*.json")):
        state = json.loads(path.read_text(encoding="utf-8"))
        index = state.get("catalogue_index")
        verdict = state.get("verdict")
        if index is not None and verdict:
            out[int(index)] = verdict
    return out


def analyse_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    """Run the structural criterion on one row; atomic and resumable."""

    path = row_state_path(index)
    _, HX, HZ = E27.parent_matrices(row)
    _, code = E27.pbb_code(row)
    n = HX.shape[1]
    k_parent = n - rank_np(HX) - rank_np(HZ)
    k_pbb = code.k
    fingerprint = E27.matrix_fingerprint(HX, HZ)
    identity = {
        "schema": SCHEMA,
        "catalogue_index": index,
        "label": E27.catalogue_label(row, index),
        "n": int(n),
        "parent_fingerprint": fingerprint,
        "pbb_fingerprint": E27.matrix_fingerprint(code.H),
    }
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if all(stored.get(key) == value for key, value in identity.items()):
            return stored

    CD = code.H[: HX.shape[0], n:]
    hint = int(row.get("d") or 2)
    block_length = int(row["ell"]) * int(row["m"])
    instance = css_side_instance(HX, HZ, "z", block_length=block_length)
    cached = cached_parent_distance(fingerprint)
    witness = None
    calls: list[dict[str, Any]] = []
    if cached is not None:
        # The pool certifies d_Z >= distance on both sides.  A SAT hit at that
        # cap is therefore automatically minimum-weight; UNSAT means the Z side
        # is strictly above and we fall through to the general climb.
        certified_lower = int(cached["distance"])
        record = decide_weight_bounded(
            instance, certified_lower, solver_name=SOLVER_NAME, conflict_budget=0
        )
        calls.append(
            {
                "reused_css_pool_lower_bound": cached["css_fingerprint"][:16],
                "certified_d_ge": certified_lower,
                "weight_cap": certified_lower,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
            }
        )
        if record["status"] == "SAT":
            witness = np.asarray(record["vector"], dtype=np.uint8)
            d_z_parent = int(record["weight"])
            if d_z_parent != certified_lower:
                raise RuntimeError(
                    f"row {index}: SAT weight {d_z_parent} below certified "
                    f"lower bound {certified_lower}"
                )
            checked = verify_witness_two_paths(instance, witness)
            if not checked["valid"]:
                raise RuntimeError(f"row {index}: pooled-path witness invalid")
    if witness is None:
        witness, d_z_parent, more = parent_min_weight_logical(
            HX, HZ, max(2, hint - 4), block_length
        )
        calls.extend(more)
    verdict = survival_test(
        HX,
        HZ,
        CD,
        witness,
        int(row["ell"]),
        int(row["m"]),
        k_parent=k_parent,
        k_pbb=k_pbb,
    )
    if not verdict.dimension_identity_holds:
        raise RuntimeError(
            f"row {index}: k_parent - k_pbb = {k_parent - k_pbb} but "
            f"dim Delta = {verdict.dim_delta}"
        )

    payload = {
        **identity,
        "k_parent": int(k_parent),
        "k_pbb": int(k_pbb),
        "dim_delta": verdict.dim_delta,
        "orbit_rank": verdict.orbit_rank,
        "counting_margin": verdict.counting_margin,
        "d_z_parent_exact": int(d_z_parent),
        "structurally_dominated": bool(verdict.dominated),
        "survivor_weight": verdict.survivor_weight,
        "survivor": (
            [int(x) for x in np.flatnonzero(verdict.survivor)]
            if verdict.survivor is not None
            else None
        ),
        "parent_distance_calls": calls,
        "detail": verdict.detail,
        "utc": E27.utc_now(),
    }
    atomic_write_json(path, payload)
    return payload


def assemble() -> dict[str, Any]:
    known = known_verdicts()
    items: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    conflicts: list[dict[str, Any]] = []
    for path in sorted(STATE_DIR.glob("row_*.json")) if STATE_DIR.exists() else []:
        state = json.loads(path.read_text(encoding="utf-8"))
        index = int(state["catalogue_index"])
        prior = known.get(index)
        dominated = bool(state["structurally_dominated"])
        bucket = (
            "STRUCTURALLY_DOMINATED" if dominated else "CRITERION_INCONCLUSIVE"
        )
        counts[bucket] = counts.get(bucket, 0) + 1
        # The falsification gate: the criterion may never claim domination on a
        # row independently certified as a reversal.
        if dominated and prior == "CERTIFIED_REVERSAL":
            conflicts.append(
                {
                    "catalogue_index": index,
                    "label": state["label"],
                    "prior_verdict": prior,
                    "survivor_weight": state["survivor_weight"],
                    "d_z_parent_exact": state["d_z_parent_exact"],
                }
            )
        items.append(
            {
                "catalogue_index": index,
                "label": state["label"],
                "n": state["n"],
                "k_parent": state["k_parent"],
                "k_pbb": state["k_pbb"],
                "dim_delta": state["dim_delta"],
                "orbit_rank": state["orbit_rank"],
                "counting_margin": state["counting_margin"],
                "d_z_parent_exact": state["d_z_parent_exact"],
                "structurally_dominated": dominated,
                "prior_exp036_verdict": prior,
                "verdict": bucket,
            }
        )
    agreed = [
        item
        for item in items
        if item["structurally_dominated"]
        and item["prior_exp036_verdict"] == "DOMINATION_PROVED"
    ]
    newly_closed = [
        item
        for item in items
        if item["structurally_dominated"]
        and item["prior_exp036_verdict"] in (None, "UNDECIDED_BUDGET")
    ]
    payload = {
        "experiment": "EXP-038",
        "schema": SCHEMA,
        "protocol": __doc__,
        "counts": counts,
        "falsification_conflicts": conflicts,
        "criterion_refuted": bool(conflicts),
        "agreed_with_exp036_dominations": len(agreed),
        "newly_closed_rows": [
            {k: item[k] for k in ("catalogue_index", "label", "n", "counting_margin")}
            for item in newly_closed
        ],
        "items": items,
        "utc": E27.utc_now(),
    }
    destination = OUT
    payload["_written_to"] = str(destination.relative_to(ROOT))
    atomic_write_json(destination, payload)
    return payload


def run(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    targets = list(enumerate(rows))
    if args.indexes:
        wanted = set(args.indexes)
        targets = [(i, r) for i, r in targets if i in wanted]
    elif args.ns:
        wanted = set(args.ns)
        targets = [
            (i, r) for i, r in targets if 2 * int(r["ell"]) * int(r["m"]) in wanted
        ]
    if args.stride > 1:
        targets = targets[args.offset :: args.stride]
    print(f"rows selected: {len(targets)}", flush=True)
    for index, row in targets:
        started = time.perf_counter()
        try:
            result = analyse_row(index, row)
            print(
                f"row {index:04d} {result['label']}: "
                f"k {result['k_parent']}->{result['k_pbb']} "
                f"dimD={result['dim_delta']} orbit_rank={result['orbit_rank']} "
                f"margin={result['counting_margin']:+d} "
                f"d_Z(parent)={result['d_z_parent_exact']} "
                f"{'DOMINATED' if result['structurally_dominated'] else 'inconclusive'}"
                f" in {time.perf_counter() - started:.1f}s",
                flush=True,
            )
        except Exception as defect:  # noqa: BLE001 - record and continue
            print(f"row {index:04d}: FAILED {defect}", flush=True)
    summary = assemble()
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))
    if summary["criterion_refuted"]:
        print("CRITERION REFUTED:", json.dumps(summary["falsification_conflicts"]))
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--indexes", type=int, nargs="*")
    runner.add_argument("--ns", type=int, nargs="*")
    runner.add_argument("--stride", type=int, default=1)
    runner.add_argument("--offset", type=int, default=0)
    sub.add_parser("assemble")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run":
        return run(args)
    summary = assemble()
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 1 if summary["criterion_refuted"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
