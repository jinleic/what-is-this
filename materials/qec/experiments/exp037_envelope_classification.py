"""EXP-037: catalogue-wide certified CSS-envelope classification.

The programme's central question is whether *any* published perturbed
bivariate-bicycle (PBB) code escapes the CSS bivariate-bicycle envelope at its
own length -- i.e. whether perturbation ever buys a genuinely better
[[n, k, d]] than every CSS BB code of the same length.

Definitions (all bounds certified, catalogue values never trusted)
-----------------------------------------------------------------
For a PBB row Q = [[n, k, d_Q]] let U(Q) be a *verified witness* upper bound:
a nontrivial logical of symplectic weight U(Q), re-checked through two
independent GF(2) paths.  For a CSS BB code C = [[n, k_C, d_C]] let L(C) be a
*certified lower* bound: both CSS sides UNSAT at cap L(C) - 1, so d_C >=
L(C), together with a verified witness of weight L(C) making d_C = L(C)
exactly.

Q is **envelope-dominated** iff some CSS BB code C at the same n satisfies

    k_C >= k        and        d_C >= U(Q).

Then d_C >= U(Q) >= d_Q and k_C >= k, so C weakly dominates Q in [[n,k,d]]
with certificates on both sides.  The asymmetry is what makes the sweep
affordable: the PBB side needs only a cheap SAT witness, while the expensive
UNSAT work is spent once per CSS code and amortised over every PBB row at
that length.

A row that is *not* dominated is a **candidate escapee** and is reported with
the exact residual obligation (which CSS lower bound is missing, or that a
PBB lower-bound proof is required).  Candidate escapees are never promoted to
a claim by this script; they are leads for a dedicated certification run.

Encoding
--------
All decisions go through ``qec_research.distance.sat_decide`` with the sound
BB translation symmetry break enabled (premise machine-checked in
``tests/test_translation_symmetry_break.py``), monolithic CNF (FR-022), and
CNF-hash-bound records.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.distance.sat_decide import (  # noqa: E402
    ENCODING_VERSION,
    css_side_instance,
    decide_weight_bounded,
    decision_cnf_digest,
    group_weight,
    symplectic_instance,
    verify_witness_two_paths,
)

_E27_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
assert _E27_SPEC and _E27_SPEC.loader
E27 = importlib.util.module_from_spec(_E27_SPEC)
_E27_SPEC.loader.exec_module(E27)

STATE_DIR = ROOT / "results" / "partial_runs" / "exp037"
PARTIAL = ROOT / "results" / "partial_runs" / "exp037_envelope_classification.json"
PROCESSED = ROOT / "results" / "processed" / "exp037_envelope_classification.json"
SOLVER_NAME = "cadical195"
SCHEMA = "exp037-envelope-classification-v1"
MAX_CLIMB = 60


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def block_length(row: dict[str, Any]) -> int:
    return int(row["ell"]) * int(row["m"])


def pbb_instance(row: dict[str, Any]):
    _, code = E27.pbb_code(row)
    return code, symplectic_instance(
        code, code.logical_basis(), block_length=block_length(row)
    )


def css_instances(row: dict[str, Any]):
    _, HX, HZ = E27.parent_matrices(row)
    block = block_length(row)
    return (HX, HZ), {
        "x": css_side_instance(HX, HZ, "x", block_length=block),
        "z": css_side_instance(HX, HZ, "z", block_length=block),
    }


def witness_upper_bound(
    instance, hint: int, *, descent_budget: int = 2_000_000
) -> dict[str, Any]:
    """Climb to the first verified witness, then descend under a budget.

    Climbing from an untrusted hint gives *some* verified witness; the descent
    then tries cap w-1 repeatedly.  A smaller upper bound makes envelope
    domination easier to establish honestly, and an ``UNSAT`` during descent
    proves the PBB distance exactly (``exact_distance`` is then set).  The
    descent is conflict-budgeted so it can never stall the sweep.
    """

    calls: list[dict[str, Any]] = []
    cap = max(1, int(hint))
    found: dict[str, Any] | None = None
    while cap <= MAX_CLIMB:
        record = decide_weight_bounded(
            instance, cap, solver_name=SOLVER_NAME, conflict_budget=0
        )
        calls.append(
            {
                "cap": cap,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
                "phase": "climb",
            }
        )
        if record["status"] == "SAT":
            found = record
            break
        cap += 1
    if found is None:
        raise RuntimeError(f"no witness up to weight {MAX_CLIMB}")
    exact = None
    while int(found["weight"]) > 1:
        target = int(found["weight"]) - 1
        record = decide_weight_bounded(
            instance,
            target,
            solver_name=SOLVER_NAME,
            conflict_budget=descent_budget,
        )
        calls.append(
            {
                "cap": target,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
                "phase": "descent",
            }
        )
        if record["status"] == "SAT":
            found = record
            continue
        if record["status"] == "UNSAT":
            # d > target and a verified witness sits at target+1: exact.
            exact = int(found["weight"])
        break
    return {
        "upper_bound": int(found["weight"]),
        "exact_distance": exact,
        "vector": found["vector"],
        "verification": found["verification"],
        "calls": calls,
    }


def certified_css_distance(instances: dict[str, Any], hint: int) -> dict[str, Any]:
    """Exact d(CSS) = min(d_X, d_Z) via witness + both-sides UNSAT below it."""

    calls: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    for side, instance in instances.items():
        found = witness_upper_bound(instance, hint)
        calls.extend({**call, "side": side} for call in found["calls"])
        if best is None or found["upper_bound"] < best["upper_bound"]:
            best = {**found, "side": side}
    assert best is not None
    upper = best["upper_bound"]
    lower_ok = True
    for side, instance in instances.items():
        record = decide_weight_bounded(
            instance, upper - 1, solver_name=SOLVER_NAME, conflict_budget=0
        )
        calls.append(
            {
                "side": side,
                "cap": upper - 1,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
            }
        )
        if record["status"] != "UNSAT":
            lower_ok = False
    return {
        "exact": bool(lower_ok),
        "distance": int(upper) if lower_ok else None,
        "certified_upper_bound": int(upper),
        "witness_side": best["side"],
        "witness_vector": best["vector"],
        "witness_verification": best["verification"],
        "calls": calls,
    }


def row_record_path(index: int) -> Path:
    return STATE_DIR / f"row_{index:04d}.json"


def css_record_path(fingerprint: str) -> Path:
    return STATE_DIR / f"css_{fingerprint[:16]}.json"


def pbb_upper_for_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    """Verified witness upper bound for one PBB row; atomic and resumable."""

    path = row_record_path(index)
    code, instance = pbb_instance(row)
    # k is domination-critical: an understated k_Q would make "k_C >= k_Q"
    # easier to satisfy.  Recompute it from the rebuilt stabilizer and require
    # the catalogue metadata to agree.
    rebuilt_k = int(code.n - code.rank())
    if rebuilt_k != int(row["k"]) or rebuilt_k != int(code.k):
        raise RuntimeError(
            f"row {index}: k disagreement - rebuilt {rebuilt_k}, "
            f"code.k {int(code.k)}, catalogue {int(row['k'])}"
        )
    identity = {
        "schema": SCHEMA,
        "catalogue_index": index,
        "label": E27.catalogue_label(row, index),
        "n": code.n,
        "k": rebuilt_k,
        "pbb_fingerprint": E27.matrix_fingerprint(code.H),
        "encoding_version": ENCODING_VERSION,
    }
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        # Records written before bounded descent carry a loose upper bound and
        # no exactness field; treat those as stale and recompute, since a
        # smaller U makes envelope domination provable for more rows.
        if "exact_distance" in stored:
            if all(stored.get(key) == value for key, value in identity.items()):
                # Re-verify the stored witness against fresh algebra.
                vector = np.asarray(stored["vector"], dtype=np.uint8)
                verification = verify_witness_two_paths(instance, vector)
                weight = group_weight(instance, vector)
                if verification["valid"] and weight == stored["upper_bound"]:
                    return stored
            raise RuntimeError(
                f"row {index}: stored PBB record failed reverification"
            )
    found = witness_upper_bound(instance, int(row["d"]))
    payload = {**identity, **found, "utc": utc_now()}
    atomic_write_json(path, payload)
    return payload


def css_exact_for_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    """Certified exact distance of one row's CSS BB parent; atomic/resumable."""

    (HX, HZ), instances = css_instances(row)
    fingerprint = E27.matrix_fingerprint(HX, HZ)
    path = css_record_path(fingerprint)
    identity = {
        "schema": SCHEMA,
        "css_fingerprint": fingerprint,
        "n": int(HX.shape[1]),
        "encoding_version": ENCODING_VERSION,
    }
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if all(stored.get(key) == value for key, value in identity.items()):
            return stored
        raise RuntimeError(f"css {fingerprint[:16]}: identity mismatch")
    result = certified_css_distance(instances, int(row["d"]))
    # CSS dimension: k = n - rank(HX) - rank(HZ).  The X- and Z-check spaces
    # live in different symplectic halves, so ranking a stacked binary matrix
    # would let overlapping row spaces INFLATE k and make domination easier to
    # claim.  Cross-checked against the independent PBB theory module.
    rx, rz = E27.rank_np(HX), E27.rank_np(HZ)
    css_k = int(HX.shape[1]) - int(rx) - int(rz)
    theory = E27.analyse_pbb(E27.pbb_code(row)[0])
    if css_k != int(theory.k_bb):
        raise RuntimeError(
            f"css {fingerprint[:16]}: k mismatch, rank formula gives {css_k} "
            f"but analyse_pbb gives {theory.k_bb}"
        )
    payload = {
        **identity,
        **result,
        "k": css_k,
        "rank_HX": int(rx),
        "rank_HZ": int(rz),
        "k_crosscheck_analyse_pbb": int(theory.k_bb),
        "source_catalogue_index": index,
        "utc": utc_now(),
    }
    atomic_write_json(path, payload)
    return payload


def pbb_lower_bound(index: int, row: dict[str, Any], cap: int) -> dict[str, Any]:
    """Try to certify d_PBB > cap by UNSAT; persist into the row record.

    This is the decisive direction for an *escape* claim: an envelope escape
    needs d_PBB strictly above the best certified same-length CSS distance at
    k_C >= k_Q, which requires a lower bound (UNSAT), never a witness.
    """

    code, instance = pbb_instance(row)
    record = decide_weight_bounded(
        instance, cap, solver_name=SOLVER_NAME, conflict_budget=0
    )
    path = row_record_path(index)
    stored = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    lowers = stored.get("lower_bound_calls", [])
    lowers.append(
        {
            "cap": int(cap),
            "status": record["status"],
            "cnf_sha256": record["cnf_sha256"],
            "encoding_version": record["encoding_version"],
            "wall_time_s": record["solver"]["wall_time_s"],
            "utc": utc_now(),
        }
    )
    stored["lower_bound_calls"] = lowers
    proven = [
        call["cap"] for call in lowers if call["status"] == "UNSAT"
    ]
    if proven:
        stored["certified_lower_gt"] = max(proven)
    if record["status"] == "SAT":
        # A lighter witness than the stored upper bound: tighten it.
        if int(record["weight"]) < int(stored.get("upper_bound", 10**9)):
            stored["upper_bound"] = int(record["weight"])
            stored["vector"] = record["vector"]
            stored["verification"] = record["verification"]
    atomic_write_json(path, stored)
    return stored


def load_css_pool() -> list[dict[str, Any]]:
    """All certified-exact CSS BB codes available as dominators."""

    pool: list[dict[str, Any]] = []
    if STATE_DIR.exists():
        for path in sorted(STATE_DIR.glob("css_*.json")):
            record = json.loads(path.read_text(encoding="utf-8"))
            if record.get("exact") and record.get("distance") is not None:
                pool.append(record)
    return pool


def classify(pool: list[dict[str, Any]]) -> dict[str, Any]:
    rows = E27.load_catalogue()
    items: list[dict[str, Any]] = []
    counts = {"DOMINATED": 0, "CANDIDATE_ESCAPEE": 0, "PENDING": 0}
    for index, row in enumerate(rows):
        path = row_record_path(index)
        if not path.exists():
            counts["PENDING"] += 1
            continue
        stored = json.loads(path.read_text(encoding="utf-8"))
        # Classification is where domination is decided, so both critical
        # dimensions are re-derived here from rebuilt algebra and the stored
        # witness is re-verified; nothing on disk is taken on trust.
        code, instance = pbb_instance(row)
        n = int(code.n)
        k = int(code.n - code.rank())
        vector = np.asarray(stored["vector"], dtype=np.uint8)
        verification = verify_witness_two_paths(instance, vector)
        upper = group_weight(instance, vector)
        if not verification["valid"]:
            raise RuntimeError(
                f"row {index}: stored witness failed reverification"
            )
        if (n, k, upper) != (
            stored["n"],
            stored["k"],
            stored["upper_bound"],
        ):
            raise RuntimeError(
                f"row {index}: rebuilt (n,k,U)=({n},{k},{upper}) disagrees with "
                f"stored ({stored['n']},{stored['k']},{stored['upper_bound']})"
            )
        dominators = [
            candidate
            for candidate in pool
            if candidate["n"] == n
            and candidate["k"] >= k
            and candidate["distance"] >= upper
        ]
        verdict = "DOMINATED" if dominators else "CANDIDATE_ESCAPEE"
        counts[verdict] += 1
        best_css = max(
            (c for c in pool if c["n"] == n and c["k"] >= k),
            key=lambda c: c["distance"],
            default=None,
        )
        items.append(
            {
                "catalogue_index": index,
                "label": stored["label"],
                "n": n,
                "k": k,
                "pbb_upper_bound": upper,
                "catalogue_d_untrusted": int(row["d"]),
                "verdict": verdict,
                "dominator_fingerprint": (
                    dominators[0]["css_fingerprint"][:16] if dominators else None
                ),
                "dominator_k_d": (
                    [dominators[0]["k"], dominators[0]["distance"]]
                    if dominators
                    else None
                ),
                "best_same_n_css_k_d": (
                    [best_css["k"], best_css["distance"]] if best_css else None
                ),
            }
        )
    return {"counts": counts, "items": items}


def assemble() -> dict[str, Any]:
    pool = load_css_pool()
    classification = classify(pool)
    payload = {
        "experiment": "EXP-037",
        "schema": SCHEMA,
        "protocol": __doc__,
        "generated_utc": utc_now(),
        "encoding_version": ENCODING_VERSION,
        "solver": SOLVER_NAME,
        "css_pool": [
            {
                "fingerprint": candidate["css_fingerprint"][:16],
                "n": candidate["n"],
                "k": candidate["k"],
                "distance": candidate["distance"],
            }
            for candidate in sorted(pool, key=lambda c: (c["n"], -c["k"]))
        ],
        "counts": classification["counts"],
        "items": classification["items"],
    }
    complete = (
        classification["counts"]["PENDING"] == 0
        and classification["counts"]["CANDIDATE_ESCAPEE"] == 0
    )
    destination = PROCESSED if complete else PARTIAL
    atomic_write_json(destination, payload)
    payload["_written_to"] = str(destination.relative_to(ROOT))
    return payload


def run(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    targets = list(enumerate(rows))
    if args.ns:
        wanted = set(args.ns)
        targets = [(i, r) for i, r in targets if 2 * block_length(r) in wanted]
    if args.indexes:
        # Preserve the caller's ordering so cheapest-first sweeps stay cheap.
        position = {index: rank for rank, index in enumerate(args.indexes)}
        targets = [(i, r) for i, r in targets if i in position]
        targets.sort(key=lambda pair: position[pair[0]])
    if args.stride > 1:
        targets = targets[args.offset :: args.stride]
    print(f"rows selected: {len(targets)}")
    for index, row in targets:
        started = time.perf_counter()
        try:
            pbb = pbb_upper_for_row(index, row)
            note = f"U(PBB)={pbb['upper_bound']}"
            if args.with_css:
                css = css_exact_for_row(index, row)
                note += (
                    f" d(CSS)={css['distance']}"
                    if css["exact"]
                    else f" d(CSS)<={css['certified_upper_bound']} (not exact)"
                )
            print(
                f"row {index:04d} {E27.catalogue_label(row, index)}: {note} "
                f"in {time.perf_counter() - started:.1f}s"
            )
        except Exception as defect:  # noqa: BLE001 - report and continue sweep
            print(f"row {index:04d}: FAILED {defect}")
    summary = assemble()
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 0


def lower(args: argparse.Namespace) -> int:
    """Certify PBB lower bounds: UNSAT at cap => d_PBB > cap."""

    rows = E27.load_catalogue()
    for index in args.indexes:
        row = rows[index]
        path = row_record_path(index)
        if not path.exists():
            pbb_upper_for_row(index, row)
        stored = json.loads(path.read_text(encoding="utf-8"))
        cap = args.cap if args.cap is not None else int(stored["upper_bound"]) - 1
        started = time.perf_counter()
        result = pbb_lower_bound(index, row, cap)
        call = result["lower_bound_calls"][-1]
        exact = (
            result.get("certified_lower_gt") is not None
            and int(result["certified_lower_gt"]) == int(result["upper_bound"]) - 1
        )
        print(
            f"row {index:04d} {E27.catalogue_label(row, index)}: cap {cap} -> "
            f"{call['status']}; U={result['upper_bound']} "
            f"lower_gt={result.get('certified_lower_gt')}"
            + ("  EXACT" if exact else "")
            + f"  {time.perf_counter() - started:.1f}s",
            flush=True,
        )
    summary = assemble()
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--ns", type=int, nargs="*")
    runner.add_argument("--indexes", type=int, nargs="*")
    runner.add_argument("--stride", type=int, default=1)
    runner.add_argument("--offset", type=int, default=0)
    runner.add_argument(
        "--with-css",
        action="store_true",
        help="also certify the row's CSS BB parent distance exactly",
    )
    lowerer = sub.add_parser("lower")
    lowerer.add_argument("--indexes", type=int, nargs="+", required=True)
    lowerer.add_argument(
        "--cap",
        type=int,
        default=None,
        help="UNSAT target; default = stored upper bound - 1 (proves exactness)",
    )
    sub.add_parser("assemble")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run":
        return run(args)
    if args.command == "lower":
        return lower(args)
    summary = assemble()
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
