"""EXP-051: exhaustive X-monotone channel queries (replaces the capped exp023 scan).

The X-M claim for a demotion-immune parent needs: no multi-row light element of
S_X (= rowspace of H_X, the CSS parent's X-check space) exists below d_P.  The
prior basis+pairs scan (``j5_syz_census.py``, ``max_w_lam=3`` + a 6000-entry
cap) is not a certificate: overlap \u2264 2 bounds only two-row sums; three
weight-6 rows can sum to weight 6.

Method (one-shot CP-SAT, ``channel_query``): variables are the physical
n-bit word; rowspace membership is encoded by even-parity constraints against a
GF(2) basis of the dual kernel, weight is bounded by 1..d_P-1, and an exclusion
clause (one differing coordinate) forbids each original translated row of H_X.
A feasible solution is a multi-row light channel: verdict CHANNEL_PRESENT
(monotonicity not decidable by light exclusion — it does NOT refute
monotonicity itself).  Status INFEASIBLE is the completeness certificate:
verdict CERTIFIED_NO_MULTIROW_LIGHT.  Every CHANNEL_PRESENT witness is
independently re-verified (``verify``): physical weight, rank identity
rank([H_X;w]) == rank(H_X), and exclusion from all original rows.

Usage: ``run [--include-open]`` (10 parents: 8 X-M targets + 2 X-U open
targets), ``verify`` (re-checks and persists witness verdicts), ``assemble
--include-open`` (strict: requires the full target set; headline
``all_xm_certified`` is computed only over the 8 X-M targets).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location(
    "exp023_light_gensets", ROOT / "experiments" / "exp023_light_gensets.py"
)
assert _SPEC and _SPEC.loader
E23 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(E23)

_SPEC39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
assert _SPEC39 and _SPEC39.loader
E39 = importlib.util.module_from_spec(_SPEC39)
sys.modules[_SPEC39.name] = E39
_SPEC39.loader.exec_module(E39)
E27 = E39.E27  # noqa: N816

SCHEMA = "exp051-xmonotone-lightscan-v1"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp051"
OUT = ROOT / "results" / "processed" / "exp051_xmonotone_lightscan.json"

# label -> certified parent distance d_P (cap = d_P - 1)
TARGETS: dict[str, int] = {
    # d_P = 8 open-window monotone incumbent
    "9_6_0175": 8,
    # phase2 d_P = 6 = wt(A)+wt(B) family
    "phase2_75": 6,
    "phase2_76": 6,
    "phase2_77": 6,
    "phase2_83": 6,
    "phase2_84": 6,
    "phase2_87": 6,
    "phase2_109": 6,
}
# The two X-U immune parents (family-closed, X-monotonicity still open):
OPEN_TARGETS: dict[str, int] = {
    "15_6_0256": 10,
    "30_6_0289": 10,
}
# Demoting opposite-number to the 9_6_0175 profile, used as a POSITIVE control:
# multi-row light below its d_P is expected to exist wt=8 rows there; omit.


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


from ortools.sat.python import cp_model  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402

SCAN_TIME_LIMIT_S = 900.0


def channel_query(label: str, HX: np.ndarray, cap: int) -> dict:
    """One-shot CP-SAT: some element of rowspace(HX) with weight in [1, cap]

    that is NOT one of the original (translated) rows of HX.

    Membership ``x in rowspace(HX)`` is encoded against a GF(2) dual basis
    (x . z = 0 for every z in the left-kernel of HX); rank-deficient
    generators are fine because only the row SPACE matters.  Exclusion clauses
    forbid every original row pattern, so a feasible solution is a multi-row
    light element and INFEASIBLE (OPTIMAL) certifies none exist.
    """

    n = HX.shape[1]
    model = cp_model.CpModel()
    x = [model.new_bool_var(f"x_{j}") for j in range(n)]
    dual = nullspace_np(HX)  # rows span (rowspace HX)^perp; x in rowspace iff x.d = 0
    dual_rows = dual.shape[0]
    # x in rowspace(HX) <=> x . z = 0 for every z spanning (rowspace)^perp
    for i in range(dual_rows):
        support = np.flatnonzero(dual[i]).tolist()
        if not support:
            continue
        # parity == 0: XOR of support literals and a 1-target trick via E23 helper
        # (XOR(inputs, 0) == 0); bind target to a model constant 0 through a var.
        zero = model.new_int_var(0, 0, f"z_{i}")  # value fixed at 0
        E23.add_xor_equality(model, [x[j] for j in support], zero)  # type: ignore[arg-type]
    model.add(sum(x) >= 1)
    model.add(sum(x) <= cap)
    # exclude every original row pattern: x != r per clause
    excluded = 0
    for idx in range(HX.shape[0]):
        supp = set(np.flatnonzero(HX[idx]).tolist())
        lits = [x[j] for j in range(n) if j not in supp] + [
            x[j].Not() for j in supp
        ]
        if lits:
            model.add_bool_or(lits)  # at least one coordinate differs
            excluded += 1
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 1
    solver.parameters.max_time_in_seconds = SCAN_TIME_LIMIT_S
    solver.parameters.random_seed = E23.derived_seed(f"xchannel:{label}", cap, 0)
    started = time.perf_counter()
    status = solver.solve(model)
    wall = time.perf_counter() - started
    status_name = solver.status_name(status)
    witness = None
    witness_weight = None
    if status_name in ("OPTIMAL", "FEASIBLE"):
        word = np.array([solver.value(v) for v in x], dtype=np.uint8)
        witness = word.tolist()
        witness_weight = int(word.sum())
    complete = status_name in ("INFEASIBLE",)
    return {
        "status": status_name,
        "complete": complete,
        "dual_ranks": dual_rows,
        "excluded_rows": excluded,
        "time_limit_s": SCAN_TIME_LIMIT_S,
        "wall_s": round(wall, 6),
        "witness": witness,
        "witness_weight": witness_weight,
    }


def scan_parent(fp: str, row: dict, label: str, d_p: int) -> dict:
    t0 = time.perf_counter()
    _, HX, _ = E27.parent_matrices(row)
    cap = d_p - 1
    result = channel_query(label, HX, cap)
    if result["complete"]:
        verdict = "CERTIFIED_NO_MULTIROW_LIGHT"
    elif result["status"] in ("OPTIMAL", "FEASIBLE"):
        verdict = "CHANNEL_PRESENT"
    else:
        verdict = "UNDECIDED"
    return {
        "schema": SCHEMA,
        "parent_fingerprint": fp,
        "label": label,
        "d_p": d_p,
        "cap": cap,
        "query_status": result["status"],
        "complete": result["complete"],
        "dual_ranks": result["dual_ranks"],
        "excluded_rows": result["excluded_rows"],
        "channel_witness_weight": result["witness_weight"],
        "channel_witness": result["witness"],
        "wall_s": round(time.perf_counter() - t0, 3),
        "verdict": verdict,
        "utc": utc_now(),
    }


def verify_record(path: Path, rows: list[dict], parents: dict) -> dict:
    """Independent verification of a CHANNEL_PRESENT witness record."""

    rec = json.loads(path.read_text(encoding="utf-8"))
    if rec.get("verdict") != "CHANNEL_PRESENT" or not rec.get("channel_witness"):
        return rec
    fp = rec["parent_fingerprint"]
    entry = parents[fp]
    row = rows[entry["members"][0]["catalogue_index"]]
    _, HX, _ = E27.parent_matrices(row)
    w = np.array(rec["channel_witness"], dtype=np.uint8)
    checks = {
        "weight": int(w.sum()),
        "declared_weight": rec["channel_witness_weight"],
        "in_rowspace": rank_np(np.vstack([HX, w])) == rank_np(HX),
        "not_original_row": not any(
            np.array_equal(HX[i].astype(np.uint8), w) for i in range(HX.shape[0])
        ),
        "weight_declared_match": int(w.sum()) == int(rec["channel_witness_weight"]),
    }
    rec["independent_verification"] = checks
    rec["independent_verified"] = all(
        checks[k] for k in ("in_rowspace", "not_original_row", "weight_declared_match")
    ) and int(w.sum()) <= rec["cap"]
    if not rec["independent_verified"]:
        raise RuntimeError(f"{rec['label']}: channel witness failed independent verification: {checks}")
    atomic_write_json(path, rec)
    return rec


def run(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    wanted = dict(TARGETS)
    if args.include_open:
        wanted.update(OPEN_TARGETS)
    for fp, entry in parents.items():
        label = entry["members"][0]["label"]
        if label not in wanted:
            continue
        path = STATE_DIR / f"scan_{fp}.json"
        if path.exists() and not args.force:
            print(f"{label}: cached", flush=True)
            continue
        payload = scan_parent(fp, rows[entry["members"][0]["catalogue_index"]], label, wanted[label])
        atomic_write_json(path, payload)
        print(
            f"{label}: d_P={payload['d_p']} cap={payload['cap']} "
            f"status={payload['query_status']} channel_wt={payload['channel_witness_weight']} "
            f"-> {payload['verdict']} ({payload['wall_s']}s)",
            flush=True,
        )
    return 0


def verify_cmd_subp(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    ok = 0
    for p in sorted(STATE_DIR.glob("scan_*.json")):
        rec = verify_record(p, rows, parents)
        if rec.get("verdict") == "CHANNEL_PRESENT":
            print(f"{rec['label']}: CHANNEL_PRESENT verified={rec['independent_verified']}", flush=True)
            ok += rec.get("independent_verified", False)
    print(f"independent verifications passed: {ok}")
    return 0


def assemble(args: argparse.Namespace) -> int:
    records = [
        json.loads(p.read_text(encoding="utf-8"))
        for p in sorted(STATE_DIR.glob("scan_*.json"))
        if json.loads(p.read_text(encoding="utf-8")).get("schema") == SCHEMA
    ]
    expected = dict(TARGETS)
    if args.include_open:
        expected.update(OPEN_TARGETS)
    labels = {r["label"] for r in records}
    missing = sorted(set(expected) - labels)
    if missing:
        print(f"missing targets: {missing}; partial assembly not written")
        return 1
    flattened = {**TARGETS, **OPEN_TARGETS}
    for r in records:
        if r["verdict"] not in (
            "CERTIFIED_NO_MULTIROW_LIGHT",
            "CHANNEL_PRESENT",
            "UNDECIDED",
        ):
            raise RuntimeError(f"{r['label']}: unknown verdict {r['verdict']}")
        if flattened.get(r["label"]) != r["d_p"]:
            raise RuntimeError(
                f"{r['label']}: d_p {r['d_p']} != expected {flattened.get(r['label'])}"
            )
        if not r.get("parent_fingerprint") or len(r["parent_fingerprint"]) != 64:
            raise RuntimeError(f"{r['label']}: missing/short parent fingerprint")
        if r["verdict"] == "CHANNEL_PRESENT" and not r.get("independent_verified"):
            raise RuntimeError(
                f"{r['label']}: CHANNEL_PRESENT record lacks independent verification"
                " (run the verify subcommand first)"
            )
    targets_records = [r for r in records if r["label"] in TARGETS]
    payload = {
        "schema": SCHEMA,
        "experiment": "exp051_xmonotone_lightscan",
        "utc": utc_now(),
        "verdicts": {r["label"]: r["verdict"] for r in records},
        "all_xm_certified": all(
            r["verdict"] == "CERTIFIED_NO_MULTIROW_LIGHT" for r in targets_records
        )
        if targets_records
        else False,
        "records": records,
    }
    atomic_write_json(OUT, payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "records"}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run")
    p.add_argument("--include-open", action="store_true")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=run)
    p = sub.add_parser("verify")
    p.set_defaults(func=verify_cmd_subp)
    p = sub.add_parser("assemble")
    p.add_argument("--include-open", action="store_true")
    p.set_defaults(func=assemble)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
