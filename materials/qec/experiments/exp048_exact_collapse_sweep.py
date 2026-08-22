"""EXP-048b — exact cap-5 census for EXP-046 D=0 demotion witnesses.

EXP-046 supplies, for each demotion-realized BB parent, a *specific* valid
sibling ``Q`` with ``C = witness.c_terms`` and ``D = 0`` plus a pure-X logical
``x0 = H_X[0]`` of weight 6.  This experiment independently rebuilds that
sibling, verifies ``x0`` with the full stabilizer code, then decides the full
symplectic distance below 6.  An UNSAT result at cap 5 together with x0
certifies ``d(Q) = d_X(Q) = 6``.  SAT results are retained as full-distance
upper bounds; they are not called X-distance certificates unless the returned
vector is pure-X.

The run is resumable and shardable.  It writes one atomic JSON record per
catalogue index under ``results/partial_runs/exp048_exact_collapse/``.  Use
``assemble`` only after all desired shards finish; no worker writes the shared
aggregate during ``run``.
"""
from __future__ import annotations

import argparse
import hashlib
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

_spec = importlib.util.spec_from_file_location("exp027_for_exp048", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

from qec_research.distance.sat_decide import (  # noqa: E402
    decide_weight_bounded,
    symplectic_instance,
)

EXP046_PATH = ROOT / "results" / "processed" / "exp046_generic_syzygy.json"
PARTIAL_DIR = ROOT / "results" / "partial_runs" / "exp048_exact_collapse"
AGGREGATE_PATH = ROOT / "results" / "processed" / "exp048_exact_collapse_sweep.json"
SCHEMA = "exp048-exact-collapse-v1"
SOLVER_NAME = "cadical195"


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def row_path(index: int) -> Path:
    return PARTIAL_DIR / f"row_{index:04d}.json"


def group_weight(vector: np.ndarray, n: int) -> int:
    vector = np.asarray(vector, dtype=np.uint8)
    return int(sum(bool(vector[j]) or bool(vector[n + j]) for j in range(n)))


def vector_sha256(vector: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(vector, dtype=np.uint8).tobytes()).hexdigest()

def expected_identity(item: dict[str, Any]) -> dict[str, Any]:
    """Rebuild identity fields used to accept or replay a row certificate."""
    index = int(item["index"])
    row = item["row"]
    parent = item["parent"]
    label = E27.catalogue_label(row, index)
    witness = parent["witness"]
    if parent.get("demotion_depth") != "D0_basis" or witness.get("d_terms"):
        raise RuntimeError(f"{label}: EXP-048 sweep accepts only D=0 witnesses")
    if not witness.get("M_zero") or not witness.get("z0_in_ker_HX"):
        raise RuntimeError(f"{label}: EXP-046 witness is not an M=0 D=0 witness")
    _, hx_parent, hz_parent = E27.parent_matrices(row)
    parent_fingerprint = E27.matrix_fingerprint(hx_parent, hz_parent)
    if parent.get("fingerprint") != parent_fingerprint:
        raise RuntimeError(f"{label}: EXP-046 parent fingerprint mismatch")
    perturbed = dict(row, C_terms=[list(t) for t in witness["c_terms"]], D_terms=[])
    _, code = E27.pbb_code(perturbed)
    delta_bar = int(witness.get("delta_bar", 0))
    k_q = int(code.n - E27.rank_np(code.H))
    expected_k_q = int(parent["k_parent"]) - delta_bar
    if k_q != expected_k_q:
        raise RuntimeError(
            f"{label}: k mismatch, rebuilt {k_q}, expected "
            f"{parent['k_parent']}-{delta_bar}={expected_k_q}"
        )
    input_sha256 = hashlib.sha256(
        json.dumps(witness["c_terms"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "schema": SCHEMA,
        "catalogue_index": index,
        "label": label,
        "n": int(code.n),
        "parent_fingerprint": parent_fingerprint,
        "pbb_fingerprint": E27.matrix_fingerprint(code.H),
        "exp046_input_sha256": input_sha256,
        "exp046_delta_bar": delta_bar,
        "expected_k_q": expected_k_q,
        "k_parent": int(parent["k_parent"]),
        "k_pbb_catalogue": int(row["k"]),
        "k_q": k_q,
        "parent_d_bound": int(parent["d_bound"]),
    }

IDENTITY_KEYS = (
    "schema",
    "catalogue_index",
    "label",
    "n",
    "parent_fingerprint",
    "pbb_fingerprint",
    "exp046_input_sha256",
    "exp046_delta_bar",
    "expected_k_q",
    "k_parent",
    "k_pbb_catalogue",
    "k_q",
    "parent_d_bound",
)


def assert_record_identity(item: dict[str, Any], record: dict[str, Any]) -> None:
    expected = expected_identity(item)
    mismatches = {
        key: (record.get(key), expected[key])
        for key in IDENTITY_KEYS
        if record.get(key) != expected[key]
    }
    if mismatches:
        raise RuntimeError(f"row identity mismatch: {mismatches}")


def load_cohort() -> tuple[list[dict[str, Any]], dict[int, dict[str, Any]]]:
    rows = E27.load_catalogue()
    exp046 = json.loads(EXP046_PATH.read_text(encoding="utf-8"))
    by_label = {E27.catalogue_label(row, i): (i, row) for i, row in enumerate(rows)}
    cohort: dict[int, dict[str, Any]] = {}
    for parent in exp046["rows"]:
        if not parent.get("demotion") or parent.get("d_bound") is None:
            continue
        label = parent["label"]
        if label not in by_label:
            raise RuntimeError(f"EXP-046 label absent from catalogue: {label}")
        index, row = by_label[label]
        witness = parent.get("witness")
        if not witness or not witness.get("c_terms"):
            raise RuntimeError(f"demotion parent lacks witness: {label}")
        cohort[index] = {"index": index, "row": row, "parent": parent}
    return rows, cohort


def certify_row(item: dict[str, Any], *, cap: int = 5) -> dict[str, Any]:
    if cap != 5:
        raise ValueError("EXP-048 exactness requires the fixed cap=5 protocol")
    index = int(item["index"])
    row = item["row"]
    parent = item["parent"]
    label = E27.catalogue_label(row, index)
    witness = parent["witness"]
    if parent.get("demotion_depth") != "D0_basis" or witness.get("d_terms"):
        raise RuntimeError(f"{label}: EXP-048 sweep accepts only D=0 witnesses")
    if not witness.get("M_zero") or not witness.get("z0_in_ker_HX"):
        raise RuntimeError(f"{label}: EXP-046 witness is not an M=0 D=0 witness")
    parent_spec, hx_parent, hz_parent = E27.parent_matrices(row)
    parent_fingerprint = E27.matrix_fingerprint(hx_parent, hz_parent)
    if parent.get("fingerprint") != parent_fingerprint:
        raise RuntimeError(f"{label}: EXP-046 parent fingerprint mismatch")
    ell = int(parent_spec.ell)
    m = int(parent_spec.m)
    seen_terms: set[tuple[int, int]] = set()
    c_terms: list[list[int]] = []
    for term in witness["c_terms"]:
        if len(term) != 2:
            raise RuntimeError(f"{label}: malformed C term {term}")
        i, j = int(term[0]), int(term[1])
        if not (0 <= i < ell and 0 <= j < m):
            raise RuntimeError(f"{label}: C term {(i, j)} outside {ell}x{m}")
        if (i, j) in seen_terms:
            raise RuntimeError(f"{label}: duplicate C term {(i, j)}")
        seen_terms.add((i, j))
        c_terms.append([i, j])
    dim = ell * m
    A = hx_parent[:, :dim]
    B = hx_parent[:, dim:]
    cvec = np.zeros(dim, dtype=np.uint8)
    for i, j in seen_terms:
        cvec[i * m + j] = 1
    if (A @ cvec % 2).any():
        raise RuntimeError(f"{label}: witness c vector is not in ker(A)")
    from qec_research.codes.bicycle import poly_matrix

    C_matrix = poly_matrix(ell, m, [tuple(t) for t in c_terms])
    M_direct = (A @ C_matrix.T) % 2
    if M_direct.any():
        raise RuntimeError(f"{label}: direct M=A C^T is nonzero")
    perturbed = dict(row, C_terms=c_terms, D_terms=[])
    _, code = E27.pbb_code(perturbed)
    n = int(code.n)
    pbb_fingerprint = E27.matrix_fingerprint(code.H)
    k_q = int(code.n - E27.rank_np(code.H))
    delta_bar = int(witness.get("delta_bar", 0))
    expected_k_q = int(parent["k_parent"]) - delta_bar
    if k_q != expected_k_q:
        raise RuntimeError(
            f"{label}: k mismatch, rebuilt {k_q}, expected "
            f"{parent['k_parent']}-{delta_bar}={expected_k_q}"
        )
    x0 = np.asarray(hx_parent[0], dtype=np.uint8)
    x0_full = np.concatenate([x0, np.zeros(n, dtype=np.uint8)])
    if int(x0.sum()) != 6:
        raise RuntimeError(f"{label}: constructor x0 has weight {int(x0.sum())}, expected 6")
    x0_valid = bool(code.is_logical(x0_full))
    if not x0_valid:
        raise RuntimeError(f"{label}: EXP-046 x0 failed independent Q logical check")

    instance = symplectic_instance(code, code.logical_basis(), block_length=n // 2)
    input_sha256 = hashlib.sha256(
        json.dumps(witness["c_terms"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    calls: list[dict[str, Any]] = []
    best_vector = x0_full
    best_weight = 6
    next_cap = int(cap)
    started = time.perf_counter()
    exact_distance: int | None = None
    terminal_status = "UNDECIDED_BUDGET"
    while next_cap >= 0:
        result = decide_weight_bounded(
            instance,
            next_cap,
            solver_name=SOLVER_NAME,
            conflict_budget=0,
        )
        call: dict[str, Any] = {
            "cap": next_cap,
            "status": result["status"],
            "solver_name": SOLVER_NAME,
            "conflict_budget": 0,
            "cnf_sha256": result["cnf_sha256"],
            "encoding_version": result.get("encoding_version"),
            "solver_wall_s": result["solver"]["wall_time_s"],
        }
        if result["status"] == "SAT":
            vector = np.asarray(result["vector"], dtype=np.uint8)
            weight = group_weight(vector, n)
            pure_x = bool(not vector[n:].any())
            call.update({
                "weight": weight,
                "pure_x": pure_x,
                "vector_sha256": vector_sha256(vector),
                "verification": result.get("verification"),
            })
            if weight >= best_weight:
                raise RuntimeError(
                    f"{label}: solver SAT model weight {weight} did not improve {best_weight}"
                )
            best_vector = vector
            best_weight = weight
            next_cap = weight - 1
            calls.append(call)
            continue
        calls.append(call)
        terminal_status = result["status"]
        if result["status"] == "UNSAT":
            exact_distance = best_weight
        break

    elapsed = time.perf_counter() - started
    best_pure_x = bool(not best_vector[n:].any())
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "experiment": "exp048_exact_collapse_sweep",
        "catalogue_index": index,
        "label": label,
        "n": n,
        "parent_fingerprint": parent_fingerprint,
        "pbb_fingerprint": pbb_fingerprint,
        "exp046_input_sha256": input_sha256,
        "exp046_delta_bar": delta_bar,
        "expected_k_q": expected_k_q,
        "k_parent": int(parent["k_parent"]),
        "k_pbb_catalogue": int(row["k"]),
        "k_q": k_q,
        "k_preserving": k_q == int(parent["k_parent"]),
        "symmetry_block_length": n // 2,
        "solver_name": SOLVER_NAME,
        "conflict_budget": 0,
        "parent_d_bound": int(parent["d_bound"]),
        "parent_d_cert_exact": parent.get("d_cert_exact"),
        "parent_d_pool_lb": parent.get("d_pool_lb"),
        "exp046_demote_depth": parent.get("demotion_depth"),
        "exp046_witness_c_terms": witness["c_terms"],
        "x0_witness": {
            "weight": 6,
            "pure_x": True,
            "vector_sha256": vector_sha256(x0_full),
            "verified_logical": x0_valid,
        },
        "calls": calls,
        "terminal_status": terminal_status,
        "best_weight": best_weight,
        "best_witness": {
            "weight": best_weight,
            "pure_x": best_pure_x,
            "vector": [int(bit) for bit in best_vector],
            "vector_sha256": vector_sha256(best_vector),
        },
        "exact_full_distance": exact_distance,
        "exact_x_distance": exact_distance if exact_distance is not None and best_pure_x else None,
        "strict_parent_drop_certified": bool(
            exact_distance is not None and exact_distance < int(parent["d_bound"])
        ),
        "wall_s": elapsed,
    }
    return payload


def run(args: argparse.Namespace) -> int:
    _, cohort = load_cohort()
    indexes = sorted(cohort)
    if args.indexes:
        wanted = set(args.indexes)
        indexes = [index for index in indexes if index in wanted]
    indexes = indexes[args.offset :: args.stride]
    print(f"rows selected: {len(indexes)}", flush=True)
    for index in indexes:
        path = row_path(index)
        if path.exists() and not args.force:
            try:
                stored = json.loads(path.read_text(encoding="utf-8"))
                if (
                    stored.get("schema") == SCHEMA
                    and stored.get("exact_full_distance") is not None
                ):
                    assert_record_identity(cohort[index], stored)
                    print(f"row {index:04d}: skip exact record", flush=True)
                    continue
            except (OSError, json.JSONDecodeError):
                pass
        started = time.perf_counter()
        try:
            payload = certify_row(cohort[index], cap=args.cap)
            atomic_json(path, payload)
            print(
                f"row {index:04d} {payload['label']}: d_bound={payload['parent_d_bound']} "
                f"best={payload['best_weight']} exact={payload['exact_full_distance']} "
                f"x_exact={payload['exact_x_distance']} in {time.perf_counter()-started:.2f}s",
                flush=True,
            )
        except Exception as exc:  # keep other shard rows running
            print(f"row {index:04d}: FAILED {type(exc).__name__}: {exc}", flush=True)
    return 0


def assemble(_: argparse.Namespace) -> int:
    _, cohort = load_cohort()
    records = []
    missing = []
    for index in sorted(cohort):
        path = row_path(index)
        if not path.exists():
            missing.append(index)
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        assert_record_identity(cohort[index], record)
        records.append(record)
    if missing:
        raise RuntimeError(f"missing row records: {missing}")
    counts: dict[str, int] = {}
    for record in records:
        key = record.get("terminal_status", "missing")
        counts[key] = counts.get(key, 0) + 1
    if any(record.get("terminal_status") != "UNSAT" for record in records):
        raise RuntimeError("not all row certificates ended UNSAT")
    exact = [r for r in records if r.get("exact_full_distance") is not None]
    strict = [r for r in exact if r.get("strict_parent_drop_certified")]
    strict_full_k = [
        r for r in strict if r.get("k_preserving") and r.get("exp046_delta_bar") == 0
    ]
    exact6 = [r for r in exact if r.get("exact_full_distance") == 6]
    payload = {
        "schema": SCHEMA,
        "experiment": "exp048_exact_collapse_sweep",
        "cohort_size": len(cohort),
        "record_count": len(records),
        "missing_indexes": missing,
        "counts_by_terminal_status": counts,
        "exact_full_distance_count": len(exact),
        "exact_distance_6_count": len(exact6),
        "strict_parent_drop_count": len(strict),
        "strict_parent_drop_full_k_count": len(strict_full_k),
        "strict_parent_drop_labels": [r["label"] for r in strict],
        "records": records,
    }
    atomic_json(AGGREGATE_PATH, payload)
    print(json.dumps({k: payload[k] for k in payload if k != "records"}, indent=2))
    print(f"wrote {AGGREGATE_PATH}")
    return 0


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--indexes", type=int, nargs="*")
    run_p.add_argument("--stride", type=int, default=1)
    run_p.add_argument("--offset", type=int, default=0)
    run_p.add_argument("--cap", type=int, default=5)
    run_p.add_argument("--force", action="store_true")
    sub.add_parser("assemble")
    return p


def main() -> int:
    args = parser().parse_args()
    if args.command == "run":
        return run(args)
    return assemble(args)


if __name__ == "__main__":
    raise SystemExit(main())
