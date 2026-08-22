"""EXP-027: certified delta>0 parent-BB domination audit.

The catalogue distance is always treated as an upper bound unless this
experiment or a persisted certificate proves both directions.  Conversely,
parent domination uses only a certified lower bound.  Thus the decisive test is

    parent_d_lower_bound >= pbb_d_upper_bound.

A strict reversal requires the opposite certified separation,

    pbb_d_lower_bound > parent_d_upper_bound,

and is accepted only after rebuilding both codes and repeating both distance
searches with a second CP-SAT seed.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import sys
import time
# Cap implicit BLAS/OpenMP fan-out before NumPy or multiprocessing imports.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_env] = "1"
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import numpy as np

# ---------------------------------------------------------------------------
# Pre-registered protocol constants.  Do not tune these after seeing outcomes.
# ---------------------------------------------------------------------------
EXPERIMENT = "EXP-027"
BASE_SEED = 20260812
PARENT_EXACT_MAX_N = 108
PBB_EXACT_MAX_N = 72
TIME_LIMIT_S = 300.0
SOLVER_WORKERS = 4
CONCURRENT_PARENT_JOBS = 2  # 2 * 4 = the shared-machine cap of 8 solver workers.
EXPECTED_CATALOGUE_ROWS = 368
EXPECTED_DELTA_POSITIVE_ROWS = 155
# An explicit non-default override exists only to verify partial-run routing.
# Such a run is clean-but-partial and MUST NOT touch any canonical EXP-027 path.
TIME_LIMIT_OVERRIDE_S = (
    float(sys.argv[sys.argv.index("--time-limit") + 1])
    if "--time-limit" in sys.argv
    else TIME_LIMIT_S
)
FULL_DECLARED_PROTOCOL = (
    TIME_LIMIT_OVERRIDE_S == TIME_LIMIT_S
    and not any(arg.startswith("--only") or arg.startswith("--subset") for arg in sys.argv)
)

ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = (
    ROOT
    / "third_party"
    / "qcode-discovery"
    / "results"
    / "campaign7_publication_merged.jsonl"
)
CERTIFICATE_DIR = ROOT / "results" / "certificates"
RAW_DIR = ROOT / "results" / "raw"
PROCESSED_PATH = ROOT / "results" / "processed" / "exp027_delta_audit.json"
REPORT_PATH = ROOT / "notes" / "agent_reports" / "exp027_delta_audit.md"
CHECKPOINT_PATH = ROOT / "results" / "partial_runs" / "exp027_checkpoint.json"

PROTOCOL: dict[str, Any] = {
    "experiment": EXPERIMENT,
    "catalogue": str(CATALOGUE.relative_to(ROOT)),
    "expected_catalogue_rows": EXPECTED_CATALOGUE_ROWS,
    "expected_delta_positive_rows": EXPECTED_DELTA_POSITIVE_ROWS,
    "scope": "every catalogue row with delta > 0",
    "parent_distance_priority": [
        "CERTIFIED_EXACT results/certificates/*.json matrix match",
        "BRAVYI_BB literature/name-label matrix match (recorded but uncertified)",
        f"exact_distance_css for n <= {PARENT_EXACT_MAX_N}",
        "UNSETTLED otherwise",
    ],
    "pbb_distance_priority": [
        "CERTIFIED_EXACT results/certificates/*.json matrix match",
        f"exact_distance_symplectic for n <= {PBB_EXACT_MAX_N}",
        "catalogue d as UPPER_BOUND only",
    ],
    "time_limit_s_per_cp_sat_sector": TIME_LIMIT_S,
    "solver_workers_per_distance_call": SOLVER_WORKERS,
    "concurrent_parent_distance_calls": CONCURRENT_PARENT_JOBS,
    "maximum_concurrent_solver_workers": CONCURRENT_PARENT_JOBS * SOLVER_WORKERS,
    "thread_environment": {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    },
    "parent_solver_search_cap": None,
    "parent_solver_mode": (
        "full exact optimization; no PBB-derived cap is passed to exact_distance_css"
    ),
    "domination_rule": (
        "certified parent_d_lower_bound >= certified pbb_d_upper_bound"
    ),
    "certified_reversal_rule": (
        "certified pbb_d_lower_bound > certified parent_d_upper_bound"
    ),
    "reversal_confirmation": (
        "rebuild from raw terms, recheck ranks, rerun both CP-SAT distances "
        "with a distinct derived seed"
    ),
    "catalogue_d_direction": "UPPER_BOUND regardless of catalogue d_is_exact",
    "bravyi_label_direction": (
        "literature_reported_uncertified; never a certified lower bound"
    ),
}

sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB,
    BBSpec,
    PBBSpec,
    bb_stabilizer,
    build_bb,
    build_pbb,
)
from qec_research.codes.pbb_theory import analyse_pbb  # noqa: E402
from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.distance import exact as exact_module  # noqa: E402
from qec_research.distance.exact import (  # noqa: E402
    exact_distance_css,
    exact_distance_symplectic,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402
from qec_research.symplectic.core import symplectic_weight  # noqa: E402


# These values are parsed from BRAVYI_BB's literature/name labels only.  They
# are not local distance certificates and MUST NOT supply a certified lower
# bound.  The four locally certified instances are matched independently by
# load_certificate_indexes() and take first priority.
REFERENCE_DISTANCE: dict[str, int] = {
    "[[72,12,6]]": 6,
    "[[90,8,10]]": 10,
    "[[108,8,10]]": 10,
    "[[144,12,12]]": 12,
    "[[288,12,18]]": 18,
    "[[360,12,<=24]]": 24,
    "[[784,24,24]]": 24,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def derived_seed(namespace: str) -> int:
    digest = hashlib.sha256(f"{BASE_SEED}:{namespace}".encode()).digest()
    return 1 + int.from_bytes(digest[:4], "big") % (2**31 - 2)


SEEDS: dict[str, Any] = {
    "base_seed": BASE_SEED,
    "derivation": "1 + uint32_be(SHA256(f'{base_seed}:{namespace}')[:4]) mod (2^31-2)",
    "primary_namespace": "exp027:primary:<matrix fingerprint>",
    "verification_namespace": "exp027:verification:<matrix fingerprint>",
}


def protocol_fingerprint() -> str:
    payload = json.dumps(
        {"protocol": PROTOCOL, "seeds": SEEDS}, sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(payload).hexdigest()


@contextmanager
def cp_sat_seed(seed: int) -> Iterator[None]:
    """Set the otherwise-unexposed CP-SAT seed used by exact.py.

    exact_distance_css and exact_distance_symplectic remain the executed APIs;
    this narrow process-local constructor wrapper only pre-sets random_seed on
    each solver they create.  Parent jobs run in separate processes, and PBB
    calls are sequential, so the temporary replacement cannot race.
    """

    original = exact_module.cp_model.CpSolver

    def seeded_solver():
        solver = original()
        solver.parameters.random_seed = int(seed)
        return solver

    exact_module.cp_model.CpSolver = seeded_solver
    try:
        yield
    finally:
        exact_module.cp_model.CpSolver = original


def matrix_fingerprint(*matrices: np.ndarray) -> str:
    h = hashlib.sha256()
    for matrix in matrices:
        a = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1)
        h.update(len(a.shape).to_bytes(1, "big"))
        for dim in a.shape:
            h.update(int(dim).to_bytes(8, "big"))
        h.update(a.tobytes())
    return h.hexdigest()


def json_safe(
    value: Any,
    *,
    _depth: int = 0,
    _seen: set[int] | None = None,
    _path: str = "$",
) -> Any:
    """Convert to JSON primitives without following cycles or unknown iterables."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if _depth > 64:
        raise TypeError(
            "json_safe depth limit exceeded "
            f"at {_path}: type={type(value).__module__}.{type(value).__qualname__}"
        )
    if isinstance(value, np.generic):
        return json_safe(
            value.item(), _depth=_depth + 1, _seen=_seen, _path=_path
        )
    if isinstance(value, Path):
        return str(value)

    if _seen is None:
        _seen = set()
    tracked = isinstance(value, (dict, list, tuple, set, np.ndarray)) or is_dataclass(value)
    if tracked:
        identity = id(value)
        if identity in _seen:
            raise TypeError(
                "json_safe cycle detected "
                f"at {_path}: type={type(value).__module__}.{type(value).__qualname__}"
            )
        _seen.add(identity)
        try:
            if isinstance(value, dict):
                return {
                    str(k): json_safe(
                        v,
                        _depth=_depth + 1,
                        _seen=_seen,
                        _path=f"{_path}.{k}",
                    )
                    for k, v in value.items()
                }
            if isinstance(value, (list, tuple)):
                return [
                    json_safe(
                        v,
                        _depth=_depth + 1,
                        _seen=_seen,
                        _path=f"{_path}[{index}]",
                    )
                    for index, v in enumerate(value)
                ]
            if isinstance(value, set):
                ordered = sorted(value, key=repr)
                return [
                    json_safe(
                        v,
                        _depth=_depth + 1,
                        _seen=_seen,
                        _path=f"{_path}[{index}]",
                    )
                    for index, v in enumerate(ordered)
                ]
            if isinstance(value, np.ndarray):
                return json_safe(
                    value.tolist(),
                    _depth=_depth + 1,
                    _seen=_seen,
                    _path=_path,
                )
            try:
                converted = asdict(value)
            except RecursionError as exc:
                raise TypeError(
                    "json_safe dataclass recursion "
                    f"at {_path}: type={type(value).__module__}.{type(value).__qualname__}"
                ) from exc
            return json_safe(
                converted,
                _depth=_depth + 1,
                _seen=_seen,
                _path=_path,
            )
        finally:
            _seen.remove(identity)

    # Never iterate an unknown object: repr is finite and makes the loss explicit.
    return repr(value)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".tmp")
    tmp.write_text(json.dumps(json_safe(payload), indent=2, sort_keys=False) + "\n")
    os.replace(tmp, path)


def terms(row: dict[str, Any], key: str) -> list[tuple[int, int]]:
    return [tuple(map(int, term)) for term in (row.get(key) or [])]


def catalogue_label(row: dict[str, Any], index: int) -> str:
    if row.get("code_id"):
        return str(row["code_id"])
    if row.get("bliss_hash"):
        return f"bliss_{row['bliss_hash']}"
    return f"catalogue_row_{index:04d}"


def load_catalogue() -> list[dict[str, Any]]:
    rows = [json.loads(line) for line in CATALOGUE.open() if line.strip()]
    if len(rows) != EXPECTED_CATALOGUE_ROWS:
        raise AssertionError(
            f"catalogue changed: expected {EXPECTED_CATALOGUE_ROWS}, found {len(rows)}"
        )
    required = {
        "ell",
        "m",
        "A_terms",
        "B_terms",
        "C_terms",
        "D_terms",
        "n",
        "k",
        "d",
    }
    for index, row in enumerate(rows):
        missing = required - row.keys()
        if missing:
            raise KeyError(f"catalogue row {index} lacks {sorted(missing)}")
    return rows


def parent_matrices(row: dict[str, Any]) -> tuple[BBSpec, np.ndarray, np.ndarray]:
    spec = BBSpec(
        ell=int(row["ell"]),
        m=int(row["m"]),
        A=terms(row, "A_terms"),
        B=terms(row, "B_terms"),
    )
    HX, HZ = build_bb(spec)
    return spec, HX, HZ


def pbb_code(row: dict[str, Any]):
    spec = PBBSpec(
        ell=int(row["ell"]),
        m=int(row["m"]),
        A=terms(row, "A_terms"),
        B=terms(row, "B_terms"),
        C=terms(row, "C_terms"),
        D=terms(row, "D_terms"),
    )
    return spec, build_pbb(spec, check=True)


def structural_item(row: dict[str, Any], index: int) -> dict[str, Any] | None:
    pspec, code = pbb_code(row)
    theory = analyse_pbb(pspec)
    if theory.delta == 0:
        return None
    parent_spec, HX, HZ = parent_matrices(row)
    rx, rz = rank_np(HX), rank_np(HZ)
    parent_k = int(HX.shape[1] - rx - rz)
    pbb_rank = rank_np(code.H)
    pbb_k_rank = int(code.n - pbb_rank)
    if code.n != int(row["n"]):
        raise AssertionError(f"{catalogue_label(row, index)} n mismatch")
    if pbb_k_rank != int(row["k"]) or theory.k_pbb != int(row["k"]):
        raise AssertionError(f"{catalogue_label(row, index)} PBB k mismatch")
    if parent_k != theory.k_bb or parent_k - int(row["k"]) != theory.delta:
        raise AssertionError(f"{catalogue_label(row, index)} parent k/delta mismatch")
    return {
        "label": catalogue_label(row, index),
        "catalogue_index": index,
        "row": row,
        "n": int(row["n"]),
        "k": int(row["k"]),
        "delta": int(theory.delta),
        "pbb_rank": int(pbb_rank),
        "parent_rank_x": int(rx),
        "parent_rank_z": int(rz),
        "parent_k": parent_k,
        "parent_fingerprint": matrix_fingerprint(HX, HZ),
        "pbb_fingerprint": matrix_fingerprint(code.H),
        "parent_spec": {
            "ell": parent_spec.ell,
            "m": parent_spec.m,
            "A": parent_spec.A,
            "B": parent_spec.B,
        },
    }


def load_certificate_indexes() -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Read the persisted schema and index only fully certified exact records."""

    parent: dict[str, dict[str, Any]] = {}
    pbb: dict[str, dict[str, Any]] = {}
    for path in sorted(CERTIFICATE_DIR.glob("*.json")):
        data = json.loads(path.read_text())
        if not data.get("CERTIFIED_EXACT"):
            continue
        A = data.get("A_terms", data.get("A"))
        B = data.get("B_terms", data.get("B"))
        if data.get("ell") is None or data.get("m") is None or A is None or B is None:
            continue
        base = {
            "ell": int(data["ell"]),
            "m": int(data["m"]),
            "A_terms": A,
            "B_terms": B,
        }
        if "C_terms" in data or "C" in data or "D_terms" in data or "D" in data:
            base.update(
                C_terms=data.get("C_terms", data.get("C", [])),
                D_terms=data.get("D_terms", data.get("D", [])),
            )
            _, code = pbb_code(base)
            fp = matrix_fingerprint(code.H)
            target = pbb
        else:
            _, HX, HZ = parent_matrices(base)
            fp = matrix_fingerprint(HX, HZ)
            target = parent
        d = int(data["d"])
        lb = int(data.get("d_lower_bound", d))
        if lb != d or not data.get("solver_reports_exact", True):
            raise AssertionError(f"inconsistent exact certificate {path}")
        target[fp] = {
            "value": d,
            "lower_bound": d,
            "upper_bound": d,
            "lower_bound_certified": True,
            "upper_bound_certified": True,
            "exact": True,
            "status": "OPTIMAL",
            "bound_direction": "EXACT",
            "provenance": f"certificate:{path.relative_to(ROOT)}",
            "certificate_file": str(path.relative_to(ROOT)),
            "solver_seed": None,
            "search_cap": data.get("search_cap"),
            "wall_time_s": float(data.get("wall_time_s", 0.0)),
            "raw_solver_result": None,
        }
    return parent, pbb


def bravyi_index() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if set(BRAVYI_BB) != set(REFERENCE_DISTANCE):
        raise AssertionError("BRAVYI_BB and registered literature table disagree")
    for name, spec in BRAVYI_BB.items():
        HX, HZ = build_bb(spec)
        reported = int(REFERENCE_DISTANCE[name])
        out[matrix_fingerprint(HX, HZ)] = {
            "value": None,
            "lower_bound": 1,
            "upper_bound": None,
            "lower_bound_certified": True,
            "upper_bound_certified": False,
            "exact": False,
            "status": "UNSETTLED_LITERATURE_REPORTED",
            "bound_direction": "TRIVIAL_LOWER_BOUND",
            "provenance": f"literature_reported_uncertified:BRAVYI_BB:{name}",
            "literature_reported_value": reported,
            "reference_name": name,
            "solver_seed": None,
            "search_cap": None,
            "wall_time_s": 0.0,
            "raw_solver_result": None,
        }
    return out


def verify_pbb_solver_witness(code, result: dict[str, Any]) -> dict[str, Any]:
    vector = result.get("witness_vector")
    if vector is None:
        return {"present": False, "valid": False}
    v = np.asarray(vector, dtype=np.uint8)
    expected = result.get("value")
    return {
        "present": True,
        "is_nontrivial_logical": bool(code.is_logical(v)),
        "weight": int(symplectic_weight(v)),
        "matches_value": expected is not None and symplectic_weight(v) == int(expected),
        "valid": bool(code.is_logical(v))
        and expected is not None
        and symplectic_weight(v) == int(expected),
    }


def verify_css_solver_witnesses(spec: BBSpec, result: dict[str, Any]) -> dict[str, Any]:
    code = bb_stabilizer(spec)
    n = code.n
    checks: dict[str, Any] = {}
    for tag in ("X", "Z"):
        support = result.get(f"d_{tag}_witness")
        if support is None:
            checks[tag] = {"present": False, "valid": False}
            continue
        v = np.zeros(2 * n, dtype=np.uint8)
        if tag == "X":
            v[np.asarray(support, dtype=int)] = 1
        else:
            v[n + np.asarray(support, dtype=int)] = 1
        expected = result.get(f"d_{tag}")
        checks[tag] = {
            "present": True,
            "is_nontrivial_logical": bool(code.is_logical(v)),
            "weight": int(symplectic_weight(v)),
            "matches_value": expected is not None
            and symplectic_weight(v) == int(expected),
            "valid": bool(code.is_logical(v))
            and expected is not None
            and symplectic_weight(v) == int(expected),
        }
    return checks


def normalize_parent_solver_result(
    spec: BBSpec,
    raw: dict[str, Any],
    *,
    seed: int,
    search_cap: int | None,
) -> dict[str, Any]:
    """Turn CSS-sector output into certified bounds for a BB parent.

    Every BB parent has d_X=d_Z by the half-swap/inversion isometry proved in
    Proposition 4.  Therefore a certified lower bound from either sector is a
    lower bound on d=min(d_X,d_Z); taking max rather than the generic CSS min is
    sound and records that structural provenance explicitly.
    """

    witnesses = verify_css_solver_witnesses(spec, raw)
    if not raw.get("d_exact"):
        # Parent lower bounds are admitted only from a complete OPTIMAL exact
        # run.  Even mathematically valid partial CP-SAT bounds are diagnostics
        # under EXP-027's strict provenance ledger.
        return {
            "value": None,
            "lower_bound": 1,
            "upper_bound": None,
            "lower_bound_certified": True,
            "upper_bound_certified": False,
            "exact": False,
            "status": "UNSETTLED_SOLVER_BOUND",
            "bound_direction": "TRIVIAL_LOWER_BOUND",
            "provenance": (
                "exact_distance_css did not return d_exact=true/OPTIMAL in all "
                "sectors; every partial bound and witness is diagnostic only"
            ),
            "solver_seed": seed,
            "search_cap": search_cap,
            "wall_time_s": float(raw.get("wall_time_s", 0.0)),
            "sector_lower_bounds_diagnostic_only": {
                "X": int(raw.get("d_X_lower_bound", 0)),
                "Z": int(raw.get("d_Z_lower_bound", 0)),
            },
            "witness_verification_diagnostic_only": witnesses,
            "raw_solver_result": raw,
        }
    distance = int(raw["d"])
    if int(raw.get("d_lower_bound", -1)) != distance:
        raise AssertionError("exact parent result has inconsistent lower bound")
    relevant_valid = any(
        check.get("valid") and check.get("weight") == distance
        for check in witnesses.values()
    )
    if not relevant_valid:
        raise AssertionError("exact parent distance lacks a valid logical witness")
    return {
        "value": distance,
        "lower_bound": distance,
        "upper_bound": distance,
        "lower_bound_certified": True,
        "upper_bound_certified": True,
        "exact": True,
        "status": "OPTIMAL",
        "bound_direction": "EXACT",
        "provenance": (
            "fresh exact_distance_css d_exact=true/OPTIMAL in all sectors + "
            "independently verified minimum-weight logical witness"
        ),
        "solver_seed": seed,
        "search_cap": search_cap,
        "wall_time_s": float(raw.get("wall_time_s", 0.0)),
        "sector_lower_bounds": {
            "X": int(raw["d_X_lower_bound"]),
            "Z": int(raw["d_Z_lower_bound"]),
        },
        "witness_verification": witnesses,
        "raw_solver_result": raw,
    }


def normalize_pbb_solver_result(
    row: dict[str, Any], raw: dict[str, Any], *, seed: int, search_cap: int
) -> dict[str, Any]:
    _, code = pbb_code(row)
    witness = verify_pbb_solver_witness(code, raw)
    solver_upper = int(raw["value"]) if raw.get("value") is not None else None
    catalogue_upper = int(row["d"])
    if solver_upper is not None and not witness["valid"]:
        raise AssertionError("PBB solver upper bound lacks a valid logical witness")
    if not raw.get("exact"):
        # A timeout is diagnostic only.  Keep the independently supplied
        # catalogue upper bound, but do not use CP-SAT's partial lower bound.
        result = catalogue_pbb_bound(row)
        result.update(
            status=f"UNSETTLED_{raw.get('status', 'PARTIAL')}_CATALOGUE_UPPER",
            provenance=(
                "exact_distance_symplectic did not certify every sector; "
                "partial solver evidence excluded; catalogue d retained only "
                "as the protocol-mandated upper bound"
            ),
            solver_seed=seed,
            search_cap=search_cap,
            wall_time_s=float(raw.get("wall_time_s", 0.0)),
            witness_verification_diagnostic_only=witness,
            raw_solver_result=raw,
        )
        return result
    lower = int(raw["lower_bound"])
    upper = min(solver_upper, catalogue_upper)
    if lower > upper:
        raise AssertionError(f"PBB solver contradicts catalogue upper bound: {lower}>{upper}")
    exact = bool(solver_upper == lower == upper and witness["valid"])
    if not exact:
        raise AssertionError("exact PBB solver result failed independent exactness checks")
    return {
        "value": upper,
        "lower_bound": lower,
        "upper_bound": upper,
        "lower_bound_certified": True,
        "upper_bound_certified": True,
        "exact": True,
        "status": "OPTIMAL",
        "bound_direction": "EXACT",
        "provenance": (
            "exact_distance_symplectic + independently verified witness; "
            "catalogue d retained as an upper-bound cross-check"
        ),
        "solver_seed": seed,
        "search_cap": search_cap,
        "wall_time_s": float(raw.get("wall_time_s", 0.0)),
        "witness_verification": witness,
        "raw_solver_result": raw,
    }


def catalogue_pbb_bound(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "value": int(row["d"]),
        "lower_bound": 1,
        "upper_bound": int(row["d"]),
        "lower_bound_certified": True,
        "upper_bound_certified": True,
        "exact": False,
        "status": "CATALOGUE_UPPER_BOUND",
        "bound_direction": "UPPER_BOUND",
        "provenance": (
            "catalogue:d (forced to UPPER_BOUND; catalogue d_is_exact ignored by protocol)"
        ),
        "catalogue_d_is_exact_ignored": bool(row.get("d_is_exact", False)),
        "catalogue_d_method": row.get("d_method"),
        "solver_seed": None,
        "search_cap": None,
        "wall_time_s": 0.0,
        "raw_solver_result": None,
    }


def unset_parent_bound() -> dict[str, Any]:
    return {
        "value": None,
        "lower_bound": 1,
        "upper_bound": None,
        "lower_bound_certified": True,
        "upper_bound_certified": False,
        "exact": False,
        "status": "UNSETTLED_SIZE_LIMIT",
        "bound_direction": "TRIVIAL_LOWER_BOUND",
        "provenance": f"n > {PARENT_EXACT_MAX_N} and no certified local distance",
        "solver_seed": None,
        "search_cap": None,
        "wall_time_s": 0.0,
        "raw_solver_result": None,
    }


def solve_pbb(row: dict[str, Any], fingerprint: str) -> dict[str, Any]:
    seed = derived_seed(f"exp027:primary:pbb:{fingerprint}")
    _, code = pbb_code(row)
    cap = int(row["d"])
    with cp_sat_seed(seed):
        raw = exact_distance_symplectic(
            code,
            time_limit_s=TIME_LIMIT_OVERRIDE_S,
            workers=SOLVER_WORKERS,
            upper_bound=cap,
        ).to_dict()
    return normalize_pbb_solver_result(row, raw, seed=seed, search_cap=cap)


def solve_parent_job(request: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    row = request["row"]
    fingerprint = request["fingerprint"]
    cap = request.get("search_cap")
    seed = int(request["seed"])
    spec, HX, HZ = parent_matrices(row)
    if matrix_fingerprint(HX, HZ) != fingerprint:
        raise AssertionError("parent fingerprint changed in worker")
    with cp_sat_seed(seed):
        raw = exact_distance_css(
            HX,
            HZ,
            time_limit_s=TIME_LIMIT_OVERRIDE_S,
            workers=SOLVER_WORKERS,
            upper_bound=None,
        )
    return fingerprint, normalize_parent_solver_result(
        spec, raw, seed=seed, search_cap=cap
    )


def checkpoint_payload(
    base: dict[str, Any], prior_wall_s: float, run_start: float, verdict: str
) -> dict[str, Any]:
    payload = dict(base)
    payload.update(
        protocol=PROTOCOL,
        seeds=SEEDS,
        protocol_fingerprint=protocol_fingerprint(),
        machine_readable_verdict=verdict,
        verdict=verdict,
        wall_s=prior_wall_s + time.monotonic() - run_start,
        updated_at=utc_now(),
    )
    return payload


def load_checkpoint() -> tuple[dict[str, Any], float]:
    if not CHECKPOINT_PATH.exists():
        return {
            "started_at": utc_now(),
            "pbb_results": {},
            "parent_results": {},
            "reversal_verifications": {},
            "per_item_results": [],
        }, 0.0
    checkpoint = json.loads(CHECKPOINT_PATH.read_text())
    if checkpoint.get("protocol_fingerprint") != protocol_fingerprint():
        # Tightened provenance/classification rules do not invalidate the raw,
        # matrix-keyed solver transcript.  Re-normalize it below before reuse.
        migrated = {
            "started_at": checkpoint.get("started_at", utc_now()),
            "pbb_results": {},
            "parent_results": {},
            "reversal_verifications": {},
            "per_item_results": [],
        }
        for fp, old in checkpoint.get("pbb_results", {}).items():
            raw = old.get("raw_solver_result")
            if raw is not None:
                migrated["pbb_results"][fp] = {
                    "_raw_for_renormalization": raw,
                    "solver_seed": old.get("solver_seed"),
                    "search_cap": old.get("search_cap"),
                }
        for fp, old in checkpoint.get("parent_results", {}).items():
            raw = old.get("raw_solver_result")
            if raw is not None:
                migrated["parent_results"][fp] = {
                    "_raw_for_renormalization": raw,
                    "solver_seed": old.get("solver_seed"),
                    "search_cap": old.get("search_cap"),
                }
        return migrated, float(checkpoint.get("wall_s", 0.0))
    return checkpoint, float(checkpoint.get("wall_s", 0.0))


def save_checkpoint(
    checkpoint: dict[str, Any], prior_wall_s: float, run_start: float, verdict: str
) -> None:
    # A live/interrupted computation is clean but incomplete evidence.  Route
    # resumable state away from results/raw/exp027_* until full coverage passes.
    route = canonical_route(clean=True, full_coverage=False)
    if route != "partial_runs":
        raise AssertionError("partial checkpoint was routed canonically")
    path = ROOT / "results" / route / "exp027_checkpoint.json"
    atomic_json(path, checkpoint_payload(checkpoint, prior_wall_s, run_start, verdict))


def reverify_reversal(
    item: dict[str, Any],
    primary_parent: dict[str, Any],
    primary_pbb: dict[str, Any],
) -> dict[str, Any]:
    """Independently rebuild and rerun a strict bound-separated reversal."""

    started = time.monotonic()
    row = item["row"]
    pspec, code = pbb_code(row)
    parent_spec, HX, HZ = parent_matrices(row)
    theory = analyse_pbb(pspec)
    parent_k = int(HX.shape[1] - rank_np(HX) - rank_np(HZ))
    pbb_k = int(code.n - rank_np(code.H))
    rank_checks = {
        "pbb_k_from_rebuilt_rank": pbb_k,
        "pbb_k_catalogue": int(row["k"]),
        "parent_k_from_rebuilt_rank": parent_k,
        "parent_k_primary": int(item["parent_k"]),
        "delta_from_rebuilt_rank": parent_k - pbb_k,
        "delta_theory": int(theory.delta),
        "all_match": bool(
            pbb_k == int(row["k"])
            and parent_k == int(item["parent_k"])
            and parent_k - pbb_k == int(theory.delta) == int(item["delta"])
        ),
    }
    rebuilt_hashes = {
        "pbb": matrix_fingerprint(code.H),
        "parent": matrix_fingerprint(HX, HZ),
        "match_primary": bool(
            matrix_fingerprint(code.H) == item["pbb_fingerprint"]
            and matrix_fingerprint(HX, HZ) == item["parent_fingerprint"]
        ),
    }

    parent_seed = derived_seed(
        f"exp027:verification:parent:{item['parent_fingerprint']}"
    )
    pbb_seed = derived_seed(f"exp027:verification:pbb:{item['pbb_fingerprint']}")
    parent_cap = int(primary_pbb["upper_bound"])
    pbb_cap = int(primary_pbb["upper_bound"])
    with cp_sat_seed(parent_seed):
        raw_parent = exact_distance_css(
            HX,
            HZ,
            time_limit_s=TIME_LIMIT_S,
            workers=SOLVER_WORKERS,
            upper_bound=parent_cap,
        )
    secondary_parent = normalize_parent_solver_result(
        parent_spec, raw_parent, seed=parent_seed, search_cap=parent_cap
    )
    with cp_sat_seed(pbb_seed):
        raw_pbb = exact_distance_symplectic(
            code,
            time_limit_s=TIME_LIMIT_S,
            workers=SOLVER_WORKERS,
            upper_bound=pbb_cap,
        ).to_dict()
    secondary_pbb = normalize_pbb_solver_result(
        row, raw_pbb, seed=pbb_seed, search_cap=pbb_cap
    )

    seeds_distinct = bool(
        parent_seed != primary_parent.get("solver_seed")
        and pbb_seed != primary_pbb.get("solver_seed")
    )
    repeated_separation = bool(
        secondary_parent["exact"]
        and secondary_pbb["exact"]
        and int(secondary_pbb["lower_bound"])
        > int(secondary_parent["upper_bound"])
        and secondary_parent["upper_bound"] == primary_parent["upper_bound"]
        and secondary_pbb["lower_bound"] == primary_pbb["lower_bound"]
    )
    double_verified = bool(
        rank_checks["all_match"]
        and rebuilt_hashes["match_primary"]
        and seeds_distinct
        and repeated_separation
    )
    result = {
        "label": item["label"],
        "n": item["n"],
        "k": item["k"],
        "delta": item["delta"],
        "terms": {
            "ell": int(row["ell"]),
            "m": int(row["m"]),
            "A_terms": row["A_terms"],
            "B_terms": row["B_terms"],
            "C_terms": row["C_terms"],
            "D_terms": row["D_terms"],
        },
        "rank_checks": rank_checks,
        "rebuilt_hashes": rebuilt_hashes,
        "primary_parent_distance": primary_parent,
        "primary_pbb_distance": primary_pbb,
        "secondary_parent_distance": secondary_parent,
        "secondary_pbb_distance": secondary_pbb,
        "seeds": {
            "base_seed": BASE_SEED,
            "primary_parent_seed": primary_parent.get("solver_seed"),
            "primary_pbb_seed": primary_pbb.get("solver_seed"),
            "verification_parent_seed": parent_seed,
            "verification_pbb_seed": pbb_seed,
            "primary_and_verification_distinct": seeds_distinct,
        },
        "double_verified": double_verified,
        "wall_s": time.monotonic() - started,
    }
    certificate_item = dict(result)
    result["certificate_payload"] = {
        "protocol": PROTOCOL,
        "seeds": result["seeds"],
        "certificate_type": "EXP027_DOUBLE_VERIFIED_DISTANCE_REVERSAL",
        "item": certificate_item,
        "summary": {
            "items": 1,
            "double_verified_reversals": int(double_verified),
        },
        "machine_readable_verdict": (
            "BREAKTHROUGH_CANDIDATE" if double_verified else "INCONCLUSIVE"
        ),
        "verdict": "BREAKTHROUGH_CANDIDATE" if double_verified else "INCONCLUSIVE",
        "wall_s": result["wall_s"],
        "created_at": utc_now(),
    }
    # Persist only after the full 155-row run reaches its canonical-write gate.
    # Until then the verification (including its certificate payload) lives in
    # the partial checkpoint and cannot be cited as a canonical experiment.
    return result


def flat_row(
    item: dict[str, Any],
    pbb: dict[str, Any],
    parent: dict[str, Any],
    verdict: str,
    verification: dict[str, Any] | None,
) -> dict[str, Any]:
    return {
        "label": item["label"],
        "catalogue_index": item["catalogue_index"],
        "n": item["n"],
        "k": item["k"],
        "delta": item["delta"],
        "pbb_d": pbb["value"],
        "pbb_d_lower_bound": pbb["lower_bound"],
        "pbb_d_upper_bound": pbb["upper_bound"],
        "pbb_d_bound_direction": pbb["bound_direction"],
        "pbb_d_provenance": pbb["provenance"],
        "pbb_d_status": pbb["status"],
        "pbb_d_lower_bound_certified": pbb["lower_bound_certified"],
        "pbb_d_upper_bound_certified": pbb["upper_bound_certified"],
        "pbb_d_exact": pbb["exact"],
        "parent_k": item["parent_k"],
        "parent_d": parent["value"],
        "parent_d_lower_bound": parent["lower_bound"],
        "parent_d_upper_bound": parent["upper_bound"],
        "parent_d_bound_direction": parent["bound_direction"],
        "parent_d_provenance": parent["provenance"],
        "parent_d_status": parent["status"],
        "parent_d_lower_bound_certified": parent["lower_bound_certified"],
        "parent_literature_reported_value": parent.get("literature_reported_value"),
        "parent_d_upper_bound_certified": parent["upper_bound_certified"],
        "parent_d_exact": parent["exact"],
        "rank_verification": {
            "pbb_rank": item["pbb_rank"],
            "pbb_k_from_rank": item["n"] - item["pbb_rank"],
            "parent_rank_x": item["parent_rank_x"],
            "parent_rank_z": item["parent_rank_z"],
            "parent_k_from_rank_formula": (
                item["n"] - item["parent_rank_x"] - item["parent_rank_z"]
            ),
            "k_relation_holds": item["parent_k"] - item["k"] == item["delta"],
        },
        "verdict": verdict,
        "double_verification": (
            {
                "performed": True,
                "double_verified": bool(verification["double_verified"]),
                "raw_certificate": (
                    "results/raw/exp027_reversal_"
                    + re.sub(r"[^A-Za-z0-9_.-]+", "_", item["label"])
                    + ".json"
                ),
            }
            if verification is not None
            else {"performed": False, "double_verified": False}
        ),
    }


def make_report(payload: dict[str, Any]) -> str:
    summary = payload["summary"]
    reversals = [r for r in payload["per_item_results"] if r["verdict"] == "CERTIFIED_REVERSAL"]
    candidates = [
        r
        for r in payload["per_item_results"]
        if r["verdict"] == "UNDECIDED_PARENT_BELOW_PBB_UPPER_BOUND"
    ]
    if reversals:
        first = (
            "BREAKTHROUGH_CANDIDATE — "
            f"{len(reversals)} delta>0 catalogue codes have double-verified "
            "certified d(PBB)>d(parent), refuting universal parent-distance domination."
        )
    else:
        first = (
            "NEGATIVE — Zero certified reversals were found; "
            f"{summary['domination_proved']} were proved dominated and "
            f"{summary['undecided']} remain undecided, consistent with (but not proving) universality."
        )

    lines = [
        first,
        "",
        "# EXP-027 delta>0 parent-domination audit",
        "",
        "## Certified outcome",
        "",
        f"All **{summary['delta_positive_rows']}** catalogue rows with `delta>0` were rebuilt and rank-checked.",
        "The mutually exclusive result counts are:",
        "",
        f"- parent domination proved: **{summary['domination_proved']}**;",
        f"- certified, independently rerun reversals: **{summary['certified_reversals']}**;",
        f"- undecided: **{summary['undecided']}**, including **{summary['uncertified_parent_below_pbb_upper']}** cases where a certified parent upper bound is below only the PBB upper bound.",
        f"- BRAVYI_BB provenance correction: **{summary['moved_from_proved_under_old_bravyi_label_rule']}** rows that would have been called proved if name-label values were incorrectly promoted are now UNSETTLED.",
        "",
        "`DOMINATION_PROVED` requires certified `parent_d_lower_bound >= pbb_d_upper_bound`; catalogue `d` was always used only in the upper-bound direction. `CERTIFIED_REVERSAL` requires certified `pbb_d_lower_bound > parent_d_upper_bound`, followed by a raw-term rebuild, fresh rank checks, and exact solver reruns with different seeds.",
        "",
    ]
    if reversals:
        lines.extend([
            "### Double-verified reversals",
            "",
            "| label | PBB | parent | k change | raw certificate |",
            "|---|---:|---:|---:|---|",
        ])
        for row in reversals:
            lines.append(
                f"| `{row['label']}` | exact(OPTIMAL) d={row['pbb_d']} | "
                f"exact(OPTIMAL) d={row['parent_d']} | {row['k']} -> {row['parent_k']} "
                f"(delta={row['delta']}) | `{row['double_verification']['raw_certificate']}` |"
            )
        lines.append("")
    if candidates:
        lines.extend([
            "### Uncertified parent-below-catalogue-upper cases",
            "",
            "These are **not** reversals: the catalogue PBB distance is only an upper bound.",
            "",
            "| label | parent upper | PBB catalogue upper | PBB certified lower |",
            "|---|---:|---:|---:|",
        ])
        for row in candidates:
            lines.append(
                f"| `{row['label']}` | {row['parent_d_upper_bound']} | "
                f"{row['pbb_d_upper_bound']} | {row['pbb_d_lower_bound']} |"
            )
        lines.append("")

    lines.extend([
        "## README/proof count changes required (not applied)",
        "",
        "No source narrative was edited by EXP-027. The exact integration changes are:",
        "",
        "1. `README.md` lines 27–28: replace the older **77 tested / 49 parent-dominates / 5 parent-weaker / 23 undecided** EXP-012 sentence with the full EXP-027 counts above. The old `parent weaker` category used only a parent witness against a PBB upper bound and must not be described as a certified reversal.",
        "2. `README.md` line 13 and lines 23–28: keep **155 delta>0** and **213/368 delta=0** unchanged, but update the universal-parent-domination status to distinguish proved domination, certified reversals, and budget-unsettled rows.",
        "3. `proofs/pbb_structure.md` lines 167–182: keep Corollary 1's **213/368** count unchanged; append the EXP-027 delta>0 audit counts and, if present above, the double-verified own-parent reversals. Lines 172–177 currently discuss only escape from each code's CSS shadow and therefore do not yet state the own-parent distance result.",
        "4. `proofs/pbb_structure.md` line 182: replace the unquantified `partially open for delta>0` ending with the exact proved/undecided/reversal split, while preserving the certified-bound caveat.",
        "",
        "## Protocol and provenance",
        "",
        f"- Parent priority: local OPTIMAL certificate, then `BRAVYI_BB` matrix match recorded as `literature_reported_uncertified` (never a lower bound), then fresh `exact_distance_css` for n <= {PARENT_EXACT_MAX_N}, otherwise unset.",
        f"- PBB priority: certificate, then `exact_distance_symplectic` for n <= {PBB_EXACT_MAX_N}, otherwise catalogue upper bound.",
        f"- Every CP-SAT subproblem used a **{TIME_LIMIT_S:.0f} s** limit and **{SOLVER_WORKERS} workers**. Two parent calls ran concurrently, for at most **{CONCURRENT_PARENT_JOBS * SOLVER_WORKERS}** solver workers.",
        f"- Seeds derive deterministically from **{BASE_SEED}** and matrix fingerprints. Reversal reruns use a distinct namespace and distinct seeds.",
        f"- Experiment wall time: **{payload['wall_s']:.2f} s**.",
        f"- Machine: `{payload['environment']['platform']}`; Python `{payload['environment']['python']}`.",
        f"- Shared-load observation: start `{payload['environment']['load_average_start']}`, end `{payload['environment']['load_average_end']}`. Solver wall times are not latency benchmarks.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd /Users/jinleic/jinleic-workspace/qec-codesign",
        "PYTHONPATH=src .venv/bin/python experiments/exp027_delta_audit.py",
        "```",
        "",
        "The partial checkpoint is resumable only when its protocol fingerprint matches. Re-running a completed experiment reuses certified cached solver records and regenerates the processed table/report.",
        "",
        "## Caveats",
        "",
        "- The catalogue enumeration is complete for the 368-row pinned file, but timeout-limited distance searches are not exhaustive proofs for undecided rows.",
        "- `time_limit_s` in the exact-distance API is per logical-sector CP-SAT problem, not a whole-code wall-clock limit.",
        "- Any parent call not returning complete `d_exact=true`/OPTIMAL is `UNSETTLED`; partial CP-SAT bounds are retained only as diagnostics and never affect classification.",
        "- `BRAVYI_BB` numeric name/comment labels are `literature_reported_uncertified`. They never supply a certified parent lower bound; the four matching local JSON certificates take precedence.",
        "- No Monte Carlo sampling is involved, so shots, failures, and Clopper–Pearson intervals are not applicable.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    run_start = time.monotonic()
    load_start = os.getloadavg()
    if not FULL_DECLARED_PROTOCOL:
        smoke = {
            "experiment": EXPERIMENT,
            "protocol": PROTOCOL,
            "requested_time_limit_s": TIME_LIMIT_OVERRIDE_S,
            "full_coverage": False,
            "clean": True,
            "machine_readable_verdict": "INCONCLUSIVE_PARTIAL",
            "verdict": "INCONCLUSIVE_PARTIAL",
            "_NOT_CANONICAL": (
                "non-default protocol override; canonical artifacts left untouched"
            ),
            "wall_s": time.monotonic() - run_start,
        }
        route = canonical_route(clean=True, full_coverage=False)
        if route != "partial_runs":
            raise AssertionError("partial smoke run was routed canonically")
        atomic_json(
            ROOT / "results" / route / "exp027_protocol_smoke.json",
            smoke,
        )
        print("PARTIAL protocol smoke -> results/partial_runs/exp027_protocol_smoke.json")
        return
    checkpoint, prior_wall_s = load_checkpoint()
    parent_certificates, pbb_certificates = load_certificate_indexes()
    references = bravyi_index()

    raw_rows = load_catalogue()
    items = [
        item
        for index, row in enumerate(raw_rows)
        if (item := structural_item(row, index)) is not None
    ]
    if len(items) != EXPECTED_DELTA_POSITIVE_ROWS:
        raise AssertionError(
            f"expected {EXPECTED_DELTA_POSITIVE_ROWS} delta>0 rows, found {len(items)}"
        )
    labels = [item["label"] for item in items]
    if len(labels) != len(set(labels)):
        raise AssertionError("catalogue fallback labels are not unique")
    print(f"EXP-027: rebuilt {len(items)} delta>0 rows; all ranks agree", flush=True)

    # Resolve the six small PBB distances before choosing shared-parent caps.
    for position, item in enumerate(items, 1):
        fp = item["pbb_fingerprint"]
        row = item["row"]
        if fp in pbb_certificates:
            result = pbb_certificates[fp]
        elif item["n"] <= PBB_EXACT_MAX_N:
            result = checkpoint["pbb_results"].get(fp)
            if result is not None and "_raw_for_renormalization" in result:
                result = normalize_pbb_solver_result(
                    row,
                    result["_raw_for_renormalization"],
                    seed=int(result["solver_seed"]),
                    search_cap=int(result["search_cap"]),
                )
                checkpoint["pbb_results"][fp] = result
            if result is None:
                print(
                    f"  PBB exact [{position}/{len(items)}] {item['label']} n={item['n']}",
                    flush=True,
                )
                result = solve_pbb(row, fp)
                checkpoint["pbb_results"][fp] = result
                save_checkpoint(
                    checkpoint, prior_wall_s, run_start, "INCONCLUSIVE_RUNNING"
                )
        else:
            result = catalogue_pbb_bound(row)
        item["pbb_distance"] = result

    # One full exact optimization per distinct parent matrix.
    parent_groups: dict[str, list[dict[str, Any]]] = {}
    for item in items:
        parent_groups.setdefault(item["parent_fingerprint"], []).append(item)
    parent_results: dict[str, dict[str, Any]] = {}
    requests: list[dict[str, Any]] = []
    for fp, group in parent_groups.items():
        representative = group[0]
        if fp in parent_certificates:
            parent_results[fp] = parent_certificates[fp]
        elif representative["n"] <= PARENT_EXACT_MAX_N:
            cached = checkpoint["parent_results"].get(fp)
            if cached is not None and "_raw_for_renormalization" in cached:
                spec, _, _ = parent_matrices(representative["row"])
                cached = normalize_parent_solver_result(
                    spec,
                    cached["_raw_for_renormalization"],
                    seed=int(cached["solver_seed"]),
                    search_cap=cached.get("search_cap"),
                )
                checkpoint["parent_results"][fp] = cached
            if cached is not None and cached.get("exact"):
                parent_results[fp] = cached
            else:
                requests.append(
                    {
                        "row": representative["row"],
                        "fingerprint": fp,
                        "search_cap": None,
                        "seed": derived_seed(f"exp027:primary:parent:{fp}"),
                    }
                )
        elif fp in references:
            # Priority (ii) is recorded, but its literature/name label is not a
            # local certificate and therefore supplies no usable lower bound.
            parent_results[fp] = references[fp]
        else:
            parent_results[fp] = unset_parent_bound()

    print(
        f"  unique parents={len(parent_groups)}; queued exact CSS calls={len(requests)}; "
        f"cached/known/unset={len(parent_groups)-len(requests)}",
        flush=True,
    )
    if requests:
        with ProcessPoolExecutor(max_workers=CONCURRENT_PARENT_JOBS) as executor:
            futures = {executor.submit(solve_parent_job, request): request for request in requests}
            for completed, future in enumerate(as_completed(futures), 1):
                request = futures[future]
                fp, result = future.result()
                parent_results[fp] = result
                checkpoint["parent_results"][fp] = result
                print(
                    f"  parent [{completed}/{len(requests)}] "
                    f"status={result['status']} bounds=[{result['lower_bound']},"
                    f"{result['upper_bound']}] cap={request['search_cap']}",
                    flush=True,
                )
                save_checkpoint(
                    checkpoint, prior_wall_s, run_start, "INCONCLUSIVE_RUNNING"
                )

    # Detect strict certified separations first, then independently rerun them.
    strict_apparent: list[dict[str, Any]] = []
    for item in items:
        parent = parent_results[item["parent_fingerprint"]]
        pbb = item["pbb_distance"]
        if (
            parent["upper_bound_certified"]
            and pbb["lower_bound_certified"]
            and parent["upper_bound"] is not None
            and int(pbb["lower_bound"]) > int(parent["upper_bound"])
        ):
            strict_apparent.append(item)

    print(f"  strict certified-bound reversal candidates={len(strict_apparent)}", flush=True)
    verifications: dict[str, dict[str, Any]] = checkpoint["reversal_verifications"]
    for item in strict_apparent:
        cached = verifications.get(item["pbb_fingerprint"])
        if cached is None or not cached.get("double_verified"):
            print(f"  double-verifying reversal {item['label']}", flush=True)
            cached = reverify_reversal(
                item,
                parent_results[item["parent_fingerprint"]],
                item["pbb_distance"],
            )
            verifications[item["pbb_fingerprint"]] = cached
            save_checkpoint(
                checkpoint, prior_wall_s, run_start, "INCONCLUSIVE_RUNNING"
            )

    final_rows: list[dict[str, Any]] = []
    for item in items:
        pbb = item["pbb_distance"]
        parent = parent_results[item["parent_fingerprint"]]
        verification = verifications.get(item["pbb_fingerprint"])
        if (
            parent["lower_bound_certified"]
            and pbb["upper_bound_certified"]
            and int(parent["lower_bound"]) >= int(pbb["upper_bound"])
        ):
            verdict = "DOMINATION_PROVED"
            verification = None
        elif (
            parent["upper_bound_certified"]
            and pbb["lower_bound_certified"]
            and parent["upper_bound"] is not None
            and int(pbb["lower_bound"]) > int(parent["upper_bound"])
            and verification is not None
            and verification.get("double_verified")
        ):
            verdict = "CERTIFIED_REVERSAL"
        elif (
            parent["upper_bound_certified"]
            and parent["upper_bound"] is not None
            and int(parent["upper_bound"]) < int(pbb["upper_bound"])
        ):
            verdict = "UNDECIDED_PARENT_BELOW_PBB_UPPER_BOUND"
        else:
            verdict = "UNDECIDED"
            verification = None
        final_rows.append(flat_row(item, pbb, parent, verdict, verification))

    counts = Counter(row["verdict"] for row in final_rows)
    certified_reversals = counts["CERTIFIED_REVERSAL"]
    machine_verdict = "BREAKTHROUGH_CANDIDATE" if certified_reversals else "NEGATIVE"
    wall_s = prior_wall_s + time.monotonic() - run_start
    summary = {
        "catalogue_rows": len(raw_rows),
        "delta_zero_rows_excluded_by_corollary_1": len(raw_rows) - len(items),
        "delta_positive_rows": len(items),
        "domination_proved": counts["DOMINATION_PROVED"],
        "certified_reversals": certified_reversals,
        "moved_from_proved_under_old_bravyi_label_rule": sum(
            1
            for row in final_rows
            if row["parent_d_status"] == "UNSETTLED_LITERATURE_REPORTED"
            and row["parent_literature_reported_value"] is not None
            and int(row["parent_literature_reported_value"])
            >= int(row["pbb_d_upper_bound"])
        ),
        "undecided": counts["UNDECIDED"]
        + counts["UNDECIDED_PARENT_BELOW_PBB_UPPER_BOUND"],
        "uncertified_parent_below_pbb_upper": counts[
            "UNDECIDED_PARENT_BELOW_PBB_UPPER_BOUND"
        ],
        "verdict_counts": dict(sorted(counts.items())),
        "by_n": {
            str(n): dict(sorted(Counter(
                row["verdict"] for row in final_rows if row["n"] == n
            ).items()))
            for n in sorted({row["n"] for row in final_rows})
        },
        "by_delta": dict(sorted(Counter(row["delta"] for row in final_rows).items())),
        "parent_status_counts": dict(
            sorted(Counter(row["parent_d_status"] for row in final_rows).items())
        ),
        "pbb_status_counts": dict(
            sorted(Counter(row["pbb_d_status"] for row in final_rows).items())
        ),
    }
    if (
        summary["domination_proved"]
        + summary["certified_reversals"]
        + summary["undecided"]
        != len(items)
    ):
        raise AssertionError("final verdict partition does not cover all rows")

    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "load_average_start": list(load_start),
        "load_average_end": list(os.getloadavg()),
        "shared_machine_caveat": (
            "other agents were active; solver wall times are operational metadata, "
            "not isolated latency measurements"
        ),
    }
    payload = {
        "experiment": EXPERIMENT,
        "protocol": PROTOCOL,
        "seeds": SEEDS,
        "protocol_fingerprint": protocol_fingerprint(),
        "per_item_results": final_rows,
        "summary": summary,
        "machine_readable_verdict": machine_verdict,
        "verdict": machine_verdict,
        "wall_s": wall_s,
        "environment": environment,
        "completed_at": utc_now(),
    }
    full_coverage = len(final_rows) == EXPECTED_DELTA_POSITIVE_ROWS
    clean = (
        full_coverage
        and all(row["rank_verification"]["k_relation_holds"] for row in final_rows)
        and summary["domination_proved"]
        + summary["certified_reversals"]
        + summary["undecided"]
        == EXPECTED_DELTA_POSITIVE_ROWS
    )
    route = canonical_route(clean=clean, full_coverage=full_coverage)
    if route != "canonical":
        failure_path = ROOT / "results" / route / "exp027_failed_full_run.json"
        payload["_NOT_CANONICAL"] = (
            "full-scope integrity checks failed; canonical artifacts left untouched"
        )
        atomic_json(failure_path, payload)
        raise RuntimeError(f"EXP-027 did not pass its canonical-write gate: {failure_path}")
    for verification in verifications.values():
        certificate_payload = verification.get("certificate_payload")
        if not certificate_payload:
            continue
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", verification["label"])
        atomic_json(
            ROOT / "results" / "raw" / f"exp027_reversal_{safe}.json",
            certificate_payload,
        )
    atomic_json(PROCESSED_PATH, payload)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(make_report(payload))

    # Promote the now-complete checkpoint only after every full-scope integrity
    # check passes.  No interrupted or partial run can reach this write.
    checkpoint["per_item_results"] = final_rows
    checkpoint["summary"] = summary
    checkpoint["completed_at"] = payload["completed_at"]
    atomic_json(
        ROOT / "results" / "raw" / "exp027_complete_run.json",
        checkpoint_payload(checkpoint, prior_wall_s, run_start, machine_verdict),
    )
    print(
        f"DONE: proved={summary['domination_proved']} "
        f"reversals={summary['certified_reversals']} "
        f"undecided={summary['undecided']} verdict={machine_verdict} "
        f"wall={wall_s:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
