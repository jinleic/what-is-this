"""EXP-070: exact connected-cluster binding of Liang's n=270 BB rows.

Each [[270,8,20]] constructor is first bound by EXP-069 from Liang et al.'s
Table III twisted torus to the canonical (15,9) rectangle.  The race-free
legacy dist-m4ri connected-cluster engine then exhausts rooted logicals through
weight 18 in both the original and exact block-swapped presentations.  Each run
is repeated.  Odd column degree excludes weight 19 and physical weight-20
witnesses close both CSS sectors by BB duality.

All searches are single-threaded.  The run command refuses to launch unless the
calling process has nice priority at least 10.

Run:
  nice -n 15 uv run python experiments/exp070_n270_connected_cluster.py prepare
  nice -n 15 uv run python experiments/exp070_n270_connected_cluster.py run \
    --target liang270a --presentation original --run-id initial \
    --solver results/partial_runs/exp067_n234_cluster/solver/dist_m4ri_arm64_darwin \
    --solver-source-commit 538d119f6e98b3782415eb78f7ce322a076f8101
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp070", "exp056_odd_distance.py")
E67 = _load("exp067_for_exp070", "exp067_n234_connected_cluster.py")
E69 = _load("exp069_for_exp070", "exp069_n270_frontier.py")

from qec_research.distance.sat_decide import css_logical_bases  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp070-n270-connected-cluster-v1"
RUN_SCHEMA = "exp070-connected-cluster-run-v1"
INPUT_SCHEMA = "exp070-connected-cluster-input-v1"
CAP = 18
DISTANCE = 20
PRESENTATIONS = ("original", "block_swapped")
RUN_IDS = ("initial", "replay")
MIN_NICE_PRIORITY = 10
INPUT_DIR = ROOT / "results" / "partial_runs" / "exp070_n270_cluster" / "inputs"
RUN_DIR = ROOT / "results" / "partial_runs" / "exp070_n270_cluster" / "runs"

TARGETS: dict[str, dict[str, Any]] = {
    key: {
        "name": source["name"],
        "source": source["source"],
        "source_url": source["source_url"],
        "ell": 15,
        "m": 9,
        "expected_k": int(source["expected_k"]),
        "expected_d": int(source["expected_d"]),
        "A": source["local_15x9"]["A"],
        "B": source["local_15x9"]["B"],
        "witness_support": source["witness_support"],
    }
    for key, source in E69.LIANG_TARGETS.items()
}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def target_payload(key: str) -> dict[str, Any]:
    if key not in TARGETS:
        raise KeyError(key)
    return {
        field: value
        for field, value in TARGETS[key].items()
        if field != "witness_support"
    }


def problem_of(key: str) -> dict[str, Any]:
    return E56.build_problem(target_payload(key))


def certificate_path(key: str) -> Path:
    return (
        ROOT
        / "results"
        / "certificates"
        / f"exp070_270_8_20_{key}_distance.json"
    )


def run_path(key: str, presentation: str, run_id: str) -> Path:
    if key not in TARGETS or presentation not in PRESENTATIONS or run_id not in RUN_IDS:
        raise ValueError("unknown EXP-070 run identity")
    return RUN_DIR / f"{key}_{presentation}_{run_id}.json"


def _presentation_matrices(
    problem: dict[str, Any], presentation: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if presentation not in PRESENTATIONS:
        raise ValueError(f"unknown presentation {presentation!r}")
    hx = np.asarray(problem["HX"], dtype=np.uint8)
    hz = np.asarray(problem["HZ"], dtype=np.uint8)
    lx, _ = css_logical_bases(hx, hz)
    if presentation == "block_swapped":
        block = int(problem["block"])
        permutation = np.r_[np.arange(block, 2 * block), np.arange(block)]
        hx, hz, lx = hx[:, permutation], hz[:, permutation], lx[:, permutation]
    return hx, hz, lx


def input_identity(key: str, presentation: str) -> dict[str, Any]:
    problem = problem_of(key)
    hx, hz, lx = _presentation_matrices(problem, presentation)
    return {
        "schema": INPUT_SCHEMA,
        "target": key,
        "presentation": presentation,
        "n": int(problem["n"]),
        "k": int(problem["k"]),
        "HX_sha256": E56.matrix_sha256(hx),
        "HZ_sha256": E56.matrix_sha256(hz),
        "LX_sha256": E56.matrix_sha256(lx),
        "HX_rank": int(rank_np(hx)),
        "HZ_rank": int(rank_np(hz)),
        "LX_rank": int(rank_np(lx)),
    }


def input_paths(key: str, presentation: str) -> dict[str, Path]:
    stem = f"{key}_{presentation}"
    return {
        "HX": INPUT_DIR / f"{stem}_HX.mtx",
        "HZ": INPUT_DIR / f"{stem}_HZ.mtx",
        "LX": INPUT_DIR / f"{stem}_LX.mtx",
    }


def _atomic_mmwrite(path: Path, matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as stream:
        mmwrite(
            stream,
            sparse.coo_matrix(matrix.astype(np.int8)),
            field="integer",
        )
    temporary.replace(path)


def prepare_inputs() -> dict[str, Any]:
    entries = []
    for key in TARGETS:
        problem = problem_of(key)
        for presentation in PRESENTATIONS:
            matrices = _presentation_matrices(problem, presentation)
            paths = input_paths(key, presentation)
            for name, matrix in zip(("HX", "HZ", "LX"), matrices):
                _atomic_mmwrite(paths[name], matrix)
            entries.append(
                {
                    **input_identity(key, presentation),
                    "paths": {
                        name: str(path.relative_to(ROOT))
                        for name, path in paths.items()
                    },
                }
            )
    manifest = {"schema": INPUT_SCHEMA, "entries": entries}
    E56.atomic_write_json(INPUT_DIR / "manifest.json", manifest)
    return manifest


def input_file_bindings(key: str, presentation: str) -> dict[str, Any]:
    problem = problem_of(key)
    expected = dict(
        zip(("HX", "HZ", "LX"), _presentation_matrices(problem, presentation))
    )
    bindings: dict[str, Any] = {}
    for name, path in input_paths(key, presentation).items():
        if not path.is_file():
            raise RuntimeError(f"missing EXP-070 input: {path}")
        loaded = np.asarray(mmread(path, spmatrix=True).toarray(), dtype=np.uint8) & 1
        if not np.array_equal(loaded, expected[name]):
            raise RuntimeError(f"EXP-070 input matrix changed: {path}")
        bindings[name] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": file_sha256(path),
        }
    return bindings


def _nice_priority() -> int:
    if not hasattr(os, "getpriority"):
        raise RuntimeError("EXP-070 requires an OS nice-priority query")
    return int(os.getpriority(os.PRIO_PROCESS, 0))


def run_search(
    key: str,
    presentation: str,
    run_id: str,
    solver: Path,
    solver_source_commit: str,
    timeout_s: float,
    *,
    force: bool,
) -> int:
    if _nice_priority() < MIN_NICE_PRIORITY:
        raise RuntimeError(
            f"EXP-070 must run at nice priority >= {MIN_NICE_PRIORITY}"
        )
    build = E67.solver_build_binding()
    if (
        solver_source_commit != E67.DIST_M4RI_COMMIT
        or not solver.is_file()
        or file_sha256(solver) != E67.DIST_M4RI_BINARY_SHA256
        or build["binary_sha256"] != E67.DIST_M4RI_BINARY_SHA256
    ):
        raise RuntimeError("solver does not match the audited legacy build")
    prepare_inputs()
    path = run_path(key, presentation, run_id)
    if path.exists() and not force:
        record = json.loads(path.read_text(encoding="utf-8"))
        validate_run_record(record, key, presentation, run_id)
        print(json.dumps(record["result"], indent=1))
        return 0

    inputs = input_paths(key, presentation)
    invocation = [
        str(solver.resolve()),
        "method=2",
        f"finH={inputs['HX'].resolve()}",
        f"finG={inputs['HZ'].resolve()}",
        "start=0",
        "dmin=1",
        f"wmax={CAP}",
        "smax=0",
        "debug=3",
        f"timeout={float(timeout_s):g}",
    ]
    running = {
        "schema": RUN_SCHEMA,
        "utc": E56.utc_now(),
        "target": key,
        "presentation": presentation,
        "run_id": run_id,
        "status": "RUNNING",
        "input": input_identity(key, presentation),
        "input_files": input_file_bindings(key, presentation),
        "invocation": invocation,
        "solver": {
            "name": "QEC-pages/dist-m4ri legacy STANDALONE method 2",
            "source": "https://github.com/QEC-pages/dist-m4ri",
            "source_commit": solver_source_commit,
            "binary_sha256": file_sha256(solver),
        },
        "protocol": {
            "engine_variant": "legacy_single_thread_do_CC_dist",
            "method": 2,
            "root": 0,
            "dmin": 1,
            "wmax": CAP,
            "smax": 0,
            "timeout_s": float(timeout_s),
            "processes": 1,
            "threads": 1,
            "parent_nice_priority": _nice_priority(),
        },
    }
    E56.atomic_write_json(path, running)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            invocation,
            cwd=solver.parent,
            capture_output=True,
            text=True,
            timeout=float(timeout_s) + 60.0,
            check=False,
        )
        result = E67._parse_result(completed.stdout)
        record = {
            **running,
            "status": "COMPLETE" if completed.returncode == 0 else "FAILED",
            "returncode": int(completed.returncode),
            "wall_time_s": time.perf_counter() - started,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
            "stderr_sha256": hashlib.sha256(completed.stderr.encode()).hexdigest(),
            "result": result,
        }
    except subprocess.TimeoutExpired as exc:
        record = {
            **running,
            "status": "UNDECIDED_BUDGET",
            "wall_time_s": time.perf_counter() - started,
            "stdout": (exc.stdout or "") if isinstance(exc.stdout, str) else "",
            "stderr": (exc.stderr or "") if isinstance(exc.stderr, str) else "",
        }
    E56.atomic_write_json(path, record)
    if record["status"] != "COMPLETE":
        return 1
    validate_run_record(record, key, presentation, run_id)
    print(json.dumps(record["result"], indent=1))
    return 0


def validate_run_record(
    record: dict[str, Any], key: str, presentation: str, run_id: str
) -> None:
    protocol = record.get("protocol", {})
    expected_invocation = [
        record.get("invocation", [None])[0],
        "method=2",
        f"finH={input_paths(key, presentation)['HX'].resolve()}",
        f"finG={input_paths(key, presentation)['HZ'].resolve()}",
        "start=0",
        "dmin=1",
        f"wmax={CAP}",
        "smax=0",
        "debug=3",
        f"timeout={float(protocol.get('timeout_s')):g}",
    ]
    try:
        result = record["result"]
        valid = bool(
            record.get("schema") == RUN_SCHEMA
            and record.get("target") == key
            and record.get("presentation") == presentation
            and record.get("run_id") == run_id
            and record.get("status") == "COMPLETE"
            and int(record.get("returncode", -1)) == 0
            and record.get("input") == input_identity(key, presentation)
            and record.get("input_files") == input_file_bindings(key, presentation)
            and record.get("invocation") == expected_invocation
            and record.get("solver", {}).get("source_commit")
            == E67.DIST_M4RI_COMMIT
            and record.get("solver", {}).get("binary_sha256")
            == E67.DIST_M4RI_BINARY_SHA256
            and protocol.get("engine_variant")
            == "legacy_single_thread_do_CC_dist"
            and protocol.get("method") == 2
            and protocol.get("root") == 0
            and protocol.get("dmin") == 1
            and protocol.get("wmax") == CAP
            and protocol.get("smax") == 0
            and protocol.get("processes") == 1
            and protocol.get("threads") == 1
            and int(protocol.get("parent_nice_priority", -1))
            >= MIN_NICE_PRIORITY
            and result.get("status") == "NO_LOGICAL_THROUGH_CAP"
            and int(result.get("exhaustive_through", -1)) == CAP
            and result.get("found_distance") is None
            and E67._parse_result(record.get("stdout", "")) == result
            and hashlib.sha256(record.get("stdout", "").encode()).hexdigest()
            == record.get("stdout_sha256")
            and hashlib.sha256(record.get("stderr", "").encode()).hexdigest()
            == record.get("stderr_sha256")
        )
    except (KeyError, TypeError, ValueError, RuntimeError):
        valid = False
    if not valid:
        raise RuntimeError(
            f"invalid EXP-070 run record: {key}/{presentation}/{run_id}"
        )


def root_reduction(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    cover = E56.build_orbit_cover(problem)
    if not cover["complete"]:
        raise RuntimeError("translation/pole cover is incomplete")
    return {
        "schema": "exp070-two-root-connected-cluster-reduction-v1",
        "translation_group_size": int(problem["block"]),
        "translations_preserve_code": True,
        "presentations": list(PRESENTATIONS),
        "root_coordinate": 0,
        "connected_component_lemma": (
            "Every Tanner-connected component of a zero-syndrome support has "
            "zero syndrome; a nontrivial word therefore has a nontrivial "
            "component no heavier than itself."
        ),
        "irreducible_subset_lemma": (
            "A minimum logical contains no nonempty proper zero-syndrome subset: "
            "a logical subset is lighter, while removing a stabilizer subset "
            "gives a lighter representative of the same logical class."
        ),
        "enumeration_argument": (
            "Every proper partial support is nonzero-syndrome, so the first "
            "unsatisfied check has another target column and the connected-cluster "
            "recursion reaches the complete minimum word."
        ),
        "coverage_argument": (
            "Translation moves any connected logical touching the first block to "
            "root 0. A logical supported only in the second block is covered after "
            "the exact block-column swap."
        ),
        "dist_m4ri_ordering_guard": (
            "start=0 is load-bearing because the engine keeps only clusters whose "
            "root is the minimum column."
        ),
        "valid": True,
    }


def upper_certificate(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    witness = np.zeros(int(problem["n"]), dtype=np.uint8)
    witness[np.asarray(TARGETS[key]["witness_support"], dtype=int)] = 1
    upper = E56.verify_upper_witness(problem, witness)
    if int(upper["weight"]) != DISTANCE:
        raise RuntimeError(f"weight-20 witness failed for {key}")
    return upper


def certificate_identity(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    return {
        "target_key": key,
        "target": target_payload(key),
        "published_transport": E69.published_transport(key),
        "n": int(problem["n"]),
        "k": int(problem["k"]),
        "HX_sha256": E56.matrix_sha256(problem["HX"]),
        "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
    }


def _run_bindings(key: str) -> dict[str, Any]:
    bindings: dict[str, Any] = {}
    identities: set[tuple[str, str]] = set()
    for presentation in PRESENTATIONS:
        bindings[presentation] = {}
        for run_id in RUN_IDS:
            path = run_path(key, presentation, run_id)
            record = json.loads(path.read_text(encoding="utf-8"))
            validate_run_record(record, key, presentation, run_id)
            bindings[presentation][run_id] = {
                "path": str(path.relative_to(ROOT)),
                "sha256": file_sha256(path),
            }
            identities.add(
                (
                    record["solver"]["source_commit"],
                    record["solver"]["binary_sha256"],
                )
            )
    build = E67.solver_build_binding()
    if identities != {(build["source_commit"], build["binary_sha256"])}:
        raise RuntimeError("EXP-070 runs do not match the solver archive")
    return bindings


def assemble(key: str) -> int:
    problem = problem_of(key)
    parity = E56.kernel_weight_parity_proof(problem, DISTANCE - 1)
    duality = E56.bb_duality_proof(problem)
    payload = {
        "schema": SCHEMA,
        "utc": E56.utc_now(),
        "identity": certificate_identity(key),
        "root_reduction": root_reduction(key),
        "solver_build": E67.solver_build_binding(),
        "protocol": {
            "algorithm": "connected-cluster exact enumeration",
            "external_engine": "QEC-pages/dist-m4ri legacy STANDALONE method 2",
            "side": "z",
            "rooted_presentations": list(PRESENTATIONS),
            "exhaustive_through_weight": CAP,
            "replay_required": True,
            "resource_policy": {
                "processes": 1,
                "threads": 1,
                "minimum_nice_priority": MIN_NICE_PRIORITY,
            },
        },
        "parity": parity,
        "duality": {
            field: value
            for field, value in duality.items()
            if field != "permutation"
        },
        "upper": upper_certificate(key),
        "lower_runs": _run_bindings(key),
        "verdict": {
            "classification": "CERTIFIED_EXACT",
            "d": DISTANCE,
            "d_X": DISTANCE,
            "d_Z": DISTANCE,
            "exact": True,
        },
    }
    if not (
        parity["valid"]
        and int(parity["effective_even_cap"]) == CAP
        and duality["valid"]
    ):
        raise RuntimeError("EXP-070 parity/duality closure failed")
    E56.atomic_write_json(certificate_path(key), payload)
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def validate_exact_certificate_payload(payload: dict[str, Any]) -> None:
    key = payload.get("identity", {}).get("target_key")
    if key not in TARGETS:
        raise RuntimeError("unknown EXP-070 target")
    problem = problem_of(key)
    parity = E56.kernel_weight_parity_proof(problem, DISTANCE - 1)
    duality = E56.bb_duality_proof(problem)
    expected_protocol = {
        "algorithm": "connected-cluster exact enumeration",
        "external_engine": "QEC-pages/dist-m4ri legacy STANDALONE method 2",
        "side": "z",
        "rooted_presentations": list(PRESENTATIONS),
        "exhaustive_through_weight": CAP,
        "replay_required": True,
        "resource_policy": {
            "processes": 1,
            "threads": 1,
            "minimum_nice_priority": MIN_NICE_PRIORITY,
        },
    }
    try:
        valid = bool(
            payload.get("schema") == SCHEMA
            and payload.get("identity") == certificate_identity(key)
            and payload.get("root_reduction") == root_reduction(key)
            and payload.get("solver_build") == E67.solver_build_binding()
            and payload.get("protocol") == expected_protocol
            and payload.get("parity") == parity
            and payload.get("duality")
            == {field: value for field, value in duality.items() if field != "permutation"}
            and payload.get("upper") == upper_certificate(key)
            and payload.get("lower_runs") == _run_bindings(key)
            and payload.get("verdict")
            == {
                "classification": "CERTIFIED_EXACT",
                "d": DISTANCE,
                "d_X": DISTANCE,
                "d_Z": DISTANCE,
                "exact": True,
            }
        )
    except (KeyError, TypeError, ValueError, RuntimeError):
        valid = False
    if not valid:
        raise RuntimeError(f"EXP-070 exact certificate failed validation: {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("prepare", "run", "assemble", "validate")
    )
    parser.add_argument("--target", choices=tuple(TARGETS))
    parser.add_argument("--presentation", choices=PRESENTATIONS)
    parser.add_argument("--run-id", choices=RUN_IDS, default="initial")
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--solver-source-commit")
    parser.add_argument("--timeout", type=float, default=7_200.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.command == "prepare":
        print(json.dumps(prepare_inputs(), indent=1))
        return 0
    if args.target is None:
        parser.error(f"{args.command} requires --target")
    if args.command == "run":
        if (
            args.presentation is None
            or args.solver is None
            or args.solver_source_commit is None
        ):
            parser.error(
                "run requires --presentation, --solver, and --solver-source-commit"
            )
        return run_search(
            args.target,
            args.presentation,
            args.run_id,
            args.solver,
            args.solver_source_commit,
            args.timeout,
            force=args.force,
        )
    if args.command == "assemble":
        return assemble(args.target)
    payload = json.loads(certificate_path(args.target).read_text(encoding="utf-8"))
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
