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
from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402

SCHEMA = "exp070-n270-connected-cluster-v1"
RUN_SCHEMA = "exp070-connected-cluster-run-v1"
INPUT_SCHEMA = "exp070-connected-cluster-input-v1"
CAP = 18
BASE_CAP = 16
DISTANCE = 20
PRESENTATIONS = ("original",)
RUN_IDS = ("initial", "replay")
MIN_NICE_PRIORITY = 10
INPUT_DIR = ROOT / "results" / "partial_runs" / "exp070_n270_cluster" / "inputs"
RUN_DIR = ROOT / "results" / "partial_runs" / "exp070_n270_cluster" / "runs"
SOLVER_BUILD = (
    ROOT / "results" / "partial_runs" / "exp070_n270_cluster" / "solver_build.json"
)
PREFIX_BINARY_SHA256 = (
    "12dd8b65ff3d6f3a7ed3609d777fb7364670689606bba000b358be639e61e769"
)
PREFIX_DIST_CC_SHA256 = (
    "14efaa86af45d92ecd8f314a5ab7b322a1ba853d3e3fe9de1697eadfa7ec3cb0"
)
PREFIX_UTIL_IO_C_SHA256 = (
    "412b2a3e63b82b7a64cbdc69cd29c1289a283a509f60e56ddaf79ce283f81e65"
)
PREFIX_UTIL_IO_H_SHA256 = (
    "7cdca39b1deac66af997f77d81879d67bc157c03a8ce4b66ee6413809512a405"
)

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


def solver_build_binding() -> dict[str, Any]:
    build = json.loads(SOLVER_BUILD.read_text(encoding="utf-8"))
    base = SOLVER_BUILD.parent
    paths = {
        "binary": base / build["solver"]["binary"],
        "upstream": base / build["solver"]["upstream_source_archive"],
        "license": base / build["solver"]["license_file"],
        "m4ri_source": base / build["m4ri"]["source_archive"],
        "m4ri_bottle": base / build["m4ri"]["bottle"],
        "m4ri_static": (
            base
            / "solver"
            / "m4ri-bottle"
            / "m4ri"
            / "20260122"
            / "lib"
            / "libm4ri.a"
        ),
        "dist_cc": base / "solver" / "dist-m4ri-src" / "src" / "dist_cc.c",
        "util_io_c": base / "solver" / "dist-m4ri-src" / "src" / "util_io.c",
        "util_io_h": base / "solver" / "dist-m4ri-src" / "src" / "util_io.h",
    }
    valid = bool(
        build.get("schema") == "exp070-prefix-dist-m4ri-build-v1"
        and build.get("solver", {}).get("upstream_commit")
        == E67.DIST_M4RI_COMMIT
        and build.get("solver", {}).get("binary_sha256")
        == PREFIX_BINARY_SHA256
        and build.get("solver", {})
        .get("prefix_partition", {})
        .get("dist_cc_c_sha256")
        == PREFIX_DIST_CC_SHA256
        and build.get("solver", {})
        .get("prefix_partition", {})
        .get("util_io_c_sha256")
        == PREFIX_UTIL_IO_C_SHA256
        and build.get("solver", {})
        .get("prefix_partition", {})
        .get("util_io_h_sha256")
        == PREFIX_UTIL_IO_H_SHA256
        and file_sha256(paths["binary"]) == PREFIX_BINARY_SHA256
        and file_sha256(paths["upstream"]) == E67.DIST_M4RI_SOURCE_SHA256
        and file_sha256(paths["license"]) == E67.DIST_M4RI_LICENSE_SHA256
        and file_sha256(paths["m4ri_source"]) == E67.M4RI_SOURCE_SHA256
        and file_sha256(paths["m4ri_bottle"])
        == build["m4ri"]["bottle_sha256"]
        and file_sha256(paths["m4ri_static"])
        == E67.M4RI_STATIC_LIBRARY_SHA256
        and file_sha256(paths["dist_cc"]) == PREFIX_DIST_CC_SHA256
        and file_sha256(paths["util_io_c"]) == PREFIX_UTIL_IO_C_SHA256
        and file_sha256(paths["util_io_h"]) == PREFIX_UTIL_IO_H_SHA256
    )
    if not valid:
        raise RuntimeError("EXP-070 prefix solver provenance is stale")
    return {
        "path": str(SOLVER_BUILD.relative_to(ROOT)),
        "sha256": file_sha256(SOLVER_BUILD),
        "binary": str(paths["binary"].relative_to(ROOT)),
        "binary_sha256": PREFIX_BINARY_SHA256,
        "source_commit": E67.DIST_M4RI_COMMIT,
        "source_archive_sha256": E67.DIST_M4RI_SOURCE_SHA256,
        "m4ri_source_archive_sha256": E67.M4RI_SOURCE_SHA256,
    }


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


def run_path(
    key: str, mode: str, run_id: str, prefix: int | None = None
) -> Path:
    if key not in TARGETS or mode not in {"base", "prefix"} or run_id not in RUN_IDS:
        raise ValueError("unknown EXP-070 run identity")
    if mode == "base":
        if prefix is not None:
            raise ValueError("base run cannot have a prefix")
        stem = "base"
    else:
        if prefix not in prefix_partition(key)["prefix_columns"]:
            raise ValueError("prefix run is outside the complete branch partition")
        stem = f"prefix{prefix}"
    return RUN_DIR / f"{key}_{stem}_{run_id}.json"


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
    mode: str,
    run_id: str,
    solver: Path,
    solver_source_commit: str,
    timeout_s: float,
    *,
    prefix: int | None,
    force: bool,
) -> int:
    if _nice_priority() < MIN_NICE_PRIORITY:
        raise RuntimeError(
            f"EXP-070 must run at nice priority >= {MIN_NICE_PRIORITY}"
        )
    build = solver_build_binding()
    if (
        solver_source_commit != E67.DIST_M4RI_COMMIT
        or not solver.is_file()
        or solver.resolve() != (ROOT / build["binary"]).resolve()
        or file_sha256(solver) != PREFIX_BINARY_SHA256
    ):
        raise RuntimeError("solver does not match the hash-bound prefix build")
    if mode == "base":
        if prefix is not None:
            raise ValueError("base run cannot have a prefix")
        dmin, wmax = 1, BASE_CAP
    elif mode == "prefix":
        if prefix not in prefix_partition(key)["prefix_columns"]:
            raise ValueError("prefix is outside the complete branch partition")
        dmin = wmax = CAP
    else:
        raise ValueError(f"unknown run mode {mode!r}")

    prepare_inputs()
    path = run_path(key, mode, run_id, prefix)
    if path.exists() and not force:
        record = json.loads(path.read_text(encoding="utf-8"))
        validate_run_record(record, key, mode, run_id, prefix)
        print(json.dumps(record["result"], indent=1))
        return 0

    inputs = input_paths(key, "original")
    invocation = [
        str(solver.resolve()),
        "method=2",
        f"finH={inputs['HX'].resolve()}",
        f"finG={inputs['HZ'].resolve()}",
        "start=0",
    ]
    if prefix is not None:
        invocation.append(f"prefix={prefix}")
    invocation.extend(
        [
            f"dmin={dmin}",
            f"wmax={wmax}",
            "smax=0",
            "debug=0",
            f"timeout={float(timeout_s):g}",
        ]
    )
    running = {
        "schema": RUN_SCHEMA,
        "utc": E56.utc_now(),
        "target": key,
        "mode": mode,
        "prefix": prefix,
        "run_id": run_id,
        "status": "RUNNING",
        "input": input_identity(key, "original"),
        "input_files": input_file_bindings(key, "original"),
        "invocation": invocation,
        "solver": {
            "name": (
                "QEC-pages/dist-m4ri legacy STANDALONE "
                "prefix-partitioned method 2"
            ),
            "source": "https://github.com/QEC-pages/dist-m4ri",
            "source_commit": solver_source_commit,
            "binary_sha256": file_sha256(solver),
            "build": build,
        },
        "protocol": {
            "engine_variant": "legacy_single_thread_do_CC_dist_prefix",
            "method": 2,
            "root": 0,
            "mode": mode,
            "prefix": prefix,
            "dmin": dmin,
            "wmax": wmax,
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
    validate_run_record(record, key, mode, run_id, prefix)
    print(json.dumps(record["result"], indent=1))
    return 0


def validate_run_record(
    record: dict[str, Any],
    key: str,
    mode: str,
    run_id: str,
    prefix: int | None,
) -> None:
    protocol = record.get("protocol", {})
    if mode == "base":
        dmin, wmax = 1, BASE_CAP
    elif mode == "prefix":
        dmin = wmax = CAP
    else:
        raise ValueError(f"unknown run mode {mode!r}")
    expected_invocation = [
        record.get("invocation", [None])[0],
        "method=2",
        f"finH={input_paths(key, 'original')['HX'].resolve()}",
        f"finG={input_paths(key, 'original')['HZ'].resolve()}",
        "start=0",
    ]
    if prefix is not None:
        expected_invocation.append(f"prefix={prefix}")
    expected_invocation.extend(
        [
            f"dmin={dmin}",
            f"wmax={wmax}",
            "smax=0",
            "debug=0",
            f"timeout={float(protocol.get('timeout_s')):g}",
        ]
    )
    try:
        result = record["result"]
        valid = bool(
            record.get("schema") == RUN_SCHEMA
            and record.get("target") == key
            and record.get("mode") == mode
            and record.get("prefix") == prefix
            and record.get("run_id") == run_id
            and record.get("status") == "COMPLETE"
            and int(record.get("returncode", -1)) == 0
            and record.get("input") == input_identity(key, "original")
            and record.get("input_files") == input_file_bindings(key, "original")
            and record.get("invocation") == expected_invocation
            and record.get("solver", {}).get("source_commit")
            == E67.DIST_M4RI_COMMIT
            and record.get("solver", {}).get("binary_sha256")
            == PREFIX_BINARY_SHA256
            and record.get("solver", {}).get("build") == solver_build_binding()
            and protocol.get("engine_variant")
            == "legacy_single_thread_do_CC_dist_prefix"
            and protocol.get("method") == 2
            and protocol.get("root") == 0
            and protocol.get("mode") == mode
            and protocol.get("prefix") == prefix
            and protocol.get("dmin") == dmin
            and protocol.get("wmax") == wmax
            and protocol.get("smax") == 0
            and protocol.get("processes") == 1
            and protocol.get("threads") == 1
            and int(protocol.get("parent_nice_priority", -1))
            >= MIN_NICE_PRIORITY
            and result.get("status") == "NO_LOGICAL_THROUGH_CAP"
            and int(result.get("exhaustive_through", -1)) == wmax
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
            f"invalid EXP-070 run record: {key}/{mode}/{prefix}/{run_id}"
        )




def prefix_partition(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    hx = np.asarray(problem["HX"], dtype=np.uint8)
    root_checks = np.flatnonzero(hx[:, 0])
    if not len(root_checks):
        raise RuntimeError("root column has zero syndrome")
    first_check = int(root_checks[0])
    incident = [int(index) for index in np.flatnonzero(hx[first_check])]
    prefixes = [index for index in incident if index > 0]
    complete = bool(
        incident[0] == 0
        and len(incident) == 6
        and len(prefixes) == 5
        and len(set(prefixes)) == 5
    )
    if not complete:
        raise RuntimeError("EXP-070 root prefix partition is incomplete")
    return {
        "schema": "exp070-root-prefix-partition-v1",
        "root_coordinate": 0,
        "first_unsatisfied_check": first_check,
        "incident_columns": incident,
        "prefix_columns": prefixes,
        "complete": True,
        "coverage_argument": (
            "A nonzero-syndrome partial cluster rooted at 0 must cancel the "
            "first unsatisfied check. Therefore every completed zero-syndrome "
            "cluster contains at least one of the five other incident columns; "
            "the five prefix runs cover the full root recursion. Overlap is safe."
        ),
    }


def pure_second_logical_enumeration(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    block = int(problem["block"])
    kernel = nullspace_np(np.asarray(problem["B"], dtype=np.uint8))
    if len(kernel) > 20:
        raise RuntimeError("pure-second kernel exceeds exhaustive enumeration cap")
    lx, _ = css_logical_bases(problem["HX"], problem["HZ"])
    pairing = np.asarray(lx[:, block:], dtype=np.uint8)
    minimum = block + 1
    minimum_vector: np.ndarray | None = None
    logical_words = 0
    batch = 1 << 15
    for start in range(1, 1 << len(kernel), batch):
        integers = np.arange(
            start, min(start + batch, 1 << len(kernel)), dtype=np.uint32
        )
        coefficients = (
            (integers[:, None] >> np.arange(len(kernel))) & 1
        ).astype(np.uint8)
        vectors = (coefficients @ kernel) % 2
        logical = np.any((vectors @ pairing.T) % 2, axis=1)
        logical_words += int(logical.sum())
        weights = vectors.sum(axis=1)
        weights[~logical] = block + 1
        index = int(weights.argmin())
        if int(weights[index]) < minimum:
            minimum = int(weights[index])
            minimum_vector = vectors[index].copy()
    if minimum_vector is None:
        raise RuntimeError("pure-second logical quotient is unexpectedly empty")
    physical = np.zeros(int(problem["n"]), dtype=np.uint8)
    physical[block:] = minimum_vector
    verification = E56.verify_upper_witness(problem, physical)
    quotient_rank = (
        rank_np(
            np.vstack(
                [
                    problem["HZ"],
                    np.hstack(
                        [
                            np.zeros((len(kernel), block), dtype=np.uint8),
                            kernel,
                        ]
                    ),
                ]
            )
        )
        - rank_np(problem["HZ"])
    )
    return {
        "schema": "exp070-pure-second-exhaustive-v1",
        "kernel_dimension": int(len(kernel)),
        "vectors_enumerated": (1 << len(kernel)) - 1,
        "logical_words": logical_words,
        "logical_quotient_rank": int(quotient_rank),
        "minimum_weight": minimum,
        "minimum_support": verification["support"],
        "kernel_basis_sha256": E56.matrix_sha256(kernel),
        "verification": verification,
        "complete": True,
    }
def root_reduction(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    cover = E56.build_orbit_cover(problem)
    if not cover["complete"]:
        raise RuntimeError("translation/pole cover is incomplete")
    return {
        "schema": "exp070-rooted-prefix-reduction-v1",
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
            "Every proper partial support is nonzero-syndrome. The recursion "
            "therefore reaches each minimum first-block-touching logical."
        ),
        "coverage_argument": (
            "Translation moves every logical touching the first block to root 0; "
            "the exact prefix partition covers that rooted recursion. Logicals "
            "supported only in the second block are exhausted independently."
        ),
        "prefix_partition_sha256": E56.canonical_json_sha256(
            prefix_partition(key)
        ),
        "pure_second_enumeration_sha256": E56.canonical_json_sha256(
            pure_second_logical_enumeration(key)
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
    bindings: dict[str, Any] = {"base": {}, "prefixes": {}}
    identities: set[tuple[str, str]] = set()
    for run_id in RUN_IDS:
        path = run_path(key, "base", run_id)
        record = json.loads(path.read_text(encoding="utf-8"))
        validate_run_record(record, key, "base", run_id, None)
        bindings["base"][run_id] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": file_sha256(path),
        }
        identities.add(
            (
                record["solver"]["source_commit"],
                record["solver"]["binary_sha256"],
            )
        )
    for prefix in prefix_partition(key)["prefix_columns"]:
        branch: dict[str, Any] = {}
        for run_id in RUN_IDS:
            path = run_path(key, "prefix", run_id, prefix)
            record = json.loads(path.read_text(encoding="utf-8"))
            validate_run_record(
                record, key, "prefix", run_id, prefix
            )
            branch[run_id] = {
                "path": str(path.relative_to(ROOT)),
                "sha256": file_sha256(path),
            }
            identities.add(
                (
                    record["solver"]["source_commit"],
                    record["solver"]["binary_sha256"],
                )
            )
        bindings["prefixes"][str(prefix)] = branch
    build = solver_build_binding()
    if identities != {(build["source_commit"], build["binary_sha256"])}:
        raise RuntimeError("EXP-070 runs do not match the prefix solver archive")
    return bindings


def certificate_protocol() -> dict[str, Any]:
    return {
        "algorithm": (
            "prefix-checkpointed connected-cluster exact enumeration plus "
            "exhaustive pure-second-block enumeration"
        ),
        "external_engine": (
            "QEC-pages/dist-m4ri legacy STANDALONE prefix-partitioned method 2"
        ),
        "side": "z",
        "rooted_presentations": list(PRESENTATIONS),
        "base_exhaustive_through_weight": BASE_CAP,
        "prefix_exact_weight": CAP,
        "replay_required": True,
        "resource_policy": {
            "processes": 1,
            "threads": 1,
            "minimum_nice_priority": MIN_NICE_PRIORITY,
        },
    }


def assemble(key: str) -> int:
    problem = problem_of(key)
    parity = E56.kernel_weight_parity_proof(problem, DISTANCE - 1)
    duality = E56.bb_duality_proof(problem)
    pure_second = pure_second_logical_enumeration(key)
    prefix = prefix_partition(key)
    payload = {
        "schema": SCHEMA,
        "utc": E56.utc_now(),
        "identity": certificate_identity(key),
        "root_reduction": root_reduction(key),
        "prefix_partition": prefix,
        "pure_second_enumeration": pure_second,
        "solver_build": solver_build_binding(),
        "protocol": certificate_protocol(),
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
        and prefix["complete"]
        and pure_second["complete"]
        and int(pure_second["minimum_weight"]) >= DISTANCE
    ):
        raise RuntimeError("EXP-070 lower-bound closure failed")
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
    pure_second = pure_second_logical_enumeration(key)
    prefix = prefix_partition(key)
    try:
        valid = bool(
            payload.get("schema") == SCHEMA
            and payload.get("identity") == certificate_identity(key)
            and payload.get("root_reduction") == root_reduction(key)
            and payload.get("prefix_partition") == prefix
            and payload.get("pure_second_enumeration") == pure_second
            and payload.get("solver_build") == solver_build_binding()
            and payload.get("protocol") == certificate_protocol()
            and payload.get("parity") == parity
            and payload.get("duality")
            == {
                field: value
                for field, value in duality.items()
                if field != "permutation"
            }
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
            and int(pure_second["minimum_weight"]) >= DISTANCE
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
    parser.add_argument("--mode", choices=("base", "prefix"))
    parser.add_argument("--prefix", type=int)
    parser.add_argument("--run-id", choices=RUN_IDS, default="initial")
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--solver-source-commit")
    parser.add_argument("--timeout", type=float, default=14_400.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.command == "prepare":
        print(json.dumps(prepare_inputs(), indent=1))
        return 0
    if args.target is None:
        parser.error(f"{args.command} requires --target")
    if args.command == "run":
        if (
            args.mode is None
            or args.solver is None
            or args.solver_source_commit is None
        ):
            parser.error(
                "run requires --mode, --solver, and --solver-source-commit"
            )
        if args.mode == "base" and args.prefix is not None:
            parser.error("base run cannot use --prefix")
        if args.mode == "prefix" and args.prefix is None:
            parser.error("prefix run requires --prefix")
        return run_search(
            args.target,
            args.mode,
            args.run_id,
            args.solver,
            args.solver_source_commit,
            args.timeout,
            prefix=args.prefix,
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
