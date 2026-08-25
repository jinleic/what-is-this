"""EXP-067: exact connected-cluster closure of the two n=234 BB bundles.

The minimum nontrivial word of an LDPC code may be chosen connected in the
variable-node Tanner graph: if its support disconnects, every component has
zero syndrome and at least one component is itself nontrivial.  Translations
are transitive inside each 117-qubit BB block.  Consequently two exhaustive
rooted searches suffice: root coordinate 0 in the original presentation, and
root coordinate 0 after swapping the two qubit blocks.

The external exact engine is QEC-pages/dist-m4ri's race-free legacy
STANDALONE method 2 (connected cluster),
pinned by source commit and binary SHA-256 in every run record.  A run is
accepted only after it exhausts every weight through 16 in both rooted
presentations.  Odd kernel weights are excluded independently, and a physical
weight-18 logical closes the distance exactly.

Commands:
  python experiments/exp067_n234_connected_cluster.py prepare
  python experiments/exp067_n234_connected_cluster.py run \
      --target bundle19 --presentation original --solver /path/to/dist_m4ri \
      --solver-source-commit COMMIT --timeout 3500 --force
  python experiments/exp067_n234_connected_cluster.py assemble --target bundle19
  python experiments/exp067_n234_connected_cluster.py validate --target bundle19
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.distance.sat_decide import css_logical_bases  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp067", "exp056_odd_distance.py")

SCHEMA = "exp067-n234-connected-cluster-v1"
RUN_SCHEMA = "exp067-connected-cluster-run-v1"
INPUT_SCHEMA = "exp067-connected-cluster-input-v1"
CAP = 16
DISTANCE = 18
DIST_M4RI_COMMIT = "538d119f6e98b3782415eb78f7ce322a076f8101"
DIST_M4RI_BINARY_SHA256 = (
    "04c051b27ae30dc353480015c01ce02b9925636d428a3a92c33a78ea5317ca66"
)
DIST_M4RI_SOURCE_SHA256 = (
    "e7fed74881ddd3a7bcd792e104dcddf37915927d68f53a333efbb3f7ab83f369"
)
DIST_M4RI_LICENSE_SHA256 = (
    "8177f97513213526df2cf6184d8ff986c675afb514d4e68a404010521b880643"
)
M4RI_SOURCE_SHA256 = (
    "986c9d4bf77b906ece535f0c405d3637dbdc8ec6bd5823c01c2f9823e3217d1c"
)
M4RI_STATIC_LIBRARY_SHA256 = (
    "b51aa6b7f720ab88d7259690545b919f34c7addb55accca35beb59d10fd267d4"
)
INPUT_DIR = ROOT / "results" / "partial_runs" / "exp067_n234_cluster" / "inputs"
RUN_DIR = ROOT / "results" / "partial_runs" / "exp067_n234_cluster" / "runs"
SOLVER_BUILD = (
    ROOT / "results" / "partial_runs" / "exp067_n234_cluster" / "solver_build.json"
)
PRESENTATIONS = ("original", "block_swapped")
RUN_IDS = ("initial", "replay")
LEGACY_RESULT_RE = re.compile(r"^-?\d+$")

TARGETS: dict[str, dict[str, Any]] = {
    "bundle19": {
        "bundle_index": 19,
        "name": "liang-twisted-torus-234-8-18-row-2",
        "source": "Liang et al. arXiv:2503.03827v3 Table III",
        "ell": 39,
        "m": 3,
        "expected_k": 8,
        "expected_d": 18,
        "A": [[0, 0], [1, 0], [5, 0]],
        "B": [[0, 0], [1, 1], [23, 2]],
        "witness_support": [
            11, 15, 17, 23, 26, 29, 32, 48, 101, 116, 117, 119, 133,
            137, 150, 166, 183, 203,
        ],
        "published": {
            "f_extra": [-1, 3],
            "g_extra": [3, -1],
            "lattice_basis": [[0, 39], [3, 6]],
            "expected_generator_images": {"x": [23, 2], "y": [34, 0]},
        },
    },
    "bundle22": {
        "bundle_index": 22,
        "name": "liang-twisted-torus-234-8-18-row-1",
        "source": "Liang et al. arXiv:2503.03827v3 Table III",
        "ell": 39,
        "m": 3,
        "expected_k": 8,
        "expected_d": 18,
        "A": [[0, 0], [1, 0], [5, 0]],
        "B": [[0, 0], [2, 2], [22, 1]],
        "witness_support": [
            14, 29, 31, 36, 51, 58, 66, 70, 82, 85, 101, 133, 146, 150,
            165, 180, 207, 215,
        ],
        "published": {
            "f_extra": [-1, -3],
            "g_extra": [3, -1],
            "lattice_basis": [[0, 39], [3, -9]],
            "expected_generator_images": {"x": [2, 2], "y": [5, 0]},
        },
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()
def solver_build_binding() -> dict[str, Any]:
    build = json.loads(SOLVER_BUILD.read_text(encoding="utf-8"))
    base = SOLVER_BUILD.parent
    artifacts = [
        (build["solver"]["binary"], DIST_M4RI_BINARY_SHA256),
        (build["solver"]["source_archive"], DIST_M4RI_SOURCE_SHA256),
        (build["solver"]["license_file"], DIST_M4RI_LICENSE_SHA256),
        (build["m4ri"]["source_archive"], M4RI_SOURCE_SHA256),
    ]
    valid = bool(
        build.get("schema") == "exp067-dist-m4ri-build-v1"
        and build.get("solver", {}).get("commit") == DIST_M4RI_COMMIT
        and build.get("solver", {}).get("binary_sha256")
        == DIST_M4RI_BINARY_SHA256
        and build.get("solver", {}).get("source_archive_sha256")
        == DIST_M4RI_SOURCE_SHA256
        and build.get("solver", {}).get("license_file_sha256")
        == DIST_M4RI_LICENSE_SHA256
        and build.get("solver", {}).get("entrypoint")
        == "legacy STANDALONE dist_rw.c -> single-thread do_CC_dist"
        and build.get("solver", {}).get("excluded_entrypoint")
        == "multithreaded dist_m4ri.c coordinator"
        and build.get("m4ri", {}).get("source_archive_sha256")
        == M4RI_SOURCE_SHA256
        and build.get("m4ri", {}).get("static_library_sha256")
        == M4RI_STATIC_LIBRARY_SHA256
        and all(file_sha256(base / path) == digest for path, digest in artifacts)
    )
    if not valid:
        raise RuntimeError("EXP-067 solver build provenance is stale")
    return {
        "path": str(SOLVER_BUILD.relative_to(ROOT)),
        "sha256": file_sha256(SOLVER_BUILD),
        "binary_sha256": DIST_M4RI_BINARY_SHA256,
        "source_commit": DIST_M4RI_COMMIT,
        "source_archive_sha256": DIST_M4RI_SOURCE_SHA256,
        "m4ri_source_archive_sha256": M4RI_SOURCE_SHA256,
    }




def target_payload(key: str) -> dict[str, Any]:
    target = TARGETS[key]
    return {
        field: value
        for field, value in target.items()
        if field not in {"bundle_index", "witness_support", "published"}
    }


def problem_of(key: str) -> dict[str, Any]:
    return E56.build_problem(target_payload(key))


def certificate_path(key: str) -> Path:
    return ROOT / "results" / "certificates" / f"exp067_234_8_18_{key}_distance.json"


def run_path(key: str, presentation: str, run_id: str = "initial") -> Path:
    if run_id not in RUN_IDS:
        raise ValueError(f"unknown run id {run_id!r}")
    suffix = "" if run_id == "initial" else f"_{run_id}"
    return RUN_DIR / f"{key}_{presentation}{suffix}.json"


def _add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (left[0] + right[0]) % 39, (left[1] + right[1]) % 3


def _scale(value: int, point: tuple[int, int]) -> tuple[int, int]:
    return (value * point[0]) % 39, (value * point[1]) % 3


def _image(
    exponent: tuple[int, int], x_image: tuple[int, int], y_image: tuple[int, int]
) -> tuple[int, int]:
    return _add(_scale(exponent[0], x_image), _scale(exponent[1], y_image))


def _canonical_support(support: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return min(
        tuple(
            sorted(
                ((x - ox) % 39, (y - oy) % 3)
                for x, y in support
            )
        )
        for ox, oy in support
    )


def _canonical_pair(
    left: list[tuple[int, int]], right: list[tuple[int, int]]
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    return tuple(sorted((_canonical_support(left), _canonical_support(right))))
def _translation_delta(
    source: list[tuple[int, int]], target: list[tuple[int, int]]
) -> tuple[int, int] | None:
    source_set, target_set = set(source), set(target)
    for delta in ((x, y) for x in range(39) for y in range(3)):
        if {_add(term, delta) for term in source_set} == target_set:
            return delta
    return None


def _rowspace_equal(left: np.ndarray, right: np.ndarray) -> bool:
    left_rank, right_rank = rank_np(left), rank_np(right)
    return bool(
        left_rank == right_rank
        and rank_np(np.vstack([left, right])) == left_rank
    )


def _published_matrix_transport(
    key: str, f_support: list[list[int]], g_support: list[list[int]]
) -> dict[str, Any]:
    source_terms = [
        [tuple(map(int, term)) for term in f_support],
        [tuple(map(int, term)) for term in g_support],
    ]
    target_terms = [
        [tuple(map(int, term)) for term in TARGETS[key]["A"]],
        [tuple(map(int, term)) for term in TARGETS[key]["B"]],
    ]
    source = E56.build_problem(
        {
            "ell": 39,
            "m": 3,
            "A": source_terms[0],
            "B": source_terms[1],
            "expected_k": 8,
        }
    )
    target = problem_of(key)
    matches = []
    for swap_blocks in (False, True):
        delta_a = _translation_delta(
            source_terms[1 if swap_blocks else 0], target_terms[0]
        )
        delta_b = _translation_delta(
            source_terms[0 if swap_blocks else 1], target_terms[1]
        )
        if delta_a is None or delta_b is None:
            continue
        mapping = np.empty(234, dtype=int)
        for coordinate in range(117):
            x, y = divmod(coordinate, 3)
            a_coordinate = 3 * ((x + delta_a[0]) % 39) + (y + delta_a[1]) % 3
            b_coordinate = 3 * ((x + delta_b[0]) % 39) + (y + delta_b[1]) % 3
            if swap_blocks:
                mapping[coordinate] = 117 + b_coordinate
                mapping[117 + coordinate] = a_coordinate
            else:
                mapping[coordinate] = a_coordinate
                mapping[117 + coordinate] = 117 + b_coordinate
        if sorted(mapping.tolist()) != list(range(234)):
            continue
        mapped_hx = np.zeros_like(source["HX"])
        mapped_hz = np.zeros_like(source["HZ"])
        mapped_hx[:, mapping] = source["HX"]
        mapped_hz[:, mapping] = source["HZ"]
        if _rowspace_equal(mapped_hx, target["HX"]) and _rowspace_equal(
            mapped_hz, target["HZ"]
        ):
            matches.append(
                {
                    "swap_blocks": swap_blocks,
                    "delta_A": list(delta_a),
                    "delta_B": list(delta_b),
                    "qubit_permutation_sha256": canonical_sha256(mapping.tolist()),
                    "HX_rowspace_equal": True,
                    "HZ_rowspace_equal": True,
                }
            )
    if len(matches) != 1:
        raise RuntimeError(f"published CSS matrix transport is not unique for {key}")
    return matches[0]




def published_transport(key: str) -> dict[str, Any]:
    """Find the unique quotient-lattice map from Liang et al. Table III."""
    target = TARGETS[key]
    row = target["published"]
    f_extra = tuple(map(int, row["f_extra"]))
    g_extra = tuple(map(int, row["g_extra"]))
    lattice = [tuple(map(int, vector)) for vector in row["lattice_basis"]]
    target_pair = _canonical_pair(
        [tuple(term) for term in target["A"]],
        [tuple(term) for term in target["B"]],
    )
    matches: list[dict[str, Any]] = []
    for x0 in range(39):
        for x1 in range(3):
            x_image = (x0, x1)
            for y0 in range(39):
                for y1 in range(3):
                    y_image = (y0, y1)
                    if any(
                        _image(relation, x_image, y_image) != (0, 0)
                        for relation in lattice
                    ):
                        continue
                    generated = {
                        _add(_scale(i, x_image), _scale(j, y_image))
                        for i in range(39)
                        for j in range(39)
                    }
                    if len(generated) != 117:
                        continue
                    f_support = [
                        (0, 0), x_image,
                        _image(f_extra, x_image, y_image),
                    ]
                    g_support = [
                        (0, 0), y_image,
                        _image(g_extra, x_image, y_image),
                    ]
                    if _canonical_pair(f_support, g_support) == target_pair:
                        matches.append(
                            {
                                "x_image": list(x_image),
                                "y_image": list(y_image),
                                "f_support": [list(term) for term in f_support],
                                "g_support": [list(term) for term in g_support],
                            }
                        )
    expected = row["expected_generator_images"]
    valid = bool(
        len(matches) == 1
        and matches[0]["x_image"] == expected["x"]
        and matches[0]["y_image"] == expected["y"]
    )
    if not valid:
        raise RuntimeError(f"published constructor transport failed for {key}")
    matrix_transport = _published_matrix_transport(
        key, matches[0]["f_support"], matches[0]["g_support"]
    )
    return {
        "schema": "exp067-liang-table-iii-transport-v1",
        "source": target["source"],
        "published_parameters": [234, 8, 18],
        "lattice_basis": [list(vector) for vector in lattice],
        "unique_maps": len(matches),
        "map": matches[0],
        "target_pair_sha256": canonical_sha256(target_pair),
        "matrix_transport": matrix_transport,
        "valid": True,
    }


def _presentation_matrices(
    problem: dict[str, Any], presentation: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if presentation not in PRESENTATIONS:
        raise ValueError(f"unknown presentation {presentation!r}")
    hx = np.asarray(problem["HX"], dtype=np.uint8)
    hz = np.asarray(problem["HZ"], dtype=np.uint8)
    lx, _lz = css_logical_bases(hx, hz)
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
def input_file_bindings(key: str, presentation: str) -> dict[str, Any]:
    problem = problem_of(key)
    expected = dict(
        zip(("HX", "HZ", "LX"), _presentation_matrices(problem, presentation))
    )
    bindings = {}
    for name, path in input_paths(key, presentation).items():
        if not path.is_file():
            raise RuntimeError(f"missing EXP-067 input file: {path}")
        loaded = np.asarray(mmread(path, spmatrix=True).toarray(), dtype=np.uint8) & 1
        if not np.array_equal(loaded, expected[name]):
            raise RuntimeError(f"EXP-067 input matrix changed: {path}")
        bindings[name] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": file_sha256(path),
        }
    return bindings
def _atomic_mmwrite(path: Path, matrix: np.ndarray) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as stream:
        mmwrite(
            stream,
            sparse.coo_matrix(matrix.astype(np.int8)),
            field="integer",
        )
    temporary.replace(path)






def prepare_inputs() -> dict[str, Any]:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    entries = []
    for key in TARGETS:
        problem = problem_of(key)
        for presentation in PRESENTATIONS:
            hx, hz, lx = _presentation_matrices(problem, presentation)
            paths = input_paths(key, presentation)
            for name, matrix in (("HX", hx), ("HZ", hz), ("LX", lx)):
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


def _parse_result(stdout: str) -> dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        token = line.strip()
        if LEGACY_RESULT_RE.fullmatch(token):
            value = int(token)
            if value < 0:
                return {
                    "status": "NO_LOGICAL_THROUGH_CAP",
                    "exhaustive_through": -value,
                    "found_distance": None,
                }
            return {
                "status": "LOGICAL_FOUND",
                "exhaustive_through": value - 1,
                "found_distance": value,
            }
    raise RuntimeError("legacy dist-m4ri output has no terminal integer result")


def run_search(
    key: str,
    presentation: str,
    solver: Path,
    solver_source_commit: str,
    timeout_s: float,
    run_id: str = "initial",
    *,
    force: bool,
) -> int:
    if not solver.is_file():
        raise RuntimeError(f"solver binary does not exist: {solver}")
    if solver_source_commit != DIST_M4RI_COMMIT:
        raise RuntimeError(
            f"solver source commit must be the pinned {DIST_M4RI_COMMIT}"
        )
    if file_sha256(solver) != DIST_M4RI_BINARY_SHA256:
        raise RuntimeError("solver binary does not match the audited legacy build")
    prepare_inputs()
    path = run_path(key, presentation, run_id)
    if path.exists() and not force:
        record = json.loads(path.read_text(encoding="utf-8"))
        validate_run_record(record, key, presentation, run_id)
        print(json.dumps(record["result"], indent=1))
        return 0
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    inputs = input_paths(key, presentation)
    args = [
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
    invocation = [str(solver.resolve()), *args]
    running = {
        "schema": RUN_SCHEMA,
        "utc": utc_now(),
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
        result = _parse_result(completed.stdout)
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
    record: dict[str, Any],
    key: str,
    presentation: str,
    run_id: str = "initial",
) -> None:
    try:
        result = record["result"]
        protocol = record["protocol"]
        valid = bool(
            record.get("schema") == RUN_SCHEMA
            and record.get("target") == key
            and record.get("presentation") == presentation
            and record.get("run_id", "initial") == run_id
            and record.get("status") == "COMPLETE"
            and int(record.get("returncode", -1)) == 0
            and record.get("input") == input_identity(key, presentation)
            and record.get("solver", {}).get("name")
            == "QEC-pages/dist-m4ri legacy STANDALONE method 2"
            and record.get("solver", {}).get("source")
            == "https://github.com/QEC-pages/dist-m4ri"
            and record.get("solver", {}).get("source_commit")
            == DIST_M4RI_COMMIT
            and record.get("solver", {}).get("binary_sha256")
            == DIST_M4RI_BINARY_SHA256
            and record.get("input_files") == input_file_bindings(key, presentation)
            and record.get("invocation")
            == [
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
            and protocol.get("method") == 2
            and protocol.get("engine_variant")
            == "legacy_single_thread_do_CC_dist"
            and protocol.get("root") == 0
            and protocol.get("dmin") == 1
            and protocol.get("wmax") == CAP
            and protocol.get("smax") == 0
            and result.get("status") == "NO_LOGICAL_THROUGH_CAP"
            and int(result.get("exhaustive_through", -1)) == CAP
            and result.get("found_distance") is None
            and _parse_result(record.get("stdout", "")) == result
            and hashlib.sha256(record.get("stdout", "").encode()).hexdigest()
            == record.get("stdout_sha256")
            and hashlib.sha256(record.get("stderr", "").encode()).hexdigest()
            == record.get("stderr_sha256")
        )
    except (KeyError, TypeError, ValueError, RuntimeError):
        valid = False
    if not valid:
        raise RuntimeError(
            f"invalid EXP-067 run record: {key}/{presentation}/{run_id}"
        )


def root_reduction(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    cover = E56.build_orbit_cover(problem)
    if not cover["complete"]:
        raise RuntimeError("translation/pole cover is incomplete")
    return {
        "schema": "exp067-two-root-connected-cluster-reduction-v1",
        "translation_group_size": int(problem["block"]),
        "translations_preserve_code": True,
        "presentations": list(PRESENTATIONS),
        "root_coordinate": 0,
        "connected_component_lemma": (
            "For a zero-syndrome support, Tanner-connected components have "
            "disjoint adjacent-check sets and therefore each has zero syndrome. "
            "Logical parities add, so a nontrivial word has a nontrivial "
            "component no heavier than itself."
        ),
        "irreducible_subset_lemma": (
            "A minimum-weight nontrivial kernel word has no nonempty proper "
            "zero-syndrome subset. If such a subset is logical it is a lighter "
            "logical; if it is a stabilizer, removing it leaves a lighter "
            "representative of the same logical class."
        ),
        "enumeration_argument": (
            "Every proper partial support of that irreducible word has nonzero "
            "syndrome. The first unsatisfied check used by dist-m4ri has another "
            "column in the remaining target support because the full syndrome "
            "is zero, so the recursive branch reaches the complete word."
        ),
        "coverage_argument": (
            "A connected logical touching the first BB block translates a "
            "supported qubit to global coordinate 0 and is enumerated in the "
            "original presentation. A logical supported only in the second "
            "block is covered after the exact block-column swap."
        ),
        "dist_m4ri_ordering_guard": (
            "start=0 is load-bearing: dist-m4ri enumerates only clusters whose "
            "root is the minimum column. start=117 is not used."
        ),
        "valid": True,
    }


def upper_certificate(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[np.asarray(TARGETS[key]["witness_support"], dtype=int)] = 1
    upper = E56.verify_upper_witness(problem, vector)
    if int(upper["weight"]) != DISTANCE:
        raise RuntimeError(f"weight-18 witness failed for {key}")
    return upper


def certificate_identity(key: str) -> dict[str, Any]:
    problem = problem_of(key)
    return {
        "target_key": key,
        "bundle_index": int(TARGETS[key]["bundle_index"]),
        "target": target_payload(key),
        "n": int(problem["n"]),
        "k": int(problem["k"]),
        "HX_sha256": E56.matrix_sha256(problem["HX"]),
        "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
    }


def _run_bindings(key: str) -> dict[str, Any]:
    bindings = {}
    solver_identities = set()
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
            solver_identities.add(
                (
                    record["solver"]["source_commit"],
                    record["solver"]["binary_sha256"],
                )
            )
    build = solver_build_binding()
    if solver_identities != {
        (build["source_commit"], build["binary_sha256"])
    }:
        raise RuntimeError("rooted runs do not match the archived solver build")
    return bindings


def assemble(key: str) -> int:
    problem = problem_of(key)
    parity = E56.kernel_weight_parity_proof(problem, DISTANCE - 1)
    duality = E56.bb_duality_proof(problem)
    runs = _run_bindings(key)
    exact = bool(
        parity["valid"]
        and int(parity["effective_even_cap"]) == CAP
        and duality["valid"]
        and upper_certificate(key)["weight"] == DISTANCE
    )
    payload = {
        "schema": SCHEMA,
        "utc": utc_now(),
        "identity": certificate_identity(key),
        "published_transport": published_transport(key),
        "root_reduction": root_reduction(key),
        "solver_build": solver_build_binding(),
        "protocol": {
            "algorithm": "connected-cluster exact enumeration",
            "external_engine": (
                "QEC-pages/dist-m4ri legacy STANDALONE method 2"
            ),
            "side": "z",
            "rooted_presentations": list(PRESENTATIONS),
            "exhaustive_through_weight": CAP,
            "replay_required": True,
        },
        "parity": parity,
        "duality": {
            field: value
            for field, value in duality.items()
            if field != "permutation"
        },
        "upper": upper_certificate(key),
        "lower_runs": runs,
        "verdict": (
            {
                "classification": "CERTIFIED_EXACT",
                "d": DISTANCE,
                "d_X": DISTANCE,
                "d_Z": DISTANCE,
                "exact": True,
            }
            if exact
            else {"classification": "INCOMPLETE", "exact": False}
        ),
    }
    path = certificate_path(key)
    E56.atomic_write_json(path, payload)
    if not exact:
        return 1
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def validate_exact_certificate_payload(payload: dict[str, Any]) -> None:
    key = payload.get("identity", {}).get("target_key")
    if key not in TARGETS:
        raise RuntimeError("unknown EXP-067 target")
    problem = problem_of(key)
    parity = E56.kernel_weight_parity_proof(problem, DISTANCE - 1)
    duality = E56.bb_duality_proof(problem)
    expected_runs = _run_bindings(key)
    try:
        valid = bool(
            payload.get("schema") == SCHEMA
            and payload.get("identity") == certificate_identity(key)
            and payload.get("published_transport") == published_transport(key)
            and payload.get("root_reduction") == root_reduction(key)
            and payload.get("solver_build") == solver_build_binding()
            and payload.get("protocol")
            == {
                "algorithm": "connected-cluster exact enumeration",
                "external_engine": (
                    "QEC-pages/dist-m4ri legacy STANDALONE method 2"
                ),
                "side": "z",
                "rooted_presentations": list(PRESENTATIONS),
                "exhaustive_through_weight": CAP,
                "replay_required": True,
            }
            and payload.get("parity") == parity
            and payload.get("duality")
            == {field: value for field, value in duality.items() if field != "permutation"}
            and payload.get("upper") == upper_certificate(key)
            and payload.get("lower_runs") == expected_runs
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
        raise RuntimeError(f"EXP-067 exact certificate failed validation: {key}")
def _refresh_exp066_summary(e55: Any, e66: Any) -> dict[str, Any]:
    shards = {}
    records = []
    initial_residual = []
    for ell, m in e66.LATTICES:
        key = f"{ell}x{m}"
        path = e55.SCREEN_DIR / f"{key}.json"
        shard = json.loads(path.read_text(encoding="utf-8"))
        shards[key] = shard
        records.extend(shard["records"])
        if (ell, m) == (39, 3):
            initial_residual = [
                record
                for record in shard["records"]
                if record.get("initial_verdict") == "solver_required"
            ]
    verdicts = Counter(record["verdict"] for record in records)
    initial_counts = Counter(
        record.get("initial_verdict", record["verdict"]) for record in records
    )
    dominated = sum(
        value for name, value in verdicts.items() if name.startswith("dominated")
    )
    survivors = verdicts.get("survivor", 0)
    undecided = verdicts.get("undecided", 0)
    with_reference = len(records) - verdicts.get("no_reference", 0)
    global_screen = json.loads(e55.SCREEN_OUT.read_text(encoding="utf-8"))
    summary = {
        "schema": e66.SCHEMA,
        "utc": utc_now(),
        "protocol": e66._screen_protocol(),
        "scope": {
            "lattices": [list(pair) for pair in e66.LATTICES],
            "n": 234,
            "k_range": [e66.K_MIN, e66.K_MAX],
            "classes": len(records),
            "represented_pairs": sum(int(record["orbit"]) for record in records),
        },
        "initial": {
            "dominated": sum(
                value
                for name, value in initial_counts.items()
                if name.startswith("dominated")
            ),
            "solver_required": initial_counts.get("solver_required", 0),
        },
        "automorphism_bundles": e66.automorphism_bundles(initial_residual),
        "routes": dict(sorted(verdicts.items())),
        "shards": {
            key: {
                "path": str(
                    (e55.SCREEN_DIR / f"{key}.json").relative_to(ROOT)
                ),
                "sha256": e55._file_sha256(e55.SCREEN_DIR / f"{key}.json"),
                "records_sha256": e66._records_identity(shard["records"]),
            }
            for key, shard in shards.items()
        },
        "verdict": {
            "complete": dominated + survivors == with_reference,
            "dominated": dominated,
            "survivors": survivors,
            "undecided": undecided,
            "all_referenced_dominated": dominated == with_reference,
            "exact_promotions": [],
            "exact_candidates_requiring_promotion": [],
        },
        "global_screen_sha256": e55._file_sha256(e55.SCREEN_OUT),
        "global_verdict": global_screen["verdict"],
    }
    E56.atomic_write_json(e66.SUMMARY, summary)
    e66.validate_summary(summary)
    return summary


def promote_reference() -> int:
    for key in TARGETS:
        certificate = json.loads(
            certificate_path(key).read_text(encoding="utf-8")
        )
        validate_exact_certificate_payload(certificate)
    e55 = _load("exp055_for_exp067_promotion", "exp055_odd_lattice_sweep.py")
    e63 = _load("exp063_for_exp067_promotion", "exp063_reference_rebind.py")
    old_screen = json.loads(e55.SCREEN_OUT.read_text(encoding="utf-8"))
    lattices = [
        f"{int(shard['ell'])}x{int(shard['m'])}"
        for shard in old_screen["lattices"]
    ]
    e63.run(lattices)
    args = argparse.Namespace(
        k_min=8,
        k_max=24,
        time_limit=120.0,
        lattices=",".join(lattices),
        min_n=18,
    )
    if e55.assemble_screen(args) != 0:
        raise RuntimeError("EXP-067 reference reassembly is incomplete")
    e55.certify_screen_survivors(args)
    e66 = _load("exp066_for_exp067_promotion", "exp066_n234_frontier.py")
    summary = _refresh_exp066_summary(e55, e66)
    print(
        json.dumps(
            {
                "exp067_exact_targets": list(TARGETS),
                "local_verdict": summary["verdict"],
                "global_verdict": summary["global_verdict"],
            },
            indent=1,
        )
    )
    return 0




def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("prepare", "run", "assemble", "validate", "promote")
    )
    parser.add_argument("--target", choices=tuple(TARGETS))
    parser.add_argument("--presentation", choices=PRESENTATIONS)
    parser.add_argument("--run-id", choices=RUN_IDS, default="initial")
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--solver-source-commit")
    parser.add_argument("--timeout", type=float, default=3500.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.command == "prepare":
        print(json.dumps(prepare_inputs(), indent=1))
        return 0
    if args.command == "promote":
        return promote_reference()
    if args.target is None:
        parser.error("--target is required")
    if args.command == "run":
        if args.presentation is None or args.solver is None or args.solver_source_commit is None:
            parser.error("run requires --presentation, --solver, and --solver-source-commit")
        return run_search(
            args.target,
            args.presentation,
            args.solver,
            args.solver_source_commit,
            args.timeout,
            args.run_id,
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
