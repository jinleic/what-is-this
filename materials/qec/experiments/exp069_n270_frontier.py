"""EXP-069: canonical n=270 odd-lattice frontier transport and screen.

The noncyclic presentations are not distinct groups.  The explicit isomorphism

    Z_15 x Z_9 -> Z_45 x Z_3,
    (x, y) |-> (CRT(y mod 9, x mod 5), x mod 3)

is audited exhaustively and used to select (15,9) as the sole screened
presentation.  Liang et al.'s two Table III [[270,8,20]] twisted-torus rows are
also transported through their quotient lattices and then through this exact
rectangular-presentation isomorphism.

Run:
  uv run python experiments/exp069_n270_frontier.py transport
  uv run python experiments/exp069_n270_frontier.py initial --lattice 15x9
  uv run python experiments/exp069_n270_frontier.py initial --lattice 27x5
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

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


E55 = _load("exp055_for_exp069", "exp055_odd_lattice_sweep.py")
E56 = _load("exp056_for_exp069", "exp056_odd_distance.py")

from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp069-n270-frontier-v1"
CANONICAL_NONCYCLIC = (15, 9)
DUPLICATE_NONCYCLIC = (45, 3)
CYCLIC = (27, 5)
K_MIN = 8
K_MAX = 24
TIME_LIMIT_S = 120.0
DISTANCE_REFERENCE_THRESHOLD = 20
PARTIAL_DIR = ROOT / "results" / "partial_runs" / "exp069_n270"

LIANG_TARGETS: dict[str, dict[str, Any]] = {
    "liang270a": {
        "name": "liang-twisted-torus-270-8-20-green",
        "source": "Liang et al. arXiv:2503.03827v3 Table III",
        "source_url": "https://arxiv.org/abs/2503.03827",
        "expected_k": 8,
        "expected_d": 20,
        "f_extra": [-1, -3],
        "g_extra": [3, -1],
        "lattice_basis": [[0, 15], [9, 6]],
        "quotient_first_coordinate": [1, 6],
        "target_45x3": {
            "A": [[0, 0], [1, 0], [26, 0]],
            "B": [[0, 0], [6, 1], [42, 2]],
        },
        "local_15x9": {
            "A": [[0, 0], [6, 1], [6, 8]],
            "B": [[0, 0], [1, 6], [2, 6]],
        },
        "witness_support": [
            0, 29, 34, 77, 82, 89, 109, 111, 116, 132,
            138, 144, 159, 172, 181, 190, 196, 199, 214, 264,
        ],
    },
    "liang270b": {
        "name": "liang-twisted-torus-270-8-20-cyan",
        "source": "Liang et al. arXiv:2503.03827v3 Table III",
        "source_url": "https://arxiv.org/abs/2503.03827",
        "expected_k": 8,
        "expected_d": 20,
        "f_extra": [-1, 3],
        "g_extra": [3, -1],
        "lattice_basis": [[0, 45], [3, -12]],
        "quotient_first_coordinate": [1, 4],
        "target_45x3": {
            "A": [[0, 0], [1, 0], [11, 0]],
            "B": [[0, 0], [4, 1], [44, 2]],
        },
        "local_15x9": {
            "A": [[0, 0], [6, 1], [6, 2]],
            "B": [[0, 0], [4, 4], [14, 8]],
        },
        "witness_support": [
            14, 23, 37, 38, 45, 70, 71, 78, 94, 101,
            102, 118, 125, 134, 141, 164, 172, 197, 204, 228,
        ],
    },
}


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def z15_z9_to_z45_z3(point: tuple[int, int]) -> tuple[int, int]:
    """The CRT isomorphism from the canonical to duplicate presentation."""
    x, y = point
    x %= 15
    y %= 9
    return (y + 9 * ((y - x) % 5)) % 45, x % 3


def z45_z3_to_z15_z9(point: tuple[int, int]) -> tuple[int, int]:
    """Inverse of :func:`z15_z9_to_z45_z3`."""
    u, v = point
    u %= 45
    v %= 3
    return (v + 3 * ((2 * (u - v)) % 5)) % 15, u % 9


def _add(
    left: tuple[int, int], right: tuple[int, int], moduli: tuple[int, int]
) -> tuple[int, int]:
    return (
        (left[0] + right[0]) % moduli[0],
        (left[1] + right[1]) % moduli[1],
    )


def coordinate_permutation() -> np.ndarray:
    """Source flat index -> target flat index for all 135 group elements."""
    return np.asarray(
        [
            3 * u + v
            for x in range(15)
            for y in range(9)
            for u, v in [z15_z9_to_z45_z3((x, y))]
        ],
        dtype=np.int16,
    )


def _transport_terms_to_45x3(
    terms: list[list[int]] | list[tuple[int, int]],
) -> list[list[int]]:
    return [list(z15_z9_to_z45_z3(tuple(term))) for term in terms]


def _transport_terms_to_15x9(
    terms: list[list[int]] | list[tuple[int, int]],
) -> list[list[int]]:
    return [list(z45_z3_to_z15_z9(tuple(term))) for term in terms]


def _rectangular_matrix_transport(
    local_a: list[list[int]], local_b: list[list[int]]
) -> dict[str, Any]:
    target_a = _transport_terms_to_45x3(local_a)
    target_b = _transport_terms_to_45x3(local_b)
    local_hx, local_hz = E55.E53.bb_from_terms(
        15, 9, local_a, local_b
    )
    target_hx, target_hz = E55.E53.bb_from_terms(
        45, 3, target_a, target_b
    )
    coordinate_map = coordinate_permutation().astype(int)
    qubit_map = np.r_[coordinate_map, 135 + coordinate_map]
    mapped_hx = np.zeros_like(local_hx)
    mapped_hz = np.zeros_like(local_hz)
    mapped_hx[np.ix_(coordinate_map, qubit_map)] = local_hx
    mapped_hz[np.ix_(coordinate_map, qubit_map)] = local_hz
    return {
        "coordinate_permutation_sha256": hashlib.sha256(
            coordinate_map.astype(np.int16).tobytes()
        ).hexdigest(),
        "qubit_permutation_bijective": bool(
            sorted(qubit_map.tolist()) == list(range(270))
        ),
        "HX_exact": bool(np.array_equal(mapped_hx, target_hx)),
        "HZ_exact": bool(np.array_equal(mapped_hz, target_hz)),
    }


def _published_quotient_image(
    key: str, point: tuple[int, int]
) -> tuple[int, int]:
    target = LIANG_TARGETS[key]
    alpha, beta = target["quotient_first_coordinate"]
    x, y = point
    return (alpha * x + beta * y) % 45, y % 3


def published_transport(key: str) -> dict[str, Any]:
    """Bind one Table III twisted torus to both rectangular presentations."""
    if key not in LIANG_TARGETS:
        raise KeyError(key)
    target = LIANG_TARGETS[key]
    f_support = [(0, 0), (1, 0), tuple(target["f_extra"])]
    g_support = [(0, 0), (0, 1), tuple(target["g_extra"])]
    target_a = [list(_published_quotient_image(key, term)) for term in f_support]
    target_b = [list(_published_quotient_image(key, term)) for term in g_support]
    if target_a != target["target_45x3"]["A"] or target_b != target["target_45x3"]["B"]:
        raise RuntimeError("published quotient transport changed")
    local_a = _transport_terms_to_15x9(target_a)
    local_b = _transport_terms_to_15x9(target_b)
    if local_a != target["local_15x9"]["A"] or local_b != target["local_15x9"]["B"]:
        raise RuntimeError("rectangular presentation transport changed")

    basis = [tuple(vector) for vector in target["lattice_basis"]]
    determinant = abs(
        basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0]
    )
    relations_vanish = all(
        _published_quotient_image(key, vector) == (0, 0) for vector in basis
    )
    quotient_images = {
        _published_quotient_image(key, (x, y))
        for x in range(45)
        for y in range(45)
    }
    target_hx, target_hz = E55.E53.bb_from_terms(
        45, 3, target_a, target_b
    )
    local_hx, local_hz = E55.E53.bb_from_terms(
        15, 9, local_a, local_b
    )
    target_k = 270 - rank_np(target_hx) - rank_np(target_hz)
    local_k = 270 - rank_np(local_hx) - rank_np(local_hz)
    if determinant != 135 or len(quotient_images) != 135:
        raise RuntimeError("published lattice quotient is not order 135")
    return {
        "source": target["source"],
        "source_url": target["source_url"],
        "published": {
            "f_support": [list(term) for term in f_support],
            "g_support": [list(term) for term in g_support],
            "lattice_basis": target["lattice_basis"],
            "reported_parameters": [270, 8, 20],
        },
        "quotient_order": determinant,
        "quotient_map": {
            "formula": (
                f"(x,y) -> (x+{target['quotient_first_coordinate'][1]}y mod 45, "
                "y mod 3)"
            ),
            "generator_images": {
                "x": list(_published_quotient_image(key, (1, 0))),
                "y": list(_published_quotient_image(key, (0, 1))),
            },
        },
        "lattice_relations_vanish": relations_vanish,
        "quotient_map_surjective": len(quotient_images) == 135,
        "target_45x3": {
            "A": target_a,
            "B": target_b,
            "k": int(target_k),
        },
        "local_15x9": {
            "A": local_a,
            "B": local_b,
            "k": int(local_k),
        },
        "matrix_transport": _rectangular_matrix_transport(local_a, local_b),
    }


@lru_cache(maxsize=1)
def presentation_transport_audit() -> dict[str, Any]:
    source = [(x, y) for x in range(15) for y in range(9)]
    coordinate_map = coordinate_permutation().astype(int)
    inverse_exact = all(
        z45_z3_to_z15_z9(z15_z9_to_z45_z3(point)) == point
        for point in source
    )
    homomorphism = all(
        z15_z9_to_z45_z3(_add(left, right, (15, 9)))
        == _add(
            z15_z9_to_z45_z3(left),
            z15_z9_to_z45_z3(right),
            (45, 3),
        )
        for left in source
        for right in source
    )
    candidates_15x9 = E55.enumerate_candidates(15, 9, K_MIN, K_MAX)
    candidates_45x3 = E55.enumerate_candidates(45, 3, K_MIN, K_MAX)
    return {
        "schema": "exp069-presentation-transport-audit-v1",
        "source_group": "Z_15 x Z_9",
        "target_group": "Z_45 x Z_3",
        "group_order": 135,
        "formula": (
            "(x,y) -> (y + 9*((y-x) mod 5) mod 45, x mod 3); "
            "inverse (u,v) -> (v + 3*(2*(u-v) mod 5) mod 15, u mod 9)"
        ),
        "all_coordinates_bijective": bool(
            sorted(coordinate_map.tolist()) == list(range(135))
        ),
        "inverse_exact": inverse_exact,
        "homomorphism_exhaustive": homomorphism,
        "coordinate_permutation_sha256": hashlib.sha256(
            coordinate_map.astype(np.int16).tobytes()
        ).hexdigest(),
        "candidate_classes": {
            "15x9": len(candidates_15x9),
            "45x3": len(candidates_45x3),
        },
        "represented_pairs": {
            "15x9": sum(int(record["orbit"]) for record in candidates_15x9),
            "45x3": sum(int(record["orbit"]) for record in candidates_45x3),
        },
        "coverage_argument": (
            "The coordinate map is a group isomorphism, hence its linear extension "
            "is a group-algebra isomorphism preserving convolution, inversion, "
            "Hamming weight, and both BB CSS rowspaces under the displayed row and "
            "qubit permutations. The (15,9) exhaustive raw weight-3 pair census "
            "therefore covers every (45,3) code; screening the latter would duplicate "
            "the same physical code family."
        ),
        "screened_presentation": list(CANONICAL_NONCYCLIC),
        "excluded_duplicate_presentation": list(DUPLICATE_NONCYCLIC),
    }


def transport_payload() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "presentation_transport": presentation_transport_audit(),
        "published_transports": {
            key: published_transport(key) for key in LIANG_TARGETS
        },
    }


def screen_protocol() -> dict[str, Any]:
    return E55._screen_protocol(
        E55._file_sha256(E55.OUT),
        E55._reference_fingerprint(),
        K_MIN,
        K_MAX,
        TIME_LIMIT_S,
    )


def verify_witness(
    record: dict[str, Any], support: list[int] | None
) -> dict[str, Any]:
    if support is None:
        raise ValueError("missing witness support")
    problem = E56.build_problem(record)
    witness = np.zeros(int(problem["n"]), dtype=np.uint8)
    witness[np.asarray(support, dtype=int)] = 1
    return E56.verify_upper_witness(problem, witness)


def _deep_seed(record: dict[str, Any]) -> int:
    return int(
        canonical_sha256(
            {
                "schema": "exp069-deep-information-set-reduction-v1",
                "ell": record["ell"],
                "m": record["m"],
                "A": record["A"],
                "B": record["B"],
                "threshold": record["threshold"],
            }
        )[:16],
        16,
    )


def deepen_record(
    record: dict[str, Any], *, tries: int = 2_000
) -> dict[str, Any]:
    """Deterministically reduce every pole class; no external solver or thread."""
    ell, m = int(record["ell"]), int(record["m"])
    problem = E56.build_problem(record)
    pole = E55.intersect(
        E55.ann_basis(record["A"], ell, m),
        E55.ann_basis(record["B"], ell, m),
    )
    physical_pole = pole[:, E55.bar_permutation(ell, m)]
    seed = _deep_seed(record)
    started = time.perf_counter()
    result = E55.reduced_witness_bound(
        physical_pole,
        problem["HZ"],
        ell,
        m,
        tries=int(tries),
        keep=E55.WITNESS_KEEP,
        seed=seed,
        target=int(record["threshold"]),
    )
    support = result.get("witness_support")
    verification = verify_witness(record, support) if support is not None else None
    if (
        verification is not None
        and int(verification["weight"]) != int(result["bound"])
    ):
        raise RuntimeError("deep-reduction weight disagrees with its support")
    return {
        "schema": "exp069-deep-information-set-reduction-v1",
        "tries": int(tries),
        "seed": seed,
        "classes_tested": int(result.get("classes_tested", 0)),
        "bound": result.get("bound"),
        "witness_support": support,
        "verification": verification,
        "wall_time_s": time.perf_counter() - started,
    }


def solve_residual_record(
    record: dict[str, Any], *, time_limit_s: float = TIME_LIMIT_S
) -> dict[str, Any]:
    """Decide one residual serially and retain any physical domination witness."""
    problem = E56.build_problem(record)
    started = time.perf_counter()
    cdcl = E55._cdcl_witness_bound(
        problem["HX"],
        problem["HZ"],
        int(problem["block"]),
        int(record["threshold"]),
    )
    output = dict(record)
    output["initial_verdict"] = record["verdict"]
    output["cdcl"] = cdcl
    if cdcl["status"] == "SAT":
        support = cdcl["witness_support"]
        verification = verify_witness(record, support)
        output.update(
            {
                "verdict": "dominated_by_cdcl_witness",
                "solver_calls": 1,
                "witness_bound": int(verification["weight"]),
                "witness_support": support,
                "wall_s": round(time.perf_counter() - started, 3),
            }
        )
        return output

    exact = E55.exact_distance_css(
        problem["HX"],
        problem["HZ"],
        time_limit_s=float(time_limit_s),
        workers=1,
        upper_bound=int(record["threshold"]),
    )
    output["exact_decision"] = exact
    decided = bool(
        exact["d_X_all_sectors_decided"]
        and exact["d_Z_all_sectors_decided"]
    )
    output.update(
        {
            "solver_calls": 2,
            "screen_decided": decided,
            "d_found": exact["d"],
            "d_exact": bool(exact["d_exact"]),
            "wall_s": round(time.perf_counter() - started, 3),
        }
    )
    if exact["d"] is None:
        output["verdict"] = "survivor" if decided else "undecided"
        return output

    threshold = int(record["threshold"])
    support: list[int] | None = None
    if exact["d_Z"] is not None and int(exact["d_Z"]) <= threshold:
        support = exact["d_Z_witness"]
    elif exact["d_X"] is not None and int(exact["d_X"]) <= threshold:
        x_vector = np.zeros(int(problem["n"]), dtype=np.uint8)
        x_vector[np.asarray(exact["d_X_witness"], dtype=int)] = 1
        permutation = np.asarray(
            E56.bb_duality_proof(problem)["permutation"], dtype=int
        )
        z_vector = (
            x_vector
            @ np.eye(int(problem["n"]), dtype=np.uint8)[permutation]
        ) % 2
        support = [
            int(index) for index in np.flatnonzero(z_vector)
        ]
    if support is None:
        output["verdict"] = "undecided"
        return output
    verification = verify_witness(record, support)
    if int(verification["weight"]) > threshold:
        raise RuntimeError("exact fallback witness exceeds the threshold")
    output.update(
        {
            "verdict": "dominated_by_exact_witness",
            "witness_bound": int(verification["weight"]),
            "witness_support": support,
            "witness_verification": verification,
        }
    )
    return output


def finalize_deep_record(
    record: dict[str, Any], deep: dict[str, Any]
) -> dict[str, Any]:
    output = dict(record)
    output["initial_verdict"] = record["verdict"]
    output["deep_reduction"] = deep
    bound = deep.get("bound")
    if bound is not None:
        output["witness_bound"] = int(bound)
        output["witness_support"] = deep["witness_support"]
    if bound is not None and int(bound) <= int(record["threshold"]):
        output["verdict"] = "dominated_by_deep_reduction"
    else:
        output["verdict"] = "undecided"
    output["solver_calls"] = 0
    output["wall_s"] = round(float(deep["wall_time_s"]), 3)
    return output


def resolve_candidate_record(
    record: dict[str, Any], *, tries: int
) -> dict[str, Any]:
    if int(record["k_parent"]) == 8:
        if int(record["threshold"]) < DISTANCE_REFERENCE_THRESHOLD:
            output = dict(record)
            output.update(
                {
                    "initial_verdict": record["verdict"],
                    "verdict": "undecided",
                    "solver_calls": 0,
                    "wall_s": 0.0,
                    "reference_pending": {
                        "required_threshold": DISTANCE_REFERENCE_THRESHOLD,
                        "reason": (
                            "The two published n=270 distance-20 rows are not "
                            "admissible until their local lower certificates pass."
                        ),
                    },
                }
            )
            return output
        return finalize_deep_record(
            record, deepen_record(record, tries=tries)
        )
    if int(record["k_parent"]) < 20:
        output = dict(record)
        output.update(
            {
                "initial_verdict": record["verdict"],
                "verdict": "undecided",
                "solver_calls": 0,
                "wall_s": 0.0,
                "exact_pending": {
                    "resource_reason": (
                        "serial per-sector fallback exceeded the bounded "
                        "screen budget"
                    ),
                    "observed_route": (
                        "CP-SAT did not finish the first k=12 comparison "
                        "class within five minutes"
                    ),
                },
            }
        )
        return output
    return solve_residual_record(record)


def _initial_path(
    lattice: tuple[int, int], protocol: dict[str, Any]
) -> Path:
    ell, m = lattice
    protocol_id = canonical_sha256(protocol)[:12]
    return PARTIAL_DIR / f"initial_{ell}x{m}_{protocol_id}.json"


def _candidate_identity(candidates: list[dict[str, Any]]) -> str:
    return canonical_sha256(
        [
            {
                "A": record["A"],
                "B": record["B"],
                "k_parent": record["k_parent"],
                "orbit": record["orbit"],
            }
            for record in candidates
        ]
    )


def screen_initial(
    lattice: tuple[int, int], *, checkpoint_every: int = 25
) -> dict[str, Any]:
    """Run the solver-free screen serially, with resumable atomic checkpoints."""
    if lattice not in {CANONICAL_NONCYCLIC, CYCLIC}:
        raise ValueError(
            f"{lattice!r} is not a canonical n=270 screen presentation"
        )
    if checkpoint_every <= 0:
        raise ValueError("checkpoint_every must be positive")
    ell, m = lattice
    protocol = screen_protocol()
    candidates = E55.enumerate_candidates(ell, m, K_MIN, K_MAX)
    candidate_sha256 = _candidate_identity(candidates)
    path = _initial_path(lattice, protocol)
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (
            payload.get("schema") != "exp069-initial-screen-v1"
            or payload.get("protocol") != protocol
            or payload.get("lattice") != [ell, m]
            or payload.get("candidate_identity_sha256") != candidate_sha256
            or int(payload.get("candidates", -1)) != len(candidates)
            or int(payload.get("next_index", -1))
            != len(payload.get("records", []))
        ):
            raise RuntimeError("stale EXP-069 initial-screen checkpoint")
        if payload.get("status") == "COMPLETE":
            if int(payload["next_index"]) != len(candidates):
                raise RuntimeError("incomplete EXP-069 checkpoint marked complete")
            return payload
        records = list(payload["records"])
        started_wall_s = float(payload.get("wall_s", 0.0))
    else:
        records = []
        started_wall_s = 0.0

    started = time.perf_counter()
    resource_policy = {
        "processes": 1,
        "threads": 1,
        "external_solvers": False,
    }

    def checkpoint(status: str) -> None:
        payload = {
            "schema": "exp069-initial-screen-v1",
            "utc": E56.utc_now(),
            "status": status,
            "protocol": protocol,
            "lattice": [ell, m],
            "n": 2 * ell * m,
            "k_range": [K_MIN, K_MAX],
            "candidate_identity_sha256": candidate_sha256,
            "candidates": len(candidates),
            "represented_pairs": sum(
                int(record["orbit"]) for record in candidates
            ),
            "next_index": len(records),
            "records": records,
            "resource_policy": resource_policy,
            "wall_s": started_wall_s + time.perf_counter() - started,
        }
        E56.atomic_write_json(path, payload)

    for index in range(len(records), len(candidates)):
        records.append(
            E55._screen_one(
                (ell, m, candidates[index], TIME_LIMIT_S, False)
            )
        )
        if len(records) % checkpoint_every == 0:
            checkpoint("RUNNING")
    checkpoint("COMPLETE")
    return json.loads(path.read_text(encoding="utf-8"))


def _resolved_path(
    lattice: tuple[int, int], protocol: dict[str, Any], tries: int
) -> Path:
    ell, m = lattice
    protocol_id = canonical_sha256(protocol)[:12]
    return PARTIAL_DIR / f"resolved_v2_{ell}x{m}_{protocol_id}_t{tries}.json"


def _screen_lattice_payload(
    lattice: tuple[int, int],
    protocol: dict[str, Any],
    records: list[dict[str, Any]],
    wall_s: float,
    deep_tries: int,
) -> dict[str, Any]:
    ell, m = lattice
    verdicts = Counter(record["verdict"] for record in records)
    return {
        "schema": "exp055-screen-v3",
        "utc": E56.utc_now(),
        "protocol": protocol,
        "exp069": {
            "schema": SCHEMA,
            "deep_tries": int(deep_tries),
            "processes": 1,
            "max_solver_workers": 1,
            "external_solvers": True,
        },
        "ell": ell,
        "m": m,
        "n": 2 * ell * m,
        "k_range": [K_MIN, K_MAX],
        "time_limit_s": TIME_LIMIT_S,
        "candidates_after_symmetry": len(records),
        "orbit_total": sum(int(record["orbit"]) for record in records),
        "verdicts": dict(sorted(verdicts.items())),
        "survivors": [
            record for record in records if record["verdict"] == "survivor"
        ],
        "no_reference": [
            record for record in records if record["verdict"] == "no_reference"
        ],
        "undecided": [
            record for record in records if record["verdict"] == "undecided"
        ],
        "solver_calls": sum(
            int(record.get("solver_calls", 0)) for record in records
        ),
        "records": records,
        "wall_s": round(float(wall_s), 3),
    }


def resolve_screen(
    lattice: tuple[int, int],
    *,
    tries: int = 2_000,
    checkpoint_every: int = 1,
) -> dict[str, Any]:
    """Resolve solver-free residuals serially and publish a validated shard."""
    if tries <= 0 or checkpoint_every <= 0:
        raise ValueError("tries and checkpoint_every must be positive")
    initial = screen_initial(lattice)
    if initial.get("status") != "COMPLETE":
        raise RuntimeError("EXP-069 initial screen is incomplete")
    protocol = initial["protocol"]
    initial_path = _initial_path(lattice, protocol)
    initial_binding = {
        "path": str(initial_path.relative_to(ROOT)),
        "sha256": E55._file_sha256(initial_path),
    }
    residual_indexes = [
        index
        for index, record in enumerate(initial["records"])
        if record["verdict"] == "solver_required"
    ]
    path = _resolved_path(lattice, protocol, tries)
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (
            payload.get("schema") != "exp069-resolved-screen-v2"
            or payload.get("protocol") != protocol
            or payload.get("lattice") != list(lattice)
            or payload.get("initial") != initial_binding
            or payload.get("residual_indexes") != residual_indexes
            or int(payload.get("tries", -1)) != tries
            or int(payload.get("next_residual", -1))
            != len(payload.get("resolved_records", []))
        ):
            raise RuntimeError("stale EXP-069 resolved-screen checkpoint")
        resolved_records = list(payload["resolved_records"])
        started_wall_s = float(payload.get("wall_s", 0.0))
    else:
        resolved_records = []
        started_wall_s = 0.0

    started = time.perf_counter()

    def checkpoint(status: str) -> None:
        E56.atomic_write_json(
            path,
            {
                "schema": "exp069-resolved-screen-v2",
                "utc": E56.utc_now(),
                "status": status,
                "protocol": protocol,
                "lattice": list(lattice),
                "initial": initial_binding,
                "residual_indexes": residual_indexes,
                "tries": int(tries),
                "next_residual": len(resolved_records),
                "resolved_records": resolved_records,
                "resource_policy": {
                    "processes": 1,
                    "max_solver_workers": 1,
                    "external_solvers": True,
                },
                "wall_s": (
                    started_wall_s + time.perf_counter() - started
                ),
            },
        )

    if (
        path.exists()
        and len(resolved_records) == len(residual_indexes)
        and json.loads(path.read_text(encoding="utf-8")).get("status")
        == "COMPLETE"
    ):
        pass
    else:
        for residual_position in range(
            len(resolved_records), len(residual_indexes)
        ):
            record = initial["records"][
                residual_indexes[residual_position]
            ]
            resolved_records.append(
                resolve_candidate_record(record, tries=tries)
            )
            if len(resolved_records) % checkpoint_every == 0:
                checkpoint("RUNNING")
        checkpoint("COMPLETE")

    records = list(initial["records"])
    for index, record in zip(residual_indexes, resolved_records):
        records[index] = record
    final_payload = _screen_lattice_payload(
        lattice,
        protocol,
        records,
        float(initial["wall_s"])
        + float(
            json.loads(path.read_text(encoding="utf-8"))["wall_s"]
        ),
        tries,
    )
    final_path = E55.SCREEN_DIR / f"{lattice[0]}x{lattice[1]}.json"
    E56.atomic_write_json(final_path, final_payload)
    E55._validate_screen_shard_aggregates(final_payload)
    E55._validate_screen_shard_records(final_payload)
    return final_payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("transport", "initial", "resolve"))
    parser.add_argument("--lattice", choices=("15x9", "27x5"))
    parser.add_argument("--tries", type=int, default=2_000)
    args = parser.parse_args()
    if args.command == "transport":
        print(json.dumps(transport_payload(), indent=1))
        return 0
    if args.lattice is None:
        parser.error(f"{args.command} requires --lattice")
    lattice = tuple(int(value) for value in args.lattice.split("x"))
    payload = (
        screen_initial(lattice)
        if args.command == "initial"
        else resolve_screen(lattice, tries=args.tries)
    )
    print(
        json.dumps(
            {
                "lattice": [int(value) for value in lattice],
                "candidates": payload.get(
                    "candidates", payload.get("candidates_after_symmetry")
                ),
                "represented_pairs": payload.get(
                    "represented_pairs", payload.get("orbit_total")
                ),
                "verdicts": (
                    dict(
                        sorted(
                            Counter(
                                record["verdict"]
                                for record in payload["records"]
                            ).items()
                        )
                    )
                ),
                "wall_s": payload["wall_s"],
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
