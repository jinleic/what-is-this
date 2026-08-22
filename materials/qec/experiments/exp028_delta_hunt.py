"""EXP-028: adversarial search for a distance-increasing PBB perturbation.

The search is deliberately outside the published PBB catalogue.  It samples
translation-deduplicated general trinomial parents, enumerates commuting
perturbations of total weight at most three, and excludes every candidate in
the safe two-block translation orbit of a catalogue row.

A hit is deliberately strict: both the CSS parent and non-CSS child distances
must be reported OPTIMAL by the repository's exact CP-SAT solvers, and the
result must survive a fresh rebuild, independent GF(2)/witness checks, and a
second solver run with seed + 1.  This is a budgeted search, not an exhaustive
statement about any lattice.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import multiprocessing as mp
import os
import platform
import queue
import sys
import time
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Sequence
# Prevent each spawned CP-SAT worker process from recursively parallelising
# BLAS/Accelerate work on the shared machine.  These must precede NumPy import.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_env] = "1"


import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec,
    PBBSpec,
    bb_stabilizer,
    build_bb,
    build_pbb,
    commutation_defect,
    monomial_matrix,
    poly_matrix,
)
from qec_research.codes.pbb_theory import analyse_pbb  # noqa: E402
from qec_research.distance.exact import (  # noqa: E402
    exact_distance_css,
    exact_distance_symplectic,
)
from qec_research.gf2.linalg import rank_bitset, rank_np, rows_to_bitsets  # noqa: E402
from qec_research.symplectic.core import symplectic_weight  # noqa: E402

# ---------------------------------------------------------------------------
# Pre-registered protocol constants.  The canonical run uses these unchanged.
EXP_ID = 28
SEED = 20260812
VERIFY_SEED = SEED + 1
LATTICES = ((6, 6), (6, 4), (4, 6))
TOTAL_PARENT_SAMPLE = 300
PARENTS_PER_LATTICE = TOTAL_PARENT_SAMPLE // len(LATTICES)
A_WEIGHT = 3
B_WEIGHT = 3
K_MIN = 4
K_MAX = 16
MAX_PERTURBATION_WEIGHT = 3
PERTURBATIONS_PER_PARENT_CAP = 200
SOLVER_TIME_LIMIT_S = 240.0
SOLVER_WALL_LIMIT_S = 240.0
SOLVER_WORKERS = 8
TOTAL_COMPUTE_BUDGET_MIN = 90.0
TOTAL_COMPUTE_BUDGET_S = 60.0 * TOTAL_COMPUTE_BUDGET_MIN
HIT_VERIFICATION_RESERVE_S = 2.0 * SOLVER_WALL_LIMIT_S
CATALOGUE = (
    ROOT
    / "third_party"
    / "qcode-discovery"
    / "results"
    / "campaign7_publication_merged.jsonl"
)
OUTPUT = ROOT / "results" / "processed" / "exp028_delta_hunt.json"

Term = tuple[int, int]
Terms = tuple[Term, ...]


def _json_terms(terms: Iterable[Sequence[int]]) -> list[list[int]]:
    return [[int(a), int(b)] for a, b in terms]


def _normal_terms(terms: Iterable[Sequence[int]], ell: int, m: int) -> Terms:
    """Order-independent polynomial support, reduced modulo the lattice."""
    return tuple(sorted((int(a) % ell, int(b) % m) for a, b in terms))


def _shift_terms(terms: Terms, ell: int, m: int, da: int, db: int) -> Terms:
    return tuple(sorted(((a + da) % ell, (b + db) % m) for a, b in terms))


def polynomial_translation_key(
    ell: int, m: int, terms: Iterable[Sequence[int]]
) -> Terms:
    """Canonical support of one polynomial under multiplication by a monomial."""
    normal = _normal_terms(terms, ell, m)
    return min(
        _shift_terms(normal, ell, m, da, db)
        for da in range(ell)
        for db in range(m)
    )


def parent_translation_key(
    ell: int,
    m: int,
    A: Iterable[Sequence[int]],
    B: Iterable[Sequence[int]],
) -> tuple[Terms, Terms]:
    """Canonical BB parent under independent translations of its two blocks."""
    return (
        polynomial_translation_key(ell, m, A),
        polynomial_translation_key(ell, m, B),
    )


def _paired_translation_key(
    ell: int,
    m: int,
    first: Iterable[Sequence[int]],
    second: Iterable[Sequence[int]],
) -> tuple[Terms, Terms]:
    first_normal = _normal_terms(first, ell, m)
    second_normal = _normal_terms(second, ell, m)
    return min(
        (
            _shift_terms(first_normal, ell, m, da, db),
            _shift_terms(second_normal, ell, m, da, db),
        )
        for da in range(ell)
        for db in range(m)
    )


def pbb_translation_key(
    ell: int,
    m: int,
    A: Iterable[Sequence[int]],
    B: Iterable[Sequence[int]],
    C: Iterable[Sequence[int]],
    D: Iterable[Sequence[int]],
) -> tuple[tuple[Terms, Terms], tuple[Terms, Terms]]:
    """Canonical PBB key under the two safe block translations.

    Translating A and C together permutes the left data block; translating B
    and D together independently permutes the right block.  The corresponding
    bottom-check permutation preserves the stabilizer code.  Excluding this
    full orbit is stronger than excluding only literal catalogue rows.
    """
    return (
        _paired_translation_key(ell, m, A, C),
        _paired_translation_key(ell, m, B, D),
    )


def canonical_trinomials(ell: int, m: int) -> list[Terms]:
    """All weight-three supports, deduplicated under lattice translation.

    This generalises EXP-018's ``canonical_A`` idea beyond its axis-restricted
    Bravyi family; that restriction has no positive-k parents on 6x4 or 4x6.
    """
    positions = tuple((a, b) for a in range(ell) for b in range(m))
    keys = {
        polynomial_translation_key(ell, m, support)
        for support in itertools.combinations(positions, A_WEIGHT)
    }
    return sorted(keys)


def canonical_A(ell: int, m: int) -> list[Terms]:
    return canonical_trinomials(ell, m)


def canonical_B(ell: int, m: int) -> list[Terms]:
    return canonical_trinomials(ell, m)


def catalogue_summary() -> tuple[dict, set[tuple]]:
    rows = [json.loads(line) for line in CATALOGUE.open() if line.strip()]
    total_hist: Counter[int] = Counter()
    joint_hist: Counter[tuple[int, int]] = Counter()
    per_lattice: dict[str, Counter[int]] = {
        f"{ell}x{m}": Counter() for ell, m in LATTICES
    }
    orbit_keys: set[tuple] = set()
    max_c = max_d = max_total = 0
    valid_term_rows = 0
    for row in rows:
        if not all(row.get(name) is not None for name in ("A_terms", "B_terms", "C_terms", "D_terms")):
            continue
        valid_term_rows += 1
        wc, wd = len(row["C_terms"]), len(row["D_terms"])
        wt = wc + wd
        total_hist[wt] += 1
        joint_hist[(wc, wd)] += 1
        max_c, max_d, max_total = max(max_c, wc), max(max_d, wd), max(max_total, wt)
        lattice = f"{int(row['ell'])}x{int(row['m'])}"
        if lattice in per_lattice:
            per_lattice[lattice][wt] += 1
        orbit_keys.add(
            (
                int(row["ell"]),
                int(row["m"]),
                pbb_translation_key(
                    int(row["ell"]),
                    int(row["m"]),
                    row["A_terms"],
                    row["B_terms"],
                    row["C_terms"],
                    row["D_terms"],
                ),
            )
        )
    summary = {
        "path": str(CATALOGUE.relative_to(ROOT)),
        "rows": len(rows),
        "rows_with_all_term_fields": valid_term_rows,
        "field_weight_definition": "len(C_terms) + len(D_terms), exactly as stored",
        "C_weight_max": max_c,
        "D_weight_max": max_d,
        "total_weight_max": max_total,
        "total_weight_histogram": {str(k): v for k, v in sorted(total_hist.items())},
        "joint_C_D_weight_histogram": {
            f"{wc},{wd}": count for (wc, wd), count in sorted(joint_hist.items())
        },
        "rows_total_weight_le_2": sum(v for w, v in total_hist.items() if w <= 2),
        "rows_total_weight_eq_3": total_hist[3],
        "rows_total_weight_gt_3": sum(v for w, v in total_hist.items() if w > 3),
        "target_lattice_total_weight_histograms": {
            lattice: {str(k): v for k, v in sorted(hist.items())}
            for lattice, hist in per_lattice.items()
        },
        "belief_total_weight_le_2_is_true": max_total <= 2,
        "weight_3_is_new_relative_to_catalogue": total_hist[3] == 0,
        "finding": (
            f"The catalogue does not stop at weight 2: {total_hist[3]} rows have total "
            f"perturbation weight 3 and the observed maximum is {max_total}."
        ),
    }
    return summary, orbit_keys


def enumerate_parent_universe(ell: int, m: int) -> tuple[list[dict], dict]:
    """Enumerate, translation-deduplicate, and filter canonical parent pairs."""
    As, Bs = canonical_A(ell, m), canonical_B(ell, m)
    seen: set[tuple] = set()
    unique_pairs: list[tuple[Terms, Terms]] = []
    for A in As:
        for B in Bs:
            key = parent_translation_key(ell, m, A, B)
            if key in seen:
                continue
            seen.add(key)
            unique_pairs.append((A, B))

    eligible: list[dict] = []
    k_hist: Counter[int] = Counter()
    connected_count = 0
    for A, B in unique_pairs:
        spec = BBSpec(ell=ell, m=m, A=list(A), B=list(B))
        code = bb_stabilizer(spec)
        rank = rank_np(code.H)
        k = code.n - rank
        k_hist[k] += 1
        connected = len(code.connected_components()) == 1
        connected_count += int(connected)
        if connected and K_MIN <= k <= K_MAX:
            eligible.append(
                {
                    "A": _json_terms(A),
                    "B": _json_terms(B),
                    "n": code.n,
                    "rank_parent": int(rank),
                    "k_parent": int(k),
                    "connected_tanner": True,
                    "translation_key_sha256": hashlib.sha256(repr(parent_translation_key(ell, m, A, B)).encode()).hexdigest(),
                }
            )
    accounting = {
        "A_family_size": len(As),
        "B_family_size": len(Bs),
        "raw_cartesian_pairs": len(As) * len(Bs),
        "translation_unique_pairs": len(unique_pairs),
        "translation_duplicates_removed": len(As) * len(Bs) - len(unique_pairs),
        "connected_pairs": connected_count,
        "eligible_connected_k_4_to_16": len(eligible),
        "k_histogram_after_translation_dedup": {str(k): v for k, v in sorted(k_hist.items())},
        "parent_family": (
            "all three-element supports for both A and B, each independently "
            "deduplicated under lattice translation; the general translation-"
            "canonical analogue of EXP-018 canonical_A"
        ),
        "connectivity_predicate": "one qubit component in BB StabilizerCode.connected_components()",
    }
    return eligible, accounting


def sample_parents() -> tuple[list[dict], dict]:
    sampled: list[dict] = []
    accounting: dict[str, dict] = {}
    for lattice_index, (ell, m) in enumerate(LATTICES):
        eligible, acc = enumerate_parent_universe(ell, m)
        lattice_seed = int(np.random.SeedSequence([SEED, lattice_index, ell, m]).generate_state(1)[0])
        rng = np.random.default_rng(lattice_seed)
        take = min(PARENTS_PER_LATTICE, len(eligible))
        chosen = rng.choice(len(eligible), size=take, replace=False).tolist() if take else []
        rows = [eligible[int(i)] for i in chosen]
        # This remains a uniform sample; ordering it by descending k only
        # spends the fixed compute budget in the delta>0-capable stratum first.
        rows.sort(key=lambda row: -int(row["k_parent"]))
        for sample_index, row in enumerate(rows):
            row.update(
                {
                    "ell": ell,
                    "m": m,
                    "lattice": f"{ell}x{m}",
                    "sample_index_within_lattice": sample_index,
                    "parent_sample_seed": lattice_seed,
                    "search_status": "sampled_not_reached",
                }
            )
            sampled.append(row)
        acc.update(
            {
                "parent_sample_seed": lattice_seed,
                "sample_target": PARENTS_PER_LATTICE,
                "parents_sampled": take,
                "sample_processing_order": "descending k within the uniform random sample",
            }
        )
        accounting[f"{ell}x{m}"] = acc
    return sampled, accounting


def _defect_signature(defect: np.ndarray) -> int:
    packed = np.packbits(np.asarray(defect, dtype=np.uint8).reshape(-1), bitorder="little")
    return int.from_bytes(packed.tobytes(), "little")


def singleton_defect_signatures(
    ell: int, m: int, A_terms: Terms, B_terms: Terms
) -> tuple[list[int], np.ndarray, np.ndarray]:
    """Linear commutation-defect signature for each C/D monomial coordinate."""
    dim = ell * m
    A = poly_matrix(ell, m, list(A_terms))
    B = poly_matrix(ell, m, list(B_terms))
    zero = np.zeros((dim, dim), dtype=np.uint8)
    signatures: list[int] = []
    for coordinate in range(2 * dim):
        monomial = monomial_matrix(ell, m, (coordinate % dim) // m, (coordinate % dim) % m)
        if coordinate < dim:
            defect = commutation_defect(A, B, monomial, zero)
        else:
            defect = commutation_defect(A, B, zero, monomial)
        signatures.append(_defect_signature(defect))
    return signatures, A, B


def combo_to_terms(combo: Sequence[int], ell: int, m: int) -> tuple[Terms, Terms]:
    dim = ell * m
    C = tuple(sorted((j // m, j % m) for j in combo if j < dim))
    D = tuple(sorted(((j - dim) // m, (j - dim) % m) for j in combo if j >= dim))
    return C, D


def enumerate_perturbations(
    parent: dict, catalogue_orbits: set[tuple], perturbation_seed: int
) -> tuple[list[dict], dict]:
    """Enumerate all supports of weight <= 3, reservoir-sampling at the cap."""
    ell, m = int(parent["ell"]), int(parent["m"])
    A_terms = _normal_terms(parent["A"], ell, m)
    B_terms = _normal_terms(parent["B"], ell, m)
    signatures, A, B = singleton_defect_signatures(ell, m, A_terms, B_terms)
    dim = ell * m
    zero = np.zeros((dim, dim), dtype=np.uint8)
    rng = np.random.default_rng(perturbation_seed)
    reservoir: list[dict] = []
    valid_hist: Counter[int] = Counter()
    outside_hist: Counter[int] = Counter()
    catalogue_excluded_hist: Counter[int] = Counter()
    duplicate_orbits = 0
    seen_orbits: set[tuple] = set()
    valid = outside = 0

    for weight in range(1, MAX_PERTURBATION_WEIGHT + 1):
        for combo in itertools.combinations(range(2 * dim), weight):
            signature = 0
            for coordinate in combo:
                signature ^= signatures[coordinate]
            if signature:
                continue
            valid += 1
            valid_hist[weight] += 1
            C_terms, D_terms = combo_to_terms(combo, ell, m)

            # Required direct predicate check, independent of the XOR signature path.
            C = poly_matrix(ell, m, list(C_terms)) if C_terms else zero
            D = poly_matrix(ell, m, list(D_terms)) if D_terms else zero
            if commutation_defect(A, B, C, D).any():
                raise AssertionError("signature accepted a noncommuting perturbation")

            orbit = pbb_translation_key(ell, m, A_terms, B_terms, C_terms, D_terms)
            if orbit in seen_orbits:
                duplicate_orbits += 1
                continue
            seen_orbits.add(orbit)
            if (ell, m, orbit) in catalogue_orbits:
                catalogue_excluded_hist[weight] += 1
                continue

            outside += 1
            outside_hist[weight] += 1
            candidate = {
                "C": _json_terms(C_terms),
                "D": _json_terms(D_terms),
                "perturbation_weight": weight,
                "support_coordinates": [int(x) for x in combo],
                "catalogue_block_translation_match": False,
                "orbit_sha256": hashlib.sha256(repr(orbit).encode()).hexdigest(),
            }
            if len(reservoir) < PERTURBATIONS_PER_PARENT_CAP:
                reservoir.append(candidate)
            else:
                slot = int(rng.integers(0, outside))
                if slot < PERTURBATIONS_PER_PARENT_CAP:
                    reservoir[slot] = candidate

    rng.shuffle(reservoir)
    accounting = {
        "perturbation_seed": perturbation_seed,
        "raw_supports_enumerated": sum(
            math.comb(2 * dim, weight)
            for weight in range(1, MAX_PERTURBATION_WEIGHT + 1)
        ),
        "commutation_valid_supports": valid,
        "commutation_valid_weight_histogram": {str(k): v for k, v in sorted(valid_hist.items())},
        "translation_duplicate_orbits_removed": duplicate_orbits,
        "catalogue_orbits_excluded": sum(catalogue_excluded_hist.values()),
        "catalogue_excluded_weight_histogram": {
            str(k): v for k, v in sorted(catalogue_excluded_hist.items())
        },
        "outside_catalogue_unique_perturbations": outside,
        "outside_catalogue_weight_histogram": {
            str(k): v for k, v in sorted(outside_hist.items())
        },
        "selected_under_cap": len(reservoir),
        "selected_weight_histogram": dict(
            sorted(Counter(int(r["perturbation_weight"]) for r in reservoir).items())
        ),
        "selection_method": (
            "uniform reservoir sample without replacement from translation-unique, "
            "commutation-valid supports outside catalogue safe block-translation orbits"
        ),
        "mixed_check_filter": (
            "verified after build: at least one generator row has nonzero support in "
            "both X and Z sectors"
        ),
    }
    return reservoir, accounting


@contextmanager
def cpsat_random_seed(seed: int):
    """Inject CP-SAT's random_seed into repository exact-solver calls."""
    original = cp_model.CpSolver

    def factory():
        solver = original()
        solver.parameters.random_seed = int(seed)
        return solver

    cp_model.CpSolver = factory
    try:
        yield
    finally:
        cp_model.CpSolver = original


def _aggregate_css_status(result: dict) -> str:
    return "OPTIMAL" if result.get("d_exact") else "PARTIAL"


def _distance_worker(task: dict, out_queue) -> None:
    """Fresh-process exact solve, enabling a hard 240 s wall-clock cap."""
    try:
        kind = task["kind"]
        ell, m = int(task["ell"]), int(task["m"])
        A, B = task["A"], task["B"]
        with cpsat_random_seed(int(task["solver_seed"])):
            if kind == "parent":
                HX, HZ = build_bb(BBSpec(ell=ell, m=m, A=A, B=B))
                result = exact_distance_css(
                    HX,
                    HZ,
                    time_limit_s=float(task["time_limit_s"]),
                    workers=int(task["workers"]),
                )
                result["status"] = _aggregate_css_status(result)
            elif kind == "pbb":
                code = build_pbb(
                    PBBSpec(
                        ell=ell,
                        m=m,
                        A=A,
                        B=B,
                        C=task["C"],
                        D=task["D"],
                    )
                )
                result = exact_distance_symplectic(
                    code,
                    time_limit_s=float(task["time_limit_s"]),
                    workers=int(task["workers"]),
                ).to_dict()
            else:
                raise ValueError(f"unknown distance task {kind}")
        out_queue.put({"ok": True, "result": result})
    except BaseException as exc:  # child must return an auditable failure record
        out_queue.put(
            {
                "ok": False,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
        )


def run_distance_task(task: dict, wall_limit_s: float) -> dict:
    """Run one exact API call with a hard outer wall cap and full status record."""
    ctx = mp.get_context("spawn")
    out_queue = ctx.Queue(maxsize=1)
    process = ctx.Process(target=_distance_worker, args=(task, out_queue))
    t0 = time.monotonic()
    process.start()
    process.join(timeout=max(0.1, wall_limit_s))
    elapsed = time.monotonic() - t0
    if process.is_alive():
        process.terminate()
        process.join(timeout=5.0)
        if process.is_alive():
            process.kill()
            process.join(timeout=5.0)
        return {
            "status": "WALL_TIMEOUT",
            "exact": False,
            "outer_wall_s": round(elapsed, 6),
            "outer_wall_limit_s": wall_limit_s,
            "solver_seed": int(task["solver_seed"]),
            "solver_time_limit_s_per_sector": float(task["time_limit_s"]),
        }
    try:
        message = out_queue.get(timeout=2.0)
    except queue.Empty:
        return {
            "status": "WORKER_NO_RESULT",
            "exact": False,
            "worker_exitcode": process.exitcode,
            "outer_wall_s": round(elapsed, 6),
            "solver_seed": int(task["solver_seed"]),
        }
    if not message["ok"]:
        return {
            "status": "WORKER_ERROR",
            "exact": False,
            "error_type": message["error_type"],
            "error": message["error"],
            "worker_exitcode": process.exitcode,
            "outer_wall_s": round(elapsed, 6),
            "solver_seed": int(task["solver_seed"]),
        }
    result = message["result"]
    result["outer_wall_s"] = round(elapsed, 6)
    result["outer_wall_limit_s"] = wall_limit_s
    result["solver_seed"] = int(task["solver_seed"])
    result["solver_time_limit_s_per_sector"] = float(task["time_limit_s"])
    result["solver_workers"] = int(task["workers"])
    return result


def _base_task(kind: str, parent: dict, seed: int, candidate: dict | None = None) -> dict:
    task = {
        "kind": kind,
        "ell": int(parent["ell"]),
        "m": int(parent["m"]),
        "A": parent["A"],
        "B": parent["B"],
        "solver_seed": int(seed),
        "time_limit_s": SOLVER_TIME_LIMIT_S,
        "workers": SOLVER_WORKERS,
    }
    if candidate is not None:
        task["C"], task["D"] = candidate["C"], candidate["D"]
    return task


def _matrix_sha256(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(matrix, dtype=np.uint8).tobytes()).hexdigest()


def verify_parent_witnesses(parent_code, result: dict) -> dict:
    n = parent_code.n
    checks = {}
    for sector in ("X", "Z"):
        support = result.get(f"d_{sector}_witness")
        if support is None:
            checks[sector] = {"present": False, "valid": False}
            continue
        vector = np.zeros(2 * n, dtype=np.uint8)
        offset = 0 if sector == "X" else n
        vector[offset + np.asarray(support, dtype=int)] = 1
        weight = symplectic_weight(vector)
        checks[sector] = {
            "present": True,
            "reported_distance": result.get(f"d_{sector}"),
            "independent_symplectic_weight": int(weight),
            "is_nontrivial_logical": bool(parent_code.is_logical(vector)),
            "valid": bool(
                parent_code.is_logical(vector)
                and weight == result.get(f"d_{sector}")
            ),
        }
    return checks


def verify_pbb_witness(code, result: dict) -> dict:
    vector_raw = result.get("witness_vector")
    if vector_raw is None:
        return {"present": False, "valid": False}
    vector = np.asarray(vector_raw, dtype=np.uint8)
    weight = symplectic_weight(vector)
    return {
        "present": True,
        "vector_length": len(vector),
        "reported_distance": result.get("value"),
        "independent_symplectic_weight": int(weight),
        "support_from_vector": np.flatnonzero(
            vector[: code.n] | vector[code.n :]
        ).astype(int).tolist(),
        "is_nontrivial_logical": bool(code.is_logical(vector)),
        "valid": bool(code.is_logical(vector) and weight == result.get("value")),
    }


def double_verify_hit(parent: dict, candidate: dict, initial_parent: dict, initial_pbb: dict) -> dict:
    """Fresh rebuild, independent arithmetic/witness checks, and seed+1 rerun."""
    ell, m = int(parent["ell"]), int(parent["m"])
    bb_spec = BBSpec(ell=ell, m=m, A=parent["A"], B=parent["B"])
    pbb_spec = PBBSpec(
        ell=ell,
        m=m,
        A=parent["A"],
        B=parent["B"],
        C=candidate["C"],
        D=candidate["D"],
    )
    parent_code = bb_stabilizer(bb_spec)
    pbb_code = build_pbb(pbb_spec, check=True)
    A = poly_matrix(ell, m, parent["A"])
    B = poly_matrix(ell, m, parent["B"])
    C = poly_matrix(ell, m, candidate["C"]) if candidate["C"] else np.zeros_like(A)
    D = poly_matrix(ell, m, candidate["D"]) if candidate["D"] else np.zeros_like(A)
    mixed = pbb_code.H[:, : pbb_code.n].any(axis=1) & pbb_code.H[:, pbb_code.n :].any(axis=1)

    parent_validation = parent_code.validate()
    pbb_validation = pbb_code.validate()
    independent = {
        "fresh_rebuild_from_terms": True,
        "parent_matrix_sha256": _matrix_sha256(parent_code.H),
        "pbb_matrix_sha256": _matrix_sha256(pbb_code.H),
        "commutation_defect_nonzero_entries": int(commutation_defect(A, B, C, D).sum()),
        "mixed_check_rows": int(mixed.sum()),
        "genuinely_non_css_mixed_predicate": bool(mixed.any()),
        "parent_validation": parent_validation,
        "pbb_validation": pbb_validation,
        "parent_rank_numpy": int(rank_np(parent_code.H)),
        "parent_rank_bitset": int(rank_bitset(rows_to_bitsets(parent_code.H), 2 * parent_code.n)),
        "pbb_rank_numpy": int(rank_np(pbb_code.H)),
        "pbb_rank_bitset": int(rank_bitset(rows_to_bitsets(pbb_code.H), 2 * pbb_code.n)),
        "initial_parent_witness_recheck": verify_parent_witnesses(parent_code, initial_parent),
        "initial_pbb_witness_recheck": verify_pbb_witness(pbb_code, initial_pbb),
    }

    rerun_parent = run_distance_task(
        _base_task("parent", parent, VERIFY_SEED), SOLVER_WALL_LIMIT_S
    )
    rerun_pbb = run_distance_task(
        _base_task("pbb", parent, VERIFY_SEED, candidate), SOLVER_WALL_LIMIT_S
    )
    rerun_parent_witness = verify_parent_witnesses(parent_code, rerun_parent)
    rerun_pbb_witness = verify_pbb_witness(pbb_code, rerun_pbb)
    rank_ok = (
        independent["parent_rank_numpy"] == independent["parent_rank_bitset"]
        and independent["pbb_rank_numpy"] == independent["pbb_rank_bitset"]
    )
    initial_witnesses_ok = (
        all(x["valid"] for x in independent["initial_parent_witness_recheck"].values())
        and independent["initial_pbb_witness_recheck"]["valid"]
    )
    rerun_witnesses_ok = (
        all(x["valid"] for x in rerun_parent_witness.values())
        and rerun_pbb_witness["valid"]
    )
    distances_agree = (
        rerun_parent.get("status") == "OPTIMAL"
        and rerun_pbb.get("status") == "OPTIMAL"
        and rerun_parent.get("d") == initial_parent.get("d")
        and rerun_pbb.get("value") == initial_pbb.get("value")
        and rerun_pbb.get("value", -1) > rerun_parent.get("d", 10**9)
    )
    verified = bool(
        independent["commutation_defect_nonzero_entries"] == 0
        and independent["genuinely_non_css_mixed_predicate"]
        and parent_validation["rank_agree"]
        and pbb_validation["rank_agree"]
        and parent_validation["commutes_numpy"]
        and parent_validation["commutes_bitset"]
        and pbb_validation["commutes_numpy"]
        and pbb_validation["commutes_bitset"]
        and rank_ok
        and initial_witnesses_ok
        and rerun_witnesses_ok
        and distances_agree
    )
    return {
        "verification_seed": VERIFY_SEED,
        "solver_seed_injection": "CP-SAT parameters.random_seed set in a fresh spawned process",
        "independent_checks": independent,
        "rerun_parent_distance": rerun_parent,
        "rerun_pbb_distance": rerun_pbb,
        "rerun_parent_witness_recheck": rerun_parent_witness,
        "rerun_pbb_witness_recheck": rerun_pbb_witness,
        "rank_paths_agree": rank_ok,
        "initial_witnesses_valid": initial_witnesses_ok,
        "rerun_witnesses_valid": rerun_witnesses_ok,
        "rerun_distances_agree_and_reverse": distances_agree,
        "verified": verified,
    }


def make_protocol() -> dict:
    return {
        "experiment": EXP_ID,
        "seed": SEED,
        "verification_seed": VERIFY_SEED,
        "lattices": [list(x) for x in LATTICES],
        "n_by_lattice": {f"{ell}x{m}": 2 * ell * m for ell, m in LATTICES},
        "total_parent_sample_target": TOTAL_PARENT_SAMPLE,
        "parents_per_lattice_target": PARENTS_PER_LATTICE,
        "parent_A_weight": A_WEIGHT,
        "parent_B_weight": B_WEIGHT,
        "parent_k_range_inclusive": [K_MIN, K_MAX],
        "parent_connected_tanner_required": True,
        "parent_translation_deduplication": "independent monomial translation orbits of A and B",
        "parent_processing_order": "descending k within each lattice's uniform sample",
        "max_perturbation_weight": MAX_PERTURBATION_WEIGHT,
        "perturbations_per_parent_cap": PERTURBATIONS_PER_PARENT_CAP,
        "commutation_predicate": "commutation_defect(A,B,C,D) is identically zero",
        "genuinely_non_css_predicate": (
            "exists a stabilizer row with nonzero X-sector support AND nonzero Z-sector support"
        ),
        "outside_catalogue_predicate": (
            "candidate (A,C)/(B,D) safe block-translation orbit absent from every catalogue row"
        ),
        "worker_thread_environment": {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        },
        "solver": "qec_research.distance.exact exact_distance_css and exact_distance_symplectic",
        "solver_status_required_for_hit": "OPTIMAL for parent and PBB",
        "solver_time_limit_s_per_sector": SOLVER_TIME_LIMIT_S,
        "solver_outer_wall_limit_s_per_api_call": SOLVER_WALL_LIMIT_S,
        "solver_workers": SOLVER_WORKERS,
        "compute_budget_minutes": TOTAL_COMPUTE_BUDGET_MIN,
        "hit_verification_reserve_s": HIT_VERIFICATION_RESERVE_S,
        "hit_rule": "both OPTIMAL and d_PBB > d_parent, followed by successful double verification",
        "early_stop": "first double-verified hit",
        "search_is_exhaustive": False,
        "candidate_priority": (
            "delta > 0 first (the only region not ruled out by Corollary 1), "
            "then larger delta, then deterministic reservoir order"
        ),
    }


def atomic_json_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=False))
    os.replace(temporary, path)


def routed_output_path(
    requested: Path,
    clean: bool,
    full_declared_scope: bool,
    stop_reason: str,
    timestamp: int | None = None,
) -> tuple[str, Path]:
    """Route noncanonical executions away from the canonical artifact."""
    route = canonical_route(clean, full_declared_scope)
    if route == "canonical" and requested == OUTPUT:
        return route, OUTPUT
    stamp = int(time.time()) if timestamp is None else int(timestamp)
    suffix = f"seed{SEED}-{stop_reason}-{stamp}.json"
    return route, ROOT / "results" / route / f"exp028_delta_hunt-{suffix}"

def format_coverage_statement(
    *,
    verified_hit: bool,
    parents_reached: int,
    perturbations_tested: int,
    per_lattice: dict[str, dict],
) -> str:
    """State aggregate coverage without multiplying or implying empty strata."""
    display = {f"{ell}x{m}": f"({ell},{m})" for ell, m in LATTICES}
    reached = [
        display[key]
        for key, counts in per_lattice.items()
        if counts["perturbations_solver_tested"] > 0
    ]
    unreached = [
        display[key]
        for key, counts in per_lattice.items()
        if counts["perturbations_solver_tested"] == 0
    ]
    outcome = "verified reversal found among" if verified_hit else "no reversal among"
    scope = f"solver-tested lattice(s) {', '.join(reached)}" if reached else "no solver-tested lattice"
    empty = (
        f"; zero solver-tested perturbations on {', '.join(unreached)}"
        if unreached
        else ""
    )
    return (
        f"{outcome} {perturbations_tested} perturbations across {parents_reached} "
        f"parents; {scope}{empty}; total perturbation weight <= 3 under budget "
        f"{TOTAL_COMPUTE_BUDGET_MIN:g} min"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help="output path; protocol constants and search budget are intentionally not overrideable",
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="check deterministic enumeration and predicates without running solvers or writing output",
    )
    return parser.parse_args()


def self_check() -> None:
    summary, catalogue_orbits = catalogue_summary()
    parents, accounting = sample_parents()
    assert summary["rows"] == 368
    assert len(parents) <= TOTAL_PARENT_SAMPLE
    assert sum(v["parents_sampled"] for v in accounting.values()) == len(parents)
    if not parents:
        raise AssertionError("parent sampler returned no eligible parent")
    parent = parents[0]
    seed = int(
        np.random.SeedSequence(
            [SEED, int(parent["ell"]), int(parent["m"]), int(parent["sample_index_within_lattice"])]
        ).generate_state(1)[0]
    )
    candidates, acc = enumerate_perturbations(parent, catalogue_orbits, seed)
    assert acc["selected_under_cap"] <= PERTURBATIONS_PER_PARENT_CAP
    for candidate in candidates:
        code = build_pbb(
            PBBSpec(
                ell=parent["ell"],
                m=parent["m"],
                A=parent["A"],
                B=parent["B"],
                C=candidate["C"],
                D=candidate["D"],
            )
        )
        mixed = code.H[:, : code.n].any(axis=1) & code.H[:, code.n :].any(axis=1)
        assert mixed.any()
    print(
        f"self-check passed: {len(parents)} sampled parents; first parent has "
        f"{acc['commutation_valid_supports']} valid supports and {len(candidates)} selected"
    )


def main() -> int:
    args = parse_args()
    if args.self_check:
        self_check()
        return 0

    started_wall = time.time()
    started_mono = time.monotonic()
    deadline = started_mono + TOTAL_COMPUTE_BUDGET_S
    protocol = make_protocol()
    catalogue, catalogue_orbits = catalogue_summary()
    parents, parent_space = sample_parents()
    evaluations: list[dict] = []
    parent_distance_results: list[dict] = []
    preliminary_hits: list[dict] = []
    verified_hits: list[dict] = []
    stop_reason = "sample_exhausted"

    print(
        f"EXP-028: {len(parents)} sampled parents on {LATTICES}; "
        f"budget {TOTAL_COMPUTE_BUDGET_MIN:.0f} min",
        flush=True,
    )
    print(f"catalogue finding: {catalogue['finding']}", flush=True)

    for parent_index, parent in enumerate(parents):
        remaining = deadline - time.monotonic()
        if remaining <= HIT_VERIFICATION_RESERVE_S + SOLVER_WALL_LIMIT_S:
            stop_reason = "compute_budget_reserve_reached_before_parent"
            break
        parent["search_status"] = "reached"
        perturbation_seed = int(
            np.random.SeedSequence(
                [
                    SEED,
                    int(parent["ell"]),
                    int(parent["m"]),
                    int(parent["sample_index_within_lattice"]),
                ]
            ).generate_state(1)[0]
        )
        candidates, perturbation_accounting = enumerate_perturbations(
            parent, catalogue_orbits, perturbation_seed
        )
        parent["perturbation_accounting"] = perturbation_accounting
        parent["search_status"] = "perturbations_enumerated"
        # Corollary 1 proves that delta=0 cannot reverse distance.  Preserve
        # those candidates in enumeration accounting, but spend solver budget
        # only on the genuinely open delta>0 region.
        candidate_meta = []
        for candidate in candidates:
            candidate_spec = PBBSpec(
                ell=int(parent["ell"]),
                m=int(parent["m"]),
                A=parent["A"],
                B=parent["B"],
                C=candidate["C"],
                D=candidate["D"],
            )
            candidate_structure = analyse_pbb(candidate_spec)
            candidate_meta.append((candidate, candidate_structure))
        candidates_delta_positive = [
            (candidate, structure)
            for candidate, structure in candidate_meta
            if structure.delta > 0 and structure.k_pbb > 0
        ]
        candidates_delta_positive.sort(key=lambda pair: -int(pair[1].delta))
        parent["selected_delta_histogram"] = {
            str(delta): count
            for delta, count in sorted(
                Counter(int(structure.delta) for _, structure in candidate_meta).items()
            )
        }
        parent["delta_zero_candidates_theory_skipped"] = sum(
            structure.delta == 0 for _, structure in candidate_meta
        )
        parent["delta_positive_candidates"] = len(candidates_delta_positive)
        parent["k_nonpositive_candidates"] = sum(
            structure.k_pbb <= 0 for _, structure in candidate_meta
        )
        if not candidates_delta_positive:
            parent["search_status"] = "no_delta_positive_positive_k_candidates"
            continue

        if not candidates:
            parent["search_status"] = "no_admissible_outside_catalogue_perturbations"
            continue

        parent_task = _base_task("parent", parent, SEED)
        parent_distance = run_distance_task(
            parent_task,
            min(SOLVER_WALL_LIMIT_S, max(0.1, deadline - time.monotonic() - HIT_VERIFICATION_RESERVE_S)),
        )
        parent_distance_id = len(parent_distance_results)
        parent_distance_results.append(
            {
                "parent_distance_id": parent_distance_id,
                "parent_index": parent_index,
                "lattice": parent["lattice"],
                "A": parent["A"],
                "B": parent["B"],
                "result": parent_distance,
            }
        )
        parent["parent_distance_id"] = parent_distance_id
        parent["parent_distance_status"] = parent_distance.get("status")
        parent["search_status"] = "parent_distance_attempted"
        if parent_distance.get("status") != "OPTIMAL":
            parent["search_status"] = "parent_distance_not_optimal"
            print(
                f"[{parent_index + 1}/{len(parents)}] {parent['lattice']} parent distance "
                f"{parent_distance.get('status')}; skipping its child solves",
                flush=True,
            )
            continue

        tested_this_parent = 0
        for candidate_index, (candidate, structure) in enumerate(candidates_delta_positive):
            remaining = deadline - time.monotonic()
            if remaining <= HIT_VERIFICATION_RESERVE_S + SOLVER_WALL_LIMIT_S:
                stop_reason = "compute_budget_reserve_reached_before_candidate"
                break

            spec = PBBSpec(
                ell=int(parent["ell"]),
                m=int(parent["m"]),
                A=parent["A"],
                B=parent["B"],
                C=candidate["C"],
                D=candidate["D"],
            )
            code = build_pbb(spec, check=True)
            mixed = code.H[:, : code.n].any(axis=1) & code.H[:, code.n :].any(axis=1)
            if not mixed.any():
                raise AssertionError("nonzero accepted perturbation failed the mixed-check predicate")
            rank = rank_np(code.H)
            k = code.n - rank
            if structure.k_pbb != k:
                raise AssertionError("analyse_pbb and direct rank disagree on k")
            record = {
                "evaluation_index": len(evaluations),
                "parent_index": parent_index,
                "parent_distance_id": parent_distance_id,
                "candidate_index_within_parent_sample": candidate_index,
                "lattice": parent["lattice"],
                "ell": int(parent["ell"]),
                "m": int(parent["m"]),
                "n": int(code.n),
                "A": parent["A"],
                "B": parent["B"],
                "C": candidate["C"],
                "D": candidate["D"],
                "perturbation_weight": int(candidate["perturbation_weight"]),
                "candidate_orbit_sha256": candidate["orbit_sha256"],
                "catalogue_block_translation_match": False,
                "commutation_defect_nonzero_entries": 0,
                "mixed_check_rows": int(mixed.sum()),
                "genuinely_non_css": True,
                "rank_pbb_numpy": int(rank),
                "k_pbb": int(k),
                "k_parent": int(parent["k_parent"]),
                "delta": int(structure.delta),
                "parent_distance_status": parent_distance.get("status"),
                "parent_distance": parent_distance.get("d"),
                "parent_distance_result_reused": True,
            }
            if k <= 0 or structure.delta <= 0:
                raise AssertionError("nonpositive-k or delta=0 candidate reached exact solver")

            pbb_distance = run_distance_task(
                _base_task("pbb", parent, SEED, candidate),
                min(
                    SOLVER_WALL_LIMIT_S,
                    max(0.1, deadline - time.monotonic() - HIT_VERIFICATION_RESERVE_S),
                ),
            )
            record["pbb_distance_result"] = pbb_distance
            record["solver_status"] = pbb_distance.get("status")
            record["pbb_distance"] = pbb_distance.get("value")
            record["both_solver_statuses_optimal"] = bool(
                parent_distance.get("status") == "OPTIMAL"
                and pbb_distance.get("status") == "OPTIMAL"
            )
            record["preliminary_reversal"] = bool(
                record["both_solver_statuses_optimal"]
                and pbb_distance.get("value") is not None
                and parent_distance.get("d") is not None
                and int(pbb_distance["value"]) > int(parent_distance["d"])
            )
            record["hit"] = False
            evaluations.append(record)
            tested_this_parent += 1

            if record["preliminary_reversal"]:
                preliminary = {
                    "parent_index": parent_index,
                    "evaluation_index": record["evaluation_index"],
                    "terms": {
                        "ell": int(parent["ell"]),
                        "m": int(parent["m"]),
                        "A": parent["A"],
                        "B": parent["B"],
                        "C": candidate["C"],
                        "D": candidate["D"],
                    },
                    "n": int(code.n),
                    "k_parent": int(parent["k_parent"]),
                    "k_pbb": int(k),
                    "delta": int(structure.delta),
                    "perturbation_weight": int(candidate["perturbation_weight"]),
                    "initial_parent_distance": parent_distance,
                    "initial_pbb_distance": pbb_distance,
                }
                verification = double_verify_hit(
                    parent, candidate, parent_distance, pbb_distance
                )
                preliminary["double_verification"] = verification
                preliminary_hits.append(preliminary)
                if verification["verified"]:
                    record["hit"] = True
                    record["hit_certificate_index"] = len(verified_hits)
                    verified_hits.append(preliminary)
                    stop_reason = "first_double_verified_hit"
                    break

        parent["perturbations_solver_tested"] = tested_this_parent
        parent["search_status"] = (
            "hit_found"
            if verified_hits
            else "delta_positive_candidates_completed"
            if tested_this_parent == len(candidates_delta_positive)
            else "candidate_search_stopped"
        )
        print(
            f"[{parent_index + 1}/{len(parents)}] {parent['lattice']} "
            f"k={parent['k_parent']} d_parent={parent_distance.get('d')} "
            f"tested={tested_this_parent}/{len(candidates_delta_positive)} delta>0 hits={len(verified_hits)}",
            flush=True,
        )
        if verified_hits or stop_reason.startswith("compute_budget"):
            break

    wall_s = time.monotonic() - started_mono
    tested_parent_indices = {
        int(row["parent_index"])
        for row in evaluations
        if row.get("solver_status") != "SKIPPED_K_NONPOSITIVE"
    }
    perturbations_tested = sum(
        1 for row in evaluations if row.get("solver_status") != "SKIPPED_K_NONPOSITIVE"
    )
    k_nonpositive_skips = sum(
        1 for row in evaluations if row.get("solver_status") == "SKIPPED_K_NONPOSITIVE"
    )
    optimal_pairs = sum(1 for row in evaluations if row.get("both_solver_statuses_optimal"))

    per_lattice: dict[str, dict] = {}
    for ell, m in LATTICES:
        lattice = f"{ell}x{m}"
        parent_indices = {i for i, p in enumerate(parents) if p["lattice"] == lattice}
        lattice_evals = [row for row in evaluations if row["lattice"] == lattice]
        per_lattice[lattice] = {
            "n": 2 * ell * m,
            "eligible_parent_universe": parent_space[lattice]["eligible_connected_k_4_to_16"],
            "parents_sampled": parent_space[lattice]["parents_sampled"],
            "parents_reached": sum(
                1 for i in parent_indices if parents[i]["search_status"] != "sampled_not_reached"
            ),
            "parents_with_distance_attempted": sum(
                1 for row in parent_distance_results if row["lattice"] == lattice
            ),
            "parents_with_child_solver_test": len(
                {
                    int(row["parent_index"])
                    for row in lattice_evals
                    if row.get("solver_status") != "SKIPPED_K_NONPOSITIVE"
                }
            ),
            "commutation_valid_supports_enumerated": sum(
                int(p.get("perturbation_accounting", {}).get("commutation_valid_supports", 0))
                for p in parents
                if p["lattice"] == lattice
            ),
            "outside_catalogue_unique_perturbations_enumerated": sum(
                int(p.get("perturbation_accounting", {}).get("outside_catalogue_unique_perturbations", 0))
                for p in parents
                if p["lattice"] == lattice
            ),
            "perturbations_selected_under_caps": sum(
                int(p.get("perturbation_accounting", {}).get("selected_under_cap", 0))
                for p in parents
                if p["lattice"] == lattice
            ),
            "perturbations_solver_tested": sum(
                1
                for row in lattice_evals
                if row.get("solver_status") != "SKIPPED_K_NONPOSITIVE"
            ),
            "k_nonpositive_skips": sum(
                1
                for row in lattice_evals
                if row.get("solver_status") == "SKIPPED_K_NONPOSITIVE"
            ),
            "both_solver_statuses_optimal": sum(
                1 for row in lattice_evals if row.get("both_solver_statuses_optimal")
            ),
            "verified_hits": sum(
                1 for hit in verified_hits if f"{hit['terms']['ell']}x{hit['terms']['m']}" == lattice
            ),
        }

    parent_count_word = len(tested_parent_indices)
    if verified_hits:
        verdict = "BREAKTHROUGH_CANDIDATE"
        verdict_reason = "distance-increasing perturbation found"
    else:
        verdict = "NEGATIVE"
        verdict_reason = "no distance-increasing perturbation found in reached budget"
    coverage_statement = format_coverage_statement(
        verified_hit=bool(verified_hits),
        parents_reached=parent_count_word,
        perturbations_tested=perturbations_tested,
        per_lattice=per_lattice,
    )

    payload = {
        "experiment": "EXP-028",
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "coverage_statement": coverage_statement,
        "protocol": protocol,
        "catalogue_perturbation_weight_audit": catalogue,
        "parent_search_space": parent_space,
        "search_accounting": {
            "parents_sampled": len(parents),
            "parents_reached": sum(p["search_status"] != "sampled_not_reached" for p in parents),
            "parent_distance_calls": len(parent_distance_results),
            "parents_with_child_solver_test": parent_count_word,
            "perturbations_solver_tested": perturbations_tested,
            "k_nonpositive_skips": k_nonpositive_skips,
            "both_solver_statuses_optimal": optimal_pairs,
            "preliminary_reversals": len(preliminary_hits),
            "double_verified_hits": len(verified_hits),
            "per_lattice": per_lattice,
            "stop_reason": stop_reason,
            "budget_s": TOTAL_COMPUTE_BUDGET_S,
            "wall_s": wall_s,
            "budget_exhausted": wall_s >= TOTAL_COMPUTE_BUDGET_S,
            "budgeted_scope_completed": stop_reason in (
                "sample_exhausted",
                "first_double_verified_hit",
                "compute_budget_reserve_reached_before_parent",
                "compute_budget_reserve_reached_before_candidate",
            ),
            "exhaustive": False,
        },
        "sampled_parents": parents,
        "parent_distance_results": parent_distance_results,
        "candidate_evaluations": evaluations,
        "preliminary_hit_certificates": preliminary_hits,
        "hits": verified_hits,
        "seeds": {
            "base": SEED,
            "verification": VERIFY_SEED,
            "parent_and_perturbation_seeds": "derived deterministically with numpy SeedSequence",
            "solver_seed_injection": "CP-SAT parameters.random_seed in fresh spawned worker",
        },
        "reproduce_command": "PYTHONPATH=src .venv/bin/python experiments/exp028_delta_hunt.py",
        "started_unix_s": started_wall,
        "finished_unix_s": time.time(),
        "wall_s": wall_s,
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
            "cpu_count": os.cpu_count(),
            "shared_load_caveat": (
                "External interactive desktop workload dominated the run queue and "
                "QoS-deprioritized the experiment. This reduces coverage achieved "
                "inside the fixed 90-minute wall budget but does not affect exact "
                "CP-SAT statuses, witnesses, ranks, or candidate counts."
            ),
        },
        "scope_caveat": (
            "Budgeted random sample from the stated canonical parent family and capped "
            "perturbation samples only; never an exhaustive lattice-wide search."
        ),
    }
    clean = all(
        row.get("status") not in ("WORKER_ERROR", "WORKER_NO_RESULT")
        for row in (
            [entry["result"] for entry in parent_distance_results]
            + [entry["pbb_distance_result"] for entry in evaluations if "pbb_distance_result" in entry]
        )
    )
    full_declared_scope = stop_reason in (
        "sample_exhausted",
        "first_double_verified_hit",
        "compute_budget_reserve_reached_before_parent",
        "compute_budget_reserve_reached_before_candidate",
    )
    route, output_path = routed_output_path(
        args.output,
        clean=clean,
        full_declared_scope=full_declared_scope,
        stop_reason=stop_reason,
    )
    payload["artifact_routing"] = {
        "clean": clean,
        "full_declared_scope": full_declared_scope,
        "route": route,
        "canonical_path": str(OUTPUT.relative_to(ROOT)),
        "actual_path": str(output_path.relative_to(ROOT)),
        "canonical_protection": (
            "canonical output is written only by a clean run that completes the "
            "pre-registered budgeted scope: sample exhaustion, budget-reserve stop, "
            "or an early double-verified hit"
        ),
    }
    atomic_json_write(output_path, payload)
    print(f"{verdict}: {coverage_statement}", flush=True)
    print(f"wrote {output_path.relative_to(ROOT)} in {wall_s:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
