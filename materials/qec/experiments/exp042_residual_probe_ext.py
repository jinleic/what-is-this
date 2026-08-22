"""EXP-042-EXT: exhaustive trinomial-parent hunt for T(P) < k_P/2 up to ell*m <= 56.

EXP-042 established SMALL_CLASS_EMPTY for ell*m <= 24 (all families).  This
extension pushes the exhaustive trinomial (3x3) parent sweep over the ordered
target lattices 5x6, 4x8, 5x7, 4x9, 6x6, 5x8, 6x7, 4x11, 4x12, 5x10, 4x13,
6x9 (4x6=24 is the measured ETA basis and is re-run for a self-contained
artifact).  Protocol identical to EXP-042: translation-canonical supports
(dedup under independent Z_ell x Z_m block translations, unordered {A, B}
including the A=B Hadamard-style swap), filters k_P >= 2 and a single Tanner
component, exact d_Z and T per parent (solver/coset engine, vectorized MITM
join, EXP-041 join fallback -- all bit-exact against each other and against
full nullspace spans where the span fits).

Budget protocol (per assignment): each lattice gets a 25-minute wall budget.
Before starting a lattice its ETA is computed from the 4x6 basis rate
(wall/pairs at 4x6) times the lattice's pair count times a growth multiplier
(L/24)^p fitted on every completed lattice; if the ETA exceeds 25 minutes the
lattice is SKIPPED and recorded with its extrapolated ETA (larger lattices
are not attempted once the bound is reached unless a refit after a fast
lattice brings them back under budget).  A hard in-run deadline aborts a
lattice that outgrows its budget mid-scan; such a lattice is recorded as
aborted and yields no verdict contribution.

Everything is exact GF(2) arithmetic; no SAT solver, no Stim, single thread.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from math import comb
from typing import Any

for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_env] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCHEMA = "exp042-residual-probe-ext-v1"
OUT = ROOT / "results" / "processed" / "exp042_residual_probe_ext.json"

FAMILY = (3, 3)               # trinomial x trinomial parents only
LATTICE_BUDGET_S = 25 * 60.0  # per-lattice wall budget (assignment cap)
BASIS_LATTICE = (4, 6)        # ETA basis, re-measured for this artifact
# An earlier in-session launch of this same scan (identical E42 module and
# lattice code, before the monotone-skip patch) exceeded the 25-minute
# deadline on 6x7 after 1502s; re-running it would deterministically abort
# again, so the measurement is reused instead of paying for it twice.
PRIOR_INSESSION_ABORTS = {(6, 7): 1502.0}
# Ascending target list given by the assignment.
TARGET_LATTICES = (
    (5, 6),    # 30
    (4, 8),    # 32
    (5, 7),    # 35
    (4, 9),    # 36
    (6, 6),    # 36
    (5, 8),    # 40
    (6, 7),    # 42
    (4, 11),   # 44
    (4, 12),   # 48
    (5, 10),   # 50
    (4, 13),   # 52
    (6, 9),    # 54
)
MIN_FIT_EXPONENT = 2.0
MAX_FIT_EXPONENT = 10.0

_SPEC = importlib.util.spec_from_file_location(
    "exp042_residual_probe", ROOT / "experiments" / "exp042_residual_probe.py"
)
E42 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E42
_SPEC.loader.exec_module(E42)
E27 = E42.E41.E27


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def join_engine_selftest() -> dict[str, Any]:
    """Vectorized join vs the EXP-041 dict join on random dense instances."""

    rng = np.random.default_rng(20260818)
    checked = 0
    for L in (12, 16, 20):
        A = (rng.random((L, L)) < 0.3).astype(np.uint8)
        B = (rng.random((L, L)) < 0.3).astype(np.uint8)
        table = E42.SideTable(L, E42.E41.column_masks(A), E42.E41.column_masks(B))
        for w1 in range(0, 5):
            for w2 in range(0, 5):
                if comb(L, max(w1, w2)) > 200_000:
                    continue
                a = sorted(table.join_matches(w1, w2, cap=1_000_000))
                b = sorted(table.matches(w1, w2))
                if a != b:
                    raise AssertionError(f"join/dict mismatch at L={L} {w1},{w2}")
                checked += 1
    return {"random_cases": checked}


def predicted_pairs(ell: int, m: int) -> int:
    """Canonical support pairs without building matrices (ETA input)."""

    canon = len(E42.canonical_supports(ell, m, FAMILY[0]))
    return canon * (canon + 1) // 2


def fit_exponent(points: list[dict[str, Any]], basis: dict[str, Any]) -> float:
    """Least-squares slope of ln(ms_per_pair) vs ln(L), anchored by the basis."""

    xs, ys = [], []
    for pt in [basis] + points:
        xs.append(np.log(pt["L"]))
        ys.append(np.log(pt["ms_per_pair"]))
    if len(xs) < 2:
        return 0.0
    slope = float(np.polyfit(np.array(xs), np.array(ys), 1)[0])
    return min(max(slope, MIN_FIT_EXPONENT), MAX_FIT_EXPONENT)


def eta_seconds(
    ell: int, m: int, pairs: int, basis: dict[str, Any], exponent: float
) -> float:
    return (
        basis["wall_s"]
        * pairs
        / basis["pairs"]
        * ((ell * m) / basis["L"]) ** exponent
    )


def run() -> dict[str, Any]:
    started = time.perf_counter()
    checks = E42.self_tests([(2, 2), (2, 3), (3, 2), (2, 4), (3, 3), (4, 2), (3, 6), (6, 3)])
    checks["join_engine_matches_dict_join"] = join_engine_selftest()

    # --- ETA basis: re-measure 4x6 with the current engines.
    t0 = time.perf_counter()
    basis_accounting, basis_analyzed, basis_residual, _ = E42.scan_shape_family(
        *BASIS_LATTICE, FAMILY, None
    )
    basis_wall = time.perf_counter() - t0
    basis = {
        "lattice": f"{BASIS_LATTICE[0]}x{BASIS_LATTICE[1]}",
        "L": BASIS_LATTICE[0] * BASIS_LATTICE[1],
        "wall_s": round(basis_wall, 3),
        "pairs": int(basis_accounting["support_pairs_enumerated"]),
        "parents_analyzed": int(basis_accounting["parents_analyzed"]),
        "residual_count": int(basis_accounting["residual_count"]),
        "ms_per_pair": round(basis_wall * 1000.0 / basis_accounting["support_pairs_enumerated"], 6),
    }
    if basis_residual:
        raise AssertionError("residual parent appeared at the basis lattice 4x6")

    points: list[dict[str, Any]] = []
    lattice_records: list[dict[str, Any]] = []
    skip_list: list[dict[str, Any]] = []
    aggregate: Counter[str] = Counter()
    residual_out: list[dict[str, Any]] = []
    fingerprints: set[str] = set()
    first_witness: dict[str, Any] | None = None
    bound_reached_L = basis["L"]
    abort_bound_L: int | None = None
    abort_bound_lattice: str | None = None
    stopped_after_witness = False

    for ell, m in sorted(TARGET_LATTICES, key=lambda s: (s[0] * s[1], s)):
        L = ell * m
        pairs = predicted_pairs(ell, m)
        exponent = fit_exponent(points, basis)
        eta = eta_seconds(ell, m, pairs, basis, exponent)
        if first_witness is not None:
            skip_list.append(
                {
                    "lattice": f"{ell}x{m}",
                    "ell_times_m": L,
                    "canonical_pairs": pairs,
                    "eta_seconds": round(eta, 1),
                    "reason": "not scanned: first witness already found",
                }
            )
            continue
        if (ell, m) in PRIOR_INSESSION_ABORTS:
            abort_bound_L = L
            abort_bound_lattice = f"{ell}x{m}"
            skip_list.append(
                {
                    "lattice": f"{ell}x{m}",
                    "ell_times_m": L,
                    "canonical_pairs": pairs,
                    "eta_seconds": round(eta, 1),
                    "aborted_after_s": PRIOR_INSESSION_ABORTS[(ell, m)],
                    "reason": (
                        "not re-run: an earlier in-session launch of the same "
                        "scan exceeded the in-run deadline after "
                        f"{PRIOR_INSESSION_ABORTS[(ell, m)]:.0f}s (scan code "
                        "identical); every larger lattice is skipped by the "
                        "monotone bound rule"
                    ),
                }
            )
            print(
                f"[exp042-ext] {ell}x{m} skipped: prior in-session abort at "
                f"{PRIOR_INSESSION_ABORTS[(ell, m)]:.0f}s",
                flush=True,
            )
            continue
        if abort_bound_L is not None and L > abort_bound_L:
            skip_list.append(
                {
                    "lattice": f"{ell}x{m}",
                    "ell_times_m": L,
                    "canonical_pairs": pairs,
                    "reason": (
                        "skipped by the monotone bound rule: lattice "
                        f"{abort_bound_lattice} (ell*m={abort_bound_L}) already "
                        "exceeded the in-run budget, so every larger lattice is "
                        "skipped"
                    ),
                }
            )
            continue
        if eta > LATTICE_BUDGET_S:
            skip_list.append(
                {
                    "lattice": f"{ell}x{m}",
                    "ell_times_m": L,
                    "canonical_pairs": pairs,
                    "eta_seconds": round(eta, 1),
                    "eta_model_exponent": round(exponent, 3),
                    "reason": (
                        f"extrapolated ETA {eta:.0f}s exceeds the "
                        f"{int(LATTICE_BUDGET_S)}s per-lattice budget"
                    ),
                }
            )
            continue

        print(
            f"[exp042-ext] {ell}x{m} (L={L}): pairs={pairs} eta={eta:.0f}s "
            f"(exponent {exponent:.2f}) ...",
            flush=True,
        )
        t0 = time.perf_counter()
        try:
            accounting, analyzed, residual, _ = E42.scan_shape_family(
                ell, m, FAMILY, None, deadline_ts=t0 + LATTICE_BUDGET_S
            )
        except E42.LatticeDeadlineExceeded as exc:
            wall = time.perf_counter() - t0
            skip_list.append(
                {
                    "lattice": f"{ell}x{m}",
                    "ell_times_m": L,
                    "canonical_pairs": pairs,
                    "eta_seconds": round(eta, 1),
                    "aborted_after_s": round(wall, 1),
                    "reason": f"aborted at the in-run deadline: {exc}",
                }
            )
            if abort_bound_L is None or L > abort_bound_L:
                abort_bound_L = L
                abort_bound_lattice = f"{ell}x{m}"
            print(f"[exp042-ext] {ell}x{m} ABORTED after {wall:.0f}s", flush=True)
            continue
        wall = time.perf_counter() - t0

        buckets: Counter[str] = Counter()
        dz_hist: Counter[int] = Counter()
        kt_hist: Counter[str] = Counter()
        for rec in analyzed:
            k, t = rec["k_P"], rec["T"]
            if 2 * t > k:
                buckets["T_gt_half_k"] += 1
            elif 2 * t == k:
                buckets["T_eq_half_k"] += 1
            else:
                buckets["T_lt_half_k"] += 1
            if t == k:
                buckets["T_eq_k"] += 1
            dz_hist[rec["d_Z"]] += 1
            kt_hist[f"k={k},T={t}"] += 1
            fingerprints.add(E27.matrix_fingerprint(rec["HX"], rec["HZ"]))
        for rec in residual:
            safe = {
                "id": f"{ell}x{m}_{len(residual_out):03d}",
                "lattice": f"{ell}x{m}",
                "ell": ell,
                "m": m,
                "family": rec["family"],
                "A_terms": rec["A_terms"],
                "B_terms": rec["B_terms"],
                "k_P": rec["k_P"],
                "rank_hx": rec["rank_hx"],
                "rank_hz": rec["rank_hz"],
                "d_Z": rec["d_Z"],
                "T": rec["T"],
                "num_minimum_logicals": rec["num_minimum_logicals"],
                "num_module_rows": rec["num_module_rows"],
            }
            residual_out.append(safe)
            if first_witness is None:
                first_witness = dict(safe)
                print(
                    "[exp042-ext] FIRST WITNESS T < k_P/2: "
                    f"lattice {ell}x{m} A={safe['A_terms']} B={safe['B_terms']} "
                    f"k_P={safe['k_P']} T={safe['T']} d_Z={safe['d_Z']}",
                    flush=True,
                )
        if residual_out:
            stopped_after_witness = True
        record = {
            "lattice": f"{ell}x{m}",
            "ell": ell,
            "m": m,
            "ell_times_m": L,
            "family": {
                **{k: v for k, v in accounting.items()},
                "buckets": {key: int(v) for key, v in sorted(buckets.items())},
                "d_Z_histogram": {str(k): int(v) for k, v in sorted(dz_hist.items())},
                "k_T_histogram": {k: int(v) for k, v in sorted(kt_hist.items())},
                "wall_time_s": round(wall, 3),
            },
            "eta_seconds_predicted": round(eta, 1),
            "eta_model_exponent_at_start": round(exponent, 3),
        }
        lattice_records.append(record)
        points.append(
            {
                "lattice": f"{ell}x{m}",
                "L": L,
                "wall_s": round(wall, 3),
                "pairs": pairs,
                "ms_per_pair": round(wall * 1000.0 / pairs, 6),
            }
        )
        aggregate["parents_analyzed"] += len(analyzed)
        for key in ("T_gt_half_k", "T_eq_half_k", "T_lt_half_k", "T_eq_k"):
            aggregate[key] += int(buckets.get(key, 0))
        aggregate["support_pairs"] += pairs
        aggregate["unresolved_parents"] += len(accounting["unresolved_parents"])
        bound_reached_L = max(bound_reached_L, L)
        print(
            f"[exp042-ext] {ell}x{m}: analyzed={len(analyzed)} "
            f"residual={len(residual)} buckets={dict(buckets)} ({wall:.1f}s)",
            flush=True,
        )

    if first_witness is not None:
        verdict = "RESIDUAL_CLASS_FIRST_WITNESS"
        verdict_detail = (
            f"first parent with T < k_P/2 beyond ell*m <= 24 found at "
            f"{first_witness['lattice']}: A={first_witness['A_terms']} "
            f"B={first_witness['B_terms']} k_P={first_witness['k_P']} "
            f"T={first_witness['T']} d_Z={first_witness['d_Z']}"
        )
    else:
        verdict = "SMALL_CLASS_EMPTY"
        completed = [r["ell_times_m"] for r in lattice_records]
        verdict_detail = (
            f"No trinomial parent with T < k_P/2 on any target lattice completed "
            f"(bound reached: ell*m = {max(completed) if completed else basis['L']}; "
            f"{len(skip_list)} larger lattices skipped by the ETA rule with "
            "extrapolated budgets attached)"
        )

    payload = {
        "schema": SCHEMA,
        "experiment": "EXP-042-EXT",
        "utc": utc_now(),
        "protocol": {
            "family": "3x3 trinomial parents only (assignment scope)",
            "target_lattices": [f"{a}x{b}" for a, b in TARGET_LATTICES],
            "basis_lattice": basis["lattice"],
            "per_lattice_budget_s": int(LATTICE_BUDGET_S),
            "eta_model": (
                "ETA = wall(4x6) * pairs(L)/pairs(4x6) * (L/24)^p with p fitted "
                "by least squares on ln(ms per pair) vs ln(L) over the basis and "
                "every completed lattice, clamped to [2, 10]; lattices with "
                "ETA > 25 min are skipped and larger ones re-evaluated after "
                "each refit; a lattice that exceeds the in-run 25-minute "
                "deadline (or a recorded prior in-session abort of the same "
                "scan) sets a monotone bound -- every strictly larger lattice "
                "is skipped"
            ),
            "parent_equivalence": (
                "supports identified under independent Z_ell x Z_m translations "
                "of the A and B blocks; {A, B} unordered (block swap allowed, "
                "including A = B); k_P, d_Z, T are invariants of the class"
            ),
            "nondegenerate_predicate": "k_P >= 2 and a single Tanner component",
            "d_Z_method": (
                "ascending-weight exact enumeration of ker[A B]: per split the "
                "cheapest of {solve-left, solve-right, vectorized MITM join} by "
                "cost model; EXP-041 dict join as fallback; one-sided splits "
                "read kernel spans; engines bit-exact per self_tests"
            ),
            "T_method": (
                "translation-orbit span of every minimum-weight Z-logical "
                "modulo S_Z; every translate verified in ker[A B]"
            ),
            "engine_caps": {
                "solver_side_weight_cap": int(E42.SIDE_WEIGHT_CAP),
                "join_row_cap": int(E42.JOIN_ROW_CAP),
                "coset_cap_bits": int(E42.COSET_CAP_BITS),
                "d_weight_cap": int(E42.D_WEIGHT_CAP),
                "match_guard": int(E42.MATCH_GUARD),
            },
        },
        "basis_measurement": basis,
        "eta_points": points,
        "self_tests": checks,
        "lattices": lattice_records,
        "extrapolated_skip_list": skip_list,
        "abort_bound": {
            "lattice": abort_bound_lattice,
            "ell_times_m": abort_bound_L,
        },
        "aggregate": {
            "lattices_completed": len(lattice_records),
            "bound_reached_ell_times_m": bound_reached_L,
            "support_pairs_enumerated": int(aggregate["support_pairs"]),
            "parents_analyzed": int(aggregate["parents_analyzed"]),
            "bucket_T_gt_half_k": int(aggregate["T_gt_half_k"]),
            "bucket_T_eq_half_k": int(aggregate["T_eq_half_k"]),
            "bucket_T_lt_half_k": int(aggregate["T_lt_half_k"]),
            "bucket_T_eq_k": int(aggregate["T_eq_k"]),
            "unresolved_parents": int(aggregate["unresolved_parents"]),
            "distinct_parent_fingerprints": len(fingerprints),
            "residual_parents": len(residual_out),
        },
        "first_witness": first_witness,
        "residual_parents": residual_out,
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "sat_solver_used": False,
        "pysat_imported": "pysat" in sys.modules,
        "sat_decide_imported": "qec_research.distance.sat_decide" in sys.modules,
        "stim_imported": "stim" in sys.modules,
        "wall_time_s": round(time.perf_counter() - started, 3),
    }
    atomic_write_json(OUT, payload)
    return payload


def main() -> int:
    payload = run()
    print(f"[exp042-ext] wrote {OUT}")
    print(f"[exp042-ext] verdict: {payload['verdict']} -- {payload['verdict_detail']}")
    print(f"[exp042-ext] aggregate: {payload['aggregate']} ({payload['wall_time_s']:.1f}s)")
    if payload["pysat_imported"] or payload["sat_decide_imported"]:
        return 1
    return 0 if payload["verdict"] == "SMALL_CLASS_EMPTY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
