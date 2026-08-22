"""EXP-023 independent artifact cross-check.

This script does not import `exp023_light_gensets.py`. It rebuilds every code
from the catalogue, recomputes structural quantities with plain NumPy GF(2)
linear algebra, and re-checks every stored witness. It additionally cross-checks
the decisive decomposition through two independent calculations:

  * brute force over all row combinations of Hamming weight <= 3 (the regime a
    hand check can cover), confirming the enumerated X-codewords on that slice;
  * a rank-test reimplementation of the coset-budget decision for every
    published light X row, replacing the producer's bitset reducer.

The exhaustive, all-coefficient-space completeness certificate remains the
producer's CP-SAT `OPTIMAL` enumeration; this script does not duplicate it.

Output: results/raw/exp023_independent_verification.json
"""

from __future__ import annotations

import os

for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_env, "1")

import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB,
    PBBSpec,
    build_bb,
    build_pbb,
)
from qec_research.codes.pbb_theory import parent_bb_matrices  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402
from qec_research.symplectic.core import symplectic_weight  # noqa: E402

PROCESSED = ROOT / "results" / "processed" / "exp023_light_gensets.json"
RAW = ROOT / "results" / "raw" / "exp023_light_elements.json"
OUT = ROOT / "results" / "raw" / "exp023_independent_verification.json"
CATALOGUE = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
BRUTE_FORCE_MAX_ROW_COMBINATION = 3
EXPECTED_N = 144
EXPECTED_RANK = 132


def gf2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return ((np.asarray(a, np.uint16) & 1) @ (np.asarray(b, np.uint16) & 1) & 1).astype(
        np.uint8
    )


def bits_to_vector(bits: str) -> np.ndarray:
    return np.frombuffer(bits.encode(), dtype=np.uint8) - ord("0")


def weight(v: np.ndarray, n: int = EXPECTED_N) -> int:
    return int(np.count_nonzero(v[:n] | v[n:]))


def in_rowspace(rows: np.ndarray, v: np.ndarray) -> bool:
    base = rank_np(rows)
    return rank_np(np.vstack([rows, v])) == base


def rebuild(label: str, spec: dict[str, Any] | None) -> np.ndarray:
    if spec is None:
        hx, hz = build_bb(BRAVYI_BB["[[144,12,12]]"])
        zero = np.zeros_like(hx)
        return np.vstack([np.hstack([hx, zero]), np.hstack([zero, hz])]).astype(np.uint8)
    return build_pbb(
        PBBSpec(
            ell=spec["ell"],
            m=spec["m"],
            A=[tuple(t) for t in spec["A_terms"]],
            B=[tuple(t) for t in spec["B_terms"]],
            C=[tuple(t) for t in spec["C_terms"]],
            D=[tuple(t) for t in spec["D_terms"]],
            name=label,
        ),
        check=True,
    ).H


def pure_z_rank(spec: dict[str, Any]) -> tuple[int, np.ndarray]:
    pbb = PBBSpec(
        ell=spec["ell"],
        m=spec["m"],
        A=[tuple(t) for t in spec["A_terms"]],
        B=[tuple(t) for t in spec["B_terms"]],
        C=[tuple(t) for t in spec["C_terms"]],
        D=[tuple(t) for t in spec["D_terms"]],
    )
    hx, hz, perturbation = parent_bb_matrices(pbb)
    kernel = np.array(
        [row for row in _left_kernel(hx)], dtype=np.uint8
    ).reshape(-1, hx.shape[0])
    extra = gf2(kernel, perturbation) if kernel.shape[0] else np.zeros((0, EXPECTED_N), np.uint8)
    z_rows = np.vstack([hz, extra]).astype(np.uint8)
    return rank_np(z_rows), z_rows


def _left_kernel(matrix: np.ndarray) -> np.ndarray:
    """{u : u M = 0}, computed by an elimination written for this script."""
    m = (np.asarray(matrix, np.uint8) & 1).copy()
    rows = m.shape[0]
    tracker = np.eye(rows, dtype=np.uint8)
    pivot = 0
    for col in range(m.shape[1]):
        nz = np.flatnonzero(m[pivot:, col])
        if not nz.size:
            continue
        p = pivot + int(nz[0])
        if p != pivot:
            m[[pivot, p]] = m[[p, pivot]]
            tracker[[pivot, p]] = tracker[[p, pivot]]
        below = np.flatnonzero(m[pivot + 1 :, col]) + pivot + 1
        if below.size:
            m[below] ^= m[pivot]
            tracker[below] ^= tracker[pivot]
        pivot += 1
        if pivot == rows:
            break
    return tracker[pivot:]


def coset_fits_rank_test(
    z0: np.ndarray, z_rows: np.ndarray, support_x: np.ndarray, budget: int
) -> bool:
    """Independent reimplementation: membership decided by rank comparison."""
    outside = np.flatnonzero(~support_x.astype(bool))
    projected = z_rows[:, outside]
    base = rank_np(projected)
    target = z0[outside]
    if rank_np(np.vstack([projected, target])) == base:
        return True
    if budget <= 0:
        return False
    for extra in itertools.combinations(range(outside.size), budget):
        probe = target.copy()
        for position in extra:
            probe[position] ^= 1
        if rank_np(np.vstack([projected, probe])) == base:
            return True
    return False


def brute_force_light_x(hx: np.ndarray, cap: int, max_rows: int) -> dict[str, Any]:
    """All codewords u H_X with |u| <= max_rows and weight <= cap."""
    rows = hx.shape[0]
    found: dict[str, int] = {}
    singles = hx.copy()
    for size in range(1, max_rows + 1):
        for combination in itertools.combinations(range(rows), size):
            word = singles[combination[0]].copy()
            for index in combination[1:]:
                word = word ^ singles[index]
            w = int(word.sum())
            if 0 < w <= cap:
                found["".join(map(str, word.tolist()))] = w
    return {
        "max_row_combination": max_rows,
        "n_distinct_codewords": len(found),
        "weight_histogram": {
            str(w): sum(1 for value in found.values() if value == w)
            for w in sorted(set(found.values()))
        },
        "codewords": set(found),
    }


def main() -> None:
    started = time.perf_counter()
    processed = json.loads(PROCESSED.read_text())
    raw = json.loads(RAW.read_text())
    raw_by_code = {row["code"]: row for row in raw["results"]}
    catalogue = {}
    for line_number, line in enumerate(CATALOGUE.read_text().splitlines(), 1):
        row = json.loads(line)
        if row.get("n") == 144 and row.get("k") == 12 and row.get("d") == 12:
            key = row.get("code_id") or f"catalogue_row_{line_number:04d}_{row['bliss_hash']}"
            catalogue[key] = row

    report: dict[str, Any] = {
        "experiment": "EXP-023-independent-crosscheck",
        "verifies": str(PROCESSED.relative_to(ROOT)),
        "protocol": {
            "brute_force_max_row_combination": BRUTE_FORCE_MAX_ROW_COMBINATION,
            "expected_n": EXPECTED_N,
            "expected_rank": EXPECTED_RANK,
            "method": (
                "independent rebuild; NumPy GF(2) ranks; bounded brute-force "
                "row-combination sweep; rank-test coset cross-check"
            ),
        },
        "seeds": [],
        "seed_note": "Deterministic exhaustive/rank calculations; no randomness used.",
        "results": [],
    }
    all_ok = True

    for item in processed["results"]:
        label = item["code"]
        spec = catalogue.get(label)
        h = rebuild(label, spec)
        raw_item = raw_by_code[label]
        basis_rows = item["basis_row_indices"]
        h_basis = h[basis_rows]
        checks: dict[str, Any] = {
            "code": label,
            "kind": item["kind"],
            "rank_H_matches": rank_np(h) == item["rank_H"] == EXPECTED_RANK,
            "recorded_basis_is_independent": rank_np(h_basis) == EXPECTED_RANK,
        }

        # Every stored element: recompute v = c H_B and re-check weight/membership.
        element_ok = True
        n_elements = 0
        for w_key, per_w in raw_item["per_w"].items():
            for element in per_w.get("elements", []):
                n_elements += 1
                c = np.array(element["c"], dtype=np.uint8)
                v = gf2(c[None, :], h_basis)[0]
                if (
                    "".join(map(str, v.tolist())) != element["v_bits_x_then_z"]
                    or weight(v) != element["symplectic_weight"]
                    or weight(v) > int(w_key)
                    or symplectic_weight(v) != weight(v)
                    or not in_rowspace(h, v)
                ):
                    element_ok = False
        checks["n_stored_elements_rechecked"] = n_elements
        checks["all_stored_elements_valid"] = element_ok

        # Recorded ranks, recomputed from the stored vectors alone.
        rank_ok = True
        for w_key, per_w in raw_item["per_w"].items():
            vectors = [bits_to_vector(e["v_bits_x_then_z"]) for e in per_w.get("elements", [])]
            stack = (
                np.vstack(vectors)
                if vectors
                else np.zeros((0, 2 * EXPECTED_N), dtype=np.uint8)
            )
            recomputed = rank_np(stack)
            claimed = per_w["rank_Vw"]
            if per_w["closed"] and claimed is not None and recomputed != claimed:
                rank_ok = False
        checks["recorded_ranks_match_stored_vectors"] = rank_ok

        if item["kind"] == "PBB catalogue member":
            r_zsub, z_rows = pure_z_rank(spec)
            checks["r_Zsub_matches"] = r_zsub == item["pure_Z_subgroup"]["r_Zsub"]
            checks["r_Zsub"] = int(r_zsub)
            hx, _, perturbation = parent_bb_matrices(
                PBBSpec(
                    ell=spec["ell"],
                    m=spec["m"],
                    A=[tuple(t) for t in spec["A_terms"]],
                    B=[tuple(t) for t in spec["B_terms"]],
                    C=[tuple(t) for t in spec["C_terms"]],
                    D=[tuple(t) for t in spec["D_terms"]],
                )
            )
            brute = brute_force_light_x(hx, 7, BRUTE_FORCE_MAX_ROW_COMBINATION)
            enumerated = {
                "".join(map(str, row.tolist()))
                for row in hx
                if 0 < int(row.sum()) <= 7
            }
            claimed_histogram = item["theorem_gate"]["decisive_decomposition"][
                "x_codeword_enumeration"
            ]["weight_histogram"]
            checks["brute_force_weight_histogram"] = brute["weight_histogram"]
            checks["brute_force_agrees_with_enumeration_histogram"] = bool(
                brute["weight_histogram"] == claimed_histogram
            )
            checks["brute_force_finds_only_published_X_rows"] = bool(
                brute["codewords"] == enumerated
            )

            # Independent stage-2: no published X-row admits a light completion.
            fits = 0
            for index in range(hx.shape[0]):
                x_part = hx[index]
                budget = 7 - int(x_part.sum())
                if budget < 0:
                    continue
                if coset_fits_rank_test(perturbation[index], z_rows, x_part, budget):
                    fits += 1
            checks["n_light_completions_found_by_rank_test"] = fits
            checks["rank_test_agrees_no_mixed_light_element"] = bool(
                fits == 0
                and item["theorem_gate"]["mixed_light_element_exists"] is False
            )
            checks["recorded_rank_V7_equals_r_Zsub"] = bool(
                item["per_w"]["7"]["rank_Vw"] == r_zsub < EXPECTED_RANK
            )

        checks["all_checks_passed"] = all(
            value is True
            for key, value in checks.items()
            if isinstance(value, bool) and key != "all_checks_passed"
        )
        all_ok = all_ok and checks["all_checks_passed"]
        report["results"].append(checks)
        print(
            f"{label:>34s} verified={checks['all_checks_passed']}",
            flush=True,
        )

    report["all_codes_crosschecked"] = all_ok
    report["verdict_under_crosscheck"] = processed["verdict"]
    report["verdict"] = "POSITIVE" if all_ok else "NEGATIVE"
    report["verdict_reason"] = (
        "Every independent artifact and bounded structural cross-check passed."
        if all_ok
        else "At least one independent artifact or bounded structural cross-check failed."
    )
    report["wall_s"] = round(time.perf_counter() - started, 6)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"all_codes_crosschecked={all_ok}; wrote {OUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
