#!/usr/bin/env python3
"""Produce the exact all-L saturation certificate for the open 2xL layer DLA.

The finite certificate has three deliberately distinct levels of evidence:

* full generator-monotone Pauli closure at L=2,...,6;
* exact Pauli-string affine-fiber production BFS at L=2,...,7;
* an independent rational Clifford grade census at L=2,...,7.

The proof of uniformity is mathematical (the bounded square move, the
no-rank-loss two-token transport lemma, and induction in L); finite runs are
checks, not substitutes for that proof.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import resource
import sys
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from e157_saturation_local_move import (
    audit_local_schemas,
    checkerboard_mask,
    grid,
    pack,
    parity,
    symplectic,
)
from e158_pair_fiber_bfs import (
    base_fiber_L2,
    extend_seed_fiber,
    pair_digest,
    pair_law,
    spread_pair_fibers,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "algebra_growth" / "alll_saturation.json"
SCRIPT = "experiments/e159_alll_saturation.py"
EXPECTED_DIMENSIONS = {
    2: 56,
    3: 1_056,
    4: 16_256,
    5: 262_656,
    6: 4_192_256,
    7: 67_117_056,
}


def rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_dimension(L: int) -> int:
    n = 2 * L
    half = 1 << (n - 1)
    return half * (half - ((-1) ** L))


def target_law(x: int, z: int, L: int) -> bool:
    n = 2 * L
    J = (1 << n) - 1
    if parity(z) != 0:
        return False
    if z == J:
        return (L & 1) == 1
    d = checkerboard_mask(L)
    return parity((J ^ z) & x) == (1 ^ parity(z & d))


def gf2_rank(vectors: Iterable[int], stop_at: int | None = None) -> int:
    basis: dict[int, int] = {}
    for value in vectors:
        row = value
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                if stop_at is not None and len(basis) == stop_at:
                    return stop_at
                break
            row ^= basis[pivot]
    return len(basis)


def nullspace_basis(rows: list[int], n: int) -> list[int]:
    matrix = [row & ((1 << n) - 1) for row in rows if row]
    pivot_columns: list[int] = []
    rank = 0
    for column in range(n):
        pivot_row = next(
            (index for index in range(rank, len(matrix)) if (matrix[index] >> column) & 1),
            None,
        )
        if pivot_row is None:
            continue
        matrix[rank], matrix[pivot_row] = matrix[pivot_row], matrix[rank]
        for index in range(len(matrix)):
            if index != rank and ((matrix[index] >> column) & 1):
                matrix[index] ^= matrix[rank]
        pivot_columns.append(column)
        rank += 1
        if rank == len(matrix):
            break
    pivot_set = set(pivot_columns)
    output: list[int] = []
    for free in range(n):
        if free in pivot_set:
            continue
        vector = 1 << free
        for index, pivot in enumerate(pivot_columns):
            if (matrix[index] >> free) & 1:
                vector ^= 1 << pivot
        if any(parity(row & vector) for row in rows):
            raise AssertionError("GF(2) nullspace construction failed")
        output.append(vector)
    return output


def pair_representative(z: int, L: int) -> int:
    n = 2 * L
    J = (1 << n) - 1
    lam = J ^ z
    rhs = 1 ^ parity(z & checkerboard_mask(L))
    if lam == 0:
        raise ValueError("the ordinary pair representative excludes z=J")
    return (lam & -lam) if rhs else 0


def adjusted_pair_point(z: int, functional: int, desired: int, L: int) -> int:
    """Choose x in the full pair fiber with functional.x = desired."""
    n = 2 * L
    J = (1 << n) - 1
    x = pair_representative(z, L)
    if parity(functional & x) == desired:
        return x
    if functional == 0 or functional & z:
        raise AssertionError("seed pairing expects a nonempty support disjoint from the next pair")
    outside_not_functional = (J ^ z) & (J ^ functional)
    if outside_not_functional == 0:
        raise AssertionError("the target z=J is the unique forbidden adaptive case")
    in_functional = functional & -functional
    outside = outside_not_functional & -outside_not_functional
    direction = in_functional ^ outside
    if parity((J ^ z) & direction) != 0 or parity(functional & direction) != 1:
        raise AssertionError("adaptive pair direction has the wrong two parities")
    x ^= direction
    if not pair_law(x, z, L) or parity(functional & x) != desired:
        raise AssertionError("adaptive pair point does not meet its two constraints")
    return x


def construct_even_z_seed(z: int, L: int, fibers: dict[int, set[int]]) -> tuple[int, int]:
    """Construct one reached label in every nonzero even z != J fiber."""
    n = 2 * L
    J = (1 << n) - 1
    if z == 0 or z == J or parity(z):
        raise ValueError("construct_even_z_seed requires nonzero even z != J")
    support = [index for index in range(n) if (z >> index) & 1]
    pairs = [
        (1 << support[index]) ^ (1 << support[index + 1])
        for index in range(0, len(support), 2)
    ]
    current_z = pairs[0]
    current_x = pair_representative(current_z, L)
    if current_x not in fibers[current_z]:
        raise AssertionError("the exact pair BFS is missing the initial seed point")
    for next_pair in pairs[1:]:
        desired = 1 ^ parity(current_x & next_pair)
        next_x = adjusted_pair_point(next_pair, current_z, desired, L)
        if next_x not in fibers[next_pair]:
            raise AssertionError("the exact pair BFS is missing an adaptive point")
        left = pack(current_x, current_z, n)
        right = pack(next_x, next_pair, n)
        if symplectic(left, right, n) != 1:
            raise AssertionError("higher-z seed parents do not anticommute")
        current_x ^= next_x
        current_z ^= next_pair
    if current_z != z or not target_law(current_x, z, L):
        raise AssertionError("higher-z telescoping produced the wrong law label")
    return current_x, current_z


def audit_universal_fibers(
    L: int, fibers: dict[int, set[int]], deadline: float
) -> dict:
    """Check every even-z row of the all-L construction, without sampling."""
    n = 2 * L
    J = (1 << n) - 1
    bonds, _ = grid(L)
    bond_masks = [(1 << u) ^ (1 << v) for u, v in bonds]

    odd_digest = hashlib.sha256()
    odd_count = 0
    odd_width = max(1, (n + 7) // 8)
    for a in range(1 << n):
        if parity(a) == 0:
            continue
        if time.process_time() > deadline:
            raise TimeoutError("NON-DECISIVE: e159 universal-fiber process-time budget expired")
        odd_count += 1
        odd_digest.update(a.to_bytes(odd_width, "little"))
        if a & (a - 1) == 0:
            continue  # X_i is a generator.
        boundary = next((edge for edge in bond_masks if parity(a & edge) == 1), None)
        if boundary is None:
            raise AssertionError("an odd proper support has no boundary edge in a connected ladder")
        anchor = min(fibers[boundary])
        mate = anchor ^ a
        if mate not in fibers[boundary]:
            raise AssertionError("a full boundary-edge fiber failed to contain the odd difference")
        if symplectic(pack(anchor, boundary, n), pack(mate, boundary, n), n) != 1:
            raise AssertionError("the pair-difference parents must anticommute")

    if odd_count != 1 << (n - 1):
        raise AssertionError("the pure-X fiber does not contain every odd vector")

    seed_digest = hashlib.sha256()
    standard_rows = 0
    direction_rank_min = n
    direction_rank_max = 0
    special_J_rank = 0
    for z in range(1 << n):
        if parity(z):
            continue
        if time.process_time() > deadline:
            raise TimeoutError("NON-DECISIVE: e159 universal-z process-time budget expired")
        if z == J:
            if L & 1:
                edge = bond_masks[0]
                predecessor = J ^ edge
                pred_x, pred_z = construct_even_z_seed(predecessor, L, fibers)
                if pred_z != predecessor or parity(pred_x & edge) != 1:
                    raise AssertionError("odd-L z=J predecessor is not forced to anticommute")
                if symplectic(pack(pred_x, predecessor, n), pack(0, edge, n), n) != 1:
                    raise AssertionError("odd-L z=J bond production is not sound")
                special_J_rank = gf2_rank([1 << index for index in range(n)])
                if special_J_rank != n:
                    raise AssertionError("singleton X toggles do not span the z=J fiber")
            else:
                # For even L, the same predecessor law fixes edge.x=0, so the
                # attempted final bond commutator vanishes exactly.
                edge = bond_masks[0]
                predecessor = J ^ edge
                lam = J ^ predecessor
                rhs = 1 ^ parity(predecessor & checkerboard_mask(L))
                if lam != edge or rhs != 0:
                    raise AssertionError("even-L z=J exclusion parity is wrong")
            continue

        standard_rows += 1
        if z == 0:
            rank = n - 1  # the already checked affine odd-x hyperplane
            seed_x = 1
        else:
            seed_x, seed_z = construct_even_z_seed(z, L, fibers)
            if seed_z != z:
                raise AssertionError("seed support mismatch")
            q0 = z & -z
            kernel = nullspace_basis([J, z], n)
            toggles = [q0] + [q0 ^ value for value in kernel]
            if not all(parity(q) == 1 and parity(q & z) == 1 for q in toggles):
                raise AssertionError("a purported pure-X toggle is outside its affine class")
            rank = gf2_rank(toggles)
            if rank != n - 1:
                raise AssertionError("pure-X toggles do not span the ordinary fiber tangent")
        direction_rank_min = min(direction_rank_min, rank)
        direction_rank_max = max(direction_rank_max, rank)
        seed_digest.update(pack(seed_x, z, n).to_bytes(max(1, (2 * n + 7) // 8), "little"))

    expected_standard_rows = (1 << (n - 1)) - 1
    if standard_rows != expected_standard_rows:
        raise AssertionError(
            f"universal audit covered {standard_rows}/{expected_standard_rows} ordinary even-z rows"
        )
    special_size = (1 << n) if (L & 1) else 0
    dimension = standard_rows * (1 << (n - 1)) + special_size
    if dimension != exact_dimension(L):
        raise AssertionError("fiber count does not equal the alternation dimension formula")
    return {
        "odd_w0_strings_checked": odd_count,
        "odd_w0_sha256": odd_digest.hexdigest(),
        "ordinary_even_z_rows_checked": standard_rows,
        "ordinary_fiber_tangent_rank_min": direction_rank_min,
        "ordinary_fiber_tangent_rank_max": direction_rank_max,
        "special_z_J_included": bool(L & 1),
        "special_z_J_tangent_rank": special_J_rank,
        "seed_labels_sha256": seed_digest.hexdigest(),
        "derived_dimension": dimension,
    }


def direct_generators(L: int) -> tuple[list[int], list[int], list[int]]:
    n = 2 * L
    x_generators = [1 << index for index in range(n)]
    bond_masks = [(1 << u) ^ (1 << v) for u, v in grid(L)[0]]
    labels = x_generators + [bond << n for bond in bond_masks]
    return labels, x_generators, bond_masks


def direct_closure_worker(L: int, cpu_budget: float, sender) -> None:
    """Isolated full Pauli closure, so every L gets an honest peak-RSS reading."""
    try:
        n = 2 * L
        mask = (1 << n) - 1
        labels, x_generators, bond_masks = direct_generators(L)
        reached = set(labels)
        queue: deque[int] = deque(labels)
        started = time.process_time()
        deadline = started + cpu_budget
        processed = 0
        while queue:
            label = queue.popleft()
            x, z = label & mask, label >> n
            for x_generator in x_generators:
                if z & x_generator:
                    output = label ^ x_generator
                    if output not in reached:
                        reached.add(output)
                        queue.append(output)
            for bond in bond_masks:
                if parity(x & bond):
                    output = label ^ (bond << n)
                    if output not in reached:
                        reached.add(output)
                        queue.append(output)
            processed += 1
            if (processed & 0x3FFFF) == 0 and time.process_time() > deadline:
                raise TimeoutError(
                    f"NON-DECISIVE: full L={L} closure exceeded {cpu_budget} process seconds"
                )

        violations = 0
        label_xor = 0
        label_sum_mod_2_64 = 0
        for label in reached:
            x, z = label & mask, label >> n
            if not target_law(x, z, L):
                violations += 1
            label_xor ^= label
            label_sum_mod_2_64 = (label_sum_mod_2_64 + label) & ((1 << 64) - 1)
        sender.send(
            {
                "L": L,
                "method": "full_generator_monotone_Pauli_closure",
                "closure_strings": len(reached),
                "set_law_violations": violations,
                "predicted_dimension": exact_dimension(L),
                "dimension_match": len(reached) == exact_dimension(L),
                "processed_queue_rows": processed,
                "label_xor": str(label_xor),
                "label_sum_mod_2_64": str(label_sum_mod_2_64),
                "process_time_seconds": time.process_time() - started,
                "peak_rss_bytes": rss_bytes(),
            }
        )
    except BaseException as exc:  # returned to the parent; producer still fails loudly
        sender.send({"L": L, "error": f"{type(exc).__name__}: {exc}"})
    finally:
        sender.close()


def isolated_direct_closure(L: int, cpu_budget: float) -> dict:
    context = mp.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(target=direct_closure_worker, args=(L, cpu_budget, sender))
    process.start()
    sender.close()
    row = receiver.recv()
    receiver.close()
    process.join()
    if process.exitcode != 0:
        raise RuntimeError(f"isolated L={L} closure exited with code {process.exitcode}")
    if "error" in row:
        raise RuntimeError(row["error"])
    return row


def clifford_grade_audit(L: int) -> dict:
    """Independent exact grade census for the rational JW/Clifford route."""
    N = 4 * L
    rung_grades = [4 * L - 4 * column + 2 for column in range(1, L)]
    available = sorted([2] + rung_grades)
    expected = list(range(2, N, 4))
    if available != expected:
        raise AssertionError("Hamiltonian-path terms do not seed every 2 mod 4 grade")
    terms = [{"grade": grade, "binomial": math.comb(N, grade)} for grade in available]
    grade_sum = sum(row["binomial"] for row in terms)
    if grade_sum != exact_dimension(L):
        raise AssertionError("Clifford grade sum disagrees with the containment dimension")
    middle = 2 * L
    return {
        "majoranas": N,
        "seeded_grades": available,
        "grade_terms": terms,
        "grade_sum": grade_sum,
        "root_of_unity_closed_form": exact_dimension(L),
        "middle_grade_present": middle in available,
        "middle_grade": middle,
        "middle_orbit_resolution": (
            "direct rational one-index-swap orbit; equivalently, after algebraic closure a "
            "basis monomial has nonzero projections to both Hodge-star summands"
            if middle in available
            else "middle grade is 0 mod 4 and is outside the target Lie class"
        ),
    }


def record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def main() -> bool:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-L", type=int, default=7)
    parser.add_argument("--full-closure-max-L", type=int, default=6)
    parser.add_argument("--cpu-budget", type=float, default=240.0)
    parser.add_argument("--closure-cpu-budget", type=float, default=180.0)
    args = parser.parse_args()
    if args.max_L < 2 or args.full_closure_max_L < 2:
        raise ValueError("L bounds must be at least 2")
    if args.full_closure_max_L > args.max_L:
        raise ValueError("full closure bound cannot exceed max-L")
    if args.cpu_budget <= 0 or args.closure_cpu_budget <= 0:
        raise ValueError("process-time budgets must be positive")

    generated_utc = datetime.now(timezone.utc).isoformat()
    main_started = time.process_time()
    deadline = main_started + args.cpu_budget
    checks: list[dict] = []
    points: list[dict] = []

    seed_z, seed_points, _ = base_fiber_L2()
    for L in range(2, args.max_L + 1):
        point_started = time.process_time()
        if L > 2:
            seed_z, seed_points, _ = extend_seed_fiber(L - 1, seed_z, seed_points)
        local = audit_local_schemas(L, deadline)
        fibers, transitions = spread_pair_fibers(L, seed_z, seed_points, deadline)
        n = 2 * L
        target_rank = n - 1
        ranks = []
        for values in fibers.values():
            anchor = min(values)
            ranks.append(gf2_rank((value ^ anchor for value in values), stop_at=target_rank))
        rank_set = sorted(set(ranks))
        universal = audit_universal_fibers(L, fibers, deadline)
        clifford = clifford_grade_audit(L)
        point_cpu = time.process_time() - point_started

        record(
            checks,
            f"L={L} bounded square schemas",
            local["triple_witnesses"] == 4 * (L - 1),
            f"{local['triple_witnesses']} exact three-corner witnesses",
        )
        record(
            checks,
            f"L={L} pair-fiber affine Pauli BFS",
            len(fibers) == math.comb(n, 2) and rank_set == [target_rank],
            f"{len(fibers)}/{math.comb(n, 2)} supports, ranks={rank_set}",
        )
        record(
            checks,
            f"L={L} universal even-z audit",
            universal["ordinary_even_z_rows_checked"] == (1 << (n - 1)) - 1,
            f"{universal['ordinary_even_z_rows_checked']} ordinary rows checked without sampling",
        )
        record(
            checks,
            f"L={L} exact dimension",
            universal["derived_dimension"] == EXPECTED_DIMENSIONS[L] == exact_dimension(L),
            f"dimension={universal['derived_dimension']}",
        )
        record(
            checks,
            f"L={L} independent Clifford census",
            clifford["grade_sum"] == universal["derived_dimension"],
            f"grades={clifford['seeded_grades']}, sum={clifford['grade_sum']}",
        )

        points.append(
            {
                "L": L,
                "n": n,
                "predicted_and_proved_dimension": exact_dimension(L),
                "local_schema": local,
                "affine_pair_bfs": {
                    "method": "exact Pauli labels; full x-set retained in every weight-two fiber",
                    "pair_supports": len(fibers),
                    "pair_fiber_size": 1 << (n - 1),
                    "pair_strings_retained": len(fibers) * (1 << (n - 1)),
                    "tangent_ranks": rank_set,
                    "transport_tree_edges": len(transitions),
                    "pair_labels_sha256": pair_digest(fibers, n),
                },
                "universal_fiber_audit": universal,
                "clifford_cross_check": clifford,
                "construction_process_time_seconds": point_cpu,
                "construction_cumulative_peak_rss_bytes": rss_bytes(),
                "full_closure": None,
            }
        )

        if L < args.max_L:
            top_left = 2 * (L - 2)
            top_right = 2 * (L - 1)
            seed_z = (1 << top_left) ^ (1 << top_right)
            seed_points = set(fibers[seed_z])
        del fibers
        gc.collect()

    full_closures: dict[int, dict] = {}
    for L in range(2, args.full_closure_max_L + 1):
        row = isolated_direct_closure(L, args.closure_cpu_budget)
        full_closures[L] = row
        record(
            checks,
            f"L={L} full exact generator closure",
            row["dimension_match"] and row["set_law_violations"] == 0,
            (
                f"reached={row['closure_strings']}, violations={row['set_law_violations']}, "
                f"cpu={row['process_time_seconds']:.6f}s, peak_rss={row['peak_rss_bytes']}"
            ),
        )

    for point in points:
        L = point["L"]
        if L in full_closures:
            point["full_closure"] = full_closures[L]
            point["evidence_strength"] = (
                "full exact generator-monotone closure + exhaustive set-law membership + "
                "affine-fiber rank certificate + Clifford census"
            )
        else:
            point["evidence_strength"] = (
                "NOT a full closure enumeration: exact affine-fiber production/rank certificate + "
                "exhaustive even-z construction audit + Clifford census"
            )

    all_passed = all(row["passed"] for row in checks)
    input_paths = [
        ROOT / "proofs" / "alternation_law.md",
        ROOT / "proofs" / "alternation_l5_char0.md",
        ROOT / "results" / "algebra" / "alternation_l5_char0.json",
        ROOT / "experiments" / "e157_saturation_local_move.py",
        ROOT / "experiments" / "e158_pair_fiber_bfs.py",
        ROOT / SCRIPT,
    ]
    payload = {
        "schema": "alll_saturation/v1",
        "meta": {
            "generated_utc": generated_utc,
            "producer": SCRIPT,
            "command": ".venv/bin/python experiments/e159_alll_saturation.py",
            "cwd": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "pythonhashseed": os.environ.get("PYTHONHASHSEED"),
            "arithmetic": "exact Python integers over GF(2) and exact binomial integers; no floating-point mathematics",
            "main_process_time_seconds": time.process_time() - main_started,
            "full_closure_child_process_time_seconds": sum(
                row["process_time_seconds"] for row in full_closures.values()
            ),
            "main_cumulative_peak_rss_bytes": rss_bytes(),
            "main_process_time_budget_seconds": args.cpu_budget,
            "per_closure_process_time_budget_seconds": args.closure_cpu_budget,
            "input_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in input_paths
            },
            "provenance": (
                "Primary proof: bounded Pauli production induction closing O1/O2.  Independent "
                "Clifford Hamiltonian-path construction (Facts 1-2) credited to the "
                "AllLQuadraticNoGo wave-18 sibling front."
            ),
        },
        "data": {
            "classification": (
                "[THEOREM] for every L>=2 over Q: g_2xL = so_m (+) so_m for even L and "
                "sp_m (+) sp_m for odd L, m=2^(2L-1)"
            ),
            "dimension_formula": "2^(n-1) * (2^(n-1) - (-1)^L), n=2L",
            "production_system": {
                "labels": "(a|b) in GF(2)^(2n)",
                "generators": "(e_v|0) and (0|e_u+e_v) for every open-ladder bond uv",
                "rule": (
                    "from p,q already reached with <p,q>=a.d+b.c=1, adjoin p+q; "
                    "[M_p,M_q]=+-2 M_(p+q) over Q"
                ),
                "target_predicate": (
                    "b even and ((b!=J and (J+b).a=1+b.d) or "
                    "(b=J and L odd)); d is either checkerboard class"
                ),
                "uniform_invariant": (
                    "at width L every weight-two b fiber is its full affine law coset, "
                    "with tangent ker((J+b).) of rank 2L-1"
                ),
                "second_independent_move": (
                    "the two diagonal length-two routes around one square produce every "
                    "pure-X three-corner string; two such strings add the two new tangent "
                    "directions in the L->L+1 step"
                ),
            },
            "points": points,
            "status": {
                "saturation": "[THEOREM] every L>=2",
                "equality": "[THEOREM] every L>=2 over Q",
                "type": "[THEOREM] D/orthogonal for even L; C/symplectic for odd L",
                "form": (
                    "[THEOREM] explicit nondegenerate signed-checkerboard form on each "
                    "prod-X half; symmetric for even L, alternating for odd L; unique up to "
                    "scale on each half"
                ),
                "scope_limit": (
                    "open 2xL ladders and characteristic zero/Q only; no periodic-boundary, "
                    "higher-row, or full 3D Ising solution claim"
                ),
            },
            "clifford_route_usage": {
                "used_by_primary_induction": False,
                "role": "independent structural and dimension cross-check",
                "middle_grade": (
                    "covered directly over Q by one-index swaps using generated grade-2 "
                    "bilinears; the Hodge-star argument is only an algebraic-closure cross-check"
                ),
            },
            "gaussian_corollary": (
                "[THEOREM] every injective Lie-bracket-preserving realization by quadratic "
                "operators on m fermionic modes obeys m(2m-1) >= dim(g_L); hence m grows "
                "exponentially in L. In particular, unitary basis changes on the physical "
                "Hilbert space have m=n=2L and are impossible for every L>=2. "
                "Non-invariant code compressions are outside this statement."
            ),
            "checks": checks,
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not all_passed:
        raise AssertionError("one or more all-L saturation checks failed")
    print(f"wrote {OUTPUT.relative_to(ROOT)}", flush=True)
    print("PASS all-L saturation and classification", flush=True)
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
