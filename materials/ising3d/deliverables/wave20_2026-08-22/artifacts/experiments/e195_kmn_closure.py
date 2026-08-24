#!/usr/bin/env python3
"""[COMPUTATION] Exact-integer Lie words with finite-field rank lower bounds on spin blocks."""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from e194_kmn_reduction import (
    GOOD_PRIMES,
    budget_tick,
    hamiltonian_complete_bipartite,
    local_term_branching_dimension,
    max_rss_bytes,
    representation_blocks,
    schur_commutant_dimension,
    sector_matrices,
    spin_labels,
    structural_upper,
)


@dataclass(frozen=True)
class ModularWord:
    blocks: tuple[np.ndarray, ...]
    parent: int
    action: int
    depth: int


def generator_blocks(
    left_size: int, right_size: int
) -> tuple[tuple[tuple[int, int], ...], tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    labels = representation_blocks(left_size, right_size)
    sectors = tuple(sector_matrices(*label) for label in labels)
    return labels, tuple(sector.a for sector in sectors), tuple(sector.b for sector in sectors)


def flatten_blocks(blocks: Sequence[np.ndarray]) -> np.ndarray:
    if not blocks:
        return np.zeros(0, dtype=np.int64)
    return np.concatenate([np.asarray(block, dtype=np.int64).reshape(-1) for block in blocks])


def generator_digest(
    labels: Sequence[tuple[int, int]],
    field_blocks: Sequence[np.ndarray],
    bond_blocks: Sequence[np.ndarray],
) -> str:
    digest = hashlib.sha256()
    for label, field, bond in zip(labels, field_blocks, bond_blocks, strict=True):
        digest.update(f"{label[0]},{label[1]};".encode("ascii"))
        digest.update(np.asarray(field, dtype="<i8").tobytes())
        digest.update(np.asarray(bond, dtype="<i8").tobytes())
    return digest.hexdigest()


def commutator_mod(generator: np.ndarray, value: np.ndarray, prime: int) -> np.ndarray:
    return (generator @ value - value @ generator) % prime


def modular_closure(left_size: int, right_size: int, prime: int) -> dict[str, object]:
    """Close under ad_A,ad_B; every accepted row is an actual integral Lie word mod p."""
    started = time.process_time()
    labels, raw_field, raw_bond = generator_blocks(left_size, right_size)
    field = tuple(np.asarray(block, dtype=np.int64) % prime for block in raw_field)
    bond = tuple(np.asarray(block, dtype=np.int64) % prime for block in raw_bond)
    generators = (field, bond)
    coordinate_dimension = sum(block.size for block in field)
    pivots: dict[int, np.ndarray] = {}
    words: list[ModularWord] = []
    frontier: list[int] = []
    depth_profile: dict[int, int] = {}

    def reduce_and_insert(
        raw_vector: np.ndarray,
        blocks: tuple[np.ndarray, ...],
        *,
        parent: int,
        action: int,
        depth: int,
    ) -> bool:
        vector = np.asarray(raw_vector, dtype=np.int64) % prime
        while True:
            nonzero = np.flatnonzero(vector)
            if not nonzero.size:
                return False
            lead = int(nonzero[0])
            existing = pivots.get(lead)
            if existing is None:
                inverse = pow(int(vector[lead]), prime - 2, prime)
                normalized = (vector * inverse) % prime
                pivots[lead] = normalized
                words.append(
                    ModularWord(
                        tuple(np.asarray(block, dtype=np.int64) % prime for block in blocks),
                        parent,
                        action,
                        depth,
                    )
                )
                depth_profile[depth] = depth_profile.get(depth, 0) + 1
                return True
            factor = int(vector[lead])
            vector = (vector - factor * existing) % prime

    for action, blocks in enumerate(generators):
        inserted = reduce_and_insert(
            flatten_blocks(blocks),
            blocks,
            parent=-1,
            action=action,
            depth=1,
        )
        if not inserted:
            raise AssertionError("the two collective seeds must be independent")
        frontier.append(len(words) - 1)

    while frontier:
        next_frontier: list[int] = []
        for parent in frontier:
            word = words[parent]
            for action, generators_for_action in enumerate(generators):
                candidate = tuple(
                    commutator_mod(generator, value, prime)
                    for generator, value in zip(
                        generators_for_action, word.blocks, strict=True
                    )
                )
                if reduce_and_insert(
                    flatten_blocks(candidate),
                    candidate,
                    parent=parent,
                    action=action,
                    depth=word.depth + 1,
                ):
                    next_frontier.append(len(words) - 1)
                    if len(words) % 64 == 0:
                        budget_tick(
                            started,
                            f"K_{{{left_size},{right_size}}} mod {prime} rank {len(words)}",
                        )
        frontier = next_frontier

    recipe_text = ";".join(
        f"{word.parent},{word.action},{word.depth},{lead}"
        for word, lead in zip(words, pivots, strict=True)
    )
    budget_tick(started, f"K_{{{left_size},{right_size}}} mod {prime} saturated")
    return {
        "tag": "[COMPUTATION]",
        "prime": prime,
        "rank_Fp": len(words),
        "rank_lower_bound_over_Q": len(words),
        "saturated_over_Fp": not frontier,
        "maximum_word_depth": max(word.depth for word in words),
        "depth_profile": {str(depth): count for depth, count in sorted(depth_profile.items())},
        "pivot_count": len(pivots),
        "coordinate_dimension": coordinate_dimension,
        "block_labels": [list(label) for label in labels],
        "block_matrix_dimensions": [(left + 1) * (right + 1) for left, right in labels],
        "generator_sha256": generator_digest(labels, raw_field, raw_bond),
        "recipe_pivot_sha256": hashlib.sha256(recipe_text.encode("ascii")).hexdigest(),
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
        "direction": (
            "rank_Fp <= dimension_Q for the integral Lie-word matrix; modular rank is "
            "a rational lower bound only"
        ),
    }


def case_closure(
    left_size: int,
    right_size: int,
    *,
    primes: Sequence[int] = GOOD_PRIMES,
) -> dict[str, object]:
    started = time.process_time()
    upper = structural_upper(left_size, right_size)
    modular = [modular_closure(left_size, right_size, prime) for prime in primes]
    ranks = [int(row["rank_Fp"]) for row in modular]
    upper_dimension = int(upper["upper_dimension"])
    exact_match = bool(ranks and all(rank == upper_dimension for rank in ranks))
    vertices = left_size + right_size
    local_dimension = local_term_branching_dimension(left_size, right_size)
    labels, field, _ = generator_blocks(left_size, right_size)
    reduced_coordinates = sum(block.size for block in field)
    full_coordinates = schur_commutant_dimension(left_size, right_size)
    if left_size != right_size and reduced_coordinates != full_coordinates:
        raise AssertionError((reduced_coordinates, full_coordinates))
    budget_tick(started, f"K_{{{left_size},{right_size}}} complete case")
    return {
        "tag": "[COMPUTATION]",
        "graph": f"K_{{{left_size},{right_size}}}",
        "part_sizes": [left_size, right_size],
        "n_vertices": vertices,
        "hamiltonian": hamiltonian_complete_bipartite(left_size, right_size),
        "maximum_degree": max(left_size, right_size),
        "spin_labels": [list(spin_labels(left_size)), list(spin_labels(right_size))],
        "representation_blocks": [list(label) for label in labels],
        "full_Sm_times_Sn_commutant_coordinates": full_coordinates,
        "reduced_faithful_coordinates": reduced_coordinates,
        "balanced_transpose_blocks_removed": left_size == right_size,
        "structural_upper": upper,
        "modular_closures": modular,
        "modular_ranks": ranks,
        "lower_equals_upper": exact_match,
        "dimension_Q": upper_dimension if exact_match else None,
        "quadratic_ceiling": vertices * (2 * vertices - 1),
        "clears_quadratic_ceiling": bool(
            exact_match and upper_dimension > vertices * (2 * vertices - 1)
        ),
        "local_term_dimension": local_dimension,
        "local_term_is_strictly_larger": bool(
            exact_match and local_dimension > upper_dimension
        ),
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def run_anchor_closures() -> tuple[dict[str, object], ...]:
    return (
        case_closure(4, 4),
        case_closure(4, 5),
    )


def main() -> int:
    rows = run_anchor_closures()
    if not all(row["lower_equals_upper"] for row in rows):
        raise AssertionError(rows)
    print(
        "PASS e195: "
        + ", ".join(f"{row['graph']}={row['dimension_Q']}" for row in rows)
        + f", rss={max_rss_bytes()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
