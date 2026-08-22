#!/usr/bin/env python3
"""Partial exact modular study of the open 4x4 Ising-layer Lie algebra.

For ``A=sum_v X_v`` and ``B=sum_(u,v) Z_u Z_v`` on the 4x4 open grid, this
producer builds literal nested Pauli-word lower minors over a 31-bit prime.
The full D4-orbit support is deliberately not preallocated: it already exceeds
an explicit support cap, whereas the sparse ordered-monomial echelon closure
can certify depth-by-depth lower bounds.  A JSON checkpoint stores the raw-word
tree and last completed frontier; a later invocation deterministically rebuilds
that completed basis and continues from it.

No characteristic-zero equality or factor type is inferred from modular rank.
The exact rational data in this script are symmetry-projector ranks and the
resulting broad D4 x global-spin-flip commutant container only.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import resource
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e97_algebra_4x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_4x4.json"
E20_PATH = ROOT / "experiments" / "e20_algebra_growth.py"
ROWS = COLS = 4
N = ROWS * COLS
GOOD_PRIME = 2_147_483_647
DEFAULT_WALL_SECONDS = 3_600.0
DEFAULT_MEMORY_LIMIT_BYTES = 48 * 1024**3
DEFAULT_SUPPORT_PROBE_CAP = 250_000
CHECKPOINT_INTERVAL_SECONDS = 60.0


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def physical_memory_bytes() -> int | None:
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (AttributeError, OSError, ValueError):
        return None


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(N))


def group_from_generators(generators: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    identity = tuple(range(N))
    group = [identity]
    seen = {identity}
    head = 0
    while head < len(group):
        current = group[head]
        head += 1
        for generator in generators:
            product = compose(generator, current)
            if product not in seen:
                seen.add(product)
                group.append(product)
    return tuple(group)


def geometric_generators() -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    def site(row: int, column: int) -> int:
        return row * COLS + column

    row_reflection = tuple(
        site(ROWS - 1 - row, column)
        for row in range(ROWS)
        for column in range(COLS)
    )
    column_reflection = tuple(
        site(row, COLS - 1 - column)
        for row in range(ROWS)
        for column in range(COLS)
    )
    transpose = tuple(
        site(column, row)
        for row in range(ROWS)
        for column in range(COLS)
    )
    return row_reflection, column_reflection, transpose


def cycles(permutation: tuple[int, ...]) -> tuple[int, ...]:
    unseen = set(range(N))
    lengths: list[int] = []
    while unseen:
        start = unseen.pop()
        current = permutation[start]
        length = 1
        while current != start:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def fixed_spin_count(permutation: tuple[int, ...], complement: int) -> int:
    lengths = cycles(permutation)
    if complement and any(length & 1 for length in lengths):
        return 0
    return 1 << len(lengths)


def c2_character_sector_dimensions() -> dict[str, int]:
    """Exact character-projector ranks for row x column x global-spin flip."""
    row_reflection, column_reflection, _ = geometric_generators()
    elements: list[tuple[int, int, int, tuple[int, ...]]] = []
    for row_bit in (0, 1):
        for column_bit in (0, 1):
            permutation = tuple(range(N))
            if row_bit:
                permutation = compose(row_reflection, permutation)
            if column_bit:
                permutation = compose(column_reflection, permutation)
            for flip_bit in (0, 1):
                elements.append((row_bit, column_bit, flip_bit, permutation))
    dimensions: dict[str, int] = {}
    for row_character in (0, 1):
        for column_character in (0, 1):
            for flip_character in (0, 1):
                numerator = sum(
                    (-1) ** (
                        row_character * row_bit
                        + column_character * column_bit
                        + flip_character * flip_bit
                    )
                    * fixed_spin_count(permutation, flip_bit)
                    for row_bit, column_bit, flip_bit, permutation in elements
                )
                if numerator % 8:
                    raise AssertionError("C2^3 character trace did not divide by eight")
                dimensions[f"{row_character}{column_character}{flip_character}"] = numerator // 8
    if sum(dimensions.values()) != 1 << N:
        raise AssertionError("C2^3 dimensions did not exhaust the physical module")
    return dimensions


def d4_class_map() -> dict[tuple[int, ...], str]:
    row_reflection, _column_reflection, transpose = geometric_generators()
    rotation = compose(transpose, row_reflection)
    identity = tuple(range(N))
    rotation_squared = compose(rotation, rotation)
    rotation_cubed = compose(rotation, rotation_squared)
    group = group_from_generators((row_reflection, transpose))
    classes: dict[tuple[int, ...], str] = {}
    for permutation in group:
        if permutation == identity:
            classes[permutation] = "identity"
        elif permutation == rotation_squared:
            classes[permutation] = "rotation_180"
        elif permutation in (rotation, rotation_cubed):
            classes[permutation] = "rotation_quarter"
        else:
            fixed_sites = sum(index == target for index, target in enumerate(permutation))
            classes[permutation] = "reflection_diagonal" if fixed_sites else "reflection_axis"
    if len(group) != 8 or set(classes.values()) != {
        "identity",
        "rotation_180",
        "rotation_quarter",
        "reflection_axis",
        "reflection_diagonal",
    }:
        raise AssertionError("D4 class construction failed")
    return classes


def d4_global_flip_multiplicity_spaces() -> dict[str, dict[str, int]]:
    """Exact D4 multiplicity-space dimensions in each global-flip eigenspace."""
    row_reflection, _column_reflection, transpose = geometric_generators()
    group = group_from_generators((row_reflection, transpose))
    classes = d4_class_map()
    character_table = {
        "A1": {"identity": 1, "rotation_180": 1, "rotation_quarter": 1, "reflection_axis": 1, "reflection_diagonal": 1},
        "A2": {"identity": 1, "rotation_180": 1, "rotation_quarter": 1, "reflection_axis": -1, "reflection_diagonal": -1},
        "B1": {"identity": 1, "rotation_180": 1, "rotation_quarter": -1, "reflection_axis": 1, "reflection_diagonal": -1},
        "B2": {"identity": 1, "rotation_180": 1, "rotation_quarter": -1, "reflection_axis": -1, "reflection_diagonal": 1},
        "E": {"identity": 2, "rotation_180": -2, "rotation_quarter": 0, "reflection_axis": 0, "reflection_diagonal": 0},
    }
    result: dict[str, dict[str, int]] = {}
    for parity_name, parity_sign in (("even", 1), ("odd", -1)):
        traces = {
            permutation: (fixed_spin_count(permutation, 0) + parity_sign * fixed_spin_count(permutation, 1)) // 2
            for permutation in group
        }
        multiplicities: dict[str, int] = {}
        for irrep, character in character_table.items():
            numerator = sum(character[classes[permutation]] * traces[permutation] for permutation in group)
            if numerator % 8:
                raise AssertionError("D4 multiplicity trace did not divide by eight")
            multiplicities[irrep] = numerator // 8
        result[parity_name] = multiplicities
    weighted_total = sum(
        result[parity][irrep] * (2 if irrep == "E" else 1)
        for parity in ("even", "odd")
        for irrep in ("A1", "A2", "B1", "B2", "E")
    )
    if weighted_total != 1 << N:
        raise AssertionError("D4 multiplicities did not exhaust the physical module")
    return result


def symmetry_decomposition() -> dict[str, Any]:
    sectors = c2_character_sector_dimensions()
    multiplicities = d4_global_flip_multiplicity_spaces()
    candidate_dimensions = [
        multiplicities[parity][irrep]
        for parity in ("even", "odd")
        for irrep in ("A1", "A2", "B1", "B2", "E")
    ]
    broad_dimension = sum(value * value for value in candidate_dimensions)
    return {
        "claim_tag": "[COMPUTATION]",
        "method": "exact finite-group character-projector traces over Q; no floating point arithmetic",
        "c2_character_sector_dimensions": sectors,
        "d4_global_flip_multiplicity_spaces": multiplicities,
        "full_special_linear_candidate_dimensions": candidate_dimensions,
        "unlinked_endomorphism_container_dimension_Q": broad_dimension,
        "one_scalar_linked_sl_candidate_dimension_Q": 1 + sum(value * value - 1 for value in candidate_dimensions),
        "interpretation": (
            "A and B commute with D4 and global spin flip, so their Q-Lie algebra acts only on these "
            "multiplicity spaces. The displayed endomorphism container is a rigorous same-basis symmetry "
            "upper container, but it is far too broad to squeeze the modular lower bound. Treating all ten "
            "multiplicity actions as linked full sl_k blocks is only a [CONJECTURE], not a factor classification."
        ),
        "comparison_to_3x3": (
            "On 3x3 the analogous symmetry container is broader than the proven 8034-dimensional linked "
            "dynamic-commutant container; hence symmetry-isotypic ranks alone do not determine 4x4 factors."
        ),
    }


def tree_checksum(recipes: list[list[int]]) -> int:
    return sum(
        (index + 1) * (parent + 2) * (action + 2) * (depth + 1)
        for index, (parent, action, depth) in enumerate(recipes)
    ) % 1_000_000_007


class SparseLiteralClosure:
    """Sparse ordered-monomial F_p echelon closure of literal nested raw words."""

    def __init__(self, e20: Any, model: Any, memory_limit_bytes: int) -> None:
        self.e20 = e20
        self.model = model
        self.worker = e20.SparseOrbitEngine(model, GOOD_PRIME, memory_limit_bytes, None)
        self.words: list[dict[int, int]] = []
        self.recipes: list[list[int]] = []
        self.word_grades: list[int] = []

    @property
    def rank(self) -> int:
        return len(self.words)

    @property
    def block_counts(self) -> list[int]:
        return [
            sum(grade == block for grade in self.word_grades)
            for block in range(4)
        ]

    def seed(self, action_index: int) -> dict[int, int]:
        terms = self.model.x_terms if action_index == 0 else self.model.zz_terms
        return self.worker.dense_seed(terms)

    def apply_raw(self, action_index: int, word: dict[int, int]) -> dict[int, int]:
        result: dict[int, int] = {}
        for source, value in word.items():
            for target, coefficient in self.worker._transitions(source, action_index).items():
                updated = (result.get(target, 0) + value * coefficient) % GOOD_PRIME
                if updated:
                    result[target] = updated
                else:
                    result.pop(target, None)
        return result

    def accept(self, raw_word: dict[int, int], recipe: list[int]) -> int | None:
        added = self.worker.reduce_add(raw_word)
        if added is None:
            return None
        if added != len(self.words):
            raise AssertionError("sparse-engine row index does not match raw word tree")
        if not raw_word:
            raise AssertionError("accepted raw word was empty")
        grades = {self.e20._pauli_grade(value, self.model.n) for value in raw_word}
        if len(grades) != 1:
            raise AssertionError("a raw Pauli word crossed a grading block")
        self.words.append(raw_word)
        self.recipes.append(recipe)
        self.word_grades.append(grades.pop())
        return added

    def restore_prefix(self, recipes: list[list[int]]) -> None:
        for index, recipe in enumerate(recipes):
            parent, action_index, depth = recipe
            if action_index not in (0, 1) or depth < 1:
                raise AssertionError("invalid saved recipe")
            if parent == -1:
                if depth != 1:
                    raise AssertionError("saved seed did not have depth one")
                raw_word = self.seed(action_index)
            else:
                if not 0 <= parent < index or depth != self.recipes[parent][2] + 1:
                    raise AssertionError("saved tree parent/depth invariant failed")
                raw_word = self.apply_raw(action_index, self.words[parent])
            if self.accept(raw_word, recipe.copy()) is None:
                raise AssertionError("saved raw word was no longer independent")


def artifact(
    *,
    symmetry: dict[str, Any],
    closure: SparseLiteralClosure,
    dimensions: list[int],
    completed_depth: int,
    completed_basis_count: int,
    completed_block_counts: list[int],
    frontier: list[int],
    partial_frontier: list[int],
    partial_depth: int | None,
    partial_candidates_processed: int,
    stop_reason: str,
    started: float,
    rss_started: int,
    memory_limit_bytes: int,
    resumed_from_checkpoint: bool,
    support_probe: dict[str, Any],
) -> dict[str, Any]:
    modular = {
        "claim_tag": "[COMPUTATION]",
        "prime": GOOD_PRIME,
        "arithmetic": (
            "exact F_p ordered-monomial echelon arithmetic; every accepted row is a literal nested Pauli "
            "word formed using [A,w]/2 or [B,w]/2, with 2 invertible modulo the stated odd prime"
        ),
        "rational_implication": (
            "each accepted echelon pivot is a nonzero minor of an integer literal-word matrix modulo p; "
            "therefore its rank is a rigorous lower bound for the corresponding Q-Lie algebra dimension"
        ),
        "closure_status": "saturated" if stop_reason == "saturated" else "partial",
        "stop_reason": stop_reason,
        "completed_depth": completed_depth,
        "dimensions_through_completed_depth": dimensions,
        "grading_block_ranks_completed_depth": completed_block_counts,
        "lower_bound_from_literal_word_minor": closure.rank,
        "grading_block_ranks_current": closure.block_counts,
        "maximum_literal_word_depth": max((recipe[2] for recipe in closure.recipes), default=0),
        "sparse_discovered_d4_orbits": len(closure.worker.discovered_orbits),
        "sparse_transition_cache_entries": len(closure.worker.transition_cache),
        "basis_estimated_bytes": closure.worker.used_basis_bytes,
        "rss_start_bytes": rss_started,
        "rss_peak_bytes": rss_bytes(),
        "memory_limit_bytes": memory_limit_bytes,
        "elapsed_seconds": f"{time.monotonic() - started:.6f}",
        "checkpoint": {
            "scheme": (
                "JSON logical checkpoint: recipes rebuild literal raw words and their same-order F_p echelon "
                "basis; resumption intentionally restarts at the last completed depth, so no opaque Python state "
                "or unverified partial reducer state is trusted"
            ),
            "basis_recipes": closure.recipes,
            "tree_checksum_mod_1000000007": tree_checksum(closure.recipes),
            "completed_basis_count": completed_basis_count,
            "completed_depth": completed_depth,
            "frontier_at_completed_depth": frontier,
            "partial_depth": partial_depth,
            "partial_frontier_accepted": partial_frontier,
            "partial_candidates_processed": partial_candidates_processed,
            "resumed_from_completed_checkpoint": resumed_from_checkpoint,
        },
    }
    checks = [
        {
            "name": "C2^3_character_projector_dimensions",
            "passed": sum(symmetry["c2_character_sector_dimensions"].values()) == 1 << N,
            "detail": "exact character traces exhaust the 65536-dimensional physical module",
        },
        {
            "name": "D4_global_flip_isotypic_multiplicities",
            "passed": sum(
                symmetry["d4_global_flip_multiplicity_spaces"][parity][irrep] * (2 if irrep == "E" else 1)
                for parity in ("even", "odd")
                for irrep in ("A1", "A2", "B1", "B2", "E")
            )
            == 1 << N,
            "detail": "exact finite-group character traces exhaust the physical module",
        },
        {
            "name": "literal_modular_lower_minor",
            "passed": closure.rank >= completed_basis_count >= 2,
            "detail": "accepted raw nested Pauli words have nonzero ordered-monomial F_p pivots",
        },
    ]
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "physical_memory_bytes": physical_memory_bytes(),
        },
        "data": {
            "layer": {
                "rows": ROWS,
                "cols": COLS,
                "boundary": "open",
                "n_sites": N,
                "n_bonds": 24,
                "generators": "A=sum_v X_v; B=sum_{(u,v) in E} Z_u Z_v",
            },
            "modular_partial_closure": modular,
            "symmetry_decomposition": symmetry,
            "support_probe": support_probe,
            "exact_Q_status": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "No dynamic commutant-centre projector calculation, linked-sector factorization, or tight "
                    "same-basis rational upper squeeze was materialized for 4x4. The D4 x global-flip symmetry "
                    "container is exact but broad and does not meet the modular lower bound."
                ),
            },
        },
        "checks": checks,
    }


def load_completed_checkpoint() -> tuple[dict[str, Any], dict[str, Any]] | None:
    if not RESULT.exists():
        return None
    try:
        stored = json.loads(RESULT.read_text(encoding="utf-8"))
        if stored.get("provenance", {}).get("script") != SCRIPT:
            return None
        modular = stored["data"]["modular_partial_closure"]
        checkpoint = modular["checkpoint"]
        if modular["prime"] != GOOD_PRIME:
            return None
        if checkpoint["completed_basis_count"] > len(checkpoint["basis_recipes"]):
            return None
        return modular, checkpoint
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def support_probe(e20: Any, model: Any, cap: int) -> dict[str, Any]:
    if cap <= 0:
        return {
            "claim_tag": "[UNRESOLVED]",
            "status": "skipped_by_cli",
            "reason": "the producer was invoked with --support-probe-cap 0",
        }
    started = time.monotonic()
    try:
        support = model.enumerate_support(cap)
    except e20.SupportLimit as error:
        return {
            "claim_tag": "[COMPUTATION]",
            "status": "cap_exceeded",
            "cap_orbits": cap,
            "observed_orbits_when_stopped": error.count,
            "wall_seconds": f"{time.monotonic() - started:.6f}",
            "interpretation": (
                "This is an exact breadth-first D4-orbit support lower bound. It excludes a full dense support "
                "allocation at this cap but does not estimate the final reachable support size."
            ),
        }
    return {
        "claim_tag": "[COMPUTATION]",
        "status": "complete_below_cap",
        "cap_orbits": cap,
        "support_orbits": len(support),
        "wall_seconds": f"{time.monotonic() - started:.6f}",
    }


def run(
    *,
    wall_seconds: float,
    memory_limit_bytes: int,
    stop_after_depth: int | None,
    resume: bool,
    support_probe_cap: int,
) -> dict[str, Any]:
    if wall_seconds <= 0:
        raise ValueError("wall_seconds must be positive")
    if memory_limit_bytes <= 0:
        raise ValueError("memory_limit_bytes must be positive")
    if stop_after_depth is not None and stop_after_depth < 1:
        raise ValueError("stop_after_depth must be positive")

    e20 = load_module("_e20_for_e97", E20_PATH)
    model = e20.OrbitModel(ROWS, COLS, False)
    symmetry = symmetry_decomposition()
    started = time.monotonic()
    rss_started = rss_bytes()
    closure = SparseLiteralClosure(e20, model, memory_limit_bytes)
    dimensions: list[int]
    completed_depth: int
    completed_basis_count: int
    completed_block_counts: list[int]
    frontier: list[int]
    resumed_from_checkpoint = False

    saved = load_completed_checkpoint() if resume else None
    if saved is not None:
        modular, checkpoint = saved
        prefix_count = checkpoint["completed_basis_count"]
        prefix_recipes = checkpoint["basis_recipes"][:prefix_count]
        closure.restore_prefix(prefix_recipes)
        dimensions = [int(value) for value in modular["dimensions_through_completed_depth"]]
        completed_depth = int(checkpoint["completed_depth"])
        completed_basis_count = prefix_count
        completed_block_counts = [int(value) for value in modular["grading_block_ranks_completed_depth"]]
        frontier = [index for index, recipe in enumerate(closure.recipes) if recipe[2] == completed_depth]
        saved_frontier = [int(index) for index in checkpoint["frontier_at_completed_depth"]]
        if frontier != saved_frontier or dimensions[-1] != closure.rank:
            raise AssertionError("saved completed frontier did not reproduce")
        resumed_from_checkpoint = True
    else:
        frontier = []
        for action_index in (0, 1):
            raw_word = closure.seed(action_index)
            if closure.accept(raw_word, [-1, action_index, 1]) is None:
                raise AssertionError("initial generator seed was dependent")
            frontier.append(closure.rank - 1)
        dimensions = [closure.rank]
        completed_depth = 1
        completed_basis_count = closure.rank
        completed_block_counts = closure.block_counts

    partial_frontier: list[int] = []
    partial_depth: int | None = None
    partial_candidates_processed = 0
    stop_reason = ""
    probe_pending = {
        "claim_tag": "[UNRESOLVED]",
        "status": "pending_until_modular_closure_stops",
    }
    last_checkpoint = time.monotonic()

    def checkpoint_now() -> None:
        snapshot = artifact(
            symmetry=symmetry,
            closure=closure,
            dimensions=dimensions,
            completed_depth=completed_depth,
            completed_basis_count=completed_basis_count,
            completed_block_counts=completed_block_counts,
            frontier=frontier,
            partial_frontier=partial_frontier,
            partial_depth=partial_depth,
            partial_candidates_processed=partial_candidates_processed,
            stop_reason=stop_reason or "running_checkpoint",
            started=started,
            rss_started=rss_started,
            memory_limit_bytes=memory_limit_bytes,
            resumed_from_checkpoint=resumed_from_checkpoint,
            support_probe=probe_pending,
        )
        write_json_atomic(RESULT, snapshot)

    checkpoint_now()
    while True:
        if stop_after_depth is not None and completed_depth >= stop_after_depth:
            stop_reason = "controlled_depth_limit"
            break
        if time.monotonic() - started >= wall_seconds:
            stop_reason = "time_wall_at_depth_boundary"
            break
        if not frontier:
            stop_reason = "saturated"
            break

        next_depth = completed_depth + 1
        partial_depth = next_depth
        partial_frontier = []
        partial_candidates_processed = 0
        wall_hit = False
        for parent in frontier:
            for action_index in (0, 1):
                if time.monotonic() - started >= wall_seconds:
                    wall_hit = True
                    break
                raw_word = closure.apply_raw(action_index, closure.words[parent])
                partial_candidates_processed += 1
                if closure.accept(raw_word, [parent, action_index, next_depth]) is not None:
                    partial_frontier.append(closure.rank - 1)
                if time.monotonic() - last_checkpoint >= CHECKPOINT_INTERVAL_SECONDS:
                    checkpoint_now()
                    last_checkpoint = time.monotonic()
            if wall_hit:
                break
        if wall_hit:
            stop_reason = "time_wall_during_depth"
            break

        completed_depth = next_depth
        frontier = partial_frontier
        dimensions.append(closure.rank)
        completed_basis_count = closure.rank
        completed_block_counts = closure.block_counts
        partial_depth = None
        partial_frontier = []
        partial_candidates_processed = 0
        checkpoint_now()
        last_checkpoint = time.monotonic()
        if not frontier:
            stop_reason = "saturated"
            break

    probe = support_probe(e20, model, support_probe_cap)
    final = artifact(
        symmetry=symmetry,
        closure=closure,
        dimensions=dimensions,
        completed_depth=completed_depth,
        completed_basis_count=completed_basis_count,
        completed_block_counts=completed_block_counts,
        frontier=frontier,
        partial_frontier=partial_frontier,
        partial_depth=partial_depth,
        partial_candidates_processed=partial_candidates_processed,
        stop_reason=stop_reason,
        started=started,
        rss_started=rss_started,
        memory_limit_bytes=memory_limit_bytes,
        resumed_from_checkpoint=resumed_from_checkpoint,
        support_probe=probe,
    )
    write_json_atomic(RESULT, final)
    return final


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wall-seconds", type=float, default=DEFAULT_WALL_SECONDS)
    parser.add_argument("--memory-limit-gib", type=float, default=DEFAULT_MEMORY_LIMIT_BYTES / 1024**3)
    parser.add_argument("--stop-after-depth", type=int)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--support-probe-cap", type=int, default=DEFAULT_SUPPORT_PROBE_CAP)
    args = parser.parse_args(argv)
    artifact_data = run(
        wall_seconds=args.wall_seconds,
        memory_limit_bytes=int(args.memory_limit_gib * 1024**3),
        stop_after_depth=args.stop_after_depth,
        resume=not args.no_resume,
        support_probe_cap=args.support_probe_cap,
    )
    modular = artifact_data["data"]["modular_partial_closure"]
    print(f"LOWER_BOUND_Fp={modular['lower_bound_from_literal_word_minor']}")
    print(f"COMPLETED_DEPTH={modular['completed_depth']}")
    print(f"STOP_REASON={modular['stop_reason']}")
    print(f"WALL_SECONDS={modular['elapsed_seconds']}")
    print(f"WROTE={RESULT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
