#!/usr/bin/env python3
"""Clean-room verifier for the open 4x4 Ising-layer partial certificate.

This verifier deliberately imports neither ``experiments/e97_algebra_4x4.py``
nor its closure helpers.  It rebuilds D4 orbit Pauli actions, exact sparse
modular elimination, C2^3 character-projector ranks, and D4 x global-flip
multiplicities from the mathematical definitions.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_structure" / "char0_4x4.json"
SCRIPT = "experiments/e97_algebra_4x4.py"
PRIME = 2_147_483_647
ROWS = COLS = 4
N = ROWS * COLS
MASK = (1 << N) - 1
BONDS = tuple(
    (row * COLS + column, (row + 1) * COLS + column)
    for row in range(ROWS - 1)
    for column in range(COLS)
) + tuple(
    (row * COLS + column, row * COLS + column + 1)
    for row in range(ROWS)
    for column in range(COLS - 1)
)


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not condition:
        raise AssertionError(name)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(N))


def group_from_generators(generators: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    identity = tuple(range(N))
    queue = [identity]
    seen = {identity}
    head = 0
    while head < len(queue):
        current = queue[head]
        head += 1
        for generator in generators:
            product = compose(generator, current)
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return tuple(queue)


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


def permute_mask(value: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for old, new in enumerate(permutation):
        if (value >> old) & 1:
            result |= 1 << new
    return result


def permute_pauli(value: int, permutation: tuple[int, ...]) -> int:
    return permute_mask(value & MASK, permutation) | (permute_mask(value >> N, permutation) << N)


def commutator_sign_half(generator: int, value: int) -> int:
    first = (((generator >> N) & (value & MASK)).bit_count()) & 1
    second = (((value >> N) & (generator & MASK)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def pauli_grade(value: int) -> int:
    x_mask = value & MASK
    z_mask = value >> N
    return ((x_mask.bit_count() & 1) << 1) | ((x_mask & z_mask).bit_count() & 1)


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


def c2_sector_dimensions() -> dict[str, int]:
    row_reflection, column_reflection, _ = geometric_generators()
    elements = []
    for row_bit in (0, 1):
        for column_bit in (0, 1):
            permutation = tuple(range(N))
            if row_bit:
                permutation = compose(row_reflection, permutation)
            if column_bit:
                permutation = compose(column_reflection, permutation)
            for flip_bit in (0, 1):
                elements.append((row_bit, column_bit, flip_bit, permutation))

    result: dict[str, int] = {}
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
                assert numerator % 8 == 0
                result[f"{row_character}{column_character}{flip_character}"] = numerator // 8
    return result


def d4_class_map() -> dict[tuple[int, ...], str]:
    row_reflection, _column_reflection, transpose = geometric_generators()
    rotation = compose(transpose, row_reflection)
    identity = tuple(range(N))
    rotation_squared = compose(rotation, rotation)
    rotation_cubed = compose(rotation, rotation_squared)
    group = group_from_generators((row_reflection, transpose))
    result: dict[tuple[int, ...], str] = {}
    for permutation in group:
        if permutation == identity:
            result[permutation] = "identity"
        elif permutation == rotation_squared:
            result[permutation] = "rotation_180"
        elif permutation in (rotation, rotation_cubed):
            result[permutation] = "rotation_quarter"
        else:
            result[permutation] = "reflection_diagonal" if sum(index == target for index, target in enumerate(permutation)) else "reflection_axis"
    assert {value for value in result.values()} == {
        "identity",
        "rotation_180",
        "rotation_quarter",
        "reflection_axis",
        "reflection_diagonal",
    }
    return result


def d4_flip_multiplicities() -> dict[str, dict[str, int]]:
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
        trace = {
            permutation: (fixed_spin_count(permutation, 0) + parity_sign * fixed_spin_count(permutation, 1)) // 2
            for permutation in group
        }
        multiplicities: dict[str, int] = {}
        for irrep, character in character_table.items():
            numerator = sum(character[classes[permutation]] * trace[permutation] for permutation in group)
            assert numerator % 8 == 0
            multiplicities[irrep] = numerator // 8
        result[parity_name] = multiplicities
    return result


class IndependentSparseClosure:
    """D4-orbit Pauli closure with a different BFS word order from the producer."""

    def __init__(self, prime: int) -> None:
        self.prime = prime
        row_reflection, column_reflection, transpose = geometric_generators()
        self.group = group_from_generators((transpose, column_reflection, row_reflection))
        self.cache: dict[int, int] = {}
        self.transition_cache: dict[tuple[int, int], dict[int, int]] = {}
        self.x_terms = tuple(1 << site for site in range(N))
        self.zz_terms = tuple(((1 << first) | (1 << second)) << N for first, second in BONDS)
        self.basis: dict[int, dict[int, int]] = {}
        self.inverse: dict[int, int] = {}
        self.accepted_words: list[dict[int, int]] = []
        self.block_counts = [0, 0, 0, 0]

    def canonical(self, value: int) -> int:
        cached = self.cache.get(value)
        if cached is None:
            cached = min(permute_pauli(value, permutation) for permutation in self.group)
            self.cache[value] = cached
        return cached

    def seed(self, action_index: int) -> dict[int, int]:
        terms = self.x_terms if action_index == 0 else self.zz_terms
        term_multiplicities: dict[int, int] = {}
        for term in terms:
            term_multiplicities[term] = term_multiplicities.get(term, 0) + 1
        coefficients: dict[int, int] = {}
        for term, multiplicity in term_multiplicities.items():
            orbit = self.canonical(term)
            previous = coefficients.setdefault(orbit, multiplicity)
            assert previous == multiplicity
        assert len({pauli_grade(orbit) for orbit in coefficients}) == 1
        return coefficients

    def transitions(self, source: int, action_index: int) -> dict[int, int]:
        key = (source, action_index)
        cached = self.transition_cache.get(key)
        if cached is not None:
            return cached
        terms = self.x_terms if action_index == 0 else self.zz_terms
        targets = {
            self.canonical(source ^ generator)
            for generator in terms
            if commutator_sign_half(generator, source)
        }
        result: dict[int, int] = {}
        for target in targets:
            coefficient = 0
            for generator in terms:
                preimage = target ^ generator
                sign = commutator_sign_half(generator, preimage)
                if sign and self.canonical(preimage) == source:
                    coefficient += sign
            if coefficient:
                result[target] = coefficient
        self.transition_cache[key] = result
        return result

    def apply_raw(self, action_index: int, word: dict[int, int]) -> dict[int, int]:
        result: dict[int, int] = {}
        for source, value in word.items():
            for target, coefficient in self.transitions(source, action_index).items():
                updated = (result.get(target, 0) + value * coefficient) % self.prime
                if updated:
                    result[target] = updated
                else:
                    result.pop(target, None)
        return result

    def reduce_add(self, raw_value: dict[int, int]) -> bool:
        value = {key: coefficient % self.prime for key, coefficient in raw_value.items() if coefficient % self.prime}
        while value:
            lead = min(value)
            existing = self.basis.get(lead)
            if existing is None:
                self.basis[lead] = value
                self.inverse[lead] = pow(value[lead], self.prime - 2, self.prime)
                self.block_counts[pauli_grade(lead)] += 1
                self.accepted_words.append(raw_value)
                return True
            factor = value[lead] * self.inverse[lead] % self.prime
            for key, coefficient in existing.items():
                updated = (value.get(key, 0) - factor * coefficient) % self.prime
                if updated:
                    value[key] = updated
                else:
                    value.pop(key, None)
        return False

    def replay_recipes(self, recipes: list[list[int]]) -> None:
        for index, recipe in enumerate(recipes):
            parent, action_index, depth = recipe
            if action_index not in (0, 1) or depth < 1:
                raise AssertionError(f"invalid recipe shape at {index}")
            if parent == -1:
                if depth != 1:
                    raise AssertionError(f"non-depth-one seed at {index}")
                raw_word = self.seed(action_index)
            else:
                if not 0 <= parent < index:
                    raise AssertionError(f"invalid parent at {index}")
                raw_word = self.apply_raw(action_index, self.accepted_words[parent])
                if depth != recipes[parent][2] + 1:
                    raise AssertionError(f"invalid recursive depth at {index}")
            if not self.reduce_add(raw_word):
                raise AssertionError(f"dependent literal word at {index}")

    def reverse_closure(self, target_depth: int) -> tuple[list[int], list[int]]:
        frontier: list[int] = []
        for action_index in (1, 0):
            raw_word = self.seed(action_index)
            check(f"reverse seed {action_index}", self.reduce_add(raw_word))
            frontier.append(len(self.accepted_words) - 1)
        dimensions = [len(self.accepted_words)]
        completed_depth = 1
        while completed_depth < target_depth:
            next_frontier: list[int] = []
            for parent in reversed(frontier):
                for action_index in (1, 0):
                    raw_word = self.apply_raw(action_index, self.accepted_words[parent])
                    if self.reduce_add(raw_word):
                        next_frontier.append(len(self.accepted_words) - 1)
            completed_depth += 1
            dimensions.append(len(self.accepted_words))
            frontier = next_frontier
            if not frontier:
                break
        return dimensions, self.block_counts


def support_probe_count(cap: int) -> int:
    """Independently replay the exact D4-orbit support BFS through ``cap``."""
    closure = IndependentSparseClosure(PRIME)
    seen = {closure.canonical(term) for term in closure.x_terms + closure.zz_terms}
    queue = sorted(seen)
    head = 0
    while head < len(queue):
        value = queue[head]
        head += 1
        for generator in closure.x_terms + closure.zz_terms:
            if commutator_sign_half(generator, value):
                target = closure.canonical(value ^ generator)
                if target not in seen:
                    seen.add(target)
                    queue.append(target)
                    if len(seen) > cap:
                        return len(seen)
    return len(seen)


def main() -> int:
    started = time.monotonic()
    data = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(data) == {"provenance", "data", "checks"}
        and data["provenance"]["script"] == SCRIPT,
    )
    payload = data["data"]
    modular = payload["modular_partial_closure"]
    checkpoint = modular["checkpoint"]
    check(
        "prime and raw-word lower-bound scope",
        modular["claim_tag"] == "[COMPUTATION]"
        and modular["prime"] == PRIME
        and "literal nested Pauli" in modular["arithmetic"],
    )
    check(
        "4x4 layer definition",
        payload["layer"]["rows"] == ROWS
        and payload["layer"]["cols"] == COLS
        and payload["layer"]["n_sites"] == N
        and payload["layer"]["n_bonds"] == len(BONDS) == 24,
    )
    support = payload["support_probe"]
    if support["status"] == "cap_exceeded":
        observed = support_probe_count(support["cap_orbits"])
        check(
            "exact D4 support-cap frontier",
            observed == support["observed_orbits_when_stopped"]
            and observed == support["cap_orbits"] + 1,
        )
    else:
        check(
            "support-probe status is explicitly nonclaiming",
            support["claim_tag"] == "[UNRESOLVED]",
        )


    expected_sectors = c2_sector_dimensions()
    expected_multiplicities = d4_flip_multiplicities()
    symmetry = payload["symmetry_decomposition"]
    check(
        "exact C2^3 projector ranks",
        symmetry["c2_character_sector_dimensions"] == expected_sectors
        and sum(expected_sectors.values()) == 1 << N,
    )
    check(
        "exact D4 x global-flip multiplicities",
        symmetry["d4_global_flip_multiplicity_spaces"] == expected_multiplicities,
    )
    candidate_dimensions = [
        expected_multiplicities[parity][irrep]
        for parity in ("even", "odd")
        for irrep in ("A1", "A2", "B1", "B2", "E")
    ]
    check(
        "symmetry-container arithmetic",
        symmetry["full_special_linear_candidate_dimensions"] == candidate_dimensions
        and symmetry["unlinked_endomorphism_container_dimension_Q"] == sum(value * value for value in candidate_dimensions)
        and symmetry["one_scalar_linked_sl_candidate_dimension_Q"]
        == 1 + sum(value * value - 1 for value in candidate_dimensions),
    )

    recipes = checkpoint["basis_recipes"]
    check(
        "checkpoint field consistency",
        checkpoint["completed_basis_count"] <= len(recipes)
        and modular["lower_bound_from_literal_word_minor"] == len(recipes)
        and checkpoint["completed_depth"] == modular["completed_depth"],
    )
    replay = IndependentSparseClosure(PRIME)
    replay.replay_recipes(recipes)
    check(
        "stored partial frontier independently replays",
        len(replay.accepted_words) == modular["lower_bound_from_literal_word_minor"]
        and replay.block_counts == modular["grading_block_ranks_current"],
    )

    reverse = IndependentSparseClosure(PRIME)
    reverse_dimensions, reverse_blocks = reverse.reverse_closure(modular["completed_depth"])
    check(
        "reverse-order completed-depth closure",
        reverse_dimensions == modular["dimensions_through_completed_depth"]
        and reverse_blocks == modular["grading_block_ranks_completed_depth"],
    )
    check(
        "completed checkpoint prefix",
        checkpoint["completed_basis_count"]
        == modular["dimensions_through_completed_depth"][-1],
    )
    print(f"WALL_SECONDS={time.monotonic() - started:.3f}")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}", file=sys.stderr)
        raise
