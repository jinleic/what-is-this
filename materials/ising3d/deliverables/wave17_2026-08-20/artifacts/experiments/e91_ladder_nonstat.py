#!/usr/bin/env python3
"""Exact nonstationary ladder-word certificates and transfer obstructions.

This is intentionally self-contained: it uses exact integer Pauli expansion, exact
Fraction arithmetic only for small characteristic-zero checks, and modular minors only
as rigorous lower-bound certificates over Q.
"""
from __future__ import annotations

import gc
import hashlib
import itertools
import json
import signal
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_nonstat.json"
SCRIPT = "experiments/e91_ladder_nonstat.py"
TIMEOUT_SECONDS = 7_200
P1 = 2_147_483_647
P2 = 2_147_483_629
HASH_WIDTH = 512
HASH_MASK = (1 << 64) - 1
HASH_SEED = 0x9E3779B97F4A7C15
HASH_MULTIPLIER_1 = 0xBF58476D1CE4E5B9
HASH_MULTIPLIER_2 = 0x94D049BB133111EB

CLASSIC_SEEDS = ("A", "BA", "ABA", "BABA")
SHORT_SEEDS = ("A", "B", "AB", "AAB")
EXPECTED_SEEDS_L2 = ("A", "B", "AB", "AAB", "BAB", "ABAB")
EXPECTED_SEEDS_L3 = (
    "A",
    "B",
    "AB",
    "AAB",
    "BAB",
    "ABAB",
    "BBAB",
    "AABAB",
    "ABBAB",
    "BABAB",
    "BBBAB",
    "AABBAB",
    "ABABAB",
    "ABBBAB",
    "BABBAB",
    "BBBBAB",
)
PATCH_REMOVED_LABELS = (
    "10111011",
    "10111101",
    "10111110",
    "11111010",
    "11111100",
    "11111101",
    "11111110",
    "11111111",
)
PATCH_EXTENSION_LABELS = (
    "000000",
    "100000",
    "010000",
    "110000",
    "001000",
    "000100",
    "100100",
    "010100",
)
THREE_PAIRS = (("ABB", "BAB"), ("ABB", "BBB"), ("BAB", "BBB"))


def generator_terms(length: int) -> dict[str, tuple[int, ...]]:
    """Exact Pauli terms of A and B for the open two-leg L-rung ladder."""
    if length < 2:
        raise ValueError("the open two-leg ladder requires length >= 2")
    n = 2 * length
    bonds: list[tuple[int, int]] = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + column, offset + column + 1) for column in range(length - 1))
    bonds.extend((column, length + column) for column in range(length))
    return {
        "A": tuple(1 << site for site in range(n)),
        "B": tuple((1 << (n + left)) | (1 << (n + right)) for left, right in bonds),
    }


@lru_cache(maxsize=None)
def cached_terms(length: int) -> dict[str, tuple[int, ...]]:
    return generator_terms(length)


def sign_half(left: int, right: int, n: int) -> int:
    """Coefficient in [Q_left,Q_right]/2: 0 or +/-1."""
    mask = (1 << n) - 1
    left_right = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    right_left = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if left_right == right_left:
        return 0
    return 1 if left_right == 0 else -1


def adjoint(terms: tuple[int, ...], vector: dict[int, int], n: int) -> dict[int, int]:
    """Apply the normalized adjoint of a generator sum using exact integers."""
    result: dict[int, int] = {}
    for term in terms:
        for pauli, coefficient in vector.items():
            sign = sign_half(term, pauli, n)
            if not sign:
                continue
            target = term ^ pauli
            updated = result.get(target, 0) + sign * coefficient
            if updated:
                result[target] = updated
            else:
                result.pop(target, None)
    return result


def vector_for_word(word: str, length: int) -> dict[int, int]:
    """Exact vector for an outermost-first generator-left-nested word."""
    terms = cached_terms(length)
    value = {term: 1 for term in terms[word[-1]]}
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


@dataclass
class ModularEchelon:
    """Sparse row echelon over a prime field, with a deterministic pivot order."""

    prime: int

    def __post_init__(self) -> None:
        self.pivots: dict[int, dict[int, int]] = {}

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def add(self, row: dict[int, int]) -> bool:
        value = {key: coefficient % self.prime for key, coefficient in row.items() if coefficient % self.prime}
        while value:
            lead = min(value)
            pivot = self.pivots.get(lead)
            if pivot is None:
                inverse = pow(value[lead], self.prime - 2, self.prime)
                self.pivots[lead] = {
                    key: coefficient * inverse % self.prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % self.prime
                }
                return True
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = (value.get(key, 0) - factor * coefficient) % self.prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
        return False


class ExactEchelon:
    """Sparse rational echelon used only where a characteristic-zero upper bound is small."""

    def __init__(self) -> None:
        self.pivots: dict[int, dict[int, Fraction]] = {}

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def add(self, row: dict[int, int]) -> bool:
        value = {key: Fraction(coefficient) for key, coefficient in row.items() if coefficient}
        while value:
            lead = min(value)
            pivot = self.pivots.get(lead)
            if pivot is None:
                inverse = Fraction(1, 1) / value[lead]
                self.pivots[lead] = {key: coefficient * inverse for key, coefficient in value.items()}
                return True
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = value.get(key, Fraction(0)) - factor * coefficient
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
        return False


def rank_mod_both(rows: list[dict[int, int]]) -> tuple[int, int]:
    first = ModularEchelon(P1)
    second = ModularEchelon(P2)
    for row in rows:
        first.add(row)
        second.add(row)
    return first.rank, second.rank


def vector_digest(vector: dict[int, int]) -> str:
    digest = hashlib.sha256()
    for code, coefficient in sorted(vector.items()):
        digest.update(f"{code}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def words_digest(words: list[str] | tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for word in words:
        digest.update(f"{word};".encode("ascii"))
    return digest.hexdigest()


def all_literal_words(max_depth: int) -> list[str]:
    return [
        "".join(letters)
        for depth in range(1, max_depth + 1)
        for letters in itertools.product("AB", repeat=depth)
    ]


def exact_seed_basis(length: int, max_depth: int) -> tuple[list[str], int]:
    echelon = ExactEchelon()
    selected: list[str] = []
    for word in all_literal_words(max_depth):
        if echelon.add(vector_for_word(word, length)):
            selected.append(word)
    return selected, echelon.rank


def expand_pair(roots: tuple[str, ...] | list[str], moves: tuple[str, str], steps: int) -> list[str]:
    words = list(roots)
    for _ in range(steps):
        words = [move + word for word in words for move in moves]
    return words


def original_recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def patch_word(extension_label: str) -> str:
    word = "BAB"
    for bit in extension_label:
        word = ("AB" if bit == "0" else "BB") + word
    return word


def patched_words_l8() -> list[str]:
    removed = set(PATCH_REMOVED_LABELS)
    original = [word for label, word in original_recursive_words(8) if label not in removed]
    words = original + [patch_word(label) for label in PATCH_EXTENSION_LABELS]
    if len(words) != 256 or max(map(len, words)) > 16:
        raise AssertionError("invalid L=8 patched family")
    return words


def p8_input_certificate(words: list[str]) -> dict[str, object]:
    started = time.monotonic()
    rows = [vector_for_word(word, 8) for word in words]
    rank_first, rank_second = rank_mod_both(rows)
    if (rank_first, rank_second) != (256, 256):
        raise AssertionError((rank_first, rank_second))
    return {
        "claim_tag": "[COMPUTATION]",
        "length": 8,
        "word_count": len(words),
        "max_depth": max(map(len, words)),
        "ranks_mod": {str(P1): rank_first, str(P2): rank_second},
        "word_digest": words_digest(words),
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def p9_transition_obstruction(words: list[str]) -> dict[str, object]:
    """Exact modular data for P8 descendants and all BAB-tail repair candidates."""
    started = time.monotonic()
    terms = cached_terms(9)
    first = ModularEchelon(P1)
    second = ModularEchelon(P2)
    seen: set[str] = set()
    max_support = 0
    base_count = 0

    def add(label: str, row: dict[int, int]) -> None:
        nonlocal max_support
        if label in seen:
            return
        seen.add(label)
        max_support = max(max_support, len(row))
        first.add(row)
        second.add(row)

    for parent in words:
        parent_vector = vector_for_word(parent, 9)
        after_b = adjoint(terms["B"], parent_vector, 18)
        add("AB" + parent, adjoint(terms["A"], after_b, 18))
        add("BB" + parent, adjoint(terms["B"], after_b, 18))
        base_count += 2
        del parent_vector, after_b
    base_ranks = (first.rank, second.rank)
    if base_ranks != (478, 478):
        raise AssertionError(f"unexpected P8 child ranks {base_ranks}")

    extra_count = 0
    for value in range(1 << 7):
        bits = f"{value:07b}"
        candidate = patch_word(bits)
        if candidate in seen:
            continue
        add(candidate, vector_for_word(candidate, 9))
        extra_count += 1
    combined_ranks = (first.rank, second.rank)
    if combined_ranks != (499, 499):
        raise AssertionError(f"unexpected BAB-repair ranks {combined_ranks}")
    return {
        "claim_tag": "[COMPUTATION]",
        "length": 9,
        "base_definition": "{AB w, BB w : w in patched P8}",
        "base_word_count": base_count,
        "base_ranks_mod": {str(P1): base_ranks[0], str(P2): base_ranks[1]},
        "repair_definition": "all unused seven-move AB/BB descendants of BAB",
        "repair_words_added": extra_count,
        "combined_word_count": base_count + extra_count,
        "combined_ranks_mod": {str(P1): combined_ranks[0], str(P2): combined_ranks[1]},
        "max_depth": 18,
        "max_full_pauli_support": max_support,
        "status": "MODULAR_DEFICIT_ONLY_NOT_A_Q_UPPER_BOUND",
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def six_seed_table(seeds: tuple[str, ...]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for length in range(2, 9):
        started = time.monotonic()
        words = expand_pair(seeds, ("AB", "BB"), length - 2)
        rows = [vector_for_word(word, length) for word in words]
        rank_first, rank_second = rank_mod_both(rows)
        records.append(
            {
                "length": length,
                "word_count": len(words),
                "max_depth": max(map(len, words)),
                "ranks_mod": {str(P1): rank_first, str(P2): rank_second},
                "target": 1 << length,
                "max_full_pauli_support": max(map(len, rows)),
                "elapsed_seconds": round(time.monotonic() - started, 6),
            }
        )
        del rows
        gc.collect()
    expected = (6, 10, 20, 40, 77, 148, 272)
    observed = tuple(record["ranks_mod"][str(P1)] for record in records)
    if observed != expected:
        raise AssertionError((observed, expected))
    return records


def alternating_words(
    roots: tuple[str, ...], pair: tuple[str, str], length: int, triple_first: bool
) -> list[str]:
    words = list(roots)
    for current_length in range(2, length):
        use_triple = (current_length % 2 == 0) if triple_first else (current_length % 2 == 1)
        moves = pair if use_triple else ("A", "B")
        words = [move + word for word in words for move in moves]
    return words


def three_letter_tables() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for orientation, roots, triple_first in (
        ("one_then_three", CLASSIC_SEEDS, False),
        ("three_then_one", SHORT_SEEDS, True),
    ):
        for pair in THREE_PAIRS:
            table: list[dict[str, object]] = []
            exact_l3: int | None = None
            for length in range(2, 9):
                started = time.monotonic()
                words = alternating_words(roots, pair, length, triple_first)
                rows = [vector_for_word(word, length) for word in words]
                rank_first, rank_second = rank_mod_both(rows)
                if length == 3:
                    exact = ExactEchelon()
                    for row in rows:
                        exact.add(row)
                    exact_l3 = exact.rank
                table.append(
                    {
                        "length": length,
                        "word_count": len(words),
                        "max_depth": max(map(len, words)),
                        "ranks_mod": {str(P1): rank_first, str(P2): rank_second},
                        "target": 1 << length,
                        "elapsed_seconds": round(time.monotonic() - started, 6),
                    }
                )
                del rows
                gc.collect()
            if exact_l3 is None or exact_l3 >= 8:
                raise AssertionError((orientation, pair, exact_l3))
            records.append(
                {
                    "claim_tag": "[FALSIFIED]",
                    "orientation": orientation,
                    "roots": list(roots),
                    "three_letter_pair": list(pair),
                    "exact_rank_Q_at_L3": exact_l3,
                    "failure_statement": "The eight L=3 words are linearly dependent over Q, so this literal period-two family cannot supply the requested all-L injection.",
                    "table": table,
                }
            )
    return records


def mix64(value: int) -> int:
    value = (value ^ HASH_SEED) & HASH_MASK
    value = ((value ^ (value >> 30)) * HASH_MULTIPLIER_1) & HASH_MASK
    value = ((value ^ (value >> 27)) * HASH_MULTIPLIER_2) & HASH_MASK
    return value ^ (value >> 31)


def hash_projection_row(vector: dict[int, int], cache: dict[int, tuple[int, int]]) -> list[int]:
    """An integer linear map from Pauli coordinates to HASH_WIDTH buckets."""
    result = [0] * HASH_WIDTH
    for code, coefficient in vector.items():
        entry = cache.get(code)
        if entry is None:
            mixed = mix64(code)
            entry = (mixed & (HASH_WIDTH - 1), 1 if ((mixed >> 9) & 1) else -1)
            cache[code] = entry
        bucket, sign = entry
        result[bucket] += sign * coefficient
    return result


class HashEchelon:
    """Dense modular echelon for a fixed 512-coordinate integer hash projection."""

    def __init__(self, prime: int) -> None:
        self.prime = prime
        self.rows: list[np.ndarray] = []
        self.pivot_row: dict[int, int] = {}

    @property
    def rank(self) -> int:
        return len(self.rows)

    def add(self, row: list[int]) -> bool:
        value = np.asarray([coefficient % self.prime for coefficient in row], dtype=np.int64)
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                return False
            lead = int(nonzero[0])
            existing = self.pivot_row.get(lead)
            if existing is None:
                inverse = pow(int(value[lead]), self.prime - 2, self.prime)
                value = value * inverse % self.prime
                self.pivot_row[lead] = len(self.rows)
                self.rows.append(value.astype(np.int32))
                return True
            factor = int(value[lead])
            value = (value - factor * self.rows[existing].astype(np.int64)) % self.prime


def hash_matrix_digest(rows: list[list[int]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(",".join(str(value) for value in row).encode("ascii"))
        digest.update(b";")
    return digest.hexdigest()


def l9_sixteen_seed_certificate(roots: tuple[str, ...]) -> dict[str, object]:
    """Find a same-word 512-row two-prime hash-minor certificate at L=9."""
    started = time.monotonic()
    terms = cached_terms(9)
    first = HashEchelon(P1)
    second = HashEchelon(P2)
    selected_second = HashEchelon(P2)
    hash_cache: dict[int, tuple[int, int]] = {}
    selected_words: list[str] = []
    selected_rows: list[list[int]] = []
    root_progress: list[dict[str, object]] = []
    leaves = 0
    max_support = 0
    stopped = False

    def visit(word: str, vector: dict[int, int], steps_remaining: int) -> None:
        nonlocal leaves, max_support, stopped
        if stopped:
            return
        if steps_remaining == 0:
            projection = hash_projection_row(vector, hash_cache)
            added_first = first.add(projection)
            second.add(projection)
            leaves += 1
            max_support = max(max_support, len(vector))
            if added_first:
                selected_words.append(word)
                selected_rows.append(projection)
                selected_second.add(projection)
            if (
                first.rank >= HASH_WIDTH
                and second.rank >= HASH_WIDTH
                and selected_second.rank >= HASH_WIDTH
            ):
                stopped = True
            return
        after_b = adjoint(terms["B"], vector, 18)
        visit("AB" + word, adjoint(terms["A"], after_b, 18), steps_remaining - 1)
        visit("BB" + word, adjoint(terms["B"], after_b, 18), steps_remaining - 1)

    for root in roots:
        if stopped:
            break
        before = first.rank
        visit(root, vector_for_word(root, 9), 6)
        root_progress.append(
            {
                "root": root,
                "rank_increment_mod_p1": first.rank - before,
                "rank_mod_p1": first.rank,
                "rank_mod_p2": second.rank,
                "leaves_visited": leaves,
            }
        )
    if not stopped or len(selected_words) != HASH_WIDTH:
        raise AssertionError((stopped, first.rank, second.rank, selected_second.rank, len(selected_words)))
    if max(map(len, selected_words)) > 18:
        raise AssertionError("depth violation")
    return {
        "claim_tag": "[COMPUTATION]",
        "length": 9,
        "target_rank": HASH_WIDTH,
        "roots_at_L3": list(roots),
        "tail_definition": "six outer AB/BB moves, depth-first AB-before-BB order",
        "leaves_visited": leaves,
        "root_progress": root_progress,
        "selected_word_count": len(selected_words),
        "selected_words": selected_words,
        "selected_words_digest": words_digest(selected_words),
        "max_depth": max(map(len, selected_words)),
        "max_full_pauli_support_seen": max_support,
        "hash_projection": {
            "width": HASH_WIDTH,
            "seed_hex": f"0x{HASH_SEED:016x}",
            "multiplier_1_hex": f"0x{HASH_MULTIPLIER_1:016x}",
            "multiplier_2_hex": f"0x{HASH_MULTIPLIER_2:016x}",
            "bucket": "mix64(pauli_code) mod 512",
            "sign": "+1 iff bit 9 of mix64(pauli_code) is one; otherwise -1",
            "selected_matrix_integer_digest": hash_matrix_digest(selected_rows),
            "rank_mod": {
                str(P1): first.rank,
                str(P2): selected_second.rank,
            },
            "all_visited_rank_mod": {str(P1): first.rank, str(P2): second.rank},
        },
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def remap_without_last_rung(mask: int, length: int) -> int:
    """Map old ladder sites after deleting the final rung, preserving row-major order."""
    result = 0
    for column in range(length - 1):
        if mask & (1 << column):
            result |= 1 << column
        if mask & (1 << (length + column)):
            result |= 1 << ((length - 1) + column)
    return result


def touches_last_rung(mask: int, length: int) -> bool:
    return bool(mask & ((1 << (length - 1)) | (1 << (2 * length - 1))))


def erase_last_rung(vector: dict[int, int], length: int) -> dict[int, int]:
    """Conditional coordinate projection to Pauli terms identity on the last rung."""
    n = 2 * length
    mask = (1 << n) - 1
    result: dict[int, int] = {}
    for code, coefficient in vector.items():
        x_mask = code & mask
        z_mask = code >> n
        if touches_last_rung(x_mask, length) or touches_last_rung(z_mask, length):
            continue
        new_code = remap_without_last_rung(x_mask, length) | (
            remap_without_last_rung(z_mask, length) << (n - 2)
        )
        result[new_code] = result.get(new_code, 0) + coefficient
    return {code: coefficient for code, coefficient in result.items() if coefficient}


def pauli_string(code: int, length: int) -> str:
    n = 2 * length
    x_mask = code & ((1 << n) - 1)
    z_mask = code >> n
    local = {(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}
    rows = []
    for row in range(2):
        rows.append(
            "".join(
                local[((x_mask >> (row * length + column)) & 1, (z_mask >> (row * length + column)) & 1)]
                for column in range(length)
            )
        )
    return "/".join(rows)


def restriction_obstruction() -> dict[str, object]:
    """Exact counterexamples to the tempting delete-a-rung transfer projection."""
    # Generic Pauli witness for old length 2 inside length 3.
    old_length = 2
    new_length = old_length + 1
    new_n = 2 * new_length
    old_last_top = old_length - 1
    new_top = old_length
    witness = (1 << old_last_top) | (1 << (new_n + new_top))  # X_v Z_u
    generic_left = erase_last_rung(adjoint(cached_terms(new_length)["B"], {witness: 1}, new_n), new_length)
    generic_right = adjoint(cached_terms(old_length)["B"], erase_last_rung({witness: 1}, new_length), 2 * old_length)
    expected = (1 << old_last_top) | (1 << ((2 * old_length) + old_last_top))
    if generic_right or generic_left != {expected: -1}:
        raise AssertionError((generic_left, generic_right, expected))

    word = "ABABAB"
    projected = erase_last_rung(vector_for_word(word, 3), 3)
    shorter = vector_for_word(word, 2)
    difference = {
        code: projected.get(code, 0) - shorter.get(code, 0)
        for code in set(projected) | set(shorter)
        if projected.get(code, 0) != shorter.get(code, 0)
    }
    if not difference:
        raise AssertionError("expected literal-word restriction failure")
    return {
        "claim_tag": "[LEMMA]",
        "generic_witness": {
            "old_length": old_length,
            "new_length": new_length,
            "P": "X_v Z_u, with v the old top boundary site and u the added top site",
            "erase_ad_B_P": {"code": expected, "pauli": pauli_string(expected, old_length), "coefficient": -1},
            "ad_B_erase_P": {},
        },
        "literal_word_counterexample": {
            "word": word,
            "from_length": 3,
            "to_length": 2,
            "projected_support_size": len(projected),
            "shorter_support_size": len(shorter),
            "difference_support_size": len(difference),
            "difference_digest": vector_digest(difference),
            "difference_terms": [
                {"code": code, "pauli": pauli_string(code, 2), "coefficient": coefficient}
                for code, coefficient in sorted(difference.items())
            ],
        },
        "conclusion": "Erase-last-rung is not an adjoint intertwiner, so it cannot by itself furnish the proposed suffix/block-triangular L-to-L+1 injection.",
    }


def make_artifact() -> dict[str, object]:
    started = time.monotonic()
    seeds_l2, rank_l2_q = exact_seed_basis(2, 4)
    seeds_l3, rank_l3_q = exact_seed_basis(3, 6)
    if tuple(seeds_l2) != EXPECTED_SEEDS_L2 or rank_l2_q != 6:
        raise AssertionError((seeds_l2, rank_l2_q))
    if tuple(seeds_l3) != EXPECTED_SEEDS_L3 or rank_l3_q != 16:
        raise AssertionError((seeds_l3, rank_l3_q))

    restriction = restriction_obstruction()
    l9_certificate = l9_sixteen_seed_certificate(tuple(seeds_l3))
    six_table = six_seed_table(tuple(seeds_l2))
    p8_words = patched_words_l8()
    p8_certificate = p8_input_certificate(p8_words)
    transition = p9_transition_obstruction(p8_words)
    triples = three_letter_tables()

    return {
        "schema_version": 1,
        "status": "FINITE_L9_CERTIFICATE_ALL_L_UNRESOLVED",
        "target": "D_(2L)(open 2xL ladder) >= 2^L",
        "claims": [
            {
                "claim_tag": "[COMPUTATION]",
                "statement": "D_18(open 2x9 ladder) >= 512 = 2^9.",
                "basis": "The selected 512 literal depth-at-most-18 words have a 512-by-512 deterministic integer hash projection minor nonzero modulo both recorded primes.",
            },
            {
                "claim_tag": "[UNRESOLVED]",
                "statement": "No all-L injective transition or unit-determinant transfer matrix is proved.",
            },
        ],
        "exact_seed_spaces": {
            "L2_depth_le_4": {"rank_Q": rank_l2_q, "basis_words": seeds_l2},
            "L3_depth_le_6": {"rank_Q": rank_l3_q, "basis_words": seeds_l3},
        },
        "restriction_obstruction": restriction,
        "l9_certificate": l9_certificate,
        "six_seed_family": {
            "claim_tag": "[COMPUTATION]",
            "roots_at_L2": seeds_l2,
            "moves": ["AB", "BB"],
            "table": six_table,
        },
        "patched_l8_certificate": p8_certificate,
        "patch_transition": transition,
        "three_letter_families": triples,
        "provenance": {
            "script": SCRIPT,
            "interpreter": ".venv/bin/python",
            "primes": [P1, P2],
            "arithmetic": "Exact integer Pauli expansion; modular nonzero minors are rational lower certificates; Fraction echelon supplies the listed small Q ranks.",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.monotonic() - started, 6),
        },
    }


def write_artifact(artifact: dict[str, object]) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = RESULT_PATH.with_suffix(".json.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(RESULT_PATH)


def on_timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s")


def main() -> int:
    signal.signal(signal.SIGALRM, on_timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        artifact = make_artifact()
        write_artifact(artifact)
    finally:
        signal.alarm(0)
    l9 = artifact["l9_certificate"]
    print(
        "PASS: L9 hash certificate "
        f"rank {l9['hash_projection']['rank_mod'][str(P1)]}/"
        f"{l9['hash_projection']['rank_mod'][str(P2)]}; "
        f"artifact={RESULT_PATH.relative_to(ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
