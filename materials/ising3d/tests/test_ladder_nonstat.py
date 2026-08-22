#!/usr/bin/env python3
"""Clean-room verifier for the ladder nonstationary research certificate.

This file intentionally imports neither experiments/e91_ladder_nonstat.py nor any
producer-only helper.  It reconstructs the open-ladder Pauli algebra from raw graph
terms and validates every claimed rank/injection in the artifact.
"""
from __future__ import annotations

import gc
import hashlib
import itertools
import json
import signal
import time
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_growth" / "ladder_nonstat.json"
P1 = 2_147_483_647
P2 = 2_147_483_629
WIDTH = 512
MASK64 = (1 << 64) - 1
SEED = 0x9E3779B97F4A7C15
MULTIPLIER_1 = 0xBF58476D1CE4E5B9
MULTIPLIER_2 = 0x94D049BB133111EB
TEST_TIMEOUT_SECONDS = 10_800

CLASSIC = ("A", "BA", "ABA", "BABA")
SHORT = ("A", "B", "AB", "AAB")
SEEDS_L2 = ("A", "B", "AB", "AAB", "BAB", "ABAB")
SEEDS_L3 = (
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
REMOVED = {
    "10111011",
    "10111101",
    "10111110",
    "11111010",
    "11111100",
    "11111101",
    "11111110",
    "11111111",
}
PATCHES = ("000000", "100000", "010000", "110000", "001000", "000100", "100100", "010100")
TRIPLES = (("ABB", "BAB"), ("ABB", "BBB"), ("BAB", "BBB"))


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not condition:
        raise AssertionError(name)


def terms(length: int) -> dict[str, tuple[int, ...]]:
    n = 2 * length
    bonds = []
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
    return terms(length)


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    first = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    second = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def ad(generator: tuple[int, ...], vector: dict[int, int], n: int) -> dict[int, int]:
    result: dict[int, int] = {}
    for term in generator:
        for code, coefficient in vector.items():
            sign = sign_half(term, code, n)
            if not sign:
                continue
            target = term ^ code
            value = result.get(target, 0) + sign * coefficient
            if value:
                result[target] = value
            else:
                result.pop(target, None)
    return result


def word_vector(word: str, length: int) -> dict[int, int]:
    generator = cached_terms(length)
    result = {term: 1 for term in generator[word[-1]]}
    for letter in reversed(word[:-1]):
        result = ad(generator[letter], result, 2 * length)
    return result


class ModEchelon:
    def __init__(self, prime: int) -> None:
        self.prime = prime
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


class QEchelon:
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


def ranks(rows: list[dict[int, int]]) -> tuple[int, int]:
    first = ModEchelon(P1)
    second = ModEchelon(P2)
    for row in rows:
        first.add(row)
        second.add(row)
    return first.rank, second.rank


def literal_words(max_depth: int) -> list[str]:
    return [
        "".join(bits)
        for depth in range(1, max_depth + 1)
        for bits in itertools.product("AB", repeat=depth)
    ]


def exact_basis(length: int, max_depth: int) -> tuple[list[str], int]:
    echelon = QEchelon()
    output = []
    for word in literal_words(max_depth):
        if echelon.add(word_vector(word, length)):
            output.append(word)
    return output, echelon.rank


def patch_word(bits: str) -> str:
    word = "BAB"
    for bit in bits:
        word = ("AB" if bit == "0" else "BB") + word
    return word


def original_l8() -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, 8):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def p8_words() -> list[str]:
    words = [word for label, word in original_l8() if label not in REMOVED]
    words.extend(patch_word(bits) for bits in PATCHES)
    check("P8 shape", len(words) == 256 and max(map(len, words)) <= 16)
    return words


def expand(roots: tuple[str, ...], moves: tuple[str, str], count: int) -> list[str]:
    words = list(roots)
    for _ in range(count):
        words = [move + word for word in words for move in moves]
    return words


def mix64(value: int) -> int:
    value = (value ^ SEED) & MASK64
    value = ((value ^ (value >> 30)) * MULTIPLIER_1) & MASK64
    value = ((value ^ (value >> 27)) * MULTIPLIER_2) & MASK64
    return value ^ (value >> 31)


def projection(vector: dict[int, int], cache: dict[int, tuple[int, int]]) -> list[int]:
    output = [0] * WIDTH
    for code, coefficient in vector.items():
        entry = cache.get(code)
        if entry is None:
            value = mix64(code)
            entry = (value & (WIDTH - 1), 1 if ((value >> 9) & 1) else -1)
            cache[code] = entry
        bucket, sign = entry
        output[bucket] += sign * coefficient
    return output


class HashRank:
    def __init__(self, prime: int) -> None:
        self.prime = prime
        self.rows: list[np.ndarray] = []
        self.pivots: dict[int, int] = {}

    @property
    def rank(self) -> int:
        return len(self.rows)

    def add(self, row: list[int]) -> bool:
        value = np.asarray([entry % self.prime for entry in row], dtype=np.int64)
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                return False
            lead = int(nonzero[0])
            old = self.pivots.get(lead)
            if old is None:
                value = value * pow(int(value[lead]), self.prime - 2, self.prime) % self.prime
                self.pivots[lead] = len(self.rows)
                self.rows.append(value.astype(np.int32))
                return True
            value = (value - int(value[lead]) * self.rows[old].astype(np.int64)) % self.prime


def matrix_digest(rows: list[list[int]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(",".join(str(value) for value in row).encode("ascii"))
        digest.update(b";")
    return digest.hexdigest()


def rebuild_l9_hash() -> dict[str, object]:
    generator = cached_terms(9)
    one = HashRank(P1)
    two = HashRank(P2)
    selected_two = HashRank(P2)
    cache: dict[int, tuple[int, int]] = {}
    selected_words: list[str] = []
    selected_rows: list[list[int]] = []
    progress: list[dict[str, object]] = []
    leaves = 0
    max_support = 0
    stopped = False

    def visit(word: str, vector: dict[int, int], remaining: int) -> None:
        nonlocal leaves, max_support, stopped
        if stopped:
            return
        if not remaining:
            row = projection(vector, cache)
            accepted = one.add(row)
            two.add(row)
            leaves += 1
            max_support = max(max_support, len(vector))
            if accepted:
                selected_words.append(word)
                selected_rows.append(row)
                selected_two.add(row)
            if one.rank >= WIDTH and two.rank >= WIDTH and selected_two.rank >= WIDTH:
                stopped = True
            return
        after_b = ad(generator["B"], vector, 18)
        visit("AB" + word, ad(generator["A"], after_b, 18), remaining - 1)
        visit("BB" + word, ad(generator["B"], after_b, 18), remaining - 1)

    for root in SEEDS_L3:
        if stopped:
            break
        before = one.rank
        visit(root, word_vector(root, 9), 6)
        progress.append(
            {
                "root": root,
                "rank_increment_mod_p1": one.rank - before,
                "rank_mod_p1": one.rank,
                "rank_mod_p2": two.rank,
                "leaves_visited": leaves,
            }
        )
    return {
        "stopped": stopped,
        "leaves": leaves,
        "rank_p1": one.rank,
        "rank_p2_all": two.rank,
        "rank_p2_selected": selected_two.rank,
        "selected_words": selected_words,
        "selected_rows": selected_rows,
        "max_support": max_support,
        "progress": progress,
    }


def remap(mask: int, length: int) -> int:
    result = 0
    for column in range(length - 1):
        if mask & (1 << column):
            result |= 1 << column
        if mask & (1 << (length + column)):
            result |= 1 << ((length - 1) + column)
    return result


def last(mask: int, length: int) -> bool:
    return bool(mask & ((1 << (length - 1)) | (1 << (2 * length - 1))))


def erase(vector: dict[int, int], length: int) -> dict[int, int]:
    n = 2 * length
    mask = (1 << n) - 1
    result: dict[int, int] = {}
    for code, coefficient in vector.items():
        x = code & mask
        z = code >> n
        if last(x, length) or last(z, length):
            continue
        target = remap(x, length) | (remap(z, length) << (n - 2))
        result[target] = result.get(target, 0) + coefficient
    return {code: coefficient for code, coefficient in result.items() if coefficient}


def verify_restriction_obstruction(artifact: dict[str, object]) -> None:
    old_length, new_length = 2, 3
    old_n, new_n = 4, 6
    old_top, new_top = 1, 2
    p = (1 << old_top) | (1 << (new_n + new_top))
    left = erase(ad(cached_terms(new_length)["B"], {p: 1}, new_n), new_length)
    right = ad(cached_terms(old_length)["B"], erase({p: 1}, new_length), old_n)
    expected = (1 << old_top) | (1 << (old_n + old_top))
    check("generic restriction witness", left == {expected: -1} and right == {})
    word = "ABABAB"
    projected = erase(word_vector(word, 3), 3)
    shorter = word_vector(word, 2)
    difference = {
        code: projected.get(code, 0) - shorter.get(code, 0)
        for code in set(projected) | set(shorter)
        if projected.get(code, 0) != shorter.get(code, 0)
    }
    record = artifact["restriction_obstruction"]["literal_word_counterexample"]
    digest = hashlib.sha256()
    for code, coefficient in sorted(difference.items()):
        digest.update(f"{code}:{coefficient};".encode("ascii"))
    check("literal restriction obstruction", bool(difference) and digest.hexdigest() == record["difference_digest"])


def verify_seed_tables(artifact: dict[str, object]) -> None:
    basis_l2, rank_l2 = exact_basis(2, 4)
    basis_l3, rank_l3 = exact_basis(3, 6)
    check("exact L2 seed basis", tuple(basis_l2) == SEEDS_L2 and rank_l2 == 6)
    check("exact L3 seed basis", tuple(basis_l3) == SEEDS_L3 and rank_l3 == 16)
    check("artifact exact seed records", artifact["exact_seed_spaces"]["L2_depth_le_4"]["rank_Q"] == 6 and artifact["exact_seed_spaces"]["L3_depth_le_6"]["rank_Q"] == 16)
    observed = []
    for length in range(2, 9):
        rows = [word_vector(word, length) for word in expand(SEEDS_L2, ("AB", "BB"), length - 2)]
        observed.append(ranks(rows))
        del rows
        gc.collect()
    expected_records = artifact["six_seed_family"]["table"]
    expected = [(record["ranks_mod"][str(P1)], record["ranks_mod"][str(P2)]) for record in expected_records]
    check("six-seed rank table", observed == expected, str(observed))


def alternating(roots: tuple[str, ...], pair: tuple[str, str], length: int, triple_first: bool) -> list[str]:
    words = list(roots)
    for current in range(2, length):
        triple = (current % 2 == 0) if triple_first else (current % 2 == 1)
        moves = pair if triple else ("A", "B")
        words = [move + word for word in words for move in moves]
    return words


def verify_three_letter(artifact: dict[str, object]) -> None:
    records = artifact["three_letter_families"]
    check("three-letter record count", len(records) == 6)
    for record in records:
        orientation = record["orientation"]
        roots = CLASSIC if orientation == "one_then_three" else SHORT
        triple_first = orientation == "three_then_one"
        pair = tuple(record["three_letter_pair"])
        observed = []
        exact_l3 = None
        for length in range(2, 9):
            rows = [word_vector(word, length) for word in alternating(roots, pair, length, triple_first)]
            observed.append(ranks(rows))
            if length == 3:
                q = QEchelon()
                for row in rows:
                    q.add(row)
                exact_l3 = q.rank
            del rows
            gc.collect()
        expected = [(item["ranks_mod"][str(P1)], item["ranks_mod"][str(P2)]) for item in record["table"]]
        check(f"three-letter {orientation} {pair}", observed == expected and exact_l3 == record["exact_rank_Q_at_L3"], str(observed))


def verify_p8_transition(artifact: dict[str, object]) -> None:
    words = p8_words()
    rows = [word_vector(word, 8) for word in words]
    p8_rank = ranks(rows)
    check("patched L8 rank", p8_rank == (256, 256))
    del rows
    gc.collect()
    generator = cached_terms(9)
    one = ModEchelon(P1)
    two = ModEchelon(P2)
    seen = set()
    for parent in words:
        v = word_vector(parent, 9)
        after_b = ad(generator["B"], v, 18)
        for label, row in (("AB" + parent, ad(generator["A"], after_b, 18)), ("BB" + parent, ad(generator["B"], after_b, 18))):
            if label not in seen:
                seen.add(label)
                one.add(row)
                two.add(row)
        del v, after_b
    transition = artifact["patch_transition"]
    check("P8 child modular ranks", (one.rank, two.rank) == (transition["base_ranks_mod"][str(P1)], transition["base_ranks_mod"][str(P2)]), f"{one.rank},{two.rank}")
    additions = 0
    for value in range(128):
        word = patch_word(f"{value:07b}")
        if word in seen:
            continue
        seen.add(word)
        row = word_vector(word, 9)
        one.add(row)
        two.add(row)
        additions += 1
    check("BAB repair count", additions == transition["repair_words_added"])
    check("P8 plus BAB modular ranks", (one.rank, two.rank) == (transition["combined_ranks_mod"][str(P1)], transition["combined_ranks_mod"][str(P2)]), f"{one.rank},{two.rank}")


def verify_l9_certificate(artifact: dict[str, object]) -> None:
    rebuilt = rebuild_l9_hash()
    record = artifact["l9_certificate"]
    check("L9 certificate terminates", rebuilt["stopped"] and rebuilt["leaves"] == record["leaves_visited"])
    check("L9 same-word two-prime rank", rebuilt["rank_p1"] == WIDTH and rebuilt["rank_p2_selected"] == WIDTH)
    words_digest = hashlib.sha256()
    for word in rebuilt["selected_words"]:
        words_digest.update(f"{word};".encode("ascii"))
    check("L9 selected words", rebuilt["selected_words"] == record["selected_words"] and words_digest.hexdigest() == record["selected_words_digest"])
    check("L9 hash matrix digest", matrix_digest(rebuilt["selected_rows"]) == record["hash_projection"]["selected_matrix_integer_digest"])
    check("L9 root progression", rebuilt["progress"] == record["root_progress"])
    check("L9 depth bound", max(map(len, rebuilt["selected_words"])) <= 18)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check("schema version", artifact["schema_version"] == 1)
    check("artifact status", artifact["status"] == "FINITE_L9_CERTIFICATE_ALL_L_UNRESOLVED")
    check("prime provenance", artifact["provenance"]["primes"] == [P1, P2])
    verify_restriction_obstruction(artifact)
    verify_seed_tables(artifact)
    verify_l9_certificate(artifact)
    verify_p8_transition(artifact)
    verify_three_letter(artifact)
    print("PASS: ladder nonstationary clean-room verifier")
    return 0


def on_timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"verifier exceeded {TEST_TIMEOUT_SECONDS}s")


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, on_timeout)
    signal.alarm(TEST_TIMEOUT_SECONDS)
    try:
        raise SystemExit(main())
    finally:
        signal.alarm(0)
