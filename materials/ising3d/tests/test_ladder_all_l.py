"""Standalone exact checks for the wave-9 ladder word-family certificate."""

from __future__ import annotations

import gc
import hashlib
import itertools
import json
import signal
import time
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_growth" / "ladder_all_l.json"
TIMEOUT_SECONDS = 3_600
P1 = 2_147_483_647
P2 = 2_147_483_629
MOVES = ("AA", "AB", "BA", "BB")
MOVE_PAIRS = tuple(itertools.combinations(MOVES, 2))
BASE_WORDS = tuple(
    "".join(letters)
    for depth in range(1, 5)
    for letters in itertools.product("AB", repeat=depth)
)
EXPECTED_BASE_HISTOGRAM = {0: 1001, 1: 8269, 2: 12708, 3: 4991, 4: 436}
EXPECTED_STAGE_HISTOGRAMS = {
    3: {3: 8, 4: 88, 5: 520, 6: 920, 7: 672, 8: 408},
    4: {14: 48, 15: 32, 16: 328},
    5: {30: 64, 31: 32, 32: 232},
    6: {61: 40, 63: 48, 64: 144},
    7: {128: 144},
}
EXPECTED_BOUNDARY_HISTOGRAMS = {
    2: {4: 436},
    3: {8: 408},
    4: {16: 328},
    5: {32: 232},
    6: {64: 144},
    7: {128: 144},
    8: {248: 144},
}
EXPECTED_REMOVED_LABELS = (
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
EXPECTED_SUPPLEMENT_WORDS = (
    "ABABABABABABBAB",
    "ABABABABABBBBAB",
    "ABABABABBBABBAB",
    "ABABABABBBBBBAB",
    "ABABABBBABABBAB",
    "ABABBBABABABBAB",
    "ABABBBABABBBBAB",
    "ABABBBABBBABBAB",
)
EXPECTED_RELATION_SECTOR_COUNTS = {(0, 0): 3, (0, 1): 0, (1, 0): 3, (1, 1): 2}
EXPECTED_SECTOR_RANKS = {(0, 0): 61, (0, 1): 64, (1, 0): 61, (1, 1): 62}

_BOUNDARY_BY_WORD: dict[str, dict[int, int]] | None = None
_NATURAL_RECORD: dict[str, object] | None = None


def generator_terms(length: int) -> dict[str, tuple[int, ...]]:
    if length < 2:
        raise ValueError("the open ladder has at least two rungs")
    n = 2 * length
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return {
        "A": tuple(1 << site for site in range(n)),
        "B": tuple((1 << (n + left)) | (1 << (n + right)) for left, right in bonds),
    }


@lru_cache(maxsize=None)
def cached_terms(length: int) -> dict[str, tuple[int, ...]]:
    return generator_terms(length)


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    left_right = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    right_left = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if left_right == right_left:
        return 0
    return 1 if left_right == 0 else -1


def adjoint(terms: tuple[int, ...], vector: dict[int, int], n: int) -> dict[int, int]:
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
    terms = cached_terms(length)
    value = {term: 1 for term in terms[word[-1]]}
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


@lru_cache(maxsize=None)
def cached_vector(length: int, word: str) -> dict[int, int]:
    return vector_for_word(word, length)


def rank_mod(rows: list[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for row in rows:
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(value[lead], prime - 2, prime)
                pivots[lead] = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                break
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def exact_rank(rows: list[dict[int, int]]) -> int:
    """Sparse rational echelon rank; used only when a modular full minor is absent."""
    pivots: dict[int, dict[int, Fraction]] = {}
    for row in rows:
        value = {key: Fraction(coefficient) for key, coefficient in row.items() if coefficient}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = 1 / value[lead]
                pivots[lead] = {key: coefficient * inverse for key, coefficient in value.items()}
                break
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = value.get(key, Fraction(0)) - factor * coefficient
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def certified_rank(rows: list[dict[int, int]], target: int) -> int:
    modular = rank_mod(rows, P1)
    if modular == target:
        return target
    return exact_rank(rows)


def family_words(base: tuple[str, ...], moves: tuple[str, str], length: int) -> list[str]:
    words = list(base)
    for _ in range(2, length):
        words = [move + word for word in words for move in moves]
    return words


def family_rows(base: tuple[str, ...], moves: tuple[str, str], length: int) -> list[dict[int, int]]:
    words = family_words(base, moves, length)
    assert len(words) == 1 << length
    assert max(map(len, words)) <= 2 * length
    return [cached_vector(length, word) for word in words]


def project_last_rungs(vector: dict[int, int], length: int, width: int) -> dict[int, int]:
    n = 2 * length
    row_mask = ((1 << width) - 1) << (length - width)
    site_mask = row_mask | (row_mask << length)
    pauli_mask = site_mask | (site_mask << n)
    return {code: coefficient for code, coefficient in vector.items() if not code & ~pauli_mask}


def vector_digest(vector: dict[int, int]) -> str:
    digest = hashlib.sha256()
    for code, coefficient in sorted(vector.items()):
        digest.update(f"{code}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def words_digest(words: list[str]) -> str:
    digest = hashlib.sha256()
    for word in words:
        digest.update(f"{word};".encode("ascii"))
    return digest.hexdigest()


def original_recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def original_word_for_label(label: str) -> str:
    return dict(original_recursive_words(8))[label]


def supplemented_word(extension_label: str) -> str:
    word = "BAB"
    for bit in extension_label:
        word = ("AB" if bit == "0" else "BB") + word
    return word


def run_natural_search() -> tuple[dict[str, object], dict[str, dict[int, int]]]:
    base_histogram: Counter[int] = Counter()
    base_sets: list[tuple[str, ...]] = []
    for base in itertools.combinations(BASE_WORDS, 4):
        rank = certified_rank([cached_vector(2, word) for word in base], 4)
        base_histogram[rank] += 1
        if rank == 4:
            base_sets.append(base)
    assert dict(base_histogram) == EXPECTED_BASE_HISTOGRAM
    assert len(base_sets) == 436

    boundary_histograms: dict[int, Counter[int]] = {2: Counter()}
    for base in base_sets:
        rows = [
            project_last_rungs(cached_vector(2, word), 2, 2)
            for word in base
        ]
        boundary_histograms[2][certified_rank(rows, 4)] += 1
    assert dict(boundary_histograms[2]) == EXPECTED_BOUNDARY_HISTOGRAMS[2]

    stage_records = []
    candidates = [(moves, base) for moves in MOVE_PAIRS for base in base_sets]
    for length, target in ((3, 8), (4, 16), (5, 32), (6, 64), (7, 128)):
        histogram: Counter[int] = Counter()
        survivors: list[tuple[tuple[str, str], tuple[str, ...]]] = []
        for moves, base in candidates:
            rank = certified_rank(family_rows(base, moves, length), target)
            histogram[rank] += 1
            if rank == target:
                survivors.append((moves, base))
        assert dict(histogram) == EXPECTED_STAGE_HISTOGRAMS[length]
        if length == 5:
            assert {moves for moves, _ in survivors} == {("AB", "BB")}
        if length == 7:
            moves, base = survivors[0]
            assert rank_mod(family_rows(base, moves, length), P2) == target

        boundary_histogram: Counter[int] = Counter()
        for moves, base in survivors:
            rows = [
                project_last_rungs(vector, length, min(3, length))
                for vector in family_rows(base, moves, length)
            ]
            boundary_histogram[certified_rank(rows, target)] += 1
        assert dict(boundary_histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[length]
        boundary_histograms[length] = boundary_histogram
        stage_records.append(
            {
                "length": length,
                "rank_histogram": dict(sorted(histogram.items())),
                "survivor_count": len(survivors),
            }
        )
        candidates = survivors

    seed_words = sorted({word for _, base in candidates for word in base})
    assert len(seed_words) == 15
    boundary_by_word: dict[str, dict[int, int]] = {}
    for seed in seed_words:
        for word in family_words((seed,), ("AB", "BB"), 8):
            boundary_by_word[word] = project_last_rungs(vector_for_word(word, 8), 8, 3)

    boundary_l8_histogram: Counter[int] = Counter()
    survivor_bases = []
    for moves, base in candidates:
        assert moves == ("AB", "BB")
        rows = [boundary_by_word[word] for word in family_words(base, moves, 8)]
        boundary_l8_histogram[exact_rank(rows)] += 1
        survivor_bases.append(base)
    assert dict(boundary_l8_histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[8]
    representative_rows = [
        boundary_by_word[word] for word in family_words(candidates[0][1], candidates[0][0], 8)
    ]
    modular_l8_representative = {
        str(prime): rank_mod(representative_rows, prime) for prime in (P1, P2)
    }
    assert all(rank == 248 for rank in modular_l8_representative.values())
    for length, histogram in boundary_histograms.items():
        assert dict(histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[length]

    record = {
        "base_rank_histogram": dict(sorted(base_histogram.items())),
        "stage_records": stage_records,
        "boundary_rank_histograms": {
            str(length): dict(sorted(histogram.items()))
            for length, histogram in sorted(boundary_histograms.items())
        },
        "boundary_rank_histogram_L8_exact": dict(sorted(boundary_l8_histogram.items())),
        "boundary_rank_L8_modular_representative": modular_l8_representative,
        "survivor_bases": [list(base) for base in survivor_bases],
        "survivor_base_digest": words_digest(["/".join(base) for base in survivor_bases]),
        "seed_words_L8": seed_words,
    }
    return record, boundary_by_word


def check_artifact_envelope() -> str:
    assert RESULT.is_file(), f"missing result artifact: {RESULT}"
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    assert set(artifact) == {"checks", "data", "provenance"}
    assert artifact["provenance"]["script"] == "experiments/e82_ladder_all_l.py"
    assert artifact["provenance"]["interpreter"] == ".venv/bin/python"
    assert artifact["data"]["status"] == "UNRESOLVED"
    assert artifact["data"]["finite_claim"].startswith("[COMPUTATION]")
    assert artifact["data"]["target"].startswith("[UNRESOLVED]")
    assert artifact["checks"] and all(check["passed"] for check in artifact["checks"])
    return f"validated envelope and {len(artifact['checks'])} passing stored checks"


def check_natural_class_obstruction() -> str:
    global _BOUNDARY_BY_WORD, _NATURAL_RECORD
    record, boundary_by_word = run_natural_search()
    _BOUNDARY_BY_WORD = boundary_by_word
    _NATURAL_RECORD = record
    cached_vector.cache_clear()
    gc.collect()
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    stored = artifact["data"]["natural_class_search"]
    assert stored == json.loads(json.dumps(record, sort_keys=True))
    return "exactly exhausted 436 depth-compatible bases and six two-letter move pairs; all 144 L=7 survivors have exact L=8 boundary rank 248"


def check_relation_structure() -> str:
    assert _BOUNDARY_BY_WORD is not None
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    analysis = artifact["data"]["relation_structure"]
    relations = analysis["exact_relations"]
    assert len(relations) == 8
    original = original_recursive_words(8)
    assert tuple(relation["pivot_label"] for relation in relations) == EXPECTED_REMOVED_LABELS
    original_rows = {label: _BOUNDARY_BY_WORD[word] for label, word in original}
    sector_rows: dict[tuple[int, int], list[dict[int, int]]] = {}
    for label, word in original:
        sector = (word.count("A") & 1, word.count("B") & 1)
        sector_rows.setdefault(sector, []).append(original_rows[label])
    assert {sector: len(rows) for sector, rows in sector_rows.items()} == {
        sector: 64 for sector in EXPECTED_RELATION_SECTOR_COUNTS
    }
    for sector, rows in sector_rows.items():
        assert rank_mod(rows, P1) == EXPECTED_SECTOR_RANKS[sector]
        assert rank_mod(rows, P2) == EXPECTED_SECTOR_RANKS[sector]

    relation_sector_counts: Counter[tuple[int, int]] = Counter()
    relation_rows_by_sector = {sector: [] for sector in EXPECTED_RELATION_SECTOR_COUNTS}
    for relation in relations:
        terms = relation["terms"]
        assert relation["term_count"] == len(terms)
        assert len({term["label"] for term in terms}) == len(terms)
        coefficients = [int(term["coefficient"]) for term in terms]
        common_factor = 0
        for coefficient in coefficients:
            common_factor = gcd(common_factor, abs(coefficient))
        assert common_factor == relation["coefficient_gcd"] == 1
        sectors = {
            (
                original_word_for_label(term["label"]).count("A") & 1,
                original_word_for_label(term["label"]).count("B") & 1,
            )
            for term in terms
        }
        assert len(sectors) == 1
        sector = sectors.pop()
        relation_sector_counts[sector] += 1
        relation_rows_by_sector[sector].append(
            {int(term["label"], 2): int(term["coefficient"]) for term in terms}
        )
        residual: dict[int, int] = {}
        for term in terms:
            vector = vector_for_word(original_word_for_label(term["label"]), 8)
            coefficient = int(term["coefficient"])
            for code, value in vector.items():
                total = residual.get(code, 0) + coefficient * value
                if total:
                    residual[code] = total
                else:
                    residual.pop(code, None)
        assert not residual
        assert relation["exact_residual_support_size"] == 0
    assert {
        sector: relation_sector_counts[sector] for sector in EXPECTED_RELATION_SECTOR_COUNTS
    } == EXPECTED_RELATION_SECTOR_COUNTS
    for sector, rows in relation_rows_by_sector.items():
        assert exact_rank(rows) == EXPECTED_RELATION_SECTOR_COUNTS[sector]

    terms = generator_terms(8)
    n = 16
    row_swap = tuple(list(range(8, 16)) + list(range(8)))
    rung_reverse = tuple(list(range(7, -1, -1)) + list(range(15, 7, -1)))
    for permutation in (row_swap, rung_reverse):
        for name, generator in terms.items():
            assert {permute_pauli(code, permutation, n) for code in generator} == set(generator), name
    assert all(
        (((left % 8) + (left // 8)) & 1) != (((right % 8) + (right // 8)) & 1)
        for left, right in ladder_bonds(8)
    )
    assert analysis["relation_sector_counts"] == {"00": 3, "01": 0, "10": 3, "11": 2}
    return "all eight full-Pauli primitive relations vanish exactly inside character sectors (00:3, 10:3, 11:2); row and rung reflection fix both generators"


def ladder_bonds(length: int) -> tuple[tuple[int, int], ...]:
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return tuple(bonds)


def permute_pauli(code: int, permutation: tuple[int, ...], n: int) -> int:
    mask = (1 << n) - 1
    x, z = code & mask, (code >> n) & mask
    mapped_x = 0
    mapped_z = 0
    for old, new in enumerate(permutation):
        mapped_x |= ((x >> old) & 1) << new
        mapped_z |= ((z >> old) & 1) << new
    return mapped_x | (mapped_z << n)


def check_patched_family() -> str:
    assert _BOUNDARY_BY_WORD is not None
    original = original_recursive_words(8)
    kept = [(label, word) for label, word in original if label not in EXPECTED_REMOVED_LABELS]
    assert len(kept) == 248
    supplements = [supplemented_word(label) for label in PATCH_EXTENSION_LABELS]
    assert tuple(supplements) == EXPECTED_SUPPLEMENT_WORDS
    assert all(len(word) == 15 for word in supplements)
    words = [word for _, word in kept] + supplements
    assert len(words) == len(set(words)) == 256
    rows = [_BOUNDARY_BY_WORD[word] for word in words]
    assert rank_mod(rows, P1) == 256
    assert rank_mod(rows, P2) == 256
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    stored = artifact["data"]["patched_family"]
    assert stored["removed_original_labels"] == list(EXPECTED_REMOVED_LABELS)
    assert [entry["extension_label"] for entry in stored["supplements"]] == list(PATCH_EXTENSION_LABELS)
    assert [entry["word"] for entry in stored["supplements"]] == supplements
    assert stored["word_sha256"] == words_digest(words)
    assert stored["boundary_row_digests_sha256"] == [vector_digest(row) for row in rows]
    assert stored["boundary_rank_mod"] == {str(P1): 256, str(P2): 256}
    return "rebuilt the 248+8 depth-<=16 patch and found boundary rank 256 modulo both recorded primes"


def run_check(name: str, function) -> bool:
    detail = function()
    print(f"PASS {name}: {detail}")
    return True


class BudgetExceeded(TimeoutError):
    pass


_CPU_START = time.process_time()


def _check_budget() -> None:
    """Per-call compute budget, in the test process's OWN CPU time.

    Elapsed wall clock measure the unrelated load on this machine; measured
    this way, a heavy co-resident job cannot change the verdict."""
    if time.process_time() - _CPU_START > TIMEOUT_SECONDS:
        raise BudgetExceeded(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s process CPU time")


def timeout(_signum: int, _frame: object) -> None:
    raise BudgetExceeded(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s CPU time")


def main() -> int:
    # SIGALRM fires in wall time under load, so it is used only as a Coarse
    # safety net at 4x the CPU budget; the decisive gate is process_time().
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(4 * TIMEOUT_SECONDS)
    try:
        checks = [
            ("artifact envelope", check_artifact_envelope),
            ("natural-class boundary obstruction", check_natural_class_obstruction),
            ("relation sector structure", check_relation_structure),
            ("patched L=8 family", check_patched_family),
        ]
        passed = []
        for name, function in checks:
            _check_budget()
            passed.append(run_check(name, function))
        print("PASS" if all(passed) else "FAIL")
        return 0 if all(passed) else 1
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
