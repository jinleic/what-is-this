"""Independent certificate checks for experiments/e74_dla_3xl.py.

This test intentionally does not import the experiment.  It rebuilds Pauli
brackets, consecutive finite-family ranks, exact minors, the product-order
obstruction, the tiny saturated 3x2 closure, and the artifact envelope.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "dla_3xl.json"
P1 = 2_147_483_647
P2 = 2_147_483_629
PRIMES = (P1, P2)
LABELS = "IXYZ"
BITS_TO_LABEL = {(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}
Vector = dict[int, int]
EXPECTED_PROFILE_3X2 = [2, 3, 5, 7, 11, 16, 26, 38, 57, 85, 129, 182, 245, 262, 263, 263]


def fail(message: str) -> int:
    print(f"FAIL: {message}")
    return 1


def grid(length: int) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    n = 3 * length
    bonds: list[tuple[int, int]] = []
    for row in range(3):
        for col in range(length):
            site = row * length + col
            if row < 2:
                bonds.append((site, site + length))
            if col + 1 < length:
                bonds.append((site, site + 1))
    a = tuple(1 << site for site in range(n))
    b = tuple((1 << (n + left)) | (1 << (n + right)) for left, right in bonds)
    return n, a, b


def half_sign(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    forward = (((left >> n) & mask) & (right & mask)).bit_count() & 1
    backward = (((right >> n) & mask) & (left & mask)).bit_count() & 1
    if forward == backward:
        return 0
    return 1 if forward == 0 else -1


def bracket(terms: Sequence[int], vector: Vector, n: int) -> Vector:
    result: Vector = {}
    for term in terms:
        for code, coefficient in vector.items():
            sign = half_sign(term, code, n)
            if sign:
                target = term ^ code
                value = result.get(target, 0) + sign * coefficient
                if value:
                    result[target] = value
                else:
                    result.pop(target, None)
    return result


def evaluate(word: str, length: int) -> Vector:
    n, a, b = grid(length)
    generators = {"A": a, "B": b}
    vector = {term: 1 for term in generators[word[-1]]}
    for letter in reversed(word[:-1]):
        vector = bracket(generators[letter], vector, n)
    return vector


def family_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def rank_mod(rows: Sequence[Vector], prime: int) -> int:
    pivots: dict[int, Vector] = {}
    for raw in rows:
        row = {column: coefficient % prime for column, coefficient in raw.items() if coefficient % prime}
        while row:
            lead = min(row)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(row[lead], prime - 2, prime)
                row = {
                    column: coefficient * inverse % prime
                    for column, coefficient in row.items()
                    if coefficient * inverse % prime
                }
                pivots[lead] = row
                break
            factor = row[lead]
            for column, coefficient in pivot.items():
                value = (row.get(column, 0) - factor * coefficient) % prime
                if value:
                    row[column] = value
                else:
                    row.pop(column, None)
    return len(pivots)


def determinant_bareiss(matrix: Sequence[Sequence[int]]) -> int:
    data = [list(map(int, row)) for row in matrix]
    size = len(data)
    sign = 1
    denominator = 1
    for column in range(size - 1):
        if data[column][column] == 0:
            swap = next((row for row in range(column + 1, size) if data[row][column]), None)
            if swap is None:
                return 0
            data[column], data[swap] = data[swap], data[column]
            sign = -sign
        pivot = data[column][column]
        for row in range(column + 1, size):
            for other in range(column + 1, size):
                data[row][other] = (
                    data[row][other] * pivot - data[row][column] * data[column][other]
                ) // denominator
            data[row][column] = 0
        denominator = pivot
    return sign * data[-1][-1]


def local_key(
    code: int,
    length: int,
    rung_reverse: bool,
    row_reverse: bool,
    ascending_order: str,
) -> tuple[int, ...]:
    n = 3 * length
    mask = (1 << n) - 1
    a, b = code & mask, (code >> n) & mask
    rank = {letter: index for index, letter in enumerate(ascending_order)}
    columns: Iterable[int] = range(length - 1, -1, -1) if rung_reverse else range(length)
    rows: Iterable[int] = range(2, -1, -1) if row_reverse else range(3)
    return tuple(
        rank[BITS_TO_LABEL[((a >> (row * length + col)) & 1, (b >> (row * length + col)) & 1)]]
        for col in columns
        for row in rows
    )


def maximum_product_leaders(vectors: Sequence[Vector], length: int) -> tuple[int, dict[int, int]]:
    histogram: dict[int, int] = {}
    for rung_reverse, row_reverse, order in itertools.product(
        (False, True), (False, True), itertools.permutations(LABELS)
    ):
        ascending = "".join(order)
        leaders = {
            max(
                vector,
                key=lambda code: local_key(code, length, rung_reverse, row_reverse, ascending),
            )
            for vector in vectors
        }
        count = len(leaders)
        histogram[count] = histogram.get(count, 0) + 1
    return max(histogram), histogram


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def symmetry_group(length: int) -> tuple[tuple[int, ...], ...]:
    n = 3 * length
    row_flip = tuple((2 - row) * length + col for row in range(3) for col in range(length))
    col_flip = tuple(row * length + length - 1 - col for row in range(3) for col in range(length))
    group = [tuple(range(n))]
    seen = set(group)
    head = 0
    while head < len(group):
        current = group[head]
        head += 1
        for generator in (row_flip, col_flip):
            product = compose(generator, current)
            if product not in seen:
                seen.add(product)
                group.append(product)
    return tuple(group)


def permute(code: int, permutation: Sequence[int], n: int) -> int:
    mask = (1 << n) - 1
    a, b = code & mask, (code >> n) & mask
    new_a = 0
    new_b = 0
    for old, new in enumerate(permutation):
        new_a |= ((a >> old) & 1) << new
        new_b |= ((b >> old) & 1) << new
    return new_a | (new_b << n)


def orbit_profile_3x2(prime: int) -> list[int]:
    """Independent tiny orbit closure and modular BFS for the complete 3x2 DLA."""

    length = 2
    n, a_terms, b_terms = grid(length)
    generators = a_terms + b_terms
    group = symmetry_group(length)
    cache: dict[int, int] = {}

    def canonical(code: int) -> int:
        if code not in cache:
            cache[code] = min(permute(code, permutation, n) for permutation in group)
        return cache[code]

    support = {canonical(term) for term in generators}
    queue = list(support)
    head = 0
    while head < len(queue):
        code = queue[head]
        head += 1
        for term in generators:
            if half_sign(term, code, n):
                target = canonical(term ^ code)
                if target not in support:
                    support.add(target)
                    queue.append(target)
    if len(support) != 304:
        raise AssertionError(f"3x2 support changed: {len(support)}")

    def seed(terms: Sequence[int]) -> Vector:
        result: Vector = {}
        for term in terms:
            key = canonical(term)
            result[key] = result.get(key, 0) + 1
        # Per-string orbit coordinates: symmetric copies have coefficient one.
        return {key: 1 for key in result}

    def action(vector: Vector, terms: Sequence[int]) -> Vector:
        # Orbit coordinates are coefficients per Pauli string.  A canonical
        # source represents its full symmetry orbit, so sum every preimage
        # contribution landing in each target orbit.
        result: Vector = {}
        for source, source_coefficient in vector.items():
            targets = {
                canonical(source ^ term)
                for term in terms
                if half_sign(term, source, n)
            }
            for target in targets:
                coefficient = 0
                for term in terms:
                    preimage = target ^ term
                    sign = half_sign(term, preimage, n)
                    if sign and canonical(preimage) == source:
                        coefficient += sign
                value = (result.get(target, 0) + coefficient * source_coefficient) % prime
                if value:
                    result[target] = value
                else:
                    result.pop(target, None)
        return result

    pivots: dict[int, Vector] = {}

    def add(raw: Vector) -> Vector | None:
        vector = {column: value % prime for column, value in raw.items() if value % prime}
        while vector:
            lead = min(vector)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(vector[lead], prime - 2, prime)
                vector = {
                    column: value * inverse % prime
                    for column, value in vector.items()
                    if value * inverse % prime
                }
                pivots[lead] = vector
                return vector
            factor = vector[lead]
            for column, value in pivot.items():
                reduced = (vector.get(column, 0) - factor * value) % prime
                if reduced:
                    vector[column] = reduced
                else:
                    vector.pop(column, None)
        return None

    frontier = [row for row in (add(seed(a_terms)), add(seed(b_terms))) if row is not None]
    dimensions = [len(pivots)]
    while frontier:
        next_frontier: list[Vector] = []
        for row in frontier:
            for terms in (a_terms, b_terms):
                accepted = add(action(row, terms))
                if accepted is not None:
                    next_frontier.append(accepted)
        dimensions.append(len(pivots))
        frontier = next_frontier
    return dimensions


def main() -> int:
    if not ARTIFACT.exists():
        return fail("missing results/algebra_growth/dla_3xl.json")
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if set(artifact) != {"provenance", "data", "checks"}:
        return fail("artifact envelope must contain provenance/data/checks")
    if not artifact["checks"] or not all(check.get("passed") is True for check in artifact["checks"]):
        return fail("artifact contains a failed or malformed check")
    data = artifact["data"]
    if data["status"] != "UNRESOLVED_ALL_L_WITH_NEW_FINITE_CERTIFICATES":
        return fail("all-L status was promoted or changed")
    if "lower bound" not in data["modular_scope"] or "cross-check only" not in data["modular_scope"]:
        return fail("modular-rank caveat is missing")

    # Independently rebuild the smallest complete closure at a different prime.
    if orbit_profile_3x2(P2) != EXPECTED_PROFILE_3X2:
        return fail("independent 3x2 orbit profile changed")

    exact_source = json.loads(
        (ROOT / "results" / "algebra_structure" / "char0_levi.json").read_text(encoding="utf-8")
    )["data"]
    exact_record = data["exact_three_by_two"]
    if not (
        exact_record["exact_dimension_Q"] == 263
        and exact_source["dimension_Q"] == 263
        and exact_source["faithful_matrix_dimension_Q"] == 263
        and exact_source["basis"]["faithful_word_count"] == 263
    ):
        return fail("paired exact-Q 3x2 dimension certificate changed")

    stored_family = {int(record["length"]): record for record in data["binary_family"]["records"]}
    for length in (4, 5, 6):
        words = family_words(length)
        vectors = [evaluate(word, length) for _, word in words]
        target = 1 << length
        if any(rank_mod(vectors, prime) != target for prime in PRIMES):
            return fail(f"independent consecutive family rank failed at L={length}")
        stored = stored_family[length]
        stored_rows = stored["rows"]
        if [
            (row["label"], row["subset_columns"], row["word"], row["support_size"])
            for row in stored_rows
        ] != [
            (
                label,
                [column + 1 for column, bit in enumerate(label) if bit == "1"],
                word,
                len(vector),
            )
            for (label, word), vector in zip(words, vectors)
        ]:
            return fail(f"stored subset-indexed words/support sizes changed at L={length}")

    # Recompute exact integer minors at consecutive sizes 4 and 5 from stored columns.
    for length in (4, 5):
        vectors = [evaluate(word, length) for _, word in family_words(length)]
        minor = stored_family[length]["integer_minor"]
        columns = [int(column) for column in minor["pauli_columns"]]
        determinant = determinant_bareiss(
            [[vector.get(column, 0) for column in columns] for vector in vectors]
        )
        if str(determinant) != minor["determinant_Z"] or determinant == 0:
            return fail(f"exact integer minor changed at L={length}")
        if any(determinant % prime == 0 for prime in PRIMES):
            return fail(f"stored minor vanished at a recorded prime at L={length}")

    obstruction = {int(record["length"]): record for record in data["product_lex_obstruction"]["records"]}
    for length, expected in ((4, 13), (5, 18)):
        vectors = [evaluate(word, length) for _, word in family_words(length)]
        maximum, histogram = maximum_product_leaders(vectors, length)
        if maximum != expected or maximum >= 1 << length:
            return fail(f"96-order obstruction changed at L={length}")
        stored_histogram = {int(key): int(value) for key, value in obstruction[length]["histogram"].items()}
        if histogram != stored_histogram or sum(histogram.values()) != 96:
            return fail(f"product-order histogram changed at L={length}")

    profile_map = {tuple(record["shape"]): record for record in data["depth_profiles"]}
    if profile_map[(3, 2)].get("exact_dimension_Q") != 263:
        return fail("3x2 exact-Q equality was dropped")
    if profile_map[(3, 3)]["lower_bound_Q"] != 8034 or not profile_map[(3, 3)]["saturated"]:
        return fail("3x3 scope/value changed")
    if profile_map[(3, 4)]["dimensions"][-2:] != [741, 1242]:
        return fail("3x4 wave-8 depth extension changed")
    if profile_map[(3, 4)]["modular_dimension"] is not None or profile_map[(3, 4)]["saturated"]:
        return fail("3x4 incomplete profile was promoted to a dimension equality")

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
