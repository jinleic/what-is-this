"""Standalone exact checks for the wave-7 two-leg-ladder leading-term search."""

from __future__ import annotations

import hashlib
import itertools
import json
import signal
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_growth" / "ladder_alll.json"
TIMEOUT_SECONDS = 1_800
P1 = 2_147_483_647
P2 = 2_147_483_629
TO_BITS = {"I": (0, 0), "X": (1, 0), "Z": (0, 1), "Y": (1, 1)}
FROM_BITS = {value: key for key, value in TO_BITS.items()}


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    lr = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    rl = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if lr == rl:
        return 0
    return 1 if lr == 0 else -1


def bracket(left: dict[int, int], right: dict[int, int], n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    for lhs, lc in left.items():
        for rhs, rc in right.items():
            sign = sign_half(lhs, rhs, n)
            if sign:
                target = lhs ^ rhs
                out[target] = out.get(target, 0) + sign * lc * rc
                if not out[target]:
                    del out[target]
    return out


def generators(length: int) -> dict[str, dict[int, int]]:
    n = 2 * length
    a = {1 << site: 1 for site in range(n)}
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    b = {(1 << (n + lhs)) | (1 << (n + rhs)): 1 for lhs, rhs in bonds}
    return {"A": a, "B": b}


def recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def vector_for_word(word: str, length: int) -> dict[int, int]:
    gens = generators(length)
    value = dict(gens[word[-1]])
    for letter in reversed(word[:-1]):
        value = bracket(gens[letter], value, 2 * length)
    return value


def recursive_vectors(length: int) -> list[tuple[str, str, dict[int, int]]]:
    return [(label, word, vector_for_word(word, length)) for label, word in recursive_words(length)]


def site_labels(pauli: int, length: int) -> tuple[str, ...]:
    n = 2 * length
    mask = (1 << n) - 1
    a, b = pauli & mask, (pauli >> n) & mask
    return tuple(FROM_BITS[((a >> site) & 1, (b >> site) & 1)] for site in range(n))


def rung_string(pauli: int, length: int) -> str:
    labels = site_labels(pauli, length)
    return " ".join(labels[col] + labels[length + col] for col in range(length))


def key_for(
    pauli: int,
    length: int,
    rung_reverse: bool = True,
    row_reverse: bool = False,
    alphabet: str = "IZYX",
) -> tuple[int, ...]:
    labels = site_labels(pauli, length)
    rank = {letter: index for index, letter in enumerate(alphabet)}
    columns = range(length - 1, -1, -1) if rung_reverse else range(length)
    key = []
    for col in columns:
        pair = (labels[col], labels[length + col])
        if row_reverse:
            pair = pair[::-1]
        key.extend(rank[letter] for letter in pair)
    return tuple(key)


def digest(vector: dict[int, int]) -> str:
    value = hashlib.sha256()
    for pauli, coefficient in sorted(vector.items()):
        value.update(f"{pauli}:{coefficient};".encode("ascii"))
    return value.hexdigest()


def rank_mod(rows: list[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for row in rows:
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        while value:
            lead = min(value)
            if lead not in pivots:
                inverse = pow(value[lead], prime - 2, prime)
                value = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                pivots[lead] = value
                break
            factor = value[lead]
            for key, coefficient in pivots[lead].items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def rational_combinations(
    rows: list[dict[int, int]],
    length: int,
    order,
) -> list[dict[int, Fraction]]:
    pivots = {}
    combinations = []
    for index, row in enumerate(rows):
        value = {pauli: Fraction(coefficient) for pauli, coefficient in row.items()}
        combination = {index: Fraction(1)}
        while value:
            lead = max(value, key=lambda pauli: order(pauli, length))
            if lead not in pivots:
                inverse = 1 / value[lead]
                value = {pauli: coefficient * inverse for pauli, coefficient in value.items()}
                combination = {row_index: coefficient * inverse for row_index, coefficient in combination.items()}
                pivots[lead] = (value, combination)
                combinations.append(combination)
                break
            basis, basis_combination = pivots[lead]
            factor = value[lead]
            for pauli, coefficient in basis.items():
                reduced = value.get(pauli, 0) - factor * coefficient
                if reduced:
                    value[pauli] = reduced
                else:
                    value.pop(pauli, None)
            for row_index, coefficient in basis_combination.items():
                reduced = combination.get(row_index, 0) - factor * coefficient
                if reduced:
                    combination[row_index] = reduced
                else:
                    combination.pop(row_index, None)
    return combinations


def combine(combination, rows):
    result = {}
    for row_index, factor in combination.items():
        for pauli, coefficient in rows[row_index].items():
            total = result.get(pauli, 0) + factor * coefficient
            if total:
                result[pauli] = total
            else:
                result.pop(pauli, None)
    return result


def lifted_step(length: int, order) -> tuple[list[tuple[int, int, dict[int, Fraction], int]], dict]:
    parents = [row[2] for row in recursive_vectors(length)]
    children = [row[2] for row in recursive_vectors(length + 1)]
    combinations = rational_combinations(parents, length, order)
    lifts = []
    for parent_index, combination in enumerate(combinations):
        for bit in (0, 1):
            value = combine({2 * row_index + bit: factor for row_index, factor in combination.items()}, children)
            leader = max(value, key=lambda pauli: order(pauli, length + 1))
            lifts.append((parent_index, bit, value, leader))
    collision = {}
    seen = {}
    for parent_index, bit, value, leader in lifts:
        if leader in seen and not collision:
            old_parent, old_bit, old_value = seen[leader]
            supports = sorted(set(old_value) | set(value))
            for left_index, left in enumerate(supports):
                for right in supports[left_index + 1 :]:
                    determinant = old_value.get(left, 0) * value.get(right, 0) - old_value.get(right, 0) * value.get(left, 0)
                    if determinant:
                        collision = {
                            "parents": (old_parent, parent_index),
                            "bits": (old_bit, bit),
                            "leader": leader,
                            "coefficients": (old_value[leader], value[leader]),
                            "support_codes": (left, right),
                            "determinant": determinant,
                        }
                        break
                if collision:
                    break
        else:
            seen[leader] = (parent_index, bit, value)
    return lifts, collision


def check_two_consecutive_lengths() -> str:
    for length in (6, 7):
        rows = recursive_vectors(length)
        assert len(rows) == 1 << length
        assert rank_mod([row[2] for row in rows], P1) == 1 << length
        assert rank_mod([row[2] for row in rows], P2) == 1 << length
        assert all(len(word) <= 2 * length for _, word, _ in rows)
    return "recomputed two-prime ranks 64 and 128 from exact L=6,7 bracket vectors"


def check_leaders_against_artifact() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    stored_by_length = {row["length"]: row for row in artifact["data"]["recursive_family"]}
    for length in (6, 7):
        actual = recursive_vectors(length)
        stored = stored_by_length[length]["rows"]
        assert len(actual) == len(stored)
        for (label, word, vector), record in zip(actual, stored):
            leader = max(vector, key=lambda pauli: key_for(pauli, length))
            assert (label, word) == (record["label"], record["word"])
            assert leader == record["leader_code"]
            assert rung_string(leader, length) == record["leader_pauli"]
            assert vector[leader] == record["leader_coefficient"]
            assert digest(vector) == record["sha256"]
    return "recomputed every fixed-order leader and digest at consecutive L=6,7"


def check_triangular_base_failure() -> str:
    order = lambda pauli, length: key_for(pauli, length)
    lifts, collision = lifted_step(2, order)
    assert len(lifts) == 8
    assert len({row[3] for row in lifts}) == 7
    assert collision == {
        "parents": (1, 3),
        "bits": (0, 1),
        "leader": 2212,
        "coefficients": (Fraction(-2), Fraction(4)),
        "support_codes": (323, 649),
        "determinant": Fraction(6),
    }
    assert rung_string(collision["leader"], 3) == "II ZI XY"
    return "exact L=2->3 lift collision has leader II ZI XY and determinant 6"


def check_all_product_orders() -> str:
    histogram = {}
    for rung_reverse, row_reverse, permutation in itertools.product(
        (False, True), (False, True), itertools.permutations("IXYZ")
    ):
        alphabet = "".join(permutation)
        order = lambda pauli, length, rr=rung_reverse, row=row_reverse, alpha=alphabet: key_for(
            pauli, length, rr, row, alpha
        )
        lifts, _ = lifted_step(2, order)
        count = len({lift[3] for lift in lifts})
        histogram[count] = histogram.get(count, 0) + 1
    assert histogram == {5: 16, 6: 36, 7: 44}
    return "all 96 product-lex orders fail: leader histogram {5:16, 6:36, 7:44}"


def check_artifact() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    assert set(artifact) == {"provenance", "data", "checks"}
    assert artifact["provenance"]["script"] == "experiments/e53_ladder_alll.py"
    assert artifact["provenance"]["interpreter"] == ".venv/bin/python"
    assert artifact["data"]["status"] == "UNRESOLVED"
    assert artifact["data"]["target"].startswith("[UNRESOLVED]")
    assert artifact["checks"] and all(check["passed"] for check in artifact["checks"])
    return f"validated envelope and {len(artifact['checks'])} passing stored checks"


def run_check(name, function) -> bool:
    detail = function()
    print(f"PASS {name}: {detail}")
    return True


def timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s")


def main() -> int:
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        checks = [
            ("two consecutive finite ranks", check_two_consecutive_lengths),
            ("leading families", check_leaders_against_artifact),
            ("triangular base failure", check_triangular_base_failure),
            ("product-order exhaustion", check_all_product_orders),
            ("artifact", check_artifact),
        ]
        passed = [run_check(name, function) for name, function in checks]
        print("PASS" if all(passed) else "FAIL")
        return 0 if all(passed) else 1
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
