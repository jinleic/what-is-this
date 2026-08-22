"""Standalone exact checks for the wave-8 ladder representation obstruction."""

from __future__ import annotations

import json
import signal
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_growth" / "ladder_representation.json"
TIMEOUT_SECONDS = 3_600
P1 = 2_147_483_647
P2 = 2_147_483_629
RECONSTRUCTION_PRIMES = (
    P1,
    P2,
    2_147_483_587,
    2_147_483_579,
    2_147_483_549,
    2_147_483_543,
    2_147_483_497,
    2_147_483_489,
)


def generator_terms(length: int) -> dict[str, tuple[int, ...]]:
    n = 2 * length
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + column, offset + column + 1) for column in range(length - 1))
    bonds.extend((column, length + column) for column in range(length))
    return {
        "A": tuple(1 << site for site in range(n)),
        "B": tuple((1 << (n + u)) | (1 << (n + v)) for u, v in bonds),
    }


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    forward = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    backward = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if forward == backward:
        return 0
    return 1 if forward == 0 else -1


def adjoint(terms: tuple[int, ...], vector: dict[int, int], n: int) -> dict[int, int]:
    result: dict[int, int] = {}
    for term in terms:
        for pauli, coefficient in vector.items():
            sign = sign_half(term, pauli, n)
            if sign:
                target = term ^ pauli
                value = result.get(target, 0) + sign * coefficient
                if value:
                    result[target] = value
                else:
                    result.pop(target, None)
    return result


def words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def vector_for_word(word: str, length: int) -> dict[int, int]:
    terms = generator_terms(length)
    value = {term: 1 for term in terms[word[-1]]}
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


def family(length: int) -> list[tuple[str, str, dict[int, int]]]:
    return [(label, word, vector_for_word(word, length)) for label, word in words(length)]


def rank_mod(rows: list[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for row in rows:
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        while value:
            lead = min(value)
            if lead not in pivots:
                inverse = pow(value[lead], prime - 2, prime)
                pivots[lead] = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                break
            factor = value[lead]
            for key, coefficient in pivots[lead].items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def project_last_three(vector: dict[int, int], length: int) -> dict[int, int]:
    n = 2 * length
    row_mask = 7 << (length - 3)
    site_mask = row_mask | (row_mask << length)
    pauli_mask = site_mask | (site_mask << n)
    return {code: coefficient for code, coefficient in vector.items() if not code & ~pauli_mask}


def gram(rows: list[dict[int, int]]) -> list[list[int]]:
    return [
        [
            sum(coefficient * right.get(code, 0) for code, coefficient in left.items())
            if len(left) <= len(right)
            else sum(coefficient * left.get(code, 0) for code, coefficient in right.items())
            for right in rows
        ]
        for left in rows
    ]


def determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            swap = next(row for row in range(column + 1, size) if work[row][column])
            work[column], work[swap] = work[swap], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = work[row][target] * pivot - work[row][column] * work[column][target]
                if column:
                    numerator, remainder = divmod(numerator, previous)
                    assert remainder == 0
                work[row][target] = numerator
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def matrix_entry(vector: dict[int, int], length: int, ket: int, bra: int) -> int:
    n = 2 * length
    mask = (1 << n) - 1
    flip = ket ^ bra
    return sum(
        coefficient * (-1 if (((pauli >> n) & mask) & ket).bit_count() & 1 else 1)
        for pauli, coefficient in vector.items()
        if (pauli & mask) == flip
    )


def evaluation_vector(vector: dict[int, int], length: int) -> dict[int, int]:
    n = 2 * length
    result = {}
    for ket in range(1 << n):
        for bra in range(1 << n):
            value = matrix_entry(vector, length, ket, bra)
            if value:
                result[(ket << n) | bra] = value
    return result


def check_base_evaluation_grams() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    stored = {row["length"]: row for row in artifact["data"]["evaluation_base_cases"]}
    for length in (2, 3):
        rows = [vector for _, _, vector in family(length)]
        coefficient_determinant = determinant(gram(rows))
        evaluations = [evaluation_vector(vector, length) for vector in rows]
        evaluation_determinant = determinant(gram(evaluations))
        scale = 1 << (2 * length)
        assert coefficient_determinant != 0
        assert evaluation_determinant == coefficient_determinant * scale ** len(rows)
        assert str(coefficient_determinant) == stored[length]["coefficient_gram_determinant"]
        assert str(evaluation_determinant) == stored[length]["evaluation_gram_determinant"]
    return "independently recomputed exact coefficient/evaluation Gram identities at L=2,3"


def check_consecutive_rank_transition() -> str:
    actual = {}
    for length in (7, 8):
        rows = [vector for _, _, vector in family(length)]
        projected = [project_last_three(row, length) for row in rows]
        actual[length] = {
            "full": [rank_mod(rows, prime) for prime in (P1, P2)],
            "projected": [rank_mod(projected, prime) for prime in (P1, P2)],
        }
    assert actual == {
        7: {"full": [128, 128], "projected": [128, 128]},
        8: {"full": [248, 248], "projected": [248, 248]},
    }
    return "recomputed exact boundary-window rank transition 128 -> 248 at L=7 -> 8"


def check_stored_exact_relations() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    rows = {label: vector for label, _, vector in family(8)}
    relations = artifact["data"]["exact_relations_L8"]
    assert len(relations) == 8
    for relation in relations:
        residual: dict[int, int] = {}
        coefficients = []
        for term in relation["terms"]:
            multiplier = int(term["coefficient"])
            coefficients.append(abs(multiplier))
            for code, coefficient in rows[term["label"]].items():
                updated = residual.get(code, 0) + multiplier * coefficient
                if updated:
                    residual[code] = updated
                else:
                    residual.pop(code, None)
        assert not residual
        assert gcd(*coefficients) == 1
        assert relation["exact_residual_support_size"] == 0
    assert sorted(relation["term_count"] for relation in relations)[:2] == [6, 7]
    return "all eight stored primitive integer relations vanish coefficientwise over Z"


def check_artifact_envelope() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    assert set(artifact) == {"provenance", "data", "checks"}
    assert artifact["provenance"]["script"] == "experiments/e65_ladder_representation.py"
    assert artifact["provenance"]["interpreter"] == ".venv/bin/python"
    assert artifact["data"]["status"] == "UNRESOLVED"
    assert artifact["data"]["target"].startswith("[UNRESOLVED]")
    assert artifact["checks"] and all(check["passed"] for check in artifact["checks"])
    reconstruction = artifact["data"]["modular_relation_reconstruction"]
    assert [record["prime"] for record in reconstruction] == list(RECONSTRUCTION_PRIMES)
    assert all(record["rank"] == 248 and record["nullity"] == 8 for record in reconstruction)
    return f"validated result envelope and {len(artifact['checks'])} passing embedded checks"


def run_check(name: str, function) -> bool:
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
            ("base evaluation Grams", check_base_evaluation_grams),
            ("consecutive rank transition", check_consecutive_rank_transition),
            ("exact null relations", check_stored_exact_relations),
            ("artifact envelope", check_artifact_envelope),
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
