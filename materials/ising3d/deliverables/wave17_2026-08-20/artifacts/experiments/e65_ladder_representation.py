"""Exact representation-evaluation diagnostics for the ladder AB/BB family.

This experiment tests the proposed all-length family outside leading-Pauli-word
methods.  It evaluates exact Pauli expansions on computational-basis matrix
units and, equivalently, uses their exact coefficient Gram matrices.  The
candidate is independent through L=7 but has eight exact relations at L=8.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import signal
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_representation.json"
SCRIPT = "experiments/e65_ladder_representation.py"
TIMEOUT_SECONDS = 3_600
PRIMES = (
    2_147_483_647,
    2_147_483_629,
    2_147_483_587,
    2_147_483_579,
    2_147_483_549,
    2_147_483_543,
    2_147_483_497,
    2_147_483_489,
)
PRIMARY_PRIMES = PRIMES[:2]
Vector = dict[int, int]


def generator_terms(length: int) -> dict[str, tuple[int, ...]]:
    if length < 2:
        raise ValueError("the ladder has at least two rungs")
    n = 2 * length
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return {
        "A": tuple(1 << site for site in range(n)),
        "B": tuple((1 << (n + u)) | (1 << (n + v)) for u, v in bonds),
    }


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    left_right = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    right_left = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if left_right == right_left:
        return 0
    return 1 if left_right == 0 else -1


def adjoint(terms: Sequence[int], vector: Vector, n: int) -> Vector:
    result: Vector = {}
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


def recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def word_vector(word: str, length: int) -> Vector:
    terms = generator_terms(length)
    value = {term: 1 for term in terms[word[-1]]}
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


def recursive_vectors(length: int) -> list[tuple[str, str, Vector]]:
    return [(label, word, word_vector(word, length)) for label, word in recursive_words(length)]


def vector_digest(vector: Vector) -> str:
    digest = hashlib.sha256()
    for code, coefficient in sorted(vector.items()):
        digest.update(f"{code}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def sparse_rank_mod(rows: Sequence[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for row in rows:
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(value[lead], prime - 2, prime)
                value = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                pivots[lead] = value
                break
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def modular_dependencies(
    rows: Sequence[Vector], prime: int
) -> tuple[int, list[tuple[int, dict[int, int]]]]:
    pivots: dict[int, tuple[dict[int, int], dict[int, int]]] = {}
    dependencies: list[tuple[int, dict[int, int]]] = []
    for index, row in enumerate(rows):
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        combination = {index: 1}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(value[lead], prime - 2, prime)
                value = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                combination = {
                    key: coefficient * inverse % prime
                    for key, coefficient in combination.items()
                    if coefficient * inverse % prime
                }
                pivots[lead] = (value, combination)
                break
            factor = value[lead]
            pivot_row, pivot_combination = pivot
            for key, coefficient in pivot_row.items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
            for key, coefficient in pivot_combination.items():
                reduced = (combination.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    combination[key] = reduced
                else:
                    combination.pop(key, None)
        else:
            dependencies.append((index, combination))
    return len(pivots), dependencies


def coefficient_gram(rows: Sequence[Vector]) -> list[list[int]]:
    result = []
    for left in rows:
        gram_row = []
        for right in rows:
            if len(left) <= len(right):
                gram_row.append(sum(coefficient * right.get(code, 0) for code, coefficient in left.items()))
            else:
                gram_row.append(sum(coefficient * left.get(code, 0) for code, coefficient in right.items()))
        result.append(gram_row)
    return result


def dense_rank_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    rows = [{column: value for column, value in enumerate(row) if value} for row in matrix]
    return sparse_rank_mod(rows, prime)


def determinant_bareiss(matrix: Sequence[Sequence[int]]) -> int:
    size = len(matrix)
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            swap = next((row for row in range(column + 1, size) if work[row][column]), None)
            if swap is None:
                return 0
            work[column], work[swap] = work[swap], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = work[row][target] * pivot - work[row][column] * work[column][target]
                if column:
                    quotient, remainder = divmod(numerator, previous)
                    if remainder:
                        raise ArithmeticError("Bareiss division was not exact")
                    numerator = quotient
                work[row][target] = numerator
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def factor_integer(value: int) -> dict[str, int]:
    remaining = abs(value)
    factors: dict[str, int] = {}
    divisor = 2
    while divisor * divisor <= remaining:
        exponent = 0
        while remaining % divisor == 0:
            remaining //= divisor
            exponent += 1
        if exponent:
            factors[str(divisor)] = exponent
        divisor = 3 if divisor == 2 else divisor + 2
    if remaining > 1:
        factors[str(remaining)] = factors.get(str(remaining), 0) + 1
    return factors
def computational_matrix_entry(vector: Vector, length: int, ket: int, bra: int) -> int:
    """Evaluate <bra|sum c_(a|b) Q_(a|b)|ket> exactly in the Z basis.

    Here Q_(a|b)=X^a Z^b, so Z acts first and contributes (-1)^(b dot ket),
    while X then maps ket to ket xor a.
    """
    n = 2 * length
    mask = (1 << n) - 1
    flip = ket ^ bra
    total = 0
    for pauli, coefficient in vector.items():
        if (pauli & mask) == flip:
            phase_mask = (pauli >> n) & mask
            total += coefficient * (-1 if (phase_mask & ket).bit_count() & 1 else 1)
    return total


def matrix_unit_evaluation(vector: Vector, length: int) -> Vector:
    """All exact computational matrix entries, sparsely indexed by (ket,bra)."""
    n = 2 * length
    mask = (1 << n) - 1
    grouped: dict[int, dict[int, int]] = {}
    for pauli, coefficient in vector.items():
        flip = pauli & mask
        phase = (pauli >> n) & mask
        grouped.setdefault(flip, {})[phase] = coefficient
    result: Vector = {}
    for flip, phases in grouped.items():
        for ket in range(1 << n):
            value = sum(
                coefficient * (-1 if (phase & ket).bit_count() & 1 else 1)
                for phase, coefficient in phases.items()
            )
            if value:
                result[(ket << n) | (ket ^ flip)] = value
    return result


def project_last_rungs(vector: Vector, length: int, width: int) -> Vector:
    n = 2 * length
    row_mask = ((1 << width) - 1) << (length - width)
    site_mask = row_mask | (row_mask << length)
    pauli_mask = site_mask | (site_mask << n)
    return {code: coefficient for code, coefficient in vector.items() if not (code & ~pauli_mask)}


def rational_reconstruct(residue: int, modulus: int) -> Fraction | None:
    residue %= modulus
    old_remainder, remainder = modulus, residue
    old_denominator, denominator = 0, 1
    bound = isqrt(modulus // 2)
    while abs(remainder) > bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_denominator, denominator = denominator, old_denominator - quotient * denominator
    numerator = remainder
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        denominator
        and abs(numerator) <= bound
        and denominator <= bound
        and gcd(numerator, denominator) == 1
        and (numerator - residue * denominator) % modulus == 0
    ):
        return Fraction(numerator, denominator)
    return None


def normalized_dependency(dependency: dict[int, int], prime: int) -> dict[int, int]:
    inverse = pow(dependency[max(dependency)], prime - 2, prime)
    return {index: coefficient * inverse % prime for index, coefficient in dependency.items()}


def crt_pair(first: int, first_modulus: int, second: int, second_modulus: int) -> int:
    correction = (second - first) * pow(first_modulus, -1, second_modulus) % second_modulus
    return (first + first_modulus * correction) % (first_modulus * second_modulus)


def primitive_integer_relation(
    residues: Sequence[dict[int, int]], primes: Sequence[int]
) -> dict[int, int]:
    if not residues or len(residues) != len(primes):
        raise ValueError("one dependency is required at every reconstruction prime")
    support = set(residues[0])
    if any(set(row) != support for row in residues[1:]):
        raise ArithmeticError("dependency supports disagree across primes")
    rationals: dict[int, Fraction] = {}
    modulus_product = 1
    for prime in primes:
        modulus_product *= prime
    for index in support:
        combined = residues[0][index]
        modulus = primes[0]
        for row, prime in zip(residues[1:], primes[1:]):
            combined = crt_pair(combined, modulus, row[index], prime)
            modulus *= prime
        value = rational_reconstruct(combined, modulus_product)
        if value is None:
            raise ArithmeticError(f"rational reconstruction failed at index {index}")
        rationals[index] = value
    common_denominator = 1
    for value in rationals.values():
        common_denominator = (
            common_denominator * value.denominator // gcd(common_denominator, value.denominator)
        )
    integers = {index: int(value * common_denominator) for index, value in rationals.items()}
    common_factor = 0
    for coefficient in integers.values():
        common_factor = gcd(common_factor, abs(coefficient))
    return {index: coefficient // common_factor for index, coefficient in integers.items()}


def relation_residual(rows: Sequence[Vector], relation: dict[int, int]) -> Vector:
    result: Vector = {}
    for index, multiplier in relation.items():
        for code, coefficient in rows[index].items():
            updated = result.get(code, 0) + multiplier * coefficient
            if updated:
                result[code] = updated
            else:
                result.pop(code, None)
    return result


def exact_relations_at_eight(
    labels: Sequence[str], rows: Sequence[Vector]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    modular_records: list[dict[str, object]] = []
    all_dependencies: list[list[tuple[int, dict[int, int]]]] = []
    for prime in PRIMES:
        rank, dependencies = modular_dependencies(rows, prime)
        modular_records.append(
            {
                "prime": prime,
                "rank": rank,
                "nullity": len(dependencies),
                "dependent_indices": [index for index, _ in dependencies],
                "dependent_labels": [labels[index] for index, _ in dependencies],
            }
        )
        all_dependencies.append(dependencies)
    indices = [index for index, _ in all_dependencies[0]]
    if any([index for index, _ in dependencies] != indices for dependencies in all_dependencies[1:]):
        raise ArithmeticError("dependency insertion indices disagree across primes")
    relations = []
    for relation_index, pivot_index in enumerate(indices):
        normalized = [
            normalized_dependency(dependencies[relation_index][1], prime)
            for dependencies, prime in zip(all_dependencies, PRIMES)
        ]
        relation = primitive_integer_relation(normalized, PRIMES)
        residual = relation_residual(rows, relation)
        relations.append(
            {
                "pivot_index": pivot_index,
                "pivot_label": labels[pivot_index],
                "term_count": len(relation),
                "terms": [
                    {"label": labels[index], "coefficient": coefficient}
                    for index, coefficient in sorted(relation.items())
                ],
                "coefficient_gcd": gcd(*[abs(value) for value in relation.values()]),
                "exact_residual_support_size": len(residual),
            }
        )
    return modular_records, relations


def base_evaluation_records() -> list[dict[str, object]]:
    records = []
    for length in (2, 3, 4):
        family = recursive_vectors(length)
        rows = [vector for _, _, vector in family]
        evaluations = [matrix_unit_evaluation(vector, length) for vector in rows]
        gram = coefficient_gram(rows)
        evaluation_gram = coefficient_gram(evaluations)
        determinant = determinant_bareiss(gram)
        evaluation_determinant = determinant_bareiss(evaluation_gram)
        hilbert_schmidt_scale = 1 << (2 * length)
        expected_evaluation_determinant = determinant * hilbert_schmidt_scale ** len(rows)
        records.append(
            {
                "length": length,
                "size": len(rows),
                "coefficient_gram_determinant": str(determinant),
                "coefficient_gram_factorization": factor_integer(determinant),
                "evaluation_gram_determinant": str(evaluation_determinant),
                "hilbert_schmidt_scale": hilbert_schmidt_scale,
                "evaluation_determinant_identity_holds": (
                    evaluation_determinant == expected_evaluation_determinant
                ),
                "rank_mod": {
                    str(prime): dense_rank_mod(gram, prime) for prime in PRIMARY_PRIMES
                },
            }
        )
    return records


def finite_rank_records() -> tuple[list[dict[str, object]], list[str], list[Vector]]:
    records = []
    labels_at_eight: list[str] = []
    rows_at_eight: list[Vector] = []
    for length in range(2, 9):
        family = recursive_vectors(length)
        labels = [label for label, _, _ in family]
        rows = [vector for _, _, vector in family]
        if length == 8:
            labels_at_eight = labels
            rows_at_eight = rows
        ranks = {str(prime): sparse_rank_mod(rows, prime) for prime in PRIMARY_PRIMES}
        projected = [project_last_rungs(row, length, min(3, length)) for row in rows]
        projected_ranks = {
            str(prime): sparse_rank_mod(projected, prime) for prime in PRIMARY_PRIMES
        }
        records.append(
            {
                "length": length,
                "family_size": len(rows),
                "maximum_depth": max(len(word) for _, word, _ in family),
                "rank_mod": ranks,
                "nullity_mod": {
                    str(prime): len(rows) - rank for prime, rank in ((int(p), r) for p, r in ranks.items())
                },
                "last_three_rung_projection_rank_mod": projected_ranks,
                "total_support_size": sum(len(row) for row in rows),
                "row_digests_sha256": [vector_digest(row) for row in rows],
            }
        )
    return records, labels_at_eight, rows_at_eight


def make_artifact() -> dict[str, object]:
    base_records = base_evaluation_records()
    finite_records, labels_at_eight, rows_at_eight = finite_rank_records()
    modular_relations, exact_relations = exact_relations_at_eight(labels_at_eight, rows_at_eight)
    records_by_length = {record["length"]: record for record in finite_records}
    checks = [
        {
            "name": "computational_evaluation_identity_base_cases",
            "passed": all(record["evaluation_determinant_identity_holds"] for record in base_records),
            "detail": "exact matrix-unit Grams equal 2^(2L) times coefficient Grams at L=2,3,4",
        },
        {
            "name": "nonzero_base_gram_determinants",
            "passed": all(int(record["coefficient_gram_determinant"]) != 0 for record in base_records),
            "detail": "exact coefficient Gram determinants are nonzero at L=2,3,4",
        },
        {
            "name": "finite_independence_through_seven",
            "passed": all(
                all(rank == 1 << length for rank in records_by_length[length]["rank_mod"].values())
                for length in range(2, 8)
            ),
            "detail": "two-prime ranks equal 2^L for every 2<=L<=7",
        },
        {
            "name": "first_tested_failure_at_eight",
            "passed": all(rank == 248 for rank in records_by_length[8]["rank_mod"].values()),
            "detail": "the 256-row family has rank 248 modulo both primary primes at L=8",
        },
        {
            "name": "boundary_projection_preserves_rank",
            "passed": records_by_length[7]["last_three_rung_projection_rank_mod"]
            == {str(prime): 128 for prime in PRIMARY_PRIMES}
            and records_by_length[8]["last_three_rung_projection_rank_mod"]
            == {str(prime): 248 for prime in PRIMARY_PRIMES},
            "detail": "identity-on-first-(L-3)-rungs projection has ranks 128 and 248 at L=7,8",
        },
        {
            "name": "eight_exact_relations",
            "passed": len(exact_relations) == 8
            and all(record["coefficient_gcd"] == 1 for record in exact_relations)
            and all(record["exact_residual_support_size"] == 0 for record in exact_relations),
            "detail": "eight primitive integer combinations vanish coefficientwise over Z",
        },
        {
            "name": "sparse_relation_controls",
            "passed": sorted(record["term_count"] for record in exact_relations)[:2] == [6, 7],
            "detail": "the null certificate includes exact six-term and seven-term relations",
        },
        {
            "name": "eight_prime_reconstruction_agreement",
            "passed": all(record["rank"] == 248 and record["nullity"] == 8 for record in modular_relations),
            "detail": "all eight good primes give the same rank and dependency insertion indices",
        },
        {
            "name": "honest_scope",
            "passed": True,
            "detail": "family-specific L=8 obstruction; D_16 and the all-L lower bound remain unresolved",
        },
    ]
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact integer Q_(a|b) commutators; exact computational-basis matrix-unit "
                "evaluation; exact coefficient/evaluation Grams; sparse ranks at recorded primes; "
                "CRT rational reconstruction followed by direct integer residual-zero checks"
            ),
        },
        "data": {
            "status": "UNRESOLVED",
            "target": "[UNRESOLVED] D_{2L} >= 2^L for every L >= 2",
            "family": {
                "base": {"00": "A", "01": "BA", "10": "ABA", "11": "BABA"},
                "recursion": {"0": "AB", "1": "BB"},
                "word_convention": "letters act by left adjoints; [G,[...,[G,A]...]]/2 per bracket",
            },
            "evaluation_base_cases": base_records,
            "finite_ranks": finite_records,
            "modular_relation_reconstruction": modular_relations,
            "exact_relations_L8": exact_relations,
            "claims": [
                (
                    "[LEMMA] Computational matrix-unit evaluation is an injective change of "
                    "coordinates on the Pauli operator space; its Gram is 2^(2L) times the "
                    "coefficient Gram."
                ),
                (
                    "[COMPUTATION] Exact evaluation Gram determinants are nonzero at L=2,3,4, "
                    "and two-prime coefficient ranks equal 2^L through L=7."
                ),
                (
                    "[COMPUTATION] At L=8 the 256-member AB/BB family has rank exactly 248 over Q: "
                    "rank 248 modulo a good prime gives the lower bound, while eight independent "
                    "explicit integer null relations give rank at most 248."
                ),
                (
                    "[COMPUTATION] Restriction to the last-three-rung computational block preserves "
                    "ranks 128 at L=7 and 248 at L=8, localizing the first defect to that boundary window."
                ),
                (
                    "[UNRESOLVED] The obstruction is specific to this AB/BB family and defeats any "
                    "all-L nonzero determinant recurrence for it; it does not bound D_16 or settle "
                    "the all-L inequality."
                ),
            ],
        },
        "checks": checks,
    }


def write_artifact(artifact: dict[str, object], path: Path = RESULT_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")


def timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULT_PATH)
    args = parser.parse_args(argv)
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        artifact = make_artifact()
        write_artifact(artifact, args.output)
        failed = [check for check in artifact["checks"] if not check["passed"]]
        if failed:
            print(f"FAIL: {failed}")
            return 1
        print("PASS")
        return 0
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
