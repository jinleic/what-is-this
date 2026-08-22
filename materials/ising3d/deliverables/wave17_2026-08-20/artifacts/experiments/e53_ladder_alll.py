"""Exact leading-term diagnostics for the two-rung ladder all-length conjecture.

The experiment recomputes depth-filtered Pauli bases, tests the natural binary
extension (AB, BB), and records an exact minimal failure of triangular leading
terms.  It deliberately reports UNRESOLVED: finite rank is not induction.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Callable, Iterable, Sequence

from ising.clifford.fast_lie import GeneratedAlgebra, P1, P2

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_alll.json"
SCRIPT = "experiments/e53_ladder_alll.py"
PRIMES = (P1, P2)
LABELS = "IXYZ"
TO_BITS = {"I": (0, 0), "X": (1, 0), "Z": (0, 1), "Y": (1, 1)}
FROM_BITS = {value: key for key, value in TO_BITS.items()}
Vector = dict[int, int]
OrderKey = Callable[[int, int], tuple[int, ...]]


def ladder_bonds(length: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    if length < 2:
        raise ValueError("the ladder has at least two rungs")
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return 2 * length, tuple(bonds)


def generator_terms(length: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    n, bonds = ladder_bonds(length)
    a_terms = tuple(1 << site for site in range(n))
    b_terms = tuple((1 << (n + lhs)) | (1 << (n + rhs)) for lhs, rhs in bonds)
    return a_terms, b_terms


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    lr = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    rl = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if lr == rl:
        return 0
    return 1 if lr == 0 else -1


def adjoint(terms: Sequence[int], vector: Vector, n: int) -> Vector:
    out: Vector = {}
    for term in terms:
        for pauli, coefficient in vector.items():
            sign = sign_half(term, pauli, n)
            if sign:
                target = term ^ pauli
                out[target] = out.get(target, 0) + sign * coefficient
                if not out[target]:
                    del out[target]
    return out


def generators(length: int) -> dict[str, Vector]:
    a_terms, b_terms = generator_terms(length)
    return {"A": {term: 1 for term in a_terms}, "B": {term: 1 for term in b_terms}}


def word_vector(word: str, length: int) -> Vector:
    gens = generators(length)
    terms = {letter: tuple(gens[letter]) for letter in "AB"}
    value = dict(gens[word[-1]])
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


def recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def recursive_vectors(length: int) -> list[tuple[str, str, Vector]]:
    return [(label, word, word_vector(word, length)) for label, word in recursive_words(length)]


def site_labels(pauli: int, length: int) -> tuple[str, ...]:
    n = 2 * length
    mask = (1 << n) - 1
    a, b = pauli & mask, (pauli >> n) & mask
    return tuple(FROM_BITS[((a >> site) & 1, (b >> site) & 1)] for site in range(n))


def rung_string(pauli: int, length: int) -> str:
    decoded = site_labels(pauli, length)
    return " ".join(decoded[col] + decoded[length + col] for col in range(length))


def monomial_key(
    pauli: int,
    length: int,
    rung_reverse: bool = True,
    row_reverse: bool = False,
    alphabet_ascending: str = "IZYX",
) -> tuple[int, ...]:
    decoded = site_labels(pauli, length)
    rank = {letter: index for index, letter in enumerate(alphabet_ascending)}
    columns: Iterable[int] = range(length - 1, -1, -1) if rung_reverse else range(length)
    key = []
    for col in columns:
        pair = (decoded[col], decoded[length + col])
        if row_reverse:
            pair = pair[::-1]
        key.extend(rank[letter] for letter in pair)
    return tuple(key)


def vector_digest(vector: Vector) -> str:
    digest = hashlib.sha256()
    for pauli, coefficient in sorted(vector.items()):
        digest.update(f"{pauli}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def rank_mod(rows: Sequence[Vector], prime: int) -> int:
    pivots: dict[int, Vector] = {}
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


def depth_basis(length: int) -> list[dict[str, object]]:
    """A modular echelon basis of all words through depth 2L, with fixed leaders."""
    key: OrderKey = lambda pauli, current: monomial_key(pauli, current)
    gens = generators(length)
    terms = {letter: tuple(gens[letter]) for letter in "AB"}
    frontier = [(letter, dict(gens[letter])) for letter in "AB"]
    pivots: dict[int, Vector] = {}
    records = []
    for depth in range(1, 2 * length + 1):
        children = []
        for word, raw in frontier:
            value = {pauli: coefficient % P1 for pauli, coefficient in raw.items() if coefficient % P1}
            while value:
                lead = max(value, key=lambda pauli: key(pauli, length))
                if lead not in pivots:
                    inverse = pow(value[lead], P1 - 2, P1)
                    value = {
                        pauli: coefficient * inverse % P1
                        for pauli, coefficient in value.items()
                        if coefficient * inverse % P1
                    }
                    pivots[lead] = value
                    records.append(
                        {
                            "word": word,
                            "depth": depth,
                            "leader_code": lead,
                            "leader_pauli": rung_string(lead, length),
                            "raw_support_size": len(raw),
                            "raw_sha256": vector_digest(raw),
                        }
                    )
                    break
                factor = value[lead]
                for pauli, coefficient in pivots[lead].items():
                    reduced = (value.get(pauli, 0) - factor * coefficient) % P1
                    if reduced:
                        value[pauli] = reduced
                    else:
                        value.pop(pauli, None)
            if depth < 2 * length:
                for outer in "AB":
                    image = adjoint(terms[outer], raw, 2 * length)
                    if image:
                        children.append((outer + word, image))
        frontier = children
    return records


def rational_echelon(rows: Sequence[Vector], length: int, key: OrderKey) -> list[dict[int, Fraction]]:
    pivots: dict[int, tuple[dict[int, Fraction], dict[int, Fraction]]] = {}
    combinations = []
    for index, row in enumerate(rows):
        value = {pauli: Fraction(coefficient) for pauli, coefficient in row.items() if coefficient}
        combination = {index: Fraction(1)}
        while value:
            lead = max(value, key=lambda pauli: key(pauli, length))
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
                reduced = value.get(pauli, Fraction(0)) - factor * coefficient
                if reduced:
                    value[pauli] = reduced
                else:
                    value.pop(pauli, None)
            for row_index, coefficient in basis_combination.items():
                reduced = combination.get(row_index, Fraction(0)) - factor * coefficient
                if reduced:
                    combination[row_index] = reduced
                else:
                    combination.pop(row_index, None)
    return combinations


def combine(combination: dict[int, Fraction], rows: Sequence[Vector]) -> dict[int, Fraction]:
    value: dict[int, Fraction] = {}
    for row_index, factor in combination.items():
        for pauli, coefficient in rows[row_index].items():
            total = value.get(pauli, Fraction(0)) + factor * coefficient
            if total:
                value[pauli] = total
            else:
                value.pop(pauli, None)
    return value


def triangular_step(length: int, key: OrderKey) -> dict[str, object]:
    parents = recursive_vectors(length)
    children = recursive_vectors(length + 1)
    combinations = rational_echelon([row[2] for row in parents], length, key)
    lifted = []
    for parent_index, combination in enumerate(combinations):
        for bit in (0, 1):
            child_combination = {
                2 * row_index + bit: factor for row_index, factor in combination.items()
            }
            value = combine(child_combination, [row[2] for row in children])
            leader = max(value, key=lambda pauli: key(pauli, length + 1))
            lifted.append((parent_index, bit, value, leader))
    first_by_leader = {}
    collision = None
    for parent_index, bit, value, leader in lifted:
        if leader in first_by_leader and collision is None:
            previous_parent, previous_bit, previous_value = first_by_leader[leader]
            support_pair = None
            determinant = Fraction(0)
            supports = sorted(set(previous_value) | set(value))
            for left_index, left in enumerate(supports):
                for right in supports[left_index + 1 :]:
                    determinant = (
                        previous_value.get(left, 0) * value.get(right, 0)
                        - previous_value.get(right, 0) * value.get(left, 0)
                    )
                    if determinant:
                        support_pair = (left, right)
                        break
                if determinant:
                    break
            collision = {
                "first_parent_index": previous_parent,
                "first_bit": previous_bit,
                "second_parent_index": parent_index,
                "second_bit": bit,
                "leader_code": leader,
                "leader_pauli": rung_string(leader, length + 1),
                "first_leader_coefficient": str(previous_value[leader]),
                "second_leader_coefficient": str(value[leader]),
                "nonproportional_minor": {
                    "support_codes": list(support_pair) if support_pair else [],
                    "support_paulis": (
                        [rung_string(code, length + 1) for code in support_pair]
                        if support_pair
                        else []
                    ),
                    "determinant": str(determinant),
                },
            }
        else:
            first_by_leader[leader] = (parent_index, bit, value)
    return {
        "from_length": length,
        "to_length": length + 1,
        "lifted_count": len(lifted),
        "distinct_leaders": len(first_by_leader),
        "collision": collision,
    }


def scan_base_step_orders() -> dict[str, object]:
    histogram: Counter[int] = Counter()
    examples = []
    for rung_reverse, row_reverse, permutation in itertools.product(
        (False, True), (False, True), itertools.permutations(LABELS)
    ):
        alphabet = "".join(permutation)
        key = (
            lambda pauli, length, rr=rung_reverse, row=row_reverse, alpha=alphabet:
            monomial_key(pauli, length, rr, row, alpha)
        )
        record = triangular_step(2, key)
        count = int(record["distinct_leaders"])
        histogram[count] += 1
        if count == 7 and len(examples) < 4:
            examples.append(
                {
                    "rung_reverse": rung_reverse,
                    "row_reverse": row_reverse,
                    "alphabet_ascending": alphabet,
                }
            )
    return {
        "order_class": "2 rung directions x 2 row directions x 24 local Pauli orders",
        "orders_checked": 96,
        "distinct_lifted_leader_histogram": {
            str(key): histogram[key] for key in sorted(histogram)
        },
        "maximum_distinct_leaders": max(histogram),
        "required": 8,
        "examples_attaining_maximum": examples,
    }


def engine_controls() -> list[dict[str, int]]:
    rows = []
    for length in (2, 3):
        n, _ = ladder_bonds(length)
        a_terms, b_terms = generator_terms(length)
        singles = list(a_terms) + list(b_terms)
        a = {term: 1 for term in a_terms}
        b = {term: 1 for term in b_terms}
        algebra = GeneratedAlgebra([a, b], singles, n, p=P1)
        dimension = algebra.closure_dim(
            [a, b],
            [list(range(len(a_terms))), list(range(len(a_terms), len(singles)))],
        )
        rows.append(
            {
                "length": length,
                "full_DLA_dimension_mod_p1": dimension,
                "support_closure_size": algebra.m,
            }
        )
    return rows


def make_artifact() -> dict[str, object]:
    depth_records = []
    expected_depth_ranks = {2: 6, 3: 16, 4: 38, 5: 85}
    for length in range(2, 6):
        basis = depth_basis(length)
        depth_records.append({"length": length, "rank_mod_p1": len(basis), "basis": basis})

    family_records = []
    for length in range(2, 8):
        vectors = recursive_vectors(length)
        rows = []
        for label, word, vector in vectors:
            leader = max(vector, key=lambda pauli: monomial_key(pauli, length))
            rows.append(
                {
                    "label": label,
                    "word": word,
                    "depth": len(word),
                    "support_size": len(vector),
                    "sha256": vector_digest(vector),
                    "leader_code": leader,
                    "leader_pauli": rung_string(leader, length),
                    "leader_coefficient": vector[leader],
                }
            )
        family_records.append(
            {
                "length": length,
                "count": len(rows),
                "rank_mod": {
                    str(prime): rank_mod([row[2] for row in vectors], prime)
                    for prime in PRIMES
                },
                "distinct_raw_leaders": len({row["leader_code"] for row in rows}),
                "rows": rows,
            }
        )

    fixed_key: OrderKey = lambda pauli, length: monomial_key(pauli, length)
    triangular = [triangular_step(length, fixed_key) for length in range(2, 7)]
    order_scan = scan_base_step_orders()
    controls = engine_controls()

    checks = []
    for record in depth_records:
        length = int(record["length"])
        checks.append(
            {
                "name": f"depth_basis_L{length}",
                "passed": record["rank_mod_p1"] == expected_depth_ranks[length],
                "detail": f"rank mod {P1}={record['rank_mod_p1']}; expected={expected_depth_ranks[length]}",
            }
        )
    for record in family_records:
        length = int(record["length"])
        target = 1 << length
        checks.append(
            {
                "name": f"recursive_family_rank_L{length}",
                "passed": all(rank == target for rank in record["rank_mod"].values()),
                "detail": f"two-prime ranks={record['rank_mod']}; target={target}",
            }
        )
    checks.extend(
        [
            {
                "name": "minimal_fixed_order_triangular_collision",
                "passed": triangular[0]["distinct_leaders"] == 7
                and triangular[0]["collision"]["nonproportional_minor"]["determinant"]
                == "6",
                "detail": (
                    "L=2->3 has seven leaders for eight lifted rows; "
                    "exact nonproportional minor determinant 6"
                ),
            },
            {
                "name": "all_product_lex_orders_fail_base_step",
                "passed": order_scan["maximum_distinct_leaders"] == 7,
                "detail": (
                    f"96-order maximum={order_scan['maximum_distinct_leaders']}; required=8"
                ),
            },
            {
                "name": "repository_engine_controls",
                "passed": controls
                == [
                    {
                        "length": 2,
                        "full_DLA_dimension_mod_p1": 11,
                        "support_closure_size": 56,
                    },
                    {
                        "length": 3,
                        "full_DLA_dimension_mod_p1": 263,
                        "support_closure_size": 1056,
                    },
                ],
                "detail": f"GeneratedAlgebra controls={controls}",
            },
            {
                "name": "honest_scope",
                "passed": True,
                "detail": "finite computations and an order-class obstruction; no all-L conclusion",
            },
        ]
    )

    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact integer Q_(a|b) commutators; fractions.Fraction triangular lifts; "
                f"modular echelon lower bounds at primes {P1},{P2}; "
                "repository GeneratedAlgebra controls"
            ),
        },
        "data": {
            "status": "UNRESOLVED",
            "target": "[UNRESOLVED] D_{2L} >= 2^L for every L >= 2",
            "fixed_order": "right rung first; top before bottom; local descending X>Y>Z>I",
            "depth_bases": depth_records,
            "recursive_family": family_records,
            "triangular_steps": triangular,
            "product_lex_base_step_scan": order_scan,
            "repository_engine_controls": controls,
            "claims": [
                (
                    "[COMPUTATION] The AB/BB family has rank 2^L modulo both recorded "
                    "primes for 2<=L<=7; each rank is a rigorous lower bound over Q."
                ),
                (
                    "[COMPUTATION] Under the fixed order its rational-echelon lifts "
                    "already have only 7 distinct leaders at L=2->3, with an exact "
                    "determinant-6 witness that the colliding rows remain independent."
                ),
                (
                    "[COMPUTATION] Every one of 96 product-lex orders fails that same "
                    "base lift; the maximum is 7 leaders, not 8."
                ),
                (
                    "[UNRESOLVED] This finite order-class obstruction neither proves "
                    "nor refutes the all-length inequality."
                ),
            ],
        },
        "checks": checks,
    }


def write_artifact(artifact: dict[str, object]) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_PATH.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    try:
        artifact = make_artifact()
        write_artifact(artifact)
        failed = [check for check in artifact["checks"] if not check["passed"]]
        if failed:
            print(f"FAIL: {failed}")
            return 1
        print("PASS")
        return 0
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
