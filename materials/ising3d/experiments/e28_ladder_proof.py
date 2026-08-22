"""Independent exact search for an all-length two-leg-ladder Lie lower bound.

This module deliberately does not import ``e20_algebra_growth``.  It expands actual
left-normed bracket words in the full Pauli basis, keeps exact integer coefficients,
and uses modular Gaussian elimination only for rank certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


PRIME_31 = 2_147_483_647
SECOND_PRIME_31 = 2_147_483_629
PRIMES = (PRIME_31, SECOND_PRIME_31)
SCRIPT = "experiments/e28_ladder_proof.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_proof.json"

Word = tuple[str, ...]  # outermost generator first; e.g. ABA = [A,[B,A]]
Vector = dict[int, int]


def ladder_bonds(length: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Return the open 2 x length ladder in row-major site order."""

    if length < 2:
        raise ValueError("a ladder must have length at least two")
    n = 2 * length
    bonds: list[tuple[int, int]] = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return n, tuple(bonds)


def chain_bonds(n: int, periodic: bool = False) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Return an open chain or a simple ring (the controls use n >= 3 for rings)."""

    if n < 2:
        raise ValueError("a chain must have at least two sites")
    bonds = [(site, site + 1) for site in range(n - 1)]
    if periodic:
        if n < 3:
            raise ValueError("the simple-ring control requires n at least three")
        bonds.append((n - 1, 0))
    return n, tuple(bonds)


def generator_terms(
    n: int, bonds: Sequence[tuple[int, int]]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Encode the Pauli terms of A and B as Q_(a|b) = X^a Z^b."""

    a_terms = tuple(1 << site for site in range(n))
    b_terms = tuple((1 << (n + left)) | (1 << (n + right)) for left, right in bonds)
    return a_terms, b_terms


def commutator_sign_half(generator: int, pauli: int, n: int) -> int:
    """Return the coefficient 0 or +/-1 in [Q_generator,Q_pauli]/2."""

    mask = (1 << n) - 1
    first = (((generator >> n) & mask) & (pauli & mask)).bit_count() & 1
    second = (((pauli >> n) & mask) & (generator & mask)).bit_count() & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def seed_vector(terms: Iterable[int]) -> Vector:
    """Build an exact-integer Pauli vector, retaining term multiplicities."""

    result: Vector = {}
    for term in terms:
        result[term] = result.get(term, 0) + 1
    return result


def adjoint(
    terms: Sequence[int], vector: Vector, n: int, modulus: int | None = None
) -> Vector:
    """Compute [sum(terms), vector]/2 exactly or modulo an odd prime."""

    result: Vector = {}
    for pauli, coefficient in vector.items():
        for generator in terms:
            sign = commutator_sign_half(generator, pauli, n)
            if not sign:
                continue
            target = generator ^ pauli
            value = result.get(target, 0) + sign * coefficient
            if modulus is not None:
                value %= modulus
            if value:
                result[target] = value
            else:
                result.pop(target, None)
    return result


def vector_for_word(
    word: Word,
    n: int,
    bonds: Sequence[tuple[int, int]],
    modulus: int | None = None,
) -> Vector:
    """Expand one actual left-normed bracket word in the Pauli basis."""

    if not word or any(letter not in ("A", "B") for letter in word):
        raise ValueError(f"invalid bracket word {word!r}")
    a_terms, b_terms = generator_terms(n, bonds)
    generators = {"A": a_terms, "B": b_terms}
    vector = seed_vector(generators[word[-1]])
    if modulus is not None:
        vector = {key: value % modulus for key, value in vector.items() if value % modulus}
    for letter in reversed(word[:-1]):
        vector = adjoint(generators[letter], vector, n, modulus)
    return vector


@dataclass
class ModularEchelon:
    """Sparse row echelon basis over one explicitly recorded prime field."""

    prime: int

    def __post_init__(self) -> None:
        self.rows: dict[int, Vector] = {}

    @property
    def rank(self) -> int:
        return len(self.rows)

    def add(self, vector: Vector) -> tuple[int, Vector] | None:
        value = {
            pauli: coefficient % self.prime
            for pauli, coefficient in vector.items()
            if coefficient % self.prime
        }
        while value:
            pivot = max(value)
            existing = self.rows.get(pivot)
            if existing is None:
                inverse = pow(value[pivot], self.prime - 2, self.prime)
                value = {
                    pauli: (coefficient * inverse) % self.prime
                    for pauli, coefficient in value.items()
                }
                self.rows[pivot] = value
                return pivot, value
            factor = value[pivot]
            for pauli, coefficient in existing.items():
                reduced = (value.get(pauli, 0) - factor * coefficient) % self.prime
                if reduced:
                    value[pauli] = reduced
                else:
                    value.pop(pauli, None)
        return None


def word_code(word: Word) -> str:
    return "".join(word)


def bracket_notation(word: Word) -> str:
    expression = word[-1]
    for letter in reversed(word[:-1]):
        expression = f"[{letter},{expression}]"
    return expression


def pauli_string(pauli: int, n: int) -> str:
    """Decode Q_(a|b), in increasing site order, using I/X/Z/Y local labels."""

    mask = (1 << n) - 1
    a_mask = pauli & mask
    b_mask = (pauli >> n) & mask
    alphabet = "IXZY"
    return "".join(
        alphabet[((a_mask >> site) & 1) + 2 * ((b_mask >> site) & 1)]
        for site in range(n)
    )


def rung_pauli_string(pauli: int, length: int) -> str:
    decoded = pauli_string(pauli, 2 * length)
    return " ".join(decoded[col] + decoded[length + col] for col in range(length))


def touched_rung_mask(pauli: int, length: int) -> int:
    decoded = pauli_string(pauli, 2 * length)
    mask = 0
    for col in range(length):
        if decoded[col] != "I" or decoded[length + col] != "I":
            mask |= 1 << col
    return mask


def boundary_lex_key(pauli: int, length: int) -> tuple[int, ...]:
    """Uniform candidate order: right rung first and X > Y > I > Z locally."""

    decoded = pauli_string(pauli, 2 * length)
    rank = {"Z": 0, "I": 1, "Y": 2, "X": 3}
    key: list[int] = []
    for col in range(length - 1, -1, -1):
        key.extend((rank[decoded[col]], rank[decoded[length + col]]))
    return tuple(key)


def vector_digest(vector: Vector) -> str:
    digest = hashlib.sha256()
    for pauli, coefficient in sorted(vector.items()):
        digest.update(f"{pauli}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def independent_record(
    length: int,
    depth: int,
    word: Word,
    vector: Vector,
    modular_pivot: int,
    include_full_support: bool,
) -> dict[str, object]:
    raw_leader = max(vector, key=lambda value: boundary_lex_key(value, length))
    rung_histogram: dict[str, int] = defaultdict(int)
    for pauli in vector:
        rung_histogram[str(touched_rung_mask(pauli, length))] += 1
    record: dict[str, object] = {
        "depth": depth,
        "word": word_code(word),
        "bracket": bracket_notation(word),
        "support_size": len(vector),
        "support_sha256": vector_digest(vector),
        "boundary_lex_leader": raw_leader,
        "boundary_lex_leader_pauli": rung_pauli_string(raw_leader, length),
        "boundary_lex_leader_rung_mask": touched_rung_mask(raw_leader, length),
        "modular_echelon_pivot": modular_pivot,
        "modular_echelon_pivot_pauli": rung_pauli_string(modular_pivot, length),
        "modular_echelon_pivot_rung_mask": touched_rung_mask(modular_pivot, length),
        "support_rung_mask_histogram": dict(sorted(rung_histogram.items(), key=lambda item: int(item[0]))),
    }
    if include_full_support:
        record["support_codes"] = sorted(vector)
    return record


def recursive_candidate_words(length: int) -> list[tuple[str, Word]]:
    """Return the tested four-base, two-child family of 2^length actual words."""

    if length < 2:
        raise ValueError("candidate family starts at length two")
    family: list[tuple[str, Word]] = [
        ("00", ("A",)),
        ("01", ("B", "A")),
        ("10", ("A", "B", "A")),
        ("11", ("B", "A", "B", "A")),
    ]
    for _ in range(2, length):
        children: list[tuple[str, Word]] = []
        for label, word in family:
            children.append((label + "0", ("A", "B") + word))
            children.append((label + "1", ("B", "B") + word))
        family = children
    return family


def depth_certificate(
    length: int,
    max_depth: int | None = None,
    include_full_supports: bool = False,
) -> dict[str, object]:
    """Expand every actual word and independently certify D_k at both primes."""

    if max_depth is None:
        max_depth = 2 * length
    n, bonds = ladder_bonds(length)
    a_terms, b_terms = generator_terms(n, bonds)
    frontier: list[tuple[Word, Vector]] = [
        (("A",), seed_vector(a_terms)),
        (("B",), seed_vector(b_terms)),
    ]
    bases = {prime: ModularEchelon(prime) for prime in PRIMES}
    dimensions = {prime: [] for prime in PRIMES}
    independent_words: list[dict[str, object]] = []
    first_by_boundary_leader: dict[int, tuple[Word, Vector]] = {}
    first_nonproportional_collision: dict[str, object] | None = None
    candidate_word_set = {word for _, word in recursive_candidate_words(length)}
    captured_candidates: dict[Word, Vector] = {}
    nonzero_words = 0

    for depth in range(1, max_depth + 1):
        next_frontier: list[tuple[Word, Vector]] = []
        for word, vector in frontier:
            if not vector:
                continue
            nonzero_words += 1
            if word in candidate_word_set:
                captured_candidates[word] = vector
            leader = max(vector, key=lambda value: boundary_lex_key(value, length))
            previous = first_by_boundary_leader.get(leader)
            if previous is None:
                first_by_boundary_leader[leader] = (word, vector)
            elif first_nonproportional_collision is None:
                previous_word, previous_vector = previous
                collision_basis = ModularEchelon(PRIME_31)
                collision_basis.add(previous_vector)
                collision_basis.add(vector)
                if collision_basis.rank == 2:
                    first_nonproportional_collision = {
                        "leader": leader,
                        "leader_pauli": rung_pauli_string(leader, length),
                        "first_word": word_code(previous_word),
                        "first_bracket": bracket_notation(previous_word),
                        "second_word": word_code(word),
                        "second_bracket": bracket_notation(word),
                        "pair_rank_mod_prime": 2,
                    }
            additions = {prime: bases[prime].add(vector) for prime in PRIMES}
            primary = additions[PRIME_31]
            if primary is not None:
                pivot, _ = primary
                independent_words.append(
                    independent_record(
                        length,
                        depth,
                        word,
                        vector,
                        pivot,
                        include_full_supports,
                    )
                )
            if depth < max_depth:
                next_frontier.append((("A",) + word, adjoint(a_terms, vector, n)))
                next_frontier.append((("B",) + word, adjoint(b_terms, vector, n)))
        for prime in PRIMES:
            dimensions[prime].append(bases[prime].rank)
        frontier = next_frontier

    candidate_rows: list[dict[str, object]] = []
    candidate_bases = {prime: ModularEchelon(prime) for prime in PRIMES}
    candidate_leaders: dict[int, str] = {}
    candidate_collision: dict[str, object] | None = None
    family_lookup = dict(recursive_candidate_words(length))
    for label, word in recursive_candidate_words(length):
        vector = captured_candidates.get(word)
        if vector is None:
            raise AssertionError(f"candidate word {word_code(word)} was not captured")
        raw_leader = max(vector, key=lambda value: boundary_lex_key(value, length))
        previous_label = candidate_leaders.setdefault(raw_leader, label)
        if previous_label != label and candidate_collision is None:
            candidate_collision = {
                "first_label": previous_label,
                "second_label": label,
                "leader": raw_leader,
                "leader_pauli": rung_pauli_string(raw_leader, length),
                "first_word": word_code(family_lookup[previous_label]),
                "second_word": word_code(word),
            }
        pivots: dict[str, int | None] = {}
        for prime in PRIMES:
            added = candidate_bases[prime].add(vector)
            pivots[str(prime)] = None if added is None else added[0]
        candidate_rows.append(
            {
                "label": label,
                "word": word_code(word),
                "depth": len(word),
                "support_size": len(vector),
                "support_sha256": vector_digest(vector),
                "boundary_lex_leader": raw_leader,
                "boundary_lex_leader_pauli": rung_pauli_string(raw_leader, length),
                "modular_pivots": pivots,
            }
        )

    return {
        "length": length,
        "n_sites": n,
        "max_depth": max_depth,
        "dimensions": {str(prime): values for prime, values in dimensions.items()},
        "D_2L": {str(prime): bases[prime].rank for prime in PRIMES},
        "nonzero_word_count": nonzero_words,
        "boundary_lex_distinct_leaders": len(first_by_boundary_leader),
        "boundary_lex_target": 1 << length,
        "boundary_lex_capacity_passed": len(first_by_boundary_leader) >= (1 << length),
        "boundary_lex_first_nonproportional_collision": first_nonproportional_collision,
        "maximal_independent_words_mod_prime": independent_words,
        "independent_support_encoding": (
            "Each support_codes list is the complete set of nonzero Pauli codes. "
            "Exact integer coefficients are bound by support_sha256, computed from sorted "
            "'code:coefficient;' records."
        ),
        "recursive_candidate": {
            "definition": "base 00:A, 01:BA, 10:ABA, 11:BABA; append bit 0 by ABw and bit 1 by BBw",
            "count": len(candidate_rows),
            "rank": {str(prime): candidate_bases[prime].rank for prime in PRIMES},
            "distinct_boundary_lex_leaders": len(candidate_leaders),
            "first_boundary_lex_collision": candidate_collision,
            "rows": candidate_rows,
        },
    }


def closure_profile(
    n: int,
    bonds: Sequence[tuple[int, int]],
    prime: int = PRIME_31,
    max_depth: int = 100,
) -> list[int]:
    """Close the two-generator span by acting only on each new echelon frontier."""

    a_terms, b_terms = generator_terms(n, bonds)
    basis = ModularEchelon(prime)
    frontier: list[Vector] = []
    for seed in (seed_vector(a_terms), seed_vector(b_terms)):
        added = basis.add(seed)
        if added is not None:
            frontier.append(added[1])
    dimensions = [basis.rank]
    for _depth in range(2, max_depth + 1):
        new_frontier: list[Vector] = []
        for vector in frontier:
            for terms in (a_terms, b_terms):
                image = adjoint(terms, vector, n, prime)
                added = basis.add(image)
                if added is not None:
                    new_frontier.append(added[1])
        dimensions.append(basis.rank)
        if not new_frontier:
            return dimensions
        frontier = new_frontier
    raise RuntimeError(f"closure did not saturate by depth {max_depth}")


def control_records() -> dict[str, object]:
    records: dict[str, object] = {"open_chains": [], "rings": []}
    for n in range(2, 7):
        _, bonds = chain_bonds(n)
        profiles = {str(prime): closure_profile(n, bonds, prime) for prime in PRIMES}
        records["open_chains"].append(
            {
                "n": n,
                "profiles": profiles,
                "dimensions": {key: values[-1] for key, values in profiles.items()},
                "expected": n * n,
            }
        )
    for n in range(3, 7):
        _, bonds = chain_bonds(n, periodic=True)
        profiles = {str(prime): closure_profile(n, bonds, prime) for prime in PRIMES}
        records["rings"].append(
            {
                "n": n,
                "profiles": profiles,
                "dimensions": {key: values[-1] for key, values in profiles.items()},
                "expected": 3 * n - 1,
            }
        )
    n, bonds = ladder_bonds(2)
    profiles = {str(prime): closure_profile(n, bonds, prime) for prime in PRIMES}
    records["ladder_2x2"] = {
        "profiles": profiles,
        "dimensions": {key: values[-1] for key, values in profiles.items()},
        "expected": 11,
    }
    return records


def make_artifact(include_full_supports: bool = True) -> dict[str, object]:
    expected = {2: 6, 3: 16, 4: 38, 5: 85, 6: 197}
    ladders = [
        depth_certificate(length, include_full_supports=include_full_supports)
        for length in range(2, 7)
    ]
    controls = control_records()
    checks: list[dict[str, object]] = []
    for record in ladders:
        length = int(record["length"])
        dimensions = record["D_2L"]
        checks.append(
            {
                "name": f"independent_depth_rank_L{length}",
                "passed": all(value == expected[length] for value in dimensions.values()),
                "detail": f"D_{{2L}}={dimensions}; expected {expected[length]}",
            }
        )
        candidate_rank = record["recursive_candidate"]["rank"]
        checks.append(
            {
                "name": f"recursive_family_rank_L{length}",
                "passed": all(value == (1 << length) for value in candidate_rank.values()),
                "detail": f"rank={candidate_rank}; target={1 << length}",
            }
        )
    for family_name in ("open_chains", "rings"):
        for row in controls[family_name]:
            checks.append(
                {
                    "name": f"control_{family_name}_{row['n']}",
                    "passed": all(value == row["expected"] for value in row["dimensions"].values()),
                    "detail": f"dimensions={row['dimensions']}; expected={row['expected']}",
                }
            )
    ladder_control = controls["ladder_2x2"]
    checks.append(
        {
            "name": "control_ladder_2x2_total",
            "passed": all(
                value == ladder_control["expected"]
                for value in ladder_control["dimensions"].values()
            ),
            "detail": (
                f"dimensions={ladder_control['dimensions']}; "
                f"expected={ladder_control['expected']}"
            ),
        }
    )
    pivot_verified = [
        int(record["length"])
        for record in ladders
        if record["boundary_lex_capacity_passed"]
    ]
    checks.append(
        {
            "name": "boundary_lex_rule_finite_counterexample",
            "passed": pivot_verified == [2, 3]
            and ladders[2]["boundary_lex_distinct_leaders"] < (1 << 4),
            "detail": (
                f"verified lengths={pivot_verified}; L=4 has "
                f"{ladders[2]['boundary_lex_distinct_leaders']} leaders for target 16"
            ),
        }
    )
    return {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": (
                "exact Python integers for Pauli coefficients; sparse elimination over "
                f"the two 31-bit primes {PRIME_31} and {SECOND_PRIME_31}"
            ),
            "primes": list(PRIMES),
            "pauli_encoding": "integer a | (b << n) represents Q_(a|b)=X^a Z^b",
            "word_encoding": "outermost generator first",
        },
        "data": {
            "claim_status": (
                "[COMPUTATION] The boundary-lex leading-term rule works for L=2..3 and "
                "has a finite capacity counterexample at L=4; the all-L inequality remains open."
            ),
            "rank_scope": (
                "[LEMMA] A nonzero modular minor is a nonzero integer minor, so every "
                "stored modular rank is a rigorous lower bound over Q. Two-prime agreement "
                "is a cross-check, not a proof of rational rank equality."
            ),
            "boundary_lex_order": (
                "rightmost rung first, top then bottom within each rung, local order X>Y>I>Z"
            ),
            "pivot_rule_verified_lengths": [2, 3],
            "pivot_rule_counterexample_length": 4,
            "ladders": ladders,
            "controls": controls,
        },
        "checks": checks,
    }


def write_artifact(artifact: dict[str, object]) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_PATH.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--compact",
        action="store_true",
        help="omit the full support-code lists (ranks and hashes are unchanged)",
    )
    args = parser.parse_args(argv)
    try:
        artifact = make_artifact(include_full_supports=not args.compact)
        write_artifact(artifact)
        failed = [check for check in artifact["checks"] if not check["passed"]]
        if failed:
            print(f"FAIL: {failed}")
            return 1
        print("PASS")
        return 0
    except Exception as exc:  # the standalone experiment must end with PASS or FAIL
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
