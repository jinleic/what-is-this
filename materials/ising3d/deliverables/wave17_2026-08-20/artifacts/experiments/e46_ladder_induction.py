"""Exact search for an all-length two-leg-ladder leading-string induction.

The experiment tests two proposed repairs of the finite leading-term certificate and
produces exact obstruction certificates when they fail.  All Pauli and bracket
arithmetic is over the integers; no floating point is used.
"""

from __future__ import annotations

import argparse
import itertools
import json
import signal
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from ising.clifford import symplectic_form

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_induction.json"
SCRIPT = "experiments/e46_ladder_induction.py"
TIMEOUT_SECONDS = 300
LOCAL_LABELS = "IXYZ"
LOCAL_BITS = {
    "I": (0, 0),
    "X": (1, 0),
    "Z": (0, 1),
    "Y": (1, 1),
}
BITS_LOCAL = {value: key for key, value in LOCAL_BITS.items()}
Vector = dict[int, int]


def ladder_bonds(length: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    if length < 2:
        raise ValueError("the ladder length must be at least two")
    bonds: list[tuple[int, int]] = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return 2 * length, tuple(bonds)


def encode_rungs(rungs: Sequence[str]) -> int:
    """Encode two-character top/bottom rung labels as Q_(a|b)=X^a Z^b."""

    length = len(rungs)
    a = 0
    b = 0
    for col, pair in enumerate(rungs):
        if len(pair) != 2 or any(label not in LOCAL_BITS for label in pair):
            raise ValueError(f"invalid rung label {pair!r}")
        for site, label in ((col, pair[0]), (length + col, pair[1])):
            abit, bbit = LOCAL_BITS[label]
            a |= abit << site
            b |= bbit << site
    return a | (b << (2 * length))


def decode_rungs(pauli: int, length: int) -> str:
    n = 2 * length
    mask = (1 << n) - 1
    a = pauli & mask
    b = (pauli >> n) & mask
    pairs: list[str] = []
    for col in range(length):
        top = BITS_LOCAL[((a >> col) & 1, (b >> col) & 1)]
        bottom = BITS_LOCAL[((a >> (length + col)) & 1, (b >> (length + col)) & 1)]
        pairs.append(top + bottom)
    return " ".join(pairs)


def commutator_sign_half(left: int, right: int, n: int) -> int:
    """Coefficient in [Q_left,Q_right]/2, in {0,+1,-1}."""

    mask = (1 << n) - 1
    first = (((left >> n) & mask) & (right & mask)).bit_count() & 1
    second = (((right >> n) & mask) & (left & mask)).bit_count() & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def commutator_half(left: Vector, right: Vector, n: int) -> Vector:
    """Compute [left,right]/2 exactly in the Q_(a|b) basis."""

    result: Vector = {}
    for lhs, lhs_coefficient in left.items():
        for rhs, rhs_coefficient in right.items():
            sign = commutator_sign_half(lhs, rhs, n)
            if not sign:
                continue
            target = lhs ^ rhs
            value = result.get(target, 0) + sign * lhs_coefficient * rhs_coefficient
            if value:
                result[target] = value
            else:
                result.pop(target, None)
    return result


def generator_vectors(length: int) -> dict[str, Vector]:
    n, bonds = ladder_bonds(length)
    a = {1 << site: 1 for site in range(n)}
    b = {(1 << (n + left)) | (1 << (n + right)): 1 for left, right in bonds}
    return {"A": a, "B": b}


def nonzero_word_vectors(length: int) -> list[tuple[str, Vector]]:
    """Every nonzero generator-left-nested word of depth at most 2L."""

    n = 2 * length
    generators = generator_vectors(length)
    frontier = [(letter, dict(generators[letter])) for letter in "AB"]
    all_words = list(frontier)
    for _depth in range(2, 2 * length + 1):
        next_frontier: list[tuple[str, Vector]] = []
        for outer in "AB":
            for word, vector in frontier:
                child = commutator_half(generators[outer], vector, n)
                if child:
                    next_frontier.append((outer + word, child))
        all_words.extend(next_frontier)
        frontier = next_frontier
    return all_words


def site_labels(pauli: int, length: int) -> tuple[str, ...]:
    n = 2 * length
    mask = (1 << n) - 1
    a = pauli & mask
    b = (pauli >> n) & mask
    return tuple(BITS_LOCAL[((a >> site) & 1, (b >> site) & 1)] for site in range(n))


def monomial_key(
    pauli: int,
    length: int,
    rung_reverse: bool,
    row_reverse: bool,
    alphabet_ascending: str,
    decoded: dict[int, tuple[str, ...]] | None = None,
) -> tuple[int, ...]:
    labels = decoded[pauli] if decoded is not None else site_labels(pauli, length)
    rank = {label: index for index, label in enumerate(alphabet_ascending)}
    columns: Iterable[int] = range(length - 1, -1, -1) if rung_reverse else range(length)
    key: list[int] = []
    for col in columns:
        pair = (labels[col], labels[length + col])
        if row_reverse:
            pair = pair[::-1]
        key.extend(rank[label] for label in pair)
    return tuple(key)


def order_id(rung_reverse: bool, row_reverse: bool, alphabet_ascending: str) -> str:
    rung_direction = "R2L" if rung_reverse else "L2R"
    row_direction = "B2T" if row_reverse else "T2B"
    descending = alphabet_ascending[::-1]
    return f"{rung_direction}:{row_direction}:local_desc={descending}"


def scan_product_lex_orders(
    word_sets: dict[int, list[tuple[str, Vector]]],
) -> dict[str, object]:
    lengths = sorted(word_sets)
    supports = {
        length: sorted({pauli for _, vector in words for pauli in vector})
        for length, words in word_sets.items()
    }
    decoded = {
        length: {pauli: site_labels(pauli, length) for pauli in supports[length]}
        for length in lengths
    }
    order_rows: list[dict[str, object]] = []
    for rung_reverse, row_reverse, alphabet in itertools.product(
        (False, True), (False, True), itertools.permutations(LOCAL_LABELS)
    ):
        alphabet_ascending = "".join(alphabet)
        capacities: dict[str, int] = {}
        for length in lengths:
            key_by_pauli = {
                pauli: monomial_key(
                    pauli,
                    length,
                    rung_reverse,
                    row_reverse,
                    alphabet_ascending,
                    decoded[length],
                )
                for pauli in supports[length]
            }
            leaders = {
                max(vector, key=key_by_pauli.__getitem__) for _, vector in word_sets[length]
            }
            capacities[str(length)] = len(leaders)
        order_rows.append(
            {
                "id": order_id(rung_reverse, row_reverse, alphabet_ascending),
                "rung_reverse": rung_reverse,
                "row_reverse": row_reverse,
                "alphabet_ascending": alphabet_ascending,
                "alphabet_descending": alphabet_ascending[::-1],
                "capacities": capacities,
                "passes_lengths": [
                    length for length in lengths if capacities[str(length)] >= (1 << length)
                ],
            }
        )

    by_length: list[dict[str, object]] = []
    for length in lengths:
        values = [int(row["capacities"][str(length)]) for row in order_rows]
        histogram = Counter(values)
        by_length.append(
            {
                "length": length,
                "depth_limit": 2 * length,
                "nonzero_bracket_words": len(word_sets[length]),
                "distinct_pauli_support_codes": len(supports[length]),
                "required": 1 << length,
                "maximum_capacity": max(values),
                "orders_reaching_maximum": sum(value == max(values) for value in values),
                "capacity_histogram": {str(key): histogram[key] for key in sorted(histogram)},
            }
        )

    passing_all = [
        row["id"]
        for row in order_rows
        if all(int(row["capacities"][str(length)]) >= (1 << length) for length in lengths)
    ]
    return {
        "definition": (
            "All 2 longitudinal directions x 2 within-rung directions x 24 strict local "
            "orders on I,X,Y,Z; the vector leader is the lexicographic maximum."
        ),
        "number_of_orders": len(order_rows),
        "by_length": by_length,
        "orders_passing_every_length_2_to_5": passing_all,
        "orders": order_rows,
    }


def get_order(
    scan: dict[str, object], rung_reverse: bool, row_reverse: bool, alphabet_ascending: str
) -> dict[str, object]:
    target = order_id(rung_reverse, row_reverse, alphabet_ascending)
    return next(row for row in scan["orders"] if row["id"] == target)


def recursive_candidate_words(length: int) -> list[tuple[str, str]]:
    family = [
        ("00", "A"),
        ("01", "BA"),
        ("10", "ABA"),
        ("11", "BABA"),
    ]
    for _ in range(2, length):
        children: list[tuple[str, str]] = []
        for label, word in family:
            children.append((label + "0", "AB" + word))
            children.append((label + "1", "BB" + word))
        family = children
    return family


def vector_for_word(word: str, length: int) -> Vector:
    generators = generator_vectors(length)
    vector = dict(generators[word[-1]])
    for letter in reversed(word[:-1]):
        vector = commutator_half(generators[letter], vector, 2 * length)
    return vector


def exact_pair_minor(first: Vector, second: Vector, length: int) -> dict[str, object]:
    supports = sorted(set(first) | set(second))
    for index, left in enumerate(supports):
        for right in supports[index + 1 :]:
            determinant = (
                first.get(left, 0) * second.get(right, 0)
                - first.get(right, 0) * second.get(left, 0)
            )
            if determinant:
                return {
                    "support_codes": [left, right],
                    "support_paulis": [decode_rungs(left, length), decode_rungs(right, length)],
                    "first_coefficients": [first.get(left, 0), first.get(right, 0)],
                    "second_coefficients": [second.get(left, 0), second.get(right, 0)],
                    "determinant": determinant,
                }
    raise ValueError("the vectors are proportional")


def recursive_collision(
    length: int,
    first_label: str,
    second_label: str,
    rung_reverse: bool,
    row_reverse: bool,
    alphabet_ascending: str,
) -> dict[str, object]:
    family = dict(recursive_candidate_words(length))
    first_word = family[first_label]
    second_word = family[second_label]
    first = vector_for_word(first_word, length)
    second = vector_for_word(second_word, length)
    decoded = {
        pauli: site_labels(pauli, length) for pauli in set(first) | set(second)
    }
    key = lambda pauli: monomial_key(  # noqa: E731
        pauli,
        length,
        rung_reverse,
        row_reverse,
        alphabet_ascending,
        decoded,
    )
    first_leader = max(first, key=key)
    second_leader = max(second, key=key)
    return {
        "length": length,
        "order": order_id(rung_reverse, row_reverse, alphabet_ascending),
        "first_label": first_label,
        "second_label": second_label,
        "first_word": first_word,
        "second_word": second_word,
        "first_leader": first_leader,
        "second_leader": second_leader,
        "leaders_equal": first_leader == second_leader,
        "leader_pauli": decode_rungs(first_leader, length),
        "leader_coefficients": [first[first_leader], second.get(first_leader, 0)],
        "exact_rank_two_minor": exact_pair_minor(first, second, length),
    }


def leader_representatives(
    words: list[tuple[str, Vector]],
    length: int,
    rung_reverse: bool,
    row_reverse: bool,
    alphabet_ascending: str,
) -> list[dict[str, object]]:
    supports = {pauli for _, vector in words for pauli in vector}
    decoded = {pauli: site_labels(pauli, length) for pauli in supports}
    key_by_pauli = {
        pauli: monomial_key(
            pauli,
            length,
            rung_reverse,
            row_reverse,
            alphabet_ascending,
            decoded,
        )
        for pauli in supports
    }
    buckets: dict[int, list[str]] = defaultdict(list)
    for word, vector in words:
        buckets[max(vector, key=key_by_pauli.__getitem__)].append(word)
    rows = []
    for leader in sorted(buckets, key=key_by_pauli.__getitem__):
        representative = min(buckets[leader], key=lambda word: (len(word), word))
        rows.append(
            {
                "leader": leader,
                "leader_pauli": decode_rungs(leader, length),
                "representative_word": representative,
                "depth": len(representative),
                "bucket_size": len(buckets[leader]),
            }
        )
    return rows


def rung_control(length: int, index: int, local_pair: str) -> int:
    rungs = ["II"] * length
    rungs[index] = local_pair
    return encode_rungs(rungs)


def nested_single_output(
    seed: int, controls: Sequence[int], subset_mask: int, n: int
) -> tuple[int, int] | None:
    pauli = seed
    coefficient = 1
    selected = [index for index in range(len(controls)) if (subset_mask >> index) & 1]
    for index in reversed(selected):
        sign = commutator_sign_half(controls[index], pauli, n)
        if not sign:
            return None
        coefficient *= sign
        pauli ^= controls[index]
    return pauli, coefficient


def finite_corner_seed_search() -> dict[str, object]:
    local_pairs = ["".join(pair) for pair in itertools.product(LOCAL_LABELS, repeat=2)]
    controls = [pair for pair in local_pairs if pair != "II"]
    templates = [
        (first, second)
        for first in local_pairs
        for second in local_pairs
        if (first, second) != ("II", "II")
    ]
    by_length: list[dict[str, object]] = []
    injective_family_ids: dict[int, set[str]] = {}
    for length in range(2, 6):
        full_injective: set[str] = set()
        max_distinct_nonzero = 0
        max_distinct_including_zero = 0
        best: tuple[str, str, str] | None = None
        for first, second in templates:
            seed = encode_rungs([first, second] + ["II"] * (length - 2))
            for control in controls:
                rung_controls = [rung_control(length, index, control) for index in range(length)]
                outputs = [
                    nested_single_output(seed, rung_controls, mask, 2 * length)
                    for mask in range(1 << length)
                ]
                nonzero_codes = {output[0] for output in outputs if output is not None}
                states = {None if output is None else output[0] for output in outputs}
                family_id = f"seed={first},{second};control={control}"
                if all(output is not None for output in outputs) and len(nonzero_codes) == (1 << length):
                    full_injective.add(family_id)
                score = (len(nonzero_codes), len(states))
                if score > (max_distinct_nonzero, max_distinct_including_zero):
                    max_distinct_nonzero, max_distinct_including_zero = score
                    best = (first, second, control)
        injective_family_ids[length] = full_injective
        by_length.append(
            {
                "length": length,
                "families_tested": len(templates) * len(controls),
                "fully_injective_families": len(full_injective),
                "maximum_distinct_nonzero_outputs": max_distinct_nonzero,
                "maximum_distinct_outputs_including_zero": max_distinct_including_zero,
                "one_maximizer": {
                    "seed_first_rung": best[0],
                    "seed_second_rung": best[1],
                    "control_rung": best[2],
                }
                if best
                else None,
            }
        )

    common = set.intersection(*(injective_family_ids[length] for length in range(2, 6)))
    canonical: list[dict[str, object]] = []
    for length in range(2, 6):
        seed = encode_rungs(["ZI", "ZI"] + ["II"] * (length - 2))
        controls_for_length = [rung_control(length, index, "XX") for index in range(length)]
        outputs = [
            nested_single_output(seed, controls_for_length, mask, 2 * length)
            for mask in range(1 << length)
        ]
        row: dict[str, object] = {
            "length": length,
            "nonzero_outputs": sum(output is not None for output in outputs),
            "distinct_nonzero_outputs": len(
                {output[0] for output in outputs if output is not None}
            ),
        }
        if length >= 3:
            first_mask = 1 << 2
            second_mask = first_mask | 1
            row["zero_collision"] = {
                "first_subset_mask": first_mask,
                "second_subset_mask": second_mask,
                "first_output": outputs[first_mask],
                "second_output": outputs[second_mask],
            }
        canonical.append(row)
    return {
        "seed_class": (
            "Every nonidentity Pauli seed on the first two rungs (255 templates) and every "
            "nonidentity translated single-rung Pauli control (15 choices)."
        ),
        "by_length": by_length,
        "families_injective_at_every_length_2_to_5": sorted(common),
        "canonical_horizontal_edge_seed": {
            "seed": "ZI ZI (top-row horizontal ZZ edge)",
            "control": "XX on each rung",
            "records": canonical,
        },
    }


def row_swap_code(pauli: int, length: int) -> int:
    n = 2 * length
    site_mask = (1 << n) - 1
    rung_mask = (1 << length) - 1
    a = pauli & site_mask
    b = (pauli >> n) & site_mask

    def swap_half(bits: int) -> int:
        return ((bits & rung_mask) << length) | ((bits >> length) & rung_mask)

    return swap_half(a) | (swap_half(b) << n)


def reflection_code(pauli: int, length: int) -> int:
    n = 2 * length
    site_mask = (1 << n) - 1
    a = pauli & site_mask
    b = (pauli >> n) & site_mask

    def reflect_sites(bits: int) -> int:
        result = 0
        for row in range(2):
            for col in range(length):
                source = row * length + col
                target = row * length + (length - 1 - col)
                result |= ((bits >> source) & 1) << target
        return result

    return reflect_sites(a) | (reflect_sites(b) << n)


def global_repeated_seed_search() -> dict[str, object]:
    pairs = [
        "".join(pair)
        for pair in itertools.product(LOCAL_LABELS, repeat=2)
        if "".join(pair) != "II"
    ]
    successful_by_length: dict[int, set[tuple[str, str]]] = {}
    for length in range(2, 6):
        successful: set[tuple[str, str]] = set()
        for seed_pair, control_pair in itertools.product(pairs, repeat=2):
            seed = encode_rungs([seed_pair] * length)
            controls = [rung_control(length, index, control_pair) for index in range(length)]
            outputs = [
                nested_single_output(seed, controls, mask, 2 * length)
                for mask in range(1 << length)
            ]
            if all(output is not None for output in outputs) and len(
                {output[0] for output in outputs if output is not None}
            ) == (1 << length):
                successful.add((seed_pair, control_pair))
        successful_by_length[length] = successful
    common = set.intersection(*(successful_by_length[length] for length in range(2, 6)))
    row_compatible = {
        pair for pair in common if pair[0][0] == pair[0][1] and pair[1][0] == pair[1][1]
    }

    witness_records: list[dict[str, object]] = []
    for length in range(2, 8):
        seed = encode_rungs(["ZI"] * length)
        controls = [rung_control(length, index, "XX") for index in range(length)]
        outputs = [
            nested_single_output(seed, controls, mask, 2 * length)
            for mask in range(1 << length)
        ]
        off_centre = [
            index
            for index, control in enumerate(controls)
            if reflection_code(control, length) != control
        ]
        witness_records.append(
            {
                "length": length,
                "outputs": len(outputs),
                "nonzero_outputs": sum(output is not None for output in outputs),
                "distinct_supports": len({output[0] for output in outputs if output is not None}),
                "normalized_coefficients": sorted(
                    {output[1] for output in outputs if output is not None}
                ),
                "full_commutator_magnitudes": sorted(
                    {1 << mask.bit_count() for mask in range(1 << length)}
                ),
                "seed_row_swap_invariant": row_swap_code(seed, length) == seed,
                "off_centre_controls_not_reflection_invariant": off_centre,
            }
        )
    return {
        "families_tested": len(pairs) ** 2,
        "injective_pair_counts": {
            str(length): len(successful_by_length[length]) for length in range(2, 6)
        },
        "same_pairs_injective_at_every_length_2_to_5": len(common),
        "same_pairs_passing_row_swap_necessary_condition": len(row_compatible),
        "canonical_combinatorial_witness": {
            "seed": "(ZI)^L",
            "controls": "C_i=XX on rung i",
            "injective_signature": "rung i is ZI if i is not in S and YX if i is in S",
            "records": witness_records,
            "membership_failure": (
                "The seed is not row-swap invariant; every off-centre C_i is not longitudinal-"
                "reflection invariant.  Hence these single strings are not elements of <A_L,B_L>."
            ),
        },
    }


def fixed_single_pauli_census() -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for length in range(2, 6):
        n = 2 * length
        fixed = [
            pauli
            for pauli in range(1 << (2 * n))
            if row_swap_code(pauli, length) == pauli
        ]
        fixed_both = [
            pauli for pauli in fixed if reflection_code(pauli, length) == pauli
        ]
        anticommuting_pairs = 0
        for index, first in enumerate(fixed):
            for second in fixed[index + 1 :]:
                anticommuting_pairs += symplectic_form(first, second, n)
        rows.append(
            {
                "length": length,
                "all_pauli_codes": 1 << (4 * length),
                "row_swap_fixed_codes": len(fixed),
                "expected_row_swap_fixed_codes": 4**length,
                "row_and_reflection_fixed_codes": len(fixed_both),
                "expected_row_and_reflection_fixed_codes": 4 ** ((length + 1) // 2),
                "anticommuting_pairs_among_row_fixed_codes": anticommuting_pairs,
            }
        )
    return {
        "statement": (
            "Every single-Pauli element of <A_L,B_L> must be fixed by row exchange.  Such a "
            "word has the same label on both sites of each rung; every two such words commute."
        ),
        "records": rows,
    }


def orbit_seed(length: int, local_label: str) -> Vector:
    codes = []
    for site in (0, length - 1, length, 2 * length - 1):
        rungs = ["II"] * length
        col = site % length
        pair = list(rungs[col])
        pair[0 if site < length else 1] = local_label
        rungs[col] = "".join(pair)
        codes.append(encode_rungs(rungs))
    return {code: 1 for code in set(codes)}


def orbit_control(length: int, index: int, local_pair: str) -> Vector:
    codes = []
    for col in {index, length - 1 - index}:
        codes.append(rung_control(length, col, local_pair))
        codes.append(rung_control(length, col, local_pair[::-1]))
    return {code: 1 for code in set(codes)}


def nested_vector_output(
    seed: Vector, controls: Sequence[Vector], subset_mask: int, n: int
) -> Vector:
    vector = dict(seed)
    selected = [index for index in range(len(controls)) if (subset_mask >> index) & 1]
    for index in reversed(selected):
        vector = commutator_half(controls[index], vector, n)
        if not vector:
            break
    return vector


def reflect_subset(mask: int, length: int) -> int:
    result = 0
    for index in range(length):
        result |= ((mask >> index) & 1) << (length - 1 - index)
    return result


def frozen_vector(vector: Vector) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(vector.items()))


def orbit_symmetrized_search() -> dict[str, object]:
    control_pairs = [
        "".join(pair)
        for pair in itertools.product(LOCAL_LABELS, repeat=2)
        if "".join(pair) != "II"
    ]
    rows: list[dict[str, object]] = []
    all_reflection_identities = True
    all_controls_commute = True
    for length in range(2, 6):
        maximum_distinct = 0
        maximizer: tuple[str, str] | None = None
        identities_tested = 0
        for seed_label in "XYZ":
            seed = orbit_seed(length, seed_label)
            for control_pair in control_pairs:
                controls = [
                    orbit_control(length, index, control_pair) for index in range(length)
                ]
                for first in range(length):
                    for second in range(first + 1, length):
                        if commutator_half(controls[first], controls[second], 2 * length):
                            all_controls_commute = False
                outputs = {
                    mask: frozen_vector(
                        nested_vector_output(seed, controls, mask, 2 * length)
                    )
                    for mask in range(1 << length)
                }
                for mask, output in outputs.items():
                    identities_tested += 1
                    if output != outputs[reflect_subset(mask, length)]:
                        all_reflection_identities = False
                distinct = len(set(outputs.values()))
                if distinct > maximum_distinct:
                    maximum_distinct = distinct
                    maximizer = (seed_label, control_pair)
        reflection_orbit_bound = ((1 << length) + (1 << ((length + 1) // 2))) // 2
        rows.append(
            {
                "length": length,
                "families_tested": 3 * len(control_pairs),
                "subset_reflection_identities_tested": identities_tested,
                "reflection_orbit_upper_bound": reflection_orbit_bound,
                "maximum_distinct_outputs_found": maximum_distinct,
                "one_maximizer": {
                    "seed_local_label": maximizer[0],
                    "control_local_pair": maximizer[1],
                }
                if maximizer
                else None,
            }
        )

    canonical: list[dict[str, object]] = []
    for length in range(2, 6):
        seed = orbit_seed(length, "X")
        controls = [orbit_control(length, index, "IY") for index in range(length)]
        first = nested_vector_output(seed, controls, 1, 2 * length)
        reflected = nested_vector_output(seed, controls, 1 << (length - 1), 2 * length)
        canonical.append(
            {
                "length": length,
                "first_subset_mask": 1,
                "reflected_subset_mask": 1 << (length - 1),
                "equal": first == reflected,
                "nonzero": bool(first),
                "output": [
                    {
                        "pauli": pauli,
                        "rung_pauli": decode_rungs(pauli, length),
                        "coefficient": coefficient,
                    }
                    for pauli, coefficient in sorted(first.items())
                ],
            }
        )
    return {
        "class": (
            "Corner seed orbit sums and row-plus-reflection orbit sums of one translated "
            "single-rung Pauli control."
        ),
        "all_control_pairs_commute": all_controls_commute,
        "all_w_S_equal_w_reflected_S": all_reflection_identities,
        "by_length": rows,
        "canonical_nonzero_collision": {
            "seed": "orbit sum of X on the four corners",
            "control": "row/reflection orbit sum of IY on rung i",
            "records": canonical,
        },
    }




def make_artifact() -> dict[str, object]:
    word_sets = {length: nonzero_word_vectors(length) for length in range(2, 6)}
    order_scan = scan_product_lex_orders(word_sets)
    old_order = get_order(order_scan, True, False, "ZIYX")
    repaired_order = get_order(order_scan, True, False, "IZYX")
    repaired_representatives = leader_representatives(
        word_sets[4], 4, True, False, "IZYX"
    )
    collision_l4_old = recursive_collision(4, "0101", "0110", True, False, "ZIYX")
    collision_l4_repaired = recursive_collision(4, "0101", "0110", True, False, "IZYX")
    collision_l5_repaired = recursive_collision(5, "01001", "01010", True, False, "IZYX")
    finite_seed = finite_corner_seed_search()
    global_seed = global_repeated_seed_search()
    fixed_single = fixed_single_pauli_census()
    orbit_search = orbit_symmetrized_search()

    data = {
        "status": "OBSTRUCTION_LEMMA",
        "target_all_L_bound": "UNRESOLVED",
        "product_lex_order_scan": order_scan,
        "old_order": old_order,
        "repaired_order": {
            **repaired_order,
            "finite_L4_distinct_leader_representatives": repaired_representatives,
        },
        "named_recursive_collisions": {
            "L4_old_order": collision_l4_old,
            "L4_repaired_order": collision_l4_repaired,
            "L5_repaired_order": collision_l5_repaired,
        },
        "fixed_corner_seed_search": finite_seed,
        "global_repeated_seed_search": global_seed,
        "single_pauli_symmetry_census": fixed_single,
        "orbit_symmetrized_search": orbit_search,
        "proved_obstructions": [
            {
                "tag": "LEMMA",
                "name": "fixed_support_seed",
                "scope": "all L and all single-rung Pauli controls",
                "statement": (
                    "Commuting disjoint rung controls cannot enlarge a Pauli seed's rung support; "
                    "selecting a control on a rung outside the seed support makes the nested "
                    "commutator zero."
                ),
            },
            {
                "tag": "THEOREM",
                "name": "single_pauli_membership_no_go",
                "scope": "the full two-generator ladder Lie algebra for every L>=2",
                "statement": (
                    "Every single-Pauli element is row-swap fixed, and all row-swap-fixed Pauli "
                    "words commute.  Therefore nested commutators of single-Pauli algebra elements "
                    "cannot furnish an injection."
                ),
            },
            {
                "tag": "LEMMA",
                "name": "orbit_reflection_collision",
                "scope": "reflection-covariant commuting depth-one rung orbit controls",
                "statement": (
                    "w_S=w_rho(S), so there are at most "
                    "(2^L+2^ceil(L/2))/2 distinct outputs."
                ),
            },
        ],
        "limitations": [
            "No all-L lower bound for D_(2L) or dim <A_L,B_L> is proved.",
            "The 96-order census covers product-lex orders only, not every abstract total order.",
            "The no-go does not exclude induction with multi-string nonlocal invariant elements.",
        ],
    }

    repaired_caps = {int(key): value for key, value in repaired_order["capacities"].items()}
    fixed_rows = fixed_single["records"]
    orbit_records = orbit_search["canonical_nonzero_collision"]["records"]
    witness_records = global_seed["canonical_combinatorial_witness"]["records"]
    checks = [
        {
            "name": "repaired_order_reaches_required_capacity_through_L4",
            "passed": all(repaired_caps[length] >= (1 << length) for length in range(2, 5)),
            "detail": f"computed capacities {repaired_caps}",
        },
        {
            "name": "repaired_order_fails_at_L5",
            "passed": repaired_caps[5] < (1 << 5),
            "detail": f"computed L=5 capacity {repaired_caps[5]} < 32",
        },
        {
            "name": "no_scanned_product_lex_order_passes_L2_to_L5",
            "passed": not order_scan["orders_passing_every_length_2_to_5"],
            "detail": (
                f"searched {order_scan['number_of_orders']} orders; passing IDs "
                f"{order_scan['orders_passing_every_length_2_to_5']}"
            ),
        },
        {
            "name": "named_leader_collisions_are_exact_and_nonproportional",
            "passed": all(
                collision["leaders_equal"]
                and collision["exact_rank_two_minor"]["determinant"] != 0
                for collision in (collision_l4_old, collision_l4_repaired, collision_l5_repaired)
            ),
            "detail": (
                "leader codes "
                f"{collision_l4_old['first_leader']}, {collision_l4_repaired['first_leader']}, "
                f"{collision_l5_repaired['first_leader']}; exact minor determinants "
                f"{collision_l4_old['exact_rank_two_minor']['determinant']}, "
                f"{collision_l4_repaired['exact_rank_two_minor']['determinant']}, "
                f"{collision_l5_repaired['exact_rank_two_minor']['determinant']}"
            ),
        },
        {
            "name": "no_fixed_two_rung_seed_family_survives_L2_to_L5",
            "passed": not finite_seed["families_injective_at_every_length_2_to_5"],
            "detail": (
                "common injective family IDs: "
                f"{finite_seed['families_injective_at_every_length_2_to_5']}"
            ),
        },
        {
            "name": "global_single_string_signature_is_combinatorially_injective",
            "passed": all(
                row["nonzero_outputs"] == (1 << row["length"])
                and row["distinct_supports"] == (1 << row["length"])
                for row in witness_records
            ),
            "detail": "exactly 2^L nonzero distinct supports for L=2,...,7",
        },
        {
            "name": "global_single_string_witness_fails_membership_symmetries",
            "passed": all(
                not row["seed_row_swap_invariant"]
                and bool(row["off_centre_controls_not_reflection_invariant"])
                for row in witness_records
            ),
            "detail": "seed fails row exchange and at least one required control fails reflection",
        },
        {
            "name": "row_fixed_single_paulis_are_exhaustively_commuting_L2_to_L5",
            "passed": all(
                row["row_swap_fixed_codes"] == row["expected_row_swap_fixed_codes"]
                and row["row_and_reflection_fixed_codes"]
                == row["expected_row_and_reflection_fixed_codes"]
                and row["anticommuting_pairs_among_row_fixed_codes"] == 0
                for row in fixed_rows
            ),
            "detail": (
                "row-fixed counts "
                f"{[row['row_swap_fixed_codes'] for row in fixed_rows]}; anticommuting counts "
                f"{[row['anticommuting_pairs_among_row_fixed_codes'] for row in fixed_rows]}"
            ),
        },
        {
            "name": "orbit_symmetrization_has_reflection_collisions",
            "passed": bool(orbit_search["all_control_pairs_commute"])
            and bool(orbit_search["all_w_S_equal_w_reflected_S"])
            and all(row["equal"] and row["nonzero"] for row in orbit_records),
            "detail": "all enumerated identities hold; canonical singleton/reflected-singleton outputs are nonzero",
        },
    ]
    return {
        "meta": {
            "title": "Two-leg ladder all-length induction obstruction",
            "schema_version": 1,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "provenance": SCRIPT,
            "arithmetic": "exact integer Pauli coefficients and GF(2) symplectic parity",
            "timeout_seconds": TIMEOUT_SECONDS,
        },
        "data": data,
        "checks": checks,
    }


def write_artifact(artifact: dict[str, object]) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_PATH.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded the {TIMEOUT_SECONDS}s hard timeout")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true", help="compute without writing JSON")
    args = parser.parse_args(argv)
    signal.signal(signal.SIGALRM, _timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        artifact = make_artifact()
        if not args.no_write:
            write_artifact(artifact)
        for check in artifact["checks"]:
            print(f"{'PASS' if check['passed'] else 'FAIL'}: {check['name']} — {check['detail']}")
        passed = all(check["passed"] for check in artifact["checks"])
        print("PASS" if passed else "FAIL")
        return 0 if passed else 1
    except TimeoutError as exc:
        print(str(exc))
        return 2
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
