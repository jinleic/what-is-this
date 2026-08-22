"""Standalone exact checks for the two-leg ladder induction obstruction.

This test deliberately reimplements the Pauli bookkeeping instead of importing the
experiment.  The JSON artifact is consulted only after the mathematical certificates
have been re-derived.
"""

from __future__ import annotations

import itertools
import json
import signal
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_growth" / "ladder_induction.json"
TIMEOUT_SECONDS = 300
LABELS = "IXYZ"
TO_BITS = {"I": (0, 0), "X": (1, 0), "Z": (0, 1), "Y": (1, 1)}
FROM_BITS = {value: key for key, value in TO_BITS.items()}


def encode(rungs: list[str]) -> int:
    length = len(rungs)
    a = b = 0
    for col, pair in enumerate(rungs):
        for site, label in ((col, pair[0]), (length + col, pair[1])):
            abit, bbit = TO_BITS[label]
            a |= abit << site
            b |= bbit << site
    return a | (b << (2 * length))


def decode(pauli: int, length: int) -> str:
    n = 2 * length
    mask = (1 << n) - 1
    a, b = pauli & mask, (pauli >> n) & mask
    return " ".join(
        FROM_BITS[((a >> col) & 1, (b >> col) & 1)]
        + FROM_BITS[((a >> (length + col)) & 1, (b >> (length + col)) & 1)]
        for col in range(length)
    )


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    lr = (((left >> n) & mask) & (right & mask)).bit_count() & 1
    rl = (((right >> n) & mask) & (left & mask)).bit_count() & 1
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


def words_through_2l(length: int) -> list[tuple[str, dict[int, int]]]:
    gens = generators(length)
    frontier = [(letter, dict(gens[letter])) for letter in "AB"]
    rows = list(frontier)
    for _ in range(2, 2 * length + 1):
        children = []
        for outer in "AB":
            for word, vector in frontier:
                value = bracket(gens[outer], vector, 2 * length)
                if value:
                    children.append((outer + word, value))
        rows.extend(children)
        frontier = children
    return rows


def site_labels(pauli: int, length: int) -> tuple[str, ...]:
    n = 2 * length
    mask = (1 << n) - 1
    a, b = pauli & mask, (pauli >> n) & mask
    return tuple(FROM_BITS[((a >> site) & 1, (b >> site) & 1)] for site in range(n))


def order_key(
    pauli: int,
    length: int,
    rung_reverse: bool,
    row_reverse: bool,
    alphabet: str,
    decoded: dict[int, tuple[str, ...]],
) -> tuple[int, ...]:
    rank = {label: index for index, label in enumerate(alphabet)}
    columns = range(length - 1, -1, -1) if rung_reverse else range(length)
    key = []
    for col in columns:
        pair = [decoded[pauli][col], decoded[pauli][length + col]]
        if row_reverse:
            pair.reverse()
        key.extend(rank[label] for label in pair)
    return tuple(key)


def order_capacities(
    word_sets: dict[int, list[tuple[str, dict[int, int]]]],
) -> dict[tuple[bool, bool, str], dict[int, int]]:
    result = {}
    for rung_reverse, row_reverse, permutation in itertools.product(
        (False, True), (False, True), itertools.permutations(LABELS)
    ):
        alphabet = "".join(permutation)
        capacities = {}
        for length, words in word_sets.items():
            supports = {pauli for _, vector in words for pauli in vector}
            decoded = {pauli: site_labels(pauli, length) for pauli in supports}
            keys = {
                pauli: order_key(
                    pauli, length, rung_reverse, row_reverse, alphabet, decoded
                )
                for pauli in supports
            }
            capacities[length] = len(
                {max(vector, key=keys.__getitem__) for _, vector in words}
            )
        result[(rung_reverse, row_reverse, alphabet)] = capacities
    return result


def recursive_words(length: int) -> dict[str, str]:
    family = {"00": "A", "01": "BA", "10": "ABA", "11": "BABA"}
    for _ in range(2, length):
        family = {
            child_label: child_word
            for label, word in family.items()
            for child_label, child_word in (
                (label + "0", "AB" + word),
                (label + "1", "BB" + word),
            )
        }
    return family


def word_vector(word: str, length: int) -> dict[int, int]:
    gens = generators(length)
    value = dict(gens[word[-1]])
    for letter in reversed(word[:-1]):
        value = bracket(gens[letter], value, 2 * length)
    return value


def exact_minor(first: dict[int, int], second: dict[int, int]) -> int:
    supports = sorted(set(first) | set(second))
    for index, lhs in enumerate(supports):
        for rhs in supports[index + 1 :]:
            determinant = (
                first.get(lhs, 0) * second.get(rhs, 0)
                - first.get(rhs, 0) * second.get(lhs, 0)
            )
            if determinant:
                return determinant
    return 0


def control(length: int, index: int, pair: str) -> int:
    rungs = ["II"] * length
    rungs[index] = pair
    return encode(rungs)


def nested_single(
    seed: int, controls: list[int], subset: int, n: int
) -> tuple[int, int] | None:
    value, coefficient = seed, 1
    for index in reversed(range(len(controls))):
        if (subset >> index) & 1:
            sign = sign_half(controls[index], value, n)
            if not sign:
                return None
            value ^= controls[index]
            coefficient *= sign
    return value, coefficient


def row_swap(pauli: int, length: int) -> int:
    n = 2 * length
    mask_n = (1 << n) - 1
    mask_l = (1 << length) - 1
    a, b = pauli & mask_n, (pauli >> n) & mask_n

    def swap(bits: int) -> int:
        return ((bits & mask_l) << length) | ((bits >> length) & mask_l)

    return swap(a) | (swap(b) << n)


def reflect(pauli: int, length: int) -> int:
    n = 2 * length
    mask = (1 << n) - 1
    a, b = pauli & mask, (pauli >> n) & mask

    def one(bits: int) -> int:
        out = 0
        for row in range(2):
            for col in range(length):
                source = row * length + col
                target = row * length + length - 1 - col
                out |= ((bits >> source) & 1) << target
        return out

    return one(a) | (one(b) << n)


def symplectic(first: int, second: int, n: int) -> int:
    mask = (1 << n) - 1
    af, bf = first & mask, (first >> n) & mask
    as_, bs = second & mask, (second >> n) & mask
    return ((bf & as_).bit_count() + (bs & af).bit_count()) & 1


def orbit_seed(length: int, label: str) -> dict[int, int]:
    codes = set()
    for site in (0, length - 1, length, 2 * length - 1):
        rungs = ["II"] * length
        col = site % length
        pair = list(rungs[col])
        pair[0 if site < length else 1] = label
        rungs[col] = "".join(pair)
        codes.add(encode(rungs))
    return {code: 1 for code in codes}


def orbit_control(length: int, index: int, pair: str) -> dict[int, int]:
    codes = set()
    for col in {index, length - 1 - index}:
        codes.add(control(length, col, pair))
        codes.add(control(length, col, pair[::-1]))
    return {code: 1 for code in codes}


def nested_vector(
    seed: dict[int, int], controls: list[dict[int, int]], subset: int, n: int
) -> dict[int, int]:
    value = dict(seed)
    for index in reversed(range(len(controls))):
        if (subset >> index) & 1:
            value = bracket(controls[index], value, n)
            if not value:
                break
    return value


def reflected_subset(mask: int, length: int) -> int:
    return sum(((mask >> index) & 1) << (length - 1 - index) for index in range(length))


def run_check(name: str, function) -> bool:
    try:
        detail = function()
        print(f"PASS: {name} — {detail}")
        return True
    except Exception as exc:
        print(f"FAIL: {name} — {type(exc).__name__}: {exc}")
        return False


def check_order_census() -> str:
    word_sets = {length: words_through_2l(length) for length in range(2, 6)}
    assert [len(word_sets[length]) for length in range(2, 6)] == [16, 64, 256, 1024]
    census = order_capacities(word_sets)
    repaired = census[(True, False, "IZYX")]
    old = census[(True, False, "ZIYX")]
    maxima = {
        length: max(capacities[length] for capacities in census.values())
        for length in range(2, 6)
    }
    assert repaired == {2: 6, 3: 9, 4: 16, 5: 26}, repaired
    assert old == {2: 5, 3: 9, 4: 15, 5: 27}, old
    assert maxima == {2: 6, 3: 10, 4: 16, 5: 27}, maxima
    passing = [
        order
        for order, capacities in census.items()
        if all(capacities[length] >= (1 << length) for length in range(2, 6))
    ]
    assert not passing, passing
    return f"96 exact orders; repaired {repaired}; maxima {maxima}; no common pass"


def check_named_collisions() -> str:
    family4 = recursive_words(4)
    first4 = word_vector(family4["0101"], 4)
    second4 = word_vector(family4["0110"], 4)
    supports4 = set(first4) | set(second4)
    decoded4 = {pauli: site_labels(pauli, 4) for pauli in supports4}
    for alphabet in ("ZIYX", "IZYX"):
        keys = {
            pauli: order_key(pauli, 4, True, False, alphabet, decoded4)
            for pauli in supports4
        }
        assert max(first4, key=keys.__getitem__) == max(second4, key=keys.__getitem__) == 33928
    assert decode(33928, 4) == "II II ZI XY"
    assert exact_minor(first4, second4) == -120

    family5 = recursive_words(5)
    first5 = word_vector(family5["01001"], 5)
    second5 = word_vector(family5["01010"], 5)
    supports5 = set(first5) | set(second5)
    decoded5 = {pauli: site_labels(pauli, 5) for pauli in supports5}
    keys5 = {
        pauli: order_key(pauli, 5, True, False, "IZYX", decoded5)
        for pauli in supports5
    }
    assert max(first5, key=keys5.__getitem__) == max(second5, key=keys5.__getitem__) == 270872
    assert decode(270872, 5) == "II II II YZ XX"
    assert exact_minor(first5, second5) == -2268
    return "L4 leader 33928/minor -120; L5 leader 270872/minor -2268"


def check_fixed_seed_search() -> str:
    pairs = ["".join(pair) for pair in itertools.product(LABELS, repeat=2)]
    templates = [(p, q) for p in pairs for q in pairs if (p, q) != ("II", "II")]
    controls_local = [pair for pair in pairs if pair != "II"]
    counts = {}
    for length in range(2, 6):
        injective = 0
        max_nonzero = 0
        for first, second in templates:
            seed = encode([first, second] + ["II"] * (length - 2))
            for local in controls_local:
                controls = [control(length, index, local) for index in range(length)]
                outputs = [
                    nested_single(seed, controls, subset, 2 * length)
                    for subset in range(1 << length)
                ]
                nonzero = {output[0] for output in outputs if output is not None}
                max_nonzero = max(max_nonzero, len(nonzero))
                if all(output is not None for output in outputs) and len(nonzero) == (1 << length):
                    injective += 1
        counts[length] = (injective, max_nonzero)
    assert counts == {2: (960, 4), 3: (0, 4), 4: (0, 4), 5: (0, 4)}, counts
    return f"3825 families/length; (injective,max-nonzero)={counts}"


def check_global_witness_and_membership() -> str:
    for length in range(2, 6):
        seed = encode(["ZI"] * length)
        controls = [control(length, index, "XX") for index in range(length)]
        outputs = [
            nested_single(seed, controls, subset, 2 * length)
            for subset in range(1 << length)
        ]
        assert all(output is not None for output in outputs)
        assert len({output[0] for output in outputs if output is not None}) == (1 << length)
        assert row_swap(seed, length) != seed
        assert any(reflect(item, length) != item for item in controls)
        for subset, output in enumerate(outputs):
            expected = ["YX" if (subset >> index) & 1 else "ZI" for index in range(length)]
            assert decode(output[0], length) == " ".join(expected)
    return "(ZI)^L and C_i=XX gives 2^L signatures, but violates ladder symmetries"


def check_single_pauli_no_go() -> str:
    records = []
    for length in range(2, 6):
        n = 2 * length
        fixed = [
            pauli for pauli in range(1 << (4 * length)) if row_swap(pauli, length) == pauli
        ]
        fixed_both = [pauli for pauli in fixed if reflect(pauli, length) == pauli]
        assert len(fixed) == 4**length
        assert len(fixed_both) == 4 ** ((length + 1) // 2)
        assert all(
            symplectic(first, second, n) == 0
            for index, first in enumerate(fixed)
            for second in fixed[index + 1 :]
        )
        records.append((len(fixed), len(fixed_both)))
    return f"(row-fixed,both-fixed) counts {records}; every pair commutes"


def check_orbit_collision() -> str:
    control_pairs = [
        "".join(pair)
        for pair in itertools.product(LABELS, repeat=2)
        if "".join(pair) != "II"
    ]
    identity_counts = []
    for length in range(2, 6):
        identities = 0
        for seed_label in "XYZ":
            seed = orbit_seed(length, seed_label)
            for pair in control_pairs:
                controls = [orbit_control(length, index, pair) for index in range(length)]
                assert all(
                    not bracket(controls[first], controls[second], 2 * length)
                    for first in range(length)
                    for second in range(first + 1, length)
                )
                outputs = {
                    subset: tuple(
                        sorted(nested_vector(seed, controls, subset, 2 * length).items())
                    )
                    for subset in range(1 << length)
                }
                for subset, output in outputs.items():
                    assert output == outputs[reflected_subset(subset, length)]
                    identities += 1
        seed = orbit_seed(length, "X")
        controls = [orbit_control(length, index, "IY") for index in range(length)]
        first = nested_vector(seed, controls, 1, 2 * length)
        last = nested_vector(seed, controls, 1 << (length - 1), 2 * length)
        assert first == last and first
        identity_counts.append(identities)
    assert identity_counts == [180, 360, 720, 1440]
    return f"exact reflection identities {identity_counts}; nonzero endpoint collision each L"


def check_artifact() -> str:
    with RESULT.open(encoding="utf-8") as handle:
        artifact = json.load(handle)
    assert set(artifact) == {"meta", "data", "checks"}
    assert artifact["meta"]["provenance"] == "experiments/e46_ladder_induction.py"
    assert artifact["data"]["status"] == "OBSTRUCTION_LEMMA"
    assert artifact["data"]["target_all_L_bound"] == "UNRESOLVED"
    assert artifact["checks"] and all(check["passed"] for check in artifact["checks"])
    return f"{len(artifact['checks'])} stored checks pass with honest UNRESOLVED target"


def timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s")


def main() -> int:
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        checks = [
            ("product-lex order census", check_order_census),
            ("named collision certificates", check_named_collisions),
            ("fixed corner seed search", check_fixed_seed_search),
            ("global signature membership failure", check_global_witness_and_membership),
            ("single-Pauli symmetry no-go", check_single_pauli_no_go),
            ("orbit-symmetrized reflection collision", check_orbit_collision),
            ("result artifact", check_artifact),
        ]
        passed = [run_check(name, function) for name, function in checks]
        print("PASS" if all(passed) else "FAIL")
        return 0 if all(passed) else 1
    except TimeoutError as exc:
        print(f"FAIL: {exc}")
        return 2
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
