"""Clean-room verifier for the sixth interlayer G/U/W correlation identity.

The test deliberately imports neither experiments.e99_c6_closed nor its
helpers.  It independently enumerates Bell(6) partitions, rebuilds the 13
connected-correlation normal form, evaluates the finite-volume algebra on a
2x2 layer, and independently re-expands the raw partition formula by exact
2D boundary-mask transfer plus planar finite-lattice inversion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "c6_closed.json"
PRIMARY = ROOT / "results" / "interlayer" / "c6_series.json"
ORDER = 12
MAX_LAYERS = 4

IntSeries = tuple[int, ...]
RatSeries = tuple[Fraction, ...]
Atom = tuple[str, tuple[int, ...]]
Mono = tuple[Atom, ...]
Poly = dict[Mono, Fraction]


# ---------------------------------------------------------------------------
# Exact truncated series and independent 2D transfer
# ---------------------------------------------------------------------------
def izero(order: int = ORDER) -> list[int]:
    return [0] * (order + 1)


def rzero(order: int = ORDER) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def sadd(dst: list[int | Fraction], src: tuple[int | Fraction, ...] | list[int | Fraction], scale: int | Fraction = 1) -> None:
    for d, value in enumerate(src):
        if value:
            dst[d] += scale * value


def smul(a: tuple[int | Fraction, ...] | list[int | Fraction], b: tuple[int | Fraction, ...] | list[int | Fraction], order: int = ORDER) -> tuple[int | Fraction, ...]:
    out: list[int | Fraction] = [0] * (order + 1)
    for i, ai in enumerate(a[: order + 1]):
        if not ai:
            continue
        for j, bj in enumerate(b[: order + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return tuple(out)


def spow(a: tuple[int | Fraction, ...], n: int, order: int = ORDER) -> tuple[int | Fraction, ...]:
    answer: tuple[int | Fraction, ...] = (1,) + (0,) * order
    for _ in range(n):
        answer = smul(answer, a, order)
    return answer


def unit_divide(numer: tuple[int, ...], denom: tuple[int, ...], order: int = ORDER) -> IntSeries:
    assert denom[0] == 1
    result = izero(order)
    for d in range(order + 1):
        value = numer[d]
        for previous in range(d):
            value -= result[previous] * denom[d - previous]
        result[d] = value
    return tuple(result)


def bonds(a: int, b: int) -> tuple[tuple[int, int], ...]:
    answer: list[tuple[int, int]] = []
    for y in range(b):
        for x in range(a):
            p = y * a + x
            if x + 1 < a:
                answer.append((p, p + 1))
            if y + 1 < b:
                answer.append((p, p + a))
    return tuple(answer)


@cache
def parity_masks(a: int, b: int, order: int = ORDER) -> tuple[IntSeries, ...]:
    """Independent mask transfer for P_S(v), retaining exact integers only."""
    table = [[0] * (order + 1) for _ in range(1 << (a * b))]
    table[0][0] = 1
    for p, q in bonds(a, b):
        flip = (1 << p) | (1 << q)
        for d in reversed(range(order)):
            for mask in range(len(table)):
                if table[mask][d]:
                    table[mask ^ flip][d + 1] += table[mask][d]
    return tuple(tuple(row) for row in table)


@cache
def moments(a: int, b: int, order: int = ORDER) -> tuple[IntSeries, ...]:
    poly = parity_masks(a, b, order)
    return tuple(unit_divide(poly[mask], poly[0], order) for mask in range(len(poly)))


# ---------------------------------------------------------------------------
# Independent Bell(6) algebra and closed G/U/W normal form
# ---------------------------------------------------------------------------
def partitions(labels: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Generate unlabeled set partitions by placing labels in growing blocks."""
    states: list[tuple[tuple[int, ...], ...]] = [()]
    for label in labels:
        next_states: list[tuple[tuple[int, ...], ...]] = []
        for state in states:
            next_states.append(state + ((label,),))
            for index in range(len(state)):
                amended = list(state)
                amended[index] = amended[index] + (label,)
                next_states.append(tuple(amended))
        states = next_states
    return tuple(states)


P6 = partitions((1, 2, 3, 4, 5, 6))


def mobius(partition: tuple[tuple[int, ...], ...]) -> int:
    return (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)


def at(kind: str, xs: tuple[int, ...] | list[int]) -> Atom:
    return kind, tuple(sorted(xs))


def mono(*xs: Atom) -> Mono:
    return tuple(sorted(xs))


def pmul(a: Poly, b: Poly) -> Poly:
    out: defaultdict[Mono, Fraction] = defaultdict(Fraction)
    for left, lc in a.items():
        for right, rc in b.items():
            out[mono(*left, *right)] += lc * rc
    return {key: value for key, value in out.items() if value}


def padd(a: Poly, b: Poly, scale: Fraction) -> Poly:
    out: defaultdict[Mono, Fraction] = defaultdict(Fraction, a)
    for key, value in b.items():
        out[key] += scale * value
    return {key: value for key, value in out.items() if value}


def pairings(xs: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], ...]:
    if not xs:
        return ((),)
    first = xs[0]
    out: list[tuple[tuple[int, int], ...]] = []
    for j in range(1, len(xs)):
        rest = xs[1:j] + xs[j + 1 :]
        for tail in pairings(rest):
            out.append(((first, xs[j]),) + tail)
    return tuple(out)


@cache
def m_to_connected(xs: tuple[int, ...]) -> Poly:
    """M_X in terms of G=K2, U=K4, W=K6, independently derived."""
    n = len(xs)
    if n == 0:
        return {(): Fraction(1)}
    if n % 2:
        return {}
    if n == 2:
        return {mono(at("G", xs)): Fraction(1)}
    if n == 4:
        answer: Poly = {mono(at("U", xs)): Fraction(1)}
        for grouping in pairings(xs):
            item = mono(*(at("G", pair) for pair in grouping))
            answer[item] = answer.get(item, Fraction(0)) + 1
        return answer
    if n == 6:
        answer = {mono(at("W", xs)): Fraction(1)}
        for pair in combinations(xs, 2):
            rest = tuple(x for x in xs if x not in pair)
            item = mono(at("G", tuple(pair)), at("U", rest))
            answer[item] = answer.get(item, Fraction(0)) + 1
        for grouping in pairings(xs):
            item = mono(*(at("G", pair) for pair in grouping))
            answer[item] = answer.get(item, Fraction(0)) + 1
        return answer
    raise AssertionError(n)


GAPS = {
    "one": {i: 0 for i in range(1, 7)},
    "two": {1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 1},
    "three": {1: 0, 2: 0, 3: 1, 4: 1, 5: 2, 6: 2},
}
COEFF = {"one": Fraction(1, 720), "two": Fraction(1, 24), "three": Fraction(1, 8)}


def layer_slot_sets(kind: str, block: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    chosen = set(block)
    labels = GAPS[kind]
    if kind == "one":
        one = tuple(sorted(chosen))
        return one, one
    a = tuple(sorted(x for x in chosen if labels[x] == 0))
    b = tuple(sorted(x for x in chosen if labels[x] == 1))
    if kind == "two":
        return a, b, tuple(sorted(a + b))
    c = tuple(sorted(x for x in chosen if labels[x] == 2))
    return a, c, tuple(sorted(a + b)), tuple(sorted(b + c))


def block_poly(kind: str, block: tuple[int, ...]) -> Poly:
    answer: Poly = {(): Fraction(1)}
    for slot_set in layer_slot_sets(kind, block):
        answer = pmul(answer, m_to_connected(slot_set))
        if not answer:
            break
    return answer


def raw_poly(kind: str) -> tuple[Poly, int]:
    answer: Poly = {}
    allowed = 0
    for partition in P6:
        current: Poly = {(): Fraction(1)}
        for block in partition:
            current = pmul(current, block_poly(kind, block))
            if not current:
                break
        if current:
            allowed += 1
            answer = padd(answer, current, COEFF[kind] * mobius(partition))
    return answer, allowed


PERMS = tuple(
    {1: 1, **dict(zip((2, 3, 4, 5, 6), perm))}
    for perm in permutations((2, 3, 4, 5, 6))
)


def remap(item: Mono, mapping: dict[int, int]) -> Mono:
    return mono(*(at(kind, tuple(mapping[x] for x in xs)) for kind, xs in item))


def canonical(item: Mono) -> Mono:
    return min(remap(item, mapping) for mapping in PERMS)


def independently_derived_formula() -> tuple[dict[Mono, Fraction], dict[str, int]]:
    merged: Poly = {}
    counts: dict[str, int] = {}
    for kind in ("one", "two", "three"):
        part, count = raw_poly(kind)
        merged = padd(merged, part, Fraction(1))
        counts[kind] = count
    answer: defaultdict[Mono, Fraction] = defaultdict(Fraction)
    for item, value in merged.items():
        answer[canonical(item)] += value
    return {item: value for item, value in answer.items() if value}, counts


# ---------------------------------------------------------------------------
# Clean-room exact q-cumulant FLM, using raw 2D layer correlations
# ---------------------------------------------------------------------------
@cache
def profile_weights(n: int) -> tuple[tuple[tuple[int, ...], int], ...]:
    values: Counter[tuple[int, ...]] = Counter()
    for partition in partitions(tuple(range(n))):
        signature = tuple(sorted((len(group) for group in partition), reverse=True))
        values[signature] += mobius(partition)
    return tuple(sorted(values.items(), reverse=True))


@cache
def qweights(vertical_bonds: int) -> tuple[tuple[Fraction, ...], ...]:
    cosh = (Fraction(1), Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1, 24), Fraction(0), Fraction(1, 720))
    tanh = (Fraction(0), Fraction(1), Fraction(0), Fraction(-1, 3), Fraction(0), Fraction(2, 15), Fraction(0))
    prefactor = spow(cosh, vertical_bonds, 6)
    answer: list[tuple[Fraction, ...]] = []
    for n in range(7):
        row: list[Fraction] = []
        for d in range(7):
            row.append(Fraction(smul(prefactor, spow(tanh, d, 6), 6)[n]) * factorial(n))
        answer.append(tuple(row))
    return tuple(answer)


def slab_columns(poly: tuple[IntSeries, ...], height: int, order: int = ORDER) -> tuple[IntSeries, ...]:
    active = tuple((mask, mask.bit_count()) for mask, value in enumerate(poly) if mask.bit_count() <= 6 and any(value))
    assert all(weight % 2 == 0 for _mask, weight in active)
    if height == 1:
        out = [tuple(izero(order)) for _ in range(7)]
        out[0] = poly[0]
        return tuple(out)
    state: dict[tuple[int, int], IntSeries] = {(mask, weight): poly[mask] for mask, weight in active}
    for _ in range(height - 2):
        update: dict[tuple[int, int], list[int]] = {}
        for (old, old_degree), carried in state.items():
            for new, new_degree in active:
                degree = old_degree + new_degree
                if degree > 6:
                    continue
                term = smul(carried, poly[old ^ new], order)
                if not any(term):
                    continue
                bucket = update.setdefault((new, degree), izero(order))
                sadd(bucket, term)
        state = {key: tuple(value) for key, value in update.items()}
    result = [izero(order) for _ in range(7)]
    for (last, degree), carried in state.items():
        sadd(result[degree], smul(carried, poly[last], order))
    return tuple(tuple(value) for value in result)


def raw_layer_moments(columns: tuple[IntSeries, ...], planar_partition: IntSeries, height: int, area: int, order: int = ORDER) -> dict[int, RatSeries]:
    denom = spow(planar_partition, height, order)
    normalized = tuple(tuple(Fraction(x) for x in unit_divide(column, denom, order)) for column in columns)
    transform = qweights(area * (height - 1))
    answer: dict[int, RatSeries] = {
        n: tuple(rzero(order)) for n in range(1, 7)
    }
    for n in (2, 4, 6):
        total = rzero(order)
        for d in range(7):
            if transform[n][d]:
                sadd(total, normalized[d], transform[n][d])
        answer[n] = tuple(total)
    return answer


def kth_cumulant(raw: dict[int, RatSeries], n: int, order: int = ORDER) -> RatSeries:
    total = rzero(order)
    for profile, sign in profile_weights(n):
        term: tuple[int | Fraction, ...] = (Fraction(1),) + (Fraction(0),) * order
        for size in profile:
            term = smul(term, raw[size], order)
        sadd(total, term, sign)
    return tuple(total)


@cache
def finite_rectangle(a: int, b: int, order: int = ORDER) -> tuple[tuple[str, RatSeries], ...]:
    a, b = sorted((a, b))
    planar = parity_masks(a, b, order)
    area = a * b
    ordinary: dict[int, dict[int, RatSeries]] = {1: {n: tuple(rzero(order)) for n in (2, 4, 6)}}
    for h in range(2, 5):
        raw = raw_layer_moments(slab_columns(planar, h, order), planar[0], h, area, order)
        ordinary[h] = {n: tuple(value / factorial(n) for value in kth_cumulant(raw, n, order)) for n in (2, 4, 6)}
    linked: dict[tuple[int, int], RatSeries] = {}
    for n in (2, 4, 6):
        for h in range(1, 5):
            value = list(ordinary[h][n])
            for subheight in range(1, h):
                sadd(value, linked[(n, subheight)], -(h - subheight + 1))
            linked[(n, h)] = tuple(value)
    return tuple(sorted((f"c{n}_h{h}", linked[(n, h)]) for n in (2, 4, 6) for h in range(1, 5)))


def shapes(order: int = ORDER) -> tuple[tuple[int, int], ...]:
    span = order // 2
    return tuple(sorted(((a, b) for a in range(1, span + 2) for b in range(1, span + 2) if a + b <= span + 2), key=lambda x: (sum(x), x)))


def independently_reexpanded_cumulants(order: int = ORDER) -> dict[str, RatSeries]:
    names = tuple(f"c{n}_h{h}" for n in (2, 4, 6) for h in range(1, 5))
    weights: dict[tuple[int, int], dict[str, RatSeries]] = {}
    total = {name: rzero(order) for name in names}
    for a, b in shapes(order):
        raw = dict(finite_rectangle(a, b, order))
        local = {name: list(raw[name]) for name in names}
        for (sa, sb), prior in weights.items():
            if sa <= a and sb <= b:
                copies = (a - sa + 1) * (b - sb + 1)
                for name in names:
                    sadd(local[name], prior[name], -copies)
        exact = {name: tuple(value) for name, value in local.items()}
        weights[(a, b)] = exact
        for name in names:
            sadd(total[name], exact[name])
    return {name: tuple(value) for name, value in total.items()}


# ---------------------------------------------------------------------------
# Direct finite-rectangle evaluation of the claimed 13-term identity
# ---------------------------------------------------------------------------
def mask_moment(moment_rows: tuple[IntSeries, ...], slots: tuple[int, ...], pos: dict[int, int]) -> IntSeries:
    mask = 0
    for slot in slots:
        mask ^= 1 << pos[slot]
    return moment_rows[mask]


def g_value(moment_rows: tuple[IntSeries, ...], slots: tuple[int, ...], pos: dict[int, int]) -> IntSeries:
    return mask_moment(moment_rows, slots, pos)


def u_value(moment_rows: tuple[IntSeries, ...], slots: tuple[int, ...], pos: dict[int, int]) -> tuple[int, ...]:
    answer = list(mask_moment(moment_rows, slots, pos))
    for grouping in pairings(slots):
        sadd(answer, smul(g_value(moment_rows, grouping[0], pos), g_value(moment_rows, grouping[1], pos)), -1)
    return tuple(answer)


def w_value(moment_rows: tuple[IntSeries, ...], slots: tuple[int, ...], pos: dict[int, int]) -> tuple[int, ...]:
    answer = list(mask_moment(moment_rows, slots, pos))
    for pair in combinations(slots, 2):
        rest = tuple(slot for slot in slots if slot not in pair)
        sadd(answer, smul(g_value(moment_rows, tuple(pair), pos), u_value(moment_rows, rest, pos)), -1)
    for grouping in pairings(slots):
        term = smul(smul(g_value(moment_rows, grouping[0], pos), g_value(moment_rows, grouping[1], pos)), g_value(moment_rows, grouping[2], pos))
        sadd(answer, term, -1)
    return tuple(answer)


def atom_value(moment_rows: tuple[IntSeries, ...], item: Atom, pos: dict[int, int]) -> tuple[int, ...]:
    kind, slots = item
    if kind == "G":
        return g_value(moment_rows, slots, pos)
    if kind == "U":
        return u_value(moment_rows, slots, pos)
    if kind == "W":
        return w_value(moment_rows, slots, pos)
    raise AssertionError(kind)


def direct_block_value(kind: str, block: tuple[int, ...], pos: dict[int, int], moment_rows: tuple[IntSeries, ...]) -> tuple[int, ...]:
    value: tuple[int | Fraction, ...] = (1,) + (0,) * ORDER
    for slots in layer_slot_sets(kind, block):
        value = smul(value, mask_moment(moment_rows, slots, pos))
    return tuple(value)


@cache
def admissible_partitions(kind: str) -> tuple[tuple[tuple[tuple[int, ...], ...], int], ...]:
    answer = []
    for partition in P6:
        if all(all(len(slot_set) % 2 == 0 for slot_set in layer_slot_sets(kind, block)) for block in partition):
            answer.append((partition, mobius(partition)))
    return tuple(answer)


def direct_profile(kind: str, pos: dict[int, int], moment_rows: tuple[IntSeries, ...]) -> tuple[int, ...]:
    answer = izero()
    for partition, sign in admissible_partitions(kind):
        term: tuple[int | Fraction, ...] = (1,) + (0,) * ORDER
        for block in partition:
            term = smul(term, direct_block_value(kind, block, pos, moment_rows))
        sadd(answer, term, sign)
    return tuple(answer)


def formula_value(formula: dict[Mono, Fraction], pos: dict[int, int], moment_rows: tuple[IntSeries, ...]) -> tuple[Fraction, ...]:
    answer = rzero()
    for item, coefficient in formula.items():
        term: tuple[int | Fraction, ...] = (1,) + (0,) * ORDER
        cache_atom: dict[Atom, tuple[int, ...]] = {}
        for one_atom in item:
            if one_atom not in cache_atom:
                cache_atom[one_atom] = atom_value(moment_rows, one_atom, pos)
            term = smul(term, cache_atom[one_atom])
        sadd(answer, term, coefficient)
    return tuple(answer)


def finite_rectangle_identity_check(formula: dict[Mono, Fraction]) -> None:
    """Exhaustively verify the S5-collected identity on one 2x2 layer."""
    moment_rows = moments(2, 2)
    left = rzero()
    right = rzero()
    for free in product(range(4), repeat=5):
        pos = {1: 0, 2: free[0], 3: free[1], 4: free[2], 5: free[3], 6: free[4]}
        one = direct_profile("one", pos, moment_rows)
        two = direct_profile("two", pos, moment_rows)
        three = direct_profile("three", pos, moment_rows)
        sadd(left, one, Fraction(1, 720))
        sadd(left, two, Fraction(1, 24))
        sadd(left, three, Fraction(1, 8))
        sadd(right, formula_value(formula, pos, moment_rows))
    assert tuple(left) == tuple(right), "2x2 finite-volume G/U/W formula mismatch"


# ---------------------------------------------------------------------------
# Artifact normalization and test entry point
# ---------------------------------------------------------------------------
def artifact_formula(payload: dict[str, object]) -> dict[Mono, Fraction]:
    entries = payload["data"]["connected_correlation_identity"]["formula"]  # type: ignore[index]
    answer: dict[Mono, Fraction] = {}
    for entry in entries:  # type: ignore[union-attr]
        atoms = tuple(
            (atom_row["kind"], tuple(atom_row["slots"]))  # type: ignore[index]
            for atom_row in entry["atoms"]  # type: ignore[index]
        )
        answer[mono(*atoms)] = Fraction(entry["coefficient"])  # type: ignore[index]
    return answer


def as_fractions(rows: list[str]) -> RatSeries:
    return tuple(Fraction(value) for value in rows)


def main() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))["data"]

    independent, counts = independently_derived_formula()
    recorded = artifact_formula(payload)
    assert len(P6) == 203
    assert counts == {"one": 31, "two": 11, "three": 5}
    assert len(independent) == 13
    assert independent == recorded, "artifact formula is not the independent Bell(6) normal form"
    assert independent[mono(at("W", (1, 2, 3, 4, 5, 6)), at("W", (1, 2, 3, 4, 5, 6)))] == Fraction(1, 720)

    # This is an exact series evaluation of the formula before any infinite-volume FLM.
    finite_rectangle_identity_check(independent)

    bulk = independently_reexpanded_cumulants()
    c2 = bulk["c2_h2"]
    c4 = tuple(bulk["c4_h2"][d] + bulk["c4_h3"][d] for d in range(ORDER + 1))
    c6 = tuple(bulk["c6_h2"][d] + bulk["c6_h3"][d] + bulk["c6_h4"][d] for d in range(ORDER + 1))
    assert c6 == as_fractions(primary["c6_direct_coupling_q"]["full_infinite_stack"])
    assert bulk["c6_h2"] == as_fractions(primary["c6_direct_coupling_q"]["same_gap_height2"])
    assert bulk["c6_h3"] == as_fractions(primary["c6_direct_coupling_q"]["two_adjacent_gaps_height3"])
    assert bulk["c6_h4"] == as_fractions(primary["c6_direct_coupling_q"]["three_adjacent_gaps_height4"])

    assert c2 == as_fractions(primary["reproduced_lower_orders"]["c2_total_q_or_w"])
    assert c4 == as_fractions(primary["reproduced_lower_orders"]["c4_direct_q"])
    total_w = tuple(c6[d] + Fraction(4, 3) * c4[d] + Fraction(23, 45) * c2[d] for d in range(ORDER + 1))
    residual_w = list(total_w)
    residual_w[0] -= Fraction(1, 6)
    assert total_w == as_fractions(primary["c6_wave_variable_w_tanh_Kz"]["total"])
    assert tuple(residual_w) == as_fractions(primary["c6_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"])
    assert all(c6[d] == 0 for d in range(1, ORDER + 1, 2))
    print("PASS test_c6_closed")


if __name__ == "__main__":
    main()
