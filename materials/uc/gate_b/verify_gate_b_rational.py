#!/usr/bin/env python3
"""Exact rational Gate B certificate for both Cartesian-power bases.

Why this file exists
--------------------
The earlier Gate B certificates prove ``A_+(B) < 0`` and ``A_+(D) < 0`` in two
ways: 256-bit Arb balls with a three-way clamp classification, and 160-bit
dyadic intervals after relaxing every feasible action interval to ``[0, 1]``.
Both are sound, but each asks a referee to trust extra machinery -- a
third-party arithmetic library in the first case, an action-set relaxation plus
atanh/exp series with hand-derived tail bounds in the second.

This module removes all of it.  It evaluates the *true* one-sided causal
Bellman optimum over the *exact* feasible action interval using ordinary Python
fractions, and it needs only three certified scalar primitives, each a dozen
lines of integer arithmetic:

``log2_lower`` / ``log2_upper``
    bit-by-bit binary logarithm.  With ``f(y) = log2 y`` and any bit ``b``,
    ``f(y) = (b + f(y^2 / 2^b)) / 2`` exactly.  Replacing the state by a
    rational upper (lower) bound and using monotonicity of ``log2`` keeps the
    running value an upper (lower) bound; the truncation tail is in ``[0, 1)``
    because the state stays in ``[1, 2]``.

``exp2_upper``
    repeated ceiling square roots of two.  ``r_j >= 2^(2^-j)`` by induction, so
    the product over the set bits of a rounded-up dyadic fraction is an upper
    bound for ``2^x``.

The remaining mathematics is elementary and stated where it is used:

*   ``h(s) = -s log2 s - (1-s) log2 (1-s)`` is upper-bounded termwise because
    both coefficients ``-s`` and ``-(1-s)`` are negative, so a *lower* bound on
    each logarithm gives an upper bound on the entropy.
*   ``h(s) + c s`` is concave with derivative ``log2((1-s)/s) + c``, so on any
    ``[lo, hi]`` its maximum is at ``lo`` when the derivative is nonpositive at
    ``lo``, at ``hi`` when it is nonnegative at ``hi``, and otherwise bounded by
    the unconstrained Fenchel value ``log2(1 + 2^c)``.  Every branch is a valid
    upper bound; the first two are exact.
*   one Bellman step is ``max_s [h(s) + sum_ac P_ac(s) V_ac]`` with
    ``P_ac(s) >= 0`` on the feasible interval and ``sum_ac P_ac(s) = 1``, hence
    it is monotone in the child values.  Replacing children by upper bounds,
    and then ``const`` and the slope ``c`` by upper bounds (legitimate because
    ``s >= 0``), keeps the result an upper bound.

Consequently the printed rational number is a rigorous upper bound for
``A_+``, obtained with no interval library, no clamp classification, and no
relaxation of the action set.  It is the most elementary of the Gate B base
certificates and the only two-sided one; its upper endpoint agrees with the Arb
upper bound to roughly 26 digits, and the Arb bound remains marginally tighter
at that scale.  Both are far below the frozen rational targets, so either
settles the base lemma.

Resumability
------------
Each fixed-order Bellman optimum is appended to a JSONL checkpoint as an exact
rational.  A rerun replays the checkpoint, recomputes a deterministic sample of
stored records to detect corruption, and reports ``new_evaluations: 0`` when the
checkpoint is already complete.  The checkpoint is a cache, not an authority:
``--fresh`` ignores it entirely, and the canonical certificate is produced by a
complete run.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb, factorial, isqrt
import os
from pathlib import Path
import sys
import time
from typing import Any, Iterable

# Precision of the certified scalar primitives.  These only control tightness;
# soundness holds for any positive values.
LOG_STEPS = 96
LOG_STATE_BITS = 128
EXP_STEPS = 96
EXP_BITS = 128
ACCUMULATOR_BITS = 96

# Checkpoint records carry the schema of the evaluator that wrote them.  A
# record with any other schema is ignored rather than trusted.
RECORD_SCHEMA = 2

ALPHA = Fraction(356_069, 10_000_000)

HERE = Path(__file__).resolve().parent
DEFAULT_CHECKPOINT = HERE / "experiments" / "rational_certificate_checkpoint.jsonl"
DEFAULT_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_rational_v1.json"
)


# ---------------------------------------------------------------------------
# certified scalar primitives
# ---------------------------------------------------------------------------


def round_up(value: Fraction, bits: int) -> Fraction:
    """Least multiple of ``2^-bits`` that is at least ``value``."""
    scale = 1 << bits
    numerator = value.numerator * scale
    denominator = value.denominator
    return Fraction(-((-numerator) // denominator), scale)


def round_down(value: Fraction, bits: int) -> Fraction:
    """Greatest multiple of ``2^-bits`` that is at most ``value``."""
    scale = 1 << bits
    return Fraction((value.numerator * scale) // value.denominator, scale)


def _log2_bits(value: Fraction, upper: bool) -> Fraction:
    """Bit-by-bit binary logarithm; upper or lower bound by rounding direction."""
    assert value > 0
    exponent = 0
    state = value
    while state >= 2:
        state /= 2
        exponent += 1
    while state < 1:
        state *= 2
        exponent -= 1

    total = Fraction(exponent)
    weight = Fraction(1)
    for _ in range(LOG_STEPS):
        weight /= 2
        squared = state * state
        state = (
            round_up(squared, LOG_STATE_BITS)
            if upper
            else round_down(squared, LOG_STATE_BITS)
        )
        if state >= 2:
            total += weight
            state /= 2
    # Invariant: 1 <= state <= 2 in lower mode and 0 < state <= 2 in upper mode,
    # so the discarded tail weight * log2(state) lies in [0, weight].
    assert state <= 2
    if upper:
        return total + weight
    assert state >= 1
    return total


@lru_cache(maxsize=None)
def log2_upper(value: Fraction) -> Fraction:
    return _log2_bits(value, True)


@lru_cache(maxsize=None)
def log2_lower(value: Fraction) -> Fraction:
    return _log2_bits(value, False)


def _sqrt_up(value: Fraction) -> Fraction:
    """Least multiple of ``2^-EXP_BITS`` that is at least ``sqrt(value)``."""
    assert value >= 1
    scale = 1 << EXP_BITS
    target = -((-value.numerator * scale * scale) // value.denominator)
    root = isqrt(target)
    if root * root < target:
        root += 1
    return Fraction(root, scale)


@lru_cache(maxsize=None)
def _two_root(index: int) -> Fraction:
    """A rational at least ``2^(2^-index)``."""
    assert index >= 0
    if index == 0:
        return Fraction(2)
    return _sqrt_up(_two_root(index - 1))


@lru_cache(maxsize=None)
def exp2_upper(value: Fraction) -> Fraction:
    """A rational at least ``2^value``."""
    integral = value.numerator // value.denominator
    fractional = value - integral
    scale = 1 << EXP_STEPS
    scaled = -((-fractional.numerator * scale) // fractional.denominator)
    if scaled == scale:
        integral += 1
        scaled = 0

    result = Fraction(1)
    for index in range(1, EXP_STEPS + 1):
        if (scaled >> (EXP_STEPS - index)) & 1:
            result = round_up(result * _two_root(index), EXP_BITS)
    if integral >= 0:
        return result * (1 << integral)
    return round_up(result / (1 << -integral), EXP_BITS)


@lru_cache(maxsize=None)
def softplus2_upper(slope: Fraction) -> Fraction:
    """A rational at least ``log2(1 + 2^slope)``."""
    return log2_upper(Fraction(1) + exp2_upper(slope))


@lru_cache(maxsize=None)
def entropy_upper(probability: Fraction) -> Fraction:
    """A rational at least the binary entropy in bits."""
    assert 0 <= probability <= 1
    if probability == 0 or probability == 1:
        return Fraction(0)
    return -probability * log2_lower(probability) - (
        1 - probability
    ) * log2_lower(1 - probability)


@lru_cache(maxsize=None)
def max_affine_entropy_upper(
    slope: Fraction, lower: Fraction, upper: Fraction
) -> Fraction:
    """A rational at least ``max_{s in [lower, upper]} h(s) + slope * s``."""
    assert 0 <= lower <= upper <= 1
    if lower == upper:
        return entropy_upper(lower) + slope * lower
    if lower > 0 and slope <= log2_lower(lower / (1 - lower)):
        # Derivative log2((1-s)/s) + slope is nonpositive at ``lower``.
        return entropy_upper(lower) + slope * lower
    if upper < 1 and slope >= log2_upper(upper / (1 - upper)):
        # Derivative is nonnegative at ``upper``.
        return entropy_upper(upper) + slope * upper
    # Unconstrained Fenchel maximum over the superset [0, 1].
    return softplus2_upper(slope)


@lru_cache(maxsize=None)
def entropy_lower(probability: Fraction) -> Fraction:
    """A rational at most the binary entropy in bits."""
    assert 0 <= probability <= 1
    if probability == 0 or probability == 1:
        return Fraction(0)
    return -probability * log2_upper(probability) - (
        1 - probability
    ) * log2_upper(1 - probability)


@lru_cache(maxsize=None)
def feasible_near_argmax(
    slope: Fraction, lower: Fraction, upper: Fraction
) -> Fraction:
    """A feasible rational action close to the maximizer of ``h(s) + slope s``.

    Only feasibility matters for soundness of a lower bound; proximity to
    ``1 / (1 + 2^-slope)`` is what makes the resulting bound tight.
    """
    assert 0 <= lower <= upper <= 1
    if lower == upper:
        return lower
    candidate = round_down(
        Fraction(1) / (Fraction(1) + exp2_upper(-slope)), ACCUMULATOR_BITS
    )
    return min(max(candidate, lower), upper)


def decimal_upper(value: Fraction, digits: int = 40) -> str:
    """Decimal string that is at least ``value``, printed to ``digits`` places."""
    scale = 10**digits
    scaled = -((-value.numerator * scale) // value.denominator)
    sign = "-" if scaled < 0 else ""
    scaled = abs(scaled)
    integer, fraction = divmod(scaled, scale)
    return f"{sign}{integer}.{fraction:0{digits}d}"


def decimal_lower(value: Fraction, digits: int = 40) -> str:
    """Decimal string that is at most ``value``, printed to ``digits`` places."""
    scale = 10**digits
    scaled = (value.numerator * scale) // value.denominator
    sign = "-" if scaled < 0 else ""
    scaled = abs(scaled)
    integer, fraction = divmod(scaled, scale)
    return f"{sign}{integer}.{fraction:0{digits}d}"


# ---------------------------------------------------------------------------
# exact family combinatorics
# ---------------------------------------------------------------------------


def exact_reimer_threshold(size: int) -> int:
    """Least integer ``r`` with ``2^(2r) >= size^size``; no floating point."""
    assert size >= 1
    target = size**size
    low, high = 0, 1
    while (1 << (2 * high)) < target:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if (1 << (2 * middle)) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def _subsets_of_weight(
    coordinates: tuple[int, ...], weight: int
) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in chosen)
        for chosen in combinations(coordinates, weight)
    )


class Base:
    """A block-cell family together with its frozen exact facts."""

    def __init__(
        self,
        name: str,
        dimension: int,
        left: tuple[int, ...],
        right: tuple[int, ...],
        cells: frozenset[tuple[int, int]],
        expected_rows: tuple[int, ...],
        expected_count: int,
        expected_missing: int,
        reimer_witness: str,
    ) -> None:
        self.name = name
        self.dimension = dimension
        self.left = left
        self.right = right
        self.cells = cells
        self.expected_rows = expected_rows
        self.expected_count = expected_count
        self.expected_missing = expected_missing
        self.reimer_witness = reimer_witness

    def reconstruct(self) -> tuple[int, ...]:
        rows = tuple(
            sorted(
                {
                    left | right
                    for left_weight, right_weight in self.cells
                    for left in _subsets_of_weight(self.left, left_weight)
                    for right in _subsets_of_weight(self.right, right_weight)
                }
            )
        )
        assert rows == self.expected_rows
        return rows

    def exact_facts(self) -> dict[str, Any]:
        family = self.reconstruct()
        size = len(family)
        assert size == len(set(family))
        counts = tuple(
            sum((row >> coordinate) & 1 for row in family)
            for coordinate in range(self.dimension)
        )
        assert counts == (self.expected_count,) * self.dimension
        incidence = sum(counts)
        threshold = exact_reimer_threshold(size)
        assert incidence >= threshold
        assert max(counts) == 2 * size // 5

        support = set(family)
        missing = sum(
            (left | right) not in support for left in family for right in family
        )
        assert missing == self.expected_missing
        defect = Fraction(missing, size**2)

        columns = tuple(
            tuple((row >> coordinate) & 1 for row in family)
            for coordinate in range(self.dimension)
        )
        assert all(any(column) for column in columns)
        assert len(set(columns)) == self.dimension

        return {
            "dimension": self.dimension,
            "family_size": size,
            "block_partition": [list(self.left), list(self.right)],
            "block_cells": [list(cell) for cell in sorted(self.cells)],
            "family_rows": list(family),
            "coordinate_counts": list(counts),
            "total_incidence": incidence,
            "reimer_threshold": threshold,
            "cap_bound": 2 * size // 5,
            "missing_ordered_join_pairs": missing,
            "closure_defect": str(defect),
            "join_success_probability": str(1 - defect),
            "active": True,
            "separating": True,
            "ordered_pairs_with_replacement": True,
            "reimer_integer_witness": self.reimer_witness,
        }


BASES: dict[str, Base] = {
    "n6": Base(
        name="n6",
        dimension=6,
        left=(0, 1, 2),
        right=(3, 4, 5),
        cells=frozenset(((0, 0), (0, 1), (1, 0), (1, 2), (2, 1))),
        expected_rows=(
            0, 1, 2, 4, 8, 11, 13, 14, 16, 19, 21, 22, 25,
            26, 28, 32, 35, 37, 38, 41, 42, 44, 49, 50, 52,
        ),
        expected_count=10,
        expected_missing=444,
        reimer_witness="25**5 < 2**24",
    ),
    "n7": Base(
        name="n7",
        dimension=7,
        left=(5, 6),
        right=(0, 1, 2, 3, 4),
        cells=frozenset(
            ((0, 1), (0, 3), (1, 0), (1, 2), (1, 5), (2, 0), (2, 1))
        ),
        expected_rows=(
            1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 32,
            35, 37, 38, 41, 42, 44, 49, 50, 52, 56, 63, 64, 67, 69, 70,
            73, 74, 76, 81, 82, 84, 88, 95, 96, 97, 98, 100, 104, 112,
        ),
        expected_count=18,
        expected_missing=1_600,
        reimer_witness="45**5 < 2**28",
    ),
}


def permute_row(row: int, permutation: tuple[int, ...]) -> int:
    image = 0
    for source, target in enumerate(permutation):
        if (row >> source) & 1:
            image |= 1 << target
    return image


def automorphism_group(
    family: tuple[int, ...], dimension: int
) -> tuple[tuple[int, ...], ...]:
    support = frozenset(family)
    return tuple(
        permutation
        for permutation in permutations(range(dimension))
        if frozenset(permute_row(row, permutation) for row in family) == support
    )


def order_orbits(
    family: tuple[int, ...], dimension: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """Canonical order representatives and the automorphism group."""
    group = automorphism_group(family, dimension)
    remaining = set(permutations(range(dimension)))
    representatives: list[tuple[int, ...]] = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(permutation[index] for index in representative)
            for permutation in group
        }
        assert len(orbit) == len(group)
        assert orbit <= remaining
        remaining -= orbit
        representatives.append(representative)
    assert len(representatives) * len(group) == factorial(dimension)
    return tuple(representatives), group


# ---------------------------------------------------------------------------
# exact rational objective
# ---------------------------------------------------------------------------


def _predecessor_law(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> tuple[tuple[Fraction, Fraction], ...]:
    groups: dict[tuple[int, ...], tuple[int, int]] = {}
    for row in family:
        prefix = tuple((row >> predecessor) & 1 for predecessor in predecessors)
        total, ones = groups.get(prefix, (0, 0))
        groups[prefix] = (total + 1, ones + ((row >> coordinate) & 1))
    merged: Counter[Fraction] = Counter()
    size = len(family)
    for total, ones in groups.values():
        merged[Fraction(ones, total)] += Fraction(total, size)
    return tuple(sorted(merged.items()))


def shapley_iid_bounds(
    family: tuple[int, ...], dimension: int
) -> tuple[Fraction, Fraction]:
    """Rational enclosure of ``Q``, the order-averaged iid OR entropy."""
    assert family
    assert all(0 <= row < 1 << dimension for row in family)
    low = Fraction(0)
    high = Fraction(0)
    for coordinate in range(dimension):
        others = tuple(index for index in range(dimension) if index != coordinate)
        for predecessor_count in range(dimension):
            weight = Fraction(
                1, dimension * comb(dimension - 1, predecessor_count)
            )
            for predecessors in combinations(others, predecessor_count):
                law = _predecessor_law(family, coordinate, predecessors)
                local_low = Fraction(0)
                local_high = Fraction(0)
                for left, left_mass in law:
                    for right, right_mass in law:
                        probability = left + right - left * right
                        mass = left_mass * right_mass
                        local_low += entropy_lower(probability) * mass
                        local_high += entropy_upper(probability) * mass
                low = round_down(low + weight * local_low, ACCUMULATOR_BITS)
                high = round_up(high + weight * local_high, ACCUMULATOR_BITS)
    assert low <= high
    return low, high


def _feasible_interval(
    left: Fraction, right: Fraction
) -> tuple[Fraction, Fraction]:
    lower = max(left, right, min(left + right, Fraction(1, 2)))
    upper = min(Fraction(1), left + right)
    return lower, upper


def bellman_bounds(
    family: tuple[int, ...],
    dimension: int,
    order: tuple[int, ...],
    canonicalize: bool = True,
) -> tuple[Fraction, Fraction, int]:
    """Rational enclosure of the true fixed-order one-sided Bellman optimum.

    The upper branch maximizes a certified entropy majorant over the exact
    feasible action interval.  The lower branch evaluates one explicitly
    feasible rational policy with pessimistic children, so it is the value of an
    admissible one-sided coupling and therefore at most the optimum.
    """
    assert tuple(sorted(order)) == tuple(range(dimension))
    assert family
    assert all(0 <= row < 1 << dimension for row in family)

    @lru_cache(maxsize=None)
    def value(
        position: int,
        left_rows: tuple[int, ...],
        right_rows: tuple[int, ...],
    ) -> tuple[Fraction, Fraction]:
        if position == dimension:
            return Fraction(0), Fraction(0)
        if canonicalize and left_rows > right_rows:
            left_rows, right_rows = right_rows, left_rows

        coordinate = order[position]
        left_zero = tuple(row for row in left_rows if not (row >> coordinate) & 1)
        left_one = tuple(row for row in left_rows if (row >> coordinate) & 1)
        right_zero = tuple(
            row for row in right_rows if not (row >> coordinate) & 1
        )
        right_one = tuple(row for row in right_rows if (row >> coordinate) & 1)
        left = Fraction(len(left_one), len(left_rows))
        right = Fraction(len(right_one), len(right_rows))

        children = {
            (left_bit, right_bit): value(position + 1, left_child, right_child)
            for left_bit, left_child in ((0, left_zero), (1, left_one))
            if left_child
            for right_bit, right_child in ((0, right_zero), (1, right_one))
            if right_child
        }
        lower, upper = _feasible_interval(left, right)

        def transitions(action: Fraction) -> tuple[tuple[int, int, Fraction], ...]:
            return (
                (0, 0, 1 - action),
                (1, 0, action - right),
                (0, 1, action - left),
                (1, 1, left + right - action),
            )

        if len(children) < 4:
            # A missing child forces left or right to be deterministic, and then
            # the feasible interval degenerates to a single action.
            assert lower == upper
            low = entropy_lower(lower)
            high = entropy_upper(lower)
            for left_bit, right_bit, probability in transitions(lower):
                assert probability >= 0
                if probability:
                    child_low, child_high = children[left_bit, right_bit]
                    low += child_low * probability
                    high += child_high * probability
            low = round_down(low, ACCUMULATOR_BITS)
            high = round_up(high, ACCUMULATOR_BITS)
            assert low <= high
            return low, high

        slope = round_up(
            -children[0, 0][1]
            + children[1, 0][1]
            + children[0, 1][1]
            - children[1, 1][1],
            ACCUMULATOR_BITS,
        )
        constant = round_up(
            children[0, 0][1]
            - children[1, 0][1] * right
            - children[0, 1][1] * left
            + children[1, 1][1] * (left + right),
            ACCUMULATOR_BITS,
        )
        high = round_up(
            constant + max_affine_entropy_upper(slope, lower, upper),
            ACCUMULATOR_BITS,
        )

        policy_slope = round_down(
            -children[0, 0][0]
            + children[1, 0][0]
            + children[0, 1][0]
            - children[1, 1][0],
            ACCUMULATOR_BITS,
        )
        action = feasible_near_argmax(policy_slope, lower, upper)
        low = entropy_lower(action)
        for left_bit, right_bit, probability in transitions(action):
            assert probability >= 0
            if probability:
                low += children[left_bit, right_bit][0] * probability
        low = round_down(low, ACCUMULATOR_BITS)
        assert low <= high
        return low, high

    root = tuple(family)
    low, high = value(0, root, root)
    return low, high, value.cache_info().currsize


# ---------------------------------------------------------------------------
# resumable driver
# ---------------------------------------------------------------------------


def _read_checkpoint(
    path: Path,
) -> tuple[
    dict[tuple[str, tuple[int, ...]], tuple[Fraction, Fraction]], int, int
]:
    records: dict[tuple[str, tuple[int, ...]], tuple[Fraction, Fraction]] = {}
    truncated = 0
    foreign = 0
    if not path.exists():
        return records, truncated, foreign
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                truncated += 1
                continue
            if payload.get("kind") != "order":
                continue
            if payload.get("schema") != RECORD_SCHEMA:
                # A record written by a different evaluator revision is a stale
                # cache entry, never a substitute for a current evaluation.
                foreign += 1
                continue
            base = payload["base"]
            order = tuple(payload["order"])
            dimension = BASES[base].dimension
            if tuple(sorted(order)) != tuple(range(dimension)):
                raise ValueError(f"corrupt order record for {base}: {order}")
            low = Fraction(payload["bellman_lower"])
            high = Fraction(payload["bellman_upper"])
            if not 0 < low <= high < dimension:
                raise ValueError(
                    f"implausible stored enclosure for {base}: {low}, {high}"
                )
            records[base, order] = (low, high)
    return records, truncated, foreign


def _append_record(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _orders(base: Base, mode: str) -> tuple[tuple[int, ...], ...]:
    if mode == "all":
        return tuple(permutations(range(base.dimension)))
    family = base.reconstruct()
    representatives, _group = order_orbits(family, base.dimension)
    return representatives


def evaluate_base(
    base: Base,
    order_mode: str,
    checkpoint: Path | None,
    sample_recheck: int,
    progress: int,
) -> dict[str, Any]:
    family = base.reconstruct()
    facts = base.exact_facts()
    orders = _orders(base, order_mode)
    stored, truncated, foreign = (
        ({}, 0, 0) if checkpoint is None else _read_checkpoint(checkpoint)
    )

    resumed = 0
    evaluated = 0
    rechecked = 0
    states = 0
    values: list[tuple[Fraction, Fraction]] = []
    started = time.monotonic()

    for index, order in enumerate(orders):
        cached = stored.get((base.name, order))
        if cached is not None:
            recheck = (
                sample_recheck > 0
                and rechecked < sample_recheck
                and index % max(1, len(orders) // sample_recheck) == 0
            )
            if recheck:
                low, high, order_states = bellman_bounds(
                    family, base.dimension, order
                )
                if (low, high) != cached:
                    raise ValueError(
                        f"checkpoint mismatch for {base.name} {order}: "
                        f"stored {cached} recomputed {(low, high)}"
                    )
                rechecked += 1
                states += order_states
            values.append(cached)
            resumed += 1
            continue

        low, high, order_states = bellman_bounds(family, base.dimension, order)
        values.append((low, high))
        evaluated += 1
        states += order_states
        if checkpoint is not None:
            _append_record(
                checkpoint,
                {
                    "kind": "order",
                    "schema": RECORD_SCHEMA,
                    "base": base.name,
                    "order": list(order),
                    "bellman_lower": str(low),
                    "bellman_upper": str(high),
                    "states": order_states,
                },
            )
        if progress and evaluated % progress == 0:
            elapsed = time.monotonic() - started
            print(
                f"[{base.name}] {evaluated} new / {len(orders)} orders "
                f"({elapsed:.1f} s)",
                file=sys.stderr,
                flush=True,
            )

    assert len(values) == len(orders)
    count = len(values)
    average_low = round_down(
        sum((low for low, _high in values), Fraction(0)) / count,
        ACCUMULATOR_BITS,
    )
    average_high = round_up(
        sum((high for _low, high in values), Fraction(0)) / count,
        ACCUMULATOR_BITS,
    )
    iid_low, iid_high = shapley_iid_bounds(family, base.dimension)
    size = Fraction(len(family))
    a_plus_low = round_down(
        (1 - ALPHA) * iid_low + ALPHA * average_low - log2_upper(size),
        ACCUMULATOR_BITS,
    )
    a_plus_high = round_up(
        (1 - ALPHA) * iid_high + ALPHA * average_high - log2_lower(size),
        ACCUMULATOR_BITS,
    )
    assert a_plus_low <= a_plus_high

    representatives, group = order_orbits(family, base.dimension)
    distinct = Counter(values)
    return {
        "base": base.name,
        "exact_base": facts,
        "order_evaluation": {
            "mode": order_mode,
            "evaluated_order_count": len(orders),
            "all_order_count": factorial(base.dimension),
            "symmetry_quotient_used": order_mode != "all",
            "automorphism_count": len(group),
            "order_orbit_count": len(representatives),
            "distinct_enclosure_count": len(distinct),
            "distinct_enclosure_multiplicities": sorted(distinct.values()),
            "bellman_state_count": states,
            "resumed_records": resumed,
            "new_evaluations": evaluated,
            "rechecked_records": rechecked,
            "truncated_checkpoint_lines": truncated,
            "foreign_schema_checkpoint_lines": foreign,
        },
        "rational_values": {
            "shapley_iid_lower": str(iid_low),
            "shapley_iid_upper": str(iid_high),
            "shapley_iid_lower_decimal": decimal_lower(iid_low),
            "shapley_iid_upper_decimal": decimal_upper(iid_high),
            "bellman_average_lower": str(average_low),
            "bellman_average_upper": str(average_high),
            "bellman_average_lower_decimal": decimal_lower(average_low),
            "bellman_average_upper_decimal": decimal_upper(average_high),
            "family_entropy_lower": str(log2_lower(size)),
            "family_entropy_upper": str(log2_upper(size)),
            "family_entropy_lower_decimal": decimal_lower(log2_lower(size)),
            "a_plus_lower": str(a_plus_low),
            "a_plus_upper": str(a_plus_high),
            "a_plus_lower_decimal": decimal_lower(a_plus_low),
            "a_plus_upper_decimal": decimal_upper(a_plus_high),
            "a_plus_enclosure_width_decimal": decimal_upper(
                a_plus_high - a_plus_low
            ),
            "minimum_order_upper": str(min(high for _low, high in values)),
            "maximum_order_upper": str(max(high for _low, high in values)),
        },
        "_a_plus": a_plus_high,
        "_a_plus_lower": a_plus_low,
        "_defect": Fraction(facts["closure_defect"]),
    }


# ---------------------------------------------------------------------------
# certificate assembly
# ---------------------------------------------------------------------------

RATIONAL_TARGETS = {
    "n6": (Fraction(-17, 1250), Fraction(-1, 80)),
    "n7": (Fraction(-7, 250), Fraction(-1, 40)),
}


def _power_consequence(
    base: Base, a_plus: Fraction, defect: Fraction, threshold: Fraction
) -> dict[str, Any]:
    success = 1 - defect
    size = len(base.expected_rows)
    ratios = {}
    for power in (1, 2, 3, 10, 100, 1_000):
        power_defect = 1 - success**power
        assert 0 < power_defect < 1
        ratio = -threshold * power / power_defect
        assert ratio > -threshold * power
        ratios[str(power)] = decimal_upper(ratio, 12)
    return {
        "family": f"{base.name}^(boxtimes k)",
        "dimension": f"{base.dimension}*k",
        "family_size": f"{size}^k",
        "coordinate_count": (
            f"{base.expected_count}*{size}^(k-1)"
        ),
        "closure_defect": f"1-({success})^k",
        "a_plus_upper": f"{threshold}*k",
        "repair_ratio_lower": (
            f"({-threshold}*k)/(1-({success})^k) > {-threshold}*k"
        ),
        "reimer_integer_witness": base.reimer_witness,
        "sampled_ratio_lower_bounds": ratios,
        "zero_defect_convention_used": False,
    }


def _asymptotic_block(results: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Growth rate of the Gate B supremum restricted to dimension at most n."""
    block: dict[str, Any] = {
        "statement": (
            "Write c(n) for the supremum of -A_+/eps_join over admissible "
            "families on at most n coordinates with A_+ < 0.  Cartesian powers "
            "of a certified base give c(n) >= |delta| * floor(n/d) with the "
            "base dimension d, hence c(n) = Omega(n).  Conversely Q >= 0 and "
            "C_+ >= 0 force -A_+ <= log2 m <= n, so for every family with "
            "eps_join >= eps0 the ratio is at most n/eps0.  Under a positive "
            "defect floor the growth is therefore exactly of order n."
        ),
        "upper_bound_reason": (
            "h >= 0 termwise gives Q >= 0; C_+ is a maximum of sums of h, so "
            "C_+ >= 0; hence A_+ >= -log2 m and -A_+ <= log2 m <= n."
        ),
        "bases": {},
    }
    for name, result in results.items():
        base = BASES[name]
        threshold = min(
            candidate
            for candidate in RATIONAL_TARGETS[name]
            if result["_a_plus"] < candidate
        )
        block["bases"][name] = {
            "base_dimension": base.dimension,
            "certified_threshold": str(threshold),
            "linear_lower_bound": (
                f"c(n) >= {-threshold} * floor(n/{base.dimension})"
            ),
            "asymptotic_slope": str(-threshold / base.dimension),
        }
    return block


def verify(
    base_names: Iterable[str] = ("n6", "n7"),
    order_mode: str = "all",
    checkpoint: Path | None = DEFAULT_CHECKPOINT,
    sample_recheck: int = 2,
    progress: int = 0,
) -> dict[str, Any]:
    results: dict[str, dict[str, Any]] = {}
    for name in base_names:
        results[name] = evaluate_base(
            BASES[name], order_mode, checkpoint, sample_recheck, progress
        )

    bases_payload: dict[str, Any] = {}
    for name, result in results.items():
        base = BASES[name]
        a_plus = result["_a_plus"]
        defect = result["_defect"]
        passed = [
            str(target)
            for target in RATIONAL_TARGETS[name]
            if a_plus < target
        ]
        assert passed, f"{name}: rational bound {a_plus} beat no frozen target"
        sharpest = min(
            (target for target in RATIONAL_TARGETS[name] if a_plus < target),
        )
        payload = {
            key: value
            for key, value in result.items()
            if not key.startswith("_")
        }
        payload["rational_consequence"] = {
            "a_plus_upper_lt": passed,
            "sharpest_frozen_target_passed": str(sharpest),
            "base_ratio_lower_gt": str(-sharpest / defect),
            "a_plus_two_sided_enclosure": [
                result["rational_values"]["a_plus_lower_decimal"],
                result["rational_values"]["a_plus_upper_decimal"],
            ],
            "enclosure_is_negative": result["_a_plus"] < 0
            and result["_a_plus_lower"] < 0,
        }
        payload["infinite_construction"] = _power_consequence(
            base, a_plus, defect, sharpest
        )
        bases_payload[name] = payload

    return {
        "schema_version": 1,
        "verdict": "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL",
        "arithmetic": "exact_python_fractions",
        "third_party_dependencies": [],
        "action_set": "exact feasible interval; no clamp-free relaxation",
        "primitives": {
            "log2_steps": LOG_STEPS,
            "log2_state_bits": LOG_STATE_BITS,
            "exp2_steps": EXP_STEPS,
            "exp2_bits": EXP_BITS,
            "accumulator_bits": ACCUMULATOR_BITS,
            "log2": (
                "bit-by-bit squaring with monotone rounding; tail in [0, 1)"
            ),
            "exp2": "repeated ceiling square roots of two",
            "entropy": (
                "termwise from log2_lower (upper bound) or log2_upper (lower "
                "bound) because both weights are negative"
            ),
            "interval_maximum": (
                "concavity case split: exact at an endpoint, otherwise the "
                "Fenchel value log2(1 + 2^slope)"
            ),
            "interval_minimum": (
                "one explicitly feasible rational action near the maximizer, "
                "evaluated with pessimistic children"
            ),
        },
        "alpha": {
            "numerator": ALPHA.numerator,
            "denominator": ALPHA.denominator,
            "dimension_independent": True,
        },
        "bases": bases_payload,
        "asymptotic": _asymptotic_block(results),
        "proof_boundary": (
            "Both finite negative bases are certified here in exact rational "
            "arithmetic against the true feasible action set.  Product "
            "additivity and all-k admissibility are the exact lemmas of "
            "PROOF.md, not finite-n extrapolation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bases", default="n6,n7", help="comma separated subset of n6,n7"
    )
    parser.add_argument(
        "--orders",
        default="all",
        choices=("all", "orbits"),
        help="evaluate every coordinate order or one per automorphism orbit",
    )
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="ignore and do not write the checkpoint",
    )
    parser.add_argument(
        "--sample-recheck",
        type=int,
        default=2,
        help="recompute this many stored records per base",
    )
    parser.add_argument("--progress", type=int, default=0)
    parser.add_argument("--write-certificate", type=Path)
    args = parser.parse_args()

    certificate = verify(
        base_names=tuple(
            name for name in args.bases.split(",") if name.strip()
        ),
        order_mode=args.orders,
        checkpoint=None if args.fresh else args.checkpoint,
        sample_recheck=args.sample_recheck,
        progress=args.progress,
    )
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write_certificate is not None:
        args.write_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.write_certificate.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
