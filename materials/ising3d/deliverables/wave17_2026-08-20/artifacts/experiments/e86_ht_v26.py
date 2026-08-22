"""Exact simple-cubic high-temperature series through v^26 by cut-capped frontiers.

Run from the repository root with::

    .venv/bin/python experiments/e86_ht_v26.py

The v^24 engine wall (proofs/ht_v24.md) was the 5x5x5 minimal-span weight.  At
order 26 three box classes exceed every spin-transfer engine, because both
engines cap the open cross-section at 22 sites:

* (5,5,5) cross-section 25 needs its degree-26 (non-minimal) weight;
* (5,5,6) cross-section 25 needs its minimal-span weight at degree 26;
* (4,6,6) cross-section 24 needs its minimal-span weight at degree 26.

This module generalizes the v^24 vertex-frontier to arbitrary boxes and to one
cut capped at four instead of two.  At degree 26 a wrapping contribution with
12 cuts crosses one cut four times and the others twice, so a single
"4-resolving" run per cut orbit extracts the degree-26 weight of (5,5,5).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from math import comb, factorial, gcd, prod
from pathlib import Path

import ising.series as flm
from ising.series import FLMSeries
from ising.transfer_matrix.crt import _descending_primes, _is_prime_32

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e86_ht_v26.py"
INTERPRETER = ".venv/bin/python"
V24_PATH = ROOT / "results" / "series" / "ht_v24.json"
PREFIX_PATH = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
RESULT_PATH = ROOT / "results" / "series" / "ht_v26.json"
PROOF_PATH = ROOT / "proofs" / "ht_v26.md"

ORDER = 26
MEMORY_BUDGET_BYTES = 40_000_000_000
WALLED_SHAPES = ((5, 5, 5), (5, 5, 6), (4, 6, 6))
# (5,5,5) cut orbits at degree 26: reflection pairs (0,3),(1,2) per axis, so
# edge planes {0,3} and middle planes {1,2} are the two orbits, six cuts each.
SPECIAL_CUTS_555 = (0, 1, 5, 9)  # x-edge, x-middle, y-middle, z-middle

# [EXTERNAL] Table I of Arisue and Fujiwara, arXiv:hep-lat/0209002 (2002),
# lists the simple-cubic HT free-energy interaction series f = 3 log cosh(beta)
# + sum_n a_n t^n, t = tanh(beta), through n = 46.  Their order-26 entry is the
# falsification target for this experiment's new coefficient.  The sha256 is of
# the exact arXiv PDF bytes retrieved and read for this comparison.
EXTERNAL_ARISUE_FUJIWARA = {
    "arxiv": "hep-lat/0209002",
    "title": "New Algorithm of the Finite Lattice Method for the High-temperature Expansion of the Ising Model in Three Dimensions",
    "table_i_sha256": "bf79ec7918694768b57e6793ed4b0c8f581aaea6a705b4273ecab42641723afc",
    "table_i_interaction": {
        2: "0",
        4: "3",
        6: "22",
        8: "375/2",
        10: "1980",
        12: "24044",
        14: "319170",
        16: "18059031/4",
        18: "201010408/3",
        20: "5162283633/5",
        22: "16397040750",
        24: "266958797382",
        26: "4437596650548",
    },
}

CRT_PRIMES = tuple(int(value) for value in _descending_primes(6))
CRT_PRODUCT = prod(CRT_PRIMES)


def _shape_sides(shape) -> tuple[int, int, int]:
    sides = tuple(int(side) for side in shape)
    if len(sides) != 3 or any(side < 1 for side in sides):
        raise ValueError("shape must be three positive side lengths")
    return sides


def _cut_inventory(shape) -> tuple[int, list[int]]:
    """Cuts and per-cut crossing-edge counts for transfer along ``shape[0]``.

    Every box edge crosses exactly one coordinate cut: the +axis edge from a
    vertex crosses the cut between that plane pair.  Digit ``k`` of the frontier
    code tracks cut ``k`` with x-cuts first, then y-cuts, then z-cuts.
    """

    a, b, c = _shape_sides(shape)
    counts = [b * c] * (a - 1) + [a * c] * (b - 1) + [a * b] * (c - 1)
    return len(counts), counts


def _digits_of(code: int, radices: list[int]) -> list[int]:
    digits = []
    rest = code
    for radix in radices:
        digits.append(rest % radix)
        rest //= radix
    if rest:
        raise AssertionError("cut code exceeds its mixed-radix range")
    return digits


def frontier_cut_polynomial(
    shape,
    modulus: int = CRT_PRIMES[0],
    max_states: int | None = None,
    special_cut: int | None = None,
):
    """Count even subgraphs with cut counters capped at two (four if special).

    The vertex frontier walks the box in transfer-direction-major order.  The
    pending register holds, in bit ``b``, the accumulated parity of all
    not-yet-consumed chosen edges ending ``b + 1`` vertices ahead, so every
    vertex degree stays even.  Cut counters are ternary (one quinary digit for
    ``special_cut``); a further crossing past the cap is discarded.
    """

    a, b, c = _shape_sides(shape)
    n_cuts, _ = _cut_inventory(shape)
    if special_cut is not None and not 0 <= int(special_cut) < n_cuts:
        raise ValueError("special cut index out of range")
    radices = [3] * n_cuts
    caps = [2] * n_cuts
    if special_cut is not None:
        radices[special_cut] = 5
        caps[special_cut] = 4
    powers = []
    accumulator = 1
    for radix in radices:
        powers.append(accumulator)
        accumulator *= radix
    width = b * c
    states: dict[int, int] = {0: 1}
    peak = 1
    timeline = []
    started = time.perf_counter()

    for x in range(a):
        for y in range(b):
            for z in range(c):
                directions = []
                if z + 1 < c:
                    directions.append((0, (a - 1) + (b - 1) + z))
                if y + 1 < b:
                    directions.append((c - 1, (a - 1) + y))
                if x + 1 < a:
                    directions.append((width - 1, x))

                options_by_parity = ([], [])
                for choice in range(1 << len(directions)):
                    parity = choice.bit_count() & 1
                    pending_delta = 0
                    increments = []
                    for index, (offset, digit) in enumerate(directions):
                        if choice & (1 << index):
                            pending_delta ^= 1 << offset
                            increments.append(
                                (powers[digit], radices[digit], caps[digit])
                            )
                    options_by_parity[parity].append(
                        (pending_delta, tuple(increments))
                    )

                updated: dict[int, int] = {}
                for packed, value in states.items():
                    pending = packed & ((1 << width) - 1)
                    code = packed >> width
                    incoming = pending & 1
                    base_pending = pending >> 1
                    for pending_delta, increments in options_by_parity[incoming]:
                        new_code = code
                        allowed = True
                        for power, radix, cap in increments:
                            if (new_code // power) % radix == cap:
                                allowed = False
                                break
                            new_code += power
                        if not allowed:
                            continue
                        key = (new_code << width) | (base_pending ^ pending_delta)
                        updated[key] = (updated.get(key, 0) + value) % modulus
                states = updated
                peak = max(peak, len(states))
                timeline.append(
                    {
                        "vertex": [x, y, z],
                        "states": len(states),
                        "elapsed_seconds": round(time.perf_counter() - started, 6),
                    }
                )
                if max_states is not None and len(states) > max_states:
                    return {
                        "completed": False,
                        "shape": [a, b, c],
                        "special_cut": special_cut,
                        "modulus": modulus,
                        "peak_states": peak,
                        "terminal_states": len(states),
                        "timeline": timeline,
                        "coefficients": None,
                    }

    coefficients: dict[int, int] = defaultdict(int)
    for packed, value in states.items():
        pending = packed & ((1 << width) - 1)
        if pending:
            raise AssertionError("nonempty parity frontier after final vertex")
        code = packed >> width
        digits = _digits_of(code, radices)
        if any(digit % 2 for digit in digits):
            raise AssertionError("an even subgraph crosses every cut evenly")
        coefficients[code] = (coefficients[code] + value) % modulus
    return {
        "completed": True,
        "shape": [a, b, c],
        "special_cut": special_cut,
        "modulus": modulus,
        "peak_states": peak,
        "terminal_states": len(states),
        "timeline": timeline,
        "coefficients": dict(coefficients),
    }


def _radix_powers(radices: list[int]) -> list[int]:
    powers = []
    value = 1
    for radix in radices:
        if radix < 2:
            raise ValueError("mixed-radix bases must be at least two")
        powers.append(value)
        value *= radix
    return powers


def _mixed_code(digits: list[int], powers: list[int]) -> int:
    return sum(digit * power for digit, power in zip(digits, powers))


def cut_polynomial_as_half_crossings(
    coefficients: dict[int, int], shape, special_cut: int | None, modulus: int
) -> tuple[tuple[int, ...], list[int]]:
    """Map terminal cut counts to ordinary half-crossing monomials.

    The exponent of a cut variable is its crossing count divided by two.  All
    ordinary cuts therefore have target exponent at most one; a distinguished
    four-crossing cut has target exponent at most two.  Unlike a square-free
    projection, this retains the product of two separate two-crossing factors.
    """

    n_cuts, _ = _cut_inventory(shape)
    terminal_radices = [3] * n_cuts
    exponent_radices = [2] * n_cuts
    if special_cut is not None:
        terminal_radices[special_cut] = 5
        exponent_radices[special_cut] = 3
    powers = _radix_powers(exponent_radices)
    result = [0] * prod(exponent_radices)
    for code, value in coefficients.items():
        terminal_digits = _digits_of(code, terminal_radices)
        if any(digit % 2 for digit in terminal_digits):
            raise AssertionError("terminal cut count is odd")
        exponents = [digit // 2 for digit in terminal_digits]
        exponent_code = _mixed_code(exponents, powers)
        result[exponent_code] = (result[exponent_code] + value) % modulus
    if result[0] != 1:
        raise AssertionError("the empty even subgraph must be unique")
    return tuple(result), exponent_radices


def truncated_mixed_log(
    coefficient: tuple[int, ...], radices: list[int], modulus: int
) -> list[int]:
    """Formal log in the ordinary mixed-radix polynomial quotient.

    If ``P = exp(L)``, coefficient comparison in
    ``partial_j P = (partial_j L) P`` gives

    ``L_a = P_a - a_j^-1 sum_{0 < g < a, g_j > 0}
       g_j L_g P_{a-g}``.

    This handles the degree-26 2+2 split of one cut across two log factors,
    which a square-free algebra necessarily discards.
    """

    expected_size = prod(radices)
    if len(coefficient) != expected_size or coefficient[0] % modulus != 1:
        raise ValueError("mixed log requires the complete unit-constant quotient")
    powers = _radix_powers(radices)
    digits_by_code = [_digits_of(code, radices) for code in range(expected_size)]
    logarithm = [0] * expected_size
    ordered_codes = sorted(
        range(1, expected_size),
        key=lambda code: (sum(digits_by_code[code]), code),
    )
    for alpha_code in ordered_codes:
        alpha = digits_by_code[alpha_code]
        axis = next(index for index, digit in enumerate(alpha) if digit)
        if gcd(alpha[axis], modulus) != 1:
            raise AssertionError("CRT modulus must invert every log denominator")
        subtotal = 0
        gamma = [0] * len(radices)

        def visit(index: int, gamma_code: int) -> None:
            nonlocal subtotal
            if index == len(radices):
                if gamma_code != alpha_code:
                    delta_code = alpha_code - gamma_code
                    subtotal += (
                        gamma[axis]
                        * logarithm[gamma_code]
                        * coefficient[delta_code]
                    )
                return
            lower = 1 if index == axis else 0
            for digit in range(lower, alpha[index] + 1):
                gamma[index] = digit
                visit(index + 1, gamma_code + digit * powers[index])

        visit(0, 0)
        logarithm[alpha_code] = (
            coefficient[alpha_code]
            - subtotal * pow(alpha[axis], -1, modulus)
        ) % modulus
    return logarithm


def ordered_partition_factor(n: int) -> int:
    """Return ``sum_k S(n,k)*(k-1)!`` by an integer Stirling recurrence."""

    stirling = [[0] * (n + 1) for _ in range(n + 1)]
    stirling[0][0] = 1
    for row in range(1, n + 1):
        for blocks in range(1, row + 1):
            stirling[row][blocks] = (
                stirling[row - 1][blocks - 1] + blocks * stirling[row - 1][blocks]
            )
    return sum(
        stirling[n][blocks] * factorial(blocks - 1)
        for blocks in range(1, n + 1)
    )


def frontier_log_bound(shape, special_cut: int | None = None) -> int:
    """Strict absolute bound for the selected ordinary-log coefficient.

    With no special cut there is one two-crossing choice per cut.  With one
    four-crossing cut, regard its two half-crossing units as labelled slots.
    Every ordinary log composition maps to at least one ordered partition of
    those slots.  Bounding each special slot by ``C(E,4)`` (also when the two
    units land in distinct factors) gives a deliberately conservative integer
    bound, sufficient for CRT uniqueness.
    """

    n_cuts, cut_edges = _cut_inventory(shape)
    if special_cut is None:
        bound = prod(comb(edges, 2) for edges in cut_edges)
        return bound * ordered_partition_factor(n_cuts)
    bound = comb(cut_edges[special_cut], 4) ** 2
    for index, edges in enumerate(cut_edges):
        if index != special_cut:
            bound *= comb(edges, 2)
    return bound * ordered_partition_factor(n_cuts + 1)


def _centered(residue: int, modulus: int) -> int:
    return residue if residue <= modulus // 2 else residue - modulus


def extract_frontier_weights(
    shape,
    modulus: int = CRT_PRODUCT,
    max_states: int | None = None,
    special_cut: int | None = None,
) -> dict:
    """Run the frontier and extract log coefficients with telemetry."""

    started = time.perf_counter()
    frontier = frontier_cut_polynomial(
        shape, modulus=modulus, max_states=max_states, special_cut=special_cut
    )
    wall = time.perf_counter() - started
    if not frontier["completed"]:
        return {
            "completed": False,
            "shape": frontier["shape"],
            "special_cut": special_cut,
            "peak_states": frontier["peak_states"],
            "terminal_states": frontier["terminal_states"],
            "timeline": frontier["timeline"],
            "wall_seconds": f"{wall:.6f}",
        }
    n_cuts, _ = _cut_inventory(shape)
    half_polynomial, exponent_radices = cut_polynomial_as_half_crossings(
        frontier["coefficients"], shape, special_cut, modulus
    )
    logarithm = truncated_mixed_log(half_polynomial, exponent_radices, modulus)
    exponent_powers = _radix_powers(exponent_radices)
    minimal_code = _mixed_code([1] * n_cuts, exponent_powers)
    payload = {
        "completed": True,
        "shape": frontier["shape"],
        "special_cut": special_cut,
        "peak_states": frontier["peak_states"],
        "terminal_states": frontier["terminal_states"],
        "timeline": frontier["timeline"],
        "wall_seconds": f"{wall:.6f}",
        "modulus": str(modulus),
        "minimal_residue": str(logarithm[minimal_code]),
        "minimal_centered": _centered(logarithm[minimal_code], modulus),
        "bound": str(frontier_log_bound(shape, special_cut)),
    }
    if special_cut is not None:
        four_cut_digits = [1] * n_cuts
        four_cut_digits[special_cut] = 2
        four_cut_code = _mixed_code(four_cut_digits, exponent_powers)
        payload["four_cut_residue"] = str(logarithm[four_cut_code])
        payload["four_cut_centered"] = _centered(
            logarithm[four_cut_code], modulus
        )
    return payload


def minimal_span_classes(order: int = ORDER) -> list[tuple[int, int, int]]:
    """Canonical classes whose minimal degree is exactly ``order``."""

    classes = []
    for total in range(3, order + 4):
        if total - 3 != order // 2:
            continue
        for a in range(1, total):
            for b in range(1, a + 1):
                c = total - a - b
                if c < 1 or c > b:
                    continue
                if 2 * (a - 1 + b - 1 + c - 1) == order:
                    classes.append((a, b, c))
    return sorted(classes)


def high_temperature_series_with_injections(
    order: int, injections: dict[tuple[int, int, int], dict[int, int]]
) -> FLMSeries:
    """Usual FLM at ``order`` with proved weights injected for walled boxes."""

    if int(order) != ORDER:
        raise ValueError("the specialized cut extraction is proved only at order 26")
    shapes = flm._ht_shapes(3, order, 0)
    weights = {}
    interaction = [Fraction(0) for _ in range(order + 1)]
    injected_seen = set()
    for shape in shapes:
        canonical = flm._canonical_shape(shape)
        minimal_degree = 2 * sum(side - 1 for side in canonical)
        if canonical in injections:
            exact_weight = [Fraction(0) for _ in range(order + 1)]
            for degree, value in injections[canonical].items():
                if degree % 2 or not minimal_degree <= degree <= order:
                    raise AssertionError(
                        "injected degree contradicts the parity/minimal-degree lemma"
                    )
                exact_weight[degree] = Fraction(value)
            injected_seen.add(canonical)
        else:
            exact_weight = list(flm._ht_box_log(canonical, order))
            for subshape, subweight in weights.items():
                if all(sub <= side for sub, side in zip(subshape, shape)):
                    placements = prod(
                        side - sub + 1 for sub, side in zip(subshape, shape)
                    )
                    for degree in range(order + 1):
                        exact_weight[degree] -= placements * subweight[degree]
        frozen_weight = tuple(exact_weight)
        weights[shape] = frozen_weight
        for degree, value in enumerate(frozen_weight):
            interaction[degree] += value
    if injected_seen != set(injections):
        raise AssertionError("an injected walled shape is absent from the box family")
    coefficients = list(interaction)
    for degree in range(2, order + 1, 2):
        coefficients[degree] += Fraction(3, degree)
    return FLMSeries(
        dimension=3,
        order=order,
        variable="v",
        coefficients=tuple(coefficients),
        interaction_coefficients=tuple(interaction),
        boxes=shapes,
        box_weights=weights,
        bound_budget=order // 2,
        bound_slack=0,
        order_bound="2*sum_i(side_i-1) <= truncation_order",
    )


def _rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _sha256_json(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _record(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}", flush=True)


def run_frontier_worker(
    shape, special_cut: int | None, timeout_seconds: int
) -> dict:
    """Run one frontier extraction in an isolated process and parse its JSON."""

    command = [
        sys.executable,
        str(ROOT / SCRIPT),
        "--worker",
        "--shape",
        ",".join(str(side) for side in shape),
    ]
    if special_cut is not None:
        command += ["--special-cut", str(special_cut)]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=ROOT,
        timeout=timeout_seconds,
    )
    parent_wall = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"frontier worker failed: {completed.stderr.strip()[-2000:]}"
        )
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    payload["parent_wall_seconds"] = f"{parent_wall:.6f}"
    payload["worker_command"] = command
    return payload


def _box_edge_count(shape) -> int:
    return flm._box_bond_count(tuple(int(side) for side in shape))


def _engine_projection(shape) -> dict:
    sides = sorted(int(side) for side in shape)
    cross = sides[0] * sides[1]
    columns = _box_edge_count(sides) + 1
    return {
        "shape": list(sides),
        "cross_section_sites": cross,
        "broken_bond_degree": columns - 1,
        "simultaneous_int64_arrays": 3,
        "projected_peak_bytes": 3 * (1 << cross) * columns * 8,
    }


def _proof_text(data: dict) -> str:
    weights = data["weights"]
    series = data["series"]
    return f"""# Exact simple-cubic HT coefficient through v^26

## Scope

[COMPUTATION] This is a finite exact coefficient calculation, not a solution of
the three-dimensional Ising model.

## Which boxes contribute at order 26

[LEMMA] Enting's finite-lattice identity is
`L(A) = sum_(R<=A) product_i(A_i-R_i+1) W(R)`, where `L(A)` is the
free-box logarithm and `W(R)` is its Moebius-inverted exact-bounding-box
weight.  Every edge of an `a*b*c` box crosses exactly one of the
`a+b+c-3` coordinate cuts, and an even subgraph crosses each cut an even
number of times.  A connected wrapping contribution of a box therefore has at
least `2*(a+b+c-3)` edges.  Thus the cancellation/order bound is
`2*sum_i(side_i-1) <= truncation_order`; every canonical class with
`a+b+c > 16` has zero degree-26 weight, while every class with
`a+b+c <= 16` is included in the inversion.  The artifact enumerates all
{data['box_class_lemma']['contributing_canonical_class_count']} canonical
classes explicitly.  The {data['box_class_lemma']['minimal_span_class_count']}
classes with `a+b+c = 16` contribute their minimal-span weight at degree 26,
with all 13 cuts crossed exactly twice.  The classes
{data['box_class_lemma']['walled']} additionally have open cross-sections of
{data['box_class_lemma']['walled_cross_sections']} sites, which both
spin-transfer engines refuse (`_MAX_CROSS_SECTION = 22`); their dense arrays
would need {data['box_class_lemma']['walled_projected_bytes']} bytes.  The
degree-26 weight of (5,5,5) (12 cuts) is not minimal-span: the cut profile is
one cut crossed four times and the other eleven twice, since
`26 = 2*11 + 4`.

## Frontier algorithm and weight extraction

[LEMMA] The v^24 vertex frontier generalizes verbatim to arbitrary boxes:
process vertices in transfer-direction-major order, keep the pending parity
register of width `b*c`, and count cut crossings in ternary digits capped at
two.  With one square-free variable per cut, the full-support coefficient of
the formal logarithm of the {0,2}-capped polynomial is exactly the box's
minimal-span weight at degree 26: any proper subbox placement misses a cut of
the box, while every wrapping degree-26 contribution crosses all 13 cuts
exactly twice, because every edge crosses one cut and each cut total is even
and positive.

[LEMMA] For the non-minimal (5,5,5) weight at degree 26, cap one
distinguished cut at four and use ordinary half-crossing variables.  The
coefficient of `x_special^2 * product_(other cuts) x` in the truncated formal
logarithm collects exactly the wrapping contributions whose total profile is
`4` on the distinguished cut and `2` on the other eleven.  Ordinary
multiplication is essential: it includes both a single four-crossing factor
and a two-plus-two split across two logarithm factors.  The 12 cuts fall into
two reflection orbits (edge planes and middle planes, six cuts each), so the
weight is `6*(c_edge + c_middle)`.

## CRT certificates and extracted weights

[COMPUTATION] The exact degree-26 finite-lattice weights are
`(4,6,6): {weights['w_664']}`, `(5,5,6): {weights['w_655']}`, and
`(5,5,5) degree-26: {weights['w_555_26']}`.

[LEMMA] Put `O_n = sum_(k=1..n) S(n,k)(k-1)!`.  At minimal span, choosing
the two crossing edges on each cut bounds every formal-log term, giving
`O_13 * product_c C(E_c,2)`.  For the four-crossing calculation, label the
two half-crossing units of the distinguished cut: every ordered-log
composition maps to at least one ordered partition of 13 slots, and bounding
each special slot by `C(25,4)` gives a strict (deliberately conservative)
bound.  The exact formulas and values used here are
{data['crt']['bound_formulas']}.

Each run uses the modulus product `{data['crt']['modulus_product']}` of six
deterministically checked 31-bit primes.  The strict absolute bounds are
{data['crt']['per_weight_bounds']}; each satisfies
`modulus_product > 2*bound` and `|weight| <= bound`, hence every centered CRT
reconstruction is unique.

[COMPUTATION] Injecting only these three weights (plus the stored certified
v^24 weight 8163299968 for (5,5,5)) into the unchanged exact finite-lattice
inversion gives

`[v^26] phi = {series['v26']}`,

with interaction coefficient `{series['interaction_v26']}`.  Every coefficient
through v^24 reproduces the stored v^24 artifact exactly.

## Independent controls

[COMPUTATION] The generalized frontier reproduces the stored minimal weights
16 (2x2x2), 9188 (3x3x3), 8655072 (4x4x4), and 8163299968 (5x5x5 cap-2 part).
The 4-resolving method is validated at lower order: for 3x3x3 it returns
`6*c_4 = {data['controls']['w333_orbit_sum']}` against the ordinary
spin-transfer FLM value `{data['controls']['w333_engine']}` at degree 14, and
for 4x4x4 `6*c_edge + 3*c_middle = {data['controls']['w444_orbit_sum']}`
against `{data['controls']['w444_engine']}` at degree 20.  Minimal-span
weights of {data['controls']['battery_agreement_count']} of the 21 sum-16
classes agree
with the unchanged engine plus Moebius inversion at degree 26; the two
exceptions are the walled classes themselves.  The (5,5,6) weight was
additionally recomputed by a second sweep direction
(`{data['controls']['alt_sweep_shape']}`, cross-section
{data['controls']['alt_sweep_cross']}) with identical result
{data['controls']['alt_sweep_weight']}.  Disk checkpoint use was zero bytes.

[EXTERNAL] Table I of Arisue and Fujiwara, arXiv:hep-lat/0209002 (2002),
tabulates the same interaction series (their `f = 3 log cosh(beta) + sum_n
a_n t^n`) through order 46 by a layer-restricted anisotropic transfer
algorithm.  All 13 tabulated entries `a_2..a_26` agree with this computation,
including the previously unpublished-in-this-repository entry
`a_26 = 4437596650548`.  Source PDF sha256:
`bf79ec7918694768b57e6793ed4b0c8f581aaea6a705b4273ecab42641723afc`.
This is an independent falsification target, not an input: no coefficient
was fit to it.

## Resource telemetry

[COMPUTATION] Pre-launch projection from v^24 telemetry (204 bytes/state at
1.7M peak states): the probed peaks
{data['frontier']['projection']} stay under the
{MEMORY_BUDGET_BYTES} byte budget, so no additional conserved label was
needed.  Measured production walls and peaks:
{data['frontier']['measured']}.
The complete coefficient run, including the unchanged engine on the other
boxes, peaked at {data['frontier']['artifact_generation_peak_rss_bytes']} bytes.
"""


def run_full_experiment(*, write_artifacts: bool = True) -> dict:
    from e18_series_extend import _hybrid_flm_engine

    checks: list[dict] = []
    v24_payload = json.loads(V24_PATH.read_text(encoding="utf-8"))
    v24_series = v24_payload["data"]["series"]
    v24_coefficients = tuple(Fraction(s) for s in v24_series["coefficients"])
    stored_w555_24 = int(v24_payload["data"]["minimal_cube"]["exact_weight"])
    prefix_payload = json.loads(PREFIX_PATH.read_text(encoding="utf-8"))
    prefix = tuple(Fraction(s) for s in prefix_payload["data"]["coefficients"])
    started = time.perf_counter()
    rss_before = _rss_bytes()

    # ---- cheap frontier controls on known minimal-span weights -------------
    cube_weights = {}
    for side in (2, 3, 4):
        payload = extract_frontier_weights((side, side, side))
        cube_weights[side] = payload["minimal_centered"]
    expected_cubes = {2: 16, 3: 9188, 4: 8_655_072}

    # ---- 4-resolving method validation at lower order ----------------------
    with _hybrid_flm_engine():
        series333 = flm.high_temperature_free_energy(3, 14)
        series444 = flm.high_temperature_free_energy(3, 20)
    w333_engine = series333.box_weights[(3, 3, 3)][14]
    w444_engine = series444.box_weights[(4, 4, 4)][20]
    cap4_333 = [extract_frontier_weights((3, 3, 3), special_cut=k) for k in (0, 1)]
    w333_orbit = 6 * cap4_333[0]["four_cut_centered"]
    cap4_444 = [extract_frontier_weights((4, 4, 4), special_cut=k) for k in (0, 1)]
    w444_orbit = (
        6 * cap4_444[0]["four_cut_centered"]
        + 3 * cap4_444[1]["four_cut_centered"]
    )

    # ---- minimal-span battery over every sum-16 class ----------------------
    classes = minimal_span_classes(ORDER)
    battery = {}
    for shape in classes:
        if shape in ((6, 5, 5), (6, 6, 4)):
            continue  # walled: computed below in isolated workers
        payload = extract_frontier_weights(shape)
        battery[shape] = payload

    # ---- memory projection before any big launch ---------------------------
    # The first axis is the transfer direction.  Every actual full-modulus
    # worker is first run modulo one prime, then projected at 220 bytes/state
    # from the v24 telemetry plus the larger six-prime value objects.
    projection = {}
    projection_specs = (
        ((6, 5, 5), None),
        ((6, 6, 4), None),
        *(((5, 5, 5), cut) for cut in SPECIAL_CUTS_555),
        ((5, 6, 5), None),
    )
    for shape, special in projection_specs:
        probe = frontier_cut_polynomial(
            shape,
            modulus=CRT_PRIMES[0],
            max_states=20_000_000,
            special_cut=special,
        )
        projected_bytes = probe["peak_states"] * 220
        projection[
            str(list(shape))
            + (f"/special={special}" if special is not None else "")
        ] = {
            "probe_completed": probe["completed"],
            "lower_bound_only": not probe["completed"],
            "probe_peak_states": probe["peak_states"],
            "projected_bytes_at_full_modulus": projected_bytes,
            "within_budget": projected_bytes < MEMORY_BUDGET_BYTES,
        }
    if not all(
        row["probe_completed"] and row["within_budget"]
        for row in projection.values()
    ):
        raise RuntimeError(
            "pre-launch frontier projection reached the state cap or memory budget; "
            "split by an additional conserved label before production"
        )

    # ---- production frontier runs in isolated workers ----------------------
    w655_payload = run_frontier_worker((6, 5, 5), None, timeout_seconds=5400)
    w664_payload = run_frontier_worker((6, 6, 4), None, timeout_seconds=5400)
    special_payloads = [
        run_frontier_worker((5, 5, 5), cut, timeout_seconds=5400)
        for cut in SPECIAL_CUTS_555
    ]
    alt655_payload = run_frontier_worker((5, 6, 5), None, timeout_seconds=7200)

    w655 = w655_payload["minimal_centered"]
    w664 = w664_payload["minimal_centered"]
    edge_coeff = special_payloads[0]["four_cut_centered"]
    middle_coeffs = [
        payload["four_cut_centered"] for payload in special_payloads[1:]
    ]
    w555_26 = 6 * (edge_coeff + middle_coeffs[0])
    w555_24_control = special_payloads[0]["minimal_centered"]

    injections = {
        (5, 5, 5): {24: stored_w555_24, 26: w555_26},
        (5, 5, 6): {26: w655},
        (4, 6, 6): {26: w664},
    }

    # ---- series assembly with the unchanged engine on all other boxes ------
    with _hybrid_flm_engine() as transfer_records:
        series = high_temperature_series_with_injections(ORDER, injections)
        independent_controls = high_temperature_series_with_injections(
            ORDER, injections
        )

    battery_agreements = {}
    for shape, payload in battery.items():
        canonical = flm._canonical_shape(shape)
        engine_value = series.box_weights[canonical][ORDER]
        battery_agreements[shape] = {
            "frontier": payload["minimal_centered"],
            "engine": engine_value,
            "agree": payload["minimal_centered"] == engine_value,
        }

    primes_ok = (
        len(set(CRT_PRIMES)) == len(CRT_PRIMES)
        and all(_is_prime_32(prime) for prime in CRT_PRIMES)
        and all(
            gcd(CRT_PRIMES[left], CRT_PRIMES[right]) == 1
            for left in range(len(CRT_PRIMES))
            for right in range(left)
        )
    )
    crt_weights = {
        "w_655": (w655, int(w655_payload["bound"])),
        "w_664": (w664, int(w664_payload["bound"])),
        "w_555_26": (w555_26, 12 * int(special_payloads[0]["bound"])),
        "w_555_26_edge_coeff": (edge_coeff, int(special_payloads[0]["bound"])),
    }
    crt_unique = all(
        CRT_PRODUCT > 2 * bound and abs(weight) <= bound
        for weight, bound in crt_weights.values()
    )

    prefix_match = series.coefficients[: len(v24_coefficients)] == v24_coefficients
    extended_match = series.coefficients[: len(prefix)] == prefix

    canonical_contributing = tuple(
        sorted(
            {
                flm._canonical_shape(shape)
                for shape in flm._ht_shapes(3, ORDER, 0)
            }
        )
    )
    _record(
        checks,
        "prelaunch_frontier_projection",
        all(
            row["probe_completed"] and row["within_budget"]
            for row in projection.values()
        ),
        f"{len(projection)} one-prime probes completed; largest projection "
        f"{max(row['projected_bytes_at_full_modulus'] for row in projection.values())} "
        f"bytes is below {MEMORY_BUDGET_BYTES}",
    )
    _record(
        checks,
        "independent_minimal_cube_controls",
        all(cube_weights[side] == value for side, value in expected_cubes.items()),
        "frontier/log weights 16, 9188, 8655072 agree with stored values",
    )
    _record(
        checks,
        "cap4_method_control_3x3x3",
        cap4_333[0]["four_cut_centered"] == cap4_333[1]["four_cut_centered"]
        and w333_orbit == w333_engine,
        f"6*c_4={w333_orbit} equals spin-transfer FLM weight {w333_engine} at v^14",
    )
    _record(
        checks,
        "cap4_method_control_4x4x4",
        w444_orbit == w444_engine,
        f"6*c_edge+3*c_middle={w444_orbit} equals FLM weight {w444_engine} at v^20",
    )
    _record(
        checks,
        "minimal_span_battery_vs_engine",
        all(row["agree"] for row in battery_agreements.values()),
        f"{sum(row['agree'] for row in battery_agreements.values())} of "
        f"{len(battery_agreements)} non-walled sum-16 classes agree with "
        "the engine+Moebius weights at degree 26",
    )
    _record(
        checks,
        "v24_minimal_weight_control",
        w555_24_control == stored_w555_24,
        f"cap-4 run reproduces the stored v^24 weight {stored_w555_24}",
    )
    _record(
        checks,
        "w26_middle_orbit_representatives",
        len(set(middle_coeffs)) == 1,
        f"middle-plane coefficients identical across axes: {middle_coeffs}",
    )
    _record(
        checks,
        "alt_sweep_direction_w655",
        alt655_payload["minimal_centered"] == w655,
        "transfer-direction charge resweep of (5,5,6) reproduces "
        f"{w655} with peak {alt655_payload['peak_states']} states",
    )
    _record(
        checks,
        "crt_moduli_are_distinct_primes",
        primes_ok and prod(CRT_PRIMES) == CRT_PRODUCT,
        "deterministic 32-bit Miller-Rabin and all pairwise gcd checks pass",
    )
    _record(
        checks,
        "crt_unique_centered_reconstruction",
        crt_unique,
        f"product {CRT_PRODUCT} exceeds twice every per-weight bound",
    )
    _record(
        checks,
        "canonical_prefix_through_v24",
        prefix_match and extended_match,
        "all 25 coefficients v^0..v^24 reproduce the stored v^24 artifact",
    )
    _record(
        checks,
        "odd_degree_vanishes",
        series.coefficients[25] == 0 and series.interaction_coefficients[25] == 0,
        "bipartite boxes force zero odd-degree coefficients",
    )
    _record(
        checks,
        "v26_coefficients",
        series.coefficients[26].denominator == 26,
        f"[v^26] phi = {series.coefficients[26]}; interaction "
        f"{series.interaction_coefficients[26]}",
    )
    _record(
        checks,
        "deterministic_repeated_injection",
        independent_controls.coefficients == series.coefficients,
        "a repeated exact finite-lattice inversion gives identical coefficients",
    )
    external_table = EXTERNAL_ARISUE_FUJIWARA["table_i_interaction"]
    external_match = all(
        series.interaction_coefficients[degree] == Fraction(value)
        for degree, value in external_table.items()
    )
    _record(
        checks,
        "external_arisue_fujiwara_table_i",
        external_match,
        "all 13 tabulated interaction coefficients v^2..v^26 of "
        "arXiv:hep-lat/0209002 Table I agree, including a_26=4437596650548",
    )

    run_usage = resource.getrusage(resource.RUSAGE_SELF)
    run_rss = _rss_bytes()
    run_elapsed = time.perf_counter() - started

    data = {
        "claim_tag": "[COMPUTATION]",
        "scope": "finite exact simple-cubic high-temperature series coefficient; not a solution of the 3D Ising model",
        "old_engine": {
            "max_supported_cross_section_sites": 22,
            "walled": {
                str(list(shape)): _engine_projection(shape)
                for shape in WALLED_SHAPES
            },
        },
        "box_class_lemma": {
            "order": ORDER,
            "contributing_classes": "a+b+c <= 16",
            "contributing_ordered_box_count": len(flm._ht_shapes(3, ORDER, 0)),
            "contributing_canonical_class_count": len(canonical_contributing),
            "contributing_canonical_classes": [
                list(shape) for shape in canonical_contributing
            ],
            "minimal_span_class_count": len(classes),
            "minimal_span_classes": [list(shape) for shape in classes],
            "walled": [list(shape) for shape in WALLED_SHAPES],
            "walled_cross_sections": [25, 25, 24],
            "walled_projected_bytes": [
                _engine_projection(shape)["projected_peak_bytes"]
                for shape in WALLED_SHAPES
            ],
        },
        "algorithm": {
            "name": "box-general cut-capped parity frontier plus ordinary half-crossing logarithm",
            "cut_count_minimal_span": 13,
            "cut_count_555_degree26": 12,
            "counter_alphabet": [0, 1, 2],
            "special_counter_alphabet": [0, 1, 2, 3, 4],
            "four_crossing_log": "ordinary mixed-radix exponents retain a 2+2 split across logarithm factors",
            "reason_exact_at_v26": "every edge crosses one cut; 13 positive even cut totals sum to 26, hence exactly two each; for (5,5,5) at degree 26 one cut totals four",
            "disk_checkpoint_bytes": 0,
        },
        "frontier": {
            "projection": projection,
            "projection_note": "v24 telemetry: 352305152 bytes / 1727071 states ~ 204 bytes per state; full-modulus values add ~16 bytes",
            "measured": {
                "w655_worker": {
                    **{
                        key: w655_payload[key]
                        for key in (
                            "peak_states",
                            "terminal_states",
                            "wall_seconds",
                            "parent_wall_seconds",
                            "max_rss_bytes",
                        )
                    },
                    "state_count_sha256": _sha256_json(
                        [row["states"] for row in w655_payload["timeline"]]
                    ),
                    "timeline": w655_payload["timeline"],
                },
                "w664_worker": {
                    **{
                        key: w664_payload[key]
                        for key in (
                            "peak_states",
                            "terminal_states",
                            "wall_seconds",
                            "parent_wall_seconds",
                            "max_rss_bytes",
                        )
                    },
                    "state_count_sha256": _sha256_json(
                        [row["states"] for row in w664_payload["timeline"]]
                    ),
                    "timeline": w664_payload["timeline"],
                },
                "special_555_workers": [
                    {
                        "special_cut": payload["special_cut"],
                        "peak_states": payload["peak_states"],
                        "terminal_states": payload["terminal_states"],
                        "wall_seconds": payload["wall_seconds"],
                        "parent_wall_seconds": payload["parent_wall_seconds"],
                        "max_rss_bytes": payload["max_rss_bytes"],
                        "state_count_sha256": _sha256_json(
                            [row["states"] for row in payload["timeline"]]
                        ),
                        "timeline": payload["timeline"],
                    }
                    for payload in special_payloads
                ],
                "alt_sweep_worker": {
                    **{
                        key: alt655_payload[key]
                        for key in (
                            "peak_states",
                            "terminal_states",
                            "wall_seconds",
                            "parent_wall_seconds",
                            "max_rss_bytes",
                        )
                    },
                    "state_count_sha256": _sha256_json(
                        [row["states"] for row in alt655_payload["timeline"]]
                    ),
                    "timeline": alt655_payload["timeline"],
                },
            },
            "artifact_generation_peak_rss_bytes": run_rss,
            "artifact_generation_rss_before_bytes": rss_before,
            "artifact_generation_elapsed_seconds": f"{run_elapsed:.6f}",
        },
        "crt": {
            "primes": list(CRT_PRIMES),
            "pairwise_coprime": primes_ok,
            "deterministic_primality_checked": primes_ok,
            "modulus_product": str(CRT_PRODUCT),
            "modulus_product_bits": CRT_PRODUCT.bit_length(),
            "centered_reconstruction_condition": "modulus_product > 2*absolute_bound",
            "bound_formulas": {
                "O_13": "sum_{k=1}^{13} S(13,k)*(k-1)!",
                "O_13_value": str(ordered_partition_factor(13)),
                "w_655": "O_13*C(25,2)^5*C(30,2)^8",
                "w_664": "O_13*C(24,2)^10*C(36,2)^3",
                "w_555_26_coefficient": "O_13*C(25,4)^2*C(25,2)^11",
                "w_555_26": "12*O_13*C(25,4)^2*C(25,2)^11",
            },
            "per_weight_bounds": {
                name: {"weight": str(weight), "bound": str(bound)}
                for name, (weight, bound) in crt_weights.items()
            },
        },
        "weights": {
            "w_664": str(w664),
            "w_655": str(w655),
            "w_555_26": str(w555_26),
            "w_555_26_edge_coefficient": str(edge_coeff),
            "w_555_26_middle_coefficient": str(middle_coeffs[0]),
            "w_555_24_control": str(w555_24_control),
            "injections": {
                str(list(shape)): {str(k): str(v) for k, v in degrees.items()}
                for shape, degrees in injections.items()
            },
            "residues": {
                "w_655": w655_payload["minimal_residue"],
                "w_664": w664_payload["minimal_residue"],
                "w_555_26_edge": special_payloads[0]["four_cut_residue"],
                "w_555_26_middle": special_payloads[1]["four_cut_residue"],
            },
        },
        "series": {
            "dimension": 3,
            "variable": "v",
            "achieved_order": ORDER,
            "normalization": "phi = log(2) + sum_n coefficients[n] v^n",
            "coefficients": [str(value) for value in series.coefficients],
            "interaction_coefficients": [
                str(value) for value in series.interaction_coefficients
            ],
            "v26": str(series.coefficients[26]),
            "interaction_v26": str(series.interaction_coefficients[26]),
            "v24_artifact_path": str(V24_PATH.relative_to(ROOT)),
            "v24_artifact_sha256": hashlib.sha256(V24_PATH.read_bytes()).hexdigest(),
            "canonical_prefix_path": str(PREFIX_PATH.relative_to(ROOT)),
            "canonical_prefix_sha256": hashlib.sha256(
                PREFIX_PATH.read_bytes()
            ).hexdigest(),
            "transfer_record_count": len(transfer_records),
            "transfer_records": [
                record
                for record in transfer_records.values()
                if record["cross_section_sites"] >= 16
            ],
        },
        "production_resource_measurement": {
            "command": ".venv/bin/python experiments/e86_ht_v26.py",
            "measurement_method": "resource.getrusage(RUSAGE_SELF) in the producing process",
            "maximum_resident_set_bytes": run_rss,
            "real_seconds": f"{run_elapsed:.6f}",
            "disk_input_blocks_before_artifact_write": int(run_usage.ru_inblock),
            "disk_output_blocks_before_artifact_write": int(run_usage.ru_oublock),
        },
        "controls": {
            "side2_weight": str(cube_weights[2]),
            "side3_weight": str(cube_weights[3]),
            "side4_weight": str(cube_weights[4]),
            "w333_engine": str(w333_engine),
            "w333_orbit_sum": str(w333_orbit),
            "w444_engine": str(w444_engine),
            "w444_orbit_sum": str(w444_orbit),
            "battery": {
                str(list(shape)): {
                    "frontier": str(row["frontier"]),
                    "engine": str(row["engine"]),
                    "agree": row["agree"],
                }
                for shape, row in battery_agreements.items()
            },
            "battery_agreement_count": sum(
                row["agree"] for row in battery_agreements.values()
            ),
            "battery_class_count": len(battery_agreements),
            "alt_sweep_shape": [5, 6, 5],
            "alt_sweep_cross": 30,
            "alt_sweep_weight": str(alt655_payload["minimal_centered"]),
            "method": "box-general vertex frontier/log extractor; cross-checked against broken-bond spin-transfer FLM",
            "external_arisue_fujiwara": {
                "tag": "[EXTERNAL]",
                "arxiv": EXTERNAL_ARISUE_FUJIWARA["arxiv"],
                "table_i_sha256": EXTERNAL_ARISUE_FUJIWARA["table_i_sha256"],
                "agreement": external_match,
                "relation": "their f = 3 log cosh(beta) + sum a_n t^n matches this interaction series exactly for n = 2..26",
            },
        },
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": INTERPRETER,
            "arithmetic": "exact Python integers/Fractions and CRT residues",
            "benchmark_kc_used": False,
        },
        "data": data,
        "checks": checks,
    }
    if write_artifacts:
        RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        PROOF_PATH.write_text(_proof_text(data), encoding="utf-8")
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more embedded checks failed")
    return payload


def _worker_main(args: argparse.Namespace) -> None:
    shape = tuple(int(side) for side in args.shape.split(","))
    payload = extract_frontier_weights(
        shape, special_cut=args.special_cut, max_states=args.max_states
    )
    payload["max_rss_bytes"] = _rss_bytes()
    print(json.dumps(payload))


def _profile_main(args: argparse.Namespace) -> None:
    shape = tuple(int(side) for side in args.shape.split(","))
    modulus = CRT_PRODUCT if args.full_modulus else CRT_PRIMES[0]
    result = frontier_cut_polynomial(
        shape,
        modulus=modulus,
        max_states=args.max_states,
        special_cut=args.special_cut,
    )
    print(
        f"shape={result['shape']} special={result['special_cut']} "
        f"completed={result['completed']} peak_states={result['peak_states']} "
        f"terminal_states={result['terminal_states']}"
    )
    if result["completed"] and args.extract_log:
        payload = extract_frontier_weights(
            shape, modulus=modulus, special_cut=args.special_cut
        )
        print(f"minimal_centered={payload['minimal_centered']}")
        if args.special_cut is not None:
            print(f"four_cut_centered={payload['four_cut_centered']}")
        print(f"bound={payload['bound']} certified={modulus > 2 * int(payload['bound'])}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--shape", type=str)
    parser.add_argument("--special-cut", type=int)
    parser.add_argument("--max-states", type=int)
    parser.add_argument("--full-modulus", action="store_true")
    parser.add_argument("--extract-log", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    if args.worker:
        _worker_main(args)
        return
    if args.profile:
        _profile_main(args)
        return
    run_full_experiment(write_artifacts=not args.no_write)
    print("PASS")


if __name__ == "__main__":
    main()
