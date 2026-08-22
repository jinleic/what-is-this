"""Independent exact/numerical falsification tests for proposed 3-D Ising solutions.

Run from the repository root with

    .venv/bin/python experiments/e06_falsify_claims.py

The combinatorial calculations use integer or :class:`fractions.Fraction` arithmetic.  All
transcendental calculations use mpmath at the precision recorded in the result artifact.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Callable, Sequence

import mpmath as mp
import numpy as np
import sympy as sp

from ising.transfer_matrix import (
    box_broken_bond_poly,
    layer_bonds,
    torus_broken_bond_poly,
)

PRECISION = 80
SCRIPT = "experiments/e06_falsify_claims.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "falsification" / "falsification.json"

mp.mp.dps = PRECISION


def _mpstr(value: mp.mpf | mp.mpc, digits: int = PRECISION - 10) -> str:
    """Stable decimal serialization without converting through binary float."""

    value = mp.mpf(value)
    if value == 0:
        return "0"
    return mp.nstr(value, digits)


def _fraction_str(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


# ---------------------------------------------------------------------------
# Critical-condition tests
# ---------------------------------------------------------------------------


def _positive_critical_root(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """Solve sinh(a beta) sinh(b beta) = 1 on the unique positive branch."""

    a, b = mp.mpf(a), mp.mpf(b)
    return mp.findroot(lambda beta: mp.sinh(a * beta) * mp.sinh(b * beta) - 1, (mp.mpf("0.1"), mp.mpf("0.5")))


def zhang_critical_data() -> dict[str, str]:
    """Recompute the 2007 claimed isotropic critical point and benchmark discrepancy."""

    root = _positive_critical_root(2, 6)
    x_golden = (mp.sqrt(5) - 1) / 2
    closed_form = -mp.log(x_golden) / 2
    benchmark = mp.mpf("0.221654626")
    benchmark_sigma = mp.mpf("0.000000005")
    discrepancy = abs(root - benchmark)
    return {
        "equation_root": _mpstr(root),
        "closed_form": _mpstr(closed_form),
        "x_c": _mpstr(x_golden),
        "equation_residual": _mpstr(mp.sinh(2 * root) * mp.sinh(6 * root) - 1),
        "root_closed_form_error": _mpstr(abs(root - closed_form)),
        "benchmark": _mpstr(benchmark),
        "benchmark_one_sigma": _mpstr(benchmark_sigma),
        "absolute_discrepancy": _mpstr(discrepancy),
        "relative_discrepancy": _mpstr(discrepancy / benchmark),
        "standard_deviations": _mpstr(discrepancy / benchmark_sigma),
    }


def degang_critical_data() -> dict[str, object]:
    """Apply the claimed anisotropic equation to cyclic permutations of fixed couplings."""

    couplings = (mp.mpf("0.5"), mp.mpf("1"), mp.mpf("1.5"))
    rows: list[dict[str, object]] = []
    roots: list[mp.mpf] = []
    for cycle in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        j, j1, j2 = (couplings[index] for index in cycle)
        root = _positive_critical_root(2 * j, 2 * (j1 + j2))
        roots.append(root)
        rows.append(
            {
                "role_assignment": {"J": _mpstr(j), "J1": _mpstr(j1), "J2": _mpstr(j2)},
                "beta_c": _mpstr(root),
                "equation_residual": _mpstr(
                    mp.sinh(2 * root * j) * mp.sinh(2 * root * (j1 + j2)) - 1
                ),
            }
        )

    isotropic = _positive_critical_root(2, 4)
    benchmark = mp.mpf("0.221654626")
    return {
        "input_couplings": [_mpstr(coupling) for coupling in couplings],
        "cyclic_permutations": rows,
        "beta_c_range": _mpstr(max(roots) - min(roots)),
        "rotational_invariance_passed": max(roots) == min(roots),
        "isotropic_beta_c": _mpstr(isotropic),
        "isotropic_benchmark": _mpstr(benchmark),
        "isotropic_absolute_discrepancy": _mpstr(abs(isotropic - benchmark)),
        "isotropic_relative_discrepancy": _mpstr(abs(isotropic - benchmark) / benchmark),
    }


# ---------------------------------------------------------------------------
# Exact low-temperature magnetisation with a transfer-matrix field axis
# ---------------------------------------------------------------------------


def _broken_in_layer(ns: int, bonds: Sequence[tuple[int, int]]) -> np.ndarray:
    states = np.arange(1 << ns, dtype=np.uint32)
    result = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for left, right in bonds:
        result += (
            ((states >> np.uint32(left)) ^ (states >> np.uint32(right))) & one
        ).astype(np.int64)
    return result


def _weighted_down_count(ns: int, weights: Sequence[int]) -> np.ndarray:
    states = np.arange(1 << ns, dtype=np.uint32)
    result = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for index, weight in enumerate(weights):
        if weight:
            result += int(weight) * (
                1 - ((states >> np.uint32(index)) & one)
            ).astype(np.int64)
    return result


def _shift_degree_rows(array: np.ndarray, shifts: np.ndarray) -> np.ndarray:
    """Multiply each state row by x**shift; the broken-bond degree is the last axis."""

    degree = array.shape[-1] - 1
    result = np.zeros_like(array)
    for shift_value in np.unique(shifts):
        shift = int(shift_value)
        selected = shifts == shift
        if shift == 0:
            result[selected] = array[selected]
        elif shift <= degree:
            result[selected, ..., shift:] = array[selected, ..., : degree + 1 - shift]
    return result


def _apply_interlayer(array: np.ndarray, ns: int, antiperiodic: bool = False) -> np.ndarray:
    """Apply all single-site interlayer bond matrices with exact polynomial shifts."""

    tail = array.shape[1:]
    degree = array.shape[-1] - 1
    value = array
    for bit in range(ns):
        high = 1 << (ns - bit - 1)
        low = 1 << bit
        value = value.reshape((high, 2, low) + tail)
        if antiperiodic:
            # Favoured states differ: [[x,1],[1,x]].
            result = value[:, ::-1].copy()
            result[..., 1:] += value[..., :degree]
        else:
            # Favoured states agree: [[1,x],[x,1]].
            result = value.copy()
            result[..., 1:] += value[:, ::-1][..., :degree]
        value = result.reshape((1 << ns,) + tail)
    return value


def _add_new_layer_down_count(array: np.ndarray, down_count: np.ndarray) -> np.ndarray:
    result = np.zeros_like(array)
    max_count = array.shape[1] - 1
    for state, count_value in enumerate(down_count):
        count = int(count_value)
        result[state, count:, :] = array[state, : max_count + 1 - count, :]
    return result


def box_broken_bond_magnetization_poly(
    shape: Sequence[int], periodic: Sequence[bool] | None = None
) -> list[list[int]]:
    """Return exact counts ``c[q][n_down]`` for a box with frozen-plus exterior.

    This is the same layer propagation as ``box_broken_bond_poly(..., plus_boundary=True)`` with
    one additional, exactly shifted axis for the number of down spins.  Summing that axis is
    asserted to reproduce the stable transfer-matrix engine.  The transfer depth must be at least
    two; callers may permute isotropic box axes to meet that condition.
    """

    shape = tuple(int(side) for side in shape)
    if len(shape) not in (2, 3) or any(side < 1 for side in shape):
        raise ValueError("shape must be a positive 2-D or 3-D box")
    dimension = len(shape)
    cross, depth = shape[:-1], shape[-1]
    if depth < 2:
        raise ValueError("transfer depth must be at least two for plus-boundary propagation")
    if periodic is None:
        periodic_cross = (False,) * (dimension - 1)
    else:
        periodic_cross = tuple(bool(value) for value in periodic)
        if len(periodic_cross) != dimension - 1:
            raise ValueError("periodic must have one entry per cross-section direction")

    ns = math.prod(cross)
    n_sites = ns * depth
    if n_sites > 62:
        raise ValueError("more than 62 sites risks int64 coefficient overflow")

    bonds = layer_bonds(cross, periodic_cross)
    physical_bonds = len(bonds) * depth + ns * (depth - 1)
    inplane_coordination = 2 * (dimension - 1)
    inplane_degrees = [0] * ns
    for left, right in bonds:
        inplane_degrees[left] += 1
        inplane_degrees[right] += 1
    ghost_weights = [inplane_coordination - degree for degree in inplane_degrees]
    if any(weight < 0 for weight in ghost_weights):
        raise AssertionError("negative ghost-bond multiplicity")
    ghost_bonds = sum(ghost_weights) * depth + 2 * ns
    max_degree = physical_bonds + ghost_bonds

    broken_layer = _broken_in_layer(ns, bonds)
    ghost_inplane = _weighted_down_count(ns, ghost_weights)
    layer_down = _weighted_down_count(ns, [1] * ns)

    vector = np.zeros((1 << ns, n_sites + 1, max_degree + 1), dtype=np.int64)
    for state, count_value in enumerate(layer_down):
        vector[state, int(count_value), 0] = 1
    vector = _shift_degree_rows(vector, broken_layer + ghost_inplane + layer_down)

    for layer_index in range(1, depth):
        vector = _apply_interlayer(vector, ns)
        vector = _add_new_layer_down_count(vector, layer_down)
        outward_ghosts = layer_down if layer_index == depth - 1 else 0
        vector = _shift_degree_rows(vector, broken_layer + ghost_inplane + outward_ghosts)

    total = vector.sum(axis=0).T
    if total.min() < 0:
        raise AssertionError("int64 overflow")
    if int(total.sum()) != 1 << n_sites:
        raise AssertionError("configuration count is not 2**N")
    result = [[int(value) for value in row] for row in total]
    while len(result) > 1 and not any(result[-1]):
        result.pop()

    marginal = [sum(row) for row in result]
    expected = box_broken_bond_poly(shape, periodic=periodic_cross, plus_boundary=True)
    while len(marginal) > 1 and marginal[-1] == 0:
        marginal.pop()
    if marginal != expected:
        raise AssertionError("magnetisation-axis marginal disagrees with transfer-matrix engine")
    return result


def _series_ratio(
    numerator: Sequence[int | Fraction], denominator: Sequence[int | Fraction], order: int
) -> list[Fraction]:
    if not denominator or denominator[0] == 0:
        raise ValueError("series denominator must have nonzero constant term")
    result = [Fraction(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        known = sum(
            Fraction(denominator[index]) * result[degree - index]
            for index in range(1, degree + 1)
            if index < len(denominator)
        )
        source = Fraction(numerator[degree]) if degree < len(numerator) else Fraction(0)
        result[degree] = (source - known) / Fraction(denominator[0])
    return result


def low_temperature_magnetization_series(max_u_order: int = 6) -> list[Fraction]:
    """Compute the infinite-volume plus-phase magnetisation through ``u**max_u_order``.

    For orders through ``x**12 = u**6``, every connected polymer or connected Mayer cluster has
    bounding-box side at most two.  The mixed backward difference at side three, using the eight
    boxes with side lengths two or three, contains all such weights and removes surface and edge
    translations.  No literature coefficient enters this calculation.
    """

    if not 0 <= max_u_order <= 6:
        raise ValueError("the proved finite-lattice exactness bound is max_u_order <= 6")
    x_order = 2 * max_u_order
    down_density_by_box: dict[tuple[int, int, int], list[Fraction]] = {}
    for shape in product((2, 3), repeat=3):
        joint = box_broken_bond_magnetization_poly(tuple(sorted(shape)))
        denominator = [sum(row) for row in joint]
        numerator = [sum(count * multiplicity for count, multiplicity in enumerate(row)) for row in joint]
        down_density_by_box[shape] = _series_ratio(numerator, denominator, x_order)

    bulk_density = [Fraction(0) for _ in range(x_order + 1)]
    for shape, density in down_density_by_box.items():
        sign = -1 if sum(side == 2 for side in shape) % 2 else 1
        for degree, coefficient in enumerate(density):
            bulk_density[degree] += sign * coefficient

    magnetization_x = [-2 * coefficient for coefficient in bulk_density]
    magnetization_x[0] += 1
    if any(magnetization_x[degree] for degree in range(1, x_order + 1, 2)):
        raise AssertionError("odd powers of x cannot be expressed as an integer u series")
    return [magnetization_x[2 * degree] for degree in range(max_u_order + 1)]


# ---------------------------------------------------------------------------
# Exact high-temperature finite-lattice references
# ---------------------------------------------------------------------------


def _broken_to_even_subgraph(
    broken: Sequence[int], n_sites: int, n_bonds: int, max_order: int | None = None
) -> tuple[int, ...]:
    """Transform the broken-bond polynomial to the exact high-T polynomial P(v)."""

    if max_order is None:
        max_order = n_bonds
    result: list[int] = []
    normalization = 1 << n_sites
    for degree in range(max_order + 1):
        total = 0
        for broken_count, multiplicity in enumerate(broken):
            if not multiplicity:
                continue
            coefficient = sum(
                (-1) ** left_degree
                * math.comb(broken_count, left_degree)
                * math.comb(n_bonds - broken_count, degree - left_degree)
                for left_degree in range(
                    max(0, degree - (n_bonds - broken_count)),
                    min(broken_count, degree) + 1,
                )
            )
            total += int(multiplicity) * coefficient
        quotient, remainder = divmod(total, normalization)
        if remainder:
            raise AssertionError("broken-bond to high-T transform was not integral")
        result.append(quotient)
    return tuple(result)


def _twisted_layer_broken(L: int, twist_x: bool, twist_y: bool) -> np.ndarray:
    ns = L * L
    states = np.arange(1 << ns, dtype=np.uint32)
    result = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for x_coord in range(L):
        for y_coord in range(L):
            left = x_coord * L + y_coord
            for axis, twisted in ((0, twist_x), (1, twist_y)):
                right_x = (x_coord + 1) % L if axis == 0 else x_coord
                right_y = (y_coord + 1) % L if axis == 1 else y_coord
                right = right_x * L + right_y
                difference = (
                    ((states >> np.uint32(left)) ^ (states >> np.uint32(right))) & one
                ).astype(np.int64)
                seam = (axis == 0 and x_coord == L - 1) or (
                    axis == 1 and y_coord == L - 1
                )
                result += 1 - difference if seam and twisted else difference
    return result


@lru_cache(maxsize=None)
def _twisted_torus_broken_poly(L: int, twists: tuple[bool, bool, bool]) -> tuple[int, ...]:
    """Exact torus broken-bond polynomial with antiferromagnetic seam twists."""

    if L != 3:
        raise ValueError("the exact winding-sector implementation is certified for L=3")
    ns = L * L
    n_sites = L**3
    n_bonds = 3 * n_sites
    dimension = 1 << ns
    block = 16
    broken_layer = _twisted_layer_broken(L, twists[0], twists[1])
    trace = np.zeros(n_bonds + 1, dtype=np.int64)
    for block_start in range(0, dimension, block):
        block_end = min(block_start + block, dimension)
        vector = np.zeros(
            (dimension, block_end - block_start, n_bonds + 1), dtype=np.int64
        )
        for state in range(block_start, block_end):
            vector[state, state - block_start, 0] = 1
        for z_coord in range(L):
            vector = _shift_degree_rows(vector, broken_layer)
            vector = _apply_interlayer(
                vector,
                ns,
                antiperiodic=bool(twists[2] and z_coord == L - 1),
            )
        for state in range(block_start, block_end):
            trace += vector[state, state - block_start]
    if trace.min() < 0 or int(trace.sum()) != 1 << n_sites:
        raise AssertionError("twisted transfer matrix overflow or state-count error")
    return tuple(int(value) for value in trace)


@lru_cache(maxsize=None)
def winding_projected_torus_polynomial(L: int = 3) -> tuple[int, ...]:
    """Return P(v) in the zero-winding sector of an exact ``L**3`` torus.

    Averaging all eight periodic/antiperiodic seam choices projects onto even subgraphs with zero
    Z2 winding.  A raw torus has a spurious winding loop at order L; after projection, the first
    finite-size graph is a pair of winding loops at order 2L.  Consequently coefficients below
    2L are the infinite-lattice coefficients.
    """

    n_sites = L**3
    n_bonds = 3 * n_sites
    untwisted = list(_twisted_torus_broken_poly(L, (False, False, False)))
    while len(untwisted) > 1 and untwisted[-1] == 0:
        untwisted.pop()
    if untwisted != torus_broken_bond_poly((L, L, L), block=16):
        raise AssertionError("custom untwisted torus propagation disagrees with stable engine")

    signed_polynomials = []
    for twists in product((False, True), repeat=3):
        signed_polynomials.append(
            _broken_to_even_subgraph(
                _twisted_torus_broken_poly(L, twists), n_sites, n_bonds
            )
        )
    projected = []
    for degree in range(n_bonds + 1):
        quotient, remainder = divmod(
            sum(polynomial[degree] for polynomial in signed_polynomials), 8
        )
        if remainder:
            raise AssertionError("winding projection was not integral")
        projected.append(quotient)
    return tuple(projected)


def _formal_log(polynomial: Sequence[int | Fraction], order: int) -> list[Fraction]:
    if not polynomial or polynomial[0] != 1:
        raise ValueError("formal logarithm requires constant term one")
    result = [Fraction(0) for _ in range(order + 1)]
    for degree in range(1, order + 1):
        correction = sum(
            Fraction(index) * result[index] * Fraction(polynomial[degree - index])
            for index in range(1, degree)
            if degree - index < len(polynomial)
        )
        source = Fraction(polynomial[degree]) if degree < len(polynomial) else Fraction(0)
        result[degree] = source - correction / Fraction(degree)
    return result


@lru_cache(maxsize=None)
def _open_box_even_polynomial(shape: tuple[int, int, int], order: int) -> tuple[int, ...]:
    oriented = tuple(sorted(shape))
    n_sites = math.prod(oriented)
    n_bonds = sum(
        (oriented[axis] - 1) * math.prod(oriented[other] for other in range(3) if other != axis)
        for axis in range(3)
    )
    broken = box_broken_bond_poly(oriented)
    return _broken_to_even_subgraph(broken, n_sites, n_bonds, order)


@lru_cache(maxsize=None)
def exact_ht_bulk_v_series(max_order: int = 12) -> tuple[Fraction, ...]:
    """Derive exact infinite simple-cubic coefficients of log P(v) through ``max_order``.

    The rectangular finite-lattice weights are obtained recursively from exact open-box transfer
    polynomials.  A connected Eulerian graph whose bounding box is (a,b,c) has at least
    ``2*((a-1)+(b-1)+(c-1))`` edges.  Thus boxes with ``a+b+c <= max_order/2 + 3`` are complete for
    the requested truncation.  At order 12 the largest volume used is 27 sites.
    """

    if max_order < 0 or max_order > 12:
        raise ValueError("this experiment certifies 0 <= max_order <= 12")
    side_limit = max_order // 2 + 1
    shapes = sorted(
        (
            (a, b, c)
            for a in range(1, side_limit + 1)
            for b in range(1, side_limit + 1)
            for c in range(1, side_limit + 1)
            if 2 * ((a - 1) + (b - 1) + (c - 1)) <= max_order
        ),
        key=lambda shape: (sum(shape), shape),
    )
    weights: dict[tuple[int, int, int], list[Fraction]] = {}
    for shape in shapes:
        weight = _formal_log(_open_box_even_polynomial(shape, max_order), max_order)
        for subshape, subweight in weights.items():
            if all(subshape[axis] <= shape[axis] for axis in range(3)):
                embeddings = math.prod(
                    shape[axis] - subshape[axis] + 1 for axis in range(3)
                )
                weight = [
                    value - embeddings * old
                    for value, old in zip(weight, subweight, strict=True)
                ]
        weights[shape] = weight
    return tuple(
        sum((weight[degree] for weight in weights.values()), Fraction(0))
        for degree in range(max_order + 1)
    )


def _k_coefficients_from_v_series(
    v_coefficients: Sequence[Fraction], max_order: int
) -> tuple[Fraction, ...]:
    symbol = sp.symbols("K")
    expression = 3 * sp.log(sp.cosh(symbol))
    expression += sum(
        sp.Rational(value.numerator, value.denominator) * sp.tanh(symbol) ** degree
        for degree, value in enumerate(v_coefficients)
        if value
    )
    expanded = sp.series(expression, symbol, 0, max_order + 1).removeO().expand()
    result = []
    for degree in range(max_order + 1):
        coefficient = sp.Rational(expanded.coeff(symbol, degree))
        result.append(Fraction(int(coefficient.p), int(coefficient.q)))
    return tuple(result)


def _torus_reference_k_coefficients(max_order: int = 5) -> tuple[Fraction, ...]:
    L = 3
    if max_order >= 2 * L:
        raise ValueError("3^3 zero-winding torus reference is exact only below order 2L=6")
    p0 = winding_projected_torus_polynomial(L)
    log_p0 = _formal_log(p0, max_order)
    v_bulk = tuple(coefficient / (L**3) for coefficient in log_p0)
    return _k_coefficients_from_v_series(v_bulk, max_order)


# ---------------------------------------------------------------------------
# General numerical Taylor falsifier
# ---------------------------------------------------------------------------


def _interpolated_taylor(
    phi_callable: Callable[[mp.mpf], mp.mpf], degree: int, step: mp.mpf
) -> tuple[mp.mpf, ...]:
    nodes = [mp.mpf(index) * step for index in range(degree + 1)]
    values = []
    imaginary_tolerance = mp.power(10, -(mp.mp.dps // 2))
    for node in nodes:
        value = mp.mpc(phi_callable(node))
        if abs(value.imag) > imaginary_tolerance:
            raise ValueError(f"candidate free energy is complex at K={node}: {value}")
        values.append(value.real)
    vandermonde = mp.matrix(
        [[node**power for power in range(degree + 1)] for node in nodes]
    )
    coefficients = mp.lu_solve(vandermonde, mp.matrix(values))
    return tuple(mp.mpf(coefficients[index]) for index in range(degree + 1))


def falsify_free_energy(
    phi_callable: Callable[[mp.mpf], mp.mpf], name: str
) -> dict[str, object]:
    """Numerically Taylor-expand a candidate and find its first exact HT-series failure.

    The reference through K**5 comes from the exact, winding-projected 3x3x3 torus.  Two
    high-precision forward-interpolation grids provide a numerical stability estimate; the
    callable must represent a real function analytic from the right at K=0.
    """

    comparison_order = 5
    fit_degree = 8
    coarse = _interpolated_taylor(phi_callable, fit_degree, mp.mpf("0.0005"))
    fine = _interpolated_taylor(phi_callable, fit_degree, mp.mpf("0.00025"))
    reference = _torus_reference_k_coefficients(comparison_order)

    candidate_strings: dict[str, str] = {}
    reference_strings: dict[str, str] = {"0": "log(2)"}
    uncertainty_strings: dict[str, str] = {}
    failures: list[dict[str, str | int]] = []
    for degree in range(comparison_order + 1):
        candidate = fine[degree]
        uncertainty = abs(fine[degree] - coarse[degree])
        exact = mp.log(2) if degree == 0 else mp.mpf(reference[degree].numerator) / reference[degree].denominator
        difference = abs(candidate - exact)
        tolerance = max(mp.mpf("1e-25"), 20 * uncertainty)
        candidate_strings[str(degree)] = _mpstr(candidate)
        uncertainty_strings[str(degree)] = _mpstr(uncertainty)
        if degree:
            reference_strings[str(degree)] = _fraction_str(reference[degree])
        if difference > tolerance:
            failures.append(
                {
                    "order": degree,
                    "candidate": _mpstr(candidate),
                    "reference": _mpstr(exact),
                    "absolute_difference": _mpstr(difference),
                    "comparison_tolerance": _mpstr(tolerance),
                }
            )

    return {
        "name": name,
        "method": "degree-8 mpmath interpolation at steps 0.0005 and 0.00025",
        "precision": mp.mp.dps,
        "reference_lattice": "3x3x3 torus, eight-twist zero-winding projection",
        "reference_exact_through_order": 5,
        "candidate_coefficients": candidate_strings,
        "reference_coefficients": reference_strings,
        "coarse_fine_differences": uncertainty_strings,
        "first_failure_order": failures[0]["order"] if failures else None,
        "failures": failures,
    }


# ---------------------------------------------------------------------------
# Degang Zhang free-energy formula and numerical comparison
# ---------------------------------------------------------------------------


def claimed_degang_xi_argument(
    omega: mp.mpf, H: mp.mpf, H1: mp.mpf, H2: mp.mpf, m: int
) -> mp.mpf:
    """The published D_m(omega), evaluated literally."""

    if m < 1:
        raise ValueError("m must be positive")
    H, H1, H2 = mp.mpf(H), mp.mpf(H1), mp.mpf(H2)
    if H <= 0:
        raise ValueError("the literal dual-coupling formula requires H > 0")
    dual = mp.atanh(mp.exp(-2 * H))
    return (
        mp.cosh(2 * H1) * mp.cosh(2 * H2) * mp.cosh(2 * dual)
        - mp.sinh(2 * H1)
        * mp.cosh(2 * H2)
        * mp.sinh(2 * dual)
        * mp.cos(omega)
        - mp.cosh(2 * H1)
        * mp.sinh(2 * H2)
        * mp.sinh(2 * dual)
        * mp.cos(m * omega)
        + mp.sinh(2 * H1)
        * mp.sinh(2 * H2)
        * mp.cosh(2 * dual)
        * mp.cos((m - 1) * omega)
    )


def claimed_degang_free_energy_m(
    H: mp.mpf, H1: mp.mpf, H2: mp.mpf, m: int
) -> mp.mpf:
    """Numerically integrate the claimed finite-m free-energy expression."""

    H, H1, H2 = mp.mpf(H), mp.mpf(H1), mp.mpf(H2)
    breakpoints = [mp.pi * index / (4 * m) for index in range(4 * m + 1)]
    integral = mp.quad(
        lambda omega: mp.acosh(claimed_degang_xi_argument(omega, H, H1, H2, m)),
        breakpoints,
    )
    value = mp.log(2 * mp.sinh(2 * H)) / 2 + integral / (2 * mp.pi)
    if abs(mp.im(value)) > mp.power(10, -(mp.mp.dps // 2)):
        raise ValueError(f"claimed free energy became complex: {value}")
    return mp.re(value)


def _projected_torus_phi(K: mp.mpf) -> mp.mpf:
    p0 = winding_projected_torus_polynomial(3)
    value = mp.fsum(mp.mpf(coefficient) * mp.tanh(K) ** degree for degree, coefficient in enumerate(p0))
    return mp.log(2) + 3 * mp.log(mp.cosh(K)) + mp.log(value) / 27


def _bulk_taylor_value(K: mp.mpf, order: int = 12) -> tuple[mp.mpf, tuple[Fraction, ...]]:
    coefficients = _k_coefficients_from_v_series(exact_ht_bulk_v_series(order), order)
    value = mp.log(2) + mp.fsum(
        mp.mpf(coefficient.numerator) / coefficient.denominator * K**degree
        for degree, coefficient in enumerate(coefficients)
        if degree
    )
    return value, coefficients


def _cumulant_tail_bound(K: mp.mpf, known_order: int = 12) -> mp.mpf:
    """Uniform high-temperature remainder majorant after an exact K series.

    At K=0 the bond variables have absolute value one and dependency-graph degree at most ten.
    The spanning-tree cumulant bound gives

        |[K**r] phi| <= 3*2**(r-1)*r**(r-2)*11**(r-1)/r! < 64**r.

    Summing the geometric majorant produces the returned rigorous bound.
    """

    ratio = 64 * abs(mp.mpf(K))
    if ratio >= 1:
        raise ValueError("the 64|K| cumulant majorant requires |K| < 1/64")
    return ratio ** (known_order + 1) / (1 - ratio)


def free_energy_comparison() -> dict[str, object]:
    """Evaluate the claimed m trend and compare it with exact finite-lattice machinery."""

    K = mp.mpf("0.001")
    m_values = (2, 4, 8, 16, 32, 64)
    trend = {m: claimed_degang_free_energy_m(K, K, K, m) for m in m_values}
    claimed_limit = trend[64]
    finite_phi = _projected_torus_phi(K)
    bulk_taylor, k_coefficients = _bulk_taylor_value(K, 12)
    thermodynamic_tail_bound = _cumulant_tail_bound(K, 12)
    finite_size_bound = abs(finite_phi - bulk_taylor) + thermodynamic_tail_bound
    discrepancy = abs(claimed_limit - finite_phi)
    trend_change = abs(trend[64] - trend[32])

    return {
        "K": _mpstr(K),
        "m_values": list(m_values),
        "phi_by_m": {str(m): _mpstr(value) for m, value in trend.items()},
        "m32_to_m64_change": _mpstr(trend_change),
        "m_to_infinity_trend_estimate": _mpstr(claimed_limit),
        "finite_lattice": "3x3x3 torus, exact eight-twist zero-winding projection",
        "finite_lattice_phi": _mpstr(finite_phi),
        "bulk_K_series_through_12": {
            str(degree): _fraction_str(coefficient)
            for degree, coefficient in enumerate(k_coefficients)
            if degree and coefficient
        },
        "bulk_taylor_value_through_12": _mpstr(bulk_taylor),
        "thermodynamic_tail_bound": _mpstr(thermodynamic_tail_bound),
        "finite_size_error_bound": _mpstr(finite_size_bound),
        "absolute_discrepancy_to_finite_lattice": _mpstr(discrepancy),
        "discrepancy_minus_finite_size_bound": _mpstr(
            max(mp.mpf(0), discrepancy - finite_size_bound)
        ),
        "discrepancy_over_finite_size_bound": _mpstr(discrepancy / finite_size_bound),
    }


def _degang_taylor_callable(K: mp.mpf) -> mp.mpf:
    if K == 0:
        return mp.log(2)
    # m=8 is already in the stable m trend, and modes depending on m enter above the tested K**5.
    return claimed_degang_free_energy_m(K, K, K, 8)


# ---------------------------------------------------------------------------
# End-to-end artifact generation
# ---------------------------------------------------------------------------


def run_analysis(write_result: bool = True) -> dict[str, object]:
    mp.mp.dps = PRECISION
    zhang = zhang_critical_data()
    magnetization = low_temperature_magnetization_series(6)
    degang = degang_critical_data()
    free_energy = free_energy_comparison()
    falsifier_report = falsify_free_energy(_degang_taylor_callable, "Degang-Zhang-2021")
    torus_p0 = winding_projected_torus_polynomial(3)
    bulk_v = exact_ht_bulk_v_series(12)

    checks = [
        {
            "name": "Zhang critical root equals claimed golden-ratio form",
            "passed": mp.mpf(zhang["root_closed_form_error"]) < mp.mpf("1e-70"),
            "detail": f"absolute difference {zhang['root_closed_form_error']}",
        },
        {
            "name": "Zhang critical point disagrees with benchmark error bar",
            "passed": mp.mpf(zhang["standard_deviations"]) > mp.mpf("1e6"),
            "detail": f"{zhang['standard_deviations']} benchmark standard deviations",
        },
        {
            "name": "exact low-temperature magnetisation coefficient through u^6",
            "passed": magnetization
            == [Fraction(1), Fraction(0), Fraction(0), Fraction(-2), Fraction(0), Fraction(-12), Fraction(14)],
            "detail": "coefficients " + ", ".join(_fraction_str(value) for value in magnetization),
        },
        {
            "name": "claimed anisotropic critical line violates cubic rotational invariance",
            "passed": not bool(degang["rotational_invariance_passed"]),
            "detail": f"cyclic beta_c range {degang['beta_c_range']}; inequality falsifies the claim",
        },
        {
            "name": "twisted torus projects away the order-L winding loop",
            "passed": torus_p0[:6] == (1, 0, 0, 0, 81, 0),
            "detail": f"P0(v) coefficients through v^5 are {list(torus_p0[:6])}",
        },
        {
            "name": "torus and open-box finite-lattice HT coefficient agree",
            "passed": Fraction(torus_p0[4], 27) == bulk_v[4] == Fraction(3),
            "detail": "[v^4] log(P)/N = 3",
        },
        {
            "name": "claimed free energy differs beyond finite-size bound",
            "passed": mp.mpf(free_energy["absolute_discrepancy_to_finite_lattice"])
            > mp.mpf(free_energy["finite_size_error_bound"]),
            "detail": (
                f"discrepancy {free_energy['absolute_discrepancy_to_finite_lattice']}, "
                f"bound {free_energy['finite_size_error_bound']}"
            ),
        },
        {
            "name": "general falsifier locates first Degang-series failure",
            "passed": falsifier_report["first_failure_order"] == 4,
            "detail": f"first failure order {falsifier_report['first_failure_order']}",
        },
    ]

    result: dict[str, object] = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": PRECISION,
        },
        "data": {
            "zhang_2007": {
                "critical_condition": zhang,
                "magnetization": {
                    "variable": "u = exp(-4K)",
                    "maximum_order": 6,
                    "coefficients_ascending": [_fraction_str(value) for value in magnetization],
                    "series": "1 - 2*u^3 - 12*u^5 + 14*u^6 + O(u^7)",
                    "published_comparison": {
                        "Perk_u6": "14",
                        "Zhang_u6": "-18",
                        "computed_u6": _fraction_str(magnetization[6]),
                        "supports": "Perk coefficient; contradicts Zhang coefficient",
                    },
                    "exact_method": (
                        "integer transfer matrix with broken-bond and n_down axes; mixed backward "
                        "difference at side 3 using all eight plus-boundary boxes with sides 2 or 3"
                    ),
                },
            },
            "degang_2021": {
                "critical_condition": degang,
                "free_energy": free_energy,
            },
            "general_falsifier": {
                "reference": {
                    "torus_size": [3, 3, 3],
                    "projection": "average over 8 seam twists, selecting zero Z2 winding",
                    "raw_torus_leading_wrap_order": 3,
                    "projected_leading_wrap_order": 6,
                    "exact_K_orders": "0 through 5",
                    "projected_P_coefficients_through_6": list(torus_p0[:7]),
                    "open_box_bulk_logP_v_coefficients_through_12": {
                        str(degree): _fraction_str(coefficient)
                        for degree, coefficient in enumerate(bulk_v)
                        if coefficient
                    },
                },
                "degang_2021": falsifier_report,
            },
        },
        "checks": checks,
    }

    if write_result:
        RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
        RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    result = run_analysis(write_result=True)
    for check in result["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"{status}: {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in result["checks"])
    print(("PASS" if passed else "FAIL") + f": wrote {RESULT_PATH.relative_to(ROOT)}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
