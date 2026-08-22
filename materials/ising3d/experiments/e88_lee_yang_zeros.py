"""Exact Lee--Yang zeros beyond the frozen cubes: 5x5x4 / 5x4x4 boxes,
4x4 and 4x5 slab ladders, the L=5 isotropic-cube resource wall, and a
four-size test of the wave-8 adversarial correction models.

Run from the repository root with

    .venv/bin/python experiments/e88_lee_yang_zeros.py

Exact roles (see proofs/lee_yang_zeros.md):

* integer modular layer transfer + CRT decides every field polynomial;
* exact Sturm isolation on the circle-restricted polynomial decides every
  zero count and every isolating interval;
* exact rational alternating cosine series decide every angle enclosure;
* exact Fraction interval arithmetic decides every separation inequality;
* numpy complex128 is used ONLY inside the fixed-field diagnostic engine
  and is validated against the exact Sturm angles before diagnostic use.

No 3D benchmark value is used anywhere.
"""

from __future__ import annotations

import gc
import hashlib
import json
import math
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from ising.lattices import hyperrect
from ising.lee_yang import (
    _CRT_PRIMES,
    field_dos_from_enumeration,
    rational_field_polynomial,
)
from ising.transfer_matrix import box_broken_bond_poly, layer_bonds

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e88_lee_yang_zeros.py"
FROZEN_SCALING = ROOT / "results" / "lee_yang" / "scaling.json"
FROZEN_PRG = ROOT / "results" / "lee_yang" / "prg.json"
OUTPUT = ROOT / "results" / "lee_yang" / "zeros_l5.json"

DPS = 120
SERIAL_DIGITS = 60
THETA_RADIUS = Fraction(1, 10**60)
COSINE_EVEN_ORDER = 80
BARRIER_EDGE = Fraction(1, 50)
BARRIER_POWERS = (2, 4)
COUPLINGS = (("x_2_over_3", 2, 3), ("x_7_over_10", 7, 10))
LADDER_LENGTHS = (2, 3, 4, 5)
DIAG_THETA_RADIUS = Fraction(1, 10**9)
FF_VALIDATION_TOL = 1e-8
MAX_STATES = 1 << 21
_T = sp.Symbol("t")

mp.mp.dps = DPS


def fraction_text(value: Fraction) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def fraction_decimal(value: Fraction, digits: int = 30) -> str:
    value = Fraction(value)
    with mp.workdps(digits + 20):
        return mp.nstr(mp.mpf(value.numerator) / value.denominator, n=digits, strip_zeros=False)


def mp_text(value: object, digits: int = SERIAL_DIGITS) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()




# ---------------------------------------------------------------------------
# Exact modular layer transfer with a field grading.
#
# Same algorithm as ising.lee_yang._modular_open_polynomial (validated there
# against brute force), generalized beyond the stored 16-site guard and with
# the grade axis truncated at the palindromic half.
# ---------------------------------------------------------------------------


def layer_statistics(ns: int, bonds) -> tuple[np.ndarray, np.ndarray]:
    states = np.arange(1 << ns, dtype=np.uint64)
    down = np.zeros(1 << ns, dtype=np.int64)
    for bit in range(ns):
        down += ((states >> np.uint64(bit)) & np.uint64(1)).astype(np.int64)
    broken = np.zeros(1 << ns, dtype=np.int64)
    for i, j in bonds:
        bi = ((states >> np.uint64(i)) & np.uint64(1)).astype(bool)
        bj = ((states >> np.uint64(j)) & np.uint64(1)).astype(bool)
        broken += (bi != bj).astype(np.int64)
    return broken, ns - down


def modular_field_residues(
    cross: tuple[int, int],
    length: int,
    numerator: int,
    denominator: int,
    modulus: int,
) -> tuple[list[int], int]:
    """One modular pass; residues of c_k for k <= n_sites//2 (palindromic half)."""

    ns = cross[0] * cross[1]
    n_states = 1 << ns
    n_sites = ns * length
    grade_cap = n_sites // 2
    in_layer = layer_bonds(cross, (False, False))
    broken, down = layer_statistics(ns, in_layer)
    weights = np.array(
        [
            pow(numerator, int(q), modulus)
            * pow(denominator, len(in_layer) - int(q), modulus)
            % modulus
            for q in broken
        ],
        dtype=np.uint64,
    )
    mod = np.uint64(modulus)
    disagree = np.uint64(numerator % modulus)
    agree = np.uint64(denominator % modulus)

    width = min(ns + 1, grade_cap + 1)
    values = np.zeros((n_states, width), dtype=np.uint64)
    keep = down < width
    values[np.arange(n_states)[keep], down[keep]] = weights[keep]
    for _ in range(1, length):
        transformed = values
        for bit in range(ns):
            hi = 1 << (ns - bit - 1)
            lo = 1 << bit
            reshaped = transformed.reshape(hi, 2, lo, width)
            updated = np.empty_like(reshaped)
            lower = reshaped[:, 0]
            upper = reshaped[:, 1]
            updated[:, 0] = (agree * lower + disagree * upper) % mod
            updated[:, 1] = (disagree * lower + agree * upper) % mod
            transformed = updated.reshape(n_states, width)
        new_width = min(width + ns, grade_cap + 1)
        shifted = np.zeros((n_states, new_width), dtype=np.uint64)
        for down_count in range(min(ns + 1, new_width)):
            selected = down == down_count
            if np.any(selected):
                available = min(width, new_width - down_count)
                shifted[selected, down_count : down_count + available] = (
                    transformed[selected, :available] * weights[selected, None]
                ) % mod
        values = shifted
        width = new_width
    if n_states > 1 << 16:
        totals = [0] * width
        chunk = 1 << 16
        for start in range(0, n_states, chunk):
            block = values[start : start + chunk].sum(axis=0, dtype=np.uint64) % mod
            for k, value in enumerate(block):
                totals[k] = (totals[k] + int(value)) % mod
    else:
        totals = [int(value) for value in values.sum(axis=0, dtype=np.uint64) % mod]
    return totals, n_sites


def open_box_bond_count(a: int, b: int, c: int) -> int:
    return (a - 1) * b * c + a * (b - 1) * c + a * b * (c - 1)


def exact_field_polynomial(
    cross: tuple[int, int], length: int, numerator: int, denominator: int
) -> dict:
    """CRT-assembled exact palindromic integer field polynomial A_k in z."""

    ns = cross[0] * cross[1]
    if (1 << ns) > MAX_STATES:
        raise MemoryError(f"cross-section {cross} has {1 << ns} states")
    n_sites = ns * length
    n_bonds = open_box_bond_count(cross[0], cross[1], length)
    bound = (1 << n_sites) * max(numerator, denominator) ** n_bonds
    half = [0] * (n_sites // 2 + 1)
    combined = 1
    primes_used = 0
    for modulus in _CRT_PRIMES:
        residues, _ = modular_field_residues(cross, length, numerator, denominator, modulus)
        inverse = pow(combined % modulus, -1, modulus)
        for k, residue in enumerate(residues):
            correction = (
                (int(residue) - (half[k] % modulus)) % modulus
            ) * inverse % modulus
            half[k] += combined * correction
        combined *= modulus
        primes_used += 1
        if combined > bound:
            break
    else:
        raise ArithmeticError("CRT prime list exhausted before the exact bound")
    for value in half:
        if not 0 <= value < bound:
            raise ArithmeticError("CRT coefficient outside the exact bound")
    if n_sites % 2 == 0:
        coefficients = half + half[-2::-1]
    else:
        coefficients = half + half[::-1]
    if coefficients != coefficients[::-1]:
        raise ArithmeticError("spin-reversal palindrome failed after CRT")
    if any(value < 0 for value in coefficients):
        raise ArithmeticError("negative field-polynomial coefficient")
    return {
        "A_k": coefficients,
        "n_sites": n_sites,
        "n_bonds": n_bonds,
        "primes_used": primes_used,
        "coefficient_bound_bits": bound.bit_length(),
    }


def modular_zero_field_total(
    cross: tuple[int, int],
    length: int,
    numerator: int,
    denominator: int,
    modulus: int,
) -> int:
    """Independent scalar zero-field transfer for the z=1 sum check."""

    ns = cross[0] * cross[1]
    n_states = 1 << ns
    in_layer = layer_bonds(cross, (False, False))
    broken, _ = layer_statistics(ns, in_layer)
    weights = np.array(
        [
            pow(numerator, int(q), modulus)
            * pow(denominator, len(in_layer) - int(q), modulus)
            % modulus
            for q in broken
        ],
        dtype=np.uint64,
    )
    mod = np.uint64(modulus)
    disagree = np.uint64(numerator % modulus)
    agree = np.uint64(denominator % modulus)
    values = weights.copy()
    for _ in range(1, length):
        transformed = values
        for bit in range(ns):
            hi = 1 << (ns - bit - 1)
            lo = 1 << bit
            reshaped = transformed.reshape(hi, 2, lo)
            updated = np.empty_like(reshaped)
            lower = reshaped[:, 0]
            upper = reshaped[:, 1]
            updated[:, 0] = (agree * lower + disagree * upper) % mod
            updated[:, 1] = (disagree * lower + agree * upper) % mod
            transformed = updated.reshape(n_states)
        values = (transformed * weights) % mod
    return int(values.sum(dtype=np.uint64) % mod)


def exact_zero_field_total(
    cross: tuple[int, int], length: int, numerator: int, denominator: int
) -> tuple[int, int]:
    """CRT reconstruction of Z(z=1), independent of the field grading."""

    n_sites = cross[0] * cross[1] * length
    n_bonds = open_box_bond_count(cross[0], cross[1], length)
    bound = (1 << n_sites) * max(numerator, denominator) ** n_bonds
    value = 0
    combined = 1
    for prime_index, modulus in enumerate(_CRT_PRIMES, start=1):
        residue = modular_zero_field_total(
            cross, length, numerator, denominator, modulus
        )
        inverse = pow(combined % modulus, -1, modulus)
        correction = ((residue - value % modulus) % modulus) * inverse % modulus
        value += combined * correction
        combined *= modulus
        if combined > bound:
            if not 0 <= value < bound:
                raise ArithmeticError("zero-field CRT value outside its exact bound")
            return value, prime_index
    raise ArithmeticError("CRT prime list exhausted in scalar zero-field check")


def sum_identity_check(
    record: dict,
    cross: tuple[int, int],
    length: int,
    numerator: int,
    denominator: int,
) -> dict:
    """Independent z=1 scalar transfer, plus a small-box broken-poly control."""

    scalar_total, scalar_primes = exact_zero_field_total(
        cross, length, numerator, denominator
    )
    field_sum = sum(record["A_k"])
    scalar_ok = field_sum == scalar_total
    broken_poly_available = record["n_sites"] <= 48 and cross[0] * cross[1] <= 16
    broken_poly_ok = None
    if broken_poly_available:
        broken_poly = box_broken_bond_poly(
            tuple(cross) + (length,), periodic=(False, False, False)
        )
        poly_total = sum(
            int(count)
            * numerator**q
            * denominator ** (record["n_bonds"] - q)
            for q, count in enumerate(broken_poly)
        )
        broken_poly_ok = (
            sum(broken_poly) == 1 << record["n_sites"] and poly_total == scalar_total
        )
    return {
        "passed": scalar_ok and (broken_poly_ok is not False),
        "field_polynomial_sum": str(field_sum),
        "independent_scalar_zero_field_total": str(scalar_total),
        "scalar_zero_field_crt_primes": scalar_primes,
        "scalar_transfer_matches_field_sum": scalar_ok,
        "built_in_broken_bond_poly_available": broken_poly_available,
        "built_in_broken_bond_poly_matches_scalar_transfer": broken_poly_ok,
    }


# ---------------------------------------------------------------------------
# Exact circle restriction and Sturm isolation (pattern of e73, kept exact).
# ---------------------------------------------------------------------------


def divide_by_z_plus_one(coefficients: list[int]) -> list[int]:
    degree = len(coefficients) - 1
    if degree % 2 != 1 or coefficients != coefficients[::-1]:
        raise ValueError("division requires an odd-degree palindrome")
    quotient = [0] * degree
    quotient[0] = coefficients[0]
    for index in range(1, degree):
        quotient[index] = coefficients[index] - quotient[index - 1]
    reconstructed = (
        [quotient[0]]
        + [quotient[index - 1] + quotient[index] for index in range(1, degree)]
        + [quotient[-1]]
    )
    if reconstructed != coefficients or quotient != quotient[::-1]:
        raise ArithmeticError("exact z+1 division failed")
    return quotient


def reduced_power_coefficients(coefficients: list[int]) -> list[int]:
    """Return Q(t), ascending in t, for the palindromic field polynomial."""

    values = (
        divide_by_z_plus_one(coefficients)
        if (len(coefficients) - 1) % 2
        else list(coefficients)
    )
    half = (len(values) - 1) // 2
    power = [values[half]]
    if half == 0:
        return power
    previous = [1]
    current = [0, 1]
    for order in range(1, half + 1):
        if order == 1:
            chebyshev = current
        else:
            chebyshev = [0] * (len(current) + 1)
            for index, value in enumerate(current):
                chebyshev[index + 1] += 2 * value
            for index, value in enumerate(previous):
                chebyshev[index] -= value
            previous, current = current, chebyshev
        if len(power) < len(chebyshev):
            power.extend([0] * (len(chebyshev) - len(power)))
        multiplier = 2 * values[half - order]
        for index, value in enumerate(chebyshev):
            power[index] += multiplier * value
    return power


def poly_sign_at_dyadic(power: list[int], mantissa: int, level: int) -> int:
    """Exact sign of Q(mantissa / 2**level), pure integer arithmetic."""

    degree = len(power) - 1
    p_num = 1
    p_den = 1 << (level * degree)
    total = 0
    for index, coefficient in enumerate(power):
        total += coefficient * p_num * p_den
        if index < degree:
            p_num *= mantissa
            p_den >>= level
    if total == 0:
        return 0
    return 1 if total > 0 else -1


def cosine_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    """Exact alternating-series bounds (lower, upper) for cos, 0 <= v < 2."""

    if not 0 <= value < 2:
        raise ValueError("cosine bound is specialized to 0 <= value < 2")
    term = Fraction(1)
    partial = term
    even_partial = None
    for order in range(1, COSINE_EVEN_ORDER + 2):
        term *= -value * value
        term /= (2 * order - 1) * (2 * order)
        partial += term
        if order == COSINE_EVEN_ORDER:
            even_partial = partial
    if even_partial is None or COSINE_EVEN_ORDER % 2:
        raise AssertionError("the configured upper partial sum must end at even order")
    odd_partial = partial
    if not odd_partial <= even_partial:
        raise ArithmeticError("alternating cosine bounds were reversed")
    return odd_partial, even_partial


def certify_first_angle(power: list[int], t_lower: Fraction, t_upper: Fraction) -> dict:
    """Exact theta_1 enclosure over a Sturm-isolated t interval."""

    with mp.workdps(220):
        t_mid = (
            mp.mpf(t_lower.numerator) / t_lower.denominator
            + mp.mpf(t_upper.numerator) / t_upper.denominator
        ) / 2
        theta_center_mp = mp.acos(t_mid)
    theta_center = Fraction(mp.nstr(theta_center_mp, 140))
    if not 0 < theta_center < 2:
        raise ArithmeticError("first zero angle outside the certified cosine branch")
    theta_lower = theta_center - THETA_RADIUS
    theta_upper = theta_center + THETA_RADIUS
    cosine_at_lower, _ = cosine_bounds(theta_lower)
    _, cosine_at_upper = cosine_bounds(theta_upper)
    enclosure = cosine_at_upper < t_lower and t_upper < cosine_at_lower
    if not enclosure:
        raise ArithmeticError("cosine series enclosure failed")
    return {
        "theta_interval": [fraction_text(theta_lower), fraction_text(theta_upper)],
        "theta_radius": fraction_text(THETA_RADIUS),
        "theta_center_decimal": mp_text(theta_center_mp, 60),
        "exact_cosine_series_enclosure": enclosure,
        "cosine_series_orders": [COSINE_EVEN_ORDER + 1, COSINE_EVEN_ORDER],
    }


def exact_zero_certificate(coefficients: list[int]) -> dict:
    """Exact Sturm isolation of every t-root and certification of theta_1."""

    power = reduced_power_coefficients(coefficients)
    degree = len(coefficients) - 1
    expected = degree // 2 if degree % 2 == 0 else (degree - 1) // 2
    sym = sp.Poly.from_list(list(reversed(power)), gens=_T, domain=sp.ZZ)
    if power[-1] == 0:
        raise ArithmeticError("reduced polynomial lost degree")
    minus1, plus1 = sp.Integer(-1), sp.Integer(1)
    count = int(sym.count_roots(minus1, plus1))
    if count != expected:
        raise ArithmeticError(f"Sturm found {count} roots in (-1,1), expected {expected}")
    isolated = sym.intervals(eps=sp.Rational(1, 10**20))
    candidates = [
        (sp.Rational(bounds[0]), sp.Rational(bounds[1]), int(multiplicity))
        for bounds, multiplicity in isolated
        if bounds[0] < plus1 and bounds[1] > minus1
    ]
    if not candidates:
        raise ArithmeticError("no isolated interval overlaps (-1,1)")
    lo_rational, hi_rational, multiplicity = max(candidates, key=lambda item: item[1])
    level = 260
    scale = 1 << level
    lower_m = -((-int(lo_rational.p) * scale) // int(lo_rational.q))
    upper_m = (int(hi_rational.p) * scale) // int(hi_rational.q)
    lower_m = max(lower_m, -(scale - 1))
    upper_m = min(upper_m, scale - 1)
    lower_sign = poly_sign_at_dyadic(power, lower_m, level)
    upper_sign = poly_sign_at_dyadic(power, upper_m, level)
    if lower_sign == 0 or upper_sign == 0 or lower_sign == upper_sign:
        raise ArithmeticError("dyadic refinement bracket does not have a sign change")
    while upper_m - lower_m > 1:
        mid = (lower_m + upper_m) // 2
        sign = poly_sign_at_dyadic(power, mid, level)
        if sign == 0:
            raise ArithmeticError("first root is exactly dyadic; unexpected branch")
        if sign == lower_sign:
            lower_m = mid
        else:
            upper_m = mid
    t_lower = Fraction(lower_m, scale)
    t_upper = Fraction(upper_m, scale)
    root_count = multiplicity
    if root_count != 1:
        raise ArithmeticError(f"first-zero isolating interval has multiplicity {root_count}")
    angle = certify_first_angle(power, t_lower, t_upper)
    return {
        "t_interval": [fraction_text(t_lower), fraction_text(t_upper)],
        "t_interval_width": fraction_text(t_upper - t_lower),
        "exact_root_count_in_t_interval": root_count,
        "sturm_root_count_open_interval": count,
        "expected_count_by_palindromic_pairing": expected,
        "all_reduced_roots_simple": all(item[2] == 1 for item in candidates),
        "first_zero_simple": multiplicity == 1,
        **angle,
    }


# ---------------------------------------------------------------------------
# Fixed-field complex128 diagnostic engine (validated against exact Sturm).
# ---------------------------------------------------------------------------

_LAYER_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}
_FIXED_BASE_WEIGHT_CACHE: dict[tuple, np.ndarray] = {}


def layer_tables(cross: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    key = tuple(cross)
    if key not in _LAYER_CACHE:
        ns = cross[0] * cross[1]
        bonds = layer_bonds(cross, (False, False))
        _LAYER_CACHE[key] = layer_statistics(ns, bonds)
    return _LAYER_CACHE[key]


def fixed_field_value(cross, length, x_num, x_den, theta: float) -> float:
    """Sign of Re(z^{-N/2} P_z(z)) via renormalized complex128 transfer."""

    ns = cross[0] * cross[1]
    broken, down = layer_tables(cross)
    x = x_num / x_den
    weight_key = (tuple(cross), x_num, x_den)
    if weight_key not in _FIXED_BASE_WEIGHT_CACHE:
        _FIXED_BASE_WEIGHT_CACHE[weight_key] = x ** broken.astype(np.float64)
    layer_weight = _FIXED_BASE_WEIGHT_CACHE[weight_key] * np.exp(
        1j * theta * down.astype(np.float64)
    )
    vector = layer_weight.copy()
    for _ in range(1, length):
        for bit in range(ns):
            hi = 1 << (ns - bit - 1)
            lo = 1 << bit
            reshaped = vector.reshape(hi, 2, lo)
            lower = reshaped[:, 0]
            upper = reshaped[:, 1]
            updated = np.empty_like(reshaped)
            updated[:, 0] = lower + x * upper
            updated[:, 1] = x * lower + upper
            vector = updated.reshape(-1)
        vector = vector * layer_weight
        scale = float(np.abs(vector).max())
        if scale > 0 and math.isfinite(scale):
            vector = vector / scale
    total = complex(vector.sum())
    n_sites = ns * length
    return float((np.exp(-0.5j * theta * n_sites) * total).real)


def fixed_field_first_zero(cross, length, x_num, x_den, guess: float) -> float:
    """Bracket and bisect the first fixed-field zero (complex128 diagnostic)."""

    def g(theta):
        return fixed_field_value(cross, length, x_num, x_den, theta)

    scan_hi = min(guess * 2.0, 1.6)
    scan_lo = max(guess * 0.4, 0.02)
    grid = [scan_lo + (scan_hi - scan_lo) * k / 23.0 for k in range(24)]
    values = [g(theta) for theta in grid]
    bracket = None
    for k in range(len(grid) - 1):
        if values[k] == 0.0:
            return grid[k]
        if values[k] * values[k + 1] < 0:
            bracket = (grid[k], grid[k + 1])
            break
    if bracket is None:
        raise ArithmeticError("no fixed-field sign change found on the scan grid")
    lo, hi = bracket
    flo = g(lo)
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        fmid = g(mid)
        if fmid == 0.0:
            return mid
        if flo * fmid < 0:
            hi = mid
        else:
            lo, flo = mid, fmid
        if hi - lo < 1e-13:
            break
    root = 0.5 * (lo + hi)
    if g(root - 1e-9) * g(root + 1e-9) > 0:
        raise ArithmeticError("fixed-field root has no confirmed sign change")
    return root


# ---------------------------------------------------------------------------
# Exact rational interval arithmetic for the four-size separation analysis.
# ---------------------------------------------------------------------------


def iv_mul(a, b):
    products = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(products), max(products))




def iv_div_exact(a, exact: Fraction):
    if exact == 0:
        raise ZeroDivisionError("divisor interval degenerate at zero")
    return iv_mul(a, (Fraction(1) / exact, Fraction(1) / exact))


def divided_difference_weights(sizes: tuple[int, ...]) -> dict[int, Fraction]:
    nodes = {size: Fraction(1, size) for size in sizes}
    weights = {}
    for size in sizes:
        product = Fraction(1)
        for other in sizes:
            if other != size:
                product *= nodes[size] - nodes[other]
        weights[size] = 1 / product
    return weights


def four_size_model(
    angle_intervals: dict[int, tuple[Fraction, Fraction]],
    sizes: tuple[int, ...],
    power: int,
    inputs_certified: bool = True,
) -> dict:
    """Exact interval feasibility of theta(L)=e+L^-p C(1/L), deg C<=2, 4 sizes."""

    ordered = tuple(sorted(sizes))
    weights = divided_difference_weights(ordered)
    a_lo = Fraction(0)
    a_hi = Fraction(0)
    for size in ordered:
        term_lo = Fraction(size**power) * (
            angle_intervals[size][0] if weights[size] > 0 else angle_intervals[size][1]
        ) * weights[size]
        term_hi = Fraction(size**power) * (
            angle_intervals[size][1] if weights[size] > 0 else angle_intervals[size][0]
        ) * weights[size]
        a_lo += term_lo
        a_hi += term_hi
    b_exact = sum((Fraction(size**power) * weights[size] for size in ordered), Fraction(0))
    if b_exact == 0:
        raise ZeroDivisionError("DD3 coefficient B(p) vanished")
    edge_interval = iv_div_exact((a_lo, a_hi), b_exact)
    frozen_value = (
        a_lo - BARRIER_EDGE * b_exact,
        a_hi - BARRIER_EDGE * b_exact,
    )
    frozen_edge_fits = frozen_value[0] <= 0 <= frozen_value[1]
    theta_min_hi = min(angle_intervals[s][1] for s in ordered)
    edge_excluded = edge_interval[1] <= 0 or edge_interval[0] >= theta_min_hi

    result = {
        "model": "theta(L)=e+L^(-p)*(c0+c1/L+c2/L^2)",
        "p": power,
        "sizes": list(ordered),
        "inputs_certified": inputs_certified,
        "b_exact": fraction_text(b_exact),
        "edge_interval_e_star": [
            fraction_text(edge_interval[0]),
            fraction_text(edge_interval[1]),
        ],
        "edge_interval_e_star_decimal": [
            fraction_decimal(edge_interval[0]),
            fraction_decimal(edge_interval[1]),
        ],
        "frozen_edge_1_over_50_dd3_interval": [
            fraction_text(frozen_value[0]),
            fraction_text(frozen_value[1]),
        ],
        "frozen_edge_1_over_50_admits_quadratic": frozen_edge_fits,
        "positive_edge_below_data": not edge_excluded,
    }
    if edge_excluded:
        result["verdict"] = "EXCLUDED_NO_ADMISSIBLE_EDGE"
        result["exclusion_reason"] = (
            "edge e*=A/B is nonpositive or fails to lie below every input angle"
        )
        return result

    fit_sizes = ordered[:3]
    matrix = [[Fraction(1), Fraction(1, s), Fraction(1, s * s)] for s in fit_sizes]
    inverse = invert_fraction_matrix(matrix)
    c_intervals = []
    for row in range(3):
        lo = Fraction(0)
        hi = Fraction(0)
        for k, size in enumerate(fit_sizes):
            rhs_lo = Fraction(size**power) * (
                angle_intervals[size][0] - edge_interval[1]
            )
            rhs_hi = Fraction(size**power) * (
                angle_intervals[size][1] - edge_interval[0]
            )
            lo += min(inverse[row][k] * rhs_lo, inverse[row][k] * rhs_hi)
            hi += max(inverse[row][k] * rhs_lo, inverse[row][k] * rhs_hi)
        c_intervals.append((lo, hi))
    c0_lo = c_intervals[0][0]
    c1_abs = max(abs(c_intervals[1][0]), abs(c_intervals[1][1]))
    c2_lo = c_intervals[2][0]
    disc_c_upper = c1_abs * c1_abs - 4 * c0_lo * c2_lo
    disc_q_upper = (power + 1) ** 2 * c1_abs * c1_abs - 4 * power * (power + 2) * c0_lo * c2_lo
    correction_positive = c0_lo > 0 and c2_lo > 0 and disc_c_upper < 0
    monotone = disc_q_upper < 0
    result.update(
        {
            "c_intervals": {
                name: [fraction_text(iv[0]), fraction_text(iv[1])]
                for name, iv in zip(("c0", "c1", "c2"), c_intervals)
            },
            "correction_discriminant_upper": fraction_text(disc_c_upper),
            "derivative_discriminant_upper": fraction_text(disc_q_upper),
            "correction_everywhere_positive": correction_positive,
            "strictly_decreasing_for_all_positive_L": monotone,
            "proof": (
                "C(u)>0 needs c0>0, c2>0, disc(C)<0; Q(L)=p*c0*L^2+(p+1)*c1*L+"
                "(p+2)*c2>0 for all L>0 needs disc(Q)<0 with positive leading "
                "coefficient; then theta'(L)=-L^(-p-3)Q(L)<0."
            ),
        }
    )
    if correction_positive and monotone:
        result["verdict"] = "ADMISSIBLE"
        predictions = {}
        for size in (6, 8):
            lo = edge_interval[0] + (
                c_intervals[0][0] + c_intervals[1][0] / size + c_intervals[2][0] / (size * size)
            ) / size**power
            hi = edge_interval[1] + (
                c_intervals[0][1] + c_intervals[1][1] / size + c_intervals[2][1] / (size * size)
            ) / size**power
            predictions[str(size)] = [fraction_text(lo), fraction_text(hi)]
        result["predictions_at_larger_sizes"] = predictions
    else:
        result["verdict"] = "EXCLUDED_BY_CORRECTION_STRUCTURE"
        reasons = []
        if not (c0_lo > 0 and c2_lo > 0):
            reasons.append("c0 or c2 not positive")
        if disc_c_upper >= 0:
            reasons.append("disc(C) interval not negative")
        if disc_q_upper >= 0:
            reasons.append("disc(Q) interval not negative")
        result["exclusion_reason"] = "; ".join(reasons)
    return result


def invert_fraction_matrix(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    work = [
        row[:] + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for pivot_index in range(size):
        pivot_row = next(row for row in range(pivot_index, size) if work[row][pivot_index] != 0)
        work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
        pivot = work[pivot_index][pivot_index]
        work[pivot_index] = [value / pivot for value in work[pivot_index]]
        for row in range(size):
            if row == pivot_index:
                continue
            factor = work[row][pivot_index]
            if factor:
                work[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(work[row], work[pivot_index])
                ]
    return [row[size:] for row in work]


def wave8_model_prediction(model: dict, edge: Fraction, size: int) -> Fraction:
    """Theta_p(size) nominal value from a stored wave-8 adversarial model."""

    coeffs = model["nominal_coefficients_exact"]
    c0 = Fraction(coeffs["c0"])
    c1 = Fraction(coeffs["c1"])
    c2 = Fraction(coeffs["c2"])
    power = int(model["p"])
    correction = c0 + c1 / size + c2 / (size * size)
    return edge + Fraction(1, size**power) * correction


# ---------------------------------------------------------------------------
# Resource-wall measurement for the isotropic 5x5x5 cube.
# ---------------------------------------------------------------------------


def resource_wall(single_prime_seconds: float, single_prime_ops: float) -> dict:
    ns = 25
    n_bonds = open_box_bond_count(5, 5, 5)
    cap = (ns * 5) // 2 + 1
    full_ops = sum(4 * ns * (1 << ns) * (ns * j + 1) for j in range(1, 5))
    capped_ops = sum(
        4 * ns * (1 << ns) * min(ns * j + 1, cap) for j in range(1, 5)
    )
    rate = single_prime_ops / single_prime_seconds
    per_coupling = {}
    for label, numerator, denominator in (
        ("x_2_over_3", 2, 3),
        ("x_7_over_10", 7, 10),
    ):
        bound = (1 << 125) * max(numerator, denominator) ** n_bonds
        combined = 1
        stored_count = 0
        stored_sufficient = False
        for modulus in _CRT_PRIMES:
            combined *= modulus
            stored_count += 1
            if combined > bound:
                stored_sufficient = True
                break
        minimum_30bit_moduli = (bound.bit_length() + 29) // 30
        estimate_moduli = stored_count if stored_sufficient else minimum_30bit_moduli
        seconds = capped_ops / rate * estimate_moduli
        per_coupling[label] = {
            "coefficient_bound_bits": bound.bit_length(),
            "stored_crt_prime_table_length": len(_CRT_PRIMES),
            "stored_crt_prime_table_sufficient": stored_sufficient,
            "crt_primes_needed_from_stored_table_exact_bound": (
                stored_count if stored_sufficient else None
            ),
            "minimum_30bit_moduli_by_bit_count": minimum_30bit_moduli,
            "estimated_seconds_truncated_at_minimum_moduli": seconds,
            "estimated_hours_truncated_at_minimum_moduli": seconds / 3600.0,
        }
    bytes_full = (1 << ns) * 126 * 8
    bytes_capped = (1 << ns) * 63 * 8
    return {
        "claim_tag": "COMPUTATION",
        "statement": (
            "The exact isotropic 5x5x5 field polynomial requires a 2^25-state "
            "graded layer transfer; measured extrapolation below, not run."
        ),
        "measured_single_prime_seconds_5x4_cross": single_prime_seconds,
        "measured_single_prime_element_ops_5x4_cross": single_prime_ops,
        "measured_rate_element_ops_per_second": rate,
        "single_prime_element_ops_full_width": full_ops,
        "single_prime_element_ops_palindrome_truncated": capped_ops,
        "final_grade_buffer_bytes_full_width": bytes_full,
        "final_grade_buffer_bytes_truncated": bytes_capped,
        "peak_working_set_bytes_two_buffers_full_width": 2 * bytes_full,
        "peak_working_set_bytes_two_buffers_truncated": 2 * bytes_capped,
        "conservative_peak_working_set_bytes_current_three_buffer_path_full_width": 3 * bytes_full,
        "conservative_peak_working_set_bytes_current_three_buffer_path_truncated": 3 * bytes_capped,
        "per_coupling": per_coupling,
        "machine_context": {
            "total_ram_bytes": 103079215104,
            "note": "96 GiB unified memory shared with sibling agents; load>60",
        },
        "scope": (
            "Linear extrapolation from the measured single-prime rate of the "
            "identical code path at a 5x4 cross-section; no 2^25-state attempt "
            "was made."
        ),
    }


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def check(checks: list, name: str, passed: bool, detail: str) -> bool:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(("PASS" if passed else "FAIL"), name, "--", detail, flush=True)
    return bool(passed)


def sorted_shape3(cross: tuple[int, int], length: int) -> tuple[int, int, int]:
    dims = sorted([cross[0], cross[1], length], reverse=True)
    return (dims[0], dims[1], dims[2])


def main() -> None:
    started = time.time()
    checks: list[dict] = []
    frozen_scaling = json.loads(FROZEN_SCALING.read_text())
    frozen_prg = json.loads(FROZEN_PRG.read_text())
    cube_data = frozen_scaling["data"]["three_dimensional"]["exact_open_cube_zero_data"]
    barrier = frozen_prg["data"]["exact_identifiability_barrier"]["couplings"]

    data: dict = {
        "exact_zero_data": {},
        "four_size_separation": {},
        "isotropic_diagnostic": {},
        "resource_wall_5x5x5": {},
    }
    all_ok = True

    families: dict[str, list[tuple[tuple[int, int], int]]] = {
        "ladder16": [((4, 4), L) for L in LADDER_LENGTHS],
        "ladder20": [((4, 5), L) for L in LADDER_LENGTHS],
        "cube3": [((3, 3), 3)],
    }
    single_prime_seconds = None
    single_prime_ops = None
    for grid_id, numerator, denominator in COUPLINGS:
        coupling_record: dict = {}
        for family, entries in families.items():
            family_record: dict = {}
            for cross, length in entries:
                t0 = time.time()
                record = exact_field_polynomial(cross, length, numerator, denominator)
                wall = time.time() - t0
                ns = cross[0] * cross[1]
                shape3 = sorted_shape3(cross, length)
                sum_certificate = sum_identity_check(
                    record, cross, length, numerator, denominator
                )
                sum_ok = sum_certificate["passed"]
                zero_record = exact_zero_certificate(record["A_k"])
                theta_center = float(zero_record["theta_center_decimal"])
                ff_theta = fixed_field_first_zero(
                    cross, length, numerator, denominator, theta_center
                )
                ff_gap = abs(ff_theta - theta_center)
                record.update(
                    {
                        "shape_sorted": list(shape3),
                        "x": f"{numerator}/{denominator}",
                        "wall_seconds": round(wall, 3),
                        "sum_identity_exact_passed": sum_ok,
                        "sum_identity_certificate": sum_certificate,
                        "first_zero_certificate": zero_record,
                        "fixed_field_theta_1_diagnostic": repr(ff_theta),
                        "fixed_field_vs_exact_gap": ff_gap,
                    }
                )
                family_record[f"{cross[0]}x{cross[1]}x{length}"] = record
                label = f"{cross[0]}x{cross[1]}x{length}"
                all_ok &= check(
                    checks,
                    f"sum_identity:{grid_id}:{label}",
                    sum_ok,
                    "exact integer identity holds",
                )
                all_ok &= check(
                    checks,
                    f"sturm_count:{grid_id}:{label}",
                    zero_record["sturm_root_count_open_interval"] == zero_record["expected_count_by_palindromic_pairing"],
                    f"count={zero_record['sturm_root_count_open_interval']}",
                )
                all_ok &= check(
                    checks,
                    f"fixed_field_validation:{grid_id}:{label}",
                    ff_gap < FF_VALIDATION_TOL,
                    f"gap={ff_gap:.2e}",
                )
                if family == "ladder16" and length == 4:
                    frozen_ak = [int(v) for v in cube_data[grid_id]["sizes"]["4"]["A_k"]]
                    all_ok &= check(
                        checks,
                        f"frozen_cube4_Ak_equality:{grid_id}",
                        record["A_k"] == frozen_ak,
                        f"degree={len(record['A_k']) - 1}",
                    )
                if family == "cube3":
                    lat = hyperrect((3, 3, 3), periodic=False)
                    enum_dos = field_dos_from_enumeration(lat)
                    enum_poly = [
                        int(v)
                        for v in rational_field_polynomial(
                            enum_dos, numerator, denominator
                        )
                    ]
                    all_ok &= check(
                        checks,
                        f"bruteforce_3cube_equality:{grid_id}",
                        record["A_k"] == enum_poly,
                        f"degree={len(record['A_k']) - 1}",
                    )
                gc.collect()
            coupling_record[family] = family_record
        data["exact_zero_data"][grid_id] = coupling_record

    # measured single-prime rate on the 5x4 cross-section (wall extrapolation)
    t0 = time.time()
    residues, _ = modular_field_residues((4, 5), 5, 2, 3, _CRT_PRIMES[0])
    single_prime_seconds = time.time() - t0
    single_prime_ops = sum(
        4 * 20 * (1 << 20) * min(20 * j + 1, 51) for j in range(1, 5)
    )
    del residues
    print(f"measured single-prime 5x4x5 pass: {single_prime_seconds:.2f}s", flush=True)

    # --- four-size exact separation on the slab ladders --------------------
    for grid_id, _, _ in COUPLINGS:
        coupling = data["exact_zero_data"][grid_id]
        analysis = {}
        for family in ("ladder16", "ladder20"):
            intervals = {}
            for key, record in coupling[family].items():
                length = int(key.split("x")[-1])
                cert = record["first_zero_certificate"]
                intervals[length] = (
                    Fraction(cert["theta_interval"][0]),
                    Fraction(cert["theta_interval"][1]),
                )
            sizes = tuple(sorted(intervals))
            models = [
                four_size_model(intervals, sizes, power) for power in BARRIER_POWERS
            ]
            analysis[family] = {
                "claim_tag": "COMPUTATION",
                "scope": (
                    "fixed-cross-section slab ladder; the isotropic relation "
                    "sigma=d/p-1 does NOT apply to slabs"
                ),
                "sizes": list(sizes),
                "models": models,
            }
        data["four_size_separation"][grid_id] = analysis

    # --- frozen isotropic protocol: exact L=5 remains intentionally unresolved
    for grid_id, _, _ in COUPLINGS:
        edge = Fraction(barrier[grid_id]["common_edge_exact"])
        predictions = {}
        for model in barrier[grid_id]["models"]:
            value = wave8_model_prediction(model, edge, 5)
            predictions[f"p={model['p']}"] = {
                "sigma_exact_d3": model["sigma_exact_d3"],
                "theta_5_prediction_exact": fraction_text(value),
                "theta_5_prediction_decimal": fraction_decimal(value),
            }
        data["isotropic_diagnostic"][grid_id] = {
            "claim_tag": "UNRESOLVED",
            "exact_isotropic_cube_sizes_available": [2, 3, 4],
            "frozen_protocol_raw_sizes_required": list(range(4, 16)),
            "protocol_instantiable": False,
            "statement": (
                "No exact theta_1(5) for the 5x5x5 cube was produced: the "
                "2^25-state resource wall below prevents an exact Sturm input. "
                "The frozen protocol is not shortened, refit, or retuned."
            ),
            "adversarial_model_verdict": (
                "Neither wave-8 model is selected or excluded for isotropic "
                "cubes because the fourth exact isotropic datum is absent."
            ),
            "wave8_model_theta5_predictions_without_comparison": predictions,
        }

    data["resource_wall_5x5x5"] = resource_wall(single_prime_seconds, single_prime_ops)

    output = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "frozen_scaling_sha256": sha256_file(FROZEN_SCALING),
            "frozen_prg_sha256": sha256_file(FROZEN_PRG),
            "exact_arithmetic": "python int / Fraction; sympy Sturm over ZZ",
            "mpmath_dps": DPS,
            "benchmarks_used": "none",
        },
        "definitions": {
            "x": "e^{-2K} bond fugacity, exact rational",
            "z": "e^{-2H_f} field fugacity",
            "A_k": "coefficients of den^n_bonds * Z in z, exact integers",
            "theta_1": "smallest positive zero angle of the field polynomial",
            "first_zero": "largest t-root in (-1,1) of the circle-restricted polynomial",
        },
        "data": data,
        "checks": checks,
    }
    OUTPUT.write_text(json.dumps(output, indent=1))
    print(f"wrote {OUTPUT} in {time.time() - started:.1f}s", flush=True)
    if not all_ok:
        raise SystemExit(1)
    print("PASS")


if __name__ == "__main__":
    main()
