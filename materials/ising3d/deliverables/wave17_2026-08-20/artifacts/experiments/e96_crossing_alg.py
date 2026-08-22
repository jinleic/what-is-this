"""Exact algebraic isolation of the four open-2x3 forced-value crossings.

This is a deliberately narrow upgrade of ``e81_spectral_interval2.py``.  It
proves that every inherited fine c-jump bracket contains exactly one forced
crossing, and represents that parameter as the t-coordinate of a unique,
nonsingular real zero of an explicit five-polynomial rational system.

It does *not* claim minimal polynomials or absence of further crossings on all
of [1/4, 7/20]: direct expansion of the first nested resultant stage was
stopped after a recorded 1,200-second wall.

All certificate decisions below use integer/rational arithmetic.  mpmath only
constructs rational Krawczyk preconditioners; the final interval inequalities
are rechecked exactly.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from fractions import Fraction
from hashlib import sha256
import argparse
import json
import math
from pathlib import Path
import sys
import time
from typing import Iterable, Sequence

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ROWS = 2
COLS = 3
SITES = ROWS * COLS
DIMENSION = 1 << SITES
MAX_T_DEGREE = 26

T, X, Z = sp.symbols("t x z")
VARIABLES = (T, sp.Symbol("a"), sp.Symbol("b"), sp.Symbol("c"), sp.Symbol("d"))

RESULT_PATH = ROOT / "results" / "spectral" / "crossing_alg.json"

# Numerical guides are not certificate data.  They were obtained from
# 400-dps sector Newton solves and are immediately converted to 45-digit
# rational box centers.  Krawczyk containment proves a zero independently.
GUIDES: dict[int, dict[str, str]] = {
    1: {
        "t": "0.2679857913395631635058572656299889509794519005568117647624297692561946581860637316692274667547047856396420923729129714437344606425254110111442210310296709311350197500669038770370717217121122466149012453190893005354366211633441749044082315557559804335456582112834312479817226871079660589153335001835541922085834469282315149982766837382115017",
        "a": "0.00184453999526173177485714039898674088258688435928789288132233",
        "b": "0.00249235820342243318571867671890244971143066702247683297476844",
        "c": "0.00572659902631792227062979406290472358706131603810974683696523",
        "d": "0.00773782953886513914727683864666871617969368028516468134432645",
    },
    2: {
        "t": "0.2909723643827632854657665045581399142883113067606008178351736899447719393030686195469742541424131631325556838437678712509697466643357174921663263480123337753463886480548357935839386271500272142275662774789689017382268503134352798016151732420886105736862417266465528022058251624328125170507023699248476054715533357540722416674468812172932571",
        "a": "0.00204719215524092525479721896249304905784185339085624106901376",
        "b": "0.00321491818903063463924435413654296981419588038130932105858597",
        "c": "0.00650205453315528717830709290400318526579196122695586004263878",
        "d": "0.0102108506674352572473980080406751076713757567453375561784765",
    },
    3: {
        "t": "0.3007641433575038941475275141584482553865814856299889959022718875415688158697567991271617905509971021515870872503594237462870279485895323606598183018432601532545789196854444967020620772326525677076871320487987355796538785304555898969246945179581410945400904022867620497188881274167587345571485161139132373698380211437730182384154361803744981",
        "a": "0.00208748312315866028014806936542013437552289931722237503978343",
        "b": "0.00351034706027657480703424855107881027329366876631544512391903",
        "c": "0.00677313300402350414382556080707359274278068711353842911331349",
        "d": "0.0113898154508476198997823864317602191688206965764807037883604",
    },
    4: {
        "t": "0.3297417801771730596824095134864066744524266453320935273835479549145863167265126261782351190404734796808322575249696446449404761722867593199079686807491021361465227076267081737207152415125144399711008292683318132679130168570750487636954956300789821410935716709756628216296796091120184144948784209485935484641879835033603782416918091012803887",
        "a": "0.00206456699119169286407664323993004843437385499655769505480318",
        "b": "0.00427765600536194469802168323833139014613160297178216945455564",
        "c": "0.00733125251611156469914253928385498355176491516158857502894029",
        "d": "0.0151899049467354898628262095477086540991757747721731970959256",
    },
}


@dataclass(frozen=True)
class SectorBasis:
    """Integer joint-character basis with diagonal Gram form."""

    character: int
    vectors: tuple[dict[int, int], ...]
    gram: tuple[int, ...]


@dataclass(frozen=True)
class CrossingSpec:
    """One inherited fine bracket and its ordered sector data."""

    number: int
    lower: Fraction
    upper: Fraction
    d_sector: int
    d_index: int
    thresholds: tuple[Fraction, ...]
    threshold_global_counts: tuple[int, ...]
    threshold_sector_counts: tuple[tuple[int, ...], ...]
    contraction_upper: Fraction


CROSSINGS = (
    CrossingSpec(
        1,
        Fraction(2248027753, 8388608000),
        Fraction(4496055507, 16777216000),
        4,
        4,
        (Fraction(1, 500), Fraction(1, 250), Fraction(3, 500), Fraction(7, 1000), Fraction(1, 125)),
        (1, 2, 3, 4, 5),
        (
            (1, 0, 0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 1),
            (2, 0, 0, 0, 0, 1, 0, 1),
            (2, 0, 0, 0, 1, 1, 0, 1),
        ),
        Fraction(1, 1000),
    ),
    CrossingSpec(
        2,
        Fraction(4881706207, 16777216000),
        Fraction(152553319, 524288000),
        5,
        5,
        (Fraction(1, 400), Fraction(1, 250), Fraction(7, 1000), Fraction(19, 2000), Fraction(21, 2000)),
        (1, 2, 3, 5, 6),
        (
            (1, 0, 0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 1),
            (2, 0, 0, 0, 1, 1, 0, 1),
            (2, 0, 0, 0, 1, 2, 0, 1),
        ),
        Fraction(1, 1000),
    ),
    CrossingSpec(
        3,
        Fraction(2522992499, 8388608000),
        Fraction(5045984999, 16777216000),
        6,
        6,
        (Fraction(1, 400), Fraction(1, 250), Fraction(7, 1000), Fraction(11, 1000), Fraction(3, 250)),
        (1, 2, 3, 6, 7),
        (
            (1, 0, 0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 1),
            (2, 0, 0, 0, 1, 2, 0, 1),
            (2, 0, 0, 0, 1, 2, 1, 1),
        ),
        Fraction(1, 1000),
    ),
    CrossingSpec(
        4,
        Fraction(553214907, 1677721600),
        Fraction(5532149071, 16777216000),
        4,
        7,
        (Fraction(3, 1000), Fraction(1, 200), Fraction(1, 125), Fraction(7, 500), Fraction(2, 125)),
        (1, 2, 3, 7, 8),
        (
            (1, 0, 0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 0),
            (1, 0, 0, 0, 0, 1, 0, 1),
            (2, 0, 0, 0, 1, 2, 1, 1),
            (2, 0, 0, 0, 2, 2, 1, 1),
        ),
        Fraction(1, 500),
    ),
)

T_RADIUS = Fraction(5, 10**11)
EIGEN_RADIUS = Fraction(2, 10**10)
INITIAL_EIGEN_RADIUS = Fraction(1, 10**10)
DERIVATIVE_TUBE_RADIUS = Fraction(5, 10**9)
ORDER_MARGIN = Fraction(1, 10**7)
DERIVATIVE_BRANCH_UPPER = Fraction(1, 20)
REFINED_BRANCH_MOTION = INITIAL_EIGEN_RADIUS + DERIVATIVE_BRANCH_UPPER * T_RADIUS
COORDINATE_CONTAINMENT_UPPER = Fraction(1, 100)


def site(row: int, col: int) -> int:
    return row * COLS + col


def layer_bonds() -> list[tuple[int, int]]:
    return [
        *( (site(row, col), site(row, col + 1)) for row in range(ROWS) for col in range(COLS - 1) ),
        *( (site(row, col), site(row + 1, col)) for row in range(ROWS - 1) for col in range(COLS) ),
    ]


def polynomial_matrix() -> tuple[list[list[list[int]]], int, int]:
    """Return coefficient lists for M(t)=(2t)^7 q(t)^4 R(t)."""

    bonds = layer_bonds()
    energies: list[int] = []
    for state in range(DIMENSION):
        spins = [1 - 2 * ((state >> (SITES - 1 - bit)) & 1) for bit in range(SITES)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    assert all((energy - parity) % 2 == 0 for energy in energies)
    exponents = [(energy - parity) // 2 for energy in energies]
    shift = -min(exponents)
    degree = max(exponent + shift for exponent in exponents)
    maximum_degree = 2 * SITES + 2 * degree
    assert (degree, shift, maximum_degree) == (7, 4, MAX_T_DEGREE)

    coefficients = [
        [[0] * (maximum_degree + 1) for _ in range(DIMENSION)]
        for _ in range(DIMENSION)
    ]
    for row in range(DIMENSION):
        for column in range(row, DIMENSION):
            entry = [0] * (maximum_degree + 1)
            for middle, exponent in enumerate(exponents):
                q_power = exponent + shift
                hamming_power = (row ^ middle).bit_count() + (middle ^ column).bit_count()
                base_power = degree - q_power + hamming_power
                integer_factor = 2 ** (degree - q_power)
                for even_power in range(q_power + 1):
                    entry[base_power + 2 * even_power] += integer_factor * math.comb(q_power, even_power)
            coefficients[row][column] = entry
            coefficients[column][row] = entry
    assert all(value >= 0 for row in coefficients for entry in row for value in entry)
    assert max(index for row in coefficients for entry in row for index, value in enumerate(entry) if value) == MAX_T_DEGREE
    return coefficients, degree, shift


def state_permutation(site_permutation: Sequence[int]) -> list[int]:
    result: list[int] = []
    for state in range(DIMENSION):
        moved = 0
        for old_site, new_site in enumerate(site_permutation):
            bit = (state >> (SITES - 1 - old_site)) & 1
            moved |= bit << (SITES - 1 - new_site)
        result.append(moved)
    return result


def make_sector_bases() -> tuple[SectorBasis, ...]:
    row_reflection = state_permutation([site(ROWS - 1 - row, col) for row in range(ROWS) for col in range(COLS)])
    column_reflection = state_permutation([site(row, COLS - 1 - col) for row in range(ROWS) for col in range(COLS)])
    spin_flip = [state ^ (DIMENSION - 1) for state in range(DIMENSION)]
    generators = (row_reflection, column_reflection, spin_flip)

    group: list[list[int]] = []
    for mask in range(1 << len(generators)):
        permutation = list(range(DIMENSION))
        for index, generator in enumerate(generators):
            if mask & (1 << index):
                permutation = [generator[state] for state in permutation]
        group.append(permutation)

    seen: set[int] = set()
    representatives: list[int] = []
    for state in range(DIMENSION):
        if state not in seen:
            orbit = {permutation[state] for permutation in group}
            representatives.append(min(orbit))
            seen.update(orbit)

    output: list[SectorBasis] = []
    for character in range(8):
        vectors: list[dict[int, int]] = []
        gram: list[int] = []
        for representative in representatives:
            vector: dict[int, int] = {}
            for group_element, permutation in enumerate(group):
                sign = -1 if (character & group_element).bit_count() % 2 else 1
                image = permutation[representative]
                vector[image] = vector.get(image, 0) + sign
            vector = {state: coefficient for state, coefficient in vector.items() if coefficient}
            if vector:
                vectors.append(vector)
                gram.append(sum(coefficient * coefficient for coefficient in vector.values()))
        output.append(SectorBasis(character, tuple(vectors), tuple(gram)))
    assert [len(basis.gram) for basis in output] == [14, 6, 6, 6, 10, 10, 6, 6]
    return tuple(output)


def reciprocal_chebyshev() -> tuple[sp.Expr, ...]:
    """Return t^k+t^-k as exact polynomials in z=t+t^-1 for k<=13."""

    values: list[sp.Expr] = [sp.Integer(2), Z]
    for _ in range(2, 14):
        values.append(sp.expand(Z * values[-1] - values[-2]))
    return tuple(values)


def entry_to_z(entry: Sequence[int], chebyshev: Sequence[sp.Expr]) -> sp.Poly:
    """Convert t^-13 entry(t) to an integer polynomial in z."""

    if len(entry) != MAX_T_DEGREE + 1:
        raise ValueError("unexpected entry degree")
    if any(entry[index] != entry[MAX_T_DEGREE - index] for index in range(MAX_T_DEGREE + 1)):
        raise AssertionError("M(t)=t^26 M(1/t) failed")
    expression: sp.Expr = sp.Integer(entry[13])
    for power in range(1, 14):
        expression += int(entry[13 + power]) * chebyshev[power]
    return sp.Poly(sp.expand(expression), Z, domain=sp.ZZ)


def sector_charpoly_z(
    coefficients: list[list[list[int]]], basis: SectorBasis, chebyshev: Sequence[sp.Expr]
) -> sp.Poly:
    """Exact generalized characteristic polynomial of t^-13 M(t) in z."""

    size = len(basis.gram)
    block: list[list[sp.Expr]] = []
    for row in range(size):
        output_row: list[sp.Expr] = []
        for column in range(size):
            entry = sp.Poly(0, Z, domain=sp.ZZ)
            for left_state, left_coefficient in basis.vectors[row].items():
                for right_state, right_coefficient in basis.vectors[column].items():
                    entry += left_coefficient * right_coefficient * entry_to_z(
                        coefficients[left_state][right_state], chebyshev
                    )
            output_row.append(entry.as_expr())
        block.append(output_row)
    determinant = sp.Poly(
        (X * sp.diag(*basis.gram) - sp.Matrix(block)).det(method="domain-ge"), X, Z, domain=sp.ZZ
    )
    _content, primitive = determinant.primitive()
    return primitive


def nontrivial_low_factor(polynomial: sp.Poly) -> tuple[sp.Poly, list[dict[str, object]]]:
    """Select the unique non-linear exact factor containing the low branch."""

    content, factors = sp.factor_list(polynomial.as_expr(), X, Z)
    metadata: list[dict[str, object]] = []
    candidates: list[sp.Poly] = []
    for factor, exponent in factors:
        poly = sp.Poly(factor, X, Z, domain=sp.ZZ)
        _factor_content, primitive = poly.primitive()
        metadata.append(
            {
                "x_degree": primitive.degree(X),
                "z_degree": primitive.degree(Z),
                "exponent": int(exponent),
                "terms": len(primitive.terms()),
                "expression_if_linear": str(primitive.as_expr()) if primitive.degree(X) == 1 else None,
            }
        )
        if primitive.degree(X) > 1:
            candidates.append(primitive)
    if content == 0 or len(candidates) != 1:
        raise AssertionError("expected exactly one non-linear low-sector factor")
    return candidates[0], metadata


def restore_t_parameter(polynomial_z: sp.Poly, eigen_degree: int) -> sp.Poly:
    """Return t^(13n) p(x/t^13,t+t^-1) in ZZ[x,t]."""

    expression: sp.Expr = sp.Integer(0)
    for (x_power, z_power), coefficient in polynomial_z.as_dict().items():
        t_power = 13 * eigen_degree - 13 * x_power - z_power
        if t_power < 0:
            raise AssertionError("reciprocal restoration produced a denominator")
        expression += int(coefficient) * X**x_power * T**t_power * (1 + T * T) ** z_power
    output = sp.Poly(sp.expand(expression), X, T, domain=sp.ZZ)
    _content, primitive = output.primitive()
    return primitive


def polynomial_digest(polynomial: sp.Poly) -> str:
    payload = json.dumps(
        [(list(exponents), str(coefficient)) for exponents, coefficient in sorted(polynomial.as_dict().items())],
        separators=(",", ":"),
    )
    return sha256(payload.encode()).hexdigest()


# ---------------------------------------------------------------------------
# Exact interval/Krawczyk arithmetic
# ---------------------------------------------------------------------------

Term = tuple[tuple[int, ...], int]


def polynomial_terms(polynomial: sp.Poly, variable_index: int) -> list[Term]:
    """Embed ZZ[x,t] into ZZ[t,a,b,c,d], placing x at variable_index."""

    output: list[Term] = []
    for (x_power, t_power), coefficient in polynomial.as_dict().items():
        exponents = [0] * 5
        exponents[0] = t_power
        exponents[variable_index] = x_power
        output.append((tuple(exponents), int(coefficient)))
    return output


def differentiate_terms(terms: Iterable[Term], variable_index: int) -> list[Term]:
    output: list[Term] = []
    for exponents, coefficient in terms:
        power = exponents[variable_index]
        if power:
            reduced = list(exponents)
            reduced[variable_index] -= 1
            output.append((tuple(reduced), coefficient * power))
    return output


def evaluate_terms(terms: Iterable[Term], values: Sequence[Fraction]) -> Fraction:
    total = Fraction(0)
    for exponents, coefficient in terms:
        term = Fraction(coefficient)
        for value, power in zip(values, exponents, strict=True):
            if power:
                term *= value**power
        total += term
    return total


def absolute_evaluate_terms(terms: Iterable[Term], upper: Sequence[Fraction]) -> Fraction:
    """Upper-bound |p| over a positive box by evaluating absolute coefficients."""

    total = Fraction(0)
    for exponents, coefficient in terms:
        term = abs(coefficient)
        for value, power in zip(upper, exponents, strict=True):
            if power:
                term *= value**power
        total += term
    return total


def matrix_product(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [
        [sum((left[row][mid] * right[mid][column] for mid in range(len(right))), Fraction(0)) for column in range(len(right[0]))]
        for row in range(len(left))
    ]


def matrix_vector_product(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    return [sum((matrix[row][column] * vector[column] for column in range(len(vector))), Fraction(0)) for row in range(len(matrix))]


def decimal_center(value: str, digits: int = 45) -> Fraction:
    """Truncate a positive decimal guide to a controlled exact rational center."""

    whole, fractional = value.split(".", 1)
    fractional = fractional[:digits].ljust(digits, "0")
    return Fraction(int(whole + fractional), 10**digits)


def mp_fraction(value: mp.mpf, digits: int = 90) -> Fraction:
    return Fraction(mp.nstr(value, digits))


def fraction_decimal(value: Fraction, digits: int = 18) -> str:
    with mp.workdps(digits + 10):
        return mp.nstr(mp.mpf(value.numerator) / value.denominator, digits)


def krawczyk_certificate(terms: Sequence[list[Term]], guide: dict[str, str], contraction_upper: Fraction) -> dict[str, object]:
    """Prove one nonsingular zero in a rational five-variable box exactly."""

    center = [decimal_center(guide[key]) for key in ("t", "a", "b", "c", "d")]
    radii = [T_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS]
    upper = [center[index] + radii[index] for index in range(5)]

    jacobian_terms = [[differentiate_terms(terms[row], column) for column in range(5)] for row in range(5)]
    hessian_terms = [
        [[differentiate_terms(jacobian_terms[row][column], direction) for direction in range(5)] for column in range(5)]
        for row in range(5)
    ]
    f_center = [evaluate_terms(row, center) for row in terms]
    jacobian_center = [[evaluate_terms(jacobian_terms[row][column], center) for column in range(5)] for row in range(5)]

    # This is merely a rational preconditioner.  All later containment tests
    # use its exact decimal-rational representation, not floating arithmetic.
    with mp.workdps(120):
        matrix = mp.matrix([[mp.mpf(value.numerator) / value.denominator for value in row] for row in jacobian_center])
        inverse = matrix**-1
        preconditioner = [[mp_fraction(inverse[row, column]) for column in range(5)] for row in range(5)]

    variation = [
        [
            sum(
                (absolute_evaluate_terms(hessian_terms[row][column][direction], upper) * radii[direction] for direction in range(5)),
                Fraction(0),
            )
            for column in range(5)
        ]
        for row in range(5)
    ]
    correction_matrix = matrix_product(preconditioner, jacobian_center)
    identity_error = [[Fraction(int(row == column)) - correction_matrix[row][column] for column in range(5)] for row in range(5)]
    absolute_preconditioner = [[abs(value) for value in row] for row in preconditioner]
    propagated_variation = matrix_product(absolute_preconditioner, variation)
    contraction = [
        [abs(identity_error[row][column]) + propagated_variation[row][column] for column in range(5)]
        for row in range(5)
    ]
    image_center = [center[index] - value for index, value in enumerate(matrix_vector_product(preconditioner, f_center))]
    containment_ratios = [
        (
            abs(image_center[row] - center[row])
            + sum((contraction[row][column] * radii[column] for column in range(5)), Fraction(0))
        )
        / radii[row]
        for row in range(5)
    ]
    contraction_norm = max(
        sum((contraction[row][column] * radii[column] / radii[row] for column in range(5)), Fraction(0))
        for row in range(5)
    )
    assert all(ratio < COORDINATE_CONTAINMENT_UPPER for ratio in containment_ratios)
    assert contraction_norm < contraction_upper < 1

    return {
        "variables": ["t", "a", "b", "c", "d"],
        "center": [str(value) for value in center],
        "radii": [str(value) for value in radii],
        "preconditioner": [[str(value) for value in row] for row in preconditioner],
        "coordinate_containment_upper_bound": str(COORDINATE_CONTAINMENT_UPPER),
        "contraction_upper_bound": str(contraction_upper),
        "coordinate_containment_ratios_decimal": [fraction_decimal(value, 14) for value in containment_ratios],
        "contraction_norm_decimal": fraction_decimal(contraction_norm, 14),
        "contraction_norm_lt_one": True,
        "nonsingular": True,
    }


# ---------------------------------------------------------------------------
# Exact inertia, ordered branches, and implicit derivative transport
# ---------------------------------------------------------------------------


def fraction_free_inertia(lower_triangle: list[list[int]]) -> int | None:
    """Exact negative-pivot count for a symmetric integer matrix."""

    size = len(lower_triangle)
    work = [row[:] for row in lower_triangle]
    previous = 1
    negatives = 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        if pivot_index + 1 == size:
            break
        for row_index in range(pivot_index + 1, size):
            left = work[row_index][pivot_index]
            for column_index in range(pivot_index + 1, row_index + 1):
                numerator = pivot * work[row_index][column_index] - left * work[column_index][pivot_index]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("nonexact Bareiss division")
                work[row_index][column_index] = quotient
        previous = pivot
    return negatives


def integer_matrix_at_fraction(coefficients: list[list[list[int]]], value: Fraction) -> tuple[list[list[int]], int]:
    """Return scale*M(value) as an integer matrix, with scale=denominator^26."""

    numerator, denominator = value.numerator, value.denominator
    numerator_powers = [numerator**power for power in range(MAX_T_DEGREE + 1)]
    denominator_powers = [denominator ** (MAX_T_DEGREE - power) for power in range(MAX_T_DEGREE + 1)]
    matrix = [
        [
            sum(coefficient * numerator_powers[power] * denominator_powers[power] for power, coefficient in enumerate(entry) if coefficient)
            for entry in row
        ]
        for row in coefficients
    ]
    return matrix, denominator**MAX_T_DEGREE


def integer_sector_blocks(matrix: list[list[int]], bases: Sequence[SectorBasis]) -> list[tuple[list[list[int]], tuple[int, ...]]]:
    output: list[tuple[list[list[int]], tuple[int, ...]]] = []
    for basis in bases:
        size = len(basis.gram)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            for column in range(row + 1):
                value = sum(
                    left * matrix[left_state][right_state] * right
                    for left_state, left in basis.vectors[row].items()
                    for right_state, right in basis.vectors[column].items()
                )
                block[row][column] = block[column][row] = value
        output.append((block, basis.gram))
    return output


def sector_count_below(block: list[list[int]], gram: Sequence[int], scale: int, shift: Fraction) -> int | None:
    size = len(gram)
    shifted = [[0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            value = shift.denominator * block[row][column]
            if row == column:
                value -= shift.numerator * scale * gram[row]
            shifted[row][column] = value
    return fraction_free_inertia(shifted)


def all_sector_counts(blocks: Sequence[tuple[list[list[int]], tuple[int, ...]]], scale: int, shift: Fraction) -> list[int]:
    counts = [sector_count_below(block, gram, scale, shift) for block, gram in blocks]
    if any(value is None for value in counts):
        raise ArithmeticError(f"singular tested inertia shift {shift}")
    return [int(value) for value in counts]


def rational_sqrt_upper(value: Fraction, denominator: int = 10**18) -> Fraction:
    target = (value.numerator * denominator * denominator + value.denominator - 1) // value.denominator
    numerator = math.isqrt(target)
    if numerator * numerator < target:
        numerator += 1
    result = Fraction(numerator, denominator)
    assert result * result >= value
    return result


def full_derivative_bound(coefficients: list[list[list[int]]], upper_t: Fraction) -> Fraction:
    frobenius_square = Fraction(0)
    row_sums: list[Fraction] = []
    for row in coefficients:
        row_sum = Fraction(0)
        for entry in row:
            derivative = sum(
                (Fraction(power * coefficient) * upper_t ** (power - 1) for power, coefficient in enumerate(entry) if power and coefficient),
                Fraction(0),
            )
            frobenius_square += derivative * derivative
            row_sum += derivative
        row_sums.append(row_sum)
    return min(rational_sqrt_upper(frobenius_square), max(row_sums))


def sector_coefficient_blocks(coefficients: list[list[list[int]]], bases: Sequence[SectorBasis]) -> list[list[list[list[int]]]]:
    output: list[list[list[list[int]]]] = []
    for basis in bases:
        size = len(basis.gram)
        block = [[[0] * (MAX_T_DEGREE + 1) for _ in range(size)] for _ in range(size)]
        for row in range(size):
            for column in range(size):
                entry = block[row][column]
                for left_state, left in basis.vectors[row].items():
                    for right_state, right in basis.vectors[column].items():
                        source = coefficients[left_state][right_state]
                        factor = left * right
                        for power, coefficient in enumerate(source):
                            entry[power] += factor * coefficient
        output.append(block)
    return output


def sector_derivative_bound(block: list[list[list[int]]], gram: Sequence[int], upper_t: Fraction) -> Fraction:
    frobenius_square = Fraction(0)
    for row in range(len(gram)):
        for column in range(len(gram)):
            derivative = sum(
                (Fraction(abs(power * coefficient)) * upper_t ** (power - 1) for power, coefficient in enumerate(block[row][column]) if power and coefficient),
                Fraction(0),
            )
            frobenius_square += derivative * derivative / (gram[row] * gram[column])
    return rational_sqrt_upper(frobenius_square)


def terms_two_variables(polynomial: sp.Poly) -> list[tuple[tuple[int, int], int]]:
    return [(tuple(map(int, exponents)), int(coefficient)) for exponents, coefficient in polynomial.as_dict().items()]


def differentiate_two(terms: Iterable[tuple[tuple[int, int], int]], variable: int) -> list[tuple[tuple[int, int], int]]:
    output: list[tuple[tuple[int, int], int]] = []
    for exponents, coefficient in terms:
        power = exponents[variable]
        if power:
            reduced = list(exponents)
            reduced[variable] -= 1
            output.append((tuple(reduced), coefficient * power))
    return output


def evaluate_two(terms: Iterable[tuple[tuple[int, int], int]], eigenvalue: Fraction, parameter: Fraction, absolute: bool = False) -> Fraction:
    total = Fraction(0)
    for (x_power, t_power), coefficient in terms:
        total += (abs(coefficient) if absolute else coefficient) * eigenvalue**x_power * parameter**t_power
    return total


def implicit_derivative_certificate(polynomial: sp.Poly, eigenvalue_center: Fraction, t_center: Fraction) -> dict[str, object]:
    """Exact bound for |lambda'|=|p_t/p_x| on a broad bootstrap tube."""

    terms = terms_two_variables(polynomial)
    p_x = differentiate_two(terms, 0)
    p_t = differentiate_two(terms, 1)
    p_xx = differentiate_two(p_x, 0)
    p_xt = differentiate_two(p_x, 1)
    p_tx = differentiate_two(p_t, 0)
    p_tt = differentiate_two(p_t, 1)
    upper_x = eigenvalue_center + DERIVATIVE_TUBE_RADIUS
    upper_t = t_center + T_RADIUS
    p_x_center = evaluate_two(p_x, eigenvalue_center, t_center)
    p_t_center = evaluate_two(p_t, eigenvalue_center, t_center)
    delta_x = evaluate_two(p_xx, upper_x, upper_t, True) * DERIVATIVE_TUBE_RADIUS + evaluate_two(p_xt, upper_x, upper_t, True) * T_RADIUS
    delta_t = evaluate_two(p_tx, upper_x, upper_t, True) * DERIVATIVE_TUBE_RADIUS + evaluate_two(p_tt, upper_x, upper_t, True) * T_RADIUS
    lower_x = abs(p_x_center) - delta_x
    upper_t_derivative = abs(p_t_center) + delta_t
    if lower_x <= 0:
        raise ArithmeticError("implicit-derivative tube reaches p_x=0")
    derivative_upper = upper_t_derivative / lower_x
    assert derivative_upper < DERIVATIVE_BRANCH_UPPER
    return {
        "tube_eigen_radius": str(DERIVATIVE_TUBE_RADIUS),
        "derivative_upper_bound": str(DERIVATIVE_BRANCH_UPPER),
        "refined_motion_upper": str(REFINED_BRANCH_MOTION),
        "actual_derivative_bound_decimal": fraction_decimal(derivative_upper, 14),
        "p_x_stays_nonzero": True,
    }


def linear_factor_excluded(linear_polynomial: sp.Poly | None, eigenvalue_center: Fraction, t_center: Fraction) -> bool:
    if linear_polynomial is None:
        return True
    # A direct center evaluation alone is not enough, so compare its rational
    # root in x at the fixed parameter to the whole initial low-eigenvalue box.
    x_poly = sp.Poly(linear_polynomial.as_expr().subs(T, sp.Rational(t_center.numerator, t_center.denominator)), X, domain=sp.QQ)
    root = -x_poly.nth(0) / x_poly.nth(1)
    return not (sp.Rational(eigenvalue_center.numerator, eigenvalue_center.denominator) - sp.Rational(INITIAL_EIGEN_RADIUS.numerator, INITIAL_EIGEN_RADIUS.denominator) <= root <= sp.Rational(eigenvalue_center.numerator, eigenvalue_center.denominator) + sp.Rational(INITIAL_EIGEN_RADIUS.numerator, INITIAL_EIGEN_RADIUS.denominator))


def ordering_certificate(
    spec: CrossingSpec,
    coefficients: list[list[list[int]]],
    bases: Sequence[SectorBasis],
    sector_blocks_coefficients: Sequence[list[list[list[int]]]],
    factor_polynomials_t: dict[int, sp.Poly],
    linear_polynomials_t: dict[int, sp.Poly | None],
) -> dict[str, object]:
    """Pin a,b,c,d to mu_0,mu_1,mu_2,mu_{r+3} throughout the fine bracket."""

    guide = GUIDES[spec.number]
    centers = {key: decimal_center(guide[key]) for key in ("t", "a", "b", "c", "d")}
    matrix, scale = integer_matrix_at_fraction(coefficients, centers["t"])
    blocks = integer_sector_blocks(matrix, bases)

    threshold_counts: list[list[int]] = []
    for threshold, expected in zip(spec.thresholds, spec.threshold_sector_counts, strict=True):
        counts = all_sector_counts(blocks, scale, threshold)
        assert counts == list(expected)
        threshold_counts.append(counts)
        # A full 10^-7 spectral gap around every separator is certified at
        # the center, then protected by the derivative bounds below.
        assert all_sector_counts(blocks, scale, threshold - ORDER_MARGIN) == counts
        assert all_sector_counts(blocks, scale, threshold + ORDER_MARGIN) == counts

    selected = (("a", 0, 0), ("b", 1, 5), ("c", 2, 7), ("d", spec.d_index, spec.d_sector))
    local_windows: list[dict[str, object]] = []
    for key, global_index, sector in selected:
        lower = centers[key] - INITIAL_EIGEN_RADIUS
        upper = centers[key] + INITIAL_EIGEN_RADIUS
        lower_counts = all_sector_counts(blocks, scale, lower)
        upper_counts = all_sector_counts(blocks, scale, upper)
        assert sum(lower_counts) == global_index
        assert sum(upper_counts) == global_index + 1
        assert upper_counts[sector] == lower_counts[sector] + 1
        assert all(upper_counts[index] == lower_counts[index] for index in range(8) if index != sector)
        assert linear_factor_excluded(linear_polynomials_t[sector], centers[key], centers["t"])
        local_windows.append(
            {
                "variable": key,
                "global_index": global_index,
                "sector": sector,
                "lower": str(lower),
                "upper": str(upper),
            }
        )

    # The selected roots lie in the corresponding threshold slots of the
    # Krawczyk box.  A forced-value interval below will show that *every*
    # forced equality in the inherited bracket uses this selected d branch.
    q0, q1, q2, qd_low, qd_high = spec.thresholds
    assert centers["a"] + EIGEN_RADIUS < q0
    assert q0 < centers["b"] - EIGEN_RADIUS < centers["b"] + EIGEN_RADIUS < q1
    assert q1 < centers["c"] - EIGEN_RADIUS < centers["c"] + EIGEN_RADIUS < q2
    assert qd_low < centers["d"] - EIGEN_RADIUS < centers["d"] + EIGEN_RADIUS < qd_high
    assert spec.lower <= centers["t"] <= spec.upper

    full_bound = full_derivative_bound(coefficients, spec.upper)
    assert full_bound * T_RADIUS < ORDER_MARGIN
    selected_sector_bounds: dict[str, Fraction] = {}
    implicit: dict[str, dict[str, object]] = {}
    for key, _global_index, sector in selected:
        sector_bound = sector_derivative_bound(sector_blocks_coefficients[sector], bases[sector].gram, spec.upper)
        assert sector_bound * T_RADIUS < ORDER_MARGIN
        selected_sector_bounds[key] = sector_bound
        implicit[key] = implicit_derivative_certificate(factor_polynomials_t[sector], centers[key], centers["t"])
        assert REFINED_BRANCH_MOTION < EIGEN_RADIUS
    # On the inherited bracket, all four selected branches are within this
    # radius of their centers.  Positivity then bounds the forced value
    # f=bc/a exactly, not just at a numerical guide point.
    assert all(centers[key] > REFINED_BRANCH_MOTION for key in ("a", "b", "c"))
    forced_lower = (centers["b"] - REFINED_BRANCH_MOTION) * (centers["c"] - REFINED_BRANCH_MOTION) / (centers["a"] + REFINED_BRANCH_MOTION)
    forced_upper = (centers["b"] + REFINED_BRANCH_MOTION) * (centers["c"] + REFINED_BRANCH_MOTION) / (centers["a"] - REFINED_BRANCH_MOTION)
    assert qd_low < forced_lower <= forced_upper < qd_high


    # The inherited exact c-count theorem gives existence inside the original
    # fine bracket.  It is contained in the Krawczyk t-box, while the bounds
    # above identify every forced equality there with the selected d branch.
    assert centers["t"] - T_RADIUS <= spec.lower <= spec.upper <= centers["t"] + T_RADIUS

    return {
        "forced_branch_identified": True,
        "thresholds": [str(value) for value in spec.thresholds],
        "global_counts": list(spec.threshold_global_counts),
        "sector_counts": threshold_counts,
        "order_margin": str(ORDER_MARGIN),
        "full_derivative_upper": str(full_bound),
        "selected_sector_derivative_uppers": {key: str(value) for key, value in selected_sector_bounds.items()},
        "local_eigenvalue_windows": local_windows,
        "initial_eigenvalue_radius": str(INITIAL_EIGEN_RADIUS),
        "implicit_branch_certificates": implicit,
        "forced_value_interval_on_fine_bracket": {
            "lower": str(forced_lower),
            "upper": str(forced_upper),
            "selected_branch_motion_upper": str(REFINED_BRANCH_MOTION),
            "target_band_lower": str(qd_low),
            "target_band_upper": str(qd_high),
            "strictly_inside_target_band": True,
        },
    }


# ---------------------------------------------------------------------------
# Result construction
# ---------------------------------------------------------------------------


def interval_json(lower: Fraction, upper: Fraction) -> dict[str, str]:
    return {
        "lower": str(lower),
        "upper": str(upper),
        "width": str(upper - lower),
        "lower_decimal": fraction_decimal(lower, 24),
        "upper_decimal": fraction_decimal(upper, 24),
    }


def build_polynomial_data(
    coefficients: list[list[list[int]]], bases: Sequence[SectorBasis]
) -> tuple[dict[int, sp.Poly], dict[int, sp.Poly | None], dict[str, object]]:
    chebyshev = reciprocal_chebyshev()
    z_characteristic = {sector: sector_charpoly_z(coefficients, bases[sector], chebyshev) for sector in (0, 4, 5, 6, 7)}
    low_factors: dict[int, sp.Poly] = {}
    linear_t: dict[int, sp.Poly | None] = {}
    factor_metadata: dict[str, object] = {}
    for sector, polynomial in z_characteristic.items():
        low_factor, factors = nontrivial_low_factor(polynomial)
        low_factors[sector] = low_factor
        linear = next((sp.Poly(item["expression_if_linear"], X, Z, domain=sp.ZZ) for item in factors if item["expression_if_linear"] is not None), None)
        linear_t[sector] = restore_t_parameter(linear, 1) if linear is not None else None
        factor_metadata[str(sector)] = {
            "characteristic_x_degree": polynomial.degree(X),
            "characteristic_z_degree": polynomial.degree(Z),
            "characteristic_terms": len(polynomial.terms()),
            "factors": factors,
            "low_factor_x_degree": low_factor.degree(X),
            "low_factor_z_degree": low_factor.degree(Z),
        }

    eigen_degrees = {0: 14, 4: 9, 5: 9, 6: 6, 7: 6}
    restored = {sector: restore_t_parameter(low_factors[sector], eigen_degrees[sector]) for sector in low_factors}
    for sector, polynomial in restored.items():
        factor_metadata[str(sector)].update(
            {
                "restored_t_x_degree": polynomial.degree(X),
                "restored_t_degree": polynomial.degree(T),
                "restored_t_terms": len(polynomial.terms()),
                "restored_t_sha256": polynomial_digest(polynomial),
            }
        )
    expected = {
        0: (14, 133, 14, 315, 2040),
        4: (9, 90, 9, 207, 880),
        5: (9, 81, 9, 198, 797),
        6: (6, 59, 6, 137, 382),
        7: (6, 55, 6, 133, 356),
    }
    for sector, values in expected.items():
        metadata = factor_metadata[str(sector)]
        assert (
            metadata["low_factor_x_degree"],
            metadata["low_factor_z_degree"],
            metadata["restored_t_x_degree"],
            metadata["restored_t_degree"],
            metadata["restored_t_terms"],
        ) == values
    return restored, linear_t, factor_metadata


def formal_resultant_wall_metadata() -> dict[str, object]:
    """An honestly recorded, non-default expansion attempt; no fake polynomial."""

    return {
        "claim_tag": "[UNRESOLVED]",
        "parameter": "z=t+t^-1",
        "first_stage": "U(z,y)=Res_b(P5(b,z), b^6 P7(y/b,z))",
        "input_x_degrees": [9, 6],
        "input_term_counts": [309, 134],
        "method": "SymPy exact multivariate resultant",
        "status": "wall",
        "wall_seconds": 1200,
        "observed": "No completed U(z,y) was produced before the 1,200-second wall; V and the final univariate R were not attempted after that wall.",
        "consequence": "No minimal polynomial or Sturm count for the full scanned region is claimed.",
    }


def build_result() -> dict[str, object]:
    started = time.perf_counter()
    coefficients, degree, shift = polynomial_matrix()
    bases = make_sector_bases()
    restored, linear_t, factor_metadata = build_polynomial_data(coefficients, bases)
    sector_blocks = sector_coefficient_blocks(coefficients, bases)

    term_cache = {(sector, variable): polynomial_terms(restored[sector], variable) for sector in restored for variable in range(1, 5)}
    crossing_records: list[dict[str, object]] = []
    for spec in CROSSINGS:
        terms = [
            term_cache[(0, 1)],
            term_cache[(5, 2)],
            term_cache[(7, 3)],
            term_cache[(spec.d_sector, 4)],
            [((0, 0, 1, 1, 0), 1), ((0, 1, 0, 0, 1), -1)],
        ]
        krawczyk = krawczyk_certificate(terms, GUIDES[spec.number], spec.contraction_upper)
        ordering = ordering_certificate(spec, coefficients, bases, sector_blocks, restored, linear_t)
        crossing_records.append(
            {
                "name": f"forced_value_crossing_{spec.number}",
                "claim_tag": "[THEOREM]",
                "fine_bracket": interval_json(spec.lower, spec.upper),
                "inherited_c_jump_existence": {
                    "claim_tag": "[THEOREM, inherited from e81]",
                    "source": "proofs/spectral_interval2.md and results/spectral/interval2.json",
                    "c_counts_left_to_right": [spec.d_index, spec.d_index + 1],
                    "use_here": "It supplies at least one forced crossing in this fine bracket; the new Krawczyk certificate proves that such a crossing is unique.",
                },
                "ordered_sectors": {"lambda0": 0, "lambda1": 5, "lambda2": 7, f"lambda{spec.d_index}": spec.d_sector},
                "exact_algebraic_representation": {
                    "kind": "unique nonsingular real zero of five integer polynomials in a rational box",
                    "variables": ["t", "a", "b", "c", "d"],
                    "system": ["P0(a,t)", "P5(b,t)", "P7(c,t)", f"P{spec.d_sector}(d,t)", "b*c-a*d"],
                    "minimal_polynomial": None,
                    "transverse_sector_product_intersection": True,
                },
                "krawczyk": krawczyk,
                "ordering": ordering,
                "numeric_guide": {
                    "method": "mpmath 400-dps sector Newton; guide only, never a certificate authority",
                    "t_decimal": GUIDES[spec.number]["t"],
                    "a_decimal": GUIDES[spec.number]["a"],
                    "b_decimal": GUIDES[spec.number]["b"],
                    "c_decimal": GUIDES[spec.number]["c"],
                    "d_decimal": GUIDES[spec.number]["d"],
                },
                "conclusion": "Exactly one forced-value crossing lies in this inherited fine bracket; its t-coordinate is the unique algebraic t-coordinate represented by the stored system-and-box certificate.",
            }
        )

    return {
        "provenance": {
            "script": "experiments/e96_crossing_alg.py",
            "generated_utc": datetime.now(UTC).isoformat(),
            "interpreter": sys.executable,
            "method": "Exact 64-state polynomial representative; reciprocal z=t+t^-1 sector characteristic polynomials; exact factorization; rational Krawczyk interval certificates; exact sector inertia and implicit polynomial derivative transport. mpmath only supplies rational preconditioners verified afterward by Fraction arithmetic.",
        },
        "data": {
            "status": "PARTIAL_THEOREM",
            "claim_tags": {
                "bracket_uniqueness_and_exact_isolation": "[THEOREM]",
                "finite_symbolic_factorizations_and_inertia": "[COMPUTATION]",
                "reciprocal_parameter_lemma": "[LEMMA]",
                "full_scan_exhaustiveness_and_minimal_polynomials": "[UNRESOLVED]",
            },
            "certified_bracket_crossing_count": 4,
            "full_scan_crossing_count": None,
            "minimal_polynomial_status": "unresolved_after_resultant_wall",
            "scan_region_t": interval_json(Fraction(1, 4), Fraction(7, 20)),
            "matrix": {
                "dimension": DIMENSION,
                "degree": MAX_T_DEGREE,
                "shift": shift,
                "scalar_degree": degree,
                "reciprocity": "M(t)=t^26 M(1/t); hence t^-13 M(t) is polynomial in z=t+t^-1",
                "sector_dimensions": [len(basis.gram) for basis in bases],
            },
            "sector_polynomials": factor_metadata,
            "crossings": crossing_records,
            "resultant_attempt": formal_resultant_wall_metadata(),
            "full_scan_exhaustiveness": {
                "claim_tag": "[UNRESOLVED]",
                "statement": "The four stored fine brackets are complete internally, but no theorem excludes extra forced crossings in the unbracketed portions of [1/4,7/20].",
                "reason": "The required global elimination polynomial R was not expanded/factored after the recorded first-resultant wall.",
            },
            "unresolved": [
                "Minimal polynomials of the four isolated t-coordinates were not computed.",
                "A Sturm root count of a global R on [1/4,7/20] was not computed.",
                "No absence theorem is claimed outside the four inherited fine brackets.",
            ],
        },
        "checks": [
            {
                "name": "clean_room_verifier",
                "command": ".venv/bin/python tests/test_crossing_alg.py",
                "status": "run separately after result emission",
            },
            {"name": "producer_elapsed_seconds", "value": time.perf_counter() - started},
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--print", action="store_true", dest="print_result", help="also print the generated JSON")
    arguments = parser.parse_args()
    result = build_result()
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if arguments.print_result:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"wrote {RESULT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
