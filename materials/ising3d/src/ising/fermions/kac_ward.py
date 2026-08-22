"""Kac--Ward control calculations and a scoped three-dimensional obstruction.

The exact two-dimensional determinant is evaluated in the cyclotomic integer
ring ``Z[zeta_8]``.  Four finite-field embeddings recover every coefficient in
the basis ``1, zeta_8, zeta_8**2, zeta_8**3``; a Hadamard bound turns the
modular reconstruction into an exact certificate rather than a probabilistic
identity test.

For the three-dimensional calculation, a scalar, translation-invariant turn
rule has weights ``U(d, d')`` for the 30 allowed ordered direction pairs.  The
strict cubic-group fixed slice has two weights: ``a`` for a straight step and
``b`` for an orthogonal turn.  The low-order trace equations prove an
obstruction only on that two-parameter slice.  Helpers for constructing the
full, generically gauge-fixed 30-weight finite-box system are provided so its
larger unresolved scope is recorded rather than silently conflated with the
proved result.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb
from typing import Iterable, Mapping, Sequence

import mpmath as mp
import sympy as sp
from sympy import GF
from sympy.ntheory import primitive_root
from sympy.polys.matrices import DomainMatrix

from ising.exact_enumeration import even_subgraph_polynomial
from ising.lattices import Lattice, cubic


DIRECTIONS_3D: tuple[tuple[int, int, int], ...] = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
)
DIRECTION_LABELS: tuple[str, ...] = ("+x", "-x", "+y", "-y", "+z", "-z")
OPPOSITE_DIRECTION: tuple[int, ...] = (1, 0, 3, 2, 5, 4)
DIRECTION_INDEX = {direction: index for index, direction in enumerate(DIRECTIONS_3D)}
ALLOWED_DIRECTION_PAIRS: tuple[tuple[int, int], ...] = tuple(
    (direction, next_direction)
    for direction in range(6)
    for next_direction in range(6)
    if next_direction != OPPOSITE_DIRECTION[direction]
)
# A directed spanning tree in the graph of direction states.  On the open
# chart where these five entries are nonzero, direction-state similarity can
# set them to one.  A common gauge factor is ineffective, hence 6 - 1 = 5.
DIRECTION_GAUGE_TREE: tuple[tuple[int, int], ...] = (
    (0, 2),
    (0, 3),
    (0, 4),
    (0, 5),
    (2, 1),
)

# Exact simple-cubic infinite-volume coefficients [v^k] log P per site,
# independently established by the finite-lattice method in
# results/series/sc_ht_free_energy.json.
BULK_SC_LOG_P: dict[int, Fraction] = {
    4: Fraction(3),
    6: Fraction(22),
    8: Fraction(375, 2),
}


@dataclass(frozen=True)
class DirectedEdge:
    """One orientation of a lattice bond."""

    tail: int
    head: int
    direction: tuple[int, ...]
    reverse: int
    bond: int


@dataclass(frozen=True)
class ExactPolynomialCertificate:
    """Certificate attached to an exact ``Z[zeta_8]`` reconstruction."""

    directed_edge_count: int
    primes: tuple[int, ...]
    modulus: int
    coefficient_bound: int
    max_nonrational_component: int


@dataclass(frozen=True)
class FullWeightSystem:
    """Finite-box equations for the translation-invariant 30-weight ansatz."""

    variables: tuple[sp.Symbol, ...]
    all_symbols: tuple[tuple[tuple[int, int], sp.Symbol], ...]
    gauge_tree: tuple[tuple[int, int], ...]
    labels: tuple[tuple[tuple[int, int, int], int], ...]
    equations: tuple[sp.Expr, ...]
    pre_gauge_unknown_count: int = 30
    post_gauge_unknown_count: int = 25


def _require_open_hypercubic(lattice: Lattice, dimension: int) -> None:
    if lattice.dim != dimension:
        raise ValueError(f"expected a {dimension}D lattice, got dimension {lattice.dim}")
    if any(lattice.periodic):
        raise ValueError("the geometric turn construction in this module requires free boundaries")


def directed_edges(lattice: Lattice) -> tuple[DirectedEdge, ...]:
    """Return both orientations of every bond, retaining bond identity.

    The reverse exclusion is by bond identity, which is also correct for the
    repository's parallel-bond convention.  Geometric routines below are
    intentionally restricted to free-boundary hypercubic lattices.
    """

    result: list[DirectedEdge] = []
    for bond_index, (tail, head) in enumerate(lattice.bonds):
        tail_coord = lattice.coord(tail)
        head_coord = lattice.coord(head)
        direction = tuple(b - a for a, b in zip(tail_coord, head_coord, strict=True))
        forward = 2 * bond_index
        reverse = forward + 1
        result.append(DirectedEdge(tail, head, direction, reverse, bond_index))
        result.append(
            DirectedEdge(head, tail, tuple(-component for component in direction), forward, bond_index)
        )
    return tuple(result)


@lru_cache(maxsize=None)
def _transition_indices(lattice: Lattice) -> tuple[tuple[int, ...], ...]:
    edges = directed_edges(lattice)
    outgoing: list[list[int]] = [[] for _ in range(lattice.n_sites)]
    for index, edge in enumerate(edges):
        outgoing[edge.tail].append(index)
    return tuple(
        tuple(next_index for next_index in outgoing[edge.head] if next_index != edge.reverse)
        for edge in edges
    )


def _quarter_turns_2d(
    direction: tuple[int, ...], next_direction: tuple[int, ...]
) -> int:
    dx, dy = direction
    ex, ey = next_direction
    dot = dx * ex + dy * ey
    cross = dx * ey - dy * ex
    if dot == 1 and cross == 0:
        return 0
    if dot == 0 and cross in (-1, 1):
        return cross
    raise ValueError(f"not an allowed square-lattice turn: {direction} -> {next_direction}")


@lru_cache(maxsize=None)
def _turn_table_2d(lattice: Lattice) -> tuple[tuple[tuple[int, int], ...], ...]:
    _require_open_hypercubic(lattice, 2)
    edges = directed_edges(lattice)
    return tuple(
        tuple(
            (next_index, _quarter_turns_2d(edge.direction, edges[next_index].direction))
            for next_index in next_indices
        )
        for edge, next_indices in zip(edges, _transition_indices(lattice), strict=True)
    )


def polynomial_square(polynomial: Sequence[int]) -> list[int]:
    """Square an integer coefficient list exactly."""

    result = [0] * (2 * len(polynomial) - 1)
    for left_degree, left in enumerate(polynomial):
        if not left:
            continue
        for right_degree, right in enumerate(polynomial):
            if right:
                result[left_degree + right_degree] += int(left) * int(right)
    return result


def formal_log_coefficients(
    polynomial: Sequence[int | Fraction], max_order: int
) -> list[Fraction]:
    """Return coefficients of ``log(polynomial)`` through ``max_order``.

    ``polynomial[0]`` must be one.  The recurrence comes from ``P' = P (log P)'``.
    """

    if not polynomial or polynomial[0] != 1:
        raise ValueError("formal logarithm requires constant coefficient one")
    source = [Fraction(value) for value in polynomial[: max_order + 1]]
    source.extend([Fraction(0)] * (max_order + 1 - len(source)))
    result = [Fraction(0)] * (max_order + 1)
    for degree in range(1, max_order + 1):
        numerator = Fraction(degree) * source[degree]
        numerator -= sum(
            Fraction(lower) * result[lower] * source[degree - lower]
            for lower in range(1, degree)
        )
        result[degree] = numerator / degree
    return result


def _modular_inverse_matrix(matrix: list[list[int]], modulus: int) -> list[list[int]]:
    size = len(matrix)
    augmented = [
        [entry % modulus for entry in row]
        + [int(row_index == column) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if augmented[row][column] % modulus
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, modulus)
        augmented[column] = [(entry * inverse) % modulus for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column] % modulus
            if factor:
                augmented[row] = [
                    (left - factor * right) % modulus
                    for left, right in zip(augmented[row], augmented[column], strict=True)
                ]
    return [row[size:] for row in augmented]


def _cyclotomic_basis_from_embeddings(
    values: Sequence[int], roots: Sequence[int], modulus: int
) -> tuple[int, int, int, int]:
    vandermonde = [
        [pow(root, exponent, modulus) for exponent in range(4)] for root in roots
    ]
    inverse = _modular_inverse_matrix(vandermonde, modulus)
    return tuple(
        sum(entry * value for entry, value in zip(row, values, strict=True)) % modulus
        for row in inverse
    )  # type: ignore[return-value]


def _next_prime_one_mod_eight(cursor: int) -> int:
    candidate = int(sp.nextprime(cursor))
    while candidate % 8 != 1:
        candidate = int(sp.nextprime(candidate))
    return candidate


def _kw_charpoly_mod(
    turn_table: tuple[tuple[tuple[int, int], ...], ...], modulus: int, root: int
) -> list[int]:
    size = len(turn_table)
    inverse_root = pow(root, -1, modulus)
    rows = [[0] * size for _ in range(size)]
    for edge, transitions in enumerate(turn_table):
        for next_edge, quarter_turn in transitions:
            if quarter_turn == 0:
                value = 1
            elif quarter_turn == 1:
                value = root
            else:
                value = inverse_root
            rows[edge][next_edge] = value
    domain = GF(modulus)
    coefficients = DomainMatrix(rows, (size, size), domain).charpoly()
    return [int(value) % modulus for value in coefficients]


def _charpoly_coefficient_bound(size: int, max_row_entries: int) -> int:
    """Integer upper bound for every cyclotomic-basis coefficient.

    Each conjugate of a degree-k characteristic coefficient is a sum of
    ``binomial(size, k)`` principal minors.  Hadamard bounds every such minor
    by ``max_row_entries**(k/2)``.  The inverse four-embedding transform has
    row l1 norm one, so the same bound applies to each ``Z[zeta_8]`` basis
    coefficient.
    """

    degree = max(1, max_row_entries)
    return max(
        comb(size, order) * degree ** ((order + 1) // 2)
        for order in range(size + 1)
    )


@lru_cache(maxsize=None)
def exact_kac_ward_polynomial(
    lattice: Lattice,
) -> tuple[tuple[int, ...], ExactPolynomialCertificate]:
    """Compute ``det(I - v U)`` exactly for a free square lattice.

    The returned coefficients are independently reconstructed in
    ``Z[zeta_8]``.  Non-rational components would raise rather than being
    rounded away.
    """

    turn_table = _turn_table_2d(lattice)
    size = len(turn_table)
    max_row_entries = max((len(row) for row in turn_table), default=0)
    bound = _charpoly_coefficient_bound(size, max_row_entries)

    residues = [[[0 for _ in range(4)] for _ in range(size + 1)]]
    # Keep a flat mutable copy; the one-element wrapper avoids special casing
    # replacement while the CRT modulus changes.
    reconstructed = residues[0]
    modulus_product = 1
    primes: list[int] = []
    cursor = 1_000_000_000
    while modulus_product <= 2 * bound:
        prime = _next_prime_one_mod_eight(cursor)
        cursor = prime
        generator = int(primitive_root(prime))
        zeta = pow(generator, (prime - 1) // 8, prime)
        roots = tuple(pow(zeta, odd, prime) for odd in (1, 3, 5, 7))
        embedded = [_kw_charpoly_mod(turn_table, prime, root) for root in roots]
        inverse_old_modulus = pow(modulus_product, -1, prime)
        for coefficient_index in range(size + 1):
            values = [embedding[coefficient_index] for embedding in embedded]
            basis_residues = _cyclotomic_basis_from_embeddings(values, roots, prime)
            for basis_index, new_residue in enumerate(basis_residues):
                old_residue = reconstructed[coefficient_index][basis_index]
                correction = (
                    (new_residue - old_residue) * inverse_old_modulus
                ) % prime
                reconstructed[coefficient_index][basis_index] = (
                    old_residue + modulus_product * correction
                )
        modulus_product *= prime
        primes.append(prime)

    signed_basis: list[tuple[int, int, int, int]] = []
    midpoint = modulus_product // 2
    for coefficient in reconstructed:
        signed_basis.append(
            tuple(value - modulus_product if value > midpoint else value for value in coefficient)
        )
    max_nonrational = max(
        (abs(value) for coefficient in signed_basis for value in coefficient[1:]),
        default=0,
    )
    if max_nonrational:
        raise ArithmeticError(
            "Kac--Ward characteristic polynomial retained a non-rational cyclotomic component"
        )
    coefficients = tuple(coefficient[0] for coefficient in signed_basis)
    certificate = ExactPolynomialCertificate(
        directed_edge_count=size,
        primes=tuple(primes),
        modulus=modulus_product,
        coefficient_bound=bound,
        max_nonrational_component=max_nonrational,
    )
    return coefficients, certificate


def kac_ward_trace_basis(lattice: Lattice, length: int) -> tuple[int, int, int, int]:
    """Return ``Tr(U**length)`` in the basis ``1,zeta_8,zeta_8^2,zeta_8^3``."""

    if length < 1:
        raise ValueError("length must be positive")
    table = _turn_table_2d(lattice)
    residue_counts = [0] * 8
    for start in range(len(table)):
        states: dict[tuple[int, int], int] = {(start, 0): 1}
        for _ in range(length):
            next_states: dict[tuple[int, int], int] = {}
            for (edge, exponent), count in states.items():
                for next_edge, quarter_turn in table[edge]:
                    key = (next_edge, (exponent + quarter_turn) % 8)
                    next_states[key] = next_states.get(key, 0) + count
            states = next_states
        for residue in range(8):
            residue_counts[residue] += states.get((start, residue), 0)
    return tuple(
        residue_counts[index] - residue_counts[index + 4] for index in range(4)
    )  # type: ignore[return-value]


def kac_ward_determinant(lattice: Lattice, value, *, dps: int = 80):
    """Evaluate the directed-edge Kac--Ward determinant with ``mpmath``."""

    table = _turn_table_2d(lattice)
    with mp.workdps(dps):
        argument = mp.mpf(value)
        zeta = mp.e ** (mp.j * mp.pi / 4)
        inverse_zeta = 1 / zeta
        matrix = mp.eye(len(table))
        for edge, transitions in enumerate(table):
            for next_edge, quarter_turn in transitions:
                phase = 1 if quarter_turn == 0 else (zeta if quarter_turn == 1 else inverse_zeta)
                matrix[edge, next_edge] -= argument * phase
        return +mp.det(matrix)


def trace_type_counts(lattice: Lattice, length: int) -> dict[int, int]:
    """Count closed non-backtracking walks by their number of straight transitions."""

    if length < 1:
        raise ValueError("length must be positive")
    _require_open_hypercubic(lattice, lattice.dim)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    counts: dict[int, int] = {}
    for start in range(len(edges)):
        states: dict[tuple[int, int], int] = {(start, 0): 1}
        for _ in range(length):
            next_states: dict[tuple[int, int], int] = {}
            for (edge_index, straight_count), multiplicity in states.items():
                direction = edges[edge_index].direction
                for next_edge in transitions[edge_index]:
                    is_straight = int(direction == edges[next_edge].direction)
                    key = (next_edge, straight_count + is_straight)
                    next_states[key] = next_states.get(key, 0) + multiplicity
            states = next_states
        for (edge_index, straight_count), multiplicity in states.items():
            if edge_index == start:
                counts[straight_count] = counts.get(straight_count, 0) + multiplicity
    return dict(sorted(counts.items()))


def closed_nonbacktracking_walk_count(lattice: Lattice, length: int) -> int:
    """Return the unweighted trace of the non-backtracking transition matrix."""

    return sum(trace_type_counts(lattice, length).values())


@lru_cache(maxsize=None)
def bulk_cubic_trace_monomials(length: int) -> dict[int, int]:
    """Count rooted closed cubic-lattice walks per site by straight transitions."""

    if length < 1:
        raise ValueError("length must be positive")
    totals: dict[int, int] = {}
    for first in range(6):
        initial_displacement = DIRECTIONS_3D[first]
        states: dict[tuple[int, tuple[int, int, int], int], int] = {
            (first, initial_displacement, 0): 1
        }
        for _ in range(1, length):
            next_states: dict[tuple[int, tuple[int, int, int], int], int] = {}
            for (last, displacement, straight_count), multiplicity in states.items():
                for next_direction in range(6):
                    if next_direction == OPPOSITE_DIRECTION[last]:
                        continue
                    step = DIRECTIONS_3D[next_direction]
                    next_displacement = tuple(
                        displacement[axis] + step[axis] for axis in range(3)
                    )
                    key = (
                        next_direction,
                        next_displacement,
                        straight_count + int(next_direction == last),
                    )
                    next_states[key] = next_states.get(key, 0) + multiplicity
            states = next_states
        for (last, displacement, straight_count), multiplicity in states.items():
            if displacement != (0, 0, 0) or first == OPPOSITE_DIRECTION[last]:
                continue
            cyclic_straight = straight_count + int(last == first)
            totals[cyclic_straight] = totals.get(cyclic_straight, 0) + multiplicity
    return dict(sorted(totals.items()))


def _trace_expression(counts: Mapping[int, int], length: int, straight, turn):
    return sp.expand(
        sum(
            multiplicity * straight**straight_count * turn ** (length - straight_count)
            for straight_count, multiplicity in counts.items()
        )
    )


def strict_cubic_equations(
    log_p_coefficients: Mapping[int, Fraction] = BULK_SC_LOG_P,
) -> tuple[sp.Symbol, sp.Symbol, dict[int, sp.Expr], dict[int, Fraction]]:
    """Build the strict cubic-covariant scalar equations.

    Proper cubic rotations are transitive on straight ordered direction pairs
    and also on orthogonal ordered pairs, leaving two scalar weights ``a,b``.
    Matching ``log det(I-vU)=2 log P`` requires
    ``Tr(U**k) = -2*k*[v**k]log(P)``.
    """

    straight, turn = sp.symbols("a b")
    equations: dict[int, sp.Expr] = {}
    required: dict[int, Fraction] = {}
    for length in sorted(log_p_coefficients):
        required[length] = -2 * length * Fraction(log_p_coefficients[length])
        trace = _trace_expression(
            bulk_cubic_trace_monomials(length), length, straight, turn
        )
        equations[length] = sp.expand(trace - required[length])
    return straight, turn, equations, required


def finite_strict_cubic_equations(
    lattice: Lattice, orders: Iterable[int] = (4, 6, 8)
) -> tuple[sp.Symbol, sp.Symbol, dict[int, sp.Expr], dict[int, Fraction]]:
    """Build the analogous exact equations for one finite free cubic box."""

    _require_open_hypercubic(lattice, 3)
    selected = tuple(sorted(set(int(order) for order in orders)))
    polynomial = even_subgraph_polynomial(lattice)
    log_p = formal_log_coefficients(polynomial, max(selected))
    straight, turn = sp.symbols("a b")
    equations: dict[int, sp.Expr] = {}
    required: dict[int, Fraction] = {}
    for length in selected:
        required[length] = -2 * length * log_p[length]
        trace = _trace_expression(trace_type_counts(lattice, length), length, straight, turn)
        equations[length] = sp.expand(trace - required[length])
    return straight, turn, equations, required


def fit_strict_cubic_candidate(lattice: Lattice) -> tuple[sp.Expr, sp.Expr]:
    """Fit the strict two-weight ansatz at orders four and six.

    The branch ``b=zeta_8`` is selected.  On ``3x3x2`` free, this gives
    ``a=sqrt((39-32*i)/7)`` and matches the determinant through degree six.
    """

    straight_symbol, turn_symbol, equations, required = finite_strict_cubic_equations(
        lattice, (4, 6)
    )
    del straight_symbol, turn_symbol
    turn = (1 + sp.I) / sp.sqrt(2)
    fourth_trace = _trace_expression(trace_type_counts(lattice, 4), 4, 0, turn)
    if sp.simplify(fourth_trace - required[4]) != 0:
        raise ValueError("the zeta_8 turn branch does not solve the fourth-order equation")
    sixth_counts = trace_type_counts(lattice, 6)
    if any(power not in (0, 2) for power in sixth_counts) or not sixth_counts.get(2):
        raise ValueError("the sixth-order equation does not determine a squared straight weight")
    constant = sixth_counts.get(0, 0) * turn**6
    coefficient = sixth_counts[2] * turn**4
    straight_squared = sp.simplify((required[6] - constant) / coefficient)
    straight = sp.sqrt(straight_squared)
    assert sp.simplify(equations[6].subs({sp.Symbol("a"): straight, sp.Symbol("b"): turn})) == 0
    return straight, turn


def strict_cubic_determinant_prefix(
    lattice: Lattice, straight, turn, max_order: int
) -> tuple[sp.Expr, ...]:
    """Compute coefficients of ``det(I-vU)`` from exact weighted traces."""

    _require_open_hypercubic(lattice, 3)
    straight = sp.sympify(straight)
    turn = sp.sympify(turn)
    traces = [sp.Integer(len(directed_edges(lattice)))]
    for length in range(1, max_order + 1):
        traces.append(
            _trace_expression(trace_type_counts(lattice, length), length, straight, turn)
        )
    coefficients: list[sp.Expr] = [sp.Integer(1)]
    for degree in range(1, max_order + 1):
        coefficient = -sum(
            traces[power] * coefficients[degree - power]
            for power in range(1, degree + 1)
        ) / degree
        coefficients.append(sp.simplify(coefficient))
    return tuple(coefficients)


def numeric_determinant_prefix(
    lattice: Lattice,
    straight,
    turn,
    max_order: int,
    *,
    dps: int = 80,
) -> tuple[mp.mpc, ...]:
    """Independently obtain a determinant prefix from dense high-precision powers."""

    _require_open_hypercubic(lattice, 3)
    edges = directed_edges(lattice)
    transitions = _transition_indices(lattice)
    with mp.workdps(dps):
        straight_value = mp.mpc(straight)
        turn_value = mp.mpc(turn)
        matrix = mp.matrix(len(edges), len(edges))
        for edge_index, next_indices in enumerate(transitions):
            for next_edge in next_indices:
                matrix[edge_index, next_edge] = (
                    straight_value
                    if edges[edge_index].direction == edges[next_edge].direction
                    else turn_value
                )
        power_matrix = mp.eye(len(edges))
        traces = [mp.mpc(len(edges))]
        for _ in range(max_order):
            power_matrix = power_matrix * matrix
            traces.append(mp.fsum(power_matrix[index, index] for index in range(len(edges))))
        coefficients = [mp.mpc(1)]
        for degree in range(1, max_order + 1):
            coefficients.append(
                -mp.fsum(
                    traces[power] * coefficients[degree - power]
                    for power in range(1, degree + 1)
                )
                / degree
            )
        return tuple(+coefficient for coefficient in coefficients)


_FULL_SYMBOLS: tuple[tuple[tuple[int, int], sp.Symbol], ...] = tuple(
    (
        pair,
        sp.Symbol(
            "u_"
            + DIRECTION_LABELS[pair[0]].replace("+", "p").replace("-", "m")
            + "_"
            + DIRECTION_LABELS[pair[1]].replace("+", "p").replace("-", "m")
        ),
    )
    for pair in ALLOWED_DIRECTION_PAIRS
)
_FULL_SYMBOL_MAP = dict(_FULL_SYMBOLS)
_FULL_SYMBOL_POSITIONS = {pair: index for index, pair in enumerate(ALLOWED_DIRECTION_PAIRS)}


@lru_cache(maxsize=None)
def _bulk_word_records(
    length: int,
) -> tuple[tuple[tuple[int, int, int], tuple[int, ...], int], ...]:
    """Aggregate cubic direction words by span and transition monomial."""

    records: dict[tuple[tuple[int, int, int], tuple[int, ...]], int] = {}

    def extend(
        first: int,
        sequence: tuple[int, ...],
        displacement: tuple[int, int, int],
        minima: tuple[int, int, int],
        maxima: tuple[int, int, int],
    ) -> None:
        if len(sequence) == length:
            if displacement != (0, 0, 0) or first == OPPOSITE_DIRECTION[sequence[-1]]:
                return
            exponents = [0] * 30
            for index, direction in enumerate(sequence):
                next_direction = sequence[(index + 1) % length]
                exponents[_FULL_SYMBOL_POSITIONS[(direction, next_direction)]] += 1
            span = tuple(maxima[axis] - minima[axis] for axis in range(3))
            key = (span, tuple(exponents))
            records[key] = records.get(key, 0) + 1
            return
        last = sequence[-1]
        for next_direction in range(6):
            if next_direction == OPPOSITE_DIRECTION[last]:
                continue
            step = DIRECTIONS_3D[next_direction]
            next_displacement = tuple(
                displacement[axis] + step[axis] for axis in range(3)
            )
            extend(
                first,
                sequence + (next_direction,),
                next_displacement,
                tuple(min(minima[axis], next_displacement[axis]) for axis in range(3)),
                tuple(max(maxima[axis], next_displacement[axis]) for axis in range(3)),
            )

    for first in range(6):
        displacement = DIRECTIONS_3D[first]
        extend(
            first,
            (first,),
            displacement,
            tuple(min(0, displacement[axis]) for axis in range(3)),
            tuple(max(0, displacement[axis]) for axis in range(3)),
        )
    return tuple(
        (span, exponents, multiplicity)
        for (span, exponents), multiplicity in sorted(records.items())
    )


def translation_invariant_trace_expression(
    shape: tuple[int, int, int], length: int, *, gauge_fix: bool = True
) -> sp.Expr:
    """Return the exact finite-box ``Tr(U**length)`` for all direction weights."""

    if any(side < 1 for side in shape):
        raise ValueError("box side lengths must be positive")
    gauge = set(DIRECTION_GAUGE_TREE) if gauge_fix else set()
    expression = sp.Integer(0)
    for span, exponents, word_multiplicity in _bulk_word_records(length):
        embeddings = word_multiplicity
        for side, width in zip(shape, span, strict=True):
            if width >= side:
                embeddings = 0
                break
            embeddings *= side - width
        if not embeddings:
            continue
        monomial = sp.Integer(embeddings)
        for position, exponent in enumerate(exponents):
            if not exponent:
                continue
            pair = ALLOWED_DIRECTION_PAIRS[position]
            if pair not in gauge:
                monomial *= _FULL_SYMBOL_MAP[pair] ** exponent
        expression += monomial
    return sp.expand(expression)


def full_weight_finite_system(
    shapes: Iterable[tuple[int, int, int]],
    orders: Iterable[int] = (4, 6, 8),
    *,
    gauge_fix: bool = True,
) -> FullWeightSystem:
    """Construct exact coefficient equations for shared 30 direction weights.

    With ``gauge_fix=True`` the five entries in ``DIRECTION_GAUGE_TREE`` are
    set to one.  This is a valid slice only on their common nonzero chart; all
    zero-entry branches remain outside the computation and must not be claimed
    as excluded.
    """

    selected_shapes = tuple(tuple(int(side) for side in shape) for shape in shapes)
    selected_orders = tuple(sorted(set(int(order) for order in orders)))
    if not selected_orders or min(selected_orders) < 1:
        raise ValueError("orders must be positive")
    equations: list[sp.Expr] = []
    labels: list[tuple[tuple[int, int, int], int]] = []
    polynomial_cache: dict[tuple[int, int, int], list[int]] = {}
    for shape in selected_shapes:
        canonical = tuple(sorted(shape))
        if canonical not in polynomial_cache:
            polynomial_cache[canonical] = even_subgraph_polynomial(
                cubic(*canonical, periodic=False)
            )
        log_p = formal_log_coefficients(polynomial_cache[canonical], max(selected_orders))
        for order in selected_orders:
            trace = translation_invariant_trace_expression(shape, order, gauge_fix=gauge_fix)
            required = -2 * order * log_p[order]
            equations.append(sp.expand(trace - required))
            labels.append((shape, order))
    gauge = set(DIRECTION_GAUGE_TREE) if gauge_fix else set()
    variables = tuple(symbol for pair, symbol in _FULL_SYMBOLS if pair not in gauge)
    return FullWeightSystem(
        variables=variables,
        all_symbols=_FULL_SYMBOLS,
        gauge_tree=DIRECTION_GAUGE_TREE if gauge_fix else (),
        labels=tuple(labels),
        equations=tuple(equations),
        pre_gauge_unknown_count=30,
        post_gauge_unknown_count=len(variables),
    )
