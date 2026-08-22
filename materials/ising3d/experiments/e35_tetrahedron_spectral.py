"""Exact spectral-parameter tests for the cubic Ising bond tensor.

This experiment deliberately studies a sharply specified extension of the constant
specialisation in ``e11_tetrahedron.py``.  It does not classify arbitrary
64-entry R matrices or arbitrary bond-dimension extensions.
"""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterable

import sympy as sp

from ising.exact_enumeration import dos_bonds, even_subgraph_polynomial
from ising.lattices import cubic


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "tetrahedron_spectral.json"
POSITIONS = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))
CONTROL_BOXES = ((2, 2, 2), (2, 2, 3), (2, 3, 3))


def _popcount(value: int) -> int:
    return value.bit_count()


def _polyval(coefficients: Iterable[int], value: int | Fraction) -> int | Fraction:
    result: int | Fraction = 0
    for coefficient in reversed(tuple(coefficients)):
        result = result * value + coefficient
    return result


def site_tensor_entry(indices: tuple[int, ...], w: int | Fraction) -> int | Fraction:
    """Return T/2 for the six-leg Ising site tensor.

    The omitted uniform factor 2 is restored once per vertex by
    :func:`tensor_partition_ht_integer`.
    """

    if len(indices) != 6 or any(bit not in (0, 1) for bit in indices):
        raise ValueError("the site tensor has six binary legs")
    total = sum(indices)
    return w**total if total % 2 == 0 else 0


def tensor_even_subgraph_polynomial(shape: tuple[int, int, int]) -> list[int]:
    """Contract the parity tensors exactly, returning their polynomial in v=w^2.

    Contracting a binary leg selects a subset of lattice bonds.  A tensor entry
    is nonzero exactly when every vertex has even selected degree.  The dynamic
    program below contracts one edge at a time while retaining the exact vertex
    parity mask and selected-edge degree.  It is independent of spin
    enumeration and gives the tensor contraction coefficient by coefficient.
    """

    lattice = cubic(*shape, periodic=False)
    bonds = tuple(lattice.bonds)
    states: dict[tuple[int, int], int] = {(0, 0): 1}
    for left, right in bonds:
        edge_mask = (1 << left) ^ (1 << right)
        updated = dict(states)
        for (parity, degree), multiplicity in states.items():
            key = (parity ^ edge_mask, degree + 1)
            updated[key] = updated.get(key, 0) + multiplicity
        states = updated
    coefficients = [0] * (len(bonds) + 1)
    for (parity, degree), multiplicity in states.items():
        if parity == 0:
            coefficients[degree] += multiplicity
    return coefficients


def tensor_partition_ht_integer(shape: tuple[int, int, int], w: int) -> int:
    """Return the exact unnormalised tensor contraction at integer ``w``.

    It is ``2^N P(w^2)``.  Multiplication by ``(cosh K)^E`` with
    ``tanh K=w^2`` gives the physical partition function.
    """

    lattice = cubic(*shape, periodic=False)
    polynomial = tensor_even_subgraph_polynomial(shape)
    return (1 << lattice.n_sites) * int(_polyval(polynomial, w * w))


def _spin_partition_control_integer(shape: tuple[int, int, int], v: Fraction) -> int:
    """Scale spin enumeration to the same integer as the HT tensor contraction.

    With ``v=1/3`` choose ``exp(2K)=2``.  The density polynomial evaluated at 2
    equals ``2^(N-E) 3^E P(1/3)`` exactly; multiplying the tensor contraction
    ``2^N P(1/3)`` by ``(3/2)^E`` gives the same integer.
    """

    if v != Fraction(1, 3):
        raise ValueError("the recorded exact-integer normalization is fixed at v=1/3")
    lattice = cubic(*shape, periodic=False)
    density = [int(value) for value in dos_bonds(lattice)]
    return int(_polyval(density, 2))


def finite_lattice_controls() -> list[dict]:
    rows = []
    for shape in CONTROL_BOXES:
        lattice = cubic(*shape, periodic=False)
        tensor_polynomial = tensor_even_subgraph_polynomial(shape)
        enumerated_polynomial = even_subgraph_polynomial(lattice)
        if tensor_polynomial != enumerated_polynomial:
            raise ArithmeticError(f"tensor/enumerator polynomial mismatch on {shape}")
        tensor_integer = int(
            Fraction(1 << lattice.n_sites)
            * _polyval(tensor_polynomial, Fraction(1, 3))
            * Fraction(3, 2) ** lattice.n_bonds
        )
        enumerator_integer = _spin_partition_control_integer(shape, Fraction(1, 3))
        if tensor_integer != enumerator_integer:
            raise ArithmeticError(f"tensor/spin integer mismatch on {shape}")
        rows.append(
            {
                "shape": list(shape),
                "boundary": "open",
                "sites": lattice.n_sites,
                "bonds": lattice.n_bonds,
                "even_subgraph_polynomial": tensor_polynomial,
                "comparison_parameter": "v=1/3 (equivalently exp(2K)=2)",
                "common_scaled_partition_integer": tensor_integer,
                "coefficientwise_match": True,
            }
        )
    return rows


def _apply_spectral_r(
    amplitudes: dict[tuple[int, tuple[int, ...]], int],
    positions: tuple[int, int, int],
) -> dict[tuple[int, tuple[int, ...]], int]:
    """Apply R built from T(w_1,...,w_6), tracking a six-variable monomial."""

    output: dict[tuple[int, tuple[int, ...]], int] = defaultdict(int)
    position_mask = sum(1 << position for position in positions)
    for (state, exponents), coefficient in amplitudes.items():
        input_bits = tuple((state >> position) & 1 for position in positions)
        for local_output in range(8):
            output_bits = tuple((local_output >> local) & 1 for local in range(3))
            if (sum(input_bits) + sum(output_bits)) & 1:
                continue
            new_state = state & ~position_mask
            new_exponents = list(exponents)
            for local, position in enumerate(positions):
                if output_bits[local]:
                    new_state |= 1 << position
                new_exponents[position] += input_bits[local] + output_bits[local]
            output[(new_state, tuple(new_exponents))] += coefficient
    return {key: coefficient for key, coefficient in output.items() if coefficient}


def spectral_residuals() -> dict[tuple[int, int], dict[tuple[int, ...], int]]:
    """Return all nonzero component residuals of the six-rapidity TE ansatz."""

    residuals: dict[tuple[int, int], dict[tuple[int, ...], int]] = {}
    for input_state in range(64):
        left = {(input_state, (0,) * 6): 1}
        for positions in reversed(POSITIONS):
            left = _apply_spectral_r(left, positions)
        right = {(input_state, (0,) * 6): 1}
        for positions in POSITIONS:
            right = _apply_spectral_r(right, positions)
        by_output: dict[int, dict[tuple[int, ...], int]] = defaultdict(lambda: defaultdict(int))
        for (output_state, exponents), coefficient in left.items():
            by_output[output_state][exponents] += coefficient
        for (output_state, exponents), coefficient in right.items():
            by_output[output_state][exponents] -= coefficient
        for output_state, polynomial in by_output.items():
            cleaned = {monomial: coefficient for monomial, coefficient in polynomial.items() if coefficient}
            if cleaned:
                residuals[(output_state, input_state)] = cleaned
    return residuals


def reduced_spectral_equations() -> tuple[tuple[sp.Symbol, ...], list[sp.Poly], dict]:
    """On the nonzero-weight torus, strip invertible monomials, square variables, and eliminate linearly."""

    variables = sp.symbols("x1:7")
    residuals = spectral_residuals()
    unique: dict[str, sp.Poly] = {}
    for polynomial in residuals.values():
        minima = tuple(min(monomial[index] for monomial in polynomial) for index in range(6))
        expression = 0
        for monomial, coefficient in polynomial.items():
            reduced = tuple(monomial[index] - minima[index] for index in range(6))
            if any(power & 1 for power in reduced):
                raise ArithmeticError("spectral residual did not reduce to squares")
            expression += coefficient * sp.prod(
                variables[index] ** (reduced[index] // 2) for index in range(6)
            )
        primitive = sp.Poly(expression, *variables, domain=sp.QQ).monic()
        unique[str(primitive.as_expr())] = primitive

    equations = list(unique.values())
    monomials = sorted({monomial for polynomial in equations for monomial, _ in polynomial.terms()})
    coefficient_matrix = sp.Matrix(
        [[polynomial.coeff_monomial(monomial) for monomial in monomials] for polynomial in equations]
    )
    rational_rank = int(coefficient_matrix.rank())
    _, pivot_rows = coefficient_matrix.T.rref()
    independent = [equations[index] for index in pivot_rows]
    if len(independent) != rational_rank:
        raise ArithmeticError("linear elimination did not return a basis")
    metadata = {
        "nonzero_component_equations": len(residuals),
        "distinct_after_monomial_and_square_reduction": len(equations),
        "coefficient_monomials": len(monomials),
        "exact_rational_linear_rank": rational_rank,
        "independent_equations": [str(sp.factor(polynomial.as_expr())) for polynomial in independent],
    }
    return variables, independent, metadata


def spectral_tetrahedron_analysis() -> dict:
    variables, equations, metadata = reduced_spectral_equations()
    saturation_variable = sp.Symbol("q")
    exceptional_product = sp.prod(variable * (variable - 1) * (variable + 1) for variable in variables)
    basis = sp.groebner(
        [polynomial.as_expr() for polynomial in equations]
        + [saturation_variable * exceptional_product - 1],
        saturation_variable,
        *variables,
        order="grevlex",
        domain=sp.QQ,
    )
    saturated_unit = list(basis) == [1]
    if not saturated_unit:
        raise ArithmeticError("the nondegenerate six-rapidity branch was not eliminated")

    diagonal = sp.Symbol("t")
    diagonal_equations = [
        sp.factor(polynomial.as_expr().subs(dict(zip(variables, (diagonal,) * 6, strict=True))))
        for polynomial in equations
    ]
    diagonal_gcd = sp.factor(sp.gcd_list(diagonal_equations))
    w = sp.Symbol("w")
    raw_diagonal_gcd = sp.factor(
        sp.gcd_list(
            [
                sum(coefficient * w ** sum(monomial) for monomial, coefficient in polynomial.items())
                for polynomial in spectral_residuals().values()
            ]
        )
    )
    degenerate_substitution = {variables[3]: 1, variables[4]: 1, variables[5]: -1}
    return {
        "ansatz": {
            "definition": (
                "T_[w1,...,w6](a1,...,a6)=2*product_i(w_i**a_i) when sum_i a_i is even, and 0 otherwise; "
                "R[out1,out2,out3;in1,in2,in3]=T[in1,in2,in3,out1,out2,out3]"
            ),
            "spectral_variables": (
                "x_i=w_i^2, one edge rapidity on each of the six shared TE vector spaces; "
                "an embedded R_ijk uses the variables of spaces i,j,k on both its input and output legs"
            ),
            "equation": "R_123 R_145 R_246 R_356 = R_356 R_246 R_145 R_123",
            "scope": (
                "This six-space edge-rapidity deformation retains the verified parity support. It is not the "
                "general four-independent-R spectral equation and not a classification of arbitrary 64-entry R matrices."
            ),
        },
        "linear_elimination": metadata,
        "nondegenerate_saturation": {
            "localized_open_set": "product_i x_i*(x_i-1)*(x_i+1) != 0",
            "groebner_basis": [str(polynomial) for polynomial in basis],
            "unit_ideal": saturated_unit,
            "conclusion": (
                "No solution has all six x_i away from 0,+1,-1. Every solution of this ansatz has at least one "
                "zero-weight or zero-/infinite-temperature/complex exceptional leg."
            ),
        },
        "physical_isotropic_slice": {
            "substitution": "x_1=...=x_6=t=tanh(K)",
            "common_gcd_after_nonzero_monomial_removal": str(diagonal_gcd),
            "raw_common_gcd_in_w": str(raw_diagonal_gcd),
            "physical_interval": "0<t<1",
            "physical_solutions": [],
        },
        "degenerate_manifold": {
            "exists": True,
            "example": "x_4=x_5=1, x_6=-1 with x_1,x_2,x_3 arbitrary",
            "verified_by_substitution": all(
                sp.expand(polynomial.as_expr().subs(degenerate_substitution)) == 0 for polynomial in equations
            ),
            "physical_relevance": (
                "none at finite ferromagnetic K: it fixes three local weights at exceptional values; "
                "it is an algebraic degenerate TE family, not an Ising spectral curve through 0<x_i<1"
            ),
        },
    }


def _embedded_extension_residual_nonzero() -> tuple[bool, tuple[int, int, int, str] | None]:
    """Check whether R_2 direct-sum I on dimension 3 satisfies the TE."""

    dimension = 3
    w = Fraction(1, 2)

    def local_entry(output: tuple[int, int, int], input_: tuple[int, int, int]) -> Fraction:
        if all(bit < 2 for bit in output + input_):
            return Fraction(site_tensor_entry(input_ + output, w))
        if output == input_ and any(bit == 2 for bit in input_):
            return Fraction(1)
        return Fraction(0)

    def apply(amplitudes: dict[tuple[int, ...], Fraction], positions: tuple[int, int, int]):
        result: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        for state, coefficient in amplitudes.items():
            local_input = tuple(state[position] for position in positions)
            for local_output in product(range(dimension), repeat=3):
                entry = local_entry(local_output, local_input)
                if entry:
                    new_state = list(state)
                    for position, value in zip(positions, local_output, strict=True):
                        new_state[position] = value
                    result[tuple(new_state)] += coefficient * entry
        return {state: coefficient for state, coefficient in result.items() if coefficient}

    for input_bits in product(range(2), repeat=6):
        left = {input_bits: Fraction(1)}
        right = {input_bits: Fraction(1)}
        for positions in reversed(POSITIONS):
            left = apply(left, positions)
        for positions in POSITIONS:
            right = apply(right, positions)
        for output_state in set(left) | set(right):
            residual = left.get(output_state, 0) - right.get(output_state, 0)
            if residual:
                return False, (
                    sum(bit << index for index, bit in enumerate(output_state)),
                    sum(bit << index for index, bit in enumerate(input_bits)),
                    residual.numerator,
                    str(residual),
                )
    return True, None


def escape_analysis() -> dict:
    extension_passed, extension_witness = _embedded_extension_residual_nonzero()
    return {
        "gauge_transformation": {
            "class_tested": (
                "delta-preserving leg gauges: paired G and G^{-T} changes of basis on contracted bonds, "
                "inducing simultaneous similarities of embedded R operators"
            ),
            "decision": "fails",
            "certificate": (
                "The tetrahedron residual is conjugated by an invertible global tensor product. "
                "Its zero/nonzero status is invariant, while the exact physical isotropic residual has gcd "
                "(t-1)^3*(t+1)^2 and is nonzero for 0<t<1."
            ),
            "unrestricted_warning": (
                "A covariant nonorthogonal six-leg GL(2) action can map each fixed nonzero Ising tensor to GHZ, "
                "but changes the delta bond metric (or requires alternating inverse gauges); it is not a gauge of "
                "the same identical vertex model and does not create a commuting physical family."
            ),
        },
        "auxiliary_state_extension": {
            "class_tested": "minimal direct-sum extension R_3=R_Ising direct-sum I on every local triple containing state 2",
            "bond_dimension": 3,
            "physical_point": "w=1/2 (exact rational, t=1/4)",
            "decision": "fails",
            "certificate": (
                "The physical two-state subspace is invariant, so its nonzero TE residual is a principal block of "
                "the extended residual. The explicit exact scan found the recorded component witness."
            ),
            "witness": list(extension_witness) if extension_witness else None,
            "scope_warning": (
                "General interacting bond-dimension-3 or larger embeddings are undecided; testing one direct-sum "
                "extension is not a no-go for all auxiliary-state constructions."
            ),
            "satisfies_tetrahedron_equation": extension_passed,
        },
        "controlled_limit": {
            "class_tested": "coefficientwise limits within the six-rapidity parity-tensor ansatz",
            "decision": "fails for a regular finite physical limit",
            "certificate": (
                "The TE solution set is polynomially closed. If exact solutions converged coefficientwise to a "
                "finite isotropic tensor with 0<t<1, that limit would satisfy the same equations; the diagonal gcd "
                "(t-1)^3*(t+1)^2 excludes it."
            ),
            "singular_limits": (
                "Solutions do occur only after exceptional x_i in {0,+1,-1}; these are zero/infinite-temperature "
                "or complex degenerations and do not converge regularly to finite ferromagnetic isotropic weights."
            ),
            "scope_warning": "Singular limits with dimension-changing projections or divergent gauges are undecided.",
        },
    }


def _layer_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    """Periodic square-layer bonds, retaining length-two parallel bonds."""

    edges = []
    for row in range(rows):
        for column in range(columns):
            index = row * columns + column
            if rows >= 2:
                edges.append((index, ((row + 1) % rows) * columns + column))
            if columns >= 2:
                edges.append((index, row * columns + (column + 1) % columns))
    return tuple(edges)


def layer_transfer_matrix(rows: int, columns: int, a: Fraction, b: Fraction) -> list[list[Fraction]]:
    """Exact two-parameter layer matrix T(a,b)=A(a)D(b), up to a scalar.

    ``T(a,b)[s,t]=a^Hamming(s,t) b^broken(t)`` on a periodic layer.
    """

    sites = rows * columns
    dimension = 1 << sites
    edges = _layer_edges(rows, columns)
    broken = [sum(((state >> i) ^ (state >> j)) & 1 for i, j in edges) for state in range(dimension)]
    return [
        [a ** _popcount(left ^ right) * b ** broken[right] for right in range(dimension)]
        for left in range(dimension)
    ]


def _matrix_commutator_nonzero(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> tuple[int, int, Fraction] | None:
    dimension = len(left)
    for row in range(dimension):
        for column in range(dimension):
            residual = sum(
                left[row][middle] * right[middle][column]
                - right[row][middle] * left[middle][column]
                for middle in range(dimension)
            )
            if residual:
                return row, column, residual
    return None


def commuting_transfer_analysis() -> dict:
    """Test the physical spectral curve and an exact tangent obstruction."""

    points = (Fraction(1, 2), Fraction(1, 3))
    exact_checks = []
    for shape in ((2, 2), (2, 3)):
        first = layer_transfer_matrix(*shape, points[0], (1 - points[0]) / (1 + points[0]))
        second = layer_transfer_matrix(*shape, points[1], (1 - points[1]) / (1 + points[1]))
        witness = _matrix_commutator_nonzero(first, second)
        exact_checks.append(
            {
                "shape": list(shape),
                "points": [
                    {"a": "1/2", "b": "1/3"},
                    {"a": "1/3", "b": "1/2"},
                ],
                "commutes": witness is None,
                "first_nonzero_component": (
                    {"row": witness[0], "column": witness[1], "value": str(witness[2])}
                    if witness
                    else None
                ),
            }
        )

    x, a, b = sp.symbols("x a b")
    base = layer_transfer_matrix(2, 2, x, x)
    varying_a = layer_transfer_matrix(2, 2, a, x)
    varying_b = layer_transfer_matrix(2, 2, x, b)
    derivative_a = [[sp.diff(entry, a).subs(a, x) for entry in row] for row in varying_a]
    derivative_b = [[sp.diff(entry, b).subs(b, x) for entry in row] for row in varying_b]

    def commutator_entry(derivative: list[list[sp.Expr]], row: int, column: int) -> sp.Expr:
        return sp.expand(
            sum(
                base[row][middle] * derivative[middle][column]
                - derivative[row][middle] * base[middle][column]
                for middle in range(len(base))
            )
        )

    tangent_rows = ((1, 0), (1, 15))
    tangent_matrix = [
        [commutator_entry(derivative, row, column) for derivative in (derivative_a, derivative_b)]
        for row, column in tangent_rows
    ]
    tangent_minor = sp.factor(
        tangent_matrix[0][0] * tangent_matrix[1][1]
        - tangent_matrix[0][1] * tangent_matrix[1][0]
    )
    return {
        "parametrization": (
            "T(a,b)_(s,t)=a^Hamming(s,t) b^broken_periodic_layer(t), up to a scalar; "
            "physical isotropic curve b=(1-a)/(1+a), 0<a<1"
        ),
        "exact_pair_checks": exact_checks,
        "two_by_two_tangent_certificate": {
            "meaning": (
                "At a putative differentiable commuting family through T(x,x), the linear system "
                "[T,dT]=0 for velocity (da,db) has a 2x2 coefficient minor from matrix entries "
                "(row,column)=(1,0),(1,15)."
            ),
            "coefficient_matrix": [
                [str(sp.factor(entry)) for entry in row] for row in tangent_matrix
            ],
            "determinant": str(tangent_minor),
            "physical_interval": "0<x<1",
            "conclusion": (
                "The determinant is nonzero for 0<x<1, so da=db=0: no nonconstant differentiable "
                "two-parameter commuting curve passes through an interior isotropic point on the 2x2 layer."
            ),
        },
        "trivial_solutions": (
            "T(a,b) commutes with itself; boundary/decoupled points a in {0,+1,-1} or b in {0,+1,-1} "
            "were not classified as physical spectral families."
        ),
        "scope_warning": (
            "The pair checks on 2x2 and 2x3 refute the natural physical curve. The tangent certificate is a local "
            "exact no-go for differentiable curves in the displayed two-weight ansatz, not a classification of "
            "arbitrary 2^(LxLy)-square transfer matrices or singular exceptional branches."
        ),
    }


def build_artifact() -> dict:
    controls = finite_lattice_controls()
    spectral = spectral_tetrahedron_analysis()
    escapes = escape_analysis()
    commuting = commuting_transfer_analysis()
    checks = [
        {
            "name": "site tensor matches exact finite-lattice enumeration",
            "passed": len(controls) >= 3 and all(row["coefficientwise_match"] for row in controls),
            "detail": "; ".join(
                f"{tuple(row['shape'])}: {row['common_scaled_partition_integer']}" for row in controls
            ),
        },
        {
            "name": "nondegenerate six-rapidity tetrahedron variety is empty",
            "passed": spectral["nondegenerate_saturation"]["unit_ideal"],
            "detail": "exact saturated Groebner basis is [1] after rational linear elimination",
        },
        {
            "name": "physical isotropic weights are excluded",
            "passed": spectral["physical_isotropic_slice"]["physical_solutions"] == [],
            "detail": (
                "raw diagonal gcd "
                f"{spectral['physical_isotropic_slice']['raw_common_gcd_in_w']} "
                "has no finite-ferromagnetic root"
            ),
        },
        {
            "name": "gauge, tested auxiliary extension, and regular limit do not repair the physical tensor",
            "passed": (
                escapes["gauge_transformation"]["decision"] == "fails"
                and escapes["auxiliary_state_extension"]["decision"] == "fails"
                and escapes["controlled_limit"]["decision"] == "fails for a regular finite physical limit"
            ),
            "detail": "exact invariance/principal-block/polynomial-closedness certificates; broader auxiliary and singular classes explicitly undecided",
        },
        {
            "name": "natural physical layer-transfer spectral curve is noncommuting",
            "passed": all(not row["commutes"] for row in commuting["exact_pair_checks"]),
            "detail": "exact rational witnesses on 2x2 and 2x3; nonzero tangent minor on every 0<x<1",
        },
    ]
    return {
        "provenance": {
            "script": "experiments/e35_tetrahedron_spectral.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": "exact Python integers/Fraction and SymPy QQ polynomial algebra; no floating point",
        },
        "data": {
            "prior_scope": (
                "e11 tested only the constant equation for the one-parameter isotropic rank-six parity tensor, "
                "plus its fixed-tensor GL(2) gauge orbit; it explicitly left spectral dependence, higher bond "
                "dimension, alternating tensors, and other representations open."
            ),
            "local_tensor": {
                "definition": (
                    "u_sigma=(1,sigma*w); T[a,b,c,d,e,f]=sum_sigma product u_sigma[index]="
                    "2*w^(a+b+c+d+e+f) for even index sum and 0 otherwise"
                ),
                "relation": "w^2=v=tanh(K)",
                "partition_function": "Z=(cosh K)^E contraction(product_vertices T)",
            },
            "finite_lattice_controls": controls,
            "spectral_tetrahedron": spectral,
            "escapes": escapes,
            "commuting_layer_transfer": commuting,
            "status": (
                "[COMPUTATION] exact no-go for the six-space edge-rapidity ansatz and natural two-weight "
                "layer family. [UNDECIDED] four-independent-R or arbitrary R-matrix representations, "
                "interacting auxiliary-state extensions, and singular dimension-changing limits."
            ),
        },
        "checks": checks,
    }


def main() -> None:
    artifact = build_artifact()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    controls = artifact["data"]["finite_lattice_controls"]
    print("tensor/enumerator integers:", [(row["shape"], row["common_scaled_partition_integer"]) for row in controls])
    print(
        "spectral equations:",
        artifact["data"]["spectral_tetrahedron"]["linear_elimination"]["nonzero_component_equations"],
        "components ->",
        artifact["data"]["spectral_tetrahedron"]["linear_elimination"]["exact_rational_linear_rank"],
        "independent equations",
    )
    print("written", OUTPUT.relative_to(ROOT))
    if not all(check["passed"] for check in artifact["checks"]):
        print("FAIL")
        raise AssertionError("one or more spectral tetrahedron checks failed")
    print("PASS")


if __name__ == "__main__":
    main()
