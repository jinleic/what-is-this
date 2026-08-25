#!/usr/bin/env python3
"""L29 reciprocal-lift classification and dyadic-compressed pullback.

Everything labelled PROVED is an exact identity, valuation argument, or finite
2-adic congruence exhaustion.  Bounded member rows are deliberately omitted.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

from h10q import FactorBudget, PrimalityBound, legendre, ramified
from l25_scaled_coupling import vp
from l26_reciprocal_tie import bridge_data, reciprocal_polynomials, reciprocal_tie
from l28_trace_field import resolvent_coefficients

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l29_reciprocal_frontier.jsonl"
REPORT = Path("/tmp/l29_reciprocal_frontier.md")


def power_of_two(exponent: int) -> F:
    return F(2) ** exponent


def unit_mod(value: F, modulus: int) -> int:
    assert value.denominator % 2
    return value.numerator * pow(value.denominator, -1, modulus) % modulus




def translated(polynomial: list[F], shift: F) -> list[F]:
    answer = [F(0)] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for out_degree in range(degree + 1):
            answer[out_degree] += (
                coefficient
                * math.comb(degree, out_degree)
                * shift ** (degree - out_degree)
            )
    return answer


def normalized_shifted_octic(a: F, Z: F) -> list[F]:
    """Return monic P(1+x), where P(b)=b^4*T(b+b^-1)."""
    p, q, r, _ = resolvent_coefficients(a, Z)
    answer = [
        r,
        4 * r,
        6 * r + q,
        4 * r + 3 * q,
        r + 3 * q + p,
        q + 2 * p,
        p,
        F(0),
        F(1),
    ]
    A = 1 + 4 * a * a
    P, _ = reciprocal_polynomials(a, Z, -4 * a**4 / A)
    direct = translated(P, F(1))
    assert answer == [coefficient / direct[-1] for coefficient in direct]
    return answer


def convolution(left: list[int], right: list[int], modulus: int) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            answer[left_degree + right_degree] = (
                answer[left_degree + right_degree]
                + left_value * right_value
            ) % modulus
    return answer


def eisenstein_quartic_factor_search(
    target: list[int], exponent: int
) -> tuple[int, list[list[int]] | None]:
    """Exhaust all possible first Eisenstein quartics modulo 2^exponent.

    The high four product coefficients determine the second monic quartic, so
    the search has 2^(4*(exponent-1)-1) rather than eight free coefficients.
    """
    modulus = 1 << exponent
    target = [value % modulus for value in target]
    even_values = range(0, modulus, 2)
    constant_values = range(2, modulus, 4)
    candidates = 0
    for f0 in constant_values:
        for f1 in even_values:
            for f2 in even_values:
                for f3 in even_values:
                    candidates += 1
                    g3 = (target[7] - f3) % modulus
                    g2 = (target[6] - f2 - f3 * g3) % modulus
                    g1 = (target[5] - f1 - f2 * g3 - f3 * g2) % modulus
                    g0 = (
                        target[4]
                        - f0
                        - f1 * g3
                        - f2 * g2
                        - f3 * g1
                    ) % modulus
                    if any(value % 2 for value in (g1, g2, g3)):
                        continue
                    if g0 % 4 != 2:
                        continue
                    left = [f0, f1, f2, f3, 1]
                    right = [g0, g1, g2, g3, 1]
                    if convolution(left, right, modulus) == target:
                        return candidates, [left, right]
    return candidates, None


def congruence_obstruction() -> dict:
    expected = {
        "t<=-4": [20, 16, 56, 16, 52, 0, 32, 0, 1],
        "t=-2": [52, 16, 56, 16, 20, 0, 32, 0, 1],
        "t=0": [52, 16, 48, 56, 60, 56, 32, 0, 1],
    }
    rows = 0
    for t in (-8, -6, -4, -2, 0):
        label = "t<=-4" if t <= -4 else f"t={t}"
        for a_residue in range(1, 64, 2):
            for z_residue in range(1, 64, 2):
                Z = power_of_two(t) * F(z_residue) ** 2
                polynomial = normalized_shifted_octic(F(a_residue), Z)
                assert all(coefficient.denominator % 2 for coefficient in polynomial)
                assert [unit_mod(coefficient, 64) for coefficient in polynomial] == expected[label]
                rows += 1

    t_zero_candidates, t_zero_factor = eisenstein_quartic_factor_search(
        [value % 16 for value in expected["t=0"]], 4
    )
    negative_candidates, negative_factor = eisenstein_quartic_factor_search(
        expected["t<=-4"], 6
    )
    assert t_zero_candidates == 2048 and t_zero_factor is None
    assert negative_candidates == 524288 and negative_factor is None

    return {
        "type": "reciprocal-lift-congruence-obstruction",
        "label": "PROVED finite 2-adic congruence lemma",
        "unit_rows": rows,
        "coefficient_order": "ascending coefficients of monic P(1+x)",
        "vectors_mod_64": expected,
        "factor_shape": (
            "the slope -1/4 forces any proper factorization to be 4+4; "
            "both monic quartics are Eisenstein"
        ),
        "t=0": {
            "modulus": 16,
            "candidate_first_factors": t_zero_candidates,
            "survivors": 0,
        },
        "t<=-4": {
            "modulus": 64,
            "candidate_first_factors": negative_candidates,
            "survivors": 0,
        },
    }


def sparse_add(
    left: dict[tuple[int, int, int], F],
    right: dict[tuple[int, int, int], F],
) -> dict[tuple[int, int, int], F]:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, F(0)) + coefficient
        if not answer[monomial]:
            del answer[monomial]
    return answer


def sparse_scale(
    polynomial: dict[tuple[int, int, int], F], scalar: F
) -> dict[tuple[int, int, int], F]:
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in polynomial.items()
        if coefficient * scalar
    }


def sparse_multiply(
    left: dict[tuple[int, int, int], F],
    right: dict[tuple[int, int, int], F],
) -> dict[tuple[int, int, int], F]:
    answer: dict[tuple[int, int, int], F] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_monomial[index] + right_monomial[index]
                for index in range(3)
            )
            answer[monomial] = (
                answer.get(monomial, F(0))
                + left_coefficient * right_coefficient
            )
    return {
        monomial: coefficient
        for monomial, coefficient in answer.items()
        if coefficient
    }


def sparse_power(
    polynomial: dict[tuple[int, int, int], F], exponent: int
) -> dict[tuple[int, int, int], F]:
    answer = {(0, 0, 0): F(1)}
    base = polynomial
    while exponent:
        if exponent & 1:
            answer = sparse_multiply(answer, base)
        base = sparse_multiply(base, base)
        exponent //= 2
    return answer


def matrix_zero(size: int) -> list[list[dict[tuple[int, int, int], F]]]:
    return [[{} for _ in range(size)] for _ in range(size)]


def matrix_identity(size: int) -> list[list[dict[tuple[int, int, int], F]]]:
    answer = matrix_zero(size)
    for index in range(size):
        answer[index][index] = {(0, 0, 0): F(1)}
    return answer


def matrix_add(
    left: list[list[dict[tuple[int, int, int], F]]],
    right: list[list[dict[tuple[int, int, int], F]]],
) -> list[list[dict[tuple[int, int, int], F]]]:
    return [
        [
            sparse_add(left[row][column], right[row][column])
            for column in range(len(left))
        ]
        for row in range(len(left))
    ]


def matrix_scale(
    matrix: list[list[dict[tuple[int, int, int], F]]], scalar: F
) -> list[list[dict[tuple[int, int, int], F]]]:
    return [
        [sparse_scale(entry, scalar) for entry in row]
        for row in matrix
    ]


def matrix_multiply(
    left: list[list[dict[tuple[int, int, int], F]]],
    right: list[list[dict[tuple[int, int, int], F]]],
) -> list[list[dict[tuple[int, int, int], F]]]:
    size = len(left)
    answer = matrix_zero(size)
    for row in range(size):
        for middle in range(size):
            if not left[row][middle]:
                continue
            for column in range(size):
                if right[middle][column]:
                    answer[row][column] = sparse_add(
                        answer[row][column],
                        sparse_multiply(
                            left[row][middle], right[middle][column]
                        ),
                    )
    return answer


def matrix_add_diagonal(
    matrix: list[list[dict[tuple[int, int, int], F]]],
    polynomial: dict[tuple[int, int, int], F],
) -> list[list[dict[tuple[int, int, int], F]]]:
    answer = [[dict(entry) for entry in row] for row in matrix]
    for index in range(len(answer)):
        answer[index][index] = sparse_add(
            answer[index][index], polynomial
        )
    return answer


def substitute_resultant_parameters(
    polynomial: dict[tuple[int, int, int], F]
) -> dict[tuple[int, int, int], F]:
    variables = [
        {(1, 0, 0): F(32)},
        {(0, 1, 0): F(128)},
        {(0, 0, 0): F(52), (0, 0, 1): F(64)},
    ]
    answer: dict[tuple[int, int, int], F] = {}
    for monomial, coefficient in polynomial.items():
        term = {(0, 0, 0): coefficient}
        for index, exponent in enumerate(monomial):
            term = sparse_multiply(
                term, sparse_power(variables[index], exponent)
            )
        answer = sparse_add(answer, term)
    return answer


def exceptional_split_stratum() -> dict:
    # Build multiplication by x in Q[p,q,r][x]/(S), then use the
    # Faddeev-LeVerrier identity for the characteristic polynomial of
    # phi(x)=x^4+6*x^3+6.  This is Res_x(S,Y-phi(x)).
    size = 8
    one = {(0, 0, 0): F(1)}
    p = {(1, 0, 0): F(1)}
    q = {(0, 1, 0): F(1)}
    r = {(0, 0, 1): F(1)}
    coefficients = [
        r,
        sparse_scale(r, 4),
        sparse_add(sparse_scale(r, 6), q),
        sparse_add(sparse_scale(r, 4), sparse_scale(q, 3)),
        sparse_add(sparse_add(r, sparse_scale(q, 3)), p),
        sparse_add(q, sparse_scale(p, 2)),
        p,
        {},
    ]
    multiplication_x = matrix_zero(size)
    for column in range(size - 1):
        multiplication_x[column + 1][column] = one
    for row, coefficient in enumerate(coefficients):
        multiplication_x[row][-1] = sparse_scale(coefficient, -1)

    x2 = matrix_multiply(multiplication_x, multiplication_x)
    x3 = matrix_multiply(x2, multiplication_x)
    x4 = matrix_multiply(x3, multiplication_x)
    multiplication_phi = matrix_add(
        matrix_add(x4, matrix_scale(x3, 6)),
        matrix_scale(matrix_identity(size), 6),
    )

    working = matrix_identity(size)
    characteristic = [one]
    for degree in range(1, size + 1):
        product = matrix_multiply(multiplication_phi, working)
        trace: dict[tuple[int, int, int], F] = {}
        for index in range(size):
            trace = sparse_add(trace, product[index][index])
        coefficient = sparse_scale(trace, F(-1, degree))
        characteristic.append(coefficient)
        working = matrix_add_diagonal(product, coefficient)

    # At t=-2: p=32*P, q=128*Q, r=52+64*R with P,Q,R in Z_2.
    substituted = [
        substitute_resultant_parameters(coefficient)
        for coefficient in characteristic
    ]
    lower_bounds_descending = [
        min(vp(coefficient, 2) for coefficient in polynomial.values())
        for polynomial in substituted
    ]
    assert lower_bounds_descending == [0, 5, 6, 11, 10, 15, 17, 19, 21]
    exact_indices = {0, 1, 2, 4, 5, 7, 8}
    for index in exact_indices:
        bound = lower_bounds_descending[index]
        odd_terms = {
            monomial
            for monomial, coefficient in substituted[index].items()
            if vp(coefficient, 2) == bound
        }
        assert odd_terms == {(0, 0, 0)}

    rows = 0
    for a_residue in range(1, 64, 2):
        for z_residue in range(1, 64, 2):
            p_value, q_value, r_value, _ = resolvent_coefficients(
                F(a_residue), F(z_residue * z_residue, 4)
            )
            assert vp(p_value, 2) == 5
            assert vp(q_value, 2) == 7
            assert vp(r_value - 52, 2) >= 6
            rows += 1

    return {
        "type": "reciprocal-lift-exceptional-split-stratum",
        "label": "PROVED by an ordinary resultant Newton polygon",
        "stratum": "v2(Z)=-2",
        "resultant": "R(Y)=Res_x(S(x),Y-(x^4+6*x^3+6))",
        "coefficient_valuation_bounds_ascending": [
            21, 19, 17, 15, 10, 11, 6, 5, 0
        ],
        "exact_vertices": [[0, 21], [4, 10], [8, 0]],
        "lower_sides": [
            {"horizontal_length": 4, "slope": "-11/4"},
            {"horizontal_length": 4, "slope": "-5/2"},
        ],
        "argument": (
            "if S were irreducible, valuation invariance under all Q_2 "
            "embeddings would give one valuation for every conjugate "
            "phi(alpha), hence a one-slope resultant; the two slopes force "
            "reducibility, and the first Newton polygon forces 4+4"
        ),
        "symbolic_resultant_terms_descending": [
            len(coefficient) for coefficient in characteristic
        ],
        "unit_rows": rows,
    }


def reciprocal_lift_classification() -> dict:
    # The positive stratum uses v2(4/w)=t+3/2>2 in K and strong Hensel at 1.
    for t in (2, 3, 4, 5, 7):
        assert F(t) + F(3, 2) > 2
    return {
        "type": "uniform-reciprocal-lift-classification",
        "label": "PROVED over Q_2",
        "hypotheses": "v2(a)=0, Z=zeta^2!=0, t=v2(Z) even, theta a root of T",
        "square_strata": ["t=-2", "t>=2"],
        "nonsquare_strata": ["t=0", "t<=-4"],
        "positive_proof": (
            "with w=theta-2, theta^2-4=w^2*(1+4/w) and "
            "v2(4/w)=t+3/2>2; strong Hensel makes 1+4/w a square"
        ),
        "negative_proof": "monic shifted-octic congruence obstruction",
        "exceptional_proof": "ordinary resultant Newton polygon at t=-2",
        "factorization_translation": (
            "P(b)=b^4*T(b+b^-1) is 4+4 over Q_2 iff "
            "theta^2-4 is square in Q_2(theta)"
        ),
    }


def compressed_Z(z: F) -> F:
    return 8 * z * z / (1 + z * z)


def compressed_pullback() -> dict:
    valuation_rows = 0
    trace_rows = 0
    for k in range(-6, 7):
        for unit in map(F, (1, 3, 5, F(1, 3), F(5, 3))):
            z = power_of_two(k) * unit
            Z = compressed_Z(z)
            expected = 3 + 2 * k if k > 0 else (2 if k == 0 else 3)
            assert vp(Z, 2) == expected >= 2
            valuation_rows += 1

    # L28's t>0 resolvent calculation uses no parity of t once t>=2.
    for t in (2, 3, 4, 5, 7):
        for a in map(F, (1, 3, 5, F(1, 3), F(5, 3))):
            for Z_unit in map(F, (1, 3, 5, F(1, 3), F(5, 3))):
                Z = power_of_two(t) * Z_unit
                p, q, _, h = resolvent_coefficients(a, Z)
                assert vp(h, 2) == 4 - 4 * t
                assert unit_mod(h / power_of_two(4 - 4 * t), 8) == 3
                assert unit_mod(q * q / power_of_two(6 - 4 * t), 8) == 1
                for y_valuation, scale, residue in (
                    (2, 6 - 4 * t, 2),
                    (2 - 2 * t, 6 - 6 * t, 4),
                ):
                    for y_unit in map(F, (1, 3, 5, 7)):
                        Y = power_of_two(y_valuation) * y_unit * y_unit
                        value = Y**3 + 2 * p * Y**2 + h * Y - q * q
                        normalized = value / power_of_two(scale)
                        assert unit_mod(normalized, 8) == residue
                        trace_rows += 1

    # Exact odd-place target equivalence, including denominator zeros as poles.
    odd_rows = 0
    for prime in (3, 5, 7, 11, 13, 17, 19):
        for k in range(-3, 4):
            if k == 0:
                samples = (F(1), F(2), F(prime - 1))
            else:
                samples = (F(prime) ** k, 2 * F(prime) ** k)
            for z in samples:
                Z = compressed_Z(z)
                assert (vp(Z, prime) > 0) == (vp(z, prime) > 0)
                odd_rows += 1

    return {
        "type": "dyadic-compressed-reciprocal-pullback",
        "label": "PROVED exact target-preserving pullback",
        "map": "Z=8*z^2/(1+z^2)",
        "guard": "1+z^2 is nonzero over Q",
        "odd_place_equivalence": "v_p(Z)>0 iff v_p(z)>0 for every odd p",
        "dyadic_valuation": {
            "v2(z)>0": "3+2*v2(z)",
            "v2(z)=0": 2,
            "v2(z)<0": 3,
        },
        "consequences": [
            "v2(Z)>=2 on every nonzero rational fibre",
            "the L28 trace irreducibility proof extends verbatim to all integer t>=2",
            "theta^2-4 is uniformly square locally, so the reciprocal octic is 4+4 over Q_2",
            "2*(theta+2) remains nonsquare by the unchanged L28 norm argument",
            "the odd target set and formal witness count are unchanged",
        ],
        "valuation_rows": valuation_rows,
        "odd_target_rows": odd_rows,
        "trace_resolvent_rows": trace_rows,
        "zero_fibre": "z=0 has c=eta=0 and the exact point y=2/a^2, r=0",
    }


def target_squareclass_slice() -> dict:
    samples = 0
    for prime in (3, 5, 7, 11, 13, 17, 19):
        for a in range(1, 2 * prime, 2):
            A = 1 + 4 * a * a
            if A % prime and legendre(A, prime) == -1:
                for numerator, denominator in ((1, 1), (3, 1), (1, 3), (5, 3)):
                    rho = F(numerator, denominator)
                    b = prime * rho * rho
                    assert (2 * b) / (2 * prime) == rho * rho
                    samples += 1
                break
        else:
            raise AssertionError(("no target residue", prime))
    return {
        "type": "target-squareclass-freezing-slice",
        "label": "PROVED reduction; member existence remains OPEN",
        "slice": "b=w*rho^2",
        "identity": "(M,2*b)_v=(M,2*w)_v at every place v",
        "target_compatibility": "v_w(b)=1 when rho is a w-unit; b remains a 2-adic unit",
        "fixed_norm_problem": "U^2-2*w*V^2=P_rec(w*rho^2)*W^2",
        "endpoint_obstruction": (
            "P_rec(0)=A*(2*a^4*Z^2)^2 and (A,2*w)_w=(A|w)=-1, "
            "so the proper conic bundle has no Q(rho)-section"
        ),
        "strict_scope": (
            "freezing removes the moving quadratic field but not the inert-prime parity problem"
        ),
        "identity_rows": samples,
    }


def fixed_squareclass_evidence() -> dict:
    examples = (
        (3, 3, 5, -5),
        (5, 5, 3, -3),
        (7, 7, 1, -5),
        (11, 11, 5, -5),
        (13, 13, 3, -1),
        (17, 17, 1, -1),
        (19, 19, 3, -5),
        (23, 23, 3, -3),
    )
    rows = []
    for prime, z_value, a_value, rho_value in examples:
        a, z, rho = F(a_value), F(z_value), F(rho_value)
        Z = z * z
        b = prime * rho * rho
        A, D, _, c = bridge_data(a, Z, b)
        eta = reciprocal_tie(Z, D, b)
        delta = -4 * a**4 / A
        M = 16 - delta * c * c - 32 * A * b * eta * eta
        assert legendre(A, prime) == -1
        assert ramified(M, 2 * b) == []
        rows.append({
            "w": prime,
            "z": z_value,
            "a": a_value,
            "rho": rho_value,
            "b": str(b),
        })
    return {
        "type": "target-squareclass-freezing-evidence",
        "label": "EVIDENCE on eight exact fibres; not a member theorem",
        "rows": rows,
        "all_global_symbols_trivial": True,
        "refusals": 0,
    }


def fixed_squareclass_parity_probe() -> dict:
    a, Z, prime = F(5), F(9), 3
    A = 1 + 4 * a * a
    counts = {0: 0, 2: 0, 4: 0}
    refusals = 0
    zero_rows = []
    for rho_value in range(1, 100, 2):
        if rho_value % prime == 0:
            continue
        rho = F(rho_value)
        b = prime * rho * rho
        _, D, _, c = bridge_data(a, Z, b)
        eta = reciprocal_tie(Z, D, b)
        delta = -4 * a**4 / A
        M = 16 - delta * c * c - 32 * A * b * eta * eta
        try:
            bad_places = ramified(M, 2 * b)
        except (FactorBudget, PrimalityBound):
            refusals += 1
            continue
        assert len(bad_places) in counts
        counts[len(bad_places)] += 1
        if not bad_places:
            zero_rows.append(rho_value)
    assert counts == {0: 11, 2: 12, 4: 1}
    assert refusals == 9
    return {
        "type": "target-squareclass-parity-probe",
        "label": "EVIDENCE ONLY; exact decided rows, refusals excluded",
        "fixed_data": {"w": prime, "a": str(a), "Z": str(Z)},
        "rho_range": "positive odd w-adic units rho<100",
        "decided_bad_place_counts": {str(key): value for key, value in counts.items()},
        "globally_good_rho": zero_rows,
        "refusals": refusals,
        "interpretation": (
            "the fixed quadratic field produces many exact members but also "
            "persistent obstruction pairs; it does not bypass parity"
        ),
    }




def literature_record() -> dict:
    return {
        "type": "recent-literature-routing",
        "label": "SOURCE-CHECKED; no theorem imported beyond its scope",
        "sources": [
            {
                "url": "https://arxiv.org/abs/1904.12845",
                "result": "Loughran-Matthiesen count rational fibres in broad fibrations",
                "first_mismatch": "the reciprocal conic bundle has a non-split fibre over a non-rational closed point",
            },
            {
                "url": "https://arxiv.org/abs/2209.08949",
                "result": "Shute treats norm-form values using the beta sieve",
                "first_mismatch": "available fixed-field factor-degree hypotheses stop below the generic high-degree slice",
            },
            {
                "url": "https://arxiv.org/abs/2506.18065",
                "result": "Diao proves averaged Chowla/Bateman-Horn and norm-form results for random binary forms",
                "first_mismatch": "the H10/Q coefficient family is fixed and low-dimensional, not an almost-all coefficient cube",
            },
            {
                "url": "https://arxiv.org/abs/2608.16108",
                "result": "Wang proves a dynamical Chowla theorem on average",
                "first_mismatch": "an averaged dynamical statement does not control this fixed sign-decorated octic sequence",
            },
            {
                "url": "https://github.com/cjdoris/ExactpAdics",
                "result": "exact p-adic factorization code informed the local diagnostic phase",
                "first_mismatch": "no external package or computer factorization is used in the theorem proof",
            },
        ],
    }


def audit_record() -> dict:
    return {
        "type": "independent-derivation-audit",
        "label": "REVISED THEN VALIDATED",
        "first_pass": (
            "flagged the original Ore-MacLane presentation as under-specified; "
            "the proof was replaced rather than defended"
        ),
        "replacement": (
            "ordinary resultant Newton polygon plus valuation invariance under "
            "Q_2-embeddings"
        ),
        "second_pass": {
            "verdict": "main irreducibility-contradiction argument is sound",
            "confidence": 0.82,
            "validated": [
                "embeddings preserve the unique local valuation",
                "the resultant roots are phi(alpha_i)",
                "two resultant slopes contradict one irreducible orbit",
                "the first Newton polygon then forces 4+4",
            ],
        },
        "modular_caveat_checked": (
            "the searches enumerate every possible first monic Eisenstein "
            "quartic modulo 16 or 64; the high coefficients determine the "
            "second, and an exact Q_2 factorization would reduce to one"
        ),
    }


def frontier_record() -> dict:
    return {
        "type": "strict-frontier",
        "label": "NO unconditional quantifier or H10/Q claim",
        "closed": [
            "the complete dyadic squareclass classification of theta^2-4 on the square-Z L26 family",
            "a target-preserving pullback forcing the locally split reciprocal stratum",
            "the exact fixed-quadratic-field reduction on b=w*rho^2",
        ],
        "open": [
            "a globally good reciprocal member",
            "the inert-prime/two-large-divisor parity estimate",
            "a rational point on the L27 sheared cover",
            "any unconditional five- or six-quantifier improvement",
            "Hilbert's tenth problem over Q",
        ],
    }


def build_report(records: list[dict], elapsed: float) -> str:
    lines = [
        "# L29 reciprocal frontier",
        "",
        "All theorem labels are scope-limited exactly as recorded in JSONL.",
        "",
    ]
    for record in records:
        lines.append(f"- **{record['type']}** — {record['label']}")
    lines.extend(["", f"Elapsed: {elapsed:.3f} seconds.", ""])
    return "\n".join(lines)


def main() -> int:
    started = time.perf_counter()
    records = [
        congruence_obstruction(),
        exceptional_split_stratum(),
        reciprocal_lift_classification(),
        compressed_pullback(),
        target_squareclass_slice(),
        fixed_squareclass_evidence(),
        fixed_squareclass_parity_probe(),
        literature_record(),
        audit_record(),
        frontier_record(),
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in records),
        encoding="utf-8",
    )
    elapsed = time.perf_counter() - started
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print(f"wrote {len(records)} records to {OUT}")
    print(f"wrote report to {REPORT}")
    print(f"elapsed {elapsed:.3f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
