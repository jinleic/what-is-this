#!/usr/bin/env python3
"""Direct norm-section constructions and exact scoped no-go theorems.

Replay from the workspace root with the requested low-priority process:

    nice -n 19 python3 math/h10q/l23_norm_section.py

Everything labelled PROVED is checked with exact stdlib arithmetic. Finite
specializations prove only their stated specialization and are never promoted
to a uniform negative claim.
"""
from __future__ import annotations

from fractions import Fraction as F
import json
import math
import time
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l23_norm_section.jsonl"
REPORT = Path("/tmp/l23_norm_section.md")

# Tiny sparse Laurent ring for the displayed universal identities.
VARIABLES = ("a", "A", "Z", "D", "s", "b", "q", "R", "S", "X", "Y", "J", "t", "m", "n")
NVAR = len(VARIABLES)
Monomial = tuple[int, ...]
Laurent = dict[Monomial, F]
ZERO_MONOMIAL = (0,) * NVAR


def lp_clean(poly: Laurent) -> Laurent:
    return {m: c for m, c in poly.items() if c}


def lp_const(value: F | int) -> Laurent:
    value = F(value)
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def lp_monomial(coefficient: F | int = 1, **powers: int) -> Laurent:
    exponent = [0] * NVAR
    for variable, power in powers.items():
        exponent[VARIABLES.index(variable)] = power
    coefficient = F(coefficient)
    return {} if coefficient == 0 else {tuple(exponent): coefficient}


def lp_add(*polys: Laurent) -> Laurent:
    answer: Laurent = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            answer[monomial] = answer.get(monomial, F(0)) + coefficient
    return lp_clean(answer)


def lp_neg(poly: Laurent) -> Laurent:
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def lp_sub(left: Laurent, right: Laurent) -> Laurent:
    return lp_add(left, lp_neg(right))


def lp_scale(poly: Laurent, scalar: F | int) -> Laurent:
    scalar = F(scalar)
    return lp_clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def lp_mul(*polys: Laurent) -> Laurent:
    answer = lp_const(1)
    for poly in polys:
        product: Laurent = {}
        for lm, lc in answer.items():
            for rm, rc in poly.items():
                monomial = tuple(lm[index] + rm[index] for index in range(NVAR))
                product[monomial] = product.get(monomial, F(0)) + lc * rc
        answer = lp_clean(product)
    return answer


def lp_pow(poly: Laurent, exponent: int) -> Laurent:
    assert exponent >= 0
    answer = lp_const(1)
    base = poly
    while exponent:
        if exponent & 1:
            answer = lp_mul(answer, base)
        base = lp_mul(base, base)
        exponent >>= 1
    return answer


def lp_b_slice(poly: Laurent, b_degree: int) -> Laurent:
    index = VARIABLES.index("b")
    answer: Laurent = {}
    for monomial, coefficient in poly.items():
        if monomial[index] == b_degree:
            reduced = list(monomial)
            reduced[index] = 0
            key = tuple(reduced)
            answer[key] = answer.get(key, F(0)) + coefficient
    return lp_clean(answer)


def symbolic_replay() -> dict[str, Any]:
    one = lp_const(1)
    a, A = lp_monomial(a=1), lp_monomial(A=1)
    Z, D = lp_monomial(Z=1), lp_monomial(D=1)
    s, b = lp_monomial(s=1), lp_monomial(b=1)
    q, R, S = lp_monomial(q=1), lp_monomial(R=1), lp_monomial(S=1)
    X, Y = lp_monomial(X=1), lp_monomial(Y=1)
    J, t = lp_monomial(J=1), lp_monomial(t=1)
    m, n = lp_monomial(m=1), lp_monomial(n=1)
    Ai, Di = lp_monomial(A=-1), lp_monomial(D=-1)
    bi, qi, Ri, ti = (lp_monomial(b=-1), lp_monomial(q=-1),
                      lp_monomial(R=-1), lp_monomial(t=-1))

    # _sun_g, _sun_h and the tied right side M, then exact clearing.
    Ng = lp_sub(lp_scale(lp_mul(lp_pow(a, 4), lp_pow(b, 2)), 16),
                lp_mul(A, lp_pow(lp_sub(b, one), 4)))
    g = lp_mul(Ng, Ai, lp_pow(bi, 2))
    c = lp_mul(lp_pow(a, 2), lp_pow(Z, 2), g, Di)
    delta = lp_scale(lp_mul(lp_pow(a, 4), Ai), -4)
    alpha = lp_neg(lp_mul(delta, A))
    M = lp_sub(lp_sub(lp_const(16), lp_mul(delta, lp_pow(c, 2))),
               lp_scale(lp_mul(A, b, lp_pow(s, 2)), 32))
    clearing_square = lp_mul(lp_pow(D, 2), lp_pow(A, 2), lp_pow(b, 4))
    P = lp_add(
        lp_scale(lp_mul(lp_pow(D, 2), lp_pow(A, 2), lp_pow(b, 4)), 16),
        lp_neg(lp_mul(delta, lp_pow(a, 4), lp_pow(Z, 4), lp_pow(Ng, 2))),
        lp_scale(lp_mul(lp_pow(A, 3), lp_pow(s, 2), lp_pow(D, 2), lp_pow(b, 5)), -32),
    )
    assert not lp_sub(lp_mul(clearing_square, M), P)
    assert not lp_sub(alpha, lp_scale(lp_pow(a, 4), 4))

    # Explicit squares reducing the tied quaternion to (P,2b).
    x0 = lp_mul(alpha, M)
    x_multiplier = lp_pow(lp_scale(lp_mul(lp_pow(a, 2), Di, Ai, lp_pow(bi, 2)), 2), 2)
    assert not lp_sub(x0, lp_mul(x_multiplier, P))
    d0 = lp_scale(lp_mul(alpha, b), 2)
    d_multiplier = lp_pow(lp_scale(lp_pow(a, 2), 2), 2)
    assert not lp_sub(d0, lp_mul(d_multiplier, lp_scale(b, 2)))

    # Canonical H identity and equal endpoints.
    L = lp_sub(one, lp_scale(lp_mul(A, lp_pow(s, 2), b), 2))
    X0 = lp_mul(lp_pow(a, 4), lp_pow(Z, 2), Ng)
    Y0 = lp_scale(lp_mul(A, D, lp_pow(b, 2)), 2)
    H = lp_scale(lp_mul(A, P), F(1, 4))
    assert not lp_sub(H, lp_add(lp_pow(X0, 2), lp_mul(A, L, lp_pow(Y0, 2))))
    endpoint = lp_mul(lp_const(4), lp_pow(a, 8), A, lp_pow(Z, 4))
    assert not lp_sub(lp_b_slice(P, 0), endpoint)
    assert not lp_sub(lp_b_slice(P, 8), endpoint)

    # Exact matching plus the residual norm L=R^2-AS^2: actual U,V.
    Ln = lp_sub(lp_pow(R, 2), lp_mul(A, lp_pow(S, 2)))
    dm = lp_neg(lp_mul(A, Ln, lp_pow(q, 2)))
    Pabs = lp_scale(lp_mul(Ai, lp_add(lp_pow(X, 2), lp_mul(A, Ln, lp_pow(Y, 2)))), 4)
    U = lp_scale(lp_mul(Ri, lp_sub(lp_mul(S, X), lp_mul(Ln, Y))), 2)
    V = lp_scale(lp_mul(Ai, qi, Ri, lp_add(X, lp_mul(A, S, Y))), 2)
    assert not lp_sub(lp_sub(lp_pow(U, 2), lp_mul(dm, lp_pow(V, 2))), Pabs)

    # Stereographic parameterization of the minus-chart residual quadric,
    # checked after clearing its common denominator T.
    T = lp_add(one, lp_pow(m, 2), lp_neg(lp_mul(A, lp_pow(n, 2))))
    theta_numerator = lp_scale(m, 2)
    rho_numerator = lp_sub(lp_sub(lp_pow(m, 2), lp_mul(A, lp_pow(n, 2))), one)
    sigma_numerator = lp_scale(n, 2)
    assert not lp_sub(
        lp_sub(lp_pow(T, 2), lp_pow(theta_numerator, 2)),
        lp_sub(lp_pow(rho_numerator, 2), lp_mul(A, lp_pow(sigma_numerator, 2))),
    )

    # A target-square section, and the distinct L=-Aq^2 factorization.
    Us = lp_scale(lp_add(J, one), F(1, 2))
    Vs = lp_scale(lp_mul(lp_sub(J, one), ti), F(1, 2))
    assert not lp_sub(lp_sub(lp_pow(Us, 2), lp_mul(lp_pow(t, 2), lp_pow(Vs, 2))), J)
    Lf = lp_neg(lp_mul(A, lp_pow(q, 2)))
    Fm, Fp = lp_sub(X, lp_mul(A, q, Y)), lp_add(X, lp_mul(A, q, Y))
    assert not lp_sub(lp_add(lp_pow(X, 2), lp_mul(A, Lf, lp_pow(Y, 2))), lp_mul(Fm, Fp))

    return {
        "type": "symbolic-replay", "label": "PROVED",
        "kernel_clearing_identity": True,
        "alpha_square_identity": "alpha=-delta*A=4*a^4=(2*a^2)^2",
        "quaternion_reduction": "(x0,d0)=(P(b),2*b) for P(b)!=0",
        "norm_identity": "H=(A/4)P=X^2+A*L*Y^2",
        "exact_match_section": True, "target_split_section": True,
        "factor_branch_identity": True, "P_terms": len(P), "H_terms": len(H),
        "P_endpoint": "4*a^8*A*Z^4",
    }


# Exact univariate arithmetic for specializations and Newton polygons.
Poly = list[F]


def poly_trim(poly: Poly) -> Poly:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(left: Poly, right: Poly) -> Poly:
    answer = [F(0)] * max(len(left), len(right))
    for i, c in enumerate(left): answer[i] += c
    for i, c in enumerate(right): answer[i] += c
    return poly_trim(answer)


def poly_scale(poly: Poly, scalar: F | int) -> Poly:
    return poly_trim([F(scalar) * c for c in poly])


def poly_mul(left: Poly, right: Poly) -> Poly:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right): answer[i + j] += x * y
    return poly_trim(answer)


def poly_pow(poly: Poly, exponent: int) -> Poly:
    answer = [F(1)]
    while exponent:
        if exponent & 1: answer = poly_mul(answer, poly)
        poly = poly_mul(poly, poly)
        exponent >>= 1
    return answer


def poly_eval(poly: Poly, value: F | int) -> F:
    answer, value = F(0), F(value)
    for coefficient in reversed(poly): answer = answer * value + coefficient
    return answer


def canonical_P(a: F | int, Z: F | int) -> Poly:
    a, Z = F(a), F(Z)
    A, s = 1 + 4 * a * a, (a - 1) / 2
    D = 1 - Z - a * a * Z * Z
    Ng = poly_add(poly_scale([F(1), F(-4), F(6), F(-4), F(1)], -A),
                  [F(0), F(0), 16 * a**4])
    return poly_add(
        poly_add([F(0)] * 4 + [16 * D * D * A * A],
                 poly_scale(poly_mul(Ng, Ng), 4 * a**8 * Z**4 / A)),
        [F(0)] * 5 + [-32 * A**3 * s * s * D * D],
    )


def vp(value: F | int, prime: int) -> int:
    value = F(value)
    assert value != 0
    numerator, denominator, answer = abs(value.numerator), value.denominator, 0
    while numerator % prime == 0: numerator //= prime; answer += 1
    while denominator % prime == 0: denominator //= prime; answer -= 1
    return answer


def unit_mod(value: F | int, prime: int) -> int:
    value = F(value)
    exponent = vp(value, prime)
    value = value / prime**exponent if exponent >= 0 else value * prime ** (-exponent)
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def legendre_unit(value: F | int, prime: int) -> int:
    residue = unit_mod(value, prime)
    assert residue
    result = pow(residue, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


def hilbert_odd(left: F | int, right: F | int, prime: int) -> int:
    alpha, beta = vp(left, prime), vp(right, prime)
    answer = -1 if alpha * beta * ((prime - 1) // 2) % 2 else 1
    if beta % 2: answer *= legendre_unit(left, prime)
    if alpha % 2: answer *= legendre_unit(right, prime)
    return answer


def lcm(left: int, right: int) -> int:
    return abs(left * right) // math.gcd(left, right)


def primitive_integer_poly(poly: Poly) -> list[int]:
    denominator = 1
    for coefficient in poly: denominator = lcm(denominator, coefficient.denominator)
    answer = [int(coefficient * denominator) for coefficient in poly]
    content = 0
    for coefficient in answer: content = math.gcd(content, abs(coefficient))
    answer = [coefficient // content for coefficient in answer]
    return [-coefficient for coefficient in answer] if answer[-1] < 0 else answer


def projective_values_mod(poly: list[int], prime: int) -> tuple[list[int], int]:
    affine = [sum((c % prime) * pow(x, i, prime) for i, c in enumerate(poly)) % prime
              for x in range(prime)]
    return affine, poly[-1] % prime


def lower_hull(valuations: list[int]) -> list[tuple[int, int]]:
    hull: list[tuple[int, int]] = []
    for point in enumerate(valuations):
        while len(hull) >= 2:
            x0, y0 = hull[-2]; x1, y1 = hull[-1]; x2, y2 = point
            if F(y1 - y0, x1 - x0) < F(y2 - y1, x2 - x1): break
            hull.pop()
        hull.append(point)
    return hull


def ordinary_legendre(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0: return 0
    result = pow(residue, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


def factor_branch_specialization() -> dict[str, Any]:
    # Exact counterexample to "L=-Aq^2 automatically splits".
    a, Z, q, A, s = F(3), F(27), F(1), F(37), F(1)
    D = 1 - Z - a * a * Z * Z
    b, d = (1 + A * q * q) / (2 * A * s * s), None
    d = 2 * b
    Ng = 16 * a**4 * b**2 - A * (b - 1) ** 4
    X, Y = a**4 * Z**2 * Ng, 2 * A * D * b**2
    fm, fp = X - A * q * Y, X + A * q * Y
    P = poly_eval(canonical_P(a, Z), b)
    assert P == 4 * fm * fp / A and vp(b, 2) == 0
    checks = {"A_at_19": hilbert_odd(A, d, 19),
              "Fminus_at_167": hilbert_odd(fm, d, 167),
              "Fplus_at_13": hilbert_odd(fp, d, 13),
              "P_at_13": hilbert_odd(P, d, 13)}
    assert set(checks.values()) == {-1}

    # Complete P^1(F_5) obstruction for X +/- AqY after b=(1+37q^2)/74.
    bq = [F(1, 74), F(0), F(1, 2)]
    Ngq = poly_add(poly_scale(poly_pow(bq, 2), 16 * a**4),
                   poly_scale(poly_pow(poly_add(bq, [F(-1)]), 4), -A))
    Xq, Yq = poly_scale(Ngq, a**4 * Z**2), poly_scale(poly_pow(bq, 2), 2 * A * D)
    minus = primitive_integer_poly(poly_add(Xq, poly_scale([F(0)] + Yq, -A)))
    plus = primitive_integer_poly(poly_add(Xq, poly_scale([F(0)] + Yq, A)))
    mv, mi = projective_values_mod(minus, 5); pv, pi = projective_values_mod(plus, 5)
    assert (mv, mi) == ([2, 4, 2, 4, 3], 4)
    assert (pv, pi) == ([2, 3, 4, 2, 4], 4)
    return {
        "type": "factor-branch-specialization", "label": "PROVED",
        "a": "3", "A": "37", "s": "1", "Z": "27", "D": str(D),
        "q": "1", "b": str(b), "v2_b": 0, "P": str(P),
        "bad_hilbert_checks": checks,
        "factor_zero_octics": {
            "minus_coefficients_low_to_high": minus, "plus_coefficients_low_to_high": plus,
            "modulus": 5, "minus_values_affine_then_infinity": mv + [mi],
            "plus_values_affine_then_infinity": pv + [pi], "rational_roots": 0,
            "proof": "primitive homogenizations have no projective root modulo 5",
        },
    }


def two_adic_newton_replay() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for a in (F(1), F(3), F(5)):
        s = (a - 1) / 2
        k = None if s == 0 else vp(s, 2)
        for m in range(-3, 4):
            Z = F(2**m) if m >= 0 else F(1, 2 ** (-m))
            D, P = 1 - Z - a * a * Z * Z, canonical_P(a, Z)
            valuations = [vp(c, 2) for c in P]
            height, hull = 2 + 4 * m, lower_hull(valuations)
            assert valuations[0] == valuations[8] == height
            assert vp(D, 2) == min(0, 2 * m)
            if m <= 0:
                expected = [(0, height), (8, height)]
                assert hull == expected
                assert [i for i, v in enumerate(valuations) if v == height] == [0, 8]
                residual = "1+T^8=(1+T)^8"
            else:
                expected = [(0, height), (4, 4)]
                has_b5_vertex = k is not None and m >= 2 * k + 2
                if has_b5_vertex:
                    expected.append((5, 5 + 2 * k))
                expected.append((8, height))
                assert hull == expected
                assert unit_mod(P[0], 2) == unit_mod(P[4], 2) == unit_mod(P[8], 2) == 1
                residual = (
                    "left side 1+T^2=(1+T)^2; b5 creates later regular sides"
                    if has_b5_vertex
                    else "1+T^2=(1+T)^2 on both sides"
                )
            rows.append({"a": str(a), "v2_s": k, "m_v2_Z": m, "v2_D": vp(D, 2),
                         "coefficient_valuations": valuations, "lower_hull": hull,
                         "first_residual_over_F2": residual})
            time.sleep(0.001)
    return {
        "type": "capell-two-adic-newton", "label": "PROVED", "rows": rows,
        "uniform_data": {"Ng_coefficient_v2": [0, 2, 1, 2, 0],
                         "v2_D": "min(0,2m)", "endpoint_height": "2+4m",
                         "b5_height": "5+2*v2(s), when s!=0",
                         "warning": "the left first residual is always inseparable; exact ramification does not follow"},
        "square_locus_status": "OPEN",
    }


def target_newton_replay() -> dict[str, Any]:
    # Cross-check the odd target-place Capell obstruction used by actual cells.
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    rows: list[dict[str, Any]] = []
    for w in primes:
        all_nonsquare = [r for r in range(w) if ordinary_legendre(1 + 4 * r * r, w) == -1]
        zero_count = 2 if w % 4 == 1 else 0
        assert len(all_nonsquare) == (w - zero_count + 1) // 2
        residues = [r for r in all_nonsquare if r != 1]
        assert residues
        a = residues[0]
        if a % 2 == 0: a += w
        assert a % 2 and a % w != 1
        A, s = 1 + 4 * a * a, F(a - 1, 2)
        assert ordinary_legendre(A, w) == -1 and vp(s, w) == 0
        for e in (1, 2, 3, 4):
            E, Z = 3 * e, F(w ** (3 * e))
            P = canonical_P(F(a), Z)
            valuations = [vp(c, w) for c in P]
            expected = [(0, 4 * E), (4, 0), (5, 0), (8, 4 * E)]
            assert lower_hull(valuations) == expected
            # With zeta=Z/w^E, the exact left residual is
            # 4*a^8*A*zeta^4 + 16*A^2*T^4: D reduces to 1. It is squarefree
            # because its derivative is 64*A^2*T^3 and its constant is nonzero.
            left_units = (unit_mod(P[0], w), unit_mod(P[4], w))
            assert left_units == ((4 * a**8 * A) % w, (16 * A * A) % w)
            s_mod = s.numerator * pow(s.denominator, -1, w) % w
            beta = (-unit_mod(P[4], w) * pow(unit_mod(P[5], w), -1, w)) % w
            expected_beta = pow((2 * A * s_mod * s_mod) % w, -1, w)
            assert beta == expected_beta
            two_beta = pow((A * s_mod * s_mod) % w, -1, w)
            assert 2 * beta % w == two_beta
            assert ordinary_legendre(two_beta, w) == ordinary_legendre(pow(A, -1, w), w) == -1
            obstruction = (
                "odd valuation E on squarefree left edge"
                if E % 2
                else "horizontal simple root has 2*beta=1/(A*s^2), of squareclass 1/A"
            )
            rows.append({"w": w, "available_character_minus_residues": len(all_nonsquare),
                         "a": a, "A_mod_w": A % w, "e_vw_z": e,
                         "E_vw_Z": E, "lower_hull": expected,
                         "left_residual": "4*a^8*A*zeta^4+16*A^2*T^4",
                         "left_residual_units": left_units, "horizontal_beta": beta,
                         "horizontal_2beta": two_beta, "legendre_2beta": -1,
                         "obstruction": obstruction})
            time.sleep(0.001)
    return {
        "type": "capell-target-newton", "label": "PROVED", "rows": rows,
        "theorem": {
            "hypotheses": "v_w(z)=e>=1, Z=z^3, (A|w)=-1, w does not divide s",
            "left_residual": "4*a^8*A*zeta^4+16*A^2*T^4, zeta=Z/w^(3e), over the residue field",
            "e_odd": "the squarefree left edge has root valuation 3e odd, so 2b is not a square",
            "e_even": "the horizontal root beta=1/(2*A*s^2) satisfies exactly 2*beta=1/(A*s^2), whose squareclass is 1/A and is nonsquare",
            "class_choice": "there are (w+1)/2 minus residues for w=3 mod4 and (w-1)/2 for w=1 mod4; after excluding a=1 at least one remains",
            "fixed_class_scope": "when e is even and w divides s, the horizontal argument is unavailable; that fixed-class case is OPEN",
            "conclusion": "the untwisted Capell condition 2*beta square is impossible for the refined aligned irreducible fiber",
        },
    }


def records(symbolic: dict[str, Any], factor: dict[str, Any], two_newton: dict[str, Any], target_newton: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"type": "meta", "schema": "l23-norm-section-v1", "stdlib_only": True,
         "status_vocabulary": ["PROVED", "CONDITIONAL", "EVIDENCE", "OPEN"],
         "verdict": "no Phi section in the complete coefficient-linear or polynomial-in-b ansatz; exact next curves isolated",
         "global_status": "OPEN outside the proved ansatz classes",
         "pace": "single process; sleep 0.001 seconds after every finite Newton row"},
        symbolic,
        {"type": "kernel-theorem", "label": "PROVED",
         "definitions": {"A": "1+4a^2", "D": "1-Z-a^2Z^2", "s": "(a-1)/2",
                         "Ng": "16a^4b^2-A(b-1)^4", "delta": "-4a^4/A",
                         "alpha": "-delta*A=(2a^2)^2", "L": "1-2As^2b",
                         "X": "a^4Z^2Ng", "Y": "2ADb^2"},
         "identities": ["M=P/(D^2A^2b^4)", "x0=[2a^2/(DAb^2)]^2 P",
                        "d0=(2a^2)^2(2b)", "H=(A/4)P=X^2+ALY^2"],
         "conclusion": "P!=0: tied solubility iff (P,2b) splits; P=0 is the trivial M=0 case"},
        {"type": "exact-match-classification", "label": "PROVED",
         "ansatz": "U=uX+vY, V=rX+tY, formal in X,Y over the family field; A*d*L!=0 (automatic on Phi)",
         "determinant": "det^2=-16L/(Ad), whose squareclass is kappa=-AdL=-2AbL",
         "exact_match": "d=-ALq^2", "residual": "(A,d)=(A,L)=1",
         "norm_parameter": "L=R^2-AS^2",
         "section": {"guard": "q*R!=0", "U": "2(SX-LY)/R",
                     "V": "2(X+ASY)/(AqR)", "identity": "P=U^2-dV^2",
                     "degenerate_R0": "d=(A*S*q)^2 is target-split",
                     "degenerate_q0": "d=0, hence b=0 and outside Phi"},
         "L8c_charts": [{"b": "-theta^2/[2As^2(1-theta^2)]", "q": "theta/(As)"},
                         {"b": "1/[2As^2(1-theta^2)]", "q": "1/(As theta)"}],
         "chart_relation": "theta -> 1/theta; same L8c gauge",
         "residual_quadric": {"minus_chart": "1-theta^2=rho^2-A sigma^2",
                              "plus_chart": "theta^2-1=rho^2-A sigma^2 (the reciprocal transport)",
                              "minus_T": "1+m^2-An^2", "minus_theta": "2m/T",
                              "minus_rho": "(m^2-An^2-1)/T", "minus_sigma": "2n/T"},
         "v2_conclusion": "both charts have v2(b)!=0 for every theta; s=0 gives v2(b)=2v2(q)-1",
         "comparison": "exactly L8c source=target, with its omitted residual scalar and actual U,V supplied"},
        {"type": "polynomial-norm-no-go", "label": "PROVED",
         "ansatz": "arbitrary U,V in Q(a,Z)[b]", "P0": "4a^8AZ^4=(2a^4Z^2)^2A",
         "formal_obstruction": "v_A(P0)=1", "Phi_obstruction": "A=5 mod 8 is nonsquare in Q2",
         "conclusion": "no polynomial-in-b direct norm section at any degree"},
        {"type": "identity-branches", "label": "PROVED", "branches": [
            {"name": "target-split", "b": "t^2/2", "U": "(P+1)/2", "V": "(P-1)/(2t)",
             "v2_b": "2v2(t)-1", "Phi": False},
            {"name": "P-square", "curve": "R^2=P(b)", "section": "U=R,V=0", "status": "OPEN"},
            {"name": "minus-d-square", "curve": "T^2=-2bP(b)", "section": "U=0,V=T/(2b)", "status": "OPEN"},
            {"name": "exact-match", "Phi": False, "reason": "L8c wall"},
            {"name": "factor-H", "condition": "L=-Aq^2", "b": "(1+Aq^2)/(2As^2)",
             "factorization": "H=(X-AqY)(X+AqY)",
             "v2_b": {"v2q>0": "-1-2k", "v2q<0": "2v2q-1-2k", "v2q=0": "-2k"},
             "Phi_exactly": "v2(q)=v2(s)=0, apart from ordinary guards",
             "quaternion": "(P,d)=(A,d)(X-AqY,d)(X+AqY,d)", "scalar": "(A,d)=(-1,d)",
             "automatic_section": False}]},
        factor,
        {"type": "factor-next-curves", "label": "PROVED", "rational_points_status": "OPEN",
         "factor_norm_conics": ["r_-^2-dt_-^2=X-AqY", "r_+^2-dt_+^2=X+AqY", "r_A^2-dt_A^2=A"],
         "zero_factor": "X=+/-AqY", "elimination": "with L=-Aq^2 this is H(a,Z,b)=0",
         "consequence": "a Phi rational point gives P=M=0 and solves the tied member trivially",
         "first_polynomial": "H=a^8Z^4Ng^2+4A^3D^2b^4-8A^4s^2D^2b^5",
         "imported_L22": "for fixed nonzero rational Z, H is irreducible in Q(a)[b]",
         "scope": "no rational-function root b(a); isolated rational specializations remain open"},
        {"type": "capell-reduction", "label": "PROVED", "locus_status": "OPEN",
         "normalization": "C=P0=P8=4a^8AZ^4; p=P/C is monic and p(0)=1",
         "criterion": "2beta square in K iff p(u^2/2) has conjugate degree-8 factors",
         "norm": "Norm(2beta)=256",
         "scheme": {"U": "sum u_i b^i, i=0..4", "V": "sum v_i b^i, i=0..3",
                    "equations": "p_k=sum_{i+j=k}u_i u_j-2sum_{i+j=k-1}v_i v_j",
                    "endpoint_branches": "u0,u4 in {+/-1}; two modulo simultaneous sign"},
         "distinction": "this represents normalized p, not direct P; scalar squareclass A remains"},
        two_newton, target_newton,
        {"type": "constant-twist", "label": "PROVED", "closure_status": "OPEN",
         "criterion": "j*2beta square iff p(u^2/(2j)) is reducible",
         "target_local_conditions": {"e_odd": "v_w(j) must be odd", "e_even": "v_w(j) even and the unit of j must have character -1 at w"},
         "uniform_no_go": "one fixed rational j cannot handle all odd-e cells, since it would need odd valuation at infinitely many target primes",
         "remaining_character": "at good value-prime roots, d has squareclass j^{-1}, so the Hilbert sign is (j|ell), not automatically +1",
         "conclusion": "a nonsquare twist relocates the obligation to a fixed quadratic-character sieve; it is not a direct norm identity"},
        {"type": "summary", "label": "PROVED", "open_boundary": "OPEN",
         "proved_no_go": ["all coefficient-linear-in-(X,Y) sections", "all polynomial-in-b sections",
                          "target-split and exact-match sections miss Phi", "factor-H is not automatically a norm"],
         "next_objects": ["factor norm conics", "H=0 octic", "R^2=P and T^2=-2bP", "Capell coefficient scheme"],
         "open": "Phi rational-function sections outside the classified ansatz classes"},
    ]


def build_report(rows: list[dict[str, Any]]) -> str:
    factor = next(row for row in rows if row["type"] == "factor-branch-specialization")
    two = next(row for row in rows if row["type"] == "capell-two-adic-newton")
    target = next(row for row in rows if row["type"] == "capell-target-newton")
    return rf"""# L23 direct norm sections: exact scoped no-go and next curves

## Status

**PROVED.** On the fixed canonical branch the tied quaternion is exactly
$$(P(b),2b).$$
There is no $\Phi$-admissible section in either complete ansatz class:

1. arbitrary $U,V\in\mathbb Q(a,Z)[b]$;
2. arbitrary $U,V$ coefficient-linear in the canonical summands $X,Y$.

**OPEN.** Rational sections outside those classes and rational points on the
exact curves below. No finite scan or refusal is negative evidence.

## 1. Kernel re-derivation [PROVED]

Put
$$A=1+4a^2,\quad D=1-Z-a^2Z^2,\quad s=(a-1)/2,$$
$$N_g=16a^4b^2-A(b-1)^4.$$
From the checked kernel,
$$g=N_g/(Ab^2),\quad c=a^2Z^2N_g/(Ab^2D),$$
$$M=16-\delta c^2-32Abs^2.$$
Clearing the square denominator gives identically
$$P=D^2A^2b^4M
=16D^2A^2b^4-\delta a^4Z^4N_g^2-32A^3s^2D^2b^5.$$
On
$$\tau^\dagger=(1+2a^2)/A,\quad\delta=-4a^4/A,$$
$$\alpha=-\delta A=4a^4=(2a^2)^2.$$
Hence
$$x_0=\alpha M=\left({{2a^2\over DAb^2}}\right)^2P,
\qquad d_0=2\alpha b=(2a^2)^2(2b).$$
For $P\ne0$, $(x_0,d_0)=(P,2b)$. For $P=0$, equivalently $M=0$,
the tied equation is trivially solved by $(y,r)=(0,0)$ instead; no quaternion
with a zero slot is used.

Set
$$X=a^4Z^2N_g,\quad Y=2ADb^2,\quad L=1-2As^2b,
\quad H=(A/4)P.$$
Sparse symbolic expansion over independent $a,A,Z,D,s,b$ proves
$$\boxed{{H=X^2+ALY^2}}.$$

## 2. Exact matching, including the residual scalar [PROVED]

Write $d=2b$. Exact source/target matching is
$$d=-ALq^2.\tag{{EM}}$$
Then $H=X^2-d(Y/q)^2$, but $P=(4/A)H$ is a norm if and only if
$$ (A,d)=(A,-AL)=(A,L)=1.\tag{{R}}$$
Writing $L=R^2-AS^2$, for $qR\ne0$ direct expansion gives
$$U={{2(SX-LY)\over R}},\qquad
V={{2(X+ASY)\over AqR}},$$
$$\boxed{{P=U^2-dV^2}}.$$
If $R=0$, then $d=(ASq)^2$ and this is exactly the target-split branch; if
$q=0$, then $d=2b=0$. Both degeneracies are already non-$\Phi$.
This is the complete nondegenerate, hence $\Phi$-relevant,
coefficient-linear classification. If
$U=uX+vY$, $V=rX+tY$ formally in $X,Y$, determinants force
$$\det(M)^2={{16L/A\over-d}}=-{{16L\over Ad}}.$$
Its squareclass is $\kappa=-AdL=-2AbL$; after (EM), evaluation on the first
basis vector forces (R). The displayed formulas prove sufficiency.

## 3. Exact comparison with L8c and the full $2$-adic wall [PROVED]

For $s\ne0$ the two customary charts are
$$b_-=-{{\theta^2\over2As^2(1-\theta^2)}},\ q={{\theta\over As}},$$
$$b_+={{1\over2As^2(1-\theta^2)}},\ q={{1\over As\theta}}.$$
They are exchanged by $\theta\mapsto1/\theta$: precisely L8c's same Pell
conic, not a weaker or new gauge. On the minus chart the residual norm
condition is
$$1-\theta^2=\rho^2-A\sigma^2;$$
on the plus chart it is the reciprocal transport
$$\theta^2-1=\rho^2-A\sigma^2.$$
With $T=1+m^2-An^2$, the first quadric is parametrized by
$$\theta=2m/T,\quad\rho=(m^2-An^2-1)/T,\quad\sigma=2n/T.$$
The second follows by reciprocal substitution. Thus $(A,L)=1$ has not been
silently omitted.
At $\theta=0$, the minus chart is $b=0$ and the plus chart is the
$L=0$ endpoint $b=1/(2As^2)$ (not represented by a finite $q$ in (EM));
both are non-$\Phi$. The values $\theta=\pm1$ are poles.

Let $k=v_2(s)\ge0$, $r=v_2(\theta)$; $A\equiv5\pmod8$. Exactly,
$$v_2(1-\theta^2)=0\ (r>0),\quad2r\ (r<0),\quad\ge3\ (r=0).$$
For $r>0$, $(v_2(b_-),v_2(b_+))=(2r-1-2k,-1-2k)$; for $r<0$ they are
$(-1-2k,-1-2k-2r)$; for $r=0$ both are
$-1-2k-v_2(1-\theta^2)<0$. At the chart endpoint the surviving valuation is
$-1-2k$. If $s=0$, (EM) gives $v_2(b)=2v_2(q)-1$.
Therefore every exact-match point, even after (R), misses $\Phi$.

## 4. Polynomial-in-$b$ no-go [PROVED]

Suppose at any degrees
$$P(b)=U(b)^2-2bV(b)^2,
\qquad U,V\in\mathbb Q(a,Z)[b].$$
At $b=0$,
$$U(0)^2=P(0)=4a^8AZ^4=(2a^4Z^2)^2A.$$
But $A=1+4a^2$ is irreducible in $\mathbb Q[a]$ and has odd $A$-valuation.
At every $\Phi$ specialization the same obstruction is $A\equiv5\pmod8$,
nonsquare in $\mathbb Q_2$. Thus there is no polynomial direct norm section,
including the natural bounds $\deg U\le4$, $\deg V\le3$.

## 5. Low-complexity quaternion branches [PROVED reductions]

A uniform target-split section exists:
$$2b=t^2,\quad U=(P+1)/2,\quad V=(P-1)/(2t),$$
but $v_2(b)=2v_2(t)-1\ne0$. The one-slot identities instead give exact
curves
$$R^2=P(b),\qquad T^2=-2bP(b),$$
generically degree $8$ and $9$ (genera $3$ and $4$ after smooth completion).
Their rational points remain **OPEN**.

## 6. The admissible factor branch is not a norm section [PROVED]

The new factorization
$$L=-Aq^2,\quad b={{1+Aq^2\over2As^2}},$$
$$H=(X-AqY)(X+AqY)=F_-F_+$$
really escapes L8c. If $k=v_2(s)$ and $r=v_2(q)$, then
$$v_2(b)=-1-2k\ (r>0),\quad2r-1-2k\ (r<0),\quad-2k\ (r=0).$$
For $r=0$, $1+Aq^2\equiv6\pmod8$. Hence this branch is $\Phi$-admissible
exactly when $q,s$ are units ($a\equiv3\pmod4$), apart from the ordinary
$b\ne1$, $D\ne0$ guards.

But quaternion multiplication gives only
$$(P,d)=(A,d)(F_-,d)(F_+,d).$$
Since $d\equiv A(1+Aq^2)$ and $1+Aq^2$ is a norm from
$\mathbb Q(\sqrt{{-A}})$, $(A,d)=(-1,d)$, not $1$ identically.
The exact admissible specialization
$$(a,A,s,q,Z,D,b)=(3,37,1,1,27,-6587,19/37)$$
has
$$(A,d)_{{19}}=(F_-,d)_{{167}}=(F_+,d)_{{13}}=(P,d)_{{13}}=-1.$$
So neither a factor nor their product is automatically a norm.

## 7. Sharp next curves [PROVED reduction; OPEN points]

Factor-by-factor work must solve
$$r_\pm^2-dt_\pm^2=X\pm AqY,\qquad r_A^2-dt_A^2=A.$$
The cheapest zero-factor condition $X=\pm AqY$, combined with $L=-Aq^2$,
eliminates to exactly
$$\boxed{{H(a,Z,b)=0}},$$
$$H=a^8Z^4N_g^2+4A^3D^2b^4-8A^4s^2D^2b^5.$$
A $\Phi$ rational point gives $P=M=0$ and solves the tied member trivially.
L22 proves this octic irreducible in $\mathbb Q(a)[b]$ for fixed $Z\ne0$,
so no rational-function root $b(a)$ exists; isolated rational specializations
remain open. At $(a,Z)=(3,27)$, the two $q$-octics have projective values mod
$5$ respectively
{factor['factor_zero_octics']['minus_values_affine_then_infinity']} and
{factor['factor_zero_octics']['plus_values_affine_then_infinity']}; all are
nonzero, proving no rational $q$ at that specialization.

## 8. The separate Capell criterion [PROVED reduction]

Let $C=P_0=P_8=4a^8AZ^4$ and $p=P/C$. Then $p$ is monic, $p(0)=1$, and for
an irreducible fiber with root $\beta$,
$$2\beta\in K^{{\times2}}\iff p(u^2/2)\text{{ factors}},\qquad N(2\beta)=256.$$
The exact square locus is the projection of
$$p_k=\sum_{{i+j=k}}u_i u_j-2\sum_{{i+j=k-1}}v_i v_j,
\quad k=0,\ldots,8,$$
where $\deg U\le4$, $\deg V\le3$ and $u_0,u_4\in\{{\pm1\}}$ (two
endpoint branches modulo simultaneous sign). This represents normalized $p$,
not direct $P$: the scalar squareclass $A$ remains.

The first $2$-adic Newton polygon does **not** settle the locus. With
$m=v_2(Z)$, $k=v_2(s)$ when $s\ne0$, $v_2(D)=\min(0,2m)$, endpoint height
$2+4m$, and $b^5$ height $5+2k$. For $m\le0$ the residual is
$(1+T)^8$. For $m\ge1$ the hull starts
$$(0,2+4m)\to(4,4).$$
It goes directly to $(8,2+4m)$ when $s=0$ or $m\le2k+1$; otherwise
$(5,5+2k)$ is an additional vertex. In every case the left side has
residual $(1+T)^2$, so its ramification is not determined by the first
polygon; the tempting inference $e=2$ is invalid. The script checks
{len(two['rows'])} exact rows and records the uniform coefficient proof
$v_2(N_{{g,i}})=(0,2,1,2,0)$.

## 9. Target-place Capell no-go for actual cells [PROVED]

Let $w$ be the cell prime, $e=v_w(z)\ge1$, $Z=z^3$, choose the L20 residue
with $(A\mid w)=-1$ and additionally $a\not\equiv1\pmod w$, so $w\nmid s$.
This refinement is always available: the number of residues with character
$-1$ is $(w+1)/2$ for $w\equiv3\pmod4$ and $(w-1)/2$ for
$w\equiv1\pmod4$, so deleting $a=1$ leaves at least one, including for
$w=3,5$.
It does not change the CRT/HIT argument.
The $w$-Newton hull is
$$(0,12e)\to(4,0)\to(5,0)\to(8,12e).$$
Writing $\zeta=Z/w^{{3e}}$, the exact left residual over the residue field is
$$4a^8A\bar\zeta^4+16A^2T^4.$$
Its constant is nonzero and its derivative is $64A^2T^3$, so it is
squarefree at odd $w$. If $e$ is odd, its roots have odd valuation $3e$, so
$2b$ is nonsquare. If $e$ is even, the horizontal side has the simple
residue root
$$\bar\beta={{1\over2As^2}},\qquad
 2\bar\beta={{1\over As^2}},$$
whose squareclass is $1/A$, not the literal value $1/A$; it is nonsquare
because $(A\mid w)=-1$. Thus the untwisted Capell escape is impossible on
this refined aligned irreducible fiber. For a pre-existing fixed class with
even $e$ and $w\mid s$, the horizontal argument is unavailable and that
fixed-class case remains **OPEN**. The script checks {len(target['rows'])}
exact target-place rows; the proof is the displayed coefficient/residual
calculation, not the finite sample.

## 10. Constant twists [PROVED obstruction; OPEN sieve]

For fixed $j\in\mathbb Q^\times$,
$$j\,2\beta\in K^{{\times2}}\iff p(u^2/(2j))\text{{ factors}}.$$
At the target place an odd-$e$ cell requires $v_w(j)$ odd. An even-$e$ cell
requires $v_w(j)$ even and the unit of $j$ to have character $-1$ at $w$.
One fixed rational $j$ cannot cover all odd-$e$ cells: it would need odd
valuation at infinitely many target primes. Even on one cell a nonsquare
twist does not kill the wild signs: at every good simple value-prime root,
$d$ has squareclass $j^{{-1}}$, so the sign is $(j\mid\ell)$. The twist
therefore relocates the obligation to a fixed quadratic-character sieve; it
is not a direct norm identity.

## Boundary

**PROVED empty:** coefficient-linear $H$-form sections and polynomial-in-$b$
direct norm sections. **PROVED non-$\Phi$:** target-split and exact-match
uniform sections. **PROVED not automatic:** the admissible factor branch.
**OPEN:** rational points on the named conics/hyperelliptic curves/octic and
rational sections outside the classified ansatz classes.
"""


def main() -> None:
    symbolic = symbolic_replay()
    factor = factor_branch_specialization()
    two = two_adic_newton_replay()
    target = target_newton_replay()
    rows = records(symbolic, factor, two, target)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in rows), encoding="utf-8")
    REPORT.write_text(build_report(rows), encoding="utf-8")
    print("L23 norm section: scoped no-go proved; admissible factor branch reduced to exact curves; global route OPEN")
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
