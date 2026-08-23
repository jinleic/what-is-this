#!/usr/bin/env python3
"""Exact replay for replacing one Schinzel polynomial by its odd factors.

The mathematical result is deliberately narrower than a uniform theorem for
``tau_dagger``.  It proves an exact factor-tuple/progression criterion, and it
proves by an admissible linear toy tuple that classical Schinzel H does not
supply the missing sign compatibility.

Replay from the repository root:

    nice -n 19 python3 math/h10q/l22_factor_tuple.py

Only standard-library arithmetic and the checked-in h10q Hilbert kernel are
used.  Finite rows are mechanism checks, never a uniform irreducibility claim.
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import h10q  # noqa: E402  (authority for P and Hilbert symbols)

OUT = HERE / "data" / "l22_factor_tuple.jsonl"
REPORT = Path("/tmp/l22_factor_tuple.md")
L21 = HERE / "data" / "l21_reducible_locus.jsonl"
PACE_EVERY = 200
PACE_SECONDS = 0.005


class Pacer:
    """Apply the requested low-duty-cycle pause to every finite search loop."""

    def __init__(self) -> None:
        self.iterations = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.iterations += 1
        if self.iterations % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()


def frac_text(value: F | int) -> str:
    q = F(value)
    return str(q.numerator) if q.denominator == 1 else str(q)


def trim(poly: list[int] | list[F]) -> list:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_eval(poly: list[int] | list[F], value: int | F) -> int | F:
    answer: int | F = 0
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def poly_mul(left: list[int] | list[F], right: list[int] | list[F]) -> list:
    out = [left[0] * 0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if x == 0:
            continue
        for j, y in enumerate(right):
            out[i + j] += x * y
    return trim(out)


def compose_linear(poly: list[int] | list[F], offset: int | F, step: int | F) -> list[F]:
    """Ascending coefficients of poly(offset + step*u), exactly."""
    out = [F(0)]
    power = [F(1)]
    linear = [F(offset), F(step)]
    for coefficient in poly:
        if len(out) < len(power):
            out.extend([F(0)] * (len(power) - len(out)))
        for i, value in enumerate(power):
            out[i] += F(coefficient) * value
        power = poly_mul(power, linear)
    return trim(out)


def primitive_integer_poly(poly: list[int] | list[F]) -> tuple[list[int], F]:
    """Return positive-leading primitive H and c with poly = c*H."""
    values = [F(value) for value in trim(poly)]
    denominator = math.lcm(*(value.denominator for value in values))
    integers = [value.numerator * (denominator // value.denominator) for value in values]
    content = math.gcd(*(abs(value) for value in integers))
    assert content > 0
    integers = [value // content for value in integers]
    scale = F(content, denominator)
    if integers[-1] < 0:
        integers = [-value for value in integers]
        scale = -scale
    assert [scale * value for value in integers] == values
    return integers, scale


def rational_sqrt(value: F) -> F | None:
    if value < 0:
        return None
    rn = math.isqrt(value.numerator)
    rd = math.isqrt(value.denominator)
    if rn * rn == value.numerator and rd * rd == value.denominator:
        return F(rn, rd)
    return None


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        PACER.tick()
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def primes_up_to(limit: int) -> Iterable[int]:
    for candidate in range(2, limit + 1):
        PACER.tick()
        if is_prime(candidate):
            yield candidate


def ff_trim(poly: list[int], prime: int) -> list[int]:
    out = [value % prime for value in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def ff_divmod(left: list[int], right: list[int], prime: int) -> tuple[list[int], list[int]]:
    remainder = ff_trim(left, prime)
    divisor = ff_trim(right, prime)
    assert divisor != [0]
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        for i, value in enumerate(divisor):
            remainder[i + shift] = (remainder[i + shift] - coefficient * value) % prime
        remainder = ff_trim(remainder, prime)
    return ff_trim(quotient, prime), remainder


def ff_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    a, b = ff_trim(left, prime), ff_trim(right, prime)
    while b != [0]:
        _, remainder = ff_divmod(a, b, prime)
        a, b = b, remainder
    inverse = pow(a[-1], -1, prime)
    return [(value * inverse) % prime for value in a]


def ff_mul_mod(left: list[int], right: list[int], modulus: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            product[i + j] = (product[i + j] + x * y) % prime
    return ff_divmod(product, modulus, prime)[1]


def ff_pow_mod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    answer = [1]
    power = ff_divmod(base, modulus, prime)[1]
    while exponent:
        if exponent & 1:
            answer = ff_mul_mod(answer, power, modulus, prime)
        power = ff_mul_mod(power, power, modulus, prime)
        exponent >>= 1
    return answer


def prime_divisors(value: int) -> list[int]:
    value = abs(value)
    answer: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        PACER.tick()
        if value % divisor == 0:
            answer.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        answer.append(value)
    return answer


def irreducible_mod(poly: list[int], prime: int) -> bool:
    reduced = ff_trim(poly, prime)
    degree = len(poly) - 1
    if len(reduced) - 1 != degree:
        return False
    inverse = pow(reduced[-1], -1, prime)
    modulus = [(value * inverse) % prime for value in reduced]
    if degree == 1:
        return True
    checks = {degree // divisor for divisor in prime_divisors(degree)}
    x = [0, 1]
    frobenius = x
    for index in range(1, degree + 1):
        frobenius = ff_pow_mod(frobenius, prime, modulus, prime)
        difference = list(frobenius)
        if len(difference) < 2:
            difference.extend([0] * (2 - len(difference)))
        difference[1] = (difference[1] - 1) % prime
        difference = ff_trim(difference, prime)
        if index in checks and len(ff_gcd(modulus, difference, prime)) != 1:
            return False
    return ff_trim(difference, prime) == [0]


def irreducibility_certificate(poly: list[int], limit: int = 2000) -> int:
    """A full-degree irreducible reduction proves irreducibility over Q."""
    if len(poly) == 2:
        return 2 if poly[-1] % 2 else 3
    for prime in primes_up_to(limit):
        if irreducible_mod(poly, prime):
            return prime
    raise AssertionError(("no modular irreducibility certificate", poly, limit))


def fixed_divisor(polys: list[list[int]]) -> int:
    """Exact product fixed divisor from degree+1 consecutive values."""
    degree = sum(len(poly) - 1 for poly in polys)
    divisor = 0
    for value in range(degree + 1):
        PACER.tick()
        product = 1
        for poly in polys:
            product *= int(poly_eval(poly, value))
        divisor = math.gcd(divisor, abs(product))
    return divisor


def squarefree_kernel(value: int) -> int:
    """Signed squarefree representative of a nonzero integer square class."""
    assert value != 0
    sign = -1 if value < 0 else 1
    residue = abs(value)
    result = 1
    divisor = 2
    while divisor * divisor <= residue:
        PACER.tick()
        parity = 0
        while residue % divisor == 0:
            residue //= divisor
            parity ^= 1
        if parity:
            result *= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if residue > 1:
        result *= residue
    return sign * result


def rational_squareclass(value: F) -> int:
    """Integer K representing value in Q*/Q*2."""
    assert value != 0
    return squarefree_kernel(value.numerator * value.denominator)


def jacobi(numerator: int, denominator: int) -> int:
    """Jacobi symbol for a positive odd denominator."""
    assert denominator > 0 and denominator % 2 == 1
    numerator %= denominator
    answer = 1
    while numerator:
        while numerator % 2 == 0:
            numerator //= 2
            if denominator % 8 in (3, 5):
                answer = -answer
        numerator, denominator = denominator, numerator
        if numerator % 4 == denominator % 4 == 3:
            answer = -answer
        numerator %= denominator
    return answer if denominator == 1 else 0


def linear_resultant(Q: list[int], factor: list[int]) -> int:
    """Res(Q,f)=N^m f(-q1/N), for Q=q1+N*t."""
    assert len(Q) == 2 and Q[1] > 0
    degree = len(factor) - 1
    value = F(Q[1]) ** degree * F(poly_eval(factor, F(-Q[0], Q[1])))
    assert value.denominator == 1 and value != 0
    return value.numerator


def signature_period(K: int, residue_classes: list[int]) -> int:
    """A safe (not necessarily minimal) period for all fixed characters."""
    moduli = [4 * abs(K), 4]
    moduli.extend(4 * abs(value) for value in residue_classes)
    return math.lcm(*moduli)


def sign_formula(
    t: int,
    Q: list[int],
    factor: list[int],
    K: int,
    resultant_squareclass: int,
) -> tuple[int, int]:
    """Return (resultant/reciprocity formula, direct Jacobi symbol)."""
    q = int(poly_eval(Q, t))
    r = int(poly_eval(factor, t))
    assert q > 0 and r > 0 and q % 2 == r % 2 == 1
    assert math.gcd(q, r) == math.gcd(K, r) == math.gcd(resultant_squareclass, q) == 1
    reciprocity = -1 if q % 4 == r % 4 == 3 else 1
    formula = jacobi(K, r) * reciprocity * jacobi(resultant_squareclass, q)
    direct = jacobi(K * q, r)
    assert formula == direct and formula in (-1, 1)
    return formula, direct


def safe_representative(
    residue: int,
    period: int,
    Q: list[int],
    factors: list[list[int]],
    K: int,
    resultant_squareclasses: list[int],
    start_shift: int = 1,
) -> int:
    """Find a positive representative avoiding the finite resultant primes."""
    for shift in range(start_shift, start_shift + 2000):
        PACER.tick()
        t = residue + shift * period
        q = int(poly_eval(Q, t))
        values = [int(poly_eval(factor, t)) for factor in factors]
        if q <= 1 or any(value <= 1 or value % 2 == 0 for value in values):
            continue
        if any(math.gcd(q, value) != 1 for value in values):
            continue
        if any(math.gcd(K, value) != 1 for value in values):
            continue
        if any(math.gcd(squareclass, q) != 1 for squareclass in resultant_squareclasses):
            continue
        return t
    raise AssertionError(("no safe representative", residue, period))


def analyze_pattern(
    *,
    name: str,
    kind: str,
    Q: list[int],
    factors: list[list[int]],
    exponents: list[int],
    d_over_Q: F,
    metadata: dict | None = None,
) -> dict:
    """Certify factor tuple, periodic signs, and every compatible L-class."""
    assert len(factors) == len(exponents) and all(exponent >= 1 for exponent in exponents)
    assert math.gcd(*[abs(value) for value in Q]) == 1 and Q[-1] > 0
    assert len({tuple(poly) for poly in factors}) == len(factors)

    factor_certificates = [irreducibility_certificate(poly) for poly in factors]
    odd_indices = [index for index, exponent in enumerate(exponents) if exponent % 2]
    odd_factors = [factors[index] for index in odd_indices]
    tuple_polys = [Q, *odd_factors]
    divisor = fixed_divisor(tuple_polys)
    admissible = divisor == 1
    K = rational_squareclass(d_over_Q)

    resultants = [linear_resultant(Q, factor) for factor in odd_factors]
    resultant_squareclasses = [
        squarefree_kernel(resultant * Q[1] ** (len(factor) - 1))
        for resultant, factor in zip(resultants, odd_factors)
    ]
    period = signature_period(K, resultant_squareclasses)
    period_primes = prime_divisors(period)

    compatible_residues: list[int] = []
    vector_counts: dict[str, int] = {}
    direct_formula_checks = 0
    period_checks = 0
    refined_checks = 0
    first_all_plus: int | None = None

    if admissible:
        for residue in range(period):
            PACER.tick()
            if any(
                int(poly_eval(poly, residue)) % prime == 0
                for prime in period_primes
                for poly in tuple_polys
            ):
                continue
            compatible_residues.append(residue)
            t = safe_representative(
                residue, period, Q, odd_factors, K, resultant_squareclasses
            )
            signs = []
            for factor, squareclass in zip(odd_factors, resultant_squareclasses):
                formula, direct = sign_formula(t, Q, factor, K, squareclass)
                assert formula == direct
                signs.append(formula)
                direct_formula_checks += 1
            key = "".join("+" if sign == 1 else "-" for sign in signs) or "vacuous"
            vector_counts[key] = vector_counts.get(key, 0) + 1
            if all(sign == 1 for sign in signs) and first_all_plus is None:
                first_all_plus = residue

            t_next = safe_representative(
                residue,
                period,
                Q,
                odd_factors,
                K,
                resultant_squareclasses,
                start_shift=(t - residue) // period + 1,
            )
            next_signs = [
                sign_formula(t_next, Q, factor, K, squareclass)[0]
                for factor, squareclass in zip(odd_factors, resultant_squareclasses)
            ]
            assert next_signs == signs
            period_checks += len(signs)

            refined = [
                primitive_integer_poly(compose_linear(poly, residue, period))
                for poly in tuple_polys
            ]
            assert all(scale == 1 for _, scale in refined)
            assert fixed_divisor([poly for poly, _ in refined]) == 1
            refined_checks += 1

    row = {
        "type": "mechanism_test",
        "label": "PROVED",
        "scope": "exact named factor pattern; not a tau_dagger uniform claim",
        "name": name,
        "kind": kind,
        "Q_coefficients_ascending": Q,
        "factor_coefficients_ascending": factors,
        "factor_exponents": exponents,
        "factor_irreducibility_certificate_primes": factor_certificates,
        "odd_factor_indices": odd_indices,
        "schinzel_tuple_size": len(tuple_polys),
        "tuple_fixed_divisor": divisor,
        "tuple_admissible": admissible,
        "d_over_Q": frac_text(d_over_Q),
        "squareclass_K": K,
        "linear_resultants": resultants,
        "resultant_formula_squareclasses": resultant_squareclasses,
        "signature_period": period,
        "signature_period_primes": period_primes,
        "compatible_residues": len(compatible_residues),
        "signature_vector_counts": vector_counts,
        "all_plus_residues": vector_counts.get("+" * len(odd_factors), 0),
        "first_all_plus_residue": first_all_plus,
        "direct_equals_resultant_formula_checks": direct_formula_checks,
        "periodicity_checks": period_checks,
        "refined_tuple_admissibility_checks": refined_checks,
        "metadata": metadata or {},
    }
    return row


def l21_rows() -> list[dict]:
    return [json.loads(line) for line in L21.read_text(encoding="utf-8").splitlines()]


def l21_pattern(
    z: F,
    tau: F,
    q1: int,
    step: int,
    rows: list[dict],
    name: str,
) -> tuple[dict, dict]:
    """Rebuild one L21a factorization and its primitive t-factors exactly."""
    matches = [
        row
        for row in rows
        if row.get("type") == "l21a"
        and F(row["z"]) == z
        and F(row["tau"]) == tau
    ]
    assert len(matches) == 1 and matches[0]["applicable"] and matches[0]["identity_holds"]
    source_row = matches[0]

    a = F(1)
    A = F(5)
    Z = z**3
    D = 1 - Z - Z * Z
    delta = 1 - A * tau * tau
    sigma = rational_sqrt(delta)
    assert sigma is not None and sigma > 0 and D != 0 and Z != 0
    ng = [F(-5), F(20), F(-14), F(20), F(-5)]

    minus = [-sigma * Z * Z * coefficient for coefficient in ng]
    plus = [sigma * Z * Z * coefficient for coefficient in ng]
    minus[2] += 4 * D * A
    plus[2] += 4 * D * A
    positive_factors_b = [minus, [-coefficient for coefficient in plus]]
    assert all(factor[-1] > 0 for factor in positive_factors_b)

    P = h10q._l10_P(a, Z, D, A, delta, F(0))
    assert poly_mul(minus, plus) == P

    Q = [q1, step]
    primitive_factors: list[list[int]] = []
    factor_scales: list[F] = []
    for factor in positive_factors_b:
        composed = compose_linear(factor, q1, step)
        primitive, scale = primitive_integer_poly(composed)
        primitive_factors.append(primitive)
        factor_scales.append(scale)

    primitive_P, P_scale = primitive_integer_poly(compose_linear(P, q1, step))
    factor_product = poly_mul(primitive_factors[0], primitive_factors[1])
    assert factor_product == primitive_P
    assert P_scale == -factor_scales[0] * factor_scales[1]

    alpha = -delta * A
    metadata = {
        "l21_source_row": {
            "z": source_row["z"],
            "tau": source_row["tau"],
            "sigma": source_row["sigma"],
            "factor_degrees": source_row["factor_degrees"],
            "identity_holds": source_row["identity_holds"],
        },
        "class_substitution": f"b={q1}+{step}t",
        "primitive_factor_scales": [frac_text(value) for value in factor_scales],
        "primitive_P_scale": frac_text(P_scale),
        "primitive_product_identity": True,
        "alpha": frac_text(alpha),
        "branch_warning": "L21a is deliberately outside tau_dagger and is only a mechanism test",
    }
    row = analyze_pattern(
        name=name,
        kind="l21a_reducible_P",
        Q=Q,
        factors=primitive_factors,
        exponents=[1, 1],
        d_over_Q=2 * alpha,
        metadata=metadata,
    )
    return row, {
        "P": P,
        "positive_factors_b": positive_factors_b,
        "primitive_P": primitive_P,
    }


def toy_hilbert_check() -> dict:
    """At t=0 the no-go tuple is (Q,R1,R2)=(5,7,23)."""
    t = 0
    Q, r1, r2 = 8 * t + 5, 8 * t + 7, 8 * t + 23
    x, d = r1 * r2, 2 * Q
    finite_places = [2, Q, r1, r2]
    symbols = {str(place): h10q.hilbert(x, d, place) for place in finite_places}
    symbols["infinity"] = h10q.hilbert(x, d, h10q.OO)
    assert symbols == {"2": 1, "5": 1, "7": -1, "23": -1, "infinity": 1}
    assert math.prod(symbols.values()) == 1
    return {
        "label": "PROVED",
        "t": t,
        "Q": Q,
        "odd_factor_primes": [r1, r2],
        "x": x,
        "d": d,
        "hilbert_symbols": symbols,
        "product": 1,
    }


def write_report(rows: list[dict], summary: dict) -> None:
    by_name = {row["name"]: row for row in rows}
    no_go = by_name["synthetic_no_go"]
    positive = by_name["synthetic_positive_with_even_factor"]
    l21_named = [row for row in rows if row["kind"] == "l21a_reducible_P"]

    table = [
        "| pattern | tuple fixed divisor | period | compatible classes | sign vectors | all-plus classes |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for row in rows:
        table.append(
            f"| `{row['name']}` | {row['tuple_fixed_divisor']} | {row['signature_period']} | "
            f"{row['compatible_residues']} | `{json.dumps(row['signature_vector_counts'], sort_keys=True)}` | "
            f"{row['all_plus_residues']} |"
        )

    l21_details: list[str] = []
    for row in l21_named:
        source = row["metadata"]["l21_source_row"]
        l21_details.extend(
            [
                f"### `{row['name']}`",
                "",
                f"This replays the L21a row `(z,tau)=({source['z']},{source['tau']})`, "
                f"substitutes `{row['metadata']['class_substitution']}`, and exactly verifies the kernel octic "
                "as the product of the two displayed positive-leading primitive quartics (the overall scale is negative).",
                f"The quartics have irreducible reductions modulo {row['factor_irreducibility_certificate_primes']}; "
                f"the odd-factor tuple has exact fixed divisor `{row['tuple_fixed_divisor']}`. "
                f"Its safe sign period is `{row['signature_period']}` and its compatible-class vectors are "
                f"`{json.dumps(row['signature_vector_counts'], sort_keys=True)}`.",
                "",
            ]
        )

    lines = [
        "# L22 — factor tuples do not automatically remove the irreducibility wall",
        "",
        "**PROVED:** the factor-tuple/progression criterion below is sufficient, and exact for progression refinement.  "
        "**CONDITIONAL:** its prime values still use classical Schinzel H.  **OPEN:** no theorem here proves the required "
        "all-plus compatibility class for every `tau_dagger` cell.  **PROVED no-go in general:** Schinzel H and Hilbert "
        "reciprocity alone cannot choose the individual signs.",
        "",
        "## 1. What classical Schinzel H gives for a reducible primitive G",
        "",
        "Fix one L19 class, write `Q(t)=q1+Nt` with `N>0`, and retain all of L19's frozen-place, moving-place, "
        "real-place, scale-squareclass, positivity, and finite-exclusion hypotheses.  Drop only the assumption that "
        "the primitive normalized polynomial `G` is irreducible.  By Gauss' lemma write",
        "",
        "    G = eta * product_i g_i^{e_i},",
        "",
        "where `eta=+/-1`, the `g_i in Z[t]` are pairwise distinct, primitive, irreducible, and chosen with positive "
        "leading coefficient.  Put `I={i:e_i is odd}`.  The exact Schinzel tuple is",
        "",
        "    T = {Q} union {g_i : i in I}.",
        "",
        "Its admissibility hypothesis is not pairwise admissibility: it is the simultaneous condition",
        "",
        "    for every prime ell there is n in Z with ell not dividing Q(n) * product_{i in I} g_i(n).",
        "",
        "Equivalently, that product has fixed divisor one.  Assume also `Res(Q,g_i)!=0`; the kernel's moving-prime "
        "unit condition supplies this in the intended application.  Classical Schinzel H applied to this entire tuple "
        "then gives infinitely many positive integers `t` for which `Q(t)` and every odd-exponent `g_i(t)` are positive "
        "primes.  After deleting finitely many `t`, these primes are pairwise distinct and avoid all fixed data: equality "
        "of two values is a root of a nonzero difference polynomial, and equality with a fixed prime is finite.",
        "",
        "Even-exponent factors are intentionally absent from the tuple.  At a prime value `R_i=g_i(t)`, other odd "
        "factors contribute valuation zero (after distinctness), while every even exponent contributes an even valuation "
        "even if its value is divisible by `R_i`.  Thus the outside odd-valuation places are exactly the `R_i`, one for "
        "each `i in I`.  If `I` is empty, no nonlinear Schinzel polynomial is needed and the sign condition below is vacuous.",
        "",
        "## 2. The individual sign and the fixed resultant character",
        "",
        "Write the square class of `d(t)/Q(t)` as a fixed squarefree integer `K`; on `tau_dagger`, "
        "`d=2*alpha*eps*f*Q` and `alpha=4a^4` is a square, so `K` represents `2*eps*f`.  At an eligible odd "
        "prime `R_i` all fixed denominators are units, `v_{R_i}(x)` is odd, and `v_{R_i}(d)=0`.  The odd-prime "
        "Hilbert formula therefore gives",
        "",
        "    (x,d)_{R_i} = (d/R_i) = (K*Q(t) / R_i).                 (1)",
        "",
        "Let `m_i=deg(g_i)` and",
        "",
        "    rho_i = Res(Q,g_i) = N^{m_i} g_i(-q1/N) != 0.",
        "",
        "Modulo the prime `Q(t)`, `g_i(t) = rho_i*N^{-m_i}`.  Quadratic reciprocity turns (1) into",
        "",
        "    sigma_i(t) = (K/R_i)",
        "                 * (-1)^(((Q(t)-1)/2)*((R_i-1)/2))",
        "                 * (rho_i/Q(t)) * (N/Q(t))^{m_i}.          (2)",
        "",
        "All exponents in the reciprocity sign are read modulo two.  If `C_i` is the signed squarefree representative "
        "of `rho_i*N^{m_i}`, the last two symbols in (2) equal `(C_i/Q(t))`.  Consequently the vector is periodic in "
        "`t`; one explicit safe period is",
        "",
        "    L = lcm(4*|K|, 4*|C_1|, ..., 4*|C_r|, 4).             (3)",
        "",
        "outside the finitely many resultant/fixed-data exceptions.  This is the obstruction: once the L19 class is "
        "fixed, these characters are not freely prescribable prime-by-prime.  They are a finite, computable function of "
        "`t mod L`.  Hilbert reciprocity fixes only their product after the frozen, moving, and infinite signs are fixed; "
        "it does not turn that product constraint into the all-plus vector.",
        "",
        "## 3. Exact progression-refinement criterion",
        "",
        "**Theorem (PROVED; conclusion CONDITIONAL on classical Schinzel H).**  Assume the factor tuple `T` above is "
        "admissible.  For `r mod L`, call `r` compatible when",
        "",
        "1. no prime `p|L` divides any value `Q(r)` or `g_i(r)` with `i in I`; and",
        "2. every value of the periodic formula (2) at `r` is `+1`.",
        "",
        "If a compatible `r` exists, classical Schinzel H for",
        "",
        "    {Q(r+Lu)} union {g_i(r+Lu): i in I}",
        "",
        "gives infinitely many L19 members for which every outside odd-valuation Hilbert symbol is `+1`.  Conversely, "
        "the refinement `t=r+Lu` cannot do this if either condition fails.  Thus existence of such an `r` is the exact "
        "extra finite compatibility premise for this method, not a consequence of Schinzel H.",
        "",
        "**Proof.** Affine substitution preserves irreducibility and positive leading coefficient.  If a prime `p` "
        "divided the content or the product of the refined tuple for every `u`, then: for `p` not dividing `L`, the map "
        "`u -> r+Lu` is a bijection modulo `p`, contradicting admissibility of `T`; for `p|L`, every refined value is "
        "its value at `r` modulo `p`, contradicting condition 1.  Hence the refined polynomials are primitive and their "
        "product has fixed divisor one.  Schinzel H applies.  Formula (2), periodicity, and condition 2 give all signs "
        "`+1`; the retained L19 argument handles all other places.  If condition 1 fails, a tuple polynomial has a fixed "
        "prime divisor on that progression.  If condition 2 fails, its corresponding odd prime has sign `-1` at every "
        "eligible Schinzel value in the progression.  This proves both directions.  □",
        "",
        "A useful inheritance observation is exact: if L20's already-proved product `Q*G` has fixed divisor one, then "
        "the smaller product `Q*product_{i in I}g_i` also has fixed divisor one.  Thus in the intended classes tuple "
        "admissibility is not a new wall; the all-plus residue is.",
        "",
        "## 4. Rigorous toy no-go: the product is plus, both entries are minus",
        "",
        "Take",
        "",
        "    Q=8t+5,   g_1=8t+7=Q+2,   g_2=8t+23=Q+18,   K=2.",
        "",
        "The three primitive linear polynomials form an admissible tuple.  At `p=2` they have no root; at `p=3` their "
        "root union omits `t=0`; and for `p>=5` three linear polynomials have at most three roots, fewer than `p`.  "
        "Suppose their values are distinct odd primes.  Since `Q=5 mod 8`, both `R_i=7 mod 8`, so `(2/R_i)=+1`; "
        "also `Q=1 mod 4`, hence reciprocity gives",
        "",
        "    (Q/R_1)=(2/Q)=-1,",
        "    (Q/R_2)=(18/Q)=(2/Q)*(3/Q)^2=-1.",
        "",
        "Therefore both Hilbert signs `(2Q/R_i)` are `-1` for every simultaneous-prime value.  Their product is `+1`, "
        "so reciprocity sees no contradiction, but the desired all-plus vector never occurs and no progression can "
        "create it.  At the literal simultaneous-prime value `t=0`, `(Q,R_1,R_2)=(5,7,23)`, `x=161`, `d=10`; the "
        "checked Hilbert vector at `(2,5,7,23,infinity)` is `(+,+,-,-,+)`.  This is a theorem about the displayed tuple, "
        "not finite evidence and not a counterexample to Schinzel H.",
        "",
        "## 5. Exact mechanism replay",
        "",
        *table,
        "",
        f"The synthetic no-go has `{no_go['all_plus_residues']}` all-plus classes out of "
        f"`{no_go['compatible_residues']}` compatible classes; its only vector is "
        f"`{json.dumps(no_go['signature_vector_counts'], sort_keys=True)}`.  The positive parity-control pattern includes "
        "an even-exponent factor, correctly omits it from the Schinzel tuple, and finds "
        f"`{positive['all_plus_residues']}` all-plus classes.",
        "",
        "The L21a rows exhibit both levers exactly.  Within every fixed admissible base class the enumerated sign "
        "vector was **constant** across all compatible refined residues, an even stronger rigidity than the proved "
        "periodicity (3): progression refinement moved nothing in these instances.  Re-choosing the base residue "
        "`q1 mod 840` for the same `z=1` factor pattern realized all four vectors `++`, `+-`, `-+`, `--` "
        "(`q1=13,17,1,11`).  In the real chain the base residue is not free: the L19/L20 residue system pins `q1` "
        "modulo `4*A*prod(S)` to freeze the finite and moving symbols, exactly as the toy's `Q=5 mod 8` constraint "
        "pins its vector at `--`.  Whether the pinned system always intersects an all-plus class is precisely the "
        "open compatibility question; these finite rows demonstrate both outcomes are possible, so no finite scan "
        "can settle it.",
        "",
        *l21_details,
        "Every table entry is an exact named computation: modular irreducibility certificates, the degree-plus-one "
        "fixed-divisor identity, the resultant formula versus the direct Jacobi symbol, a full safe-period enumeration, "
        "and an exact fixed-divisor replay for every locally compatible refined tuple.  These finite mechanism tests are "
        "not used as evidence for `tau_dagger`.",
        "",
        "## 6. Scope of the route",
        "",
        "**PROVED:** whole-polynomial irreducibility is stronger than necessary.  A reducible `G` works if the odd-factor "
        "tuple has a compatible all-plus progression, and repeated even factors need no Schinzel prime value.  In that "
        "precise conditional theorem, the old premise is genuinely weakened.",
        "",
        "**OPEN uniformly:** no argument here supplies a compatible all-plus residue for the factors of every arbitrary "
        "`tau_dagger` specialization.  The toy tuple proves that classical Schinzel H plus progression refinement cannot "
        "supply such a residue in general.  Consequently this route does not remove the last non-Schinzel obligation in "
        "L19; absent a new `tau_dagger`-specific compatibility theorem, it relocates the gap from irreducibility of `G` "
        "to factorization-dependent sign compatibility.  No stronger prime-value conjecture has been introduced.",
        "",
        f"Artifact: `data/l22_factor_tuple.jsonl`.  Pacing: {summary['pacing']['iterations']} finite-loop iterations, "
        f"{summary['pacing']['sleeps']} sleeps of {PACE_SECONDS} seconds.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    started = time.perf_counter()
    rows: list[dict] = []

    no_go = analyze_pattern(
        name="synthetic_no_go",
        kind="synthetic_reciprocity_obstruction",
        Q=[5, 8],
        factors=[[7, 8], [23, 8]],
        exponents=[1, 1],
        d_over_Q=F(2),
        metadata={
            "uniform_symbol_proof": [
                "(2/(8t+7))=+1 and ((8t+5)/(8t+7))=(2/(8t+5))=-1",
                "(2/(8t+23))=+1 and ((8t+5)/(8t+23))=(18/(8t+5))=-1",
            ],
            "all_plus_impossible_for_every_simultaneous_prime_value": True,
        },
    )
    assert no_go["tuple_admissible"]
    assert no_go["signature_vector_counts"] == {"--": no_go["compatible_residues"]}
    assert no_go["all_plus_residues"] == 0
    no_go["metadata"]["literal_hilbert_replay"] = toy_hilbert_check()
    rows.append(no_go)

    positive = analyze_pattern(
        name="synthetic_positive_with_even_factor",
        kind="synthetic_parity_control",
        Q=[1, 8],
        factors=[[9, 8], [17, 8], [33, 8]],
        exponents=[1, 2, 3],
        d_over_Q=F(2),
        metadata={
            "factor_relations": ["g0=Q+8", "g1=Q+16", "g2=Q+32"],
            "even_factor_index_omitted_from_schinzel_tuple": 1,
        },
    )
    assert positive["tuple_admissible"] and positive["odd_factor_indices"] == [0, 2]
    assert positive["signature_vector_counts"] == {"++": positive["compatible_residues"]}
    rows.append(positive)

    source_rows = l21_rows()
    l21_specs = [
        (F(1), F(0), 1, 8, "l21a_z1_tau0_step8_inadmissible_control", None),
        (F(1), F(0), 1, 840, "l21a_z1_tau0_step840_q1", "-+"),
        (F(1), F(0), 11, 840, "l21a_z1_tau0_step840_q11", "--"),
        (F(1), F(0), 13, 840, "l21a_z1_tau0_step840_q13", "++"),
        (F(1), F(0), 17, 840, "l21a_z1_tau0_step840_q17", "+-"),
        (F(2), F(0), 1, 840, "l21a_z2_tau0_step840_q1", "-+"),
    ]
    l21_auxiliary = []
    for z, tau, q1, step, name, expected_vector in l21_specs:
        row, auxiliary = l21_pattern(z, tau, q1, step, source_rows, name)
        rows.append(row)
        l21_auxiliary.append(auxiliary)
        if expected_vector is None:
            assert row["tuple_fixed_divisor"] > 1 and not row["tuple_admissible"]
            assert row["signature_vector_counts"] == {}
        else:
            assert row["tuple_admissible"]
            assert row["signature_vector_counts"] == {
                expected_vector: row["compatible_residues"]
            }
    observed_z1_vectors = {
        next(iter(row["signature_vector_counts"]))
        for row in rows
        if row["kind"] == "l21a_reducible_P" and row["tuple_admissible"]
        and row["name"].startswith("l21a_z1")
    }
    assert observed_z1_vectors == {"++", "+-", "-+", "--"}

    summary = {
        "type": "summary",
        "label": "PROVED",
        "conditional_input": "classical Schinzel H for each compatible refined irreducible tuple",
        "uniform_tau_dagger_status": "OPEN",
        "generalized_implication": (
            "odd-factor tuple admissibility plus one locally nonvanishing all-plus residue modulo the explicit "
            "character period suffices under classical Schinzel H"
        ),
        "no_go_theorem": (
            "classical Schinzel H and Hilbert reciprocity alone do not force individual all-plus signs; "
            "the admissible tuple (8t+5,8t+7,8t+23) has the constant odd-factor vector (-1,-1)"
        ),
        "effect_on_irreducibility_gap": (
            "whole-G irreducibility is weakened for compatible reducible factorizations, but the uniform L19 gap "
            "is relocated to an unproved tau_dagger-specific all-plus compatibility theorem"
        ),
        "rows": len(rows),
        "synthetic_rows": sum(row["kind"].startswith("synthetic") for row in rows),
        "l21a_rows": sum(row["kind"] == "l21a_reducible_P" for row in rows),
        "l21a_admissible_rows": sum(
            row["kind"] == "l21a_reducible_P" and row["tuple_admissible"] for row in rows
        ),
        "l21a_all_plus_available_rows": sum(
            row["kind"] == "l21a_reducible_P" and row["all_plus_residues"] > 0 for row in rows
        ),
        "l21a_z1_distinct_base_residue_vectors": sorted(observed_z1_vectors),
        "l21a_vector_constant_within_each_progression": all(
            len(row["signature_vector_counts"]) == 1
            for row in rows
            if row["kind"] == "l21a_reducible_P" and row["tuple_admissible"]
        ),
        "refusals": 0,
        "refusals_used_as_evidence": False,
        "finite_tests_used_as_uniform_evidence": False,
        "pacing": {
            "iterations": PACER.iterations,
            "sleeps": PACER.sleeps,
            "every": PACE_EVERY,
            "seconds": PACE_SECONDS,
        },
        "wall_seconds": time.perf_counter() - started,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.write(json.dumps(summary, sort_keys=True) + "\n")

    write_report(rows, summary)
    print(json.dumps(summary, sort_keys=True))
    print(f"wrote {OUT} ({len(rows) + 1} rows)")
    print(f"wrote {REPORT}")
    print("VERDICT: ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
