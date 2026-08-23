#!/usr/bin/env python3
"""Close the a=1 reciprocal branch, including the cube specialization Z=z^3.

The proof is stronger than the requested cube statement in Z, but is confined
to the constructed branch a=1, A=5, s=0, tau=3/5=tau_dagger,
delta=-4/5. On that branch, for every rational Z != 0, the trace quartic and
its reciprocal degree-eight lift are irreducible. Exact finite searches and
modular scans are emitted separately as EVIDENCE.
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import h10q  # noqa: E402  (the stipulated arithmetic authority)

OUT = HERE / "data" / "l22_reciprocal_cube.jsonl"
REPORT = Path("/tmp/l22_reciprocal_cube.md")
SEARCH_HEIGHT = 500
PACE_EVERY = 200
PACE_SECONDS = 0.01
CELL_PRIME_LIMIT = 47
BRANCH = {
    "a": 1,
    "A": 5,
    "s": "0",
    "tau": "3/5=tau_dagger",
    "delta": "-4/5",
}


class Pacer:
    """Keep finite scans at a deliberately low duty cycle."""

    def __init__(self) -> None:
        self.candidates = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.candidates += 1
        if self.candidates % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def is_square_int(value: int) -> bool:
    if value < 0:
        return False
    root = math.isqrt(value)
    return root * root == value


def quadratic_character(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0:
        return 0
    return 1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1


def primes_upto(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in primes
               if prime * prime <= candidate):
            primes.append(candidate)
    return primes


def poly_mul(left: list[F], right: list[F]) -> list[F]:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            answer[i + j] += left_coefficient * right_coefficient
    return answer


def poly_scale(coefficients: list[F], scalar: F | int) -> list[F]:
    return [F(scalar) * coefficient for coefficient in coefficients]


def trace_polynomial(octic: list[F]) -> list[F]:
    """Return T with octic(b)=b^4 T(b+b^-1)."""
    assert len(octic) == 9
    assert all(octic[index] == octic[8 - index] for index in range(9))
    trace = [
        octic[4] + 2 * octic[0] - 2 * octic[2],
        octic[3] - 3 * octic[1],
        octic[2] - 4 * octic[0],
        octic[1],
        octic[0],
    ]
    lifted = [
        trace[4], trace[3], trace[2] + 4 * trace[4],
        trace[1] + 3 * trace[3],
        trace[0] + 2 * trace[2] + 6 * trace[4],
        trace[1] + 3 * trace[3], trace[2] + 4 * trace[4],
        trace[3], trace[4],
    ]
    assert lifted == octic
    return trace


# Sparse Q[Z,u] polynomials, keyed by (Z-degree, u-degree).
Sparse = dict[tuple[int, int], F]


def sparse_clean(poly: Sparse) -> Sparse:
    return {monomial: coefficient for monomial, coefficient in poly.items()
            if coefficient}


def sparse_add(left: Sparse, right: Sparse) -> Sparse:
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, F(0)) + coefficient
    return sparse_clean(answer)


def sparse_scale(poly: Sparse, scalar: F | int) -> Sparse:
    return sparse_clean({monomial: F(scalar) * coefficient
                         for monomial, coefficient in poly.items()})


def sparse_mul(left: Sparse, right: Sparse) -> Sparse:
    answer: Sparse = {}
    for (z1, u1), c1 in left.items():
        for (z2, u2), c2 in right.items():
            monomial = (z1 + z2, u1 + u2)
            answer[monomial] = answer.get(monomial, F(0)) + c1 * c2
    return sparse_clean(answer)


def sparse_pow(poly: Sparse, exponent: int) -> Sparse:
    answer: Sparse = {(0, 0): F(1)}
    base = poly
    while exponent:
        if exponent & 1:
            answer = sparse_mul(answer, base)
        base = sparse_mul(base, base)
        exponent >>= 1
    return answer


def symbolic_record() -> dict[str, Any]:
    one: Sparse = {(0, 0): F(1)}
    z: Sparse = {(1, 0): F(1)}
    u: Sparse = {(0, 1): F(1)}
    D = sparse_add(sparse_add(one, sparse_scale(z, -1)),
                   sparse_scale(sparse_pow(z, 2), -1))
    x = sparse_add(u, sparse_scale(one, -2))
    X = sparse_add(sparse_scale(one, 16),
                   sparse_scale(sparse_pow(x, 2), -5))
    trace = sparse_add(sparse_scale(sparse_pow(D, 2), 400),
                       sparse_scale(sparse_mul(sparse_pow(z, 4),
                                               sparse_pow(X, 2)), F(4, 5)))
    rational_part = sparse_mul(sparse_pow(z, 2), X)
    t_part = sparse_scale(D, 10)
    norm_product = sparse_add(sparse_pow(rational_part, 2),
                              sparse_scale(sparse_pow(t_part, 2), 5))
    assert norm_product == sparse_scale(trace, F(5, 4))

    samples = [F(-2), F(-2, 3), F(1, 2), F(1), F(8, 27), F(27, 8)]
    checked: list[dict[str, str]] = []
    H = [F(-4), F(20), F(-5)]  # 16-5(u-2)^2
    for Z in samples:
        D_value = 1 - Z - Z * Z
        assert Z and D_value
        octic = h10q._l10_P(1, Z, D_value, 5, F(-4, 5), F(0))
        trace_value = trace_polynomial(octic)
        expected = poly_scale(poly_mul(H, H), F(4, 5) * Z**4)
        expected[0] += 400 * D_value**2
        assert trace_value == expected
        checked.append({"Z": frac_text(Z), "D": frac_text(D_value)})
    return {
        "type": "symbolic_identity", "label": "PROVED",
        "D": "1-Z-Z^2",
        "trace_formula": "T(u)=400*D^2+(4/5)*Z^4*(16-5*(u-2)^2)^2",
        "factorization": "(5/4)T=(Z^2*(16-5*(u-2)^2)-10*D*t)*(Z^2*(16-5*(u-2)^2)+10*D*t), t^2=-5",
        "sparse_identity_terms": len(norm_product),
        "kernel_specializations": checked,
    }


def trace_parity_record() -> dict[str, Any]:
    residue_rows = []
    squares_mod_8 = {0, 1, 4}
    for p_mod_2, q_mod_2 in ((0, 1), (1, 0), (1, 1)):
        d_mod_2 = (q_mod_2 * q_mod_2 - p_mod_2 * q_mod_2
                   - p_mod_2 * p_mod_2) % 2
        numerator_mod_8 = (125 * d_mod_2 * d_mod_2
                           + 64 * pow(p_mod_2, 4)) % 8
        assert d_mod_2 == 1 and numerator_mod_8 == 5
        assert numerator_mod_8 not in squares_mod_8
        residue_rows.append({
            "p_mod_2": p_mod_2, "q_mod_2": q_mod_2,
            "d_mod_2": d_mod_2,
            "F1_numerator_mod_8": numerator_mod_8,
        })
    return {
        "type": "trace_reducibility", "label": "PROVED",
        "equivalent_condition": "T reducible over Q iff 80*Z^2+50*D*t is a square in Q(t), t^2=-5",
        "equivalent_rational_system": [
            "r^2-5*s^2=80*Z^2", "r*s=25*D",
        ],
        "necessary_not_asserted_sufficient_condition": "125*D^2+64*Z^4 is a rational square",
        "norm_identity": "Norm(80*Z^2+50*D*t)=100*(125*D^2+64*Z^4)",
        "uniform_obstruction": "For Z=p/q reduced, d=q^2-p*q-p^2 is odd and q^4*(125*D^2+64*Z^4)=125*d^2+64*p^4=5 mod 8",
        "primitive_parity_fibers": residue_rows,
        "conclusion": "T is irreducible for every rational Z!=0",
    }


def curve_map_record() -> dict[str, Any]:
    """Verify both curve maps as identities in two-variable polynomial rings."""
    one: Sparse = {(0, 0): F(1)}
    first: Sparse = {(1, 0): F(1)}
    second: Sparse = {(0, 1): F(1)}

    # Here first=h and second=v.
    h2 = sparse_pow(first, 2)
    curve_rhs = sparse_mul(
        sparse_add(sparse_scale(h2, 5), one),
        sparse_add(sparse_scale(h2, 5), sparse_scale(one, 16)),
    )
    curve_relation = sparse_add(sparse_pow(second, 2),
                                sparse_scale(curve_rhs, -1))
    x = sparse_scale(h2, 25)
    y = sparse_scale(sparse_mul(first, second), 25)
    elliptic_rhs = sparse_mul(
        x,
        sparse_mul(sparse_add(x, sparse_scale(one, 5)),
                   sparse_add(x, sparse_scale(one, 80))),
    )
    forward_difference = sparse_add(sparse_pow(y, 2),
                                    sparse_scale(elliptic_rhs, -1))
    assert forward_difference == sparse_scale(
        sparse_mul(h2, curve_relation), 625
    )

    # Now first=h and second=X=(v+4)/h^2.
    X2_minus_25 = sparse_add(sparse_pow(second, 2),
                             sparse_scale(one, -25))
    x_prime = sparse_add(sparse_scale(second, 8), sparse_scale(one, 85))
    y_prime = sparse_scale(sparse_mul(first, X2_minus_25), 8)
    e_prime_rhs = sparse_mul(
        x_prime,
        sparse_mul(sparse_add(x_prime, sparse_scale(one, -45)),
                   sparse_add(x_prime, sparse_scale(one, -125))),
    )
    birational_difference = sparse_add(
        sparse_pow(y_prime, 2), sparse_scale(e_prime_rhs, -1)
    )
    defining_relation = sparse_add(
        sparse_mul(X2_minus_25, h2), sparse_scale(x_prime, -1)
    )
    assert birational_difference == sparse_scale(
        sparse_mul(X2_minus_25, defining_relation), 64
    )
    return {
        "type": "curve_maps",
        "label": "PROVED",
        "two_covering_identity": "y^2-x*(x+5)*(x+80)=625*h^2*(v^2-(5*h^2+1)*(5*h^2+16))",
        "birational_identity": "y'^2-x'*(x'-45)*(x'-125)=64*(X^2-25)*((X^2-25)*h^2-(8*X+85))",
        "birational_forward": "X=(v+4)/h^2, Y=h*(X^2-25), x'=8*X+85, y'=8*Y",
        "birational_inverse": "h=8*y'/((x'-85)^2-1600), v=((x'-85)/8)*h^2-4",
    }


def primitive_cover_count(a: int, b: int, d: int, modulus: int,
                          prime: int) -> int:
    assert b % d == 0
    squares = {value * value % modulus for value in range(modulus)}
    count = 0
    for M in range(modulus):
        for e in range(modulus):
            PACER.tick()
            if M % prime == 0 and e % prime == 0:
                continue
            rhs = (d * pow(M, 4, modulus)
                   + a * M * M * e * e
                   + (b // d) * pow(e, 4, modulus)) % modulus
            if rhs in squares:
                count += 1
    return count


def descent_records() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    E_cases = {
        1: ("realized", "O"), -1: ("modular", "mod 3"),
        2: ("modular", "mod 16"), -2: ("modular", "mod 16"),
        5: ("modular", "mod 3"), -5: ("realized", "(-5,0)"),
        10: ("modular", "mod 16"), -10: ("modular", "mod 16"),
    }
    E_prime_cases = {
        1: ("realized", "O"), -1: ("modular", "mod 16"),
        3: ("modular", "mod 16"), -3: ("real", "negative definite"),
        5: ("realized", "(45,0)"), -5: ("modular", "mod 16"),
        15: ("modular", "mod 16"), -15: ("real", "negative definite"),
    }
    for curve, a, b, cases in (
        ("E: y^2=x^3+85*x^2+400*x", 85, 400, E_cases),
        ("E': y^2=x^3-170*x^2+5625*x", -170, 5625,
         E_prime_cases),
    ):
        for d, (kind, witness) in cases.items():
            record: dict[str, Any] = {
                "type": "isogeny_cover", "label": "PROVED",
                "curve": curve, "square_class": d,
                "cover": f"N^2={d}*M^4+({a})*M^2*e^2+({b // d})*e^4",
            }
            if kind == "modular":
                modulus = 3 if witness == "mod 3" else 16
                prime = 3 if modulus == 3 else 2
                count = primitive_cover_count(a, b, d, modulus, prime)
                assert count == 0
                record.update({"status": "excluded", "obstruction": witness,
                               "primitive_residue_solutions": count})
            elif kind == "real":
                assert d < 0 and b // d < 0 and a < 0
                record.update({"status": "excluded", "obstruction": witness})
            else:
                record.update({"status": "realized", "witness": witness})
            records.append(record)
    rank_record = {
        "type": "two_isogeny_descent", "label": "PROVED",
        "E": "y^2=x*(x+5)*(x+80)",
        "E_prime": "y^2=x*(x-45)*(x-125)",
        "alpha_E": [1, -5], "alpha_E_prime": [1, 5],
        "rank_formula": "2^rank=|alpha(E)|*|alpha(E')|/4",
        "rank": 0,
    }
    return records, rank_record


def curve_point_count(prime: int) -> int:
    count = 1
    for x in range(prime):
        rhs = x * (x + 5) * (x + 80) % prime
        if rhs == 0:
            count += 1
        elif quadratic_character(rhs, prime) == 1:
            count += 2
    return count


def torsion_record() -> dict[str, Any]:
    counts = {prime: curve_point_count(prime) for prime in (7, 11, 13)}
    assert counts == {7: 8, 11: 16, 13: 12}
    assert math.gcd(*counts.values()) == 4
    return {
        "type": "elliptic_torsion", "label": "PROVED",
        "good_reduction_point_counts": {str(p): n for p, n in counts.items()},
        "gcd": 4,
        "rational_points": ["O", "(0,0)", "(-5,0)", "(-80,0)"],
        "conclusion": "E(Q)=E(Q)[2]",
    }


def normalized_curve_fibers() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for prime in primes_upto(CELL_PRIME_LIMIT):
        squares = {value * value % prime for value in range(prime)}
        count = 0
        for numerator in range(prime):
            for denominator in range(prime):
                PACER.tick()
                if numerator == denominator == 0:
                    continue
                rhs = (25 * pow(numerator, 4, prime)
                       + 85 * numerator * numerator * denominator * denominator
                       + 16 * pow(denominator, 4, prime)) % prime
                if rhs in squares:
                    count += 1
        assert count > 0
        records.append({
            "type": "normalized_curve_modular_fiber", "label": "EVIDENCE",
            "prime": prime, "primitive_pairs_with_square_rhs": count,
            "verdict": "locally populated; not a uniform obstruction",
        })
    return records


def vp_integer(value: int, prime: int) -> int:
    assert value
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def hensel_square_root(value: int, prime: int, exponent: int,
                       root_mod_prime: int) -> int:
    root = root_mod_prime % prime
    modulus = prime
    assert (root * root - value) % modulus == 0 and root % prime
    for _ in range(1, exponent):
        quotient = (value - root * root) // modulus
        correction = quotient * pow(2 * root, -1, prime) % prime
        root += correction * modulus
        modulus *= prime
        assert (root * root - value) % modulus == 0
    return root % modulus


def xi_integer(Z: int) -> tuple[int, int, int, int]:
    D = 1 - Z - Z * Z
    factor_64 = 125 * D * D + 64 * Z**4
    factor_1024 = 125 * D * D + 1024 * Z**4
    return factor_64 * factor_1024, D, factor_64, factor_1024


def cell_prime_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for prime in primes_upto(CELL_PRIME_LIMIT):
        if prime == 2:
            continue
        lifted = 0
        for unit in range(1, prime):
            PACER.tick()
            z = prime * unit
            Z = z**3
            Xi, D, factor_64, _ = xi_integer(Z)
            assert D % prime
            if prime == 5:
                assert vp_integer(factor_64, prime) == 3
                assert vp_integer(Xi, prime) == 6
                unit_Xi = Xi // prime**6
                root = hensel_square_root(unit_Xi, prime, 6,
                                          (D * D) % prime)
                assert (root * root - unit_Xi) % prime**6 == 0
            else:
                expected_character = quadratic_character(5, prime)
                assert quadratic_character(factor_64, prime) == expected_character
                root = hensel_square_root(Xi, prime, 6,
                                          (125 * D * D) % prime)
                assert (root * root - Xi) % prime**6 == 0
            lifted += 1
        records.append({
            "type": "cell_prime_scan", "label": "EVIDENCE", "w": prime,
            "unit_residues_checked_for_z_over_w": lifted,
            "Xi_square_root_lifted_mod_w_power": 6,
            "trace_F1_local_square_class": (
                "5 (valuation 3)" if prime == 5
                else f"5, Legendre={quadratic_character(5, prime)}"
            ),
            "a_equals_1_L20_admissible": (
                prime != 5 and quadratic_character(5, prime) == -1),
            "verdict": "Xi is a w-adic square; w supplies no lift obstruction",
        })
    return records


def bounded_search() -> dict[str, Any]:
    checked = 0
    trace_norm_square_hits: list[list[int]] = []
    xi_square_hits: list[list[int]] = []
    for denominator in range(1, SEARCH_HEIGHT + 1):
        for numerator in range(-SEARCH_HEIGHT, SEARCH_HEIGHT + 1):
            if numerator == 0 or math.gcd(abs(numerator), denominator) != 1:
                continue
            PACER.tick()
            d = (denominator**6 - numerator**3 * denominator**3
                 - numerator**6)
            assert d
            z_fourth_numerator = numerator**12
            factor_64 = 125 * d * d + 64 * z_fourth_numerator
            factor_1024 = 125 * d * d + 1024 * z_fourth_numerator
            assert factor_64 % 8 == 5
            checked += 1
            if is_square_int(factor_64):
                trace_norm_square_hits.append([numerator, denominator])
            if is_square_int(factor_64 * factor_1024):
                xi_square_hits.append([numerator, denominator])
    assert checked == 304_462
    assert not trace_norm_square_hits and not xi_square_hits
    return {
        "type": "bounded_rational_search", "label": "EVIDENCE",
        "variable": "z=p/q in lowest terms; Z=z^3",
        "box": {"numerator": f"0<|p|<={SEARCH_HEIGHT}",
                "denominator": f"1<=q<={SEARCH_HEIGHT}",
                "gcd": "gcd(|p|,q)=1"},
        "inputs_checked": checked,
        "trace_norm_square_hits": trace_norm_square_hits,
        "Xi_square_hits": xi_square_hits, "counterexamples": 0,
        "role": "reproducible evidence only; the uniform proof is independent",
    }


def build_report(payload: list[dict[str, Any]], elapsed: float) -> str:
    search = next(row for row in payload if row["type"] == "bounded_rational_search")
    fiber_rows = [row for row in payload
                  if row["type"] == "normalized_curve_modular_fiber"]
    cell_rows = [row for row in payload if row["type"] == "cell_prime_scan"]
    cover_rows = [row for row in payload if row["type"] == "isogeny_cover"]
    lines = [
        "# L22 reciprocal-cube closure", "", "## Status", "",
        "**PROVED.** Fix the constructed reciprocal branch",
        "`a=1, A=5, s=0, tau=3/5=tau_dagger, delta=-4/5`. Let `Z` be any",
        "nonzero rational number and put `D=1-Z-Z^2`. The normalized",
        "degree-eight polynomial evaluated on this branch is irreducible over",
        "`Q`. In particular this holds for every nonzero rational `z` after",
        "`Z=z^3`. There are no residual curves on this branch.", "",
        "The branch qualification is indispensable. No statement here applies",
        "to another `tau` or `delta`: the L21a square-delta reducibility lemma",
        "gives reducible loci outside this branch, including `tau=0` and",
        "`tau=1/3`. The theorem is stronger only in its rational-`Z` scope.",
        "Every finite scan below is labelled **EVIDENCE** and is not used to",
        "infer the theorem.", "",
        "## 1. Reciprocal trace and exact reducibility equations (PROVED)", "",
        "For `u=b+b^-1`, direct expansion against `h10q._l10_P` gives", "",
        "```text", "P(b)=b^4 T(u),",
        "T(u)=400D^2+(4/5)Z^4(16-5(u-2)^2)^2.", "```", "",
        "In `E=Q(t)`, `t^2=-5`, there is the exact factorization", "",
        "```text", "(5/4)T = Q_-(u) Q_+(u),",
        "Q_±(u)=Z^2(16-5(u-2)^2) ± 10Dt.", "```", "",
        "Because `D!=0`, the conjugate quadratics are distinct and coprime.",
        "If `Q_+` is irreducible over `E`, any factor of `T` over `Q` becomes",
        "over `E` a product of irreducible factors chosen from `{Q_+,Q_-}`.",
        "Galois conjugation swaps those two factors, so no proper nonempty",
        "choice descends to `Q`. Conversely, if `Q_+` splits over `E`, the",
        "norm of either linear factor is a proper quadratic factor over `Q`.",
        "Thus the following are **equivalent**:", "",
        "1. `T` is reducible over `Q`;",
        "2. `Q_+` (equivalently `Q_-`) is reducible over `E`;",
        "3. `C=80Z^2+50Dt` is a square in `E`;",
        "4. there are `r,s in Q` satisfying",
        "   `r^2-5s^2=80Z^2` and `rs=25D`.", "",
        "The discriminant calculation behind (3) is", "", "```text",
        "disc(Q_+)=4Z^2(80Z^2+50Dt).", "```", "",
        "Taking the norm in (3) gives the **necessary** condition", "",
        "```text", "125D^2+64Z^4 is a rational square,",
        "N_E/Q(C)=100(125D^2+64Z^4).", "```", "",
        "This norm condition is not asserted sufficient. Write `Z=p/q` in",
        "lowest terms and `d=q^2-pq-p^2`. The three primitive parity classes",
        "all make `d` odd, so", "", "```text",
        "q^4(125D^2+64Z^4)=125d^2+64p^4 = 5 (mod 8),", "```", "",
        "which cannot be an integer square. Therefore `T` is irreducible for",
        "every rational `Z!=0`. This replaces the earlier `(5|w)=-1` local",
        "argument by a uniform 2-adic proof.", "",
        "## 2. Reciprocal lift and the norm direction (PROVED)", "",
        "Let `u` be a root of the irreducible quartic `T` and `K=Q(u)`.",
        "Since `b` satisfies `b^2-ub+1=0`, the following are **equivalent**:",
        "", "```text",
        "P irreducible over Q  <=>  u^2-4 is not a square in K.", "```", "",
        "The exact norm is", "", "```text",
        "N_K/Q(u^2-4)=T(2)T(-2)/lc(T)^2",
        "             =(4/(25Z^4))^2 Xi(Z),",
        "Xi(Z)=(125D^2+64Z^4)(125D^2+1024Z^4).", "```", "",
        "Thus `Xi(Z)` nonsquare is **sufficient** for irreducibility of the",
        "lift. No converse is used or claimed: a square norm need not make an",
        "element a square in `K`.", "",
        "## 3. Xi is never a rational square (PROVED)", "",
        "Put `h=5D/(8Z^2)`. Since `64Z^4` is a rational square, `Xi(Z)`",
        "is a square exactly when the explicit genus-one curve", "", "```text",
        "C: v^2=(5h^2+1)(5h^2+16)", "```", "",
        "has the corresponding rational point. It has the following exact",
        "birational Weierstrass model. For `h!=0`, put", "",
        "```text",
        "X=(v+4)/h^2,       Y=h(X^2-25),",
        "x'=8X+85,          y'=8Y.",
        "```", "",
        "Then `y'^2=x'(x'-45)(x'-125)`, defining `E'`. Away from the",
        "finitely many projective exceptional points, the inverse is", "",
        "```text",
        "h=8y'/((x'-85)^2-1600),",
        "v=((x'-85)/8)h^2-4.",
        "```", "",
        "The identity follows directly from",
        "`(X^2-25)h^2=8X+85`. There is also the particularly simple",
        "2-covering map to the 2-isogenous curve", "",
        "```text",
        "C -> E: y^2=x(x+5)(x+80),",
        "x=25h^2, y=25hv.",
        "```", "",
        "Only the latter forward implication is needed: an `Xi` square",
        "produces an affine point of `E` with positive square `x` when",
        "`h!=0`. The birational formulas record the reverse scope exactly.", "",
        "A complete 2-isogeny descent proves `E(Q)` has rank zero. For",
        "`E_{a,b}: y^2=x^3+a x^2+b x`, the descent map has square classes",
        "represented by squarefree divisors `d|b`; its covering is", "",
        "```text", "N^2=dM^4+aM^2e^2+(b/d)e^4, gcd(M,e)=1.", "```", "",
        "For `E`, `b=400`, so the complete candidate list is",
        "`d in {±1,±2,±5,±10}`. The 2-isogenous curve is",
        "`E': y^2=x(x-45)(x-125)`; there `b'=5625`, giving exactly",
        "`d in {±1,±3,±5,±15}`. The complete square-class table is:", "",
        "| curve | d | result | obstruction/witness |",
        "|---|---:|---|---|",
    ]
    for row in cover_rows:
        curve = "E'" if row["curve"].startswith("E'") else "E"
        detail = row.get("obstruction", row.get("witness", ""))
        lines.append(f"| {curve} | {row['square_class']} | {row['status']} | {detail} |")
    lines.extend([
        "", "Here are the congruence exclusions without a black-box rank call.",
        "For `E` and `d=-1,5`, the right side of the covering is `2 mod 3`",
        "for each of the three primitive parity patterns modulo 3. For even",
        "`d`, its value modulo 16 lies in",
        "`{2,3,6,7,8,10,11,12,14,15}`; none is in the square set",
        "`{0,1,4,9}`. For `E'` and `d=-1,3,-5,15`, the analogous primitive",
        "values modulo 16 lie in `{3,7,11,12,15}`. Finally, for",
        "`d=-3,-15` on `E'`, all three coefficients are negative, so the",
        "form is negative definite over `R`. Thus no excluded class has a",
        "primitive rational covering point.", "",
        "Consequently `alpha(E)={1,-5}` and `alpha(E')={1,5}`. The formula",
        "", "```text", "2^rank = |alpha(E)| |alpha(E')| / 4", "```", "",
        "gives `rank E(Q)=0`. At good primes 7, 11, and 13 the point counts",
        "are 8, 16, and 12; their gcd is 4. The four visible 2-torsion points",
        "therefore exhaust the group:", "", "```text",
        "E(Q)={O,(0,0),(-5,0),(-80,0)}.", "```", "",
        "If `D,Z` are nonzero, then `h!=0`, while the image point has",
        "`x=25h^2>0`. No point in this group has positive `x`. Therefore",
        "`Xi(Z)` is not a rational square.", "",
        "For `z=p/q`, the requested cube square curve is", "", "```text",
        "Y^2=(125d^2+64p^12)(125d^2+1024p^12),",
        "d=q^6-p^3q^3-p^6,", "```", "",
        "and it pulls back to `C` under `h=5d/(8p^6)`. The descent leaves no",
        "residual curve.", "",
        "## 4. What the cell prime w can and cannot do (PROVED)", "",
        "If `v_w(z)>=1`, then `v_w(Z)>=3` and `D` is a `w`-adic unit. In",
        "`Q_w`, for every odd `w`,", "", "```text", "Xi=(125D^2)^2",
        " *(1+64Z^4/(125D^2))*(1+1024Z^4/(125D^2)).", "```", "",
        "For `w!=5` both parenthesized units lie in `1+w Z_w`; for `w=5`",
        "they lie in `1+5^9 Z_5`. Such units are squares at an odd prime.",
        "Hence `Xi` is always a square in `Q_w`: the cell prime supplies no",
        "lift obstruction. This is a proved negative result, not a refusal.", "",
        "For the trace norm, when `w!=5` the first factor has local square",
        "class `5`. Thus `(5|w)=-1` recovers the old local trace obstruction,",
        "while mod 8 also covers `(5|w)=+1`. At `w=5` it has valuation 3.",
        "The `a=1` L20 class is available exactly when `w!=5` and",
        "`(5|w)=-1`.", "",
        "## 5. Reproducible finite evidence (EVIDENCE)", "",
        f"The exact search used `0<|p|<={SEARCH_HEIGHT}`, `1<=q<={SEARCH_HEIGHT}`,",
        f"and `gcd(|p|,q)=1`, with `z=p/q`, `Z=z^3`. It checked",
        f"**{search['inputs_checked']}** inputs and found zero trace-norm square",
        "cases, zero `Xi` square cases, and zero counterexamples. These counts",
        "are evidence only.", "",
        f"The normalized curve was scanned over all {len(fiber_rows)} prime",
        f"fibers through {CELL_PRIME_LIMIT}; every fiber was populated, so no",
        "single-prime obstruction was promoted. The successful congruence",
        "covering is the complete isogeny table (mod 3, mod 16, real place).",
        "", f"For {len(cell_rows)} odd cell primes through {CELL_PRIME_LIMIT},",
        "every nonzero residue `z/w` was checked and a square root of `Xi` was",
        "Hensel-lifted through `w^6`, corroborating the proved absence of a",
        "cell-prime lift obstruction.", "", "## Final conclusion", "",
        "**PROVED:** on the fixed branch",
        "`a=1, A=5, s=0, tau=3/5=tau_dagger, delta=-4/5`, the normalized",
        "degree-eight polynomial at `Z=z^3` is irreducible over `Q` for every",
        "rational `z!=0` with `D=1-z^3-z^6!=0`. In fact the same branch",
        "polynomial is irreducible for every rational `Z!=0`. The trace",
        "criterion is an equivalence; the norm-nonsquare tests are used only",
        "in their valid sufficient directions. No residual curve remains.",
        "No claim is made for the L21a square-delta branches.", "",
        f"Replay wall-clock: {elapsed:.3f} s. Pacing: {PACER.candidates}",
        f"finite-loop iterations, {PACER.sleeps} sleeps of {PACE_SECONDS} s",
        f"every {PACE_EVERY} iterations.",
    ])
    return "\n".join(lines) + "\n"


def main() -> None:
    started = time.perf_counter()
    symbolic = symbolic_record()
    trace = trace_parity_record()
    maps = curve_map_record()
    covers, rank = descent_records()
    torsion = torsion_record()
    fibers = normalized_curve_fibers()
    cell_scans = cell_prime_records()
    search = bounded_search()
    assert rank["rank"] == 0
    assert torsion["rational_points"] == [
        "O", "(0,0)", "(-5,0)", "(-80,0)"
    ]
    elapsed = time.perf_counter() - started
    lift = {
        "type": "reciprocal_lift", "label": "PROVED",
        "equivalence_after_trace_irreducibility": "P irreducible iff u^2-4 is nonsquare in Q(u)",
        "norm_identity": "Norm(u^2-4)=(4/(25*Z^4))^2*Xi(Z)",
        "Xi": "(125*D^2+64*Z^4)*(125*D^2+1024*Z^4)",
        "logical_direction": "Xi nonsquare is sufficient; converse not claimed",
    }
    xi_theorem = {
        "type": "Xi_theorem", "label": "PROVED",
        "normalization": "h=5*D/(8*Z^2), v^2=(5*h^2+1)*(5*h^2+16)",
        "birational_model": "X=(v+4)/h^2, Y=h*(X^2-25), x'=8*X+85, y'=8*Y on E': y'^2=x'*(x'-45)*(x'-125)",
        "birational_inverse": "h=8*y'/((x'-85)^2-1600), v=((x'-85)/8)*h^2-4",
        "two_covering_map": "x=25*h^2, y=25*h*v on E: y^2=x*(x+5)*(x+80)",
        "elliptic_group": ["O", "(0,0)", "(-5,0)", "(-80,0)"],
        "conclusion": "Xi(Z) is nonsquare for every rational Z,D with Z*D!=0",
    }
    cell_theorem = {
        "type": "cell_prime_theorem", "label": "PROVED",
        "hypothesis": "w odd prime, v_w(z)>=1, Z=z^3",
        "Xi_local_behavior": "Xi is a square in Q_w; no cell-prime lift obstruction",
        "trace_local_behavior": "F1 has square class 5 for w!=5 and valuation 3 for w=5",
        "a_equals_1_admissibility": "w!=5 and Legendre(5,w)=-1",
    }
    summary = {
        "type": "summary", "label": "PROVED", "uniform_verdict": "PROVED",
        "branch": BRANCH,
        "theorem": "On a=1,A=5,s=0,tau=3/5=tau_dagger,delta=-4/5, the normalized P is irreducible over Q for every rational Z!=0 with D=1-Z-Z^2!=0",
        "cube_corollary": "On a=1,A=5,s=0,tau=3/5=tau_dagger,delta=-4/5, the normalized P at Z=z^3 is irreducible for every rational z!=0 with D!=0",
        "branch_scope_note": "No claim for other tau or delta; L21a proves square-delta reducibility outside this branch, including tau=0 and tau=1/3",
        "residual_curves": [], "counterexamples": 0,
        "finite_search_inputs": search["inputs_checked"],
        "wall_seconds": elapsed,
        "pacing": {"iterations": PACER.candidates, "every": PACE_EVERY,
                   "sleep_seconds": PACE_SECONDS, "sleeps": PACER.sleeps},
    }
    meta = {
        "type": "meta", "label": "PROVED", "schema": "l22-reciprocal-cube-v1",
        "source_script": "math/h10q/l22_reciprocal_cube.py",
        "kernel_authority": "math/h10q/h10q.py",
        "refusals_are_evidence": False,
        "finite_scans_imply_uniform_theorem": False,
    }
    payload = [meta, symbolic, trace, lift, xi_theorem, maps, rank, torsion,
               cell_theorem, *covers, *fibers, *cell_scans, search, summary]
    for row in payload:
        row.setdefault("branch", BRANCH)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(
        json.dumps(row, sort_keys=True, separators=(",", ":"))
        for row in payload) + "\n", encoding="utf-8")
    REPORT.write_text(build_report(payload, elapsed), encoding="utf-8")
    print(f"L22 reciprocal cube: PROVED on tau_dagger=3/5, delta=-4/5 "
          f"for all rational Z!=0; search={search['inputs_checked']}, "
          f"counterexamples=0; {elapsed:.3f}s")


if __name__ == "__main__":
    main()
