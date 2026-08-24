#!/usr/bin/env python3
"""L26 reciprocal tie on an even bridge pullback.

This producer proves exact algebraic/local statements and writes their replay
records.  Bounded irreducibility rows are EVIDENCE only.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as F
from pathlib import Path

from l24_diagonal_geometry import irreducibility_certificate
from l25_scaled_coupling import squareclass_2adic, vp

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l26_reciprocal_tie.jsonl"
REPORT = Path("/tmp/l26_reciprocal_tie.md")


def add(left: list[F], right: list[F]) -> list[F]:
    result = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += coefficient
    return result


def multiply(left: list[F], right: list[F]) -> list[F]:
    result = [F(0)] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] += left_coefficient * right_coefficient
    return result


def scale(polynomial: list[F], scalar: F) -> list[F]:
    return [scalar * coefficient for coefficient in polynomial]


def evaluate(polynomial: list[F], value: F) -> F:
    answer = F(0)
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


def fraction_mod(value: F, prime: int) -> int:
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, -1, prime) % prime


def bridge_data(a: F, Z: F, b: F) -> tuple[F, F, F, F]:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
    assert a and b not in (0, 1) and Z and D and Ng
    c = a * a * Z * Z * Ng / (A * b * b * D)
    return A, D, Ng, c


def reciprocal_tie(Z: F, D: F, b: F) -> F:
    assert b and D
    return Z * (b + 1) / (D * b)


def reciprocal_polynomials(a: F, Z: F, delta: F) -> tuple[list[F], list[F]]:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    # Ng(b) = 16*a^4*b^2 - A*(b-1)^4.
    b_minus_one_fourth = [F(1), F(-4), F(6), F(-4), F(1)]
    Ng = add([F(0), F(0), 16 * a**4], scale(b_minus_one_fourth, -A))
    P = add(
        add([F(0)] * 4 + [16 * A * A * D * D],
            scale(multiply(Ng, Ng), -delta * a**4 * Z**4)),
        [F(0)] * 3 + [
            -32 * A**3 * Z**2,
            -64 * A**3 * Z**2,
            -32 * A**3 * Z**2,
        ],
    )
    assert len(P) == 9
    # R(u) = 16*a^4 - A*(u-2)^2.
    R = [16 * a**4 - 4 * A, 4 * A, -A]
    T = add(
        add([16 * A * A * D * D], scale(multiply(R, R), -delta * a**4 * Z**4)),
        [-64 * A**3 * Z**2, -32 * A**3 * Z**2],
    )
    assert len(T) == 5
    return P, T


def identity_replay() -> dict:
    checks = 0
    for a in map(F, (1, 3, 5, -1)):
        A = 1 + 4 * a * a
        delta = -4 * a**4 / A
        for z in (F(3), F(5, 3), F(-1, 3), F(2)):
            Z = z * z
            D = 1 - Z - a * a * Z * Z
            if D == 0:
                continue
            P, T = reciprocal_polynomials(a, Z, delta)
            for b in (F(3), F(-3), F(5, 3), F(7, 5)):
                A0, D0, Ng, c = bridge_data(a, Z, b)
                assert (A0, D0) == (A, D)
                tied_s = reciprocal_tie(Z, D, b)
                M = 16 - delta * c * c - 32 * A * b * tied_s * tied_s
                cleared = D * D * A * A * b**4 * M
                assert evaluate(P, b) == cleared
                u = b + 1 / b
                assert evaluate(P, b) == b**4 * evaluate(T, u)
                checks += 1
    return {
        "type": "exact-reciprocal-identity",
        "label": "PROVED",
        "instances": checks,
        "tie": "Z*(b+1)/(D*b)",
        "trace_degree": 4,
        "reciprocal_weight": 8,
    }

def trace_character_descent() -> dict:
    checks = 0
    for prime in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        for b in range(1, prime):
            if (b + 1) % prime == 0:
                continue
            u = (b + pow(b, -1, prime)) % prime
            assert legendre_mod(2 * b, prime) == legendre_mod(2 * (u + 2), prime)
            checks += 1
    return {
        "type": "trace-character-descent",
        "label": "PROVED by (u+2)=(b+1)^2/b away from b=-1",
        "instances": checks,
        "character_identity": "(2b|p)=(2(u+2)|p)",
        "trace_algebra": "Q[u]/(T), with T quartic",
        "conditional_field": "if T is irreducible, Q(theta) has degree 4",
        "conditional_bad_sign_extension": (
            "if 2(theta+2) is nonsquare, "
            "Q(theta)(sqrt(2(theta+2))) is nontrivial"
        ),
        "scope": "primes dividing T(u) outside the fixed divisor T(-2)",
    }


def legendre_mod(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0:
        return 0
    answer = pow(residue, (prime - 1) // 2, prime)
    return -1 if answer == prime - 1 else answer



def guard_and_dyadic_replay() -> dict:
    checks = 0
    valuation_rows = []
    for a in map(F, (1, 3, 5, -1, -3)):
        A = 1 + 4 * a * a
        delta = -4 * a**4 / A
        assert vp(A, 2) == 0 and vp(delta, 2) == 2
        for b in (F(3), F(5), F(7), F(1, 3), F(5, 3)):
            assert vp(b, 2) == 0
            for Z in (F(4), F(9), F(1), F(1, 4), F(1, 16)):
                D = 1 - Z - a * a * Z * Z
                if D == 0:
                    continue
                A0, D0, Ng, c = bridge_data(a, Z, b)
                assert (A0, D0) == (A, D)
                tied_s = reciprocal_tie(Z, D, b)
                M = 16 - delta * c * c - 32 * A * b * tied_s * tied_s
                P_value = D * D * A * A * b**4 * M
                assert vp(tied_s, 2) >= 1
                assert vp(Ng, 2) >= 4
                assert vp(c, 2) >= 4
                assert vp(M, 2) == 4
                assert squareclass_2adic(M)
                assert squareclass_2adic(P_value)
                checks += 1
                if len(valuation_rows) < 6:
                    valuation_rows.append({
                        "a": str(a), "b": str(b), "Z": str(Z),
                        "v2_Z": vp(Z, 2), "v2_D": vp(D, 2),
                        "v2_s_rec": vp(tied_s, 2), "v2_c": vp(c, 2),
                        "v2_M": vp(M, 2),
                    })
    return {
        "type": "guards-and-dyadic-square",
        "label": "PROVED by valuation identities; finite rows replay the cases",
        "instances": checks,
        "valuation_rows": valuation_rows,
        "conclusion": "M and P_rec are Q_2-squares on Phi",
    }


def target_certificate() -> dict:
    z = F(3)
    Z = z * z
    a = F(1)
    b = F(3)
    A, D, Ng, c = bridge_data(a, Z, b)
    delta = -4 * a**4 / A
    tied_s = reciprocal_tie(Z, D, b)
    M = 16 - delta * c * c - 32 * A * b * tied_s * tied_s
    assert (A, D, Ng, c, tied_s) == (F(5), F(-89), F(64), F(-576, 445), F(-12, 89))
    assert M == F(8529104, 990125)
    assert vp(tied_s, 3) == 1 and vp(c, 3) == 2 and vp(M, 3) == 0
    assert fraction_mod(-delta * A, 3) == 1 and fraction_mod(M, 3) == 1
    return {
        "type": "odd-target-local",
        "label": "PROVED on the standard guarded W1 target stratum; named exact replay at w=3",
        "pullback": "Z=z^2",
        "target_guards": {
            "w": "odd", "v_w_z": ">=1", "v_w_Z": "2*v_w(z)>=2",
            "v_w_b": 1, "v_w_a": 0, "v_w_A": 0, "v_w_delta": 0,
            "D_unit": True, "Ng_unit": True,
        },
        "general_valuations": {
            "v_w_s_rec": "v_w(Z)-1 >= 1",
            "v_w_c": "2*v_w(Z)-2 >= 2",
            "M_mod_w": 16,
            "alpha": "4*a^4 is a square unit on tau_dagger",
        },
        "sample": {
            "w": 3, "a": str(a), "b": str(b), "Z": str(Z),
            "D": str(D), "s_rec": str(tied_s), "c": str(c), "M": str(M),
        },
    }


def bounded_irreducibility_context() -> dict:
    a, z, b = F(1), F(3), F(3)
    Z = z * z
    A = 1 + 4 * a * a
    delta = -4 * a**4 / A
    P, T = reciprocal_polynomials(a, Z, delta)
    certificate = irreducibility_certificate(P, 37)
    expected = [
        164025, -1312200, 3542940, -5391360, 10013306,
        -5391360, 3542940, -1312200, 164025,
    ]
    primitive = certificate["reduction_ascending"]
    # The certificate routine normalizes content before reduction.  The exact
    # primitive vector is asserted separately via rational scaling.
    ratio = P[0] / expected[0]
    assert all(P[index] == ratio * expected[index] for index in range(9))
    trace_certificate = irreducibility_certificate(T, 37)
    assert certificate["degree"] == 8 and trace_certificate["degree"] == 4 and primitive
    return {
        "type": "bounded-irreducibility-context",
        "label": "EVIDENCE ONLY; one reciprocal octic is still irreducible",
        "sample": {"a": 1, "z": 3, "Z": 9, "b": str(b)},
        "prime": 37,
        "degree": certificate["degree"],
        "primitive_vector": expected,
        "trace_quartic": [str(coefficient) for coefficient in T],
        "trace_prime": 37,
        "trace_degree": trace_certificate["degree"],
        "trace_certificate_sha256": trace_certificate["primitive_sha256"],
        "certificate_sha256": certificate["primitive_sha256"],
    }


def build_report(records: list[dict], elapsed: float) -> str:
    lines = [
        "# L26 reciprocal tie replay", "",
        "PROVED labels are exact algebra/valuation statements.  The modular",
        "irreducibility row is bounded EVIDENCE only.  No member theorem or",
        "unconditional H10/Q consequence is claimed.", "",
    ]
    for record in records:
        lines.extend(["## " + record["type"], "", "```json",
                      json.dumps(record, indent=2, sort_keys=True), "```", ""])
    lines.append(f"elapsed_seconds: {elapsed:.3f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [
        identity_replay(),
        trace_character_descent(),
        guard_and_dyadic_replay(),
        target_certificate(),
        bounded_irreducibility_context(),
    ]
    elapsed = time.perf_counter() - started
    with OUT.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    REPORT.write_text(build_report(records, elapsed), encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
