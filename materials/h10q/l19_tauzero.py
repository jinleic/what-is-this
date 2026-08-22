#!/usr/bin/env python3
"""Prove and audit the tau=0 fixed-factor escape from the L10b wall.

The new algebra is factor-by-factor and therefore covers composite
A = 1 + 4*a^2; no h10q kernel certificate is modified or generalized here.
The bounded audit independently predicts every symbol returned by
l13_filter.aligned for the stated eligible cases.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proved arithmetic kernel is authoritative)
import l13_filter  # noqa: E402

OUT = ROOT / "data" / "l19_tauzero.jsonl"
REPORT = Path("/tmp/l19_tauzero.md")
SOURCE = "l19_tauzero.py"
A_VALUES = tuple(range(1, 26, 2))
W_MIN = 101
W_MAX = 400
Q1_MAX = 200
PACE_EVERY = 200
PACE_SECONDS = 0.01
SMALL_CONTROLLED = (3, 5, 7)

HYPOTHESIS = (
    "H_tau0: a is a nonzero odd integer; A=1+4*a^2; z=-w and f=w; "
    "w and q1 are distinct odd proven primes; w=3 (mod 4); "
    "gcd(A,2*z*D)=1 for D=1-z^3-a^2*z^6; "
    "q1 is outside {3,5,7} and gcd(q1,2*A*a*z*D)=1; "
    "tau=0, delta=1, alpha=-A, eps in {+1,-1}, and b=eps*f*q1."
)


class Pacer:
    def __init__(self) -> None:
        self.candidates = 0

    def tick(self) -> None:
        self.candidates += 1
        if self.candidates % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)


def prime_character(value: int, prime: int) -> int:
    """Legendre character, with zero represented explicitly."""
    residue = value % prime
    if residue == 0:
        return 0
    return 1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1


@lru_cache(maxsize=None)
def proven_factorization(value: int) -> tuple[tuple[int, int], ...]:
    factors = h10q.factorint(abs(value))
    assert math.prod(prime**exponent for prime, exponent in factors.items()) == abs(value)
    assert all(h10q._is_prime(prime) for prime in factors)
    return tuple(sorted((int(prime), int(exponent)) for prime, exponent in factors.items()))


def factor_dict(value: int) -> dict[int, int]:
    return dict(proven_factorization(value))


def factored_jacobi(value: int, factors: dict[int, int]) -> int:
    """Jacobi(value | product p^e), using kernel-proved prime factors."""
    result = 1
    for prime, exponent in factors.items():
        character = prime_character(value, prime)
        if character == 0:
            return 0
        if exponent % 2:
            result *= character
    return result


def valuation(value: F, prime: int) -> int:
    value = F(value)
    assert value
    numerator = value.numerator
    denominator = value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def unit_residue(value: F, prime: int, known_valuation: int | None = None) -> int:
    value = F(value)
    exponent = valuation(value, prime) if known_valuation is None else known_valuation
    unit = value / F(prime) ** exponent
    return (unit.numerator % prime) * pow(unit.denominator % prime, -1, prime) % prime


def direct_hilbert(value_x: F, value_d: F, place: int | str) -> int:
    """Independent literal Hilbert formula used only as the audit predictor."""
    value_x = F(value_x)
    value_d = F(value_d)
    if place == "oo":
        return -1 if value_x < 0 and value_d < 0 else 1

    prime = int(place)
    alpha = valuation(value_x, prime)
    beta = valuation(value_d, prime)
    unit_x = unit_residue(value_x, prime, alpha)
    unit_d = unit_residue(value_d, prime, beta)
    if prime == 2:
        epsilon = lambda residue: ((residue - 1) // 2) % 2
        omega = lambda residue: ((residue * residue - 1) // 8) % 2
        parity = (
            epsilon(unit_x) * epsilon(unit_d)
            + alpha * omega(unit_d)
            + beta * omega(unit_x)
        ) % 2
        return -1 if parity else 1

    parity = (alpha * beta * ((prime - 1) // 2)) % 2
    if beta % 2 and prime_character(unit_x, prime) == -1:
        parity ^= 1
    if alpha % 2 and prime_character(unit_d, prime) == -1:
        parity ^= 1
    return -1 if parity else 1


def symbol_map(symbols: dict) -> dict[str, int]:
    return {str(place): int(value) for place, value in symbols.items()}


def base_failure(A: int, w: int, D: int) -> str | None:
    if math.gcd(A, 2 * w * abs(D)) != 1:
        return "gcd(A,2*z*D)!=1"
    return None


def q1_failure(a: int, A: int, w: int, D: int, q1: int) -> str | None:
    if q1 == 2:
        return "q1-even"
    if q1 in SMALL_CONTROLLED:
        return "q1-in-{3,5,7}"
    if q1 == w:
        return "q1=w"
    if any(value % q1 == 0 for value in (A, a, w, D)):
        return "q1-not-unit-in-A*a*z*D"
    return None


def arithmetic_data(a: int, w: int, eps: int, q1: int) -> dict:
    z = F(-w)
    Z = z**3
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    assert D.denominator == 1
    b = F(eps * w * q1)
    c = h10q._sun_h(F(a), b, Z)
    assert c is not None
    s = F(a - 1, 2)
    M = F(16) - c * c - 32 * A * b * s * s
    assert M
    return {
        "a": a,
        "A": A,
        "z": z,
        "Z": Z,
        "D": int(D),
        "b": b,
        "c": c,
        "M": M,
        "x": -A * M,
        "d": -2 * A * b,
    }


def specialized_prediction(data: dict, eps: int, w: int, q1: int) -> tuple[dict, list[str]]:
    """Predict the complete l13_filter.aligned map without calling hilbert()."""
    a = data["a"]
    A = data["A"]
    b = data["b"]
    c = data["c"]
    M = data["M"]
    x = data["x"]
    d = data["d"]
    factors = factor_dict(A)
    failures: list[str] = []
    predicted: dict[int | str, int] = {}

    # The factor-by-factor extension of the L10b A-place computation.
    for prime, exponent in factors.items():
        A_prime = A // prime**exponent
        expected_unit_x = (
            (A_prime % prime)
            * pow(unit_residue(F(prime) ** exponent * c, prime), 2, prime)
        ) % prime
        checks = {
            "v_p(c)=-e": valuation(c, prime) == -exponent,
            "v_p(M)=-2e": valuation(M, prime) == -2 * exponent,
            "v_p(x)=-e": valuation(x, prime) == -exponent,
            "v_p(d)=e": valuation(d, prime) == exponent,
            "u_x=A'*(p^e*c)^2": unit_residue(x, prime) == expected_unit_x,
            "u_d=-2*A'*b": unit_residue(d, prime)
            == (-2 * A_prime * b.numerator * pow(b.denominator, -1, prime)) % prime,
        }
        failures.extend(
            f"A-place-{prime}:{name}" for name, passed in checks.items() if not passed
        )
        predicted[prime] = (
            1 if exponent % 2 == 0 else prime_character(-2 * eps * w * q1, prime)
        )

    # At f=w, z=-w changes c from a pole to a fourth-order zero.
    w_checks = {
        "v_w(c)=4": valuation(c, w) == 4,
        "v_w(M)=0": valuation(M, w) == 0,
        "v_w(x)=0": valuation(x, w) == 0,
        "v_w(d)=1": valuation(d, w) == 1,
        "u_x=-16*A": unit_residue(x, w) == (-16 * A) % w,
    }
    failures.extend(f"f-place:{name}" for name, passed in w_checks.items() if not passed)
    predicted[w] = prime_character(-A, w)

    # At the distinct wild prime q1, x has even valuation -4.
    expected_wild_unit = (
        (A % q1) * pow(unit_residue(F(q1) ** 2 * c, q1), 2, q1)
    ) % q1
    q1_checks = {
        "v_q1(c)=-2": valuation(c, q1) == -2,
        "v_q1(M)=-4": valuation(M, q1) == -4,
        "v_q1(x)=-4": valuation(x, q1) == -4,
        "v_q1(d)=1": valuation(d, q1) == 1,
        "u_x=A*(q1^2*c)^2": unit_residue(x, q1) == expected_wild_unit,
    }
    failures.extend(
        f"wild-place:{name}" for name, passed in q1_checks.items() if not passed
    )
    predicted[q1] = prime_character(A, q1)

    # For odd a,z,b, c is in 16 Z_2, v_2(M)=4, and the exact 2-adic
    # formula reduces to (x,d)_2=+1 iff b=a (mod 4).
    two_checks = {
        "v2(c)>=4": valuation(c, 2) >= 4,
        "v2(M)=4": valuation(M, 2) == 4,
        "v2(x)=4": valuation(x, 2) == 4,
        "v2(d)=1": valuation(d, 2) == 1,
    }
    failures.extend(f"2-place:{name}" for name, passed in two_checks.items() if not passed)
    predicted[2] = 1 if b.numerator % 4 == a % 4 else -1

    # These fixed places are not part of the new reciprocity identity.  The
    # direct local formula predicts them independently; the existence proof
    # below instead makes d a square at each one not dividing A.
    for prime in SMALL_CONTROLLED:
        if prime not in factors:
            predicted[prime] = direct_hilbert(x, d, prime)

    predicted["oo"] = direct_hilbert(x, d, "oo")
    return predicted, failures


def empty_character_branch() -> dict:
    return {
        "eligible": 0,
        "anti_correlation_breaks": 0,
        "key_predicted_aligned": 0,
        "key_actual_aligned": 0,
        "full_predicted_aligned": 0,
        "full_actual_aligned": 0,
        "mismatches": 0,
    }


def audit_branch(a: int, w: int, eps: int, q1_values: list[int], pacer: Pacer) -> tuple[dict, list[dict]]:
    A = 1 + 4 * a * a
    factors = factor_dict(A)
    z = F(-w)
    Z = z**3
    D_fraction = 1 - Z - a * a * Z * Z
    assert D_fraction.denominator == 1
    D = int(D_fraction)
    fixed_failure = base_failure(A, w, D)
    exclusions: Counter[str] = Counter()
    branches = {str(character): empty_character_branch() for character in (-1, 1)}
    mismatch_rows: list[dict] = []
    eligible = 0
    actual_computed = 0
    actual_none = 0
    anti_break_count = 0
    key_predicted_count = 0
    key_actual_count = 0
    full_predicted_count = 0
    full_actual_count = 0
    first_key = None
    first_full = None

    for q1 in q1_values:
        pacer.tick()
        assert h10q._is_prime(q1)
        b = F(eps * w * q1)
        actual = l13_filter.aligned(F(a), z, F(0), b)
        if actual is None:
            actual_none += 1
            exclusions["aligned-returned-none"] += 1
            continue
        actual_computed += 1

        failure = fixed_failure or q1_failure(a, A, w, D, q1)
        if failure is not None:
            exclusions[failure] += 1
            continue

        eligible += 1
        data = arithmetic_data(a, w, eps, q1)
        predicted, derivation_failures = specialized_prediction(data, eps, w, q1)
        q1_character = predicted[q1]
        qbranch = branches[str(q1_character)]
        qbranch["eligible"] += 1

        A_product_actual = math.prod(actual[prime] for prime in factors)
        A_product_predicted = math.prod(predicted[prime] for prime in factors)
        expected_pair_product = factored_jacobi(-2 * eps * w, factors)
        actual_pair_product = A_product_actual * actual[q1]
        predicted_pair_product = A_product_predicted * predicted[q1]
        anti_break_predicted = expected_pair_product == 1
        anti_break_actual = actual_pair_product == 1
        if anti_break_predicted:
            anti_break_count += 1
            qbranch["anti_correlation_breaks"] += 1

        key_places = set(factors) | {w, q1}
        key_predicted = all(predicted[place] == 1 for place in key_places)
        key_actual = all(actual[place] == 1 for place in key_places)
        full_predicted = all(value == 1 for value in predicted.values())
        full_actual = all(value == 1 for value in actual.values())
        key_predicted_count += int(key_predicted)
        key_actual_count += int(key_actual)
        full_predicted_count += int(full_predicted)
        full_actual_count += int(full_actual)
        qbranch["key_predicted_aligned"] += int(key_predicted)
        qbranch["key_actual_aligned"] += int(key_actual)
        qbranch["full_predicted_aligned"] += int(full_predicted)
        qbranch["full_actual_aligned"] += int(full_actual)

        if key_actual and first_key is None:
            first_key = {"q1": q1, "symbols": symbol_map(actual)}
        if full_actual and first_full is None:
            first_full = {"q1": q1, "symbols": symbol_map(actual)}

        symbol_differences = {
            str(place): {"predicted": predicted.get(place), "actual": actual.get(place)}
            for place in set(predicted) | set(actual)
            if predicted.get(place) != actual.get(place)
        }
        failed_checks = list(derivation_failures)
        if predicted_pair_product != expected_pair_product:
            failed_checks.append("predicted-pair-product-identity")
        if actual_pair_product != expected_pair_product:
            failed_checks.append("actual-pair-product-identity")
        if anti_break_predicted != anti_break_actual:
            failed_checks.append("anti-correlation-status")
        if key_predicted != key_actual:
            failed_checks.append("key-alignment-status")
        if full_predicted != full_actual:
            failed_checks.append("full-alignment-status")
        if symbol_differences or failed_checks:
            qbranch["mismatches"] += 1
            mismatch_rows.append(
                {
                    "type": "audit-mismatch",
                    "label": "OPEN",
                    "hypothesis_clause": HYPOTHESIS,
                    "a": a,
                    "A": A,
                    "w": w,
                    "eps": eps,
                    "q1": q1,
                    "failed_checks": failed_checks,
                    "symbol_differences": symbol_differences,
                    "expected_pair_product": expected_pair_product,
                    "predicted_pair_product": predicted_pair_product,
                    "actual_pair_product": actual_pair_product,
                }
            )

    minus_two_character = factored_jacobi(-2 * eps, factors)
    f_character = factored_jacobi(w, factors)
    character_A_mod_w = prime_character(A, w)
    criterion_breaks = (
        fixed_failure is None
        and minus_two_character != 0
        and f_character != 0
        and minus_two_character * f_character == 1
    )
    assert A % 8 == 5 and minus_two_character == -1
    if fixed_failure is None:
        assert f_character == character_A_mod_w
    assert anti_break_count == (eligible if criterion_breaks else 0)

    record = {
        "type": "audit-branch",
        "label": "PROVED",
        "bounded_scan_label": "EVIDENCE",
        "hypothesis_clause": HYPOTHESIS,
        "cell": [w, [-1, 1]],
        "a": a,
        "A": A,
        "A_mod_8": A % 8,
        "A_factorization": {str(prime): exponent for prime, exponent in factors.items()},
        "A_prime": len(factors) == 1 and next(iter(factors.values())) == 1,
        "eps": eps,
        "f": w,
        "tau": "0",
        "base_hypothesis_met": fixed_failure is None,
        "base_failure": fixed_failure,
        "characters": {
            "(f|A)": f_character,
            "(A|w)": character_A_mod_w,
            "(-2eps|A)": minus_two_character,
            "pair_product": minus_two_character * f_character,
        },
        "criterion_breaks_anti_correlation": criterion_breaks,
        "q1_total": len(q1_values),
        "actual_aligned_calls_completed": actual_computed,
        "actual_none": actual_none,
        "eligible": eligible,
        "exclusion_counts": dict(sorted(exclusions.items())),
        "q1_character_branches_(A|q1)": branches,
        "anti_correlation_break_count": anti_break_count,
        "key_predicted_aligned_count": key_predicted_count,
        "key_actual_aligned_count": key_actual_count,
        "full_predicted_aligned_count": full_predicted_count,
        "full_actual_aligned_count": full_actual_count,
        "first_key_actual_aligned": first_key,
        "first_full_actual_aligned": first_full,
        "mismatch_count": len(mismatch_rows),
        "refusals_used_as_evidence": False,
    }
    return record, mismatch_rows


def character_sum_witness(w: int) -> dict:
    """Construct an odd a satisfying both (A|w)=-1 and gcd(A,D)=1."""
    residues = [
        residue
        for residue in range(w)
        if prime_character(1 + 4 * residue * residue, w) == -1
    ]
    assert len(residues) == (w + 1) // 2
    chosen_residue = residues[0]

    # If p divides both A and D, then 4D=(Z-2)^2 mod p for Z=-w^3,
    # hence p divides w^3+2.  Force a=0 at every such possible p, so A=1.
    avoidance_factors = factor_dict(w**3 + 2)
    avoidance_modulus = math.prod(avoidance_factors)
    assert math.gcd(avoidance_modulus, w) == 1
    multiplier = (chosen_residue * pow(avoidance_modulus, -1, w)) % w
    a = avoidance_modulus * multiplier
    if a % 2 == 0:
        a += avoidance_modulus * w
    assert a % 2 == 1 and a % w == chosen_residue and a % avoidance_modulus == 0

    A = 1 + 4 * a * a
    z = -w
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    assert prime_character(A, w) == -1
    assert math.gcd(A, 2 * z * D) == 1
    return {
        "type": "uniform-witness",
        "label": "PROVED",
        "hypothesis_clause": (
            "w is an odd proven prime with w=3 (mod 4); choose r with "
            "(1+4r^2|w)=-1, then impose a=r (mod w), a=0 modulo "
            "rad(w^3+2), and a odd."
        ),
        "w": w,
        "w_mod_20": w % 20,
        "nonsquare_residue_count": len(residues),
        "expected_nonsquare_residue_count": (w + 1) // 2,
        "chosen_residue_mod_w": chosen_residue,
        "avoidance_rad_w3_plus_2": avoidance_modulus,
        "avoidance_factorization": {
            str(prime): exponent for prime, exponent in avoidance_factors.items()
        },
        "constructed_a": a,
        "constructed_A": A,
        "A_mod_8": A % 8,
        "(A|w)": prime_character(A, w),
        "gcd(A,2*z*D)": math.gcd(A, 2 * z * D),
        "criterion_breaks": True,
    }


def theorem_records() -> list[dict]:
    return [
        {
            "type": "criterion-theorem",
            "label": "PROVED",
            "hypothesis_clause": HYPOTHESIS,
            "scope": (
                "Factor-by-factor algebraic extension of L10b to integral composite A; "
                "no h10q or l13_filter function was changed or generalized."
            ),
            "at_p_power_exactly_dividing_A": {
                "notation": "p^e || A, A'=A/p^e",
                "valuations": {"v_p(c)": "-e", "v_p(M)": "-2e", "v_p(x)": "-e", "v_p(d)": "e"},
                "unit": "u_x = A'*(p^e*c)^2 (mod p)",
                "symbol": "(x,d)_p = 1 if e is even; (-2*eps*f*q1|p) if e is odd",
                "prime_A_specialization": "v_A(x)=-1, v_A(d)=1, u_x=(A*c)^2 mod A",
            },
            "wild_symbol": "(x,d)_q1=(A|q1)",
            "f_symbol": "for f=w prime and z=-w, (x,d)_w=(-A|w)",
            "aggregate_identity": (
                "product_{p|A}(x,d)_p * (A|q1) = "
                "(-2*eps|A)*(f|A)"
            ),
            "exact_break_criterion": "(-2*eps|A)*(f|A)=+1",
            "A_mod_8_cases": {
                "1": "(-2eps|A)=+1, so break iff (f|A)=+1",
                "5": "(-2eps|A)=-1, so break iff (f|A)=-1",
            },
            "admissible_odd_a_specialization": (
                "A=5 (mod 8), hence the L10b anti-correlation breaks iff (f|A)=-1. "
                "For prime A, the A-place and wild place are both +1 exactly when "
                "(f|A)=-1 and (A|q1)=+1.  For composite A, also require "
                "(-2*eps*f*q1|p)=+1 at every p of odd exponent in A."
            ),
        },
        {
            "type": "always-breakable-theorem",
            "label": "PROVED",
            "hypothesis_clause": "w is a proven prime with w=11 or 19 (mod 20), z=-w, f=w.",
            "verdict": "PROVED-always-breakable",
            "scope": (
                "The frozen/class-side canonical 5-wall is always broken by some tau=0 "
                "fixed-factor family with variable odd a and infinitely many prime q1. "
                "This does not prove an emergent-free member or an unconditional closure."
            ),
            "proof_steps": [
                "For w=3 mod 4, sum_{r mod w}(1+4r^2|w)=-1 and there are no zeros; therefore exactly (w+1)/2 residues have character -1.",
                "Choose such r.  CRT with a=0 modulo rad(w^3+2), followed by a parity lift, gives odd a with (A|w)=-1 and gcd(A,2*z*D)=1.",
                "Every p|A is 1 mod 4.  The criterion gives (w|A)=(A|w)=-1, aligns the f=w place, and removes the A/wild anti-correlation.",
                "Choose q1 locally so (-2*eps*w*q1|p)=+1 for every odd-exponent p|A, b=a mod 4, and d=-2*A*b is a square at each of 3,5,7 not dividing A.",
                "CRT gives a reduced residue class for q1.  Dirichlet gives infinitely many primes in it; discard finitely many data divisors and take q1 large enough that x>0 at infinity.",
                "Then every frozen place 2,3,5,7,p|A,w,q1,infinity has symbol +1; the existing L10 Taylor-exponent lemma, with every prime of fixed f frozen, produces an aligned fixed-factor progression.",
            ],
        },
    ]


def build_wall_rows(ws: list[int], audit_rows: list[dict], witnesses: list[dict]) -> list[dict]:
    witness_by_w = {row["w"]: row for row in witnesses}
    rows = []
    for w in ws:
        local = [row for row in audit_rows if row["cell"][0] == w]
        break_as = sorted(
            {
                row["a"]
                for row in local
                if row["base_hypothesis_met"] and row["criterion_breaks_anti_correlation"]
            }
        )
        full_candidates = [
            row for row in local if row["first_full_actual_aligned"] is not None
        ]
        first_full = None
        if full_candidates:
            selected = min(
                full_candidates,
                key=lambda row: (
                    row["first_full_actual_aligned"]["q1"],
                    row["a"],
                    0 if row["eps"] == 1 else 1,
                ),
            )
            first_full = {
                "a": selected["a"],
                "eps": selected["eps"],
                **selected["first_full_actual_aligned"],
            }
        a7 = next(row for row in local if row["a"] == 7 and row["eps"] == 1)
        rows.append(
            {
                "type": "wall-summary",
                "label": "PROVED",
                "bounded_scan_label": "EVIDENCE",
                "hypothesis_clause": HYPOTHESIS,
                "cell": [w, [-1, 1]],
                "finite_a_pool": list(A_VALUES),
                "finite_pool_break_a": break_as,
                "finite_pool_has_break": bool(break_as),
                "finite_pool_first_full_q1_le_200": first_full,
                "a7_fixed_family_verdict": (
                    "breakable" if a7["criterion_breaks_anti_correlation"] else "obstructed"
                ),
                "a7_(w|197)": a7["characters"]["(f|A)"],
                "uniform_constructed_a": witness_by_w[w]["constructed_a"],
                "uniform_theorem_applies": witness_by_w[w]["criterion_breaks"],
                "per_branch_counts": {
                    "audit_rows": len(local),
                    "eligible": sum(row["eligible"] for row in local),
                    "key_actual_aligned": sum(row["key_actual_aligned_count"] for row in local),
                    "full_actual_aligned": sum(row["full_actual_aligned_count"] for row in local),
                    "mismatches": sum(row["mismatch_count"] for row in local),
                },
            }
        )
    return rows


def build_report(
    audit_rows: list[dict],
    wall_rows: list[dict],
    summary: dict,
    verdict: dict,
) -> str:
    lines = [
        "# L19 tau=0 fixed-factor wall criterion",
        "",
        "## Verdict",
        "",
        "**PROVED-always-breakable (frozen/class-side scope).**  For every prime "
        "$w\\equiv11,19\\pmod{20}$, some variable odd $a$ and infinitely many prime "
        "$q_1$ give an aligned $\\tau=0$, $f=w$ fixed-factor class.  This does "
        "not prove that any member is emergent-free, so unconditional member closure remains OPEN.",
        "",
        "The fixed $a=7$ subfamily is not universal: it is **PROVED-breakable-iff** "
        "$(w\\mid197)=-1$.  The three closures at 131, 139, and 151 all lie on that "
        "character branch; they do not prove a universal $a=7$ closure theorem.",
        "",
        f"Exact tested hypothesis clause: `{HYPOTHESIS}`",
        "",
        "## Algebraic criterion (PROVED)",
        "",
        "Put $A=1+4a^2$, $b=\\varepsilon f q_1$, $\\delta=1$, "
        "$\\alpha=-A$, $x=-AM$, and $d=-2Ab$.  Let $p^e\\Vert A$ and "
        "$A'=A/p^e$.  Every $p\\mid A$ satisfies $p\\equiv1\\pmod4$.  Under "
        "the displayed hypothesis, the exact valuations and unit are",
        "",
        "$$v_p(c)=-e,\\quad v_p(M)=-2e,\\quad v_p(x)=-e,\\quad "
        "v_p(d)=e,\\qquad u_x\\equiv A'(p^ec)^2\\pmod p.$$",
        "",
        "Consequently",
        "",
        "$$ (x,d)_p=\\begin{cases}1,&e\\text{ even},\\\\"
        "(-2\\varepsilon f q_1\\mid p),&e\\text{ odd}.\\end{cases} $$",
        "",
        "For prime $A$ this specializes exactly to $v_A(x)=-1$, $v_A(d)=1$, "
        "$u_x=-A^2M\\equiv(Ac)^2\\pmod A$, and",
        "",
        "$$ (x,d)_A=(-2\\varepsilon\\mid A)(f\\mid A)(q_1\\mid A). $$",
        "",
        "At the wild prime, $(x,d)_{q_1}=(A\\mid q_1)$.  Reciprocity (all "
        "$p\\mid A$ are $1$ mod 4) therefore gives the factor-by-factor identity",
        "",
        "$$\\prod_{p\\mid A}(x,d)_p\\,(A\\mid q_1)="
        "(-2\\varepsilon\\mid A)(f\\mid A).$$",
        "",
        "Thus anti-correlation breaks exactly when "
        "$(-2\\varepsilon\\mid A)(f\\mid A)=+1$.  If $A\\equiv1\\pmod8$ "
        "this means $(f\\mid A)=+1$; if $A\\equiv5\\pmod8$ it means "
        "$(f\\mid A)=-1$.  For admissible odd $a$, always $A\\equiv5\\pmod8$ "
        "and $(-2\\varepsilon\\mid A)=-1$ for both signs.",
        "",
        "For $f=w$ and $z=-w$, one has $v_w(c)=4$, $v_w(M)=v_w(x)=0$, "
        "$v_w(d)=1$, and $(x,d)_w=(-A\\mid w)$.  Since these wall primes are "
        "$3$ mod 4 and $(w\\mid A)=(A\\mid w)$, the same condition "
        "$(w\\mid A)=-1$ also makes the $w$-symbol $+1$.",
        "",
        "For composite $A$, aggregate correlation is not enough to align each "
        "$A$-place: additionally require $(-2\\varepsilon wq_1\\mid p)=+1$ "
        "for every $p$ occurring to odd exponent.  This is a PROVED "
        "factor-by-factor extension of the algebra, not a generalized kernel certificate; "
        "no kernel function was changed.",
        "",
        "## Why every canonical 5-wall is breakable (PROVED)",
        "",
        "For prime $w\\equiv3\\pmod4$, the standard quadratic character sum is "
        "$\\sum_{r\\bmod w}(1+4r^2\\mid w)=-1$.  There are no zero terms, so "
        "exactly $(w+1)/2$ residues have character $-1$.  Choose one.  To retain "
        "the L10 unit hypothesis, impose by CRT $a\\equiv0$ modulo "
        "$\\operatorname{rad}(w^3+2)$ and use an odd lift.  Indeed, for $p\\mid A$",
        "",
        "$$4D\\equiv(Z-2)^2\\pmod p,\\qquad Z=-w^3,$$",
        "",
        "so this auxiliary congruence forces $\\gcd(A,D)=1$ while preserving "
        "$(A\\mid w)=-1$.",
        "",
        "For the remaining places choose a reduced residue class for $q_1$ by CRT: "
        "make $(-2\\varepsilon wq_1\\mid p)=+1$ at each odd-exponent "
        "$p\\mid A$; impose $\\varepsilon wq_1\\equiv a\\pmod4$ (the exact "
        "2-adic criterion); and make $d=-2A\\varepsilon wq_1$ a square at each "
        "of $3,5,7$ not dividing $A$.  Dirichlet supplies infinitely many prime "
        "$q_1$ in the resulting class.  After discarding finitely many data divisors "
        "and taking $q_1$ large, the real symbol is also $+1$.  The existing L10 "
        "Taylor-exponent lemma (freezing every prime of fixed $f$) then turns this "
        "base into an aligned fixed-factor progression.",
        "",
        "## Mechanical audit",
        "",
        f"Measured wall-clock for kernel self-test plus theorem construction and audit "
        f"(before artifact serialization): **{summary['audit_wall_seconds']:.6f} s**.",
        "",
        f"The script called `l13_filter.aligned` on all {summary['candidate_count']} "
        f"$(a,w,\\varepsilon,q_1)$ tuples.  {summary['eligible_count']} met the exact "
        "hypothesis.  It independently predicted every returned symbol from the "
        "valuation/unit formulas above (and the literal local formula at unrelated "
        "small places).",
        "",
        f"- symbol-map mismatches: **{summary['symbol_map_mismatch_count']}**",
        f"- valuation/unit or pair-identity mismatches: **{summary['derivation_mismatch_count']}**",
        f"- full aligned/obstructed status mismatches: **{summary['full_status_mismatch_count']}**",
        f"- full aligned eligible bases in the bounded scan: **{summary['full_actual_aligned_count']}**",
        f"- refusals used as evidence: **no**",
        "",
        "### Per-(a,w,eps) audit table",
        "",
        "| a | w | eps | A | (w|A) | base H | eligible/46 | key +/+ | full aligned | mismatches |",
        "|---:|---:|---:|---:|---:|:---:|---:|---:|---:|---:|",
    ]
    for row in audit_rows:
        lines.append(
            "| {a} | {w} | {eps:+d} | {A} | {character:+d} | {base} | "
            "{eligible}/46 | {key} | {full} | {mismatch} |".format(
                a=row["a"],
                w=row["cell"][0],
                eps=row["eps"],
                A=row["A"],
                character=row["characters"]["(f|A)"],
                base="yes" if row["base_hypothesis_met"] else "no",
                eligible=row["eligible"],
                key=row["key_actual_aligned_count"],
                full=row["full_actual_aligned_count"],
                mismatch=row["mismatch_count"],
            )
        )
    lines.extend(
        [
            "",
            "### Per-wall result",
            "",
            "| w | finite-pool break a values | first full q1<=200 witness | a=7 | uniform theorem |",
            "|---:|:---|:---|:---:|:---:|",
        ]
    )
    for row in wall_rows:
        witness = row["finite_pool_first_full_q1_le_200"]
        witness_text = (
            "none"
            if witness is None
            else f"a={witness['a']}, eps={witness['eps']:+d}, q1={witness['q1']}"
        )
        lines.append(
            f"| {row['cell'][0]} | {','.join(map(str, row['finite_pool_break_a']))} | "
            f"{witness_text} | {row['a7_fixed_family_verdict']} | PROVED |"
        )
    lines.extend(
        [
            "",
            "## Scope ledger",
            "",
            "- **PROVED:** exact factor-by-factor criterion; fixed $a$ breakability iff "
            "the stated character condition; existence of an aligned variable-$a$ "
            "tau=0 fixed-factor class for every wall prime.",
            "- **EVIDENCE:** the bounded $a\\le25$, $q_1\\le200$ audit and the previously "
            "observed k=0 member closures.",
            "- **OPEN:** an unconditional emergent-free member (hence full closure) for every $w$.",
            "",
            f"Primary machine verdict: `{verdict['verdict']}`; member-level verdict: "
            f"`{verdict['member_level_verdict']}`.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    started = time.perf_counter()
    h10q._selftest()
    q1_values = list(h10q.primerange(2, Q1_MAX + 1))
    ws = [w for w in h10q.primerange(W_MIN, W_MAX + 1) if w % 20 in (11, 19)]
    assert all(h10q._is_prime(w) and w % 4 == 3 for w in ws)

    records = [
        {
            "type": "metadata",
            "label": "PROVED",
            "schema": "l19-tauzero-v1",
            "source": SOURCE,
            "kernel_authority": ["h10q.py", "l13_filter.py"],
            "kernel_selftest_first": True,
            "a_values": list(A_VALUES),
            "w_range_inclusive": [W_MIN, W_MAX],
            "w_clause": "proven prime and w=11 or 19 (mod 20)",
            "w_values": ws,
            "eps_values": [1, -1],
            "q1_range_inclusive": [2, Q1_MAX],
            "q1_values": q1_values,
            "pace": {"sleep_seconds": PACE_SECONDS, "every_candidates": PACE_EVERY},
            "hypothesis_clause": HYPOTHESIS,
            "refusals_are_evidence": False,
        }
    ]
    records.extend(theorem_records())

    uniform_witnesses = [character_sum_witness(w) for w in ws]
    records.extend(uniform_witnesses)

    pacer = Pacer()
    audit_rows: list[dict] = []
    mismatches: list[dict] = []
    for a in A_VALUES:
        for w in ws:
            for eps in (1, -1):
                row, local_mismatches = audit_branch(a, w, eps, q1_values, pacer)
                audit_rows.append(row)
                mismatches.extend(local_mismatches)
    records.extend(audit_rows)
    records.extend(mismatches)

    wall_rows = build_wall_rows(ws, audit_rows, uniform_witnesses)
    records.extend(wall_rows)
    audit_seconds = time.perf_counter() - started

    mismatch_case_count = len(mismatches)
    symbol_map_mismatch_count = sum(bool(row["symbol_differences"]) for row in mismatches)
    derivation_mismatch_count = sum(bool(row["failed_checks"]) for row in mismatches)
    full_status_mismatch_count = sum(
        "full-alignment-status" in row["failed_checks"] for row in mismatches
    )
    summary = {
        "type": "audit-summary",
        "label": "PROVED" if mismatch_case_count == 0 else "OPEN",
        "bounded_scan_label": "EVIDENCE",
        "hypothesis_clause": HYPOTHESIS,
        "audit_wall_seconds": audit_seconds,
        "w_count": len(ws),
        "audit_branch_count": len(audit_rows),
        "candidate_count": pacer.candidates,
        "actual_aligned_call_count": sum(
            row["actual_aligned_calls_completed"] for row in audit_rows
        ),
        "eligible_count": sum(row["eligible"] for row in audit_rows),
        "prime_A_eligible_count": sum(
            row["eligible"] for row in audit_rows if row["A_prime"]
        ),
        "composite_A_eligible_count": sum(
            row["eligible"] for row in audit_rows if not row["A_prime"]
        ),
        "anti_correlation_break_count": sum(
            row["anti_correlation_break_count"] for row in audit_rows
        ),
        "key_actual_aligned_count": sum(
            row["key_actual_aligned_count"] for row in audit_rows
        ),
        "full_actual_aligned_count": sum(
            row["full_actual_aligned_count"] for row in audit_rows
        ),
        "mismatch_case_count": mismatch_case_count,
        "symbol_map_mismatch_count": symbol_map_mismatch_count,
        "derivation_mismatch_count": derivation_mismatch_count,
        "full_status_mismatch_count": full_status_mismatch_count,
        "mismatches": mismatches,
        "finite_a_pool_breaks_every_bounded_wall": all(
            row["finite_pool_has_break"] for row in wall_rows
        ),
        "finite_a_pool_has_full_q1_le_200_for_every_bounded_wall": all(
            row["finite_pool_first_full_q1_le_200"] is not None for row in wall_rows
        ),
        "uniform_character_sum_witness_for_every_wall": all(
            row["uniform_theorem_applies"] for row in wall_rows
        ),
        "refusals_used_as_evidence": False,
    }
    records.append(summary)

    a7_break_w = [
        row["cell"][0]
        for row in wall_rows
        if row["a7_fixed_family_verdict"] == "breakable"
    ]
    a7_obstructed_w = [
        row["cell"][0]
        for row in wall_rows
        if row["a7_fixed_family_verdict"] == "obstructed"
    ]
    verdict = {
        "type": "verdict",
        "label": "PROVED",
        "hypothesis_clause": "w is prime and w=11 or 19 (mod 20), z=-w, f=w.",
        "verdict": "PROVED-always-breakable",
        "scope": "variable-a tau=0 frozen/class-side fixed-factor family",
        "fixed_a7_verdict": "PROVED-breakable-iff-(w|197)=-1",
        "fixed_a7_bounded_break_w": a7_break_w,
        "fixed_a7_bounded_obstructed_w": a7_obstructed_w,
        "member_level_verdict": "OPEN-unconditional-emergent-free-member",
        "three_k0_closures_label": "EVIDENCE-not-a-universal-member-theorem",
        "refusals_used_as_evidence": False,
    }
    records.append(verdict)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
    REPORT.write_text(build_report(audit_rows, wall_rows, summary, verdict), encoding="utf-8")

    # One-pass replay verification of both requested artifacts.
    loaded = [json.loads(line) for line in OUT.read_text(encoding="utf-8").splitlines()]
    loaded_summary = next(row for row in loaded if row.get("type") == "audit-summary")
    loaded_verdict = next(row for row in loaded if row.get("type") == "verdict")
    assert loaded_summary["candidate_count"] == len(A_VALUES) * len(ws) * 2 * len(q1_values)
    assert loaded_summary["actual_aligned_call_count"] == loaded_summary["candidate_count"]
    assert loaded_summary["mismatch_case_count"] == 0
    assert loaded_summary["symbol_map_mismatch_count"] == 0
    assert loaded_summary["derivation_mismatch_count"] == 0
    assert loaded_summary["full_status_mismatch_count"] == 0
    assert loaded_summary["uniform_character_sum_witness_for_every_wall"]
    assert loaded_verdict["verdict"] == "PROVED-always-breakable"
    report_text = REPORT.read_text(encoding="utf-8")
    assert "PROVED-always-breakable" in report_text
    assert "Measured wall-clock" in report_text
    assert "unconditional member closure remains OPEN" in report_text

    print(f"WROTE {OUT}")
    print(f"WROTE {REPORT}")
    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
    print(json.dumps(verdict, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
