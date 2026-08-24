#!/usr/bin/env python3
"""Audit Schinzel admissibility for the uniform aligned-class construction.

Replay from the workspace repository root:

    nice -n 19 python3 math/h10q/l20_admissible.py

The arithmetic verdicts come from h10q.py.  Existing JSON artifacts are used
only as inputs to independently replay their parameters/certificates.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h10q  # noqa: E402  (the stipulated exact-arithmetic authority)


OUT = HERE / "data" / "l20_admissible.jsonl"
REPORT = Path("/tmp/l20_admissible.md")
L19 = HERE / "data" / "l19_classexist.jsonl"
L13_FACTORS = HERE / "data" / "l12_factorP.json"
A_POOL = (1, 3, 5, 7)
IRRED_PRIME_LIMIT = 1000
FROZEN_SAMPLE_T = (0, 1, 2, 5)
PACE_EVERY = 200
PACE_SECONDS = 0.005


class Pacer:
    """Keep the requested low-duty-cycle pause in candidate/certificate loops."""

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


def gcd_many(values: list[int]) -> int:
    answer = 0
    for value in values:
        answer = math.gcd(answer, abs(value))
    return answer


def evaluate(coefficients: list[int] | list[F], value: int | F) -> int | F:
    answer: int | F = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def compose(P: list[F], b0: F, step: F) -> list[F]:
    """Ascending coefficients of P(b0 + step*t), exactly."""
    out = [F(0)] * len(P)
    for degree, coefficient in enumerate(P):
        if not coefficient:
            continue
        for t_degree in range(degree + 1):
            out[t_degree] += (
                coefficient
                * math.comb(degree, t_degree)
                * b0 ** (degree - t_degree)
                * step**t_degree
            )
    return out


def primitive_integer_poly(coefficients: list[F]) -> tuple[list[int], F]:
    """Return primitive G and positive c such that polynomial = c*G."""
    denominator = 1
    for coefficient in coefficients:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in coefficients]
    content = gcd_many(integers)
    assert content > 0
    return [value // content for value in integers], F(content, denominator)


def quadratic_character(value: int, prime: int) -> int:
    """Ordinary Legendre character, including its zero value."""
    residue = value % prime
    if residue == 0:
        return 0
    return 1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1


def remove_support(value: F, support: list[int] | set[int]) -> F:
    answer = F(value)
    for prime in support:
        answer /= F(prime) ** h10q.vp(answer, prime)
    return answer


def supported_integer(value: int, support: list[int] | set[int]) -> bool:
    remainder = abs(value)
    for prime in support:
        while remainder % prime == 0:
            remainder //= prime
    return remainder == 1


def choose_a(prime: int) -> int:
    for a in A_POOL:
        if quadratic_character(1 + 4 * a * a, prime) == -1:
            return a
    raise AssertionError((prime, "the frozen grid's a-pool is empty"))


def trace_polynomial(P: list[F]) -> list[F]:
    """Q with P(b)=b^4*Q(b+b^-1), for a palindromic octic P."""
    assert len(P) == 9 and all(P[index] == P[8 - index] for index in range(9))
    Q = [
        P[4] + 2 * P[0] - 2 * P[2],
        P[3] - 3 * P[1],
        P[2] - 4 * P[0],
        P[1],
        P[0],
    ]
    lifted = [
        Q[4],
        Q[3],
        Q[2] + 4 * Q[4],
        Q[1] + 3 * Q[3],
        Q[0] + 2 * Q[2] + 6 * Q[4],
        Q[1] + 3 * Q[3],
        Q[2] + 4 * Q[4],
        Q[3],
        Q[4],
    ]
    assert lifted == P
    return Q


def class_parameters(
    w: int,
    unit: tuple[int, int],
    a_value: int,
    f_value: int | None = None,
) -> dict[str, Any]:
    a = F(a_value)
    f = w if f_value is None else f_value
    z = F(w) * F(*unit)
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    assert a.denominator == A.denominator == s.denominator == 1
    assert delta == -4 * a**4 / A < 0
    assert alpha == (2 * a * a) ** 2
    assert z and D
    assert h10q.vp(z, w) >= 1 and h10q.vp(z, f) > 0
    assert quadratic_character(A.numerator, f) == -1
    S = sorted(
        {2, 3, 5, 7, f}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
    )
    P = h10q._l10_P(a, Z, D, A, delta, s)
    symbolic_lead = 4 * a**8 * A * Z**4
    assert len(P) == 9 and P[0] == P[-1] == symbolic_lead > 0
    perturbation = 32 * A**3 * s**2 * D**2
    palindromic_core = list(P)
    palindromic_core[5] += perturbation
    assert all(
        palindromic_core[index] == palindromic_core[8 - index]
        for index in range(9)
    )
    assert all(
        P[index] - P[8 - index]
        == (
            perturbation
            if index == 3
            else -perturbation
            if index == 5
            else 0
        )
        for index in range(9)
    )
    trace_Q = trace_polynomial(P) if s == 0 else None
    assert (trace_Q is not None) == (a == 1)
    return {
        "w": w,
        "unit": unit,
        "a": a,
        "z": z,
        "A": A,
        "tau": tau,
        "delta": delta,
        "alpha": alpha,
        "Z": Z,
        "D": D,
        "s": s,
        "S": S,
        "P": P,
        "f": f,
        "palindromic_core": True,
        "palindromic_perturbation": perturbation,
        "trace_Q": trace_Q,
    }


def construct_class(
    w: int,
    unit: tuple[int, int],
    a_value: int,
    f_value: int | None = None,
) -> dict[str, Any]:
    """Rebuild the first clean class in L19's explicit residue subsystem."""
    data = class_parameters(w, unit, a_value, f_value)
    f = data["f"]
    a, z = data["a"], data["z"]
    A, delta, alpha = data["A"], data["delta"], data["alpha"]
    Z, D, s, S, P = (
        data["Z"],
        data["D"],
        data["s"],
        data["S"],
        data["P"],
    )
    attempts = 0
    for q1 in range(3, 10000, 2):
        PACER.tick()
        attempts += 1
        if q1 in S or math.gcd(q1, math.prod(S)) != 1:
            continue
        if any(
            quadratic_character(2 * f * q1, prime) != 1
            for prime in S
            if prime not in (2, f)
        ):
            continue
        if not h10q._is_prime(q1):
            continue
        if h10q.vp(z, q1) != 0 or h10q.vp(D, q1) != 0:
            continue

        b0 = F(f * q1)
        Pb0 = F(evaluate(P, b0))
        if Pb0 == 0:
            continue
        ks = {prime: h10q._l10_exponent(P, b0, prime) for prime in S}
        N = 8
        for prime, exponent in ks.items():
            N = math.lcm(N, prime**exponent)
        N = math.lcm(N, 4 * A.numerator * A.denominator)
        if math.gcd(q1, N) != 1:
            continue

        c0 = h10q._sun_h(a, b0, Z)
        if c0 is None:
            continue
        M0 = 16 - delta * c0 * c0 - 32 * A * b0 * s * s
        assert M0 == Pb0 / (b0**4 * D**2 * A**2)
        x0, d0 = alpha * M0, 2 * alpha * b0
        symbols = {prime: h10q.hilbert(x0, d0, prime) for prime in S}
        if not all(value == 1 for value in symbols.values()):
            continue
        if h10q.legendre(A, q1) != 1:
            continue

        assert math.gcd(q1, 4 * A.numerator * math.prod(S)) == 1
        assert supported_integer(N, S)
        assert all(N % prime == 0 for prime in S)
        data.update(
            {
                "f": f,
                "eps": 1,
                "q1": q1,
                "b0": b0,
                "Pb0": Pb0,
                "ks": ks,
                "N": N,
                "symbols": symbols,
                "x0": x0,
                "d0": d0,
                "attempts": attempts,
            }
        )
        return data
    raise AssertionError((w, unit, a_value, f, "no clean q1 below 10000"))


def expected_content(data: dict[str, Any]) -> F:
    """The exact content formula, including denominator primes outside S."""
    answer = F(1)
    for prime in data["S"]:
        answer *= F(prime) ** h10q.vp(data["Pb0"], prime)
    outside_denominator = h10q._l10_supp(data["z"].denominator) - set(data["S"])
    for prime in outside_denominator:
        assert h10q.vp(data["z"], prime) < 0
        answer *= F(prime) ** (12 * h10q.vp(data["z"], prime))
    return answer


def outside_support_is_square(value: F, support: list[int]) -> bool:
    remainder = remove_support(value, support)
    numerator_root = math.isqrt(abs(remainder.numerator))
    denominator_root = math.isqrt(remainder.denominator)
    return (
        remainder > 0
        and numerator_root * numerator_root == remainder.numerator
        and denominator_root * denominator_root == remainder.denominator
    )


def irred_certificate(G: list[int], support: list[int]) -> int:
    for prime in h10q.primerange(2, IRRED_PRIME_LIMIT + 1):
        PACER.tick()
        if prime in support:
            continue
        if h10q._l13_irred8(G, prime):
            return prime
    raise AssertionError("no degree-8 Frobenius certificate <= 1000")


def audit_class(data: dict[str, Any], row_type: str) -> tuple[dict[str, Any], dict[str, int]]:
    P, S = data["P"], data["S"]
    f, q1, N = data["f"], data["q1"], data["N"]
    b0, step = data["b0"], F(f * N)
    F_coefficients = compose(P, b0, step)
    G, content = primitive_integer_poly(F_coefficients)
    content_formula = expected_content(data)
    assert content == content_formula
    assert len(G) == 9 and gcd_many(G) == 1 and G[-1] > 0
    assert F_coefficients[-1] == P[-1] * step**8 == content * G[-1]

    # The Taylor exponent makes the constant coefficient uniquely minimal at
    # every p in S: nonconstant coefficients gain p (and gain 8 at p=2).
    unit_residues: dict[str, dict[str, int]] = {}
    for prime in S:
        modulus = 8 if prime == 2 else prime
        assert G[0] % prime != 0
        assert all(coefficient % modulus == 0 for coefficient in G[1:])
        assert q1 % prime != 0 and N % modulus == 0
        for t in range(prime):
            assert (q1 + N * t) % prime != 0
            assert int(evaluate(G, t)) % prime != 0
        unit_residues[str(prime)] = {
            "q1_mod_p": q1 % prime,
            "G0_mod_p": G[0] % prime,
            "nonconstant_modulus": modulus,
        }

    # Degree nine: the gcd of the first ten product values is the exact fixed
    # divisor.  This mechanically backs the uniform root-counting proof.
    pair_product_gcd = 0
    for t in range(10):
        pair_product_gcd = math.gcd(
            pair_product_gcd,
            abs((q1 + N * t) * int(evaluate(G, t))),
        )
    assert pair_product_gcd == 1

    frozen_checks = 0
    for t in FROZEN_SAMPLE_T:
        Q = q1 + N * t
        b = F(f * Q)
        Pb = F(evaluate(P, b))
        assert Pb == content * int(evaluate(G, t))
        x = data["alpha"] * Pb / (b**4 * data["D"] ** 2 * data["A"] ** 2)
        d = 2 * data["alpha"] * b
        for prime in S:
            assert h10q.vp(Pb, prime) == h10q.vp(data["Pb0"], prime)
            assert h10q.hilbert(x, d, prime) == data["symbols"][prime] == 1
            frozen_checks += 1

    certificate_prime = irred_certificate(G, S)
    primitive_P, _P_scale = primitive_integer_poly(P)
    # The same prime certifies P: outside S the affine step and content are
    # units, so reduction commutes with the invertible affine substitution.
    assert h10q._l13_irred8(primitive_P, certificate_prime)

    outside_part = remove_support(content, S)
    row = {
        "type": row_type,
        "label": "PROVED",
        "scope": "exact per-row cross-check; not a uniform irreducibility theorem",
        "cell": [data["w"], list(data["unit"])],
        "z": frac_text(data["z"]),
        "a": frac_text(data["a"]),
        "A": frac_text(data["A"]),
        "tau": frac_text(data["tau"]),
        "delta": frac_text(data["delta"]),
        "eps": data["eps"],
        "f": data["f"],
        "q1": data["q1"],
        "N": str(data["N"]),
        "S": S,
        "ks": {str(prime): exponent for prime, exponent in data["ks"].items()},
        "q1_search_candidates": data["attempts"],
        "F_definition": "h10q._l10_P(a,Z,D,A,delta,s) evaluated at f*(q1+N*t)",
        "P_palindromic_core": data["palindromic_core"],
        "P_palindromic": data["trace_Q"] is not None,
        "P_reciprocal_perturbation": frac_text(data["palindromic_perturbation"]),
        "P_trace_Q_coefficients_ascending": (
            [frac_text(coefficient) for coefficient in data["trace_Q"]]
            if data["trace_Q"] is not None
            else None
        ),
        "F_scale_c": frac_text(content),
        "F_scale_formula": frac_text(content_formula),
        "F_scale_outside_S": frac_text(outside_part),
        "F_scale_S_supported": outside_part == 1,
        "F_scale_square_outside_S": outside_support_is_square(content, S),
        "G_coefficients_ascending": [str(coefficient) for coefficient in G],
        "G_degree": 8,
        "G_primitive": True,
        "G_positive_leading_coefficient": True,
        "G_irreducibility_certificate_prime": certificate_prime,
        "G_irreducible_over_Q": True,
        "gcd_q1_N": math.gcd(q1, N),
        "support_N_equals_S": supported_integer(N, S) and all(N % p == 0 for p in S),
        "pair_product_fixed_divisor": pair_product_gcd,
        "pair_product_gcd_values_t": [0, 9],
        "S_unit_congruences": unit_residues,
        "frozen_symbol_sample_t": list(FROZEN_SAMPLE_T),
        "frozen_symbol_checks": frozen_checks,
        "frozen_symbol_mismatches": 0,
    }
    counts = {
        "content_mismatches": 0,
        "fixed_divisor_mismatches": 0,
        "symbol_mismatches": 0,
        "irreducibility_mismatches": 0,
    }
    return row, counts


def l19_materialized() -> tuple[dict[tuple[int, tuple[int, int]], dict[str, Any]], dict[str, Any]]:
    rows = [json.loads(line) for line in L19.read_text(encoding="utf-8").splitlines()]
    assert len(rows) == 165 and rows[0]["type"] == "meta" and rows[-1]["type"] == "summary"
    materialized: dict[tuple[int, tuple[int, int]], dict[str, Any]] = {}
    for row in rows[1:-1]:
        key = (int(row["cell"][0]), tuple(row["cell"][1]))
        if row["type"] == "canonical-escape-cell":
            cert = row["constructed_clean_class"]
        elif row["type"] == "canonical-collision-cell":
            cert = row["constructed_class"]
        else:
            raise AssertionError(row["type"])
        materialized[key] = cert
    assert len(materialized) == 163
    return materialized, rows[-1]


def compare_l19(data: dict[str, Any], prior: dict[str, Any]) -> None:
    assert prior["a"] == frac_text(data["a"])
    assert prior["A"] == frac_text(data["A"])
    assert prior["tau"] == frac_text(data["tau"])
    assert int(prior["eps"]) == data["eps"]
    assert int(prior["f"]) == data["f"]
    assert int(prior["q1"]) == data["q1"]
    assert int(prior["N"]) == data["N"]
    assert list(prior["S"]) == data["S"]
    assert {int(p): int(k) for p, k in prior["ks"].items()} == data["ks"]
    assert prior["ok"] is True


def replay_l13c() -> dict[str, Any]:
    """Replay all 706 old grid irreducibility claims through h10q._l13_irred8."""
    payload = json.loads(L13_FACTORS.read_text(encoding="utf-8"))
    records = payload["records"]
    assert payload["meta"]["n_records"] == len(records) == 706
    branch_counts = Counter()
    certificate_counts = Counter()
    for record in records:
        PACER.tick()
        assert record["degP"] == 8 and record["factor_degrees"] == [[8, 1]]
        a = F(1)
        A = F(5)
        z = F(record["w"]) * F(*record["u"])
        Z = z**3
        D = 1 - Z - Z**2
        s = F(0)
        branch = record["branch"]
        if branch == "tau=(1+2a^2)/A":
            delta = F(-4, 5)
        elif branch == "tau=2a/A":
            delta = F(1, 5)
        else:
            raise AssertionError(branch)
        P = h10q._l10_P(a, Z, D, A, delta, s)
        trace_polynomial(P)
        primitive_P, _scale = primitive_integer_poly(P)
        certificate_primes = sorted(
            int(prime)
            for prime, cycle_type in record["frobenius"].items()
            if cycle_type == [8]
        )
        assert certificate_primes
        certificate_prime = certificate_primes[0]
        assert h10q._is_prime(certificate_prime)
        assert h10q._l13_irred8(primitive_P, certificate_prime)
        branch_counts[branch] += 1
        certificate_counts[certificate_prime] += 1
    return {
        "label": "PROVED",
        "rows": len(records),
        "replayed_irreducible": len(records),
        "mismatches": 0,
        "original_engine": "Sympy 1.14 factor_list over QQ",
        "replay_engine": "h10q._l13_irred8 exact Frobenius certificate modulo p",
        "reciprocal_structure": "P is palindromic on all 706 a=1 rows",
        "branches": dict(sorted(branch_counts.items())),
        "maximum_certificate_prime": max(certificate_counts),
        "certificate_prime_counts": {
            str(prime): count for prime, count in sorted(certificate_counts.items())
        },
        "uniform_conclusion": "OPEN; these are 706 finite specializations at a=1",
    }


def content_counterexample() -> tuple[dict[str, Any], dict[str, int]]:
    """A genuine all-rational-cell construction where c is not S-supported."""
    # w=3, z=3/11 is represented by unit 1/11 and has v_3(z)=1.
    data = construct_class(3, (1, 11), choose_a(3))
    row, counts = audit_class(data, "content-counterexample")
    assert row["a"] == "1" and row["q1"] == 19
    assert row["N"] == "1905120" and row["S"] == [2, 3, 5, 7]
    assert row["F_scale_c"] == "9072/15692141883605"
    assert row["F_scale_outside_S"] == "1/3138428376721"  # 11^-12
    assert row["F_scale_S_supported"] is False
    assert row["F_scale_square_outside_S"] is True
    row["literal_claim_verdict"] = "REFUTED: c is not S-supported"
    row["corrected_claim_verdict"] = "PROVED: c is S-supported modulo Q^x squares"
    return row, counts


def write_report(summary: dict[str, Any], counterexample: dict[str, Any]) -> None:
    cross = summary["mechanical_crosscheck"]
    legacy = summary["legacy_L13c_grid"]
    conditions = summary["conditions"]
    lines = [
        "# L20 — Schinzel admissibility of the uniform classes",
        "",
        "**Overall verdict — OPEN uniformly.** Conditions (b), (c), and (d), and positivity in (a), follow uniformly from the construction. Irreducibility in (a) is proved for every finite row replayed here but has no uniform proof for arbitrary `(w,z)`. Therefore this report does **not** close the chain and does **not** claim that Schinzel H implies the six-quantifier definition.",
        "",
        "There is also a literal correction to the normalization premise: for arbitrary rational `z`, the primitive scale `c` need not be supported on `S`. Its part outside `S` is always a twelfth power, which is sufficient for the parity/square-class argument, but it is not literally absent.",
        "",
        "## 1. Symbolic polynomial and exact content",
        "",
        "For the uniform square-branch construction put",
        "",
        "    A=1+4a^2,  tau=(1+2a^2)/A,  delta=-4a^4/A,",
        "    alpha=4a^4,  Z=z^3,  D=1-Z-a^2Z^2,  s=(a-1)/2,",
        "    B(t)=f(q1+Nt),  Ng(B)=16a^4B^2-A(B-1)^4.",
        "",
        "The exact `h10q._l10_P` specialization is",
        "",
        "    F(t)=P(B(t))",
        "        =16D^2A^2B(t)^4 + (4a^8Z^4/A)Ng(B(t))^2",
        "         -32A^3s^2D^2B(t)^5.",
        "",
        "If `P(b)=sum_i P_i b^i`, `b0=f*q1`, and",
        "",
        "    T_j=sum_{i=j}^8 binom(i,j) P_i b0^(i-j),",
        "",
        "then `F(t)=sum_j T_j(fN)^j t^j`. Let `F=cG`, where `G` is primitive integral and `c>0`. The exact valuation formula proved and replayed here is",
        "",
        "    c = (prod_{p in S} p^{v_p(P(b0))}) * d_out(z)^(-12),",
        "",
        "where `d_out(z)=prod_{p notin S, v_p(z)<0} p^{-v_p(z)}` is the part of the denominator of `z` outside `S`.",
        "",
        "At `p in S`, `_l10_exponent` and `p^k|N` make the constant Taylor coefficient uniquely minimal: every nonconstant coefficient has valuation at least `v_p(P(b0))+1`, and at `2` at least `v_2(P(b0))+3`. Hence `v_p(c)=v_p(P(b0))`. At `p notin S`, `fN` is a unit, so affine substitution preserves the coefficient ideal. If `v_p(z)>=0`, a coefficient of `P` is a unit; if `v_p(z)<0`, then `v_p(D)=2v_p(Z)` and the minimum coefficient valuation is `4v_p(Z)=12v_p(z)`. This proves the formula.",
        "",
        f"**Literal S-support is REFUTED for arbitrary rational cells.** The exact constructed class `(w,z)=(3,3/11)` has `a=1`, `q1=19`, `N=1905120`, `S={{2,3,5,7}}`, and `c={counterexample['F_scale_c']} = 2^4*3^4*7*5^(-1)*11^(-12)`. Thus the outside part is `11^-12`, not `1`. The corrected statement is uniform and sufficient: `c` is S-supported modulo a rational square (indeed modulo a twelfth power). On the canonical grid, denominators use only `3,5,7`, so literal S-support holds on {cross['classes']}/{cross['classes']} rows.",
        "",
        "Why the correction is sufficient downstream: at simultaneous Schinzel values `G(t)=R`, every valuation of `c` outside `S` is even. Since `alpha` is a square and the denominator `D^2A^2B^4` is a square, every outside place other than `R` and the moving prime `Q=q1+Nt` has even `v_p(x)` and `v_p(d)=0`, hence Hilbert symbol `+1`; `R` is the unique additional odd-valuation place. Thus the L19 parity/reciprocity implication needs the square class of `c` to be S-supported, not literal support.",
        "",
        "## 2. Condition (a): positive lead and irreducibility",
        "",
        "The constant and leading coefficients of `P` are both",
        "",
        "    P_0=P_8=4a^8 A Z^4 > 0.",
        "",
        "Consequently `F_8=P_8(fN)^8>0`, and division by `c>0` gives `G_8>0`: **positive leading coefficient is PROVED uniformly**.",
        "",
        "The substitution `b -> f(q1+Nt)` is an invertible affine change over `Q`, with inverse `t=(b-fq1)/(fN)`. A factorization of `F` pulls back to one of `P`, and multiplying by `c` changes none. Thus",
        "",
        "    G irreducible over Q  <=>  P irreducible over Q.",
        "",
        "There is an exact reciprocal reduction. At the L20 checkpoint it stopped short of the missing theorem; L22 later closes that theorem on the fixed canonical branch. Writing `C=32A^3s^2D^2`, direct coefficient comparison gives",
        "",
        "    P(b)=R(b)-C*b^5,       R_i=R_{8-i},",
        "    P(b)-b^8P(1/b)=C*(b^3-b^5).",
        "",
        f"This identity was exact on all {cross['palindromic_core_classes']} reconstructed classes. When `a=1` (equivalently `s=0`), `C=0`, so `P` is palindromic; this covers {cross['palindromic_classes']} constructed rows and all {legacy['rows']} legacy rows. Then `P(b)=b^4Q(b+1/b)` with",
        "",
        "    Q(u)=P0*u^4+P1*u^3+(P2-4P0)*u^2+(P3-3P1)*u",
        "         +(P4+2P0-2P2).",
        "",
        "For a root `u` of `Q`, `P` is irreducible exactly when `Q` is irreducible and `u^2-4` is not a square in `Q(u)`. At the L20 checkpoint neither assertion was proved uniformly in `z`; L22b now proves both on `a=1`, `tau=3/5=tau_dagger`, `delta=-4/5`, and L22a independently proves the full fixed-fiber theorem. For `s!=0`, near-symmetry alone still does not prove irreducibility; the L22a endpoint elimination does.",
        "",
        "The reciprocal pairing forces the `a=1` Galois group into the hyperoctahedral group `C2 wr S4` of order 384. The recorded `D8 wr C2`-consistent 2-group has order 128 and its observed cycle types do not contradict that containment; the artifact supplies no generic group certificate. Thus the structural and recorded Galois pictures are compatible, not a uniform proof.",
        "",
        f"**Constructed-branch scope (L21, superseded frontier).** L21 gives an exact reducible locus on other branches when `s=0` and `delta` is a rational square, and proves that the canonical branch avoids it because `delta=-4a^4/A<0`. The present replay found **0 reducible instances in {cross['tau_dagger_classes_checked']} exact constructed-branch classes** (the 353-cell grid plus `(3,3/11)`). Avoiding the L21 locus was necessary but not sufficient at that checkpoint; L22a now proves uniform fixed-fiber irreducibility on `tau_dagger`.",
        "",
        "The obvious construction prime does not itself prove irreducibility: since `f|Z` and `D=1 (mod f)`,",
        "",
        "    P(b) = 16D^2A^2 b^4 (1-2As^2b)  (mod f),",
        "",
        "which is reducible and has degree at most five. Moreover `v_f(P_0)=v_f(P_8)=4v_f(Z)>0` while `v_f(P_4)=0`, so the Newton polygon is not the single Eisenstein/Dumas segment needed for an octic irreducibility conclusion. The Galois-pattern observations in `data/l12_factorP.json` are finite-specialization evidence, not a generic Galois-group proof.",
        "",
        f"Under the explicit extra condition 'the primitive polynomial has an irreducible degree-8 reduction modulo some prime', irreducibility is PROVED. The kernel's exact Frobenius test supplied such a prime `<= {cross['maximum_irreducibility_certificate_prime']}` for **{cross['irreducible']}/{cross['classes']}** independently reconstructed uniform grid classes (`a` chosen from {list(A_POOL)}). It also replayed **{legacy['replayed_irreducible']}/{legacy['rows']}** L13c `(cell,branch)` rows whose original engine was {legacy['original_engine']}. These are theorem-level per-row certificates, but a finite grid is not the uniform theorem.",
        "",
        f"**Condition (a) current verdict:** positivity **PROVED uniformly** here; irreducibility **PROVED uniformly by L22a**. At the historical L20 checkpoint it was OPEN, while modular certificates already proved all {cross['classes']} constructed grid classes and {legacy['rows']} legacy grid rows.",
        "",
        "## 3. Condition (b): gcd(q1,N)=1",
        "",
        "Every prime divisor of `A` occurs in `supp(delta)` because `gcd(a,A)=1` and `delta=-4a^4/A`. Hence the construction",
        "",
        "    N=lcm(8, {p^{k_p}:p in S}, 4A)",
        "",
        "has prime support exactly `S`. The CRT residue for `q1` is a unit modulo `4*A*prod(S)`, so it is a unit at every prime dividing `N`. Therefore `gcd(q1,N)=1` uniformly.",
        "",
        f"**Condition (b) verdict: PROVED uniformly.** Mechanical check: {cross['gcd_q1_N_one']}/{cross['classes']}, 0 mismatches.",
        "",
        "## 4. Condition (c): no fixed prime divisor",
        "",
        "Let `L(t)=q1+Nt`. If `p in S`, then `p|N`, `q1` is a unit mod `p`, and the Taylor congruence below gives `G(t)=G(0) != 0 (mod p)`. Thus neither factor ever vanishes modulo `p`. This decides every prime in `S`, including `2,3,5,7`.",
        "",
        "If `p notin S`, then `p` does not divide `N`, so `L` has exactly one root modulo `p`. Were `L(t)G(t)` zero for every residue, the nonzero polynomial `G mod p` (nonzero because `G` is primitive) would vanish at the other `p-1` residues. Since `deg G=8`, this forces `p-1<=8`, hence `p<=9`. But all primes `2,3,5,7` are in `S`, a contradiction. This is a uniform proof, not an instance scan.",
        "",
        f"As an independent exact check, the degree-nine product's values at `t=0,...,9` have gcd 1 on **{cross['fixed_divisor_one']}/{cross['classes']}** rows. Ten values are exact by the finite-difference lemma for an integral degree-nine polynomial.",
        "",
        "**Condition (c) verdict: PROVED uniformly.**",
        "",
        "## 5. Condition (d): S-units and frozen symbols",
        "",
        "Here 'G is an S-unit' follows the project terminology: `G(t)` is a p-adic unit for every `p in S`; it does not mean the conventional global group of rational S-units.",
        "",
        "The Taylor-exponent inequalities and the content formula give the polynomial congruences",
        "",
        "    G(t) = G(0) != 0 (mod p)       for odd p in S,",
        "    G(t) = G(0) != 0 (mod 8)       for p=2,",
        "    q1+Nt = q1 != 0 (mod p), and mod 8 at p=2.",
        "",
        "They hold coefficientwise, hence for every integer `t`. Equivalently, `P(B(t))/P(b0)` and `B(t)/b0` lie in `1+pZ_p` (in `1+8Z_2` at 2), so both are local squares. Since",
        "",
        "    x(t)/x(0) = [P(B(t))/P(b0)] * [b0/B(t)]^4,",
        "    d(t)/d(0) = B(t)/b0,",
        "",
        "the square classes and all frozen Hilbert symbols are constant. This is exactly the L10/L12 Taylor-exponent lemma; no sampling is needed for the universal quantifier in `t`.",
        "",
        f"**Condition (d) verdict: PROVED uniformly.** The replay checked the coefficient congruences on {cross['S_unit_congruence_classes']}/{cross['classes']} classes and additionally made {cross['frozen_symbol_checks']} exact sampled-member symbol comparisons, with 0 mismatches.",
        "",
        "## 6. Mechanical scope and final chain verdict",
        "",
        f"The script independently reconstructed **{cross['classes']}** classes across **{cross['cells']}** canonical cells and **{cross['distinct_w']}** target primes, choosing `a` by `(1+4a^2|w)=-1` and searching the explicit CRT subsystem. It matched all **{cross['l19_materialized_matches']}/{cross['l19_materialized_rows']}** materialized L19 clean-class records and the L19 353-cell aggregate totals. Content, normalized polynomial, modular irreducibility, gcd, fixed divisor, S-unit congruences, and frozen symbols had **{cross['total_mismatches']} mismatches**. The explicit `(3,3/11)` content counterexample was audited through the same path.",
        "",
        "**Honest final sentence (current).** L20 proves positivity and conditions (b)–(d), and its historical irreducibility gap is now closed by the fixed-fiber theorem L22a (`l22_elimination.py`, `data/l22_elimination.jsonl`). Thus classical Schinzel H implies the six-quantifier definition through L19–L22. The conclusion remains CONDITIONAL because Schinzel H is unproved. The literal premise 'c is S-supported' must still be replaced by the proved sufficient statement 'the square class of c is S-supported' when `z` has denominator primes outside `S`.",
        "",
        f"Labels: (a) `{conditions['a_irreducibility']['label']}` uniformly; (b) `{conditions['b_gcd']['label']}`; (c) `{conditions['c_no_fixed_divisor']['label']}`; (d) `{conditions['d_S_units_and_symbols']['label']}`. Refusals: {summary['refusals']} (none used as evidence). Verification wall-clock: **{summary['wall_seconds']:.3f} s**. Pacing: {summary['pacing']['candidates']} loop steps, {summary['pacing']['sleeps']} sleeps of {PACE_SECONDS} s.",
        "",
        "Artifacts: `math/h10q/l20_admissible.py`, `math/h10q/data/l20_admissible.jsonl`, `/tmp/l20_admissible.md`.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    materialized, l19_summary = l19_materialized()
    grid = sorted(h10q._l9_grid())
    assert len(grid) == 353

    rows: list[dict[str, Any]] = []
    total_counts = Counter()
    a_counts = Counter()
    certificate_counts = Counter()
    l19_matches = 0
    attempts = 0
    maximum_q1 = 0
    for w, unit in grid:
        a = choose_a(w)
        data = construct_class(w, unit, a)
        row, counts = audit_class(data, "constructed-grid-class")
        rows.append(row)
        total_counts.update(counts)
        a_counts[a] += 1
        certificate_counts[row["G_irreducibility_certificate_prime"]] += 1
        attempts += data["attempts"]
        maximum_q1 = max(maximum_q1, data["q1"])

    # The 103 canonical materializations may use a different numerator prime
    # f than the uniform f=w grid replay.  Rebuild each materialized class with
    # its own (a,f), rather than conflating those two valid constructions.
    for (w, unit), prior in sorted(materialized.items()):
        materialized_data = construct_class(
            w,
            unit,
            int(prior["a"]),
            int(prior["f"]),
        )
        compare_l19(materialized_data, prior)
        l19_matches += 1

    assert l19_matches == len(materialized) == 163
    uniform_l19 = l19_summary["uniform_theorem"]
    assert uniform_l19["grid_cells"] == len(rows)
    assert uniform_l19["grid_certificates_ok"] == len(rows)
    assert {int(a): count for a, count in uniform_l19["bounded_a_counts_by_cell"].items()} == dict(a_counts)
    assert uniform_l19["maximum_constructed_clean_q1"] == maximum_q1
    assert uniform_l19["constructed_q1_search_candidates"] == attempts
    assert all(row["F_scale_S_supported"] for row in rows)

    counterexample, counter_counts = content_counterexample()
    total_counts.update(counter_counts)
    legacy = replay_l13c()

    elapsed = time.perf_counter() - started
    frozen_checks = sum(row["frozen_symbol_checks"] for row in rows) + counterexample["frozen_symbol_checks"]
    summary = {
        "type": "summary",
        "label": "PROVED",
        "verdict": "Schinzel admissibility conditions (a)-(d) are proved: L20 proves positivity and (b)-(d), while L22a proves uniform fixed-fiber irreducibility; literal S-support of c is refuted but S-support modulo squares is proved and sufficient",
        "conditions": {
            "a_irreducibility": {
                "label": "PROVED",
                "l20_checkpoint_label": "OPEN",
                "positive_leading_coefficient": "PROVED_UNIFORMLY",
                "affine_equivalence": "G irreducible iff P irreducible",
                "reciprocal_structure": "P=R-32*A^3*s^2*D^2*b^5 with R palindromic",
                "a_equals_1_reduction": "L22b proves P irreducible for every rational nonzero Z on a=1,tau=3/5=tau_dagger,delta=-4/5",
                "palindromic_constructed_classes": a_counts[1],
                "palindromic_legacy_classes": legacy["rows"],
                "L21_known_reducible_locus": "s=0 and delta a rational square (other branches)",
                "constructed_branch_avoids_L21_locus": "PROVED because delta=-4*a^4/A<0",
                "tau_dagger_classes_checked": len(rows) + 1,
                "tau_dagger_reducible_instances": 0,
                "uniform_irreducibility": "PROVED_BY_L22A_FIXED_FIBER_THEOREM",
                "uniform_irreducibility_producer": "l22_elimination.py",
                "uniform_irreducibility_artifact": "data/l22_elimination.jsonl",
                "explicit_extra_condition": "primitive P has irreducible degree-8 reduction modulo some prime",
                "constructed_grid_under_extra_condition": len(rows),
                "constructed_grid_total": len(rows),
                "legacy_L13c_under_extra_condition": legacy["replayed_irreducible"],
                "legacy_L13c_total": legacy["rows"],
            },
            "b_gcd": {
                "label": "PROVED",
                "verdict": "PROVED_UNIFORMLY",
                "proof": "supp(N)=S and q1 is a unit modulo 4*A*prod(S)",
            },
            "c_no_fixed_divisor": {
                "label": "PROVED",
                "verdict": "PROVED_UNIFORMLY",
                "proof": "p in S: both factors are units; p outside S: one linear root plus at most eight G-roots forces p<=9, but 2,3,5,7 lie in S",
            },
            "d_S_units_and_symbols": {
                "label": "PROVED",
                "verdict": "PROVED_UNIFORMLY",
                "proof": "Taylor exponent gives coefficientwise mod-p/mod-8 congruences and constant local square classes",
                "terminology": "S-unit means a p-adic unit at every p in S",
            },
        },
        "normalization": {
            "label": "PROVED",
            "exact_formula": "c=(prod_{p in S} p^vp(P(b0),p))*d_out(z)^(-12)",
            "literal_S_support_uniform": "REFUTED",
            "literal_S_support_on_canonical_grid": len(rows),
            "canonical_grid_rows": len(rows),
            "square_class_S_supported_uniformly": "PROVED",
            "counterexample": {
                "cell": counterexample["cell"],
                "z": counterexample["z"],
                "a": counterexample["a"],
                "q1": counterexample["q1"],
                "N": counterexample["N"],
                "S": counterexample["S"],
                "c": counterexample["F_scale_c"],
                "outside_S": counterexample["F_scale_outside_S"],
            },
        },
        "mechanical_crosscheck": {
            "label": "PROVED",
            "classes": len(rows),
            "cells": len(rows),
            "distinct_w": len({row["cell"][0] for row in rows}),
            "a_counts": {str(a): count for a, count in sorted(a_counts.items())},
            "q1_search_candidates": attempts,
            "maximum_q1": maximum_q1,
            "l19_materialized_rows": len(materialized),
            "l19_materialized_matches": l19_matches,
            "l19_aggregate_match": True,
            "content_formula_matches": len(rows),
            "literal_scale_S_supported": len(rows),
            "palindromic_core_classes": len(rows),
            "palindromic_classes": a_counts[1],
            "tau_dagger_classes_checked": len(rows) + 1,
            "tau_dagger_reducible_instances": 0,
            "irreducible": len(rows),
            "maximum_irreducibility_certificate_prime": max(certificate_counts),
            "irreducibility_certificate_counts": {
                str(prime): count for prime, count in sorted(certificate_counts.items())
            },
            "gcd_q1_N_one": len(rows),
            "fixed_divisor_one": len(rows),
            "S_unit_congruence_classes": len(rows),
            "frozen_symbol_checks": frozen_checks,
            "content_counterexample_checked": 1,
            "total_mismatches": sum(total_counts.values()),
        },
        "legacy_L13c_grid": legacy,
        "chain_consequence": {
            "label": "PROVED",
            "schinzel_H_implies_six_quantifier_definition": True,
            "conditional_on": "classical Schinzel H",
            "residual_gap": None,
            "unconditional_open_input": "classical Schinzel H is unproved",
            "normalization_correction": "replace c S-supported by squareclass(c) S-supported for arbitrary rational z",
        },
        "refusals": 0,
        "refusals_are_evidence": False,
        "wall_seconds": elapsed,
        "pacing": {
            "candidates": PACER.candidates,
            "every": PACE_EVERY,
            "sleep_seconds": PACE_SECONDS,
            "sleeps": PACER.sleeps,
        },
        "sources": {
            "kernel_authority": "math/h10q/h10q.py",
            "construction_comparison": "math/h10q/data/l19_classexist.jsonl",
            "legacy_irreducibility_input": "math/h10q/data/l12_factorP.json",
            "script": "math/h10q/l20_admissible.py",
        },
    }
    assert summary["mechanical_crosscheck"]["classes"] >= 100
    assert summary["mechanical_crosscheck"]["total_mismatches"] == 0

    meta = {
        "type": "meta",
        "label": "PROVED",
        "schema": "l20-admissible-v2",
        "source_script": "math/h10q/l20_admissible.py",
        "kernel_authority": "math/h10q/h10q.py",
        "refusals_are_evidence": False,
        "symbolic_F": "16*D^2*A^2*B^4+(4*a^8*Z^4/A)*(16*a^4*B^2-A*(B-1)^4)^2-32*A^3*s^2*D^2*B^5; B=f*(q1+N*t)",
        "content_formula": "c=(prod_{p in S}p^vp(P(b0),p))*d_out(z)^(-12)",
    }
    output_rows = [meta, *rows, counterexample, summary]
    with OUT.open("w", encoding="utf-8") as handle:
        for row in output_rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(summary, counterexample)

    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
    print(f"wrote {OUT} ({len(output_rows)} rows)")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
