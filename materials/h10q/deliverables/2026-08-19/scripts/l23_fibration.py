#!/usr/bin/env python3
"""Exact fixed-cell conic-bundle audit on the canonical square branch.

This script is deliberately theorem-matching rather than a search.  It proves
all algebraic identities and local residue calculations with standard-library
exact arithmetic, replays the 293 checked Capell obstructions, and writes a
source-ready report.  It does not promote a finite replay to a uniform theorem.

Replay from the workspace root:

    nice -n 19 python3 math/h10q/l23_fibration.py
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from fractions import Fraction as F
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l23_fibration.jsonl"
REPORT = Path("/tmp/l23_fibration.md")
STEP_II = HERE / "data" / "l17_stepii.jsonl"
PACE_EVERY = 64
PACE_SECONDS = 0.003

Poly = list[F]
Monomial = tuple[int, int, int]  # exponents of (a,Z,b)
MPoly = dict[Monomial, F]


def trim(poly: Poly) -> Poly:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer or [F(0)]


def p_add(*polys: Poly) -> Poly:
    size = max((len(poly) for poly in polys), default=1)
    answer = [F(0)] * size
    for poly in polys:
        for index, coefficient in enumerate(poly):
            answer[index] += coefficient
    return trim(answer)


def p_scale(poly: Poly, scalar: int | F) -> Poly:
    scalar = F(scalar)
    return trim([scalar * coefficient for coefficient in poly])


def p_mul(left: Poly, right: Poly) -> Poly:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            answer[i + j] += first * second
    return trim(answer)


def p_pow(poly: Poly, exponent: int) -> Poly:
    assert exponent >= 0
    answer = [F(1)]
    base = poly
    power = exponent
    while power:
        if power & 1:
            answer = p_mul(answer, base)
        base = p_mul(base, base)
        power //= 2
    return answer


def p_derivative(poly: Poly) -> Poly:
    if len(poly) == 1:
        return [F(0)]
    return trim([F(index) * poly[index] for index in range(1, len(poly))])


def m_const(value: int | F) -> MPoly:
    value = F(value)
    return {(0, 0, 0): value} if value else {}


def m_var(index: int) -> MPoly:
    exponent = [0, 0, 0]
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def m_add(*polys: MPoly) -> MPoly:
    answer: MPoly = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            answer[monomial] = answer.get(monomial, F(0)) + coefficient
            if not answer[monomial]:
                del answer[monomial]
    return answer


def m_scale(poly: MPoly, scalar: int | F) -> MPoly:
    scalar = F(scalar)
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in poly.items()
        if scalar * coefficient
    }


def m_mul(left: MPoly, right: MPoly) -> MPoly:
    answer: MPoly = {}
    for first_monomial, first_coefficient in left.items():
        for second_monomial, second_coefficient in right.items():
            monomial = tuple(
                first_monomial[index] + second_monomial[index]
                for index in range(3)
            )
            answer[monomial] = (
                answer.get(monomial, F(0))
                + first_coefficient * second_coefficient
            )
            if not answer[monomial]:
                del answer[monomial]
    return answer


def m_pow(poly: MPoly, exponent: int) -> MPoly:
    assert exponent >= 0
    answer = m_const(1)
    base = poly
    power = exponent
    while power:
        if power & 1:
            answer = m_mul(answer, base)
        base = m_mul(base, base)
        power //= 2
    return answer


def frac_text(value: int | F) -> str:
    value = F(value)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def mod_fraction(value: F, prime: int) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    assert denominator
    return numerator * pow(denominator, -1, prime) % prime


def eval_mod(poly: Poly, value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(poly):
        answer = (
            answer * value + mod_fraction(coefficient, prime)
        ) % prime
    return answer


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def legendre(value: int, prime: int) -> int:
    value %= prime
    if not value:
        return 0
    result = pow(value, (prime - 1) // 2, prime)
    return -1 if result == prime - 1 else result


def canonical_polynomial(a: F, z: F) -> dict[str, Any]:
    """Return the exact canonical-branch kernel polynomial in b."""
    A = 1 + 4 * a * a
    s = (a - 1) / 2
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    tau = (A + 1) / (2 * A)
    delta = 1 - A * tau * tau
    alpha = -delta * A
    assert delta == -4 * a**4 / A
    assert alpha == 4 * a**4 == (2 * a * a) ** 2

    b_minus_one = [F(-1), F(1)]
    Ng = p_add(
        p_scale(p_pow(b_minus_one, 4), -A),
        [F(0), F(0), 16 * a**4],
    )
    P = p_add(
        [F(0)] * 4 + [16 * D * D * A * A],
        p_scale(p_mul(Ng, Ng), -delta * a**4 * Z**4),
        [F(0)] * 5 + [-32 * A**3 * s * s * D * D],
    )
    assert len(P) == 9

    H = p_scale(P, A / 4)
    X = p_scale(Ng, a**4 * Z**2)
    Y = [F(0), F(0), 2 * A * D]
    L = [F(1), -2 * A * s * s]
    assert H == p_add(p_mul(X, X), p_scale(p_mul(L, p_mul(Y, Y)), A))

    endpoint = 4 * A * a**8 * Z**4
    assert P[0] == P[8] == endpoint
    assert H[0] == H[8] == (A * a**4 * Z**2) ** 2
    return {
        "a": a,
        "z": z,
        "A": A,
        "s": s,
        "Z": Z,
        "D": D,
        "tau": tau,
        "delta": delta,
        "alpha": alpha,
        "Ng": Ng,
        "P": P,
        "H": H,
    }


def symbolic_identity_record() -> dict[str, Any]:
    """Verify A*P=4H and H=X^2+A*L*Y^2 in Q[a,Z,b]."""
    one = m_const(1)
    a = m_var(0)
    Z = m_var(1)
    b = m_var(2)
    A = m_add(one, m_scale(m_pow(a, 2), 4))
    s = m_scale(m_add(a, m_scale(one, -1)), F(1, 2))
    D = m_add(
        one,
        m_scale(Z, -1),
        m_scale(m_mul(m_pow(a, 2), m_pow(Z, 2)), -1),
    )
    Ng = m_add(
        m_scale(m_mul(m_pow(a, 4), m_pow(b, 2)), 16),
        m_scale(m_mul(A, m_pow(m_add(b, m_scale(one, -1)), 4)), -1),
    )
    X = m_mul(m_mul(m_pow(a, 4), m_pow(Z, 2)), Ng)
    Y = m_scale(m_mul(m_mul(A, D), m_pow(b, 2)), 2)
    L = m_add(
        one,
        m_scale(m_mul(m_mul(A, m_pow(s, 2)), b), -2),
    )
    H = m_add(m_pow(X, 2), m_mul(m_mul(A, L), m_pow(Y, 2)))
    direct_H = m_add(
        m_mul(m_mul(m_pow(a, 8), m_pow(Z, 4)), m_pow(Ng, 2)),
        m_scale(m_mul(m_mul(m_pow(A, 3), m_pow(D, 2)), m_pow(b, 4)), 4),
        m_scale(
            m_mul(
                m_mul(m_mul(m_pow(A, 4), m_pow(s, 2)), m_pow(D, 2)),
                m_pow(b, 5),
            ),
            -8,
        ),
    )
    assert H == direct_H
    P_bar = m_scale(H, 4)
    return {
        "type": "canonical-birational-model",
        "label": "PROVED",
        "canonical_branch": (
            "tau_dagger=(A+1)/(2A), delta=-4a^4/A, "
            "alpha=-delta*A=4a^4=(2a^2)^2"
        ),
        "surface": "U^2-2b*V^2=P_{a,Z}(b)*W^2",
        "threefold_cleared": "A*(U^2-2b*V^2)=4H(a,Z,b)*W^2",
        "function_field_quaternion": "(P_{a,Z}(b),2b)",
        "definitions": {
            "A": "1+4a^2",
            "D": "1-Z-a^2Z^2",
            "s": "(a-1)/2",
            "Ng": "16a^4b^2-A(b-1)^4",
            "P": (
                "16D^2A^2b^4+(4a^8/A)Z^4Ng^2"
                "-32A^3s^2D^2b^5"
            ),
            "H": "(A/4)P=X^2+A*L*Y^2",
            "X": "a^4Z^2Ng",
            "Y": "2ADb^2",
            "L": "1-2As^2b",
        },
        "symbolic_monomials": {
            "H": len(H),
            "A_times_P": len(P_bar),
        },
        "dimension": {
            "fixed_Z_variable_a": 3,
            "fixed_a_and_Z": 2,
        },
    }


def geometry_record() -> dict[str, Any]:
    return {
        "type": "conic-bundle-geometry",
        "label": "PROVED",
        "premises": [
            "a*Z*D != 0",
            "P is irreducible of degree 8 (the proved L22 premise)",
        ],
        "discriminant_squareclass": "2b*P(b)",
        "closed_discriminant_points": [
            {"point": "b=0", "degree": 1, "splitting_class": "A"},
            {
                "point": "P(b)=0",
                "degree": 8,
                "residue_field": "K=Q(theta)",
                "splitting_class": "2theta",
            },
            {"point": "b=infinity", "degree": 1, "splitting_class": "A"},
        ],
        "geometric_degenerate_fibres": 10,
        "closed_factor_degrees": [1, 8, 1],
        "natural_model_K_square": -2,
        "K_square_formula": "K_X^2=8-r with r=10 on the natural regular model",
        "relative_minimal_actual_nonsplit_case": True,
        "split_case_note": (
            "if 2theta is square, the eight split fibres remain on the natural "
            "model but their chosen components are contractible"
        ),
        "endpoint_identities": {
            "P0_equals_P8": "4*A*a^8*Z^4=A*(2a^4Z^2)^2",
            "infinity_residue": (
                "(2b)^8/P(b) at infinity =256/P8, squareclass A"
            ),
            "root_residue_norm": (
                "Norm_K/Q(2theta)=2^8*P0/P8=256=16^2"
            ),
        },
        "rank": {
            "if_2theta_nonsquare_in_K": 10,
            "if_2theta_square_in_K": 2,
            "definition": "degree of the non-split locus, not number of geometric singular fibres",
        },
    }


def brauer_record() -> dict[str, Any]:
    return {
        "type": "vertical-brauer-group",
        "label": "PROVED",
        "method": (
            "Faddeev residues at {0,P,infinity}, the conic-kernel class beta=(P,2b), "
            "and purity for a relatively minimal conic bundle"
        ),
        "residue_vectors": {
            "coordinates": ["0", "P", "infinity"],
            "beta_when_2theta_nonsquare": [1, 1, 1],
            "beta_when_2theta_square": [1, 0, 1],
            "A_b": [1, 0, 1],
            "faddeev_relation": "e_0=e_infinity; Norm(2theta)=256 is square",
        },
        "A_nonsquare_reason": "a in 1+2Z_2 gives A=1+4a^2 ==5 mod 8",
        "answers": {
            "2theta_nonsquare": {
                "BrX_mod_BrQ": "Z/2",
                "generator": "(A,b)",
                "dimension_F2": 1,
            },
            "2theta_square": {
                "BrX_mod_BrQ": "0",
                "reason": "(A,b) has the same residue vector as beta, which vanishes on X",
                "dimension_F2": 0,
            },
        },
        "scope": "for the smooth proper relatively minimal surface X_{a,Z}",
    }


def theorem_match_record() -> dict[str, Any]:
    return {
        "type": "unconditional-theorem-match",
        "label": "OPEN",
        "strongest_exact_implication": {
            "theorem": (
                "Harpaz-Wei-Wittenberg, Rational points on fibrations with few "
                "non-split fibres, Theorem 1.3(i)"
            ),
            "source": "https://arxiv.org/abs/2109.03547",
            "statement_used": (
                "for rank at most 2, rational points are dense in the full "
                "Brauer-Manin set"
            ),
            "match_if": "2theta is a square in K=Q[b]/(P)",
            "then": (
                "the octic fibre is split, rank=2 and BrX/BrQ=0; the aligned "
                "local certificates give an adelic point, and weak approximation "
                "can impose b in Z_2^times"
            ),
            "logical_status": "PROVED implication",
        },
        "proved_rank_10_scopes": {
            "odd_target_valuation": (
                "for v_w(z) odd, every admissible a: the left Newton edge "
                "forces a place with odd valuation of 2theta"
            ),
            "all_cells_refined": (
                "after the free choice a!=1 mod w and the proved L22 HIT "
                "selection inside that progression"
            ),
            "BrX_mod_BrQ": "Z/2 generated by (A,b)",
            "first_failed_hypothesis": "rank<=2 (the rank is 10)",
            "rank_3_clauses": "inapplicable for the same numerical reason",
            "verdict": "no unconditional Hasse/weak-approximation theorem found that decides these fixed surfaces",
        },
        "unresolved_preexisting_scope": (
            "when v_w(z) is even and w divides s, the target-prime argument is "
            "OPEN unless an exact per-row residue certificate is supplied"
        ),
        "conditional_fallback": {
            "theorem": "Harpaz-Wei-Wittenberg Theorem 1.4",
            "cyclic_splitting_fields": True,
            "first_unproved_input": (
                "homogeneous Schinzel (HH_1) for the irreducible octic closed point P"
            ),
            "status": "CONDITIONAL",
        },
    }


def comparison_records() -> list[dict[str, Any]]:
    return [
        {
            "type": "comparison-chatelet",
            "label": "PROVED",
            "source": "https://arxiv.org/abs/2209.08949",
            "theorem_scope": (
                "classical Chatelet: N_{L/Q}(x)=f(t) with one fixed quadratic "
                "field L and separable f of degree 3 or 4"
            ),
            "first_failed_hypothesis": (
                "the quadratic algebra Q(b)(sqrt(2b)) varies with b; it is not a "
                "constant number-field extension"
            ),
            "additional_obstruction": (
                "the relatively minimal bundle has 10, not at most 4, geometric "
                "degenerate fibres (K_X^2=-2)"
            ),
            "non_birational_base_change": (
                "b=u^2/2 is a degree-2 base change that splits the generic conic, "
                "not a birational change of the Q-base"
            ),
            "applies": False,
        },
        {
            "type": "comparison-linear-split-fibres",
            "label": "PROVED",
            "source": "https://arxiv.org/abs/1304.3333",
            "theorem_scope": (
                "Harpaz-Skorobogatov-Wittenberg/Green-Tao-Ziegler applies when all "
                "degenerate geometric fibres are defined over Q"
            ),
            "first_failed_hypothesis": (
                "the eight roots of irreducible P form one degree-8 closed point, "
                "so those geometric fibres are not Q-rational"
            ),
            "applies": False,
        },
        {
            "type": "comparison-beta-sieve",
            "label": "PROVED",
            "source": "https://arxiv.org/abs/2209.08949",
            "theorem_scope": (
                "Shute Theorem 1.4 proves the Harpaz-Wittenberg conjecture for "
                "residue fields k_i/Q of degree at most 2, or at most 3 under "
                "extra sieve hypotheses"
            ),
            "first_failed_hypothesis": "[K:Q]=deg(P)=8>3",
            "fixed_norm_form_mismatch": (
                "the direct norm-form theorems also require a fixed number field; "
                "sqrt(2b) is not fixed"
            ),
            "applies": False,
        },
        {
            "type": "comparison-browning-schindler",
            "label": "PROVED",
            "source": "https://arxiv.org/abs/1509.07744",
            "theorem_scope": (
                "Browning-Schindler Theorem 1.1 treats rank at most 3 over Q, "
                "with at least one non-split fibre above a rational point"
            ),
            "first_failed_hypothesis": "rank=10>3",
            "applies": False,
        },
        {
            "type": "comparison-random-conic-bundles",
            "label": "PROVED",
            "source": "https://arxiv.org/abs/2604.07047",
            "theorem_scope": (
                "Frei-Sofos Theorem 1.1 is a density-one statement as coefficients "
                "vary, for arbitrary fixed degree patterns"
            ),
            "first_failed_hypothesis": (
                "it supplies no pointwise conclusion for this fixed two-parameter "
                "coefficient family or any specified member"
            ),
            "applies": False,
        },
    ]


def phi_record() -> dict[str, Any]:
    return {
        "type": "phi-local-openness",
        "label": "PROVED",
        "identity": (
            "Phi_1^{2} on this base is a in 1+2Z_2 and b in Z_2^times"
        ),
        "fixed_a_condition": "already satisfied by the chosen odd a",
        "b_condition": "Z_2^times is an open subset of Q_2",
        "weak_approximation_compatible": True,
        "not_Zariski_open": True,
        "important_distinction": (
            "b in Phi is a local open condition; demanding b=eps*f*Q with Q prime "
            "in a fixed progression is not.  A surface weak-approximation theorem "
            "would only need the former, which is the original existential condition."
        ),
    }


def target_prime_replay() -> dict[str, Any]:
    """Finite replay of the target-prime lemma; the proof itself is symbolic."""
    primes = [prime for prime in range(3, 100) if is_prime(prime)]
    checked = 0
    qualifying = 0
    for w in primes:
        good = [
            residue
            for residue in range(w)
            if legendre(1 + 4 * residue * residue, w) == -1
        ]
        expected = (w - legendre(-1, w)) // 2
        assert len(good) == expected
        assert good and all(residue for residue in good)
        choices = [residue for residue in good if residue != 1]
        assert choices
        a_residue = choices[0]
        # Choose the odd CRT lift of a_residue and take z=w.  The theorem only
        # uses the displayed residue classes, not this finite representative.
        a_value = a_residue if a_residue % 2 else a_residue + w
        data = canonical_polynomial(F(a_value), F(w))
        A = mod_fraction(data["A"], w)
        s = mod_fraction(data["s"], w)
        assert A and s and legendre(A, w) == -1
        beta = pow(2 * A * s * s % w, -1, w)
        P = data["P"]
        assert eval_mod(P, beta, w) == 0
        derivative = eval_mod(p_derivative(P), beta, w)
        assert derivative
        assert legendre(2 * beta, w) == -1
        # Directly replay P mod w =16*A^2*b^4*(1-2*A*s^2*b).
        for b_value in range(w):
            expected_value = (
                16 * A * A * pow(b_value, 4, w)
                * (1 - 2 * A * s * s * b_value)
            ) % w
            assert eval_mod(P, b_value, w) == expected_value
            checked += 1
        qualifying += len(good)
        if checked % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
    return {
        "type": "target-prime-capell-no-go",
        "label": "PROVED",
        "hypotheses": [
            "w is an odd prime",
            "v_w(Z)>0",
            "a,A=1+4a^2,s=(a-1)/2 are w-units",
            "(A|w)=-1",
        ],
        "reduction": "P(b)=16*A^2*b^4*(1-2*A*s^2*b) mod w",
        "simple_root": "beta=(2*A*s^2)^(-1)",
        "root_squareclass": "2*beta=1/(A*s^2), hence (2*beta|w)=-1",
        "conclusion": (
            "the degree-one w-adic place theta->beta makes 2theta nonsquare; "
            "therefore 2theta is not a square in K"
        ),
        "choice_of_a": (
            "the residues with (1+4a^2|w)=-1 occur in nonzero +/- pairs; "
            "choose one not congruent to 1 mod w, then CRT with a odd"
        ),
        "qualifying_residue_count": "(w-(-1|w))/2",
        "odd_target_valuation_case": {
            "hypothesis": "e=v_w(z) is odd; no condition on s",
            "newton_edge": "(0,12e) to (4,0), slope -3e",
            "residual": "endpoint quartic is separable because w is odd",
            "conclusion": (
                "a place above w has ord(theta)=3e odd, so 2theta is nonsquare"
            ),
        },
        "uniform_scope": (
            "every fixed cell after refining the existing L20 choice by "
            "a not congruent to 1 mod w and applying the proved L22 HIT "
            "selection inside that progression; no new branch parameter"
        ),
        "residual_scope": (
            "for a pre-existing class with v_w(z) even and w|s this target-prime "
            "proof is OPEN; only the separately replayed per-row certificates apply"
        ),
        "finite_replay": {
            "label": "EVIDENCE",
            "primes": len(primes),
            "qualifying_residues": qualifying,
            "point_evaluations": checked,
            "maximum_prime": max(primes),
        },
    }


def capell_record() -> dict[str, Any]:
    return {
        "type": "capell-square-locus",
        "label": "PROVED",
        "premise": "P is irreducible of degree 8 and theta=b mod P",
        "equivalences": [
            "2theta is a square in K=Q(theta)",
            "P(u^2/2) factors over Q as two degree-8 factors exchanged by u->-u",
            (
                "H=(A/4)P has a representation H=E(b)^2-2b*O(b)^2 with "
                "E in Q[b], deg E<=4, and O in Q[b], deg O<=3"
            ),
        ],
        "exact_arithmetic_locus": (
            "for H=sum(h_k b^k), it is the projection of h_k="
            "sum_{i+j=k}e_i e_j-2 sum_{i+j=k-1}o_i o_j (0<=k<=8), "
            "with E=sum_{0..4}e_i b^i and O=sum_{0..3}o_i b^i"
        ),
        "endpoint_normalization": (
            "H_0=H_8=(A*a^4*Z^2)^2, so the scalar in a conjugate factorization "
            "is a rational square and can be absorbed into E,O"
        ),
        "intersection_with_refined_L20_locus": "EMPTY by the target-prime lemma",
        "warning": (
            "theta is the algebraic root of P; the rational Phi variable b is not "
            "being substituted into this field criterion"
        ),
    }


def twist_record() -> dict[str, Any]:
    return {
        "type": "constant-twist-audit",
        "label": "PROVED",
        "twisted_capell": (
            "2theta/j square in K iff P(j*u^2/2) has the conjugate degree-8 factorization"
        ),
        "effect_at_simple_value_primes": (
            "it changes the odd-valuation Hilbert sign to (j|R); it does not make "
            "the octic fibre split and does not force that sign to +1"
        ),
        "norm_identity_no_go": (
            "if also P=c*(E^2-j*O^2) with c in Q^times, then at theta, "
            "E(theta)/O(theta)=sqrt(j) lies in K (O(theta)!=0 by deg O<8). "
            "Thus j is a square in K, and 2theta/j square would imply 2theta square, "
            "contradicting the target-prime lemma"
        ),
        "without_norm_identity": (
            "requiring every odd value-prime R to satisfy (j|R)=+1 is the same "
            "Chebotarev half-sieve condition, not a Dirichlet-only replacement"
        ),
    }


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def replay_step_ii() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Replay every simple residue witness proving 2theta nonsquare per row."""
    source_rows = load_jsonl(STEP_II)
    class_rows = [row for row in source_rows if row.get("type") == "class"]
    source_summary = next(
        row for row in source_rows if row.get("type") == "summary"
    )
    assert len(class_rows) == source_summary["n_classes"] == 293
    proof_rows: list[dict[str, Any]] = []
    primes: list[int] = []
    for index, row in enumerate(class_rows, start=1):
        w, unit = row["cell"]
        params = row["params"]
        a = F(params["a"])
        z = F(w) * F(unit[0], unit[1])
        data = canonical_polynomial(a, z)
        P = data["P"]
        alpha = data["alpha"]
        witness = row["proof"]["simple_bad_root_witness"]
        prime = int(witness["p"])
        residue = int(witness["r"])
        assert is_prime(prime) and witness["prime_proved"]
        eps = int(params["eps"])
        f = int(params["f"])
        q1 = int(params["q1"])
        modulus = int(params["N"])
        base = eps * f * q1
        step = eps * f * modulus
        b_mod = (base + step * residue) % prime
        value = eval_mod(P, b_mod, prime)
        derivative = (
            eval_mod(p_derivative(P), b_mod, prime) * (step % prime)
        ) % prime
        bad_arg = mod_fraction(2 * alpha, prime) * b_mod % prime
        assert value == int(witness["F_mod_p"]) == 0
        assert derivative == int(witness["dFdk_mod_p"]) != 0
        assert bad_arg == int(witness["bad_arg_mod_p"]) != 0
        assert legendre(bad_arg, prime) == -1
        assert row["proof"]["bad_unit_nonsquare_in_Qtheta"]
        primes.append(prime)
        proof_rows.append({
            "type": "recorded-capell-obstruction",
            "label": "PROVED_PER_ROW",
            "cell": row["cell"],
            "family": row["family"],
            "a": frac_text(a),
            "witness_prime": prime,
            "affine_parameter_residue": residue,
            "theta_residue": b_mod,
            "P_value_mod_prime": value,
            "simple_derivative_mod_prime": derivative,
            "two_theta_squareclass_argument": bad_arg,
            "legendre": -1,
            "deduction": "2theta is nonsquare in Q[b]/(P)",
        })
        if index % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
    summary = {
        "type": "recorded-capell-obstruction-summary",
        "label": "PROVED_PER_ROW",
        "source": str(STEP_II.relative_to(HERE)),
        "rows": len(proof_rows),
        "families": {
            family: sum(row["family"] == family for row in proof_rows)
            for family in ("L11", "ESC")
        },
        "witness_prime_min": min(primes),
        "witness_prime_max": max(primes),
        "meaning": (
            "each simple degree-one residue place has 2theta nonsquare; these are "
            "exact rowwise obstructions, not evidence for unrecorded a-values"
        ),
    }
    assert summary["families"] == {"L11": 190, "ESC": 103}
    return summary, proof_rows


def source_record() -> dict[str, Any]:
    return {
        "type": "source-map",
        "label": "PROVED",
        "sources": [
            {
                "citation": (
                    "Y. Harpaz, D. Wei, O. Wittenberg, Rational points on "
                    "fibrations with few non-split fibres"
                ),
                "url": "https://arxiv.org/abs/2109.03547",
                "used": "Theorems 1.3 and 1.4; definition of rank",
            },
            {
                "citation": (
                    "Y. Harpaz, A. Skorobogatov, O. Wittenberg, The "
                    "Hardy-Littlewood conjecture and rational points"
                ),
                "url": "https://arxiv.org/abs/1304.3333",
                "used": "unconditional theorem when geometric bad fibres are Q-rational",
            },
            {
                "citation": (
                    "T. Browning, D. Schindler, Strong approximation and a "
                    "conjecture of Harpaz and Wittenberg"
                ),
                "url": "https://arxiv.org/abs/1509.07744",
                "used": "Theorem 1.1 rank <=3 case",
            },
            {
                "citation": "A. Shute, Polynomials represented by norm forms via the beta sieve",
                "url": "https://arxiv.org/abs/2209.08949",
                "used": "Theorem 1.4 degree <=2/3 cases and Chatelet comparison",
            },
            {
                "citation": "C. Frei, E. Sofos, Random conic bundle surfaces satisfy the Hasse principle",
                "url": "https://arxiv.org/abs/2604.07047",
                "used": "Theorem 1.1 is statistical, not pointwise",
            },
            {
                "citation": (
                    "J.-L. Colliot-Thelene, A. Skorobogatov, "
                    "The Brauer-Grothendieck Group"
                ),
                "url": "https://link.springer.com/book/10.1007/978-3-030-74248-5",
                "used": "Faddeev/purity framework and conic-bundle Brauer generators",
            },
        ],
        "withdrawn_sources_used": False,
    }


def write_report(records: list[dict[str, Any]]) -> None:
    replay = next(
        row for row in records
        if row["type"] == "recorded-capell-obstruction-summary"
    )
    target = next(
        row for row in records if row["type"] == "target-prime-capell-no-go"
    )
    lines = [
        "# L23: fixed-cell conic bundle and unconditional fibration theorems",
        "",
        "## Status",
        "",
        "- **PROVED:** on `tau_dagger`, the member problem is the rational-point problem on",
        "  `X_{a,Z}: U^2-2bV^2=P_{a,Z}(b)W^2`.",
        "- **PROVED:** its natural conic-bundle model has ten geometric degenerate fibres,",
        "  in closed degrees `1+8+1`; the splitting classes are `A, 2theta, A`.",
        "- **PROVED:** if `2theta` is nonsquare, the non-split rank is 10 and",
        "  `Br(X)/Br(Q)=Z/2`, generated by `(A,b)`.  If it is square, the rank is 2",
        "  and the Brauer quotient is zero.",
        "- **PROVED rank-10 scopes:** if `v_w(z)` is odd, the target Newton edge works",
        "  for every admissible `a`.  For every cell, one can instead make the free",
        "  refinement `a!=1 (mod w)` and then use the proved L22 HIT selection in that",
        "  progression; the simple target-prime root makes `2theta` nonsquare.",
        "- **OPEN residual scope:** for a pre-existing class with even `v_w(z)` and",
        "  `w|s`, this uniform target argument is undecided; the 293 stored classes",
        "  below are settled only by their exact per-row certificates.",
        "- **OPEN:** no checked unconditional Hasse/weak-approximation theorem decides",
        "  the proved fixed rank-10 scopes.  The first failed Harpaz--Wei--Wittenberg",
        "  hypothesis is exactly `rank<=2`.",
        "",
        "Strict labels are used below.  Finite scans and theorem refusals are never negative evidence.",
        "",
        "## 1. Birational model",
        "",
        "Put",
        "",
        "```text",
        "A=1+4a^2,  s=(a-1)/2,  D=1-Z-a^2Z^2,",
        "Ng=16a^4b^2-A(b-1)^4,",
        "P=16D^2A^2b^4+(4a^8/A)Z^4Ng^2-32A^3s^2D^2b^5.",
        "```",
        "",
        "On the canonical branch `tau_dagger=(A+1)/(2A)`,",
        "`delta=-4a^4/A` and `alpha=-delta*A=4a^4=(2a^2)^2`.  Square-class",
        "invariance therefore turns the tied quaternion into `(P(b),2b)`.  Its conic is",
        "",
        "```text",
        "U^2-2bV^2=P(b)W^2.                                      (1)",
        "```",
        "",
        "For fixed `Z` and variable `a`, (1) is the generic conic of a threefold over",
        "the `(a,b)`-plane.  Clearing the only displayed denominator gives",
        "",
        "```text",
        "A(U^2-2bV^2)=4H(a,Z,b)W^2,",
        "H=(A/4)P=X^2+A*L*Y^2,",
        "X=a^4Z^2Ng,  Y=2ADb^2,  L=1-2As^2b.",
        "```",
        "",
        "The script verifies this identity in `Q[a,Z,b]`, not by specialization.  Fixing",
        "`a` gives the surface to which one-variable fibration theorems apply.",
        "",
        "## 2. Degenerate fibres and splitting fields",
        "",
        "The diagonal determinant of (1) has square class `2bP(b)`.  Since the proved",
        "L22 premise says that `P` is irreducible of degree eight, and characteristic zero",
        "gives separability, the closed discriminant points are",
        "",
        "| point | degree | component/splitting square class |",
        "|---|---:|---|",
        "| `b=0` | 1 | `P(0) ~ A` |",
        "| `P(b)=0`, `theta=b mod P` | 8 | `2theta` in `K=Q(theta)` |",
        "| `b=infinity` | 1 | `A` |",
        "",
        "Indeed `P_0=P_8=4Aa^8Z^4=A(2a^4Z^2)^2`, and the infinity residue is",
        "`256/P_8 ~ A`.  Thus the natural regular model has ten geometric degenerate",
        "fibres and `K_X^2=8-10=-2`.  Also",
        "",
        "```text",
        "Norm_K/Q(2theta)=2^8 P_0/P_8=256=16^2.",
        "```",
        "",
        "This norm equality is necessary but not sufficient for `2theta` to be a square.",
        "",
        "## 3. Vertical Brauer group",
        "",
        "Faddeev residues at `(0,P,infinity)` take values in the spans of",
        "`(A,2theta,A)`.  Corestriction gives the sole relation `e_0=e_infinity`,",
        "because `A` is nonsquare (`A=5 mod 8` in `Q_2`) while `Norm(2theta)` is a",
        "square.  Quotienting by the generic conic class `beta=(P,2b)`, which vanishes",
        "on `X`, gives exactly",
        "",
        "```text",
        "2theta nonsquare: Br(X)/Br(Q) = Z/2, generator (A,b); rank=10.",
        "2theta square:    Br(X)/Br(Q) = 0;                    rank=2.",
        "```",
        "",
        "Here rank is the degree of the *non-split* locus.  It is not the number of",
        "geometric singular fibres on the natural model, which stays ten; in the",
        "hypothetical square case the eight split fibres have contractible components.",
        "The conic-bundle generator framework is the one in Colliot--Thélène--",
        "Skorobogatov, *The Brauer--Grothendieck Group*; the displayed three-entry",
        "residue calculation is carried out here in full.",
        "",
        "",
        "## 4. The exact Capell/Fable locus and its target-prime obstruction",
        "",
        "For irreducible `P`, the following are equivalent (**PROVED**):",
        "",
        "1. `2theta` is a square in `K`;",
        "2. `P(u^2/2)` factors into conjugate degree-eight factors over `Q`;",
        "3. `H=E(b)^2-2b O(b)^2` for `deg E<=4`, `deg O<=3`.",
        "",
        "The last line is an exact finite cover of the `(a,Z)`-locus: if",
        "`H=sum h_k b^k`, equate each coefficient with",
        "`sum_{i+j=k}e_i e_j-2 sum_{i+j=k-1}o_i o_j`.  It is not the historical",
        "L8c substitution of a *rational* Phi-value `b`; `theta` is the algebraic root.",
        "",
        "There are two **PROVED** target-prime no-go scopes.  First put",
        "`e=v_w(z)`.  If `e` is odd, no condition on `s` is needed: the left",
        "Newton edge of `P` joins `(0,12e)` to `(4,0)`, has slope `-3e`, and its",
        "endpoint quartic residual polynomial is separable because `w` is odd.",
        "Thus a place above `w` has `ord(theta)=3e`, which is odd; since `2` is a",
        "`w`-unit, `2theta` is nonsquare.  Hence the rank is 10 for every admissible",
        "`a` on every odd-`v_w(z)` cell.",
        "",
        "For arbitrary `e>=1`, refine L20's existing choice of `a` so that",
        "`(A|w)=-1` and `a!=1 (mod w)`; the good nonzero residues occur in `+/-`",
        "pairs, so this is free.  Apply the proved L22 HIT selection inside that",
        "refined progression.  Then `s` is a `w`-unit, `D=1 mod w`, and",
        "",
        "```text",
        "P(b)=16A^2b^4(1-2As^2b) (mod w).",
        "beta=(2As^2)^(-1) is a simple root,",
        "2beta=1/(As^2) is a nonsquare mod w.",
        "```",
        "",
        "Hensel gives a degree-one `w`-adic place of `K` at which `2theta` is a",
        "nonsquare.  This proves rank 10 on every cell after the free refinement plus",
        "L22 HIT.  It does **not** prove a uniform claim for a pre-existing class with",
        "even `v_w(z)` and `w|s`; that residual case remains **OPEN** unless certified",
        "row by row.  The script replayed such a simple-root certificate on all",
        f"{replay['rows']} stored aligned classes ({replay['families']['L11']} L11,",
        f"{replay['families']['ESC']} ESC; witness primes {replay['witness_prime_min']}..",
        f"{replay['witness_prime_max']}) as **PROVED per row**.",
        "",
        "The small modular check inside this script is only **EVIDENCE**:",
        f"{target['finite_replay']['point_evaluations']} evaluations through prime",
        f"{target['finite_replay']['maximum_prime']}.  It is not part of the proof.",
        "",
        "## 5. Exact theorem comparison",
        "",
        "### Harpaz--Wei--Wittenberg",
        "",
        "The aligned local certificates give points at the controlled places; at every",
        "other odd place one may choose a unit `b` with `P(b)` also a unit, and at infinity",
        "choose `b>0`.  Hence the surface has an adelic point.  In the hypothetical",
        "rank-two case the Brauer quotient is zero, so this adelic point lies in the",
        "full Brauer--Manin set.",
        "",
        "Theorem 1.3(i) of `arXiv:2109.03547` is the strongest exact match: rank at",
        "most two implies density of rational points in the full Brauer--Manin set.",
        "If `2theta` were square, the rank would be two and the Brauer quotient zero,",
        "so this would be an unconditional weak-approximation route.  The target-prime",
        "lemmas prove that its first hypothesis fails on every odd-valuation cell for",
        "every admissible `a`, and on every cell after the free refinement plus L22 HIT:",
        "the rank is 10.  The rank-three clauses fail numerically as well.  Theorem 1.4",
        "instead asks for homogeneous Schinzel `(HH_1)` at the irreducible octic point;",
        "that is conditional.",
        "",
        "### Chatelet and fixed norm forms",
        "",
        "A classical Chatelet surface has `N_{L/Q}(x)=f(t)` for one fixed quadratic",
        "field and `deg f=3` or `4`.  Here the quadratic algebra is",
        "`Q(b)(sqrt(2b))`, which varies with the base, and the bundle has ten bad fibres.",
        "The substitution `b=u^2/2` is a degree-two base change splitting the generic",
        "conic, not a birational change over `Q`.  Thus Chatelet/CTSSD does not apply.",
        "",
        "### Green--Tao/Matthiesen/Browning and beta-sieve cases",
        "",
        "`arXiv:1304.3333` requires all geometric bad fibres to be defined over `Q`;",
        "the irreducible octic supplies eight non-rational ones.  Browning--Schindler",
        "Theorem 1.1 (`arXiv:1509.07744`) first fails at `rank<=3`, since this rank is 10.",
        "Shute's Theorem 1.4 (`arXiv:2209.08949`) treats residue fields of degree at",
        "most two, or at most three with extra hypotheses; its first failure is",
        "`[K:Q]=8`.  Its direct norm-form theorems also use a fixed number field.",
        "Frei--Sofos (`arXiv:2604.07047`) is a density-one theorem as coefficients vary",
        "and gives no result for this fixed surface.",
        "",
        "## 6. Constant twists do not repair the exact section",
        "",
        "A fixed twist `j` with `2theta/j` square only changes a simple value-prime",
        "sign to `(j|R)`; it does not force `+1`.  If one also asks for a polynomial norm",
        "identity `P=c(E^2-jO^2)`, then at `theta` the ratio `E(theta)/O(theta)` puts",
        "`sqrt(j)` in `K`.  Therefore `j` is already a square in `K`, and the twisted",
        "criterion implies the forbidden untwisted criterion.  Without the norm identity,",
        "forcing `(j|R)=+1` at every odd value-prime is precisely the dimension-`1/2`",
        "Chebotarev sieve problem, not Dirichlet for one linear prime.",
        "",
        "## 7. The Phi restriction",
        "",
        "For this base, `Phi_1^{2}` is exactly",
        "",
        "```text",
        "a in 1+2Z_2,   b in Z_2^times.",
        "```",
        "",
        "After fixing odd `a`, the condition on `b` is an open local condition and weak",
        "approximation can impose it.  It is not Zariski open.  Requiring the proof-side",
        "shape `b=eps*f*Q` with `Q` prime in an arithmetic progression is not local-open,",
        "but it is also not part of the original Phi predicate; a genuine surface theorem",
        "would bypass that auxiliary shape.",
        "",
        "## Verdict",
        "",
        "The birational comparison isolates a precise obstruction rather than an analogy.",
        "On odd-valuation cells for every admissible `a`, and on every cell after the",
        "free `a!=1 mod w` refinement plus L22 HIT, the fixed surface is a rank-10",
        "conic bundle with one degree-eight non-split closed point and a nontrivial",
        "vertical Brauer class.  The only known unconditional few-fibre theorem that",
        "would close the member step requires that degree-eight fibre to split.  For",
        "pre-existing even-valuation classes with `w|s`, only the 293 exact rowwise",
        "certificates are asserted.  Removing Schinzel still requires the half-dimensional",
        "locally-split-value sieve (or a genuinely new fixed rank-10 theorem).  **OPEN.**",
        "",
        "## Sources",
        "",
        "- Y. Harpaz, D. Wei, O. Wittenberg, *Rational points on fibrations with few",
        "  non-split fibres*, arXiv:2109.03547, Theorems 1.3--1.4.",
        "- Y. Harpaz, A. Skorobogatov, O. Wittenberg, *The Hardy--Littlewood",
        "  conjecture and rational points*, arXiv:1304.3333.",
        "- T. Browning, D. Schindler, *Strong approximation and a conjecture of",
        "  Harpaz and Wittenberg*, arXiv:1509.07744, Theorem 1.1.",
        "- A. Shute, *Polynomials represented by norm forms via the beta sieve*,",
        "  arXiv:2209.08949, Theorem 1.4.",
        "- C. Frei, E. Sofos, *Random conic bundle surfaces satisfy the Hasse principle*,",
        "  arXiv:2604.07047, Theorem 1.1 (statistical scope only).",
        "- J.-L. Colliot--Thélène, A. Skorobogatov, *The Brauer--Grothendieck",
        "  Group*, Springer, 2021 (Faddeev/purity and conic-bundle Brauer groups).",
        "",
        "Machine verdict: **OPEN**.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    started = time.monotonic()
    records: list[dict[str, Any]] = [
        source_record(),
        symbolic_identity_record(),
        geometry_record(),
        brauer_record(),
        phi_record(),
        capell_record(),
        target_prime_replay(),
        twist_record(),
        theorem_match_record(),
        *comparison_records(),
    ]
    replay_summary, replay_rows = replay_step_ii()
    records.append(replay_summary)
    records.extend(replay_rows)
    elapsed = time.monotonic() - started
    records.append({
        "type": "summary",
        "label": "OPEN",
        "birational_model": "PROVED",
        "degenerate_fibres": "PROVED:10 with closed degrees 1,8,1",
        "vertical_brauer": (
            "PROVED:Z/2 on odd-v_w(z) cells for every admissible a, and on all "
            "cells after a!=1 mod w plus L22 HIT"
        ),
        "capell_rank2_route": (
            "PROVED impossible on those scopes; OPEN uniformly for pre-existing "
            "even-valuation w|s classes absent a row certificate"
        ),
        "recorded_capell_obstructions": f"PROVED_PER_ROW:{replay_summary['rows']}",
        "unconditional_member_existence": "OPEN",
        "first_failed_named_hypothesis": (
            "Harpaz-Wei-Wittenberg rank<=2; rank=10 on the proved scopes"
        ),
        "records": len(records) + 1,
        "wall_seconds": round(elapsed, 6),
    })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = "\n".join(
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records
    ) + "\n"
    OUT.write_text(payload, encoding="utf-8")
    write_report(records)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    print(
        "L23 fibration: model PROVED; fibres 1+8+1; rank=10 and Br/Br0=Z/2 "
        "on odd-v_w(z) cells and after refined L20+L22 HIT; residual scope OPEN; "
        "unconditional verdict OPEN"
    )
    print(f"replayed {replay_summary['rows']} exact recorded Capell obstructions")
    print(f"wrote {OUT} sha256={digest}")
    print(f"wrote {REPORT}; wall={elapsed:.3f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
