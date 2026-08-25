#!/usr/bin/env python3
"""L31 exact reciprocal member and sharpened remaining frontiers.

This producer has four deliberately different outcomes.

* It promotes one L29 fixed-field reciprocal row to a proved global member by
  an explicit factor/primality certificate and a complete Hilbert-symbol table.
* It isolates an exact rational point on the L30 constant-two quartic before
  the required cube pullback, and proves precisely why that point is not an
  admissible target-cell point.
* It rewrites the L23 two-large-bad-divisor problem as a root-class dispersion
  error beyond the Bombieri--Vinogradov range.  This is a proved reduction, not
  the missing estimate.
* It identifies AP1 (R_bad <= 1 in one selected aligned class per cell) as the
  exact parity-based replacement for classical Schinzel H.  AP1 remains open.

No finite search is promoted to a uniform theorem.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

from h10q import OO, factorint, hilbert, legendre, ramified, vp
from l26_reciprocal_tie import (
    bridge_data,
    evaluate,
    reciprocal_polynomials,
    reciprocal_tie,
)
from l30_quartic_frontier import constant_two_data

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l31_frontier_push.jsonl"
REPORT = Path("/tmp/l31_frontier_push.md")


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def unit_mod(value: F, modulus: int) -> int:
    value = F(value)
    assert math.gcd(value.denominator, modulus) == 1
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def rational_cube(value: F) -> F | None:
    value = F(value)
    sign = -1 if value < 0 else 1
    value = abs(value)

    def integer_cube(number: int) -> int | None:
        if number == 0:
            return 0
        root = 1 << ((number.bit_length() + 2) // 3)
        while True:
            candidate = (2 * root + number // (root * root)) // 3
            if candidate >= root:
                break
            root = candidate
        for candidate in range(max(0, root - 2), root + 3):
            if candidate**3 == number:
                return candidate
        return None

    numerator = integer_cube(value.numerator)
    denominator = integer_cube(value.denominator)
    if numerator is None or denominator is None:
        return None
    return F(sign * numerator, denominator)


def h2_value(a: F, b: F, c: F, m: F) -> F:
    A = 1 + 4 * a * a
    B = 2 * b
    s = (a - 1) / 2
    q = (
        (s + 1) ** 2 * m * m
        - 2 * A * (s * s + 1) * m
        + (s - 1) ** 2 * A * A
    )
    return -(c * c - 4 * A) * (m - A) ** 2 - 16 * A * B * q - 64 * A * m


def quartic_canonical_sections() -> dict:
    """The four canonical hyperbola parameters are not L30 roots on Phi."""
    checks = 0
    for a_residue in range(1, 64, 2):
        a = F(a_residue)
        A = 1 + 4 * a * a
        s = (a - 1) / 2
        T = s * (1 + 2 * a * a) - a * a * (s * s + 1)
        assert T.denominator == 1 and T.numerator % 2
        for b_residue in range(1, 64, 2):
            b = F(b_residue)
            for c_quotient in range(8):
                c = F(16 * c_quotient)
                for m in (F(1), A * A):
                    C2 = h2_value(a, b, c, m) / (4 * m)
                    assert vp(C2, 2) == 7
                    assert unit_mod(C2, 256) == 128
                    checks += 1
    return {
        "type": "constant-two-canonical-section-no-go",
        "label": "PROVED uniformly on Phi",
        "sections": ["lambda=1", "lambda=-1", "lambda=A", "lambda=-A"],
        "congruence": "C2 == 128 (mod 256)",
        "proof_m_1": (
            "H2(1)/A=64*(a^4-1)+512*b*a^2*T-16*a^4*c^2/A, "
            "T=s*(1+2*a^2)-a^2*(s^2+1) odd"
        ),
        "proof_m_A2": (
            "Q_s(A^2)/A^2 == 16 (mod 32); the same division by 4*m "
            "leaves the unit class 128 (mod 256)"
        ),
        "exact_residue_checks": checks,
        "scope": "excludes only the four canonical hyperbola parameters, not the full quartic",
    }


def constant_c_bridge_section() -> list[dict]:
    """A rational H2 root just outside the mandatory cube pullback."""
    a = F(1)
    A = F(5)
    c = F(-64, 25)
    lam = F(3)
    m = lam * lam
    b = F(2101, 25000) - 2 * m / (m - 5) ** 2
    assert b == F(-3253, 3125)

    u = (b - 1) ** 2 / b
    bridge_q = -F(5, 4) * u
    Z = 2 / (bridge_q + 1)
    assert bridge_q == F(10169721, 2033125)
    assert Z == F(2033125, 6101423)

    A0, D, Ng, bridge_c = bridge_data(a, Z, b)
    assert A0 == A and bridge_c == c
    A1, c1, X, rho, C1, C2, m1, H = constant_two_data(a, b, Z, lam)
    assert (A1, c1, m1) == (A, c, m)
    assert C1 == C2 == H == 0
    assert X == F(7, 3) and rho == F(2, 3)
    assert vp(a, 2) == vp(b, 2) == 0

    numerator_factors = factorint(Z.numerator)
    denominator_factors = factorint(Z.denominator)
    assert numerator_factors == {5: 4, 3253: 1}
    assert denominator_factors == {1009: 1, 6047: 1}
    assert rational_cube(Z) is None
    assert vp(b, 3253) == vp(Z, 3253) == 1
    assert vp(c, 3253) == 0

    identity_checks = 0
    for test_b in (F(-7), F(-3), F(-1, 3), F(3), F(7, 5)):
        test_u = (test_b - 1) ** 2 / test_b
        Ng_over_b2 = 16 - 5 * test_u * test_u
        discriminant = 5 + 4 * Ng_over_b2 / (5 * c)
        assert discriminant == (F(5, 4) * test_u) ** 2
        identity_checks += 1

    cube_checks = 0
    for test_z in map(F, (-3, -1, F(-1, 2), F(1, 2), 1, 2, 3)):
        test_t = test_z**3
        plus = -64 * test_t * test_t + 96 * test_t + 64
        minus = 96 * test_t * test_t - 224 * test_t + 64
        assert plus == -32 * (2 * test_t + 1) * (test_t - 2)
        assert minus == 32 * (3 * test_t - 1) * (test_t - 2)
        cube_checks += 1

    section = {
        "type": "constant-c-bridge-section",
        "label": "PROVED exact rational section before the cube condition",
        "specialization": {"a": "1", "c": "-64/25"},
        "bridge_discriminant": (
            "q^2=5+4*(16-5*u^2)/(5*c)=(5*u/4)^2, "
            "u=(b-1)^2/b"
        ),
        "uniqueness": (
            "for finite constant c, 5+64/(5*c)-4*u^2/c is a polynomial square "
            "in Q[u] iff c=-64/25"
        ),
        "H2_root_family": "b=2101/25000-2*lambda^2/(lambda^2-5)^2",
        "bridge_maps": [
            "Z_plus=8*b/(5*b^2-6*b+5)",
            "Z_minus=8*b/(-5*b^2+14*b-5)",
        ],
        "exact_identity_checks": identity_checks,
    }
    point = {
        "type": "constant-two-off-cube-rational-point",
        "label": "PROVED exact near miss; NOT an admissible L30 target point",
        "parameters": {
            "a": frac_text(a),
            "b": frac_text(b),
            "Z": frac_text(Z),
            "lambda": frac_text(lam),
            "c": frac_text(c),
            "X": frac_text(X),
            "rho": frac_text(rho),
        },
        "guards": {"D_nonzero": D != 0, "Ng_nonzero": Ng != 0, "Phi_at_2": True},
        "Z_factorization": {
            "numerator": {str(p): e for p, e in numerator_factors.items()},
            "denominator": {str(p): e for p, e in denominator_factors.items()},
        },
        "why_not_L30": [
            "Z is not a rational cube",
            "at w=3253, v_w(b)=v_w(Z)=1 but v_w(c)=0",
            "the standard cube-pullback target stratum requires v_w(Z)=3*v_w(z) and v_w(c)>=4",
        ],
    }
    cube_reduction = {
        "type": "constant-c-cube-pullback-reduction",
        "label": "PROVED reduction; EXTERNAL-CAS COMPLETE on this section",
        "plus_curve": "d^2=-64*z^6+96*z^3+64",
        "minus_curve": "d^2=96*z^6-224*z^3+64",
        "factorizations": [
            "-32*(2*z^3+1)*(z^3-2)",
            "32*(3*z^3-1)*(z^3-2)",
        ],
        "genus": 2,
        "squarefree_reason": (
            "the t=z^3 roots are {-1/2,2} and {1/3,2}; all are distinct "
            "and nonzero, so both degree-6 polynomials are squarefree"
        ),
        "exact_factorization_checks": cube_checks,
        "external_audit": {
            "artifact": "data/l31_magma_genus2.json",
            "engine": "Magma V2.29-9 RationalPointsGenus2",
            "complete_rational_points_both_curves": ["(0:-8:1)", "(0:8:1)"],
            "stdlib_replayed": False,
        },
        "derivation": (
            "substitute Z=z^3 into the two bridge maps and require the quadratic "
            "in b to have square discriminant"
        ),
        "strict_scope": (
            "the complete external CAS calculation closes only the unique constant-c "
            "section; it is not a no-point theorem for the full L30 quartic"
        ),
    }
    return [section, point, cube_reduction]


def pocklington_step(n: int, certified_factor: int, base: int) -> dict:
    assert (n - 1) % certified_factor == 0
    assert certified_factor * certified_factor > n
    fermat = pow(base, n - 1, n)
    residue = pow(base, (n - 1) // certified_factor, n)
    gcd_value = math.gcd(residue - 1, n)
    assert fermat == 1 and gcd_value == 1
    return {
        "n": str(n),
        "certified_factor": str(certified_factor),
        "cofactor": str((n - 1) // certified_factor),
        "base": base,
        "F_squared_gt_n": True,
        "fermat_residue": fermat,
        "factor_residue": str(residue),
        "gcd_residue_minus_one": gcd_value,
    }


def reciprocal_member() -> dict:
    """Promote the w=13 L29 row with a complete exact certificate."""
    w = 13
    a, z, Z, q0 = F(3), F(13), F(169), F(1)
    b = F(w) * q0 * q0
    rho = q0 + 1 / (w * q0)
    lam = w * q0 - 1 / q0
    u = b + 1 / b
    assert (b, rho, lam, u) == (F(13), F(14, 13), F(12), F(170, 13))
    assert lam * lam == w * (w * rho * rho - 4)

    A, D, Ng, c = bridge_data(a, Z, b)
    eta = reciprocal_tie(Z, D, b)
    delta = -4 * a**4 / A
    M = 16 - delta * c * c - 32 * A * b * eta * eta
    assert A == 37
    assert D == -257217 == -3 * 83 * 1033
    assert Ng == -548208 == -(2**4) * (3**6) * 47
    assert c == F(277941456, 3172343)
    assert eta == F(-182, 257217)
    assert delta == F(-324, 37)

    Q = 14082426920623718389
    denominator = 3**2 * 37**3 * 83**2 * 1033**2
    assert M == F(2**4 * Q, denominator)
    assert M == F(225318830729979494224, 3351232116513117)

    r1 = 1173535576718643199
    r2 = 6150025556911
    r3 = 46202581
    assert Q - 1 == 2**2 * 3 * r1
    assert r1 - 1 == 2 * 3**2 * 10601 * r2
    assert r2 - 1 == 2 * 3**3 * 5 * 17 * 29 * r3
    assert r3 - 1 == 2**2 * 3**2 * 5 * 283 * 907
    assert factorint(r3) == {r3: 1}
    primality_chain = [
        pocklington_step(Q, r1, 2),
        pocklington_step(r1, r2, 3),
        pocklington_step(r2, r3, 3),
    ]
    assert factorint(Q) == {Q: 1}

    expected_places = (2, 3, 13, 37, 83, 1033, Q, OO)
    symbol_rows = []
    for place in expected_places:
        symbol = hilbert(M, 2 * w, place)
        assert symbol == 1
        symbol_rows.append({"place": "infinity" if place == OO else str(place), "symbol": symbol})
    assert ramified(M, 2 * w) == []
    assert unit_mod(M / 16, 8) == 1
    assert unit_mod(M, 13) == 3 and pow(4, 2, 13) == 3
    assert legendre(F(26), 37) == 1
    assert Q % 8 == 5 and Q % 13 == 5
    assert legendre(F(2), Q) == legendre(F(13), Q) == -1

    P, T = reciprocal_polynomials(a, Z, delta)
    P_value = evaluate(P, b)
    T_value = evaluate(T, u)
    assert P_value == (D * A * b * b) ** 2 * M
    assert T_value == (D * A) ** 2 * M == F(2**4 * Q, 37)

    return {
        "type": "fixed-field-reciprocal-global-member",
        "label": "PROVED exact global member",
        "target": w,
        "parameters": {
            "a": frac_text(a),
            "z": frac_text(z),
            "Z": frac_text(Z),
            "q0": frac_text(q0),
            "b": frac_text(b),
            "rho": frac_text(rho),
            "lambda": frac_text(lam),
            "u": frac_text(u),
        },
        "bridge": {
            "A": frac_text(A),
            "D": frac_text(D),
            "Ng": frac_text(Ng),
            "c": frac_text(c),
            "eta": frac_text(eta),
            "delta": frac_text(delta),
        },
        "norm_value_M": frac_text(M),
        "M_factorization": "2^4*Q/(3^2*37^3*83^2*1033^2)",
        "Q": str(Q),
        "Q_primality": {"method": "recursive Pocklington", "chain": primality_chain},
        "hilbert_symbols": symbol_rows,
        "ramified_places": [],
        "reciprocal_lift": "lambda^2=w*(w*rho^2-4)=144",
        "trace_value": frac_text(T_value),
        "conclusion": (
            "(M,26) splits at every place; Hasse-Minkowski gives a rational point "
            "on U^2-26*V^2=M*W^2, and W is nonzero because 26 is nonsquare"
        ),
        "strict_scope": "one target w=13; not a per-target or uniform-in-w member theorem",
    }


def reciprocal_slice_scope() -> dict:
    checks = 0
    for w in (3, 5, 7, 11, 13, 17, 19):
        for numerator, denominator in ((1, 1), (3, 1), (1, 3), (5, 3)):
            q0 = F(numerator, denominator)
            b = w * q0 * q0
            rho = q0 + 1 / (w * q0)
            lam = w * q0 - 1 / q0
            assert b + 1 / b == w * rho * rho - 2
            assert lam * lam == w * (w * rho * rho - 4)
            checks += 1
    return {
        "type": "fixed-field-reciprocal-scope",
        "label": "PROVED reduction; uniform member theorem OPEN",
        "mandatory_slice": "b=w*q0^2",
        "lift": "rho=q0+1/(w*q0), lambda=w*q0-1/q0",
        "homogeneous_form": "G_w(r,s)=s^16*P_rec(w*r^2/s^2)",
        "degree": 16,
        "member_condition": "(G_w(r,s),2*w)_v=+1 at every place v",
        "q_support": (
            "away from fixed support, primes of q0 have even M-valuation "
            "(-8*v_p(q0) at zeros and the reciprocal even analogue at poles)"
        ),
        "A_support": (
            "if p|A has odd multiplicity under the standard unit guards, the fixed "
            "local condition is (2*w|p)=+1"
        ),
        "identity_checks": checks,
        "consequence": (
            "the exact w=13 member does not prove a density theorem; the general slice "
            "is a sign-decorated degree-16 binary-form norm problem"
        ),
    }


def crt_pair(left: int, p: int, right: int, q: int) -> int:
    return (left + p * (((right - left) * pow(p, -1, q)) % q)) % (p * q)


def dispersion_reduction() -> dict:
    p, q = 5, 7
    bad_p = (1, 3)
    bad_q = (2, 4, 6)
    classes = {crt_pair(rp, p, rq, q) for rp in bad_p for rq in bad_q}
    assert len(classes) == len(bad_p) * len(bad_q)
    for residue in range(p * q):
        assert (residue in classes) == (residue % p in bad_p and residue % q in bad_q)

    exact_odd_classes = []
    root = 2
    for residue in range(p * p):
        if residue % p == root and (residue - root) % (p * p) != 0:
            exact_odd_classes.append(residue)
    assert len(exact_odd_classes) == p - 1

    root_signs_p = (1, -1, -1, 1)
    root_signs_q = (-1, 1, -1)
    direct_bad_pairs = sum(
        left == right == -1 for left in root_signs_p for right in root_signs_q
    )
    Rp, Rq = len(root_signs_p), len(root_signs_q)
    Cp, Cq = sum(root_signs_p), sum(root_signs_q)
    character_expansion = (Rp * Rq - Rp * Cq - Cp * Rq + Cp * Cq) // 4
    assert direct_bad_pairs == character_expansion

    return {
        "type": "two-large-bad-divisor-dispersion-reduction",
        "label": "PROVED exact reduction; required estimate OPEN",
        "root_classes": "C_p^-={r mod p:G(r)=0 and chi_p(u(r))=-1}",
        "pair_classes": "C_pq^-=CRT(C_p^- x C_q^-)",
        "Lambda_sum": (
            "S_2^root=sum_{p<q, p,q>=z, D<pq<=X^(8+o(1))} "
            "sum_{c in C_pq^-} sum_{X<t<=2X, t=c mod pq} Lambda(Q(t))"
        ),
        "AP_decomposition": (
            "after dyadic truncation in the finite divisor range, "
            "A_{pq,c}=N*X/phi(N*p*q)+E_{pq,c}; the positive main is "
            "C_N*X*sum g_root(p)g_root(q)"
        ),
        "range": {
            "z": "X^(49/100)",
            "D_BV": "X^(1/2)/(log X)^B",
            "first_pair_product": "p*q>=X^(98/100)>D_BV",
            "largest_relevant_modulus": "X^(8+o(1))",
        },
        "prime_count_odd_sector": (
            "B_2^odd(X)=#{t in the L23a small-prime-clean set: there exist "
            "distinct p,q>=z with odd G-valuations and bad signs}"
        ),
        "overstrong_sufficient_bound": (
            "B_2^odd(X)=o(X/(log X)^(3/2)) would suffice, but is not proved "
            "and is not the natural expected asymptotic"
        ),
        "Lambda_weighted_analogue": (
            "the corresponding small-clean von Mangoldt bound is "
            "o(X/sqrt(log X)); this is distinct from the raw S_2^root sum above"
        ),
        "natural_missing_theorem": (
            "a fixed-family Buchstab/Hilbert-detector asymptotic with a positive "
            "zero-bad main term C*X/(log X)^(3/2); the two-large sector is its "
            "first term beyond the available distribution level"
        ),
        "positive_main_warning": (
            "Chebotarev plus partial summation gives "
            "sum_{X^0.49<=p<=X^8} g_root(p)=(1/2)*log(8/0.49)+o(1); "
            "the raw pair main is order X in Lambda weight (X/log X in prime count) "
            "and cannot be discarded by signed-character cancellation or AP errors alone"
        ),
        "why_existing_tools_stop": [
            "Bombieri-Vinogradov controls squarefree p*q only through X^(1/2-o(1))",
            "exact v_p(G)=1 uses p^2 classes and leaves product level X^(1/4-o(1))",
            "for p*q>X, the elementary X/(p*q)+1 bound is dominated by +1",
            "character cancellation controls signed terms but leaves the positive R_p*R_q/4 term",
            "Frei-Sofos is L2 over coefficient boxes, not pointwise for this fixed octic/AP",
        ],
        "finite_CRT_checks": p * q,
        "finite_exact_odd_classes": len(exact_odd_classes),
        "character_identity_bad_pairs": direct_bad_pairs,
    }


def ap1_equivalence() -> dict:
    truth_rows = []
    for bad_count in range(17):
        if bad_count % 2:
            continue
        implication = (bad_count <= 1) == (bad_count == 0)
        assert implication
        truth_rows.append({"R_bad": bad_count, "R_bad_le_1_iff_zero": implication})

    Q0, g1, g2 = 5, 7, 23
    signs = (legendre(F(2 * Q0), g1), legendre(F(2 * Q0), g2))
    assert signs == (-1, -1) and signs[0] * signs[1] == 1

    return {
        "type": "AP1-Schinzel-replacement",
        "label": "PROVED equivalence inside the L19-L22 protocol; AP1 OPEN",
        "AP1": (
            "for every cell, one selected L20/L22 aligned class has a Q(t)-prime "
            "member with R_bad(t)<=1"
        ),
        "parity": (
            "alignment fixes every place in S union {Q,infinity} to +1; all unsupported "
            "even valuations are +1; Hilbert reciprocity makes R_bad even"
        ),
        "equivalence": "AP1 iff intermediate per-cell H within this protocol",
        "strict_weakness": (
            "R_bad<=1 is the weakest parity-based threshold; total P2 is insufficient "
            "because two bad odd-valuation primes have product sign +1"
        ),
        "L23_relation": (
            "L23a removes bad root primes below X^0.49 but permits R_bad in "
            "{0,2,...,16}; it does not imply AP1"
        ),
        "toy_factor_tuple": {
            "polynomials_at_t_0": [Q0, g1, g2],
            "individual_signs": list(signs),
            "product_sign": signs[0] * signs[1],
        },
        "parity_truth_rows": truth_rows,
        "chain_if_AP1_were_proved": "AP1 => H => conditional Theorem C becomes unconditional",
        "current_chain": "classical Schinzel H => H => Theorem C",
    }


def frontier_record() -> dict:
    return {
        "type": "L31-strict-frontier",
        "label": "MIXED: one exact construction; three requested uniform/global theorems OPEN",
        "outcomes": {
            "global_H2_cube_pullback_root": "OPEN; exact off-cube point proved and its unique constant-c section externally closed",
            "fixed_field_reciprocal_member": "PROVED for the exact target w=13; uniform-in-w theorem OPEN",
            "two_large_bad_divisor_estimate": "OPEN; positive main and exact beyond-BV detector target isolated",
            "remove_classical_Schinzel_H": "OPEN; AP1 is the exact parity substitute and is not implied by L23",
        },
        "unchanged": [
            "classical Schinzel H remains the sole conjectural input to the six-universal result",
            "no unconditional quantifier record changes",
            "H10(Q) remains open",
        ],
    }


def build_report(records: list[dict], elapsed: float) -> str:
    member = next(row for row in records if row["type"] == "fixed-field-reciprocal-global-member")
    off_cube = next(row for row in records if row["type"] == "constant-two-off-cube-rational-point")
    return "\n".join(
        [
            "# L31 frontier push",
            "",
            "## Exact construction",
            "",
            f"The fixed-field reciprocal row `w={member['target']}` is a PROVED global member.",
            f"Its norm value is `{member['norm_value_M']}` and every Hilbert symbol is +1.",
            "The mandatory reciprocal lift is checked directly; no trace-base shortcut is used.",
            "",
            "## Quartic near miss",
            "",
            f"The exact point `{off_cube['parameters']}` satisfies H2(lambda^2)=0 and Phi at 2,",
            "but Z is not a cube and its shared target prime leaves c a unit.  It is not an L30 member.",
            "The unique constant-c bridge section reduces cube compatibility to two genus-2 curves;",
            "an external complete Magma audit finds only z=0 and closes that section.",
            "",
            "## Analytic boundary",
            "",
            "The two-large-bad-divisor CRT expansion begins at modulus X^0.98, beyond BV.",
            "Its unsigned positive main survives signed cancellation; the required fixed-family",
            "Buchstab/Hilbert-detector asymptotic remains OPEN.",
            "AP1 (R_bad<=1 per cell) is equivalent to intermediate H by even Hilbert parity,",
            "but L23a allows 0,2,...,16 large bad primes and does not prove AP1.",
            "",
            "## Scope",
            "",
            "One exact reciprocal member is not a uniform member theorem.  No classical Schinzel-H",
            "removal, unconditional count change, or H10(Q) conclusion is claimed.",
            "",
            f"Wall time: {elapsed:.3f} s.",
        ]
    ) + "\n"


def main() -> int:
    started = time.perf_counter()
    records: list[dict] = [quartic_canonical_sections()]
    records.extend(constant_c_bridge_section())
    records.extend(
        [
            reciprocal_member(),
            reciprocal_slice_scope(),
            dispersion_reduction(),
            ap1_equivalence(),
            frontier_record(),
        ]
    )
    elapsed = time.perf_counter() - started
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in records))
    REPORT.write_text(build_report(records, elapsed))
    print(
        "l31_frontier_push: exact reciprocal member w=13; "
        "off-cube H2 point proved; dispersion/AP1 remain open; "
        f"rows={len(records)} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
