#!/usr/bin/env python3
"""Audit the half-dimensional bad-prime sieve and its large-prime barrier.

This is an exact replay of checked-in local data plus theorem bookkeeping. It
is not a numerical experiment purporting to prove an analytic distribution
statement. Analytic output records the hypotheses and conclusions of named
sieve theorems.

Replay from the workspace root:
    nice -n 19 python3 math/h10q/l23_half_sieve.py
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
STEP = DATA / "l17_stepii.jsonl"
BADROOTS = DATA / "l17_badroots_closures.jsonl"
L20 = DATA / "l20_admissible.jsonl"
OUT = DATA / "l23_half_sieve.jsonl"
REPORT = Path("/tmp/l23_half_sieve.md")
PACE_EVERY = 25
PACE_SECONDS = 0.005
PRIME_CERTIFICATE_LIMIT = 300

PolyQ = list[F]
PolyP = list[int]


def frac_text(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def qtrim(poly: PolyQ) -> PolyQ:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def qadd(left: PolyQ, right: PolyQ) -> PolyQ:
    answer = [F(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return qtrim(answer)


def qscale(poly: PolyQ, scalar: F | int) -> PolyQ:
    return qtrim([F(scalar) * value for value in poly])


def qmul(left: PolyQ, right: PolyQ) -> PolyQ:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return qtrim(answer)


def qpow(poly: PolyQ, exponent: int) -> PolyQ:
    assert exponent >= 0
    answer = [F(1)]
    base = qtrim(poly)
    power = exponent
    while power:
        if power & 1:
            answer = qmul(answer, base)
        power >>= 1
        if power:
            base = qmul(base, base)
    return answer


def qcompose_affine(poly: PolyQ, origin: F, step: F) -> PolyQ:
    answer = [F(0)]
    affine = [origin, step]
    power = [F(1)]
    for coefficient in poly:
        answer = qadd(answer, qscale(power, coefficient))
        power = qmul(power, affine)
    return qtrim(answer)


def primitive_integer_poly(poly: PolyQ) -> list[int]:
    denominator = 1
    for coefficient in poly:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [coefficient.numerator * (denominator // coefficient.denominator) for coefficient in poly]
    content = math.gcd(*[abs(value) for value in integers if value])
    integers = [value // content for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def canonical_P(a_int: int, z: F) -> tuple[PolyQ, dict[str, F]]:
    """Return P(b) on tau_dagger, in ascending powers of b."""
    a = F(a_int)
    A = 1 + 4 * a * a
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    delta = -4 * a**4 / A
    b = [F(0), F(1)]
    bm1 = [F(-1), F(1)]
    Ng = qadd(qscale(qpow(b, 2), 16 * a**4), qscale(qpow(bm1, 4), -A))
    first = qscale(qpow(b, 4), 16 * D * D * A * A)
    second = qscale(qpow(Ng, 2), -delta * a**4 * Z**4)
    third = qscale(qpow(b, 5), -32 * A**3 * s * s * D * D)
    P = qadd(qadd(first, second), third)
    assert len(P) == 9 and P[0] == P[8] != 0
    return P, {"A": A, "Z": Z, "D": D, "s": s, "delta": delta}


def capell_composition(P: PolyQ, twist: F = F(1)) -> PolyQ:
    """Return P(u^2/(2*twist)); u^2=2*twist*b at a root."""
    answer = [F(0)] * (2 * (len(P) - 1) + 1)
    for degree, coefficient in enumerate(P):
        answer[2 * degree] = coefficient / ((2 * twist) ** degree)
    return qtrim(answer)


def ptrim(poly: PolyP, prime: int) -> PolyP:
    answer = [value % prime for value in poly]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def padd(left: PolyP, right: PolyP, prime: int) -> PolyP:
    answer = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] = (answer[index] + value) % prime
    for index, value in enumerate(right):
        answer[index] = (answer[index] + value) % prime
    return ptrim(answer, prime)


def psub(left: PolyP, right: PolyP, prime: int) -> PolyP:
    return padd(left, [-value for value in right], prime)


def pmul(left: PolyP, right: PolyP, prime: int) -> PolyP:
    answer = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] = (answer[i + j] + left_value * right_value) % prime
    return ptrim(answer, prime)


def pdivmod(dividend: PolyP, divisor: PolyP, prime: int) -> tuple[PolyP, PolyP]:
    dividend = ptrim(dividend, prime)
    divisor = ptrim(divisor, prime)
    assert divisor != [0]
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    inverse = pow(divisor[-1], prime - 2, prime)
    while dividend != [0] and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        coefficient = dividend[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            dividend[index + shift] = (dividend[index + shift] - coefficient * value) % prime
        dividend = ptrim(dividend, prime)
    return ptrim(quotient, prime), dividend


def pmod(poly: PolyP, modulus: PolyP, prime: int) -> PolyP:
    return pdivmod(poly, modulus, prime)[1]


def pgcd(left: PolyP, right: PolyP, prime: int) -> PolyP:
    left, right = ptrim(left, prime), ptrim(right, prime)
    while right != [0]:
        left, right = right, pmod(left, right, prime)
    inverse = pow(left[-1], prime - 2, prime)
    return ptrim([(value * inverse) % prime for value in left], prime)


def ppowmod(base: PolyP, exponent: int, modulus: PolyP, prime: int) -> PolyP:
    answer = [1]
    base = pmod(base, modulus, prime)
    power = exponent
    while power:
        if power & 1:
            answer = pmod(pmul(answer, base, prime), modulus, prime)
        power >>= 1
        if power:
            base = pmod(pmul(base, base, prime), modulus, prime)
    return answer


def primes_upto(limit: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, limit + 1):
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
    return primes


def is_irreducible_mod(poly: list[int], prime: int) -> bool:
    f = ptrim(poly, prime)
    degree = len(f) - 1
    if degree != 16:
        return False
    inverse = pow(f[-1], prime - 2, prime)
    f = [(coefficient * inverse) % prime for coefficient in f]
    x = [0, 1]
    power = x
    assert degree == 16
    # Rabin: 16 has only the prime divisor 2.
    for exponent in range(1, degree + 1):
        power = ppowmod(power, prime, f, prime)
        if exponent == degree // 2 and pgcd(psub(power, x, prime), f, prime) != [1]:
            return False
    return pmod(psub(power, x, prime), f, prime) == [0]


def v2_int(value: int) -> int:
    assert value
    value = abs(value)
    answer = 0
    while value % 2 == 0:
        value //= 2
        answer += 1
    return answer


def v2(value: F) -> int:
    assert value
    return v2_int(value.numerator) - v2_int(value.denominator)


def vp(value: F, prime: int) -> int:
    assert value
    numerator, denominator = abs(value.numerator), value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def lower_newton_hull(poly: PolyQ, prime: int = 2) -> list[tuple[int, int]]:
    points = [(degree, vp(coefficient, prime)) for degree, coefficient in enumerate(poly) if coefficient]
    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if F(y1 - y0, x1 - x0) >= F(y2 - y1, x2 - x1):
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def is_prime_trial(value: int) -> bool:
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


def legendre(value: F | int, prime: int) -> int:
    value = F(value)
    residue = value.numerator * pow(value.denominator, prime - 2, prime) % prime
    if residue == 0:
        return 0
    symbol = pow(residue, (prime - 1) // 2, prime)
    return -1 if symbol == prime - 1 else symbol


def replay_local_masses() -> dict[str, Any]:
    step = load_jsonl(STEP)
    badroots = load_jsonl(BADROOTS)
    assert len(step) == 295 and step[0]["type"] == "meta" and step[-1]["type"] == "summary"
    assert len(badroots) == 295 and badroots[0]["type"] == "meta" and badroots[-1]["type"] == "summary"
    classes = step[1:-1]
    root_classes = badroots[1:-1]
    assert len(classes) == len(root_classes) == 293
    assert all(row["proof"]["F_degree"] == 8 for row in classes)
    assert all(row["proof"]["bad_unit_nonsquare_in_Qtheta"] for row in classes)
    assert all(row["proof"]["fixed_small_prime_bad_masses"] == {"3": 0, "5": 0, "7": 0} for row in classes)

    local_pairs = exact_pairs = excluded_pairs = simple_roots = step_roots = 0
    exact_sum_mismatches = bound_violations = 0
    max_ratio = F(0)
    max_ratio_at: dict[str, Any] | None = None
    for class_row in root_classes:
        for local in class_row["bad_root_lists"]:
            p = int(local["p"])
            known = F(0)
            residual = F(0)
            for root in local["roots"]:
                if root["kind"] == "simple":
                    simple_roots += 1
                    assert F(root["mass_num"], root["mass_den"]) == F(1, p + 1)
                else:
                    step_roots += 1
                if root.get("mass_num") is not None:
                    mass = F(root["mass_num"], root["mass_den"])
                else:
                    mass = F(root["resolved_mass_num"], root["resolved_mass_den"])
                known += mass
                residual += F(root.get("residual_mass_upper_bound", "0"))
            local_pairs += 1
            if local["excluded"]:
                excluded_pairs += 1
            else:
                exact_pairs += 1
                stored = F(local["local_mass_num"], local["local_mass_den"])
                if stored != known:
                    exact_sum_mismatches += 1
                assert stored == known
            upper = known + residual
            ratio = upper / F(8, p + 1)
            if ratio > max_ratio:
                max_ratio = ratio
                max_ratio_at = {"cell": class_row["cell"], "p": p, "upper": frac_text(upper)}
            if upper > F(8, p + 1):
                bound_violations += 1
    summary = step[-1]["proof_replay"]
    assert summary["positive_local_pairs_p_le_10000_checked_against_bound"] == 102439
    assert summary["bound_violations"] == 0
    assert excluded_pairs == summary["excluded_pairs_p_le_10000"] == 67
    assert badroots[-1]["n_excluded_pairs"] == 67
    assert badroots[-1]["reconciliation"]["missing"] == 0
    assert exact_sum_mismatches == bound_violations == 0
    return {
        "type": "local_mass_replay",
        "label": "PROVED",
        "scope": "exact checked-in p-adic rows; Chebotarev dimension is theorem-based, not inferred from this finite window",
        "classes": len(classes),
        "class_degree_8": sum(row["proof"]["F_degree"] == 8 for row in classes),
        "class_nonsquare_certificates_replayed_by_L18": sum(row["proof"]["bad_unit_nonsquare_in_Qtheta"] for row in classes),
        "fixed_small_masses": {"3": 0, "5": 0, "7": 0},
        "bad_root_local_pairs_stored": local_pairs,
        "exact_local_pairs": exact_pairs,
        "excluded_singular_pairs_with_certified_interval": excluded_pairs,
        "simple_root_entries": simple_roots,
        "recursive_step_root_entries": step_roots,
        "exact_mass_sum_mismatches": exact_sum_mismatches,
        "upper_bound_violations": bound_violations,
        "max_ratio_to_8_over_p_plus_1_in_stored_lists": frac_text(max_ratio),
        "max_ratio_location": max_ratio_at,
        "L18_independent_pairs_p_le_10000": summary["positive_local_pairs_p_le_10000_checked_against_bound"],
        "L18_independent_bound_violations": summary["bound_violations"],
        "observed_occurrence_reconciliation": badroots[-1]["reconciliation"],
        "theoretical_chebotarev_dimension": 0.5,
        "finite_window_product_exponent_evidence": step[-1]["finite_root_window_evidence"]["product_c_lower"],
    }


def cell_prime(cell: list[Any]) -> int:
    w = int(cell[0])
    assert is_prime_trial(w) and w % 2 == 1
    return w


def w_adic_check(P: PolyQ, params: dict[str, F], cell: list[Any], a: int) -> dict[str, Any]:
    """Check the two rigorous local Capell obstructions when applicable."""
    w = cell_prime(cell)
    m = vp(params["Z"], w)
    A, s = params["A"], params["s"]
    assert legendre(A, w) == -1
    hull = lower_newton_hull(P, w)
    conclusion = False
    mechanism = "not-applicable"
    if m > 0 and m % 2:
        # The left edge is regular: its endpoint residual c0+c4*T^4 is
        # squarefree because w is odd and its roots are nonzero.
        assert hull[0] == (0, 4 * m) and hull[1] == (4, 0)
        conclusion = True
        mechanism = "odd-left-edge-valuation"
    elif m > 0 and m % 2 == 0 and s and vp(s, w) == 0:
        assert (4, 0) in hull and (5, 0) in hull
        beta = F(1, 1) / (2 * A * s * s)
        assert legendre(2 * beta, w) == -1
        conclusion = True
        mechanism = "horizontal-simple-root-nonsquare-residue"
    return {
        "w": w,
        "v_w_Z": m,
        "a_mod_w": a % w,
        "w_divides_s": s == 0 or vp(s, w) > 0,
        "mechanism": mechanism,
        "proves_2b_nonsquare": conclusion,
        "P_w_adic_lower_hull": [[x, y] for x, y in hull],
    }


def reciprocal_a_one_check(params: dict[str, F], a: int) -> dict[str, Any]:
    """Prove the reciprocal a=1 Capell obstruction by two trace norms."""
    if a != 1:
        return {"applicable": False, "proves_2b_nonsquare": False}
    Z, D = params["Z"], params["D"]
    p, q = Z.numerator, Z.denominator
    d = q * q - p * q - p * p
    assert D == F(d, q * q) and d % 2
    norm_minus_squareclass = 125 * d * d + 64 * p**4
    norm_plus_squareclass = 125 * d * d + 1024 * p**4
    assert norm_minus_squareclass % 8 == norm_plus_squareclass % 8 == 5
    assert math.isqrt(norm_minus_squareclass) ** 2 != norm_minus_squareclass
    assert math.isqrt(norm_plus_squareclass) ** 2 != norm_plus_squareclass
    return {
        "applicable": True,
        "proves_2b_nonsquare": True,
        "trace_coordinate": "u=b+b^-1",
        "square_implication": "2*b square forces one of 2*(u-2), 2*(u+2) square in Q(u)",
        "norm_squareclasses": {
            "2*(u-2)": str(norm_minus_squareclass),
            "2*(u+2)": str(norm_plus_squareclass),
        },
        "both_cleared_norms_mod_8": 5,
    }


def capell_grid_audit() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    l20 = load_jsonl(L20)
    rows = [row for row in l20 if row["type"] == "constructed-grid-class"]
    summary = next(row for row in l20 if row["type"] == "summary")
    assert len(rows) == summary["mechanical_crosscheck"]["classes"] == 353
    prime_list = [p for p in primes_upto(PRIME_CERTIFICATE_LIMIT) if p != 2]
    records: list[dict[str, Any]] = []
    certificate_counts: Counter[int] = Counter()
    hull_counts: Counter[str] = Counter()
    w_mechanisms: Counter[str] = Counter()
    failures: list[dict[str, Any]] = []
    affine_mismatches = norm_mismatches = 0
    for index, row in enumerate(rows):
        a = int(row["a"])
        z = F(row["z"])
        assert a % 2 == 1 and row["G_irreducible_over_Q"]
        P, params = canonical_P(a, z)
        origin = F(row["eps"]) * F(row["f"]) * F(row["q1"])
        step = F(row["eps"]) * F(row["f"]) * F(row["N"])
        primitive_composed = primitive_integer_poly(qcompose_affine(P, origin, step))
        expected_G = [int(value) for value in row["G_coefficients_ascending"]]
        if primitive_composed != expected_G:
            affine_mismatches += 1
        assert primitive_composed == expected_G
        norm = F(2**8) * P[0] / P[8]
        if norm != 256:
            norm_mismatches += 1
        assert norm == 256
        hull = lower_newton_hull(P)
        hull_text = ";".join(f"{x}:{y}" for x, y in hull)
        hull_counts[hull_text] += 1
        w_check = w_adic_check(P, params, row["cell"], a)
        reciprocal_check = reciprocal_a_one_check(params, a)
        w_mechanisms[w_check["mechanism"]] += 1
        certificate = None
        if not w_check["proves_2b_nonsquare"] and not reciprocal_check["proves_2b_nonsquare"]:
            C = primitive_integer_poly(capell_composition(P))
            certificate = next((p for p in prime_list if is_irreducible_mod(C, p)), None)
        capell_proved = (
            certificate is not None
            or w_check["proves_2b_nonsquare"]
            or reciprocal_check["proves_2b_nonsquare"]
        )
        record = {
            "type": "capell_grid_class",
            "label": "PROVED" if capell_proved else "OPEN",
            "scope": "this exact L20 class only; absence of a certificate/local obstruction is a refusal, never negative evidence",
            "cell": row["cell"],
            "a": str(a),
            "z": str(z),
            "Z": frac_text(params["Z"]),
            "v2_Z": v2(params["Z"]),
            "P_irreducible_authority": True,
            "P_u2_over_2_irreducibility_certificate_prime": certificate,
            "w_adic_obstruction": w_check,
            "reciprocal_a_one_obstruction": reciprocal_check,
            "capell_conclusion_2b_nonsquare_in_Qb_mod_P": capell_proved,
            "norm_2b": frac_text(norm),
            "norm_is_square_but_not_sufficient": True,
            "P_2adic_lower_hull": [[x, y] for x, y in hull],
        }
        records.append(record)
        if not capell_proved:
            failures.append({"cell": row["cell"], "a": a, "z": str(z), "hull": hull_text})
        if certificate is not None:
            certificate_counts[certificate] += 1
        if (index + 1) % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
    capell_summary = {
        "type": "capell_grid_summary",
        "label": "PROVED" if not failures else "OPEN",
        "scope": "exact 353-row L20 grid plus the uniform reciprocal a=1 trace-norm lemma",
        "classes": len(rows),
        "affine_reconstruction_mismatches": affine_mismatches,
        "norm_identity_mismatches": norm_mismatches,
        "composition_irreducibility_certificates": sum(row["P_u2_over_2_irreducibility_certificate_prime"] is not None for row in records),
        "w_adic_obstructions": sum(row["w_adic_obstruction"]["proves_2b_nonsquare"] for row in records),
        "reciprocal_a_one_obstructions": sum(row["reciprocal_a_one_obstruction"]["proves_2b_nonsquare"] for row in records),
        "total_rows_proved_nonsquare": len(rows) - len(failures),
        "certificate_refusals_after_local_checks": len(failures),
        "refusals_are_negative_evidence": False,
        "certificate_prime_limit": PRIME_CERTIFICATE_LIMIT,
        "certificate_prime_counts": {str(p): count for p, count in sorted(certificate_counts.items())},
        "w_adic_mechanism_counts": dict(sorted(w_mechanisms.items())),
        "two_adic_hull_counts": dict(sorted(hull_counts.items())),
        "uncertified_rows": failures,
        "capell_equivalence": "for irreducible P, 2*b is a square in Q[b]/(P) iff P(u^2/2) is reducible over Q",
        "norm_warning": "Norm(2*b)=2^8*P(0)/lc(P)=256 is a necessary square norm only, not a square certificate",
        "twist_warning": "2*j*b square only makes the root sign (2*b|p) equal (j|p); it does not make every sign +1 or give a fixed norm identity",
    }
    return records, capell_summary


def analytic_records() -> list[dict[str, Any]]:
    return [
        {
            "type": "prime_sequence_sieve",
            "label": "PROVED",
            "hypotheses": [
                "Q(t)=q1+N*t with N>0 and gcd(q1,N)=1",
                "G primitive irreducible degree 8 and Q*G fixed-divisor-free",
                "outside a finite set, G-roots are simple and Q is nonzero there",
                "u(t)=2*alpha*eps*f*Q(t); frozen, moving, and infinite Hilbert symbols are +1",
                "u(theta) nonsquare in K=Q(theta), G(theta)=0",
                "each finite exceptional prime has a CRT-compatible residue with Q(t)*G(t) nonzero, supplied by fixed-divisor-freeness",
            ],
            "definitions": {
                "bad_root_count": "r_-(p)=#{r mod p:G(r)=0 and Legendre(u(r),p)=-1}",
                "exact_Haar_odd_mass": "m_p=r_-(p)/(p+1) at good primes",
                "prime_sequence_exact_odd_density": "g_odd(p)=p*m_p/(p-1)=p*r_-(p)/(p^2-1)",
                "root_class_oversieve_density": "g_root(p)=r_-(p)/(p-1)=(p+1)*m_p/(p-1)",
                "oversieve_difference": "g_root(p)-g_odd(p)=r_-(p)/(p^2-1), a convergent correction",
            },
            "dimension": "kappa=1/2",
            "dimension_formula": "sum_{p<z} g_root(p)*log(p)=(1/2)*log(z)+O(1)",
            "source_of_half": "Chebotarev in K(sqrt(u(theta)))/Q plus nonsquareness; not a finite-data fit",
            "Euler_product": "prod_{p<z}(1-g_root(p))~C*(log z)^(-1/2), C>0",
            "finite_exception_handling": "choose one Q*G-nonzero residue at each exceptional prime and absorb the CRT class into the fixed modulus before applying Dirichlet/Bombieri-Vinogradov",
        },
        {
            "type": "beta_sieve_parameters",
            "label": "PROVED",
            "unconditional_distribution": {
                "theorem": "Bombieri-Vinogradov for primes in AP, fixed N",
                "level": "D=X^(1/2)/(log X)^B",
                "theta": 0.5,
                "weighted_root_classes": "bounded multiplicative root weights cost fixed divisor-function/log factors",
            },
            "semilinear_sieve": {
                "dimension_kappa": 0.5,
                "lower_sifting_limit_beta_kappa": 1,
                "parameter": "s=log(D)/log(z); positive lower main term requires s>1",
                "concrete_choice": "z=X^(49/100), s tends to 50/49>1",
                "general_choice": "z=X^(1/2-epsilon), fixed epsilon>0",
                "lower_bound": "#{X<t<=2X:Q(t) prime and no bad-sign root p<z} >> X/(log X)^(3/2)",
            },
            "exact_odd_cylinder_variant": {
                "warning": "distinguishing v_p(G)=1 from v_p(G)>=2 uses classes mod p^2",
                "plain_BV_base_level": "underlying product d<=X^(1/4-o(1)) when AP modulus is d^2",
                "best_simple_route": "use the mod-p root oversieve; its convergent excess preserves dimension",
            },
            "large_prime_arithmetic": {
                "degree": 8,
                "G_size": "|G(t)|=X^(8+o(1))",
                "at_z_X_49_over_100": "at most floor(8/(49/100))=16 factors >=z, with multiplicity",
                "reciprocity": "odd-valuation bad-sign count is even",
                "remaining_possibilities": [0, 2, 4, 6, 8, 10, 12, 14, 16],
                "cutoff_needed_for_parity_alone": "z>|G(t)|^(1/2)=X^(4+o(1)); conventional level D>X^(4+o(1))",
                "Elliott_Halberstam_still_insufficient": "theta=1 only permits z<X^(1-o(1)), leaving up to 8 factors",
            },
        },
        {
            "type": "almost_prime_threshold",
            "label": "PROVED",
            "parity_implication": "Omega(G(t))<=R forces zero bad primes from reciprocity only when R<=1",
            "exact_required_total_almost_prime_threshold": 1,
            "why_R_2_fails": "two distinct odd-valuation bad primes obey even parity but obstruct at both places",
            "current_general_theorems": [
                {"theorem": "Kao Theorem 1/Table 1", "degree": 8, "threshold_R": 12, "scope": "even granting fixed-AP adaptation, above R=1"},
                {"theorem": "Irving Theorem 1.1/Table 1", "degree": 8, "threshold_R": 14, "status": "superseded by Kao R=12"},
                {"theorem": "Chen P2 for a linear companion to a prime", "degree": 1, "threshold_R": 2, "status": "still above R=1"},
            ],
            "first_missing_estimate": "parity-sensitive saving for Q(t) prime and two large bad odd-valuation divisors p1,p2 of G(t), p1,p2>=X^(1/2-epsilon), p1*p2>D",
            "equivalent_missing_input": "an R_bad<=1 theorem for sign-decorated divisors of a degree-8 polynomial at prime arguments",
            "verdict": "standard half-dimensional, weighted, Chen, and known prime-argument almost-prime sieves do not remove Schinzel H",
        },
        {
            "type": "capell_uniform_local_lemma",
            "label": "PROVED",
            "hypotheses": "cell prime w odd; m=v_w(Z)>0; (A|w)=-1; P on tau_dagger",
            "odd_m": "left Newton edge (0,4m)->(4,0) is regular; e=1 and ord(2b)=m odd, so 2b is nonsquare",
            "even_m": "if w does not divide s, horizontal b4-b5 residual has simple beta=(2*A*s^2)^(-1), and (2*beta|w)=(A|w)=-1",
            "residue_choice": "#{a mod w:(1+4a^2|w)=-1}=(w+1)/2 for w=3 mod4 and (w-1)/2 for w=1 mod4; one can avoid a=1 (for w=5 use a=+/-2)",
            "reciprocal_a_one": "with P=b^4*T(u), T=400*D^2+(4/5)*Z^4*(16-5*(u-2)^2)^2, 2*b square forces one of 2*(u-2),2*(u+2) square; their norm squareclasses are 125*D^2+64*Z^4 and 125*D^2+1024*Z^4, whose cleared numerators are 5 mod 8",
            "scope_warning": "full compatibility of that refined residue choice with every global class condition is a separate class-construction statement",
            "nonsquare_twist": "2*j*b square implies (2*b|p)=(j|p) at simple value roots, not +1 uniformly; twisted Capell yields a variable norm U(b)^2-2*j*b*V(b)^2, not a fixed norm U^2-j*V^2",
        },
    ]


def build_report(local: dict[str, Any], capell: dict[str, Any]) -> str:
    local_w = capell["w_adic_obstructions"]
    reciprocal = capell["reciprocal_a_one_obstructions"]
    total = capell["total_rows_proved_nonsquare"]
    refused = capell["certificate_refusals_after_local_checks"]
    lines = [
        "# L23 — half-dimensional bad-prime sieve",
        "",
        "## Status",
        "",
        "**PROVED (barrier).** The dimension-`1/2` sieve is genuine and standard semilinear lower-bound technology gives many prime members free of every bad root prime below `X^(1/2-epsilon)`. It does **not** eliminate bad primes above that cutoff. Reciprocity says their number is even, but does not distinguish zero from two. The first missing input is a parity-sensitive estimate for two large bad divisors of a degree-eight value at a prime argument. No cited unconditional theorem supplies it; Schinzel H is not removed by this route.",
        "",
        "Finite local and Capell rows below are exact. Analytic distribution statements are invoked only with named hypotheses; finite scans are never analytic evidence.",
        "",
        "## 1. Sequence and exact local factor",
        "",
        "Fix an aligned class:",
        "",
        "```text",
        "Q(t)=q1+N*t,  b(t)=eps*f*Q(t),",
        "G(t)=primitive degree-8 part of P(b(t)),",
        "u(t)=2*alpha*eps*f*Q(t).",
        "```",
        "",
        "Take `X<t<=2X` with `Q(t)` prime and put fixed-data, discriminant, and resultant primes in a finite exceptional set. Fixed-divisor-freeness supplies, at each such prime, a residue on which `Q(t)G(t)` is nonzero. Choose these finitely many residues by CRT and absorb them into the fixed modulus. Thus exceptional primes are clean rather than silently discarded.",
        "",
        "At every remaining `p` define",
        "",
        "```text",
        "r_-(p)=#{r mod p:G(r)=0 and (u(r)|p)=-1}.",
        "```",
        "",
        "Roots are simple and `Q(r)!=0`. One simple root has Haar odd-valuation mass `1/(p+1)`. Relative to the prime sequence,",
        "",
        "```text",
        "m_p       = r_-(p)/(p+1),",
        "g_odd(p)  = p*m_p/(p-1) = p*r_-(p)/(p^2-1),",
        "g_root(p) = r_-(p)/(p-1) = (p+1)*m_p/(p-1).",
        "```",
        "",
        "`g_root` deliberately removes a bad-sign root even when the later valuation is even. This harmless oversieve keeps modulus `p`; its excess `r_-(p)/(p^2-1)` is summable. Exact valuation one versus at least two would require modulus `p^2`.",
        "",
        f"**Local replay (PROVED).** This script checked {local['exact_local_pairs']} exact local sums and {local['excluded_singular_pairs_with_certified_interval']} singular interval rows across {local['classes']} classes; every simple contribution is `1/(p+1)`, with {local['upper_bound_violations']} violations of `m_p<=8/(p+1)`. L18 independently records {local['L18_independent_pairs_p_le_10000']} positive pairs through `10^4`, also with zero violations. These are replays, not the source of Chebotarev.",
        "",
        "For a root `theta`, nonsquareness of `u(theta)` makes `K(sqrt(u(theta)))/K` nontrivial. Chebotarev gives",
        "",
        "```text",
        "sum_{p<z} g_root(p)log p=(1/2)log z+O(1),",
        "prod_{p<z}(1-g_root(p))~C/(log z)^(1/2).",
        "```",
        "",
        "Thus `kappa=1/2` exactly. A square `u(theta)` would collapse the tail; §5 treats the proposed `2b` square escape.",
        "",
        "## 2. Distribution and beta parameters",
        "",
        "CRT makes a squarefree product `d` of bad root conditions a bounded-multiplicity union of reduced classes for `Q(t)` modulo `Nd`. Bombieri–Vinogradov (fixed `N`, standard divisor/log weights) gives",
        "",
        "```text",
        "D=X^(1/2)/(log X)^B.",
        "```",
        "",
        "The dimension-`1/2` beta sieve is the **semilinear sieve**, with lower sifting limit `beta_(1/2)=1`. If `s=log D/log z`, positivity requires `s>1`. Hence",
        "",
        "```text",
        "z=X^(49/100), s->50/49>1,",
        "#{X<t<=2X:Q(t) prime, no bad root p<z} >> X/(log X)^(3/2).",
        "```",
        "",
        "Any fixed `z=X^(1/2-epsilon)` works. This is an unconditional **small-prime-clean** theorem, not an all-prime theorem. Directly encoding odd valuations uses AP moduli divisible by `p^2` and reduces the underlying product level to `X^(1/4-o(1))`; the mod-`p` oversieve avoids that loss.",
        "",
        "## 3. Large-prime barrier",
        "",
        "For `t~X`, `|G(t)|=X^(8+o(1))`. At `z=X^(49/100)` it can retain `floor(8/(49/100))=16` prime factors at least `z`. Reciprocity only narrows the bad odd-valuation count to `0,2,...,16`. Sifting to `z` is not eliminating all factors.",
        "",
        "Parity alone would force zero only for `z>|G(t)|^(1/2)=X^(4+o(1))`, which at semilinear sifting limit one demands conventional level `D>X^(4+o(1))`, versus BV's `X^(1/2-o(1))`. Even Elliott–Halberstam (`D=X^(1-o(1))`) leaves up to eight large factors.",
        "",
        "Every failing small-clean member has at least two bad primes `p1,p2>=z`, already with `p1*p2>D`. The first unavailable estimate is a saving bound for",
        "",
        "```text",
        "sum_{p1,p2>=z, p1*p2>D} #{t~X:Q(t) prime,",
        " p1*p2|G(t), both signs bad, both valuations odd}.",
        "```",
        "",
        "This parity-sensitive two-large-divisor/dispersion sector is outside BV and the beta-sieve fundamental lemma.",
        "",
        "## 4. Weighted and almost-prime route",
        "",
        "If `G(t) in P_R`, reciprocity forces zero bad factors only at the exact threshold `R<=1`. `R=2` fails: two distinct bad odd-valuation factors have product sign `+1` but obstruct at both places. Thus even Chen's `P_2` benchmark for a linear companion is insufficient.",
        "",
        "Kao's Theorem 1/Table 1 gives only `P_12` for an irreducible octic at prime arguments (improving Irving's `P_14`). Even granting the fixed-AP adaptation to `Q(t)`, 12 is far above 1. No cited Chen/Iwaniec-style theorem gives `R_bad<=1` here.",
        "",
        "## 5. Capell escape and twists",
        "",
        "For irreducible `P`, `K=Q[b]/(P)`, Capell gives",
        "",
        "```text",
        "2b square in K <=> P(u^2/2) reducible over Q.",
        "```",
        "",
        "The automatic identity `Norm(2b)=2^8P(0)/lc(P)=256` is necessary, not sufficient.",
        "",
        f"This script reconstructed all 353 L20 grid polynomials. The cell-prime Newton/residue obstruction applies to {local_w} rows. The reciprocal `a=1` trace-norm lemma applies to {reciprocal} rows (overlapping the local proof except on three even-valuation rows). Together they prove `2b` nonsquare on {total}/353 exact rows, with {refused} refusals.",
        "",
        "For completeness, on `a=1` write `P(b)=b^4T(v)` with `v=b+b^-1` and `T(v)=400D^2+(4/5)Z^4(16-5(v-2)^2)^2`. If `y^2=2b`, the involution `b<->b^-1` sends `y` to `+/-2/y`; accordingly one of `y+2/y` and `y-2/y` lies in `Q(v)`, so one of `2(v+2)` and `2(v-2)` is a square there. Their norms have squareclasses `125D^2+1024Z^4` and `125D^2+64Z^4`. For reduced `Z=p/q`, put `d=q^2-pq-p^2`; `d` is odd, and the cleared norms `125d^2+1024p^4` and `125d^2+64p^4` are both `5 mod 8`, hence nonsquares. This proves the reciprocal no-go for every rational nonzero `Z`.",
        "",
        "There is also a uniform local no-go after a compatible residue choice. Let `m=v_w(Z)>0` and `(A|w)=-1`. If `m` is odd, the regular left edge `(0,4m)->(4,0)` gives an unramified place with `ord(2b)=m` odd. If `m` is even and `w` does not divide `s`, the horizontal `b4-b5` residual has simple root `beta=(2As^2)^(-1)` and `(2beta|w)=(A|w)=-1`. The number of residues with `(1+4a^2|w)=-1` is `(w+1)/2` for `w=3 mod4` and `(w-1)/2` for `w=1 mod4`, so `a=1` can be avoided locally (at `w=5`, use `a=+/-2`). Compatibility with the full class construction is a separate structural statement.",
        "",
        "A nonsquare constant twist does not by itself rescue the argument. If `2jb` is square at the algebraic root, then at simple value primes `(2b|p)=(j|p)`, not uniformly `+1`. Twisted Capell yields a **variable** norm shape `U(b)^2-2jbV(b)^2`; it does not supply the separate fixed norm `U^2-jV^2` that would force inert `j`-primes to even valuation. Both conditions would be needed.",
        "",
        "A universal 2-adic shortcut is also unsafe: the nonreciprocal `b^5` coefficient can change the hull, and denominator-two residuals in characteristic two can be inseparable. No such claim is made.",
        "",
        "## 6. Verdict",
        "",
        "**The proved dimension-`1/2` sieve does not remove Schinzel H.** It supplies many members clean below `X^(1/2-epsilon)`. The first unavailable estimate is the two-large-bad-divisor parity sector, equivalently an `R_bad<=1` theorem for an octic at prime arguments. Current general almost-prime technology gives total `P_12` only.",
        "",
        "## Sources",
        "",
        "- H. Iwaniec, *The half dimensional sieve*, Acta Arith. 29 (1976), 69–95: https://eudml.org/doc/205410",
        "- J. Friedlander and H. Iwaniec, *Opera de Cribro*, AMS Colloquium Publications 57 (2010), Theorem 11.12 and p. 225 (semilinear sifting limit).",
        "- K. Nath and L. Xie, Lemma 3.6 and Remark 3.7 (`beta_(1/2)=1`), arXiv:2501.16723: https://arxiv.org/abs/2501.16723",
        "- P.-H. Kao, Theorem 1 and Table 1 (`r(8)=12`), arXiv:1606.03505: https://arxiv.org/abs/1606.03505",
        "- A. J. Irving, Theorem 1.1 and Table 1 (`r(8)=14`), arXiv:1410.3333: https://arxiv.org/abs/1410.3333",
        "- E. Bombieri, *On the large sieve*, Mathematika 12 (1965), 201–225; A. I. Vinogradov, *The density hypothesis for Dirichlet L-series*, Izv. Akad. Nauk SSSR 29 (1965), 903–934.",
        "- J. R. Chen, *On the representation of a large even integer as the sum of a prime and the product of at most two primes*, Sci. Sinica 16 (1973), 157–176.",
        "",
        f"Machine records: `{OUT}`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    started = time.monotonic()
    local = replay_local_masses()
    capell_rows, capell_summary = capell_grid_audit()
    records: list[dict[str, Any]] = [
        {
            "type": "meta",
            "artifact": "l23_half_sieve",
            "version": 1,
            "labels": {"exact_replays": "PROVED", "analytic_verdict": "PROVED", "capell_rows": "PROVED"},
            "inputs": {str(path.relative_to(HERE.parent.parent)): sha256(path) for path in (STEP, BADROOTS, L20)},
            "sources": [
                "Iwaniec 1976 half-dimensional sieve",
                "Bombieri-Vinogradov theorem",
                "Nath-Xie arXiv:2501.16723 Lemma 3.6/Remark 3.7",
                "Kao arXiv:1606.03505 Theorem 1/Table 1",
                "Irving arXiv:1410.3333 Theorem 1.1/Table 1",
                "Chen 1973 P2 theorem",
            ],
        },
        local,
        *analytic_records(),
        *capell_rows,
        capell_summary,
    ]
    OUT.write_text("".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in records), encoding="utf-8")
    REPORT.write_text(build_report(local, capell_summary), encoding="utf-8")
    elapsed = time.monotonic() - started
    assert OUT.exists() and REPORT.exists()
    print(
        "L23 half-sieve: PROVED barrier; "
        f"local_pairs={local['bad_root_local_pairs_stored']}, "
        f"capell={capell_summary['total_rows_proved_nonsquare']}/353, "
        f"refusals={capell_summary['certificate_refusals_after_local_checks']}, {elapsed:.3f}s"
    )


if __name__ == "__main__":
    main()
