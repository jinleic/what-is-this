#!/usr/bin/env python3
"""Prove and replay aligned-class existence for the L12b escape family.

The table-backed l12_class.l12_class_cert is generalized here only in its
inputs: (a, z, f, eps, q1) are supplied directly.  Its mathematics is
unchanged and remains the h10q kernel's Taylor-exponent certificate.  Every
JSON row produced from that generalization is labelled explicitly.

Replay from the workspace repository root:

    nice -n 19 python3 math/h10q/l19_classexist.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h10q  # noqa: E402  (the stipulated arithmetic authority)


SAMPLE = HERE / "data" / "l12b_class_sample.json"
OUT = HERE / "data" / "l19_classexist.jsonl"
REPORT = Path("/tmp/l19_classexist.md")
A_POOL = (1, 3, 5, 7)
PACE_EVERY = 200
PACE_SECONDS = 0.01
CANONICAL_NONEMPTY_HYPOTHESIS = (
    "constructed clean class only: a=1, A=5, tau=3/5, eps=+1, "
    "f prime with v_f(z)>0 and (f|5)=-1; q1 is a unit modulo "
    "M=4*A*prod(S), and (2*f*q1|p)=+1 for every odd p in S\\{f}"
)
CANONICAL_COLLISION_HYPOTHESIS = (
    "canonical wall test: a=1, A=5, tau=3/5, eps=+1, f prime "
    "with v_f(z)>0 and (f|5)=+1; split every reduced q1 residue "
    "by (q1|5)=+1/-1 and test Hilb_5 and wild_q1 simultaneously"
)
RECORDED_CLASS_HYPOTHESIS = (
    "recorded _L12_ESCAPE base: exact h10q Hilbert symbols at its "
    "recorded prime q1 plus the Taylor-exponent class modulus; the "
    "recorded q1 need not satisfy the stronger square-d system (R)"
)
UNIFORM_HYPOTHESIS = (
    "f is an odd numerator prime of z; a is odd with "
    "(1+4*a^2|f)=-1; eps=+1; q1 is a unit modulo "
    "M=4*A*prod(S), and (2*f*q1|p)=+1 for every odd "
    "p in S\\{f}"
)


class Pacer:
    """Apply the requested low-duty-cycle pause to candidate loops."""

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


def place_text(place) -> str:
    return "inf" if place == h10q.OO else str(place)


def symbol_map(symbols: dict) -> dict[str, int]:
    return {place_text(place): int(value) for place, value in symbols.items()}


def evaluate(coefficients: list[F], value: F) -> F:
    return sum(
        coefficient * value**degree
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    )


def quadratic_character(value: int, prime: int) -> int:
    """Legendre character with an explicit zero branch.

    h10q.legendre intentionally takes the character of a p-adic unit part, so
    it has no zero value.  Residue enumeration needs the ordinary character.
    """
    residue = value % prime
    if residue == 0:
        return 0
    return h10q.legendre(F(residue), prime)


def euler_phi(value: int) -> int:
    answer = value
    for prime in h10q.factorint(value):
        answer = answer // prime * (prime - 1)
    return answer


def controlled_set(a: int, f: int) -> tuple[F, F, F, list[int]]:
    aa = F(a)
    A = 1 + 4 * aa * aa
    tau = (1 + 2 * aa * aa) / A
    delta = 1 - A * tau * tau
    alpha = -delta * A
    S = sorted(
        {2, 3, 5, 7, f}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
    )
    return A, delta, alpha, S


def residue_modulus(A: F, S: list[int]) -> int:
    assert A.denominator == 1
    return 4 * A.numerator * math.prod(S)


def finite_excluded_support(a: F, z: F, delta: F, D: F, f: int) -> list[int]:
    """Exact finite exclusion set, used only for the 103 authority comparisons."""
    return sorted(
        h10q._l10_supp(D)
        | h10q._l10_supp(z.numerator)
        | h10q._l10_supp(delta)
        | h10q._l10_supp(a)
        | {f}
    )


def generalized_escape_cert(
    a_value: int | F,
    z_value: F,
    f: int,
    eps: int,
    q1: int,
) -> dict | None:
    """Direct-input generalization of the L12 exponent-lemma certificate.

    GENERALIZATION LABEL: the table lookup and the restrictions a=1 and a
    pre-registered f are removed.  The certificate, frozen set, Taylor
    exponents, moving character and asymptotic-real calculation are otherwise
    exactly the h10q.py::_l10_class_cert / l12_class construction.
    """
    a = F(a_value)
    z = F(z_value)
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    if delta == 0:
        return None
    alpha = -delta * A
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    if a == 0 or z == 0 or D == 0 or eps not in (1, -1):
        return None
    if f == 2 or not h10q._is_prime(f) or not h10q._is_prime(q1):
        return None
    if h10q.vp(a, 2) != 0:
        return None

    S = sorted(
        {2, 3, 5, 7, f}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
    )
    if q1 in S:
        return None
    if any(h10q.vp(value, q1) != 0 for value in (delta, A, a, F(f))):
        return None
    if h10q.vp(z, f) <= 0:
        return None
    if any(h10q.vp(value, f) != 0 for value in (delta, A, a, D)):
        return None

    s = (a - 1) / 2
    P = h10q._l10_P(a, Z, D, A, delta, s)
    b0 = F(eps * f * q1)
    ks = {prime: h10q._l10_exponent(P, b0, prime) for prime in S}
    N = 8
    for prime, exponent in ks.items():
        N = math.lcm(N, prime**exponent)
    N = math.lcm(N, 4 * A.numerator * A.denominator)
    if math.gcd(q1, N) != 1:
        return None

    c0 = h10q._sun_h(a, b0, Z)
    if c0 is None:
        return None
    M0 = 16 - delta * c0 * c0 - 32 * A * b0 * s * s
    if M0 == 0:
        return None
    Pb0 = evaluate(P, b0)
    assert M0 == Pb0 / (b0**4 * D * D * A * A)
    x0, d0 = alpha * M0, alpha * 2 * b0
    syms = {prime: h10q.hilbert(x0, d0, prime) for prime in S}

    nonzero = [index for index, coefficient in enumerate(P) if coefficient]
    lead = P[nonzero[-1]]
    Q0 = 1 + max(
        (abs(coefficient / lead) for coefficient in P[: nonzero[-1]]),
        default=F(0),
    )
    sign_x = (
        (1 if alpha > 0 else -1)
        * (1 if lead > 0 else -1)
        * eps ** nonzero[-1]
    )
    sign_d = 1 if alpha * eps > 0 else -1
    syms[h10q.OO] = -1 if sign_x < 0 and sign_d < 0 else 1
    syms["q1"] = h10q.legendre(A, q1)
    return {
        "ok": all(value == 1 for value in syms.values()),
        "N": N,
        "S": S,
        "ks": ks,
        "syms": syms,
        "Q0": Q0,
        "f": f,
        "eps": eps,
        "q": q1,
        "a": a,
        "A": A,
        "tau": tau,
        "delta": delta,
        "alpha": alpha,
        "D": D,
        "kernel_identity": True,
        "generalized_kernel_cert": True,
    }


def cert_record(cert: dict) -> dict:
    return {
        "a": frac_text(cert["a"]),
        "A": frac_text(cert["A"]),
        "tau": frac_text(cert["tau"]),
        "eps": int(cert["eps"]),
        "f": int(cert["f"]),
        "q1": int(cert["q"]),
        "N": str(cert["N"]),
        "S": list(cert["S"]),
        "ks": {str(prime): int(exponent) for prime, exponent in cert["ks"].items()},
        "syms": symbol_map(cert["syms"]),
        "Q0": frac_text(cert["Q0"]),
        "ok": bool(cert["ok"]),
        "kernel_identity": bool(cert["kernel_identity"]),
        "generalized_kernel_cert": True,
    }


def forcing_condition(a: int, f: int, q_residue: int, S: list[int]) -> bool:
    return all(
        quadratic_character(2 * f * q_residue, prime) == 1
        for prime in S
        if prime not in (2, f)
    )


def moving_character_from_residue(A: int, residue: int) -> int:
    """Compute (A|q) from q residues at p|A using reciprocity.

    Every p|A is 1 mod 4 in the construction, so (p|q)=(q|p).
    """
    value = 1
    for prime, exponent in h10q.factorint(A).items():
        assert prime % 4 == 1
        if exponent % 2:
            value *= quadratic_character(residue, prime)
    return value


def canonical_profile(f: int, direct: bool) -> dict:
    """Count the canonical a=1 residue system, directly when requested."""
    A, _delta, _alpha, S = controlled_set(1, f)
    assert A == 5
    M = residue_modulus(A, S)
    phi = euler_phi(M)
    f_character = quadratic_character(f, 5)
    assert f_character in (-1, 1)
    imposed_places = [prime for prime in S if prime not in (2, f)]
    forcing_count = phi // (2 ** len(imposed_places))
    forcing_branch = -f_character
    branches = {
        "+1": {
            "q1_character_mod_5": 1,
            "unit_residues": phi // 2,
            "frozen_forcing_residues": forcing_count if forcing_branch == 1 else 0,
            "simultaneous_all_plus_residues": (
                forcing_count if forcing_branch == 1 else 0
            ),
            "hilbert_5": -f_character,
            "wild_q1": 1,
        },
        "-1": {
            "q1_character_mod_5": -1,
            "unit_residues": phi // 2,
            "frozen_forcing_residues": forcing_count if forcing_branch == -1 else 0,
            "simultaneous_all_plus_residues": 0,
            "hilbert_5": f_character,
            "wild_q1": -1,
        },
    }
    examples: list[int] = []
    direct_allowed = None
    direct_units = None
    direct_branches = None
    if direct:
        direct_allowed = 0
        direct_units = 0
        direct_branches = {
            "+1": {"unit_residues": 0, "frozen_forcing_residues": 0,
                    "simultaneous_all_plus_residues": 0},
            "-1": {"unit_residues": 0, "frozen_forcing_residues": 0,
                    "simultaneous_all_plus_residues": 0},
        }
        for residue in range(1, M + 1):
            PACER.tick()
            if math.gcd(residue, M) != 1:
                continue
            direct_units += 1
            q_character = quadratic_character(residue, 5)
            branch = "+1" if q_character == 1 else "-1"
            direct_branches[branch]["unit_residues"] += 1
            forced = forcing_condition(1, f, residue, S)
            if forced:
                direct_branches[branch]["frozen_forcing_residues"] += 1
            moving = moving_character_from_residue(5, residue)
            simultaneous = forced and moving == 1
            if simultaneous:
                direct_allowed += 1
                direct_branches[branch]["simultaneous_all_plus_residues"] += 1
                if len(examples) < 12:
                    examples.append(residue % M)
        assert direct_units == phi
        assert direct_allowed == (
            forcing_count if f_character == -1 else 0
        )
        for branch in branches:
            for field in (
                "unit_residues",
                "frozen_forcing_residues",
                "simultaneous_all_plus_residues",
            ):
                assert direct_branches[branch][field] == branches[branch][field]

    return {
        "exact_hypothesis_clause": (
            CANONICAL_NONEMPTY_HYPOTHESIS
            if f_character == -1
            else CANONICAL_COLLISION_HYPOTHESIS
        ),
        "a": 1,
        "A": 5,
        "eps": 1,
        "f": f,
        "S": S,
        "modulus_definition": "M=4*A*prod(S), with S a set of distinct primes",
        "M": M,
        "M_factorization": {
            str(prime): exponent for prime, exponent in h10q.factorint(M).items()
        },
        "unit_residues": phi,
        "f_character_mod_5": f_character,
        "criterion": "nonempty iff (f|5)=-1",
        "criterion_label": "PROVED",
        "branches_by_q1_character_mod_5": branches,
        "simultaneous_all_plus_residues": (
            forcing_count if f_character == -1 else 0
        ),
        "direct_enumeration": bool(direct),
        "direct_residues_visited": M if direct else 0,
        "direct_allowed_examples": examples,
    }


def parameter_character_check(f: int) -> dict:
    counts = Counter()
    suitable: list[int] = []
    character_sum = 0
    for residue in range(f):
        PACER.tick()
        A_mod_f = (1 + 4 * residue * residue) % f
        character = quadratic_character(A_mod_f, f)
        character_sum += character
        counts[character] += 1
        if character == -1:
            suitable.append(residue)
    zero_expected = 2 if f % 4 == 1 else 0
    negative_expected = (f - zero_expected + 1) // 2
    assert character_sum == -1
    assert counts[0] == zero_expected
    assert counts[-1] == negative_expected and suitable
    return {
        "f": f,
        "sum_character_1_plus_4r2": character_sum,
        "counts": {
            "nonresidue": counts[-1],
            "zero": counts[0],
            "residue": counts[1],
        },
        "expected_nonresidue_count": (
            "(f-1)/2" if f % 4 == 1 else "(f+1)/2"
        ),
        "first_nonresidue_residue": suitable[0],
        "label": "PROVED_BY_EXACT_ENUMERATION_FOR_THIS_f",
    }


def two_adic_residue_check() -> list[dict]:
    """Exhaust the 16 odd (a,b) residue pairs in the 2-adic proof."""
    rows = []
    for a_mod_8 in (1, 3, 5, 7):
        s = (a_mod_8 - 1) // 2
        A_mod_8 = (1 + 4 * a_mod_8 * a_mod_8) % 8
        assert A_mod_8 == 5
        for b_mod_8 in (1, 3, 5, 7):
            unit_x = (1 - 2 * A_mod_8 * s * s * b_mod_8) % 8
            symbol = h10q.hilbert(64 * unit_x, 8 * b_mod_8, 2)
            assert symbol == 1
            rows.append(
                {
                    "a_mod_8": a_mod_8,
                    "b_mod_8": b_mod_8,
                    "unit_x_mod_8": unit_x,
                    "unit_d_mod_8": b_mod_8,
                    "hilbert_2": symbol,
                }
            )
    assert len(rows) == 16
    return rows


def choose_pool_a(f: int) -> int:
    for a in A_POOL:
        A = 1 + 4 * a * a
        if A % f and quadratic_character(A, f) == -1:
            return a
    raise AssertionError((f, "bounded a pool empty"))


def find_clean_safe_cert(a: int, z: F, f: int) -> tuple[dict, int]:
    """Find one small kernel-proved q1 in the theorem's nonempty CRT system."""
    A, _delta, _alpha, S = controlled_set(a, f)
    assert A.denominator == 1
    Z = z**3
    D = 1 - Z - F(a * a) * Z * Z
    attempts = 0
    for q1 in range(3, 10000, 2):
        PACER.tick()
        attempts += 1
        if q1 in S or math.gcd(q1, math.prod(S)) != 1:
            continue
        if not forcing_condition(a, f, q1, S):
            continue
        if not h10q._is_prime(q1):
            continue
        if h10q.vp(z, q1) != 0 or h10q.vp(D, q1) != 0:
            continue
        assert h10q.legendre(A, q1) == 1
        cert = generalized_escape_cert(a, z, f, 1, q1)
        if cert is not None and cert["ok"]:
            return cert, attempts
    raise AssertionError((a, z, f, "no safe prime below 10000"))


def sample_comparison(sample_row: dict, cert: dict, excluded: list[int]) -> dict:
    checks = {
        "ok": bool(sample_row["ok_cert"]) == bool(cert["ok"]),
        "f": int(sample_row["f"]) == int(cert["f"]),
        "eps": int(sample_row["eps"]) == int(cert["eps"]),
        "q1": int(sample_row["q"]) == int(cert["q"]),
        "N": int(sample_row["N"]) == int(cert["N"]),
        "S": list(sample_row["S"]) == list(cert["S"]),
        "ks": {
            int(prime): int(exponent) for prime, exponent in sample_row["ks"].items()
        }
        == cert["ks"],
        "Q0": F(sample_row["Q0"]) == cert["Q0"],
        "excluded": list(sample_row["excluded"]) == excluded,
        "sample_failures_empty": not sample_row["failures"],
    }
    assert all(checks.values()), checks
    return {"all": True, "checks": checks}


def canonical_collision_cells() -> list[tuple[int, tuple[int, int]]]:
    cells = sorted(h10q._L13_ESCAPE2)
    assert len(cells) == 60
    for w, unit in cells:
        z = F(w) * F(*unit)
        odd_numerator = [
            prime
            for prime in h10q.factorint(abs(z.numerator))
            if prime != 2
        ]
        assert odd_numerator == [w]
        assert quadratic_character(w, 5) == 1
    return cells


def write_report(summary: dict, target_rows: list[dict], collision_rows: list[dict]) -> None:
    canonical = summary["canonical_103"]
    uniform = summary["uniform_theorem"]
    collision_groups: dict[int, list[list[int]]] = {}
    for row in collision_rows:
        collision_groups.setdefault(row["cell"][0], []).append(row["cell"][1])

    lines = [
        "# L19 — class existence for the escape family",
        "",
        "**Class-existence verdict — PROVED, uniformly in the cell.  Member-existence verdict — still CONDITIONAL on Schinzel H.**  Class existence for the escape family reduces to Dirichlet + nonemptiness of an explicit residue system.  The residue system is always nonempty after choosing the square-branch parameter `a` by the character-sum lemma below.  Thus clause (ii), existence of one aligned class, is removed entirely; clause (i), a qualifying member of that class, is not.",
        "",
        "This result generalizes the *inputs* of the table-backed L12 certificate: `a`, `f`, and `q1` are accepted directly rather than through `_L12_ESCAPE`.  Every affected JSON row is labelled `generalized_kernel_cert=true`.  The Taylor exponent, frozen set, Hilbert symbols, moving character, real sign, and kernel identity are still computed only by `h10q.py`.",
        "",
        "## 1. Exact construction",
        "",
        "Let `(w,z)` be a cell and choose any odd prime `f` with `v_f(z)>0` (one may always take `f=w`).  Choose an odd integer `a` such that",
        "",
        "    A = 1 + 4 a^2,                 (A|f) = -1,",
        "",
        "where `(A|f)` is the ordinary Legendre symbol, so in particular `f` divides neither `a` nor `A`.  Take the square branch",
        "",
        "    tau=(1+2a^2)/A,  delta=-4a^4/A,  alpha=(2a^2)^2,  eps=+1.",
        "Taking `f=w` never fails: the cell hypothesis gives `v_w(z)>=1`; `(A|w)=-1` gives `w` coprime to `a*A` and hence to `delta`; and `v_w(Z)>=3` gives `D=1-Z-a^2Z^2=1 (mod w)`.  Thus every fixed-factor guard holds.  There is no denominator, gcd, or branch exception, including `w=5` (take `a=3`, `A=37`).",
        "",
        "",
        "Put",
        "",
        "    S={2,3,5,7} union supp(alpha) union supp(delta) union {f},",
        "    M=4*A*prod_{p in S} p.",
        "",
        "The explicit residue system for `q1` is",
        "",
        "    q1 in (Z/MZ)^x,",
        "    (2*f*q1|p)=+1  for every odd p in S with p != f.          (R)",
        "",
        "Each line of (R) is a nonempty half of the units at one distinct prime.  CRT therefore gives exactly",
        "",
        "    phi(M) / 2^{#{odd p in S: p != f}}",
        "",
        "residue classes modulo `M`; in particular the system is nonempty.  Dirichlet gives infinitely many positive primes `q1` in any one of them.  Discarding primes dividing the fixed cell data and the finitely many roots/degenerate values removes only finitely many choices.",
        "",
        "For such a prime, set `b0=f*q1`, form the kernel polynomial `P`, put",
        "",
        "    k_p = _l10_exponent(P,b0,p),",
        "    N = lcm(8, {p^k_p : p in S}, 4*A).",
        "",
        "Then `gcd(q1,N)=1`, so Dirichlet also guarantees infinitely many prime members of the progression `Q=q1 (mod N)` (the simultaneous degree-8 prime value is still the Schinzel input).",
        "",
        "## 2. Why every required symbol is +1",
        "",
        "Write `Z=z^3`, `D=1-Z-a^2 Z^2`, `s=(a-1)/2`, and use the kernel polynomial",
        "",
        "    P(b)=16D^2A^2b^4-delta*a^4*Z^4*Ng(b)^2-32A^3s^2D^2b^5.",
        "",
        "Because `alpha`, `D^2`, `A^2`, and `b^4` are squares, `x` has the square class of `P(b)`.",
        "",
        "* At `f`, `v_f(b)=1`, `v_f(Z)>=3`, `D` is a unit, and `Ng(b)=-A (mod f)`.  The three displayed terms have valuations `4`, at least `12`, and at least `5`; the unique first term is a square.  Hence `x` is a square in `Q_f`.",
        "* At `2`, the first term has uniquely least valuation, but for `s` odd its unit need not be a square; the symbol is nevertheless always `+1`.  Exactly, `v_2(x)=6`, `v_2(d)=3`, `u_d=b (mod 8)`, and `u_x=1-2*A*s^2*b (mod 8)`.  If `s` is even then `u_x=1`; if `s` is odd, `A=5 (mod 8)` and `u_x=1-2b`.  Substitution of `b=1,3,5,7 (mod 8)` in the kernel's 2-adic Hilbert formula gives exponent `0` in all four cases.",
        "* At every other odd frozen `p`, (R) makes `2*f*q1` a square unit; because `alpha` is itself a square, `d=2*alpha*f*q1` is a square (its valuation need not be zero when `p|a`).  Thus `(x,d)_p=+1` without any condition on the valuation of `P`.",
        "* At infinity, `d>0`, so the real Hilbert symbol is `+1`.",
        "* At the moving prime, the kernel gives `(x,d)_q1=(A|q1)`.  Every prime divisor `p|A` is `1 mod 4`, since `(2a)^2=-1 (mod p)`.  Reciprocity and (R) give",
        "",
        "      (A|q1) = product_{p^e||A}(q1|p)^e",
        "             = (2f|A) = (2|A)(f|A) = (-1)(-1) = +1,",
        "",
        "  because `A=5 mod 8` and `(f|A)=(A|f)=-1`.  This is the compatibility step that rules out a hidden wall.",
        "",
        "The exponent lemma then freezes the finite symbols at `S` throughout the `N`-class, while `4A|N` freezes the moving character.  Emergent places outside `S` are not frozen; those remain on the member/Schinzel side.  This proves class existence.",
        "",
        "## 3. Why a always exists",
        "",
        "For every odd prime `f`, the standard quadratic-character sum is",
        "",
        "    sum_{r mod f} (1+4r^2|f) = -1.",
        "",
        "There are two zero residues when `f=1 mod 4` and none when `f=3 mod 4`.  Consequently the number of nonresidue choices is `(f-1)/2` in the first case and `(f+1)/2` in the second, always positive.  An integer lift can be made odd because adding the odd number `f` flips parity.  Thus the construction is uniform in the cell and does not require primes of the form `1+4a^2`; composite `A` is allowed and its actual prime support is frozen.",
        "This is the same uniform **parameter-selection lemma** used by the sibling `tau=0` construction in its `w=3 mod 4` scope: both choose a moving `A=1+4a^2` with `(A|w)=-1` from the character sum.  The branch-specific local arguments are different.  The `tau=0` proof additionally chooses the odd lift `a=0 (mod rad(w^3+2))` and uses its own product identity; the square-branch proof here instead uses the explicit square-`d` residue system (R).",
        "",
        "The replay also checks a genuinely composite value `A=325=5^2*13` (`a=9`, `f=z=7`): the generalized certificate is `ok`, so the theorem is not relying on the prime-`A` grid examples.",
        "",
        "",
        "## 4. Canonical a=1 criterion and the 103 rows",
        "",
        "For `a=1`, `A=5`, the construction specializes to",
        "",
        "    Hilb_5=-(f|5)(q1|5),       wild=(5|q1)=(q1|5).",
        "",
        "Hence the canonical system is nonempty iff `(f|5)=-1`; then the compatible branch is `(q1|5)=+1`.  If `(f|5)=+1`, the two symbols are opposite for every admissible `q1`: this is the exact 5-wall collision, not a bounded no-hit.",
        "",
        f"The replay checked all **{canonical['rows']}/{canonical['rows']}** `_L12_ESCAPE` rows.  It directly enumerated every residue modulo `M=4*A*prod(S)` for each of the {canonical['distinct_f']} distinct canonical `f` values ({canonical['direct_residues_visited']} candidates, paced every {PACE_EVERY}); every system was nonempty.  It independently rebuilt the generalized certificate at the recorded `q1` and matched `N,S,k_p,Q0`, the finite exclusion set, and `ok` against `data/l12b_class_sample.json`: **{canonical['recorded_certificate_matches']}/{canonical['rows']}**.  It also constructed a clean prime base in the stronger square-`d` subsystem at **{canonical['constructed_clean_certificates']}/{canonical['rows']}** rows.  The recorded base itself lies in that deliberately stronger subsystem in {canonical['recorded_in_forcing_subsystem']}/{canonical['rows']} rows; the other recorded rows remain valid by their exact Hilbert certificate.",
        "",
        "Per-`f` direct counts:",
        "",
        "| f | M | units | `(q1|5)=+1`: system-(R) all-plus classes | `(q1|5)=-1`: system-(R) all-plus classes |",
        "|---:|---:|---:|---:|---:|",
    ]
    seen_f: dict[int, dict] = {}
    for row in target_rows:
        seen_f.setdefault(row["f"], row["residue_system"])
    for f, profile in sorted(seen_f.items()):
        branches = profile["branches_by_q1_character_mod_5"]
        lines.append(
            f"| {f} | {profile['M']} | {profile['unit_residues']} | "
            f"{branches['+1']['simultaneous_all_plus_residues']} | "
            f"{branches['-1']['simultaneous_all_plus_residues']} |"
        )

    lines.extend(
        [
            "",
            "**103/103 verdict: PROVED.  Canonical collisions among these recorded choices: none.**",
            "",
            "## 5. Exact canonical collisions and their repair",
            "",
            "The canonical `a=1` system is empty at the following 60 walled grid cells.  Each has only the odd numerator prime `f=w` and `(w|5)=+1`.  This list is an exact theorem-level canonical obstruction:",
            "",
            "| w | exact unit-coordinate cells u | canonical system-(R) all-plus counts (+/-) | repaired a | A |",
            "|---:|:---|:---|---:|---:|",
        ]
    )
    collision_by_w: dict[int, dict] = {}
    for row in collision_rows:
        collision_by_w.setdefault(row["cell"][0], row)
    for w in sorted(collision_groups):
        first = collision_by_w[w]
        profile = first["canonical_residue_system"]
        plus_count = profile["branches_by_q1_character_mod_5"]["+1"][
            "simultaneous_all_plus_residues"
        ]
        minus_count = profile["branches_by_q1_character_mod_5"]["-1"][
            "simultaneous_all_plus_residues"
        ]
        units = ", ".join(str(tuple(unit)) for unit in collision_groups[w])
        lines.append(
            f"| {w} | `{units}` | `{plus_count}/{minus_count}` | "
            f"{first['constructed_class']['a']} | {first['constructed_class']['A']} |"
        )

    lines.extend(
        [
            "",
            f"The uniform construction repairs all **{len(collision_rows)}/{len(collision_rows)}**: choose `a=5` for `w=11`, `a=3` for `w in {{19,31,59,79}}`, and `a=7` for `w=71`, then apply (R).  The exact generalized certificate is `ok` at every row.  Therefore **no cell is empty for every admissible `a` in the bounded pool {A_POOL}**; the empty-all-pool list is `[]`.  The canonical collision is real, but it disappears when `a` changes because `A` changes.",
            "",
            "## 6. Mechanical totals and conditional consequence",
            "",
            f"The script additionally applied the uniform construction with `f=w` to the entire frozen grid: **{uniform['grid_certificates_ok']}/{uniform['grid_cells']}** exact generalized certificates are `ok`, with every chosen base prime outside the finite moving-prime exclusions.  The grid needs only `a in {list(A_POOL)}`.  All {uniform['grid_prime_character_sum_checks']} distinct grid-prime character tables satisfy the exact `-1` sum and predicted nonresidue count.",
            "",
            "**Precise upgrade.** The aligned-class-existence clause is no longer a hypothesis for the escape family; since every cell has the odd numerator prime `w`, it is no longer a hypothesis for any cell.  What remains conditional is Schinzel's Hypothesis H for the chosen pair `{q1+Nt,G(t)}` (plus verification of the pair's irreducibility/admissibility when invoking Schinzel).  This report does **not** prove a prime value of the degree-8 polynomial, an emergent-free member, the L6 assembly lemma, or Hilbert's Tenth over `Q`.",
            "",
            f"Refusals: {summary['refusals']} (none were used as evidence).  Verification wall-clock: **{summary['wall_seconds']:.3f} s**.  Candidate loop: {summary['pacing']['candidates']} candidates, {summary['pacing']['sleeps']} sleeps of {PACE_SECONDS} s.",
            "",
            "Artifacts: `math/h10q/l19_classexist.py`, `math/h10q/data/l19_classexist.jsonl`, `/tmp/l19_classexist.md`.",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    sample_payload = json.loads(SAMPLE.read_text(encoding="utf-8"))
    assert len(sample_payload) == 103
    target_items = sorted(h10q._L12_ESCAPE.items())
    assert len(target_items) == 103
    target_f = sorted({row[0] for _cell, row in target_items})

    direct_profiles = {f: canonical_profile(f, direct=True) for f in target_f}
    target_rows: list[dict] = []
    recorded_in_forcing = 0
    direct_residues_visited = sum(
        profile["direct_residues_visited"] for profile in direct_profiles.values()
    )
    for index, (cell, table_row) in enumerate(target_items):
        PACER.tick()
        w, unit = cell
        f, eps, recorded_q1 = table_row
        assert eps == 1
        z = F(w) * F(*unit)
        profile = direct_profiles[f]
        assert profile["simultaneous_all_plus_residues"] > 0
        assert quadratic_character(f, 5) == -1

        recorded_cert = generalized_escape_cert(1, z, f, eps, recorded_q1)
        assert recorded_cert is not None and recorded_cert["ok"]
        excluded = finite_excluded_support(
            F(1), z, recorded_cert["delta"], recorded_cert["D"], f
        )
        sample_key = str((w, unit))
        assert sample_key in sample_payload
        comparison = sample_comparison(
            sample_payload[sample_key], recorded_cert, excluded
        )
        in_forcing = forcing_condition(1, f, recorded_q1, recorded_cert["S"])
        recorded_in_forcing += int(in_forcing)

        clean_cert, search_attempts = find_clean_safe_cert(1, z, f)
        assert clean_cert["ok"]
        assert forcing_condition(1, f, clean_cert["q"], clean_cert["S"])
        assert h10q.vp(z, clean_cert["q"]) == 0
        assert h10q.vp(clean_cert["D"], clean_cert["q"]) == 0

        target_rows.append(
            {
                "type": "canonical-escape-cell",
                "label": "PROVED",
                "family": "ESC",
                "cell": [w, list(unit)],
                "z": frac_text(z),
                "a": "1",
                "A": "5",
                "tau": "3/5",
                "eps": eps,
                "f": f,
                "exact_hypothesis_clauses": {
                    "constructed_clean_class": CANONICAL_NONEMPTY_HYPOTHESIS,
                    "recorded_class": RECORDED_CLASS_HYPOTHESIS,
                },
                "residue_system": profile,
                "recorded_class": cert_record(recorded_cert),
                "recorded_class_sample_comparison": comparison,
                "recorded_q1_in_stronger_forcing_subsystem": in_forcing,
                "constructed_clean_class": cert_record(clean_cert),
                "constructed_q1_search_candidates": search_attempts,
                "canonical_verdict": "PROVED_NONEMPTY",
                "collision": None,
                "source": "generalized h10q exponent-lemma certificate; compared to data/l12b_class_sample.json",
                "generalized_kernel_cert": True,
            }
        )

    assert len(target_rows) == 103
    assert all(row["canonical_verdict"] == "PROVED_NONEMPTY" for row in target_rows)

    collision_cells = canonical_collision_cells()
    collision_f = sorted({w for w, _unit in collision_cells})
    collision_profiles = {f: canonical_profile(f, direct=False) for f in collision_f}

    grid = sorted(h10q._l9_grid())
    assert len(grid) == 353
    grid_primes = sorted({w for w, _unit in grid})
    character_rows = [parameter_character_check(f) for f in grid_primes]
    pool_choice = {f: choose_pool_a(f) for f in grid_primes}
    assert set(pool_choice.values()) <= set(A_POOL)

    grid_certs: dict[tuple[int, tuple[int, int]], dict] = {}
    grid_search_attempts = 0
    grid_q_max = 0
    grid_a_counts = Counter()
    for cell in grid:
        PACER.tick()
        w, unit = cell
        a = pool_choice[w]
        z = F(w) * F(*unit)
        cert, attempts = find_clean_safe_cert(a, z, w)
        assert cert["ok"]
        assert h10q.vp(z, cert["q"]) == 0
        assert h10q.vp(cert["D"], cert["q"]) == 0
        grid_certs[cell] = cert
        grid_search_attempts += attempts
        grid_q_max = max(grid_q_max, cert["q"])
        grid_a_counts[a] += 1
    assert len(grid_certs) == 353

    collision_rows: list[dict] = []
    for cell in collision_cells:
        PACER.tick()
        w, unit = cell
        canonical = collision_profiles[w]
        assert canonical["simultaneous_all_plus_residues"] == 0
        cert = grid_certs[cell]
        assert cert["a"] == pool_choice[w]
        collision_rows.append(
            {
                "type": "canonical-collision-cell",
                "label": "PROVED",
                "family": "ESC",
                "cell": [w, list(unit)],
                "z": frac_text(F(w) * F(*unit)),
                "f": w,
                "canonical_exact_hypothesis_clause": CANONICAL_COLLISION_HYPOTHESIS,
                "canonical_residue_system": canonical,
                "canonical_verdict": "PROVED_EMPTY_BY_FIVE_WALL",
                "overall_class_existence_verdict": "PROVED_BY_ALTERNATE_A",
                "uniform_exact_hypothesis_clause": UNIFORM_HYPOTHESIS,
                "constructed_class": cert_record(cert),
                "bounded_a_pool": list(A_POOL),
                "empty_for_all_a_in_bounded_pool": False,
                "source": "generalized h10q exponent-lemma certificate",
                "generalized_kernel_cert": True,
            }
        )
    assert len(collision_rows) == 60

    two_adic_rows = two_adic_residue_check()
    composite_cert, composite_attempts = find_clean_safe_cert(9, F(7), 7)
    assert composite_cert["A"] == 325
    assert h10q.factorint(325) == {5: 2, 13: 1}
    assert composite_cert["ok"]


    elapsed = time.perf_counter() - started
    summary = {
        "type": "summary",
        "label": "PROVED",
        "verdict": "class existence for the escape family reduces to Dirichlet + nonemptiness of an explicit residue system; the system is uniformly nonempty after choosing a",
        "class_existence_verdict": "PROVED_UNIFORMLY_FOR_EVERY_CELL",
        "member_existence_verdict": "CONDITIONAL_ON_SCHINZEL_H",
        "generalized_certificate": {
            "enabled": True,
            "label_on_every_affected_row": "generalized_kernel_cert=true",
            "change": "direct (a,z,f,eps,q1) inputs replace the _L12_ESCAPE table lookup; h10q Taylor/Hilbert mathematics unchanged",
        },
        "canonical_103": {
            "rows": 103,
            "proved_nonempty": 103,
            "collisions": [],
            "distinct_f": len(target_f),
            "f_values": target_f,
            "direct_residue_profiles": len(direct_profiles),
            "direct_residues_visited": direct_residues_visited,
            "recorded_certificate_matches": 103,
            "recorded_in_forcing_subsystem": recorded_in_forcing,
            "constructed_clean_certificates": 103,
        },
        "canonical_collision_shapes": {
            "count": 60,
            "cells": [[w, list(unit)] for w, unit in collision_cells],
            "canonical_verdict": "PROVED_EMPTY_BY_FIVE_WALL",
            "overall_verdict": "PROVED_BY_ALTERNATE_A",
        },
        "uniform_theorem": {
            "label": "PROVED",
            "exact_hypothesis_clause": UNIFORM_HYPOTHESIS,
            "parameter_lemma": "sum_r (1+4r^2|f)=-1, so a nonresidue A exists for every odd prime f",
            "f_equals_w_always_admissible": True,
            "f_equals_w_proof": "v_w(z)>0; (A|w)=-1 gives w coprime to a*A*delta; D=1 mod w",
            "shared_parameter_lemma": "same character-sum selection of A as the tau=0 sibling; branch-specific local-symbol systems differ",
            "residue_system": "q1 unit mod M=4*A*prod(S); (2*f*q1|p)=+1 for odd p in S\\{f}",
            "residue_count": "phi(M)/2^#{odd p in S, p!=f}",
            "dirichlet_input": "one nonempty reduced residue class modulo M",
            "moving_symbol_identity": "(A|q1)=(2f|A)=(2|A)(f|A)=+1",
            "two_adic_residue_checks": len(two_adic_rows),
            "two_adic_rows": two_adic_rows,
            "composite_A_check": {
                "label": "PROVED",
                "a": 9,
                "A": 325,
                "A_factorization": {"5": 2, "13": 1},
                "f": 7,
                "z": "7",
                "q1_search_candidates": composite_attempts,
                "certificate": cert_record(composite_cert),
            },
            "grid_cells": len(grid),
            "grid_certificates_ok": len(grid_certs),
            "grid_prime_character_sum_checks": len(character_rows),
            "grid_character_rows": character_rows,
            "bounded_a_pool": list(A_POOL),
            "bounded_a_counts_by_cell": {
                str(a): count for a, count in sorted(grid_a_counts.items())
            },
            "empty_for_all_a_in_bounded_pool": 0,
            "empty_cells": [],
            "maximum_constructed_clean_q1": grid_q_max,
            "constructed_q1_search_candidates": grid_search_attempts,
        },
        "conditional_consequence": {
            "removed": "existence of one verified aligned class per cell",
            "remains": "Schinzel H for the chosen linear/degree-8 pair and the hypotheses required to invoke it",
            "not_claimed": [
                "an unconditional prime value of the degree-8 polynomial",
                "an emergent-free member without Schinzel H",
                "the L6 assembly lemma",
                "Hilbert's Tenth Problem over Q",
            ],
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
            "kernel": "math/h10q/h10q.py",
            "recorded_escape_table": "h10q.py::_L12_ESCAPE",
            "recorded_sample": "math/h10q/data/l12b_class_sample.json",
            "script": "math/h10q/l19_classexist.py",
        },
    }

    rows = [
        {
            "type": "meta",
            "label": "PROVED",
            "schema": "l19-classexist-v1",
            "source_script": "math/h10q/l19_classexist.py",
            "kernel_authority": "math/h10q/h10q.py",
            "exact_hypothesis_clause": UNIFORM_HYPOTHESIS,
            "refusals_are_evidence": False,
        },
        *target_rows,
        *collision_rows,
        summary,
    ]
    with OUT.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(summary, target_rows, collision_rows)

    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
    print(f"wrote {OUT} ({len(rows)} rows)")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
