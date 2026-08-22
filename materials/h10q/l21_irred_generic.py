#!/usr/bin/env python3
"""Generic and cellwise irreducibility audit for the constructed L21 branch.

Replay from the workspace repository root:

    nice -n 19 python3 math/h10q/l21_irred_generic.py

The theorem-level argument written to /tmp/l21_irred_generic.md uses an exact
irreducible specialization to prove irreducibility over Q(a,Z), then repairs
the fixed-fiber issue in a two-variable HIT argument by certifying the a=1
fiber separately at every one of the 353 target cells.
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


OUT = HERE / "data" / "l21_irred_generic.jsonl"
REPORT = Path("/tmp/l21_irred_generic.md")
IRRED_PRIME_LIMIT = 1000
PACE_EVERY = 100
PACE_SECONDS = 0.002
GRID_SAMPLE_SIZE = 32
OFFGRID_PRIMES = (
    101,
    103,
    107,
    109,
    113,
    127,
    131,
    137,
    139,
    149,
    151,
    157,
    163,
    167,
    173,
    179,
)


class Pacer:
    """Sleep-pace the finite-field certificate loops."""

    def __init__(self) -> None:
        self.steps = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.steps += 1
        if self.steps % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()
PRIMES = tuple(h10q.primerange(2, IRRED_PRIME_LIMIT + 1))
CERT_CACHE: dict[tuple[F, F, F], tuple[int | None, tuple[int, ...]]] = {}


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def quadratic_character(value: int, prime: int) -> int:
    """Ordinary Legendre character, including the zero value."""
    residue = value % prime
    if residue == 0:
        return 0
    return 1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1


def primitive_integer_poly(coefficients: list[F]) -> list[int]:
    """Primitive ascending integer coefficients with positive leading term."""
    denominator = 1
    for coefficient in coefficients:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in coefficients]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, coefficient)
    assert content > 0
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return integers


def polynomial_degree(coefficients: list[F]) -> int:
    degree = len(coefficients) - 1
    while degree >= 0 and coefficients[degree] == 0:
        degree -= 1
    return degree


def poly_mul(left: list[F], right: list[F]) -> list[F]:
    product = [F(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if x:
            for j, y in enumerate(right):
                if y:
                    product[i + j] += x * y
    return product


def branch_parameters(a_value: int | F, z_value: int | F, tau: F) -> dict[str, Any]:
    a = F(a_value)
    z = F(z_value)
    A = 1 + 4 * a * a
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    delta = 1 - A * tau * tau
    s = (a - 1) / 2
    P = h10q._l10_P(a, Z, D, A, delta, s)
    return {
        "a": a,
        "z": z,
        "A": A,
        "Z": Z,
        "D": D,
        "tau": tau,
        "delta": delta,
        "s": s,
        "P": P,
    }


def constructed_parameters(a_value: int | F, z_value: int | F) -> dict[str, Any]:
    a = F(a_value)
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    data = branch_parameters(a, z_value, tau)
    assert data["delta"] == -4 * a**4 / A
    assert polynomial_degree(data["P"]) == 8
    assert data["P"][0] == data["P"][8] == 4 * a**8 * A * data["Z"] ** 4
    return data


def irreducibility_certificate(
    a_value: int | F,
    z_value: int | F,
    tau: F,
) -> tuple[int | None, tuple[int, ...]]:
    """Return an exact irreducible-reduction prime, if one occurs <= limit."""
    key = (F(a_value), F(z_value), F(tau))
    cached = CERT_CACHE.get(key)
    if cached is not None:
        return cached
    data = branch_parameters(*key)
    primitive = tuple(primitive_integer_poly(data["P"]))
    certificate = None
    if len(primitive) == 9:
        for prime in PRIMES:
            PACER.tick()
            if h10q._l13_irred8(list(primitive), prime):
                certificate = prime
                break
    answer = (certificate, primitive)
    CERT_CACHE[key] = answer
    return answer


def rational_sqrt(value: F) -> F | None:
    if value < 0:
        return None
    numerator_root = math.isqrt(value.numerator)
    denominator_root = math.isqrt(value.denominator)
    if (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    ):
        return F(numerator_root, denominator_root)
    return None


def square_delta_factors(data: dict[str, Any], sigma: F) -> tuple[list[F], list[F]]:
    """The lead L21 square-delta factorization, represented ascending."""
    a, A, Z, D = data["a"], data["A"], data["Z"], data["D"]
    q = [F(1)]
    for _ in range(4):
        q = poly_mul(q, [F(-1), F(1)])
    Ng = [F(0)] * 5
    for index, coefficient in enumerate(q):
        Ng[index] -= A * coefficient
    Ng[2] += 16 * a**4
    center = [F(0), F(0), 4 * D * A]
    twist = [sigma * a * a * Z * Z * coefficient for coefficient in Ng]
    width = max(len(center), len(twist))
    center += [F(0)] * (width - len(center))
    twist += [F(0)] * (width - len(twist))
    minus = [center[i] - twist[i] for i in range(width)]
    plus = [center[i] + twist[i] for i in range(width)]
    assert poly_mul(minus, plus) == data["P"]
    assert polynomial_degree(minus) == polynomial_degree(plus) == 4
    return minus, plus


def generic_specialization() -> dict[str, Any]:
    """Exact (a,Z)=(1,27) certificate used in the generic proof."""
    data = constructed_parameters(1, 3)
    certificate, primitive = irreducibility_certificate(1, 3, data["tau"])
    expected = (
        13286025,
        -106288200,
        286978140,
        -403895160,
        840899786,
        -403895160,
        286978140,
        -106288200,
        13286025,
    )
    assert primitive == expected
    assert certificate == 17
    reduction = [coefficient % 17 for coefficient in primitive]
    x_p4 = h10q._l13_xpk(reduction, 17, 4)
    x_p8 = h10q._l13_xpk(reduction, 17, 8)
    gcd_certificate = h10q._l13_pgcd(
        [
            (coefficient - x_coefficient) % 17
            for coefficient, x_coefficient in zip(
                x_p4 + [0] * 2,
                [0, 1] + [0] * len(x_p4),
            )
        ],
        reduction,
        17,
    )
    assert reduction == [15, 16, 1, 3, 5, 3, 1, 16, 15]
    assert x_p4 == [8, 9, 10, 11, 10, 9, 8, 16]
    assert x_p8 == [0, 1]
    assert gcd_certificate == [11]
    return {
        "type": "generic-specialization",
        "label": "PROVED",
        "branch": "tau_dagger=(1+2*a^2)/(1+4*a^2)",
        "a": "1",
        "z": "3",
        "Z": "27",
        "A": "5",
        "delta": "-4/5",
        "primitive_P_coefficients_ascending": list(primitive),
        "prime": 17,
        "reduction_coefficients_ascending": reduction,
        "x_pow_17_pow_4_mod_P_ascending": x_p4,
        "x_pow_17_pow_8_mod_P_ascending": x_p8,
        "gcd_x_pow_17_pow_4_minus_x_with_P": gcd_certificate,
        "finite_field_criterion": "x^(17^8)=x mod P and gcd(x^(17^4)-x,P)=1",
        "generic_verdict": "P is irreducible of degree 8 in Q(a,Z)[b] on tau_dagger",
        "integral_model": "H=(A/4)P=a^8*Z^4*N_g^2+4*A^3*D^2*b^4-2*A^4*(a-1)^2*D^2*b^5",
        "leading_coefficient_H": "a^8*A^2*Z^4",
        "newton_places": {
            "Z=0": {
                "coefficient_valuations_degrees_0_to_8": [4, 4, 4, 4, 0, 0, 4, 4, 4],
                "lower_slopes": ["-1", "0", "4/3"],
            },
            "a=0": {
                "coefficient_valuations_degrees_0_to_8": [8, 8, 8, 8, 0, 0, 8, 8, 8],
                "lower_slopes": ["-2", "0", "8/3"],
            },
            "A=0": {
                "coefficient_valuations_degrees_0_to_8": [2, 2, 1, 1, 0, 1, 1, 2, 2],
                "lower_slopes": ["-1/2", "1/2"],
            },
            "Z=infinity": {
                "coefficient_valuations_degrees_0_to_8": [-4] * 9,
                "lower_slopes": ["0"],
            },
        },
    }


def grid_fiber_certificates() -> list[dict[str, Any]]:
    """Certify P(1,z) at every target z, proving each Q(a)-fiber generic."""
    rows: list[dict[str, Any]] = []
    for w, unit in sorted(h10q._l9_grid()):
        z = F(w) * F(*unit)
        data = constructed_parameters(1, z)
        certificate, primitive = irreducibility_certificate(1, z, data["tau"])
        assert certificate is not None
        assert h10q._l13_irred8(list(primitive), certificate)
        rows.append(
            {
                "type": "grid-fiber-certificate",
                "label": "PROVED",
                "cell": [w, list(unit)],
                "z": frac_text(z),
                "Z": frac_text(z**3),
                "specialization_a": 1,
                "specialization_need_not_be_character_admissible": True,
                "degree": 8,
                "irreducibility_certificate_prime": certificate,
                "vertical_conclusion": "P(a,z) irreducible in Q(a)[b]",
            }
        )
    assert len(rows) == 353
    return rows


def evenly_spaced(items: list[Any], count: int) -> list[Any]:
    assert 1 < count <= len(items)
    indices = [round(index * (len(items) - 1) / (count - 1)) for index in range(count)]
    assert len(set(indices)) == count
    return [items[index] for index in indices]


def sampled_cells() -> list[tuple[str, int, tuple[int, int]]]:
    grid = sorted(h10q._l9_grid())
    cells = [("grid", w, unit) for w, unit in evenly_spaced(grid, GRID_SAMPLE_SIZE)]
    units = list(h10q._L9_U_POOL)
    for index, w in enumerate(OFFGRID_PRIMES):
        assert h10q._is_prime(w) and w >= 100
        unit = units[(5 * index + 2) % len(units)]
        cells.append(("off-grid", w, (unit.numerator, unit.denominator)))
    return cells


def first_irreducible_admissible_rows() -> list[dict[str, Any]]:
    """Enumerate character-admissible positive odd a in increasing order."""
    rows: list[dict[str, Any]] = []
    for scope, w, unit in sampled_cells():
        z = F(w) * F(*unit)
        attempted: list[dict[str, Any]] = []
        for a in range(1, 20001, 2):
            A = 1 + 4 * a * a
            if quadratic_character(A, w) != -1:
                continue
            data = constructed_parameters(a, z)
            if data["D"] == 0:
                continue
            certificate, primitive = irreducibility_certificate(a, z, data["tau"])
            attempted.append(
                {
                    "index": len(attempted) + 1,
                    "a": a,
                    "A": A,
                    "certificate_prime": certificate,
                }
            )
            if certificate is not None:
                assert h10q._l13_irred8(list(primitive), certificate)
                rows.append(
                    {
                        "type": "first-irreducible-admissible-a",
                        "label": "PROVED_PER_ROW",
                        "scope": scope,
                        "cell": [w, list(unit)],
                        "z": frac_text(z),
                        "admissibility": "a positive odd, D!=0, and (1+4*a^2|w)=-1",
                        "first_irreducible_index": len(attempted),
                        "first_irreducible_a": a,
                        "A": A,
                        "irreducibility_certificate_prime": certificate,
                        "attempted_admissible_values": attempted,
                    }
                )
                break
            if len(attempted) >= 64:
                raise AssertionError((scope, w, unit, "64 admissible a without certificate"))
        else:
            raise AssertionError((scope, w, unit, "no admissible a below 20000"))
    assert len(rows) == GRID_SAMPLE_SIZE + len(OFFGRID_PRIMES)
    return rows


def hunt_tau_values(a: F) -> list[tuple[str, F]]:
    A = 1 + 4 * a * a
    return [
        ("tau_dagger", (1 + 2 * a * a) / A),
        ("tau_alt=2a/A", 2 * a / A),
        ("tau=0", F(0)),
        ("tau=1/3", F(1, 3)),
        ("tau=1/2", F(1, 2)),
        ("tau=-1/2", F(-1, 2)),
        ("tau=1", F(1)),
    ]


def reducible_hunt() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Scan a rational (a,z,tau) box, proving every decided row exactly."""
    a_values = tuple(F(a) for a in range(1, 32, 2))
    z_values = tuple(
        sorted(
            {
                F(numerator, denominator)
                for denominator in range(1, 5)
                for numerator in range(-12, 13)
                if numerator
            }
        )
    )
    assert len(a_values) == 16 and len(z_values) == 64
    rows: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()
    branch_counts: dict[str, Counter[str]] = {}
    for a in a_values:
        for z in z_values:
            for tau_label, tau in hunt_tau_values(a):
                data = branch_parameters(a, z, tau)
                degree = polynomial_degree(data["P"])
                sigma = rational_sqrt(data["delta"])
                row: dict[str, Any] = {
                    "type": "reducible-hunt",
                    "label": "EVIDENCE_BOX_WITH_EXACT_ROW_CERTIFICATES",
                    "a": frac_text(a),
                    "z": frac_text(z),
                    "tau_label": tau_label,
                    "tau": frac_text(tau),
                    "delta": frac_text(data["delta"]),
                    "degree": degree,
                }
                if degree != 8:
                    status = "DEGREE_DROP"
                    row["status"] = status
                elif a == 1 and sigma is not None:
                    minus, plus = square_delta_factors(data, sigma)
                    status = "PROVED_REDUCIBLE_EXPLICIT_4x4"
                    row.update(
                        {
                            "status": status,
                            "sigma": frac_text(sigma),
                            "factor_minus_coefficients_ascending": [frac_text(x) for x in minus],
                            "factor_plus_coefficients_ascending": [frac_text(x) for x in plus],
                            "factor_identity": "P=(4*D*A*b^2-sigma*a^2*Z^2*N_g)*(4*D*A*b^2+sigma*a^2*Z^2*N_g)",
                        }
                    )
                else:
                    certificate, primitive = irreducibility_certificate(a, z, tau)
                    if certificate is None:
                        status = "OPEN_NO_IRREDUCIBLE_REDUCTION_AT_P_LE_1000"
                        row["status"] = status
                    else:
                        assert h10q._l13_irred8(list(primitive), certificate)
                        status = "PROVED_IRREDUCIBLE_MOD_P"
                        row.update(
                            {
                                "status": status,
                                "irreducibility_certificate_prime": certificate,
                            }
                        )
                status_counts[status] += 1
                branch_counts.setdefault(tau_label, Counter())[status] += 1
                rows.append(row)
    expected_rows = len(a_values) * len(z_values) * 7
    assert len(rows) == expected_rows == 7168
    dagger = [row for row in rows if row["tau_label"] == "tau_dagger"]
    assert len(dagger) == 1024
    assert all(row["status"] == "PROVED_IRREDUCIBLE_MOD_P" for row in dagger)
    summary = {
        "type": "hunt-summary",
        "label": "EVIDENCE_BOX_WITH_EXACT_ROW_CERTIFICATES",
        "box": {
            "a": "positive odd integers 1<=a<=31",
            "z": "all distinct reduced n/d from -12<=n<=12, n!=0, 1<=d<=4",
            "tau": [label for label, _tau in hunt_tau_values(F(1))],
            "a_values": len(a_values),
            "z_values": len(z_values),
            "rows": len(rows),
        },
        "status_counts": dict(sorted(status_counts.items())),
        "branch_status_counts": {
            branch: dict(sorted(counts.items()))
            for branch, counts in branch_counts.items()
        },
        "constructed_branch": {
            "rows": len(dagger),
            "proved_irreducible_mod_p": sum(
                row["status"] == "PROVED_IRREDUCIBLE_MOD_P" for row in dagger
            ),
            "reducible_instances": sum(
                row["status"] == "PROVED_REDUCIBLE_EXPLICIT_4x4" for row in dagger
            ),
            "unresolved": sum(row["status"].startswith("OPEN") for row in dagger),
        },
        "off_branch_reducible_locus_replayed": "a=1 and delta=sigma^2; tau=0 (sigma=1) and tau=1/3 (sigma=2/3) in this box",
    }
    return rows, summary


def write_report(
    generic: dict[str, Any],
    fiber_rows: list[dict[str, Any]],
    first_rows: list[dict[str, Any]],
    hunt_summary: dict[str, Any],
    summary: dict[str, Any],
) -> None:
    first_distribution = Counter(row["first_irreducible_index"] for row in first_rows)
    fiber_primes = Counter(row["irreducibility_certificate_prime"] for row in fiber_rows)
    hunt = hunt_summary
    lines = [
        "# L21 generic irreducibility and Hilbert specialization",
        "",
        "## Verdict",
        "",
        "**PROVED on the constructed branch** `tau_dagger=(1+2a^2)/(1+4a^2)`: the octic `P` is irreducible in `Q(a,Z)[b]`.  For each of the 353 frozen target cells, `P(a,z)` is also irreducible in `Q(a)[b]`; quantitative Hilbert irreducibility then proves that the infinite character-admissible progression contains an `a` with irreducible specialized `P`.  Thus the generic-plus-specialization route supplies SOME admissible irreducible `a` in every target cell.",
        "",
        "The two-variable generic theorem by itself would not justify that last sentence: a thin subset of `A^2(Q)` may contain an entire vertical line `Z=Z_0`.  The 353 exact vertical-fiber certificates below are the necessary repair.  Beyond the frozen grid, the generic theorem and the finite off-grid table are not mislabeled as a theorem for every rational `Z`.",
        "",
        "## 1. Integral model and generic proof",
        "",
        "Put `A=1+4a^2`, `D=1-Z-a^2Z^2`, `s=(a-1)/2`, and use the constructed identity `delta=-4a^4/A`.  Multiplication by the nonzero element `A/4` does not change reducibility over `K=Q(a,Z)`, and gives",
        "",
        "    H=(A/4)P",
        "     =a^8 Z^4 N_g(b)^2 + 4 A^3 D^2 b^4",
        "       - 2 A^4 (a-1)^2 D^2 b^5  in Q[a,Z][b].",
        "",
        "Its degree and constant coefficients are both `L=a^8 A^2 Z^4`; in particular `deg_b H=8`.  Let `R=Q[a,Z]` and `S=R[1/L]`.  The ring `S` is a localization of a UFD, hence integrally closed, and `H/L` is monic in `S[b]`.",
        "",
        "If `H` factored in `K[b]`, take both factors monic after dividing by `L`.  A monic divisor of a monic polynomial over an integrally closed domain belongs to that domain (its coefficients are integral over `S` and lie in `Frac(S)`).  Thus the factors lie in `S[b]`.  Any specialization with `L!=0` preserves their positive degrees.",
        "",
        "Specialize `(a,Z)=(1,27)` (that is, `z=3`).  The primitive ascending coefficient vector is",
        "",
        f"    {generic['primitive_P_coefficients_ascending']}",
        "",
        "Modulo 17 it is",
        "",
        f"    q={generic['reduction_coefficients_ascending']}.",
        "",
        f"Exact Frobenius arithmetic gives `x^(17^8) mod q={generic['x_pow_17_pow_8_mod_P_ascending']}`, `x^(17^4) mod q={generic['x_pow_17_pow_4_mod_P_ascending']}`, and `gcd(x^(17^4)-x,q)={generic['gcd_x_pow_17_pow_4_minus_x_with_P']}`.  For degree 8, the only prime divisor of 8 is 2, so these are exactly the finite-field irreducibility conditions.  Hence `q` is irreducible over `F_17`, the specialization is irreducible over `Q`, and the supposed generic factorization is impossible.",
        "",
        "**Generic irreducibility verdict: PROVED, degree 8, over `Q(a,Z)` on `tau_dagger`.**",
        "",
        "### Cheap Newton polygons checked first",
        "",
        "For the integral model `H`, coefficient valuations in degrees 0 through 8 and lower slopes are:",
        "",
        "| place | valuations | lower slopes |",
        "|:--|:--|:--|",
    ]
    for place, record in generic["newton_places"].items():
        lines.append(
            f"| `{place}` | `{record['coefficient_valuations_degrees_0_to_8']}` | "
            f"`{record['lower_slopes']}` |"
        )
    lines.extend(
        [
            "",
            "None is a pure Eisenstein polygon or a single segment with slope denominator 8: the finite places have multiple segments with denominators at most 3, and `Z=infinity` is horizontal.  They therefore do not prove octic irreducibility; the exact specialization above is decisive.",
            "",
            "## 2. Fixed fibers and the actual HIT hypotheses",
            "",
            f"For every frozen cell `z=w*u`, the script specialized only `a=1`, primitive-normalized `P(1,z)`, and found an irreducible degree-8 reduction modulo a prime at most {max(fiber_primes)}.  This succeeded on **{len(fiber_rows)}/{len(fiber_rows)}** cells.  The witness `a=1` need not satisfy the cell's character condition: it is used only in the same localization lemma with `R_z=Q[a]` and `L_z=a^8(1+4a^2)^2Z^4`.  Consequently `P(a,z)` is irreducible in `Q(a)[b]` separately for every one of the 353 vertical fibers.",
            "",
            f"Fiber certificate-prime distribution: `{dict(sorted(fiber_primes.items()))}`.",
            "",
            "Characteristic zero makes each irreducible octic separable.  Hilbert's irreducibility theorem therefore applies to each `P(a,z)` and says that the rational `a` where the specialization drops degree or becomes reducible form a thin subset `E_z` of `A^1(Q)`.  In the integral quantitative form of Cohen--Serre (Serre, *Lectures on the Mordell--Weil Theorem*, section 13, Theorems 1--2),",
            "",
            "    #{n in E_z intersect Z : |n|<=B} = O_z(B^(1/2) log B).",
            "",
            "This is the theorem used; no density assertion is being inferred merely from the word `thin`.",
            "",
            "For a cell prime `w`, the exact character sum `sum_{r mod w}(1+4r^2|w)=-1` supplies a residue `r` with character `-1`.  CRT with parity supplies an odd `a_0 mod 2w` above `r`.  Every positive member",
            "",
            "    a=a_0+2w*k,  k>=0,",
            "",
            "has the required character.  The class-existence theorem then solves the `q1` side conditions after `a` is chosen; they do not remove this `a`-progression.  The condition `D=1-Z-a^2Z^2!=0` removes at most two values.  Up to `B`, the progression has `B/(2w)+O(1)` members, whereas the thin exceptional set has `O_z(sqrt(B) log B)`.  Since the former dominates the latter, the progression cannot be contained in `E_z`; for all sufficiently large `B` it contains a good `a`.  This is the full thin-set-dodging argument.",
            "",
            "## 3. First irreducible admissible `a`: mechanical table",
            "",
            f"The table has {len(first_rows)} cells: {sum(row['scope']=='grid' for row in first_rows)} spread across the frozen grid and {sum(row['scope']=='off-grid' for row in first_rows)} with `w>100`.  Indexing is among positive odd values satisfying `(1+4a^2|w)=-1` and `D!=0`, in increasing order.",
            "",
            "| scope | cell `(w,u)` | z | first index | a | A | certificate p |",
            "|:--|:--|--:|--:|--:|--:|--:|",
        ]
    )
    for row in first_rows:
        lines.append(
            f"| {row['scope']} | `{tuple([row['cell'][0], tuple(row['cell'][1])])}` | "
            f"{row['z']} | {row['first_irreducible_index']} | {row['first_irreducible_a']} | "
            f"{row['A']} | {row['irreducibility_certificate_prime']} |"
        )
    lines.extend(
        [
            "",
            f"First-index distribution: **`{dict(sorted(first_distribution.items()))}`**.  Every row is theorem-level per row because its displayed prime is replayed by the kernel's exact Frobenius criterion; the off-grid collection is finite evidence about generality, not the all-`Z` theorem.",
            "",
            "## 4. Reducible-instance hunt",
            "",
            f"The active box contains **{hunt['box']['rows']}** octics: {hunt['box']['a']}; {hunt['box']['z']}; and seven tau choices `{hunt['box']['tau']}`.  Status totals are `{hunt['status_counts']}`.",
            "",
            f"On the constructed branch `tau_dagger`, **{hunt['constructed_branch']['proved_irreducible_mod_p']}/{hunt['constructed_branch']['rows']}** rows have exact irreducible reductions, with **{hunt['constructed_branch']['reducible_instances']} reducible and {hunt['constructed_branch']['unresolved']} unresolved**.  Thus the hunt found no reducible constructed-branch instance.",
            "",
            "The hunt did find a genuine off-branch reducible family, replaying the lead L21 square-delta lemma: when `a=1` and `delta=sigma^2`,",
            "",
            "    P=(4DAb^2-sigma*a^2*Z^2*N_g)(4DAb^2+sigma*a^2*Z^2*N_g).",
            "",
            "Both factors have degree 4 for the recorded nonzero `z`.  In this box these are the `tau=0` (`sigma=1`) and `tau=1/3` (`sigma=2/3`) rows.  This refutes any blanket all-`tau` irreducibility claim, but it cannot occur on `tau_dagger`, where `delta=-4a^4/A<0` for nonzero rational `a`.",
            "",
            "Branch-by-branch hunt counts:",
            "",
            "| tau branch | exact status counts |",
            "|:--|:--|",
        ]
    )
    for branch, counts in hunt["branch_status_counts"].items():
        lines.append(f"| `{branch}` | `{counts}` |")
    lines.extend(
        [
            "",
            "## 5. Replay",
            "",
            "From the repository root:",
            "",
            "    nice -n 19 python3 math/h10q/l21_irred_generic.py",
            "",
            f"Wall-clock: **{summary['wall_seconds']:.3f} s**.  Pacing: {summary['pacing']['steps']} finite-field loop steps and {summary['pacing']['sleeps']} sleeps of {PACE_SECONDS} s.  Refusals/unresolved hunt rows: **{summary['hunt']['unresolved']}**; none were used as evidence.",
            "",
            "Artifacts: `math/h10q/l21_irred_generic.py`, `math/h10q/data/l21_irred_generic.jsonl`, `/tmp/l21_irred_generic.md`.",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    generic = generic_specialization()
    fiber_rows = grid_fiber_certificates()
    first_rows = first_irreducible_admissible_rows()
    hunt_rows, hunt_summary = reducible_hunt()

    first_distribution = Counter(row["first_irreducible_index"] for row in first_rows)
    unresolved = sum(
        row["status"].startswith("OPEN") or row["status"] == "DEGREE_DROP"
        for row in hunt_rows
    )
    elapsed = time.perf_counter() - started
    summary = {
        "type": "summary",
        "label": "PROVED",
        "scope": "constructed tau_dagger branch; all 353 frozen target cells",
        "generic_irreducibility_Q_a_Z": "PROVED",
        "generic_degree": 8,
        "generic_specialization": {"a": 1, "Z": 27, "prime": 17},
        "vertical_fibers": {
            "proved_irreducible_over_Q_a": len(fiber_rows),
            "target_cells": 353,
            "method": "a=1 exact irreducible reduction separately at each z",
        },
        "HIT": {
            "hypotheses_verified": [
                "P(a,z) irreducible in Q(a)[b] for each frozen z",
                "degree 8",
                "separable in characteristic zero",
                "leading coefficient nonzero on positive odd a and z!=0",
            ],
            "exceptional_set": "thin in A^1(Q) for each fixed z",
            "integer_count": "O_z(B^(1/2) log B)",
            "admissible_progression_count": "B/(2w)+O(1)",
            "dodge": "PROVED for every frozen target cell",
        },
        "mechanical_first_a": {
            "rows": len(first_rows),
            "grid_rows": sum(row["scope"] == "grid" for row in first_rows),
            "off_grid_rows": sum(row["scope"] == "off-grid" for row in first_rows),
            "first_index_distribution": dict(sorted(first_distribution.items())),
        },
        "hunt": {
            "rows": len(hunt_rows),
            "status_counts": hunt_summary["status_counts"],
            "constructed_branch": hunt_summary["constructed_branch"],
            "unresolved": unresolved,
        },
        "off_grid_uniform_claim": "OPEN; finite rows only",
        "wall_seconds": elapsed,
        "pacing": {"steps": PACER.steps, "sleeps": PACER.sleeps},
        "refusals": unresolved,
    }
    rows: list[dict[str, Any]] = [
        {
            "type": "meta",
            "artifact": "l21_irred_generic",
            "version": 1,
            "labels": {
                "generic": "PROVED",
                "frozen_grid_HIT_dodge": "PROVED",
                "finite_hunt": "EVIDENCE_WITH_EXACT_ROW_CERTIFICATES",
                "off_grid_uniformity": "OPEN",
            },
            "irred_engine": "h10q._l13_irred8 exact Frobenius criterion",
            "prime_limit": IRRED_PRIME_LIMIT,
        },
        generic,
        *fiber_rows,
        *first_rows,
        *hunt_rows,
        hunt_summary,
        summary,
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(generic, fiber_rows, first_rows, hunt_summary, summary)
    print(
        "L21 generic irreducibility: "
        f"fibers={len(fiber_rows)}, first-a={len(first_rows)}, hunt={len(hunt_rows)}, "
        f"unresolved={unresolved}, wall={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
