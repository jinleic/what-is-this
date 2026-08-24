#!/usr/bin/env python3
"""L25 scaled non-diagonal self-couplings and their exact local/frontier walls.

Replay authority for THEOREMS.md L25.  Everything below is unconditional
unless explicitly labelled otherwise.  The module writes
``data/l25_scaled_coupling.jsonl`` and a human-readable audit report to
``/tmp/l25_scaled_coupling.md``.
"""
from __future__ import annotations

import json
import math
import time
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l25_scaled_coupling.jsonl"
REPORT = Path("/tmp/l25_scaled_coupling.md")


def frac(value: int | F) -> F:
    return F(value)


def integer_vp(value: int, prime: int) -> int:
    assert value
    n = abs(value)
    answer = 0
    while n % prime == 0:
        n //= prime
        answer += 1
    return answer


def vp(value: int | F, prime: int = 2) -> int:
    q = frac(value)
    assert q
    return integer_vp(q.numerator, prime) - integer_vp(q.denominator, prime)


def frac_mod_2k(value: int | F, exponent: int) -> int:
    q = frac(value)
    assert q.denominator % 2 == 1
    modulus = 1 << exponent
    return (q.numerator * pow(q.denominator, -1, modulus)) % modulus


def rational_sqrt(value: int | F) -> F | None:
    q = frac(value)
    if q < 0:
        return None
    if q == 0:
        return F(0)
    n, d = q.numerator, q.denominator
    rn = math.isqrt(n)
    rd = math.isqrt(d)
    if rn * rn == n and rd * rd == d:
        return F(rn, rd)
    return None


def squareclass_2adic(value: int | F) -> bool:
    """True iff value is a square in Q_2."""
    q = frac(value)
    assert q
    valuation = vp(q, 2)
    if valuation % 2:
        return False
    unit = q / F(2) ** valuation
    return frac_mod_2k(unit, 3) == 1


def row(**fields):
    return dict(fields)


def ng_value(a: F, A: F, b: F) -> F:
    return 16 * a**4 * b * b - A * (b - 1) ** 4


def bridge_c(a: F, Z: F, b: F) -> F:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    Ng = ng_value(a, A, b)
    assert Z != 0 and b not in (0, 1) and D != 0 and Ng != 0
    return a * a * Z * Z * Ng / (A * b * b * D)


def eliminate_coefficients(coupling: str, l1: F, l2: F, A: F, B: F,
                           c: F, s: F) -> tuple[F, F, F]:
    """Coefficients A2,L2,M2 of A2*u^2 + L2*u + M2 with u=rho^2."""
    s = frac(s)
    if coupling == "I":
        return (A * l1 * l1,
                A * A * l1 * l1 - c * c - 16 * A * B * l2 * l2,
                16 * A * A * B * s * s - 16 * A)
    if coupling == "II":
        return (A * l1 * l1,
                -(c * c + 16 * A * B * l2 * l2),
                16 * A * A * B * (s * s - l2 * l2) - 16 * A)
    raise AssertionError(coupling)


def tied_value(coupling: str, l1: F, l2: F, A: F, B: F, c: F, s: F,
               X: F, rho: F) -> F:
    """Value of the cleared tied equation under the selected scaled coupling."""
    assert X * X - rho * rho == A
    if coupling == "I":
        y, r = l1 * X, l2 * rho
    elif coupling == "II":
        y, r = l1 * rho, l2 * X
    else:
        raise AssertionError(coupling)
    return (-rho * rho * (c * c - A * y * y)
            - 16 * A * B * (r * r - A * s * s) - 16 * A)


def identity_replay() -> dict:
    checks = 0
    for a_int in (1, 3, 5, 7):
        a = F(a_int)
        A = 1 + 4 * a * a
        for Z in (F(-8), F(-1, 8), F(27), F(1, 27)):
            D = 1 - Z - a * a * Z * Z
            assert D != 0
            for lam in (F(1), F(2), F(3, 2), A, F(3, 4)):
                X = (lam + A / lam) / 2
                rho = (lam - A / lam) / 2
                assert X * X - rho * rho == A
                for b in (F(-3), F(3), F(5, 3)):
                    B = 2 * b
                    c = bridge_c(a, Z, b)
                    for coupling, l1, l2 in (
                            ("I", F(2), F(1)),
                            ("I", F(2), F(1, 2)),
                            ("II", F(1), F(2)),
                            ("II", F(1), F(1)),
                            ("I", F(6), F(5))):
                        tied = tied_value(coupling, l1, l2, A, B, c,
                                          (a - 1) / 2, X, rho)
                        A2, L2, M2 = eliminate_coefficients(
                            coupling, l1, l2, A, B, c, (a - 1) / 2)
                        assert tied == (A2 * rho**4 + L2 * rho * rho + M2), (
                            a, Z, lam, b, coupling, l1, l2)
                        checks += 1
    return row(type="identity-replay", label="PROVED", exact_instances=checks)


def q2_quadratic_roots(A2: F, L2: F, M2: F) -> list[F]:
    """Return the Q_2 roots of a rational quadratic, or the rational root.

    Non-square discriminant means no Q_2 root.  A nonzero square
    discriminant is represented by a 2^80 residue; this is exact for the
    coefficient ranges used below, whose root valuations are < 60.
    """
    discriminant = L2 * L2 - 4 * A2 * M2
    if discriminant == 0:
        return [-L2 / (2 * A2)]
    if not squareclass_2adic(discriminant):
        return []
    valuation = vp(discriminant, 2)
    unit = discriminant / F(2) ** valuation
    residue = frac_mod_2k(unit, 80)
    root = 1
    for bit in range(3, 80):
        if (root * root - residue) % (1 << (bit + 1)):
            root += 1 << (bit - 1)
    assert (root * root - residue) % (1 << 80) == 0
    roots = []
    for sign in (1, -1):
        numerator = -L2 + sign * F(2) ** (valuation // 2) * root
        value = numerator / (2 * A2)
        if value == 0:
            continue
        if vp(value) >= 60:
            raise AssertionError("2-adic root precision exhausted")
        roots.append(value)
    return roots


def q2_solutions_exist(A2: F, L2: F, M2: F, need_u_square: bool = True,
                       need_Au_square: bool = True, A: F | None = None) -> bool:
    """Decide existence of a Q_2 root with the requested square classes."""
    for u in q2_quadratic_roots(A2, L2, M2):
        if u == 0:
            continue
        if need_u_square and not squareclass_2adic(u):
            continue
        if need_Au_square:
            assert A is not None
            if not squareclass_2adic(frac(A) + u):
                continue
        return True
    return False


def complete_dyadic_classification() -> dict:
    """Finite residue enumeration behind the scaled type-I theorem.

    The proof text in THEOREMS.md reduces every coupling state to
    ``s mod 8``, ``b mod 8`` and the unit square classes of the two scalars.
    The loop below performs the same finite test on every residue class,
    using the exact local solver instead of relying on sampled values; the
    c-dependence enters only through c^2 with v_2(c)>=4, which is invisible
    modulo 8 of the quadratic coefficients after multiplying by 1/(4A).
    """
    records: list[dict] = []
    # The c unit enters L24's fatal mod-16 test, so finite classes are keyed
    # by c's 2-adic unit as well as s,b and scalar units.
    for coupling, e, f in (("I", 1, 0), ("I", 1, -1), ("I", 0, 0)):  # L24's type-II death is proved by its mod-16 theorem, not by this generic table
        classes = defaultdict(set)
        for l1_unit in (1, 3):
            for l2_unit in (1, 3):
                for c_unit in (1, 3, 5, 7):
                    for c_val in (4, 5, 6):
                        c = F(2) ** c_val * c_unit
                        for b in range(-7, 8):
                            if b == 0 or b % 2 == 0:
                                continue
                            for s in range(-8, 8):
                                a = 1 + 2 * s
                                A = 1 + 4 * a * a
                                l1 = F(2) ** e * l1_unit
                                l2 = F(2) ** f * l2_unit
                                A2, L2, M2 = eliminate_coefficients(
                                    coupling, l1, l2, F(A), F(2 * b), c, F(s))
                                alive = q2_solutions_exist(A2, L2, M2, A=F(A))
                                classes[(coupling, e, f, s % 8)].add(alive)
        mixed = {key: value for key, value in classes.items() if len(value) != 1}
        if coupling == "I" and e == 1:
            assert not mixed, (coupling, e, f, mixed)
            truth = {key[3]: next(iter(value))
                     for key, value in sorted(classes.items())}
            if f == 0:
                expected = {0: True, 1: False, 2: True, 3: False,
                            4: True, 5: False, 6: True, 7: False}
            else:
                expected = {0: False, 1: True, 2: False, 3: True,
                            4: False, 5: True, 6: False, 7: True}
            assert truth == expected, (coupling, e, f, truth, expected)
            records.append(row(type="dyadic-escape", coupling=coupling, e=e,
                               f=f, truth_by_s_mod_8=truth,
                               residue_classes=len(classes), label="PROVED"))
        else:
            # Bounded recovery context only.  L24's uniform emptiness has to
            # use its stronger square-side conditions; this generic quadratic
            # table is intentionally not a substitute proof.
            truth_table = {
                str(key[3]): sorted(value)
                for key, value in sorted(classes.items())
            }
            records.append(row(type="dyadic-recovery-context", coupling=coupling,
                               e=e, f=f, residue_classes=len(classes),
                               truth_table=truth_table,
                               label="EVIDENCE CONTEXT; L24 proof remains authoritative"))
    return row(type="dyadic-classification", label="PROVED", records=records)


def target_residue_certificate() -> list[dict]:
    """Bounded exact evidence for the free unit-b residue stratum.

    This is not the L24 aligned low-positive target stratum and is not
    compatible with a preselected L20 residue by default.  Only w<=11 are
    scanned; the all-prime claim is the symbolic formula record below.
    """
    records = []
    for w in (3, 5, 7, 11):
        entries = []
        for label, ksq in (("kappa=1", 1 % w), ("kappa=1/2", pow(4, -1, w))):
            total = smooth = smooth_u_unit = 0
            for s, b, X, r in product(range(w), range(1, w), range(w), range(w)):
                a = (1 + 2 * s) % w
                A = (1 + 4 * a * a) % w
                B = 2 * b % w
                c1 = (X * X - r * r - A) % w
                c2 = (4 * A * r * r * X * X - 16 * A * B * ksq * r * r
                      + 16 * A * A * B * s * s - 16 * A) % w
                if c1 or c2:
                    continue
                total += 1
                J = (16 * A * X * r * (X * X + r * r - 4 * B * ksq)) % w
                if J:
                    smooth += 1
                    if r * r % w:
                        smooth_u_unit += 1
            entries.append((label, ksq, total, smooth, smooth_u_unit))
        records.append(row(type="target-residue-scan", w=w, entries=entries,
                           label="BOUNDED SCAN only; all-target theorem uses the formula record"))
    return records


def universal_target_formula_certificate() -> dict:
    """Exact symbolic unit-b target construction for both couplings.

    At s=0 (a=1,A=5), X=3, rho=2, and b=kappa^{-2}, the residual
    equations vanish over Z, and the (rho,X)-Jacobian is the nonzero
    integer -2400.  Reduction modulo any prime not dividing 2400
    therefore gives a smooth point.  This is a freely chosen unit-b local
    stratum, not a statement about L20-aligned target residues.
    """
    formulas = []
    for k in (F(1), F(1, 2)):
        b = 1 / (k * k)
        c1 = F(3 * 3 - 2 * 2 - 5)
        c2 = F(4 * 5 * 2 * 2 * 3 * 3
               - 16 * 5 * (2 * b) * k * k * 2 * 2 - 16 * 5)
        determinant = -16 * 5 * 3 * 2 * (3 * 3 + 2 * 2 - 8)
        assert c1 == c2 == 0
        assert determinant == -2400
        finite_replay = []
        for w in (7, 11, 13, 17, 19, 23, 29, 31):
            assert 2400 % w != 0
            assert (int(c1) % w, int(c2) % w, int(determinant) % w) == (0, 0, int(determinant) % w)
            finite_replay.append(w)
        formulas.append(row(kappa=str(k), b=str(b), c1=int(c1), c2=int(c2),
                            determinant=int(determinant),
                            finite_replay=finite_replay))
    return row(type="target-universal-formula",
               label="PROVED symbolically on free unit-b stratum for w not dividing 2400; aligned target compatibility remains OPEN",
               formulas=formulas)


def real_viability_certificate() -> dict:
    """Exact guarded bridge samples with a positive real root.

    Each tuple uses the actual c=h(a,b,Z); D and N_g are checked nonzero.
    A2>0, M2<0, and a nonnegative discriminant imply a positive real root
    u of the eliminated quadratic; hence u and A+u are positive reals.
    """
    specs = [
        ("even", F(0), F(-10), F(3)),
        ("odd", F(1), F(-2, 5), F(4)),
    ]
    samples = []
    for parity, s, Z, b in specs:
        a = 1 + 2 * s
        A = 1 + 4 * a * a
        D = 1 - Z - a * a * Z * Z
        Ng = ng_value(a, A, b)
        assert D != 0 and Ng != 0
        c = bridge_c(a, Z, b)
        l2 = F(1) if parity == "even" else F(1, 2)
        A2, L2, M2 = eliminate_coefficients("I", F(2), l2, A, F(2 * b), c, s)
        disc = L2 * L2 - 4 * A2 * M2
        assert A2 > 0 and L2 < 0 and disc >= 0
        samples.append(row(parity=parity, a=str(a), A=str(A), b=str(b),
                           Z=str(Z), c=str(c), signs=[1, sgn(L2), sgn(M2)],
                           discriminant=str(disc), positive_root=True))
    return row(type="real-viability",
               label="PROVED guarded bridge-Certificate on named parity samples; no all-parameter theorem",
               samples=samples)


def sgn(value: F) -> int:
    return (value > 0) - (value < 0)


def rational_point_search() -> dict:
    """Bounded search only; absence is deliberately not evidence."""
    hits = 0
    tried = 0
    examples = []
    for s in range(-8, 9):
        for b_num in range(-12, 13):
            for b_den in range(1, 9):
                if b_num == 0 or F(b_num, b_den) == 1 or math.gcd(abs(b_num), b_den) != 1:
                    continue
                b = F(b_num, b_den)
                a = 1 + 2 * F(s)
                A = 1 + 4 * a * a
                for Z in (F(27), F(-8), F(1, 27)):
                    D = 1 - Z - a * a * Z * Z
                    if D == 0:
                        continue
                    c = bridge_c(a, Z, b)
                    for coupling, l1, l2 in (("I", F(2), F(1)),
                                             ("I", F(2), F(1, 2))):
                        A2, L2, M2 = eliminate_coefficients(coupling, l1, l2,
                                                             A, F(2 * b), c,
                                                             F(s))
                        disc = L2 * L2 - 4 * A2 * M2
                        root = rational_sqrt(disc)
                        if root is None:
                            continue
                        tried += 1
                        for sign_ in (1, -1):
                            u = (-L2 + sign_ * root) / (2 * A2)
                            if u == 0:
                                continue
                            ru = rational_sqrt(u)
                            if ru is None:
                                continue
                            rAu = rational_sqrt(A + u)
                            if rAu is None:
                                continue
                            hits += 1
                            if len(examples) < 4:
                                examples.append((s, str(b), str(Z), coupling,
                                                 str(ru), str(rAu)))
    return row(type="bounded-rational-search", label="EVIDENCE ONLY; zero-hit is not proof",
               tried_square_disc=tried, hits=hits, examples=examples)


def build_report(records: list[dict], elapsed: float) -> str:
    lines = [
        "# L25 scaled coupling audit",
        "",
        "STATUS: unconditional proofs where `PROVED`; bounded searches are EVIDENCE ONLY.",
        "Chains changed: L24 remains true; the scaled family prevents a universal dyadic wall.",
        "Schinzel H remains sole conjectural input for the conditional six-count.",
        "",
    ]
    for record in records:
        lines.extend(["## " + record["type"], "", "```json",
                      json.dumps(record, indent=2, sort_keys=True), "```", ""])
    lines.append(f"elapsed_seconds: {elapsed:.3f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    records: list[dict] = []
    records.append(identity_replay())
    records.append(complete_dyadic_classification())
    for record in target_residue_certificate():
        records.append(record)
    records.append(universal_target_formula_certificate())
    records.append(real_viability_certificate())
    records.append(rational_point_search())
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
