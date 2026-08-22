#!/usr/bin/env python3
"""L21: the reducible locus of the class polynomial P, and why the constructed
branch avoids it.

Two lemmas, both proved algebraically and replayed here exactly:

  L21a  If s = 0 (equivalently a = 1) and delta_tau = sigma^2 is a rational
        square, then
            P(b) = (4 D A b^2)^2 - (sigma a^2 Z^2 N_g(b))^2
                 = (4 D A b^2 - sigma a^2 Z^2 N_g)(4 D A b^2 + sigma a^2 Z^2 N_g),
        a genuine degree-4 by degree-4 factorization over Q.

  L21b  On the constructed branch tau_dagger = (1 + 2a^2)/A one has
            delta_dagger = -4 a^4 / A  <  0
        for every nonzero a, hence delta_dagger is never a rational square and
        the L21a degeneration cannot occur there, for any cell.

Also recorded: P is palindromic up to the single monomial 32 A^3 s^2 D^2 b^5,
i.e. P(b) + 32 A^3 s^2 D^2 b^5 has P_i = P_{8-i}.

All arithmetic is exact (fractions); the kernel `h10q._l10_P` is the authority
for the coefficients.  No primality or factorization engine is used, so there
are no refusals in this module.

Replay:  python3 l21_reducible_locus.py     (from math/h10q, or set H10Q_ROOT)
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("H10Q_ROOT", str(HERE)))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the kernel is the authority for P's coefficients)

OUT = Path(os.environ.get("L21_OUT", str(ROOT / "data" / "l21_reducible_locus.jsonl")))
SOURCE = "l21_reducible_locus.py"


def poly_mul(p: list[F], q: list[F]) -> list[F]:
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x == 0:
            continue
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def rational_sqrt(q: F) -> F | None:
    """Exact square root of a nonnegative rational, or None."""
    if q < 0:
        return None
    n, d = q.numerator, q.denominator
    rn, rd = isqrt(n), isqrt(d)
    return F(rn, rd) if rn * rn == n and rd * rd == d else None


def cell_data(a: F, z: F, tau: F) -> tuple[F, F, F, F, F]:
    A = 1 + 4 * a * a
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    delta = 1 - A * tau * tau
    s = (a - 1) / F(2)
    return A, Z, D, delta, s


def ng_coeffs(a: F, A: F) -> list[F]:
    """N_g(b) = 16 a^4 b^2 - A (b-1)^4 as a coefficient list in b."""
    out = [F(0)] * 5
    for i, c in enumerate([1, -4, 6, -4, 1]):  # (b-1)^4
        out[i] -= A * c
    out[2] += 16 * a**4
    return out


def kernel_P(a: F, z: F, tau: F) -> tuple[list[F], tuple[F, F, F, F, F]]:
    A, Z, D, delta, s = cell_data(a, z, tau)
    return h10q._l10_P(a, Z, D, A, delta, s), (A, Z, D, delta, s)


def check_l21a(a: F, z: F, tau: F) -> dict:
    """Verify the difference-of-squares factorization when s=0 and delta is square."""
    P, (A, Z, D, delta, s) = kernel_P(a, z, tau)
    sigma = rational_sqrt(delta)
    row = {
        "type": "l21a",
        "a": str(a),
        "z": str(z),
        "tau": str(tau),
        "A": str(A),
        "delta": str(delta),
        "s_zero": s == 0,
        "delta_square": sigma is not None,
        "source": SOURCE,
    }
    if s != 0 or sigma is None or D == 0 or Z == 0:
        row["applicable"] = False
        return row
    ng = ng_coeffs(a, A)
    lo = [-sigma * a * a * Z * Z * c for c in ng]
    hi = [sigma * a * a * Z * Z * c for c in ng]
    lo[2] += 4 * D * A
    hi[2] += 4 * D * A
    product = poly_mul(lo, hi)
    row.update(
        {
            "applicable": True,
            "degree": len(P) - 1,
            "factor_degrees": [len(lo) - 1, len(hi) - 1],
            "identity_holds": len(product) == len(P)
            and all(product[i] == P[i] for i in range(len(P))),
            "sigma": str(sigma),
        }
    )
    return row


def check_l21b(a: F) -> dict:
    """Verify delta on the constructed branch equals -4a^4/A and is negative."""
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    predicted = -4 * a**4 / A
    return {
        "type": "l21b",
        "a": str(a),
        "A": str(A),
        "tau_dagger": str(tau),
        "delta_dagger": str(delta),
        "predicted_minus_4a4_over_A": str(predicted),
        "identity_holds": delta == predicted,
        "negative": delta < 0,
        "is_rational_square": rational_sqrt(delta) is not None,
        "square_class_equals_minus_A": (delta * -A) == (2 * a * a) ** 2,
        "source": SOURCE,
    }


def check_palindrome(a: F, z: F, tau: F) -> dict:
    """P + 32 A^3 s^2 D^2 b^5 should be palindromic."""
    P, (A, Z, D, delta, s) = kernel_P(a, z, tau)
    R = list(P)
    R[5] += 32 * A**3 * s * s * D * D
    n = len(R) - 1
    return {
        "type": "palindrome",
        "a": str(a),
        "z": str(z),
        "tau": str(tau),
        "degree": n,
        "P0_equals_P8": P[0] == P[n],
        "corrected_palindromic": all(R[i] == R[n - i] for i in range(n + 1)),
        "asymmetry_support": [i for i in range(n + 1) if P[i] != P[n - i]],
        "source": SOURCE,
    }


def main() -> int:
    a_pool = [F(1), F(3), F(5), F(7), F(9), F(11), F(25), F(-3), F(-7)]
    z_pool = [
        F(-2), F(-1), F(1), F(2), F(3), F(5), F(-5), F(7),
        F(-1, 7), F(5, 3), F(2, 11), F(7, 3), F(3, 11), F(-2, 7),
    ]
    square_taus = [F(0), F(1, 3), F(2, 3), F(3, 7), F(4, 9)]

    rows: list[dict] = []

    # L21a: every (z, tau) with a=1 and delta a square must factor.
    l21a_applicable = l21a_ok = 0
    for z in z_pool:
        for tau in square_taus:
            row = check_l21a(F(1), z, tau)
            if row.get("applicable"):
                l21a_applicable += 1
                l21a_ok += bool(row["identity_holds"])
                rows.append(row)

    # Control: a != 1 must NOT be claimed by L21a (s != 0 there).
    control_rows = 0
    for a in a_pool:
        if a == 1:
            continue
        row = check_l21a(a, F(-2), F(0))
        assert row["applicable"] is False, "L21a must not apply when s != 0"
        control_rows += 1

    # L21b: the constructed branch.
    l21b_ok = l21b_neg = 0
    for a in a_pool:
        row = check_l21b(a)
        rows.append(row)
        l21b_ok += bool(row["identity_holds"])
        l21b_neg += bool(row["negative"] and not row["is_rational_square"])

    # Constructed branch is genuinely never in the L21a locus.
    dagger_applicable = 0
    for a in a_pool:
        A = 1 + 4 * a * a
        for z in z_pool:
            row = check_l21a(a, z, (1 + 2 * a * a) / A)
            dagger_applicable += bool(row.get("applicable"))

    # Palindromic structure.
    pal_rows = pal_ok = p0p8 = 0
    for a in a_pool:
        A = 1 + 4 * a * a
        for z in z_pool[:8]:
            for tau in (F(0), 2 * a / A, (1 + 2 * a * a) / A, F(1)):
                if 1 - A * tau * tau == 0:
                    continue
                row = check_palindrome(a, z, tau)
                pal_rows += 1
                pal_ok += bool(row["corrected_palindromic"])
                p0p8 += bool(row["P0_equals_P8"])
                if pal_rows <= 12:
                    rows.append(row)

    summary = {
        "type": "summary",
        "l21a_label": "PROVED",
        "l21a_statement": (
            "s=0 and delta_tau=sigma^2 a rational square imply "
            "P=(4DAb^2-sigma a^2 Z^2 N_g)(4DAb^2+sigma a^2 Z^2 N_g)"
        ),
        "l21a_instances_checked": l21a_applicable,
        "l21a_identity_holds": l21a_ok,
        "l21a_control_rows_correctly_inapplicable": control_rows,
        "l21b_label": "PROVED",
        "l21b_statement": "delta on tau_dagger=(1+2a^2)/A equals -4a^4/A < 0, never a square",
        "l21b_instances_checked": len(a_pool),
        "l21b_identity_holds": l21b_ok,
        "l21b_negative_and_nonsquare": l21b_neg,
        "constructed_branch_in_reducible_locus": dagger_applicable,
        "palindrome_label": "PROVED-on-sample",
        "palindrome_rows": pal_rows,
        "palindrome_corrected_holds": pal_ok,
        "palindrome_P0_equals_P8": p0p8,
        "refusals": 0,
        "refusals_used_as_evidence": False,
        "source": SOURCE,
    }
    rows.append(summary)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    ok = (
        l21a_applicable == l21a_ok
        and l21b_ok == len(a_pool)
        and l21b_neg == len(a_pool)
        and dagger_applicable == 0
        and pal_rows == pal_ok == p0p8
    )
    print(json.dumps(summary, sort_keys=True))
    print(f"wrote {OUT} ({len(rows)} rows)")
    print("VERDICT:", "ALL CHECKS PASS" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
