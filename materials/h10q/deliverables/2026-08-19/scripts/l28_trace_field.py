#!/usr/bin/env python3
"""L28 local trace-field closure for the reciprocal even pullback.

The producer replays the exact centered quartic, its Ferrari resolvent,
the uniform dyadic congruences, and the norm squareclass certificate.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

from l24_diagonal_geometry import irreducibility_certificate
from l25_scaled_coupling import squareclass_2adic, vp
from l26_reciprocal_tie import evaluate, multiply, reciprocal_polynomials

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l28_trace_field.jsonl"
REPORT = Path("/tmp/l28_trace_field.md")


def power_of_two(exponent: int) -> F:
    return F(2) ** exponent


def unit_mod_8(value: F) -> int:
    assert vp(value, 2) == 0
    return value.numerator * pow(value.denominator, -1, 8) % 8


def centered_quartic(a: F, Z: F) -> tuple[F, F, F, F, F, list[F]]:
    A = 1 + 4 * a * a
    D = 1 - Z - a * a * Z * Z
    assert a and Z and D
    c4 = 4 * a**8 * A * Z**4
    c2 = -128 * a**12 * Z**4
    c1 = -32 * A**3 * Z**2
    c0 = (16 * A * A * D * D + 1024 * a**16 * Z**4 / A
          - 128 * A**3 * Z**2)
    centered = [c0, c1, c2, F(0), c4]

    delta = -4 * a**4 / A
    _, trace = reciprocal_polynomials(a, Z, delta)
    translated = [F(0)] * 5
    for degree, coefficient in enumerate(trace):
        for w_degree in range(degree + 1):
            translated[w_degree] += (
                coefficient * math.comb(degree, w_degree)
                * F(2) ** (degree - w_degree)
            )
    assert translated == centered
    return A, D, c4, c2, c1, centered


def resolvent_coefficients(a: F, Z: F) -> tuple[F, F, F, F]:
    A, D, c4, c2, c1, centered = centered_quartic(a, Z)
    p = c2 / c4
    q = c1 / c4
    r = centered[0] / c4
    h = p * p - 4 * r
    assert p == -32 * a**4 / A
    assert q == -8 * A**2 / (a**8 * Z**2)
    assert r == (
        4 * A * D**2 / (a**8 * Z**4)
        + 256 * a**8 / A**2
        - 32 * A**2 / (a**8 * Z**2)
    )
    return p, q, r, h


def resolvent_value(p: F, q: F, h: F, Y: F) -> F:
    return Y**3 + 2 * p * Y**2 + h * Y - q**2


def centered_identity_replay() -> dict:
    checks = 0
    for a in map(F, (1, 3, 5, -1)):
        for z in (F(1), F(3), F(2), F(1, 2), F(5, 4)):
            Z = z * z
            A, D, c4, c2, c1, centered = centered_quartic(a, Z)
            assert centered[-1] == c4 and centered[2] == c2
            assert centered[1] == c1 and centered[3] == 0
            assert A % 2 and D
            checks += 1
    return {
        "type": "centered-trace-quartic",
        "label": "PROVED exact identity",
        "instances": checks,
        "shift": "w=u-2",
        "leading": "4*a^8*A*Z^4",
        "coefficients_descending": [
            "4*a^8*A*Z^4",
            "0",
            "-128*a^12*Z^4",
            "-32*A^3*Z^2",
            "16*A^2*D^2+1024*a^16*Z^4/A-128*A^3*Z^2",
        ],
    }


def local_irreducibility_replay() -> dict:
    rows = 0
    subcases: dict[str, int] = {}
    for t in (-6, -4, -2, 0, 2, 4, 6):
        for a in (F(1), F(3), F(5), F(1, 3), F(5, 3)):
            assert vp(a, 2) == 0
            for z_unit in (F(1), F(3), F(5), F(1, 3), F(5, 3)):
                Z = power_of_two(t) * z_unit * z_unit
                A, D, c4, c2, c1, centered = centered_quartic(a, Z)
                assert vp(A, 2) == 0 and unit_mod_8(A) == 5
                assert vp(D, 2) == (0 if t >= 0 else 2 * t)
                expected_c0 = 4 if t >= 0 else 4 + 4 * t
                assert [vp(centered[index], 2) for index in (0, 1, 2, 4)] == [
                    expected_c0, 5 + 2 * t, 7 + 4 * t, 2 + 4 * t,
                ]

                p, q, r, h = resolvent_coefficients(a, Z)
                assert vp(p, 2) == 5
                assert vp(q, 2) == 3 - 2 * t
                assert vp(r, 2) == (2 - 4 * t if t >= 0 else 2)
                h_scale = 4 - 4 * t if t >= 0 else 4
                assert vp(h, 2) == h_scale
                assert unit_mod_8(h / power_of_two(h_scale)) == 3
                q2_scale = 6 - 4 * t
                assert unit_mod_8(q * q / power_of_two(q2_scale)) == 1

                if t == 0:
                    cases = [("t=0,vY=2", 2, 6, 3)]
                elif t > 0:
                    cases = [
                        ("t>0,vY=2", 2, 6 - 4 * t, 2),
                        ("t>0,vY=2-2t", 2 - 2 * t, 6 - 6 * t, 4),
                    ]
                else:
                    cases = [
                        ("t<0,vY=2", 2, 6, 4),
                        ("t<0,vY=2-4t", 2 - 4 * t, 6 - 4 * t, 2),
                    ]
                for name, y_valuation, scale, expected_residue in cases:
                    for y_unit in map(F, (1, 3, 5, 7, 9, 11, 13, 15)):
                        Y = power_of_two(y_valuation) * y_unit * y_unit
                        normalized = resolvent_value(p, q, h, Y) / power_of_two(scale)
                        assert unit_mod_8(normalized / power_of_two(vp(normalized, 2))) in (1, 3, 5, 7)
                        assert normalized.numerator * pow(normalized.denominator, -1, 8) % 8 == expected_residue
                        rows += 1
                        subcases[name] = subcases.get(name, 0) + 1
    return {
        "type": "uniform-trace-irreducibility",
        "label": "PROVED over Q_2, hence over Q",
        "finite_unit_rows": rows,
        "finite_unit_subcases": subcases,
        "quartic_newton_slopes": {
            "t>=0": "t-1/2",
            "t<0": "-1/2",
        },
        "linear_factor_exclusion": "all root valuations are half-integral",
        "ferrari_resolvent": "Y^3+2*p*Y^2+(p^2-4*r)*Y-q^2",
        "quadratic_factor_exclusion": "every possible square Y has normalized residue 2, 3, or 4 mod 8",
        "guards": "v2(a)=0, Z!=0, t=v2(Z) even; D!=0 follows from A nonsquare",
    }


def norm_squareclass_replay() -> dict:
    rows = 0
    minimum_correction_valuation = 10**9
    sample = None
    for a in (F(1), F(3), F(5), F(1, 3), F(5, 3)):
        for t in (-6, -4, -2, 0, 2, 4, 6):
            for z_unit in (F(1), F(3), F(5), F(1, 3), F(5, 3)):
                Z = power_of_two(t) * z_unit * z_unit
                A = 1 + 4 * a * a
                D = 1 - Z - a * a * Z * Z
                delta = -4 * a**4 / A
                _, trace = reciprocal_polynomials(a, Z, delta)
                leading = trace[-1]
                trace_at_minus_two = evaluate(trace, F(-2))
                norm = 16 * trace_at_minus_two / leading
                q0 = a**4 - A
                E = A**3 * D**2 + 64 * a**8 * Z**4 * q0**2
                expected = 64 * E / (a**8 * Z**4 * A**2)
                correction = E / (A**3 * D**2) - 1
                assert norm == expected
                assert vp(q0, 2) == 2
                assert vp(correction, 2) >= 10
                assert squareclass_2adic(norm / A)
                assert not squareclass_2adic(norm)
                minimum_correction_valuation = min(
                    minimum_correction_valuation, vp(correction, 2)
                )
                rows += 1
                if sample is None and a == 1 and Z == 9:
                    sample = {
                        "a": str(a), "Z": str(Z), "A": str(A),
                        "D": str(D), "norm": str(norm),
                        "v2_norm": vp(norm, 2),
                    }
    assert sample is not None
    return {
        "type": "bad-character-norm-squareclass",
        "label": "PROVED: Norm(2*(theta+2)) lies in A*Q_2^{x2}",
        "identity": (
            "Norm=64*(A^3*D^2+64*a^8*Z^4*(a^4-A)^2)"
            "/(a^8*Z^4*A^2)"
        ),
        "correction": "E/(A^3*D^2) is in 1+2^10*Z_2",
        "minimum_replayed_correction_valuation": minimum_correction_valuation,
        "finite_unit_rows": rows,
        "sample": sample,
    }


def capell_certificate() -> dict:
    a, Z = F(1), F(9)
    A = 1 + 4 * a * a
    _, trace = reciprocal_polynomials(a, Z, -4 * a**4 / A)
    base = [F(-2), F(0), F(1, 2)]
    composition = [F(0)]
    power = [F(1)]
    for coefficient in trace:
        if len(composition) < len(power):
            composition.extend([F(0)] * (len(power) - len(composition)))
        for index, value in enumerate(power):
            composition[index] += coefficient * value
        power = multiply(power, base)
    while len(composition) > 1 and composition[-1] == 0:
        composition.pop()
    certificate = irreducibility_certificate(composition, 41)
    assert certificate["degree"] == 8
    return {
        "type": "capell-bad-sign-cover",
        "label": "PROVED uniformly by L28 irreducibility + norm; exact mod-41 sample replayed",
        "polynomial": "T(v^2/2-2)",
        "uniform_degree": 8,
        "sample": {"a": 1, "Z": 9, "prime": 41},
        "sample_certificate_sha256": certificate["primitive_sha256"],
        "consequence": "the bad-sign quadratic extension is always nontrivial on the canonical L26 domain",
    }


def frontier_record() -> dict:
    return {
        "type": "strict-frontier",
        "label": "OPEN globally",
        "closed": [
            "the trace quartic is irreducible over Q_2 for every canonical L26 parameter",
            "2*(theta+2) is uniformly nonsquare in the trace field",
            "the Capell bad-sign cover has degree 8",
        ],
        "open": [
            "a globally good reciprocal member b",
            "the reciprocal lift squareclass theta^2-4",
            "the two-large-bad-divisor parity estimate",
            "any unconditional H10/Q or quantifier improvement",
        ],
        "chain_changed": False,
    }


def build_report(records: list[dict], elapsed: float) -> str:
    lines = [
        "# L28 trace-field replay", "",
        "PROVED labels are exact identities, Newton-polygon deductions, or",
        "finite residue replays of the stated uniform congruence lemmas.",
        "No global member theorem or unconditional H10/Q consequence is claimed.", "",
    ]
    for record in records:
        lines.extend(["## " + record["type"], "", "```json",
                      json.dumps(record, indent=2, sort_keys=True), "```", ""])
    lines.append(f"elapsed_seconds: {elapsed:.3f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [
        centered_identity_replay(),
        local_irreducibility_replay(),
        norm_squareclass_replay(),
        capell_certificate(),
        frontier_record(),
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
