#!/usr/bin/env python3
"""Replay the established off-grid closure protocol at (w,z)=(131,-131).

The L11 root route is exhausted exactly.  The only escape factor is f=131;
for its canonical a=1, tau=3/5 class, the frozen symbol at 5 and the wild
Legendre symbol are opposite for every admissible prime q1.  Thus this script
emits a precise negative protocol report, not a closure row and not a claim
that the cell is impossible by every future method.  All primality decisions
are delegated to math/h10q/h10q.py.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path


def find_root() -> Path:
    candidates = []
    if os.environ.get("H10Q_ROOT"):
        candidates.append(Path(os.environ["H10Q_ROOT"]).expanduser())
    here = Path(__file__).resolve().parent
    cwd = Path.cwd()
    candidates.extend((cwd, cwd / "math" / "h10q", here))
    for candidate in candidates:
        if (candidate / "h10q.py").is_file():
            return candidate.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l12_class  # noqa: E402


CELL = (131, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
A = F(1)
TAU = F(3, 5)
F_ESC = W
PROBE_Q_MAX = 1000
OUT = ROOT / "data" / "l17_horizon131.jsonl"
REPORT = Path(os.environ.get("HORIZON131_REPORT", "/tmp/l17_horizon131.report"))


def emit(handle, record: dict) -> None:
    # Integer prime keys and the string ``oo`` coexist in symbol maps, so do
    # not ask json.dumps to recursively sort dictionary keys.
    handle.write(json.dumps(record) + "\n")
    handle.flush()


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def cert_summary(cert: dict | None) -> dict:
    if cert is None:
        return {}
    return {
        "ok": bool(cert["ok"]),
        "N": str(cert["N"]),
        "S": cert["S"],
        "ks": cert["ks"],
        "syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "excluded": cert.get("excluded", []),
    }


def fresh_l12_cert(eps: int, q1: int):
    """Run the canonical L12 certificate with one ephemeral off-grid row."""
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(CELL, sentinel)
    l12_class._L12_ESCAPE[CELL] = (F_ESC, eps, q1)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(CELL, None)
        else:
            l12_class._L12_ESCAPE[CELL] = previous


def discover_l11() -> tuple[list[int], list[int], dict]:
    """Exhaust the prescribed roots r mod w and their odd lifted a-values."""
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_pool = sorted(
        {
            abs(r + sign * W * lift)
            for r in roots
            for sign in (1, -1)
            for lift in range(3)
            if r + sign * W * lift != 0
        }
    )
    a_pool = [a for a in a_pool if a % 2 == 1]
    minus_one_character = h10q.legendre(F(-1), W)
    assert h10q._is_prime(W)
    assert W % 4 == 3 and minus_one_character == -1
    assert roots == [] and a_pool == []
    proof = {
        "prime_w_proven": True,
        "tested_r": "1<=r<131",
        "equation": "1+4*r^2 == 0 (mod 131)",
        "roots": roots,
        "odd_a_pool": a_pool,
        "legendre_minus_one_mod_w": minus_one_character,
        "reason": "a root would make -1 a square modulo the proven prime 131",
    }
    return roots, a_pool, proof


def character_mod_5(value: int) -> int:
    """Quadratic character of a nonzero residue modulo 5."""
    residue = value % 5
    assert residue != 0
    return 1 if pow(residue, 2, 5) == 1 else -1


def escape_obstruction_proof() -> dict:
    """Give the exact valuation/residue proof blocking every ESC class.

    Put b=eps*131*q1.  For q1 != 5, b is a 5-adic unit.  With Z3=z^3,
    D=1-Z3-Z3^2, and t=(b-1)^2/b, the Sun kernel has

        g=(16-5*t^2)/5,  c=Z3^2*g/D.

    Here Z3=4 and D=1 modulo 5, while 16-5*t^2=1 modulo 5.  Consequently
    v5(c)=-1, v5(M)=-3, v5(x)=-3, and v5(d)=0.  The Hilbert symbol at 5 is
    therefore the quadratic character of the unit d=8*eps*131*q1.  Since
    eps and 131 are squares modulo 5 but 8 is not, it equals -(q1|5).
    Quadratic reciprocity gives (q1|5)=(5|q1), the wild class symbol.
    """
    z3 = Z**3
    D = 1 - z3 - z3 * z3
    z3_mod_5 = h10q.unit_mod(z3, 5)
    D_mod_5 = h10q.unit_mod(D, 5)
    g_numerator_mod_5 = 16 % 5
    assert z3_mod_5 == 4
    assert D_mod_5 == 1
    assert g_numerator_mod_5 == 1
    assert h10q.vp(z3, 5) == 0 and h10q.vp(D, 5) == 0

    residue_table = []
    for eps in (1, -1):
        assert character_mod_5(eps) == 1
        for q1_mod_5 in (1, 2, 3, 4):
            d_unit_mod_5 = (8 * eps * W * q1_mod_5) % 5
            hilbert_5 = character_mod_5(d_unit_mod_5)
            wild_q1 = character_mod_5(q1_mod_5)
            assert hilbert_5 == -wild_q1
            residue_table.append(
                {
                    "eps": eps,
                    "q1_mod_5": q1_mod_5,
                    "d_unit_mod_5": d_unit_mod_5,
                    "hilbert_5": hilbert_5,
                    "wild_q1": wild_q1,
                    "product": hilbert_5 * wild_q1,
                }
            )

    return {
        "scope": "every prime q1 admissible to the canonical ESC class",
        "canonical_data": {
            "a": str(A),
            "tau": str(TAU),
            "f": F_ESC,
            "eps": [1, -1],
            "b": "eps*131*q1",
        },
        "excluded_q1": [2, 3, 5, 7, 131],
        "valuation_proof": {
            "Z3_mod_5": z3_mod_5,
            "D_mod_5": D_mod_5,
            "g_numerator_mod_5": g_numerator_mod_5,
            "v5_c": -1,
            "v5_M": -3,
            "v5_x": -3,
            "v5_d": 0,
        },
        "reciprocity_identity": "hilbert_5 = -(q1|5) = -(5|q1) = -wild_q1",
        "residue_table": residue_table,
        "conclusion": "the frozen p=5 and wild q1 symbols cannot both be +1",
    }


def probe_escape_classes() -> tuple[list[dict], list[dict]]:
    """Replay all proven-prime attempts through q1<1001.

    This bounded log illustrates the exact obstruction; the universal residue
    proof above, not the absence of a hit in this list, justifies the negative
    protocol conclusion.
    """
    attempts = []
    aligned = []
    S = {2, 3, 5, 7, F_ESC}
    for eps in (1, -1):
        for q1 in h10q.primerange(2, PROBE_Q_MAX + 1):
            time.sleep(0.1)
            assert h10q._is_prime(q1)
            record = {
                "eps": eps,
                "f": F_ESC,
                "q1": q1,
                "q1_proven_prime": True,
            }
            if q1 in S:
                cert = fresh_l12_cert(eps, q1)
                assert cert is None
                record.update(
                    {
                        "status": "ineligible:q1-in-controlled-set",
                        "controlled_set": sorted(S),
                    }
                )
                attempts.append(record)
                continue
            try:
                cert = fresh_l12_cert(eps, q1)
            except Exception as exc:
                record["status"] = refusal(exc)
                attempts.append(record)
                continue
            record.update(cert_summary(cert))
            if cert is None:
                record["status"] = "refused:none"
                attempts.append(record)
                continue

            b = F(eps * F_ESC * q1)
            z3 = Z**3
            D = 1 - z3 - z3 * z3
            c0 = h10q._sun_h(A, b, z3)
            assert c0 is not None
            A0 = 1 + 4 * A * A
            delta = 1 - A0 * TAU * TAU
            alpha = -delta * A0
            M0 = 16 - delta * c0 * c0
            x0, d0 = alpha * M0, alpha * 2 * b
            assert h10q.vp(c0, 5) == -1
            assert h10q.vp(M0, 5) == -3
            assert h10q.vp(x0, 5) == -3
            assert h10q.vp(d0, 5) == 0
            assert h10q.vp(D, 5) == 0

            hilbert_5 = h10q.hilbert(x0, d0, 5)
            wild_q1 = h10q.legendre(A0, q1)
            assert cert["syms"][5] == hilbert_5
            assert cert["syms"]["q1"] == wild_q1
            assert hilbert_5 == -wild_q1
            assert hilbert_5 * wild_q1 == -1
            assert not cert["ok"]
            record.update(
                {
                    "status": "obstructed:p5-versus-wild",
                    "q1_mod_5": q1 % 5,
                    "v5_c": h10q.vp(c0, 5),
                    "v5_M": h10q.vp(M0, 5),
                    "v5_x": h10q.vp(x0, 5),
                    "v5_d": h10q.vp(d0, 5),
                    "hilbert_5": hilbert_5,
                    "wild_q1": wild_q1,
                    "symbol_product": hilbert_5 * wild_q1,
                }
            )
            attempts.append(record)
            if cert["ok"]:
                aligned.append(
                    {
                        "family": "ESC-fresh",
                        "a": str(A),
                        "eps": eps,
                        "f": F_ESC,
                        "q1": q1,
                        "tau": str(TAU),
                        "cert": cert,
                    }
                )
    return attempts, aligned


def main() -> None:
    h10q._selftest()
    assert Z == F(-131)
    assert h10q.factorint(abs(Z.numerator)) == {131: 1}
    assert h10q._l11_escape_primes(Z) == []
    assert CELL not in h10q._L11_CLASSES
    assert CELL not in l12_class._L12_ESCAPE

    roots, a_pool, l11_proof = discover_l11()
    obstruction = escape_obstruction_proof()
    attempts, aligned = probe_escape_classes()

    refusals = [record for record in attempts if record["status"].startswith("refused:")]
    eligible = [
        record
        for record in attempts
        if record["status"] == "obstructed:p5-versus-wild"
    ]
    ineligible = [
        record
        for record in attempts
        if record["status"] == "ineligible:q1-in-controlled-set"
    ]
    assert not refusals
    assert eligible
    assert not aligned
    assert all(record["symbol_product"] == -1 for record in eligible)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        emit(
            handle,
            {
                "type": "meta",
                "schema": "l17-off-grid-negative-v1",
                "task": "Off-grid H-class closure attempt",
                "engine": "math/h10q/h10q.py proven-primality only",
                "cell": [W, list(UT)],
                "z": str(Z),
                "source_script": "math/h10q/l17_horizon131.py",
                "closure_row_emitted": False,
                "protocol_status": "PROVED-OBSTRUCTED",
                "cell_status": "OPEN",
                "probe_q_max": PROBE_Q_MAX,
                "refusals_are_not_evidence": True,
            },
        )
        emit(
            handle,
            {
                "type": "l11-discovery",
                "family": "L11-fresh",
                "cell": [W, list(UT)],
                **l11_proof,
            },
        )
        emit(
            handle,
            {
                "type": "protocol-obstruction",
                "family": "ESC-fresh",
                "cell": [W, list(UT)],
                **obstruction,
            },
        )
        for record in attempts:
            emit(
                handle,
                {
                    "type": "class-attempt",
                    "family": "ESC-fresh",
                    "cell": [W, list(UT)],
                    "a": str(A),
                    "tau": str(TAU),
                    **record,
                },
            )
        emit(
            handle,
            {
                "type": "member-phase",
                "family": "ESC-fresh",
                "cell": [W, list(UT)],
                "status": "not-entered",
                "reason": "the member ladder requires an aligned class; none exists in this protocol",
                "aligned_classes": 0,
            },
        )
        emit(
            handle,
            {
                "type": "summary",
                "cell": [W, list(UT)],
                "status": "OPEN",
                "protocol_result": "PROVED-OBSTRUCTED",
                "closure_row_emitted": False,
                "l11_roots": roots,
                "l11_a_pool": a_pool,
                "class_attempts": len(attempts),
                "eligible_attempts": len(eligible),
                "ineligible_attempts": len(ineligible),
                "refusals": len(refusals),
                "aligned_classes": len(aligned),
                "reason": "no L11 root; canonical ESC alignment is universally blocked by hilbert_5=-wild_q1",
                "scope_note": "this is not a proof that the cell is impossible by every future family",
            },
        )

    report_lines = [
        "STATUS: OPEN (no closure row)",
        "LABEL: PROVED protocol obstruction; the cell itself is not proved impossible",
        "TARGET: cell [131,[-1,1]], z=-131",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        "PROVED L11: " + json.dumps(l11_proof, sort_keys=True),
        "PROVED ESC: " + json.dumps(obstruction, sort_keys=True),
        "ATTEMPT LOG: "
        + json.dumps(
            {
                "q1_range": f"proven primes 2..{PROBE_Q_MAX}",
                "eps": [1, -1],
                "attempts": len(attempts),
                "eligible": len(eligible),
                "ineligible": len(ineligible),
                "refusals": len(refusals),
                "aligned_classes": len(aligned),
            },
            sort_keys=True,
        ),
        "MEMBER LADDER: not entered because the protocol produced no aligned class",
        "RESULT: no VERIFIED closure at [131,[-1,1]] by the established L11/ESC recipe",
        "SCOPE: the target remains OPEN; no claim is made against other future constructions",
        "ARTIFACT: " + str(OUT),
        "REPLAY: python3 math/h10q/l17_horizon131.py (run from repo or set H10Q_ROOT)",
    ]
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("STATUS OPEN; PROVED protocol obstruction")
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
