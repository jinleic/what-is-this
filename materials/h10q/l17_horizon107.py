#!/usr/bin/env python3
"""Replay and serialize the off-grid H-class closure at z = 107*(-1).

All arithmetic and primality decisions come from math/h10q/h10q.py.  The
L12 escape row is injected only in this process because the checked-in table
is intentionally restricted to the on-grid cells.
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
        candidates.append(Path(os.environ["H10Q_ROOT"]))
    cwd = Path.cwd()
    candidates.extend((cwd, cwd / "math" / "h10q"))
    for p in candidates:
        if (p / "h10q.py").is_file():
            return p.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l12_class  # noqa: E402
import l13_filter  # noqa: E402


TARGET = (107, (-1, 1))
W, UT = TARGET
Z = F(W) * F(*UT)
A = F(1)
TAU = F(3, 5)
F_ESC = 107
EPS = 1
Q1 = 31


def replay() -> tuple[dict, dict]:
    """Run L11, L12, kernel, and L13 checks, returning row and diagnostics."""
    assert Z == F(-107)
    # L11h says this target is walled exactly when no numerator prime is 1 mod 4.
    z_factorization = h10q.factorint(abs(Z.numerator))
    l11_escape_primes = h10q._l11_escape_primes(Z)
    assert z_factorization == {107: 1}
    assert l11_escape_primes == []
    assert TARGET not in h10q._L11_CLASSES
    time.sleep(0.1)

    # l12_class_cert is deliberately table-backed.  Install this off-grid
    # escape in-process, replay the complete alignment certificate, then undo
    # the injection before returning.
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(TARGET, sentinel)
    l12_class._L12_ESCAPE[TARGET] = (F_ESC, EPS, Q1)
    try:
        cert = l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(TARGET, None)
        else:
            l12_class._L12_ESCAPE[TARGET] = previous
    assert cert is not None and cert["ok"]
    assert cert["f"] == F_ESC and cert["eps"] == EPS and cert["q"] == Q1
    assert cert["N"] == 141377214454560
    assert cert["S"] == [2, 3, 5, 7, 107]
    assert all(v == 1 for v in cert["syms"].values())
    assert cert["excluded"] == [2, 5, 59, 71, 107, 1259, 56911]
    assert h10q._is_prime(F_ESC)
    assert h10q._is_prime(Q1)
    assert h10q.math.gcd(Q1, cert["N"]) == 1
    time.sleep(0.1)

    # Exact kernel witness for b = eps*f*Q with Q=q1+kN and k=0.
    k_zero = 0
    Q = Q1 + k_zero * cert["N"]
    assert Q == Q1 and h10q._is_prime(Q)
    b = F(EPS * F_ESC * Q)
    A0 = 1 + 4 * A * A
    delta = 1 - A0 * TAU * TAU
    alpha = -delta * A0
    Z3 = Z ** 3
    D = 1 - Z3 - A * A * Z3 * Z3
    s = (A - 1) / 2
    c0 = h10q._sun_h(A, b, Z3)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A0 * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    P = h10q._l10_P(A, Z3, D, A0, delta, s)
    Pb = sum(c * b ** i for i, c in enumerate(P))
    assert M0 == Pb / (b ** 4 * D * D * A0 * A0)
    member_syms = h10q._l12_syms(A, Z, TAU, b)
    assert member_syms is not None and all(v == 1 for v in member_syms.values())
    aligned_syms = l13_filter.aligned(A, Z, TAU, b)
    assert aligned_syms is not None and all(v == 1 for v in aligned_syms.values())
    time.sleep(0.1)

    # L13f: no small emergent place, then a proven-prime cofactor rung.
    smooth_status, smooth_detail = l13_filter.smooth_emergent(A, Z, TAU, b)
    assert smooth_status == "cofactor-big"
    emergent, (rn, rd), smooth = smooth_detail
    assert emergent == [] and not smooth and rd == 1
    verdict, info = l13_filter.cofactor_decide(
        A, Z, TAU, b, detail=smooth_detail
    )
    assert verdict == "zero"
    assert info["rung"] == "prime"
    assert info["parts"][0]["kind"] == "prime"
    assert info["places"] and info["places"][0][1] == 1
    cofactor_prime = int(info["places"][0][0])
    assert h10q._is_prime(cofactor_prime)
    tied = h10q._l7_tied_status(A, b, Z, TAU)
    assert tied is True
    ramified = h10q.ramified(x0, d0)
    assert ramified == []

    row = {
        "family": "ESC",
        "cell": [W, list(UT)],
        "a": "1",
        "eps": EPS,
        "f": F_ESC,
        "q1": Q1,
        "N": str(cert["N"]),
        "k_zero": k_zero,
        "Q": str(Q),
        "rung": info["rung"],
        "tied": str(tied),
        "ramified_empty": not ramified,
        "source": "l17_horizon107.py",
    }
    diagnostics = {
        "status": "VERIFIED",
        "target_z": str(Z),
        "l11": {
            "z_factorization": z_factorization,
            "escape_primes_1_mod_4": l11_escape_primes,
            "table_row_present": TARGET in h10q._L11_CLASSES,
        },
        "l12": {
            "f": F_ESC,
            "eps": EPS,
            "q1": Q1,
            "N": cert["N"],
            "S": cert["S"],
            "ks": cert["ks"],
            "syms": cert["syms"],
            "excluded": cert["excluded"],
        },
        "kernel": {
            "b": str(b),
            "c0": str(c0),
            "M0": str(M0),
            "x0": str(x0),
            "d0": str(d0),
            "member_syms": member_syms,
            "aligned_syms": aligned_syms,
            "ramified": ramified,
            "tied": tied,
        },
        "l13": {
            "smooth_status": smooth_status,
            "small_emergent": emergent,
            "rn": str(rn),
            "rd": str(rd),
            "cofactor_verdict": verdict,
            "cofactor_rung": info["rung"],
            "cofactor_info": info,
            "cofactor_prime_proven": True,
        },
        "closure_row": row,
    }
    return row, diagnostics


def main() -> None:
    row, diagnostics = replay()
    data_path = ROOT / "data" / "l17_horizon107.jsonl"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "type": "meta",
        "schema": "l17-off-grid-closure-v1",
        "task": "Off-grid H-class closure",
        "engine": "math/h10q/h10q.py proven-primality only",
        "cell": [W, list(UT)],
        "family": "ESC",
        "verified": True,
        "replay": "full in-script",
        "source_script": "math/h10q/l17_horizon107.py",
    }
    with data_path.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps(meta, sort_keys=True) + "\n")
        fh.write(json.dumps(row, sort_keys=True) + "\n")

    report_path = Path(os.environ.get("HORIZON107_REPORT", "/tmp/l17_horizon107.report"))
    report_lines = [
        "STATUS: VERIFIED",
        "TARGET: cell [107,[-1,1]], z=-107, family=ESC",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        "L11: z factorization = " + json.dumps(diagnostics["l11"]["z_factorization"], sort_keys=True),
        "L11: numerator escape primes (1 mod 4) = " + json.dumps(diagnostics["l11"]["escape_primes_1_mod_4"]),
        "L11: checked-in class row present = " + str(diagnostics["l11"]["table_row_present"]),
        "L12: " + json.dumps(diagnostics["l12"], sort_keys=False),
        "KERNEL: " + json.dumps(diagnostics["kernel"], sort_keys=False),
        "L13: " + json.dumps(diagnostics["l13"], sort_keys=False),
        "CLOSURE ROW: " + json.dumps(row, sort_keys=True),
        "ARTIFACT: " + str(data_path),
        "REPLAY: python math/h10q/l17_horizon107.py (run from repo or set H10Q_ROOT)",
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps(row, sort_keys=True))
    print(f"WROTE {data_path}")
    print(f"WROTE {report_path}")


if __name__ == "__main__":
    main()
