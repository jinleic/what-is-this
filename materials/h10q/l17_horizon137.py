#!/usr/bin/env python3
"""Fresh L11/L12 discovery and proven L13 closure for (137,-137).

The target lies beyond the frozen w<=97 grid.  This replay discovers the L11
root classes directly, probes the L12 escape route with an ephemeral table row,
and verifies one member using only the proven-primality h10q engine and the L13
cofactor ladder.  Refusals are logged but are never treated as evidence.
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
    candidates.extend(
        (
            here,
            cwd,
            cwd / "math" / "h10q",
            Path("/Users/jinleic/jinleic-workspace/math/h10q"),
        )
    )
    for candidate in candidates:
        if (candidate / "h10q.py").is_file():
            return candidate.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402
import l12_class  # noqa: E402
import l13_filter  # noqa: E402

CELL = (137, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
OUT = ROOT / "data" / "l17_horizon137.jsonl"
REPORT = Path(os.environ.get("HORIZON137_REPORT", "/tmp/l17_horizon137.report"))
CLASS_Q_MAX = int(os.environ.get("H137_CLASS_Q_MAX", "149"))
MAX_K = int(os.environ.get("H137_MAX_K", "0"))
A_T_MAX = int(os.environ.get("H137_A_T_MAX", "2"))
MEMBER_SECONDS = float(os.environ.get("H137_MEMBER_SECONDS", "900"))


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def cert_summary(cert):
    if cert is None:
        return {}
    return {
        "N": str(cert["N"]),
        "S": cert["S"],
        "ks": cert["ks"],
        "syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "excluded": cert.get("excluded", []),
    }


def fresh_l12_cert(eps: int, q1: int):
    """Call the canonical L12 cert with one temporary off-grid escape row."""
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(CELL, sentinel)
    l12_class._L12_ESCAPE[CELL] = (W, eps, q1)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(CELL, None)
        else:
            l12_class._L12_ESCAPE[CELL] = previous


def discover_l11():
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_pool = sorted(
        {
            abs(r + sign * W * t)
            for r in roots
            for sign in (1, -1)
            for t in range(A_T_MAX + 1)
            if r + sign * W * t != 0
        }
    )
    a_pool = [a for a in a_pool if a % 2 == 1]
    attempts = []
    rows = []
    for a_int in a_pool:
        a = F(a_int)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        for eps in (1, -1):
            for q1 in h10q.primerange(3, CLASS_Q_MAX + 1):
                time.sleep(0.1)
                try:
                    cert = h10q._l10_class_cert(a, Z, tau, eps, q1)
                    status = (
                        "aligned"
                        if cert is not None and cert.get("ok")
                        else "cert-not-ok"
                        if cert is not None
                        else "refused:none"
                    )
                except Exception as exc:
                    cert = None
                    status = refusal(exc)
                attempt = {
                    "a": a_int,
                    "eps": eps,
                    "q1": q1,
                    "status": status,
                }
                attempt.update(cert_summary(cert))
                attempts.append(attempt)
                if cert is not None and cert.get("ok"):
                    rows.append(
                        {
                            "family": "L11-fresh",
                            "a": str(a),
                            "eps": eps,
                            "f": 1,
                            "q1": q1,
                            "tau": str(tau),
                            "cert": cert,
                        }
                    )
    return roots, a_pool, attempts, rows


def discover_l12():
    attempts = []
    rows = []
    for eps in (1, -1):
        for q1 in h10q.primerange(3, CLASS_Q_MAX + 1):
            time.sleep(0.1)
            try:
                cert = fresh_l12_cert(eps, q1)
                status = (
                    "aligned"
                    if cert is not None and cert.get("ok")
                    else "cert-not-ok"
                    if cert is not None
                    else "refused:none"
                )
            except Exception as exc:
                cert = None
                status = refusal(exc)
            attempt = {"f": W, "eps": eps, "q1": q1, "status": status}
            attempt.update(cert_summary(cert))
            attempts.append(attempt)
            if cert is not None and cert.get("ok"):
                rows.append(
                    {
                        "family": "ESC-fresh",
                        "a": "1",
                        "eps": eps,
                        "f": W,
                        "q1": q1,
                        "tau": str(F(3, 5)),
                        "cert": cert,
                    }
                )
    return attempts, rows


def member_decide(row, k: int, deadline: float):
    cert = row["cert"]
    q1 = int(row["q1"])
    N = int(cert["N"])
    eps = int(row["eps"])
    f = int(row["f"])
    a = F(row["a"])
    tau = F(row["tau"])
    Q = q1 + k * N
    result = {"k": k, "Q": str(Q), "q1": q1, "N": str(N)}
    if Q in set(cert.get("excluded", [])):
        result["verdict"] = "excluded"
        return result, None
    try:
        prime = h10q._is_prime(Q)
    except Exception as exc:
        result["verdict"] = refusal(exc)
        return result, None
    if not prime:
        result["verdict"] = "nonprime"
        return result, None
    b = F(eps * f * Q)
    try:
        smooth_status, detail = l13_filter.smooth_emergent(a, Z, tau, b)
    except Exception as exc:
        result["verdict"] = refusal(exc)
        return result, None
    result["smooth_status"] = smooth_status
    result["smooth_detail"] = detail
    if smooth_status == "zero":
        result.update({"verdict": "zero", "rung": "square"})
        return result, (k, Q, "square", b, None)
    if smooth_status == "cofactor-big":
        try:
            verdict, info = l13_filter.cofactor_decide(
                a, Z, tau, b, detail=detail, deadline=deadline
            )
        except Exception as exc:
            result["verdict"] = refusal(exc)
            return result, None
        result["verdict"] = verdict
        result["rung"] = info.get("rung") if isinstance(info, dict) else None
        result["cofactor_info"] = info
        if verdict == "zero":
            return result, (k, Q, result["rung"], b, info)
        return result, None
    result["verdict"] = smooth_status
    return result, None


def close_certificate(row, hit):
    """Independently replay every fact used by the emitted closure row."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])
    tau = F(row["tau"])
    cert = row["cert"]

    assert cert["ok"] and all(value == 1 for value in cert["syms"].values())
    assert Q == q1 + k * int(cert["N"])
    assert b == F(eps * f * Q)
    assert h10q._is_prime(Q)
    if f != 1:
        assert h10q._is_prime(f)
    if row["family"] == "L11-fresh":
        cert_replay = h10q._l10_class_cert(a, Z, tau, eps, q1)
    else:
        cert_replay = fresh_l12_cert(eps, q1)
    assert cert_replay == cert and cert_replay["ok"]

    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z3 = Z**3
    D = 1 - Z3 - a * a * Z3 * Z3
    s = (a - 1) / 2
    c0 = h10q._sun_h(a, b, Z3)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    P = h10q._l10_P(a, Z3, D, A, delta, s)
    Pb = sum(coefficient * b**i for i, coefficient in enumerate(P))
    assert M0 == Pb / (b**4 * D * D * A * A)

    member_syms = h10q._l12_syms(a, Z, tau, b)
    aligned_syms = l13_filter.aligned(a, Z, tau, b)
    assert member_syms is not None and all(value == 1 for value in member_syms.values())
    assert aligned_syms is not None and all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, Z, tau, b)
    assert smooth_status == "cofactor-big"
    emergent, (rn, rd), smooth = smooth_detail
    assert emergent == [] and not smooth and rd == 1
    verdict, info = l13_filter.cofactor_decide(
        a, Z, tau, b, detail=smooth_detail, deadline=time.time() + 600
    )
    assert verdict == "zero" and info["rung"] == rung
    assert prior_info is not None and prior_info["rn"] == info["rn"]
    assert info["places"]
    assert all(place[1] == 1 and place[2] == "proved" for place in info["places"])
    if info["rung"] == "prime":
        assert len(info["places"]) == 1
        assert h10q._is_prime(int(info["places"][0][0]))

    tied = h10q._l7_tied_status(a, b, Z, tau)
    assert tied is True
    ramified = h10q.ramified(x0, d0)
    assert ramified == []

    closure = {
        "family": "L11" if row["family"] == "L11-fresh" else "ESC",
        "cell": [W, list(UT)],
        "a": str(a),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": str(tied),
        "ramified_empty": not ramified,
        "source": "l17_horizon137.py",
    }
    replay = {
        "type": "kernel-replay",
        "cell": [W, list(UT)],
        "family": row["family"],
        "a": str(a),
        "b": str(b),
        "tau": str(tau),
        "z": str(Z),
        "q1": q1,
        "Q": str(Q),
        "k": k,
        "N": str(cert["N"]),
        "S": cert["S"],
        "ks": cert["ks"],
        "class_syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "A": str(A),
        "delta": str(delta),
        "alpha": str(alpha),
        "D": str(D),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "P_degree": len(P) - 1,
        "kernel_identity": True,
        "member_syms": member_syms,
        "aligned_syms": aligned_syms,
        "small_emergent": emergent,
        "stripped_rn": str(rn),
        "stripped_rd": str(rd),
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "cofactor_prime_proven": info["rung"] == "prime",
        "tied": tied,
        "ramified": ramified,
        "ramified_empty": not ramified,
    }
    return closure, replay


def emit(out, record):
    # Mixed integer/string Hilbert-place keys preclude sort_keys=True.
    out.write(json.dumps(record) + "\n")


def write_artifact(
    roots,
    a_pool,
    l11_attempts,
    l12_attempts,
    rows,
    member_records,
    closure,
    replay,
):
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        emit(
            out,
            {
                "type": "meta",
                "schema": "l17-off-grid-closure-v1",
                "task": "Off-grid H-class closure",
                "engine": "math/h10q/h10q.py proven-primality only",
                "proven_primality_only": True,
                "selftest": True,
                "cell": [W, list(UT)],
                "z": str(Z),
                "verified": closure is not None,
                "replay": "full in-script",
                "source_script": "math/h10q/l17_horizon137.py",
                "l11_roots": roots,
                "l11_a_pool": a_pool,
                "class_q_max": CLASS_Q_MAX,
                "max_k": MAX_K,
            },
        )
        for attempt in l11_attempts:
            emit(
                out,
                {
                    "type": "class-attempt",
                    "family": "L11-fresh",
                    "cell": [W, list(UT)],
                    **attempt,
                },
            )
        for attempt in l12_attempts:
            emit(
                out,
                {
                    "type": "class-attempt",
                    "family": "ESC-fresh",
                    "cell": [W, list(UT)],
                    **attempt,
                },
            )
        for row in rows:
            cert = row["cert"]
            emit(
                out,
                {
                    "type": "aligned-class",
                    "family": row["family"],
                    "cell": [W, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    "tau": row["tau"],
                    "N": str(cert["N"]),
                    "S": cert["S"],
                    "ks": cert["ks"],
                    "syms": cert["syms"],
                    "Q0": str(cert["Q0"]),
                    "excluded": cert.get("excluded", []),
                    "ok": cert["ok"],
                },
            )
        for record in member_records:
            emit(out, record)
        if closure is None:
            emit(
                out,
                {
                    "type": "summary",
                    "status": "OPEN",
                    "cell": [W, list(UT)],
                    "reason": "no proved zero member in the bounded logged search",
                },
            )
        else:
            emit(out, replay)
            emit(out, closure)
            emit(
                out,
                {
                    "type": "summary",
                    "status": "VERIFIED",
                    "cell": [W, list(UT)],
                    "family": closure["family"],
                    "k_zero": closure["k_zero"],
                    "Q": closure["Q"],
                },
            )


def write_report(
    roots,
    a_pool,
    l11_attempts,
    l12_attempts,
    l11_rows,
    l12_rows,
    member_records,
    closure,
    replay,
):
    status = "VERIFIED" if closure is not None else "OPEN"
    lines = [
        "H137 off-grid closure discovery and replay",
        "STATUS: " + status,
        "TARGET: cell [137,[-1,1]], z=-137",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        "L11 DISCOVERY: "
        + json.dumps(
            {
                "roots": roots,
                "a_pool": a_pool,
                "attempts": len(l11_attempts),
                "aligned_rows": len(l11_rows),
                "q_max": CLASS_Q_MAX,
            }
        ),
        "L12 PROBE: "
        + json.dumps(
            {
                "temporary_f": W,
                "attempts": len(l12_attempts),
                "aligned_rows": len(l12_rows),
                "q_max": CLASS_Q_MAX,
            }
        ),
        "MEMBER ATTEMPTS: " + json.dumps(member_records),
    ]
    if closure is None:
        lines.append(
            "RESULT: OPEN after the bounded search; every attempt is recorded in the JSONL artifact."
        )
    else:
        lines.extend(
            (
                "RESULT: VERIFIED closure row " + json.dumps(closure),
                "L11 CLASS: "
                + json.dumps(
                    {
                        "a": replay["a"],
                        "eps": closure["eps"],
                        "q1": closure["q1"],
                        "N": closure["N"],
                        "S": replay["S"],
                        "ks": replay["ks"],
                        "syms": replay["class_syms"],
                        "Q0": replay["Q0"],
                    }
                ),
                "KERNEL: "
                + json.dumps(
                    {
                        "b": replay["b"],
                        "tau": replay["tau"],
                        "A": replay["A"],
                        "delta": replay["delta"],
                        "alpha": replay["alpha"],
                        "D": replay["D"],
                        "c0": replay["c0"],
                        "M0": replay["M0"],
                        "x0": replay["x0"],
                        "d0": replay["d0"],
                        "P_degree": replay["P_degree"],
                        "kernel_identity": replay["kernel_identity"],
                        "aligned_syms": replay["aligned_syms"],
                    }
                ),
                "L13: "
                + json.dumps(
                    {
                        "small_emergent": replay["small_emergent"],
                        "rn": replay["stripped_rn"],
                        "rd": replay["stripped_rd"],
                        "verdict": replay["emergent_verdict"],
                        "cofactor_info": replay["cofactor_info"],
                        "cofactor_prime_proven": replay["cofactor_prime_proven"],
                    }
                ),
                "CROSS-CHECKS: "
                + json.dumps(
                    {
                        "tied": replay["tied"],
                        "ramified_empty": replay["ramified_empty"],
                        "closure_basis": "L11 class alignment plus L13 zero-bad ladder",
                    }
                ),
            )
        )
    lines.extend(
        (
            "ARTIFACT: " + str(OUT),
            "REPLAY: python3 math/h10q/l17_horizon137.py (from workspace or set H10Q_ROOT)",
        )
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    h10q._selftest()
    assert Z == F(-137)
    assert h10q._is_prime(W)
    assert h10q.factorint(abs(Z.numerator)) == {137: 1}
    assert h10q._l11_escape_primes(Z) == [137]
    assert CELL not in h10q._L11_CLASSES

    roots, a_pool, l11_attempts, l11_rows = discover_l11()
    l12_attempts, l12_rows = discover_l12()
    assert roots == [50, 87]
    if A_T_MAX == 2:
        assert a_pool == [87, 187, 361]

    rows = l11_rows + l12_rows
    rows.sort(
        key=lambda row: (
            0
            if (
                row["family"] == "L11-fresh"
                and row["a"] == "87"
                and row["eps"] == 1
                and row["q1"] == 149
            )
            else 1,
            int(row["q1"]),
            abs(int(row["a"])),
            int(row["eps"]),
            row["family"],
        )
    )

    closure = None
    replay = None
    member_records = []
    for row in rows:
        deadline = time.time() + MEMBER_SECONDS
        for k in range(MAX_K + 1):
            time.sleep(0.1)
            result, hit = member_decide(row, k, deadline)
            member_records.append(
                {
                    "type": "member",
                    "family": row["family"],
                    "cell": [W, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    "N": str(row["cert"]["N"]),
                    **result,
                }
            )
            if hit is not None:
                closure, replay = close_certificate(row, hit)
                break
            if time.time() >= deadline:
                break
        if closure is not None:
            break

    write_artifact(
        roots,
        a_pool,
        l11_attempts,
        l12_attempts,
        rows,
        member_records,
        closure,
        replay,
    )
    write_report(
        roots,
        a_pool,
        l11_attempts,
        l12_attempts,
        l11_rows,
        l12_rows,
        member_records,
        closure,
        replay,
    )

    records = [json.loads(line) for line in OUT.read_text(encoding="utf-8").splitlines()]
    assert records[0]["type"] == "meta"
    assert records[-1]["type"] == "summary"
    if closure is None:
        assert records[-1]["status"] == "OPEN"
        raise RuntimeError("no proved L13 zero member found; see the full attempt log")

    assert closure == {
        "family": "L11",
        "cell": [137, [-1, 1]],
        "a": "87",
        "eps": 1,
        "f": 1,
        "q1": 149,
        "N": "8850572640",
        "k_zero": 0,
        "Q": "149",
        "rung": "prime",
        "tied": "True",
        "ramified_empty": True,
        "source": "l17_horizon137.py",
    }
    assert records[-1]["status"] == "VERIFIED"
    assert any(record == closure for record in records)
    assert REPORT.is_file() and "STATUS: VERIFIED" in REPORT.read_text(encoding="utf-8")
    print(json.dumps(closure, sort_keys=True))
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
