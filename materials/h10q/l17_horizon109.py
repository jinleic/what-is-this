#!/usr/bin/env python3
"""Fresh L11/L12 discovery and L13 closure replay for (109,-109).

The checked-in class tables stop at w<100.  This script discovers the
off-grid L11 class from the root of 1+4*a^2 modulo 109, probes the fresh L12
escape path without persisting a table edit, and replays one member through
the proven-primality L13 cofactor ladder.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path


def find_root() -> Path:
    env = os.environ.get("H10Q_ROOT")
    if env:
        p = Path(env).expanduser().resolve()
        if (p / "h10q.py").exists():
            return p
    here = Path(__file__).resolve().parent
    candidates = [
        Path.cwd() / "math" / "h10q",
        here.parent / "math" / "h10q",
        Path("/Users/jinleic/jinleic-workspace/math/h10q"),
    ]
    for p in candidates:
        if (p / "h10q.py").exists():
            return p.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402
import l12_class  # noqa: E402
import l13_filter  # noqa: E402

CELL = (109, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
CLASS_Q_MAX = int(os.environ.get("H109_CLASS_Q_MAX", "200"))
MAX_K = int(os.environ.get("H109_MAX_K", "0"))


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def cert_summary(cert):
    if cert is None:
        return {}
    return {
        "N": str(cert["N"]),
        "S": cert["S"],
        "syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "excluded": cert.get("excluded", []),
    }


def fresh_l12_cert(eps: int, q: int):
    """Run the table-backed L12 cert with one ephemeral off-grid row."""
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(CELL, sentinel)
    l12_class._L12_ESCAPE[CELL] = (W, eps, q)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(CELL, None)
        else:
            l12_class._L12_ESCAPE[CELL] = previous


def discover_l11():
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_t_max = int(os.environ.get("H109_A_T_MAX", "2"))
    a_pool = sorted(
        {
            abs(r + s * W * t)
            for r in roots
            for s in (1, -1)
            for t in range(a_t_max + 1)
            if r + s * W * t != 0
        }
    )
    a_pool = [a for a in a_pool if a % 2 == 1]
    rows = []
    attempts = []
    for a_int in a_pool:
        a = F(a_int)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        for eps in (1, -1):
            for q in h10q.primerange(3, CLASS_Q_MAX + 1):
                time.sleep(0.1)
                try:
                    cert = h10q._l10_class_cert(a, Z, tau, eps, q)
                    status = (
                        "aligned" if cert is not None and cert.get("ok")
                        else "cert-not-ok" if cert is not None else "refused:none"
                    )
                except Exception as exc:
                    cert = None
                    status = refusal(exc)
                rec = {"a": a_int, "eps": eps, "q1": q, "status": status}
                rec.update(cert_summary(cert))
                attempts.append(rec)
                if cert is not None and cert.get("ok"):
                    rows.append(
                        {
                            "family": "L11-fresh",
                            "a": str(a),
                            "eps": eps,
                            "f": 1,
                            "q1": q,
                            "tau": str(tau),
                            "cert": cert,
                        }
                    )
    return roots, a_pool, attempts, rows


def discover_l12():
    attempts = []
    rows = []
    for eps in (1, -1):
        for q in h10q.primerange(3, CLASS_Q_MAX + 1):
            time.sleep(0.1)
            try:
                cert = fresh_l12_cert(eps, q)
                status = (
                    "aligned" if cert is not None and cert.get("ok")
                    else "cert-not-ok" if cert is not None else "refused:none"
                )
            except Exception as exc:
                cert = None
                status = refusal(exc)
            rec = {"f": W, "eps": eps, "q1": q, "status": status}
            rec.update(cert_summary(cert))
            attempts.append(rec)
            if cert is not None and cert.get("ok"):
                rows.append(
                    {
                        "family": "ESC-fresh",
                        "a": "1",
                        "eps": eps,
                        "f": W,
                        "q1": q,
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
    k, Q, rung, b, info = hit
    a = F(row["a"])
    eps = int(row["eps"])
    f = int(row["f"])
    tau = F(row["tau"])
    tied = None
    try:
        tied = str(h10q._l7_tied_status(a, b, Z, tau))
    except Exception as exc:
        tied = refusal(exc)
    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    try:
        c0 = h10q._sun_h(a, b, Z**3)
        M0 = 16 - delta * c0 * c0 - 32 * A * b * ((a - 1) / 2) ** 2
        x0, d0 = alpha * M0, alpha * 2 * b
    except Exception as exc:
        c0 = refusal(exc)
        M0 = None
        x0 = d0 = None
    try:
        ram_empty = ramified = (h10q.ramified(x0, d0) == []) if x0 is not None else refusal(RuntimeError())
    except Exception as exc:
        ramified = refusal(exc)
        ram_empty = ramified
    closure = {
        "family": "L11" if row["family"] == "L11-fresh" else "ESC",
        "cell": [W, list(UT)],
        "a": str(a),
        "eps": eps,
        "f": f,
        "q1": int(row["q1"]),
        "N": str(row["cert"]["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": tied,
        "ramified_empty": ram_empty,
        "source": "l17_horizon109.py",
    }
    replay = {
        "type": "kernel-replay",
        "cell": [W, list(UT)],
        "family": row["family"],
        "a": str(a),
        "b": str(b),
        "tau": str(tau),
        "z": str(Z),
        "q1": int(row["q1"]),
        "Q": str(Q),
        "k": k,
        "A": str(A),
        "delta": str(delta),
        "alpha": str(alpha),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "tied": tied,
        "ramified": ramified,
        "ramified_empty": ram_empty,
        "emergent_verdict": "zero",
        "cofactor_info": info,
    }
    return closure, replay


def main():
    # Explicit self-tests are required proof hygiene for every replay.
    h10q._selftest()
    # The target is off-grid but has a numerator escape prime 109 == 1 (mod 4),
    # so the fresh L11 route is the relevant branch.
    assert Z == F(-109)
    assert h10q.factorint(abs(Z.numerator)) == {109: 1}
    assert h10q._l11_escape_primes(Z) == [109]
    assert CELL not in h10q._L11_CLASSES
    roots, a_pool, l11_attempts, l11_rows = discover_l11()
    l12_attempts, l12_rows = discover_l12()
    assert roots == [38, 71] and 71 in a_pool
    assert not l12_rows
    rows = l11_rows + l12_rows
    assert any(
        r["family"] == "L11-fresh" and r["a"] == "71"
        and r["eps"] == 1 and r["q1"] == 193 and r["cert"]["ok"]
        for r in rows
    )
    # The first known zero is q1=193 on the minimal odd root a=71.  Sorting
    # this row first avoids spending a factor budget on earlier large-
    # cofactor classes while retaining every discovery attempt in the JSONL.
    rows.sort(
        key=lambda r: (
            0
            if (r["family"] == "L11-fresh" and r["a"] == "71"
                and r["eps"] == 1 and r["q1"] == 193)
            else 1,
            int(r["q1"]),
            abs(int(r["a"])),
            int(r["eps"]),
        )
    )
    out_path = ROOT / "data" / "l17_horizon109.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    closure = None
    replay = None
    member_records = []
    with out_path.open("w", encoding="utf-8") as out:
        meta = {
            "type": "meta",
            "schema": "l17-off-grid-closure-v1",
            "task": "Off-grid H-class closure",
            "engine": "math/h10q/h10q.py proven-primality only",
            "cell": [W, list(UT)],
            "family": "L11",
            "verified": True,
            "replay": "full in-script",
            "source_script": "math/h10q/l17_horizon109.py",
            "l11_roots": roots,
            "l11_a_pool": a_pool,
            "class_q_max": CLASS_Q_MAX,
            "max_k": MAX_K,
        }
        out.write(json.dumps(meta, sort_keys=True) + "\n")
        for rec in l11_attempts:
            out.write(json.dumps({"type": "class-attempt", "family": "L11-fresh", **rec}) + "\n")
        for rec in l12_attempts:
            out.write(json.dumps({"type": "class-attempt", "family": "ESC-fresh", **rec}) + "\n")
        for row in rows:
            cert = row["cert"]
            out.write(json.dumps({
                "type": "aligned-class",
                "family": row["family"],
                "cell": [W, list(UT)],
                "a": row["a"], "eps": row["eps"], "f": row["f"],
                "q1": row["q1"], "tau": row["tau"], "N": str(cert["N"]),
                "S": cert["S"], "ks": cert["ks"], "syms": cert["syms"],
                "Q0": str(cert["Q0"]), "excluded": cert.get("excluded", []),
                "ok": cert["ok"],
            }) + "\n")
        for row in rows:
            deadline = time.time() + float(os.environ.get("H109_MEMBER_SECONDS", "900"))
            for k in range(MAX_K + 1):
                time.sleep(0.1)
                result, hit = member_decide(row, k, deadline)
                member_record = {
                    "type": "member",
                    "family": row["family"],
                    "cell": [W, list(UT)],
                    "a": row["a"], "eps": row["eps"], "f": row["f"],
                    "q1": row["q1"], "N": str(row["cert"]["N"]),
                    **result,
                }
                out.write(json.dumps(member_record, sort_keys=True) + "\n")
                member_records.append(member_record)
                if hit is not None:
                    closure, replay = close_certificate(row, hit)
                    out.write(json.dumps(replay, sort_keys=True) + "\n")
                    out.write(json.dumps(closure, sort_keys=True) + "\n")
                    break
                if time.time() >= deadline:
                    break
            if closure is not None:
                break
        if closure is None:
            out.write(json.dumps({
                "type": "summary", "status": "OPEN",
                "cell": [W, list(UT)],
                "reason": "no proved zero member in searched classes",
            }, sort_keys=True) + "\n")
        else:
            out.write(json.dumps({
                "type": "summary", "status": "VERIFIED",
                "cell": [W, list(UT)], "family": closure["family"],
                "k_zero": closure["k_zero"], "Q": closure["Q"],
            }, sort_keys=True) + "\n")
    if closure is None or replay is None:
        raise RuntimeError("no L13 zero-bad member found")
    assert closure["family"] == "L11"
    assert closure["cell"] == [109, [-1, 1]]
    assert closure["a"] == "71" and closure["eps"] == 1
    assert closure["f"] == 1 and closure["q1"] == 193
    assert closure["N"] == "962112480"
    assert closure["k_zero"] == 0 and closure["Q"] == "193"
    assert closure["rung"] == "prime"

    # Independently reconstruct the decisive kernel/L13 facts for the report.
    row = next(r for r in rows if r["family"] == "L11-fresh" and r["a"] == "71" and r["eps"] == 1 and r["q1"] == 193)
    cert = row["cert"]
    a = F(row["a"])
    tau = F(row["tau"])
    Q = int(closure["Q"])
    b = F(closure["eps"] * closure["f"] * Q)
    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    c0 = h10q._sun_h(a, b, Z**3)
    M0 = 16 - delta * c0 * c0 - 32 * A * b * ((a - 1) / 2) ** 2
    D = 1 - Z**3 - a * a * Z**6
    P = h10q._l10_P(a, Z**3, D, A, delta, (a - 1) / 2)
    Pb = sum(c * b**i for i, c in enumerate(P))
    assert M0 == Pb / (b**4 * D * D * A * A)
    member_syms = l13_filter.aligned(a, Z, tau, b)
    assert member_syms is not None and all(v == 1 for v in member_syms.values())
    info = next(m.get("cofactor_info") for m in member_records if m.get("verdict") == "zero")
    cofactor_prime = int(info["places"][0][0])
    assert h10q._is_prime(cofactor_prime)
    assert info["places"][0][1] == 1
    diagnostics = {
        "status": "VERIFIED",
        "target_z": str(Z),
        "l11": {
            "z_factorization": h10q.factorint(abs(Z.numerator)),
            "escape_primes_1_mod_4": h10q._l11_escape_primes(Z),
            "table_row_present": CELL in h10q._L11_CLASSES,
            "roots": roots,
            "a_pool": a_pool,
            "aligned_rows": len(l11_rows),
        },
        "l12": {
            "f": W,
            "aligned_rows": len(l12_rows),
            "attempts": len(l12_attempts),
            "note": "fresh table injection only; no aligned off-grid L12 row",
        },
        "class": {
            "a": str(a), "eps": closure["eps"], "f": closure["f"],
            "q1": closure["q1"], "N": closure["N"], "S": cert["S"],
            "ks": cert["ks"], "syms": cert["syms"],
            "Q0": str(cert["Q0"]), "excluded": cert.get("excluded", []),
        },
        "kernel": {
            "b": str(b), "tau": str(tau), "A": str(A),
            "delta": str(delta), "alpha": str(alpha), "D": str(D),
            "c0": str(c0), "M0": str(M0),
            "x0": replay["x0"], "d0": replay["d0"],
            "P_degree": len(P) - 1, "P_identity": True,
            "aligned_symbols": member_syms,
        },
        "l13": {
            "smooth_status": next(m["smooth_status"] for m in member_records if m.get("verdict") == "zero"),
            "small_emergent": [], "rn": info["rn"], "rd": info["rd"],
            "cofactor_verdict": "zero", "cofactor_rung": info["rung"],
            "cofactor_prime": str(cofactor_prime),
            "cofactor_prime_proven": True,
            "cofactor_info": info,
        },
        "cross_checks": {
            "tied": closure["tied"],
            "ramified_empty": closure["ramified_empty"],
            "refusals_are_not_evidence": True,
            "closure_basis": "L11 class alignment + L13 zero-bad ladder",
        },
        "closure_row": closure,
    }
    report_path = Path(os.environ.get("HORIZON109_REPORT", "/tmp/l17_horizon109.report"))
    report_lines = [
        "STATUS: VERIFIED",
        "TARGET: cell [109,[-1,1]], z=-109, family=L11",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        "DISCOVERY: " + json.dumps(diagnostics["l11"], sort_keys=True),
        "L12 PROBE: " + json.dumps(diagnostics["l12"], sort_keys=True),
        "L11 CLASS: " + json.dumps(diagnostics["class"]),
        "KERNEL: " + json.dumps(diagnostics["kernel"]),
        "L13: " + json.dumps(diagnostics["l13"], sort_keys=True),
        "CROSS-CHECKS: " + json.dumps(diagnostics["cross_checks"], sort_keys=True),
        "CLOSURE ROW: " + json.dumps(closure, sort_keys=True),
        "ARTIFACT: " + str(out_path),
        "REPLAY: python math/h10q/l17_horizon109.py (from workspace or with H10Q_ROOT)",
        "NOTE: tied-status and ramified cross-checks may refuse FactorBudget for this large-a witness; the verified closure is the L13 zero-bad certificate, and no refusal is used as evidence.",
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps(closure, sort_keys=True))
    print("WROTE", out_path)
    print("WROTE", report_path)


if __name__ == "__main__":
    main()
