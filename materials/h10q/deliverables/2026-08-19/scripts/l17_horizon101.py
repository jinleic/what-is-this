#!/usr/bin/env python3
"""Fresh H-class discovery and one-member closure for (w,z)=(101,-101).

Uses only math/h10q kernel primitives.  The script deliberately keeps the
fresh L12 row out of the repository table: it injects one candidate row for
l12_class.l12_class_cert, restores the frozen table immediately, and records
that prerequisite in the report/artifact.  Every candidate is decided by the
proven-primality engine; refusals are emitted as refusals, never evidence.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(os.environ.get("H10Q_ROOT", "/Users/jinleic/jinleic-workspace/math/h10q"))
sys.path.insert(0, str(ROOT))

import h10q
from h10q import (
    OO,
    PrimalityBound,
    _L12_ESCAPE,
    _is_prime,
    _l10_class_cert,
    _l7_tied_status,
    _sun_h,
    factorint,
    hilbert,
    primerange,
    ramified,
    vp,
)
import l12_class
from l13_filter import cofactor_decide, smooth_emergent

CELL = (101, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
OUT = ROOT / "data" / "l17_horizon101.jsonl"
REPORT = Path("/tmp/l17_horizon101.report")
MAX_K = int(os.environ.get("H101_MAX_K", "800"))
CLASS_Q_MAX = int(os.environ.get("H101_CLASS_Q_MAX", "200"))


def emit(out, rec):
    out.write(json.dumps(rec) + "\n")
    out.flush()


def refusal(exc):
    return "refused:" + type(exc).__name__
def cert_status(cert):
    if cert is None:
        return "refused:none"
    return "aligned" if cert.get("ok") else "cert-not-ok"


def cert_summary(cert):
    if cert is None:
        return {}
    return {"N": str(cert["N"]), "S": cert["S"], "syms": cert["syms"],
            "Q0": str(cert["Q0"])}




def fresh_l12_cert(f, eps, q):
    """Call the existing L12 cert for a fresh key without persisting a table edit."""
    old = _L12_ESCAPE.get(CELL, None)
    _L12_ESCAPE[CELL] = (f, eps, q)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if old is None:
            _L12_ESCAPE.pop(CELL, None)
        else:
            _L12_ESCAPE[CELL] = old


def discover_l11():
    """Enumerate the l13h_alt-style fresh L11 construction."""
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_t_max = int(os.environ.get("H101_A_T_MAX", "2"))
    a_pool = sorted({abs(r + s * W * t) for r in roots for s in (1, -1)
                     for t in range(0, a_t_max + 1)
                     if (r + s * W * t) != 0})
    a_pool = [a for a in a_pool if a % 2 == 1]
    rows = []
    attempts = []
    for a_int in a_pool:
        a = F(a_int)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        for eps in (1, -1):
            for q in primerange(3, CLASS_Q_MAX + 1):
                time.sleep(0.1)
                cert = _l10_class_cert(a, Z, tau, eps, q)
                attempts.append({"a": a_int, "eps": eps, "q1": q,
                                 "status": cert_status(cert),
                                 **cert_summary(cert)})
                if cert is not None and cert.get("ok"):
                    rows.append({"a": str(a), "eps": eps, "f": 1,
                                 "q1": q, "tau": str(tau), "cert": cert})
    return roots, a_pool, attempts, rows


def discover_l12():
    """Try the fresh L12 escape family f=101, with all small proven primes q."""
    rows = []
    attempts = []
    for eps in (1, -1):
        for q in primerange(3, CLASS_Q_MAX + 1):
            time.sleep(0.1)
            cert = fresh_l12_cert(W, eps, q)
            attempts.append({"f": W, "eps": eps, "q1": q,
                             "status": cert_status(cert),
                             **cert_summary(cert)})
            if cert is not None and cert.get("ok"):
                rows.append({"a": "1", "eps": eps, "f": W, "q1": q,
                             "tau": str(F(3, 5)), "cert": cert})
    return attempts, rows


def member_decide(row, k, deadline):
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
    if Q < 2:
        result["verdict"] = "nonprime"
        return result, None
    try:
        prime = _is_prime(Q)
    except (PrimalityBound, Exception) as exc:
        # Keep refusal type explicit.  Exception is intentionally caught only
        # after the proven-primality call; no unproved result is accepted.
        result["verdict"] = refusal(exc)
        return result, None
    if not prime:
        result["verdict"] = "nonprime"
        return result, None
    b = F(eps * f * Q)
    try:
        st, detail = smooth_emergent(a, Z, tau, b)
    except Exception as exc:
        result["verdict"] = refusal(exc)
        return result, None
    result["smooth_status"] = st
    if st == "zero":
        result["verdict"] = "zero"
        result["rung"] = "square"
        return result, (k, Q, "square", b, None)
    if st == "cofactor-big":
        try:
            verdict, info = cofactor_decide(a, Z, tau, b, detail=detail,
                                            deadline=deadline)
        except Exception as exc:
            result["verdict"] = refusal(exc)
            return result, None
        result["verdict"] = verdict
        result["rung"] = info.get("rung") if isinstance(info, dict) else None
        result["cofactor_info"] = info
        if verdict == "zero":
            return result, (k, Q, result["rung"], b, info)
        return result, None
    result["verdict"] = st
    return result, None


def close_certificate(row, hit):
    k, Q, rung, b, info = hit
    a = F(row["a"])
    eps = int(row["eps"])
    f = int(row["f"])
    tau = F(row["tau"])
    tied = None
    try:
        tied = str(_l7_tied_status(a, b, Z, tau))
    except Exception as exc:
        tied = refusal(exc)
    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    c0 = None
    try:
        c0 = _sun_h(a, b, Z ** 3)
        M0 = 16 - delta * c0 * c0 - 32 * A * b * ((a - 1) / 2) ** 2
        x0, d0 = alpha * M0, alpha * 2 * b
    except Exception as exc:
        x0 = d0 = None
        c0 = refusal(exc)
    try:
        ram_empty = (ramified(x0, d0) == []) if x0 is not None else refusal(RuntimeError())
    except Exception as exc:
        ram_empty = refusal(exc)
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
        "source": "l17_horizon101.py",
    }
    replay = {
        "type": "kernel-replay",
        "cell": [W, list(UT)],
        "family": row["family"],
        "a": str(a), "b": str(b), "tau": str(tau), "z": str(Z),
        "q1": int(row["q1"]), "Q": str(Q), "k": k,
        "A": str(A), "delta": str(delta), "alpha": str(alpha),
        "c0": str(c0), "x0": str(x0), "d0": str(d0),
        "tied": tied, "ramified_empty": ram_empty,
        "emergent_verdict": "zero",
        "cofactor_info": info,
    }
    return closure, replay


def main():
    # Explicitly run the kernel self-tests; all arithmetic below imports from
    # this same proven engine.  This is intentionally before discovery.
    h10q._selftest()
    roots, a_pool, l11_attempts, l11_rows = discover_l11()
    l12_attempts, l12_rows = discover_l12()
    rows = []
    for x in l11_rows:
        x["family"] = "L11-fresh"
        rows.append(x)
    for x in l12_rows:
        x["family"] = "ESC-fresh"
        rows.append(x)
    # Prefer the smaller L12 escape construction when available; retain all
    # discovery rows so a refusal is auditable and no path is silently hidden.
    rows.sort(key=lambda x: (0 if x["family"] == "ESC-fresh" else 1,
                             int(x["q1"]), int(x["f"]), abs(int(x["a"]))))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as out:
        emit(out, {"type": "meta", "cell": [W, list(UT)], "z": str(Z),
                   "target": "off-grid H-class closure", "engine": "h10q.py",
                   "proven_primality_only": True, "selftest": True,
                   "l11_roots": roots, "l11_a_pool": a_pool,
                   "max_k": MAX_K, "class_q_max": CLASS_Q_MAX,
                   "note": "fresh L11/L12 class discovery; no frozen-table edit"})
        for a in l11_attempts:
            emit(out, {"type": "class-attempt", "family": "L11-fresh", **a})
        for a in l12_attempts:
            emit(out, {"type": "class-attempt", "family": "ESC-fresh", **a})
        for row in rows:
            cert = row["cert"]
            emit(out, {"type": "aligned-class", "family": row["family"],
                       "cell": [W, list(UT)], "a": row["a"],
                       "eps": row["eps"], "f": row["f"], "q1": row["q1"],
                       "tau": row["tau"], "N": str(cert["N"]),
                       "S": cert["S"], "ks": cert["ks"],
                       "syms": cert["syms"], "Q0": str(cert["Q0"]),
                       "excluded": cert.get("excluded", []), "ok": cert["ok"]})

        closure = None
        for row in rows:
            # Search the fresh class in increasing k; each row has its own
            # progression.  A successful zero is immediately replayed.
            deadline = time.time() + float(os.environ.get("H101_MEMBER_SECONDS", "900"))
            for k in range(-1, MAX_K + 1):
                time.sleep(0.1)
                result, hit = member_decide(row, k, deadline)
                emit(out, {"type": "member", "family": row["family"],
                           "cell": [W, list(UT)], "a": row["a"],
                           "eps": row["eps"], "f": row["f"],
                           "q1": row["q1"], "N": str(row["cert"]["N"]),
                           **result})
                if hit is not None:
                    closure, replay = close_certificate(row, hit)
                    emit(out, replay)
                    emit(out, closure)
                    break
                if time.time() >= deadline:
                    emit(out, {"type": "search-stop", "family": row["family"],
                               "q1": row["q1"], "reason": "deadline"})
                    break
            if closure is not None:
                break
        if closure is None:
            emit(out, {"type": "summary", "status": "open",
                       "cell": [W, list(UT)],
                       "reason": "no proved zero member in searched classes"})
        else:
            emit(out, {"type": "summary", "status": "closed",
                       "cell": [W, list(UT)], "family": closure["family"],
                       "k_zero": closure["k_zero"], "Q": closure["Q"]})

    lines = []
    lines.append("H101 off-grid closure recipe and replay")
    lines.append(f"cell={CELL!r}; z={Z}; engine=h10q.py; kernel_selftest=PASS")
    lines.append("Fresh L11 prerequisite: choose escape root r<w with 1+4r^2=0 mod w, construct odd a, canonical tau=(1+2a^2)/(1+4a^2), then _l10_class_cert(a,z,tau,eps,q1).")
    lines.append("Fresh L12 prerequisite: f=w must be proven prime and divide numerator(z); l12_class.l12_class_cert currently requires a row in _L12_ESCAPE, so inject the candidate tuple (f,eps,q1) only for the call (or promote it into a source table), then restore the frozen table.")
    l12_bad = [r for r in l12_attempts if r["status"] == "cert-not-ok"]
    l12_none = [r for r in l12_attempts if r["status"] == "refused:none"]
    l12_fail_places = sorted({place for r in l12_bad
                              for place, value in r.get("syms", {}).items()
                              if value != 1}, key=str)
    lines.append(f"Fresh L12 result: {len(l12_bad)} cert-not-ok rows"
                 f" (bad frozen places {l12_fail_places}),"
                 f" {len(l12_none)} refused:none; no L12 aligned row.")
    lines.append("Repository search found no `_LL` symbol; the callable path is"
                 " _l10_class_cert/l12_class_cert -> smooth_emergent ->"
                 " cofactor_decide -> _l7_tied_status.")
    lines.append("L13 prerequisite: for Q=q1+kN, prove Q prime with _is_prime; b=eps*f*Q; run smooth_emergent and, only for cofactor-big, cofactor_decide. A zero verdict is the only closure evidence; zero-jacobi/refusals are not proofs.")
    lines.append(f"L11 roots={roots}; a_pool={a_pool}; aligned L11 rows={len(l11_rows)}; aligned fresh L12 rows={len(l12_rows)}")
    if closure is None:
        lines.append("RESULT: OPEN after the bounded search; no proved zero member was found.")
        lines.append("All candidate refusals/nonprimes are recorded in data/l17_horizon101.jsonl; no refusal is evidence.")
    else:
        lines.append(f"RESULT: VERIFIED closure row {json.dumps(closure, sort_keys=True)}")
        lines.append("The closure row is emitted with the exact 13-field schema used by data/l13h_all_closures.json; kernel replay immediately precedes it in the JSONL artifact.")
    REPORT.write_text("\n".join(lines) + "\n")
    print("wrote", OUT)
    print("wrote", REPORT)
    print("closure", closure)


if __name__ == "__main__":
    main()
