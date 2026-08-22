#!/usr/bin/env python3
"""Fresh ESC-class discovery and one-member closure for (w,z)=(127,-127).

The target is off the frozen w<=97 grid.  Because 127 is 3 modulo 4, the L11
root condition has no solution; this replay therefore uses the ESC family with
f=127.  The L12 escape tuple is installed only around each certificate call,
and the checked-in tables are never modified.  Every primality and
factorization decision is delegated to math/h10q/h10q.py.
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
    for path in candidates:
        if (path / "h10q.py").is_file():
            return path.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l12_class  # noqa: E402
import l13_filter  # noqa: E402


CELL = (127, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
TAU = F(3, 5)
F_ESC = W
OUT = ROOT / "data" / "l17_horizon127.jsonl"
REPORT = Path(os.environ.get("HORIZON127_REPORT", "/tmp/l17_horizon127.report"))
CLASS_Q_MAX = int(os.environ.get("H127_CLASS_Q_MAX", "149"))
MAX_K = int(os.environ.get("H127_MAX_K", "0"))


def emit(out, record):
    """Write and immediately flush one replay record."""
    out.write(json.dumps(record) + "\n")
    out.flush()


def refusal(exc):
    """Record an engine refusal without treating it as mathematical evidence."""
    return "refused:" + type(exc).__name__


def fresh_l12_cert(eps, q):
    """Run canonical l12_class_cert with one temporary off-grid ESC row."""
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(CELL, sentinel)
    l12_class._L12_ESCAPE[CELL] = (F_ESC, eps, q)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(CELL, None)
        else:
            l12_class._L12_ESCAPE[CELL] = previous


def discover_l11():
    """Exhaust the L11 congruence and its canonical odd-a lift pool."""
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_pool = sorted({
        abs(r + sign * W * shift)
        for r in roots
        for sign in (1, -1)
        for shift in range(0, 3)
        if r + sign * W * shift != 0
    })
    return roots, [a for a in a_pool if a % 2 == 1]


def discover_l12():
    """Try the fresh ESC family f=w for every proven prime in the bound."""
    attempts = []
    rows = []
    for eps in (1, -1):
        for q in h10q.primerange(3, CLASS_Q_MAX + 1):
            time.sleep(0.1)
            try:
                cert = fresh_l12_cert(eps, q)
            except Exception as exc:
                attempts.append({
                    "f": F_ESC,
                    "eps": eps,
                    "q1": q,
                    "status": refusal(exc),
                })
                continue
            ok = cert is not None and cert.get("ok")
            attempts.append({
                "f": F_ESC,
                "eps": eps,
                "q1": q,
                "status": "aligned" if ok else "cert-not-ok",
            })
            if ok:
                rows.append({
                    "family": "ESC-fresh",
                    "a": "1",
                    "eps": eps,
                    "f": F_ESC,
                    "q1": q,
                    "tau": str(TAU),
                    "cert": cert,
                })
    rows.sort(key=lambda row: (
        int(row["q1"]),
        0 if int(row["eps"]) == 1 else 1,
    ))
    return attempts, rows


def member_decide(row, k, deadline):
    """Decide one prime-ladder member through the canonical L13 ladder."""
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
        result["verdict"] = "zero"
        result["rung"] = "square"
        return result, (k, Q, "square", b, None)
    if smooth_status != "cofactor-big":
        result["verdict"] = smooth_status
        return result, None

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


def close_certificate(row, hit):
    """Re-run the L12, kernel, and L13 checks before emitting a closure."""
    k, Q, rung, b, _prior_info = hit
    a = F(row["a"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])
    tau = F(row["tau"])
    cert = row["cert"]

    replayed_cert = fresh_l12_cert(eps, q1)
    assert replayed_cert is not None and replayed_cert.get("ok")
    assert int(replayed_cert["N"]) == int(cert["N"])
    assert replayed_cert["syms"] == cert["syms"]
    assert Q == q1 + k * int(cert["N"])
    assert h10q._is_prime(Q)
    assert h10q._is_prime(f)

    A0 = 1 + 4 * a * a
    delta = 1 - A0 * tau * tau
    alpha = -delta * A0
    Z3 = Z ** 3
    D = 1 - Z3 - a * a * Z3 * Z3
    s = (a - 1) / 2
    c0 = h10q._sun_h(a, b, Z3)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A0 * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    polynomial = h10q._l10_P(a, Z3, D, A0, delta, s)
    evaluated = sum(coefficient * b ** degree
                    for degree, coefficient in enumerate(polynomial))
    assert M0 == evaluated / (b ** 4 * D * D * A0 * A0)

    member_syms = h10q._l12_syms(a, Z, tau, b)
    assert member_syms is not None and all(
        value == 1 for value in member_syms.values()
    )
    aligned_syms = l13_filter.aligned(a, Z, tau, b)
    assert aligned_syms is not None and all(
        value == 1 for value in aligned_syms.values()
    )

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, Z, tau, b)
    assert smooth_status == "cofactor-big"
    emergent, (_rn, rd), smooth = smooth_detail
    assert emergent == [] and not smooth and rd == 1
    verdict, info = l13_filter.cofactor_decide(
        a, Z, tau, b, detail=smooth_detail, deadline=time.time() + 600
    )
    assert verdict == "zero" and info["rung"] == rung
    assert info["places"] and all(
        place[1] == 1 and place[2] == "proved" for place in info["places"]
    )
    if info["rung"] == "prime":
        assert info["parts"][0]["kind"] == "prime"
        assert h10q._is_prime(int(info["places"][0][0]))

    cofactor_polynomial = h10q._l13_make_P(Z, "sq")
    irreducibility_prime = None
    for prime in h10q.primerange(7, 1000):
        time.sleep(0.1)
        if prime == 5 or cofactor_polynomial[-1] % prime == 0:
            continue
        if h10q._l13_irred8(cofactor_polynomial, prime):
            irreducibility_prime = prime
            break
    assert len(cofactor_polynomial) - 1 == 8
    assert irreducibility_prime == 17

    tied = h10q._l7_tied_status(a, b, Z, tau)
    assert tied is True
    ramified = h10q.ramified(x0, d0)
    assert ramified == []

    closure = {
        "family": "ESC",
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
        "source": "l17_horizon127.py",
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
        "A": str(A0),
        "delta": str(delta),
        "alpha": str(alpha),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
        "l12_replayed": True,
        "member_syms": member_syms,
        "aligned_syms": aligned_syms,
        "tied": tied,
        "ramified_empty": not ramified,
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "small_emergent": emergent,
        "cofactor_rung": info["rung"],
        "cofactor_prime_proven": info["rung"] == "prime",
        "P_degree": len(cofactor_polynomial) - 1,
        "P_irreducibility_prime": irreducibility_prime,
        "P_irreducible_over_Q": True,
        "P_coefficients": cofactor_polynomial,
    }
    return closure, replay


def main():
    """Discover the off-grid class, search it, and write a full replay log."""
    h10q._selftest()
    assert Z == F(-127)
    assert W % 4 == 3
    assert h10q._is_prime(W)
    z_factorization = h10q.factorint(abs(Z.numerator))
    l11_escape_primes = h10q._l11_escape_primes(Z)
    roots, a_pool = discover_l11()
    assert z_factorization == {127: 1}
    assert roots == [] and a_pool == []
    assert l11_escape_primes == []
    assert CELL not in h10q._L11_CLASSES

    l12_attempts, l12_rows = discover_l12()
    assert l12_rows

    OUT.parent.mkdir(parents=True, exist_ok=True)
    closure = None
    with OUT.open("w", encoding="utf-8") as out:
        emit(out, {
            "type": "meta",
            "schema": "l17-off-grid-closure-v1",
            "task": "Off-grid H-class closure",
            "engine": "h10q.py",
            "proven_primality_only": True,
            "selftest": True,
            "cell": [W, list(UT)],
            "z": str(Z),
            "family": "ESC",
            "l11_roots": roots,
            "l11_a_pool": a_pool,
            "l11_escape_primes_1_mod_4": l11_escape_primes,
            "structural_finding": (
                "w=127 is 3 mod 4 and the exhaustive L11 root set is empty; "
                "the working route is ESC with f=w"
            ),
            "probe_space": {
                "a": ["1"],
                "eps": [1, -1],
                "f": [F_ESC],
                "q1": {
                    "kind": "proven-primes",
                    "min": 3,
                    "max": CLASS_Q_MAX,
                },
                "k": {"min": 0, "max": MAX_K},
            },
            "stop_policy": "first member whose zero verdict passes the full replay",
            "note": "fresh discovery; each L12 ESC row is temporary",
        })
        for attempt in l12_attempts:
            emit(out, {
                "type": "class-attempt",
                "family": "ESC-fresh",
                "cell": [W, list(UT)],
                **attempt,
            })
        for row in l12_rows:
            cert = row["cert"]
            emit(out, {
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
            })

        for row in l12_rows:
            deadline = time.time() + float(
                os.environ.get("H127_MEMBER_SECONDS", "900")
            )
            for k in range(0, MAX_K + 1):
                time.sleep(0.1)
                result, hit = member_decide(row, k, deadline)
                emit(out, {
                    "type": "member",
                    "family": row["family"],
                    "cell": [W, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    "N": str(row["cert"]["N"]),
                    **result,
                })
                if hit is not None:
                    closure, replay = close_certificate(row, hit)
                    emit(out, replay)
                    emit(out, closure)
                    break
            if closure is not None:
                break

        if closure is None:
            emit(out, {
                "type": "summary",
                "status": "open",
                "cell": [W, list(UT)],
                "reason": "no proved zero member in searched ESC classes",
            })
        else:
            emit(out, {
                "type": "summary",
                "status": "closed",
                "label": "PROVED",
                "cell": [W, list(UT)],
                "family": closure["family"],
                "k_zero": closure["k_zero"],
                "Q": closure["Q"],
            })

    report_lines = [
        "H127 off-grid closure recipe and replay",
        f"STATUS: {'PROVED (VERIFIED replay)' if closure is not None else 'OPEN'}",
        f"TARGET: cell [127,[-1,1]], z={Z}, family=ESC",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        "STRUCTURAL FINDING (PROVED): 127 is 3 modulo 4; exhaustive residue checking "
        "gives no solution to 1+4r^2=0 modulo 127, so the L11 path is unavailable. "
        "The working path is ESC with f=127.",
        f"L11 REPLAY (PROVED): z factorization={json.dumps(z_factorization, sort_keys=True)}; "
        f"escape primes 1 mod 4={l11_escape_primes}; roots={roots}; a_pool={a_pool}",
        f"SEARCH LOG: ESC attempts q<= {CLASS_Q_MAX}; aligned classes={len(l12_rows)}; "
        f"max_k={MAX_K}. Every attempt and member verdict is in the JSONL artifact.",
    ]
    if closure is None:
        report_lines.append(
            "RESULT: OPEN after the bounded search; refusals and nonzero verdicts do not close the cell."
        )
    else:
        report_lines.extend([
            "L12 (PROVED): f=127, temporary off-grid tuple, canonical class alignment replayed.",
            "L13 (PROVED): smooth_emergent and cofactor_decide replayed; zero is the closure verdict.",
            f"RESULT: PROVED VERIFIED closure row {json.dumps(closure, sort_keys=True)}",
            "SCHEMA: the closure row has the exact 13 fields of data/l13h_all_closures.json; "
            "a kernel-replay record immediately precedes it.",
        ])
    report_lines.extend([
        f"ARTIFACT: {OUT}",
        "REPLAY: python3 math/h10q/l17_horizon127.py "
        "(run from repository root or set H10Q_ROOT)",
    ])
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("wrote", OUT)
    print("wrote", REPORT)
    print("closure", closure)


if __name__ == "__main__":
    main()
