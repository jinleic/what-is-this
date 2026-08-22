#!/usr/bin/env python3
"""Fresh H-class discovery and one-member closure for (w,z)=(103,-103).

The target is off the frozen w<=97 grid.  The L12 escape tuple is injected only
for each in-process certificate call; the checked-in tables are not modified.
Every primality/factorization decision is delegated to math/h10q/h10q.py.
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


CELL = (103, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
A = F(1)
TAU = F(3, 5)
F_ESC = 103
OUT = ROOT / "data" / "l17_horizon103.jsonl"
REPORT = Path(os.environ.get("HORIZON103_REPORT", "/tmp/l17_horizon103.report"))
CLASS_Q_MAX = int(os.environ.get("H103_CLASS_Q_MAX", "41"))
MAX_K = int(os.environ.get("H103_MAX_K", "0"))


def emit(out, rec):
    # Hilbert-symbol maps mix integer prime keys with the string ``oo``.
    out.write(json.dumps(rec) + "\n")
    out.flush()


def refusal(exc):
    return "refused:" + type(exc).__name__


def fresh_l12_cert(eps, q):
    """Run canonical l12_class_cert with a temporary off-grid table row."""
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
    """Replay the fresh L11 root/a construction used by l13h_alt."""
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_pool = sorted({abs(r + s * W * t) for r in roots for s in (1, -1)
                     for t in range(0, 3) if (r + s * W * t) != 0})
    a_pool = [a for a in a_pool if a % 2 == 1]
    return roots, a_pool


def discover_l12():
    """Try the fresh escape family f=w with small proven-prime q1 values."""
    rows = []
    attempts = []
    for eps in (1, -1):
        for q in h10q.primerange(3, CLASS_Q_MAX + 1):
            time.sleep(0.1)
            cert = fresh_l12_cert(eps, q)
            ok = cert is not None and cert.get("ok")
            attempts.append({"f": F_ESC, "eps": eps, "q1": q,
                             "status": "ok" if ok else "refused:class-cert"})
            if ok:
                rows.append({"family": "ESC-fresh", "a": "1", "eps": eps,
                             "f": F_ESC, "q1": q, "tau": str(TAU),
                             "cert": cert})
    rows.sort(key=lambda row: (int(row["q1"]),
                               0 if int(row["eps"]) == 1 else 1))
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
    if smooth_status == "cofactor-big":
        try:
            verdict, info = l13_filter.cofactor_decide(
                a, Z, tau, b, detail=detail, deadline=deadline)
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
    """Re-run kernel/L12/L13 checks immediately before emitting the row."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    eps = int(row["eps"])
    f = int(row["f"])
    tau = F(row["tau"])
    assert Q == int(row["q1"]) + k * int(row["cert"]["N"])
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
    P = h10q._l10_P(a, Z3, D, A0, delta, s)
    Pb = sum(c * b ** i for i, c in enumerate(P))
    assert M0 == Pb / (b ** 4 * D * D * A0 * A0)

    member_syms = h10q._l12_syms(a, Z, tau, b)
    assert member_syms is not None and all(v == 1 for v in member_syms.values())
    aligned_syms = l13_filter.aligned(a, Z, tau, b)
    assert aligned_syms is not None and all(v == 1 for v in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, Z, tau, b)
    assert smooth_status == "cofactor-big"
    emergent, (rn, rd), smooth = smooth_detail
    assert emergent == [] and not smooth and rd == 1
    verdict, info = l13_filter.cofactor_decide(
        a, Z, tau, b, detail=smooth_detail, deadline=time.time() + 600)
    assert verdict == "zero" and info["rung"] == rung
    assert info["places"] and all(place[1] == 1 and place[2] == "proved"
                                    for place in info["places"])
    if info["rung"] == "prime":
        assert info["parts"][0]["kind"] == "prime"
        assert h10q._is_prime(int(info["places"][0][0]))

    # L13's degree-eight cofactor object is irreducible over Q by the
    # checked-in Frobenius certificate machinery; p=19 is the first witness.
    Pi = h10q._l13_make_P(Z, "sq")
    irreducibility_prime = None
    for p in h10q.primerange(7, 1000):
        if p == 5 or Pi[-1] % p == 0:
            continue
        if h10q._l13_irred8(Pi, p):
            irreducibility_prime = p
            break
        time.sleep(0.1)
    assert len(Pi) - 1 == 8 and irreducibility_prime == 19

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
        "q1": int(row["q1"]),
        "N": str(row["cert"]["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": str(tied),
        "ramified_empty": not ramified,
        "source": "l17_horizon103.py",
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
        "A": str(A0),
        "delta": str(delta),
        "alpha": str(alpha),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
        "member_syms": member_syms,
        "aligned_syms": aligned_syms,
        "tied": tied,
        "ramified_empty": not ramified,
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "small_emergent": emergent,
        "cofactor_rung": info["rung"],
        "cofactor_prime_proven": info["rung"] == "prime",
        "P_degree": len(Pi) - 1,
        "P_irreducibility_prime": irreducibility_prime,
        "P_irreducible_over_Q": True,
        "P_coefficients": Pi,
    }
    return closure, replay


def main():
    # Explicit self-tests are part of the same proven kernel used below.
    h10q._selftest()
    assert Z == F(-103)
    z_factorization = h10q.factorint(abs(Z.numerator))
    l11_escape_primes = h10q._l11_escape_primes(Z)
    roots, a_pool = discover_l11()
    assert z_factorization == {103: 1}
    assert roots == [] and a_pool == []
    assert l11_escape_primes == []
    assert CELL not in h10q._L11_CLASSES

    l12_attempts, l12_rows = discover_l12()
    assert l12_rows

    OUT.parent.mkdir(parents=True, exist_ok=True)
    closure = None
    with OUT.open("w", encoding="utf-8") as out:
        emit(out, {"type": "meta", "schema": "l17-off-grid-closure-v1",
                   "task": "Off-grid H-class closure", "engine": "h10q.py",
                   "proven_primality_only": True, "selftest": True,
                   "cell": [W, list(UT)], "z": str(Z), "family": "ESC",
                   "l11_roots": roots, "l11_a_pool": a_pool,
                   "l11_escape_primes_1_mod_4": l11_escape_primes,
                   "max_k": MAX_K, "class_q_max": CLASS_Q_MAX,
                   "note": "fresh L11/L12 discovery; temporary L12 row only"})
        for attempt in l12_attempts:
            emit(out, {"type": "class-attempt", "family": "ESC-fresh",
                       "cell": [W, list(UT)], **attempt})
        for row in l12_rows:
            cert = row["cert"]
            emit(out, {"type": "aligned-class", "family": row["family"],
                       "cell": [W, list(UT)], "a": row["a"],
                       "eps": row["eps"], "f": row["f"], "q1": row["q1"],
                       "tau": row["tau"], "N": str(cert["N"]),
                       "S": cert["S"], "ks": cert["ks"],
                       "syms": cert["syms"], "Q0": str(cert["Q0"]),
                       "excluded": cert.get("excluded", []), "ok": cert["ok"]})

        # Search discovered classes in q1 order.  The target closes at k=0
        # for q1=41, eps=-1; all refusals/nonzero verdicts remain evidence-free.
        for row in l12_rows:
            deadline = time.time() + float(os.environ.get("H103_MEMBER_SECONDS", "900"))
            for k in range(0, MAX_K + 1):
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

    report_lines = [
        "H103 off-grid closure recipe and replay",
        f"STATUS: {'VERIFIED' if closure is not None else 'OPEN'}",
        f"TARGET: cell [103,[-1,1]], z={Z}, family=ESC",
        "ENGINE: math/h10q/h10q.py proven-primality engine; refusals are not evidence",
        f"L11: z factorization={json.dumps(z_factorization, sort_keys=True)}; "
        f"escape primes 1 mod 4={l11_escape_primes}; roots={roots}; a_pool={a_pool}",
        f"DISCOVERY: ESC attempts q<= {CLASS_Q_MAX}; aligned classes={len(l12_rows)}; max_k={MAX_K}",
    ]
    if closure is None:
        report_lines.append("RESULT: OPEN after the bounded search; no proved zero member was found.")
    else:
        report_lines.extend([
            "L12: f=103, temporary off-grid tuple, canonical class alignment replayed",
            "L13: smooth_emergent + cofactor_decide replayed; zero is the only closure verdict",
            f"RESULT: VERIFIED closure row {json.dumps(closure, sort_keys=True)}",
            "The closure row uses the exact 13-field schema of data/l13h_all_closures.json;"
            " kernel-replay immediately precedes it in the JSONL artifact.",
        ])
    report_lines.extend([
        f"ARTIFACT: {OUT}",
        f"REPLAY: python math/h10q/l17_horizon103.py (run from repo or set H10Q_ROOT)",
    ])
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("wrote", OUT)
    print("wrote", REPORT)
    print("closure", closure)


if __name__ == "__main__":
    main()
