#!/usr/bin/env python3
"""Deterministic off-grid class discovery and closure replay for (113, -113).

The script enumerates the fresh L11 construction, probes the canonical L12
escape construction through an ephemeral table entry, and sends aligned class
members through the proven-primality L13 ladder.  It never persists a class
table mutation and never treats a refusal as evidence.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l12_class  # noqa: E402
import l13_filter  # noqa: E402

CELL = (113, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
OUT = ROOT / "data" / "l17_horizon113.jsonl"
REPORT = ROOT / "data" / "l17_horizon113.report"

A_T_MAX = 2
CLASS_Q_MIN = 3
CLASS_Q_MAX = 200
MAX_K = 50
PACE_SECONDS = 0.1
CLOSURE_FIELDS = (
    "family",
    "cell",
    "a",
    "eps",
    "f",
    "q1",
    "N",
    "k_zero",
    "Q",
    "rung",
    "tied",
    "ramified_empty",
    "source",
)


def emit(out, record):
    """Write one deterministic JSONL row without sorting mixed symbol keys."""
    out.write(json.dumps(record, separators=(",", ":")) + "\n")
    out.flush()


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def cert_status(cert) -> str:
    if cert is not None and cert.get("ok"):
        return "aligned"
    return "cert-not-ok"


def cert_summary(cert) -> dict:
    if cert is None:
        return {}
    summary = {}
    for key in ("N", "S", "ks", "syms", "Q0", "excluded", "ok"):
        if key not in cert:
            continue
        value = cert[key]
        summary[key] = str(value) if key in ("N", "Q0") else value
    return summary


def fresh_l12_cert(eps: int, q1: int):
    """Call the canonical L12 cert with one temporary off-grid ESC tuple."""
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


def discover_l11(q1_candidates: list[int]):
    """Enumerate L11 roots, the odd canonical a-pool, and class attempts."""
    roots = [r for r in range(1, W) if (1 + 4 * r * r) % W == 0]
    a_pool = sorted(
        {
            abs(r + sign * W * t)
            for r in roots
            for sign in (1, -1)
            for t in range(A_T_MAX + 1)
            if r + sign * W * t != 0
            and abs(r + sign * W * t) % 2 == 1
        }
    )
    attempts = []
    rows = []
    for a_int in a_pool:
        a = F(a_int)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        for eps in (1, -1):
            for q1 in q1_candidates:
                time.sleep(PACE_SECONDS)
                try:
                    cert = h10q._l10_class_cert(a, Z, tau, eps, q1)
                    status = cert_status(cert)
                except Exception as exc:
                    cert = None
                    status = refusal(exc)
                attempt = {
                    "a": a_int,
                    "tau": str(tau),
                    "eps": eps,
                    "q1": q1,
                    "status": status,
                    **cert_summary(cert),
                }
                attempts.append(attempt)
                if status == "aligned":
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


def discover_esc(q1_candidates: list[int]):
    """Probe the f=w ESC path without persisting the off-grid tuple."""
    attempts = []
    rows = []
    for eps in (1, -1):
        for q1 in q1_candidates:
            time.sleep(PACE_SECONDS)
            try:
                cert = fresh_l12_cert(eps, q1)
                status = cert_status(cert)
            except Exception as exc:
                cert = None
                status = refusal(exc)
            attempt = {
                "f": W,
                "eps": eps,
                "q1": q1,
                "status": status,
                **cert_summary(cert),
            }
            attempts.append(attempt)
            if status == "aligned":
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


def member_decide(row: dict, k: int):
    """Run one Q=q1+kN member through primality and the complete L13 ladder."""
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
        Q_is_prime = h10q._is_prime(Q)
    except Exception as exc:
        result["verdict"] = refusal(exc)
        return result, None
    result["Q_prime_proved"] = bool(Q_is_prime)
    if not Q_is_prime:
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
    if smooth_status != "cofactor-big":
        result["verdict"] = smooth_status
        return result, None

    try:
        verdict, info = l13_filter.cofactor_decide(
            a, Z, tau, b, detail=detail, deadline=None
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


def replay_class(row: dict):
    """Re-run the decisive aligned class certificate from its defining path."""
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    q1 = int(row["q1"])
    if row["family"] == "L11-fresh":
        cert = h10q._l10_class_cert(a, Z, tau, eps, q1)
    else:
        cert = fresh_l12_cert(eps, q1)
    assert cert is not None and cert.get("ok")
    assert cert == row["cert"]
    assert all(value == 1 for value in cert["syms"].values())
    return cert


def close_certificate(row: dict, hit):
    """Replay class, kernel identity, member alignment, and the zero ladder."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])

    replayed_cert = replay_class(row)
    assert Q == q1 + k * int(replayed_cert["N"])
    assert h10q._is_prime(Q)
    assert f == 1 or h10q._is_prime(f)
    assert b == F(eps * f * Q)

    member_syms = h10q._l12_syms(a, Z, tau, b)
    aligned_syms = l13_filter.aligned(a, Z, tau, b)
    assert member_syms is not None
    assert aligned_syms is not None
    assert all(value == 1 for value in member_syms.values())
    assert all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, Z, tau, b)
    if rung == "square":
        assert smooth_status == "zero"
        replayed_verdict = "zero"
        replayed_info = None
    else:
        assert smooth_status == "cofactor-big"
        replayed_verdict, replayed_info = l13_filter.cofactor_decide(
            a, Z, tau, b, detail=smooth_detail, deadline=None
        )
        assert replayed_verdict == "zero"
        assert replayed_info["rung"] == rung
        assert replayed_info == prior_info
        assert all(
            place[1] == 1 and place[2] == "proved"
            for place in replayed_info["places"]
        )

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

    try:
        tied = str(h10q._l7_tied_status(a, b, Z, tau))
    except Exception as exc:
        tied = refusal(exc)
    try:
        ramified_empty = h10q.ramified(x0, d0) == []
    except Exception as exc:
        ramified_empty = refusal(exc)

    closure = {
        "family": "L11" if row["family"] == "L11-fresh" else "ESC",
        "cell": [W, list(UT)],
        "a": str(a),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(replayed_cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": tied,
        "ramified_empty": ramified_empty,
        "source": "l17_horizon113.py",
    }
    assert tuple(closure) == CLOSURE_FIELDS

    replay = {
        "type": "kernel-replay",
        "proof_label": "PROVED",
        "cell": [W, list(UT)],
        "family": row["family"],
        "a": str(a),
        "b": str(b),
        "tau": str(tau),
        "z": str(Z),
        "q1": q1,
        "Q": str(Q),
        "k": k,
        "N": str(replayed_cert["N"]),
        "class_cert_replayed": True,
        "class_syms": replayed_cert["syms"],
        "member_syms": member_syms,
        "aligned_syms": aligned_syms,
        "member_alignment_all_plus_one": True,
        "smooth_status": smooth_status,
        "emergent_verdict": replayed_verdict,
        "cofactor_info": replayed_info,
        "cofactor_replay_matches_initial": replayed_info == prior_info,
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
        "tied": tied,
        "ramified_empty": ramified_empty,
        "l13_replayed": True,
        "refusals_used_as_evidence": False,
    }
    return closure, replay


def member_order_key(row: dict):
    """Put the exact discovered ESC witness first, then use a stable order."""
    witness = (
        row["family"] == "ESC-fresh"
        and int(row["eps"]) == 1
        and int(row["q1"]) == 41
    )
    return (
        0 if witness else 1,
        0 if row["family"] == "ESC-fresh" else 1,
        int(row["q1"]),
        0 if int(row["eps"]) == 1 else 1,
        abs(int(row["a"])),
    )


def attempt_status_counts(attempts: list[dict]) -> dict[str, int]:
    return dict(sorted(Counter(record["status"] for record in attempts).items()))


def main():
    h10q._selftest()

    assert Z == F(-113)
    assert h10q._is_prime(W)
    assert h10q.factorint(abs(Z.numerator)) == {113: 1}
    assert CELL not in h10q._L11_CLASSES

    q1_candidates = list(h10q.primerange(CLASS_Q_MIN, CLASS_Q_MAX + 1))
    roots, a_pool, l11_attempts, l11_rows = discover_l11(q1_candidates)
    esc_attempts, esc_rows = discover_esc(q1_candidates)
    assert roots == [49, 64]
    assert a_pool == [49, 177, 275]
    assert l11_rows
    assert esc_rows
    assert any(
        row["family"] == "ESC-fresh"
        and row["eps"] == 1
        and row["q1"] == 41
        and row["cert"]["ok"]
        for row in esc_rows
    )

    all_class_attempts = [
        {"type": "class-attempt", "family": "L11-fresh", "cell": [W, list(UT)], **record}
        for record in l11_attempts
    ] + [
        {"type": "class-attempt", "family": "ESC-fresh", "cell": [W, list(UT)], **record}
        for record in esc_attempts
    ]
    assert all(
        record["status"] in ("aligned", "cert-not-ok")
        or record["status"].startswith("refused:")
        for record in all_class_attempts
    )

    rows = l11_rows + esc_rows
    rows.sort(key=member_order_key)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    closure = None
    replay = None
    member_attempts = []

    with OUT.open("w", encoding="utf-8") as out:
        meta = {
            "type": "meta",
            "schema": "l17-off-grid-closure-v1",
            "task": "Off-grid closure at (w,z)=(113,-113)",
            "engine": "math/h10q/h10q.py proven-primality only",
            "proven_primality_only": True,
            "selftest_first": True,
            "cell": [W, list(UT)],
            "z": str(Z),
            "source_script": "math/h10q/l17_horizon113.py",
            "l11_roots": roots,
            "l11_a_pool": a_pool,
            "probe_space": {
                "root_r_range": [1, W - 1],
                "root_congruence": "1+4*r^2 == 0 (mod w)",
                "a_t_range": [0, A_T_MAX],
                "a_signs": [1, -1],
                "a_parity": "odd",
                "eps": [1, -1],
                "q1_range_inclusive": [CLASS_Q_MIN, CLASS_Q_MAX],
                "q1_candidates": q1_candidates,
                "families": ["L11-fresh", "ESC-fresh"],
                "esc_f": W,
                "member_k_range_inclusive": [0, MAX_K],
                "member_order": "discovered ESC (eps=1,q1=41) first, then ESC/L11 by q1,eps,a",
                "stop_rule": "first zero verdict that passes the full in-script replay",
            },
            "refusals_are_evidence": False,
        }
        emit(out, meta)
        for record in all_class_attempts:
            emit(out, record)
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

        for row in rows:
            for k in range(MAX_K + 1):
                time.sleep(PACE_SECONDS)
                result, hit = member_decide(row, k)
                member_record = {
                    "type": "member",
                    "family": row["family"],
                    "cell": [W, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    **result,
                }
                emit(out, member_record)
                member_attempts.append(member_record)
                if hit is not None:
                    closure, replay = close_certificate(row, hit)
                    emit(out, replay)
                    emit(out, closure)
                    break
            if closure is not None:
                break

        if closure is None:
            emit(
                out,
                {
                    "type": "summary",
                    "status": "OPEN",
                    "proof_label": "EVIDENCE",
                    "cell": [W, list(UT)],
                    "reason": "no proved zero member in the complete bounded scan",
                    "refusals_used_as_evidence": False,
                },
            )
        else:
            emit(
                out,
                {
                    "type": "summary",
                    "status": "VERIFIED",
                    "proof_label": "PROVED",
                    "cell": [W, list(UT)],
                    "family": closure["family"],
                    "k_zero": closure["k_zero"],
                    "Q": closure["Q"],
                    "replay": "full in-script",
                },
            )

    artifact_records = [
        json.loads(line)
        for line in OUT.read_text(encoding="utf-8").splitlines()
        if line
    ]
    assert artifact_records[0]["type"] == "meta"
    closure_records = [
        record for record in artifact_records if tuple(record) == CLOSURE_FIELDS
    ]
    if closure is None:
        assert closure_records == []
    else:
        assert closure_records == [closure]
        assert closure["cell"] == [113, [-1, 1]]
        assert closure["family"] == "ESC"
        assert closure["a"] == "1"
        assert closure["eps"] == 1 and closure["f"] == 113
        assert closure["q1"] == 41 and closure["k_zero"] == 0
        assert closure["Q"] == "41" and closure["rung"] == "prime"
        assert replay is not None and replay["l13_replayed"]

    l11_counts = attempt_status_counts(l11_attempts)
    esc_counts = attempt_status_counts(esc_attempts)
    report_lines = [
        "H113 off-grid closure discovery and replay",
        "STATUS: " + ("VERIFIED" if closure is not None else "OPEN"),
        "LABEL: "
        + (
            "PROVED (exact machine verification, replayed)"
            if closure is not None
            else "EVIDENCE (bounded negative scan only)"
        ),
        "TARGET: cell [113,[-1,1]], z=-113",
        "ENGINE: math/h10q/h10q.py proven-primality engine only",
        "RULE: refusals are logged and are never evidence",
        "PROBE SPACE: " + json.dumps(meta["probe_space"], separators=(",", ":")),
        "L11 ROOTS: " + json.dumps(roots),
        "L11 A-POOL: " + json.dumps(a_pool),
        "L11 STATUS COUNTS: " + json.dumps(l11_counts, separators=(",", ":")),
        "ESC STATUS COUNTS: " + json.dumps(esc_counts, separators=(",", ":")),
        "ALIGNED CLASSES: L11=" + str(len(l11_rows)) + ", ESC=" + str(len(esc_rows)),
        "CLASS ATTEMPT LOG:",
    ]
    report_lines.extend(
        json.dumps(record, separators=(",", ":")) for record in all_class_attempts
    )
    report_lines.append("MEMBER ATTEMPT LOG:")
    report_lines.extend(
        json.dumps(record, separators=(",", ":")) for record in member_attempts
    )
    if closure is None:
        report_lines.append(
            "RESULT: OPEN after the complete bounded probe; no refusal contributes evidence."
        )
    else:
        report_lines.extend(
            [
                "RESULT: VERIFIED closure row "
                + json.dumps(closure, separators=(",", ":")),
                "REPLAY: aligned class certificate, all controlled member symbols, "
                "kernel identity, smooth_emergent, and cofactor_decide zero were rerun.",
                "SCHEMA: closure row has exactly the 13 fields of data/l13h_all_closures.json.",
            ]
        )
    report_lines.extend(
        [
            "ARTIFACT: " + str(OUT),
            "REPORT: " + str(REPORT),
            "REPLAY COMMAND: python3 math/h10q/l17_horizon113.py",
        ]
    )
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(json.dumps(closure, separators=(",", ":")))
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
