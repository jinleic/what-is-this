#!/usr/bin/env python3
"""Bounded off-grid closure sweep for w in {163, 167, 173, 179}.

Each cell is z=-w (the stored cell key is ``[w,[-1,1]]``).  The script
reconstructs the canonical L11 and ESC classes without modifying the frozen
class tables, searches the fixed k=0 member rung, and fully replays every
closure through the proven-primality engine in h10q.py.  The negative w=179
result is a proved obstruction for this L11/ESC protocol, not a claim of
nonexistence outside the protocol.
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


W_VALUES = (163, 167, 173, 179)
UT = (-1, 1)
CLASS_Q_MIN = 3
CLASS_Q_MAX = 199
A_T_MAX = 2
MAX_K = 0
PACE_SECONDS = 0.1
MEMBER_SECONDS = 300.0
OUT = ROOT / "data" / "l18_horizon_sweep_b.jsonl"
REPORT = Path("/tmp/l18_horizon_sweep_b.report")

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

EXPECTED_CLOSURES = {
    163: {
        "family": "ESC",
        "cell": [163, [-1, 1]],
        "a": "1",
        "eps": -1,
        "f": 163,
        "q1": 31,
        "N": "8118888818554080",
        "k_zero": 0,
        "Q": "31",
        "rung": "prime",
        "tied": "True",
        "ramified_empty": True,
        "source": "l18_horizon_sweep_b.py",
    },
    167: {
        "family": "ESC",
        "cell": [167, [-1, 1]],
        "a": "1",
        "eps": -1,
        "f": 167,
        "q1": 31,
        "N": "3927933644755680",
        "k_zero": 0,
        "Q": "31",
        "rung": "prime",
        "tied": "True",
        "ramified_empty": True,
        "source": "l18_horizon_sweep_b.py",
    },
    173: {
        "family": "ESC",
        "cell": [173, [-1, 1]],
        "a": "1",
        "eps": 1,
        "f": 173,
        "q1": 11,
        "N": "4686108096892320",
        "k_zero": 0,
        "Q": "11",
        "rung": "prime",
        "tied": "True",
        "ramified_empty": True,
        "source": "l18_horizon_sweep_b.py",
    },
}


def emit(handle, record: dict) -> None:
    """Write one deterministic JSONL record; symbol maps have mixed key types."""
    handle.write(json.dumps(record, separators=(",", ":")) + "\n")
    handle.flush()


def refusal(exc: BaseException) -> str:
    """Log an engine refusal without converting it into negative evidence."""
    return "refused:" + type(exc).__name__


def cert_status(cert: dict | None) -> str:
    if cert is None:
        return "ineligible"
    return "aligned" if cert.get("ok") else "cert-not-ok"


def cert_summary(cert: dict | None) -> dict:
    if cert is None:
        return {"cert_present": False}
    return {
        "cert_present": True,
        "ok": bool(cert.get("ok")),
        "N": str(cert["N"]),
        "S": cert["S"],
        "ks": cert["ks"],
        "syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "excluded": cert.get("excluded", []),
    }


def five_character(q1: int) -> int:
    """Return (5|q1), with 0 at q1=5."""
    return 0 if q1 == 5 else h10q.legendre(F(5), q1)


def fresh_l12_cert(w: int, eps: int, q1: int):
    """Call l12_class_cert with one temporary canonical ESC tuple f=w."""
    cell = (w, UT)
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(cell, sentinel)
    l12_class._L12_ESCAPE[cell] = (w, eps, q1)
    try:
        return l12_class.l12_class_cert(w, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(cell, None)
        else:
            l12_class._L12_ESCAPE[cell] = previous


def discover_l11(w: int, q1_candidates: list[int]):
    """Exhaust the prescribed roots and the canonical odd lifted a-pool."""
    z = F(-w)
    roots = [r for r in range(1, w) if (1 + 4 * r * r) % w == 0]
    a_pool = sorted(
        {
            abs(r + sign * w * shift)
            for r in roots
            for sign in (1, -1)
            for shift in range(A_T_MAX + 1)
            if r + sign * w * shift != 0
            and abs(r + sign * w * shift) % 2 == 1
        }
    )
    attempts: list[dict] = []
    rows: list[dict] = []
    for a_int in a_pool:
        a = F(a_int)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        for eps in (1, -1):
            for q1 in q1_candidates:
                time.sleep(PACE_SECONDS)
                try:
                    cert = h10q._l10_class_cert(a, z, tau, eps, q1)
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


def discover_esc(w: int, q1_candidates: list[int]):
    """Scan canonical f=w ESC classes in both (5|q1) residue branches."""
    attempts: list[dict] = []
    rows: list[dict] = []
    for eps in (1, -1):
        for q1 in q1_candidates:
            time.sleep(PACE_SECONDS)
            branch = five_character(q1)
            try:
                cert = fresh_l12_cert(w, eps, q1)
                status = cert_status(cert)
            except Exception as exc:
                cert = None
                status = refusal(exc)
            summary = cert_summary(cert)
            syms = None if cert is None else cert["syms"]
            attempt = {
                "f": w,
                "eps": eps,
                "q1": q1,
                "five_character": branch,
                "five_character_name": "(5|q1)",
                "status": status,
                "sym_5": None if syms is None else syms.get(5),
                "sym_q1": None if syms is None else syms.get("q1"),
                **summary,
            }
            attempts.append(attempt)
            if status == "aligned":
                rows.append(
                    {
                        "family": "ESC-fresh",
                        "a": "1",
                        "eps": eps,
                        "f": w,
                        "q1": q1,
                        "tau": str(F(3, 5)),
                        "cert": cert,
                    }
                )
    return attempts, rows


def status_counts(records: list[dict]) -> dict[str, int]:
    return dict(sorted(Counter(record["status"] for record in records).items()))


def summarize_attempt_group(records: list[dict]) -> dict:
    pair_counts = Counter(
        f"{record['sym_5']},{record['sym_q1']}"
        for record in records
        if record.get("cert_present")
    )
    return {
        "q1_count": len({record["q1"] for record in records}),
        "attempt_count": len(records),
        "certificate_count": sum(bool(record.get("cert_present")) for record in records),
        "aligned_count": sum(record["status"] == "aligned" for record in records),
        "refusal_count": sum(record["status"].startswith("refused:") for record in records),
        "status_counts": status_counts(records),
        "symbol_pair_counts_sym5_symq1": dict(sorted(pair_counts.items())),
    }


def esc_branch_counts(attempts: list[dict]) -> dict:
    """Count the +1 and -1 character branches separately, including eps."""
    result = {}
    for branch in (-1, 0, 1):
        branch_records = [
            record for record in attempts if record["five_character"] == branch
        ]
        result[str(branch)] = {
            **summarize_attempt_group(branch_records),
            "by_eps": {
                str(eps): summarize_attempt_group(
                    [record for record in branch_records if record["eps"] == eps]
                )
                for eps in (1, -1)
            },
        }
    assert result["-1"]["attempt_count"] > 0
    assert result["1"]["attempt_count"] > 0
    return result


def member_decide(w: int, row: dict, k: int, deadline: float):
    """Run one Q=q1+kN member through the complete canonical L13 ladder."""
    cert = row["cert"]
    q1 = int(row["q1"])
    modulus = int(cert["N"])
    eps = int(row["eps"])
    f = int(row["f"])
    a = F(row["a"])
    tau = F(row["tau"])
    z = F(-w)
    Q = q1 + k * modulus
    result = {
        "k": k,
        "Q": str(Q),
        "q1": q1,
        "N": str(modulus),
    }

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
        smooth_status, detail = l13_filter.smooth_emergent(a, z, tau, b)
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
            a, z, tau, b, detail=detail, deadline=deadline
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


def replay_class(w: int, row: dict):
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    q1 = int(row["q1"])
    if row["family"] == "L11-fresh":
        cert = h10q._l10_class_cert(a, F(-w), tau, eps, q1)
    else:
        cert = fresh_l12_cert(w, eps, q1)
    assert cert is not None and cert.get("ok")
    assert cert == row["cert"]
    assert all(value == 1 for value in cert["syms"].values())
    return cert


def close_certificate(w: int, row: dict, hit):
    """Replay class, kernel identity, alignment, and proven zero verdict."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])
    z = F(-w)

    cert = replay_class(w, row)
    assert Q == q1 + k * int(cert["N"])
    assert h10q._is_prime(Q)
    assert f == 1 or h10q._is_prime(f)
    assert b == F(eps * f * Q)

    member_syms = h10q._l12_syms(a, z, tau, b)
    aligned_syms = l13_filter.aligned(a, z, tau, b)
    assert member_syms is not None and all(value == 1 for value in member_syms.values())
    assert aligned_syms is not None and all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, z, tau, b)
    if rung == "square":
        assert smooth_status == "zero"
        replayed_verdict = "zero"
        replayed_info = None
    else:
        assert smooth_status == "cofactor-big"
        replayed_verdict, replayed_info = l13_filter.cofactor_decide(
            a, z, tau, b, detail=smooth_detail, deadline=time.time() + MEMBER_SECONDS
        )
        assert replayed_verdict == "zero"
        assert replayed_info == prior_info
        assert replayed_info["rung"] == rung
        assert replayed_info["places"]
        assert all(
            place[1] == 1 and place[2] == "proved"
            for place in replayed_info["places"]
        )
        if rung == "prime":
            assert len(replayed_info["places"]) == 1
            assert h10q._is_prime(int(replayed_info["places"][0][0]))

    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    c0 = h10q._sun_h(a, b, Z)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    polynomial = h10q._l10_P(a, Z, D, A, delta, s)
    evaluated = sum(
        coefficient * b**degree
        for degree, coefficient in enumerate(polynomial)
    )
    assert M0 == evaluated / (b**4 * D * D * A * A)

    tied = h10q._l7_tied_status(a, b, z, tau)
    ramified = h10q.ramified(x0, d0)
    assert tied is True
    assert ramified == []

    closure = {
        "family": "L11" if row["family"] == "L11-fresh" else "ESC",
        "cell": [w, list(UT)],
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
        "source": "l18_horizon_sweep_b.py",
    }
    assert tuple(closure) == CLOSURE_FIELDS

    replay = {
        "type": "kernel-replay",
        "proof_label": "PROVED",
        "cell": [w, list(UT)],
        "family": row["family"],
        "a": str(a),
        "b": str(b),
        "tau": str(tau),
        "z": str(z),
        "q1": q1,
        "Q": str(Q),
        "k": k,
        "N": str(cert["N"]),
        "class_cert_replayed": True,
        "class_syms": cert["syms"],
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
        "P_degree": len(polynomial) - 1,
        "kernel_identity": True,
        "tied": tied,
        "ramified": ramified,
        "ramified_empty": not ramified,
        "l13_replayed": True,
        "refusals_used_as_evidence": False,
    }
    return closure, replay


def prove_179_protocol_obstruction(
    roots: list[int], esc_attempts: list[dict], branch_counts: dict
) -> dict:
    """Prove the L11/canonical-ESC protocol wall at z=-179."""
    w = 179
    assert roots == [] and w % 4 == 3
    assert h10q.legendre(F(-1), w) == -1
    assert h10q.legendre(F(8 * w), 5) == -1

    refusals = [
        record for record in esc_attempts if record["status"].startswith("refused:")
    ]
    assert refusals == []
    eligible = [record for record in esc_attempts if record["cert_present"]]
    assert eligible
    assert {record["five_character"] for record in eligible} == {-1, 1}
    for record in eligible:
        character = record["five_character"]
        assert record["sym_5"] == -character
        assert record["sym_q1"] == character
        assert record["status"] == "cert-not-ok"
    assert not any(record["status"] == "aligned" for record in esc_attempts)

    ineligible_q1 = sorted(
        {record["q1"] for record in esc_attempts if not record["cert_present"]}
    )
    assert ineligible_q1 == [3, 5, 7, 179]
    assert branch_counts["-1"]["certificate_count"] > 0
    assert branch_counts["1"]["certificate_count"] > 0

    proof_steps = [
        "L11 requires 1+4r^2=0 mod 179, hence -1 must be a square; this is impossible because 179=3 mod 4.",
        "For canonical ESC, a=1, tau=3/5, alpha=4, b=eps*179*q1, and d0=8b.",
        "At p=5 and eligible q1, b is a 5-adic unit. In g=16/5-((b-1)^2/b)^2 the first term has valuation -1 and the second has valuation at least 0, so v5(g)=-1.",
        "Writing Z=(-179)^3, one has Z=1 mod 5 and D=1-Z-Z^2=-1 mod 5; hence c0=Z^2*g/D has v5(c0)=-1.",
        "Since delta=-4/5 and s=0, M0=16+4*c0^2/5 has v5(M0)=-3, so v5(x0)=v5(4*M0)=-3 is odd while v5(d0)=0.",
        "The odd-prime Hilbert formula gives (x0,d0)_5=(d0|5)=(8*eps*179*q1|5). Because eps in {+1,-1} is a square modulo 5 and (8*179|5)=-1, this equals -(q1|5)=-(5|q1).",
        "The wild q1 symbol in the same L12 certificate is (A|q1)=(5|q1). Therefore the 5-adic and wild symbols are opposite for every eligible q1 and cannot both equal +1.",
    ]
    return {
        "type": "protocol-obstruction",
        "status": "OBSTRUCTED",
        "proof_label": "PROVED",
        "scope": "canonical L11 plus canonical ESC (f=w) protocol only",
        "cell": [179, list(UT)],
        "z": "-179",
        "name": "five-adic/wild reciprocity wall",
        "statement": "For every eps in {+1,-1} and every eligible odd prime q1, Hilb_5=-(5|q1) while Hilb_q1=(5|q1); simultaneous acceptance is impossible.",
        "l11_roots": roots,
        "esc_ineligible_q1": ineligible_q1,
        "esc_branch_counts": branch_counts,
        "proof_steps": proof_steps,
        "bounded_scan_role": "discriminating replay of both residue branches; the obstruction follows from the displayed identity, not from the finite scan",
        "refusals_used_as_evidence": False,
    }


def aligned_class_record(w: int, row: dict) -> dict:
    cert = row["cert"]
    return {
        "type": "aligned-class",
        "proof_label": "PROVED",
        "family": row["family"],
        "cell": [w, list(UT)],
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
    }


def member_order_key(row: dict):
    return (
        0 if row["family"] == "ESC-fresh" else 1,
        int(row["q1"]),
        0 if int(row["eps"]) == 1 else 1,
        abs(int(row["a"])),
    )


def run_cell(w: int, q1_candidates: list[int], handle) -> dict:
    cell = [w, list(UT)]
    z = F(-w)
    assert h10q._is_prime(w)
    assert h10q.factorint(abs(z.numerator)) == {w: 1}
    assert (w, UT) not in h10q._L11_CLASSES

    roots, a_pool, l11_attempts, l11_rows = discover_l11(w, q1_candidates)
    esc_attempts, esc_rows = discover_esc(w, q1_candidates)
    branch_counts = esc_branch_counts(esc_attempts)
    escape_primes = h10q._l11_escape_primes(z)

    cell_meta = {
        "type": "cell-meta",
        "proof_label": "EVIDENCE",
        "cell": cell,
        "z": str(z),
        "l11_roots": roots,
        "l11_a_pool": a_pool,
        "l11_escape_primes_1_mod_4": escape_primes,
        "probe_space": {
            "root_r_range_inclusive": [1, w - 1],
            "root_congruence": "1+4*r^2 == 0 (mod w)",
            "a_shift_range_inclusive": [0, A_T_MAX],
            "a_signs": [1, -1],
            "a_parity": "odd",
            "families": ["L11-fresh", "ESC-fresh"],
            "esc_f": w,
            "eps": [1, -1],
            "q1_range_inclusive": [CLASS_Q_MIN, CLASS_Q_MAX],
            "q1_candidates": q1_candidates,
            "esc_residue_split": "(5|q1)=-1,0,+1 counted separately",
            "member_k_range_inclusive": [0, MAX_K],
            "stop_rule": "first zero verdict passing full in-script replay",
        },
        "refusals_are_evidence": False,
    }
    emit(handle, cell_meta)

    for attempt in l11_attempts:
        emit(
            handle,
            {
                "type": "class-attempt",
                "proof_label": "EVIDENCE",
                "family": "L11-fresh",
                "cell": cell,
                "refusal_counts_as_evidence": False,
                **attempt,
            },
        )
    for attempt in esc_attempts:
        emit(
            handle,
            {
                "type": "class-attempt",
                "proof_label": "EVIDENCE",
                "family": "ESC-fresh",
                "cell": cell,
                "refusal_counts_as_evidence": False,
                **attempt,
            },
        )
    for row in sorted(l11_rows + esc_rows, key=member_order_key):
        emit(handle, aligned_class_record(w, row))

    closure = None
    replay = None
    member_attempts: list[dict] = []
    if w != 179:
        rows = sorted(esc_rows + l11_rows, key=member_order_key)
        for row in rows:
            deadline = time.time() + MEMBER_SECONDS
            for k in range(MAX_K + 1):
                time.sleep(PACE_SECONDS)
                result, hit = member_decide(w, row, k, deadline)
                member_record = {
                    "type": "member",
                    "proof_label": "EVIDENCE",
                    "family": row["family"],
                    "cell": cell,
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    "refusals_used_as_evidence": False,
                    **result,
                }
                member_attempts.append(member_record)
                emit(handle, member_record)
                if hit is not None:
                    closure, replay = close_certificate(w, row, hit)
                    emit(handle, replay)
                    emit(handle, closure)
                    break
            if closure is not None:
                break

    obstruction = None
    if w == 179:
        assert l11_rows == [] and esc_rows == []
        obstruction = prove_179_protocol_obstruction(
            roots, esc_attempts, branch_counts
        )
        emit(handle, obstruction)
    else:
        assert closure == EXPECTED_CLOSURES[w]

    if closure is not None:
        result_record = {
            "type": "cell-result",
            "status": "VERIFIED",
            "proof_label": "PROVED",
            "cell": cell,
            "family": closure["family"],
            "q1": closure["q1"],
            "k_zero": closure["k_zero"],
            "Q": closure["Q"],
            "rung": closure["rung"],
            "esc_branch_counts": branch_counts,
            "refusals_used_as_evidence": False,
        }
    else:
        assert obstruction is not None
        result_record = {
            "type": "cell-result",
            "status": "OBSTRUCTED",
            "proof_label": "PROVED",
            "scope": obstruction["scope"],
            "cell": cell,
            "name": obstruction["name"],
            "esc_branch_counts": branch_counts,
            "refusals_used_as_evidence": False,
        }
    emit(handle, result_record)

    return {
        "w": w,
        "cell": cell,
        "z": str(z),
        "status": result_record["status"],
        "proof_label": result_record["proof_label"],
        "roots": roots,
        "a_pool": a_pool,
        "l11_status_counts": status_counts(l11_attempts),
        "esc_status_counts": status_counts(esc_attempts),
        "l11_aligned_count": len(l11_rows),
        "esc_aligned_count": len(esc_rows),
        "esc_branch_counts": branch_counts,
        "member_status_counts": dict(
            sorted(Counter(record["verdict"] for record in member_attempts).items())
        ),
        "member_attempt_count": len(member_attempts),
        "closure": closure,
        "obstruction": obstruction,
    }


def write_report(results: list[dict], q1_candidates: list[int]) -> None:
    lines = [
        "L18 horizon sweep B: w=163,167,173,179 at u=-1",
        "ENGINE: math/h10q/h10q.py proven-primality engine only",
        "RULE: refusals are logged and never used as evidence",
        "BOUNDS: q1 proven primes in [3,199], L11 shifts t in [0,2], member k=0",
        "RESIDUE EXPERIMENT: every ESC q1 is counted separately by (5|q1)=-1,0,+1 and eps=+1,-1",
        "Q1 CANDIDATES: " + json.dumps(q1_candidates, separators=(",", ":")),
    ]
    for result in results:
        lines.extend(
            [
                "",
                "CELL: " + json.dumps(result["cell"], separators=(",", ":")),
                "STATUS: " + result["status"],
                "LABEL: " + result["proof_label"],
                "L11 ROOTS: " + json.dumps(result["roots"]),
                "L11 A-POOL: " + json.dumps(result["a_pool"]),
                "L11 STATUS COUNTS: "
                + json.dumps(result["l11_status_counts"], separators=(",", ":")),
                "ESC STATUS COUNTS: "
                + json.dumps(result["esc_status_counts"], separators=(",", ":")),
                "ALIGNED CLASSES: L11="
                + str(result["l11_aligned_count"])
                + ", ESC="
                + str(result["esc_aligned_count"]),
                "ESC (5|q1) BRANCH COUNTS: "
                + json.dumps(result["esc_branch_counts"], separators=(",", ":")),
                "MEMBER ATTEMPT COUNT: " + str(result["member_attempt_count"]),
                "MEMBER STATUS COUNTS: "
                + json.dumps(result["member_status_counts"], separators=(",", ":")),
            ]
        )
        if result["closure"] is not None:
            lines.extend(
                [
                    "RESULT: VERIFIED closure row "
                    + json.dumps(result["closure"], separators=(",", ":")),
                    "REPLAY: class certificate, controlled symbols, kernel identity, tied status, ramification, smooth/cofactor ladder, and proven prime cofactor were rerun.",
                ]
            )
        else:
            obstruction = result["obstruction"]
            lines.extend(
                [
                    "RESULT: PROVED protocol obstruction (not an unrestricted nonexistence claim).",
                    "THEOREM: " + obstruction["statement"],
                    "PROOF:",
                ]
            )
            lines.extend(
                f"  {index}. {step}"
                for index, step in enumerate(obstruction["proof_steps"], start=1)
            )
            lines.append(
                "SCAN ROLE: both residue branches were replayed as a discriminating experiment; the theorem, not the finite scan, proves the protocol wall."
            )
    lines.extend(
        [
            "",
            "SWEEP SUMMARY: VERIFIED=3, PROVED_PROTOCOL_OBSTRUCTION=1",
            "ARTIFACT: " + str(OUT),
            "REPORT: " + str(REPORT),
            "REPLAY COMMAND: python3 math/h10q/l18_horizon_sweep_b.py",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def verify_artifact(results: list[dict]) -> None:
    records = [
        json.loads(line)
        for line in OUT.read_text(encoding="utf-8").splitlines()
        if line
    ]
    assert records[0]["type"] == "meta"
    cell_results = [record for record in records if record.get("type") == "cell-result"]
    assert [record["cell"][0] for record in cell_results] == list(W_VALUES)
    assert [record["status"] for record in cell_results] == [
        "VERIFIED",
        "VERIFIED",
        "VERIFIED",
        "OBSTRUCTED",
    ]
    closures = [record for record in records if tuple(record) == CLOSURE_FIELDS]
    assert closures == [EXPECTED_CLOSURES[w] for w in (163, 167, 173)]
    obstructions = [
        record for record in records if record.get("type") == "protocol-obstruction"
    ]
    assert len(obstructions) == 1
    assert obstructions[0]["cell"] == [179, [-1, 1]]
    assert obstructions[0]["proof_label"] == "PROVED"
    assert obstructions[0]["esc_branch_counts"]["-1"]["attempt_count"] > 0
    assert obstructions[0]["esc_branch_counts"]["1"]["attempt_count"] > 0
    assert [result["status"] for result in results] == [
        "VERIFIED",
        "VERIFIED",
        "VERIFIED",
        "OBSTRUCTED",
    ]


def main() -> None:
    h10q._selftest()
    q1_candidates = list(
        h10q.primerange(CLASS_Q_MIN, CLASS_Q_MAX + 1)
    )
    assert q1_candidates[0] == 3 and q1_candidates[-1] == 199

    OUT.parent.mkdir(parents=True, exist_ok=True)
    results = []
    with OUT.open("w", encoding="utf-8") as handle:
        emit(
            handle,
            {
                "type": "meta",
                "schema": "l18-off-grid-horizon-sweep-v1",
                "task": "bounded off-grid closure sweep B",
                "proof_label": "EVIDENCE",
                "engine": "math/h10q/h10q.py proven-primality only",
                "proven_primality_only": True,
                "selftest_first": True,
                "w_values": list(W_VALUES),
                "u": list(UT),
                "source_script": "math/h10q/l18_horizon_sweep_b.py",
                "refusals_are_evidence": False,
            },
        )
        for w in W_VALUES:
            results.append(run_cell(w, q1_candidates, handle))
        emit(
            handle,
            {
                "type": "sweep-result",
                "proof_label": "PROVED",
                "statuses": {
                    str(result["w"]): result["status"] for result in results
                },
                "verified_closures": sum(
                    result["status"] == "VERIFIED" for result in results
                ),
                "proved_protocol_obstructions": sum(
                    result["status"] == "OBSTRUCTED" for result in results
                ),
                "refusals_used_as_evidence": False,
            },
        )

    write_report(results, q1_candidates)
    verify_artifact(results)
    assert REPORT.is_file()
    report_text = REPORT.read_text(encoding="utf-8")
    assert "SWEEP SUMMARY: VERIFIED=3, PROVED_PROTOCOL_OBSTRUCTION=1" in report_text
    assert "Hilb_5=-(5|q1)" in report_text
    print(json.dumps({str(result["w"]): result["status"] for result in results}))
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
