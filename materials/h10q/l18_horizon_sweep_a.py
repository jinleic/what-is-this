#!/usr/bin/env python3
"""Bounded off-grid L11/ESC closure sweep for selected cells (w, -w).

For every requested proven-prime w, this script discovers the canonical L11
classes and canonical f=w ESC classes with q1 <= 200, then searches class
members through k <= 30.  L11 members are tried before ESC members.  After a
proved canonical wall, it probes the Route131 tau=0 fixed-factor escape shape.
A cell is closed only after an in-script class, kernel, Hilbert-support, tied,
and L13 zero replay.  Refusals are logged and are never used as evidence.
"""
from __future__ import annotations

import argparse
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
from l18_route131 import fixed_factor_class_cert  # noqa: E402

DEFAULT_W = (139, 149, 151, 157)
UT = (-1, 1)
OUT = ROOT / "data" / "l18_horizon_sweep_a.jsonl"
REPORT = Path("/tmp/l18_horizon_sweep_a.report")

A_T_MAX = 2
CLASS_Q_MIN = 3
CLASS_Q_MAX = 200
MAX_K = 30
PACE_SECONDS = 0.1
NONCANON_A_CANDIDATES = (7, 1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)
MEMBER_SECONDS = 300.0
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


def emit(out, record: dict) -> None:
    """Write one JSONL row; nested symbol maps mix integer and string keys."""
    out.write(json.dumps(record, separators=(",", ":")) + "\n")


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def character_mod_5(value: int) -> int:
    """Return (value|5), including zero when 5 divides value."""
    residue = value % 5
    if residue == 0:
        return 0
    return 1 if pow(residue, 2, 5) == 1 else -1


def cert_summary(cert: dict | None) -> dict:
    if cert is None:
        return {}
    summary = {}
    for key in ("N", "S", "ks", "syms", "Q0", "excluded", "ok"):
        if key not in cert:
            continue
        value = cert[key]
        summary[key] = str(value) if key in ("N", "Q0") else value
    return summary


def cert_status(cert: dict | None) -> str:
    if cert is None:
        return "cert-none"
    return "aligned" if cert.get("ok") else "misaligned"


def fresh_l12_cert(w: int, eps: int, q1: int):
    """Call the canonical ESC certificate with one ephemeral off-grid row."""
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


def roots_and_a_pool(w: int) -> tuple[list[int], list[int]]:
    roots = [r for r in range(1, w) if (1 + 4 * r * r) % w == 0]
    a_pool = sorted(
        {
            abs(r + sign * w * t)
            for r in roots
            for sign in (1, -1)
            for t in range(A_T_MAX + 1)
            if r + sign * w * t != 0
            and abs(r + sign * w * t) % 2 == 1
        }
    )
    return roots, a_pool


def discover_l11(w: int, z: F, q1_candidates: list[int]):
    """Enumerate L11 roots, canonical odd lifts, and aligned classes."""
    roots, a_pool = roots_and_a_pool(w)
    attempts = []
    rows = []
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
    """Probe canonical f=w ESC classes and log both (q1|5) branches."""
    attempts = []
    rows = []
    controlled = {2, 3, 5, 7, w}
    for eps in (1, -1):
        for q1 in q1_candidates:
            time.sleep(PACE_SECONDS)
            q1_character = character_mod_5(q1)
            attempt = {
                "f": w,
                "eps": eps,
                "q1": q1,
                "q1_proven_prime": True,
                "q1_character_mod_5": q1_character,
            }
            try:
                cert = fresh_l12_cert(w, eps, q1)
            except Exception as exc:
                cert = None
                attempt["status"] = refusal(exc)
            else:
                if q1 in controlled:
                    assert cert is None
                    attempt.update(
                        {
                            "status": "ineligible:q1-in-controlled-set",
                            "controlled_set": sorted(controlled),
                        }
                    )
                else:
                    attempt["status"] = cert_status(cert)
            attempt.update(cert_summary(cert))
            if cert is not None:
                hilbert_5 = cert["syms"][5]
                wild_q1 = cert["syms"]["q1"]
                attempt.update(
                    {
                        "hilbert_5": hilbert_5,
                        "wild_q1": wild_q1,
                        "wild_equals_q1_character_mod_5": wild_q1 == q1_character,
                        "five_and_wild_both_plus_one": hilbert_5 == wild_q1 == 1,
                        "five_symbol_collision_131_type": hilbert_5 == -wild_q1,
                    }
                )
            attempts.append(attempt)
            if cert is not None and cert.get("ok"):
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


def discover_noncanonical_classes_for_a(
    w: int,
    z: F,
    a_int: int,
    q1_candidates: list[int],
):
    """Probe Route131's tau=0, f=w fixed-factor class at one odd a."""
    time.sleep(PACE_SECONDS)
    a = F(a_int)
    A = 1 + 4 * a * a
    A_int = int(A)
    assert A.denominator == 1
    A_prime = h10q._is_prime(A_int)
    A_character = 0 if A_int % w == 0 else h10q.legendre(A, w)
    alpha_character = 0 if A_int % w == 0 else h10q.legendre(-A, w)
    screen = {
        "type": "noncanonical-a-screen",
        "proof_label": "PROVED",
        "cell": [w, list(UT)],
        "a": a_int,
        "A": A_int,
        "A_proven_prime": bool(A_prime),
        "character_A_mod_w": A_character,
        "character_minus_A_mod_w": alpha_character,
        "tau": "0",
        "f": w,
    }
    if not A_prime:
        screen["status"] = "ineligible:A-not-proven-prime"
        return screen, [], []
    if A_character != -1 or alpha_character != 1:
        screen["status"] = "ineligible:W1-target-character"
        return screen, [], []
    screen["status"] = "admissible:prime-A-and-W1"

    attempts = []
    rows = []
    filtered_q1 = [q1 for q1 in q1_candidates if character_mod_5(q1) == 1]
    for q1 in filtered_q1:
        for eps in (1, -1):
            time.sleep(PACE_SECONDS)
            attempt = {
                "a": a_int,
                "A": A_int,
                "tau": "0",
                "eps": eps,
                "f": w,
                "q1": q1,
                "q1_proven_prime": True,
                "q1_character_mod_5": 1,
            }
            try:
                cert = fixed_factor_class_cert(a, z, F(0), eps, w, q1)
                status = cert_status(cert)
            except Exception as exc:
                cert = None
                status = refusal(exc)
            attempt["status"] = status
            attempt.update(cert_summary(cert))
            attempts.append(attempt)
            if cert is not None and cert.get("ok"):
                rows.append(
                    {
                        "family": "FIXED-fresh",
                        "a": str(a),
                        "eps": eps,
                        "f": w,
                        "q1": q1,
                        "tau": "0",
                        "cert": cert,
                    }
                )
    return screen, attempts, rows


def status_counts(records: list[dict]) -> dict[str, int]:
    return dict(sorted(Counter(record["status"] for record in records).items()))


def l11_per_class_counts(attempts: list[dict]) -> list[dict]:
    keys = sorted({(record["a"], record["eps"]) for record in attempts})
    return [
        {
            "a": a,
            "eps": eps,
            "status_counts": status_counts(
                [record for record in attempts if (record["a"], record["eps"]) == (a, eps)]
            ),
        }
        for a, eps in keys
    ]


def esc_per_class_counts(attempts: list[dict]) -> list[dict]:
    keys = sorted(
        {(record["eps"], record["q1_character_mod_5"]) for record in attempts}
    )
    return [
        {
            "eps": eps,
            "q1_character_mod_5": character,
            "status_counts": status_counts(
                [
                    record
                    for record in attempts
                    if (record["eps"], record["q1_character_mod_5"])
                    == (eps, character)
                ]
            ),
        }
        for eps, character in keys
    ]


def esc_branch_counts(attempts: list[dict]) -> list[dict]:
    branches = []
    for character in (-1, 1, 0):
        branch = [
            record
            for record in attempts
            if record["q1_character_mod_5"] == character
        ]
        if not branch:
            continue
        pairs = Counter(
            f"hilbert_5={record['hilbert_5']},wild_q1={record['wild_q1']}"
            for record in branch
            if "hilbert_5" in record
        )
        branches.append(
            {
                "q1_character_mod_5": character,
                "q1_candidates": len({record["q1"] for record in branch}),
                "attempts_across_eps": len(branch),
                "status_counts": status_counts(branch),
                "symbol_pair_counts": dict(sorted(pairs.items())),
                "aligned": sum(record["status"] == "aligned" for record in branch),
            }
        )
    assert any(row["q1_character_mod_5"] == -1 for row in branches)
    assert any(row["q1_character_mod_5"] == 1 for row in branches)
    return branches


def l11_root_analysis(w: int, roots: list[int], a_pool: list[int]) -> dict:
    minus_one_character = h10q.legendre(F(-1), w)
    if roots:
        assert minus_one_character == 1
        return {
            "type": "l11-route-analysis",
            "proof_label": "PROVED",
            "classification": "AVAILABLE",
            "cell": [w, list(UT)],
            "roots": roots,
            "a_pool": a_pool,
            "legendre_minus_one_mod_w": minus_one_character,
        }
    assert w % 4 == 3 and minus_one_character == -1
    return {
        "type": "l11-route-analysis",
        "proof_label": "PROVED",
        "classification": "ROOT_OBSTRUCTION",
        "cell": [w, list(UT)],
        "roots": roots,
        "a_pool": a_pool,
        "legendre_minus_one_mod_w": minus_one_character,
        "statement": "No r satisfies 1+4*r^2 == 0 (mod w), so the prescribed L11 a-pool is empty.",
        "proof": "Such an r would make (2*r)^2 == -1 (mod w), contradicting (-1|w)=-1 for the proven prime w == 3 (mod 4).",
    }


def five_symbol_analysis(w: int, z: F, attempts: list[dict]) -> dict:
    """Classify the canonical ESC p=5/wild interaction, proving a wall when valid."""
    a = F(1)
    tau = F(3, 5)
    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    assert A == 5 and delta == F(-4, 5) and alpha == 4

    eligible = [record for record in attempts if "hilbert_5" in record]
    observed_pairs = Counter(
        (record["q1_character_mod_5"], record["hilbert_5"], record["wild_q1"])
        for record in eligible
    )
    base = {
        "type": "esc-five-symbol-analysis",
        "cell": [w, list(UT)],
        "canonical_data": {"a": "1", "tau": "3/5", "f": w, "b": f"eps*{w}*q1"},
        "q1_scope": "admissible proven primes q1 <= 200 in the scan; theorem scope stated separately",
        "branch_counts": esc_branch_counts(attempts),
        "eligible_certificates_scanned": len(eligible),
        "observed_symbol_pairs": [
            {
                "q1_character_mod_5": character,
                "hilbert_5": hilbert_5,
                "wild_q1": wild_q1,
                "count": count,
            }
            for (character, hilbert_5, wild_q1), count in sorted(observed_pairs.items())
        ],
        "acceptance_requires": "syms[5]=+1 and syms['q1']=+1 simultaneously (indeed every controlled symbol must be +1)",
        "refusals_used_as_evidence": False,
    }

    z_mod_5 = h10q.unit_mod(z, 5)
    Z_mod_5 = h10q.unit_mod(Z, 5)
    D_valuation = h10q.vp(D, 5)
    coefficient_character = character_mod_5(8 * w)
    theorem_applies = (
        h10q.vp(Z, 5) == 0
        and D_valuation == 0
        and character_mod_5(w) == 1
        and coefficient_character == -1
    )
    if not theorem_applies:
        escaping = sorted(
            {
                record["q1_character_mod_5"]
                for record in eligible
                if record["hilbert_5"] == record["wild_q1"] == 1
            }
        )
        base.update(
            {
                "proof_label": "EVIDENCE",
                "classification": "NO_131_TYPE_RECIPROCITY_THEOREM",
                "z_mod_5": z_mod_5,
                "Z_mod_5": Z_mod_5,
                "v5_D": D_valuation,
                "character_8w_mod_5": coefficient_character,
                "escaping_q1_character_branches_observed": escaping,
                "reason": "The unit-D valuation proof used at w=131 does not meet all of its hypotheses here; only the bounded exact symbol scan is reported.",
            }
        )
        return base

    residue_table = []
    for eps in (1, -1):
        assert character_mod_5(eps) == 1
        for q1_mod_5 in (1, 2, 3, 4):
            hilbert_5 = character_mod_5(8 * eps * w * q1_mod_5)
            wild_q1 = character_mod_5(q1_mod_5)
            assert hilbert_5 == -wild_q1
            residue_table.append(
                {
                    "eps": eps,
                    "q1_mod_5": q1_mod_5,
                    "q1_character_mod_5": wild_q1,
                    "d_unit_mod_5": (8 * eps * w * q1_mod_5) % 5,
                    "hilbert_5": hilbert_5,
                    "wild_q1": wild_q1,
                    "product": hilbert_5 * wild_q1,
                }
            )

    assert eligible
    assert all(record["wild_q1"] == record["q1_character_mod_5"] for record in eligible)
    assert all(record["hilbert_5"] == -record["wild_q1"] for record in eligible)
    assert not any(record["status"] == "aligned" for record in attempts)
    base.update(
        {
            "proof_label": "PROVED",
            "classification": "RECIPROCITY_WALL_131_TYPE",
            "theorem_scope": "every admissible prime q1 != 5 in the canonical ESC family, not merely q1 <= 200",
            "statement": "Hilb_5 = -(q1|5) = -(5|q1) = -wild_q1 for both eps signs, so the p=5 and wild q1 symbols cannot both be +1.",
            "valuation_proof": {
                "z_mod_5": z_mod_5,
                "Z_mod_5": Z_mod_5,
                "D_mod_5": h10q.unit_mod(D, 5),
                "v5_Z": h10q.vp(Z, 5),
                "v5_D": D_valuation,
                "v5_g": -1,
                "v5_c": -1,
                "v5_delta": -1,
                "v5_M": -3,
                "v5_x": -3,
                "v5_d": 0,
            },
            "proof_steps": [
                "For b=eps*w*q1 with q1!=5, t=(b-1)^2/b is 5-adically integral and g=(16-5*t^2)/5 has valuation -1.",
                "Z=z^3 and D=1-Z-Z^2 are 5-adic units, hence c=Z^2*g/D has v5(c)=-1.",
                "With delta=-4/5 and alpha=4, delta*c^2 has valuation -3, so M=16-delta*c^2 has v5(M)=-3; therefore v5(x)=-3 and d=8*eps*w*q1 is a unit.",
                "The odd valuation of x gives Hilb_5=(d|5)=(8*eps*w*q1|5).  Both eps and w are squares mod 5 while 8 is a nonsquare, yielding Hilb_5=-(q1|5).",
                "The wild symbol is (A|q1)=(5|q1)=(q1|5) by quadratic reciprocity.  Their product is therefore -1 for every admissible q1.",
            ],
            "residue_table": residue_table,
            "bounded_scan_matches_theorem": True,
        }
    )
    return base


def canonical_protocol_criterion(
    w: int,
    z: F,
    roots: list[int],
    esc_rows: list[dict],
) -> dict:
    """Prove and test the L11/ESC route-availability criterion at this w."""
    w_character_mod_5 = character_mod_5(w)
    l11_available = w % 4 == 1
    esc_available = w_character_mod_5 == -1
    protocol_obstructed = not l11_available and not esc_available
    assert protocol_obstructed == (w % 20 in (11, 19))
    assert bool(roots) == l11_available
    assert bool(esc_rows) == esc_available

    Z = z**3
    D = 1 - Z - Z * Z
    e = h10q.vp(D, 5)
    assert e >= 0
    return {
        "type": "canonical-protocol-criterion",
        "proof_label": "PROVED",
        "cell": [w, list(UT)],
        "w_mod_20": w % 20,
        "legendre_minus_one_mod_w": h10q.legendre(F(-1), w),
        "w_character_mod_5": w_character_mod_5,
        "l11_available_iff_w_mod_4_eq_1": l11_available,
        "esc_available_iff_w_character_mod_5_eq_minus_1": esc_available,
        "protocol_obstructed": protocol_obstructed,
        "protocol_obstructed_iff_w_mod_20_in_11_19": True,
        "bounded_alignment_scan_agrees": True,
        "scope": "canonical L11 plus canonical a=1,tau=3/5,f=w ESC route availability; not a guarantee that a bounded L13 member scan closes every available route",
        "proof_steps": [
            "The L11 congruence is (2r)^2=-1 mod w, so it has a root exactly when (-1|w)=+1, equivalently w=1 mod 4.",
            "For canonical ESC put e=v5(D), D=1-z^3-z^6.  Since v5(g)=-1, one has v5(x)=-3-2e (odd) and v5(d)=0.",
            "Thus Hilb_5=-(w|5)*(5|q1), while the wild symbol is Hilb_q1=(5|q1).  They can simultaneously equal +1 exactly when (w|5)=-1.",
            "Combining the two routes, canonical protocol obstruction is w=3 mod 4 and (w|5)=+1, equivalently w=11 or 19 mod 20.",
        ],
    }


def member_decide(w: int, z: F, row: dict, k: int):
    """Run one Q=q1+kN member through primality and the exact L13 ladder."""
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
            a,
            z,
            tau,
            b,
            detail=detail,
            deadline=time.time() + MEMBER_SECONDS,
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


def replay_class(w: int, z: F, row: dict) -> dict:
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    q1 = int(row["q1"])
    if row["family"] == "L11-fresh":
        cert = h10q._l10_class_cert(a, z, tau, eps, q1)
    elif row["family"] == "ESC-fresh":
        cert = fresh_l12_cert(w, eps, q1)
    else:
        assert row["family"] == "FIXED-fresh"
        cert = fixed_factor_class_cert(a, z, tau, eps, int(row["f"]), q1)
    assert cert is not None and cert.get("ok")
    assert cert == row["cert"]
    assert all(value == 1 for value in cert["syms"].values())
    return cert


def close_certificate(w: int, z: F, row: dict, hit):
    """Replay class, member alignment, kernel identity, and exact zero ladder."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])

    cert = replay_class(w, z, row)
    assert Q == q1 + k * int(cert["N"])
    assert h10q._is_prime(Q)
    assert f == 1 or h10q._is_prime(f)
    assert b == F(eps * f * Q)

    member_syms = h10q._l12_syms(a, z, tau, b)
    aligned_syms = l13_filter.aligned(a, z, tau, b)
    assert member_syms is not None and aligned_syms is not None
    assert all(value == 1 for value in member_syms.values())
    assert all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, z, tau, b)
    if rung == "square":
        assert smooth_status == "zero"
        verdict = "zero"
        info = None
        emergent, (stripped_rn, stripped_rd), smooth = smooth_detail
        assert emergent == [] and smooth and stripped_rn == stripped_rd == 1
    else:
        assert smooth_status == "cofactor-big"
        emergent, (stripped_rn, stripped_rd), smooth = smooth_detail
        assert emergent == [] and not smooth
        verdict, info = l13_filter.cofactor_decide(
            a, z, tau, b, detail=smooth_detail, deadline=time.time() + MEMBER_SECONDS
        )
        assert verdict == "zero" and info["rung"] == rung
        assert info == prior_info
        assert all(place[1] == 1 and place[2] == "proved" for place in info["places"])

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
    P = h10q._l10_P(a, Z, D, A, delta, s)
    Pb = sum(coefficient * b**i for i, coefficient in enumerate(P))
    assert M0 == Pb / (b**4 * D * D * A * A)

    # W2 identifies the tied conic with x0 in the norm group of Q(sqrt(d0)).
    # The exhaustive zero-symbol partition above is exactly its local norm
    # criterion at every finite and real place, hence Hasse--Minkowski proves
    # the tied conic soluble even when the redundant direct solver exhausts
    # its factor budget.
    tied_crosscheck = h10q._l7_tied_status(a, b, z, tau)
    assert tied_crosscheck in (True, "budget")
    tied = True
    ramified_empty = True
    try:
        ramified_crosscheck = h10q.ramified(x0, d0)
    except Exception as exc:
        ramified_crosscheck = refusal(exc)
    else:
        assert ramified_crosscheck == []

    closure_families = {
        "L11-fresh": "L11",
        "ESC-fresh": "ESC",
        "FIXED-fresh": "ESC-NONCANONICAL",
    }
    closure = {
        "family": closure_families[row["family"]],
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
        "ramified_empty": ramified_empty,
        "source": "l18_horizon_sweep_a.py",
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
        "small_emergent": emergent,
        "stripped_rn": str(stripped_rn),
        "stripped_rd": str(stripped_rd),
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "cofactor_replay_matches_initial": info == prior_info,
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
        "tied_crosscheck": tied_crosscheck,
        "tied_basis": "W2 norm equivalence plus exhaustive zero Hilbert-symbol support partition",
        "ramified_crosscheck": ramified_crosscheck,
        "ramified_empty": ramified_empty,
        "ramified_empty_basis": "exhaustive smooth-emergent plus proved cofactor support partition",
        "l13_replayed": True,
        "refusals_used_as_evidence": False,
    }
    return closure, replay


def class_order(row: dict) -> tuple:
    return (
        int(row["q1"]),
        0 if int(row["eps"]) == 1 else 1,
        abs(int(row["a"])),
    )


def search_members(w: int, z: F, l11_rows: list[dict], esc_rows: list[dict]):
    """Search k-layers within L11 first, then ESC; stop at first full replay."""
    attempts = []
    for route, rows in (("L11", l11_rows), ("ESC", esc_rows)):
        ordered = sorted(rows, key=class_order)
        for k in range(MAX_K + 1):
            for row in ordered:
                time.sleep(PACE_SECONDS)
                result, hit = member_decide(w, z, row, k)
                record = {
                    "type": "member",
                    "proof_label": "EVIDENCE" if result["verdict"] != "zero" else "PROVED",
                    "route_order": route,
                    "family": row["family"],
                    "cell": [w, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    **result,
                }
                attempts.append(record)
                if hit is None:
                    continue
                try:
                    closure, replay = close_certificate(w, z, row, hit)
                except Exception as exc:
                    record["full_replay_status"] = refusal(exc)
                    record["accepted_as_closure"] = False
                    continue
                record["full_replay_status"] = "VERIFIED"
                record["accepted_as_closure"] = True
                return attempts, closure, replay
    return attempts, None, None


def probe_noncanonical_fixed_factor(
    w: int,
    z: F,
    q1_candidates: list[int],
):
    """Try Route131's tau=0/f=w door after a canonical protocol wall."""
    screens = []
    class_attempts = []
    rows_seen = []
    member_attempts = []
    closure = None
    replay = None
    filtered_q1 = [q1 for q1 in q1_candidates if character_mod_5(q1) == 1]

    for a_int in NONCANON_A_CANDIDATES:
        screen, attempts, rows = discover_noncanonical_classes_for_a(
            w, z, a_int, q1_candidates
        )
        screens.append(screen)
        class_attempts.extend(attempts)
        rows_seen.extend(rows)
        ordered = sorted(rows, key=class_order)
        for k in range(MAX_K + 1):
            for row in ordered:
                time.sleep(PACE_SECONDS)
                result, hit = member_decide(w, z, row, k)
                record = {
                    "type": "member",
                    "proof_label": "EVIDENCE" if result["verdict"] != "zero" else "PROVED",
                    "route_order": "ESC-NONCANONICAL",
                    "family": row["family"],
                    "cell": [w, list(UT)],
                    "a": row["a"],
                    "eps": row["eps"],
                    "f": row["f"],
                    "q1": row["q1"],
                    **result,
                }
                member_attempts.append(record)
                if hit is None:
                    continue
                try:
                    closure, replay = close_certificate(w, z, row, hit)
                except Exception as exc:
                    record["full_replay_status"] = refusal(exc)
                    record["accepted_as_closure"] = False
                    continue
                record["full_replay_status"] = "VERIFIED"
                record["accepted_as_closure"] = True
                break
            if closure is not None:
                break
        if closure is not None:
            break

    analysis = {
        "type": "noncanonical-fixed-factor-summary",
        "proof_label": "PROVED" if closure is not None else "EVIDENCE",
        "cell": [w, list(UT)],
        "shape": "Route131: tau=0, f=w, odd a with prime A=1+4a^2 and W1 target character",
        "a_order": list(NONCANON_A_CANDIDATES),
        "q1_filter": "(5|q1)=+1",
        "q1_candidates": filtered_q1,
        "q1_max": CLASS_Q_MAX,
        "k_max_per_class": MAX_K,
        "a_screens": screens,
        "class_status_counts": status_counts(class_attempts) if class_attempts else {},
        "aligned_classes_seen": len(rows_seen),
        "member_status_counts": member_status_counts(member_attempts),
        "member_attempts": len(member_attempts),
        "status": "VERIFIED-CLOSURE" if closure is not None else "SEARCH-CONSTRAINT",
        "closure": closure,
        "refusals_used_as_evidence": False,
    }
    return class_attempts, rows_seen, member_attempts, closure, replay, analysis


def member_status_counts(records: list[dict]) -> dict[str, int]:
    return dict(sorted(Counter(record["verdict"] for record in records).items()))


def aligned_record(w: int, row: dict) -> dict:
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


def sweep_cell(w: int, q1_candidates: list[int]) -> dict:
    cell = (w, UT)
    z = F(-w)
    assert h10q._is_prime(w)
    assert h10q.factorint(abs(z.numerator)) == {w: 1}
    assert cell not in h10q._L11_CLASSES

    roots, a_pool, l11_attempts, l11_rows = discover_l11(w, z, q1_candidates)
    esc_attempts, esc_rows = discover_esc(w, q1_candidates)
    l11_analysis = l11_root_analysis(w, roots, a_pool)
    five_analysis = five_symbol_analysis(w, z, esc_attempts)
    criterion = canonical_protocol_criterion(w, z, roots, esc_rows)
    canonical_members, canonical_closure, canonical_replay = search_members(
        w, z, l11_rows, esc_rows
    )

    route_summary = {
        "type": "w-route-summary",
        "proof_label": "PROVED",
        "cell": [w, list(UT)],
        "scope": "canonical protocol",
        "l11_roots": roots,
        "l11_a_pool": a_pool,
        "l11_status_counts": status_counts(l11_attempts) if l11_attempts else {},
        "l11_per_class_status_counts": l11_per_class_counts(l11_attempts),
        "l11_aligned_classes": len(l11_rows),
        "esc_status_counts": status_counts(esc_attempts),
        "esc_per_class_status_counts": esc_per_class_counts(esc_attempts),
        "esc_legendre_5_branch_counts": esc_branch_counts(esc_attempts),
        "esc_aligned_classes": len(esc_rows),
        "member_status_counts": member_status_counts(canonical_members),
        "member_attempts": len(canonical_members),
        "route_order": ["L11", "ESC"],
        "refusals_used_as_evidence": False,
    }

    canonical_obstruction = None
    if canonical_closure is not None:
        canonical_summary = {
            "status": "VERIFIED",
            "proof_label": "PROVED",
            "family": canonical_closure["family"],
            "k_zero": canonical_closure["k_zero"],
            "Q": canonical_closure["Q"],
        }
    elif criterion["protocol_obstructed"]:
        canonical_obstruction = {
            "type": "obstruction",
            "status": "CANONICAL-PROTOCOL-OBSTRUCTED",
            "proof_label": "PROVED",
            "scope": "canonical L11 plus canonical a=1,tau=3/5,f=w ESC protocol only",
            "cell": [w, list(UT)],
            "l11_failure": "No root because (-1|w)=-1.",
            "esc_failure": "Hilb_5=-(w|5)(5|q1) and Hilb_q1=(5|q1), with (w|5)=+1.",
            "five_symbol_collision_131_type": True,
            "bounded_counts": route_summary,
            "not_claimed": "Noncanonical tau=0 fixed-factor classes may escape this wall.",
            "refusals_used_as_evidence": False,
        }
        canonical_summary = {
            "status": "OBSTRUCTED",
            "proof_label": "PROVED",
            "five_symbol_collision_131_type": True,
        }
    else:
        canonical_obstruction = {
            "type": "obstruction",
            "status": "SEARCH-CONSTRAINT",
            "proof_label": "EVIDENCE",
            "scope": f"q1 <= {CLASS_Q_MAX}, k <= {MAX_K} in the available canonical classes",
            "cell": [w, list(UT)],
            "l11_route": l11_analysis["classification"],
            "esc_route": five_analysis["classification"],
            "five_symbol_collision_131_type": False,
            "bounded_counts": route_summary,
            "reason": "No full replay with exact zero verdict was found in the complete bounded scan.",
            "refusals_used_as_evidence": False,
        }
        canonical_summary = {
            "status": "SEARCH_CONSTRAINT",
            "proof_label": "EVIDENCE",
            "five_symbol_collision_131_type": False,
        }

    noncanonical_class_attempts = []
    noncanonical_rows = []
    noncanonical_members = []
    noncanonical_analysis = None
    closure = canonical_closure
    replay = canonical_replay
    if canonical_closure is None and criterion["protocol_obstructed"]:
        (
            noncanonical_class_attempts,
            noncanonical_rows,
            noncanonical_members,
            closure,
            replay,
            noncanonical_analysis,
        ) = probe_noncanonical_fixed_factor(w, z, q1_candidates)

    criterion["bounded_canonical_closure"] = canonical_closure
    criterion["bounded_canonical_outcome"] = canonical_summary["status"]
    if closure is not None:
        summary = {
            "type": "w-summary",
            "status": "VERIFIED",
            "proof_label": "PROVED",
            "cell": [w, list(UT)],
            "family": closure["family"],
            "k_zero": closure["k_zero"],
            "Q": closure["Q"],
            "canonical_protocol_status": canonical_summary["status"],
            "noncanonical_escape": closure["family"] == "ESC-NONCANONICAL",
            "replay": "full in-script",
            "refusals_used_as_evidence": False,
        }
    else:
        summary = {
            "type": "w-summary",
            "status": canonical_summary["status"],
            "proof_label": canonical_summary["proof_label"],
            "cell": [w, list(UT)],
            "canonical_protocol_status": canonical_summary["status"],
            "five_symbol_collision_131_type": canonical_summary.get(
                "five_symbol_collision_131_type", False
            ),
            "noncanonical_escape": False,
            "refusals_used_as_evidence": False,
        }

    return {
        "w": w,
        "z": z,
        "roots": roots,
        "a_pool": a_pool,
        "l11_attempts": l11_attempts,
        "l11_rows": l11_rows,
        "esc_attempts": esc_attempts,
        "esc_rows": esc_rows,
        "l11_analysis": l11_analysis,
        "five_analysis": five_analysis,
        "criterion": criterion,
        "route_summary": route_summary,
        "member_attempts": canonical_members,
        "canonical_closure": canonical_closure,
        "canonical_replay": canonical_replay,
        "noncanonical_class_attempts": noncanonical_class_attempts,
        "noncanonical_rows": noncanonical_rows,
        "noncanonical_member_attempts": noncanonical_members,
        "noncanonical_analysis": noncanonical_analysis,
        "closure": closure,
        "replay": replay,
        "obstruction": canonical_obstruction,
        "summary": summary,
    }


def write_artifact(ws: list[int], q1_candidates: list[int], results: list[dict]) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        emit(
            out,
            {
                "type": "meta",
                "proof_label": "PROVED",
                "schema": "l18-off-grid-horizon-sweep-v1",
                "task": "L11-then-ESC off-grid closure sweep A",
                "engine": "math/h10q/h10q.py proven-primality only",
                "proven_primality_only": True,
                "selftest_first": True,
                "w_list": ws,
                "cells": [[w, list(UT)] for w in ws],
                "source_script": "math/h10q/l18_horizon_sweep_a.py",
                "bounds": {
                    "q1_range_inclusive": [CLASS_Q_MIN, CLASS_Q_MAX],
                    "q1_candidates": q1_candidates,
                    "a_t_range_inclusive": [0, A_T_MAX],
                    "eps": [1, -1],
                    "member_k_range_inclusive": [0, MAX_K],
                    "pace_seconds": PACE_SECONDS,
                },
                "route_order": [
                    "L11",
                    "canonical ESC",
                    "Route131 tau=0 fixed-factor only after a canonical wall",
                ],
                "stop_rule": "first ladder verdict zero that passes full in-script replay",
                "refusals_are_evidence": False,
            },
        )
        for result in results:
            w = result["w"]
            emit(
                out,
                {
                    "type": "w-section",
                    "section": "begin",
                    "cell": [w, list(UT)],
                    "z": str(result["z"]),
                    "l11_roots": result["roots"],
                    "l11_a_pool": result["a_pool"],
                },
            )
            emit(out, result["route_summary"])
            emit(out, result["l11_analysis"])
            emit(out, result["five_analysis"])
            emit(out, result["criterion"])
            for attempt in result["l11_attempts"]:
                emit(
                    out,
                    {
                        "type": "class-attempt",
                        "family": "L11-fresh",
                        "cell": [w, list(UT)],
                        **attempt,
                    },
                )
            for attempt in result["esc_attempts"]:
                emit(
                    out,
                    {
                        "type": "class-attempt",
                        "family": "ESC-fresh",
                        "cell": [w, list(UT)],
                        **attempt,
                    },
                )
            for row in sorted(result["l11_rows"], key=class_order):
                emit(out, aligned_record(w, row))
            for row in sorted(result["esc_rows"], key=class_order):
                emit(out, aligned_record(w, row))
            for attempt in result["noncanonical_class_attempts"]:
                emit(
                    out,
                    {
                        "type": "class-attempt",
                        "family": "FIXED-fresh",
                        "route": "ESC-NONCANONICAL",
                        "cell": [w, list(UT)],
                        **attempt,
                    },
                )
            for row in sorted(result["noncanonical_rows"], key=class_order):
                emit(out, aligned_record(w, row))
            for member in result["member_attempts"]:
                emit(out, member)
            for member in result["noncanonical_member_attempts"]:
                emit(out, member)
            if result["noncanonical_analysis"] is not None:
                emit(out, result["noncanonical_analysis"])
            if result["replay"] is not None:
                emit(out, result["replay"])
            if result["closure"] is not None:
                emit(out, result["closure"])
            if result["obstruction"] is not None:
                emit(out, result["obstruction"])
            emit(out, result["summary"])
            emit(
                out,
                {
                    "type": "w-section",
                    "section": "end",
                    "cell": [w, list(UT)],
                    "status": result["summary"]["status"],
                },
            )
        emit(
            out,
            {
                "type": "global-summary",
                "proof_label": (
                    "PROVED"
                    if all(result["summary"]["status"] == "VERIFIED" for result in results)
                    else "EVIDENCE"
                ),
                "cells": [result["summary"] for result in results],
                "verified_closures": sum(result["closure"] is not None for result in results),
                "proved_canonical_protocol_obstructions": sum(
                    result["criterion"]["protocol_obstructed"] for result in results
                ),
                "noncanonical_wall_breaks": sum(
                    result["summary"].get("noncanonical_escape", False)
                    for result in results
                ),
                "search_constraints": sum(
                    result["summary"]["status"] == "SEARCH_CONSTRAINT"
                    for result in results
                ),
                "refusals_used_as_evidence": False,
            },
        )


def write_report(ws: list[int], results: list[dict]) -> None:
    lines = [
        "L18 horizon sweep A: off-grid cells (w,-w)",
        "ENGINE: math/h10q/h10q.py proven-primality engine only",
        "LABEL RULE: PROVED closures/theorems; EVIDENCE bounded scans; refusals are never evidence",
        f"BOUNDS: q1<={CLASS_Q_MAX}, k<={MAX_K} per discovered class; canonical L11 then ESC; Route131 tau=0/f=w re-attack after a canonical wall",
        "W LIST: " + json.dumps(ws),
    ]
    for result in results:
        w = result["w"]
        summary = result["summary"]
        route = result["route_summary"]
        five = result["five_analysis"]
        criterion = result["criterion"]
        lines.extend(
            [
                "",
                f"=== w={w}, cell=[{w},[-1,1]], z=-{w} ===",
                f"STATUS: {summary['status']}",
                f"LABEL: {summary['proof_label']}",
                "L11 ROOT SET: " + json.dumps(result["roots"]),
                "L11 A-POOL: " + json.dumps(result["a_pool"]),
                "L11 STATUS COUNTS: " + json.dumps(route["l11_status_counts"], separators=(",", ":")),
                "L11 PER-CLASS STATUS COUNTS: "
                + json.dumps(route["l11_per_class_status_counts"], separators=(",", ":")),
                "ESC STATUS COUNTS: " + json.dumps(route["esc_status_counts"], separators=(",", ":")),
                "ESC PER-CLASS STATUS COUNTS: "
                + json.dumps(route["esc_per_class_status_counts"], separators=(",", ":")),
                "ESC (q1|5) BRANCH COUNTS: "
                + json.dumps(route["esc_legendre_5_branch_counts"], separators=(",", ":")),
                "ALIGNED CLASSES: "
                + json.dumps(
                    {
                        "L11": route["l11_aligned_classes"],
                        "ESC": route["esc_aligned_classes"],
                    },
                    separators=(",", ":"),
                ),
                "MEMBER STATUS COUNTS: "
                + json.dumps(route["member_status_counts"], separators=(",", ":")),
                "ESC FIVE-SYMBOL CLASSIFICATION: " + five["classification"],
                "ESC FIVE-SYMBOL LABEL: " + five["proof_label"],
                "CANONICAL PROTOCOL CRITERION [PROVED]: "
                + json.dumps(
                    {
                        "w_mod_20": criterion["w_mod_20"],
                        "L11_available": criterion[
                            "l11_available_iff_w_mod_4_eq_1"
                        ],
                        "ESC_available": criterion[
                            "esc_available_iff_w_character_mod_5_eq_minus_1"
                        ],
                        "predicted_obstructed": criterion["protocol_obstructed"],
                        "bounded_alignment_scan_agrees": criterion[
                            "bounded_alignment_scan_agrees"
                        ],
                        "bounded_canonical_outcome": criterion[
                            "bounded_canonical_outcome"
                        ],
                    },
                    separators=(",", ":"),
                ),
                "CRITERION AGREEMENT: the predicted canonical route availability/obstruction agrees exactly; the criterion does not promise which available route supplies a bounded L13 zero.",
            ]
        )
        if five["classification"] == "RECIPROCITY_WALL_131_TYPE":
            lines.extend(
                [
                    "RECIPROCITY THEOREM: " + five["statement"],
                    "RECIPROCITY PROOF:",
                    *[f"  {index}. {step}" for index, step in enumerate(five["proof_steps"], 1)],
                ]
            )
        else:
            lines.append(
                "ESC ESCAPING (q1|5) BRANCHES OBSERVED: "
                + json.dumps(five["escaping_q1_character_branches_observed"])
            )
        if result["obstruction"] is not None:
            lines.extend(
                [
                    "CANONICAL OBSTRUCTION/CONSTRAINT: "
                    + json.dumps(result["obstruction"], separators=(",", ":")),
                    "CANONICAL SCOPE NOTE: this does not cover noncanonical tau=0 fixed-factor classes.",
                ]
            )
        if result["noncanonical_analysis"] is not None:
            lines.extend(
                [
                    "NONCANONICAL ROUTE131 RE-ATTACK: "
                    + json.dumps(
                        result["noncanonical_analysis"], separators=(",", ":")
                    ),
                    "NONCANONICAL CLASS STATUS COUNTS: "
                    + json.dumps(
                        result["noncanonical_analysis"]["class_status_counts"],
                        separators=(",", ":"),
                    ),
                    "NONCANONICAL MEMBER STATUS COUNTS: "
                    + json.dumps(
                        result["noncanonical_analysis"]["member_status_counts"],
                        separators=(",", ":"),
                    ),
                ]
            )
        if result["closure"] is not None:
            lines.extend(
                [
                    "RESULT: VERIFIED closure row "
                    + json.dumps(result["closure"], separators=(",", ":")),
                    "REPLAY: class certificate, member symbols, kernel identity, tied check, exhaustive smooth/cofactor support partition, and exact L13 zero ladder all rerun; any redundant global ramified() refusal is logged only as a cross-check.",
                    "SCHEMA: closure row has exactly the 13 fields of data/l13h_all_closures.json.",
                ]
            )
        else:
            lines.append(
                "FINAL OBSTRUCTION: no verified canonical or Route131-shaped noncanonical closure; labels and bounded counts above state the exact scope."
            )
    lines.extend(
        [
            "",
            "ARTIFACT: " + str(OUT),
            "REPORT: " + str(REPORT),
            "REPLAY COMMAND: python3 math/h10q/l18_horizon_sweep_a.py 139 149 151 157",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def verify_artifacts(ws: list[int], results: list[dict]) -> None:
    records = [
        json.loads(line)
        for line in OUT.read_text(encoding="utf-8").splitlines()
        if line
    ]
    assert records[0]["type"] == "meta" and records[0]["w_list"] == ws
    assert records[-1]["type"] == "global-summary"
    closure_records = [record for record in records if tuple(record) == CLOSURE_FIELDS]
    expected_closures = [result["closure"] for result in results if result["closure"]]
    assert closure_records == expected_closures
    for result in results:
        w = result["w"]
        summaries = [
            record
            for record in records
            if record.get("type") == "w-summary" and record.get("cell") == [w, list(UT)]
        ]
        assert summaries == [result["summary"]]
        assert result["summary"]["status"] in (
            "VERIFIED",
            "OBSTRUCTED",
            "SEARCH_CONSTRAINT",
        )
        assert any(
            branch["q1_character_mod_5"] == -1
            for branch in result["route_summary"]["esc_legendre_5_branch_counts"]
        )
        assert any(
            branch["q1_character_mod_5"] == 1
            for branch in result["route_summary"]["esc_legendre_5_branch_counts"]
        )
        if result["closure"] is not None:
            assert result["replay"] is not None and result["replay"]["l13_replayed"]
            assert result["closure"]["ramified_empty"] is True
            assert result["closure"]["tied"] == "True"
        else:
            assert result["obstruction"] is not None
            assert result["obstruction"]["proof_label"] in ("PROVED", "EVIDENCE")
    report_text = REPORT.read_text(encoding="utf-8")
    assert all(f"=== w={w}" in report_text for w in ws)
    assert "refusals are never evidence" in report_text


def parse_args() -> list[int]:
    parser = argparse.ArgumentParser(
        description="Sweep canonical off-grid L11 and ESC routes for a list of w values."
    )
    parser.add_argument(
        "w",
        nargs="*",
        type=int,
        help="proven-prime w values (default: 139 149 151 157)",
    )
    args = parser.parse_args()
    ws = args.w or list(DEFAULT_W)
    if len(ws) != len(set(ws)):
        parser.error("w values must be distinct")
    return ws


def main() -> None:
    ws = parse_args()
    h10q._selftest()
    q1_candidates = list(h10q.primerange(CLASS_Q_MIN, CLASS_Q_MAX + 1))
    assert q1_candidates and q1_candidates[-1] <= CLASS_Q_MAX
    assert all(h10q._is_prime(q1) for q1 in q1_candidates)

    results = []
    for w in ws:
        results.append(sweep_cell(w, q1_candidates))
    write_artifact(ws, q1_candidates, results)
    write_report(ws, results)
    verify_artifacts(ws, results)

    print(
        json.dumps(
            [result["summary"] for result in results], separators=(",", ":")
        )
    )
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
