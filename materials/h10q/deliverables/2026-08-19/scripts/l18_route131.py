#!/usr/bin/env python3
"""Break the canonical w=131 wall and replay the resulting exact closure.

The four requested extension families are run in order.  Every bounded probe
is written to data/l18_route131.jsonl; refusals are recorded but never used as
mathematical evidence.  The decisive route keeps the non-coprime factor
f=131 but replaces the canonical (a,tau)=(1,3/5) by the named L11 branch
(a,tau)=(7,0).  All primality and factorization decisions use h10q.py.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l13_filter  # noqa: E402

W = 131
TARGET_UT = (-1, 1)
TARGET_CELL = (W, TARGET_UT)
TARGET_Z = F(W) * F(*TARGET_UT)
CANONICAL_A = F(1)
CANONICAL_TAU = F(3, 5)
Q_PROBES = tuple(h10q.primerange(11, 98))
COMPOSITE_MULTIPLIERS = (3, 5, 7, 15, 21, 35, 105, 1181, 4279348769)
EXTRA_FROZEN_PRIMES = (11, 13, 1181, 4279348769)
COLUMN_UTS = ((-3, 1), (-2, 1), (-1, 3), (1, 3), (1, 1), (2, 1), (3, 1), (5, 1))
COLUMN_Q_PROBES = (19, 29, 41)
PACE_SECONDS = 0.1
OUT = ROOT / "data" / "l18_route131.jsonl"
REPORT = Path("/tmp/l18_route131.report")
SOURCE = "l18_route131.py"

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


def pace() -> None:
    time.sleep(PACE_SECONDS)


def emit(handle, record: dict) -> None:
    """Write deterministic compact JSON without sorting mixed place keys."""
    handle.write(json.dumps(record, separators=(",", ":")) + "\n")
    handle.flush()


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def character_mod_5(value: int) -> int:
    residue = value % 5
    assert residue != 0
    return 1 if pow(residue, 2, 5) == 1 else -1


@lru_cache(maxsize=None)
def integer_support(value: int) -> tuple[int, ...]:
    factors = h10q.factorint(abs(value)) if abs(value) != 1 else {}
    for prime in factors:
        assert h10q._is_prime(prime)
    return tuple(sorted(factors))


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
        "excluded": cert["excluded"],
    }


def detail_summary(status: str, detail) -> dict | None:
    if detail is None:
        return None
    if status == "cofactor-big" or status == "bad":
        emergent, (rn, rd), smooth = detail
        return {
            "small_bad_places": emergent,
            "remainder_numerator": str(rn),
            "remainder_denominator": str(rd),
            "smooth": smooth,
        }
    return {"detail": str(detail)}


def fixed_factor_class_cert(
    a: F,
    z: F,
    tau: F,
    eps: int,
    f: int,
    q1: int,
    extra_frozen: tuple[int, ...] = (),
) -> dict | None:
    """Exponent-lemma class certificate for b=eps*f*Q, fixed f arbitrary.

    This is l12_class.l12_class_cert with the two restrictions irrelevant to
    this experiment removed: a need not be 1 and f need not be prime.  Every
    prime of f is frozen.  The proof is unchanged: Taylor exponents freeze
    x and d at S, Q == q1 mod 4A freezes (A|Q), and sufficiently large class
    members have the recorded real sign.
    """
    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    if delta == 0:
        return None
    alpha = -delta * A
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    if a == 0 or z == 0 or D == 0 or f <= 0 or eps not in (1, -1):
        return None
    if h10q.vp(a, 2) != 0 or not h10q._is_prime(q1):
        return None

    f_support = set(integer_support(f))
    for prime in extra_frozen:
        assert h10q._is_prime(prime)
    S = sorted(
        {2, 3, 5, 7}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
        | f_support
        | set(extra_frozen)
    )
    if q1 in S:
        return None
    if any(h10q.vp(value, q1) != 0 for value in (delta, A, a, F(f))):
        return None

    s = (a - 1) / 2
    P = h10q._l10_P(a, Z, D, A, delta, s)
    b0 = F(eps * f * q1)
    ks = {prime: h10q._l10_exponent(P, b0, prime) for prime in S}
    N = 8
    for prime, exponent in ks.items():
        modulus = prime**exponent
        N = N * modulus // math.gcd(N, modulus)
    character_modulus = 4 * A.numerator * A.denominator
    N = N * character_modulus // math.gcd(N, character_modulus)
    if math.gcd(q1, N) != 1:
        return None

    c0 = h10q._sun_h(a, b0, Z)
    if c0 is None:
        return None
    M0 = 16 - delta * c0 * c0 - 32 * A * b0 * s * s
    if M0 == 0:
        return None
    x0, d0 = alpha * M0, alpha * 2 * b0
    syms = {prime: h10q.hilbert(x0, d0, prime) for prime in S}

    nonzero = [i for i, coefficient in enumerate(P) if coefficient != 0]
    lead = P[nonzero[-1]]
    Q0 = 1 + max(
        (abs(coefficient / lead) for coefficient in P[: nonzero[-1]]),
        default=F(0),
    )
    sign_x = (
        (1 if alpha > 0 else -1)
        * (1 if lead > 0 else -1)
        * eps ** nonzero[-1]
    )
    sign_d = 1 if alpha * eps > 0 else -1
    syms[h10q.OO] = -1 if sign_x < 0 and sign_d < 0 else 1
    syms["q1"] = h10q.legendre(A, q1)

    excluded = sorted(
        h10q._l10_supp(D)
        | h10q._l10_supp(z.numerator)
        | h10q._l10_supp(delta)
        | h10q._l10_supp(a)
        | f_support
    )
    return {
        "ok": all(value == 1 for value in syms.values()),
        "N": N,
        "S": S,
        "ks": ks,
        "syms": syms,
        "Q0": Q0,
        "f": f,
        "eps": eps,
        "q": q1,
        "excluded": excluded,
    }


def fixed_class_attempt(
    extension: str,
    a: F,
    z: F,
    tau: F,
    eps: int,
    f: int,
    q1: int,
    extra_frozen: tuple[int, ...] = (),
) -> tuple[dict, dict | None]:
    """Run and exactly label one fixed-factor class attempt."""
    pace()
    controlled = (
        {2, 3, 5, 7}
        | h10q._l10_supp(-(1 - (1 + 4 * a * a) * tau * tau) * (1 + 4 * a * a))
        | h10q._l10_supp(1 - (1 + 4 * a * a) * tau * tau)
        | set(integer_support(f))
        | set(extra_frozen)
    )
    record = {
        "type": "class-attempt",
        "label": "PROVED",
        "extension": extension,
        "family": "fixed-factor",
        "cell": [W, list(TARGET_UT)] if z == TARGET_Z else None,
        "z": str(z),
        "a": str(a),
        "tau": str(tau),
        "f": f,
        "eps": eps,
        "q1": q1,
        "q1_mod_5": q1 % 5,
        "character_5_q1": h10q.legendre(F(5), q1),
        "extra_frozen": list(extra_frozen),
    }
    if q1 in controlled:
        record.update(
            {
                "status": "ineligible:q1-in-controlled-set",
                "controlled_set": sorted(controlled),
                "mathematical_evidence": False,
            }
        )
        return record, None
    try:
        cert = fixed_factor_class_cert(a, z, tau, eps, f, q1, extra_frozen)
    except Exception as exc:
        record.update(
            {
                "status": refusal(exc),
                "mathematical_evidence": False,
            }
        )
        return record, None
    if cert is None:
        record.update(
            {
                "status": "refused:none",
                "mathematical_evidence": False,
            }
        )
        return record, None
    bad_places = [str(place) for place, value in cert["syms"].items() if value != 1]
    record.update(cert_summary(cert))
    record.update(
        {
            "status": "aligned" if cert["ok"] else "obstructed:symbol-minus-one",
            "bad_places": bad_places,
            "mathematical_evidence": True,
        }
    )
    return record, cert


def member_attempt(
    extension: str,
    family: str,
    cell: tuple[int, tuple[int, int]],
    a: F,
    z: F,
    tau: F,
    eps: int,
    f: int,
    q1: int,
    cert: dict,
) -> tuple[dict, tuple | None]:
    """Decide the k=0 member with the complete L13 cofactor ladder."""
    pace()
    Q = q1
    b = F(eps * f * Q)
    record = {
        "type": "member-attempt",
        "label": "PROVED",
        "extension": extension,
        "family": family,
        "cell": [cell[0], list(cell[1])],
        "a": str(a),
        "tau": str(tau),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(cert["N"]),
        "k": 0,
        "Q": str(Q),
        "b": str(b),
    }
    try:
        assert h10q._is_prime(Q)
        smooth_status, smooth_detail = l13_filter.smooth_emergent(a, z, tau, b)
    except Exception as exc:
        record.update(
            {
                "status": refusal(exc),
                "mathematical_evidence": False,
            }
        )
        return record, None

    record["smooth_status"] = smooth_status
    record["smooth_detail"] = detail_summary(smooth_status, smooth_detail)
    if smooth_status == "zero":
        record.update(
            {
                "status": "zero",
                "rung": "square",
                "mathematical_evidence": True,
            }
        )
        return record, (0, Q, "square", b, None)
    if smooth_status != "cofactor-big":
        record.update(
            {
                "status": smooth_status,
                "mathematical_evidence": not smooth_status.startswith("refused:"),
            }
        )
        return record, None

    try:
        verdict, info = l13_filter.cofactor_decide(
            a,
            z,
            tau,
            b,
            detail=smooth_detail,
            deadline=time.time() + 600,
        )
    except Exception as exc:
        verdict, info = refusal(exc), None
    record["status"] = verdict
    record["cofactor_info"] = info
    record["rung"] = info.get("rung") if isinstance(info, dict) else None
    record["mathematical_evidence"] = not verdict.startswith("refused:")
    if verdict == "zero":
        return record, (0, Q, record["rung"], b, info)
    return record, None


def canonical_proof() -> dict:
    """Universal p=5/wild reciprocity identity for the canonical protocol."""
    z3 = TARGET_Z**3
    D = 1 - z3 - z3 * z3
    assert h10q.vp(z3, 5) == 0
    assert h10q.vp(D, 5) == 0
    assert h10q.unit_mod(z3, 5) == 4
    assert h10q.unit_mod(D, 5) == 1
    w_character = h10q.legendre(F(W), 5)
    assert w_character == 1
    residue_table = []
    for eps in (1, -1):
        assert character_mod_5(eps) == 1
        for residue in (1, 2, 3, 4):
            q_character = character_mod_5(residue)
            hilbert_5 = -w_character * q_character
            wild = q_character
            residue_table.append(
                {
                    "eps": eps,
                    "q1_mod_5": residue,
                    "character_5_q1": q_character,
                    "hilbert_5": hilbert_5,
                    "hilbert_q1": wild,
                    "product": hilbert_5 * wild,
                }
            )
            assert hilbert_5 * wild == -1
    return {
        "type": "canonical-reciprocity-theorem",
        "label": "PROVED",
        "extension": "canonical-baseline",
        "cell": [W, list(TARGET_UT)],
        "a": "1",
        "tau": "3/5",
        "f": W,
        "valuation_proof": {
            "v5_c": -1,
            "v5_M": -3,
            "v5_x": -3,
            "v5_d": 0,
            "z_cubed_mod_5": 4,
            "D_mod_5": 1,
        },
        "identity": "Hilb_5=-(131|5)(5|q1)=-(5|q1); Hilb_q1=(5|q1)",
        "canonical_general_criterion": "simultaneous +1 iff (w|5)=-1 and (5|q1)=+1",
        "residue_table": residue_table,
        "conclusion": "the canonical a=1,tau=3/5,f=131 protocol is reciprocity-obstructed",
    }


def run_canonical_baseline(handle) -> dict:
    emit(handle, canonical_proof())
    branch_counts: Counter[int] = Counter()
    pair_counts: Counter[str] = Counter()
    attempts = 0
    for eps in (1, -1):
        for q1 in Q_PROBES:
            record, cert = fixed_class_attempt(
                "canonical-baseline",
                CANONICAL_A,
                TARGET_Z,
                CANONICAL_TAU,
                eps,
                W,
                q1,
            )
            assert cert is not None
            q_character = h10q.legendre(F(5), q1)
            assert cert["syms"][5] == -q_character
            assert cert["syms"]["q1"] == q_character
            assert not cert["ok"]
            record["status"] = "obstructed:p5-versus-wild"
            record["hilbert_5"] = cert["syms"][5]
            record["hilbert_q1"] = cert["syms"]["q1"]
            record["symbol_product"] = -1
            emit(handle, record)
            attempts += 1
            branch_counts[q_character] += 1
            pair_counts[f"{cert['syms'][5]},{cert['syms']['q1']}"] += 1
    assert set(branch_counts) == {-1, 1}
    summary = {
        "type": "family-summary",
        "label": "PROVED",
        "extension": "canonical-baseline",
        "status": "universally-obstructed",
        "attempts": attempts,
        "q1_range": [Q_PROBES[0], Q_PROBES[-1]],
        "character_branch_attempts": {str(k): v for k, v in sorted(branch_counts.items())},
        "symbol_pair_attempts": dict(sorted(pair_counts.items())),
        "universal_reason": "Hilb_5*Hilb_q1=-1 for every admissible q1",
    }
    emit(handle, summary)
    return summary


def run_composite_f(handle) -> dict:
    attempts = 0
    aligned = 0
    refusals = 0
    per_multiplier: dict[str, Counter[str]] = {}
    for multiplier in COMPOSITE_MULTIPLIERS:
        f = W * multiplier
        per_multiplier[str(multiplier)] = Counter()
        for eps in (1, -1):
            for q1 in Q_PROBES:
                record, cert = fixed_class_attempt(
                    "a-composite-f",
                    CANONICAL_A,
                    TARGET_Z,
                    CANONICAL_TAU,
                    eps,
                    f,
                    q1,
                )
                record["multiplier"] = multiplier
                emit(handle, record)
                attempts += 1
                per_multiplier[str(multiplier)][record["status"]] += 1
                if record["status"] == "aligned":
                    aligned += 1
                if record["status"].startswith("refused:"):
                    refusals += 1
    assert aligned == 0
    summary = {
        "type": "family-summary",
        "label": "EVIDENCE",
        "extension": "a-composite-f",
        "status": "no-aligned-class-in-bounded-grid",
        "multipliers": list(COMPOSITE_MULTIPLIERS),
        "q1_range": [Q_PROBES[0], Q_PROBES[-1]],
        "eps": [1, -1],
        "attempts": attempts,
        "aligned": aligned,
        "refusals": refusals,
        "per_multiplier_status": {
            key: dict(sorted(counts.items())) for key, counts in per_multiplier.items()
        },
        "scope": "finite exact scan; not a universal negative",
        "refusals_used_as_evidence": False,
    }
    emit(handle, summary)
    return summary


def run_enlarged_frozen(handle) -> dict:
    emit(
        handle,
        {
            "type": "enlarged-set-theorem",
            "label": "PROVED",
            "extension": "b-enlarged-frozen-set",
            "cell": [W, list(TARGET_UT)],
            "statement": "No superset of {2,3,5,7,131} can align the canonical family.",
            "proof": "Enlarging S only restricts Q and freezes more unchanged pointwise symbols; Hilb_5*Hilb_q1=-1 remains an existing pair in every subclass.",
        },
    )
    attempts = 0
    eligible = 0
    ineligible = 0
    refusals = 0
    for prime in EXTRA_FROZEN_PRIMES:
        assert h10q._is_prime(prime)
        for eps in (1, -1):
            for q1 in Q_PROBES:
                record, cert = fixed_class_attempt(
                    "b-enlarged-frozen-set",
                    CANONICAL_A,
                    TARGET_Z,
                    CANONICAL_TAU,
                    eps,
                    W,
                    q1,
                    (prime,),
                )
                record["added_prime"] = prime
                if cert is not None:
                    q_character = h10q.legendre(F(5), q1)
                    assert cert["syms"][5] == -q_character
                    assert cert["syms"]["q1"] == q_character
                    assert not cert["ok"]
                    record["status"] = "obstructed:p5-versus-wild"
                    record["hilbert_5"] = cert["syms"][5]
                    record["hilbert_q1"] = cert["syms"]["q1"]
                    eligible += 1
                elif record["status"].startswith("ineligible:"):
                    ineligible += 1
                else:
                    refusals += 1
                emit(handle, record)
                attempts += 1
    summary = {
        "type": "family-summary",
        "label": "PROVED",
        "extension": "b-enlarged-frozen-set",
        "status": "universally-obstructed",
        "added_primes_probed": list(EXTRA_FROZEN_PRIMES),
        "attempts": attempts,
        "eligible": eligible,
        "ineligible": ineligible,
        "refusals": refusals,
        "reason": "a stricter congruence class cannot alter the canonical p=5/wild reciprocity identity",
        "refusals_used_as_evidence": False,
    }
    emit(handle, summary)
    return summary


def close_target(row: dict, cert: dict, hit: tuple) -> tuple[dict, dict]:
    """Independently replay the exact noncanonical target closure."""
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    f = int(row["f"])
    q1 = int(row["q1"])
    replayed_cert = fixed_factor_class_cert(a, TARGET_Z, tau, eps, f, q1)
    assert replayed_cert == cert and replayed_cert["ok"]
    assert Q == q1 + k * replayed_cert["N"]
    assert h10q._is_prime(Q) and h10q._is_prime(f)
    assert b == F(eps * f * Q)

    member_syms = h10q._l12_syms(a, TARGET_Z, tau, b)
    aligned_syms = l13_filter.aligned(a, TARGET_Z, tau, b)
    assert member_syms is not None and all(value == 1 for value in member_syms.values())
    assert aligned_syms is not None and all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, TARGET_Z, tau, b)
    assert smooth_status == "cofactor-big"
    verdict, info = l13_filter.cofactor_decide(
        a, TARGET_Z, tau, b, detail=smooth_detail, deadline=None
    )
    assert verdict == "zero" and info == prior_info and info["rung"] == rung
    assert len(info["places"]) == 1
    assert info["places"][0][1:] == [1, "proved", "prime"]
    cofactor_prime = int(info["places"][0][0])
    assert h10q._is_prime(cofactor_prime)

    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z3 = TARGET_Z**3
    D = 1 - Z3 - a * a * Z3 * Z3
    s = (a - 1) / 2
    c0 = h10q._sun_h(a, b, Z3)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    P = h10q._l10_P(a, Z3, D, A, delta, s)
    Pb = sum(coefficient * b**i for i, coefficient in enumerate(P))
    assert M0 == Pb / (b**4 * D * D * A * A)

    tied = h10q._l7_tied_status(a, b, TARGET_Z, tau)
    ramified = h10q.ramified(x0, d0)
    assert tied is True and ramified == []
    closure = {
        "family": "ESC",
        "cell": [W, list(TARGET_UT)],
        "a": str(a),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(replayed_cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": str(tied),
        "ramified_empty": not ramified,
        "source": SOURCE,
    }
    assert tuple(closure) == CLOSURE_FIELDS
    replay = {
        "type": "kernel-replay",
        "label": "PROVED",
        "extension": "c-noncanonical-tau",
        "cell": [W, list(TARGET_UT)],
        "z": str(TARGET_Z),
        "a": str(a),
        "A": str(A),
        "tau": str(tau),
        "delta": str(delta),
        "alpha": str(alpha),
        "f": f,
        "eps": eps,
        "q1": q1,
        "N": str(replayed_cert["N"]),
        "k": k,
        "Q": str(Q),
        "b": str(b),
        "class_cert_replayed": True,
        "class_syms": replayed_cert["syms"],
        "member_syms": member_syms,
        "aligned_syms": aligned_syms,
        "smooth_status": smooth_status,
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "cofactor_prime": str(cofactor_prime),
        "cofactor_prime_proven": True,
        "D": str(D),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
        "tied": tied,
        "ramified": ramified,
        "ramified_empty": True,
        "refusals_used_as_evidence": False,
    }
    return closure, replay


def run_noncanonical_tau(handle) -> tuple[dict, dict, dict]:
    roots = [residue for residue in range(1, W) if (1 + 4 * residue * residue) % W == 0]
    minus_one_character = h10q.legendre(F(-1), W)
    assert roots == [] and minus_one_character == -1
    emit(
        handle,
        {
            "type": "hensel-condition-theorem",
            "label": "PROVED",
            "extension": "c-noncanonical-tau",
            "cell": [W, list(TARGET_UT)],
            "canonical_A_zero_roots": roots,
            "character_minus_one": minus_one_character,
            "canonical_reason": "A=1+4a^2 cannot vanish mod 131 because -1 is a nonsquare",
            "replacement_condition": "W1: (A|131)=-1 and (-delta*A|131)=+1",
            "named_branches": ["tau=0", "tau=2a/A", "tau=1", "tau=(1+2a^2)/A"],
        },
    )

    alternate_pool = []
    for residue in range(W):
        pace()
        a_int = residue if residue % 2 else residue + W
        A = 1 + 4 * a_int * a_int
        character_A = h10q.legendre(F(A), W)
        accepted = character_A == -1
        record = {
            "type": "hensel-residue-attempt",
            "label": "PROVED",
            "extension": "c-noncanonical-tau",
            "residue_mod_131": residue,
            "odd_lift_a": a_int,
            "A": A,
            "character_A": character_A,
            "status": "admissible:W1-A-nonsquare" if accepted else "rejected:A-square",
        }
        emit(handle, record)
        if accepted:
            alternate_pool.append(a_int)
    alternate_pool = sorted(set(alternate_pool))
    assert alternate_pool

    q1 = 41
    closure = None
    replay = None
    class_attempts = 0
    member_attempts = 0
    refusals = 0
    for a_int in alternate_pool:
        a = F(a_int)
        A = 1 + 4 * a * a
        branches = (
            ("tau0", F(0)),
            ("canonical", 2 * a / A),
            ("tau1", F(1)),
            ("square", (1 + 2 * a * a) / A),
        )
        for branch, tau in branches:
            alpha = -(1 - A * tau * tau) * A
            target_character = (
                0 if h10q.vp(alpha, W) != 0 else h10q.legendre(alpha, W)
            )
            for eps in (1, -1):
                if target_character != 1:
                    pace()
                    record = {
                        "type": "class-attempt",
                        "label": "PROVED",
                        "extension": "c-noncanonical-tau",
                        "family": branch,
                        "cell": [W, list(TARGET_UT)],
                        "z": str(TARGET_Z),
                        "a": str(a),
                        "A": str(A),
                        "tau": str(tau),
                        "alpha": str(alpha),
                        "target_character": target_character,
                        "f": W,
                        "eps": eps,
                        "q1": q1,
                        "status": "ineligible:target-character",
                        "mathematical_evidence": False,
                    }
                    emit(handle, record)
                    class_attempts += 1
                    continue
                record, cert = fixed_class_attempt(
                    "c-noncanonical-tau", a, TARGET_Z, tau, eps, W, q1
                )
                record["family"] = branch
                record["A"] = str(A)
                record["alpha"] = str(alpha)
                record["target_character"] = target_character
                emit(handle, record)
                class_attempts += 1
                if record["status"].startswith("refused:"):
                    refusals += 1
                if cert is None or not cert["ok"]:
                    continue
                member_record, hit = member_attempt(
                    "c-noncanonical-tau",
                    branch,
                    TARGET_CELL,
                    a,
                    TARGET_Z,
                    tau,
                    eps,
                    W,
                    q1,
                    cert,
                )
                emit(handle, member_record)
                member_attempts += 1
                if member_record["status"].startswith("refused:"):
                    refusals += 1
                if hit is None:
                    continue
                closure, replay = close_target(record, cert, hit)
                emit(handle, replay)
                emit(handle, closure)
                break
            if closure is not None:
                break
        if closure is not None:
            break
    assert closure is not None and replay is not None
    assert closure["cell"] == [W, list(TARGET_UT)]
    assert closure["a"] == "7" and closure["q1"] == 41
    summary = {
        "type": "family-summary",
        "label": "PROVED",
        "extension": "c-noncanonical-tau",
        "status": "VERIFIED-target-closure",
        "canonical_A_zero_roots": roots,
        "alternate_W1_pool_size": len(alternate_pool),
        "alternate_W1_pool_first": alternate_pool[:12],
        "class_attempts_before_closure": class_attempts,
        "member_attempts_before_closure": member_attempts,
        "refusals": refusals,
        "closure": closure,
        "reason_no_L10_conflict": "b=131*41 shares 131 with z, so L10/L12 coprime-b walls do not apply",
        "refusals_used_as_evidence": False,
    }
    emit(handle, summary)
    return closure, replay, summary


def builtin_class_attempt(
    extension: str,
    cell: tuple[int, tuple[int, int]],
    a: F,
    z: F,
    tau: F,
    eps: int,
    q1: int,
) -> tuple[dict, dict | None]:
    pace()
    record = {
        "type": "class-attempt",
        "label": "PROVED",
        "extension": extension,
        "family": "L11",
        "cell": [cell[0], list(cell[1])],
        "z": str(z),
        "a": str(a),
        "tau": str(tau),
        "eps": eps,
        "f": 1,
        "q1": q1,
    }
    try:
        cert = h10q._l10_class_cert(a, z, tau, eps, q1)
    except Exception as exc:
        record.update({"status": refusal(exc), "mathematical_evidence": False})
        return record, None
    if cert is None:
        record.update({"status": "refused:none", "mathematical_evidence": False})
        return record, None
    record.update(cert_summary({**cert, "excluded": []}))
    record["status"] = "aligned" if cert["ok"] else "obstructed:symbol-minus-one"
    record["bad_places"] = [
        str(place) for place, value in cert["syms"].items() if value != 1
    ]
    record["mathematical_evidence"] = True
    return record, cert


def close_column_hit(
    cell: tuple[int, tuple[int, int]],
    z: F,
    row: dict,
    cert: dict,
    hit: tuple,
) -> tuple[dict, dict]:
    k, Q, rung, b, prior_info = hit
    a = F(row["a"])
    tau = F(row["tau"])
    eps = int(row["eps"])
    q1 = int(row["q1"])
    replayed_cert = h10q._l10_class_cert(a, z, tau, eps, q1)
    assert replayed_cert == cert and replayed_cert["ok"]
    assert k == 0 and Q == q1 and h10q._is_prime(Q) and b == F(eps * Q)
    aligned_syms = l13_filter.aligned(a, z, tau, b)
    assert aligned_syms is not None and all(value == 1 for value in aligned_syms.values())
    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, z, tau, b)
    assert smooth_status == "cofactor-big"
    verdict, info = l13_filter.cofactor_decide(a, z, tau, b, detail=smooth_detail)
    assert verdict == "zero" and info == prior_info and info["rung"] == rung
    assert info["places"] and all(
        place[1] == 1 and place[2] == "proved" for place in info["places"]
    )

    A = 1 + 4 * a * a
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z3 = z**3
    D = 1 - Z3 - a * a * Z3 * Z3
    s = (a - 1) / 2
    c0 = h10q._sun_h(a, b, Z3)
    M0 = 16 - delta * c0 * c0 - 32 * A * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    P = h10q._l10_P(a, Z3, D, A, delta, s)
    assert M0 == sum(c * b**i for i, c in enumerate(P)) / (
        b**4 * D * D * A * A
    )
    try:
        tied = str(h10q._l7_tied_status(a, b, z, tau))
    except Exception as exc:
        tied = refusal(exc)
    try:
        ramified_empty: bool | str = h10q.ramified(x0, d0) == []
    except Exception as exc:
        ramified_empty = refusal(exc)

    closure = {
        "family": "L11",
        "cell": [cell[0], list(cell[1])],
        "a": str(a),
        "eps": eps,
        "f": 1,
        "q1": q1,
        "N": str(cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": tied,
        "ramified_empty": ramified_empty,
        "source": SOURCE,
    }
    assert tuple(closure) == CLOSURE_FIELDS
    replay = {
        "type": "kernel-replay",
        "label": "PROVED",
        "extension": "d-other-u",
        "cell": [cell[0], list(cell[1])],
        "z": str(z),
        "a": str(a),
        "tau": str(tau),
        "b": str(b),
        "q1": q1,
        "N": str(cert["N"]),
        "class_cert_replayed": True,
        "aligned_syms": aligned_syms,
        "smooth_status": smooth_status,
        "emergent_verdict": verdict,
        "cofactor_info": info,
        "complete_l13_zero_certificate": True,
        "D": str(D),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
        "tied_crosscheck": tied,
        "ramified_crosscheck": ramified_empty,
        "crosscheck_refusals_used_as_evidence": False,
    }
    return closure, replay


def run_other_u(handle) -> tuple[dict | None, dict]:
    attempts = 0
    class_attempts = 0
    member_attempts = 0
    partial_closure = None
    for ut in COLUMN_UTS:
        pace()
        cell = (W, ut)
        z = F(W) * F(*ut)
        numerator_factors = h10q.factorint(abs(z.numerator))
        for prime in numerator_factors:
            assert h10q._is_prime(prime)
        escape_primes = h10q._l11_escape_primes(z)
        record = {
            "type": "column-cell-attempt",
            "label": "PROVED",
            "extension": "d-other-u",
            "cell": [W, list(ut)],
            "z": str(z),
            "numerator_factorization": numerator_factors,
            "escape_primes_1_mod_4": escape_primes,
            "status": "L11-candidate" if escape_primes else "proved-no-L11-root",
        }
        emit(handle, record)
        attempts += 1
        if not escape_primes:
            continue

        prime = min(escape_primes)
        roots = [r for r in range(1, prime) if (1 + 4 * r * r) % prime == 0]
        assert roots
        a_pool = sorted({r if r % 2 else r + prime for r in roots})
        a = F(a_pool[0])
        A = 1 + 4 * a * a
        tau = (1 + 2 * a * a) / A
        for eps in (1, -1):
            for q1 in COLUMN_Q_PROBES:
                class_record, cert = builtin_class_attempt(
                    "d-other-u", cell, a, z, tau, eps, q1
                )
                class_record["root_prime"] = prime
                class_record["roots"] = roots
                emit(handle, class_record)
                class_attempts += 1
                if cert is None or not cert["ok"]:
                    continue
                member_record, hit = member_attempt(
                    "d-other-u", "L11", cell, a, z, tau, eps, 1, q1, cert
                )
                emit(handle, member_record)
                member_attempts += 1
                if hit is None:
                    continue
                partial_closure, replay = close_column_hit(cell, z, class_record, cert, hit)
                emit(handle, replay)
                emit(handle, partial_closure)
                break
            if partial_closure is not None:
                break
        if partial_closure is not None:
            break
    assert partial_closure is not None
    assert partial_closure["cell"] == [W, [5, 1]]
    summary = {
        "type": "family-summary",
        "label": "PROVED",
        "extension": "d-other-u",
        "status": "VERIFIED-partial-column-extension",
        "cells_attempted": attempts,
        "class_attempts": class_attempts,
        "member_attempts": member_attempts,
        "closure": partial_closure,
        "closure_basis": "complete L13 zero certificate; independent crosscheck refusals, if any, are not evidence",
    }
    emit(handle, summary)
    return partial_closure, summary


def run_w179_generalization(handle) -> dict:
    """Bounded extra probe requested by the lead; no negative inference."""
    w = 179
    cell = (w, (-1, 1))
    z = F(-w)
    a = F(7)
    tau = F(0)
    q_candidates = tuple(q for q in Q_PROBES if h10q.legendre(F(5), q) == 1)
    attempts = 0
    aligned = 0
    member_attempts = 0
    refusals = 0
    hit = None
    for eps in (1, -1):
        for q1 in q_candidates:
            record, cert = fixed_class_attempt(
                "bounded-generalization-w179", a, z, tau, eps, w, q1
            )
            record["cell"] = [w, [-1, 1]]
            emit(handle, record)
            attempts += 1
            if cert is None or not cert["ok"]:
                if record["status"].startswith("refused:"):
                    refusals += 1
                continue
            aligned += 1
            member_record, member_hit = member_attempt(
                "bounded-generalization-w179",
                "tau0-fixed-factor",
                cell,
                a,
                z,
                tau,
                eps,
                w,
                q1,
                cert,
            )
            emit(handle, member_record)
            member_attempts += 1
            if member_record["status"].startswith("refused:"):
                refusals += 1
            if member_hit is not None:
                hit = member_hit
                break
        if hit is not None:
            break
    assert hit is None
    summary = {
        "type": "family-summary",
        "label": "EVIDENCE",
        "extension": "bounded-generalization-w179",
        "status": "no-closure-in-bounded-probe",
        "cell": [w, [-1, 1]],
        "a": "7",
        "tau": "0",
        "f": w,
        "q1_range": [Q_PROBES[0], Q_PROBES[-1]],
        "q1_filter": "(5|q1)=+1",
        "eps": [1, -1],
        "class_attempts": attempts,
        "aligned_classes": aligned,
        "member_attempts": member_attempts,
        "refusals": refusals,
        "scope": "bounded k=0 evidence only; refusals carry no evidentiary weight",
        "refusals_used_as_evidence": False,
    }
    emit(handle, summary)
    return summary


def write_report(
    canonical: dict,
    composite: dict,
    enlarged: dict,
    target_closure: dict,
    noncanonical: dict,
    partial_closure: dict | None,
    other_u: dict,
    generalization: dict,
) -> None:
    lines = [
        "STATUS: VERIFIED exact target closure",
        "LABEL: PROVED",
        "TARGET: cell [131,[-1,1]], z=-131",
        "ENGINE: math/h10q/h10q.py proven-primality engine only; refusals are never evidence",
        "CANONICAL PROTOCOL — PROVED: Hilb_5=-(5|q1), Hilb_q1=(5|q1); the two required symbols are anti-correlated.",
        "CANONICAL BRANCH COUNTS — PROVED: "
        + json.dumps(canonical, sort_keys=True),
        "(a) COMPOSITE f — EVIDENCE: " + json.dumps(composite, sort_keys=True),
        "(b) ENLARGED S — PROVED negative: " + json.dumps(enlarged, sort_keys=True),
        "(c) NONCANONICAL tau — PROVED exact closure: "
        + json.dumps(target_closure, sort_keys=True),
        "ROUTE: a=7, A=197, tau=0, f=131, eps=+1, q1=Q=41, b=5371; class N=45046449133189418880.",
        "CONSISTENCY: b shares 131 with z, so this uses the L12b non-coprimality door and does not contradict the coprime-b tau=0 wall.",
        "NONCANONICAL SEARCH — PROVED: " + json.dumps(noncanonical, sort_keys=True),
        "(d) OTHER u — PROVED partial column extension: "
        + json.dumps(other_u, sort_keys=True),
        "PARTIAL CLOSURE ROW: " + json.dumps(partial_closure, sort_keys=True),
        "w=179 GENERALIZATION — EVIDENCE only: "
        + json.dumps(generalization, sort_keys=True),
        "CONCLUSION: the w=131 obstruction is a theorem about the canonical a=1,tau=3/5 protocol, not an intrinsic obstruction of the cell.",
        "ARTIFACT: " + str(OUT),
        "REPLAY: python3 math/h10q/l18_route131.py",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    h10q._selftest()
    assert TARGET_Z == F(-131)
    assert h10q._is_prime(W)
    assert h10q.factorint(abs(TARGET_Z.numerator)) == {131: 1}
    target_D = 1 - TARGET_Z**3 - TARGET_Z**6
    assert h10q.factorint(abs(target_D.numerator)) == {1181: 1, 4279348769: 1}
    assert set(h10q.legendre(F(5), q) for q in Q_PROBES) == {-1, 1}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as handle:
        emit(
            handle,
            {
                "type": "meta",
                "label": "PROVED",
                "schema": "l18-route131-v1",
                "task": "four-family route beyond the canonical w=131 protocol wall",
                "engine": "math/h10q/h10q.py proven-primality only",
                "source_script": "math/h10q/l18_route131.py",
                "target_cell": [W, list(TARGET_UT)],
                "z": str(TARGET_Z),
                "pace_seconds": PACE_SECONDS,
                "refusals_are_not_evidence": True,
            },
        )
        canonical = run_canonical_baseline(handle)
        composite = run_composite_f(handle)
        enlarged = run_enlarged_frozen(handle)
        target_closure, _, noncanonical = run_noncanonical_tau(handle)
        partial_closure, other_u = run_other_u(handle)
        generalization = run_w179_generalization(handle)
        emit(
            handle,
            {
                "type": "summary",
                "label": "PROVED",
                "status": "VERIFIED",
                "target_cell": [W, list(TARGET_UT)],
                "target_closure": target_closure,
                "partial_column_closure": partial_closure,
                "canonical_protocol_status": "PROVED-OBSTRUCTED",
                "intrinsic_cell_status": "PROVED-SOLUBLE",
                "conclusion": "canonical obstruction is protocol-specific, not intrinsic",
                "refusals_used_as_evidence": False,
            },
        )

    write_report(
        canonical,
        composite,
        enlarged,
        target_closure,
        noncanonical,
        partial_closure,
        other_u,
        generalization,
    )
    print("STATUS VERIFIED")
    print("TARGET CLOSURE", json.dumps(target_closure, sort_keys=True))
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
