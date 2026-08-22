#!/usr/bin/env python3
"""Search and exactly replay the noncanonical fixed-factor route at w=179.

The class constructor is the fixed-factor exponent-lemma generalization already
used by l18_route131.py: compared with h10q._l10_class_cert, odd ``a`` is not
fixed to 1 and ``b = eps*f*Q`` has every prime of the fixed factor ``f`` added
to the controlled set.  Every emitted aligned-class row names that
generalization explicitly.  Member decisions use only the proven h10q kernel
and the L13 smooth/cofactor ladder.  Refusals are recorded, never interpreted
as negative evidence.
"""
from __future__ import annotations

import json
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import h10q
import l13_filter
from l18_route131 import fixed_factor_class_cert

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "l19_cell179.jsonl"
REPORT = Path("/tmp/l19_cell179.md")
SOURCE = "l19_cell179.py"

W = 179
UT = (-1, 1)
CELL = (W, UT)
Z = F(-W)
TAU = F(0)
FIXED_FACTOR = W
A_VALUES = (7, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25)
EPS_VALUES = (1, -1)
Q1_MAX = 400
K_MAX = 60
PACE_EVERY = 200
PACE_SECONDS = 0.01

CERT_FUNCTION = "l18_route131.fixed_factor_class_cert"
CERT_GENERALIZATION = (
    "local exponent-lemma generalization of h10q._l10_class_cert: arbitrary "
    "odd a and fixed b-factor f, with every prime of f frozen"
)
HYPOTHESIS_CLAUSE = (
    "For tau=0 and b=eps*179*Q, the generalized fixed-factor class "
    "certificate has every controlled, moving-prime, and asymptotic-real "
    "symbol +1; an eligible proven prime Q=q1+kN (0<=k<=60) has "
    "l13_filter.aligned all +1 and the complete smooth_emergent plus "
    "cofactor_decide ladder returns verdict zero; acceptance additionally "
    "requires h10q._l7_tied_status is True and h10q.ramified is empty."
)

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


class Pacer:
    """Pace the two candidate loops without sleeping on every cheap row."""

    def __init__(self) -> None:
        self.candidates = 0

    def tick(self) -> None:
        self.candidates += 1
        if self.candidates % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)


def emit(handle, record: dict) -> None:
    handle.write(json.dumps(record, separators=(",", ":")) + "\n")
    handle.flush()


def refusal_name(exc: BaseException) -> str:
    return type(exc).__name__


def cert_payload(cert: dict) -> dict:
    return {
        "N": str(cert["N"]),
        "S": cert["S"],
        "ks": cert["ks"],
        "syms": cert["syms"],
        "Q0": str(cert["Q0"]),
        "excluded": cert["excluded"],
        "all_symbols_plus_one": bool(cert["ok"]),
    }


def controlled_set(a: F) -> set[int]:
    A = 1 + 4 * a * a
    delta = 1 - A * TAU * TAU
    alpha = -delta * A
    return (
        {2, 3, 5, 7}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
        | set(h10q.factorint(FIXED_FACTOR))
    )


def class_none_reason(a: F, q1: int) -> str:
    controlled = controlled_set(a)
    if q1 in controlled:
        return "q1-in-controlled-set"
    A = 1 + 4 * a * a
    if any(h10q.vp(value, q1) != 0 for value in (F(1), A, a, F(FIXED_FACTOR))):
        return "q1-divides-fixed-data"
    return "cert-none"


def enumerate_classes(handle, q1_candidates: tuple[int, ...], pacer: Pacer):
    aligned = []
    branch_counts: dict[tuple[int, int], Counter] = defaultdict(Counter)
    total_counts: Counter = Counter()
    refused = []

    for a_int in A_VALUES:
        a = F(a_int)
        for eps in EPS_VALUES:
            branch = branch_counts[(a_int, eps)]
            for q1 in q1_candidates:
                pacer.tick()
                base = {
                    "type": "class-attempt",
                    "label": "PROVED",
                    "cell": [W, list(UT)],
                    "family": "ESC-NONCANONICAL",
                    "a": str(a),
                    "tau": str(TAU),
                    "eps": eps,
                    "f": FIXED_FACTOR,
                    "q1": q1,
                    "q1_proven_prime": True,
                    "q1_max_clause": q1 <= Q1_MAX,
                    "hypothesis_clause": HYPOTHESIS_CLAUSE,
                    "certificate_function": CERT_FUNCTION,
                    "certificate_generalization": CERT_GENERALIZATION,
                }
                try:
                    cert = fixed_factor_class_cert(
                        a, Z, TAU, eps, FIXED_FACTOR, q1
                    )
                except (h10q.FactorBudget, h10q.PrimalityBound) as exc:
                    reason = refusal_name(exc)
                    base.update(
                        {
                            "label": "OPEN",
                            "status": "refused:" + reason,
                            "reason_code": reason,
                            "mathematical_evidence": False,
                        }
                    )
                    refused.append(
                        {"phase": "class", "a": a_int, "eps": eps, "q1": q1,
                         "reason": reason}
                    )
                    branch[reason] += 1
                    total_counts[reason] += 1
                    emit(handle, base)
                    continue

                if cert is None:
                    reason = class_none_reason(a, q1)
                    base.update(
                        {
                            "status": "ineligible:" + reason,
                            "reason_code": reason,
                            "mathematical_evidence": True,
                        }
                    )
                    branch[reason] += 1
                    total_counts[reason] += 1
                    emit(handle, base)
                    continue

                base.update(cert_payload(cert))
                if not cert["ok"]:
                    bad_places = [
                        str(place)
                        for place, symbol in cert["syms"].items()
                        if symbol != 1
                    ]
                    base.update(
                        {
                            "status": "obstructed:symbol-minus-one",
                            "reason_code": "symbol-minus-one",
                            "bad_places": bad_places,
                            "mathematical_evidence": True,
                        }
                    )
                    branch["symbol-minus-one"] += 1
                    total_counts["symbol-minus-one"] += 1
                    emit(handle, base)
                    continue

                base.update(
                    {
                        "status": "aligned",
                        "reason_code": "all-symbols-plus-one",
                        "mathematical_evidence": True,
                    }
                )
                branch["aligned"] += 1
                total_counts["aligned"] += 1
                emit(handle, base)
                aligned.append(
                    {
                        "a": a,
                        "eps": eps,
                        "f": FIXED_FACTOR,
                        "q1": q1,
                        "tau": TAU,
                        "cert": cert,
                    }
                )

    aligned.sort(key=lambda row: (int(row["cert"]["N"]), int(row["a"]),
                                  -int(row["eps"]), int(row["q1"])))
    for class_index, row in enumerate(aligned):
        row["class_order_index"] = class_index
        emit(
            handle,
            {
                "type": "aligned-class",
                "label": "PROVED",
                "cell": [W, list(UT)],
                "family": "ESC-NONCANONICAL",
                "class_order_index": class_index,
                "class_order_key": [str(row["cert"]["N"]), str(row["a"]),
                                    row["eps"], row["q1"]],
                "a": str(row["a"]),
                "tau": str(TAU),
                "eps": row["eps"],
                "f": row["f"],
                "q1": row["q1"],
                "certificate_function": CERT_FUNCTION,
                "certificate_generalization": CERT_GENERALIZATION,
                "hypothesis_clause": HYPOTHESIS_CLAUSE,
                **cert_payload(row["cert"]),
            },
        )

    branch_rows = []
    for a_int in A_VALUES:
        for eps in EPS_VALUES:
            counts = branch_counts[(a_int, eps)]
            branch_rows.append(
                {
                    "a": a_int,
                    "eps": eps,
                    "q1_candidates": len(q1_candidates),
                    "counts": dict(sorted(counts.items())),
                    "refused": sum(
                        counts[reason] for reason in ("FactorBudget", "PrimalityBound")
                    ),
                    "decided_fraction": (
                        len(q1_candidates)
                        - sum(counts[reason] for reason in ("FactorBudget", "PrimalityBound"))
                    )
                    / len(q1_candidates),
                }
            )
    emit(
        handle,
        {
            "type": "class-pool-summary",
            "label": "PROVED" if not refused else "OPEN",
            "cell": [W, list(UT)],
            "hypothesis_clause": HYPOTHESIS_CLAUSE,
            "a_values": list(A_VALUES),
            "tau": str(TAU),
            "eps_values": list(EPS_VALUES),
            "f": FIXED_FACTOR,
            "q1_bound": [2, Q1_MAX],
            "q1_candidates": list(q1_candidates),
            "q1_count": len(q1_candidates),
            "attempts": len(A_VALUES) * len(EPS_VALUES) * len(q1_candidates),
            "aligned_classes": len(aligned),
            "reason_counts": dict(sorted(total_counts.items())),
            "per_branch_counts": branch_rows,
            "class_refusals": refused,
            "refusals_used_as_evidence": False,
        },
    )
    return aligned, branch_rows, total_counts, refused


def member_base(row: dict, k: int, Q: int) -> dict:
    cert = row["cert"]
    return {
        "type": "member-attempt",
        "label": "PROVED",
        "cell": [W, list(UT)],
        "family": "ESC-NONCANONICAL",
        "class_order_index": row["class_order_index"],
        "a": str(row["a"]),
        "tau": str(row["tau"]),
        "eps": row["eps"],
        "f": row["f"],
        "q1": row["q1"],
        "N": str(cert["N"]),
        "k": k,
        "Q": str(Q),
        "hypothesis_clause": HYPOTHESIS_CLAUSE,
    }


def detail_digits(detail) -> int:
    _emergent, (rn, rd), _smooth = detail
    return max(len(str(rn)), len(str(rd)))


def reason_from_verdict(verdict: str) -> str:
    if verdict.startswith("refused:"):
        suffix = verdict.removeprefix("refused:")
        if suffix.startswith("cofactor-big-digits-"):
            return "cofactor-big"
        return suffix
    if verdict in ("zero-jacobi", "bad-jacobi"):
        return "PrimalityBound"
    return verdict


def exact_xd(a: F, b: F, tau: F) -> tuple[F, F, dict]:
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
    Pb = sum(coefficient * b**degree for degree, coefficient in enumerate(P))
    assert M0 == Pb / (b**4 * D * D * A * A)
    return x0, d0, {
        "A": str(A),
        "delta": str(delta),
        "alpha": str(alpha),
        "D": str(D),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
    }


def replay_closure(row: dict, job: dict, verdict: str, info):
    """Require a second full ladder replay plus literal tied/ramified checks."""
    a = row["a"]
    eps = row["eps"]
    f = row["f"]
    q1 = row["q1"]
    cert = row["cert"]
    k = job["k"]
    Q = job["Q"]
    b = job["b"]
    rung = job["rung"]

    replayed_cert = fixed_factor_class_cert(a, Z, TAU, eps, f, q1)
    assert replayed_cert == cert and replayed_cert is not None and replayed_cert["ok"]
    assert tuple(replayed_cert["syms"].values()).count(1) == len(replayed_cert["syms"])
    assert Q == q1 + k * int(replayed_cert["N"])
    assert Q not in set(replayed_cert["excluded"])
    assert h10q._is_prime(Q) is True
    assert b == F(eps * f * Q)

    aligned_syms = l13_filter.aligned(a, Z, TAU, b)
    assert aligned_syms is not None and all(symbol == 1 for symbol in aligned_syms.values())
    smooth_status, smooth_detail = l13_filter.smooth_emergent(a, Z, TAU, b)
    if rung == "square":
        assert verdict == "zero" and smooth_status == "zero"
        replay_verdict, replay_info = "zero", None
    else:
        assert smooth_status == "cofactor-big"
        replay_verdict, replay_info = l13_filter.cofactor_decide(
            a, Z, TAU, b, detail=smooth_detail, deadline=None
        )
        assert replay_verdict == "zero" and replay_info == info
        assert replay_info["rung"] == rung
        assert all(
            place[1] == 1 and place[2] == "proved"
            for place in replay_info["places"]
        )

    tied = h10q._l7_tied_status(a, b, Z, TAU)
    if tied is not True:
        return None, {
            "status": "rejected:tied-not-literal-True",
            "tied": str(tied),
            "reason_code": "FactorBudget" if tied == "budget" else "tied-false",
        }

    x0, d0, identity = exact_xd(a, b, TAU)
    try:
        ramified = h10q.ramified(x0, d0)
    except (h10q.FactorBudget, h10q.PrimalityBound) as exc:
        return None, {
            "status": "refused:" + refusal_name(exc),
            "tied": "True",
            "reason_code": refusal_name(exc),
            **identity,
        }
    if ramified:
        return None, {
            "status": "rejected:ramified-nonempty",
            "tied": "True",
            "reason_code": "symbol-minus-one",
            "ramified": [str(place) for place in ramified],
            **identity,
        }

    closure = {
        "family": "ESC-NONCANONICAL",
        "cell": [W, list(UT)],
        "a": str(a),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": rung,
        "tied": "True",
        "ramified_empty": True,
        "source": SOURCE,
    }
    assert tuple(closure) == CLOSURE_FIELDS
    replay = {
        "type": "kernel-replay",
        "label": "PROVED",
        "status": "VERIFIED",
        "cell": [W, list(UT)],
        "family": "ESC-NONCANONICAL",
        "a": str(a),
        "tau": str(TAU),
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": str(cert["N"]),
        "k": k,
        "Q": str(Q),
        "b": str(b),
        "certificate_function": CERT_FUNCTION,
        "certificate_generalization": CERT_GENERALIZATION,
        "class_cert_replayed": True,
        "class_syms": cert["syms"],
        "aligned_syms": aligned_syms,
        "member_alignment_all_plus_one": True,
        "smooth_status": smooth_status,
        "smooth_detail": [
            smooth_detail[0],
            [str(smooth_detail[1][0]), str(smooth_detail[1][1])],
            smooth_detail[2],
        ],
        "emergent_verdict": replay_verdict,
        "cofactor_info": replay_info,
        "cofactor_replay_matches_initial": replay_info == info,
        "rung": rung,
        "tied": True,
        "ramified": [],
        "ramified_empty": True,
        "refusals_used_as_evidence": False,
        **identity,
    }
    return closure, replay


def scan_class(handle, row: dict, pacer: Pacer, accept_closure: bool):
    """Scan every k, deciding cofactor-big rows in ascending cofactor size."""
    started = time.monotonic()
    counts: Counter = Counter()
    refusal_counts: Counter = Counter()
    pending = []
    unexplored = []
    refused_positions: set[int] = set()
    closure_checks_not_needed = 0
    positions_decided = 0
    prime_members = 0
    member_decisions = 0
    largest_digits = 0
    closure = replay = None

    cert = row["cert"]
    N = int(cert["N"])
    excluded = set(cert["excluded"])
    for k in range(K_MAX + 1):
        pacer.tick()
        Q = row["q1"] + k * N
        record = member_base(row, k, Q)
        if Q in excluded:
            record.update(
                {
                    "verdict": "excluded",
                    "reason_code": "excluded",
                    "position_decided": True,
                    "mathematical_evidence": True,
                }
            )
            counts["excluded"] += 1
            positions_decided += 1
            emit(handle, record)
            continue

        try:
            is_prime = h10q._is_prime(Q)
        except (h10q.FactorBudget, h10q.PrimalityBound) as exc:
            reason = refusal_name(exc)
            record.update(
                {
                    "label": "OPEN",
                    "verdict": "refused:" + reason,
                    "reason_code": reason,
                    "position_decided": False,
                    "mathematical_evidence": False,
                }
            )
            counts["refused:" + reason] += 1
            refusal_counts[reason] += 1
            refused_positions.add(k)
            unexplored.append(
                {
                    "class_order_index": row["class_order_index"],
                    "a": str(row["a"]),
                    "eps": row["eps"],
                    "q1": row["q1"],
                    "N": str(N),
                    "k": k,
                    "Q": str(Q),
                    "phase": "Q-primality",
                    "reason": reason,
                }
            )
            emit(handle, record)
            continue
        if not is_prime:
            record.update(
                {
                    "Q_prime_proved": False,
                    "verdict": "composite-proved",
                    "reason_code": "composite-proved",
                    "position_decided": True,
                    "mathematical_evidence": True,
                }
            )
            counts["composite-proved"] += 1
            positions_decided += 1
            emit(handle, record)
            continue

        prime_members += 1
        b = F(row["eps"] * row["f"] * Q)
        record.update({"Q_prime_proved": True, "b": str(b)})
        aligned_syms = l13_filter.aligned(row["a"], Z, TAU, b)
        if aligned_syms is None:
            record.update(
                {
                    "label": "OPEN",
                    "verdict": "refused:aligned-none",
                    "reason_code": "aligned-none",
                    "position_decided": False,
                    "mathematical_evidence": False,
                }
            )
            counts["refused:aligned-none"] += 1
            refusal_counts["aligned-none"] += 1
            refused_positions.add(k)
            unexplored.append(
                {
                    "class_order_index": row["class_order_index"],
                    "a": str(row["a"]), "eps": row["eps"], "q1": row["q1"],
                    "N": str(N), "k": k, "Q": str(Q),
                    "phase": "member-alignment", "reason": "aligned-none",
                }
            )
            emit(handle, record)
            continue
        bad_places = [str(place) for place, symbol in aligned_syms.items() if symbol != 1]
        if bad_places:
            record.update(
                {
                    "aligned_syms": aligned_syms,
                    "verdict": "symbol-minus-one",
                    "reason_code": "symbol-minus-one",
                    "bad_places": bad_places,
                    "position_decided": True,
                    "member_decided": True,
                    "mathematical_evidence": True,
                }
            )
            counts["symbol-minus-one"] += 1
            positions_decided += 1
            member_decisions += 1
            emit(handle, record)
            continue

        smooth_status, smooth_detail = l13_filter.smooth_emergent(
            row["a"], Z, TAU, b
        )
        record.update(
            {
                "aligned_syms": aligned_syms,
                "member_alignment_all_plus_one": True,
                "smooth_status": smooth_status,
            }
        )
        if smooth_status == "bad":
            emergent, (rn, rd), smooth = smooth_detail
            record.update(
                {
                    "verdict": "bad",
                    "reason_code": "bad",
                    "small_bad_places": emergent,
                    "stripped_rn": str(rn),
                    "stripped_rd": str(rd),
                    "smooth": smooth,
                    "position_decided": True,
                    "member_decided": True,
                    "mathematical_evidence": True,
                }
            )
            counts["bad"] += 1
            positions_decided += 1
            member_decisions += 1
            emit(handle, record)
            continue
        if smooth_status not in ("zero", "cofactor-big"):
            reason = smooth_status
            record.update(
                {
                    "label": "OPEN",
                    "verdict": "refused:" + reason,
                    "reason_code": reason,
                    "position_decided": False,
                    "mathematical_evidence": False,
                }
            )
            counts["refused:" + reason] += 1
            refusal_counts[reason] += 1
            refused_positions.add(k)
            unexplored.append(
                {
                    "class_order_index": row["class_order_index"],
                    "a": str(row["a"]), "eps": row["eps"], "q1": row["q1"],
                    "N": str(N), "k": k, "Q": str(Q),
                    "phase": "smooth-emergent", "reason": reason,
                }
            )
            emit(handle, record)
            continue

        digits = 0 if smooth_status == "zero" else detail_digits(smooth_detail)
        largest_digits = max(largest_digits, digits)
        if smooth_status == "cofactor-big":
            counts["cofactor-big"] += 1
        pending.append(
            {
                "record": record,
                "k": k,
                "Q": Q,
                "b": b,
                "smooth_status": smooth_status,
                "smooth_detail": smooth_detail,
                "cofactor_digits": digits,
            }
        )

    pending.sort(key=lambda job: (job["cofactor_digits"], job["k"]))
    for job in pending:
        record = job["record"]
        record["cofactor_order_key"] = [job["cofactor_digits"], job["k"]]
        if job["smooth_status"] == "zero":
            verdict, info, rung = "zero", None, "square"
        else:
            verdict, info = l13_filter.cofactor_decide(
                row["a"], Z, TAU, job["b"], detail=job["smooth_detail"],
                deadline=None
            )
            rung = info.get("rung") if isinstance(info, dict) else None
        job["rung"] = rung
        record.update(
            {
                "verdict": verdict,
                "reason_code": reason_from_verdict(verdict),
                "rung": rung,
                "cofactor_digits": job["cofactor_digits"],
                "cofactor_info": info,
            }
        )
        counts[verdict] += 1

        if verdict == "INCONSISTENT":
            emit(handle, record)
            raise AssertionError("L13 parity inconsistency")
        if verdict.startswith("refused:") or verdict in ("zero-jacobi", "bad-jacobi"):
            reason = reason_from_verdict(verdict)
            record.update(
                {
                    "label": "OPEN",
                    "position_decided": False,
                    "member_decided": False,
                    "mathematical_evidence": False,
                }
            )
            refusal_counts[reason] += 1
            refused_positions.add(job["k"])
            unexplored.append(
                {
                    "class_order_index": row["class_order_index"],
                    "a": str(row["a"]), "eps": row["eps"], "q1": row["q1"],
                    "N": str(N), "k": job["k"], "Q": str(job["Q"]),
                    "phase": "cofactor", "reason": reason,
                    "raw_verdict": verdict,
                    "cofactor_digits": job["cofactor_digits"],
                }
            )
            emit(handle, record)
            continue

        record.update(
            {
                "position_decided": verdict != "zero",
                "member_decided": True,
                "mathematical_evidence": True,
            }
        )
        member_decisions += 1
        if verdict != "zero":
            positions_decided += 1
            emit(handle, record)
            continue

        if not accept_closure or closure is not None:
            record.update(
                {
                    "position_decided": True,
                    "accepted_as_closure": False,
                    "closure_replay_status": "not-needed:cell-already-closed",
                }
            )
            positions_decided += 1
            closure_checks_not_needed += 1
            emit(handle, record)
            continue
        candidate_closure, candidate_replay = replay_closure(row, job, verdict, info)
        if candidate_closure is None:
            reason = candidate_replay["reason_code"]
            record.update(
                {
                    "accepted_as_closure": False,
                    "closure_replay": candidate_replay,
                }
            )
            refusal_counts[reason] += 1
            refused_positions.add(job["k"])
            unexplored.append(
                {
                    "class_order_index": row["class_order_index"],
                    "a": str(row["a"]), "eps": row["eps"], "q1": row["q1"],
                    "N": str(N), "k": job["k"], "Q": str(job["Q"]),
                    "phase": "closure-crosscheck", "reason": reason,
                    "cofactor_digits": job["cofactor_digits"],
                }
            )
            emit(handle, record)
            continue

        record.update(
            {
                "position_decided": True,
                "accepted_as_closure": True,
                "full_replay_status": "VERIFIED",
            }
        )
        positions_decided += 1
        emit(handle, record)
        closure, replay = candidate_closure, candidate_replay
        continue

    total_positions = K_MAX + 1
    summary = {
        "type": "class-member-summary",
        "label": "PROVED" if closure is not None or not refusal_counts else "OPEN",
        "cell": [W, list(UT)],
        "class_order_index": row["class_order_index"],
        "a": str(row["a"]),
        "tau": str(TAU),
        "eps": row["eps"],
        "f": row["f"],
        "q1": row["q1"],
        "N": str(N),
        "k_protocol": [0, K_MAX],
        "positions_total": total_positions,
        "positions_decided": positions_decided,
        "positions_refused": len(refused_positions),
        "refusal_events": sum(refusal_counts.values()),
        "closure_checks_not_needed": closure_checks_not_needed,
        "decided_fraction": positions_decided / total_positions,
        "prime_members_proved": prime_members,
        "prime_members_decided": member_decisions,
        "reason_counts": dict(sorted(counts.items())),
        "refusal_reason_counts": dict(sorted(refusal_counts.items())),
        "cofactor_big_count": counts["cofactor-big"],
        "largest_cofactor_digits_attempted": largest_digits,
        "cofactor_order": "ascending stripped-cofactor digit count, then k",
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "closure": closure,
        "unexplored": unexplored,
        "refusals_used_as_evidence": False,
        "position_partition_verified": (
            positions_decided + len(refused_positions) == total_positions
        ),
    }
    emit(handle, summary)
    return summary, closure, replay


def write_report(result: dict) -> None:
    closure = result["closure"]
    verdict = "VERIFIED" if closure is not None else "OPEN"
    lines = [
        "# L19 cell 179",
        "",
        f"- **Verdict:** **{verdict}** for `[179,[-1,1]]`.",
        f"- Wall-clock time: **{result['wall_seconds']:.6f} s** (single process).",
        "- CPU protocol: `nice -n 19`; hot loops sleep 0.01 s every 200 candidates.",
        f"- Largest stripped-cofactor digit count attempted: **{result['largest_cofactor_digits']}**.",
        f"- Proven base primes: **{result['q1_count']}**, all primes `<= 400`.",
        f"- Class attempts: **{result['class_attempts']}**; aligned classes recorded: **{result['aligned_classes']}**.",
        f"- Member positions reached: **{result['positions_reached']}**; decided: **{result['positions_decided']}**; decided fraction: **{result['decided_fraction']:.6f}**.",
        f"- Refused positions: **{result['positions_refused']}** ({result['refusal_events']} refusal events); zero-member closure checks not needed after the verified hit: **{result['closure_checks_not_needed']}**. Refusals are never evidence.",
        "",
        "## Exact hypothesis clause",
        "",
        HYPOTHESIS_CLAUSE,
        "",
        "## Ordering and authority",
        "",
        "Aligned classes are ordered by increasing `N` (then `a`, `eps`, `q1`). Within each class, `k=0..60` is generated exactly; rows reaching the cofactor ladder are attempted by increasing stripped-cofactor digit count, then `k`. Class rows explicitly label the fixed-factor certificate as a local generalization. Member primality is accepted only from `h10q._is_prime`; member alignment and decisions use `l13_filter.aligned`, `smooth_emergent`, and `cofactor_decide`.",
        "",
        "## Per-class results reached",
        "",
        "| order | a | eps | q1 | N | decided/61 | refused | checks skipped | cofactor-big | max digits |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for summary in result["class_summaries"]:
        lines.append(
            "| {class_order_index} | {a} | {eps} | {q1} | {N} | "
            "{positions_decided}/61 | {positions_refused} | "
            "{closure_checks_not_needed} | {cofactor_big_count} | "
            "{largest_cofactor_digits_attempted} |".format(**summary)
        )
    lines += [
        "",
        "Each class's exact verdict and refusal reason maps are stored in its `class-member-summary` JSONL row.",
        "",
        "## Class-pool reason counts",
        "",
        "```json",
        json.dumps(result["class_reason_counts"], indent=2, sort_keys=True),
        "```",
        "",
        "## Member verdict counts",
        "",
        "```json",
        json.dumps(result["member_verdict_counts"], indent=2, sort_keys=True),
        "```",
        "",
        "## Member refusal reason counts",
        "",
        "```json",
        json.dumps(result["member_refusal_reason_counts"], indent=2, sort_keys=True),
        "```",
        "",
    ]
    if closure is not None:
        lines += [
            "## Verified closure row",
            "",
            "```json",
            json.dumps(closure, separators=(",", ":")),
            "```",
            "",
            "The row was accepted only after a second class/member/ladder replay, literal `_l7_tied_status is True`, and a successful `ramified(...) == []` computation.",
            "",
            "Replay from the h10q repository root:",
            "",
            "```sh",
            "nice -n 19 python3 l19_cell179.py",
            "```",
        ]
    else:
        lines += [
            "## Resumable cursor",
            "",
            f"The exact unresolved remainder is the `resume-cursor` row in `{OUT.name}` ({len(result['resume_cursor'])} entries):",
            "",
            "```json",
            json.dumps(result["resume_cursor"], indent=2),
            "```",
        ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def verify_artifacts(result: dict) -> None:
    records = [json.loads(line) for line in OUT.read_text(encoding="utf-8").splitlines()]
    aligned_rows = [record for record in records if record.get("type") == "aligned-class"]
    assert len(aligned_rows) == result["aligned_classes"]
    assert all(record["all_symbols_plus_one"] for record in aligned_rows)
    assert [int(record["N"]) for record in aligned_rows] == sorted(
        int(record["N"]) for record in aligned_rows
    )
    summaries = [record for record in records if record.get("type") == "class-member-summary"]
    assert all(record["position_partition_verified"] for record in summaries)
    closure_rows = [record for record in records if tuple(record) == CLOSURE_FIELDS]
    if result["closure"] is None:
        assert closure_rows == []
        assert any(record.get("type") == "resume-cursor" for record in records)
        assert "**OPEN**" in REPORT.read_text(encoding="utf-8")
    else:
        assert closure_rows == [result["closure"]]
        replay_rows = [record for record in records if record.get("type") == "kernel-replay"]
        assert len(replay_rows) == 1 and replay_rows[0]["status"] == "VERIFIED"
        assert replay_rows[0]["tied"] is True and replay_rows[0]["ramified"] == []
        assert "**VERIFIED**" in REPORT.read_text(encoding="utf-8")


def main() -> None:
    started = time.monotonic()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    q1_candidates = tuple(h10q.primerange(2, Q1_MAX + 1))
    assert q1_candidates and q1_candidates[-1] <= Q1_MAX
    assert all(h10q._is_prime(q1) is True for q1 in q1_candidates)
    pacer = Pacer()

    with OUT.open("w", encoding="utf-8") as handle:
        emit(
            handle,
            {
                "type": "meta",
                "label": "PROVED-protocol",
                "cell": [W, list(UT)],
                "z": str(Z),
                "a_values": list(A_VALUES),
                "tau": str(TAU),
                "eps_values": list(EPS_VALUES),
                "f": FIXED_FACTOR,
                "q1_bound": [2, Q1_MAX],
                "k_bound": [0, K_MAX],
                "hypothesis_clause": HYPOTHESIS_CLAUSE,
                "certificate_function": CERT_FUNCTION,
                "certificate_generalization": CERT_GENERALIZATION,
                "class_order": "increasing N, then a, eps(+ before -), q1",
                "cofactor_order": "within class: increasing stripped-cofactor digits, then k",
                "primality_authority": "h10q._is_prime (deterministic MR/Pocklington only)",
                "member_authority": [
                    "l13_filter.aligned",
                    "l13_filter.smooth_emergent",
                    "l13_filter.cofactor_decide",
                    "h10q._l7_tied_status",
                    "h10q.ramified",
                ],
                "refusals_are_evidence": False,
            },
        )
        aligned, branch_rows, class_reason_counts, class_refusals = enumerate_classes(
            handle, q1_candidates, pacer
        )

        class_summaries = []
        closure = replay = None
        resume_cursor = list(class_refusals)
        for row in aligned:
            summary, candidate_closure, candidate_replay = scan_class(
                handle, row, pacer, closure is None
            )
            class_summaries.append(summary)
            resume_cursor.extend(summary["unexplored"])
            if closure is None and candidate_closure is not None:
                closure, replay = candidate_closure, candidate_replay

        if closure is not None:
            emit(handle, replay)
            emit(handle, closure)
        else:
            for row in aligned[len(class_summaries):]:
                resume_cursor.append(
                    {
                        "class_order_index": row["class_order_index"],
                        "a": str(row["a"]),
                        "eps": row["eps"],
                        "q1": row["q1"],
                        "N": str(row["cert"]["N"]),
                        "k_next": 0,
                        "k_max": K_MAX,
                        "phase": "class-not-started",
                        "reason": "prior-work-exhausted",
                    }
                )
            emit(
                handle,
                {
                    "type": "resume-cursor",
                    "label": "OPEN",
                    "cell": [W, list(UT)],
                    "protocol": "retry exactly these class/member positions; no refusal is a tested negative",
                    "entries": resume_cursor,
                },
            )

        positions_reached = len(class_summaries) * (K_MAX + 1)
        positions_decided = sum(s["positions_decided"] for s in class_summaries)
        positions_refused = sum(s["positions_refused"] for s in class_summaries)
        refusal_events = sum(s["refusal_events"] for s in class_summaries)
        closure_checks_not_needed = sum(
            s["closure_checks_not_needed"] for s in class_summaries
        )
        assert positions_decided + positions_refused == positions_reached
        largest_digits = max(
            (s["largest_cofactor_digits_attempted"] for s in class_summaries),
            default=0,
        )
        member_verdict_counts: Counter = Counter()
        member_refusal_reason_counts: Counter = Counter()
        for summary in class_summaries:
            member_verdict_counts.update(summary["reason_counts"])
            member_refusal_reason_counts.update(summary["refusal_reason_counts"])
        wall_seconds = time.monotonic() - started
        result = {
            "closure": closure,
            "wall_seconds": wall_seconds,
            "largest_cofactor_digits": largest_digits,
            "q1_count": len(q1_candidates),
            "class_attempts": len(A_VALUES) * len(EPS_VALUES) * len(q1_candidates),
            "aligned_classes": len(aligned),
            "positions_reached": positions_reached,
            "positions_decided": positions_decided,
            "positions_refused": positions_refused,
            "refusal_events": refusal_events,
            "closure_checks_not_needed": closure_checks_not_needed,
            "decided_fraction": positions_decided / positions_reached if positions_reached else 0.0,
            "class_summaries": class_summaries,
            "class_reason_counts": dict(sorted(class_reason_counts.items())),
            "member_verdict_counts": dict(sorted(member_verdict_counts.items())),
            "member_refusal_reason_counts": dict(
                sorted(member_refusal_reason_counts.items())
            ),
            "resume_cursor": resume_cursor,
            "branch_rows": branch_rows,
        }
        emit(
            handle,
            {
                "type": "summary",
                "label": "PROVED" if closure is not None else "OPEN",
                "status": "VERIFIED-CLOSURE" if closure is not None else "OPEN",
                "cell": [W, list(UT)],
                "closure": closure,
                "wall_seconds": round(wall_seconds, 6),
                "largest_cofactor_digits_attempted": largest_digits,
                "class_attempts": result["class_attempts"],
                "aligned_classes": len(aligned),
                "aligned_classes_searched": len(class_summaries),
                "positions_reached": positions_reached,
                "positions_decided": positions_decided,
                "positions_refused": positions_refused,
                "refusal_events": refusal_events,
                "closure_checks_not_needed": closure_checks_not_needed,
                "decided_fraction": result["decided_fraction"],
                "class_reason_counts": result["class_reason_counts"],
                "member_verdict_counts": result["member_verdict_counts"],
                "member_refusal_reason_counts": result["member_refusal_reason_counts"],
                "resumable_cursor_entries": 0 if closure is not None else len(resume_cursor),
                "refusals_used_as_evidence": False,
                "source": SOURCE,
            },
        )

    write_report(result)
    verify_artifacts(result)
    print(json.dumps({
        "status": "VERIFIED-CLOSURE" if closure is not None else "OPEN",
        "closure": closure,
        "aligned_classes": len(aligned),
        "aligned_classes_searched": len(class_summaries),
        "positions_decided": positions_decided,
        "positions_reached": positions_reached,
        "decided_fraction": result["decided_fraction"],
        "largest_cofactor_digits_attempted": largest_digits,
        "wall_seconds": round(wall_seconds, 6),
        "data": str(OUT),
        "report": str(REPORT),
    }, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    main()
