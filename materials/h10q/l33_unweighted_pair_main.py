#!/usr/bin/env python3
"""L33 unconditional unweighted two-large bad-root pair main term.

The theorem is analytic: L23's Chebotarev law plus elementary CRT counting
shows that the unsigned two-large root-pair sector already has a positive
order-X main term on any prime window X^alpha..X^beta with
0<alpha<beta<1/2.  This producer replays the exact finite identities and keeps
the von Mangoldt/prime-weighted step explicitly outside its license.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l33_unweighted_pair_main.jsonl"
REPORT = Path("/tmp/l33_unweighted_pair_main.md")


def crt(left: int, p: int, right: int, q: int) -> int:
    return (left + p * (((right - left) * pow(p, -1, q)) % q)) % (p * q)


def interval_class_count(lower: int, upper: int, residue: int, modulus: int) -> int:
    return (upper - residue) // modulus - (lower - residue) // modulus


def exact_crt_replay() -> dict:
    bad_sets = {5: (1, 3), 7: (2, 4, 6), 11: (1, 5, 8)}
    rows = 0
    max_error = 0.0
    samples = []
    for p_index, p in enumerate(bad_sets):
        for q in tuple(bad_sets)[p_index + 1 :]:
            classes = {
                crt(left, p, right, q)
                for left in bad_sets[p]
                for right in bad_sets[q]
            }
            assert len(classes) == len(bad_sets[p]) * len(bad_sets[q])
            for residue in range(p * q):
                assert (residue in classes) == (
                    residue % p in bad_sets[p] and residue % q in bad_sets[q]
                )
            for X in (101, 1000, 10007, 100000):
                actual = sum(interval_class_count(X, 2 * X, c, p * q) for c in classes)
                expected = X * len(classes) / (p * q)
                error = abs(actual - expected)
                assert error <= len(classes)
                max_error = max(max_error, error)
                rows += 1
                if len(samples) < 6:
                    samples.append(
                        {
                            "p": p,
                            "q": q,
                            "X": X,
                            "classes": len(classes),
                            "actual": actual,
                            "expected": expected,
                            "error": error,
                        }
                    )
    return {
        "type": "unweighted-pair-CRT-replay",
        "label": "PROVED exact finite identity",
        "statement": "count in each CRT class is X/(p*q)+O(1), with total O(r_-(p)r_-(q))",
        "rows": rows,
        "maximum_observed_error": max_error,
        "samples": samples,
    }


def asymptotic_theorem() -> dict:
    alpha = 0.49
    beta = 0.499
    mu = 0.5 * math.log(beta / alpha)
    constant = 0.5 * mu * mu
    assert 0 < alpha < beta < 0.5
    assert constant > 0
    return {
        "type": "unweighted-two-large-root-pair-main",
        "label": "PROVED from L23 Chebotarev plus elementary CRT counting",
        "definition": (
            "U_2(X;alpha,beta)=sum_{X^alpha<=p<q<=X^beta} "
            "sum_{c in C_pq^-} #{X<t<=2X:t=c mod p*q}"
        ),
        "theorem": (
            "U_2(X;alpha,beta)~(1/8)*log(beta/alpha)^2*X "
            "for every fixed 0<alpha<beta<1/2"
        ),
        "proof": [
            "L23 Chebotarev and partial summation give sum r_-(p)/p=(1/2)*log(beta/alpha)+o(1)",
            "CRT gives r_-(p)r_-(q) classes and X/(p*q)+O(1) points per class",
            "the diagonal sum sum r_-(p)^2/p^2 is o(1)",
            "the total O(1)-per-class error is O(X^(2*beta)/log(X)^2)=o(X) because beta<1/2",
        ],
        "sample_window": {"alpha": alpha, "beta": beta},
        "sample_mu": mu,
        "sample_positive_constant": constant,
        "consequence": (
            "the unsigned root-pair sector is rigorously order X even inside p*q<X; "
            "no cancellation-only dispersion estimate can make it negligible"
        ),
        "strict_scope": (
            "unweighted t, mod-p root oversieve; this does not evaluate Lambda(Q(t)), "
            "impose exact odd valuations, or condition on L23 small-prime cleanliness"
        ),
    }


def prime_weight_boundary() -> dict:
    return {
        "type": "prime-weight-boundary",
        "label": "OPEN",
        "missing_transfer": (
            "replace the elementary class count by sum Lambda(Q(t)) uniformly on the "
            "factorable moduli p*q in [X^(2*alpha),X^(2*beta)]"
        ),
        "why_BV_stops": "2*alpha=0.98>1/2, so every modulus lies beyond Bombieri-Vinogradov",
        "exact_odd_boundary": "v_p(G)=v_q(G)=1 uses p^2*q^2 and is strictly harder",
        "needed_result": "a pointwise fixed-family Buchstab/analytic-Hilbert-detector asymptotic",
    }


def build_report(records: list[dict], elapsed: float) -> str:
    theorem = next(row for row in records if row["type"] == "unweighted-two-large-root-pair-main")
    return "\n".join(
        [
            "# L33 unweighted pair main",
            "",
            theorem["theorem"],
            "",
            "This is an unconditional positive main term in the mod-p root oversieve.",
            "It proves that signed cancellation cannot dispose of the unsigned pair sector.",
            "The Lambda(Q(t))-weighted transfer and exact odd valuations remain open.",
            "",
            f"Wall time: {elapsed:.3f} s.",
        ]
    ) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [exact_crt_replay(), asymptotic_theorem(), prime_weight_boundary()]
    elapsed = time.perf_counter() - started
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in records))
    REPORT.write_text(build_report(records, elapsed))
    print(
        "l33_unweighted_pair_main: positive order-X unweighted pair main proved; "
        "Lambda/exact-odd transfer open; "
        f"rows={len(records)} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
