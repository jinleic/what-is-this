#!/usr/bin/env python3
"""L34 sharp pair-moment criterion for AP1.

Hilbert reciprocity makes the remaining bad-place count even.  Therefore AP1
does not require the two-large sector to be negligible: it is enough that its
factorial second moment have relative constant strictly below one.  The
unrestricted Chebotarev local-mass envelope is below one precisely for
z=X^theta with theta>8*exp(-2*sqrt(2)); L23's theta=0.49 lies in that
window.  The product constraint on actual divisor pairs can only lower this
envelope.  The pointwise prime-weighted moment estimate itself remains open.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l34_ap1_pair_threshold.jsonl"
REPORT = Path("/tmp/l34_ap1_pair_threshold.md")


def parity_pair_criterion() -> dict:
    rows = []
    for bad_count in range(0, 17, 2):
        pairs = bad_count * (bad_count - 1) // 2
        zero_indicator = int(bad_count == 0)
        assert zero_indicator >= 1 - pairs
        rows.append(
            {
                "R_bad": bad_count,
                "pair_moment": pairs,
                "one_minus_pair_moment": 1 - pairs,
                "zero_indicator": zero_indicator,
            }
        )
    return {
        "type": "AP1-factorial-second-moment-criterion",
        "label": "PROVED exact combinatorial inequality",
        "pointwise_inequality": "1_{R_bad=0} >= 1-binomial(R_bad,2) for even R_bad",
        "summed_inequality": "N_good(X)>=#A_z(X)-P_2(X)>=#A_z(X)-P_2^root(X)",
        "definition": (
            "P_2=sum binomial(R_bad,2); P_2^root uses bad-sign root primes "
            "p|G and oversieves even valuations"
        ),
        "consequence": (
            "P_2^root(X)<#A_z(X) implies a zero-bad member and hence AP1; "
            "exact odd valuations are unnecessary for this upper-bound route"
        ),
        "rows": rows,
    }


def local_constant_threshold() -> dict:
    degree = 8.0
    theta = 0.49
    mu = 0.5 * math.log(degree / theta)
    pair_constant = 0.5 * mu * mu
    threshold = degree * math.exp(-2 * math.sqrt(2))
    endpoint_mu = 0.5 * math.log(degree / 0.5)
    endpoint_pair = 0.5 * endpoint_mu * endpoint_mu
    assert threshold < theta < 0.5
    assert pair_constant < 1
    assert endpoint_pair < pair_constant
    return {
        "type": "octic-pair-main-threshold",
        "label": "PROVED algebraic threshold for the unrestricted local-mass envelope",
        "degree": 8,
        "mu(theta)": "(1/2)*log(8/theta)",
        "unrestricted_relative_pair_envelope": "mu(theta)^2/2",
        "envelope_below_one_iff": "theta>8*exp(-2*sqrt(2))",
        "threshold": threshold,
        "L23_theta": theta,
        "L23_mu": mu,
        "L23_pair_constant": pair_constant,
        "L23_margin_below_one": 1 - pair_constant,
        "BV_endpoint": {
            "theta": 0.5,
            "mu": endpoint_mu,
            "pair_constant": endpoint_pair,
            "margin_below_one": 1 - endpoint_pair,
        },
        "interpretation": (
            "a conditioned prime-weighted pair-moment upper bound by this "
            "unrestricted envelope would imply AP1; the divisor-product constraint "
            "can only reduce the true pair region, and an o-bound is unnecessary"
        ),
    }


def analytic_boundary() -> dict:
    return {
        "type": "sharp-pair-moment-analytic-boundary",
        "label": "OPEN estimate, exact sufficient form proved",
        "sufficient_estimate": (
            "P_2^root(X)<=(mu(0.49)^2/2+o(1))*#A_z(X), with matching singular series/local conditioning"
        ),
        "envelope_scope": (
            "mu^2/2 integrates the unrestricted square of prime log-sizes; "
            "actual pairs also satisfy log_X(p)+log_X(q)<=8, so this is a safe "
            "target envelope rather than a claimed sharp actual constant"
        ),
        "why_it_closes": (
            "P_2<=P_2^root, mu(0.49)^2/2<1, and Hilbert parity makes every "
            "failure contribute at least one exact bad pair"
        ),
        "why_current_sieve_stops": [
            "the lower semilinear sieve is used at s just above its limit 1 and does not give a sharp relative constant",
            "every root-pair modulus p*q is at least X^0.98, beyond Bombieri-Vinogradov",
            "L33 proves the unweighted pair main is positive, so cancellation to o(X) is the wrong objective",
        ],
        "Schinzel_status": (
            "unconditional seven-variable untied result remains separate; the six-variable chain still needs this moment theorem/AP1"
        ),
    }


def build_report(records: list[dict], elapsed: float) -> str:
    constants = next(row for row in records if row["type"] == "octic-pair-main-threshold")
    return "\n".join(
        [
            "# L34 AP1 pair threshold",
            "",
            "For even R_bad, 1_{R_bad=0} >= 1-binomial(R_bad,2).",
            f"At theta=0.49 the unrestricted pair envelope is {constants['L23_pair_constant']:.12f}<1.",
            "Thus a conditioned root-oversieve pair upper bound would imply AP1; an o-bound is not needed.",
            "The prime-weighted, small-prime-clean root-pair moment estimate remains open.",
            "",
            f"Wall time: {elapsed:.3f} s.",
        ]
    ) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [parity_pair_criterion(), local_constant_threshold(), analytic_boundary()]
    elapsed = time.perf_counter() - started
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in records))
    REPORT.write_text(build_report(records, elapsed))
    print(
        "l34_ap1_pair_threshold: sharp relative pair criterion proved; "
        "prime-weighted moment/AP1 open; "
        f"rows={len(records)} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
