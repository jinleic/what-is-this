#!/usr/bin/env python3
"""Exact size-biased block reduction for the remaining Liu H2 obstruction."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import solve_equation_parameters
from liu9_boundary_layer import gap_mp

OUTPUT = HERE / "verification/results/liu9-block-kernel.json"
QONE_SHA256 = "eb645792526578d115da8356f7d9378812164273900a60607691a8f32661f19a"


def _digest(payload) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _h(value):
    if value == 0 or value == 1:
        return mpmath.mpf(0)
    return -value * mpmath.log(value) - (1 - value) * mpmath.log(1 - value)


def _mu(value):
    if value == 0:
        return mpmath.mpf(1)
    if value == 1:
        return mpmath.mpf(0)
    return -(1 - value) * mpmath.log(1 - value) / value


def _pi(left, right):
    return left * right * (1 + (1 - left) * (1 - right))


def _c_kernel(left, right, mean, beta):
    return (
        (2 * mean - 1) / 2 * (
            _h(left) / left + _h(right) / right)
        - mean * (_mu(left) + _mu(right) - _mu(left * right))
        + beta * mean * (
            _h(_pi(left, right)) - _h(left * right)) / (left * right)
    )


def build_report() -> dict:
    parameters = solve_equation_parameters(100)
    with mpmath.workdps(100):
        x = mpmath.mpf("0.063")
        q = mpmath.mpf("0.95")
        qbar = 1 - q
        mean = qbar * x + q
        q_weights = (qbar * x / mean, q / mean)
        supports = (x, mpmath.mpf(1))
        frozen = parameters.mean * sum(
            q_weights[i] * q_weights[j] * _c_kernel(
                supports[i], supports[j], parameters.mean, parameters.beta)
            for i in range(2) for j in range(2)
        )
        signed_protocol = _h(_pi(x, x)) - 2 * _h(x)
        surrogate = frozen + parameters.beta * q * qbar * signed_protocol
        values = (
            mpmath.mpf(0), mpmath.mpf(0), q,
            mpmath.mpf("0.5"), mpmath.mpf("0.5"), x,
            mpmath.mpf("0.5"), mpmath.mpf("0.5"), mpmath.mpf(1),
        )
        raw_gap = gap_mp(values, parameters.beta)
    if not surrogate < 0 or not raw_gap > 0 or not mean > parameters.mean:
        raise AssertionError("shortcut discriminator lost its signs")

    report = {
        "tool": "liu9_block_kernel.py",
        "claim_status": "OPEN_REDUCTION",
        "qone_kernel_report_sha256": QONE_SHA256,
        "block_kernel": {
            "B00": "((1-q)^2/M) C_M + beta*q*(1-q)*k_pi",
            "B11": "(q^2/M) C_M + beta*q*(1-q)*k_pi",
            "B01": "(q*(1-q)/M) C_M - beta*q*(1-q)*k_pi",
            "measures": "alpha_i(dx)=x P_i(dx)",
        },
        "equivalence": (
            "The full raw gap equals <alpha0^2,B00>+<alpha1^2,B11>+"
            "2<alpha0 alpha1,B01>; H2 is exactly block copositivity on the "
            "positive paired measures with shared masses."
        ),
        "false_shortcut": {
            "P0": "delta_0.063",
            "P1": "delta_1",
            "q": "0.95",
            "mean": mpmath.nstr(mean, 80),
            "m_frozen_surrogate": mpmath.nstr(surrogate, 80),
            "raw_gap_gap_mp": mpmath.nstr(raw_gap, 80),
            "interpretation": (
                "Refutes freezing C_M to C_m inside the full correction; "
                "does not refute H2 because the true raw gap is positive."
            ),
        },
        "refuted_routes": [
            "k_pi is PSD",
            "cross product kernel dominates both component diagonals",
            "rank-one square completion",
            "m-frozen global correction is nonnegative",
        ],
        "remaining_obligation": (
            "Prove the displayed 2x2 block kernel copositive for every "
            "q in [0,1], mixture mean M>=m, and positive paired measures, "
            "or exhibit a negative raw-gap configuration."
        ),
        "negative_mean_feasible_raw_gap_found": False,
    }
    report["report_sha256"] = _digest(report)
    return report


def main() -> int:
    report = build_report()
    OUTPUT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 EXACT BLOCK-KERNEL REDUCTION")
    print("OPEN block copositivity; false shortcut %s; raw gap %s"
          % (report["false_shortcut"]["m_frozen_surrogate"],
             report["false_shortcut"]["raw_gap_gap_mp"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
