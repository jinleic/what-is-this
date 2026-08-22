#!/usr/bin/env python3
"""Standalone verifier for e66: the Simon--Lieb finite-family obstruction.

Independently recomputes a sample of the exact integer sign certificates from
the library (never importing the experiment module), re-derives the control
polynomial identity, re-does the directed rounding comparisons, and validates
the stored artifact envelope and conclusions.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from decimal import Decimal, ROUND_FLOOR, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.rigorous_bounds.simon_lieb import (  # noqa: E402
    atanh_rational_interval,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
)

RESULT = ROOT / "results" / "bounds" / "simon_lieb_next.json"
INCUMBENT = ROOT / "results" / "bounds" / "improved_bounds.json"
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(encoded).hexdigest(),
    }


def floor_decimal(value: str, places: int = 40) -> str:
    with localcontext() as context:
        context.prec = max(110, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == "experiments/e66_simon_lieb_next.py"
        and bool(artifact["checks"])
        and all(item["passed"] for item in artifact["checks"]),
    )
    data = artifact["data"]
    records = {tuple(record["shape"]): record for record in data["certificates"]}

    # 1. Independent recomputation of exact sign certificates (library only).
    for shape in ((2, 2, 3), (3, 3, 4)):
        record = records[shape]
        p, q = map(int, record["safe_exact_t"].split("/"))
        residual = exact_criterion_residual(shape, p, q)
        check(
            f"safe sign recomputed on {'x'.join(map(str, shape))}",
            residual < 0 and fingerprint(residual) == record["safe_residual_tQ_minus_P"],
            f"exact tQ-P < 0 at t={p}/{q}, fingerprint matches",
        )
    for shape in ((2, 2, 2), (3, 3, 4)):
        record = records[shape]
        cp, cq = map(int, record["challenge_exact_t"].split("/"))
        check("challenge point is 23/110 on " + "x".join(map(str, shape)), (cp, cq) == (23, 110))
        residual = exact_criterion_residual(shape, cp, cq)
        check(
            f"challenge rejection recomputed on {'x'.join(map(str, shape))}",
            residual >= 0
            and fingerprint(residual) == record["challenge_residual_tQ_minus_P"]
            and record["challenge_passes"] is False,
            "exact tQ-P >= 0 at t=23/110, fingerprint matches",
        )

    # 2. Control: complete spin enumeration versus parity transfer at t=1/5.
    p_poly, q_poly = enumerate_pq_polynomials((2, 2, 2))
    stored_controls = data["independent_exact_controls"]
    transfer_p, transfer_q = exact_pq_at_rational((2, 2, 2), 1, 5)
    degree_p = len(p_poly) - 1
    degree_q = len(q_poly) - 1
    scaled_p = sum(c * 5 ** (degree_p - i) for i, c in enumerate(p_poly))
    scaled_q = sum(c * 5 ** (degree_q - i) for i, c in enumerate(q_poly))
    check(
        "control enumeration vs transfer",
        transfer_p == scaled_p
        and transfer_q == scaled_q
        and list(p_poly) == stored_controls["P_coefficients"]
        and list(q_poly) == stored_controls["Q_coefficients"],
    )

    # 3. Directed rounding and the ordering that defines the negative result.
    stored_incumbent = data["final_lower_bound"]["decimal_lower"]
    challenge = data["challenge_beyond_incumbent"]
    interval = atanh_rational_interval(23, 110, dps=90)
    check(
        "challenge interval reproduced",
        floor_decimal(interval[0]) == challenge["K_decimal_lower"]
        and Decimal(challenge["K_decimal_lower"]) > Decimal(stored_incumbent),
        "atanh(23/110) directed lower end matches and strictly exceeds the incumbent",
    )
    best = data["best_completed_exact_certificate"]
    bp, bq = map(int, best["safe_exact_t"].split("/"))
    best_interval = atanh_rational_interval(bp, bq, dps=90)
    check(
        "best endpoint reproduced and below incumbent",
        floor_decimal(best_interval[0]) == best["safe_K_decimal_lower"]
        and Decimal(best["safe_K_decimal_lower"]) < Decimal(stored_incumbent),
        f"best completed exact endpoint {best['safe_K_decimal_lower']}",
    )
    check(
        "conclusion fields honest",
        data["final_lower_bound"]["improves_incumbent"] is False
        and "[UNRESOLVED]" in data["classification"]
        and "not_a_global_no_go" in data["resource_obstruction"],
    )
    check(
        "every stored shape rejects the challenge",
        all(record["challenge_passes"] is False for record in records.values())
        and len(records) == 20,
    )

    # 4. Search-order policy is benchmark-free.
    shapes = [tuple(shape) for shape in data["search_policy"]["shapes"]]
    check(
        "search order is size then lexicographic",
        shapes == sorted(shapes, key=lambda s: (math.prod(s), s)),
    )

    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {FAILURES}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
