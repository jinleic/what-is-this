"""Audit the proposed finite Kesten SAW-ratio bound and compute its consequence.

The finite one-sided inequality is deliberately kept separate from Kesten's
proved ratio limit.  Exact integer brackets and directed interval arithmetic
record what the proposed inequality *would* imply, while the certified output
remains the incumbent bound unless an auditable one-sided theorem is present.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path

import mpmath as mp


SCRIPT = "experiments/e55_kesten_bound.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "kesten_bound.json"
DPS = 90
BRACKET_DIGITS = 90
REPORT_DECIMAL_PLACES = 40
TASK_INCUMBENT_DISPLAY = Decimal("0.21221190116616784")
INCUMBENT_OUTWARD_LOWER = Decimal(
    "0.2122119011661678393310862783954278184914"
)

COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: 2941370856334701726560670,
}
CANDIDATE_N = (31, 32, 33, 34)

SOURCE_SPECS = {
    "schram2011_exact_saw": (
        "sources/fulltext/schram_barkema_bisseling2011.pdf",
        "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12",
    ),
    "bauerschmidt2012_saw_lectures": (
        "sources/fulltext/bauerschmidt_duminil_copin_goodman_slade2012_saw_lectures.pdf",
        "f41468f516eea8d4ba5866d8c444636fa579ca57bcb27210204718515a2530cd",
    ),
    "madras_slade1996_ch7_preview": (
        "sources/fulltext/madras_slade1996_ch7_preview.pdf",
        "4f3db221eb1a4fb8162b4f1fed325c3f41ad9d7912550deb51ee8a0200bc8135",
    ),
    "grimmett2020_kesten_work": (
        "sources/fulltext/grimmett2020_kesten_work.pdf",
        "1ab2dc69831fcdf3ca26a05506cb6766dffd6d24cbcb8fc427e9555b0f4729b7",
    ),
    "schram2017_bcc_fcc": (
        "sources/fulltext/schram_barkema_bisseling_clisby2017_bcc_fcc.pdf",
        "002b53281c491a0589b43139918db1cd22d0c35edfa4f9f9060e744175ebfc2f",
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _interval_endpoints(value: object) -> tuple[str, str]:
    text = str(value).strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError(f"unexpected interval representation: {text}")
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return lower.strip(), upper.strip()


def _floor_decimal(value: str, places: int = REPORT_DECIMAL_PLACES) -> str:
    with localcontext() as context:
        context.prec = 200
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _ceil_decimal(value: str, places: int = REPORT_DECIMAL_PLACES) -> str:
    with localcontext() as context:
        context.prec = 200
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_CEILING), "f")


def _sqrt_ratio_bracket(
    numerator: int, denominator: int, *, digits: int = BRACKET_DIGITS
) -> tuple[int, int, int]:
    """Bracket sqrt(numerator/denominator) by adjacent scaled integers."""

    if numerator <= 0 or denominator <= 0:
        raise ValueError("ratio must be positive")
    scale = 10**digits
    target = numerator * scale * scale
    lo = 0
    hi = scale
    while denominator * hi * hi <= target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if denominator * mid * mid <= target:
            lo = mid
        else:
            hi = mid
    if not (
        denominator * lo * lo <= target < denominator * hi * hi
        and hi == lo + 1
    ):
        raise AssertionError("failed exact square-root bracket")
    return lo, hi, scale


def _inverse_count_root_bracket(
    count: int, steps: int, *, digits: int = BRACKET_DIGITS
) -> tuple[int, int, int]:
    """Bracket count**(-1/steps) by adjacent scaled integers."""

    scale = 10**digits
    target = scale**steps
    lo, hi = 0, scale
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if count * mid**steps <= target:
            lo = mid
        else:
            hi = mid
    if not (count * lo**steps <= target < count * hi**steps and hi == lo + 1):
        raise AssertionError("failed exact inverse-root bracket")
    return lo, hi, scale


def _rational_decimal_interval(lo: int, hi: int, scale: int) -> tuple[str, str]:
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        lower = _interval_endpoints(mp.iv.mpf(lo) / scale)[0]
        upper = _interval_endpoints(mp.iv.mpf(hi) / scale)[1]
        return lower, upper
    finally:
        mp.iv.dps = previous


def _atanh_interval(lo: int, hi: int, scale: int) -> tuple[str, str]:
    """Directed enclosure of atanh(v) for lo/scale <= v < hi/scale."""

    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        v_lo = mp.iv.mpf(lo) / scale
        v_hi = mp.iv.mpf(hi) / scale
        k_lo = mp.iv.log((1 + v_lo) / (1 - v_lo)) / 2
        k_hi = mp.iv.log((1 + v_hi) / (1 - v_hi)) / 2
        return _interval_endpoints(k_lo)[0], _interval_endpoints(k_hi)[1]
    finally:
        mp.iv.dps = previous


def _bracket_payload(lo: int, hi: int, scale: int) -> dict[str, object]:
    decimal_lo, decimal_hi = _rational_decimal_interval(lo, hi, scale)
    return {
        "lower_numerator": str(lo),
        "upper_numerator": str(hi),
        "denominator": str(scale),
        "exact_relation": (
            "[COMPUTATION] lower_numerator/denominator <= target "
            "< upper_numerator/denominator; numerators are adjacent integers"
        ),
        "directed_decimal_interval": [decimal_lo, decimal_hi],
    }


def _record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def main() -> None:
    checks: list[dict[str, object]] = []

    sources: dict[str, object] = {}
    for key, (relative_path, expected_hash) in SOURCE_SPECS.items():
        path = ROOT / relative_path
        actual_hash = _sha256(path) if path.is_file() else None
        verified = actual_hash == expected_hash
        sources[key] = {
            "local_path": relative_path,
            "expected_sha256": expected_hash,
            "actual_sha256": actual_hash,
            "verified": verified,
        }
        _record(
            checks,
            f"source_sha256_{key}",
            verified,
            f"[COMPUTATION] {relative_path} SHA-256 equals the audited manifest value",
        )

    candidates: list[dict[str, object]] = []
    brackets_valid = True
    for n in CANDIDATE_N:
        c_n = COUNTS[n]
        c_next = COUNTS[n + 2]
        mu_lo, mu_hi, mu_scale = _sqrt_ratio_bracket(c_next, c_n)
        v_lo, v_hi, v_scale = _sqrt_ratio_bracket(c_n, c_next)
        k_lo, k_hi = _atanh_interval(v_lo, v_hi, v_scale)
        brackets_valid &= (
            c_n * mu_lo * mu_lo <= c_next * mu_scale * mu_scale
            < c_n * mu_hi * mu_hi
        )
        brackets_valid &= (
            c_next * v_lo * v_lo <= c_n * v_scale * v_scale
            < c_next * v_hi * v_hi
        )
        candidates.append(
            {
                "n": n,
                "parity": "odd" if n % 2 else "even",
                "counts": {str(n): str(c_n), str(n + 2): str(c_next)},
                "status": (
                    "[UNRESOLVED] This is a conditional arithmetic consequence only; "
                    "the required finite one-sided ratio theorem was not audited"
                ),
                "unproved_assumption": (
                    f"[UNRESOLVED] mu^2 <= c_{n + 2}/c_{n} for this finite n"
                ),
                "conditional_mu_upper_exact": f"sqrt({c_next}/{c_n})",
                "conditional_mu_upper_bracket": _bracket_payload(
                    mu_lo, mu_hi, mu_scale
                ),
                "conditional_v_lower_exact": f"sqrt({c_n}/{c_next})",
                "conditional_v_bracket": _bracket_payload(v_lo, v_hi, v_scale),
                "conditional_k_exact": f"atanh(sqrt({c_n}/{c_next}))",
                "conditional_k_directed_interval": [k_lo, k_hi],
                "conditional_k_outward_rounded_40": [
                    _floor_decimal(k_lo),
                    _ceil_decimal(k_hi),
                ],
                "comparison_with_task_incumbent": (
                    "[COMPUTATION] conditional lower endpoint is strictly greater "
                    "than 0.21221190116616784"
                    if Decimal(k_lo) > TASK_INCUMBENT_DISPLAY
                    else "[COMPUTATION] conditional lower endpoint is not greater than 0.21221190116616784"
                ),
            }
        )

    best_n = max(CANDIDATE_N, key=lambda n: Fraction(COUNTS[n], COUNTS[n + 2]))
    best = next(candidate for candidate in candidates if candidate["n"] == best_n)

    even_ratio_decreases = COUNTS[36] * COUNTS[32] < COUNTS[34] * COUNTS[34]
    odd_ratio_decreases = COUNTS[35] * COUNTS[31] < COUNTS[33] * COUNTS[33]
    _record(
        checks,
        "all_exact_sqrt_brackets",
        brackets_valid,
        "[COMPUTATION] every square-root bracket passes its cleared-integer inequalities",
    )
    _record(
        checks,
        "even_ratio_stability",
        even_ratio_decreases,
        "[COMPUTATION] c_36/c_34 < c_34/c_32 by exact cross multiplication",
    )
    _record(
        checks,
        "odd_ratio_stability",
        odd_ratio_decreases,
        "[COMPUTATION] c_35/c_33 < c_33/c_31 by exact cross multiplication",
    )
    _record(
        checks,
        "best_conditional_candidate",
        best_n == 34,
        "[COMPUTATION] n=34 maximizes exact c_n/c_(n+2) among n=31,32,33,34",
    )
    _record(
        checks,
        "conditional_candidate_beats_task_display",
        Decimal(best["conditional_k_directed_interval"][0])
        > TASK_INCUMBENT_DISPLAY,
        "[COMPUTATION] the n=34 conditional lower enclosure exceeds the task's incumbent display",
    )

    incumbent_v_lo, incumbent_v_hi, incumbent_scale = _inverse_count_root_bracket(
        COUNTS[36], 36
    )
    incumbent_k_lo, incumbent_k_hi = _atanh_interval(
        incumbent_v_lo, incumbent_v_hi, incumbent_scale
    )
    incumbent_outward = [
        _floor_decimal(incumbent_k_lo),
        _ceil_decimal(incumbent_k_hi),
    ]
    incumbent_bracket_valid = (
        COUNTS[36] * incumbent_v_lo**36 <= incumbent_scale**36
        < COUNTS[36] * incumbent_v_hi**36
    )
    _record(
        checks,
        "incumbent_inverse_root_bracket",
        incumbent_bracket_valid,
        "[COMPUTATION] c_36^(-1/36) bracket passes exact cleared-integer inequalities",
    )
    _record(
        checks,
        "incumbent_outward_lower_reproduced",
        Decimal(incumbent_outward[0]) == INCUMBENT_OUTWARD_LOWER,
        "[COMPUTATION] independent 90-dps directed rounding reproduces the stored 40-place lower endpoint",
    )

    theorem_audit = {
        "proposed_finite_statement": (
            "[UNRESOLVED] For every n on Z^d, c_(n+2) >= mu^2 c_n"
        ),
        "verdict": (
            "[UNRESOLVED] rejected as a repository certificate: no obtained source "
            "states or proves the finite one-sided inequality"
        ),
        "primary_source": (
            "[EXTERNAL] Kesten (1963), DOI 10.1063/1.1704022: the Crossref/AIP "
            "abstract states c_(n+2)/c_n - mu^2 tends to zero; full text was not obtainable"
        ),
        "madras_slade": (
            "[EXTERNAL] Madras--Slade, Chapter 7, Theorem 7.3.4(a), p.248, states "
            "lim_(n->infinity) c_(n+2)/c_n = mu^2. The alleged 'Theorem 1.2.5' "
            "is not the ratio theorem; (1.2.5) is an equation in the subadditivity section"
        ),
        "open_lecture_notes": (
            "[EXTERNAL] Bauerschmidt--Duminil-Copin--Goodman--Slade, Section 1.3, "
            "equation (1.17), states only the same ratio limit"
        ),
        "quantitative_form": (
            "[EXTERNAL] Grimmett (2020), Theorem 5.1(a), records "
            "abs(c_(n+2)/c_n - mu^2) <= A n^(-1/3) for an unspecified constant A; "
            "this supplies neither the proposed sign nor a usable numerical A"
        ),
        "logical_gap": (
            "[LEMMA] Convergence of a sequence to mu^2 does not determine which side "
            "a finite term lies on; therefore a single exact finite ratio cannot be "
            "promoted to an upper bound on mu from the audited ratio-limit theorem"
        ),
        "missing_for_upgrade": (
            "[UNRESOLVED] A page-level theorem and reproducible proof of the finite "
            "one-sided sign, with hypotheses covering nearest-neighbour SAWs on Z^3, "
            "or a complete independent proof"
        ),
    }
    _record(
        checks,
        "finite_theorem_not_promoted",
        theorem_audit["verdict"].startswith("[UNRESOLVED]"),
        "[COMPUTATION] result classification keeps the unaudited finite inequality out of the certified endpoint",
    )

    data = {
        "model": "nearest-neighbour ferromagnetic Ising model on the simple cubic lattice; K=beta*J and v=tanh(K)",
        "classification": (
            "[UNRESOLVED] No improved certified lower bound: the proposed finite "
            "Kesten inequality could not be audited to repository standard"
        ),
        "theorem_audit": theorem_audit,
        "source_certificates": sources,
        "exact_saw_counts": {
            "statement": (
                "[EXTERNAL COMPUTATION] Schram--Barkema--Bisseling Table I gives "
                "these exact rooted simple-cubic SAW counts"
            ),
            "counts": {str(n): str(COUNTS[n]) for n in sorted(COUNTS)},
        },
        "conditional_ratio_candidates": candidates,
        "conditional_best": {
            "statement": (
                "[COMPUTATION] Conditional on the unaudited finite inequality, n=34 "
                "is the strongest tested parity-resolved candidate"
            ),
            "n": best_n,
            "exact": best["conditional_k_exact"],
            "directed_interval": best["conditional_k_directed_interval"],
            "outward_rounded_40": best["conditional_k_outward_rounded_40"],
            "certified_for_ising": False,
        },
        "ratio_stability": {
            "even": (
                "[COMPUTATION] c_36/c_34 < c_34/c_32 exactly; this is a sanity "
                "observation and is not used as a monotonicity assumption"
            ),
            "odd": (
                "[COMPUTATION] c_35/c_33 < c_33/c_31 exactly; this is a sanity "
                "observation and is not used as a monotonicity assumption"
            ),
        },
        "longer_count_audit": {
            "statement": (
                "[EXTERNAL] Schram--Barkema--Bisseling--Clisby (2017), Introduction, "
                "calls N_max=36 the record simple-cubic series as of 2017"
            ),
            "outcome": (
                "[UNRESOLVED] No freely auditable post-2011 exact simple-cubic c_n "
                "with n>36 was located; this search result is not a proof of nonexistence"
            ),
        },
        "final_certified_lower_bound": {
            "status": (
                "[THEOREM from EXTERNAL EXACT COMPUTATION] incumbent unchanged; "
                "uses SAW domination, submultiplicativity, and exact c_36"
            ),
            "exact": f"atanh({COUNTS[36]}^(-1/36))",
            "v_exact_bracket": _bracket_payload(
                incumbent_v_lo, incumbent_v_hi, incumbent_scale
            ),
            "directed_interval": [incumbent_k_lo, incumbent_k_hi],
            "outward_rounded_40": incumbent_outward,
            "improves_incumbent": False,
            "reason": (
                "[UNRESOLVED] The only numerically stronger endpoint depends on a "
                "finite one-sided ratio inequality not present in the audited sources"
            ),
        },
        "comparison_only": {
            "task_incumbent_display": str(TASK_INCUMBENT_DISPLAY),
            "statement": (
                "[COMPUTATION] The conditional n=34 value exceeds the task display, "
                "but the certified final value remains the prior outward-rounded endpoint"
            ),
        },
    }

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact integer SAW data; adjacent rational brackets by integer bisection; "
                "mpmath.iv directed rounding at 90 dps; source SHA-256 audit"
            ),
            "mpmath_iv_dps": DPS,
            "rational_bracket_digits": BRACKET_DIGITS,
        },
        "data": data,
        "checks": checks,
    }
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise SystemExit(f"embedded checks failed: {failed}")
    print("PASS")


if __name__ == "__main__":
    main()
