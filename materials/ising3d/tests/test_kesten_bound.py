"""Independent standalone verification of the Kesten-ratio audit artifact."""

from __future__ import annotations

import hashlib
import json
import math
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "kesten_bound.json"
MANIFEST = ROOT / "sources" / "manifest.yaml"
DPS = 100
DIGITS = 90
SCALE = 10**DIGITS
EXPECTED_INCUMBENT_LOWER = Decimal(
    "0.2122119011661678393310862783954278184914"
)
TASK_INCUMBENT_DISPLAY = Decimal("0.21221190116616784")

EXPECTED_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: 2941370856334701726560670,
}
EXPECTED_SOURCES = {
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


def _interval_endpoints(value: object) -> tuple[Decimal, Decimal]:
    text = str(value).strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise AssertionError(f"unexpected interval representation: {text}")
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return Decimal(lower.strip()), Decimal(upper.strip())


def _sqrt_ratio_bracket_isqrt(
    numerator: int, denominator: int
) -> tuple[int, int]:
    """Independent adjacent bracket using math.isqrt, not experiment bisection."""

    target = numerator * SCALE * SCALE
    lo = math.isqrt(target // denominator)
    while denominator * (lo + 1) * (lo + 1) <= target:
        lo += 1
    while denominator * lo * lo > target:
        lo -= 1
    hi = lo + 1
    assert denominator * lo * lo <= target < denominator * hi * hi
    return lo, hi


def _inverse_count_root_bracket(count: int, steps: int) -> tuple[int, int]:
    """Independent exact bisection for the incumbent inverse 36th root."""

    target = SCALE**steps
    low, high = 0, SCALE
    while low + 1 < high:
        middle = (low + high) >> 1
        if count * pow(middle, steps) <= target:
            low = middle
        else:
            high = middle
    assert count * pow(low, steps) <= target < count * pow(high, steps)
    return low, high


def _directed_atanh(lo: int, hi: int) -> tuple[Decimal, Decimal]:
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        v_lo = mp.iv.mpf(lo) / SCALE
        v_hi = mp.iv.mpf(hi) / SCALE
        lower_value = mp.iv.log((1 + v_lo) / (1 - v_lo)) / 2
        upper_value = mp.iv.log((1 + v_hi) / (1 - v_hi)) / 2
        return _interval_endpoints(lower_value)[0], _interval_endpoints(upper_value)[1]
    finally:
        mp.iv.dps = previous


def _verify_stored_bracket(
    stored: dict[str, object], expected_lo: int, expected_hi: int
) -> None:
    assert int(stored["lower_numerator"]) == expected_lo
    assert int(stored["upper_numerator"]) == expected_hi
    assert int(stored["denominator"]) == SCALE
    assert expected_hi == expected_lo + 1


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}
    provenance = payload["provenance"]
    assert provenance["script"] == "experiments/e55_kesten_bound.py"
    assert provenance["interpreter"] == ".venv/bin/python"
    assert provenance["mpmath_iv_dps"] == 90
    assert provenance["rational_bracket_digits"] == DIGITS
    assert provenance["generated_utc"].endswith("+00:00")

    manifest_text = MANIFEST.read_text(encoding="utf-8")
    source_payload = payload["data"]["source_certificates"]
    for key, (relative_path, expected_hash) in EXPECTED_SOURCES.items():
        path = ROOT / relative_path
        assert _sha256(path) == expected_hash
        assert source_payload[key]["actual_sha256"] == expected_hash
        assert source_payload[key]["verified"] is True
        assert key in manifest_text
        assert expected_hash in manifest_text

    stored_counts = {
        int(n): int(value)
        for n, value in payload["data"]["exact_saw_counts"]["counts"].items()
    }
    assert stored_counts == EXPECTED_COUNTS

    candidates = payload["data"]["conditional_ratio_candidates"]
    assert [candidate["n"] for candidate in candidates] == [31, 32, 33, 34]
    for candidate in candidates:
        n = candidate["n"]
        c_n = EXPECTED_COUNTS[n]
        c_next = EXPECTED_COUNTS[n + 2]

        mu_lo, mu_hi = _sqrt_ratio_bracket_isqrt(c_next, c_n)
        v_lo, v_hi = _sqrt_ratio_bracket_isqrt(c_n, c_next)
        _verify_stored_bracket(
            candidate["conditional_mu_upper_bracket"], mu_lo, mu_hi
        )
        _verify_stored_bracket(candidate["conditional_v_bracket"], v_lo, v_hi)

        assert c_n * mu_lo * mu_lo <= c_next * SCALE * SCALE
        assert c_next * SCALE * SCALE < c_n * mu_hi * mu_hi
        assert c_next * v_lo * v_lo <= c_n * SCALE * SCALE
        assert c_n * SCALE * SCALE < c_next * v_hi * v_hi

        recomputed_lo, recomputed_hi = _directed_atanh(v_lo, v_hi)
        stored_lo = Decimal(candidate["conditional_k_directed_interval"][0])
        stored_hi = Decimal(candidate["conditional_k_directed_interval"][1])
        assert stored_lo <= recomputed_lo <= recomputed_hi <= stored_hi

        rounded_lo = Decimal(candidate["conditional_k_outward_rounded_40"][0])
        rounded_hi = Decimal(candidate["conditional_k_outward_rounded_40"][1])
        assert rounded_lo <= stored_lo
        assert stored_hi <= rounded_hi
        assert recomputed_lo > TASK_INCUMBENT_DISPLAY
        assert candidate["status"].startswith("[UNRESOLVED]")

    best_n = max(
        (31, 32, 33, 34),
        key=lambda n: Fraction(EXPECTED_COUNTS[n], EXPECTED_COUNTS[n + 2]),
    )
    assert best_n == 34
    conditional_best = payload["data"]["conditional_best"]
    assert conditional_best["n"] == best_n
    assert conditional_best["certified_for_ising"] is False
    assert Decimal(conditional_best["directed_interval"][0]) > TASK_INCUMBENT_DISPLAY

    assert (
        EXPECTED_COUNTS[36] * EXPECTED_COUNTS[32]
        < EXPECTED_COUNTS[34] * EXPECTED_COUNTS[34]
    )
    assert (
        EXPECTED_COUNTS[35] * EXPECTED_COUNTS[31]
        < EXPECTED_COUNTS[33] * EXPECTED_COUNTS[33]
    )

    audit = payload["data"]["theorem_audit"]
    assert audit["verdict"].startswith("[UNRESOLVED]")
    assert "Theorem 7.3.4(a)" in audit["madras_slade"]
    assert "equation (1.17)" in audit["open_lecture_notes"]
    assert "abs(c_(n+2)/c_n - mu^2)" in audit["quantitative_form"]
    assert audit["missing_for_upgrade"].startswith("[UNRESOLVED]")

    final_bound = payload["data"]["final_certified_lower_bound"]
    incumbent_lo, incumbent_hi = _inverse_count_root_bracket(
        EXPECTED_COUNTS[36], 36
    )
    _verify_stored_bracket(final_bound["v_exact_bracket"], incumbent_lo, incumbent_hi)
    recomputed_lo, recomputed_hi = _directed_atanh(incumbent_lo, incumbent_hi)
    stored_lo = Decimal(final_bound["directed_interval"][0])
    stored_hi = Decimal(final_bound["directed_interval"][1])
    assert stored_lo <= recomputed_lo <= recomputed_hi <= stored_hi
    outward_lo = Decimal(final_bound["outward_rounded_40"][0])
    outward_hi = Decimal(final_bound["outward_rounded_40"][1])
    assert outward_lo == EXPECTED_INCUMBENT_LOWER
    assert outward_lo <= stored_lo <= stored_hi <= outward_hi
    assert final_bound["improves_incumbent"] is False

    assert all(check["passed"] is True for check in payload["checks"])
    print("PASS")


if __name__ == "__main__":
    main()
