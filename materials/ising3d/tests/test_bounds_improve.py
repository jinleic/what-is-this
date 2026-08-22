"""Standalone independent checks for the improved critical-coupling bounds."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "improved_bounds.json"
INCUMBENT = ROOT / "results" / "bounds" / "kc_bounds.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
SAW_SOURCE_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
PEIERLS_SOURCE = ROOT / "sources" / "fulltext" / "bonati2014_peierls.pdf"
PEIERLS_SOURCE_SHA256 = "2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36"
C36 = 2941370856334701726560670
DPS = 80
BENCHMARK = Decimal("0.221654626")  # sanity check only


def _interval_endpoints(value: object) -> tuple[Decimal, Decimal]:
    text = str(value).strip()
    assert text.startswith("[") and text.endswith("]"), text
    lo, hi = text[1:-1].split(",", maxsplit=1)
    return Decimal(lo.strip()), Decimal(hi.strip())


def _recompute_saw_certificate() -> tuple[str, str, Decimal, Decimal]:
    """Rebuild the exact root bracket and interval without experiment helpers."""

    scale = 10**90
    target = scale**36
    lo, hi = 0, scale
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if C36 * mid**36 < target:
            lo = mid
        else:
            hi = mid
    assert C36 * lo**36 < target <= C36 * hi**36

    old = mp.iv.dps
    try:
        mp.iv.dps = DPS
        vlo = mp.iv.mpf(lo) / scale
        vhi = mp.iv.mpf(hi) / scale
        klo = mp.iv.log((1 + vlo) / (1 - vlo)) / 2
        khi = mp.iv.log((1 + vhi) / (1 - vhi)) / 2
        k_lower = _interval_endpoints(klo)[0]
        k_upper = _interval_endpoints(khi)[1]
    finally:
        mp.iv.dps = old
    return f"{lo}/{scale}", f"{hi}/{scale}", k_lower, k_upper

def _recompute_mu_interval() -> tuple[Decimal, Decimal]:
    old = mp.iv.dps
    try:
        mp.iv.dps = DPS
        enclosure = mp.iv.exp(mp.iv.log(mp.iv.mpf(C36)) / 36)
        return _interval_endpoints(enclosure)
    finally:
        mp.iv.dps = old


def _peierls_sign(p: int, q: int) -> int:
    """Sign of F(p/q)-1/2 after clearing the positive denominator."""

    return 3 * p**3 * (9 * q * q - 11 * p * q + 4 * p * p) - 8 * (
        q - p
    ) ** 3 * q * q


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    incumbent = json.loads(INCUMBENT.read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e33_bounds_improve.py"
    assert payload["provenance"]["precision"]["mpmath_iv_dps"] == DPS
    assert payload["checks"] and all(item["passed"] for item in payload["checks"])

    candidates = {item["id"]: item for item in payload["data"]["candidates"]}
    required = {
        "published_exact_saw_c36",
        "macdonald_mu_4_7114_source_audit",
        "finite_count_saw_c14",
        "simon_lieb_4x4x32",
        "griffiths_simon_finite_region",
        "infrared_watson",
        "peierls_surface_count",
    }
    assert required <= candidates.keys()
    for candidate in candidates.values():
        assert candidate["theorem"]
        assert "certificate" in candidate
        assert isinstance(candidate["improves_incumbent"], bool)
    print("candidate schema and theorem certificates: PASS")

    saw = candidates["published_exact_saw_c36"]
    assert saw["certificate"]["published_exact_count"] == {
        "steps": 36,
        "c_n": C36,
    }
    assert hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest() == SAW_SOURCE_SHA256
    assert hashlib.sha256(PEIERLS_SOURCE.read_bytes()).hexdigest() == PEIERLS_SOURCE_SHA256
    mu_lo, mu_hi = _recompute_mu_interval()
    stored_mu_lo, stored_mu_hi = map(
        Decimal, saw["certificate"]["mu_upper_interval_directed_rounding"]
    )
    assert stored_mu_lo <= mu_lo <= mu_hi <= stored_mu_hi
    assert stored_mu_hi < Decimal(
        incumbent["data"]["lower_bound_saw"]["mu_upper_interval"][0]
    )
    vlo, vhi, lo, hi = _recompute_saw_certificate()
    assert saw["certificate"]["v_exact_rational_bracket"] == [vlo, vhi]
    stored_lo, stored_hi = map(
        Decimal, saw["certificate"]["K_interval_directed_rounding"]
    )
    reported_lower = Decimal(saw["certificate"]["certified_decimal_lower"])
    assert stored_lo <= lo <= hi <= stored_hi
    assert reported_lower <= stored_lo
    assert stored_lo - reported_lower < Decimal("1e-40")
    rejected = candidates["macdonald_mu_4_7114_source_audit"]
    assert rejected["certificate"]["decision"] == "not used in any certified endpoint"
    assert rejected["improves_incumbent"] is False
    print("published exact c36 endpoint independently certified: PASS")

    old_lower, old_upper = map(
        Decimal,
        incumbent["data"]["final_certified_interval"]["decimal_outward_rounded"],
    )
    final_lower, final_upper = map(
        Decimal,
        payload["data"]["final_certified_interval"]["decimal_outward_rounded"],
    )
    assert final_lower == reported_lower > old_lower
    assert final_upper == old_upper
    interval_record = payload["data"]["final_certified_interval"]
    assert Decimal(interval_record["incumbent_width"]) == old_upper - old_lower
    assert Decimal(interval_record["final_width"]) == final_upper - final_lower
    assert (
        Decimal(interval_record["lower_endpoint_improvement"])
        == final_lower - old_lower
    )
    assert saw["improves_incumbent"] is True
    assert candidates["infrared_watson"]["improves_incumbent"] is False
    print("lower endpoint improves and upper endpoint remains certified: PASS")

    peierls = candidates["peierls_surface_count"]["certificate"]
    xlo_text, xhi_text = peierls["x_root_exact_rational_bracket"]
    xlo_p, xlo_q = map(int, xlo_text.split("/"))
    xhi_p, xhi_q = map(int, xhi_text.split("/"))
    assert xlo_q == xhi_q and xhi_p == xlo_p + 1
    assert _peierls_sign(xlo_p, xlo_q) < 0 <= _peierls_sign(xhi_p, xhi_q)
    peierls_upper = Decimal(peierls["certified_decimal_upper"])
    assert peierls_upper > final_upper
    prefix = peierls["exact_surface_prefix_from_repo"][
        "coefficients_area_0_to_28_nonzero"
    ]
    assert prefix == {
        "0": 1,
        "6": 36,
        "10": 75,
        "12": 555,
        "14": 250,
        "16": 2102,
        "18": 5697,
        "20": 8412,
        "22": 27806,
        "24": 58192,
        "26": 122144,
        "28": 276078,
    }
    print("Peierls rational root bracket and exact finite prefix: PASS")

    diagnostic = candidates["simon_lieb_4x4x32"]["diagnosis"][
        "boundary_orbit_diagnostic_at_certified_v"
    ]
    kappa = Decimal(diagnostic["kappa"])
    lateral = Decimal(diagnostic["lateral_faces_contribution"])
    caps = Decimal(diagnostic["end_caps_contribution"])
    assert abs(kappa - lateral - caps) < Decimal("2e-15")
    assert caps / kappa < Decimal("0.000004")
    assert Decimal(diagnostic["ten_largest_orbits_fraction"]) > Decimal("0.45")
    print("Simon--Lieb lateral-boundary saturation diagnosis: PASS")

    # Required comparison-only sanity check.  It does not feed any endpoint.
    assert final_lower < BENCHMARK < final_upper
    assert payload["data"]["final_certified_interval"]["benchmark_role"].startswith(
        "comparison-only"
    )
    print("benchmark containment (comparison only): PASS")
    print("PASS: improved rigorous-bounds certificate verified")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}")
        raise
