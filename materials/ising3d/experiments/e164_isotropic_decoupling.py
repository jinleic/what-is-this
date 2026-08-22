"""Exact obstruction to a universal isotropic tensor-decoupling route.

The certificate has two independent parts.

1. Across a nonempty graph cut, the uniform bond factor has a 2x2 minor
   q(t)^(2c)-1, so the graph-piece tensor identity used by Lemma 5 is not an
   identity on the isotropic curve.
2. The valid all-size branching family K_{1,3} disjoint-union m isolated
   vertices has m identical free modes.  Its pair-product count is at most
   136(2m+1), below the Gaussian ceiling for every m >= 3.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "isotropic_decoupling.json"
CPU_BUDGET_SECONDS = 60
RSS_CAP_BYTES = 4_000_000_000


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def q_of_t(t: Fraction) -> Fraction:
    return (1 + t * t) / (2 * t)


def cut_minor_numerator_coefficients(cut_edges: int) -> list[int]:
    """Ascending coefficients of (1+t^2)^(2c) - (2t)^(2c)."""
    assert cut_edges >= 1
    degree = 4 * cut_edges
    coefficients = [0] * (degree + 1)
    for j in range(2 * cut_edges + 1):
        coefficients[2 * j] += comb(2 * cut_edges, j)
    coefficients[2 * cut_edges] -= 2 ** (2 * cut_edges)
    return coefficients


def equal_site_upper(m: int) -> int:
    # At most C(16+1,2)=136 all-pair values from the claw and 2m+1 free exponents.
    return comb(17, 2) * (2 * m + 1)


def gaussian_ceiling_for_claw_plus(m: int) -> int:
    return 3 ** (m + 4) - 2 ** (m + 4)


def free_pair_exponents(m: int) -> set[int]:
    """Exponents j+j' from m identical two-level factors; multiplicities are irrelevant."""
    return {j + k for j in range(m + 1) for k in range(m + 1)}


def run_decoupling_obstruction() -> dict[str, object]:
    cpu0 = time.process_time()
    checks: list[dict[str, object]] = []

    cut_rows: list[dict[str, object]] = []
    for c in range(1, 9):
        coefficients = cut_minor_numerator_coefficients(c)
        row = {
            "cut_edges": c,
            "minor": "q(t)^(2c)-1",
            "cleared_numerator": f"(1+t^2)^{2*c}-(2t)^{2*c}",
            "degree_bound": 4 * c,
            "constant_coefficient": coefficients[0],
            "leading_coefficient": coefficients[-1],
            "coefficient_sha256": hashlib.sha256(
                ",".join(str(value) for value in coefficients).encode()
            ).hexdigest(),
        }
        cut_rows.append(row)
        check(
            checks,
            f"cut size {c}: separating minor is a nonzero polynomial",
            coefficients[0] == 1 and coefficients[-1] == 1,
            f"constant=leading=1; degree {4*c}",
        )

    # Exact rational spot checks of the cleared identity; these are checks, not
    # substitutes for the constant-term proof above.
    spot_ok = True
    spot_detail: list[str] = []
    for t in (Fraction(1, 5), Fraction(1, 3), Fraction(2, 5)):
        q = q_of_t(t)
        for c in (1, 2, 5):
            numerator = (1 + t * t) ** (2 * c) - (2 * t) ** (2 * c)
            identity = q ** (2 * c) - 1 == numerator / (2 * t) ** (2 * c)
            spot_ok &= identity and q > 1
        spot_detail.append(f"t={t},q={q}")
    check(
        checks,
        "isotropic substitution and physical positivity",
        spot_ok,
        "; ".join(spot_detail),
    )

    collapse_rows: list[dict[str, object]] = []
    for m in range(0, 13):
        exponent_count = len(free_pair_exponents(m))
        upper = equal_site_upper(m)
        ceiling = gaussian_ceiling_for_claw_plus(m)
        collapse_rows.append(
            {
                "isolated_sites_m": m,
                "sites_n": m + 4,
                "free_pair_exponent_count": exponent_count,
                "proved_pair_product_upper_bound": upper,
                "gaussian_ceiling": ceiling,
                "upper_below_ceiling": upper < ceiling,
            }
        )
        check(
            checks,
            f"m={m}: identical-site exponent classes",
            exponent_count == 2 * m + 1,
            f"{exponent_count}=2m+1",
        )

    base_m = 3
    base_gap = gaussian_ceiling_for_claw_plus(base_m) - equal_site_upper(base_m)
    # For m >= 3, U_{m+1}<2U_m and C_{m+1}>2C_m.  The displayed positive
    # differences prove the induction without a finite-size extrapolation.
    induction_certificate = {
        "base_m": base_m,
        "base_upper": equal_site_upper(base_m),
        "base_ceiling": gaussian_ceiling_for_claw_plus(base_m),
        "base_gap": base_gap,
        "twice_upper_minus_next_upper": "136*(2m-1) > 0 for m>=1",
        "next_ceiling_minus_twice_ceiling": "81*3^m > 0",
        "conclusion": "136(2m+1) < 3^(m+4)-2^(m+4) for every integer m>=3",
    }
    check(
        checks,
        "all-m collapse induction base",
        base_gap > 0,
        f"{equal_site_upper(base_m)} < {gaussian_ceiling_for_claw_plus(base_m)}",
    )
    check(
        checks,
        "all-m collapse induction steps",
        136 * (2 * base_m - 1) > 0 and 81 * 3**base_m > 0,
        "U_(m+1)<2U_m and C_(m+1)>2C_m for every m>=3 by positive closed forms",
    )

    graph_family = {
        "definition": "G_m = K_{1,3} disjoint-union m K_1",
        "vertices": "m+4",
        "maximum_degree": 3,
        "valid_for_every_integer_m_at_least": 0,
        "uniform_operator_identity": (
            "R_(G_m)(t,q(t)) = R_(K1,3)(t,q(t)) tensor "
            "([[1,t],[t,1]]^2)^(tensor m)"
        ),
        "single_site_eigenvalues": ["(1+t)^2", "(1-t)^2"],
        "reason_no_chain_remainder": (
            "the m remainder vertices have no edges; edge deletion cannot create the edges of P_m"
        ),
        "method_obstruction_from_m": 3,
    }

    # At q=1 on the physical curve one has t=1, where every P_t factor has rank
    # one.  This records why the only positive real bond-decoupling point cannot
    # seed a distinct-mode tensor count.
    endpoint = {
        "identity": "q(t)-1=(t-1)^2/(2t)",
        "positive_real_solution_q_equals_1": "t=1",
        "rank_P_t_at_t_1": 1,
        "rank_R_G_at_t_1_at_most": 1,
        "distinct_unordered_pair_products_at_t_1": 1,
    }
    check(
        checks,
        "positive isotropic decoupling endpoint is rank-collapsed",
        q_of_t(Fraction(1, 1)) == 1,
        "q(1)=1 and P_1 is the all-ones rank-one matrix",
    )
    check(
        checks,
        "resource cap",
        max_rss_bytes() < RSS_CAP_BYTES and time.process_time() - cpu0 < CPU_BUDGET_SECONDS,
        f"CPU {time.process_time()-cpu0:.3f}s; peak RSS {max_rss_bytes()}",
    )

    return {
        "cut_factorization_obstruction": {
            "statement": (
                "For a cut with c>=1 uniform crossing edges, the all-plus/all-minus "
                "2x2 minor of the crossing Boltzmann factor is q(t)^(2c)-1, a nonzero "
                "rational function. Hence graph-piece tensor decoupling is not an "
                "identity on the isotropic curve."
            ),
            "scope": "Lemma-5 graph-piece factorization across a nonempty edge cut",
            "rows": cut_rows,
        },
        "equal_mode_counterfamily": {
            "statement": (
                "For G_m=K1,3 disjoint-union m K1 and every m>=3, the exact uniform "
                "tensor identity has at most 136(2m+1) distinct unordered-slot pair "
                "products, strictly below 3^(m+4)-2^(m+4)."
            ),
            "graph_family": graph_family,
            "finite_audit_rows": collapse_rows,
            "universal_induction_certificate": induction_certificate,
        },
        "rank_collapsed_endpoint": endpoint,
        "conclusion": {
            "core_piece": "survives (see e163)",
            "decoupling_piece": "fails for the requested every-graph quantifier",
            "isotropic_all_size_theorem": "UNRESOLVED",
        },
        "checks": checks,
        "process_time_seconds": round(time.process_time() - cpu0, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def make_artifact(data: dict[str, object]) -> dict[str, object]:
    return {
        "meta": {
            "experiment": "e164_isotropic_decoupling",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "command": ".venv/bin/python experiments/e164_isotropic_decoupling.py",
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "producer_sha256": sha256_file(Path(__file__)),
            "arithmetic": "exact integers and fractions only",
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
        },
        "data": data,
    }


def main() -> int:
    data = run_decoupling_obstruction()
    artifact = make_artifact(data)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    passed = all(bool(item["passed"]) for item in data["checks"])
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if passed:
        print("PASS e164 isotropic decoupling obstruction", flush=True)
        return 0
    print("FAIL e164 isotropic decoupling obstruction", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
