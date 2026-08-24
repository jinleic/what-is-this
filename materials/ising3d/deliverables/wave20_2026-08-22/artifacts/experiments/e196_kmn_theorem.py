#!/usr/bin/env python3
"""[THEOREM][COMPUTATION] Integrate K_{m,n} reduction, exact anchors, and scope limits."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from e194_kmn_reduction import (
    CPU_BUDGET_SECONDS,
    GOOD_PRIMES,
    RSS_CAP_BYTES,
    balanced_family_upper,
    budget_tick,
    is_prime,
    max_rss_bytes,
    opposite_parity_family_upper,
    reduction_audit,
)
from e195_kmn_closure import case_closure

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "kmn_twogen.json"
CASE_SPECS = ((2, 3), (3, 3), (3, 4), (4, 4), (4, 5))
REGRESSION_DIMENSIONS = {"K_{2,3}": 44, "K_{3,3}": 63, "K_{3,4}": 167}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_check(name: str, passed: object, detail: object) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def largest_opposite_parity_sector_dimension(m: int) -> int:
    return ((m + 1) * (m + 2) // 2) ** 2


def threshold_polynomial_numerator(m: int) -> int:
    """Four times q_max^2-(2m+1)(4m+1)."""
    return m * (m**3 + 6 * m**2 - 19 * m - 12)


def family_tables(maximum_m: int = 8) -> dict[str, object]:
    opposite_rows = []
    balanced_rows = []
    for m in range(1, maximum_m + 1):
        opposite = opposite_parity_family_upper(m, m + 1)
        vertices = 2 * m + 1
        ceiling = vertices * (2 * vertices - 1)
        largest = largest_opposite_parity_sector_dimension(m)
        opposite_rows.append(
            {
                "tag": "[THEOREM]",
                "graph": f"K_{{{m},{m + 1}}}",
                "m": m,
                "upper_dimension": opposite["upper_dimension"],
                "quadratic_ceiling": ceiling,
                "largest_spin_sector_container_dimension": largest,
                "upper_is_below_ceiling": opposite["upper_dimension"] < ceiling,
                "largest_sector_already_exceeds_ceiling": largest > ceiling,
                "four_times_largest_minus_ceiling": threshold_polynomial_numerator(m),
            }
        )
        if m >= 2:
            balanced = balanced_family_upper(m)
            balanced_vertices = 2 * m
            balanced_ceiling = balanced_vertices * (2 * balanced_vertices - 1)
            balanced_rows.append(
                {
                    "tag": "[THEOREM]",
                    "graph": f"K_{{{m},{m}}}",
                    "m": m,
                    "upper_dimension": balanced["upper_dimension"],
                    "quadratic_ceiling": balanced_ceiling,
                    "upper_is_below_ceiling": balanced["upper_dimension"]
                    < balanced_ceiling,
                    "semisimple_upper": balanced["semisimple_upper"],
                    "global_center_upper": balanced["global_center_upper"],
                }
            )
    return {
        "tag": "[THEOREM]",
        "opposite_parity_closed_form": (
            "dim h_{m,n} <= 1 + S_m*S_n/4 - a_m*a_n for m,n of opposite "
            "parity, S_t=C(t+3,3)-1_{t even}, a_t=ceil(t/2)"
        ),
        "hamiltonian_near_balanced_rows": opposite_rows,
        "balanced_representation_rows": balanced_rows,
        "threshold_method_limit": {
            "tag": "[THEOREM]",
            "statement": (
                "For K_{m,m+1}, m>=3, even the largest-spin gl(q) container has "
                "q^2>(2m+1)(4m+1); therefore this collective-symmetry upper-bound "
                "method cannot force the two-sum dimension below the quadratic ceiling."
            ),
            "q": "(m+1)(m+2)/2",
            "four_times_difference": "m(m^3+6m^2-19m-12)>0 for m>=3",
            "scope": (
                "This is a limitation of the proved upper container, not an all-size "
                "lower bound and not a decision of the universal-threshold question."
            ),
        },
    }


def run_all() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    reduction = reduction_audit()
    cases = [case_closure(left, right) for left, right in CASE_SPECS]
    by_graph = {str(row["graph"]): row for row in cases}
    families = family_tables()
    new_anchors = [by_graph["K_{4,4}"], by_graph["K_{4,5}"]]
    regression_values = {
        graph: by_graph[graph]["dimension_Q"] for graph in REGRESSION_DIMENSIONS
    }
    opposite_rows = families["hamiltonian_near_balanced_rows"]
    polynomial_positive = all(
        threshold_polynomial_numerator(m) > 0
        and threshold_polynomial_numerator(m + 1)
        > threshold_polynomial_numerator(m)
        for m in range(3, 64)
    )
    checks = [
        make_check(
            "collective-spin matrix identities",
            reduction["all_square_sums_pass"]
            and reduction["all_sector_identities_pass"],
            {
                "square_sum_rows": len(reduction["square_sum_rows"]),
                "sector_rows": len(reduction["sector_identity_rows"]),
            },
        ),
        make_check(
            "finite-field moduli are prime",
            all(is_prime(prime) for prime in GOOD_PRIMES),
            list(GOOD_PRIMES),
        ),
        make_check(
            "all audited complete-bipartite graphs are Hamiltonian and branching",
            all(row["hamiltonian"] and row["maximum_degree"] >= 3 for row in cases),
            [row["graph"] for row in cases],
        ),
        make_check(
            "two modular lower bounds meet the representation upper in every case",
            all(
                len(row["modular_ranks"]) == len(GOOD_PRIMES)
                and len(set(row["modular_ranks"])) == 1
                and row["lower_equals_upper"]
                and row["dimension_Q"] == row["structural_upper"]["upper_dimension"]
                for row in cases
            ),
            {row["graph"]: row["modular_ranks"] for row in cases},
        ),
        make_check(
            "prior exact complete-bipartite anchors regress",
            regression_values == REGRESSION_DIMENSIONS,
            regression_values,
        ),
        make_check(
            "new K4,4 and K4,5 anchors are exact",
            {row["graph"]: row["dimension_Q"] for row in new_anchors}
            == {"K_{4,4}": 137, "K_{4,5}": 471},
            {
                row["graph"]: {
                    "lower": row["modular_ranks"],
                    "upper": row["structural_upper"]["upper_dimension"],
                }
                for row in new_anchors
            },
        ),
        make_check(
            "both new anchors clear the physical quadratic ceiling",
            all(row["clears_quadratic_ceiling"] for row in new_anchors),
            {
                row["graph"]: [row["dimension_Q"], row["quadratic_ceiling"]]
                for row in new_anchors
            },
        ),
        make_check(
            "K4,5 bypasses the legacy 256-coordinate rational wall",
            by_graph["K_{4,5}"]["reduced_faithful_coordinates"] > 256
            and all(
                closure["saturated_over_Fp"]
                for closure in by_graph["K_{4,5}"]["modular_closures"]
            )
            and by_graph["K_{4,5}"]["lower_equals_upper"],
            {
                "legacy_wall": 256,
                "representation_coordinates": by_graph["K_{4,5}"][
                    "reduced_faithful_coordinates"
                ],
                "rank": by_graph["K_{4,5}"]["dimension_Q"],
            },
        ),
        make_check(
            "modular directions are lower-bound directions",
            all(
                closure["rank_lower_bound_over_Q"] <= row["dimension_Q"]
                for row in cases
                for closure in row["modular_closures"]
            ),
            "rank_Fp <= rank_Q; equality is obtained only after the independent upper bound",
        ),
        make_check(
            "local-term and two-sum algebras are strictly distinguished",
            all(row["local_term_is_strictly_larger"] for row in cases),
            {
                row["graph"]: [row["dimension_Q"], row["local_term_dimension"]]
                for row in cases
            },
        ),
        make_check(
            "near-balanced symmetry upper stops lying below the ceiling after m=2",
            [
                row["m"]
                for row in opposite_rows
                if row["upper_is_below_ceiling"]
            ]
            == [1, 2]
            and all(
                row["largest_sector_already_exceeds_ceiling"]
                for row in opposite_rows
                if row["m"] >= 3
            )
            and polynomial_positive,
            {
                row["graph"]: [row["upper_dimension"], row["quadratic_ceiling"]]
                for row in opposite_rows
            },
        ),
        make_check(
            "resource cap",
            max_rss_bytes() < RSS_CAP_BYTES,
            {
                "peak_rss_bytes": max_rss_bytes(),
                "rss_cap_bytes": RSS_CAP_BYTES,
                "process_time_seconds": round(time.process_time() - started, 6),
                "process_time_budget_seconds": CPU_BUDGET_SECONDS,
            },
        ),
    ]
    if not all(bool(row["passed"]) for row in checks):
        raise AssertionError(checks)
    budget_tick(started, "integrated K_mn theorem audit")
    data = {
        "claim_tags": [
            "[THEOREM]",
            "[LEMMA]",
            "[COMPUTATION]",
            "[CONJECTURE]",
            "[EXTERNAL]",
            "[UNRESOLVED]",
        ],
        "headline": {
            "tag": "[THEOREM][COMPUTATION]",
            "statement": (
                "Schur-Weyl spin blocks plus parity, bilinear-form, and part-swap "
                "constraints give characteristic-zero upper bounds met by modular Lie "
                "words: dim h(K4,4)=137 and dim h(K4,5)=471."
            ),
        },
        "representation_reduction": reduction,
        "complete_bipartite_cases": cases,
        "family_bounds": families,
        "local_term_contrast": {
            "tag": "[THEOREM][EXTERNAL][COMPUTATION]",
            "statement": (
                "The collective two-sum algebra is polynomially symmetry-reduced; the "
                "independently generated local-term algebra is the branching Clifford "
                "grade algebra and is exponentially larger in every audited row."
            ),
            "rows": [
                {
                    "tag": "[COMPUTATION]",
                    "graph": row["graph"],
                    "two_sum_dimension": row["dimension_Q"],
                    "local_term_dimension": row["local_term_dimension"],
                }
                for row in cases
            ],
        },
        "universal_threshold_status": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "Complete-bipartite symmetry supplies no eventual counterexample: the "
                "proved container already exceeds the quadratic ceiling for K_{m,m+1}, "
                "m>=3, and the exact K3,4, K4,4, K4,5 anchors clear it. No all-size "
                "lower bound or universal n0 follows."
            ),
        },
        "scope": {
            "tag": "[THEOREM]",
            "statement": (
                "These are finite-dimensional Lie-algebra statements for A=sum X and "
                "B=sum ZZ. They do not compute thermodynamic free energy, critical "
                "coupling, or solve the three-dimensional Ising model."
            ),
        },
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }
    return data, checks


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]]
) -> dict[str, object]:
    sources = (
        "experiments/e194_kmn_reduction.py",
        "experiments/e195_kmn_closure.py",
        "experiments/e196_kmn_theorem.py",
        "tests/test_kmn_twogen.py",
        "proofs/kmn_twogen.md",
    )
    return {
        "meta": {
            "experiment": "e196_kmn_theorem",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e196_kmn_theorem.py",
            "working_directory": str(ROOT),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact integral Sym^r matrices; exact representation-theoretic upper "
                "dimensions; integral Lie words ranked over two explicitly checked primes"
            ),
            "modular_direction": (
                "each finite-field rank is only a lower bound on rational rank; exactness "
                "is asserted only when that lower bound meets the proved characteristic-zero upper"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "single_process": True,
            },
            "source_sha256": {
                source: sha256_file(ROOT / source) for source in sources
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks = run_all()
    artifact = make_artifact(data, checks)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    exact = {
        row["graph"]: row["dimension_Q"]
        for row in data["complete_bipartite_cases"]
        if row["graph"] in {"K_{4,4}", "K_{4,5}"}
    }
    print(
        f"PASS e196: {exact}, {len(checks)} checks, "
        f"cpu={data['process_time_seconds']}s, rss={data['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
