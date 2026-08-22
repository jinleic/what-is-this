"""Produce the all-L quadratic no-go certificate artifact.

The decisive data are exact nested-bracket evaluations from e160 and exact
integer dimension-to-mode arithmetic from e161.  The producer covers L=2..8;
the proof note supplies the uniform induction behind the compact schemas.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for candidate in (ROOT / "src", ROOT / "experiments"):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import e160_hamiltonian_clifford_family as family  # noqa: E402
import e161_quadratic_mode_bounds as bounds  # noqa: E402

OUTPUT = ROOT / "results" / "algebra_growth" / "alll_quadratic.json"
CPU_BUDGET_SECONDS = 90.0
EXACT_DIMENSIONS = {2: 56, 3: 1056, 4: 16256, 5: 262656}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def source_hashes() -> dict[str, str]:
    paths = [
        ROOT / "experiments" / "e160_hamiltonian_clifford_family.py",
        ROOT / "experiments" / "e161_quadratic_mode_bounds.py",
        ROOT / "experiments" / "e162_alll_quadratic.py",
        ROOT / "src" / "ising" / "clifford" / "__init__.py",
        ROOT / "proofs" / "allsize_gaussian.md",
        ROOT / "proofs" / "algebraic_obstruction.md",
        ROOT / "proofs" / "alternation_law.md",
    ]
    return {str(path.relative_to(ROOT)): sha256_file(path) for path in paths}


def build_payload(started: float) -> dict:
    deadline = started + CPU_BUDGET_SECONDS
    checks: list[dict] = []
    sizes: list[dict] = []

    for L in range(2, 9):
        certificate = family.certify_size(L, deadline)
        row = bounds.size_row(L)
        ceiling = row["physical_gaussian_ceiling"]

        check(
            checks,
            f"L{L}_local_term_inventory",
            certificate["local_term_count"] == 5 * L - 2,
            f"{certificate['local_term_count']} terms = 2L fields + (3L-2) bonds",
        )
        check(
            checks,
            f"L{L}_path_brackets",
            certificate["path_bilinears"]["count"] == ceiling,
            f"all {ceiling}=C(4L,2) path-bilinear words evaluated to nonzero distinct Paulis",
        )
        check(
            checks,
            f"L{L}_polynomial_no_go_family",
            certificate["polynomial_family_count"] == row["polynomial_family_bound"] == ceiling + 1,
            f"{certificate['polynomial_family_count']} > physical so(4L) ceiling {ceiling}",
        )
        check(
            checks,
            f"L{L}_ring_orbits",
            certificate["ring_family_count"] == row["ring_family_bound"] == 2 * ceiling,
            "grade-2 and grade-(4L-2) nested-bracket families are disjoint",
        )
        check(
            checks,
            f"L{L}_materialised_hypercubes",
            certificate["central_hypercube_count"] == row["materialised_hypercube_count"]
            and certificate["materialised_family_count"] == row["materialised_family_bound"],
            (
                f"every one of {certificate['central_hypercube_count']} central-hypercube members "
                "was evaluated; the stored materialised union is duplicate-free"
            ),
        )
        compact = certificate["compact_central_orbit"]
        check(
            checks,
            f"L{L}_compact_orbit_schema",
            compact["grade"] == row["compact_central_grade"]
            and compact["claimed_count"] == row["compact_central_orbit_bound"],
            (
                f"grade {compact['grade']}, claimed count {compact['claimed_count']}; "
                f"{compact['single_swap_schema_checks']} exact seed-swap checks"
            ),
        )
        check(
            checks,
            f"L{L}_strong_dimension_bound",
            row["strongest_dimension_bound"]
            == (row["ring_family_bound"] if L == 2 else row["ring_family_bound"] + compact["claimed_count"]),
            f"dim(g_L) >= {row['strongest_dimension_bound']}",
        )
        modes = row["minimum_modes_from_strongest_bound"]
        check(
            checks,
            f"L{L}_mode_minimality",
            row["predecessor_gaussian_dimension"] < row["strongest_dimension_bound"]
            <= row["attained_gaussian_dimension"],
            (
                f"m={modes - 1} gives {row['predecessor_gaussian_dimension']} < bound; "
                f"m={modes} gives {row['attained_gaussian_dimension']}"
            ),
        )
        check(
            checks,
            f"L{L}_binomial_growth_inequality",
            5 * (4 * L + 1) * compact["claimed_count"] >= (1 << (4 * L + 1)),
            "5(4L+1) C_L >= 2^(4L+1), exact integers",
        )

        certificate.update(row)
        sizes.append(certificate)

    finite_refinements: list[dict] = []
    for L, dimension in EXACT_DIMENSIONS.items():
        m = bounds.minimum_modes(dimension)
        finite_refinements.append(
            {
                "L": L,
                "exact_characteristic_zero_dimension": dimension,
                "minimum_modes_from_exact_dimension": m,
                "predecessor_gaussian_dimension": bounds.gaussian_lie_dimension(m - 1),
                "attained_gaussian_dimension": bounds.gaussian_lie_dimension(m),
                "source": "proofs/alternation_law.md (equality certified only at L=2..5)",
            }
        )
        check(
            checks,
            f"L{L}_certified_exact_refinement_consistent",
            dimension >= sizes[L - 2]["strongest_dimension_bound"]
            and bounds.gaussian_lie_dimension(m - 1) < dimension <= bounds.gaussian_lie_dimension(m),
            f"exact dim {dimension} implies m >= {m}; not used for the all-L theorem",
        )

    payload = {
        "meta": {
            "schema": "ising3d.alll_quadratic.v1",
            "status": "PASS" if all(item["passed"] for item in checks) else "FAIL",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "provenance": {
                "producer": "experiments/e162_alll_quadratic.py",
                "family_engine": "experiments/e160_hamiltonian_clifford_family.py",
                "bound_engine": "experiments/e161_quadratic_mode_bounds.py",
                "standalone_verifier": "tests/test_alll_quadratic.py",
                "source_sha256": source_hashes(),
            },
            "environment": {
                "interpreter_requested": ".venv/bin/python",
                "executable": sys.executable,
                "python_version": platform.python_version(),
                "python_implementation": platform.python_implementation(),
                "platform": platform.platform(),
                "machine": platform.machine(),
                "cwd": os.getcwd(),
                "pythonhashseed": os.environ.get("PYTHONHASHSEED"),
                "arithmetic": "exact Python integers and F_2 packed Pauli vectors; no floating-point decisions",
                "process_time_budget_seconds": CPU_BUDGET_SECONDS,
                "process_time_seconds": None,
            },
        },
        "data": {
            "statements": [
                {
                    "tag": "THEOREM",
                    "name": "dimension necessary condition",
                    "statement": (
                        "If all local terms embed into the quadratic algebra of m fermionic modes, "
                        "then dim(g_L) <= m(2m-1), hence m >= ceil((1+sqrt(1+8 dim(g_L)))/4)."
                    ),
                },
                {
                    "tag": "THEOREM",
                    "name": "all-L polynomial witness",
                    "statement": (
                        "For every L>=2 and n=2L, dim(g_L) >= n(2n-1)+1 > n(2n-1); "
                        "there is no same-physical-space quadratic realization."
                    ),
                },
                {
                    "tag": "THEOREM",
                    "name": "all-L quantitative central-rung bound",
                    "statement": (
                        "A central rung gives C_L=binom(4L,2L)/2 for odd L and "
                        "C_L=binom(4L,2L+2) for even L; for L>=3 this grade is disjoint from "
                        "the two perimeter grades, so dim(g_L)>=2(2L)(4L-1)+C_L."
                    ),
                },
                {
                    "tag": "THEOREM",
                    "name": "exponential mode growth",
                    "statement": (
                        "For every L>=3, any invariant-code quadratic realization needs at least "
                        "minimum_modes(strongest_dimension_bound) modes and in particular "
                        "m > 2^(2L)/sqrt(5(4L+1))."
                    ),
                },
                {
                    "tag": "UNRESOLVED",
                    "name": "larger-space sufficiency",
                    "statement": (
                        "The lower bound does not construct or exclude a realization once m meets it; "
                        "mere non-invariant compression into a larger Fock space is outside the Lie-embedding hypothesis."
                    ),
                },
            ],
            "bracket_recipes": {
                "path_bilinear": "H[a,a+1]=G[a]; H[a,b]=[G[b-1],H[a,b-1]].",
                "dual_bilinear": (
                    "Start at the closing rung with holes {1,4L}; bracket with H[1,a] if a>1, "
                    "then with H[b,4L] if b<4L."
                ),
                "johnson_target": (
                    "For seed support S0 and target S, pair sorted(S0\\S) with sorted(S\\S0); "
                    "successively bracket with the corresponding H[min,max]."
                ),
                "middle_grade_selection": (
                    "At odd L, retain one lexicographic representative from each pair {S,S^c}; "
                    "this gives exactly half the central binomial without assuming middle-grade saturation."
                ),
            },
            "finite_verification_range": [2, 8],
            "sizes": sizes,
            "certified_exact_finite_refinements": finite_refinements,
            "checks": checks,
        },
    }
    payload["meta"]["environment"]["process_time_seconds"] = round(time.process_time() - started, 9)
    return payload


def main() -> None:
    started = time.process_time()
    try:
        payload = build_payload(started)
    except family.NonDecisiveBudget as error:
        print(str(error))
        raise SystemExit(2) from error
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checks = payload["data"]["checks"]
    failed = [item["name"] for item in checks if not item["passed"]]
    if failed:
        print(f"FAIL: {', '.join(failed)}")
        raise SystemExit(1)
    print(f"artifact={OUTPUT.relative_to(ROOT)} checks={len(checks)}")
    print("PASS")


if __name__ == "__main__":
    main()
