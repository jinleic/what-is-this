"""Integrate the exact Clifford-grade refutation and corrected trichotomy."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from math import isqrt
from pathlib import Path

from e178_clifford_setup import (
    RSS_CAP_BYTES,
    grade_two_mod_four_dimension,
    max_rss_bytes,
    run_setup,
)
from e179_grade_generation import run_generation_audit
from e180_graph_verification import run_graph_audit

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "clifford_grade.json"
CPU_BUDGET_SECONDS = 240.0


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimum_quadratic_modes(dimension: int) -> int:
    """Smallest integer m satisfying dimension <= m(2m-1)."""
    if dimension < 0:
        raise ValueError(dimension)
    estimate = max(0, (1 + isqrt(1 + 8 * dimension)) // 4 - 1)
    while estimate * (2 * estimate - 1) < dimension:
        estimate += 1
    while estimate and (estimate - 1) * (2 * estimate - 3) >= dimension:
        estimate -= 1
    return estimate


def run_all() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    setup = run_setup()
    generation = run_generation_audit()
    graph_audit = run_graph_audit()
    exact_by_label = {row["label"]: row for row in graph_audit["exact_rows"]}
    characterized_by_label = {
        row["label"]: row for row in graph_audit["characterization_rows"]
    }

    c6 = exact_by_label["cycle_C6_counterexample"]
    c8 = exact_by_label["cycle_C8_counterexample"]
    c4 = exact_by_label["known_grid_2x2_C4"]
    nonbip = exact_by_label["nonbipartite_C5_plus_chord_0_3"]
    grid_3x4 = characterized_by_label["new_grid_3x4"]

    verdict = {
        "tag": "[THEOREM]",
        "target_statement": "REFUTED",
        "load_bearing_counterexample": {
            "graph": "C6",
            "hypotheses": (
                "simple, bipartite, Hamiltonian; its Hamiltonian path has one non-path edge"
            ),
            "chord_grade": 10,
            "N_majoranas": 12,
            "actual_grade_set": [2, 10],
            "actual_dimension": 132,
            "proposed_grade_set": [2, 6, 10],
            "proposed_dimension": 1056,
            "exact_reason": (
                "two 10-subsets of a 12-set overlap in at least 8 places; only odd overlap 9 "
                "can give a nonzero bracket, and its output has grade 2"
            ),
            "raw_generator_computation_sha256": c6["closure_sha256"],
        },
        "corrected_trichotomy": [
            {
                "branch": "path",
                "graph_condition": "Delta(Gamma)<=2 and no chord",
                "grade_set": [2],
                "dimension": "n(2n-1)=C(2n,2)",
                "interpretation": "open-chain quadratic algebra so(2n)",
            },
            {
                "branch": "even_cycle",
                "graph_condition": "Delta(Gamma)=2 and the sole chord joins path endpoints",
                "grade_set": [2, "2n-2"],
                "dimension": "2*C(2n,2)=2n(2n-1)",
                "interpretation": "ring/parity-doubled quadratic-size algebra",
            },
            {
                "branch": "branching",
                "graph_condition": "Delta(Gamma)>=3, equivalently some chord is not the endpoint pair",
                "grade_set": "every noncentral grade k=2 mod 4",
                "dimension_even_n": "2^(2n-2)-(-1)^(n/2)2^(n-1)",
                "dimension_odd_n": "2^(2n-2)-1",
                "interpretation": "exponential Clifford-reversal algebra",
            },
        ],
        "threshold_equivalence": (
            "For a graph carrying a Hamiltonian path, no chord gives the path; only the endpoint "
            "chord gives a cycle; any other chord meets an interior path vertex and forces "
            "Delta>=3. Conversely Delta<=2 forces path or cycle."
        ),
    }

    obligations = [
        {
            "number": 1,
            "tag": "[LEMMA]",
            "status": "PROVED",
            "result": (
                "With gamma_(2r)=X_<r Z_r and gamma_(2r+1)=X_<r Y_r, "
                "X_r=i gamma_(2r)gamma_(2r+1) and Z_pZ_q=i^(q-p) "
                "gamma_(2p+1)...gamma_(2q), of grade 2(q-p)."
            ),
        },
        {
            "number": 2,
            "tag": "[LEMMA]+[COMPUTATION]",
            "status": "PROVED_AND_COMPUTED",
            "result": (
                "Bipartite chords have odd path distance and grade 2 mod 4. A non-bipartite graph "
                "has a same-path-parity chord. C5 plus chord (0,3) closes exactly to dimension "
                "510, versus 255 for the 2-mod-4 class and 510 for the full noncentral even class."
            ),
        },
        {
            "number": 3,
            "tag": "[THEOREM]",
            "status": "REFUTED_AS_STATED; CORRECTED_VERSION_PROVED",
            "result": (
                "The claim fails at k=N-2. For 6<=k<=N-4, overlap k-3 of two k-monomials "
                "produces grade 6, and overlap 1 with grade 6 raises r to r+4."
            ),
        },
        {
            "number": 4,
            "tag": "[LEMMA]",
            "status": "PROVED",
            "result": (
                "Grade-2 brackets perform single-index swaps, so the Johnson graph connects all "
                "coordinate k-monomials. In middle grade the seed and its disjoint Hodge complement "
                "are connected; equivalently both nonisomorphic Hodge projections are nonzero."
            ),
        },
        {
            "number": 5,
            "tag": "[THEOREM]",
            "status": "REQUESTED_SCALAR_RATIONALE_REFUTED; CORRECTION_PROVED",
            "result": (
                "The volume is not scalar on the full 2^n-dimensional spinor space: it is a phase "
                "times global X parity and is scalar only on each half-spin sector. It is excluded "
                "because it is central in Cl^even, absent from the generators, and cannot be a "
                "nonzero even-monomial commutator."
            ),
        },
        {
            "number": 6,
            "tag": "[LEMMA]",
            "status": "PROVED",
            "result": (
                "The noncentral 2-mod-4 span is bracket-closed and contains every bipartite local "
                "generator; this gives containment independently of the generation inclusion."
            ),
        },
        {
            "number": 7,
            "tag": "[COMPUTATION]",
            "status": "COMPLETED",
            "result": (
                "Raw-generator closures reproduce 56,1056,16256,262656,65535; new 3x3-minus-edge "
                "also gives 65535; 3x4 and C10-plus-chord are certified non-materialized predictions; "
                "P6, C6, C8, and the non-bipartite control are exact closures."
            ),
        },
        {
            "number": 8,
            "tag": "[THEOREM]",
            "status": "UNIVERSAL_EXPONENTIAL_CLAIM_REFUTED; BRANCHWISE CONSEQUENCE_PROVED",
            "result": (
                "Branching graphs require exponentially many quadratic modes by the dimension bound. "
                "Cycles have only quadratic dimension and hence only a linear mode lower bound; they "
                "are the ring/free-fermion branch. Exact dimensions still give the n-mode full-space "
                "no-go for every cycle."
            ),
        },
    ]

    selected_dimensions = [
        ("P6", 6, exact_by_label["path_P6_control"]["dimension"], "path"),
        ("C4=2x2", 4, c4["dimension"], "even_cycle"),
        ("C6", 6, c6["dimension"], "even_cycle"),
        ("C8", 8, c8["dimension"], "even_cycle"),
        ("2x3", 6, exact_by_label["known_grid_2x3"]["dimension"], "branching"),
        ("3x3", 9, exact_by_label["known_grid_3x3"]["dimension"], "branching"),
        ("2x5", 10, exact_by_label["known_grid_2x5"]["dimension"], "branching"),
        ("3x4", 12, grid_3x4["dimension"], "branching"),
    ]
    mode_bounds = [
        {
            "graph": label,
            "n": n,
            "branch": branch,
            "dimension": dimension,
            "n_mode_gaussian_bound": n * (2 * n - 1),
            "minimum_modes_from_dimension": minimum_quadratic_modes(dimension),
        }
        for label, n, dimension, branch in selected_dimensions
    ]

    consequences = {
        "tag": "[THEOREM]",
        "dimension_bound": "a Lie algebra of m-mode Majorana bilinears has dimension at most m(2m-1)",
        "exact_mode_formula": "m_min(D)=ceil((1+sqrt(1+8D))/4)",
        "mode_bounds": mode_bounds,
        "branching_asymptotic": (
            "D=Theta(4^n), hence m_min=Theta(2^n) for Delta(Gamma)>=3"
        ),
        "cycle_asymptotic": (
            "D=2n(2n-1)=Theta(n^2), hence the dimension argument gives only m_min=Theta(n)"
        ),
        "physical_space_scope": (
            "A unitary change of basis on the n-mode physical Fock space would keep m=n. "
            "Every cycle has D=2n(2n-1)>n(2n-1), and every branching case is larger still, "
            "so neither is globally n-mode quadratic; the ring remains solvable after parity-sector "
            "resolution or a linear-size enlargement."
        ),
        "claw_comparison": (
            "The induced-claw/frustration theorem already gives a yes/no term-wise obstruction at "
            "Delta>=3. The grade route adds exact exponential dimensions and quantitative mode counts. "
            "At Delta<=2 it computes the path and cycle algebras, including C4 (dimension 56) where "
            "the claw obstruction does not apply."
        ),
    }

    top_checks = [
        {
            "name": "all Jordan-Wigner phase checks",
            "passed": all(bool(row["passed"]) for row in setup["checks"]),
            "detail": f"{len(setup['checks'])} exact checks",
        },
        {
            "name": "all grade-generation checks",
            "passed": all(bool(row["passed"]) for row in generation["checks"]),
            "detail": f"{len(generation['checks'])} exact checks",
        },
        {
            "name": "all graph checks",
            "passed": all(bool(row["passed"]) for row in graph_audit["checks"]),
            "detail": f"{len(graph_audit['checks'])} exact checks",
        },
        {
            "name": "load-bearing C6 counterexample",
            "passed": c6["dimension"] == 132 and grade_two_mod_four_dimension(6) == 1056,
            "detail": "132 != 1056",
        },
        {
            "name": "C4 is cycle branch with accidental full grade class",
            "passed": c4["branch"] == "even_cycle" and c4["dimension"] == 56,
            "detail": f"grade histogram={c4['grade_histogram']}",
        },
        {
            "name": "non-bipartite necessity test",
            "passed": nonbip["dimension"] == 510,
            "detail": "C5 plus chord: 510 vs 255",
        },
        {
            "name": "new n=12 prediction",
            "passed": grid_3x4["dimension"] == 4_192_256,
            "detail": "non-materialized exact grade characterization",
        },
        {
            "name": "resource cap",
            "passed": max_rss_bytes() < RSS_CAP_BYTES,
            "detail": f"peak_rss_bytes={max_rss_bytes()}",
        },
    ]
    if not all(bool(row["passed"]) for row in top_checks):
        raise AssertionError("integrated Clifford-grade checks failed")

    data = {
        "verdict": verdict,
        "obligations": obligations,
        "jordan_wigner_setup": setup,
        "grade_generation": generation,
        "graph_verification": graph_audit,
        "centre_and_representation_convention": {
            "tag": "[THEOREM]",
            "matrix_space": "full 2^n-dimensional physical spinor space",
            "dimension_convention": (
                "number of linearly independent anti-Hermitian Pauli strings in the generated Lie "
                "algebra; identity and any nongenerated central parity string are excluded"
            ),
            "volume_identity": "P_X=product_r X_r=i^n gamma_0...gamma_(2n-1)",
            "full_space_scalar": False,
            "half_spin_scalar": True,
        },
        "nonbipartite_test": nonbip,
        "consequences": consequences,
        "integrated_process_time_seconds": round(time.process_time() - started, 6),
        "integrated_peak_rss_bytes": max_rss_bytes(),
    }
    return data, top_checks


def make_artifact(data: dict[str, object], checks: list[dict[str, object]]) -> dict[str, object]:
    sources = [
        ROOT / "experiments" / "e178_clifford_setup.py",
        ROOT / "experiments" / "e179_grade_generation.py",
        ROOT / "experiments" / "e180_graph_verification.py",
        Path(__file__),
    ]
    return {
        "meta": {
            "experiment": "e181_clifford_grade",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e181_clifford_grade.py",
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact Pauli GF(2) support closure, exact powers-of-i word multiplication, integer "
                "binomial counts, and finite Clifford-support overlap calculus"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in sources
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
    print(f"wrote {ARTIFACT.relative_to(ROOT)}")
    print(
        json.dumps(
            {
                "checks": len(checks),
                "process_time_seconds": data["integrated_process_time_seconds"],
                "peak_rss_bytes": data["integrated_peak_rss_bytes"],
            }
        )
    )
    print("PASS e181 Clifford-grade refutation and trichotomy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
