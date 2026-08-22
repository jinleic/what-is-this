"""Exact bond-dimension-two tetrahedron-equation test for the 3D Ising site tensor."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ising.tensor_networks import ising_tetrahedron_analysis


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "tetrahedron.json"


def main() -> None:
    analysis = ising_tetrahedron_analysis()
    isotropic = analysis["isotropic_substitution"]
    gauge = analysis["arbitrary_leg_gl2_gauge"]
    metric_gauge = analysis["delta_preserving_identical_tensor_gauge"]
    checks = [
        {
            "name": "constant tetrahedron equation was specialized nontrivially",
            "passed": isotropic["nonzero_equations"] > 0,
            "detail": (
                f"{isotropic['nonzero_equations']} of 4096 component equations are nonzero, "
                f"with {isotropic['distinct_nonzero_polynomials']} distinct residual polynomials"
            ),
        },
        {
            "name": "finite-temperature ferromagnetic isotropic tensor fails the raw equation",
            "passed": isotropic["physical_interval_solutions"] == [],
            "detail": (
                f"gcd of all nonzero residuals is {isotropic['common_gcd']}; "
                "it has no root in 0<w<1"
            ),
        },
        {
            "name": "explicit unrestricted six-leg GL(2) gauge reaches a TE solution",
            "passed": gauge["satisfies_tetrahedron_equation"] is True,
            "detail": "the displayed invertible gauge maps the rank-two Ising tensor to the GHZ diagonal projector",
        },
        {
            "name": "the exhibited gauge is not an identical delta-bond-preserving gauge",
            "passed": metric_gauge["possible_for_finite_positive_K"] is False,
            "detail": metric_gauge["reason"],
        },
    ]
    artifact = {
        "provenance": {
            "script": "experiments/e11_tetrahedron.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "arithmetic": "exact Python integers and SymPy polynomials over Q[w]",
            "sympy_role": "compute the reduced univariate Groebner basis, factor it, and solve its generator",
        },
        "data": analysis,
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    print("=== bond-dimension-two 3D Ising tensor ===")
    print("general TE system: 64 variables, 4096 quartic equations")
    print(
        "isotropic residual:",
        isotropic["nonzero_equations"],
        "nonzero components; gcd =",
        isotropic["common_gcd"],
    )
    print("physical roots in 0<w<1:", isotropic["physical_interval_solutions"])
    print("unrestricted leg gauge -> GHZ TE solution:", gauge["satisfies_tetrahedron_equation"])
    print("written", OUTPUT.relative_to(ROOT))
    if not all(check["passed"] for check in checks):
        print("FAIL")
        raise AssertionError("one or more tetrahedron checks failed")
    print("PASS")


if __name__ == "__main__":
    main()
