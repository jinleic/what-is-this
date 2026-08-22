"""Exact formalization checks for Gaussian spin spectra and the physical curve.

This is the lightweight first producer for the layer-group spectral front.  It
records only identities that can be checked symbolically or combinatorially;
the group-theoretic implication directions are stated in the accompanying
proof note and integrated by e190.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _spin_weight_rows(max_modes: int = 8) -> list[dict[str, object]]:
    """Check the full-spin torus weights against Boolean subset exponents."""
    rows: list[dict[str, object]] = []
    for modes in range(1, max_modes + 1):
        signs = list(itertools.product((-1, 1), repeat=modes))
        normalized = {tuple((sign + 1) // 2 for sign in row) for row in signs}
        subsets = set(itertools.product((0, 1), repeat=modes))
        rows.append(
            {
                "modes": modes,
                "spin_weight_count": len(signs),
                "normalized_exponent_count": len(normalized),
                "subset_count": len(subsets),
                "bijection": normalized == subsets and len(signs) == 1 << modes,
            }
        )
    return rows


def run_formalization() -> dict[str, object]:
    checks: list[dict[str, object]] = []
    weight_rows = _spin_weight_rows()
    _check(
        checks,
        "full-spin torus weights biject with Boolean subsets",
        all(bool(row["bijection"]) for row in weight_rows),
        "checked exactly for n=1..8; the proof is epsilon_i -> (epsilon_i+1)/2",
    )

    t = sp.symbols("t", nonzero=True)
    x = (1 + t) / (1 - t)
    y = (1 + t**2) / (2 * t)
    curve = (y - 1) * x**2 - (y + 1)
    inverse_t = (x - 1) / (x + 1)
    curve_zero = sp.cancel(curve) == 0
    inverse_zero = sp.cancel(inverse_t - t) == 0
    _check(
        checks,
        "physical multiplicative coordinates obey the curve equation",
        curve_zero,
        "x=exp(a)=(1+t)/(1-t), y=exp(2b)=(1+t^2)/(2t)",
    )
    _check(
        checks,
        "the physical curve parametrization is birational away from endpoints",
        inverse_zero,
        "t=(x-1)/(x+1)",
    )

    # The curve polynomial is visibly nonzero in Q[x,y]; checking its coefficient
    # support prevents an accidental simplification of the ambient two-variable
    # obstruction into the zero polynomial.
    X, Y = sp.symbols("x y")
    physical_polynomial = sp.Poly((Y - 1) * X**2 - (Y + 1), X, Y)
    support = sorted((monomial, int(coefficient)) for monomial, coefficient in physical_polynomial.terms())
    expected_support = [((0, 0), -1), ((0, 1), -1), ((2, 0), -1), ((2, 1), 1)]
    _check(
        checks,
        "physical curve is a proper ambient hypersurface",
        support == expected_support,
        f"coefficient support={support}",
    )

    data = {
        "tag": "[LEMMA]",
        "scalar_spin_group": (
            "GSpinSpec_n := C^times rho(Spin(2n,C)) acting on the full 2^n-dimensional "
            "spinor/Fock space"
        ),
        "nonzero_subset_product_multiset": (
            "Sigma_n(a,u)={a product_{i in S} u_i : S subseteq [n]}, with a,u_i in C^times, "
            "including slot multiplicity"
        ),
        "implication_directions": [
            {
                "tag": "[LEMMA]",
                "direction": "conjugate into scalar-spin => nonzero subset-product spectrum",
                "hypotheses": "no diagonalizability hypothesis is needed",
            },
            {
                "tag": "[THEOREM]",
                "direction": "nonzero subset-product spectrum + diagonalizable => conjugate into scalar-spin",
                "hypotheses": "both nonzero spectrum and diagonalizability are necessary for this converse",
            },
            {
                "tag": "[THEOREM]",
                "direction": "for real physical V(a,b), the two formulations are equivalent",
                "hypotheses": (
                    "V is similar to exp(aA/2) exp(bB) exp(aA/2), which is Hermitian positive "
                    "definite, hence invertible and diagonalizable"
                ),
            },
        ],
        "torus_weight_conversion": {
            "spin_weights": "c product_i r_i^(epsilon_i), epsilon_i in {-1,+1}",
            "subset_parameters": "a=c product_i r_i^(-1), u_i=r_i^2",
            "converse_choice": "choose r_i^2=u_i and c=a product_i r_i",
            "finite_checks": weight_rows,
        },
        "nonconverse_without_hypotheses": {
            "tag": "[THEOREM]",
            "witness": "the nontrivial Jordan block J_2(1) has the one-mode multiset {1,1}",
            "failure": (
                "it is not conjugate to the scalar Spin(2,C) torus because that torus is "
                "diagonalizable; zero subset factors likewise cannot represent a group element"
            ),
        },
        "physical_curve": {
            "tag": "[LEMMA]",
            "parameter": "t=tanh(a/2)",
            "x": "exp(a)=(1+t)/(1-t)",
            "y": "exp(2b)=(1+t^2)/(2t)",
            "equation": "F_phys(x,y)=(y-1)x^2-(y+1)=0",
            "proper_hypersurface_support": [
                {"x_degree": monomial[0], "y_degree": monomial[1], "coefficient": coefficient}
                for monomial, coefficient in support
            ],
        },
    }
    return {"data": data, "checks": checks}


def main() -> int:
    result = run_formalization()
    for item in result["checks"]:
        print(f"[{'PASS' if item['passed'] else 'FAIL'}] {item['name']}: {item['detail']}")
    if all(bool(item["passed"]) for item in result["checks"]):
        print("PASS e188 Gaussian-group formalization")
        return 0
    print("FAIL e188 Gaussian-group formalization")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
