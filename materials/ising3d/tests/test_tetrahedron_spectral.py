"""Standalone exact checks for the spectral tetrahedron investigation."""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.lattices import cubic


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "e35_tetrahedron_spectral.py"


def _load_experiment():
    spec = importlib.util.spec_from_file_location("e35_tetrahedron_spectral", EXPERIMENT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {EXPERIMENT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    experiment = _load_experiment()

    expected = {
        (2, 2, 2): (
            36450,
            [1, 0, 0, 0, 6, 0, 16, 0, 9, 0, 0, 0, 0],
        ),
        (2, 2, 3): (
            16394562,
            [1, 0, 0, 0, 11, 0, 36, 0, 99, 0, 184, 0, 157, 0, 20, 0, 4, 0, 0, 0, 0],
        ),
        (2, 3, 3): (
            246853161090,
            [
                1, 0, 0, 0, 20, 0, 78, 0, 402, 0, 1640, 0, 4909, 0,
                11636, 0, 17685, 0, 16840, 0, 9374, 0, 2542, 0, 376,
                0, 32, 0, 1, 0, 0, 0, 0, 0,
            ],
        ),
    }
    controls = experiment.finite_lattice_controls()
    assert len(controls) == 3
    for row in controls:
        shape = tuple(row["shape"])
        integer, polynomial = expected[shape]
        assert row["common_scaled_partition_integer"] == integer
        assert row["even_subgraph_polynomial"] == polynomial
        # Independent repository spin enumerator is the control; this assertion
        # does not reuse the tensor contraction implementation.
        assert even_subgraph_polynomial(cubic(*shape, periodic=False)) == polynomial
        assert row["coefficientwise_match"] is True
    print("tensor/enumerator finite boxes: PASS", [row["common_scaled_partition_integer"] for row in controls])

    spectral = experiment.spectral_tetrahedron_analysis()
    reduction = spectral["linear_elimination"]
    assert reduction["nonzero_component_equations"] == 1792
    assert reduction["distinct_after_monomial_and_square_reduction"] == 28
    assert reduction["exact_rational_linear_rank"] == 7
    assert spectral["nondegenerate_saturation"]["groebner_basis"] == ["1"]
    assert spectral["nondegenerate_saturation"]["unit_ideal"] is True
    assert spectral["physical_isotropic_slice"]["raw_common_gcd_in_w"] == (
        "w**2*(w - 1)**3*(w + 1)**3*(w**2 + 1)**2"
    )
    assert spectral["physical_isotropic_slice"]["physical_solutions"] == []
    assert spectral["degenerate_manifold"]["verified_by_substitution"] is True
    print("six-rapidity TE elimination: PASS [1] on the nonexceptional open set")

    escapes = experiment.escape_analysis()
    assert escapes["gauge_transformation"]["decision"] == "fails"
    auxiliary = escapes["auxiliary_state_extension"]
    assert auxiliary["decision"] == "fails"
    assert auxiliary["satisfies_tetrahedron_equation"] is False
    assert auxiliary["witness"] == [20, 0, -675, "-675/4096"]
    assert "undecided" in auxiliary["scope_warning"]
    assert escapes["controlled_limit"]["decision"] == "fails for a regular finite physical limit"
    assert "undecided" in escapes["controlled_limit"]["scope_warning"]
    print("gauge/direct-sum/regular-limit escape checks: PASS")

    transfer = experiment.commuting_transfer_analysis()
    assert [row["shape"] for row in transfer["exact_pair_checks"]] == [[2, 2], [2, 3]]
    assert all(row["commutes"] is False for row in transfer["exact_pair_checks"])
    assert [row["first_nonzero_component"]["value"] for row in transfer["exact_pair_checks"]] == [
        "-1213745/30233088",
        "-127393295/2176782336",
    ]
    x = sp.Symbol("x")
    assert sp.factor(sp.sympify(transfer["two_by_two_tangent_certificate"]["determinant"]) - (
        4 * x**6 * (x - 1) ** 4 * (x + 1) ** 4 * (x**2 + 1) ** 2 * (x**4 + 3)
    )) == 0
    # Direct exact matrix smoke check at an additional untuned rational pair.
    first = experiment.layer_transfer_matrix(2, 2, Fraction(2, 5), Fraction(3, 7))
    second = experiment.layer_transfer_matrix(2, 2, Fraction(1, 4), Fraction(5, 9))
    assert experiment._matrix_commutator_nonzero(first, second) is not None
    print("periodic-layer commuting-family obstruction: PASS")
    print("PASS")


if __name__ == "__main__":
    main()
