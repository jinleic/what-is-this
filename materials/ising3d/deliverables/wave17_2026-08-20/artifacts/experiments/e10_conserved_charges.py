"""Local conserved-charge and transfer-matrix commutant tests for Ising layers.

All local-charge ranks are exact over two prime fields.  Small transfer commutants use exact
integer characteristic polynomials.  The largest layers use the same exact integer matrices but
resolve eigenvalue multiplicities with symmetry-blocked mpmath diagonalization at 50 and 80
working decimal digits.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ising.tensor_networks import (
    exact_transfer_commutant,
    scan_local_conserved_charges,
    stable_transfer_commutant,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "conserved_charges.json"
PRIMES = (2_147_483_647, 2_147_483_629)
COUPLINGS = ((1, 2), (1, 3))


def _exact_commutant_rows() -> list[dict]:
    rows = []
    for shape in ((4,), (6,), (2, 2), (2, 3)):
        for numerator, denominator in COUPLINGS:
            rows.append(
                exact_transfer_commutant(
                    shape,
                    x_numerator=numerator,
                    x_denominator=denominator,
                )
            )
    return rows


def _multiprecision_commutant_rows() -> list[dict]:
    rows = []
    for shape in ((8,), (3, 3)):
        for numerator, denominator in COUPLINGS:
            rows.append(
                stable_transfer_commutant(
                    shape,
                    x_numerator=numerator,
                    x_denominator=denominator,
                    precisions=(50, 80),
                    tolerance_exponents=(30, 45),
                )
            )
    return rows


def main() -> None:
    chain = scan_local_conserved_charges(
        dimension=1,
        max_support=5,
        linear_size=11,
        field_numerator=2,
        primes=PRIMES,
    )
    grid = scan_local_conserved_charges(
        dimension=2,
        max_support=4,
        linear_size=9,
        field_numerator=2,
        primes=PRIMES,
    )
    chain_g3 = scan_local_conserved_charges(
        dimension=1,
        max_support=5,
        linear_size=11,
        field_numerator=3,
        primes=PRIMES,
    )
    grid_g3 = scan_local_conserved_charges(
        dimension=2,
        max_support=4,
        linear_size=9,
        field_numerator=3,
        primes=PRIMES,
    )
    exact_commutants = _exact_commutant_rows()
    multiprecision_commutants = _multiprecision_commutant_rows()

    chain_dimensions = [row["kernel_dimension"] for row in chain["rows"]]
    grid_dimensions = [row["kernel_dimension"] for row in grid["rows"]]
    checks = [
        {
            "name": "two prime fields give identical local-charge nullities",
            "passed": all(
                len(set(row["kernel_dimension_by_prime"])) == 1
                for scan in (chain, grid, chain_g3, grid_g3)
                for row in scan["rows"]
            ),
            "detail": f"prime fields {PRIMES}",
        },
        {
            "name": "1D free-fermion control reproduces a growing local-charge tower",
            "passed": chain_dimensions == [3, 5, 7, 9],
            "detail": (
                "total kernel dimensions at r=2,3,4,5 are "
                f"{chain_dimensions}; subtracting span{{I,H}} gives "
                f"{[row['nontrivial_dimension'] for row in chain['rows']]}"
            ),
        },
        {
            "name": "2D square-layer search finds only the defined trivial local subspace",
            "passed": grid_dimensions == [2, 2, 2],
            "detail": (
                "kernel dimensions at r=2,3,4 are "
                f"{grid_dimensions} in ansatz spaces of dimensions "
                f"{[row['ansatz_dimension'] for row in grid['rows']]}; since rank_Fp <= rank_Q "
                "and I,H are rational kernel vectors, this certifies kernel_Q = span{I,H}"
            ),
        },
        {
            "name": "local-charge dimensions are stable at two generic rational fields",
            "passed": (
                [row["kernel_dimension"] for row in chain_g3["rows"]] == chain_dimensions
                and [row["kernel_dimension"] for row in grid_g3["rows"]] == grid_dimensions
            ),
            "detail": "the same dimensions were obtained at g=2 and g=3",
        },
        {
            "name": "exact transfer commutants are stable across two generic rational couplings",
            "passed": all(
                len({row["commutant_dimension"] for row in exact_commutants if tuple(row["shape"]) == shape}) == 1
                for shape in ((4,), (6,), (2, 2), (2, 3))
            ),
            "detail": "x=e^(-2K) was 1/2 and 1/3; all matrices were integer-scaled exactly",
        },
        {
            "name": "large-layer multiplicities are stable across 50 and 80 digits",
            "passed": all(row["precision_stability_passed"] for row in multiprecision_commutants),
            "detail": "relative clustering tolerances were 1e-30 and 1e-45 respectively",
        },
        {
            "name": "open 1D n=8 transfer matrix has a simple generic spectrum",
            "passed": all(
                row["commutant_dimension"] == row["matrix_dimension"] == 256
                for row in multiprecision_commutants
                if row["shape"] == [8]
            ),
            "detail": (
                "local commuting charges still exist: for a finite simple-spectrum matrix every "
                "commuting matrix is a polynomial in V, so the unrestricted commutant is not a locality diagnostic"
            ),
        },
    ]

    artifact = {
        "provenance": {
            "script": "experiments/e10_conserved_charges.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "local_charge_arithmetic": "exact sparse Gaussian elimination over two listed prime fields",
            "small_commutant_arithmetic": "exact integer matrices and Q[lambda] square-free decomposition",
            "large_commutant_arithmetic": (
                "exact integer matrix construction; mpmath symmetry-sector eigenspectra at 50 and 80 decimal digits"
            ),
            "rank_tolerances": ["1e-30 at 50 dps", "1e-45 at 80 dps"],
        },
        "data": {
            "local_conserved_charges": {
                "chain": chain,
                "square_layer": grid,
                "square_layer_rational_certificate": (
                    "For an integer matrix rank_Fp <= rank_Q.  The finite-field nullity is 2 and "
                    "the explicit rational vectors I,H give rational nullity at least 2; therefore "
                    "the rational kernel is exactly span_Q{I,H} at every reported r."
                ),
            },
            "generic_field_stability": {
                "chain_g3": chain_g3,
                "square_layer_g3": grid_g3,
            },
            "transfer_matrix_commutants": {
                "definition": "dimension of the full complex matrix algebra {M : M V = V M}",
                "dimension_formula": (
                    "V is diagonalizable, hence dim Comm(V)=sum over distinct eigenvalues of multiplicity squared"
                ),
                "exact_rows": exact_commutants,
                "multiprecision_rows": multiprecision_commutants,
                "interpretation_warning": (
                    "The unrestricted finite-matrix commutant always contains all polynomials in V "
                    "and has dimension at least 2^n.  It does not distinguish local conserved charges."
                ),
            },
        },
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    print("=== local conserved charges ===")
    print(
        "1D chain (r, ansatz, kernel, nontrivial):",
        [
            (
                row["max_support"],
                row["ansatz_dimension"],
                row["kernel_dimension"],
                row["nontrivial_dimension"],
            )
            for row in chain["rows"]
        ],
    )
    print(
        "2D square layer (r, ansatz, kernel, nontrivial):",
        [
            (
                row["max_support"],
                row["ansatz_dimension"],
                row["kernel_dimension"],
                row["nontrivial_dimension"],
            )
            for row in grid["rows"]
        ],
    )
    print("=== full transfer-matrix commutants ===")
    for row in exact_commutants + multiprecision_commutants:
        print(
            f"shape={tuple(row['shape'])} x={row['x']} D={row['matrix_dimension']} "
            f"dim Comm(V)={row['commutant_dimension']} method={row['method']}"
        )
    print("written", OUTPUT.relative_to(ROOT))
    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        print("FAIL", failed)
        raise AssertionError(f"checks failed: {failed}")
    print("PASS")


if __name__ == "__main__":
    main()
