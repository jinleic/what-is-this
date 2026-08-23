"""Produce the exact finite corner-transfer signature certificate.

The producer first constructs the graph-defined d=2 matrix and d=3 tensor from
e206, then applies the logarithm-free characteristic tests from e207.  The
radius-two 2D block is the smallest fixed-corner block where a two-mode
subset-product relation is non-vacuous.  The radius-one 3D octant is the
smallest object with a summed bulk spin and three nontrivial face legs.
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
from itertools import permutations
from math import comb
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from e206_ctm_definition import (  # noqa: E402
    PLUS,
    contract_third_leg,
    corner_matrix_2d_polynomial,
    corner_tensor_3d_r1_scaled,
    evaluate_polynomial_matrix_with_common_scale,
    evaluate_polynomial_with_common_scale,
    orthant_edges,
    permute_face_state,
    radius_one_face_covectors,
    shared_edge_sector_r1,
    square_face_symmetry_maps,
    tensor_nonzero_count,
)
from e207_ctm_spectrum import (  # noqa: E402
    characteristic_polynomial,
    exact_rank,
    integral_fraction_polynomial,
    leading_principal_minors,
    matrix_power_traces,
    pair_product_characteristic_polynomial,
    polynomial_multiply,
    scale_invariant_second_moment,
    squarefree_gcd_descending,
    two_dimensional_symbolic_certificate,
)
from ising.exact_enumeration import dos_bonds  # noqa: E402
from ising.lattices import cubic, square  # noqa: E402


ARTIFACT = ROOT / "results" / "spectral" / "ctm_signature.json"
CPU_BUDGET_SECONDS = 60.0
RSS_CAP_BYTES = 2_000_000_000
Q_NUMERATOR = 1
Q_DENOMINATOR = 2


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class Budget:
    def __init__(self) -> None:
        self.started = time.process_time()
        self.stage = "initialization"

    def used(self) -> float:
        return time.process_time() - self.started

    def check(self, stage: str) -> None:
        self.stage = stage
        used = self.used()
        rss = max_rss_bytes()
        if used >= CPU_BUDGET_SECONDS:
            raise RuntimeError(
                f"process-time budget exceeded at {stage}: "
                f"{used:.3f}>={CPU_BUDGET_SECONDS}"
            )
        if rss >= RSS_CAP_BYTES:
            raise RuntimeError(f"RSS cap exceeded at {stage}: {rss}>={RSS_CAP_BYTES}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_record(value: Fraction) -> dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
    }


def record_check(
    checks: list[dict[str, object]],
    name: str,
    passed: bool,
    detail: str,
    tag: str = "[COMPUTATION]",
) -> None:
    checks.append({"tag": tag, "name": name, "passed": bool(passed), "detail": detail})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def is_symmetric(matrix: Sequence[Sequence[int]]) -> bool:
    return all(
        int(matrix[row][column]) == int(matrix[column][row])
        for row in range(len(matrix))
        for column in range(len(matrix))
    )


def tensor_is_coordinate_covariant(tensor: Sequence[Sequence[Sequence[int]]]) -> bool:
    """Audit all S3 coordinate permutations, including induced face-bit reorderings."""
    for permutation in permutations(range(3)):
        transported: dict[tuple[int, int], tuple[int, int]] = {}
        for axis in range(3):
            for state in range(16):
                transported[(axis, state)] = permute_face_state(
                    3, 1, axis, state, permutation
                )
        for first in range(16):
            for second in range(16):
                for third in range(16):
                    new_states = [0, 0, 0]
                    for old_axis, state in enumerate((first, second, third)):
                        new_axis, new_state = transported[(old_axis, state)]
                        new_states[new_axis] = new_state
                    if tensor[first][second][third] != tensor[new_states[0]][new_states[1]][new_states[2]]:
                        return False
    return True


def covector_is_invariant(covector: Sequence[int]) -> bool:
    spin_flip = all(covector[state] == covector[state ^ 15] for state in range(16))
    square_invariant = all(
        covector[state] == covector[mapping[state]]
        for mapping in square_face_symmetry_maps()
        for state in range(16)
    )
    return spin_flip and square_invariant


def covectors_are_proportional(left: Sequence[int], right: Sequence[int]) -> bool:
    for first in range(len(left)):
        for second in range(first + 1, len(left)):
            if left[first] * right[second] != left[second] * right[first]:
                return False
    return True

def scaled_partition_from_satisfied_dos(lattice: object) -> int:
    bond_count = int(lattice.n_bonds)
    counts = dos_bonds(lattice)
    return sum(
        int(counts[satisfied])
        * Q_NUMERATOR ** (bond_count - satisfied)
        * Q_DENOMINATOR**satisfied
        for satisfied in range(bond_count + 1)
    )



def run_certificate() -> tuple[dict[str, object], list[dict[str, object]], Budget]:
    budget = Budget()
    checks: list[dict[str, object]] = []

    edge_count_2d = len(orthant_edges(2, 2))
    matrix_2d_polynomial = corner_matrix_2d_polynomial(radius=2, fixed_corner_spin=PLUS)
    matrix_2d = evaluate_polynomial_matrix_with_common_scale(
        matrix_2d_polynomial,
        Q_NUMERATOR,
        Q_DENOMINATOR,
        edge_count_2d,
    )
    scale_2d = Q_DENOMINATOR**edge_count_2d
    budget.check("2D finite corner construction")

    symbolic_2d = two_dimensional_symbolic_certificate(matrix_2d_polynomial)
    determinant_root = symbolic_2d["determinant_root_coefficients_ascending"]
    expected_root = [0] * 6
    for power in range(7):
        degree = 6 + 2 * power
        if len(expected_root) <= degree:
            expected_root.extend([0] * (degree + 1 - len(expected_root)))
        expected_root[degree] = (-1) ** power * comb(6, power)
    while len(expected_root) > 1 and expected_root[-1] == 0:
        expected_root.pop()
    scaled_repeated_product = evaluate_polynomial_with_common_scale(
        determinant_root,
        Q_NUMERATOR,
        Q_DENOMINATOR,
        2 * edge_count_2d,
    )
    budget.check("2D symbolic factorization")

    charpoly_2d = characteristic_polynomial(matrix_2d)
    pair_charpoly_2d = pair_product_characteristic_polynomial(matrix_2d)
    char_gcd_2d = integral_fraction_polynomial(squarefree_gcd_descending(charpoly_2d))
    pair_gcd_2d = integral_fraction_polynomial(
        squarefree_gcd_descending(pair_charpoly_2d)
    )
    minors_2d = leading_principal_minors(matrix_2d)
    pair_slots_2d = len(matrix_2d) * (len(matrix_2d) - 1) // 2
    distinct_pairs_2d = pair_slots_2d - (len(pair_gcd_2d) - 1)
    enumerator_2d = scaled_partition_from_satisfied_dos(
        square(3, 3, periodic=False)
    )
    kernel_sum_2d = 2 * sum(sum(row) for row in matrix_2d)

    record_check(
        checks,
        "2D corner uses the open-grid bond count",
        edge_count_2d == 2 * 2 * 3,
        f"computed {edge_count_2d} bonds",
        "[LEMMA]",
    )
    record_check(
        checks,
        "2D fixed-corner kernel is a symmetric 4x4 matrix",
        len(matrix_2d) == 4 and is_symmetric(matrix_2d),
        f"dimension={len(matrix_2d)}",
        "[LEMMA]",
    )
    record_check(
        checks,
        "2D corner sum matches the repository exact enumerator",
        kernel_sum_2d == enumerator_2d,
        f"fixed-corner blocks total={kernel_sum_2d}; DOS total={enumerator_2d}",
        "[COMPUTATION]",
    )
    record_check(
        checks,
        "2D determinant square is derived symbolically",
        polynomial_multiply(determinant_root, determinant_root)
        == symbolic_2d["determinant_coefficients_ascending"],
        "det(A_2(q)) is an exact square in Z[q]",
        "[THEOREM]",
    )
    record_check(
        checks,
        "2D determinant root has the closed form q^6(1-q^2)^6",
        determinant_root == expected_root,
        f"root coefficients={determinant_root}",
        "[THEOREM]",
    )
    record_check(
        checks,
        "2D pair polynomial contains the complementary-product factor twice",
        symbolic_2d["first_remainder_coefficients_ascending"] == [0]
        and symbolic_2d["second_remainder_coefficients_ascending"] == [0],
        "two exact synthetic-division remainders vanish in Z[q]",
        "[THEOREM]",
    )
    record_check(
        checks,
        "2D symbolic factor is not a manufactured third repetition",
        symbolic_2d["third_remainder_coefficients_ascending"] != [0],
        "the third exact synthetic-division remainder is nonzero",
        "[THEOREM]",
    )
    record_check(
        checks,
        "2D q=1/2 spectrum is positive and simple",
        all(value > 0 for value in minors_2d) and char_gcd_2d == [1],
        f"leading minors={minors_2d}; char gcd={char_gcd_2d}",
        "[COMPUTATION]",
    )
    record_check(
        checks,
        "2D q=1/2 pair-product relation is non-vacuous and exact",
        pair_gcd_2d == [1, -scaled_repeated_product]
        and distinct_pairs_2d == pair_slots_2d - 1,
        (
            f"pair gcd={pair_gcd_2d}; repeated product={scaled_repeated_product}; "
            f"distinct pair products={distinct_pairs_2d}/{pair_slots_2d}"
        ),
        "[THEOREM]",
    )

    scale_3d, tensor_3d = corner_tensor_3d_r1_scaled(
        Q_NUMERATOR, Q_DENOMINATOR
    )
    edge_count_3d = len(orthant_edges(3, 1))
    nonzero_tensor_entries = tensor_nonzero_count(tensor_3d)
    enumerator_3d = scaled_partition_from_satisfied_dos(
        cubic(2, 2, 2, periodic=False)
    )
    tensor_sum_3d = sum(
        value for plane in tensor_3d for row in plane for value in row
    )
    budget.check("3D radius-one tensor construction")

    covectors = radius_one_face_covectors()
    free_matrix = contract_third_leg(tensor_3d, covectors["free"])
    aligned_matrix = contract_third_leg(tensor_3d, covectors["aligned"])
    biased_matrix = contract_third_leg(tensor_3d, covectors["biased"])
    free_sector = shared_edge_sector_r1(free_matrix, fixed_spin=PLUS)
    biased_sector = shared_edge_sector_r1(biased_matrix, fixed_spin=PLUS)

    charpoly_3d = characteristic_polynomial(free_sector)
    pair_charpoly_3d = pair_product_characteristic_polynomial(free_sector)
    char_gcd_3d = integral_fraction_polynomial(squarefree_gcd_descending(charpoly_3d))
    pair_gcd_3d = integral_fraction_polynomial(
        squarefree_gcd_descending(pair_charpoly_3d)
    )
    minors_3d = leading_principal_minors(free_sector)
    pair_slots_3d = len(free_sector) * (len(free_sector) - 1) // 2
    distinct_pairs_3d = pair_slots_3d - (len(pair_gcd_3d) - 1)

    free_shape = scale_invariant_second_moment(free_matrix)
    biased_shape = scale_invariant_second_moment(biased_matrix)
    free_sector_shape = scale_invariant_second_moment(free_sector)
    biased_sector_shape = scale_invariant_second_moment(biased_sector)
    free_rank = exact_rank(free_matrix)
    aligned_rank = exact_rank(aligned_matrix)
    biased_rank = exact_rank(biased_matrix)
    budget.check("3D exact spectral and ambiguity calculations")

    tensor_covariant = tensor_is_coordinate_covariant(tensor_3d)
    record_check(
        checks,
        "3D octant uses the open-cube bond count and common scale",
        edge_count_3d == 3 * 1 * 4 and scale_3d == Q_DENOMINATOR**edge_count_3d,
        f"bonds={edge_count_3d}; scale={scale_3d}",
        "[LEMMA]",
    )
    record_check(
        checks,
        "3D honest analogue is a nontrivial 16x16x16 boundary tensor",
        len(tensor_3d) == 16 and nonzero_tensor_entries == 1 << 7,
        f"nonzero compatible entries={nonzero_tensor_entries}",
        "[LEMMA]",
    )
    record_check(
        checks,
        "3D tensor sum matches the repository exact enumerator",
        tensor_sum_3d == enumerator_3d,
        f"tensor total={tensor_sum_3d}; DOS total={enumerator_3d}",
        "[COMPUTATION]",
    )
    record_check(
        checks,
        "3D boundary tensor is S3 coordinate-covariant",
        tensor_covariant,
        "all six coordinate permutations and induced face-bit maps agree",
        "[LEMMA]",
    )
    record_check(
        checks,
        "declared free-face plus-edge spectralization is symmetric 4x4",
        len(free_sector) == 4 and is_symmetric(free_sector),
        f"dimension={len(free_sector)}",
        "[COMPUTATION]",
    )
    record_check(
        checks,
        "3D free-face plus-edge spectrum is positive and simple",
        all(value > 0 for value in minors_3d) and char_gcd_3d == [1],
        f"leading minors={minors_3d}; char gcd={char_gcd_3d}",
        "[COMPUTATION]",
    )
    record_check(
        checks,
        "3D free-face plus-edge pair polynomial is squarefree",
        pair_gcd_3d == [1] and distinct_pairs_3d == pair_slots_3d,
        f"pair gcd={pair_gcd_3d}; distinct pair products={distinct_pairs_3d}/{pair_slots_3d}",
        "[THEOREM]",
    )
    record_check(
        checks,
        "3D free-face plus-edge block fails the two-mode additive/multiplicative signature",
        charpoly_3d[-1] != 0
        and char_gcd_3d == [1]
        and pair_gcd_3d == [1]
        and distinct_pairs_3d == 6,
        "a nonzero simple four-slot subset-product spectrum would force a repeated complementary product",
        "[THEOREM]",
    )

    free_and_biased_symmetric = covector_is_invariant(
        covectors["free"]
    ) and covector_is_invariant(covectors["biased"])
    both_strictly_positive = all(value > 0 for value in covectors["free"]) and all(
        value > 0 for value in covectors["biased"]
    )
    nonproportional_covectors = not covectors_are_proportional(
        covectors["free"], covectors["biased"]
    )
    record_check(
        checks,
        "two strictly positive invariant third-face covectors exist",
        free_and_biased_symmetric
        and both_strictly_positive
        and nonproportional_covectors,
        "h_free=1 and h_biased=1+1_{all face spins aligned}",
        "[THEOREM]",
    )
    record_check(
        checks,
        "invariant positive contractions have non-scale-equivalent spectra",
        free_shape != biased_shape
        and free_sector_shape != biased_sector_shape
        and free_rank == biased_rank == 16,
        (
            f"full I2: {free_shape} versus {biased_shape}; sector I2: "
            f"{free_sector_shape} versus {biased_sector_shape}; ranks={free_rank},{biased_rank}"
        ),
        "[THEOREM]",
    )
    record_check(
        checks,
        "aligned-face contraction exposes collision-test vacuity",
        aligned_rank < free_rank and aligned_rank == 8,
        f"rank(h_aligned contraction)={aligned_rank} versus rank(h_free)={free_rank}",
        "[COMPUTATION]",
    )

    budget.check("final resource audit")
    record_check(
        checks,
        "process-time and peak-RSS budgets",
        budget.used() < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        (
            f"CPU={budget.used():.3f}s<{CPU_BUDGET_SECONDS}; "
            f"RSS={max_rss_bytes()}<{RSS_CAP_BYTES}"
        ),
        "[COMPUTATION]",
    )

    data: dict[str, object] = {
        "finite_kernel_definition": {
            "tag": "[LEMMA]",
            "statement": (
                "For Q_d(r)={0,...,r}^d with open nearest-neighbour bonds and faces "
                "F_i={x_i=0}, C_{d,r}(q)[alpha_1,...,alpha_d] is zero for incompatible "
                "face assignments and otherwise sums q^(broken bonds) over spins outside "
                "the union of the faces. Every graph bond occurs exactly once."
            ),
            "spin_bit_convention": "bit 1 is spin +1",
            "edge_weight": "1 for equal spins and q=e^(-2K) for unequal spins; local q is repository x",
            "two_dimensional_type": "symmetric matrix after fixing the shared corner spin",
            "three_dimensional_type": "S3-covariant rank-three tensor with overlapping face legs",
        },
        "repository_enumerator_cross_checks": {
            "tag": "[COMPUTATION]",
            "two_dimensional_open_3x3_scaled_partition": enumerator_2d,
            "two_dimensional_sum_of_fixed_corner_blocks": kernel_sum_2d,
            "three_dimensional_open_2x2x2_scaled_partition": enumerator_3d,
            "three_dimensional_boundary_tensor_sum": tensor_sum_3d,
        },
        "two_dimensional_integrable_control": {
            "tag": "[THEOREM]",
            "scope": "radius r=2 fixed-plus block; exact symbolic q and exact point q=1/2",
            "radius": 2,
            "dimension": len(matrix_2d),
            "edge_count": edge_count_2d,
            "q": f"{Q_NUMERATOR}/{Q_DENOMINATOR}",
            "common_matrix_scale": scale_2d,
            "scaled_matrix": matrix_2d,
            "characteristic_polynomial_descending": charpoly_2d,
            "characteristic_gcd_derivative_descending": char_gcd_2d,
            "pair_product_polynomial_descending": pair_charpoly_2d,
            "pair_product_gcd_derivative_descending": pair_gcd_2d,
            "leading_principal_minors": minors_2d,
            "distinct_pair_products": distinct_pairs_2d,
            "pair_product_slots": pair_slots_2d,
            "symbolic": {
                "tag": "[THEOREM]",
                "determinant_root": "q^6(1-q^2)^6",
                **symbolic_2d,
            },
            "verdict": (
                "The simple positive four-slot spectrum is exactly relabellable as "
                "c*{1,x,y,xy}; this is a non-vacuous two-mode additive-in-exponents "
                "signature certified only by products, with no logarithms."
            ),
        },
        "three_dimensional_honest_analogue": {
            "tag": "[LEMMA]",
            "scope": "radius r=1 octant (one summed bulk spin) at exact q=1/2",
            "radius": 1,
            "leg_dimension": 16,
            "edge_count": edge_count_3d,
            "common_tensor_scale": scale_3d,
            "nonzero_compatible_entries": nonzero_tensor_entries,
            "intrinsic_characteristic_polynomial": None,
            "reason": (
                "A rank-three tensor is not an endomorphism. A face covector or another "
                "contraction rule is required before a characteristic polynomial exists."
            ),
        },
        "declared_free_face_spectralization": {
            "tag": "[THEOREM]",
            "definition": (
                "Contract the third face with h_free(gamma)=1 and fix the shared edge of "
                "the two remaining faces to spin +1; the remaining two spins on each face "
                "index a 4x4 matrix."
            ),
            "scaled_matrix": free_sector,
            "characteristic_polynomial_descending": charpoly_3d,
            "characteristic_gcd_derivative_descending": char_gcd_3d,
            "pair_product_polynomial_descending": pair_charpoly_3d,
            "pair_product_gcd_derivative_descending": pair_gcd_3d,
            "leading_principal_minors": minors_3d,
            "distinct_pair_products": distinct_pairs_3d,
            "pair_product_slots": pair_slots_3d,
            "verdict": (
                "This precisely declared finite regularization has six distinct slot-pair "
                "products and therefore is not a nonzero two-mode subset-product spectrum."
            ),
        },
        "definition_level_obstruction": {
            "tag": "[THEOREM]",
            "statement": (
                "The graph-defined 3D corner data, positivity, global spin flip, and square-face "
                "symmetry do not select a unique matrix spectrum. The strictly positive invariant "
                "covectors h_free and h_biased are nonproportional and their contractions have "
                "different scale-invariant second spectral moments."
            ),
            "covectors": {
                "free": covectors["free"],
                "biased": covectors["biased"],
                "aligned_counterexample": covectors["aligned"],
            },
            "full_contraction_ranks": {
                "free": free_rank,
                "biased": biased_rank,
                "aligned_counterexample": aligned_rank,
            },
            "full_scale_invariant_second_moments": {
                "free": fraction_record(free_shape),
                "biased": fraction_record(biased_shape),
            },
            "plus_edge_sector_scale_invariant_second_moments": {
                "free": fraction_record(free_sector_shape),
                "biased": fraction_record(biased_sector_shape),
            },
            "biased_plus_edge_characteristic_polynomial_descending": characteristic_polynomial(
                biased_sector
            ),
            "consequence": (
                "The free-face failure is an exact finite boundary-condition result, not an "
                "intrinsic spectrum of the trilinear octant tensor. Choosing a contraction to "
                "obtain a favorable collision pattern would add unproved boundary data."
            ),
        },
        "conclusion": {
            "tag": "[UNRESOLVED]",
            "standard_3d_ising_solved": False,
            "infinite_volume_ctm_signature_decided": False,
            "finite_result": (
                "The radius-two 2D control has the exact Baxter-style two-mode product relation; "
                "the radius-one 3D free-face plus-edge regularization fails it exactly."
            ),
            "route_limitation": (
                "The honest 3D analogue is tensor-valued and admits symmetry-preserving, positive, "
                "spectrally inequivalent contractions, so this finite computation cannot be "
                "promoted to a canonical 3D CTM spectral no-go without an additional theorem "
                "selecting and controlling a contraction through the thermodynamic limit."
            ),
        },
    }
    return data, checks, budget


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]], budget: Budget
) -> dict[str, object]:
    source_paths = [
        ROOT / "experiments" / "e206_ctm_definition.py",
        ROOT / "experiments" / "e207_ctm_spectrum.py",
        Path(__file__),
    ]
    return {
        "meta": {
            "tag": "[COMPUTATION]",
            "experiment": "e208_ctm_certificate",
            "provenance": "experiments/e208_ctm_certificate.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e208_ctm_certificate.py",
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact Z[q], integer matrices, Fraction Euclidean gcds, exact Bareiss minors, "
                "and exact spin enumeration; no floating point, eigenvalue approximation, or logarithm"
            ),
            "modular_arithmetic_used": False,
            "benchmarks_used_to_select_or_fit": False,
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "observed_process_time_seconds": round(budget.used(), 6),
                "observed_peak_rss_bytes": max_rss_bytes(),
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in source_paths
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks, budget = run_certificate()
    artifact = make_artifact(data, checks, budget)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if all(bool(item["passed"]) for item in checks):
        print("PASS e208 finite CTM signature certificate", flush=True)
        return 0
    print("FAIL e208 finite CTM signature certificate", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
