#!/usr/bin/env python3
"""Behavioral and adversarial checks for Gate B definitions and product proof."""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path
from random import Random
import sys
import unittest

import sympy
from flint import arb

HERE = Path(__file__).resolve().parent
UC_ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(UC_ROOT))

import search_n8_block_symmetric as n8
import verify_gate_b_dyadic as dyadic
import verify_gate_b_n6_dyadic as n6_dyadic
import verify_gate_b as verifier
import shapley_n7_block_symmetric as n7
from shapley_n6_shared_bellman import one_sided_costs
from shapley_entropy import ALPHA as PROJECT_ALPHA
from shapley_global_coupling import local_order_parts


CANDIDATE = HERE / "candidates" / "n7_block_extremizer.json"
N8_RESULT = HERE / "candidates" / "n8_k1_census.json"
AUTHORITATIVE_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_arb_v4.json"
)
INDEPENDENT_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_dyadic_v1.json"
)
SECOND_BASE_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_n6_dyadic_v1.json"
)
TENSORIZATION_AUDIT = HERE / "candidates" / "tensorization_audit.json"
TENSORIZATION_CHECKPOINT = (
    HERE / "experiments" / "tensorization_audit_checkpoint.jsonl"
)
POWER_AUDIT = HERE / "candidates" / "power_admissibility_audit.json"
SQUARE_AUDIT = HERE / "candidates" / "square_direct_audit.json"


def average_certified_cost(family: tuple[int, ...], dimension: int) -> arb:
    orders = tuple(permutations(range(dimension)))
    roots, _states, _ambiguous = verifier.certified_one_sided_costs(
        family, dimension, orders
    )
    return sum(roots, verifier.ZERO) / len(roots)


def cartesian_product(
    left: tuple[int, ...], left_dimension: int, right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        left_row | (right_row << left_dimension)
        for left_row in left
        for right_row in right
    )


def direct_defect(family: tuple[int, ...]) -> Fraction:
    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    return Fraction(missing, len(family) ** 2)


class ExactDefinitionTests(unittest.TestCase):
    def test_reimer_threshold_is_integer_exact(self) -> None:
        expected = {3: 3, 5: 6, 25: 59, 45: 124, 70: 215, 128: 448}
        for size, threshold in expected.items():
            self.assertEqual(verifier.exact_reimer_threshold(size), threshold)
            self.assertEqual(n8.exact_reimer_threshold(size), threshold)
            if size <= n7.SUBSET_COUNT:
                self.assertEqual(n7.exact_reimer_threshold(size), threshold)
            self.assertGreaterEqual(1 << (2 * threshold), size**size)
            if threshold:
                self.assertLess(1 << (2 * (threshold - 1)), size**size)

    def test_base_constraints_and_normalization(self) -> None:
        data = json.loads(CANDIDATE.read_text())
        base = verifier.exact_family_checks(data)
        self.assertEqual(base["defect"], Fraction(64, 81))
        self.assertEqual(base["counts"], (18,) * 7)
        self.assertEqual(base["incidence"], 126)
        self.assertEqual(base["reimer_threshold"], 124)
        self.assertTrue(base["active"])
        self.assertTrue(base["separating"])
        self.assertTrue(base["cell_union_verified"])

    def test_alpha_and_power_certificate_are_consistent(self) -> None:
        candidate = json.loads(CANDIDATE.read_text())
        alpha = verifier.fraction_from_json(candidate["alpha"])
        self.assertEqual(alpha, Fraction(356069, 10_000_000))
        self.assertEqual(Fraction(str(PROJECT_ALPHA)), alpha)

        certificate = json.loads(AUTHORITATIVE_CERTIFICATE.read_text())
        self.assertEqual(certificate["verdict"], "PROVED_GATE_B_UNBOUNDED")
        self.assertEqual(certificate["schema_version"], 2)
        self.assertTrue(certificate["exact_base"]["cell_union_verified"])
        self.assertTrue(
            certificate["exact_base"]["ordered_pairs_with_replacement"]
        )
        self.assertEqual(certificate["alpha"]["numerator"], alpha.numerator)
        self.assertEqual(certificate["alpha"]["denominator"], alpha.denominator)
        self.assertTrue(certificate["alpha"]["dimension_independent"])
        powers = certificate["infinite_construction"]["checked_power_arithmetic"]
        self.assertEqual([entry["power"] for entry in powers], list(range(1, 9)))
        self.assertTrue(
            all(entry["cap_holds_with_equality"] for entry in powers)
        )
        self.assertTrue(all(entry["reimer_holds_strictly"] for entry in powers))
        self.assertTrue(
            all("certified_ratio_lower_bound" in entry for entry in powers)
        )
        repair_ratio_lower = arb(
            certificate["arb"]["repair_ratio_lower_certificate"]["lower"]
        )
        self.assertGreater(repair_ratio_lower, arb("0.0362"))
        self.assertFalse(
            certificate["infinite_construction"][
                "zero_defect_extended_value_convention_used"
            ]
        )

        checks = certificate["tensorization_checks"][
            "actual_square_fixed_orders"
        ]
        for check in checks.values():
            residual = arb(check["residual"]["ball"])
            self.assertTrue(residual.contains(0))
        square_a_plus = arb(
            certificate["derived_from_proved_identities"]["square_a_plus"]["ball"]
        )
        self.assertLess(square_a_plus.upper(), arb("-0.057"))

    def test_direct_square_matches_symbolic_formula(self) -> None:
        data = json.loads(CANDIDATE.read_text())
        base = verifier.exact_family_checks(data)
        square = verifier.direct_square_checks(base)
        self.assertEqual(square["family_size"], 2025)
        self.assertEqual(square["coordinate_count"], 810)
        self.assertEqual(square["closure_defect"], Fraction(6272, 6561))
        self.assertEqual(square["missing_ordered_join_pairs"], 3_920_000)

    def test_ordered_pairs_with_replacement_convention(self) -> None:
        family = (0, 1, 2)
        self.assertEqual(direct_defect(family), Fraction(2, 9))
        self.assertTrue(all((row | row) in family for row in family))


class IndependentDyadicCertificateTests(unittest.TestCase):
    def test_elementary_function_enclosures(self) -> None:
        half_entropy = dyadic.binary_entropy(Fraction(1, 2))
        self.assertTrue(half_entropy.contains_fraction(Fraction(1)))
        self.assertTrue(
            dyadic.natural_log(dyadic.Interval.from_int(2)).contains_interval(
                dyadic.LN2
            )
        )
        self.assertTrue(
            dyadic.exponential(dyadic.Interval.zero()).contains_fraction(
                Fraction(1)
            )
        )

    def test_dyadic_enclosures_dominate_arb_controls(self) -> None:
        def lower(interval: dyadic.Interval) -> arb:
            return verifier.aq(Fraction(interval.lo, dyadic.SCALE))

        def upper(interval: dyadic.Interval) -> arb:
            return verifier.aq(Fraction(interval.hi, dyadic.SCALE))

        for value in (
            Fraction(1, 45),
            Fraction(1, 3),
            Fraction(1, 2),
            Fraction(2),
            Fraction(45),
        ):
            dyadic_log = dyadic.natural_log(dyadic.Interval.from_fraction(value))
            arb_log = verifier.aq(value).log()
            self.assertLessEqual(lower(dyadic_log), arb_log.lower())
            self.assertGreaterEqual(upper(dyadic_log), arb_log.upper())

        for value in (
            Fraction(-5),
            Fraction(-1),
            Fraction(0),
            Fraction(1, 3),
            Fraction(5),
        ):
            dyadic_exp = dyadic.exponential(dyadic.Interval.from_fraction(value))
            arb_exp = verifier.aq(value).exp()
            self.assertLessEqual(lower(dyadic_exp), arb_exp.lower())
            self.assertGreaterEqual(upper(dyadic_exp), arb_exp.upper())

        family = dyadic.reconstruct_family()
        orders, _automorphisms = dyadic.order_representatives(family)
        actual_roots, _states, _ambiguous = verifier.certified_one_sided_costs(
            family, dyadic.DIMENSION, orders
        )
        for order, actual in zip(orders, actual_roots):
            relaxed, _states = dyadic.relaxed_bellman_cost(family, order)
            self.assertGreater(upper(relaxed), actual.upper())

        dyadic_iid = dyadic.certified_shapley_iid(family)
        arb_iid = verifier.certified_shapley_iid(family, dyadic.DIMENSION)
        self.assertLess(lower(dyadic_iid), arb_iid.lower())
        self.assertGreater(upper(dyadic_iid), arb_iid.upper())

    def test_order_orbit_reduction_is_value_invariant(self) -> None:
        family = dyadic.reconstruct_family()
        representatives, automorphism_count = dyadic.order_representatives(family)
        self.assertEqual(automorphism_count, 240)
        self.assertEqual(len(representatives), 21)

        support = frozenset(family)
        automorphisms = [
            permutation
            for permutation in permutations(range(dyadic.DIMENSION))
            if frozenset(
                verifier.permute_subset(row, permutation) for row in family
            )
            == support
        ]
        self.assertEqual(len(automorphisms), automorphism_count)

        random = Random(20260826)
        for representative in random.sample(representatives, 3):
            baseline, _states = dyadic.relaxed_bellman_cost(
                family, representative
            )
            for automorphism in random.sample(automorphisms, 3):
                relabeled = tuple(
                    automorphism[coordinate] for coordinate in representative
                )
                value, _states = dyadic.relaxed_bellman_cost(family, relabeled)
                self.assertLess(
                    abs(value.lo - baseline.lo) / dyadic.SCALE, 1e-40
                )
                self.assertLess(
                    abs(value.hi - baseline.hi) / dyadic.SCALE, 1e-40
                )

    def test_independent_relaxation_certificate(self) -> None:
        result = dyadic.verify()
        self.assertEqual(result["verdict"], "PROVED_GATE_B_UNBOUNDED")
        self.assertEqual(result["arithmetic"], "pure_python_dyadic_intervals")
        self.assertEqual(result["exact_base"]["closure_defect"], "64/81")
        self.assertEqual(result["exact_base"]["join_success_probability"], "17/81")
        self.assertTrue(result["relaxation"]["clamp_free_upper_bound"])
        self.assertEqual(result["rational_consequence"]["a_plus_upper_lt"], "-1/40")
        self.assertEqual(
            result["rational_consequence"]["base_ratio_lower_gt"], "81/2560"
        )
        self.assertEqual(
            json.loads(INDEPENDENT_CERTIFICATE.read_text()),
            result,
        )

    def test_second_base_all_order_relaxation_certificate(self) -> None:
        result = n6_dyadic.verify()
        self.assertEqual(
            result["verdict"],
            "PROVED_GATE_B_UNBOUNDED_SECOND_BASE",
        )
        self.assertEqual(result["exact_base"]["closure_defect"], "444/625")
        self.assertEqual(
            result["exact_base"]["join_success_probability"],
            "181/625",
        )
        self.assertFalse(
            result["all_order_evaluation"]["symmetry_quotient_used"]
        )
        self.assertEqual(
            result["all_order_evaluation"]["evaluated_order_count"],
            720,
        )
        self.assertEqual(
            result["rational_consequence"]["a_plus_upper_lt"],
            "-1/80",
        )
        self.assertEqual(
            json.loads(SECOND_BASE_CERTIFICATE.read_text()),
            result,
        )


class ProductIdentityTests(unittest.TestCase):
    def assert_product_additive(
        self,
        left: tuple[int, ...],
        left_dimension: int,
        right: tuple[int, ...],
        right_dimension: int,
    ) -> None:
        product = cartesian_product(left, left_dimension, right)
        product_dimension = left_dimension + right_dimension

        q_left = verifier.certified_shapley_iid(left, left_dimension)
        q_right = verifier.certified_shapley_iid(right, right_dimension)
        q_product = verifier.certified_shapley_iid(product, product_dimension)
        self.assertTrue((q_product - q_left - q_right).contains(0))

        c_left = average_certified_cost(left, left_dimension)
        c_right = average_certified_cost(right, right_dimension)
        c_product = average_certified_cost(product, product_dimension)
        self.assertTrue((c_product - c_left - c_right).contains(0))

        epsilon_expected = 1 - (1 - direct_defect(left)) * (
            1 - direct_defect(right)
        )
        self.assertEqual(direct_defect(product), epsilon_expected)

    def test_hand_derived_product(self) -> None:
        self.assert_product_additive((0, 1, 2), 2, (0, 1), 1)

    def test_seeded_random_products(self) -> None:
        random = Random(20260825)
        nonempty_families = [
            tuple(row for row in range(4) if mask >> row & 1)
            for mask in range(1, 1 << 4)
        ]
        for _ in range(4):
            left = random.choice(nonempty_families)
            right = random.choice(((0,), (1,), (0, 1)))
            self.assert_product_additive(left, 2, right, 1)

    def test_exhaustive_n2_products_attack_tensorization(self) -> None:
        families = tuple(
            tuple(row for row in range(4) if mask >> row & 1)
            for mask in range(1, 1 << 4)
        )
        factor_orders = tuple(permutations(range(2)))
        product_orders = tuple(permutations(range(4)))
        factor_costs = {}
        factor_iid = {}
        for family in families:
            roots, _states = one_sided_costs(
                family, factor_orders, dimension=2
            )
            factor_costs[family] = dict(zip(factor_orders, roots))
            factor_iid[family] = {
                order: local_order_parts(family, order)[0]
                for order in factor_orders
            }

        for left in families:
            for right in families:
                product = cartesian_product(left, 2, right)
                roots, _states = one_sided_costs(
                    product, product_orders, dimension=4
                )
                for global_order, product_cost in zip(product_orders, roots):
                    left_order = tuple(
                        coordinate
                        for coordinate in global_order
                        if coordinate < 2
                    )
                    right_order = tuple(
                        coordinate - 2
                        for coordinate in global_order
                        if coordinate >= 2
                    )
                    self.assertAlmostEqual(
                        product_cost,
                        factor_costs[left][left_order]
                        + factor_costs[right][right_order],
                        places=12,
                    )
                    self.assertAlmostEqual(
                        local_order_parts(product, global_order)[0],
                        factor_iid[left][left_order]
                        + factor_iid[right][right_order],
                        places=12,
                    )
                expected_defect = 1 - (1 - direct_defect(left)) * (
                    1 - direct_defect(right)
                )
                self.assertEqual(direct_defect(product), expected_defect)

    def test_shared_float_evaluator_matches_standalone_arb(self) -> None:
        family = (0, 1, 2, 5, 7)
        orders = tuple(permutations(range(3)))
        float_roots, _states = one_sided_costs(family, orders, dimension=3)
        arb_roots, _states, _ambiguous = verifier.certified_one_sided_costs(
            family, 3, orders
        )
        for float_value, arb_value in zip(float_roots, arb_roots):
            self.assertLess(abs(float_value - float(arb_value)), 2e-14)


class SymmetryAndSearchTests(unittest.TestCase):
    def test_n8_block_order_representatives_partition_all_orders(self) -> None:
        representatives = n8.block_order_representatives(1)
        internal = tuple(
            (0,) + permutation
            for permutation in permutations(range(1, n8.DIMENSION))
        )
        covered: set[tuple[int, ...]] = set()
        for representative in representatives:
            orbit = {
                tuple(automorphism[coordinate] for coordinate in representative)
                for automorphism in internal
            }
            self.assertEqual(len(orbit), len(internal))
            self.assertTrue(covered.isdisjoint(orbit))
            covered.update(orbit)
        self.assertEqual(covered, set(permutations(range(n8.DIMENSION))))

    def test_n8_canonicalization_matches_full_orbit(self) -> None:
        families, cell_count, raw_accepted = n8.enumerate_canonical(1)
        self.assertEqual(cell_count, 16)
        self.assertEqual(raw_accepted, 136)
        family = families[37]
        block_canonical = n8.canonical_block_family(family, 1)
        full_canonical = min(
            n8.family_mask(
                tuple(n8.permute_subset(row, permutation) for row in family)
            )
            for permutation in permutations(range(n8.DIMENSION))
        )
        self.assertEqual(block_canonical, full_canonical)

    def test_completed_n8_result_has_consistent_scope(self) -> None:
        result = json.loads(N8_RESULT.read_text())
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["raw_mask_count"], 65_535)
        self.assertEqual(result["raw_accepted_count"], 136)
        self.assertEqual(result["canonical_family_count"], 136)
        self.assertEqual(result["normalized_family_count"], 130)
        self.assertEqual(result["negative_a_plus_count"], 2)
        self.assertTrue(result["largest_repair_ratio"]["normalized"])
        self.assertLess(result["largest_repair_ratio"]["repair_ratio"], 0.035)


class SymbolicBellmanAlgebraTests(unittest.TestCase):
    """Machine-check the algebra that both certificates and the proof reuse."""

    def test_step_reduces_to_affine_plus_entropy(self) -> None:
        s, p, r = sympy.symbols("s p r", real=True)
        v00, v10, v01, v11, vg = sympy.symbols(
            "v00 v10 v01 v11 vg", real=True
        )
        weights = ((1 - s, v00), (s - r, v10), (s - p, v01), (p + r - s, v11))
        self.assertEqual(
            sympy.simplify(sum(weight for weight, _value in weights)), 1
        )

        entropy = -(s * sympy.log(s) + (1 - s) * sympy.log(1 - s)) / sympy.log(2)
        objective = entropy + sum(weight * value for weight, value in weights)
        slope = -v00 + v10 + v01 - v11
        constant = v00 - r * v10 - p * v01 + (p + r) * v11
        self.assertEqual(
            sympy.simplify(objective - (constant + s * slope + entropy)), 0
        )

        shifted_weights = (
            (1 - s, v00 + vg),
            (s - r, v10 + vg),
            (s - p, v01 + vg),
            (p + r - s, v11 + vg),
        )
        shifted = entropy + sum(
            weight * value for weight, value in shifted_weights
        )
        shifted_slope = -(v00 + vg) + (v10 + vg) + (v01 + vg) - (v11 + vg)
        self.assertEqual(sympy.simplify(shifted - objective - vg), 0)
        self.assertEqual(sympy.simplify(shifted_slope - slope), 0)

    def test_unconstrained_maximum_is_softplus(self) -> None:
        s = sympy.symbols("s", positive=True)
        slope = sympy.symbols("slope", real=True)
        entropy = -(s * sympy.log(s) + (1 - s) * sympy.log(1 - s)) / sympy.log(2)
        objective = s * slope + entropy

        roots = sympy.solve(sympy.diff(objective, s), s)
        self.assertEqual(len(roots), 1)
        logistic = 1 / (1 + 2 ** (-slope))
        self.assertEqual(sympy.simplify(roots[0] - logistic), 0)

        maximum = sympy.simplify(
            objective.subs(s, logistic)
            - sympy.log(1 + 2**slope) / sympy.log(2)
        )
        self.assertEqual(sympy.simplify(sympy.expand_log(maximum, force=True)), 0)

        second = sympy.simplify(sympy.diff(objective, s, 2))
        self.assertEqual(second, 1 / (s * (s - 1) * sympy.log(2)))
        self.assertLess(second.subs(s, sympy.Rational(1, 3)), 0)



class AuditArtifactTests(unittest.TestCase):
    def test_tensorization_audit_found_no_counterexample(self) -> None:
        summary = json.loads(TENSORIZATION_AUDIT.read_text())
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["counterexample_count"], 0)
        self.assertEqual(summary["exact_defect_identity_failures"], 0)
        self.assertGreaterEqual(summary["case_count"], 292)
        self.assertGreater(summary["unequal_dimension_case_count"], 100)
        self.assertLess(summary["worst_fixed_order_bellman_gap"], 1e-12)
        self.assertLess(summary["worst_fixed_order_iid_gap"], 1e-12)
        self.assertLess(summary["worst_order_averaged_iid_gap"], 1e-12)

        records = [
            json.loads(line)
            for line in TENSORIZATION_CHECKPOINT.read_text().splitlines()
            if line.strip()
        ]
        self.assertEqual(len(records), summary["case_count"])
        self.assertEqual(
            len({record["case_id"] for record in records}), len(records)
        )
        self.assertTrue(all(record["defect_identity_exact"] for record in records))

    def test_direct_square_audit_confirms_product_identities(self) -> None:
        summary = json.loads(SQUARE_AUDIT.read_text())
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["checked_order_count"], 5)
        self.assertTrue(summary["all_within_tolerance"])
        self.assertLess(summary["worst_bellman_gap"], 1e-12)
        self.assertLess(summary["worst_iid_gap"], 1e-12)

        exact = summary["exact_square"]
        self.assertEqual(exact["family_size"], 2025)
        self.assertEqual(exact["dimension"], 14)
        self.assertEqual(exact["coordinate_counts"], [810] * 14)
        self.assertEqual(exact["cap_bound"], 810)
        self.assertTrue(exact["cap_holds"])
        self.assertTrue(exact["reimer_holds_strictly"])
        self.assertEqual(exact["missing_ordered_join_pairs"], 3_920_000)
        self.assertEqual(exact["closure_defect"], "6272/6561")
        self.assertTrue(exact["closure_defect_matches_product"])
        self.assertTrue(exact["active"] and exact["separating"])
        for order in summary["orders"].values():
            self.assertAlmostEqual(
                order["direct_square_bellman"],
                order["base_block_bellman_sum"],
                places=12,
            )

    def test_constructed_powers_are_admissible_and_match_defects(self) -> None:
        summary = json.loads(POWER_AUDIT.read_text())
        self.assertEqual(summary["status"], "complete")
        self.assertTrue(summary["all_admissible"])
        self.assertTrue(summary["all_defects_match_product_formula"])
        self.assertTrue(summary["all_counts_uniform"])
        self.assertTrue(summary["all_incidences_match"])

        cube = summary["powers"]["3"]
        self.assertEqual(cube["family_size"], 91_125)
        self.assertEqual(cube["dimension"], 21)
        self.assertEqual(cube["distinct_coordinate_counts"], [36_450])
        self.assertEqual(cube["cap_bound"], 36_450)
        self.assertTrue(cube["cap_holds_with_equality"])
        self.assertEqual(cube["total_incidence"], 765_450)
        self.assertTrue(cube["reimer_holds_strictly"])
        self.assertEqual(cube["brute_force_pair_count"], 8_303_765_625)
        self.assertEqual(cube["missing_ordered_join_pairs"], 8_227_000_000)
        self.assertEqual(cube["closure_defect"], "526528/531441")
        self.assertEqual(
            Fraction(cube["closure_defect"]),
            1 - Fraction(17, 81) ** 3,
        )


class GeneralImplicationTests(unittest.TestCase):
    """The finite base plus exact integer facts must give the all-k statement."""

    def test_cap_is_met_with_exact_integer_equality_for_every_power(self) -> None:
        for power in range(1, 65):
            size = 45**power
            count = 18 * 45 ** (power - 1)
            self.assertEqual(2 * size % 5, 0)
            self.assertEqual(2 * size // 5, count)
            self.assertGreaterEqual(size, 3)

    def test_reimer_strictness_reduces_to_one_integer_witness(self) -> None:
        # 2^(2I) > m^m with m = 45^k and I = 126k45^(k-1) is, after writing both
        # sides as an exponent times k45^(k-1), exactly 2^252 > 45^45.
        self.assertLess(45**5, 2**28)
        self.assertLess(45**45, 2**252)

        for power in range(1, 65):
            size = 45**power
            incidence = 126 * power * 45 ** (power - 1)
            common = power * 45 ** (power - 1)
            self.assertEqual(2 * incidence, 252 * common)
            self.assertEqual(power * size, 45 * common)

        # Direct exact m^m comparisons are only tractable for the small powers.
        for power in (1, 2):
            size = 45**power
            incidence = 126 * power * 45 ** (power - 1)
            self.assertGreater(1 << (2 * incidence), size**size)
            self.assertGreater(incidence, verifier.exact_reimer_threshold(size))

    def test_second_base_power_implication_is_exact(self) -> None:
        self.assertLess(25**5, 2**24)
        self.assertLess(25**25, 2**120)
        success = Fraction(181, 625)
        delta = Fraction(1, 80)
        for power in range(1, 65):
            size = 25**power
            count = 10 * 25 ** (power - 1)
            incidence = 60 * power * 25 ** (power - 1)
            common = power * 25 ** (power - 1)
            self.assertEqual(count, 2 * size // 5)
            self.assertEqual(2 * incidence, 120 * common)
            self.assertEqual(power * size, 25 * common)
            defect = 1 - success**power
            self.assertGreater(defect, 0)
            self.assertLess(defect, 1)
            self.assertGreater(delta * power / defect, delta * power)

        for power in (1, 2):
            size = 25**power
            incidence = 60 * power * 25 ** (power - 1)
            self.assertGreater(1 << (2 * incidence), size**size)

    def test_ratio_diverges_from_the_certified_bounds(self) -> None:
        success = Fraction(17, 81)
        for delta in (Fraction(7, 250), Fraction(1, 40)):
            previous = Fraction(0)
            for power in range(1, 200):
                defect = 1 - success**power
                self.assertGreater(defect, 0)
                self.assertLess(defect, 1)
                ratio = delta * power / defect
                self.assertGreater(ratio, delta * power)
                self.assertGreater(ratio, previous)
                previous = ratio
            self.assertGreater(delta * 100_000 / (1 - success**100_000), 1_000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
