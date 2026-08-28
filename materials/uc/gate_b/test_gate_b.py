#!/usr/bin/env python3
"""Behavioral and adversarial checks for Gate B definitions and product proof."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
import json
from pathlib import Path
from random import Random
import sys
from tempfile import TemporaryDirectory
import unittest

import sympy
from flint import arb

HERE = Path(__file__).resolve().parent
UC_ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(UC_ROOT))

import search_n8_block_symmetric as n8
import audit_n6_square_direct as n6_square
import verify_gate_b_dyadic as dyadic
import verify_gate_b_n6_dyadic as n6_dyadic
import verify_gate_b_n6_arb as n6_arb
import verify_gate_b_rational as rational
import audit_tensorization_exact as tensor_exact
import verify_gate_b as verifier
import shapley_n7_block_symmetric as n7
import search_local_defect as localdefect
import hunt_local_ratio as localratio
import bound_local_regime as bounds
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
SECOND_BASE_ARB_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_n6_arb_v1.json"
)
TENSORIZATION_AUDIT = HERE / "candidates" / "tensorization_audit.json"
TENSORIZATION_CHECKPOINT = (
    HERE / "experiments" / "tensorization_audit_checkpoint.jsonl"
)
POWER_AUDIT = HERE / "candidates" / "power_admissibility_audit.json"
SQUARE_AUDIT = HERE / "candidates" / "square_direct_audit.json"
N6_SQUARE_AUDIT = HERE / "candidates" / "n6_square_direct_audit.json"
RATIONAL_CERTIFICATE = (
    HERE / "certificates" / "gate_b_unbounded_rational_v1.json"
)
TENSORIZATION_EXACT_AUDIT = (
    HERE / "candidates" / "tensorization_exact_audit.json"
)
LOWDEFECT_CERTIFICATE = (
    HERE / "certificates" / "gate_b_lowdefect_rational_v1.json"
)
N8_EXACT_FRONTIER = HERE / "candidates" / "n8_k2_exact_frontier.json"
N8_EXACT_CHECKPOINT = (
    HERE / "experiments" / "n8_k2_exact_checkpoint.jsonl"
)
N8_CERTIFICATE = (
    HERE / "certificates" / "gate_b_n8_rational_v1.json"
)
CLONE_CERTIFICATE = (
    HERE / "certificates" / "gate_b_n8_clone_rational_v1.json"
)
CLONE_POWER_AUDIT = HERE / "candidates" / "clone_power_audit.json"
N8MAX_CERTIFICATE = (
    HERE / "certificates" / "gate_b_n8_max_rational_v1.json"
)
N8_K3_FRONTIER = HERE / "candidates" / "n8_k3_exact_frontier.json"
N8_K4_FRONTIER = HERE / "candidates" / "n8_k4_exact_frontier.json"
N8_K4_CERTIFICATE = (
    HERE / "certificates" / "gate_b_n8_k4_rational_v1.json"
)
CLONE_SATURATION_AUDIT = HERE / "candidates" / "clone_saturation_audit.json"


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

    def test_second_base_arb_certificate_and_dyadic_dominance(self) -> None:
        result = n6_arb.verify()
        self.assertEqual(
            result["verdict"],
            "PROVED_GATE_B_UNBOUNDED_SECOND_BASE_ARB",
        )
        self.assertEqual(
            result["rational_consequence"]["a_plus_upper_lt"],
            "-17/1250",
        )
        self.assertEqual(
            result["rational_consequence"]["base_ratio_lower_gt"],
            "17/888",
        )
        self.assertEqual(
            json.loads(SECOND_BASE_ARB_CERTIFICATE.read_text()),
            result,
        )

        actual_upper = arb(result["arb"]["a_plus"]["upper"])
        self.assertLess(actual_upper, arb("-0.0136"))
        dyadic_certificate = json.loads(SECOND_BASE_CERTIFICATE.read_text())
        relaxed_upper = arb(
            dyadic_certificate["intervals"]["a_plus_upper_relaxation"][
                "upper_decimal"
            ]
        )
        self.assertGreater(relaxed_upper, actual_upper)


class ExactRationalCertificateTests(unittest.TestCase):
    """The exact rational route must bracket the true objective on both bases."""

    # An Arb ball comparison is only decidable when the two do not overlap, and
    # at points such as log2(1/2) or h(1/2) the rational primitive is exactly
    # right.  These cross-checks therefore allow a slack far below the width of
    # any bound the certificate relies on; soundness itself is argued in
    # PROOF.md, not by agreement with Arb.
    SLACK = Fraction(1, 10**30)

    def assert_at_most(self, value: Fraction, reference: arb) -> None:
        self.assertTrue(verifier.aq(value - self.SLACK) < reference)

    def assert_at_least(self, value: Fraction, reference: arb) -> None:
        self.assertTrue(verifier.aq(value + self.SLACK) > reference)

    def test_certified_log_and_exp_primitives_bracket_arb(self) -> None:
        for value in (
            Fraction(1, 625),
            Fraction(1, 45),
            Fraction(1, 3),
            Fraction(1, 2),
            Fraction(2),
            Fraction(25),
            Fraction(45),
            Fraction(625),
        ):
            reference = verifier.aq(value).log() / verifier.LOG2
            self.assert_at_most(rational.log2_lower(value), reference)
            self.assert_at_least(rational.log2_upper(value), reference)
            width = rational.log2_upper(value) - rational.log2_lower(value)
            self.assertLess(width, Fraction(1, 10**24))

        for exponent in (
            Fraction(-9, 2),
            Fraction(-1),
            Fraction(0),
            Fraction(1, 3),
            Fraction(7, 2),
        ):
            reference = (verifier.aq(exponent) * verifier.LOG2).exp()
            self.assert_at_least(rational.exp2_upper(exponent), reference)

    def test_entropy_and_fenchel_bounds_bracket_arb(self) -> None:
        for value in (
            Fraction(1, 625),
            Fraction(3, 25),
            Fraction(2, 5),
            Fraction(1, 2),
            Fraction(17, 18),
        ):
            reference = verifier.binary_entropy(value)
            self.assert_at_least(rational.entropy_upper(value), reference)
            self.assert_at_most(rational.entropy_lower(value), reference)
            self.assertLess(
                rational.entropy_upper(value) - rational.entropy_lower(value),
                Fraction(1, 10**24),
            )

        for slope in (Fraction(-4), Fraction(0), Fraction(3, 7), Fraction(6)):
            reference = (
                verifier.ONE + (verifier.aq(slope) * verifier.LOG2).exp()
            ).log() / verifier.LOG2
            self.assert_at_least(rational.softplus2_upper(slope), reference)

    def test_interval_maximum_dominates_every_feasible_action(self) -> None:
        cases = (
            (Fraction(2), Fraction(1, 5), Fraction(4, 5)),
            (Fraction(-3), Fraction(1, 25), Fraction(1, 2)),
            (Fraction(0), Fraction(0), Fraction(1)),
            (Fraction(11, 3), Fraction(1, 2), Fraction(1)),
            (Fraction(-7, 2), Fraction(0), Fraction(3, 10)),
        )
        for slope, lower, upper in cases:
            bound = rational.max_affine_entropy_upper(slope, lower, upper)
            for step in range(65):
                action = lower + (upper - lower) * Fraction(step, 64)
                actual = verifier.binary_entropy(action) + verifier.aq(
                    slope * action
                )
                self.assert_at_least(bound, actual)
            feasible = rational.feasible_near_argmax(slope, lower, upper)
            self.assertGreaterEqual(feasible, lower)
            self.assertLessEqual(feasible, upper)

    def test_rational_enclosure_contains_arb_optimum_for_every_orbit(self) -> None:
        for name in ("n6", "n7"):
            base = rational.BASES[name]
            family = base.reconstruct()
            representatives, group = rational.order_orbits(
                family, base.dimension
            )
            arb_roots, _states, _ambiguous = verifier.certified_one_sided_costs(
                family, base.dimension, representatives
            )
            for order, reference in zip(representatives, arb_roots):
                low, high, _order_states = rational.bellman_bounds(
                    family, base.dimension, order
                )
                # The rational pair encloses the true optimum; the Arb value is
                # itself an upper bound for it, so it must sit above the lower
                # endpoint and agree with the upper endpoint to high precision.
                self.assertLess(low, high)
                self.assertLess(high - low, Fraction(1, 10**24))
                self.assertLessEqual(verifier.aq(low), reference.lower())
                self.assertLessEqual(
                    reference.upper(),
                    verifier.aq(high + Fraction(1, 10**20)),
                )
            self.assertEqual(
                len(representatives) * len(group), factorial(base.dimension)
            )

    def test_rational_upper_bound_is_sharper_than_dyadic_relaxation(self) -> None:
        base = rational.BASES["n6"]
        family = base.reconstruct()
        representatives, _group = rational.order_orbits(family, base.dimension)
        for order in representatives:
            _low, high, _states = rational.bellman_bounds(
                family, base.dimension, order
            )
            relaxed, _relaxed_states = dyadic.relaxed_bellman_cost(family, order)
            self.assertLess(high, Fraction(relaxed.hi, dyadic.SCALE))

    def test_state_canonicalisation_does_not_change_the_value(self) -> None:
        base = rational.BASES["n6"]
        family = base.reconstruct()
        order = (3, 0, 4, 1, 5, 2)
        canonical = rational.bellman_bounds(family, base.dimension, order)
        direct = rational.bellman_bounds(
            family, base.dimension, order, canonicalize=False
        )
        self.assertEqual(canonical[:2], direct[:2])
        self.assertGreaterEqual(direct[2], canonical[2])

    def test_orbit_quotient_reproduces_the_all_order_certificate(self) -> None:
        certificate = json.loads(RATIONAL_CERTIFICATE.read_text())
        recomputed = rational.verify(
            base_names=("n6", "n7"),
            order_mode="orbits",
            checkpoint=None,
            sample_recheck=0,
        )
        for name in ("n6", "n7"):
            stored = certificate["bases"][name]
            fresh = recomputed["bases"][name]
            self.assertEqual(
                stored["order_evaluation"]["evaluated_order_count"],
                factorial(rational.BASES[name].dimension),
            )
            self.assertFalse(
                stored["order_evaluation"]["symmetry_quotient_used"]
            )
            self.assertEqual(
                stored["order_evaluation"]["distinct_enclosure_multiplicities"],
                [stored["order_evaluation"]["automorphism_count"]]
                * stored["order_evaluation"]["order_orbit_count"],
            )
            for key in ("a_plus_lower", "a_plus_upper", "shapley_iid_upper"):
                self.assertEqual(
                    stored["rational_values"][key],
                    fresh["rational_values"][key],
                )

    def test_certificate_is_internally_consistent(self) -> None:
        certificate = json.loads(RATIONAL_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL"
        )
        self.assertEqual(certificate["third_party_dependencies"], [])
        for name, targets in (
            ("n6", ["-17/1250", "-1/80"]),
            ("n7", ["-7/250", "-1/40"]),
        ):
            payload = certificate["bases"][name]
            values = payload["rational_values"]
            low = Fraction(values["a_plus_lower"])
            high = Fraction(values["a_plus_upper"])
            self.assertLess(low, high)
            self.assertLess(high, 0)
            self.assertEqual(
                payload["rational_consequence"]["a_plus_upper_lt"], targets
            )
            self.assertTrue(
                payload["rational_consequence"]["enclosure_is_negative"]
            )
            for target in targets:
                self.assertLess(high, Fraction(target))

            iid = Fraction(values["shapley_iid_upper"])
            bellman = Fraction(values["bellman_average_upper"])
            entropy = Fraction(values["family_entropy_lower"])
            alpha = Fraction(
                certificate["alpha"]["numerator"],
                certificate["alpha"]["denominator"],
            )
            self.assertEqual(
                high,
                rational.round_up(
                    (1 - alpha) * iid + alpha * bellman - entropy,
                    rational.ACCUMULATOR_BITS,
                ),
            )
            self.assertEqual(
                Fraction(payload["exact_base"]["closure_defect"]),
                1 - Fraction(payload["exact_base"]["join_success_probability"]),
            )

    def test_checkpoint_rejects_foreign_and_truncated_records(self) -> None:
        base = rational.BASES["n6"]
        with TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "records.jsonl"
            first = rational.evaluate_base(
                base, "orbits", checkpoint, sample_recheck=0, progress=0
            )
            self.assertEqual(
                first["order_evaluation"]["new_evaluations"],
                first["order_evaluation"]["evaluated_order_count"],
            )

            resumed = rational.evaluate_base(
                base, "orbits", checkpoint, sample_recheck=2, progress=0
            )
            self.assertEqual(resumed["order_evaluation"]["new_evaluations"], 0)
            self.assertEqual(resumed["order_evaluation"]["rechecked_records"], 2)
            self.assertEqual(
                resumed["rational_values"]["a_plus_upper"],
                first["rational_values"]["a_plus_upper"],
            )

            lines = checkpoint.read_text().splitlines()
            foreign = json.loads(lines[0])
            foreign["schema"] = rational.RECORD_SCHEMA + 1
            with checkpoint.open("a") as handle:
                handle.write(json.dumps(foreign, sort_keys=True) + "\n")
                handle.write('{"kind": "order", "base": "n6"')
            partial = rational.evaluate_base(
                base, "orbits", checkpoint, sample_recheck=0, progress=0
            )
            self.assertEqual(partial["order_evaluation"]["new_evaluations"], 0)
            self.assertEqual(
                partial["order_evaluation"]["foreign_schema_checkpoint_lines"], 1
            )
            self.assertEqual(
                partial["order_evaluation"]["truncated_checkpoint_lines"], 1
            )
            self.assertEqual(
                partial["rational_values"]["a_plus_lower"],
                first["rational_values"]["a_plus_lower"],
            )

    def test_linear_growth_corollary_is_exact(self) -> None:
        certificate = json.loads(RATIONAL_CERTIFICATE.read_text())
        for name, dimension, delta in (
            ("n6", 6, Fraction(17, 1250)),
            ("n7", 7, Fraction(7, 250)),
        ):
            payload = certificate["bases"][name]
            success = Fraction(
                payload["exact_base"]["join_success_probability"]
            )
            self.assertEqual(
                certificate["asymptotic"]["bases"][name]["asymptotic_slope"],
                str(delta / dimension),
            )
            for coordinates in range(dimension, 40 * dimension + 1, dimension):
                power = coordinates // dimension
                defect = 1 - success**power
                self.assertGreater(defect, 0)
                self.assertLess(defect, 1)
                ratio = delta * power / defect
                self.assertGreater(ratio, delta * power)
                # Linear lower bound in the number of coordinates.
                self.assertGreater(ratio, delta * coordinates / dimension - 1)
            # Matching upper bound under a defect floor: -A_+ <= log2 m <= n.
            floor = Fraction(1, 2)
            self.assertLess(
                delta * 1000 / (1 - success**1000),
                Fraction(dimension * 1000) / floor,
            )

    def test_exact_tensorization_audit_found_no_counterexample(self) -> None:
        summary = json.loads(TENSORIZATION_EXACT_AUDIT.read_text())
        self.assertEqual(summary["verdict"], "NO_EXACT_COUNTEREXAMPLE")
        self.assertEqual(summary["arithmetic"], "exact_python_fractions")
        self.assertEqual(summary["third_party_dependencies"], [])
        self.assertEqual(summary["bellman_intersection_failures"], 0)
        self.assertEqual(summary["defect_identity_failures"], 0)
        self.assertEqual(summary["q_intersection_failures"], 0)
        self.assertGreaterEqual(summary["total_orders_checked"], 15_000)
        self.assertGreaterEqual(summary["unequal_dimension_cases"], 20)
        self.assertEqual(
            Fraction(summary["worst_bellman_gap_decimal"]), Fraction(0)
        )
        self.assertLess(
            Fraction(summary["worst_q_gap_decimal"]), Fraction(1, 10**24)
        )
        for case in summary["cases"]:
            self.assertTrue(case["defect_identity_exact"])
            self.assertTrue(case["q_enclosures_intersect"])
            self.assertEqual(case["bellman_intersection_failures"], 0)
            self.assertEqual(
                case["orders_checked"],
                factorial(case["left_dimension"] + case["right_dimension"]),
            )

    def test_exact_audit_replays_one_case_from_scratch(self) -> None:
        summary = json.loads(TENSORIZATION_EXACT_AUDIT.read_text())
        target = next(
            case
            for case in summary["cases"]
            if case["name"] == "triple_antichain_x_triple_chain"
        )
        replay = tensor_exact.evaluate_case(
            {
                "name": target["name"],
                "left": tuple(target["left_rows"]),
                "left_dimension": target["left_dimension"],
                "right": tuple(target["right_rows"]),
                "right_dimension": target["right_dimension"],
            }
        )
        self.assertEqual(replay, target)
        self.assertEqual(replay["orders_checked"], 720)

    def test_induced_order_extraction_is_faithful(self) -> None:
        order = (4, 0, 3, 1, 5, 2)
        left, right = tensor_exact.induced_orders(order, 3)
        self.assertEqual(left, (0, 1, 2))
        self.assertEqual(right, (1, 0, 2))
        self.assertEqual(tensor_exact._relabel(right), (1, 0, 2))
        self.assertEqual(tensor_exact._relabel((4, 0, 3)), (2, 0, 1))


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

    def test_second_base_direct_square_audit(self) -> None:
        summary = json.loads(N6_SQUARE_AUDIT.read_text())
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["checked_order_count"], 5)
        self.assertTrue(summary["all_within_tolerance"])
        self.assertLess(summary["worst_bellman_gap"], 1e-12)
        self.assertLess(summary["worst_iid_gap"], 1e-12)

        exact = summary["exact_square"]
        self.assertEqual(exact["family_size"], 625)
        self.assertEqual(exact["dimension"], 12)
        self.assertEqual(exact["coordinate_counts"], [250] * 12)
        self.assertEqual(exact["cap_bound"], 250)
        self.assertTrue(exact["cap_holds_with_equality"])
        self.assertTrue(exact["reimer_holds_strictly"])
        self.assertEqual(exact["missing_ordered_join_pairs"], 357_864)
        self.assertEqual(exact["closure_defect"], "357864/390625")
        self.assertTrue(exact["closure_defect_matches_product"])
        self.assertTrue(exact["active"] and exact["separating"])
        for order in summary["orders"].values():
            self.assertAlmostEqual(
                order["direct_square_bellman"],
                order["base_block_bellman_sum"],
                places=12,
            )

    def test_second_base_square_checkpoint_tolerates_truncated_tail(self) -> None:
        with TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "checkpoint.jsonl"
            checkpoint.write_text(
                '{"order_name": "complete", "value": 1}\n'
                '{"order_name": "truncated"'
            )
            original = n6_square.CHECKPOINT
            try:
                n6_square.CHECKPOINT = checkpoint
                records = n6_square.load_records()
            finally:
                n6_square.CHECKPOINT = original
        self.assertEqual(
            records,
            {"complete": {"order_name": "complete", "value": 1}},
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


class LocalRegimeTests(unittest.TestCase):
    """The local-regime lemmas and the instruments that measure the regime."""

    def test_reimer_threshold_agrees_across_implementations(self) -> None:
        for size in range(3, 129):
            self.assertEqual(
                localdefect.reimer_threshold(size),
                n7.exact_reimer_threshold(size),
            )

    def test_defect_is_zero_exactly_on_union_closed_families(self) -> None:
        closed = (0, 1, 2, 3)  # a full sublattice on two coordinates
        self.assertEqual(localdefect.defect_count(closed), 0)
        open_family = (1, 2)
        self.assertEqual(localdefect.defect_count(open_family), 2)

    def test_lemma_4_failures_come_in_pairs(self) -> None:
        random = Random(20260827)
        checked = 0
        for _trial in range(400):
            size = random.randrange(3, 8)
            rows = tuple(sorted(random.sample(range(16), size)))
            bad = localdefect.defect_count(rows)
            self.assertEqual(bad % 2, 0)
            if bad:
                self.assertGreaterEqual(bad, 2)
                self.assertGreaterEqual(
                    Fraction(bad, size ** 2), Fraction(2, size ** 2))
                checked += 1
        self.assertGreater(checked, 100)

    def test_lemma_5_defect_multiplies_through_success(self) -> None:
        left_families = ((1, 2), (0, 1, 2), (1, 2, 3))
        right_families = ((1, 2), (0, 1, 3), (0, 1, 2, 3))
        for left in left_families:
            for right in right_families:
                product = tuple(sorted(
                    a | (b << 2) for a in left for b in right
                ))
                joint = Fraction(
                    localdefect.defect_count(product), len(product) ** 2)
                factor_left = Fraction(
                    localdefect.defect_count(left), len(left) ** 2)
                factor_right = Fraction(
                    localdefect.defect_count(right), len(right) ** 2)
                self.assertEqual(
                    1 - joint, (1 - factor_left) * (1 - factor_right))
                self.assertGreaterEqual(joint, factor_left)
                self.assertGreaterEqual(joint, factor_right)

    def test_certified_bases_have_the_published_defects(self) -> None:
        base_six = localratio.known_base(6)
        self.assertEqual(len(base_six), 25)
        self.assertEqual(
            Fraction(localdefect.defect_count(base_six), 625),
            Fraction(444, 625),
        )
        base_seven = localratio.known_base(7)
        self.assertEqual(len(base_seven), 45)
        self.assertEqual(
            Fraction(localdefect.defect_count(base_seven), 2025),
            Fraction(64, 81),
        )
        for rows, dimension in ((base_six, 6), (base_seven, 7)):
            self.assertTrue(localdefect.admissible(rows, dimension))
            self.assertTrue(localdefect.normalized(rows, dimension))

    def test_calibration_ratios_match_the_reported_constants(self) -> None:
        for dimension, defect, expected in (
            (6, Fraction(444, 625), 0.0192456471),
            (7, Fraction(64, 81), 0.0362591270),
        ):
            rows = localratio.known_base(dimension)
            a_plus, _q, _c = localratio.full_objective(rows, dimension)
            self.assertLess(a_plus, 0.0)
            self.assertAlmostEqual(
                localratio.ratio(a_plus, defect), expected, places=7)

    def test_exact_enumeration_reproduces_the_small_dimension_minima(self) -> None:
        result = localdefect.exact_minimum(4)
        self.assertTrue(result["complete"])
        expected = {"3": "2/9", "5": "6/25", "8": "1/4"}
        for size, defect in expected.items():
            self.assertEqual(result["sizes"][size]["min_defect"], defect)
        for entry in result["sizes"].values():
            rows = tuple(entry["witness"])
            self.assertTrue(localdefect.admissible(rows, 4))
            self.assertEqual(
                localdefect.defect_count(rows), entry["min_defect_count"])

    def test_witness_reader_accepts_foreign_schema_but_rejects_bad_families(
        self,
    ) -> None:
        """Schema governs work skipping; only re-verification governs evidence."""
        good = localratio.known_base(6)
        inadmissible = tuple(range(25))  # violates the coordinate cap
        self.assertFalse(localdefect.admissible(inadmissible, 6))
        with TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.jsonl"
            path.write_text("\n".join(json.dumps(payload, sort_keys=True) for payload in (
                {
                    "schema": localdefect.RECORD_SCHEMA,
                    "dimension": 6, "size": 25, "restart": 0,
                    "best_family": list(inadmissible),
                    "best_defect_count": 0,
                },
                {
                    "schema": localdefect.RECORD_SCHEMA - 1,
                    "dimension": 6, "size": 25, "restart": 1,
                    "best_family": list(good),
                    "best_defect_count": 444,
                },
            )) + "\n{ truncated")

            done, foreign, truncated = localdefect.read_checkpoint(path)
            self.assertEqual(len(done), 1)
            self.assertEqual(foreign, 1)
            self.assertEqual(truncated, 1)

            witnesses = localdefect.read_witnesses(path, 6)
            self.assertIn(25, witnesses)
            self.assertEqual(witnesses[25]["family"], list(good))
            self.assertEqual(witnesses[25]["defect_count"], 444)
            self.assertEqual(
                witnesses[25]["from_schema"], localdefect.RECORD_SCHEMA - 1)

    def test_seeder_reaches_every_feasible_size(self) -> None:
        for dimension in (6, 7):
            for size in localdefect.feasible_sizes(dimension):
                generator = Random(localdefect.restart_seed(1, dimension, size, 0))
                degrees = localdefect.degree_sequences(size, dimension, generator)
                rows = localdefect.seed_family(size, dimension, degrees, generator)
                self.assertIsNotNone(rows, (dimension, size))
                self.assertTrue(localdefect.admissible(rows, dimension))

    def test_sweep_respects_the_hard_defect_cap(self) -> None:
        budget = Fraction(3, 4)
        outcome = localratio.sweep(
            6, 25, budget, localratio.known_base(6),
            Random(4242), 25, 12,
        )
        rows = tuple(outcome["family"])
        self.assertTrue(localdefect.admissible(rows, 6))
        self.assertLessEqual(Fraction(outcome["defect"]), budget)
        self.assertEqual(
            Fraction(localdefect.defect_count(rows), 625),
            Fraction(outcome["defect"]),
        )

    def test_lemma_6_ratio_bound_holds_on_every_recorded_witness(self) -> None:
        checkpoint = HERE / "experiments" / "local_defect_checkpoint.jsonl"
        if not checkpoint.exists():
            self.skipTest("local defect checkpoint not present")
        for dimension in (6, 7):
            for size, row in localdefect.read_witnesses(checkpoint, dimension).items():
                defect = Fraction(row["defect_count"], size ** 2)
                self.assertGreaterEqual(defect, Fraction(2, size ** 2))
                self.assertLessEqual(
                    Fraction(1) / defect,
                    Fraction(size ** 2, 2),
                )

    def test_lowdefect_certificate_is_negative_at_the_claimed_defect(self) -> None:
        """The third base must certify negativity at defect 1336/2025."""
        certificate = json.loads(LOWDEFECT_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL")
        base = certificate["bases"]["n7lo"]
        facts = base["exact_base"]
        rows = tuple(facts["family_rows"])

        self.assertEqual(len(rows), 45)
        self.assertTrue(localdefect.admissible(rows, 7))
        self.assertTrue(localdefect.normalized(rows, 7))
        self.assertEqual(facts["coordinate_counts"], [18] * 7)
        self.assertEqual(facts["total_incidence"], 126)
        self.assertEqual(facts["reimer_threshold"], 124)
        self.assertEqual(facts["cap_bound"], 18)

        defect = Fraction(1336, 2025)
        self.assertEqual(Fraction(facts["closure_defect"]), defect)
        self.assertEqual(localdefect.defect_count(rows), 1336)

        # The whole point of this base is the defect, not the sharpness: it must
        # sit strictly below both published bases.
        self.assertLess(defect, Fraction(444, 625))
        self.assertLess(defect, Fraction(64, 81))

        consequence = base["rational_consequence"]
        self.assertTrue(consequence["enclosure_is_negative"])
        self.assertEqual(consequence["sharpest_frozen_target_passed"], "-1/1200")
        lower, upper = (
            Fraction(value)
            for value in consequence["a_plus_two_sided_enclosure"]
        )
        self.assertLess(lower, upper)
        self.assertLess(upper, Fraction(-1, 1200))
        self.assertLess(upper - lower, Fraction(1, 10 ** 26))

        evaluation = base["order_evaluation"]
        self.assertEqual(evaluation["all_order_count"], 5040)
        self.assertEqual(evaluation["automorphism_count"], 1)
        self.assertEqual(evaluation["distinct_enclosure_count"], 5040)
        self.assertFalse(evaluation["symmetry_quotient_used"])

    def test_lowdefect_base_reconstructs_without_a_cell_pattern(self) -> None:
        base = rational.BASES["n7lo"]
        self.assertIsInstance(base, rational.ExplicitBase)
        facts = base.exact_facts()
        self.assertIsNone(facts["block_cells"])
        self.assertIsNone(facts["block_partition"])
        self.assertEqual(facts["missing_ordered_join_pairs"], 1336)
        self.assertEqual(base.reconstruct(), tuple(sorted(base.expected_rows)))


class SizeCeilingAndDefectFloorTests(unittest.TestCase):
    """Lemma 9, Corollary 10, Lemma 11, Lemma 13 and Proposition 14."""

    def test_size_ceiling_holds_for_every_admissible_size(self) -> None:
        """log2 m <= 4n/5, as the exact integer test m^5 <= 2^(4n)."""
        for dimension in range(3, 13):
            for size in bounds.feasible_sizes(dimension):
                self.assertTrue(
                    bounds.size_ceiling_holds(size, dimension),
                    f"ceiling failed at n={dimension} m={size}",
                )
            largest = bounds.max_admissible_size(dimension)
            # Sharpness: one size above the ceiling must fail the integer test.
            self.assertFalse(bounds.size_ceiling_holds(1 << dimension, dimension))
            self.assertTrue(bounds.size_ceiling_holds(largest, dimension))

    def test_size_ceiling_matches_the_reimer_witness_strings(self) -> None:
        """Each base's reimer_witness is its instance of Lemma 9."""
        for name, expected in (
            ("n6", (25, 6)), ("n7", (45, 7)), ("n7lo", (45, 7)),
            ("n8lo", (75, 8)), ("n8hi", (75, 8)),
        ):
            size, dimension = expected
            base = rational.BASES[name]
            self.assertEqual(len(base.reconstruct()), size)
            self.assertEqual(base.dimension, dimension)
            self.assertEqual(base.reimer_witness, f"{size}**5 < 2**{4 * dimension}")
            self.assertLess(size ** 5, 1 << (4 * dimension))

    def test_defect_floors_hold_over_every_small_admissible_family(self) -> None:
        """Lemma 11 and Lemma 13 never exceed the true defect, and 13 is tight."""
        tight = False
        for dimension in (3, 4):
            for size in bounds.feasible_sizes(dimension):
                for rows in combinations(range(1 << dimension), size):
                    if not bounds.admissible(rows, dimension):
                        continue
                    observed = bounds.defect(rows)
                    growth = bounds.union_growth_floor(rows, dimension)
                    if growth is not None:
                        self.assertGreaterEqual(observed, growth)
                    capacity = bounds.capacity_floor(rows)
                    self.assertGreaterEqual(observed, capacity)
                    if observed == capacity:
                        tight = True
        self.assertTrue(tight, "Lemma 13 should be attained with equality")

    def test_certified_bases_admit_no_addition(self) -> None:
        """Proposition 14, checked over every candidate set."""
        for name in ("n6", "n7", "n7lo"):
            base = rational.BASES[name]
            rows = base.reconstruct()
            self.assertEqual(len(rows), bounds.max_admissible_size(base.dimension))
            self.assertEqual(
                set(bounds.degrees(rows, base.dimension)),
                {bounds.cap(len(rows))},
            )
            self.assertEqual(bounds.admissible_additions(rows, base.dimension), ())

    def test_dominant_set_is_present_exactly_when_the_floor_is_vacuous(self) -> None:
        """Corollary 12: a low-defect family needs a set of size >= (8/5) sbar."""
        for name in ("n6", "n7", "n7lo", "n8lo", "n8hi"):
            base = rational.BASES[name]
            rows = base.reconstruct()
            largest = max(bin(row).count("1") for row in rows)
            threshold = Fraction(8, 5) * Fraction(
                bounds.incidence(rows), len(rows)
            )
            floor_value = bounds.union_growth_floor(rows, base.dimension)
            self.assertEqual(floor_value is None, largest >= threshold)

    def test_bound_primitives_agree_with_the_search_module(self) -> None:
        """The standalone checker must not drift from the search implementation."""
        for name in ("n6", "n7", "n7lo"):
            base = rational.BASES[name]
            rows = base.reconstruct()
            self.assertEqual(
                bounds.defect_count(rows), localdefect.defect_count(rows)
            )
            self.assertEqual(bounds.defect(rows), localdefect.defect(rows))
            self.assertEqual(
                bounds.admissible(rows, base.dimension),
                localdefect.admissible(rows, base.dimension),
            )
        for size in range(3, 90):
            self.assertEqual(
                bounds.reimer_threshold(size), localdefect.reimer_threshold(size)
            )


class EightCoordinateFrontierTests(unittest.TestCase):
    """Proposition 15: the n=8 witnesses and the exact re-ranking."""

    def test_n8_bases_are_admissible_normalized_and_at_the_claimed_defect(self) -> None:
        for name, defect_value, missing in (
            ("n8lo", Fraction(1144, 1875), 3_432),
            ("n8hi", Fraction(468, 625), 4_212),
        ):
            base = rational.BASES[name]
            rows = base.reconstruct()
            facts = base.exact_facts()
            self.assertEqual(len(rows), 75)
            self.assertTrue(bounds.admissible(rows, 8))
            self.assertTrue(localdefect.normalized(rows, 8))
            self.assertEqual(facts["coordinate_counts"], [30] * 8)
            self.assertEqual(facts["total_incidence"], 240)
            self.assertEqual(facts["reimer_threshold"], 234)
            self.assertEqual(Fraction(facts["closure_defect"]), defect_value)
            self.assertEqual(bounds.defect_count(rows), missing)
            self.assertEqual(facts["missing_ordered_join_pairs"], missing)

    def test_size_seventy_five_is_infeasible_at_seven_coordinates(self) -> None:
        """The eighth coordinate is what admits these sizes at all."""
        for size in (70, 75):
            self.assertGreater(bounds.reimer_threshold(size), 7 * bounds.cap(size))
            self.assertLessEqual(bounds.reimer_threshold(size), 8 * bounds.cap(size))
            self.assertNotIn(size, bounds.feasible_sizes(7))
            self.assertIn(size, bounds.feasible_sizes(8))

    def test_exact_frontier_report_improves_both_published_numbers(self) -> None:
        report = json.loads(N8_EXACT_FRONTIER.read_text())
        self.assertEqual(report["candidates_considered"], 36)
        self.assertEqual(report["float_negative_count"], 21)
        self.assertEqual(report["certified_negative_count"], 20)
        self.assertEqual(report["foreign_schema_lines"], 0)
        self.assertEqual(report["truncated_lines"], 0)
        frontier = report["lowest_defect_certified_negative"]
        self.assertTrue(frontier["certified_negative"])
        self.assertTrue(frontier["beats_published_frontier_defect"])
        self.assertLess(
            Fraction(frontier["exact_facts"]["closure_defect"]), Fraction(1336, 2025)
        )

        sharpest = report["highest_ratio_certified_negative"]
        self.assertTrue(sharpest["beats_published_ratio"])
        self.assertGreater(
            Fraction(sharpest["exact_ratio_lower"]),
            Fraction(286_491, 10_000_000) / Fraction(64, 81),
        )

    def test_one_float_negative_is_exactly_positive(self) -> None:
        """The float64 census produced a false negative; the exact route caught it.

        The mirror direction matters just as much: no family that reads
        non-negative in float64 is exactly negative, so the frontier is not
        hiding behind a rounding error in the other direction.
        """
        report = json.loads(N8_EXACT_FRONTIER.read_text())
        self.assertEqual(report["float_negative_but_exactly_nonnegative"], ["2553/3200"])
        self.assertEqual(report["float_nonnegative_but_exactly_negative"], [])
        self.assertGreater(
            Fraction(report["largest_float_exact_gap_decimal"]), Fraction(1, 1000)
        )

    def test_every_certified_record_carries_an_exact_two_sided_enclosure(self) -> None:
        records = [
            json.loads(line)
            for line in N8_EXACT_CHECKPOINT.read_text().splitlines()
            if line.strip()
        ]
        self.assertEqual(len(records), 36)
        for record in records:
            lower = Fraction(record["a_plus_lower"])
            upper = Fraction(record["a_plus_upper"])
            self.assertLessEqual(lower, upper)
            self.assertLess(upper - lower, Fraction(1, 10 ** 24))
            self.assertEqual(record["certified_negative"], upper < 0)
            facts = record["exact_facts"]
            self.assertTrue(facts["admissible"])
            self.assertTrue(facts["size_ceiling_holds"])
            rows = tuple(record["family_rows"])
            self.assertEqual(Fraction(facts["closure_defect"]), bounds.defect(rows))

    def test_n8_certificate_beats_its_frozen_targets(self) -> None:
        certificate = json.loads(N8_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL"
        )
        self.assertEqual(certificate["third_party_dependencies"], [])
        expected = {
            "n8lo": (Fraction(1144, 1875), "-1/500", Fraction(1336, 2025)),
            "n8hi": (Fraction(468, 625), "-693/25000", None),
        }
        for name, (defect_value, target, frontier) in expected.items():
            base = certificate["bases"][name]
            self.assertEqual(
                Fraction(base["exact_base"]["closure_defect"]), defect_value
            )
            consequence = base["rational_consequence"]
            self.assertTrue(consequence["enclosure_is_negative"])
            self.assertEqual(consequence["sharpest_frozen_target_passed"], target)
            lower, upper = (
                Fraction(value)
                for value in consequence["a_plus_two_sided_enclosure"]
            )
            self.assertLess(lower, upper)
            self.assertLess(upper, Fraction(target))
            self.assertLess(upper - lower, Fraction(1, 10 ** 24))
            if frontier is not None:
                self.assertLess(defect_value, frontier)

    def test_n8hi_certified_ratio_beats_the_published_base(self) -> None:
        """The rational ratio bound, not just the decimal, must clear 0.0362591."""
        certificate = json.loads(N8_CERTIFICATE.read_text())
        ratio = Fraction(
            certificate["bases"]["n8hi"]["rational_consequence"]["base_ratio_lower_gt"]
        )
        published = Fraction(286_491, 10_000_000) / Fraction(64, 81)
        self.assertGreater(ratio, published)
        self.assertEqual(ratio, Fraction(77, 2080))

    def test_the_symmetry_quotient_is_exact(self) -> None:
        """Orbits all have size |Aut| and the Bellman value is constant on them.

        The n=8 certificate is the first here to use a nontrivial quotient, so
        the two facts that make a uniform average over representatives equal to
        the average over all n! orders are checked rather than assumed.
        """
        for name in ("n8lo", "n8hi"):
            base = rational.BASES[name]
            rows = base.reconstruct()
            representatives, group = rational.order_orbits(rows, base.dimension)
            self.assertEqual(
                len(representatives) * len(group), factorial(base.dimension)
            )
            self.assertEqual(len(group), 1440)
            random = Random(20260827)
            for representative in representatives[:2]:
                reference = rational.bellman_bounds(
                    rows, base.dimension, representative
                )[:2]
                for _ in range(3):
                    element = group[random.randrange(len(group))]
                    mate = tuple(element[index] for index in representative)
                    self.assertEqual(
                        rational.bellman_bounds(rows, base.dimension, mate)[:2],
                        reference,
                    )

    def test_n8max_sits_at_the_maximum_admissible_size(self) -> None:
        """The best separating ratio is attained at m=80 with maximum incidence."""
        base = rational.BASES["n8max"]
        rows = base.reconstruct()
        facts = base.exact_facts()
        self.assertEqual(len(rows), 80)
        self.assertEqual(len(rows), bounds.max_admissible_size(8))
        self.assertTrue(bounds.admissible(rows, 8))
        self.assertTrue(facts["separating"])
        self.assertEqual(facts["coordinate_counts"], [32] * 8)
        self.assertEqual(facts["cap_bound"], 32)
        # 8 * 32 = 256 is the largest incidence any n=8 family can carry.
        self.assertEqual(facts["total_incidence"], 256)
        self.assertEqual(facts["total_incidence"], 8 * bounds.cap(80))
        self.assertEqual(facts["reimer_threshold"], 253)
        self.assertEqual(Fraction(facts["closure_defect"]), Fraction(2343, 3200))

    def test_n8max_has_the_best_separating_ratio(self) -> None:
        certificate = json.loads(N8MAX_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL"
        )
        consequence = certificate["bases"]["n8max"]["rational_consequence"]
        ratio = Fraction(consequence["base_ratio_lower_gt"])
        self.assertEqual(ratio, Fraction(2432, 58575))
        # beats the published base and the earlier separating record n8hi
        self.assertGreater(ratio, Fraction(286_491, 10_000_000) / Fraction(64, 81))
        self.assertGreater(ratio, Fraction(693, 25_000) / Fraction(468, 625))
        # but does not improve the asymptotic slope
        slope = Fraction(
            certificate["asymptotic"]["bases"]["n8max"]["asymptotic_slope"]
        )
        self.assertEqual(slope, Fraction(19, 5000))
        self.assertLess(slope, Fraction(7, 250) / 7)
        self.assertLess(slope, Fraction(3, 640))

    def test_k3_class_is_exactly_re_ranked(self) -> None:
        """k=3 exposes the float evaluator failing in *both* directions."""
        report = json.loads(N8_K3_FRONTIER.read_text())
        self.assertEqual(report["float_negative_count"], 54)
        self.assertEqual(report["certified_negative_count"], 52)
        # Three float negatives are exactly non-negative, so only 51 of the 54
        # survive; the 52nd certified negative is a family the census called
        # non-negative.
        self.assertEqual(len(report["float_negative_but_exactly_nonnegative"]), 3)
        # ... and one family the float census called non-negative is exactly
        # negative.  That is a *missed* witness, not a false alarm, and it is why
        # the re-ranking threshold must sit above zero.
        self.assertEqual(
            report["float_nonnegative_but_exactly_negative"], ["94/125"]
        )
        self.assertGreater(
            Fraction(report["largest_float_exact_gap_decimal"]), Fraction(3, 1000)
        )
        # k=3 improves neither the defect frontier nor the slope
        frontier = report["lowest_defect_certified_negative"]
        self.assertFalse(frontier["beats_published_frontier_defect"])
        self.assertGreater(
            Fraction(frontier["exact_facts"]["closure_defect"]), Fraction(139, 245)
        )

    def test_census_block_orders_are_only_a_screen(self) -> None:
        """The census order set is not an automorphism-orbit system.

        This pins the diagnosis rather than the symptom.  The census averages
        `C_+` over one order per `S_k x S_(8-k)` pattern, which is a valid
        representative system only when the family's automorphism group *is* that
        block group.  For `n8lo` the group has the same order as `S_2 x S_6` but
        is not it, so the block orders badly mis-weight the true orbits.  If this
        test ever starts passing trivially, the census objective has changed and
        the screening threshold has to be revisited.
        """
        base = rational.BASES["n8lo"]
        rows = base.reconstruct()
        group = rational.automorphism_group(rows, 8)
        self.assertEqual(len(group), 1440)

        preserving = [
            g for g in group
            if set(g[0:2]) == {0, 1} and set(g[2:8]) == set(range(2, 8))
        ]
        # same order as |S_2 x S_6| = 1440, but not that group
        self.assertEqual(len(preserving), 120)
        self.assertLess(len(preserving), len(group))

        def canonical(order):
            return min(tuple(g[i] for i in order) for g in group)

        block_orders = n8.block_order_representatives(2)
        orbit_orders, _group = rational.order_orbits(rows, 8)
        self.assertEqual(len(block_orders), len(orbit_orders))
        # the honest orbit system covers every orbit exactly once ...
        self.assertEqual(len({canonical(o) for o in orbit_orders}), 28)
        # ... the census one covers only 12 of the 28
        self.assertEqual(len({canonical(tuple(o)) for o in block_orders}), 12)

    def test_screening_threshold_exceeds_the_census_error_bound(self) -> None:
        """Any subset average is within alpha*(max C_+ - min C_+) of the truth."""
        for name in ("n8lo", "n8hi", "n8clone_hi", "n8max"):
            rows = rational.BASES[name].reconstruct()
            representatives, _group = rational.order_orbits(rows, 8)
            values = [
                rational.bellman_bounds(rows, 8, order)[1]
                for order in representatives
            ]
            bound = PROJECT_ALPHA * (max(values) - min(values))
            self.assertLess(bound, Fraction(1, 100))


class ClonedCoordinateGrowthTests(unittest.TestCase):
    """Proposition 16: separation is unused, and the growth constant improves."""

    def test_cloned_bases_are_admissible_and_not_separating(self) -> None:
        for name, defect_value, missing in (
            ("n8clone_lo", Fraction(139, 245), 2_780),
            ("n8clone_hi", Fraction(317, 490), 3_170),
        ):
            base = rational.BASES[name]
            rows = base.reconstruct()
            facts = base.exact_facts()
            self.assertEqual(len(rows), 70)
            self.assertTrue(bounds.admissible(rows, 8))
            self.assertFalse(facts["separating"])
            self.assertEqual(facts["duplicate_column_pairs"], [[0, 1]])
            self.assertTrue(facts["active"])
            self.assertEqual(facts["coordinate_counts"], [28] * 8)
            self.assertEqual(facts["cap_bound"], 28)
            self.assertEqual(facts["total_incidence"], 224)
            self.assertEqual(facts["reimer_threshold"], 215)
            self.assertEqual(Fraction(facts["closure_defect"]), defect_value)
            self.assertEqual(bounds.defect_count(rows), missing)

    def test_a_declared_separating_base_still_asserts_separation(self) -> None:
        """Dropping the check must not silently weaken the normalized bases."""
        self.assertTrue(rational.Base.requires_separation)
        self.assertTrue(rational.ExplicitBase.requires_separation)
        self.assertFalse(rational.ClonedCoordinateBase.requires_separation)
        for name in ("n6", "n7", "n7lo", "n8lo", "n8hi"):
            self.assertTrue(rational.BASES[name].exact_facts()["separating"])
            self.assertEqual(
                rational.BASES[name].exact_facts()["duplicate_column_pairs"], []
            )
        # A cloned family declared as separating must fail loudly.
        rows = rational.BASES["n8clone_hi"].reconstruct()
        impostor = rational.ExplicitBase(
            name="impostor", dimension=8, rows=rows,
            expected_count=28, expected_missing=3_170,
            reimer_witness="70**5 < 2**32",
        )
        with self.assertRaises(AssertionError):
            impostor.exact_facts()

    def test_reimer_for_the_base_certifies_every_power(self) -> None:
        """m^m <= 2^(2I) is exactly I >= R_m, and it is what powers need."""
        for name in ("n6", "n7", "n7lo", "n8lo", "n8hi", "n8clone_lo", "n8clone_hi"):
            base = rational.BASES[name]
            rows = base.reconstruct()
            size = len(rows)
            incidence = bounds.incidence(rows)
            self.assertGreaterEqual(incidence, bounds.reimer_threshold(size))
            self.assertLessEqual(size ** size, 1 << (2 * incidence))
            # and the two forms agree
            self.assertEqual(
                incidence >= bounds.reimer_threshold(size),
                size ** size <= 1 << (2 * incidence),
            )

    def test_clone_power_audit_admits_every_power(self) -> None:
        report = json.loads(CLONE_POWER_AUDIT.read_text())
        self.assertEqual(report["verdict"], "EVERY_POWER_ADMISSIBLE")
        self.assertEqual(report["base"], "n8clone_hi")
        self.assertTrue(report["every_power_degree_equals_cap"])
        self.assertTrue(report["all_integer_powers_admissible"])
        self.assertTrue(report["all_instantiated_admissible"])
        square = next(r for r in report["instantiated"] if r["power"] == 2)
        self.assertEqual(square["size"], 4_900)
        self.assertTrue(square["defect_identity_holds"])
        self.assertEqual(
            Fraction(square["closure_defect"]),
            1 - (1 - Fraction(317, 490)) ** 2,
        )

    def test_growth_constant_improves_on_the_seven_coordinate_base(self) -> None:
        certificate = json.loads(CLONE_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL"
        )
        clone = certificate["asymptotic"]["bases"]["n8clone_hi"]
        self.assertEqual(Fraction(clone["asymptotic_slope"]), Fraction(3, 640))
        self.assertEqual(Fraction(clone["certified_threshold"]), Fraction(-3, 80))
        self.assertEqual(clone["base_dimension"], 8)
        published_slope = Fraction(7, 250) / 7
        self.assertGreater(Fraction(clone["asymptotic_slope"]), published_slope)
        self.assertEqual(
            Fraction(clone["asymptotic_slope"]) / published_slope, Fraction(75, 64)
        )
        consequence = certificate["bases"]["n8clone_hi"]["rational_consequence"]
        upper = Fraction(consequence["a_plus_two_sided_enclosure"][1])
        self.assertLess(upper, Fraction(-3, 80))
        self.assertGreater(
            Fraction(consequence["base_ratio_lower_gt"]),
            Fraction(286_491, 10_000_000) / Fraction(64, 81),
        )


class FourFourClassTests(unittest.TestCase):
    """The S_4 x S_4 class closes block-symmetric coverage at n=8."""

    def test_k4_class_is_exactly_re_ranked_with_no_screen_error(self) -> None:
        report = json.loads(N8_K4_FRONTIER.read_text())
        self.assertEqual(report["float_negative_count"], 73)
        self.assertEqual(report["certified_negative_count"], 73)
        self.assertEqual(report["float_negative_but_exactly_nonnegative"], [])
        self.assertEqual(report["float_nonnegative_but_exactly_negative"], [])

    def test_lowest_defect_anywhere_is_four_ninths(self) -> None:
        base = rational.BASES["n8tiny"]
        rows = base.reconstruct()
        facts = base.exact_facts()
        self.assertEqual(len(rows), 15)
        self.assertTrue(bounds.admissible(rows, 8))
        self.assertEqual(Fraction(facts["closure_defect"]), Fraction(4, 9))
        self.assertEqual(bounds.defect_count(rows), 100)
        self.assertFalse(facts["separating"])
        self.assertEqual(
            facts["duplicate_column_pairs"],
            [[2, 3], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]],
        )
        # improves every previously certified defect
        for previous in (Fraction(1336, 2025), Fraction(139, 245), Fraction(1144, 1875)):
            self.assertLess(Fraction(4, 9), previous)
        # Corollary 12 requires a dominant set, and it has one: [8] itself.
        self.assertEqual(max(bin(row).count("1") for row in rows), 8)
        self.assertIsNone(bounds.union_growth_floor(rows, 8))

    def test_separating_base_improves_the_growth_constant(self) -> None:
        """n8best raises the slope above 1/250 without dropping separation."""
        certificate = json.loads(N8_K4_CERTIFICATE.read_text())
        self.assertEqual(
            certificate["verdict"], "PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL"
        )
        base = certificate["bases"]["n8best"]
        self.assertTrue(base["exact_base"]["separating"])
        self.assertEqual(
            Fraction(base["exact_base"]["closure_defect"]), Fraction(284, 375)
        )
        slope = Fraction(
            certificate["asymptotic"]["bases"]["n8best"]["asymptotic_slope"]
        )
        self.assertEqual(slope, Fraction(177, 40_000))
        self.assertGreater(slope, Fraction(7, 250) / 7)   # beats the n=7 base
        self.assertLess(slope, Fraction(3, 640))          # but not the clone
        ratio = Fraction(base["rational_consequence"]["base_ratio_lower_gt"])
        self.assertEqual(ratio, Fraction(531, 11_360))
        self.assertGreater(ratio, Fraction(2432, 58_575))  # beats n8max
        self.assertGreater(ratio, Fraction(286_491, 10_000_000) / Fraction(64, 81))

    def test_every_registered_base_is_admissible_and_certified(self) -> None:
        """Ten bases; each admissible, each with a frozen target it beats."""
        self.assertEqual(len(rational.BASES), 10)
        for name, base in rational.BASES.items():
            rows = base.reconstruct()
            self.assertTrue(
                bounds.admissible(rows, base.dimension), f"{name} inadmissible"
            )
            self.assertIn(name, rational.RATIONAL_TARGETS)
            self.assertTrue(bounds.size_ceiling_holds(len(rows), base.dimension))
            facts = base.exact_facts()
            self.assertTrue(facts["active"])


class CloneSaturationTests(unittest.TestCase):
    """Proposition 18: cloning is defect-free, lowers A_+, and saturates."""

    def test_cloning_preserves_size_defect_and_admissibility(self) -> None:
        report = json.loads(CLONE_SATURATION_AUDIT.read_text())
        self.assertEqual(report["verdict"], "CLONING_IS_DEFECT_FREE_AND_SATURATES")
        self.assertTrue(report["size_is_invariant"])
        self.assertTrue(report["defect_is_invariant"])
        self.assertTrue(report["admissibility_preserved_throughout"])
        for row in report["ladder"]:
            self.assertEqual(row["closure_defect"], "4/9")
            self.assertEqual(row["family_size"], 15)
            self.assertTrue(row["admissible"])

    def test_the_collapsed_core_is_admissible_but_positive(self) -> None:
        """n8tiny's clones buy negativity, not feasibility."""
        report = json.loads(CLONE_SATURATION_AUDIT.read_text())
        self.assertEqual(report["collapsed_dimension"], 5)
        core = tuple(report["collapsed_rows"])
        self.assertTrue(bounds.admissible(core, 5))
        self.assertEqual(bounds.defect(core), Fraction(4, 9))
        base = report["ladder"][0]
        self.assertEqual(base["clones"], 0)
        self.assertFalse(base["certified_negative"])
        self.assertGreater(Fraction(base["a_plus_upper"]), 0)
        # and the sign flips only once enough clones are added
        self.assertEqual(report["sign_flips_at_clone_count"], 3)

    def test_clones_of_the_other_two_bases_buy_feasibility_instead(self) -> None:
        """n8clone_lo/hi collapse to inadmissible 7-coordinate families."""
        import audit_clone_saturation as saturation
        for name in ("n8clone_lo", "n8clone_hi"):
            rows = rational.BASES[name].reconstruct()
            core, kept = saturation.collapse(rows, 8)
            self.assertEqual(len(kept), 7)
            self.assertEqual(bounds.defect(core), bounds.defect(rows))
            self.assertFalse(bounds.admissible(core, 7))
            # precisely: Reimer fails, the cap does not
            self.assertEqual(bounds.incidence(core), 196)
            self.assertEqual(bounds.reimer_threshold(len(core)), 215)
            self.assertLessEqual(max(bounds.degrees(core, 7)), bounds.cap(len(core)))

    def test_the_cloning_gain_saturates(self) -> None:
        report = json.loads(CLONE_SATURATION_AUDIT.read_text())
        self.assertTrue(report["improvement_is_shrinking"])
        ratios = [Fraction(r) for r in report["delta_ratios_decimal"]]
        self.assertTrue(ratios)
        for ratio in ratios:
            self.assertLess(ratio, 1)
        deltas = [
            Fraction(row["delta_from_previous_decimal"])
            for row in report["ladder"]
            if row["delta_from_previous_decimal"] is not None
        ]
        for delta in deltas:
            self.assertLess(delta, 0)          # every clone helps
        for earlier, later in zip(deltas, deltas[1:]):
            self.assertLess(abs(later), abs(earlier))   # but by less each time


if __name__ == "__main__":
    unittest.main(verbosity=2)
