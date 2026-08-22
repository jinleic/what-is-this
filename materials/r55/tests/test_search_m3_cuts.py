"""Focused tests for exact affine certificates and primal witness validation.

Tests cover the Certificate and PrimalWitness dataclasses, exact affine bound
computation, search-record status derivation, and the aggregate verdict over
the four persisted records. All tests use tiny synthetic states and no mocks
of solver internals. Objective registry is frozen and verifiers reject
unknown IDs, noncanonical rationals, wrong signs, and every constraint
mismatch.
"""

import pathlib
import sys
import unittest
from dataclasses import replace
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from m3_deficiency_cone import State
from search_m3_cuts import (
    ACCEPTED, ACCEPTED_CUT, Certificate, NO_CUT, PrimalWitness, REJECTED,
    ROUTE_IDS, UNRESOLVED, _terminal_status, _unavailable_record,
    _verify_searches, exact_affine_bound, verify_certificate,
    verify_primal_witness, verify_search_record,
)


class CertificateStructureTests(unittest.TestCase):
    """Test Certificate dataclass immutability and field structure."""

    def test_certificate_is_frozen(self):
        """Certificate instances cannot be mutated."""
        cert = Certificate("total_deficiency", "1", "0", "0", "0")
        with self.assertRaises((AttributeError, TypeError)):
            cert.alpha = "2"

    def test_certificate_has_required_fields(self):
        """Certificate must have exactly five fields."""
        cert = Certificate("total_deficiency", "1", "2", "3", "4")
        self.assertEqual(cert.objective_id, "total_deficiency")
        self.assertEqual(cert.alpha, "1")
        self.assertEqual(cert.beta, "2")
        self.assertEqual(cert.gamma, "3")
        self.assertEqual(cert.delta, "4")

    def test_certificate_rationals_are_canonical_reduced_strings(self):
        """Rational coefficients serialize as reduced fraction strings."""
        cert = Certificate("total_deficiency", "1/2", "-3/4", "5/6", "0")
        self.assertEqual(cert.alpha, "1/2")
        self.assertEqual(cert.beta, "-3/4")
        self.assertEqual(cert.gamma, "5/6")
        self.assertEqual(cert.delta, "0")


class PrimalWitnessStructureTests(unittest.TestCase):
    """Test PrimalWitness dataclass immutability and field structure."""

    def test_primal_witness_is_frozen(self):
        """PrimalWitness instances cannot be mutated."""
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises((AttributeError, TypeError)):
            witness.objective_id = "other"

    def test_primal_witness_has_required_fields(self):
        """PrimalWitness must have objective_id and sorted weights."""
        witness = PrimalWitness("degree20_count", ((0, "45"), (1, "0")))
        self.assertEqual(witness.objective_id, "degree20_count")
        self.assertEqual(witness.weights, ((0, "45"), (1, "0")))

    def test_primal_witness_weights_are_sorted_tuples(self):
        """Weights are immutable sorted (state_index, rational_string) tuples."""
        witness = PrimalWitness("degree20_count", ((0, "1/3"), (1, "2/3")))
        self.assertEqual(len(witness.weights), 2)
        self.assertEqual(witness.weights[0], (0, "1/3"))
        self.assertEqual(witness.weights[1], (1, "2/3"))


class CertificateValidationTests(unittest.TestCase):
    """Test verify_certificate enforces constraints and rejects violations."""

    def test_valid_certificate_accepts_within_constraint(self):
        """Certificate satisfying gamma,delta >= 0 and all state inequalities passes."""
        states = [
            State(20, 0, 0, 0, -12, -4, 2),
            State(24, 0, 0, 0, -12, -2, 4),
        ]
        # Construct a valid certificate: alpha=1, beta=0, gamma=0, delta=0
        # objective(state) = state.deficiency (always 0 in our test)
        # constraint: 0 <= 1 + 0*balance + 0*g_lo - 0*g_hi = 1 ✓
        cert = Certificate("total_deficiency", "1", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)
    def test_exact_affine_bound_returns_tight_bound_and_rejects_lower_alpha(self):
        """The exact affine search finds alpha=15 and bound 45*15=675."""
        states = [
            State(20, 0, 0, 10, -1, 0, 0),
            State(21, 0, 0, 20, 1, 0, 0),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        self.assertEqual(verify_certificate(states, cert), Fraction(675))
        tampered = replace(cert, alpha=str(Fraction(cert.alpha) - 1))
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)


    def test_certificate_rejects_unknown_objective_id(self):
        """Verifier rejects objectives not in frozen registry."""
        states = [State(20, 0, 0, 0, -12, -4, 2)]
        cert = Certificate("unknown_objective", "1", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_negative_gamma(self):
        """Gamma must be >= 0; negative gamma is rejected."""
        states = [State(20, 0, 0, 0, 0, -1, 0)]
        cert = Certificate("total_deficiency", "1", "0", "-1", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_negative_delta(self):
        """Delta must be >= 0; negative delta is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 1)]
        cert = Certificate("total_deficiency", "1", "0", "0", "-1")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_tampered_alpha(self):
        """If alpha is too small, some state inequality fails."""
        states = [State(20, 0, 0, 5, 0, 0, 0)]  # deficiency=5
        # Constraint: 5 <= alpha + 0*balance + 0*g_lo - 0*g_hi = alpha
        # alpha must be >= 5; alpha < 5 fails
        cert = Certificate("total_deficiency", "4", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_violation_with_negative_g_lo(self):
        """When g_lo < 0 and gamma > 0, must satisfy: objective <= alpha + gamma*g_lo."""
        states = [State(20, 0, 0, 10, 0, -5, 0)]  # deficiency=10, g_lo=-5
        # Constraint: 10 <= alpha + gamma*(-5)
        # With alpha=100, gamma=20: 100 + 20*(-5) = 0 < 10, fails
        cert = Certificate("total_deficiency", "100", "0", "20", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_violation_with_positive_g_hi(self):
        """When g_hi > 0 and delta > 0, must satisfy: objective <= alpha - delta*g_hi."""
        states = [State(20, 0, 0, 10, 0, 0, 5)]  # deficiency=10, g_hi=5
        # Constraint: 10 <= alpha - delta*5
        # With alpha=100, delta=20: 100 - 20*5 = 0 < 10, fails
        cert = Certificate("total_deficiency", "100", "0", "0", "20")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_noncanonical_rational_alpha(self):
        """Noncanonical rationals (e.g., 2/4 instead of 1/2) are rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        cert = Certificate("total_deficiency", "2/4", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_accepts_fractional_coefficients(self):
        """Fractional coefficients in canonical form are accepted."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        cert = Certificate("total_deficiency", "1/2", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_degree20_count_objective(self):
        """degree20_count objective: 1 if d==20, 0 otherwise."""
        states = [
            State(20, 0, 0, 0, 0, 0, 0),   # degree20_count=1
            State(21, 0, 0, 0, 0, 0, 0),   # degree20_count=0
        ]
        # Upper bound 1 requires: alpha >= 1 (from first state)
        cert = Certificate("degree20_count", "1", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_deficiency_ge8_count_objective(self):
        """deficiency_ge8_count objective: 1 if deficiency>=8, 0 otherwise."""
        states = [
            State(20, 0, 0, 8, 0, 0, 0),   # deficiency_ge8_count=1
            State(20, 0, 0, 7, 0, 0, 0),   # deficiency_ge8_count=0
        ]
        cert = Certificate("deficiency_ge8_count", "1", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_certificate_checked_against_all_states(self):
        """Certificate is validated against every state in the list."""
        states = [
            State(20, 0, 0, 0, 0, 0, 0),
            State(21, 0, 0, 10, 0, 0, 0),  # deficiency=10
            State(22, 0, 0, 5, 0, 0, 0),
        ]
        # Must satisfy all three: alpha >= 0, alpha >= 10, alpha >= 5
        cert = Certificate("total_deficiency", "10", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)


class PrimalWitnessValidationTests(unittest.TestCase):
    """Test verify_primal_witness enforces constraints and rejects violations."""

    def test_valid_single_state_witness(self):
        """A single state with weight 45 satisfies balance=0 and returns objective value."""
        states = [State(20, 0, 0, 0, 0, -1, 1)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        result = verify_primal_witness(states, witness)
        # degree20_count(state 0) = 1 (d==20), total = 45*1 = 45
        self.assertEqual(result, Fraction(45))

    def test_witness_rejects_unknown_objective_id(self):
        """Witness with unknown objective_id is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_negative_weight(self):
        """Weights must be nonnegative; negative weight is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("total_deficiency", ((0, "-1"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_wrong_total_weight(self):
        """Sum of weights must equal 45; wrong total is rejected."""
        states = [
            State(20, 0, 0, 0, 0, 0, 0),
            State(21, 0, 0, 0, 0, 0, 0),
        ]
        witness = PrimalWitness("total_deficiency", ((0, "30"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_nonzero_balance(self):
        """Sum of excess_balance*weight must equal 0; nonzero is rejected."""
        states = [
            State(20, 0, 0, 0, 10, 0, 0),  # excess_balance=10
            State(21, 0, 0, 0, -5, 0, 0),  # excess_balance=-5
        ]
        # weights (0, "30"), (1, "15"): balance = 30*10 + 15*(-5) = 300-75 = 225 ≠ 0
        witness = PrimalWitness("total_deficiency", ((0, "30"), (1, "15")))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_positive_aggregate_g_lo(self):
        """Sum of g_lo*weight must be <= 0; positive is rejected."""
        states = [
            State(20, 0, 0, 0, 0, 1, 0),  # g_lo=1
            State(21, 0, 0, 0, 0, 0, 0),
        ]
        # weights (0, "45"): g_lo_sum = 45*1 = 45 > 0, rejected
        witness = PrimalWitness("total_deficiency", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_negative_aggregate_g_hi(self):
        """Sum of g_hi*weight must be >= 0; negative is rejected."""
        states = [
            State(20, 0, 0, 0, 0, 0, -2),  # g_hi=-2
            State(21, 0, 0, 0, 0, 0, 0),
        ]
        # weights (0, "45"): g_hi_sum = 45*(-2) = -90 < 0, rejected
        witness = PrimalWitness("total_deficiency", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_out_of_range_state_index(self):
        """State index must be in [0, len(states)); out-of-range is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("total_deficiency", ((5, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_duplicate_state_index(self):
        """Each state index may appear at most once; duplicates are rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0), State(21, 0, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("total_deficiency", ((0, "25"), (0, "20")))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_unsorted_indices(self):
        """Weights must be sorted by state index; unsorted is rejected."""
        states = [
            State(20, 0, 0, 0, 0, 0, 0),
            State(21, 0, 0, 0, 0, 0, 0),
            State(22, 0, 0, 0, 0, 0, 0),
        ]
        # Indices (1, 0, 2) are not sorted
        witness = PrimalWitness("total_deficiency", ((1, "20"), (0, "25")))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_rejects_noncanonical_rational_weight(self):
        """An unreduced weight string is rejected before the total check."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        # 90/2 is exactly 45, so the total is right and only the spelling
        # is wrong: the canonicality guard is what must reject this.
        witness = PrimalWitness("total_deficiency", ((0, "90/2"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_witness_computes_objective_value(self):
        """verify_primal_witness returns sum of objective(s)*weight."""
        states = [
            State(20, 0, 0, 5, 0, 0, 0),   # deficiency=5
            State(21, 0, 0, 10, 0, 0, 0),  # deficiency=10
        ]
        # weights (0, "30"), (1, "15"): sum_deficiency = 30*5 + 15*10 = 150+150 = 300
        witness = PrimalWitness("total_deficiency", ((0, "30"), (1, "15")))
        result = verify_primal_witness(states, witness)
        self.assertEqual(result, Fraction(300))

    def test_witness_validates_degree20_count(self):
        """degree20_count objective is 1 only for d==20."""
        states = [
            State(20, 0, 0, 0, 0, 0, 0),   # degree20_count=1
            State(21, 0, 0, 0, 0, 0, 0),   # degree20_count=0
            State(22, 0, 0, 0, 0, 0, 0),   # degree20_count=0
        ]
        # weights (0, "45"): sum = 45*1 = 45
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        result = verify_primal_witness(states, witness)
        self.assertEqual(result, Fraction(45))

    def test_witness_validates_deficiency_ge8_count(self):
        """deficiency_ge8_count objective is 1 only for deficiency>=8."""
        states = [
            State(20, 0, 0, 7, 0, 0, 0),   # deficiency_ge8_count=0
            State(21, 0, 0, 8, 0, 0, 0),   # deficiency_ge8_count=1
            State(22, 0, 0, 9, 0, 0, 0),   # deficiency_ge8_count=1
        ]
        # weights (1, "20"), (2, "25"): sum = 20*1 + 25*1 = 45
        witness = PrimalWitness("deficiency_ge8_count", ((1, "20"), (2, "25")))
        result = verify_primal_witness(states, witness)
        self.assertEqual(result, Fraction(45))


class SearchRecordValidationTests(unittest.TestCase):
    """Test verify_search_record enforces schema and rejects violations."""

    def test_valid_search_record_accepted_status(self):
        """Valid search record with ACCEPTED_EXACT_CUT status is accepted."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": {
                "objective_id": "total_deficiency",
                "alpha": "1",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "45",
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        result = verify_search_record(states, record)
        self.assertIsNotNone(result)
    def test_total_deficiency_accepts_exact_threshold_bound_315(self):
        """The total-deficiency route accepts an exact upper bound of 315."""
        states = [State(21, 0, 0, 7, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": {
                "objective_id": "total_deficiency",
                "alpha": "7",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "315",
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        self.assertIsNotNone(verify_search_record(states, record))

    def test_degree20_exact_bound_one_is_unresolved_not_accepted(self):
        """A count bound exactly one is unresolved, never an accepted cut."""
        states = [State(21, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": {
                "objective_id": "degree20_count",
                "alpha": "1/45",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "1",
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)
        record["route_status"] = "CERTIFICATION_UNRESOLVED"
        self.assertIsNotNone(verify_search_record(states, record))


    def test_valid_search_record_rejected_status(self):
        """Valid search record with REJECTED_BY_EXACT_WITNESS status is accepted."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": {
                "objective_id": "degree20_count",
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": "45",
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        result = verify_search_record(states, record)
        self.assertIsNotNone(result)

    def test_valid_search_record_unresolved_status(self):
        """Valid search record with CERTIFICATION_UNRESOLVED status is accepted."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "CERTIFICATION_UNRESOLVED",
        }
        result = verify_search_record(states, record)
        self.assertIsNotNone(result)

    def test_record_rejects_extra_field(self):
        """Search record with extra 'accepted' field is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "CERTIFICATION_UNRESOLVED",
            "accepted": True,  # extra field
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_missing_required_field(self):
        """Search record missing a required field is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            # missing exact_witness_value and route_status
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_invalid_route_status(self):
        """Search record with invalid route_status is rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "INVALID_STATUS",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_mismatched_certificate_objective(self):
        """Certificate objective_id must match record objective_id."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": {
                "objective_id": "degree20_count",  # mismatch
                "alpha": "1",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "45",
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_mismatched_witness_objective(self):
        """Witness objective_id must match record objective_id."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": {
                "objective_id": "total_deficiency",  # mismatch
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": "45",
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_wrong_exact_upper_bound(self):
        """exact_upper_bound must match 45*alpha from certificate."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": {
                "objective_id": "total_deficiency",
                "alpha": "1",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "50",  # wrong: should be 45*1=45
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_wrong_exact_witness_value(self):
        """exact_witness_value must match computed objective value from witness."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": {
                "objective_id": "degree20_count",
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": "50",  # wrong: should be 45*1=45
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_must_have_certificate_for_accepted_status(self):
        """ACCEPTED_EXACT_CUT requires non-null certificate."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_must_have_witness_for_rejected_status(self):
        """REJECTED_BY_EXACT_WITNESS requires non-null primal_witness."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_inconsistent_null_certificate_and_bound(self):
        """If certificate is not null, exact_upper_bound must not be null."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "total_deficiency",
            "certificate": {
                "objective_id": "total_deficiency",
                "alpha": "1",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": None,  # inconsistent
            "primal_witness": None,
            "exact_witness_value": None,
            "route_status": "ACCEPTED_EXACT_CUT",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_inconsistent_null_witness_and_value(self):
        """If primal_witness is not null, exact_witness_value must not be null."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": {
                "objective_id": "degree20_count",
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": None,  # inconsistent
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_requires_object_form_weight_entries(self):
        """A persisted weight is an object; other spellings are rejected."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        record = {
            "objective_id": "degree20_count",
            "certificate": None,
            "exact_upper_bound": None,
            "primal_witness": {
                "objective_id": "degree20_count",
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": "45",
            "route_status": "REJECTED_BY_EXACT_WITNESS",
        }
        self.assertEqual(
            verify_search_record(states, record),
            "REJECTED_BY_EXACT_WITNESS",
        )
        # A JSON array pair, a Python pair, an extra key, a missing key.
        spellings = (
            [[0, "45"]],
            [(0, "45")],
            [{"state_index": 0, "weight": "45", "w": "1"}],
            [{"state_index": 0}],
        )
        for spelling in spellings:
            record["primal_witness"]["weights"] = spelling
            with self.assertRaises(ValueError):
                verify_search_record(states, record)


class BoundaryAndEdgeCaseTests(unittest.TestCase):
    """Test edge cases: fractional coefficients, multiple states, boundary conditions."""

    def test_fractional_alpha_coefficient(self):
        """Alpha can be fractional in canonical form."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        cert = Certificate("total_deficiency", "1/2", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_fractional_beta_coefficient(self):
        """Beta can be fractional and negative."""
        states = [State(20, 0, 0, 0, 5, 0, 0)]  # excess_balance=5
        # objective <= alpha + beta*5; alpha >= 0, beta = -1/2
        # 0 <= alpha - 5/2, so alpha >= 5/2
        cert = Certificate("total_deficiency", "3", "-1/2", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_zero_coefficients(self):
        """All-zero coefficients (except possibly alpha) are valid."""
        states = [State(20, 0, 0, 0, 0, 0, 0)]
        cert = Certificate("total_deficiency", "0", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_witness_with_partial_state_coverage(self):
        """Witness can use a subset of states if weights sum to 45 and constraints hold."""
        states = [
            State(20, 0, 0, 10, 0, 0, 0),  # deficiency=10
            State(21, 0, 0, 8, 0, 0, 0),   # deficiency=8
            State(22, 0, 0, 12, 0, 0, 0),  # deficiency=12
        ]
        # Use only states 0 and 2: weights (0, "20"), (2, "25")
        # balance: 20*0 + 25*0 = 0 ✓, deficiency_sum = 20*10 + 25*12 = 200+300 = 500
        witness = PrimalWitness("total_deficiency", ((0, "20"), (2, "25")))
        result = verify_primal_witness(states, witness)
        self.assertEqual(result, Fraction(500))

    def test_large_number_of_states(self):
        """Certificate validation scales to many states."""
        states = [State(20 + i, 0, 0, i, 0, 0, 0) for i in range(10)]
        # All deficiencies 0..9; max=9
        cert = Certificate("total_deficiency", "10", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_witness_with_many_states(self):
        """Witness validation scales to many states."""
        states = [State(20, 0, 0, 1, -i, 0, 0) for i in range(5)]
        # excess_balance: 0, -1, -2, -3, -4
        # To get balance=0: weights must satisfy w0*0 + w1*(-1) + w2*(-2) + w3*(-3) + w4*(-4) = 0
        # and sum(weights) = 45
        # Example: w0=45 (and others=0): 0 + 0 + 0 + 0 + 0 = 0 ✓
        witness = PrimalWitness("total_deficiency", ((0, "45"),))
        result = verify_primal_witness(states, witness)
        self.assertEqual(result, Fraction(45))

    def test_certificate_with_mixed_sign_beta(self):
        """Beta can be positive or negative; no sign restriction."""
        states = [State(20, 0, 0, 0, 10, 0, 0)]  # excess_balance=10
        # Positive beta: objective <= alpha + beta*balance (more restrictive for positive balance)
        cert = Certificate("total_deficiency", "100", "5", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)

    def test_certificate_zero_bounds_g_lo_and_g_hi(self):
        """States with g_lo=0 and g_hi=0 don't depend on gamma/delta."""
        states = [State(20, 0, 0, 5, 0, 0, 0)]  # g_lo=0, g_hi=0
        cert = Certificate("total_deficiency", "5", "0", "0", "0")
        result = verify_certificate(states, cert)
        self.assertIsNotNone(result)


class TerminalStatusTests(unittest.TestCase):
    """Test _terminal_status maps route statuses to the campaign verdict."""

    def test_all_rejected_routes_give_no_cut(self):
        """Three rejected routes are the only path to NO_CUT_IN_FROZEN_CONE."""
        self.assertEqual(
            _terminal_status((REJECTED, REJECTED, REJECTED)), (NO_CUT, 0),
        )

    def test_one_accepted_route_takes_precedence(self):
        """An accepted route outranks rejected and unresolved ones."""
        self.assertEqual(
            _terminal_status((REJECTED, ACCEPTED, UNRESOLVED)),
            (ACCEPTED_CUT, 1),
        )
        self.assertEqual(
            _terminal_status((ACCEPTED, ACCEPTED, REJECTED)),
            (ACCEPTED_CUT, 2),
        )

    def test_unresolved_route_blocks_no_cut(self):
        """Anything short of all rejected and none accepted is unresolved."""
        self.assertEqual(
            _terminal_status((REJECTED, REJECTED, UNRESOLVED)),
            (UNRESOLVED, 0),
        )
        self.assertEqual(
            _terminal_status((UNRESOLVED, UNRESOLVED, UNRESOLVED)),
            (UNRESOLVED, 0),
        )


class SearchesVerificationTests(unittest.TestCase):
    """Test _verify_searches and _terminal_status over the four records."""

    # One state that rejects all three frozen routes at weight 45: the
    # total deficiency is 45*8 = 360 > 315, and both counts reach 45,
    # never below the exclusive edge of 1.
    STATES = (State(20, 4, 4, 8, 0, 0, 0),)
    VALUES = {
        "total_deficiency": "360",
        "degree20_count": "45",
        "deficiency_ge8_count": "45",
    }

    def _records(self):
        """Three exactly rejected routes in frozen order, then route four."""
        records = [
            {
                "objective_id": route_id,
                "certificate": None,
                "exact_upper_bound": None,
                "primal_witness": {
                    "objective_id": route_id,
                    "weights": [{"state_index": 0, "weight": "45"}],
                },
                "exact_witness_value": self.VALUES[route_id],
                "route_status": REJECTED,
            }
            for route_id in ROUTE_IDS
        ]
        records.append(_unavailable_record())
        return records

    def test_verify_searches_returns_the_frozen_route_statuses(self):
        """All three routes verify as rejected, so the campaign is NO_CUT."""
        statuses = _verify_searches(self.STATES, self._records())
        self.assertEqual(statuses, (REJECTED, REJECTED, REJECTED))
        self.assertEqual(_terminal_status(statuses), (NO_CUT, 0))

    def test_verify_searches_rejects_permuted_route_order(self):
        """Records must sit in the frozen route order, not just be present."""
        records = self._records()
        records[0], records[1] = records[1], records[0]
        with self.assertRaises(ValueError):
            _verify_searches(self.STATES, records)

    def test_verify_searches_rejects_wrong_record_count(self):
        """Exactly four records: three routes plus the unavailable fourth."""
        records = self._records()
        with self.assertRaises(ValueError):
            _verify_searches(self.STATES, records[:3])
        with self.assertRaises(ValueError):
            _verify_searches(self.STATES, records + [_unavailable_record()])

    def test_verify_searches_rejects_unavailable_record_with_evidence(self):
        """The fourth record is pinned: any evidence or status edit fails."""
        evidence = {
            "certificate": {
                "objective_id": "total_deficiency",
                "alpha": "8",
                "beta": "0",
                "gamma": "0",
                "delta": "0",
            },
            "exact_upper_bound": "360",
            "primal_witness": {
                "objective_id": "total_deficiency",
                "weights": [{"state_index": 0, "weight": "45"}],
            },
            "exact_witness_value": "360",
            "route_status": REJECTED,
        }
        for field, value in evidence.items():
            records = self._records()
            records[-1][field] = value
            with self.assertRaises(ValueError):
                _verify_searches(self.STATES, records)

    def test_verify_searches_rejects_array_form_weight(self):
        """A persisted weight stays an object inside the searches array."""
        for spelling in ([[0, "45"]], [(0, "45")]):
            records = self._records()
            records[0]["primal_witness"]["weights"] = spelling
            with self.assertRaises(ValueError):
                _verify_searches(self.STATES, records)


if __name__ == "__main__":
    unittest.main()
