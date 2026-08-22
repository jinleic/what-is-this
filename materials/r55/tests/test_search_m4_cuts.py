"""Task 3 Step 1 (TDD RED): exact six-coefficient certificates, m=4 rows.

Mirrors the hardened m=3 search suite upgraded to six coefficients
(alpha, beta, gamma, delta, epsilon, zeta) and the extra m=4 row: every
state carries h_lo/h_hi on top of the m=3 g endpoints, certificates gain
+epsilon*h_lo - zeta*h_hi, and witnesses add aggregate h_lo <= 0 <= h_hi.
The plan Task 3 Step 1 sketch appears verbatim in ExactCutTests. All
states are tiny synthetic nine-field cells; no catalogs. Red until
math/r55/src/search_m4_cuts.py exists (State is re-exported there).
"""

import pathlib
import sys
import unittest
from dataclasses import replace
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from search_m4_cuts import (  # expected to fail before Task 3 implementation
    ACCEPTED, ACCEPTED_CUT, Certificate, FALLBACK_SUPPORT, MAX_SUPPORT,
    NO_CUT, PrimalWitness, REJECTED, ROUTE_IDS, State, UNRESOLVED,
    _terminal_status, _unavailable_record, _verify_searches,
    exact_affine_bound, verify_certificate, verify_primal_witness,
    verify_search_record,
)

Z = State(20, 0, 0, 0, 0, 0, 0, 0, 0)  # cone-zero degree-20 cell


def _cert(oid, alpha, **kw):
    """Persisted certificate entry: seven keys, zero-sign defaults."""
    entry = {"objective_id": oid, "alpha": alpha, "beta": "0", "gamma": "0",
             "delta": "0", "epsilon": "0", "zeta": "0"}
    entry.update(kw)
    return entry


def _wit(oid, *weights):
    """Persisted witness entry; weights spelled as exact-key objects."""
    return {"objective_id": oid, "weights": [
        {"state_index": index, "weight": weight} for index, weight in weights]}


def _rec(oid="total_deficiency", certificate=None, exact_upper_bound=None,
         primal_witness=None, exact_witness_value=None,
         route_status=UNRESOLVED, **extra):
    """One search record; the overrides spell every stored field."""
    record = {"objective_id": oid, "certificate": certificate,
              "exact_upper_bound": exact_upper_bound,
              "primal_witness": primal_witness,
              "exact_witness_value": exact_witness_value,
              "route_status": route_status}
    record.update(extra)
    return record


class CertificateStructureTests(unittest.TestCase):
    """Certificate dataclass immutability and field structure."""

    def test_certificate_is_frozen(self):
        cert = Certificate("total_deficiency", "1", "0", "0", "0", "0", "0")
        with self.assertRaises((AttributeError, TypeError)):
            cert.alpha = "2"

    def test_certificate_has_required_fields(self):
        """objective_id plus exactly six coefficients, in order."""
        cert = Certificate("total_deficiency", "1", "2", "3", "4", "5", "6")
        self.assertEqual(cert.objective_id, "total_deficiency")
        self.assertEqual((cert.alpha, cert.beta, cert.gamma, cert.delta,
                          cert.epsilon, cert.zeta), ("1", "2", "3", "4", "5", "6"))

    def test_certificate_rationals_are_canonical_reduced_strings(self):
        """Rational coefficients serialize as reduced fraction strings."""
        cert = Certificate("total_deficiency", "1/2", "-3/4", "5/6", "0",
                           "7/8", "-9/10")
        self.assertEqual((cert.alpha, cert.beta, cert.gamma, cert.delta,
                          cert.epsilon, cert.zeta),
                         ("1/2", "-3/4", "5/6", "0", "7/8", "-9/10"))


class PrimalWitnessStructureTests(unittest.TestCase):
    """PrimalWitness dataclass immutability and field structure."""

    def test_primal_witness_is_frozen(self):
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises((AttributeError, TypeError)):
            witness.objective_id = "other"

    def test_primal_witness_has_required_fields(self):
        witness = PrimalWitness("degree20_count", ((0, "45"), (1, "0")))
        self.assertEqual(witness.objective_id, "degree20_count")
        self.assertEqual(witness.weights, ((0, "45"), (1, "0")))

    def test_primal_witness_weights_are_sorted_tuples(self):
        """Immutable sorted (state_index, rational_string) tuples."""
        witness = PrimalWitness("degree20_count", ((0, "1/3"), (1, "2/3")))
        self.assertEqual(witness.weights[0], (0, "1/3"))
        self.assertEqual(witness.weights[1], (1, "2/3"))


class ExactCutTests(unittest.TestCase):
    """The plan Task 3 Step 1 sketch, verbatim."""

    def test_rationalized_certificate_checks_every_state(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2, -8, 1),
            State(24, 0, 0, 0, -12, -2, 4, -6, 3),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        verify_certificate(states, cert)

    def test_tampered_alpha_is_rejected(self):
        states = [
            State(20, 0, 0, 0, -12, -4, 2, -8, 1),
            State(24, 0, 0, 0, -12, -2, 4, -6, 3),
        ]
        cert = exact_affine_bound(states, "total_deficiency")
        tampered = replace(cert, alpha=str(Fraction(cert.alpha) - 1))
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_negative_epsilon_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -1, 1)]
        cert = exact_affine_bound(states, "deficiency_ge8_count")
        tampered = replace(cert, epsilon="-1")
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_exact_primal_rejection_witness_with_m4_rows(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -2, 2)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(45))

    def test_primal_witness_violating_m4_row_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, 1, 2)]  # h_lo > 0
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_wrong_objective_id_is_rejected(self):
        states = [State(20, 0, 0, 0, 0, -1, 1, -1, 1)]
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)


class CertificateValidationTests(unittest.TestCase):
    """verify_certificate: signs, canonicality, every state, tight LP."""

    def test_valid_certificate_accepts_within_constraint(self):
        """Sign coefficients >= 0 and slack at every state passes."""
        states = [State(20, 0, 0, 0, -12, -4, 2, -8, 1),
                  State(24, 0, 0, 0, -12, -2, 4, -6, 3)]
        # objective = deficiency = 0 <= 1 at both states
        cert = Certificate("total_deficiency", "1", "0", "0", "0", "0", "0")
        self.assertIsNotNone(verify_certificate(states, cert))

    def test_certificate_rejects_unknown_objective_id(self):
        cert = Certificate("unknown_objective", "1", "0", "0", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate([Z], cert)

    def test_certificate_rejects_negative_sign_coefficients(self):
        """gamma, delta, epsilon, zeta each must be >= 0 in isolation."""
        for index in range(4):
            signs = ["0", "0", "0", "0"]
            signs[index] = "-1"
            cert = Certificate("total_deficiency", "1", "0", *signs)
            with self.assertRaises(ValueError, msg=repr(cert)):
                verify_certificate([Z], cert)

    def test_certificate_rejects_tampered_alpha(self):
        """If alpha is too small, some state inequality fails."""
        states = [State(20, 0, 0, 5, 0, 0, 0, 0, 0)]  # deficiency=5
        cert = Certificate("total_deficiency", "4", "0", "0", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_certificate_rejects_per_row_shrunken_right_sides(self):
        """Each interval row can sink the right side below the objective
        (100 -> 0 < 10): +gamma*g_lo, -delta*g_hi, +epsilon*h_lo, -zeta*h_hi."""
        rows = ((State(20, 0, 0, 10, 0, -5, 0, 0, 0), ("20", "0", "0", "0")),
                (State(20, 0, 0, 10, 0, 0, 5, 0, 0), ("0", "20", "0", "0")),
                (State(20, 0, 0, 10, 0, 0, 0, -5, 0), ("0", "0", "20", "0")),
                (State(20, 0, 0, 10, 0, 0, 0, 0, 5), ("0", "0", "0", "20")))
        for state, signs in rows:
            cert = Certificate("total_deficiency", "100", "0", *signs)
            with self.assertRaises(ValueError, msg=str(state)):
                verify_certificate((state,), cert)

    def test_positive_h_lo_slack_is_required_not_decorative(self):
        """h_lo > 0 with epsilon carrying the bound: 10 <= 4 + 3*2 is
        exactly tight; collapsing epsilon tightens the right side to 4."""
        states = [State(20, 0, 0, 10, 0, 0, 0, 2, 4)]  # h_lo=2, h_hi=4
        tight = Certificate("total_deficiency", "4", "0", "0", "0", "3", "0")
        self.assertEqual(verify_certificate(states, tight), Fraction(180))
        with self.assertRaises(ValueError):
            verify_certificate(states, replace(tight, epsilon="0"))

    def test_mixed_g_h_straddle_certificate_is_tight_and_tamper_evident(self):
        """One state straddling both rows uses all four sign terms:
        10 <= 35/2 + 1*(-2) - 1*3 + (1/2)*(-4) - (1/10)*5 = 10 exactly."""
        states = [State(20, 0, 0, 10, 0, -2, 3, -4, 5)]
        cert = Certificate("total_deficiency", "35/2", "0", "1", "1", "1/2", "1/10")
        self.assertEqual(verify_certificate(states, cert), Fraction(1575, 2))
        with self.assertRaises(ValueError):
            verify_certificate(states, replace(cert, alpha="17"))

    def test_certificate_rejects_noncanonical_rational_alpha(self):
        """Noncanonical rationals (2/4 instead of 1/2) are rejected."""
        cert = Certificate("total_deficiency", "2/4", "0", "0", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate([Z], cert)

    def test_count_objectives_pin_degree_and_deficiency_windows(self):
        states = [State(20, 0, 0, 0, 0, 0, 0, 0, 0),   # degree20=1, ge8=0
                  State(21, 0, 0, 8, 0, 0, 0, 0, 0)]   # degree20=0, ge8=1
        for oid in ("degree20_count", "deficiency_ge8_count"):
            cert = Certificate(oid, "1", "0", "0", "0", "0", "0")
            self.assertIsNotNone(verify_certificate(states, cert))

    def test_certificate_checked_against_all_states(self):
        """The certificate is validated against every state in the list."""
        states = [State(20, 0, 0, 0, 0, 0, 0, 0, 0),
                  State(21, 0, 0, 10, 0, 0, 0, 0, 0),  # deficiency=10
                  State(22, 0, 0, 5, 0, 0, 0, 0, 0)]
        cert = Certificate("total_deficiency", "10", "0", "0", "0", "0", "0")
        self.assertIsNotNone(verify_certificate(states, cert))

    def test_every_state_check_includes_the_h_row(self):
        """A later state's h endpoint must fail the sweep:
        state 0: 4 <= 3 + 1 = 4; state 1: 4 <= 3 - 8*1 = -5."""
        states = [State(20, 0, 0, 4, 0, 0, 0, 1, 2),
                  State(24, 0, 0, 4, 0, 0, 0, -8, 0)]
        cert = Certificate("total_deficiency", "3", "0", "0", "0", "1", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, cert)

    def test_exact_affine_bound_tight_and_beats_trivial_candidate(self):
        """Hand-computed h-row LP: alpha >= max(10+2e, 20-2e) minimized at
        epsilon=5/2 gives alpha=15, bound 675; the trivial (all-zero)
        candidate has bound 45*20=900 and must not shadow the better one."""
        states = [State(20, 0, 0, 10, 0, 0, 0, -2, 0),  # >= 10 + 2*epsilon
                  State(21, 0, 0, 20, 0, 0, 0, 2, 2),   # >= 20 - 2*epsilon
                  State(22, 0, 0, 12, 0, -4, 1, 0, 0)]  # >= 12 always
        cert = exact_affine_bound(states, "total_deficiency")
        self.assertEqual((cert.alpha, cert.epsilon), ("15", "5/2"))
        bound = verify_certificate(states, cert)
        self.assertEqual(bound, Fraction(675))
        self.assertLess(bound, Fraction(900))
        with self.assertRaises(ValueError):
            verify_certificate(states, replace(cert, alpha="14"))


class PrimalWitnessValidationTests(unittest.TestCase):
    """verify_primal_witness: every relation, both h aggregates."""

    def test_valid_single_state_witness(self):
        """Weight 45 on a 0-straddling state; value is the count."""
        states = [State(20, 0, 0, 0, 0, -1, 1, -2, 2)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(45))

    def test_witness_rejects_unknown_objective_id(self):
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness([Z], witness)

    def test_witness_rejects_bad_weights(self):
        """Negative weight, wrong total, noncanonical 90/2 == 45."""
        for weights in (((0, "-1"),), ((0, "30"),), ((0, "90/2"),)):
            witness = PrimalWitness("total_deficiency", weights)
            with self.assertRaises(ValueError, msg=repr(weights)):
                verify_primal_witness([Z], witness)

    def test_witness_rejects_each_violated_aggregate(self):
        """Zero balance; aggregates g_lo <= 0, g_hi >= 0, and the two new
        m=4 aggregates h_lo <= 0, h_hi >= 0 — each in isolation."""
        rows = ((State(20, 0, 0, 0, 10, 0, 0, 0, 0), "balance"),
                (State(20, 0, 0, 0, 0, 1, 0, 0, 0), "g_lo > 0"),
                (State(20, 0, 0, 0, 0, 0, -1, 0, 0), "g_hi < 0"),
                (State(20, 0, 0, 0, 0, 0, 0, 1, 0), "h_lo > 0"),
                (State(20, 0, 0, 0, 0, 0, 0, 0, -1), "h_hi < 0"))
        for state, message in rows:
            witness = PrimalWitness("total_deficiency", ((0, "45"),))
            with self.assertRaises(ValueError, msg=message):
                verify_primal_witness((state,), witness)

    def test_witness_accepts_straddling_g_and_h_aggregates(self):
        """Exact weights (225/7, 90/7): balance 0; g = (-495, 225)/7;
        h = (-540, 855)/7 — both h aggregates strictly correct-signed."""
        states = [State(20, 0, 0, 0, -2, -1, 1, -2, 3),
                  State(21, 1, 0, 0, 5, -3, 0, -1, 2)]
        witness = PrimalWitness("total_deficiency", ((0, "225/7"), (1, "90/7")))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(0))

    def test_witness_rejects_bad_index_structure(self):
        """Out-of-range, duplicate, and unsorted indices are rejected."""
        states = [State(d, 0, 0, 0, 0, 0, 0, 0, 0) for d in (20, 21, 22)]
        for weights in (((5, "45"),), ((0, "25"), (0, "20")),
                        ((1, "20"), (0, "25"))):
            witness = PrimalWitness("total_deficiency", weights)
            with self.assertRaises(ValueError, msg=repr(weights)):
                verify_primal_witness(states, witness)

    def test_witness_computes_objective_value(self):
        """The returned value is sum of objective*weight: 30*5+15*10=300."""
        states = [State(20, 0, 0, 5, 0, 0, 0, 0, 0),
                  State(21, 0, 0, 10, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("total_deficiency", ((0, "30"), (1, "15")))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(300))

    def test_witness_validates_count_objectives(self):
        states = [State(20, 0, 0, 7, 0, 0, 0, 0, 0),   # degree20=1, ge8=0
                  State(21, 0, 0, 8, 0, 0, 0, 0, 0)]   # degree20=0, ge8=1
        first = PrimalWitness("degree20_count", ((0, "20"), (1, "25")))
        self.assertEqual(verify_primal_witness(states, first), Fraction(20))
        second = PrimalWitness("deficiency_ge8_count", ((0, "20"), (1, "25")))
        self.assertEqual(verify_primal_witness(states, second), Fraction(25))


class SearchRecordValidationTests(unittest.TestCase):
    """verify_search_record: strict schema, mismatches, stored == derived."""

    def test_valid_search_record_accepted_status(self):
        record = _rec(certificate=_cert("total_deficiency", "1"),
                      exact_upper_bound="45", route_status=ACCEPTED)
        self.assertEqual(verify_search_record([Z], record), ACCEPTED)

    def test_total_deficiency_accepts_exact_threshold_bound_315(self):
        """The total-deficiency route accepts an exact upper bound of 315."""
        states = [State(21, 0, 0, 7, 0, 0, 0, 0, 0)]
        record = _rec(certificate=_cert("total_deficiency", "7"),
                      exact_upper_bound="315", route_status=ACCEPTED)
        self.assertEqual(verify_search_record(states, record), ACCEPTED)

    def test_count_bound_exactly_one_is_unresolved_not_accepted(self):
        """A count bound exactly one is unresolved, never an accepted cut."""
        states = [State(21, 0, 0, 0, 0, 0, 0, 0, 0)]
        for oid in ("degree20_count", "deficiency_ge8_count"):
            record = _rec(oid=oid, certificate=_cert(oid, "1/45"),
                          exact_upper_bound="1", route_status=ACCEPTED)
            with self.assertRaises(ValueError):
                verify_search_record(states, record)
            record["route_status"] = UNRESOLVED
            self.assertEqual(verify_search_record(states, record), UNRESOLVED)

    def test_valid_search_record_rejected_status(self):
        record = _rec(oid="degree20_count",
                      primal_witness=_wit("degree20_count", (0, "45")),
                      exact_witness_value="45", route_status=REJECTED)
        self.assertEqual(verify_search_record([Z], record), REJECTED)

    def test_valid_search_record_unresolved_status(self):
        self.assertEqual(verify_search_record([Z], _rec()), UNRESOLVED)

    def test_record_rejects_extra_field(self):
        """An extra 'accepted' field is rejected."""
        with self.assertRaises(ValueError):
            verify_search_record([Z], _rec(accepted=True))

    def test_record_rejects_missing_required_field(self):
        record = _rec()
        del record["exact_witness_value"], record["route_status"]
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_invalid_route_status(self):
        with self.assertRaises(ValueError):
            verify_search_record([Z], _rec(route_status="INVALID_STATUS"))

    def test_record_rejects_stored_status_not_derived(self):
        """The stored status decides nothing: value 0 does not clear >315,
        so a record cannot be stored as REJECTED."""
        states = [State(21, 0, 0, 0, 0, 0, 0, 0, 0)]  # witness value 0
        record = _rec(primal_witness=_wit("total_deficiency", (0, "45")),
                      exact_witness_value="0", route_status=REJECTED)
        with self.assertRaises(ValueError):
            verify_search_record(states, record)

    def test_record_rejects_mismatched_certificate_objective(self):
        record = _rec(certificate=_cert("degree20_count", "1"),
                      exact_upper_bound="45", route_status=ACCEPTED)
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_mismatched_witness_objective(self):
        record = _rec(oid="degree20_count",
                      primal_witness=_wit("total_deficiency", (0, "45")),
                      exact_witness_value="45", route_status=REJECTED)
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_wrong_exact_upper_bound(self):
        record = _rec(certificate=_cert("total_deficiency", "1"),
                      exact_upper_bound="50",  # must be 45*1=45
                      route_status=ACCEPTED)
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_wrong_exact_witness_value(self):
        record = _rec(oid="degree20_count",
                      primal_witness=_wit("degree20_count", (0, "45")),
                      exact_witness_value="50",  # must be 45*1=45
                      route_status=REJECTED)
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_must_have_certificate_for_accepted_status(self):
        with self.assertRaises(ValueError):
            verify_search_record([Z], _rec(route_status=ACCEPTED))

    def test_record_must_have_witness_for_rejected_status(self):
        record = _rec(oid="degree20_count", route_status=REJECTED)
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_inconsistent_null_certificate_and_bound(self):
        record = _rec(certificate=_cert("total_deficiency", "1"),
                      route_status=ACCEPTED)  # bound stays None
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_rejects_inconsistent_null_witness_and_value(self):
        record = _rec(oid="degree20_count",
                      primal_witness=_wit("degree20_count", (0, "45")),
                      route_status=REJECTED)  # value stays None
        with self.assertRaises(ValueError):
            verify_search_record([Z], record)

    def test_record_requires_object_form_weight_entries(self):
        """Exact-key objects only; arrays, tuples, extra and missing keys
        are all rejected."""
        record = _rec(oid="degree20_count",
                      primal_witness=_wit("degree20_count", (0, "45")),
                      exact_witness_value="45", route_status=REJECTED)
        self.assertEqual(verify_search_record([Z], record), REJECTED)
        for spelling in ([[0, "45"]], [(0, "45")],
                         [{"state_index": 0, "weight": "45", "w": "1"}],
                         [{"state_index": 0}]):
            record["primal_witness"]["weights"] = spelling
            with self.assertRaises(ValueError):
                verify_search_record([Z], record)


class BoundaryAndEdgeCaseTests(unittest.TestCase):
    """Fractional coefficients, tight six-term certificates, scale."""

    def test_fractional_alpha_coefficient(self):
        """Alpha can be fractional in canonical form."""
        cert = Certificate("total_deficiency", "1/2", "0", "0", "0", "0", "0")
        self.assertIsNotNone(verify_certificate([Z], cert))

    def test_all_six_coefficients_can_vary_and_stay_tight(self):
        """0 <= 1 + (-1/5)*5 + 0 + 0 + (1/3)*3 - (1/2)*2 = 0 exactly."""
        states = [State(20, 0, 0, 0, 5, 0, 0, 3, 2)]
        cert = Certificate("total_deficiency", "1", "-1/5", "0", "0", "1/3", "1/2")
        self.assertEqual(verify_certificate(states, cert), Fraction(45))
        with self.assertRaises(ValueError):
            verify_certificate(states, replace(cert, epsilon="0"))

    def test_zero_coefficients(self):
        """All-zero coefficients (except possibly alpha) are valid."""
        cert = Certificate("total_deficiency", "0", "0", "0", "0", "0", "0")
        self.assertIsNotNone(verify_certificate([Z], cert))

    def test_witness_with_partial_state_coverage(self):
        """A strict subset whose weights sum to 45: 20*10 + 25*12 = 500."""
        states = [State(20, 0, 0, 10, 0, 0, 0, 0, 0),
                  State(21, 0, 0, 8, 0, 0, 0, 0, 0),
                  State(22, 0, 0, 12, 0, 0, 0, 0, 0)]
        witness = PrimalWitness("total_deficiency", ((0, "20"), (2, "25")))
        self.assertEqual(verify_primal_witness(states, witness), Fraction(500))

    def test_large_number_of_states(self):
        """Certificate validation scales to many states with h windows."""
        states = [State(20 + i, 0, 0, i, 0, 0, 0, -i, i) for i in range(10)]
        cert = Certificate("total_deficiency", "10", "0", "0", "0", "0", "0")
        self.assertIsNotNone(verify_certificate(states, cert))


class FrozenSearchConstantsTests(unittest.TestCase):
    """The frozen m=4 registry and basis-reconstruction knobs."""

    def test_route_ids_are_the_frozen_m3_routes_in_order(self):
        self.assertEqual(
            ROUTE_IDS,
            ("total_deficiency", "degree20_count", "deficiency_ge8_count"))

    def test_support_knobs_are_the_m4_frozen_values(self):
        """FALLBACK_SUPPORT=24 and MAX_SUPPORT=6 differ from m=3's 32/4."""
        self.assertEqual(FALLBACK_SUPPORT, 24)
        self.assertEqual(MAX_SUPPORT, 6)


class TerminalStatusTests(unittest.TestCase):
    """_terminal_status maps route statuses to the campaign verdict."""

    def test_all_rejected_routes_give_no_cut(self):
        """Three rejected routes are the only path to NO_CUT."""
        self.assertEqual(
            _terminal_status((REJECTED, REJECTED, REJECTED)), (NO_CUT, 0))

    def test_one_accepted_route_takes_precedence(self):
        """An accepted route outranks rejected and unresolved ones."""
        self.assertEqual(_terminal_status((REJECTED, ACCEPTED, UNRESOLVED)),
                         (ACCEPTED_CUT, 1))
        self.assertEqual(_terminal_status((ACCEPTED, ACCEPTED, REJECTED)),
                         (ACCEPTED_CUT, 2))

    def test_unresolved_route_blocks_no_cut(self):
        self.assertEqual(_terminal_status((REJECTED, REJECTED, UNRESOLVED)),
                         (UNRESOLVED, 0))
        self.assertEqual(_terminal_status((UNRESOLVED, UNRESOLVED, UNRESOLVED)),
                         (UNRESOLVED, 0))


class SearchesVerificationTests(unittest.TestCase):
    """_verify_searches and _terminal_status over the four records."""

    # One state that rejects all three frozen routes at weight 45: the
    # total deficiency is 45*8 = 360 > 315, and both counts reach 45,
    # never below the exclusive edge of 1.
    STATES = (State(20, 4, 4, 8, 0, 0, 0, 0, 0),)
    VALUES = {
        "total_deficiency": "360",
        "degree20_count": "45",
        "deficiency_ge8_count": "45",
    }

    def _records(self):
        """Three exactly rejected routes in frozen order, then route four."""
        records = [
            _rec(oid=route_id, primal_witness=_wit(route_id, (0, "45")),
                 exact_witness_value=self.VALUES[route_id],
                 route_status=REJECTED)
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
        """The fourth record equals the pinned template exactly: any
        evidence field or status edit fails."""
        evidence = {
            "certificate": _cert("total_deficiency", "8"),
            "exact_upper_bound": "360",
            "primal_witness": _wit("total_deficiency", (0, "45")),
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
