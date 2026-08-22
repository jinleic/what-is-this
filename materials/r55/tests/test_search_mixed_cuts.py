#!/usr/bin/env python3
"""Frozen eight-coefficient cut-search contracts over the mixed n=45 cone.

Every observable contract below is specified in
notes/engstrom_identity_mixed_plan_2026-08-18.md, Task 3.
"""

import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mixed_deficiency_cone import State  # noqa: E402
from search_mixed_cuts import (  # noqa: E402  (RED before Task 3)
    Certificate,
    PrimalWitness,
    _exact_alpha,
    exact_affine_bound,
    verify_certificate,
    verify_primal_witness,
)


# d=20, deficiency/excess_balance zero: mirror the m4 synthetic convention
# so a single support state carries aggregate balance 45*0 = 0.
def S(g_lo, g_hi, h_lo, h_hi, f_lo, f_hi, deficiency=0):
    return State(
        20, 0, 0, deficiency, 0, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi)


class SearchMixedCutTests(unittest.TestCase):
    def test_rationalized_certificate_checks_every_state(self):
        states = [S(-4, 2, -8, 1, -3, 5), S(-2, 4, -6, 3, -7, 1)]
        certificate = exact_affine_bound(states, "total_deficiency")
        self.assertIsInstance(verify_certificate(states, certificate), object)

    def test_tampered_alpha_is_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5), S(-2, 4, -6, 3, -7, 1)]
        certificate = exact_affine_bound(states, "total_deficiency")
        fields = certificate.__dict__
        tampered = Certificate(
            fields["objective_id"],
            str(__import__("fractions").Fraction(fields["alpha"]) - 1),
            fields["beta"], fields["gamma"], fields["delta"],
            fields["epsilon"], fields["zeta"], fields["eta"], fields["theta"],
        )
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_tampered_alpha_plus_is_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5, deficiency=1), S(-2, 4, -6, 3, -7, 1)]
        certificate = exact_affine_bound(states, "total_deficiency")
        fields = certificate.__dict__
        tampered = Certificate(
            fields["objective_id"],
            str(__import__("fractions").Fraction(fields["alpha"]) + 1),
            fields["beta"], fields["gamma"], fields["delta"],
            fields["epsilon"], fields["zeta"], fields["eta"], fields["theta"],
        )
        with self.assertRaises(ValueError):
            verify_certificate(states, tampered)

    def test_negative_eta_is_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5)]
        certificate = Certificate(
            "total_deficiency", "0", "0", "0", "0", "0", "0", "-1", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, certificate)

    def test_negative_theta_is_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5)]
        certificate = Certificate(
            "total_deficiency", "0", "0", "0", "0", "0", "0", "0", "-1")
        with self.assertRaises(ValueError):
            verify_certificate(states, certificate)

    def test_negative_gamma_delta_epsilon_zeta_still_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5)]
        for name, index in (("gamma", 3), ("delta", 4), ("epsilon", 5),
                            ("zeta", 6)):
            with self.subTest(sign_gate=name):
                parts = ["total_deficiency", "0", "0", "0", "0", "0", "0",
                         "0", "0"]
                parts[index] = "-1"
                with self.assertRaises(ValueError):
                    verify_certificate(states, Certificate(*parts))

    def test_wrong_certificate_arity_is_rejected(self):
        states = [S(-4, 2, -8, 1, -3, 5)]
        for name, index in (("gamma", 3), ("eta", 7), ("theta", 8)):
            with self.subTest(sign_gate=name):
                parts = ["total_deficiency", "0", "0", "0", "0", "0", "0",
                         "0", "0"]
                parts[index] = "-1/2"
                with self.assertRaises(ValueError):
                    verify_certificate(states, Certificate(*parts))
        with self.assertRaises(ValueError):
            verify_certificate(
                states, {"objective_id": "total_deficiency", "alpha": "0"})
        with self.assertRaises(ValueError):
            verify_certificate(
                states, Certificate("total_deficiency", "x", "0", "0", "0",
                                    "0", "0", "0", "0"))
        with self.assertRaises(ValueError):
            verify_certificate(
                states, Certificate("total_deficiency", "1/2", "0", "0", "0",
                                    "0", "0", "0", "2/4"))

    def test_exact_primal_rejection_witness_with_mixed_rows(self):
        states = [S(-1, 1, -2, 2, -1, 1, deficiency=0)]
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        self.assertEqual(verify_primal_witness(states, witness), 45)

    def test_primal_witness_violating_f_row_is_rejected(self):
        states = [S(-1, 1, -2, 2, 1, 2, deficiency=0)]  # f_lo > 0
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_primal_witness_violating_f_hi_row_is_rejected(self):
        states = [S(-1, 1, -2, 2, -2, -1, deficiency=0)]  # f_hi < 0
        witness = PrimalWitness("degree20_count", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)

    def test_wrong_objective_id_is_rejected(self):
        states = [S(-1, 1, -2, 2, -1, 1)]
        witness = PrimalWitness("not_registered", ((0, "45"),))
        with self.assertRaises(ValueError):
            verify_primal_witness(states, witness)
        certificate = Certificate(
            "not_registered", "0", "0", "0", "0", "0", "0", "0", "0")
        with self.assertRaises(ValueError):
            verify_certificate(states, certificate)

    def test_perturbed_f_window_moves_alpha(self):
        # A cone whose verified certificate carries eta > 0: the binding
        # state's positive f_lo lets eta lower alpha, the anchor's negative
        # f_lo pins eta from above, and every other interval coefficient
        # strictly raises the requirement, so the optimum uses eta alone.
        binding = State(20, 0, 0, 0, 0, -1, 1, -1, 1, 4, 6)
        anchor = State(20, 0, 0, -2, 0, -1, 1, -1, 1, -2, 3)
        states = [binding, anchor]
        certificate = exact_affine_bound(states, "total_deficiency")
        eta = __import__("fractions").Fraction(certificate.eta)
        self.assertGreater(eta, 0)
        tightened = [
            State(20, 0, 0, 0, 0, -1, 1, -1, 1, 3, 6),
            anchor,
        ]
        # Exact moved bound: tightening f_lo by one raises least alpha by eta
        moved = _exact_alpha(
            tightened, "total_deficiency",
            __import__("fractions").Fraction(certificate.beta),
            __import__("fractions").Fraction(certificate.gamma),
            __import__("fractions").Fraction(certificate.delta),
            __import__("fractions").Fraction(certificate.epsilon),
            __import__("fractions").Fraction(certificate.zeta),
            eta,
            __import__("fractions").Fraction(certificate.theta),
        )
        self.assertEqual(
            moved, __import__("fractions").Fraction(certificate.alpha) + eta)
        with self.assertRaises(ValueError):
            verify_certificate(tightened, certificate)

    def test_m4_witnesses_remain_valid_mixed_witnesses(self):
        # The frozen m4 exact primal witnesses against the REAL artifact
        # states satisfy sum f_lo*w <= 0 <= sum f_hi*w: the mixed row does
        # not reopen any rejected route.  This is the decisive carryover.
        import json
        from fractions import Fraction
        m4 = json.loads(
            (ROOT / "data" / "higher_identity_m4.json").read_text(
                encoding="ascii"))
        eng = json.loads(
            (ROOT / "data" / "engstrom_identity.json").read_text(
                encoding="ascii"))
        m4_keys = [(s["d"], s["a"], s["b"]) for s in m4["states"]]
        eng_index = {
            (s["d"], s["a"], s["b"]): s for s in eng["states"]
        }
        self.assertEqual(sorted(eng_index), sorted(set(m4_keys)))
        count = 0
        for record in m4["searches"]:
            if not record["objective_id"].endswith("count") and (
                    record["objective_id"] != "total_deficiency"):
                continue
            witness = record["primal_witness"]
            if witness is None:
                continue
            low = Fraction(0)
            high = Fraction(0)
            for pair in witness["weights"]:
                state = eng_index[m4_keys[pair["state_index"]]]
                weight = Fraction(pair["weight"])
                low += weight * state["f_lo"]
                high += weight * state["f_hi"]
            self.assertLessEqual(low, 0, record["objective_id"])
            self.assertGreaterEqual(high, 0, record["objective_id"])
            count += 1
        self.assertEqual(count, 3)


if __name__ == "__main__":
    unittest.main()
