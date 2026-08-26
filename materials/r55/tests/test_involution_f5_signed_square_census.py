"""Regression tests for the complete fixed-five W/T/Ramsey census."""

from __future__ import annotations

import hashlib
import json
import sys
import unittest
from collections import Counter
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import involution_f5_signed_completion as completion  # noqa: E402
import involution_f5_t_square_census as t_census  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_DATA = _ROOT / "data"
_W_PATH = _DATA / "involution_f5_w_square_census.json"
_SIGNED_PATH = _DATA / "involution_f5_signed_square_census.json"
_RAMSEY_PATH = _DATA / "involution_f5_ramsey_square_census.json"
_PROOF_PATH = _DATA / "involution_f5_w_dfs_certificate_source0.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class TestExactSquareFrontiers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.w = json.loads(_W_PATH.read_text())
        cls.signed = json.loads(_SIGNED_PATH.read_text())
        cls.ramsey = json.loads(_RAMSEY_PATH.read_text())

    def test_artifact_commitments_and_exact_domain_count(self):
        self.assertEqual(_sha256(_W_PATH),
                         "e99963ed0494b90a7cd3cb117113b317132b711e0592867d27353e76e8a11a6a")
        self.assertEqual(_sha256(_SIGNED_PATH),
                         "01a805f5a783d572966785cedb5c31c89b69ff986ed6bdc77bf3a5dc4c3b8150")
        self.assertEqual(_sha256(_RAMSEY_PATH),
                         "6069740d3c5339dc1d734e10f70d07be53f7046ddcbf07c37974fab322f359c1")
        census = self.w["census"]
        self.assertEqual(census["input_support_orbits"], 705)
        self.assertEqual(census["input_labelled_supports"], 6627)
        self.assertEqual(census["total_row_domains"], 48_865_656)
        self.assertEqual(census["w_square_unsatisfiable_support_orbits"], 550)
        self.assertEqual(census["w_square_satisfiable_support_orbits"], 155)

    def test_w_square_survivors_have_exact_d5_weight(self):
        survivors = [record for record in self.w["census"]["candidates"]
                     if record["w_square_satisfiable"]]
        self.assertEqual(len(survivors), 155)
        self.assertEqual(sum(record["orbit_size"] for record in survivors), 1401)
        self.assertEqual(Counter(record["orbit_size"] for record in survivors),
                         {1: 1, 5: 28, 10: 126})

    def test_complete_signed_square_classification(self):
        census = self.signed["census"]
        self.assertEqual(census["w_square_survivor_support_orbits"], 155)
        self.assertEqual(census["signed_square_unsatisfiable_support_orbits"], 47)
        self.assertEqual(census["signed_square_satisfiable_support_orbits"], 108)
        self.assertTrue(self.signed["claim"]["complete_over_all_w_solutions"])
        self.assertTrue(
            self.signed["claim"]["complete_t_existence_for_each_w_support"])
        survivors = [record for record in census["candidates"]
                     if record["signed_square_satisfiable"]]
        rejected = [record for record in census["candidates"]
                    if not record["signed_square_satisfiable"]]
        self.assertEqual(sum(record["orbit_size"] for record in survivors), 986)
        self.assertEqual(sum(record["orbit_size"] for record in rejected), 415)

    def test_every_committed_signed_witness_reconstructs_an_srg(self):
        checked = 0
        for record in self.signed["census"]["candidates"]:
            if not record["signed_square_satisfiable"]:
                continue
            instance = completion.load_instance(record["source_index"])
            states = t_census.relation_states(
                record["w_half_trits"], record["t_half_trits"])
            result = completion.verify_state_assignment(instance, states)
            self.assertTrue(result["srg_45_22_10_11"])
            self.assertFalse(result["ramsey_good"])
            self.assertIsNotNone(result["first_K5"])
            self.assertIsNotNone(result["first_I5"])
            checked += 1
        self.assertEqual(checked, 108)

    def test_all_signed_switching_classes_have_both_obstructions(self):
        census = self.ramsey["census"]
        self.assertEqual(census["signed_square_survivor_support_orbits"], 108)
        self.assertEqual(census["ramsey_good_support_orbits"], 0)
        self.assertEqual(census["ramsey_bad_support_orbits"], 108)
        self.assertEqual(census["signed_completions_tested"], 73_336)
        self.assertEqual(census["signed_completions_with_k5"], 73_336)
        self.assertEqual(census["signed_completions_with_i5"], 73_336)
        for record in census["candidates"]:
            self.assertFalse(record["ramsey_square_satisfiable"])
            self.assertEqual(record["signed_completions_with_k5"],
                             record["signed_completions_tested"])
            self.assertEqual(record["signed_completions_with_i5"],
                             record["signed_completions_tested"])
        self.assertEqual(
            self.ramsey["claim"]["negative_certificate_status"],
            "PENDING_VERIPB_CAKEPB")
        self.assertFalse(self.ramsey["claim"]["general_ramsey_bound_claimed"])


class TestPinnedProofControl(unittest.TestCase):
    def test_source_zero_has_a_real_veripb_cakepb_certificate(self):
        proof = json.loads(_PROOF_PATH.read_text())
        self.assertEqual(proof["source_index"], 0)
        self.assertEqual(proof["checkers"]["veripb"]["verified_conclusion"],
                         "UNSAT")
        self.assertEqual(proof["checkers"]["cakepb"]["verified_conclusion"],
                         "UNSAT")
        self.assertTrue(proof["checkers"]["veripb"]["force_checked_deletion"])
        self.assertFalse(
            proof["artifacts"]["veripb_proof"]["unchecked_deletion_present"])
        self.assertFalse(
            proof["artifacts"]["cakepb_kernel_proof"][
                "unchecked_deletion_present"])
        self.assertFalse(proof["coverage"]["campaign_complete"])
        for key in ("formula", "veripb_proof", "cakepb_kernel_proof"):
            record = proof["artifacts"][key]
            path = _ROOT / record["relative_path"].removeprefix("r55/")
            self.assertEqual(path.stat().st_size, record["gzip_bytes"])
            self.assertEqual(_sha256(path), record["gzip_sha256"])


if __name__ == "__main__":
    unittest.main()
