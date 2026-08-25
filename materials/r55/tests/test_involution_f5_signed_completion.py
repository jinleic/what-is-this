"""Tests for the exact fixed-five signed-completion model.

The model is scoped to the 705 Ramsey-support survivors.  Solver output is not
accepted as a negative result unless its proof is checked by the pinned VeriPB
and CakePB toolchain.
"""

from __future__ import annotations

import dataclasses
import itertools
import json
import sys
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import certify_involution_f5_w_square as certificate  # noqa: E402
import involution_f5_signed_completion as completion  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_FRONTIER = _ROOT / "data" / "involution_f5_ramsey_filter.json"
_SUPPORT = _ROOT / "data" / "involution_f5_support_census.json"
_UNFILTERED_CONTROL_W = (
    "101120111011221211012010111121101120210111122021201111121112210102101111"
    "211211201012112111020112011121121012201121112110010220211111200110210112"
    "0100101111101101101110011011111112211122101121"
)


class TestSignedRelationModel(unittest.TestCase):
    def test_four_one_hot_states_recover_p_q_w_t(self):
        self.assertEqual(
            {state: completion.relation_values(state)
             for state in completion.RELATION_STATES},
            {
                "00": (0, 0, 1, 0),
                "01": (0, 1, 0, 1),
                "10": (1, 0, 0, -1),
                "11": (1, 1, -1, 0),
            },
        )

    def test_unordered_pair_has_one_shared_state_variable(self):
        instance = completion.load_instance(0)
        registry = completion.VariableRegistry()
        self.assertEqual(
            completion.state_variable(registry, 2, 7, "01"),
            completion.state_variable(registry, 7, 2, "01"),
        )
        self.assertEqual(instance.source_index, 0)

    def test_w_row_domains_are_exact_for_first_survivor(self):
        instance = completion.load_instance(0)
        domains = completion.enumerate_w_row_domains(instance, 0)
        self.assertEqual(len(domains), 1779)
        self.assertEqual(len(set(domains)), len(domains))
        for domain in domains:
            completion.verify_w_row_domain(instance, 0, domain)


class TestGroupActions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.instances = completion.load_all_instances()

    def test_explicit_d5_expansion_has_6627_labelled_supports(self):
        signatures = []
        for instance in self.instances:
            images = completion.explicit_d5_images(instance)
            self.assertEqual(len(images), instance.orbit_size)
            self.assertEqual(len(images), len(set(images)))
            signatures.extend(images)
        self.assertEqual(len(signatures), 6627)
        self.assertEqual(len(signatures), len(set(signatures)))

    def test_residual_action_is_equivariant_on_all_auxiliary_families(self):
        instance = max(
            self.instances,
            key=lambda item: len(completion.residual_generators(item)),
        )
        generators = completion.residual_generators(instance)
        self.assertTrue(generators)
        registry = completion.build_variable_registry(instance)
        semantic_names = tuple(registry.semantic_names())
        for permutation in generators:
            image_names = tuple(
                completion.permute_semantic_variable(name, permutation)
                for name in semantic_names
            )
            self.assertEqual(set(image_names), set(semantic_names))
            for name, image in zip(semantic_names, image_names, strict=True):
                restored = completion.permute_semantic_variable(
                    image, completion.inverse_permutation(permutation))
                self.assertEqual(restored, name)
    def test_residual_action_preserves_the_full_pb_formula(self):
        instance = completion.load_instance(2)
        formula, registry = completion.build_signed_formula(
            instance, include_ramsey=False)
        baseline = {
            (constraint.terms, constraint.comparator, constraint.rhs)
            for constraint in formula.constraints
        }
        for permutation in completion.residual_generators(instance):
            image = {
                completion.permute_constraint_signature(
                    constraint, registry, permutation)
                for constraint in formula.constraints
            }
            self.assertEqual(image, baseline)




class TestFormulaAndGraphChecks(unittest.TestCase):
    def test_full_formula_contains_both_square_layers(self):
        instance = completion.load_instance(0)
        formula, registry = completion.build_signed_formula(
            instance, include_ramsey=False)
        self.assertGreater(formula.constraint_count, 0)
        self.assertTrue(any(name[0] == "wprod"
                            for name in registry.semantic_names()))
        self.assertTrue(any(name[0] == "tprod"
                            for name in registry.semantic_names()))
        opb = formula.to_opb(registry)
        self.assertTrue(opb.startswith("* #variable="))
        self.assertNotIn("del id", opb)
        self.assertNotIn("delc", opb)
    def test_w_square_stage_omits_every_t_product(self):
        instance = completion.load_instance(0)
        formula, registry = completion.build_signed_formula(
            instance, include_t_square=False, include_ramsey=False)
        self.assertTrue(any(name[0] == "wprod"
                            for name in registry.semantic_names()))
        self.assertFalse(any(name[0] == "tprod"
                             for name in registry.semantic_names()))
        self.assertFalse(any(constraint.label.startswith("t-square:")
                             for constraint in formula.constraints))
    def test_complete_w_row_tables_are_equivariant_extensions(self):
        instance = completion.load_instance(0)
        formula, registry = completion.build_signed_formula(
            instance,
            include_t_square=False,
            include_w_domain_tables=True,
            include_ramsey=False,
        )
        row_variables = {
            name for name in registry.semantic_names() if name[0] == "wrow"
        }
        self.assertEqual(len(row_variables), 62463)
        self.assertTrue(any(constraint.label.startswith("wrow-onehot:")
                            for constraint in formula.constraints))
        for permutation in completion.residual_generators(instance):
            self.assertEqual(
                {completion.permute_semantic_variable(name, permutation)
                 for name in row_variables},
                row_variables,
            )
    def test_unfiltered_control_satisfies_every_w_table_constraint(self):
        instance = dataclasses.replace(
            completion.load_instance(393), restrictions=())
        formula, registry = completion.build_signed_formula(
            instance,
            include_t_square=False,
            include_w_domain_tables=True,
            include_ramsey=False,
        )
        pair_states = {}
        rows = [[0] * 20 for _ in range(20)]
        position = 0
        for left, right in itertools.combinations(range(20), 2):
            value = int(_UNFILTERED_CONTROL_W[position]) - 1
            position += 1
            rows[left][right] = rows[right][left] = value
            pair_states[left, right] = {
                1: "00", 0: "01", -1: "11",
            }[value]
        assignment = {}
        for variable, name in enumerate(
                registry.semantic_names(), start=1):
            if name[0] == "state":
                _, left, right, state = name
                assignment[variable] = int(
                    pair_states[left, right] == state)
            elif name[0] == "wprod":
                _, left, right, centre, left_state, right_state = name
                assignment[variable] = int(
                    pair_states[tuple(sorted((left, centre)))] == left_state
                    and pair_states[
                        tuple(sorted((right, centre)))] == right_state
                )
            elif name[0] == "wrow":
                _, row, code = name
                assignment[variable] = int(
                    completion.decode_w_domain(code) == tuple(rows[row]))
            else:
                self.fail(f"unexpected W-only semantic variable {name}")
        failed = [constraint.label for constraint in formula.constraints
                  if not constraint.is_satisfied(assignment)]
        self.assertEqual(failed, [])







    def test_reconstruction_uses_parallel_and_crossed_bits(self):
        instance = completion.load_instance(0)
        states = {(left, right): "00"
                  for left in range(20) for right in range(left + 1, 20)}
        graph = completion.reconstruct_graph(instance, states)
        x0, x1 = completion.orbit_vertices(0)
        y0, y1 = completion.orbit_vertices(1)
        self.assertFalse(completion.has_edge(graph, x0, y0))
        self.assertFalse(completion.has_edge(graph, x1, y1))
        self.assertFalse(completion.has_edge(graph, x0, y1))
        self.assertFalse(completion.has_edge(graph, x1, y0))
        states[0, 1] = "11"
        graph = completion.reconstruct_graph(instance, states)
        self.assertTrue(completion.has_edge(graph, x0, y0))
        self.assertTrue(completion.has_edge(graph, x1, y1))
        self.assertTrue(completion.has_edge(graph, x0, y1))
        self.assertTrue(completion.has_edge(graph, x1, y0))

    def test_independent_ramsey_check_finds_both_obstruction_types(self):
        complete = tuple(((1 << 45) - 1) ^ (1 << vertex)
                         for vertex in range(45))
        empty = (0,) * 45
        self.assertIsNotNone(completion.first_clique(complete, 5))
        self.assertIsNotNone(completion.first_independent_set(empty, 5))
        self.assertFalse(completion.is_srg_45_22_10_11(complete))
        self.assertFalse(completion.is_srg_45_22_10_11(empty))

    def test_certificate_toolchain_is_commit_pinned(self):
        self.assertEqual(
            completion.VERIPB_COMMIT,
            "6d38dab246af9c321b8f17cb5a187f2fbb9e491d",
        )
        self.assertEqual(
            completion.CAKEPB_COMMIT,
            "a7593ef22de2fc0b47a688f2d4f08e6b742735af",
        )
        self.assertTrue(completion.CHECKED_DELETION_REQUIRED)
        self.assertFalse(completion.PROJECTED_COUNTS_ARE_ORBIT_COUNTS)
    def test_certificate_policy_rejects_every_unchecked_deletion(self):
        scratch = (
            _ROOT.parents[1] / "scratch/r55-involution-f5/certificate-tests"
        )
        scratch.mkdir(parents=True, exist_ok=True)
        unchecked = scratch / "unchecked.pbp"
        checked = scratch / "checked.pbp"
        clean = scratch / "clean.pbp"
        for payload in (
                "del id 1 2\n",
                "del\trange 1 2\n",
                "deld 3\n",
                "d 4\n",
                "wiplvl 0\n"):
            unchecked.write_text(payload)
            with self.assertRaises(certificate.CertificateViolation):
                certificate._assert_checked_deletion_policy(
                    unchecked, allow_checked_deletion=True)
        checked.write_text("delc 1 : : subproof\n")
        clean.write_text("pseudo-Boolean proof version 2.0\n")
        with self.assertRaises(certificate.CertificateViolation):
            certificate._assert_checked_deletion_policy(
                checked, allow_checked_deletion=False)
        certificate._assert_checked_deletion_policy(
            checked, allow_checked_deletion=True)
        certificate._assert_checked_deletion_policy(
            clean, allow_checked_deletion=False)

    def test_certificate_requires_exact_unsat_conclusion_and_tool_hashes(self):
        scratch = (
            _ROOT.parents[1] / "scratch/r55-involution-f5/certificate-tests"
        )
        no_conclusion = scratch / "no-conclusion.pbp"
        unsat = scratch / "unsat.pbp"
        substitute = scratch / "substitute-checker"
        no_conclusion.write_text(
            "pseudo-Boolean proof version 2.0\nconclusion NONE\n")
        unsat.write_text(
            "pseudo-Boolean proof version 2.0\nconclusion UNSAT : 1\n")
        substitute.write_text("not a pinned checker\n")
        with self.assertRaises(certificate.CertificateViolation):
            certificate._require_unsat_conclusion(
                no_conclusion, "no-conclusion control")
        certificate._require_unsat_conclusion(unsat, "UNSAT control")
        with self.assertRaises(certificate.CertificateViolation):
            certificate._require_pinned_file(
                substitute, certificate.PINNED_VERIPB_SHA256, "VeriPB")
        for path, digest, label in (
                (certificate.DEFAULT_VERIPB,
                 certificate.PINNED_VERIPB_SHA256, "VeriPB"),
                (certificate.DEFAULT_CAKEPB,
                 certificate.PINNED_CAKEPB_SHA256, "CakePB"),
                (certificate.DEFAULT_ROUNDINGSAT,
                 certificate.PINNED_ROUNDINGSAT_SHA256, "RoundingSat"),
                (certificate.ROUNDINGSAT_LOGGER_SOURCE,
                 certificate.PINNED_ROUNDINGSAT_LOGGER_SHA256,
                 "deletion-free Logger.cpp")):
            certificate._require_pinned_file(path, digest, label)




class TestPinnedInputs(unittest.TestCase):
    def test_frontier_and_support_inputs_are_canonical(self):
        frontier = json.loads(_FRONTIER.read_text())
        support = json.loads(_SUPPORT.read_text())
        self.assertEqual(frontier["census"]["surviving_d5_orbits"], 705)
        self.assertEqual(support["support_census"]["dihedral_orbits"], 844)
        self.assertEqual(
            completion.FRONTIER_SHA256,
            "9dfc89b79b36041a71c113d1f1d4c4a76478b95a4131fb5b12666c16be681f02",
        )


if __name__ == "__main__":
    unittest.main()
