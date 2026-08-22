"""Focused tests for the disjoint mixed-tier certificate checker."""

import ast
import contextlib
import copy
import io
import json
import pathlib
import re
import sys
import tempfile
import unittest
from fractions import Fraction

R55_ROOT = pathlib.Path(__file__).resolve().parents[1]
MATH_ROOT = R55_ROOT.parent
SOURCE = R55_ROOT / "src" / "check_mixed_certificate.py"
ARTIFACT = R55_ROOT / "data" / "engstrom_identity.json"
sys.path.insert(0, str(R55_ROOT / "src"))

import check_mixed_certificate as checker  # noqa: E402


class DisjointImportTests(unittest.TestCase):
    def test_source_has_only_the_permitted_project_import(self):
        source = SOURCE.read_text(encoding="utf-8")
        forbidden = (
            "mixed_deficiency_cone",
            "m4_deficiency_cone",
            "m3_deficiency_cone",
            "search_mixed_cuts",
            "verify_mixed_independent",
            "mixed_subgraph_identities",
            "m4_subgraph_identities",
            "subgraph_identities",
        )
        import_grep = re.compile(
            r"(?m)^\s*(?:from|import)\s+(?:"
            + "|".join(re.escape(name) for name in forbidden)
            + r")\b"
        )
        self.assertIsNone(import_grep.search(source))

        tree = ast.parse(source, filename=str(SOURCE))
        standard_library = {
            "argparse", "hashlib", "itertools", "json", "multiprocessing",
            "os", "sys", "fractions", "pathlib",
        }
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    self.assertIn(alias.name.split(".", 1)[0], standard_library)
            elif isinstance(node, ast.ImportFrom):
                if node.module == "check_ramsey":
                    self.assertEqual(
                        {alias.name for alias in node.names},
                        {"parse_graph6_line", "popcount"},
                    )
                else:
                    self.assertIn(node.module.split(".", 1)[0], standard_library)


class MotifKernelTests(unittest.TestCase):
    def test_all_labeled_graphs_through_order_six_match_brute_force(self):
        self.assertEqual(checker.cross_validate_motif_kernel(6), 33_867)


class ArtifactVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(ARTIFACT.read_text(encoding="ascii"))

    def _temporary_document(self, document):
        directory = tempfile.TemporaryDirectory()
        path = pathlib.Path(directory.name) / "engstrom_identity.json"
        path.write_text(
            json.dumps(document, sort_keys=True, ensure_ascii=True),
            encoding="ascii",
        )
        return directory, path

    def _assert_rejected(self, mutate, token, *, input_hashes=False):
        candidate = copy.deepcopy(self.document)
        mutate(candidate)
        directory, path = self._temporary_document(candidate)
        try:
            loaded = checker.load_analysis(path)
            with self.assertRaises(checker.CheckViolation) as caught:
                checker.verify_analysis(
                    loaded,
                    MATH_ROOT,
                    check_input_hashes=input_hashes,
                )
        finally:
            directory.cleanup()
        self.assertIn(token, str(caught.exception))

    def test_untouched_artifact_passes(self):
        candidate = copy.deepcopy(self.document)
        directory, path = self._temporary_document(candidate)
        try:
            loaded = checker.load_analysis(path)
            self.assertEqual(
                checker.verify_analysis(loaded, MATH_ROOT),
                checker.DISPOSITION_NO_CUT,
            )
        finally:
            directory.cleanup()

    def test_no_input_hashes_cli_prints_unproved_provenance_notice(self):
        candidate = copy.deepcopy(self.document)
        directory, path = self._temporary_document(candidate)
        stdout = io.StringIO()
        stderr = io.StringIO()
        try:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                code = checker.main([
                    "--root", str(MATH_ROOT),
                    "--no-input-hashes",
                    str(path),
                ])
        finally:
            directory.cleanup()
        self.assertEqual(code, 0, stderr.getvalue())
        self.assertIn("MIXED PROVENANCE UNPROVED", stdout.getvalue())
        self.assertTrue(
            stdout.getvalue().rstrip().endswith(
                "MIXED EVIDENCE VERIFIED: MIXED_NO_CUT_IN_FROZEN_BASIS"
            )
        )

    def test_rejects_changed_state_f_lo(self):
        def mutate(document):
            document["states"][0]["f_lo"] += 1

        self._assert_rejected(mutate, "f_lo")

    def test_rejects_lowered_alpha(self):
        def mutate(document):
            document["searches"][0]["certificate"]["alpha"] = "7"

        self._assert_rejected(mutate, "alpha")

    def test_rejects_padded_alpha(self):
        def mutate(document):
            document["searches"][0]["certificate"]["alpha"] = "9"

        self._assert_rejected(mutate, "least admissible")

    def test_rejects_changed_primal_weight(self):
        def mutate(document):
            pair = document["searches"][0]["primal_witness"]["weights"][0]
            pair["weight"] = str(Fraction(pair["weight"]) + 1)

        self._assert_rejected(mutate, "weights")

    def test_rejects_changed_trust_root_url(self):
        def mutate(document):
            document["trust_roots"][0]["url"] += "#drift"

        self._assert_rejected(mutate, "trust_roots[0]")

    def test_rejects_changed_g_endpoint_without_consulting_frozen_m4_states(self):
        def mutate(document):
            document["states"][0]["g_lo"] += 1

        self._assert_rejected(mutate, "g_lo")

    def test_rejects_changed_histogram_count(self):
        def mutate(document):
            histogram = document["catalog_motif_histograms"][0]["histograms"]["q"]
            histogram[0]["count"] += 1

        self._assert_rejected(mutate, "counts sum")

    def test_rejects_changed_n49_field(self):
        def mutate(document):
            document["n49_replay"]["row_values"][0] += 1

        self._assert_rejected(mutate, "n49_replay")

    def test_rejects_changed_route_status(self):
        def mutate(document):
            document["searches"][0]["route_status"] = checker.ACCEPTED

        self._assert_rejected(mutate, "route_status")

    def test_rejects_deleted_state(self):
        def mutate(document):
            document["states"].pop()

        self._assert_rejected(mutate, "states")

    def test_rejects_changed_schema_version(self):
        def mutate(document):
            document["schema_version"] = 2

        self._assert_rejected(mutate, "schema_version")

    def test_rejects_swapped_f_interval(self):
        def mutate(document):
            state = next(
                state for state in document["states"]
                if state["f_lo"] < state["f_hi"]
            )
            state["f_lo"], state["f_hi"] = state["f_hi"], state["f_lo"]

        self._assert_rejected(mutate, "inverted")

    def test_rejects_changed_recorded_input_sha256(self):
        def mutate(document):
            digest = document["inputs"][0]["sha256"]
            document["inputs"][0]["sha256"] = (
                ("1" if digest[0] == "0" else "0") + digest[1:]
            )

        self._assert_rejected(mutate, "sha256", input_hashes=True)


if __name__ == "__main__":
    unittest.main()
