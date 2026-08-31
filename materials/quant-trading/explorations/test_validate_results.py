#!/usr/bin/env python3
"""Regression tests for exploration artifact return applicability."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from validate_results import EXPECTED_RETURNS_HEADER, validate_directory


class ReturnsApplicabilityTests(unittest.TestCase):
    def make_artifact(self, returns_applicable: object = ...) -> Path:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name) / "screen"
        directory.mkdir()
        strategy = {
            "schema_version": 1,
            "id": "screen",
            "status": "falsified",
            "holdout_accessed": False,
            "trial_count": 0,
        }
        if returns_applicable is not ...:
            strategy["returns_applicable"] = returns_applicable
        (directory / "strategy.json").write_text(json.dumps(strategy), encoding="utf-8")
        (directory / "trials.json").write_text('{"trial_count": 0}\n', encoding="utf-8")
        (directory / "run.py").write_text("", encoding="utf-8")
        return directory

    def test_explicit_non_return_screen_may_be_falsified_without_returns(self) -> None:
        result = validate_directory(self.make_artifact(False), [])
        self.assertEqual(result["errors"], [])
        self.assertEqual(result["return_rows"], 0)

    def test_falsified_artifact_requires_returns_by_default(self) -> None:
        result = validate_directory(self.make_artifact(), [])
        self.assertIn("empirical result is missing returns.csv", result["errors"])

    def test_returns_applicable_must_be_boolean(self) -> None:
        result = validate_directory(self.make_artifact("false"), [])
        self.assertIn("returns_applicable must be boolean when present", result["errors"])

    def test_non_return_screen_rejects_returns_file(self) -> None:
        directory = self.make_artifact(False)
        row = ["2026-08-30", "0", "0", "0", "0", "0"]
        content = ",".join(EXPECTED_RETURNS_HEADER) + "\n" + ",".join(row) + "\n"
        (directory / "returns.csv").write_text(content, encoding="utf-8")
        result = validate_directory(directory, [])
        self.assertIn("returns_applicable=false but returns.csv exists", result["errors"])


if __name__ == "__main__":
    unittest.main()
