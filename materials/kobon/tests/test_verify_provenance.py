#!/usr/bin/env python3
"""Regression tests for the n=12 provenance verifier."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
DECIDE = HERE.parent / "n12_decide"
sys.path.insert(0, str(DECIDE))
import verify_provenance  # noqa: E402


class VerifyProvenanceTests(unittest.TestCase):
    def test_rejects_stale_manifest_before_regeneration(self) -> None:
        current = {"schema": "current"}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stale_manifest = root / "manifest.json"
            stale_manifest.write_text(
                json.dumps({"schema": "stale"}, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with (
                mock.patch.object(verify_provenance, "MANIFEST", stale_manifest),
                mock.patch.object(
                    verify_provenance.make_manifest, "build", return_value=current
                ),
                mock.patch.object(verify_provenance.engine, "dump_instances") as dump,
            ):
                with self.assertRaisesRegex(RuntimeError, "manifest.json is stale"):
                    verify_provenance.verify(root / "generated")
                dump.assert_not_called()


if __name__ == "__main__":
    unittest.main()
