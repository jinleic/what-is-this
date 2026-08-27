from __future__ import annotations

import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest.mock import patch


SYNC_PATH = Path(__file__).resolve().parents[1] / "tools" / "sync.py"
SPEC = importlib.util.spec_from_file_location("what_is_this_sync", SYNC_PATH)
assert SPEC is not None and SPEC.loader is not None
sync = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sync)


class ProgressSelectionTest(unittest.TestCase):
    def test_current_uc_anchor_beats_legacy_campaign(self) -> None:
        detail = "Verified detail. " * 30
        progress_text = textwrap.dedent(
            f"""\
            # Progress Ledger

            ### UC-GATE-B-INDEPENDENCE (2026-08-26)
            {detail}

            ### BOLD FRONTIER CAMPAIGN IV (2026-08-25)
            {detail}
            """
        )
        project = next(
            project for project in sync.PROJECTS
            if project["slug"] == "uc"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            progress = Path(tmpdir) / "PROGRESS.md"
            progress.write_text(progress_text)
            with patch.object(sync, "PROGRESS", progress):
                header, _ = sync.newest_progress_block(project)

        self.assertEqual(
            header,
            "### UC-GATE-B-INDEPENDENCE (2026-08-26)",
        )

    def test_renamed_ising_wave_uses_current_anchor(self) -> None:
        detail = "Verified detail. " * 30
        progress_text = textwrap.dedent(
            f"""\
            # Progress Ledger

            ### ISING-W26 (2026-08-24)
            {detail}

            ### ISING3D WAVE 22 (2026-08-23)
            {detail}
            """
        )
        project = next(
            project for project in sync.PROJECTS
            if project["slug"] == "ising3d"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            progress = Path(tmpdir) / "PROGRESS.md"
            progress.write_text(progress_text)
            with patch.object(sync, "PROGRESS", progress):
                header, _ = sync.newest_progress_block(project)

        self.assertEqual(header, "### ISING-W26 (2026-08-24)")

    def test_new_campaigns_are_registered(self) -> None:
        registered = {project["slug"] for project in sync.PROJECTS}
        self.assertTrue({"e389", "liu_h1"} <= registered)


if __name__ == "__main__":
    unittest.main()
