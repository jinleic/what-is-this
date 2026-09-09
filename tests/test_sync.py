from __future__ import annotations

import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path


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
            header, _ = sync.newest_progress_block(project, progress)

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
            header, _ = sync.newest_progress_block(project, progress)

        self.assertEqual(header, "### ISING-W26 (2026-08-24)")

    def test_new_campaigns_are_registered(self) -> None:
        registered = {project["slug"] for project in sync.PROJECTS}
        self.assertTrue({"e389", "liu_h1"} <= registered)

    def test_math_domain_mirror_includes_research_status(self) -> None:
        self.assertIn("RESEARCH_STATUS.md", sync.DOMAINS["math"]["shared_files"])


class DomainSyncTest(unittest.TestCase):
    def test_research_domains_and_targets_are_registered(self) -> None:
        self.assertEqual(
            set(sync.DOMAINS),
            {"math", "physics", "cs", "quant-trading"},
        )
        registered = {
            (project.get("domain", "math"), project["slug"])
            for project in sync.PROJECTS
        }
        self.assertTrue(
            {
                ("physics", "qldpc-dec"),
                ("physics", "msd"),
                ("physics", "shadows"),
                ("physics", "na-compiler"),
                ("physics", "fss-bb"),
                ("physics", "qlops"),
                ("cs", "mceliece"),
                ("cs", "kg"),
                ("cs", "omega"),
                ("cs", "delcap"),
                ("cs", "oct-rank"),
                ("cs", "rs-pe3d"),
                ("cs", "mm3"),
                ("quant-trading", "overview"),
            }
            <= registered
        )

    def test_math_urls_stay_stable_and_new_domains_are_namespaced(self) -> None:
        math_project = next(
            project for project in sync.PROJECTS
            if project["slug"] == "uc"
        )
        physics_project = next(
            project for project in sync.PROJECTS
            if project.get("domain") == "physics"
            and project["slug"] == "qldpc-dec"
        )
        quant_project = next(
            project for project in sync.PROJECTS
            if project.get("domain") == "quant-trading"
        )

        self.assertEqual(sync.project_page(math_project), "problems/uc.html")
        self.assertEqual(sync.project_material(math_project), "materials/uc")
        self.assertEqual(
            sync.project_page(physics_project),
            "physics/qldpc-dec.html",
        )
        self.assertEqual(
            sync.project_material(physics_project),
            "materials/physics/qldpc-dec",
        )
        self.assertEqual(
            sync.project_page(quant_project),
            "quant-trading/index.html",
        )
        self.assertEqual(
            sync.project_material(quant_project),
            "materials/quant-trading",
        )

    def test_progress_selection_uses_the_requested_domain_ledger(self) -> None:
        detail = "Verified domain detail. " * 30
        physics_text = textwrap.dedent(
            f"""\
            # Progress Ledger

            ## 2026-08-31 — qldpc-dec gate closed
            {detail}
            """
        )
        project = {
            "domain": "physics",
            "slug": "qldpc-dec",
            "keywords": ("qldpc-dec",),
            "anchor": None,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            progress = Path(tmpdir) / "PROGRESS.md"
            progress.write_text(physics_text)
            header, _ = sync.newest_progress_block(project, progress)

        self.assertEqual(
            header,
            "## 2026-08-31 — qldpc-dec gate closed",
        )

    def test_progress_anchor_is_literal(self) -> None:
        detail = "Verified anchor detail. " * 30
        progress_text = textwrap.dedent(
            f"""\
            # Progress Ledger

            ### KOBON CURRENT
            {detail}

            ### KOBON (2026-08-21/22)
            {detail}
            """
        )
        project = {
            "slug": "kobon",
            "keywords": ("KOBON",),
            "anchor": "KOBON (2026-08-21/22)",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            progress = Path(tmpdir) / "PROGRESS.md"
            progress.write_text(progress_text)
            header, _ = sync.newest_progress_block(project, progress)

        self.assertEqual(header, "### KOBON (2026-08-21/22)")

    def test_progress_links_distinguish_domain_root_documents(self) -> None:
        entry = {
            "domain": "physics",
            "name": "Target",
            "source": "physics/target",
            "material": "materials/physics/target",
            "session": "",
            "updated": "2026-08-31T00:00:00Z",
            "shared_files": ["RESULTS.md"],
            "progress": "## Result\n\nVerified detail.",
            "readme_head": "",
            "files": [],
        }

        _, _, links = sync.build_problem(entry)

        self.assertEqual(
            links["md-progress"],
            {
                "base": "../materials/physics/",
                "root": "../materials/domains/physics/",
            },
        )

    def test_relative_path_exclusion_keeps_derived_data(self) -> None:
        project = {
            "papers": (),
            "exclude_paths": ("data/raw",),
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir)
            (source / "data" / "raw").mkdir(parents=True)
            (source / "data" / "derived").mkdir(parents=True)
            (source / "data" / "raw" / "ticks.json").write_text("{}")
            (source / "data" / "derived" / "summary.json").write_text("{}")

            files = sync.collect_sources(project, source)

        self.assertNotIn("data/raw/ticks.json", files)
        self.assertIn("data/derived/summary.json", files)

    def test_high_confidence_secret_patterns_are_blocked(self) -> None:
        self.assertTrue(sync.sensitive_hit(b"token=sk-" + b"x" * 32))
        self.assertFalse(sync.sensitive_hit(b"token loaded from environment"))


class DriftGuardTest(unittest.TestCase):
    def _run(self, root: Path, projects: list[dict]) -> dict[str, list[str]]:
        domains = {"math": {"root": root}}
        return sync.drift_guard(domains=domains, projects=projects)

    def test_unregistered_state_json_is_flagged(self) -> None:
        # dirname/slug case mismatch is the trap: registration matches the
        # on-disk dirname (LIU_H1), not the slug (liu_h1).
        projects = [dict(slug="liu_h1", dirname="LIU_H1")]
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for name in ("LIU_H1", "newtarget", "scratch", ".hidden"):
                d = root / name
                d.mkdir()
                (d / "state.json").write_text("{}")
            (root / "nostate").mkdir()

            drift = self._run(root, projects)

        self.assertEqual(drift["unregistered"], ["math/newtarget"])

    def test_registered_but_missing_source_is_flagged(self) -> None:
        projects = [
            dict(slug="overview", domain="quant-trading", dirname="."),
            dict(slug="gone", dirname="gone"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            drift = self._run(Path(tmpdir), projects)

        self.assertEqual(drift["missing"], ["math/gone"])
        self.assertEqual(drift["unregistered"], [])

    def test_ignored_dirs_are_not_flagged(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            smoke = root / "smoke"
            smoke.mkdir()
            (smoke / "state.json").write_text("{}")

            drift = self._run(root, [])

        self.assertEqual(drift["unregistered"], [])


if __name__ == "__main__":
    unittest.main()
