#!/usr/bin/env python3
"""Behavioral tests for the frozen 21:00-23:00 UTC net-premium runner."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_utc_2100_2300_net_premium as runner


DAY_MS = 86_400_000
HOUR_MS = 3_600_000
FEE = 0.0005


def synthetic_bars(day_ms: int, p_in: float = 100.0, p_out: float = 101.0,
                   p_0: float = 100.0, p_24: float = 100.0) -> dict[int, dict]:
    """Twenty-four hourly bars with the four frozen price anchors set."""
    bars = {}
    for hour in range(24):
        open_ms = day_ms + hour * HOUR_MS
        bars[open_ms] = {"open": 100.0, "close": 100.0}
    bars[day_ms]["open"] = p_0
    bars[day_ms + 21 * HOUR_MS]["open"] = p_in
    bars[day_ms + 23 * HOUR_MS]["open"] = p_out
    bars[day_ms + 23 * HOUR_MS]["close"] = p_24
    return bars


def frozen_test_contract() -> dict:
    contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
    contract["windows"]["primary"] = {
        "start_utc_inclusive": "2020-01-01T00:00:00Z",
        "end_utc_exclusive": "2020-04-01T00:00:00Z",
        "expected_dates": 91,
        "expected_pooled_observations": 182,
    }
    contract["minimum_sample"]["primary"] = {
        "paired_dates": 60,
        "per_asset": 60,
        "pooled": 120,
        "included_quarters": 1,
    }
    return contract


def synthetic_days(start_ms: int, p_out: float, days: int = 91,
                   p_24: float = 100.0,
                   duplicated: frozenset[int] = frozenset()
                   ) -> tuple[dict[int, dict[int, dict]], set[int]]:
    return ({
        start_ms + index * DAY_MS: synthetic_bars(
            start_ms + index * DAY_MS, p_out=p_out, p_24=p_24)
        for index in range(days)
    }, set(duplicated))


class DateGridTests(unittest.TestCase):
    def test_expected_dates_cover_frozen_windows_exactly(self) -> None:
        primary = runner.expected_dates(
            runner.parse_iso_ms("2020-01-01T00:00:00Z"),
            runner.parse_iso_ms("2022-01-01T00:00:00Z"))
        replication = runner.expected_dates(
            runner.parse_iso_ms("2022-04-01T00:00:00Z"),
            runner.parse_iso_ms("2024-01-01T00:00:00Z"))
        self.assertEqual(len(primary), 731)
        self.assertEqual(len(replication), 640)
        self.assertEqual(primary[0], runner.parse_iso_ms("2020-01-01T00:00:00Z"))
        self.assertEqual(primary[-1], runner.parse_iso_ms("2021-12-31T00:00:00Z"))
        self.assertEqual(replication[-1],
                         runner.parse_iso_ms("2023-12-31T00:00:00Z"))
        self.assertTrue(all(day % DAY_MS == 0 for day in primary))

    def test_quarter_labels_follow_calendar_quarters(self) -> None:
        self.assertEqual(runner.quarter_label("2020-01-01"), "2020Q1")
        self.assertEqual(runner.quarter_label("2020-03-31"), "2020Q1")
        self.assertEqual(runner.quarter_label("2020-04-01"), "2020Q2")
        self.assertEqual(runner.quarter_label("2023-12-31"), "2023Q4")


class DailyObservationTests(unittest.TestCase):
    def test_net_return_charges_both_taker_fees_exactly(self) -> None:
        day_ms = runner.parse_iso_ms("2020-01-01T00:00:00Z")
        flat = runner.daily_observation(
            "BTCUSDT", day_ms, synthetic_bars(day_ms, p_out=100.0))
        self.assertAlmostEqual(flat["net"], -2 * FEE, places=15)
        self.assertEqual(flat["x"], 0.0)
        up = runner.daily_observation(
            "BTCUSDT", day_ms, synthetic_bars(day_ms, p_out=101.0))
        self.assertAlmostEqual(up["net"], 1.01 * (1 - FEE) - (1 + FEE), places=15)
        self.assertEqual(up["date"], "2020-01-01")
        self.assertEqual(up["quarter"], "2020Q1")

    def test_concentration_subtracts_pro_rata_day_return(self) -> None:
        day_ms = runner.parse_iso_ms("2020-06-15T00:00:00Z")
        bars = synthetic_bars(day_ms, p_in=100.0, p_out=102.0, p_0=100.0,
                              p_24=112.0)
        observation = runner.daily_observation("ETHUSDT", day_ms, bars)
        expected = math.log(1.02) - (2 / 24) * math.log(1.12)
        self.assertAlmostEqual(observation["x"], expected, places=15)
        self.assertEqual(observation["quarter"], "2020Q2")

    def test_missing_duplicate_or_nonpositive_bar_rejects_date(self) -> None:
        day_ms = runner.parse_iso_ms("2020-01-01T00:00:00Z")
        missing = synthetic_bars(day_ms)
        del missing[day_ms + 5 * HOUR_MS]
        with self.assertRaisesRegex(ValueError, "complete hourly grid"):
            runner.daily_observation("BTCUSDT", day_ms, missing)
        nonpositive = synthetic_bars(day_ms)
        nonpositive[day_ms + 21 * HOUR_MS]["open"] = 0.0
        with self.assertRaisesRegex(ValueError, "nonpositive price"):
            runner.daily_observation("BTCUSDT", day_ms, nonpositive)
        self.assertEqual(
            runner._rejection_key(ValueError("complete hourly grid missing")),
            "incomplete-hourly-grid")
        with self.assertRaises(ValueError):
            runner._rejection_key(ValueError("internal invariant"))


class GateTests(unittest.TestCase):
    def test_all_five_gates_must_pass_for_promotion(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        with mock.patch.object(
                runner, "load_asset_days",
                side_effect=lambda *_args: synthetic_days(start_ms, 101.0)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "pass-on-this-sample")
        self.assertEqual(set(result["gate_results"]), {
            "G1_net_magnitude", "G2_net_asset_consistency",
            "G3_quarter_stability", "G4_concentration",
            "G5_dependence_robustness",
        })
        self.assertTrue(all(
            gate["pass"] for gate in result["gate_results"].values()))
        self.assertEqual(result["sample"]["paired_dates"], 91)
        self.assertEqual(result["sample"]["pooled_observations"], 182)
        self.assertEqual(result["gate_results"]["G5_dependence_robustness"]
                         ["one_sided_block_bootstrap_p_net"], 1 / 20001)

    def test_one_failed_asset_gate_falsifies_window(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])

        def loader(_inventory, asset, *_args):
            return synthetic_days(start_ms, 101.0 if asset == "BTCUSDT" else 100.1)

        with mock.patch.object(runner, "load_asset_days", side_effect=loader):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items()
                  if not gate["pass"]}
        self.assertEqual(failed, {"G2_net_asset_consistency"})

    def test_fee_sized_gross_return_fails_magnitude_gate(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        # Gross +0.12% per date nets +0.02% (below the 3 bps floor).
        with mock.patch.object(
                runner, "load_asset_days",
                side_effect=lambda *_args: synthetic_days(start_ms, 100.12)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items()
                  if not gate["pass"]}
        self.assertEqual(failed, {"G1_net_magnitude"})

    def test_sub_pro_rata_window_fails_concentration_gate(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        # 21:00-23:00 rises 1% net-positive, but the whole day rises 13x that
        # in log terms, so the two-hour window carries less than its pro-rata
        # share and the drift-neutral concentration gate must fail.
        p_24 = 100.0 * math.exp(13 * math.log(1.01))
        with mock.patch.object(
                runner, "load_asset_days",
                side_effect=lambda *_args: synthetic_days(
                    start_ms, 101.0, p_24=p_24)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items()
                  if not gate["pass"]}
        self.assertEqual(failed, {"G4_concentration", "G5_dependence_robustness"})

    def test_insufficient_sample_is_inconclusive_without_gates(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        with mock.patch.object(
                runner, "load_asset_days",
                side_effect=lambda *_args: synthetic_days(start_ms, 101.0, days=40)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "inconclusive")
        self.assertIsNone(result["gate_results"])
        self.assertEqual(result["sample"]["paired_dates"], 40)
        self.assertEqual(
            result["sample"]["rejected_asset_dates"]["BTCUSDT"]
            ["incomplete-hourly-grid"], 51)

    def test_quarter_stability_counts_included_quarters(self) -> None:
        rows = []
        for date, net in (("2020-01-05", 0.01), ("2020-01-06", 0.01),
                          ("2020-04-05", -0.01), ("2020-04-06", -0.01),
                          ("2020-07-05", 0.01), ("2020-07-06", 0.02)):
            rows.append({"date": date, "quarter": runner.quarter_label(date),
                         "net": net, "x": 0.0})
        stability = runner._quarter_stability(rows, "net")
        self.assertEqual(stability["included_quarters"], 3)
        self.assertEqual(stability["positive_quarters"], 2)
        self.assertAlmostEqual(stability["positive_fraction"], 2 / 3)

    def test_duplicate_timestamp_rejects_only_that_asset_date(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        duplicated = frozenset({start_ms + 3 * DAY_MS, start_ms + 40 * DAY_MS})

        def loader(_inventory, asset, *_args):
            return synthetic_days(
                start_ms, 101.0,
                duplicated=duplicated if asset == "ETHUSDT" else frozenset())

        with mock.patch.object(runner, "load_asset_days", side_effect=loader):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "pass-on-this-sample")
        self.assertEqual(result["sample"]["paired_dates"], 89)
        self.assertEqual(result["sample"]["rejected_asset_dates"]["ETHUSDT"],
                         {"duplicate-timestamp-bar": 2})
        self.assertEqual(result["sample"]["rejected_asset_dates"]["BTCUSDT"], {})
        self.assertEqual(result["sample"]["unpaired_complete_dates"]["BTCUSDT"], 2)

    def test_quarter_instability_alone_falsifies_window(self) -> None:
        contract = frozen_test_contract()
        contract["windows"]["primary"]["end_utc_exclusive"] = "2020-07-01T00:00:00Z"
        contract["windows"]["primary"]["expected_dates"] = 182
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])

        def loader(*_args):
            first, _ = synthetic_days(start_ms, 102.0, days=91)
            second, _ = synthetic_days(start_ms + 91 * DAY_MS, 100.09, days=91)
            return {**first, **second}, set()

        with mock.patch.object(runner, "load_asset_days", side_effect=loader):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items()
                  if not gate["pass"]}
        self.assertEqual(failed, {"G3_quarter_stability"})
        stability = result["gate_results"]["G3_quarter_stability"]["net"]
        self.assertEqual(stability["included_quarters"], 2)
        self.assertEqual(stability["positive_quarters"], 1)

    def test_quarter_stability_uses_exact_three_fifths_rule(self) -> None:
        contract = frozen_test_contract()
        contract["windows"]["primary"]["end_utc_exclusive"] = "2021-04-01T00:00:00Z"
        contract["windows"]["primary"]["expected_dates"] = 456
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        quarter_starts = ["2020-01-01", "2020-04-01", "2020-07-01",
                          "2020-10-01", "2021-01-01"]
        signs = [102.0, 102.0, 102.0, 100.09, 100.09]  # 3 of 5 positive

        def loader(*_args):
            days: dict[int, dict[int, dict]] = {}
            for index, quarter in enumerate(quarter_starts):
                q_start = runner.parse_iso_ms(f"{quarter}T00:00:00Z")
                q_end = (runner.parse_iso_ms(f"{quarter_starts[index + 1]}T00:00:00Z")
                         if index + 1 < len(quarter_starts)
                         else runner.parse_iso_ms("2021-04-01T00:00:00Z"))
                bars, _ = synthetic_days(q_start, signs[index],
                                         days=(q_end - q_start) // DAY_MS)
                days.update(bars)
            self.assertEqual(len(days), 456)
            return days, set()

        with mock.patch.object(runner, "load_asset_days", side_effect=loader):
            result = runner.run_window(contract, "primary", {})
        gate = result["gate_results"]["G3_quarter_stability"]
        self.assertEqual(gate["net"]["included_quarters"], 5)
        self.assertEqual(gate["net"]["positive_quarters"], 3)
        self.assertTrue(gate["pass"])


class BootstrapTests(unittest.TestCase):
    def test_block_bootstrap_is_deterministic_and_date_paired(self) -> None:
        rows = []
        for index in range(30):
            date = f"2020-01-{index + 1:02d}"
            for asset, net in (("BTCUSDT", 0.001 * (index % 5 - 1)),
                               ("ETHUSDT", 0.0005)):
                rows.append({"asset": asset, "date": date, "net": net, "x": 0.0})
        first = runner.block_bootstrap_p(rows, "net", replicates=200)
        second = runner.block_bootstrap_p(rows, "net", replicates=200)
        self.assertEqual(first, second)
        self.assertGreater(first, 0.0)
        self.assertLessEqual(first, 1.0)
        unpaired = rows[:-1]
        with self.assertRaisesRegex(ValueError, "asset-paired"):
            runner.block_bootstrap_p(unpaired, "net", replicates=10)

    def test_draw_truncates_circular_blocks_to_sample_size(self) -> None:
        dates = [f"d{i}" for i in range(10)]
        picked = runner._draw_with_rng(dates, runner.random.Random(1), 7)
        self.assertEqual(len(picked), 10)
        self.assertTrue(set(picked) <= set(dates))


def _write_archive(directory: Path, name: str, rows: list[list[str]],
                   header: bool) -> Path:
    archive = directory / name
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    if header:
        writer.writerow(runner.EXPECTED_HEADER)
    writer.writerows(rows)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(name.replace(".zip", ".csv"), buffer.getvalue())
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    sidecar = archive.with_name(archive.name + ".CHECKSUM")
    sidecar.write_text(f"{digest}  {archive.name}\n", encoding="utf-8")
    return archive


def _entry(root: Path, archive: Path, header_present: bool) -> dict:
    sidecar = archive.with_name(archive.name + ".CHECKSUM")
    return {
        "path": str(archive.relative_to(root)),
        "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "checksum_sidecar": str(sidecar.relative_to(root)),
        "checksum_sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
        "size_bytes": archive.stat().st_size,
        "zip_member": archive.name.replace(".zip", ".csv"),
        "header": ",".join(runner.EXPECTED_HEADER),
        "header_present": header_present,
    }


def _bar_row(open_ms: int, open_price: float = 100.0,
             close_price: float = 100.0) -> list[str]:
    return [str(open_ms), f"{open_price:.2f}", "100.00", "100.00",
            f"{close_price:.2f}", "1.0", str(open_ms + HOUR_MS - 1), "100.0",
            "10", "0.5", "50.0", "0"]


class ArchiveVerificationTests(unittest.TestCase):
    def test_archive_form_must_match_inventory_header_flag(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            day_ms = runner.parse_iso_ms("2020-01-01T00:00:00Z")
            rows = [_bar_row(day_ms + hour * HOUR_MS) for hour in range(24)]
            headerless = _write_archive(root, "BTCUSDT-1h-2020-01.zip", rows, False)
            headed = _write_archive(root, "BTCUSDT-1h-2022-01.zip", rows, True)
            self.assertEqual(
                runner.verify_archive(_entry(root, headerless, False), root),
                "BTCUSDT-1h-2020-01.csv")
            self.assertEqual(
                runner.verify_archive(_entry(root, headed, True), root),
                "BTCUSDT-1h-2022-01.csv")
            with self.assertRaisesRegex(ValueError, "first-row form"):
                runner.verify_archive(_entry(root, headerless, True), root)
            with self.assertRaisesRegex(ValueError, "first-row form"):
                runner.verify_archive(_entry(root, headed, False), root)

    def test_sidecar_and_size_drift_are_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            day_ms = runner.parse_iso_ms("2020-01-01T00:00:00Z")
            rows = [_bar_row(day_ms + hour * HOUR_MS) for hour in range(24)]
            archive = _write_archive(root, "ETHUSDT-1h-2020-01.zip", rows, False)
            entry = _entry(root, archive, False)
            drifted = dict(entry, size_bytes=entry["size_bytes"] + 1)
            with self.assertRaisesRegex(ValueError, "byte size drift"):
                runner.verify_archive(drifted, root)
            sidecar = archive.with_name(archive.name + ".CHECKSUM")
            sidecar.write_text("0" * 64 + f"  {archive.name}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "sidecar file hash drift"):
                runner.verify_archive(entry, root)

    def test_parse_archive_reads_both_forms_and_flags_duplicate_days(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            day_ms = runner.parse_iso_ms("2020-01-01T00:00:00Z")
            rows = [_bar_row(day_ms + hour * HOUR_MS, 100.0 + hour, 100.0 + hour)
                    for hour in range(24)]
            for header in (False, True):
                archive = _write_archive(
                    root, f"BTCUSDT-1h-20{20 + int(header)}-01.zip", rows, header)
                entry = _entry(root, archive, header)
                days: dict[int, dict[int, dict]] = {}
                duplicated_days: set[int] = set()
                runner._parse_archive(entry, day_ms, day_ms + DAY_MS,
                                      {day_ms}, set(), days, duplicated_days, root)
                self.assertEqual(len(days[day_ms]), 24)
                self.assertEqual(days[day_ms][day_ms + 21 * HOUR_MS]["open"], 121.0)
                self.assertEqual(duplicated_days, set())
            duplicated = _write_archive(
                root, "ETHUSDT-1h-2020-01.zip", rows + rows[:1], False)
            days, duplicated_days = {}, set()
            runner._parse_archive(_entry(root, duplicated, False), day_ms,
                                  day_ms + DAY_MS, {day_ms}, set(), days,
                                  duplicated_days, root)
            self.assertEqual(duplicated_days, {day_ms})
            self.assertEqual(len(days[day_ms]), 24)


class LifecycleTests(unittest.TestCase):
    def test_write_json_new_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            runner.write_json_new(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                runner.write_json_new(path, {"a": 2})
            self.assertEqual(json.loads(path.read_text()), {"a": 1})

    def test_existing_output_is_refused_before_any_work(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "exists.json"
            output.write_text("{}", encoding="utf-8")
            with mock.patch.object(runner, "verify_static_inputs") as verify:
                with self.assertRaisesRegex(SystemExit, "already exists"):
                    runner.main(["--contract-sha256", "unused",
                                 "--window", "primary",
                                 "--output", str(output)])
            verify.assert_not_called()

    def test_contract_hash_mismatch_prevents_outcome_read(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "never.json"
            with mock.patch.object(runner, "verify_static_inputs") as verify:
                with self.assertRaisesRegex(SystemExit, "contract hash mismatch"):
                    runner.main(["--contract-sha256", "0" * 64,
                                 "--window", "primary",
                                 "--output", str(output)])
            verify.assert_not_called()
            self.assertFalse(output.exists())

    def test_replication_requires_committed_primary_pass(self) -> None:
        expected = {"contract_sha256": "c", "input_inventory_sha256": "i",
                    "prior_window_audit_sha256": "a", "loop_policy_sha256": "p",
                    "runner_sha256": "r"}
        primary = {**expected, "role": "primary",
                   "result": {"status": "falsified"}}
        with self.assertRaisesRegex(ValueError, "committed primary pass"):
            runner.validate_primary_artifact(primary, expected)
        primary["result"]["status"] = "pass-on-this-sample"
        runner.validate_primary_artifact(primary, expected)
        drifted = dict(primary, runner_sha256="other")
        with self.assertRaisesRegex(ValueError, "runner_sha256 mismatch"):
            runner.validate_primary_artifact(drifted, expected)

    def test_static_integrity_failure_writes_no_artifact(self) -> None:
        contract_sha = hashlib.sha256(
            runner.CONTRACT_PATH.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "never.json"
            with mock.patch.object(
                    runner, "verify_static_inputs",
                    side_effect=ValueError("archive hash drift: x")):
                with mock.patch.object(runner, "run_window") as run:
                    with self.assertRaisesRegex(SystemExit, "static integrity"):
                        runner.main(["--contract-sha256", contract_sha,
                                     "--window", "primary",
                                     "--output", str(output)])
            run.assert_not_called()
            self.assertFalse(output.exists())

    def test_main_primary_then_replication_end_to_end(self) -> None:
        contract_sha = hashlib.sha256(
            runner.CONTRACT_PATH.read_bytes()).hexdigest()
        integrity = {"input_inventory_sha256": "i", "prior_window_audit_sha256": "a",
                     "loop_policy_sha256": "p", "archives_verified": 0}
        passing = {"window": {}, "status": "pass-on-this-sample", "reason": None,
                   "sample": {"paired_dates": 600, "pooled_observations": 1200},
                   "gate_results": {}}
        # Artifacts must live inside the workspace root: the runner records
        # the primary artifact path relative to it and refuses outsiders.
        with tempfile.TemporaryDirectory(dir=runner.EXPLORATIONS,
                                         prefix=".tmp-utc2123-") as directory:
            primary = Path(directory) / "primary.json"
            replication = Path(directory) / "replication.json"
            with mock.patch.object(runner, "verify_static_inputs",
                                   return_value=({}, integrity)), \
                 mock.patch.object(runner, "run_window", return_value=passing):
                self.assertEqual(runner.main([
                    "--contract-sha256", contract_sha, "--window", "primary",
                    "--output", str(primary)]), 0)
                artifact = json.loads(primary.read_text(encoding="utf-8"))
                self.assertEqual(artifact["role"], "primary")
                self.assertEqual(artifact["contract_sha256"], contract_sha)
                primary_sha = hashlib.sha256(primary.read_bytes()).hexdigest()
                with self.assertRaisesRegex(SystemExit, "requires --primary-artifact"):
                    runner.main(["--contract-sha256", contract_sha,
                                 "--window", "replication",
                                 "--output", str(replication)])
                with self.assertRaisesRegex(SystemExit, "does not match committed"):
                    runner.main(["--contract-sha256", contract_sha,
                                 "--window", "replication",
                                 "--output", str(replication),
                                 "--primary-artifact", str(primary),
                                 "--primary-artifact-sha256", "0" * 64])
                self.assertFalse(replication.exists())
                self.assertEqual(runner.main([
                    "--contract-sha256", contract_sha, "--window", "replication",
                    "--output", str(replication),
                    "--primary-artifact", str(primary),
                    "--primary-artifact-sha256", primary_sha]), 0)
                replicated = json.loads(replication.read_text(encoding="utf-8"))
                self.assertEqual(replicated["role"], "replication")
                self.assertEqual(replicated["primary_artifact_sha256"], primary_sha)
                with self.assertRaisesRegex(SystemExit, "already exists"):
                    runner.main(["--contract-sha256", contract_sha,
                                 "--window", "replication",
                                 "--output", str(replication),
                                 "--primary-artifact", str(primary),
                                 "--primary-artifact-sha256", primary_sha])

    def test_main_replication_refuses_falsified_primary(self) -> None:
        contract_sha = hashlib.sha256(
            runner.CONTRACT_PATH.read_bytes()).hexdigest()
        integrity = {"input_inventory_sha256": "i", "prior_window_audit_sha256": "a",
                     "loop_policy_sha256": "p", "archives_verified": 0}
        falsified = {"window": {}, "status": "falsified", "reason": "gate false",
                     "sample": {"paired_dates": 600, "pooled_observations": 1200},
                     "gate_results": {}}
        with tempfile.TemporaryDirectory() as directory:
            primary = Path(directory) / "primary.json"
            replication = Path(directory) / "replication.json"
            with mock.patch.object(runner, "verify_static_inputs",
                                   return_value=({}, integrity)), \
                 mock.patch.object(runner, "run_window", return_value=falsified) as run:
                runner.main(["--contract-sha256", contract_sha, "--window",
                             "primary", "--output", str(primary)])
                primary_sha = hashlib.sha256(primary.read_bytes()).hexdigest()
                with self.assertRaisesRegex(SystemExit, "committed primary pass"):
                    runner.main(["--contract-sha256", contract_sha,
                                 "--window", "replication",
                                 "--output", str(replication),
                                 "--primary-artifact", str(primary),
                                 "--primary-artifact-sha256", primary_sha])
                self.assertEqual(run.call_count, 1)
            self.assertFalse(replication.exists())


if __name__ == "__main__":
    unittest.main()
