#!/usr/bin/env python3
"""Behavioral tests for the frozen weekday participation runner."""

from __future__ import annotations

import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_weekday_participation as runner


DAY_MS = 86_400_000


def synthetic_accumulator(weekday_count: float = 120.0,
                          weekend_count: float = 100.0,
                          weekday_quote: float = 1200.0,
                          weekend_quote: float = 1000.0) -> dict:
    return {
        "total_bars": 10_080,
        "weekday_bars": 7_200,
        "weekend_bars": 2_880,
        "weekday_count_sum": weekday_count * 7_200,
        "weekend_count_sum": weekend_count * 2_880,
        "weekday_quote_sum": weekday_quote * 7_200,
        "weekend_quote_sum": weekend_quote * 2_880,
        "zero_trade_bars": 0,
        "zero_quote_bars": 0,
    }


def frozen_test_contract() -> dict:
    contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
    contract["windows"]["primary"] = {
        "start_utc_inclusive": "2024-01-01T00:00:00Z",
        "end_utc_exclusive": "2024-02-26T00:00:00Z",
        "expected_complete_weeks": 8,
        "expected_pooled_observations": 16,
    }
    contract["minimum_sample"]["primary"] = {
        "paired_weeks": 8,
        "per_asset": 8,
        "pooled": 16,
        "included_monday_months": 1,
    }
    return contract


def synthetic_weeks(start_ms: int, weekday_count: float,
                    weeks: int = 8) -> dict[int, dict]:
    return {
        start_ms + index * 7 * DAY_MS:
            synthetic_accumulator(weekday_count=weekday_count)
        for index in range(weeks)
    }


class WeekdayParticipationTests(unittest.TestCase):
    def test_complete_week_starts_are_monday_anchored(self) -> None:
        starts = runner.complete_week_starts(
            runner.parse_iso_ms("2024-07-01T00:00:00Z"),
            runner.parse_iso_ms("2025-01-01T00:00:00Z"))
        self.assertEqual(len(starts), 26)
        self.assertEqual(starts[0], runner.parse_iso_ms("2024-07-01T00:00:00Z"))
        self.assertEqual(starts[-1], runner.parse_iso_ms("2024-12-23T00:00:00Z"))

    def test_group_assignment_and_cardinality_are_frozen(self) -> None:
        week_ms = runner.parse_iso_ms("2024-01-01T00:00:00Z")
        accumulator = runner.new_accumulator()
        for offset in range(10_080):
            runner.add_minute(
                accumulator, week_ms + offset * runner.MINUTE_MS,
                quote_volume=1000.0, trades=100)
        self.assertEqual(accumulator["weekday_bars"], 7_200)
        self.assertEqual(accumulator["weekend_bars"], 2_880)
        observation = runner.weekly_observation(
            "BTCUSDT", week_ms, accumulator)
        self.assertEqual(observation["x_count"], 0.0)
        self.assertEqual(observation["x_quote"], 0.0)

    def test_weekly_observation_uses_weekday_weekend_ratios(self) -> None:
        week_ms = runner.parse_iso_ms("2024-01-01T00:00:00Z")
        observation = runner.weekly_observation(
            "ETHUSDT", week_ms, synthetic_accumulator())
        self.assertAlmostEqual(observation["x_count"], runner.math.log(1.2))
        self.assertAlmostEqual(observation["x_quote"], runner.math.log(1.2))

    def test_missing_or_zero_bar_rejects_asset_week(self) -> None:
        week_ms = runner.parse_iso_ms("2024-01-01T00:00:00Z")
        missing = synthetic_accumulator()
        missing["total_bars"] -= 1
        with self.assertRaisesRegex(ValueError, "complete minute grid"):
            runner.weekly_observation("BTCUSDT", week_ms, missing)
        zero = synthetic_accumulator()
        zero["zero_trade_bars"] = 1
        with self.assertRaisesRegex(ValueError, "zero trade"):
            runner.weekly_observation("BTCUSDT", week_ms, zero)

    def test_all_six_gates_must_pass_for_promotion(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        with mock.patch.object(
                runner, "load_asset_weeks",
                side_effect=lambda *_args: synthetic_weeks(start_ms, 120.0)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "pass-on-this-sample")
        self.assertEqual(set(result["gate_results"]), {
            "G1_count_magnitude", "G2_count_asset_consistency",
            "G3_quote_magnitude", "G4_quote_asset_consistency",
            "G5_month_stability", "G6_dependence_robustness",
        })
        self.assertTrue(all(
            gate["pass"] for gate in result["gate_results"].values()))

    def test_one_failed_asset_gate_falsifies_window(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        counts = {"BTCUSDT": 140.0, "ETHUSDT": 100.0}

        def load(_inventory: dict, asset: str, _start: int,
                 _end: int) -> dict:
            return synthetic_weeks(start_ms, counts[asset])

        with mock.patch.object(runner, "load_asset_weeks", side_effect=load):
            result = runner.run_window(contract, "primary", {})
        failed = {name for name, gate in result["gate_results"].items()
                  if not gate["pass"]}
        self.assertEqual(result["status"], "falsified")
        self.assertEqual(failed, {"G2_count_asset_consistency"})

    def test_insufficient_sample_is_inconclusive_without_gates(self) -> None:
        contract = frozen_test_contract()
        contract["minimum_sample"]["primary"].update({
            "paired_weeks": 9, "per_asset": 9, "pooled": 18})
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        with mock.patch.object(
                runner, "load_asset_weeks",
                side_effect=lambda *_args: synthetic_weeks(start_ms, 120.0)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "inconclusive")
        self.assertIsNone(result["gate_results"])

    def test_bootstrap_formula_tail_and_exact_two_asset_length(self) -> None:
        rows = []
        for index in range(8):
            week = f"2024-01-{index * 7 + 1:02d}"
            rows.extend([
                {"week": week, "asset": "BTCUSDT", "x_count": 0.2},
                {"week": week, "asset": "ETHUSDT", "x_count": 0.1},
            ])
        first = runner.block_bootstrap_p(
            rows, "x_count", replicates=100, seed=9, block_weeks=2)
        second = runner.block_bootstrap_p(
            rows, "x_count", replicates=100, seed=9, block_weeks=2)
        self.assertEqual(first, second)
        self.assertEqual(first, 1 / 101)
        zeros = [{**row, "x_count": 0.0} for row in rows]
        self.assertEqual(runner.block_bootstrap_p(
            zeros, "x_count", replicates=10, seed=9, block_weeks=2), 1.0)
        mixed = [{**row, "x_count": 1.0 if row["asset"] == "BTCUSDT" else -2.0}
                 for row in rows]
        self.assertEqual(runner.block_bootstrap_p(
            mixed, "x_count", replicates=10, seed=9, block_weeks=2), 1.0)
        lengths: list[int] = []

        def capture(values) -> float:
            values = list(values)
            lengths.append(len(values))
            return 1.0

        with mock.patch.object(runner, "median", side_effect=capture):
            p_value = runner.block_bootstrap_p(
                rows, "x_count", replicates=3, seed=9, block_weeks=2)
        self.assertEqual(lengths, [16, 16, 16])
        self.assertEqual(p_value, 1 / 4)

    def test_timestamp_normalization_and_parser_pin_both_units(self) -> None:
        minute_ms = runner.parse_iso_ms("2024-01-01T00:01:00Z")
        minute_us = minute_ms * 1000
        self.assertEqual(runner.normalize_epoch(str(minute_ms)), minute_ms)
        self.assertEqual(runner.normalize_epoch(str(minute_us)), minute_ms)
        boundary = runner.MICROSECOND_THRESHOLD
        self.assertEqual(runner.normalize_epoch(str(boundary)), boundary)
        self.assertEqual(
            runner.normalize_epoch(str(boundary + 1)), (boundary + 1) // 1000)

        week_ms = runner.parse_iso_ms("2024-01-01T00:00:00Z")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, raw_timestamp in (
                    ("milliseconds.zip", minute_ms),
                    ("microseconds.zip", minute_us)):
                row = [str(raw_timestamp), "100", "101", "99", "100", "1",
                       str(raw_timestamp), "1000", "10", "1", "1", "0"]
                with zipfile.ZipFile(root / name, "w") as archive:
                    archive.writestr(
                        "bars.csv", ",".join(runner.EXPECTED_HEADER) + "\n"
                        + ",".join(row) + "\n")
                weeks: dict[int, dict] = {}
                with mock.patch.object(runner, "WORKSPACE_ROOT", root):
                    runner._parse_archive(
                        {"path": name, "zip_member": "bars.csv"},
                        week_ms, week_ms + 7 * DAY_MS, {week_ms}, set(), weeks)
                self.assertEqual(weeks[week_ms]["total_bars"], 1)
                self.assertEqual(weeks[week_ms]["weekday_bars"], 1)

    def test_archive_verification_checks_sidecar_crc_member_and_header(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "sample.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr(
                    "sample.csv", ",".join(runner.EXPECTED_HEADER) + "\n")
            digest = runner.sha256_file(archive)
            sidecar = root / "sample.zip.CHECKSUM"
            sidecar.write_text(f"{digest}  sample.zip\n", encoding="utf-8")
            entry = {
                "path": "sample.zip", "sha256": digest,
                "checksum_sidecar": "sample.zip.CHECKSUM",
                "checksum_sidecar_sha256": runner.sha256_file(sidecar),
                "size_bytes": archive.stat().st_size,
                "zip_member": "sample.csv",
                "header": ",".join(runner.EXPECTED_HEADER),
            }
            self.assertEqual(runner.verify_archive(entry, root), "sample.csv")
            sidecar.write_text(f"{'0' * 64}  sample.zip\n", encoding="utf-8")
            entry["checksum_sidecar_sha256"] = runner.sha256_file(sidecar)
            with self.assertRaisesRegex(ValueError, "upstream sidecar"):
                runner.verify_archive(entry, root)

    def test_replication_requires_same_hashes_and_primary_pass(self) -> None:
        expected = {
            "contract_sha256": "c", "input_inventory_sha256": "i",
            "prior_window_audit_sha256": "a", "loop_policy_sha256": "p",
            "runner_sha256": "r",
        }
        primary = {**expected, "role": "primary",
                   "result": {"status": "falsified"}}
        with self.assertRaisesRegex(ValueError, "primary pass"):
            runner.validate_primary_artifact(primary, expected)
        primary["result"]["status"] = "pass-on-this-sample"
        runner.validate_primary_artifact(primary, expected)
        primary["runner_sha256"] = "different"
        with self.assertRaisesRegex(ValueError, "runner_sha256"):
            runner.validate_primary_artifact(primary, expected)

    def test_write_json_new_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            runner.write_json_new(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                runner.write_json_new(path, {"a": 2})
            self.assertEqual(json.loads(path.read_text()), {"a": 1})

    def test_existing_output_is_refused_before_any_read(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "exists.json"
            output.write_text("{}\n", encoding="utf-8")
            with mock.patch.object(runner, "verify_static_inputs") as verify, \
                    mock.patch.object(runner, "run_window") as run:
                with self.assertRaisesRegex(SystemExit, "output already exists"):
                    runner.main([
                        "--contract-sha256", "unused", "--window", "primary",
                        "--output", str(output),
                    ])
                verify.assert_not_called()
                run.assert_not_called()

    def test_integrity_or_schema_failure_writes_no_artifact(self) -> None:
        contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            with mock.patch.object(
                    runner, "verify_static_inputs",
                    side_effect=ValueError("drift")), \
                    mock.patch.object(runner, "run_window") as run:
                with self.assertRaisesRegex(SystemExit, "static integrity"):
                    runner.main([
                        "--contract-sha256", contract_sha, "--window", "primary",
                        "--output", str(output),
                    ])
                run.assert_not_called()
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
