#!/usr/bin/env python3
"""Behavioral tests for the frozen quarter-hour clock-phase runner."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_quarter_hour_algorithmic_burst as runner


DAY_MS = 86_400_000
MINUTE_MS = 60_000


def required_bars(day_ms: int, q_trades: int = 110,
                  c_trades: int = 100) -> dict[int, dict]:
    bars: dict[int, dict] = {}
    for hour in range(24):
        for minute in runner.QUARTER_MINUTES | runner.CONTROL_MINUTES:
            is_quarter = minute in runner.QUARTER_MINUTES
            absolute_return = 0.002 if is_quarter else 0.001
            bars[day_ms + (hour * 60 + minute) * MINUTE_MS] = {
                "open": 100.0,
                "close": 100.0 * math.exp(absolute_return),
                "number_of_trades": q_trades if is_quarter else c_trades,
            }
    return bars


def frozen_test_contract() -> dict:
    contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
    contract["windows"]["primary"] = {
        "start_utc_inclusive": "2025-01-01T00:00:00Z",
        "end_utc_exclusive": "2025-01-22T00:00:00Z",
        "expected_dates": 21,
        "expected_pooled_observations": 42,
    }
    contract["minimum_sample"]["primary"] = {
        "paired_dates": 21,
        "per_asset": 21,
        "pooled": 42,
        "included_calendar_months": 1,
    }
    return contract


def synthetic_days(start_ms: int, q_trades: int,
                   dates: int = 21) -> dict[int, dict[int, dict]]:
    return {
        start_ms + index * DAY_MS:
            required_bars(start_ms + index * DAY_MS, q_trades=q_trades)
        for index in range(dates)
    }


class QuarterHourRunnerTests(unittest.TestCase):
    def test_frozen_clock_groups_exclude_minute_zero(self) -> None:
        self.assertEqual(runner.QUARTER_MINUTES, {15, 30, 45})
        self.assertEqual(
            runner.CONTROL_MINUTES, {5, 10, 20, 25, 35, 40, 50, 55})
        self.assertNotIn(0, runner.QUARTER_MINUTES | runner.CONTROL_MINUTES)
        self.assertEqual(len(runner.required_offsets()), 264)

    def test_daily_observation_uses_count_and_within_bar_return_ratios(self) -> None:
        day_ms = 1_767_225_600_000
        observation = runner.daily_observation(
            "BTCUSDT", day_ms, required_bars(day_ms))
        self.assertEqual(observation["quarter_bars"], 72)
        self.assertEqual(observation["control_bars"], 192)
        self.assertAlmostEqual(observation["x_count"], math.log(1.10), places=12)
        self.assertAlmostEqual(observation["x_absret"], math.log(2.0), places=10)

    def test_top_of_hour_bar_cannot_change_observation(self) -> None:
        day_ms = 1_767_225_600_000
        bars = required_bars(day_ms)
        baseline = runner.daily_observation("BTCUSDT", day_ms, bars)
        for hour in range(24):
            bars[day_ms + hour * 60 * MINUTE_MS] = {
                "open": 1.0, "close": 1_000_000.0,
                "number_of_trades": 1_000_000_000,
            }
        observed = runner.daily_observation("BTCUSDT", day_ms, bars)
        self.assertEqual(observed, baseline)

    def test_missing_required_bar_fails_closed(self) -> None:
        day_ms = 1_767_225_600_000
        bars = required_bars(day_ms)
        del bars[min(bars)]
        with self.assertRaisesRegex(ValueError, "required bars"):
            runner.daily_observation("ETHUSDT", day_ms, bars)

    def test_zero_trade_required_bar_rejects_asset_date(self) -> None:
        day_ms = 1_767_225_600_000
        bars = required_bars(day_ms)
        bars[min(bars)]["number_of_trades"] = 0
        with self.assertRaisesRegex(ValueError, "zero trades"):
            runner.daily_observation("BTCUSDT", day_ms, bars)

    def test_all_six_gates_must_pass_for_window_promotion(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])

        def load(_inventory: dict, asset: str, _start: int,
                 _end: int) -> dict:
            return synthetic_days(start_ms, 120 if asset in runner.ASSETS else 0)

        with mock.patch.object(runner, "load_asset_days", side_effect=load):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "pass-on-this-sample")
        self.assertEqual(set(result["gate_results"]), {
            "G1_count_magnitude", "G2_count_asset_consistency",
            "G3_absret_magnitude", "G4_absret_asset_consistency",
            "G5_month_stability", "G6_dependence_robustness",
        })
        self.assertTrue(all(
            gate["pass"] for gate in result["gate_results"].values()))
        self.assertEqual(
            result["gate_results"]["G6_dependence_robustness"]
            ["one_sided_block_bootstrap_p_count"],
            1 / (runner.BOOTSTRAP_REPLICATES + 1),
        )

    def test_one_failed_gate_falsifies_window(self) -> None:
        contract = frozen_test_contract()
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])
        q_trades = {"BTCUSDT": 130, "ETHUSDT": 100}

        def load(_inventory: dict, asset: str, _start: int,
                 _end: int) -> dict:
            return synthetic_days(start_ms, q_trades[asset])

        with mock.patch.object(runner, "load_asset_days", side_effect=load):
            result = runner.run_window(contract, "primary", {})
        failed = {
            name for name, gate in result["gate_results"].items()
            if not gate["pass"]
        }
        self.assertEqual(result["status"], "falsified")
        self.assertEqual(failed, {"G2_count_asset_consistency"})

    def test_insufficient_sample_is_inconclusive_without_gate_evaluation(
            self) -> None:
        contract = frozen_test_contract()
        contract["minimum_sample"]["primary"].update({
            "paired_dates": 22, "per_asset": 22, "pooled": 44})
        start_ms = runner.parse_iso_ms(
            contract["windows"]["primary"]["start_utc_inclusive"])

        with mock.patch.object(
                runner, "load_asset_days",
                side_effect=lambda *_args: synthetic_days(start_ms, 120)):
            result = runner.run_window(contract, "primary", {})
        self.assertEqual(result["status"], "inconclusive")
        self.assertIsNone(result["gate_results"])

    def test_unknown_daily_rejection_is_an_integrity_error(self) -> None:
        with self.assertRaisesRegex(ValueError, "internal invariant"):
            runner._rejection_key(ValueError("internal invariant"))

    def test_block_bootstrap_is_deterministic_and_date_paired(self) -> None:
        rows = []
        for index in range(28):
            date = f"2025-01-{index + 1:02d}"
            rows.extend([
                {"date": date, "asset": "BTCUSDT", "x_count": 0.2,
                 "x_absret": 0.1},
                {"date": date, "asset": "ETHUSDT", "x_count": 0.1,
                 "x_absret": 0.05},
            ])
        first = runner.block_bootstrap_p(
            rows, "x_count", replicates=100, seed=7, block_dates=7)
        second = runner.block_bootstrap_p(
            rows, "x_count", replicates=100, seed=7, block_dates=7)
        self.assertEqual(first, second)
        picked = runner.draw_block_dates(
            sorted({row["date"] for row in rows}), seed=11, block_dates=7)
        self.assertEqual(len(picked), 28)
        self.assertTrue(all(
            {row["asset"] for row in rows if row["date"] == date}
            == set(runner.ASSETS) for date in picked))
        self.assertEqual(first, 1 / 101)
        zero_rows = [{**row, "x_count": 0.0} for row in rows]
        self.assertEqual(
            runner.block_bootstrap_p(
                zero_rows, "x_count", replicates=100, seed=7, block_dates=7),
            1.0,
        )
        mixed_rows = [
            {**row, "x_count": 1.0 if row["asset"] == "BTCUSDT" else -2.0}
            for row in rows
        ]
        self.assertEqual(
            runner.block_bootstrap_p(
                mixed_rows, "x_count", replicates=10, seed=7, block_dates=7),
            1.0,
        )
        sample_lengths: list[int] = []

        def capture(values) -> float:
            values = list(values)
            sample_lengths.append(len(values))
            return 1.0

        with mock.patch.object(runner, "median", side_effect=capture):
            observed = runner.block_bootstrap_p(
                rows, "x_count", replicates=3, seed=7, block_dates=7)
        self.assertEqual(sample_lengths, [56, 56, 56])
        self.assertEqual(observed, 1 / 4)

    def test_timestamp_normalization_and_parser_pin_both_units(self) -> None:
        realistic_ms = 1_735_690_500_000
        realistic_us = realistic_ms * 1000
        self.assertEqual(runner.normalize_epoch(str(realistic_ms)), realistic_ms)
        self.assertEqual(runner.normalize_epoch(str(realistic_us)), realistic_ms)
        boundary = runner.MICROSECOND_THRESHOLD
        self.assertEqual(runner.normalize_epoch(str(boundary)), boundary)
        self.assertEqual(
            runner.normalize_epoch(str(boundary + 1)), (boundary + 1) // 1000)

        day_ms = 1_735_689_600_000
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, raw_timestamp in (
                    ("milliseconds.zip", realistic_ms),
                    ("microseconds.zip", realistic_us)):
                archive = root / name
                row = [
                    str(raw_timestamp), "100", "101", "99", "101", "1",
                    str(raw_timestamp), "1", "10", "1", "1", "0",
                ]
                with zipfile.ZipFile(archive, "w") as handle:
                    handle.writestr(
                        "bars.csv",
                        ",".join(runner.EXPECTED_HEADER) + "\n"
                        + ",".join(row) + "\n")
                days: dict[int, dict[int, dict]] = {}
                with mock.patch.object(runner, "WORKSPACE_ROOT", root):
                    runner._parse_archive(
                        {"path": name, "zip_member": "bars.csv"},
                        day_ms, day_ms + DAY_MS, set(), days)
                self.assertIn(realistic_ms, days[day_ms])


    def test_archive_verification_checks_sidecar_crc_and_member(self) -> None:
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
                "zip_member": "sample.csv", "header_present": True,
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
        }
        primary = {**expected, "role": "primary",
                   "result": {"status": "falsified"}}
        with self.assertRaisesRegex(ValueError, "primary pass"):
            runner.validate_primary_artifact(primary, expected)
        primary["result"]["status"] = "pass-on-this-sample"
        runner.validate_primary_artifact(primary, expected)
        primary["loop_policy_sha256"] = "different"
        with self.assertRaisesRegex(ValueError, "loop_policy_sha256"):
            runner.validate_primary_artifact(primary, expected)

    def test_write_json_new_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            runner.write_json_new(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                runner.write_json_new(path, {"a": 2})
            self.assertEqual(json.loads(path.read_text()), {"a": 1})

    def test_static_integrity_failure_prevents_outcome_read_and_artifact(
            self) -> None:
        contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            with mock.patch.object(
                    runner, "verify_static_inputs",
                    side_effect=ValueError("drift")) as verify, \
                    mock.patch.object(runner, "run_window") as run:
                with self.assertRaisesRegex(SystemExit, "static integrity"):
                    runner.main([
                        "--contract-sha256", contract_sha, "--window", "primary",
                        "--output", str(output),
                    ])
                verify.assert_called_once()
                run.assert_not_called()
                self.assertFalse(output.exists())

    def test_schema_failure_during_outcome_read_writes_no_artifact(self) -> None:
        contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        integrity = {
            "input_inventory_sha256": "i",
            "prior_window_audit_sha256": "a",
            "loop_policy_sha256": "p",
            "archives_verified": 20,
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            with mock.patch.object(
                    runner, "verify_static_inputs",
                    return_value=({}, integrity)), \
                    mock.patch.object(
                        runner, "run_window",
                        side_effect=ValueError("duplicate timestamp")):
                with self.assertRaisesRegex(SystemExit, "outcome input/schema"):
                    runner.main([
                        "--contract-sha256", contract_sha, "--window", "primary",
                        "--output", str(output),
                    ])
                self.assertFalse(output.exists())

    def test_existing_output_is_refused_before_integrity_or_outcome_reads(self) -> None:
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


if __name__ == "__main__":
    unittest.main()
