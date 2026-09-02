#!/usr/bin/env python3
"""Behavioral tests for the frozen activity/volatility runner."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
import zipfile
from pathlib import Path

import run_perp_spot_activity_volatility as runner


class ActivityVolatilityTests(unittest.TestCase):
    def test_epoch_normalization_handles_ms_and_us(self) -> None:
        self.assertEqual(runner.normalize_epoch("1735689600000"), 1735689600000)
        self.assertEqual(runner.normalize_epoch("1735689600000000"), 1735689600000)

    def test_daily_observation_uses_only_frozen_hours(self) -> None:
        hour = 3_600_000
        decision = 1_800_000_000_000
        times = [decision + k * hour for k in range(-25, 24)]
        spot = {
            t: {"close": 100.0 * math.exp(0.001 * i), "quote_volume": 10.0}
            for i, t in enumerate(times)
        }
        perp = {
            t: {"close": 100.0, "quote_volume": 20.0}
            for t in times
        }
        obs = runner.daily_observation("BTCUSDT", decision, spot, perp)
        self.assertAlmostEqual(obs["x"], math.log(2.0), places=12)
        self.assertEqual(obs["signal_first_ms"], decision - 24 * hour)
        self.assertEqual(obs["signal_last_ms"], decision - hour)
        self.assertEqual(obs["outcome_first_return_ms"], decision)
        self.assertEqual(obs["outcome_last_return_ms"], decision + 23 * hour)
        self.assertEqual(obs["trailing_return_count"], 24)
        self.assertEqual(obs["future_return_count"], 24)

    def test_missing_required_hour_fails_closed(self) -> None:
        hour = 3_600_000
        decision = 1_800_000_000_000
        times = [decision + k * hour for k in range(-25, 24)]
        spot = {t: {"close": 100.0, "quote_volume": 10.0} for t in times}
        perp = {t: {"close": 100.0, "quote_volume": 20.0} for t in times}
        del spot[decision + 7 * hour]
        with self.assertRaisesRegex(ValueError, "missing required spot hour"):
            runner.daily_observation("BTCUSDT", decision, spot, perp)

    def test_average_rank_spearman_handles_ties(self) -> None:
        self.assertAlmostEqual(runner.spearman([1, 1, 3], [2, 2, 5]), 1.0)
        self.assertAlmostEqual(runner.spearman([1, 2, 3], [3, 2, 1]), -1.0)

    def test_date_block_bootstrap_is_deterministic_and_paired(self) -> None:
        rows = []
        for day in range(20):
            for asset_i, asset in enumerate(("BTCUSDT", "ETHUSDT")):
                x = float(day * 2 + asset_i)
                rows.append({"date": f"d{day:02d}", "asset": asset, "x": x, "y": x})
        first = runner.block_bootstrap_p(rows, replicates=200, seed=7, block_dates=3)
        second = runner.block_bootstrap_p(rows, replicates=200, seed=7, block_dates=3)
        self.assertEqual(first, second)
        self.assertLess(first, 0.05)

    def test_write_json_new_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            runner.write_json_new(path, {"a": 1})
            self.assertEqual(json.loads(path.read_text()), {"a": 1})
            with self.assertRaises(FileExistsError):
                runner.write_json_new(path, {"a": 2})

    def test_bootstrap_replicate_length_equals_observed_dates(self) -> None:
        seen: list[int] = []
        real_pearson = runner.pearson

        def spy(xs, ys):
            seen.append(len(xs))
            return real_pearson(xs, ys)

        rows = [
            {"date": f"d{day:02d}", "asset": asset, "x": float(day), "y": float(day)}
            for day in range(30)
            for asset in ("BTCUSDT", "ETHUSDT")
        ]
        runner.pearson = spy
        try:
            runner.block_bootstrap_p(rows, replicates=25, seed=3, block_dates=7)
        finally:
            runner.pearson = real_pearson
        # 30 dates x 2 paired assets; 7-date blocks must not inflate n past 60
        self.assertTrue(seen)
        self.assertEqual(set(seen), {60})

    def test_replication_refused_without_committed_primary_pass(self) -> None:
        contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        with tempfile.TemporaryDirectory() as directory:
            primary = Path(directory) / "primary.json"
            primary.write_text(json.dumps({
                "role": "primary",
                "contract_sha256": contract_sha,
                "result": {"status": "falsified"},
            }), encoding="utf-8")
            with self.assertRaises(SystemExit) as ctx:
                runner.main([
                    "--contract-sha256", contract_sha,
                    "--window", "replication",
                    "--output", str(Path(directory) / "rep.json"),
                    "--primary-artifact", str(primary),
                ])
            self.assertIn("only after a committed primary pass", str(ctx.exception))

    def test_replication_refused_when_primary_artifact_absent(self) -> None:
        contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit) as ctx:
                runner.main([
                    "--contract-sha256", contract_sha,
                    "--window", "replication",
                    "--output", str(Path(directory) / "rep.json"),
                    "--primary-artifact", str(Path(directory) / "absent.json"),
                ])
            self.assertIn("committed primary artifact absent", str(ctx.exception))

    def test_single_member_archive_is_enforced(self) -> None:
        source = Path(
            "../data/raw/binance/spot-1h/BTCUSDT/BTCUSDT-1h-2025-07.zip"
        ).resolve()
        self.assertTrue(source.is_file())
        with zipfile.ZipFile(source) as archive:
            self.assertEqual(len(archive.namelist()), 1)


if __name__ == "__main__":
    unittest.main()
