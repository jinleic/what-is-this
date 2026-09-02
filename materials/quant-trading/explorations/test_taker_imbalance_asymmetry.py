#!/usr/bin/env python3
"""Behavioral tests for the taker-imbalance / semivariance-asymmetry runner."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import run_taker_imbalance_asymmetry as runner


def synthetic(decision: int, spot_taker: float, perp_taker: float,
              forward: list[float], anchor_taker: float | None = None,
              omit_signal_hour: int | None = None) -> tuple[dict, dict]:
    """Trailing bars t-25h..t-2h carry the signal; t-1h is the anchor only."""
    hour = runner.HOUR_MS
    spot: dict[int, dict] = {}
    perp: dict[int, dict] = {}
    for k in range(2, 26):  # signal bars t-25h .. t-2h
        t = decision - k * hour
        if omit_signal_hour is not None and k == omit_signal_hour:
            continue
        spot[t] = {"close": 100.0, "quote_volume": 10.0,
                   "taker_buy_quote_volume": 10.0 * spot_taker}
        perp[t] = {"close": 100.0, "quote_volume": 20.0,
                   "taker_buy_quote_volume": 20.0 * perp_taker}
    # anchor bar t-1h: close feeds y; its volumes must never feed x
    anchor_share = spot_taker if anchor_taker is None else anchor_taker
    spot[decision - hour] = {"close": 100.0, "quote_volume": 1_000_000.0,
                             "taker_buy_quote_volume": 1_000_000.0 * anchor_share}
    perp[decision - hour] = {"close": 100.0, "quote_volume": 1_000_000.0,
                             "taker_buy_quote_volume": 1_000_000.0 * anchor_share}
    price = 100.0
    for k, ret in enumerate(forward):  # forward bars t .. t+23h
        price *= math.exp(ret)
        spot[decision + k * hour] = {"close": price, "quote_volume": 10.0,
                                     "taker_buy_quote_volume": 5.0}
    return spot, perp


FORWARD = ([0.01, -0.02] * 6) + ([0.005, -0.005] * 6)


class TakerImbalanceTests(unittest.TestCase):
    def test_signal_is_signed_share_difference(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.40, 0.60, FORWARD)
        obs = runner.daily_observation("BTCUSDT", decision, spot, perp)
        self.assertAlmostEqual(obs["x"], 0.20, places=12)
        self.assertAlmostEqual(obs["perp_taker_share"], 0.60, places=12)
        self.assertAlmostEqual(obs["spot_taker_share"], 0.40, places=12)

    def test_outcome_uses_only_forward_bars_and_is_a_ratio_of_semivariances(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.5, 0.5, FORWARD)
        obs = runner.daily_observation("BTCUSDT", decision, spot, perp)
        down = sum(r * r for r in FORWARD if r < 0)
        up = sum(r * r for r in FORWARD if r > 0)
        self.assertAlmostEqual(obs["y"], math.log(down / up), places=12)
        self.assertEqual(obs["forward_return_count"], 24)
        self.assertEqual(obs["outcome_first_return_ms"], decision)
        self.assertEqual(obs["outcome_last_return_ms"], decision + 23 * runner.HOUR_MS)

    def test_one_sided_forward_window_is_rejected(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.5, 0.5, [0.001] * 24)
        with self.assertRaisesRegex(ValueError, "insufficient signed returns"):
            runner.daily_observation("BTCUSDT", decision, spot, perp)

    def test_missing_perp_signal_hour_fails_closed(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.5, 0.5, FORWARD)
        del perp[decision - 5 * runner.HOUR_MS]
        with self.assertRaisesRegex(ValueError, "missing required perp hour"):
            runner.daily_observation("BTCUSDT", decision, spot, perp)

    def test_anchor_bar_volumes_never_enter_the_signal(self) -> None:
        """Embargo hour t-1h carries a million-unit lopsided taker print; if it
        leaked into x, the share difference would move away from 0.20."""
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.40, 0.60, FORWARD, anchor_taker=1.0)
        obs = runner.daily_observation("BTCUSDT", decision, spot, perp)
        self.assertAlmostEqual(obs["x"], 0.20, places=12)
        self.assertEqual(obs["signal_last_ms"], decision - 2 * runner.HOUR_MS)
        self.assertEqual(obs["signal_first_ms"], decision - 25 * runner.HOUR_MS)

    def test_missing_anchor_bar_fails_closed(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.5, 0.5, FORWARD)
        del spot[decision - runner.HOUR_MS]
        with self.assertRaisesRegex(ValueError, "missing required spot hour"):
            runner.daily_observation("BTCUSDT", decision, spot, perp)

    def test_earliest_signal_bar_is_required(self) -> None:
        decision = 1_800_000_000_000
        spot, perp = synthetic(decision, 0.5, 0.5, FORWARD, omit_signal_hour=25)
        with self.assertRaisesRegex(ValueError, "missing required spot hour"):
            runner.daily_observation("BTCUSDT", decision, spot, perp)

    def test_quintile_tie_break_is_deterministic(self) -> None:
        rows = [{"x": 0.0, "y": float(i), "date": f"d{i:02d}",
                 "asset": "BTCUSDT" if i % 2 else "ETHUSDT"} for i in range(20)]
        first = runner.quintile_gap(rows)
        second = runner.quintile_gap(list(reversed(rows)))
        self.assertEqual(first, second)

    def test_taker_volume_exceeding_quote_volume_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "taker_buy_quote_volume"):
            runner.validate_bar("x.zip", 0, close=100.0, quote_volume=10.0,
                                taker_buy_quote_volume=10.5)
        runner.validate_bar("x.zip", 0, close=100.0, quote_volume=10.0,
                            taker_buy_quote_volume=10.0)

    def test_bootstrap_is_deterministic_and_exact_length(self) -> None:
        seen: list[int] = []
        real = runner.pearson

        def spy(xs, ys):
            seen.append(len(xs))
            return real(xs, ys)

        rows = [
            {"date": f"d{day:02d}", "asset": asset, "x": float(day), "y": float(day)}
            for day in range(30)
            for asset in ("BTCUSDT", "ETHUSDT")
        ]
        runner.pearson = spy
        try:
            first = runner.block_bootstrap_p(rows, replicates=50, seed=5, block_dates=7)
            second = runner.block_bootstrap_p(rows, replicates=50, seed=5, block_dates=7)
        finally:
            runner.pearson = real
        self.assertEqual(first, second)
        self.assertEqual(set(seen), {60})

    def test_replication_requires_committed_primary_pass(self) -> None:
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

    def test_write_json_new_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r.json"
            runner.write_json_new(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                runner.write_json_new(path, {"a": 2})


if __name__ == "__main__":
    unittest.main()
