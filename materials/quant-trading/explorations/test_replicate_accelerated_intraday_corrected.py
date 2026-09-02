#!/usr/bin/env python3
"""Behavioral tests for the feasibility-corrected chronological runner."""

from __future__ import annotations

import unittest

from replicate_accelerated_intraday_corrected import (
    alignment_provenance,
    instrument_clock_provenance,
    validate_common_cutoffs,
    validate_window_feasibility,
)


BLOCK_MS = 15 * 60 * 1000


class FeasibilitySafeguardTests(unittest.TestCase):
    def test_ninety_minutes_is_rejected_and_105_minutes_is_reachable(self) -> None:
        start_ms = 0
        with self.assertRaisesRegex(ValueError, "structurally unreachable"):
            validate_window_feasibility(start_ms, 6 * BLOCK_MS, 6, 12)

        report = validate_window_feasibility(start_ms, 7 * BLOCK_MS, 6, 12)
        self.assertEqual(report["window_blocks"], 7)
        self.assertEqual(report["max_complete_liquidation_onset_blocks"], 6)
        self.assertEqual(report["max_complete_deribit_asset_block_transitions"], 12)
        self.assertTrue(report["all_registered_stops_reachable"])

    def test_window_must_not_overlap_prior_replication(self) -> None:
        with self.assertRaisesRegex(ValueError, "overlaps"):
            validate_window_feasibility(
                0, 7 * BLOCK_MS, 6, 12, not_before_ms=BLOCK_MS
            )

        report = validate_window_feasibility(
            BLOCK_MS, 8 * BLOCK_MS, 6, 12, not_before_ms=BLOCK_MS
        )
        self.assertTrue(report["nonoverlap_with_prior_window"])
        self.assertEqual(report["not_before_ms"], BLOCK_MS)

    def test_window_must_be_aligned_to_parent_blocks(self) -> None:
        with self.assertRaisesRegex(ValueError, "aligned"):
            validate_window_feasibility(1, 7 * BLOCK_MS + 1, 6, 12)


class CoverageSafeguardTests(unittest.TestCase):
    def setUp(self) -> None:
        self.end_ms = 10_000
        self.records = {
            "liquidation-capture.jsonl": [
                {"kind": "event", "payload": {"data": {"E": self.end_ms + 7}}}
            ],
            "trades.jsonl": [
                {"trade": {"timestamp": self.end_ms + 11}}
            ],
            "index.jsonl": [
                {"data": {"index_name": "btc_usd", "timestamp": self.end_ms}},
                {"data": {"index_name": "eth_usd", "timestamp": self.end_ms + 3}},
            ],
        }

    def test_common_cutoff_requires_every_stream_through_end(self) -> None:
        report = validate_common_cutoffs(self.records, self.end_ms)
        self.assertTrue(report["all_required_streams_cover_end"])
        self.assertEqual(report["common_cutoff_ms"], self.end_ms)
        self.assertEqual(report["streams"]["btc_usd_index"]["first_at_or_after_end_ms"], self.end_ms)

    def test_missing_or_truncated_stream_fails_closed(self) -> None:
        missing_eth = {name: list(rows) for name, rows in self.records.items()}
        missing_eth["index.jsonl"] = missing_eth["index.jsonl"][:1]
        with self.assertRaisesRegex(ValueError, "eth_usd_index"):
            validate_common_cutoffs(missing_eth, self.end_ms)

        truncated_trade = {name: list(rows) for name, rows in self.records.items()}
        truncated_trade["trades.jsonl"] = [{"trade": {"timestamp": self.end_ms - 1}}]
        with self.assertRaisesRegex(ValueError, "deribit_option_trade"):
            validate_common_cutoffs(truncated_trade, self.end_ms)


class ProvenanceSafeguardTests(unittest.TestCase):
    def test_alignment_provenance_preserves_exact_missing_reason(self) -> None:
        outcomes = [
            {
                "event_id": "event-1",
                "symbol": "BTCUSDT",
                "side": "BUY",
                "onset_ms": 100,
                "cohort_complete": False,
                "P_0": {"observed_ts_ms": 99, "price": 10.0},
                "P_1m": {"target_ms": 160, "observed_ts_ms": 161, "price": 10.1},
                "P_5m": {"target_ms": 400, "observed_ts_ms": 401, "price": 10.2},
                "P_15m": {"target_ms": 1000, "observed": None, "missing": "beyond_common_cutoff"},
            }
        ]
        report = alignment_provenance(outcomes)
        self.assertEqual(report["incomplete_alignment_counts"]["horizon_15m"], 1)
        self.assertEqual(report["incomplete_alignment_counts"]["initial"], 0)
        self.assertEqual(
            report["incomplete_onsets"][0]["alignments"]["P_15m"]["missing"],
            "beyond_common_cutoff",
        )
        self.assertNotIn("missing", report["incomplete_onsets"][0]["alignments"]["P_5m"])

    def test_instrument_exchange_and_receipt_clocks_remain_separate(self) -> None:
        report = instrument_clock_provenance(
            [
                {
                    "received_ts_ms": 500,
                    "instrument": {"creation_timestamp": 100},
                },
                {
                    "received_ts_ms": 700,
                    "instrument": {"creation_timestamp": 300},
                },
            ]
        )
        self.assertEqual(report["exchange_clock"]["field"], "instrument.creation_timestamp")
        self.assertEqual(report["exchange_clock"]["latest_ms"], 300)
        self.assertEqual(report["capture_availability_clock"]["field"], "received_ts_ms")
        self.assertEqual(report["capture_availability_clock"]["latest_ms"], 700)


if __name__ == "__main__":
    unittest.main()
