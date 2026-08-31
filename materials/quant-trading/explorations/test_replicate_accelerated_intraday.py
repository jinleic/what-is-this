#!/usr/bin/env python3
"""Boundary and integrity tests for the frozen chronological replication."""

from __future__ import annotations

import unittest

import replicate_accelerated_intraday as replication


class ReplicationWindowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.start = 1_000_000
        self.end = 1_900_000

    def test_half_open_window_includes_start_and_excludes_end(self) -> None:
        self.assertFalse(replication.in_half_open(self.start - 1, self.start, self.end))
        self.assertTrue(replication.in_half_open(self.start, self.start, self.end))
        self.assertTrue(replication.in_half_open(self.end - 1, self.start, self.end))
        self.assertFalse(replication.in_half_open(self.end, self.start, self.end))

    def test_liquidation_filter_keeps_pre_window_context_but_not_post_window(self) -> None:
        rows = [
            {"kind": "control", "ts_utc": "1970-01-01T00:15:00.000Z", "control": {"type": "connect"}},
            {"kind": "event", "payload": {"E": self.start - 1}},
            {"kind": "event", "payload": {"E": self.start}},
            {"kind": "event", "payload": {"E": self.end}},
        ]
        filtered = replication.filter_liquidation_context(rows, self.end)
        self.assertEqual(filtered, rows[:3])

    def test_deribit_filter_applies_window_before_metric_code(self) -> None:
        records = {
            "trades": [
                {"trade": {"timestamp": self.start - 1}},
                {"trade": {"timestamp": self.start}},
                {"trade": {"timestamp": self.end}},
            ],
            "index": [
                {"data": {"timestamp": self.start - 5_000}},
                {"data": {"timestamp": self.start}},
                {"data": {"timestamp": self.end}},
            ],
            "instruments": [
                {"received_ts_ms": self.end - 1},
                {"received_ts_ms": self.end},
            ],
        }
        filtered, counts = replication.filter_deribit_records(
            records, self.start, self.end
        )
        self.assertEqual([r["trade"]["timestamp"] for r in filtered["trades"]], [self.start])
        self.assertEqual(
            [r["data"]["timestamp"] for r in filtered["index"]],
            [self.start - 5_000, self.start],
        )
        self.assertEqual(len(filtered["instruments"]), 1)
        self.assertEqual(counts["excluded_at_or_after_end"]["trades"], 1)
        self.assertEqual(counts["excluded_at_or_after_end"]["index"], 1)
        self.assertEqual(counts["excluded_at_or_after_end"]["instruments"], 1)

    def test_liquidation_onset_membership_is_half_open(self) -> None:
        onsets = [
            {"onset_ms": self.start - 1},
            {"onset_ms": self.start},
            {"onset_ms": self.end - 1},
            {"onset_ms": self.end},
        ]
        self.assertEqual(
            replication.filter_liquidation_onsets(onsets, self.start, self.end),
            onsets[1:3],
        )


if __name__ == "__main__":
    unittest.main()
