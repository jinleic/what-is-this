#!/usr/bin/env python3
"""Behavioral tests for the event-count-stopped Deribit taker-GEX runner."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import replicate_accelerated_intraday as BASE
import run_deribit_gex_eventstop as RUNNER
from run_deribit_gex_eventstop import (
    BLOCK_MS,
    derive_window_end,
    floor_block,
    validate_feasibility,
    write_bytes_exclusive,
)


PARENT_CONTRACT = BASE.EXPLORATIONS / "accelerated-intraday-contract.json"
CORRECTED_CONTRACT = BASE.EXPLORATIONS / "accelerated-intraday-feasibility-corrected-contract.json"


def _trade(ts: int, **extra: object) -> dict:
    return {"record_type": "option_trade", "received_ts_ms": ts, "dedup_id": f"t{ts}", "trade": {"timestamp": ts}, **extra}


def _index(name: str, ts: int) -> dict:
    return {"record_type": "index_price", "received_ts_ms": ts, "dedup_id": f"{name}:{ts}", "data": {"index_name": name, "timestamp": ts}}


def _pin(path: Path) -> dict:
    return {"path": str(path), "sha256": BASE.sha256_file(path)}


class WindowDerivationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.start_ms = 4 * BLOCK_MS
        self.records = {
            "a/trades.jsonl": [_trade(self.start_ms + 10), _trade(9 * BLOCK_MS + 500)],
            "a/index.jsonl": [
                _index("btc_usd", self.start_ms),
                _index("btc_usd", 10 * BLOCK_MS + 1),
                _index("eth_usd", self.start_ms),
                _index("eth_usd", 8 * BLOCK_MS + BLOCK_MS - 1),
            ],
        }

    def test_end_is_last_block_boundary_covered_by_slowest_stream(self) -> None:
        report = derive_window_end(self.records, self.start_ms)
        self.assertEqual(report["common_cutoff_ms"], 9 * BLOCK_MS - 1)
        self.assertEqual(report["end_ms"], 8 * BLOCK_MS)
        self.assertTrue(report["all_required_streams_cover_end"])
        self.assertEqual(report["streams"]["eth_usd_index"]["first_at_or_after_end_ms"], 9 * BLOCK_MS - 1)

    def test_missing_stream_fails_closed(self) -> None:
        records = {
            "a/trades.jsonl": self.records["a/trades.jsonl"],
            "a/index.jsonl": [row for row in self.records["a/index.jsonl"] if row["data"]["index_name"] == "btc_usd"],
        }
        with self.assertRaisesRegex(ValueError, "eth_usd_index"):
            derive_window_end(records, self.start_ms)

    def test_prefix_ending_inside_first_block_fails_closed(self) -> None:
        records = {
            "a/trades.jsonl": [_trade(self.start_ms + 1)],
            "a/index.jsonl": [_index("btc_usd", self.start_ms + 2), _index("eth_usd", self.start_ms + 3)],
        }
        with self.assertRaisesRegex(ValueError, "before the first complete block"):
            derive_window_end(records, self.start_ms)

    def test_floor_block_is_idempotent_on_boundaries(self) -> None:
        self.assertEqual(floor_block(7 * BLOCK_MS), 7 * BLOCK_MS)
        self.assertEqual(floor_block(7 * BLOCK_MS + BLOCK_MS - 1), 7 * BLOCK_MS)


class FeasibilityTests(unittest.TestCase):
    def test_minimum_seven_blocks_reach_twelve_transitions(self) -> None:
        with self.assertRaisesRegex(ValueError, "unreachable"):
            validate_feasibility(0, 6 * BLOCK_MS, 12, not_before_ms=0)
        report = validate_feasibility(0, 7 * BLOCK_MS, 12, not_before_ms=0)
        self.assertEqual(report["max_complete_asset_block_transitions"], 12)
        self.assertTrue(report["registered_stop_reachable"])

    def test_start_must_not_precede_prior_window_end(self) -> None:
        with self.assertRaisesRegex(ValueError, "overlaps"):
            validate_feasibility(BLOCK_MS, 20 * BLOCK_MS, 12, not_before_ms=2 * BLOCK_MS)

    def test_boundaries_must_be_block_aligned(self) -> None:
        with self.assertRaisesRegex(ValueError, "aligned"):
            validate_feasibility(1, 20 * BLOCK_MS, 12, not_before_ms=0)


class ParentStopIsOutcomeBlindTests(unittest.TestCase):
    """The stop may depend on counts, signs, and RV availability, never RV magnitude."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.module = BASE.load_module("eventstop_test_parent_deribit", BASE.DERIBIT_MODULE)
        cls.sample_stop = {
            "minimum_gex_computable_trades": 4,
            "minimum_complete_asset_block_transitions": 4,
            "minimum_transitions_per_asset": 2,
            "minimum_positive_dealer_gex_blocks": 2,
            "minimum_negative_dealer_gex_blocks": 2,
        }

    def _transitions(self, rvs: list[float]) -> list[dict]:
        rows = []
        for block, ((underlying, gex), rv) in enumerate(
            zip([("BTC", 1.0), ("ETH", -1.0), ("BTC", -2.0), ("ETH", 3.0), ("BTC", 1.5), ("ETH", -0.5)], rvs)
        ):
            rows.append(
                {
                    "underlying": underlying,
                    "signal_block": block,
                    "outcome_end_ms": (block + 2) * BLOCK_MS,
                    "dealer_gex": gex,
                    "computable_trades_in_signal_block": 1,
                    "next_rv": rv,
                }
            )
        return rows

    def test_stop_is_invariant_to_rv_permutation_that_flips_every_gate(self) -> None:
        deadline = 10 * BLOCK_MS
        increasing = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        decreasing = increasing[::-1]
        gate_cfg = {
            "spearman_dealer_gex_vs_next_rv_lte": -0.25,
            "negative_vs_positive_dealer_gex_next_rv_ratio_gte": 1.1,
            "both_asset_spearman_lt_zero": True,
            "leave_one_block_out_nonpositive_fraction_gte": 0.7,
        }
        gates_up, _ = self.module.evaluate_gates(self._transitions(increasing)[:4], gate_cfg)
        gates_down, _ = self.module.evaluate_gates(self._transitions(decreasing)[:4], gate_cfg)
        self.assertNotEqual(
            gates_up["spearman_dealer_gex_vs_next_rv_lte"]["value"],
            gates_down["spearman_dealer_gex_vs_next_rv_lte"]["value"],
        )
        stop_up, trail_up = self.module.find_earliest_stop(self._transitions(increasing), self.sample_stop, deadline)
        stop_down, trail_down = self.module.find_earliest_stop(self._transitions(decreasing), self.sample_stop, deadline)
        self.assertEqual(stop_up, 5 * BLOCK_MS)
        self.assertEqual(stop_up, stop_down)
        self.assertEqual(trail_up, trail_down)

    def test_stop_respects_deadline(self) -> None:
        stop, trail = self.module.find_earliest_stop(self._transitions([1.0] * 6), self.sample_stop, 4 * BLOCK_MS)
        self.assertIsNone(stop)
        self.assertEqual(len(trail), 3)


class MultiDayReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = BASE.load_module("eventstop_test_parent_deribit_replay", BASE.DERIBIT_MODULE)

    def test_cross_file_duplicate_survives_once_in_receipt_order(self) -> None:
        early = _trade(100)
        late = dict(_trade(100), received_ts_ms=200)
        records = {"x/day2/trades.jsonl": [late], "x/day1/trades.jsonl": [early]}
        combined = BASE.records_by_filename(records, "trades.jsonl")
        self.assertEqual([row["received_ts_ms"] for row in combined], [100, 200])
        deduped, drops = self.module.dedupe(self.module.replay_order(combined))
        self.assertEqual(drops, 1)
        self.assertEqual(len(deduped), 1)
        self.assertEqual(deduped[0]["received_ts_ms"], 100)

    def test_metadata_received_at_end_is_ineligible(self) -> None:
        end_ms = 10 * BLOCK_MS
        rows = {
            "trades": [],
            "index": [],
            "instruments": [
                {"received_ts_ms": end_ms - 1, "instrument": {"instrument_name": "A"}},
                {"received_ts_ms": end_ms, "instrument": {"instrument_name": "B"}},
            ],
        }
        kept, counts = BASE.filter_deribit_records(rows, 0, end_ms)
        self.assertEqual([row["instrument"]["instrument_name"] for row in kept["instruments"]], ["A"])
        self.assertEqual(counts["excluded_at_or_after_end"]["instruments"], 1)


class PublicationTests(unittest.TestCase):
    def test_exclusive_write_never_replaces(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "out.json"
            digest = write_bytes_exclusive(target, b"{}\n")
            self.assertEqual(digest, BASE.sha256_file(target))
            with self.assertRaises(FileExistsError):
                write_bytes_exclusive(target, b"[]\n")
            self.assertEqual(target.read_bytes(), b"{}\n")


class FailClosedLifecycleTests(unittest.TestCase):
    """A precondition failure must leave every trial output absent."""

    def _write_jsonl(self, path: Path, rows: list[dict]) -> dict:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = "".join(json.dumps(row) + "\n" for row in rows).encode("utf-8")
        path.write_bytes(payload)
        return {"path": str(path), "prefix_bytes": len(payload), "sha256": BASE.sha256_bytes(payload)}

    def _fixture(self, tmp: Path, *, with_eth: bool) -> tuple[Path, Path]:
        parent = json.loads(PARENT_CONTRACT.read_text(encoding="utf-8"))
        corrected = json.loads(CORRECTED_CONTRACT.read_text(encoding="utf-8"))
        start_ms = BASE.parse_iso_ms(corrected["observation_window"]["end_utc_exclusive"])
        index_rows = [_index("btc_usd", start_ms + 1), _index("btc_usd", start_ms + 30 * BLOCK_MS)]
        if with_eth:
            index_rows += [_index("eth_usd", start_ms + 1), _index("eth_usd", start_ms + 30 * BLOCK_MS)]
        prefixes = [
            self._write_jsonl(tmp / "day1" / "trades.jsonl", [_trade(start_ms + 5)]),
            self._write_jsonl(tmp / "day1" / "index.jsonl", index_rows),
            self._write_jsonl(
                tmp / "day1" / "instruments.jsonl",
                [{"record_type": "instrument", "received_ts_ms": start_ms, "dedup_id": "i1", "instrument": {"instrument_name": "X", "creation_timestamp": 1}}],
            ),
        ]
        run_dir = tmp / "run"
        run_dir.mkdir()
        contract = {
            "status": RUNNER.CONTRACT_STATUS,
            "created_at_utc": "2026-09-02T00:00:00.000Z",
            "classification": "test",
            "claim_boundary": "test",
            "runner": _pin(RUNNER.RUNNER_PATH),
            "lineage": {"parent_contract": _pin(PARENT_CONTRACT), "corrected_contract": _pin(CORRECTED_CONTRACT)},
            "governance": {},
            "stopping_rule": {},
            "screen_definition": parent["screens"][RUNNER.SCREEN],
            "observation_window": {"start_utc_inclusive": corrected["observation_window"]["end_utc_exclusive"]},
            "frozen_input_prefixes": prefixes,
            "frozen_rule_files": {},
            "baseline_artifact_hashes": {},
            "governing_authority": {},
            "outputs": {
                "snapshot_filename": "snapshot.json",
                "results_filename": "results.json",
                "owner_dir_copy": str(tmp / "owner" / "copy.json"),
            },
        }
        contract_path = tmp / "contract.json"
        contract_path.write_bytes(RUNNER.encode_json(contract))
        (run_dir / "manifest.json").write_text(
            json.dumps({"run_id": "test", "prereg_sha256": BASE.sha256_file(contract_path), "status": "RUNNING"}),
            encoding="utf-8",
        )
        return contract_path, run_dir

    def _assert_no_outputs(self, tmp: Path, run_dir: Path) -> None:
        self.assertFalse((run_dir / "snapshot.json").exists())
        self.assertFalse((run_dir / "results.json").exists())
        self.assertFalse((tmp / "owner" / "copy.json").exists())

    def test_missing_stream_aborts_before_any_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            contract_path, run_dir = self._fixture(tmp, with_eth=False)
            with self.assertRaisesRegex(ValueError, "eth_usd_index"):
                RUNNER.main(["--contract", str(contract_path), "--run-dir", str(run_dir)])
            self._assert_no_outputs(tmp, run_dir)

    def test_existing_snapshot_refuses_before_reading_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            contract_path, run_dir = self._fixture(tmp, with_eth=True)
            (run_dir / "snapshot.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "CRASHED"):
                RUNNER.main(["--contract", str(contract_path), "--run-dir", str(run_dir)])
            self.assertFalse((run_dir / "results.json").exists())

    def test_manifest_must_bind_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            contract_path, run_dir = self._fixture(tmp, with_eth=True)
            (run_dir / "manifest.json").write_text(
                json.dumps({"run_id": "test", "prereg_sha256": "0" * 64, "status": "RUNNING"}), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "prereg_sha256"):
                RUNNER.main(["--contract", str(contract_path), "--run-dir", str(run_dir)])
            self._assert_no_outputs(tmp, run_dir)

    def test_changed_prefix_aborts_before_any_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            contract_path, run_dir = self._fixture(tmp, with_eth=True)
            with (tmp / "day1" / "trades.jsonl").open("r+b") as handle:
                handle.seek(0)
                handle.write(b" ")
            with self.assertRaisesRegex(ValueError, "prefix hash mismatch"):
                RUNNER.main(["--contract", str(contract_path), "--run-dir", str(run_dir)])
            self._assert_no_outputs(tmp, run_dir)


if __name__ == "__main__":
    unittest.main()
