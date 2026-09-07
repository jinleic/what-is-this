#!/usr/bin/env python3
"""Regression tests for the live public-data collectors."""

from __future__ import annotations

import asyncio
from collections import Counter, deque
import importlib.util
import io
import json
import os
import subprocess
import signal
import socket
import tempfile
import sys
import time
import unittest
import urllib.error
from unittest import mock
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


deribit = load_module("collector_integrity_deribit", HERE / "focused-deribit-gex" / "capture.py")
liquidation = load_module(
    "collector_integrity_liquidation",
    HERE / "focused-liquidation-overlay" / "capture.py",
)

def start_output_lock_holder(lock_path: Path) -> subprocess.Popen:
    script = (
        "import fcntl, pathlib, sys, time;"
        "path = pathlib.Path(sys.argv[1]);"
        "path.parent.mkdir(parents=True, exist_ok=True);"
        "handle = path.open('a+');"
        "fcntl.flock(handle.fileno(), fcntl.LOCK_EX);"
        "print('LOCKED', flush=True);"
        "time.sleep(60)"
    )
    holder = subprocess.Popen(
        [sys.executable, "-c", script, str(lock_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if holder.stdout is None or holder.stdout.readline().strip() != "LOCKED":
        stderr = holder.stderr.read() if holder.stderr is not None else ""
        holder.wait(timeout=5)
        if holder.stdout is not None:
            holder.stdout.close()
        if holder.stderr is not None:
            holder.stderr.close()
        raise RuntimeError(f"lock holder failed: {stderr}")
    return holder


def stop_process(process: subprocess.Popen) -> None:
    if process.poll() is None:
        process.terminate()
        process.wait(timeout=5)
    if process.stdout is not None:
        process.stdout.close()
    if process.stderr is not None:
        process.stderr.close()


class ScriptedDeribitWs:
    def __init__(self, messages: list[dict]):
        self.messages = list(messages)
        self.sent: list[dict] = []

    def send_text(self, text: str) -> None:
        self.sent.append(json.loads(text))

    def recv_message(self, timeout: float):
        del timeout
        if not self.messages:
            raise socket.timeout()
        return deribit.OP_TEXT, json.dumps(self.messages.pop(0)).encode("utf-8")

    def close(self) -> None:
        pass


class FakeHttpResponse:
    """Stateful HTTP response with realistic incremental reads."""

    def __init__(
        self,
        chunks: list[bytes],
        *,
        status: int = 200,
        delay: float = 0.0,
    ) -> None:
        self.chunks = deque(chunks)
        self.status = status
        self.delay = delay
        self.reads = 0

    def read(self, size: int = -1) -> bytes:
        self.reads += 1
        if size is None or size < 0:
            body = bytearray()
            while self.chunks:
                if self.delay:
                    time.sleep(self.delay)
                body += self.chunks.popleft()
            return bytes(body)
        if not self.chunks:
            return b""
        if self.delay:
            time.sleep(self.delay)
        chunk = self.chunks.popleft()
        if len(chunk) > size:
            self.chunks.appendleft(chunk[size:])
            return chunk[:size]
        return chunk

    def __enter__(self) -> "FakeHttpResponse":
        return self

    def __exit__(self, _exc_type, _exc, _tb) -> bool:
        return False


class DeribitTickerIntegrityTests(unittest.TestCase):
    instrument_name = "BTC-1JAN27-90000-C"

    def make_collector(self):
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        # Tests of metadata provide explicit responses; other fixtures stay offline.
        fetch = mock.patch.object(
            deribit,
            "fetch_instruments",
            side_effect=ConnectionError("offline fixture has no metadata response"),
        )
        fetch.start()
        self.addCleanup(fetch.stop)

        def finish_refresh() -> None:
            thread = collector._chain_refresh_thread
            if thread is not None:
                thread.join(timeout=5.0)
                self.assertFalse(thread.is_alive(), "fixture refresh outlived its test")

        self.addCleanup(finish_refresh)
        collector.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        return collector, out

    def ticker_payload(self, timestamp: int = 1_800_000_000_000) -> dict:
        return {
            "timestamp": timestamp,
            "state": "open",
            "instrument_name": self.instrument_name,
            "index_price": 90_000.0,
            "underlying_price": 90_100.0,
            "open_interest": 12.5,
            "mark_price": 0.01,
            "mark_iv": 55.0,
            "greeks": {
                "delta": 0.5,
                "gamma": 0.00001,
                "theta": -10.0,
                "vega": 20.0,
                "rho": 5.0,
            },
        }
    def trade_payload(
        self,
        sequence: int,
        trade_id: str | None = None,
        received_ms: int | None = None,
    ) -> dict:
        return {
            "timestamp": (
                1_800_000_000_000 + sequence
                if received_ms is None
                else received_ms - 1
            ),
            "instrument_name": self.instrument_name,
            "trade_id": trade_id or f"trade-{sequence}",
            "trade_seq": sequence,
            "direction": "buy",
            "amount": 1.0,
            "price": 0.01,
            "index_price": 90_000.0,
            "iv": 55.0,
        }

    def instrument_payload(
        self,
        index: int = 0,
        *,
        currency: str = "btc",
        name: str | None = None,
        expiration_timestamp: int | None = None,
        state: str = "open",
        is_active: bool = True,
    ) -> dict:
        upper = currency.upper()
        instrument_name = (
            name
            or f"{upper}-30DEC30-{90_000 + index}-"
            f"{'C' if index % 2 == 0 else 'P'}"
        )
        if expiration_timestamp is None:
            identity = deribit.Collector._option_instrument_identity(
                instrument_name
            )
            if identity is None:
                raise ValueError(
                    f"test fixture name is not an option: {instrument_name}"
                )
            expiration_timestamp = int(
                identity.expiry.timestamp() * 1000
            ) + 8 * 60 * 60 * 1000
        return {
            "instrument_name": instrument_name,
            "instrument_id": index + 1,
            "kind": "option",
            "creation_timestamp": 1_700_000_000_000 + index,
            "expiration_timestamp": expiration_timestamp,
            "strike": float(90_000 + index),
            "option_type": "call" if index % 2 == 0 else "put",
            "base_currency": upper,
            "state": state,
            "is_active": is_active,
            "contract_size": 1.0,
            "min_trade_amount": 0.1 if upper == "BTC" else 1.0,
            "lot_size": 1 if upper == "BTC" else 10,
            "tick_size": 0.0001,
            "settlement_currency": upper,
            "price_index": f"{currency.lower()}_usd",
        }

    def restart_cursor(
        self,
        root: Path,
        path: Path,
        offset: int | None = None,
    ) -> dict:
        offset = path.stat().st_size if offset is None else offset
        anchor_start = max(0, offset - 64)
        with path.open("rb") as fh:
            fh.seek(anchor_start)
            anchor = fh.read(offset - anchor_start)
        return {
            "path": path.relative_to(root).as_posix(),
            "offset": offset,
            "anchor_sha256": deribit.hashlib.sha256(anchor).hexdigest(),
        }

    def write_restart_checkpoint(
        self,
        root: Path,
        *,
        sequences: dict[str, int],
        ticker_buckets: dict[str, int] | None = None,
        cursors: dict[str, dict] | None = None,
        schema_version: int | None = None,
        ticker_sample_ms: int | None = None,
        required_source_timestamps: dict[str, int] | None = None,
        trade_source_timestamps: dict[str, int] | None = None,
    ) -> None:
        schema_version = (
            deribit.TRADE_SEQUENCE_STATE_SCHEMA_VERSION
            if schema_version is None
            else schema_version
        )
        payload = {
            "schema_version": schema_version,
            "sequences": sequences,
            "ticker_buckets": ticker_buckets or {},
            "cursors": cursors or {},
        }
        if schema_version >= 3:
            payload["required_source_timestamps"] = (
                required_source_timestamps or {}
            )
        if schema_version >= 4:
            payload["ticker_sample_ms"] = (
                deribit.TICKER_SAMPLE_SECONDS * 1000
                if ticker_sample_ms is None
                else ticker_sample_ms
            )
        if schema_version >= 5:
            payload["trade_source_timestamps"] = (
                {
                    instrument_name: 1_800_000_000_000
                    for instrument_name in sequences
                }
                if trade_source_timestamps is None
                else trade_source_timestamps
            )
        if schema_version >= 6:
            payload["ticker_cap_state"] = {
                "day_utc": None,
                "total_exhausted": False,
            }
        (root / deribit.TRADE_SEQUENCE_STATE_FILE).write_text(
            json.dumps(payload),
            encoding="utf-8",
        )

    def trade_record(
        self,
        sequence: int,
        received_ms: int,
        *,
        trade: dict | None = None,
        dedup_id: str | None = None,
    ) -> dict:
        trade = (
            self.trade_payload(sequence, received_ms=received_ms)
            if trade is None
            else trade
        )
        canonical_id = deribit.sha256_id(
            f"trade|{trade['instrument_name']}|{trade['trade_id']}|"
            f"{trade['trade_seq']}|{trade['timestamp']}"
        )
        return {
            "record_type": "option_trade",
            "dedup_id": canonical_id if dedup_id is None else dedup_id,
            "received_ts_ms": received_ms,
            "trade": trade,
        }

    def ticker_record(
        self,
        timestamp: int,
        *,
        ticker: dict | None = None,
        dedup_id: str | None = None,
        sample_interval_ms: int | None = None,
        received_ms: int | None = None,
    ) -> dict:
        ticker = (
            self.ticker_payload(timestamp)
            if ticker is None
            else ticker
        )
        return {
            "record_type": "option_ticker",
            "dedup_id": (
                deribit.sha256_id(
                    f"ticker|{self.instrument_name}|{timestamp}"
                )
                if dedup_id is None
                else dedup_id
            ),
            "received_ts_ms": (
                timestamp if received_ms is None else received_ms
            ),
            "sample_interval_ms": (
                deribit.TICKER_SAMPLE_SECONDS * 1000
                if sample_interval_ms is None
                else sample_interval_ms
            ),
            "instrument_name": self.instrument_name,
            "ticker": ticker,
        }


    def test_ticker_channels_are_required_by_default(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        self.assertTrue(collector.capture_tickers)
        self.assertIn(
            f"ticker.{self.instrument_name}.agg2",
            collector.default_channels(),
        )

    def test_ticker_cap_exhaustion_preserves_critical_capture(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-cap-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            ticker_day_cap=1,
            ticker_total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        collector.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        channel = f"ticker.{self.instrument_name}.agg2"
        received_ms = 1_800_000_000_000
        self.assertTrue(
            collector.ingest_ticker(
                self.instrument_name,
                self.ticker_payload(received_ms),
                received_ms,
            )
        )
        next_ms = received_ms + collector.ticker_sample_ms
        self.assertTrue(
            collector.ingest_ticker(
                self.instrument_name,
                self.ticker_payload(next_ms),
                next_ms,
            )
        )
        trade_received_ms = next_ms + 1_000
        collector.ingest_trade(
            self.trade_payload(1, received_ms=trade_received_ms),
            trade_received_ms,
        )
        index_received_ms = trade_received_ms + 1_000
        self.assertTrue(
            collector.ingest_index(
                "btc_usd",
                {
                    "timestamp": index_received_ms,
                    "price": 90_000.0,
                    "index_name": "btc_usd",
                },
                index_received_ms,
            )
        )

        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        self.assertEqual((day_dir / "ticker.jsonl").stat().st_size, 0)
        self.assertEqual(
            len((day_dir / "trades.jsonl").read_text().splitlines()),
            1,
        )
        self.assertEqual(
            len((day_dir / "index.jsonl").read_text().splitlines()),
            1,
        )
        self.assertEqual(
            collector._required_last_source_timestamp_ms[channel],
            next_ms,
        )
        self.assertEqual(collector.stats["ticker_cap_exhaustions"], 1)
        self.assertEqual(collector.stats["ticker_records_dropped_cap"], 2)
        self.assertFalse(collector.stats.get("cap_caught", False))
        self.assertTrue(collector.capture_tickers)
        self.assertIn(channel, collector.default_channels())
        events = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "events.jsonl").read_text().splitlines()
        ]
        cap_events = [
            event
            for event in events
            if event.get("event") == "ticker_cap_exceeded"
        ]
        self.assertEqual(len(cap_events), 1)
        self.assertEqual(cap_events[0]["detail"]["scope"], "daily")
        self.assertEqual(
            cap_events[0]["detail"]["day_cap_bytes"],
            1,
        )

    def test_ticker_budget_accounting_is_independent(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-budget-"
        ))
        sink = deribit.Sink(
            out,
            day_cap_bytes=10_000_000,
            total_cap_bytes=100_000_000,
            ticker_day_cap_bytes=10_000_000,
            ticker_total_cap_bytes=100_000_000,
        )
        self.addCleanup(sink.close)
        received_ms = 1_800_000_000_000
        ticker_record = {
            "record_type": "option_ticker",
            "dedup_id": "ticker-budget-1",
            "received_ts_ms": received_ms,
        }
        sink.write("ticker", ticker_record, received_ms)
        ticker_bytes = sink.ticker_day_bytes
        critical_bytes = sink.day_bytes
        self.assertGreater(ticker_bytes, 0)
        self.assertGreater(critical_bytes, 0)
        sink.ticker_day_cap = ticker_bytes
        with self.assertRaises(deribit.TickerCapError):
            sink.write(
                "ticker",
                {**ticker_record, "dedup_id": "ticker-budget-2"},
                received_ms + 1,
            )
        sink.write(
            "trades",
            {
                "record_type": "option_trade",
                "dedup_id": "critical-after-ticker-cap",
                "received_ts_ms": received_ms + 1,
            },
            received_ms + 1,
        )
        self.assertGreater(sink.day_bytes, critical_bytes)
        self.assertEqual(sink.ticker_day_bytes, ticker_bytes)

    def test_ticker_total_cap_reconciles_separately_on_restart(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-total-cap-"
        ))
        received_ms = 1_800_000_000_000
        first = deribit.Sink(
            out,
            day_cap_bytes=10_000_000,
            total_cap_bytes=100_000_000,
            ticker_day_cap_bytes=10_000_000,
            ticker_total_cap_bytes=100_000_000,
        )
        try:
            first.write(
                "ticker",
                {
                    "record_type": "option_ticker",
                    "dedup_id": "ticker-before-restart",
                    "received_ts_ms": received_ms,
                },
                received_ms,
            )
            ticker_bytes = first.ticker_day_bytes
            critical_bytes = first.day_bytes
        finally:
            first.close()

        restarted = deribit.Sink(
            out,
            day_cap_bytes=10_000_000,
            total_cap_bytes=100_000_000,
            ticker_day_cap_bytes=10_000_000,
            ticker_total_cap_bytes=ticker_bytes,
        )
        self.addCleanup(restarted.close)
        self.assertEqual(restarted.prior_bytes, critical_bytes)
        self.assertEqual(restarted.ticker_prior_bytes, ticker_bytes)
        restarted.write(
            "events",
            {
                "record_type": "capture_event",
                "dedup_id": "critical-after-restart",
                "event": "critical-after-restart",
                "received_ts_ms": received_ms,
            },
            received_ms,
        )
        with self.assertRaises(deribit.TickerCapError) as cap_error:
            restarted.write(
                "ticker",
                {
                    "record_type": "option_ticker",
                    "dedup_id": "ticker-after-restart",
                    "received_ts_ms": received_ms,
                },
                received_ms,
            )
        self.assertEqual(cap_error.exception.scope, "cumulative")
        restarted.write(
            "index",
            {
                "record_type": "index_price",
                "dedup_id": "critical-index-after-restart",
                "received_ts_ms": received_ms,
            },
            received_ms,
        )


    def test_ticker_cap_pause_survives_restart(self) -> None:
        received_ms = (
            int(time.time() * 1000) // 86_400_000 * 86_400_000
            + 43_200_000
        )
        for scope, restarted_sample_seconds in (
            ("daily", 1),
            ("cumulative", 1),
            ("daily", 2),
            ("cumulative", 2),
        ):
            with self.subTest(
                scope=scope,
                restarted_sample_seconds=restarted_sample_seconds,
            ):
                out = Path(tempfile.mkdtemp(
                    prefix=f"collector-integrity-ticker-{scope}-restart-"
                ))
                first = deribit.Collector(
                    out,
                    currencies=("btc",),
                    day_cap=10_000_000,
                    total_cap=100_000_000,
                    heartbeat=False,
                    ticker_sample_seconds=1,
                    ticker_day_cap=(
                        1 if scope == "daily" else 10_000_000
                    ),
                    ticker_total_cap=(
                        1 if scope == "cumulative" else 100_000_000
                    ),
                )
                first.instruments = {
                    self.instrument_name: {
                        "instrument_name": self.instrument_name,
                        "kind": "option",
                    }
                }
                self.assertTrue(
                    first.ingest_ticker(
                        self.instrument_name,
                        self.ticker_payload(received_ms),
                        received_ms,
                    )
                )
                self.assertEqual(first.stats["ticker_records"], 0)
                first.sink.close()

                restarted = deribit.Collector(
                    out,
                    currencies=("btc",),
                    day_cap=10_000_000,
                    total_cap=100_000_000,
                    heartbeat=False,
                    ticker_sample_seconds=restarted_sample_seconds,
                    ticker_day_cap=10_000_000,
                    ticker_total_cap=100_000_000,
                )
                try:
                    restarted.instruments = {
                        self.instrument_name: {
                            "instrument_name": self.instrument_name,
                            "kind": "option",
                        }
                    }
                    restarted.rehydrate_state()
                    next_ms = received_ms + 1_000
                    self.assertTrue(
                        restarted.ingest_ticker(
                            self.instrument_name,
                            self.ticker_payload(next_ms),
                            next_ms,
                        )
                    )
                    self.assertEqual(restarted.stats["ticker_records"], 0)
                    self.assertEqual(
                        restarted.stats["ticker_records_dropped_cap"],
                        1,
                    )
                    next_day_ms = received_ms + 86_400_000
                    self.assertTrue(
                        restarted.ingest_ticker(
                            self.instrument_name,
                            self.ticker_payload(next_day_ms),
                            next_day_ms,
                        )
                    )
                    self.assertEqual(
                        restarted.stats["ticker_records"],
                        1 if scope == "daily" else 0,
                    )
                    self.assertTrue(
                        restarted.ingest_index(
                            "btc_usd",
                            {"timestamp": next_day_ms, "price": 90_000.0},
                            next_day_ms,
                        )
                    )
                    day_dir = out / f"day{deribit.utc_day(next_day_ms)}"
                    self.assertEqual(
                        len((day_dir / "ticker.jsonl").read_text().splitlines()),
                        1 if scope == "daily" else 0,
                    )
                    index_records = [
                        deribit.strict_json_loads(line)
                        for line in (day_dir / "index.jsonl").read_text().splitlines()
                    ]
                    self.assertEqual(index_records[0]["data"]["price"], 90_000.0)
                finally:
                    restarted.sink.close()
    def test_daily_ticker_cap_recovers_on_next_utc_day(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        self.assertTrue(
            collector.ingest_ticker(
                self.instrument_name,
                self.ticker_payload(received_ms),
                received_ms,
            )
        )
        collector.sink.ticker_day_cap = collector.sink.ticker_day_bytes
        blocked_ms = received_ms + collector.ticker_sample_ms
        self.assertTrue(
            collector.ingest_ticker(
                self.instrument_name,
                self.ticker_payload(blocked_ms),
                blocked_ms,
            )
        )
        blocked_day = deribit.utc_day(blocked_ms)
        self.assertEqual(collector._ticker_cap_day, blocked_day)

        resumed_ms = received_ms + 86_400_000
        self.assertTrue(
            collector.ingest_ticker(
                self.instrument_name,
                self.ticker_payload(resumed_ms),
                resumed_ms,
            )
        )
        self.assertIsNone(collector._ticker_cap_day)
        self.assertEqual(collector.stats["ticker_records"], 2)
        self.assertEqual(collector.stats["ticker_records_dropped_cap"], 1)
        resumed_path = (
            out
            / f"day{deribit.utc_day(resumed_ms)}"
            / "ticker.jsonl"
        )
        resumed = [
            deribit.strict_json_loads(line)
            for line in resumed_path.read_text().splitlines()
        ]
        self.assertEqual(len(resumed), 1)
        self.assertEqual(
            resumed[0]["sample_interval_ms"],
            collector.ticker_sample_ms,
        )
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text().splitlines()
        ]
        self.assertEqual(
            sum(
                event.get("event") == "ticker_cap_recovered"
                for event in events
            ),
            1,
        )

    def test_ticker_partial_tail_cap_never_blocks_critical_sink(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-tail-cap-"
        ))
        received_ms = 1_800_000_000_000
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True)
        partial = b'{"partial":true}'
        ticker_path = day_dir / "ticker.jsonl"
        ticker_path.write_bytes(partial)
        sink = deribit.Sink(
            out,
            day_cap_bytes=10_000_000,
            total_cap_bytes=100_000_000,
            ticker_day_cap_bytes=len(partial),
            ticker_total_cap_bytes=len(partial),
        )
        self.addCleanup(sink.close)
        sink.write(
            "events",
            {
                "record_type": "capture_event",
                "dedup_id": "critical-with-blocked-ticker-tail",
                "event": "critical-probe",
                "received_ts_ms": received_ms,
            },
            received_ms,
        )
        self.assertEqual(ticker_path.read_bytes(), partial)
        with self.assertRaises(deribit.TickerCapError):
            sink.write(
                "ticker",
                {
                    "record_type": "option_ticker",
                    "dedup_id": "blocked-after-partial-tail",
                    "received_ts_ms": received_ms,
                },
                received_ms,
            )
        sink.write(
            "trades",
            {
                "record_type": "option_trade",
                "dedup_id": "critical-still-open",
                "received_ts_ms": received_ms,
            },
            received_ms,
        )
        self.assertIsNone(sink._io_error)
        events = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "events.jsonl").read_text().splitlines()
        ]
        self.assertTrue(events[0]["ticker_tail_recovery_blocked"])

    def test_cli_exposes_independent_ticker_caps(self) -> None:
        args = deribit.parse_args([
            "--max-ticker-day-bytes",
            "123",
            "--max-ticker-total-bytes",
            "456",
        ])
        self.assertEqual(args.max_ticker_day_bytes, 123)
        self.assertEqual(args.max_ticker_total_bytes, 456)

    def test_ticker_sample_interval_below_one_millisecond_is_rejected(
        self,
    ) -> None:
        rejected = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-precision-"
        )) / "capture"
        for seconds in (0.0001, 0.0009):
            with self.subTest(seconds=seconds), self.assertRaises(ValueError):
                deribit.Collector(
                    rejected,
                    currencies=("btc",),
                    heartbeat=False,
                    ticker_sample_seconds=seconds,
                )
        self.assertFalse(rejected.exists())
        for seconds in (0.001, 0.0015):
            with self.subTest(seconds=seconds):
                out = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-ticker-precision-"
                )) / "capture"
                collector = deribit.Collector(
                    out,
                    currencies=("btc",),
                    heartbeat=False,
                    ticker_sample_seconds=seconds,
                )
                self.addCleanup(collector.sink.close)
                self.assertEqual(collector.ticker_sample_ms, 1)

    def test_run_records_capture_configuration_once_across_retry(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_minutes = 1.0
        deribit._STOP_REQUESTED.clear()
        self.addCleanup(deribit._STOP_REQUESTED.clear)
        attempts = 0

        def rehydrate():
            nonlocal attempts
            attempts += 1
            if attempts == 1:
                raise ConnectionError("injected startup retry")
            return {"sequence_checkpoint_loaded": False}

        def seed_then_stop():
            collector.stop = True

        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                collector,
                "rehydrate_state",
                side_effect=rehydrate,
            ),
            mock.patch.object(
                collector,
                "seed_instruments",
                side_effect=seed_then_stop,
            ),
            mock.patch.object(deribit, "stop_aware_sleep"),
        ):
            rc = collector.run()
        self.assertEqual(rc, 0)
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text().splitlines()
        ]
        configurations = [
            event
            for event in events
            if event.get("event") == "collector_configuration"
        ]
        self.assertEqual(len(configurations), 1)
        self.assertEqual(
            configurations[0]["detail"],
            {
                "currencies": ["btc"],
                "capture_tickers": True,
                "ticker_sample_ms": collector.ticker_sample_ms,
                "critical_day_cap_bytes": 10_000_000,
                "critical_total_cap_bytes": 100_000_000,
                "ticker_day_cap_bytes": 268_435_456,
                "ticker_total_cap_bytes": 10_737_418_240,
                "heartbeat": False,
                "max_minutes": 1.0,
            },
        )

    def test_invalid_jsonrpc_values_are_classified_in_every_receive_phase(
        self,
    ) -> None:
        phases = (
            "subscribe_ack",
            "heartbeat_ack",
            "required_coverage",
            "steady_state",
        )
        for phase in phases:
            for value in ([], None):
                with self.subTest(phase=phase, value=value):
                    collector, out = self.make_collector()
                    collector._current_connection_id = (
                        f"{collector.run_id}:c1"
                    )
                    collector._deadline_monotonic = (
                        deribit.time.monotonic() + 1.0
                    )
                    collector.ws = ScriptedDeribitWs([value])
                    if phase == "heartbeat_ack":
                        collector.heartbeat = True
                        collector.active_channels = set(
                            collector.default_channels()
                        )
                    elif phase == "required_coverage":
                        collector.active_channels = {
                            "deribit_price_index.btc_usd"
                        }
                    elif phase == "steady_state":
                        now_ms = int(deribit.time.time() * 1000)
                        collector.last_utc_day = deribit.utc_day(now_ms)
                        collector._next_chain_refresh_attempt = float("inf")
                        collector.last_message_monotonic = (
                            deribit.time.monotonic()
                        )
                    try:
                        with self.assertRaises(
                            deribit.DeribitProtocolFrameError
                        ):
                            if phase == "subscribe_ack":
                                collector._await_ack(
                                    1,
                                    {"deribit_price_index.btc_usd"},
                                    timeout_s=0.1,
                                )
                            elif phase == "heartbeat_ack":
                                with mock.patch.object(
                                    collector,
                                    "_await_ack",
                                    return_value=True,
                                ):
                                    collector.subscribe_all()
                            elif phase == "required_coverage":
                                collector.await_required_coverage(
                                    timeout_s=0.1
                                )
                            else:
                                collector._loop()
                    finally:
                        collector.sink.close()
                    events = [
                        json.loads(line)
                        for path in out.glob("day*/events.jsonl")
                        for line in path.read_text(
                            encoding="utf-8"
                        ).splitlines()
                    ]
                    shape_errors = [
                        row
                        for row in events
                        if row.get("event") == "invalid_jsonrpc_shape"
                    ]
                    self.assertEqual(len(shape_errors), 1)
                    self.assertEqual(
                        shape_errors[0]["detail"]["phase"],
                        phase,
                    )
                    self.assertLess(
                        len(json.dumps(shape_errors[0])),
                        1_000,
                    )
                    self.assertEqual(
                        collector._required_last_received_monotonic,
                        {},
                    )
                    self.assertEqual(
                        collector.stats["trades_records"]
                        + collector.stats["ticker_records"]
                        + collector.stats["index_records"],
                        0,
                    )

    def test_invalid_subscription_params_do_not_advance_message_health(
        self,
    ) -> None:
        for params in ([], None, "not-an-object"):
            with self.subTest(params=params):
                collector, out = self.make_collector()
                collector._current_connection_id = (
                    f"{collector.run_id}:c1"
                )
                before = (
                    collector.msg_seq,
                    collector.last_message_ms,
                    collector.last_message_monotonic,
                )
                try:
                    with self.assertRaises(
                        deribit.DeribitProtocolFrameError
                    ):
                        collector._handle_notification(
                            {
                                "method": "subscription",
                                "params": params,
                            },
                            1_800_000_000_000,
                            phase="steady_state",
                        )
                finally:
                    collector.sink.close()
                self.assertEqual(
                    (
                        collector.msg_seq,
                        collector.last_message_ms,
                        collector.last_message_monotonic,
                    ),
                    before,
                )
                events = [
                    json.loads(line)
                    for path in out.glob("day*/events.jsonl")
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                self.assertEqual(
                    sum(
                        row.get("event") == "invalid_jsonrpc_shape"
                        for row in events
                    ),
                    1,
                )

    def test_invalid_jsonrpc_ack_results_are_classified(self) -> None:
        collector, _ = self.make_collector()
        collector._current_connection_id = f"{collector.run_id}:c1"
        collector._deadline_monotonic = deribit.time.monotonic() + 1.0
        collector.ws = ScriptedDeribitWs([
            {"id": 1, "result": 123},
        ])
        with self.assertRaises(deribit.DeribitProtocolFrameError):
            collector._await_ack(
                1,
                {"deribit_price_index.btc_usd"},
                timeout_s=0.1,
            )
        collector.sink.close()

        heartbeat_collector, _ = self.make_collector()
        heartbeat_collector.heartbeat = True
        heartbeat_collector._current_connection_id = (
            f"{heartbeat_collector.run_id}:c1"
        )
        heartbeat_collector._deadline_monotonic = (
            deribit.time.monotonic() + 1.0
        )
        heartbeat_collector.active_channels = set(
            heartbeat_collector.default_channels()
        )
        heartbeat_collector.ws = ScriptedDeribitWs([
            {"id": 2, "result": []},
        ])
        try:
            with (
                mock.patch.object(
                    heartbeat_collector,
                    "_await_ack",
                    return_value=True,
                ),
                self.assertRaises(deribit.DeribitProtocolFrameError),
            ):
                heartbeat_collector.subscribe_all()
        finally:
            heartbeat_collector.sink.close()

    def test_post_ready_invalid_jsonrpc_shape_reconnects_once(self) -> None:
        collector, out = self.make_collector()
        collector.max_minutes = 1.0
        attempts = 0

        class ProtocolWebSocket:
            def __init__(self, attempt: int):
                self.attempt = attempt
                self.calls = 0

            def recv_message(self, timeout: float):
                del timeout
                self.calls += 1
                if self.attempt == 1:
                    return deribit.OP_TEXT, b"[]"
                if self.calls == 1:
                    notification = {
                        "method": "subscription",
                        "params": {
                            "channel": "trades.option.BTC.100ms",
                            "data": [self_trade],
                        },
                    }
                    return (
                        deribit.OP_TEXT,
                        json.dumps(notification).encode("utf-8"),
                    )
                collector.stop = True
                raise socket.timeout()

            def close(self) -> None:
                return None

        self_trade = self.trade_payload(
            1,
            received_ms=int(deribit.time.time() * 1000),
        )

        def connect():
            nonlocal attempts
            attempts += 1
            return ProtocolWebSocket(attempts)

        stdout = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(collector, "seed_instruments"),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                side_effect=connect,
            ),
            mock.patch.object(
                collector,
                "subscribe_all",
                side_effect=lambda: bool(
                    collector.active_channels.update(
                        {"trades.option.BTC.100ms"}
                    )
                    or True
                ),
            ),
            mock.patch.object(
                collector,
                "await_required_coverage",
                return_value=True,
            ),
            mock.patch.object(deribit, "stop_aware_sleep"),
            mock.patch.object(deribit.sys, "stdout", stdout),
        ):
            rc = collector.run()
        self.assertEqual(rc, 0)
        self.assertEqual(attempts, 2)
        self.assertEqual(collector.stats["reconnects"], 1)
        self.assertEqual(
            stdout.getvalue().count(deribit.READY_BANNER),
            1,
        )
        day_dir = out / f"day{deribit.utc_day(int(deribit.time.time() * 1000))}"
        trades = [
            json.loads(line)
            for line in (day_dir / "trades.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        self.assertEqual(len(trades), 1)
        events = [
            json.loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        self.assertEqual(
            sum(
                row.get("event") == "invalid_jsonrpc_shape"
                for row in events
            ),
            1,
        )
        protocol_disconnects = [
            row
            for row in events
            if row.get("event") == "connection_disconnect"
            and "DeribitProtocolFrameError"
            in row.get("detail", {}).get("reason", "")
        ]
        self.assertEqual(len(protocol_disconnects), 1)

    def test_deribit_rejects_non_rfc_constants_before_routing(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_500
        constants = (
            float("nan"),
            float("inf"),
            float("-inf"),
        )
        for context, constant in (
            (context, constant)
            for context in ("trade", "ticker", "index")
            for constant in constants
        ):
            with self.subTest(context=context, constant=constant):
                if context == "trade":
                    trade = self.trade_payload(101)
                    trade["mark_price"] = constant
                    channel = "trades.option.BTC.100ms"
                    data = [trade]
                elif context == "ticker":
                    ticker = self.ticker_payload(received_ms)
                    ticker["greeks"]["delta"] = constant
                    channel = f"ticker.{self.instrument_name}.agg2"
                    data = ticker
                else:
                    channel = "deribit_price_index.btc_usd"
                    data = {
                        "timestamp": received_ms,
                        "price": 90_000.0,
                        "index_name": "btc_usd",
                        "extra": constant,
                    }
                frame = json.dumps({
                    "method": "subscription",
                    "params": {
                        "channel": channel,
                        "data": data,
                    },
                }).encode("utf-8")
                with self.assertRaises(
                    deribit.DeribitProtocolFrameError
                ):
                    collector._decode_jsonrpc_object(
                        frame,
                        "steady_state",
                        received_ms,
                    )

        self.assertEqual(collector._last_trade_seq, {})
        self.assertEqual(collector.deduper.seen, set())
        self.assertEqual(collector._ticker_last_bucket, {})
        self.assertEqual(
            collector._required_last_source_timestamp_ms,
            {},
        )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        for kind in ("trades", "ticker", "index"):
            path = day_dir / f"{kind}.jsonl"
            self.assertEqual(
                path.read_text(encoding="utf-8").splitlines()
                if path.exists()
                else [],
                [],
            )
        events = [
            json.loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        audits = [
            row
            for row in events
            if row.get("event") == "invalid_jsonrpc_shape"
            and row.get("detail", {}).get("reason")
            == "non_finite_json_constant"
        ]
        self.assertEqual(len(audits), 9)
        self.assertTrue(
            all(
                set(row["detail"])
                == {
                    "phase",
                    "reason",
                    "value_type",
                    "frame_bytes",
                }
                for row in audits
            )
        )

    def test_deribit_sink_rejects_non_finite_record_before_append(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_500
        with self.assertRaises(deribit.JsonSerializationIntegrityError):
            collector.sink.write(
                "trades",
                {
                    "record_type": "option_trade",
                    "trade": {"extra": float("nan")},
                },
                now_ms=received_ms,
            )
        self.assertEqual(collector.sink.day_bytes, 0)
        self.assertEqual(list(out.glob("day*")), [])

    def test_deribit_non_rfc_history_never_rehydrates_state(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-non-rfc-"))
        received_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True)
        trade = self.trade_payload(101)
        trade["mark_price"] = float("nan")
        (day_dir / "trades.jsonl").write_text(
            json.dumps({
                "record_type": "option_trade",
                "dedup_id": "poisoned-trade-id",
                "received_ts_ms": received_ms,
                "trade": trade,
            })
            + "\n",
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        stats = collector.rehydrate_state()
        self.assertGreaterEqual(stats["malformed"], 1)
        self.assertEqual(collector._last_trade_seq, {})
        self.assertNotIn("poisoned-trade-id", collector.deduper.seen)

    def test_deribit_non_rfc_checkpoint_fails_closed_explicitly(
        self,
    ) -> None:
        collector, out = self.make_collector()
        received_ms = 1_800_000_000_500
        collector.ingest_trade(self.trade_payload(101), received_ms)
        self.assertTrue(collector._save_sequence_state())
        collector.sink.close()
        state_path = out / deribit.TRADE_SEQUENCE_STATE_FILE
        payload = json.loads(state_path.read_text(encoding="utf-8"))
        payload["sequences"][self.instrument_name] = float("nan")
        state_path.write_text(
            json.dumps(payload),
            encoding="utf-8",
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        with self.assertRaises(OSError) as failure:
            restarted.rehydrate_state()
        self.assertIn(
            "NonFiniteJsonConstantError",
            str(failure.exception),
        )

    def test_pre_ready_invalid_jsonrpc_shape_never_announces_ready(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_minutes = 1.0

        class ConnectedWebSocket:
            def close(self) -> None:
                return None

        def reject_subscription_shape() -> bool:
            collector.stop = True
            collector._decode_jsonrpc_object(
                b"null",
                "subscribe_ack",
                int(deribit.time.time() * 1000),
            )
            raise AssertionError("unreachable")

        stdout = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(collector, "seed_instruments"),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                return_value=ConnectedWebSocket(),
            ),
            mock.patch.object(
                collector,
                "subscribe_all",
                side_effect=reject_subscription_shape,
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
        ):
            rc = collector.run()
        self.assertEqual(rc, 0)
        self.assertNotIn(deribit.READY_BANNER, stdout.getvalue())
        events = [
            json.loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            sum(
                row.get("event") == "invalid_jsonrpc_shape"
                for row in events
            ),
            1,
        )
        self.assertFalse(
            any(
                row.get("event") == "connection_ready"
                for row in events
            )
        )


    def test_ticker_notification_is_routed_and_downsampled(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        received_ms = 1_800_000_000_100
        collector._route(channel, self.ticker_payload(), received_ms)
        collector._route(channel, self.ticker_payload(timestamp=1_800_000_000_500), received_ms + 400)

        ticker_file = out / f"day{deribit.utc_day(received_ms)}" / "ticker.jsonl"
        rows = [json.loads(line) for line in ticker_file.read_text(encoding="utf-8").splitlines()]
        self.assertEqual(collector.stats["ticker_pushes"], 2)
        self.assertEqual(collector.stats["ticker_records"], 1)
        self.assertEqual(collector.stats["ticker_samples_suppressed"], 1)
        self.assertEqual(rows[0]["record_type"], "option_ticker")
        self.assertEqual(
            rows[0]["sample_interval_ms"],
            collector.ticker_sample_ms,
        )
        self.assertEqual(rows[0]["ticker"]["open_interest"], 12.5)
        self.assertEqual(rows[0]["ticker"]["greeks"]["gamma"], 0.00001)

    def test_incomplete_ticker_never_poisons_same_bucket(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        received_ms = 1_800_000_000_000
        missing = object()
        positive_invalid = (
            missing,
            None,
            False,
            0,
            -1,
            "1",
            float("nan"),
            float("inf"),
            float("-inf"),
        )
        nonnegative_invalid = (
            missing,
            None,
            False,
            -1,
            "1",
            float("nan"),
            float("inf"),
            float("-inf"),
        )
        finite_invalid = (
            missing,
            None,
            False,
            "1",
            float("nan"),
            float("inf"),
            float("-inf"),
        )
        cases = [
            (field, value)
            for field in ("index_price", "underlying_price")
            for value in positive_invalid
        ]
        cases.extend(
            (field, value)
            for field in ("open_interest", "mark_iv")
            for value in nonnegative_invalid
        )
        cases.extend(
            ("greeks", value)
            for value in (missing, None, [], "greeks")
        )
        cases.extend(
            (f"greeks.{field}", value)
            for field in ("delta", "gamma", "theta", "vega", "rho")
            for value in finite_invalid
        )
        cases.append(("greeks.gamma", -1))

        for index, (field, value) in enumerate(cases):
            with self.subTest(field=field, value=repr(value)):
                payload = self.ticker_payload(timestamp=received_ms)
                if field.startswith("greeks."):
                    target = payload["greeks"]
                    key = field.removeprefix("greeks.")
                else:
                    target = payload
                    key = field
                if value is missing:
                    del target[key]
                else:
                    target[key] = value
                collector._route(
                    channel,
                    payload,
                    received_ms,
                    received_monotonic=100.0 + index,
                )

        ticker_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "ticker.jsonl"
        )
        self.assertEqual(ticker_path.stat().st_size, 0)
        self.assertEqual(
            collector._required_last_received_monotonic,
            {},
        )
        self.assertEqual(
            collector._required_last_source_timestamp_ms,
            {},
        )
        self.assertEqual(collector._ticker_last_bucket, {})
        self.assertEqual(collector.deduper.seen, set())
        self.assertFalse(collector._ready_printed)
        self.assertEqual(
            collector.stats["ticker_schema_errors"],
            len(cases),
        )
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        schema_events = [
            event
            for event in events
            if event.get("event") == "ticker_schema_error"
        ]
        self.assertEqual(len(schema_events), len(cases))
        self.assertEqual(
            Counter(
                event["detail"]["field"]
                for event in schema_events
            ),
            Counter(field for field, _ in cases),
        )
        for event in schema_events:
            self.assertEqual(event["detail"]["channel"], channel)
            self.assertEqual(event["detail"]["body_type"], "dict")
            self.assertNotIn("body", event["detail"])
            self.assertNotIn("payload", event["detail"])

        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms),
            received_ms,
            received_monotonic=999.0,
        )
        rows = [
            deribit.strict_json_loads(line)
            for line in ticker_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(len(rows), 1)
        self.assertEqual(
            set(rows[0]["ticker"]["greeks"]),
            {"delta", "gamma", "theta", "vega", "rho"},
        )
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            999.0,
        )
        self.assertEqual(
            collector._required_last_source_timestamp_ms[channel],
            received_ms,
        )
        self.assertEqual(
            collector._ticker_last_bucket[self.instrument_name],
            received_ms // collector.ticker_sample_ms,
        )
        self.assertEqual(len(collector.deduper.seen), 1)

    def test_invalid_numeric_inputs_cannot_satisfy_readiness(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        ticker_channel = f"ticker.{self.instrument_name}.agg2"
        index_channel = "deribit_price_index.btc_usd"
        bad_ticker = self.ticker_payload()
        bad_ticker["open_interest"] = "12.5"
        collector._route(ticker_channel, bad_ticker, received_ms)
        collector._route(
            index_channel,
            {
                "timestamp": received_ms,
                "price": 90_000.0,
                "index_name": "eth_usd",
            },
            received_ms,
        )
        collector._route(
            index_channel,
            {
                "timestamp": received_ms + 1,
                "price": float("nan"),
                "index_name": "btc_usd",
            },
            received_ms + 1,
        )
        self.assertEqual(
            collector._required_last_received_monotonic,
            {},
        )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        self.assertEqual((day_dir / "ticker.jsonl").stat().st_size, 0)
        self.assertEqual((day_dir / "index.jsonl").stat().st_size, 0)
        self.assertEqual(collector.stats["ticker_schema_errors"], 1)
        self.assertEqual(collector.stats["index_schema_errors"], 2)

    def test_future_ticker_cannot_poison_health_or_sampling(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        channel = f"ticker.{self.instrument_name}.agg2"
        future_ms = received_ms + 10 * 365 * 24 * 60 * 60 * 1000
        with self.assertRaises(deribit.SourceTimestampIntegrityError):
            collector._route(
                channel,
                self.ticker_payload(timestamp=future_ms),
                received_ms,
                received_monotonic=100.0,
            )
        self.assertNotIn(channel, collector._required_last_received_monotonic)
        self.assertNotIn(
            channel,
            collector._required_last_source_timestamp_ms,
        )
        self.assertNotIn(
            self.instrument_name,
            collector._ticker_last_bucket,
        )

        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms + 1),
            received_ms + 1,
            received_monotonic=101.0,
        )
        ticker_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "ticker.jsonl"
        )
        self.assertEqual(len(ticker_path.read_text().splitlines()), 1)
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            101.0,
        )

    def test_ticker_liveness_tracks_repeats_not_regressions(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        channel = f"ticker.{self.instrument_name}.agg2"
        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms),
            received_ms,
            received_monotonic=100.0,
        )
        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms),
            received_ms + 1,
            received_monotonic=200.0,
        )
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            200.0,
        )
        collector._check_required_continuity(now_monotonic=240.0)
        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms - 1),
            received_ms + 2,
            received_monotonic=300.0,
        )
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            200.0,
        )
        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms + 3),
            received_ms + 3,
            received_monotonic=400.0,
        )
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            400.0,
        )
        ticker_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "ticker.jsonl"
        )
        self.assertEqual(len(ticker_path.read_text().splitlines()), 1)

    def test_index_liveness_tracks_repeats_not_regressions(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        channel = "deribit_price_index.btc_usd"

        def route(timestamp: int, receipt: int, monotonic: float) -> None:
            collector._route(
                channel,
                {
                    "timestamp": timestamp,
                    "price": 90_000.0,
                    "index_name": "btc_usd",
                },
                receipt,
                received_monotonic=monotonic,
            )

        route(received_ms, received_ms, 100.0)
        route(received_ms, received_ms + 1, 200.0)
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            200.0,
        )
        route(received_ms - 1, received_ms + 2, 300.0)
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            200.0,
        )
        with self.assertRaises(deribit.SourceTimestampIntegrityError):
            route(
                received_ms + 10 * 365 * 24 * 60 * 60 * 1000,
                received_ms + 3,
                400.0,
            )
        route(received_ms + 4, received_ms + 4, 500.0)
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            500.0,
        )
        index_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "index.jsonl"
        )
        self.assertEqual(len(index_path.read_text().splitlines()), 2)

    def test_source_timestamps_checkpoint_and_reject_legacy_poison(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        ticker_channel = f"ticker.{self.instrument_name}.agg2"
        index_channel = "deribit_price_index.btc_usd"
        collector._route(
            ticker_channel,
            self.ticker_payload(timestamp=received_ms),
            received_ms,
            received_monotonic=100.0,
        )
        collector._route(
            index_channel,
            {
                "timestamp": received_ms,
                "price": 90_000.0,
                "index_name": "btc_usd",
            },
            received_ms,
            received_monotonic=100.0,
        )
        collector._sequence_state_dirty = True
        self.assertTrue(collector._save_sequence_state())
        expected = dict(collector._required_last_source_timestamp_ms)
        collector.sink.close()

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        restarted.rehydrate_state()
        self.assertEqual(
            restarted._required_last_source_timestamp_ms,
            expected,
        )

        poison_out = Path(tempfile.mkdtemp(prefix="collector-integrity-timestamp-poison-"))
        day_dir = poison_out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True)
        ticker_path = day_dir / "ticker.jsonl"
        future_ms = received_ms + 10 * 365 * 24 * 60 * 60 * 1000
        ticker_path.write_text(
            json.dumps({
                "record_type": "option_ticker",
                "received_ts_ms": received_ms,
                "instrument_name": self.instrument_name,
                "ticker": {"timestamp": future_ms},
            }) + "\n",
            encoding="utf-8",
        )
        self.write_restart_checkpoint(
            poison_out,
            sequences={},
            ticker_buckets={
                self.instrument_name: future_ms
                // deribit.TICKER_SAMPLE_SECONDS
                // 1000,
            },
            cursors={
                "ticker": self.restart_cursor(
                    poison_out,
                    ticker_path,
                    offset=0,
                ),
            },
        )
        poison = deribit.Collector(
            poison_out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(poison.sink.close)
        stats = poison.rehydrate_state()
        self.assertNotIn(
            self.instrument_name,
            poison._ticker_last_bucket,
        )
        self.assertNotIn(
            ticker_channel,
            poison._required_last_source_timestamp_ms,
        )
        self.assertGreater(stats["source_timestamp_rejected"], 0)

    def test_readiness_requires_ticker_and_index_payloads(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        ticker_channel = f"ticker.{self.instrument_name}.agg2"
        index_channel = "deribit_price_index.btc_usd"
        now_ms = int(deribit.time.time() * 1000)
        collector.active_channels = set(collector.default_channels())
        collector.ws = ScriptedDeribitWs([
            {
                "jsonrpc": "2.0",
                "method": "subscription",
                "params": {
                    "channel": ticker_channel,
                    "data": self.ticker_payload(timestamp=now_ms),
                },
            },
            {
                "jsonrpc": "2.0",
                "method": "subscription",
                "params": {
                    "channel": index_channel,
                    "data": {
                        "timestamp": now_ms,
                        "price": 90_000.0,
                        "index_name": "btc_usd",
                    },
                },
            },
        ])
        self.assertTrue(collector.await_required_coverage(timeout_s=0.1))

        missing, _ = self.make_collector()
        self.addCleanup(missing.sink.close)
        missing.active_channels = set(missing.default_channels())
        missing.ws = ScriptedDeribitWs([
            {
                "jsonrpc": "2.0",
                "method": "subscription",
                "params": {
                    "channel": ticker_channel,
                    "data": self.ticker_payload(timestamp=now_ms),
                },
            },
        ])
        self.assertFalse(missing.await_required_coverage(timeout_s=0.01))

    def test_ack_wait_preserves_interleaved_ticker_payload(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        channels = collector.default_channels()
        now_ms = int(deribit.time.time() * 1000)
        collector.ws = ScriptedDeribitWs([
            {
                "jsonrpc": "2.0",
                "method": "subscription",
                "params": {
                    "channel": channel,
                    "data": self.ticker_payload(timestamp=now_ms),
                },
            },
            {
                "jsonrpc": "2.0",
                "method": "subscription",
                "params": {
                    "channel": "deribit_price_index.btc_usd",
                    "data": {
                        "timestamp": now_ms,
                        "price": 90_000.0,
                        "index_name": "btc_usd",
                    },
                },
            },
            {"jsonrpc": "2.0", "id": 1, "result": channels},
        ])
        self.assertTrue(collector.subscribe_all())
        self.assertTrue(collector.await_required_coverage(timeout_s=0.01))
        ticker_file = (
            out
            / f"day{deribit.utc_day(int(deribit.time.time() * 1000))}"
            / "ticker.jsonl"
        )
        self.assertEqual(len(ticker_file.read_text(encoding="utf-8").splitlines()), 1)

    def test_readiness_rejects_stale_pre_ack_channel_observations(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        expected = {
            f"ticker.{self.instrument_name}.agg2",
            "deribit_price_index.btc_usd",
        }
        collector.active_channels = set(expected)
        stale_monotonic = (
            deribit.time.monotonic()
            - deribit.TICKER_STALE_SECONDS
            - 1
        )
        collector._required_last_received_monotonic = {
            channel: stale_monotonic
            for channel in expected
        }

        class SilentWebSocket:
            def recv_message(self, timeout: float):
                del timeout
                raise socket.timeout()

        collector.ws = SilentWebSocket()
        self.assertFalse(
            collector.await_required_coverage(
                expected,
                timeout_s=0.001,
            )
        )

    def test_dynamic_subscription_mismatch_fails_closed(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        collector.ws = ScriptedDeribitWs([
            {"jsonrpc": "2.0", "id": 1, "result": []},
        ])
        with self.assertRaises(ConnectionError):
            collector._modify([channel], "public/subscribe")
        self.assertNotIn(channel, collector.active_channels)

    def test_session_and_required_channel_watchdogs_fail_closed(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        ticker_channel = f"ticker.{self.instrument_name}.agg2"
        index_channel = "deribit_price_index.btc_usd"
        collector.active_channels = {ticker_channel, index_channel}
        collector.last_message_monotonic = 1_000.0
        with self.assertRaises(ConnectionError):
            collector._check_connection_liveness(
                1_000.0 + deribit.SESSION_IDLE_TIMEOUT_SECONDS + 1.0
            )

        now_monotonic = 10_000.0
        collector._next_required_health_check = 0.0
        collector._required_last_received_monotonic[ticker_channel] = (
            now_monotonic
        )
        collector._required_last_received_monotonic[index_channel] = (
            now_monotonic - deribit.TICKER_STALE_SECONDS - 1.0
        )
        with self.assertRaises(ConnectionError):
            collector._check_required_continuity(now_monotonic)

    def test_trade_gap_advances_after_durable_boundary_and_reconnects_once(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_100
        collector.ingest_trade(self.trade_payload(100), received_ms)
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(
                self.trade_payload(99, trade_id="late-99"),
                received_ms + 1,
            )
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(
                self.trade_payload(102),
                received_ms + 2,
            )
        collector.ingest_trade(
            self.trade_payload(103),
            received_ms + 3,
        )
        trade_file = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        self.assertEqual(
            [
                json.loads(line)["trade"]["trade_seq"]
                for line in trade_file.read_text(
                    encoding="utf-8"
                ).splitlines()
            ],
            [100, 102, 103],
        )
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            103,
        )
        event_file = trade_file.with_name("events.jsonl")
        events = [
            json.loads(line)
            for line in event_file.read_text(encoding="utf-8").splitlines()
        ]
        gaps = [
            row
            for row in events
            if row.get("event") == "trade_continuity_gap"
        ]
        self.assertEqual(len(gaps), 1)
        self.assertEqual(gaps[0]["detail"]["expected"], 101)
        self.assertEqual(gaps[0]["detail"]["observed"], 102)

    def test_out_of_order_trade_batch_does_not_create_false_gap(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        collector.ingest_trade(self.trade_payload(100, received_ms=received_ms), received_ms)
        collector._route(
            "trades.option.BTC.100ms",
            [
                self.trade_payload(102, received_ms=received_ms),
                self.trade_payload(101, received_ms=received_ms),
            ],
            received_ms + 1,
        )
        self.assertEqual(collector.stats["trade_sequence_gaps"], 0)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            102,
        )
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        self.assertEqual(
            len(trade_path.read_text(encoding="utf-8").splitlines()),
            3,
        )

    def test_malformed_trade_does_not_drop_valid_batch_siblings(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        collector.ingest_trade(self.trade_payload(100, received_ms=received_ms), received_ms)
        with self.assertRaises(ConnectionError):
            collector._route(
                "trades.option.BTC.100ms",
                [
                    {"instrument_name": self.instrument_name},
                    self.trade_payload(102, received_ms=received_ms),
                    self.trade_payload(101, received_ms=received_ms),
                ],
                received_ms + 1,
            )
        self.assertEqual(collector.stats["trade_schema_errors"], 1)
        self.assertEqual(collector.stats["trade_sequence_gaps"], 0)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            102,
        )
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        sequences = [
            json.loads(line)["trade"]["trade_seq"]
            for line in trade_path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(sequences, [100, 101, 102])
    def test_invalid_option_trade_fields_never_mutate_state(self) -> None:
        missing = object()
        invalid_by_field = {
            "instrument_name": (missing, "", 1, "BTC-NOT-CURRENT-C"),
            "trade_id": (missing, "", 1),
            "trade_seq": (missing, True, 0, -1, 1.5),
            "timestamp": (missing, True, 0, -1, 1.5),
            "direction": (missing, "", "hold", 1),
            "amount": (
                missing,
                True,
                0,
                -1,
                "1",
                float("nan"),
                float("inf"),
            ),
            "price": (
                missing,
                True,
                0,
                -1,
                "1",
                float("nan"),
                float("inf"),
            ),
            "index_price": (
                missing,
                True,
                0,
                -1,
                "1",
                float("nan"),
                float("inf"),
            ),
            "iv": (
                missing,
                True,
                -1,
                "1",
                float("nan"),
                float("inf"),
            ),
        }
        cases = [
            ("body", [], "body", False),
            (
                "non-option-instrument",
                self.trade_payload(101),
                "instrument_name",
                True,
            ),
        ]
        for field, invalid_values in invalid_by_field.items():
            for index, value in enumerate(invalid_values):
                trade = self.trade_payload(101)
                if value is missing:
                    del trade[field]
                else:
                    trade[field] = value
                cases.append((f"{field}-{index}", trade, field, False))

        received_ms = 1_800_000_000_500
        for name, trade, expected_field, non_option in cases:
            with self.subTest(name=name):
                collector, out = self.make_collector()
                self.addCleanup(collector.sink.close)
                if non_option:
                    collector.instruments[self.instrument_name]["kind"] = (
                        "future"
                    )
                with self.assertRaises(ConnectionError):
                    collector.ingest_trade(trade, received_ms)
                self.assertEqual(collector.stats["trade_schema_errors"], 1)
                self.assertEqual(collector._last_trade_seq, {})
                self.assertEqual(collector.deduper.seen, set())
                self.assertFalse(collector._sequence_state_dirty)
                day_dir = out / f"day{deribit.utc_day(received_ms)}"
                trade_path = day_dir / "trades.jsonl"
                trade_lines = (
                    trade_path.read_text(encoding="utf-8").splitlines()
                    if trade_path.exists()
                    else []
                )
                self.assertEqual(trade_lines, [])
                self.assertFalse(
                    (out / deribit.TRADE_SEQUENCE_STATE_FILE).exists()
                )
                audits = [
                    row
                    for row in (
                        json.loads(line)
                        for line in (
                            day_dir / "events.jsonl"
                        ).read_text(encoding="utf-8").splitlines()
                    )
                    if row.get("event") == "trade_schema_error"
                ]
                self.assertEqual(len(audits), 1)
                self.assertEqual(
                    audits[0]["detail"],
                    {
                        "field": expected_field,
                        "body_type": type(trade).__name__,
                    },
                )

    def test_corrected_trade_after_schema_rejection_is_accepted(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_500
        invalid = self.trade_payload(101)
        invalid["price"] = float("nan")
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(invalid, received_ms)
        self.assertEqual(collector._last_trade_seq, {})
        self.assertEqual(collector.deduper.seen, set())

        collector.ingest_trade(self.trade_payload(101), received_ms + 1)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            101,
        )
        self.assertEqual(len(collector.deduper.seen), 1)
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        rows = [
            json.loads(line)
            for line in trade_path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["trade"]["trade_seq"], 101)


    def test_unknown_new_listing_trade_is_durable_before_reconnect(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector._current_connection_id = f"{collector.run_id}:c1"
        channel = "trades.option.BTC.100ms"
        received_ms = 1_800_000_000_500
        instrument_name = "BTC-30DEC30-95000-C"
        trigger = {
            **self.trade_payload(101, "new-listing-101"),
            "instrument_name": instrument_name,
        }

        with self.assertRaises(ConnectionError) as raised:
            collector._route(channel, [trigger], received_ms)

        self.assertIsInstance(
            raised.exception,
            deribit.UnresolvedInstrumentTradeError,
        )
        self.assertEqual(
            collector._last_trade_seq[instrument_name],
            101,
        )
        self.assertEqual(len(collector.deduper.seen), 1)
        self.assertEqual(collector.stats["trades_unresolved"], 1)
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        first = [
            deribit.strict_json_loads(line)
            for line in trade_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(len(first), 1)
        self.assertEqual(first[0]["record_type"], "option_trade")
        self.assertEqual(
            first[0]["instrument_metadata_status"],
            "unresolved_at_receipt",
        )
        self.assertEqual(first[0]["source_channel"], channel)
        self.assertEqual(first[0]["received_ts_ms"], received_ms)
        self.assertEqual(
            first[0]["connection_id"],
            collector._current_connection_id,
        )
        self.assertEqual(first[0]["trade"], trigger)

        collector.instruments[instrument_name] = {
            "instrument_name": instrument_name,
            "kind": "option",
        }
        following = {
            **self.trade_payload(102, "new-listing-102"),
            "instrument_name": instrument_name,
        }
        collector._route(channel, [following], received_ms + 1)
        rows = [
            deribit.strict_json_loads(line)
            for line in trade_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            [row["trade"]["trade_seq"] for row in rows],
            [101, 102],
        )
        self.assertNotIn("instrument_metadata_status", rows[1])
        self.assertEqual(collector.stats["trade_sequence_gaps"], 0)

    def test_late_trade_after_authoritative_removal_is_retained(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = "trades.option.BTC.100ms"
        received_ms = 1_800_000_000_500
        instrument_name = "BTC-30DEC30-96000-P"
        collector.instruments[instrument_name] = {
            "instrument_name": instrument_name,
            "kind": "option",
        }
        del collector.instruments[instrument_name]
        trade = {
            **self.trade_payload(201, "late-201"),
            "instrument_name": instrument_name,
        }

        with self.assertRaises(ConnectionError):
            collector._route(channel, [trade], received_ms)

        path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        rows = [
            deribit.strict_json_loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["trade"], trade)
        self.assertEqual(
            rows[0]["instrument_metadata_status"],
            "unresolved_at_receipt",
        )

    def test_unknown_trade_validation_rejects_malformed_or_foreign_names(
        self,
    ) -> None:
        invalid_names = (
            "BTC-NOT-AN-OPTION",
            "SOL-30DEC30-90000-C",
            "ETH-30DEC30-90000-C",
            "BTC-30FOO30-90000-C",
            "BTC-30DEC30-zero-C",
        )
        received_ms = 1_800_000_000_500
        for instrument_name in invalid_names:
            with self.subTest(instrument_name=instrument_name):
                collector, out = self.make_collector()
                self.addCleanup(collector.sink.close)
                trade = {
                    **self.trade_payload(101),
                    "instrument_name": instrument_name,
                }
                with self.assertRaises(ConnectionError):
                    collector._route(
                        "trades.option.BTC.100ms",
                        [trade],
                        received_ms,
                    )
                trade_paths = list(out.glob("day*/trades.jsonl"))
                self.assertTrue(
                    not trade_paths
                    or all(path.stat().st_size == 0 for path in trade_paths)
                )
                self.assertEqual(collector._last_trade_seq, {})
                self.assertEqual(collector.deduper.seen, set())

    def test_unknown_trade_batch_and_sink_commit_order(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = "trades.option.BTC.100ms"
        received_ms = 1_800_000_000_500
        unknown_name = "BTC-30DEC30-97000-C"
        known_name = "BTC-31DEC30-98000-C"
        collector.instruments[known_name] = {
            "instrument_name": known_name,
            "kind": "option",
        }
        unknown = {
            **self.trade_payload(1, "unknown-1"),
            "instrument_name": unknown_name,
        }
        known = {
            **self.trade_payload(1, "known-1"),
            "instrument_name": known_name,
        }

        with self.assertRaises(ConnectionError):
            collector._route(
                channel,
                [unknown, known],
                received_ms,
            )
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        rows = [
            deribit.strict_json_loads(line)
            for line in trade_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(
            {row["trade"]["instrument_name"] for row in rows},
            {unknown_name, known_name},
        )
        self.assertEqual(
            set(collector._last_trade_seq),
            {unknown_name, known_name},
        )

        failed, _ = self.make_collector()
        self.addCleanup(failed.sink.close)
        with (
            mock.patch.object(
                failed.sink,
                "write",
                side_effect=OSError("injected unresolved write failure"),
            ),
            self.assertRaises(OSError),
        ):
            failed._route(channel, [unknown], received_ms)
        self.assertEqual(failed._last_trade_seq, {})
        self.assertEqual(failed.deduper.seen, set())
        self.assertEqual(failed.stats["trades_unresolved"], 0)
        self.assertFalse(failed._sequence_state_dirty)

    def test_unresolved_trade_dedup_rehydrates_on_restart(self) -> None:
        collector, out = self.make_collector()
        received_ms = int(deribit.time.time() * 1000)
        channel = "trades.option.BTC.100ms"
        instrument_name = "BTC-30DEC30-99000-C"
        trade = {
            **self.trade_payload(101, "restart-unresolved-101", received_ms=received_ms),
            "instrument_name": instrument_name,
        }
        with self.assertRaises(ConnectionError):
            collector._route(channel, [trade], received_ms)
        trade_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "trades.jsonl"
        )
        before = trade_path.read_bytes()
        collector.sink.close()

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state(recent_trade_ids=10)
        self.assertEqual(stats["trade_ids"], 1)
        self.assertEqual(
            restarted._last_trade_seq[instrument_name],
            101,
        )
        restarted._route(channel, [trade], received_ms + 1)
        self.assertEqual(trade_path.read_bytes(), before)
        self.assertEqual(restarted.stats["dup_trades"], 1)

    def test_batched_trade_gaps_have_distinct_capture_event_ids(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        other_instrument = "BTC-1JAN27-91000-C"
        collector.instruments[other_instrument] = {
            "instrument_name": other_instrument,
            "kind": "option",
        }
        first_a = self.trade_payload(100, received_ms=received_ms)
        first_b = {
            **self.trade_payload(200, received_ms=received_ms),
            "instrument_name": other_instrument,
            "trade_id": "other-200",
        }
        collector.ingest_trade(first_a, received_ms - 1)
        collector.ingest_trade(first_b, received_ms - 1)
        gap_a = self.trade_payload(102, received_ms=received_ms)
        gap_b = {
            **self.trade_payload(202, received_ms=received_ms),
            "instrument_name": other_instrument,
            "trade_id": "other-202",
        }
        with self.assertRaises(deribit.TradeContinuityError):
            collector._route(
                "trades.option.BTC.100ms",
                [gap_a, gap_b],
                received_ms,
            )
        event_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "events.jsonl"
        )
        gaps = [
            row
            for row in (
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("event") == "trade_continuity_gap"
        ]
        self.assertEqual(len(gaps), 2)
        self.assertEqual(
            {row["detail"]["instrument_name"] for row in gaps},
            {self.instrument_name, other_instrument},
        )
        self.assertEqual(len({row["dedup_id"] for row in gaps}), 2)
        self.assertEqual(
            len({row["capture_event_seq"] for row in gaps}),
            2,
        )

    def test_deduper_retains_exact_configured_fifo_capacity(self) -> None:
        deduper = deribit.Deduper(max_ids=3)
        for dedup_id in ("a", "b", "c", "d"):
            deduper.commit(dedup_id)
        self.assertEqual(list(deduper.order), ["b", "c", "d"])
        self.assertEqual(deduper.seen, {"b", "c", "d"})
        self.assertFalse(deduper.is_duplicate("a"))
        self.assertTrue(deduper.is_duplicate("b"))

    def test_restart_rejects_invalid_trade_payloads_and_forged_ids(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-deribit-replay-trade-schema-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True)
        valid = self.trade_record(100, received_ms)
        invalid_trade = self.trade_payload(101, received_ms=received_ms)
        invalid_trade["amount"] = 0
        poisoned = self.trade_record(
            101,
            received_ms,
            trade=invalid_trade,
        )
        forged = self.trade_record(
            102,
            received_ms,
            dedup_id="f" * 64,
        )
        trade_path = day_dir / "trades.jsonl"
        trade_path.write_text(
            "".join(
                json.dumps(row) + "\n"
                for row in (valid, poisoned, forged)
            ),
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)

        stats = collector.rehydrate_state(recent_trade_ids=10)

        self.assertEqual(
            collector._last_trade_seq,
            {self.instrument_name: 100},
        )
        self.assertEqual(
            collector.deduper.seen,
            {valid["dedup_id"]},
        )
        self.assertEqual(stats["trade_sequence_replay_rejected"], 2)
        self.assertEqual(stats["trade_tail_replay_rejected"], 2)
        collector.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        collector.ingest_trade(
            self.trade_payload(101, received_ms=received_ms),
            received_ms + 1,
        )
        self.assertEqual(collector.stats["trade_sequence_gaps"], 0)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            101,
        )
        self.assertIn(poisoned["dedup_id"], collector.deduper.seen)

    def test_restart_ticker_replay_requires_current_semantic_record(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-deribit-replay-ticker-schema-"
        ))
        timestamp = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(timestamp)}"
        day_dir.mkdir(parents=True)
        incomplete = self.ticker_payload(timestamp)
        del incomplete["greeks"]["gamma"]
        negative_gamma = self.ticker_payload(timestamp)
        negative_gamma["greeks"]["gamma"] = -0.1
        invalid_rows = (
            self.ticker_record(timestamp, ticker=incomplete),
            self.ticker_record(timestamp, ticker=negative_gamma),
            self.ticker_record(timestamp, dedup_id="f" * 64),
            self.ticker_record(
                timestamp,
                sample_interval_ms=(
                    deribit.TICKER_SAMPLE_SECONDS * 1000 + 1
                ),
            ),
        )
        ticker_path = day_dir / "ticker.jsonl"
        ticker_path.write_text(
            "".join(
                json.dumps(row) + "\n" for row in invalid_rows
            ),
            encoding="utf-8",
        )
        rejected = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        stats = rejected.rehydrate_state()
        self.assertEqual(rejected._ticker_last_bucket, {})
        self.assertEqual(
            rejected._required_last_source_timestamp_ms,
            {},
        )
        self.assertEqual(stats["ticker_replay_rejected"], 4)
        rejected.sink.close()

        valid = self.ticker_record(timestamp)
        with ticker_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(valid) + "\n")
        accepted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(accepted.sink.close)
        accepted.rehydrate_state()
        self.assertEqual(
            accepted._ticker_last_bucket[self.instrument_name],
            timestamp // accepted.ticker_sample_ms,
        )
        self.assertEqual(
            accepted._required_last_source_timestamp_ms[
                f"ticker.{self.instrument_name}.agg2"
            ],
            timestamp,
        )

    def test_v3_checkpoint_rebuilds_from_canonical_raw_history(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-deribit-v3-migration-"
        ))
        timestamp = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(timestamp)}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        ticker_path = day_dir / "ticker.jsonl"
        trade_path.write_text(
            json.dumps(self.trade_record(100, timestamp)) + "\n",
            encoding="utf-8",
        )
        ticker_path.write_text(
            json.dumps(self.ticker_record(timestamp)) + "\n",
            encoding="utf-8",
        )
        self.write_restart_checkpoint(
            out,
            schema_version=3,
            sequences={self.instrument_name: 999},
            ticker_buckets={self.instrument_name: 1},
            required_source_timestamps={
                f"ticker.{self.instrument_name}.agg2": timestamp - 1,
            },
            cursors={
                "trades": self.restart_cursor(out, trade_path),
                "ticker": self.restart_cursor(out, ticker_path),
            },
        )
        migrated = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        stats = migrated.rehydrate_state(recent_trade_ids=1)
        self.assertFalse(stats["sequence_checkpoint_loaded"])
        self.assertTrue(stats["sequence_checkpoint_rebuilt"])
        self.assertEqual(
            stats["sequence_checkpoint_rebuild_reason"],
            "legacy_schema_v3",
        )
        self.assertTrue(stats["sequence_checkpoint_written"])
        self.assertEqual(
            migrated._last_trade_seq[self.instrument_name],
            100,
        )
        self.assertEqual(
            migrated._ticker_last_bucket[self.instrument_name],
            timestamp // migrated.ticker_sample_ms,
        )
        migrated.sink.close()
        checkpoint = deribit.strict_json_loads(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            checkpoint["schema_version"],
            deribit.TRADE_SEQUENCE_STATE_SCHEMA_VERSION,
        )
        self.assertEqual(
            checkpoint["ticker_sample_ms"],
            deribit.TICKER_SAMPLE_SECONDS * 1000,
        )

        trusted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(trusted.sink.close)
        trusted_stats = trusted.rehydrate_state(recent_trade_ids=1)
        self.assertTrue(trusted_stats["sequence_checkpoint_loaded"])
        self.assertFalse(trusted_stats["sequence_checkpoint_rebuilt"])
        self.assertEqual(trusted_stats["cursor_replay_bytes"], 0)
        self.assertEqual(
            trusted._last_trade_seq[self.instrument_name],
            100,
        )

    def test_checkpoint_sampling_interval_mismatch_rebuilds_raw_state(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-deribit-interval-migration-"
        ))
        timestamp = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(timestamp)}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        ticker_path = day_dir / "ticker.jsonl"
        trade_path.write_text(
            json.dumps(self.trade_record(100, timestamp)) + "\n",
            encoding="utf-8",
        )
        ticker_path.write_text(
            json.dumps(self.ticker_record(timestamp)) + "\n",
            encoding="utf-8",
        )
        self.write_restart_checkpoint(
            out,
            schema_version=deribit.TRADE_SEQUENCE_STATE_SCHEMA_VERSION,
            ticker_sample_ms=deribit.TICKER_SAMPLE_SECONDS * 1000 + 1,
            sequences={self.instrument_name: 999},
            ticker_buckets={self.instrument_name: 1},
            cursors={
                "trades": self.restart_cursor(out, trade_path),
                "ticker": self.restart_cursor(out, ticker_path),
            },
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)

        stats = collector.rehydrate_state()

        self.assertFalse(stats["sequence_checkpoint_loaded"])
        self.assertTrue(stats["sequence_checkpoint_rebuilt"])
        self.assertEqual(
            stats["sequence_checkpoint_rebuild_reason"],
            "ticker_sample_ms_mismatch",
        )
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            100,
        )
        self.assertEqual(
            collector._ticker_last_bucket[self.instrument_name],
            timestamp // collector.ticker_sample_ms,
        )

    def test_sampling_migration_rejects_corrupt_ticker_cap_state(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        self.write_restart_checkpoint(
            out,
            sequences={},
            ticker_sample_ms=collector.ticker_sample_ms + 1,
        )
        state_path = out / deribit.TRADE_SEQUENCE_STATE_FILE
        payload = json.loads(state_path.read_text(encoding="utf-8"))
        payload["ticker_cap_state"]["total_exhausted"] = "false"
        original = json.dumps(payload)
        state_path.write_text(original, encoding="utf-8")

        with self.assertRaises(OSError):
            collector.rehydrate_state()
        self.assertEqual(state_path.read_text(encoding="utf-8"), original)

    def test_active_trade_rehydrate_restores_only_requested_dedup_tail(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-tail-"))
        received_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True)
        rows = []
        for sequence in range(250):
            rows.append(json.dumps(
                self.trade_record(sequence, received_ms)
            ))
        (day_dir / "trades.jsonl").write_text(
            "\n".join(rows) + "\n",
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        stats = collector.rehydrate_state(recent_trade_ids=20)
        self.assertEqual(stats["trade_ids"], 20)
        self.assertEqual(
            collector.deduper.seen,
            {
                self.trade_record(sequence, received_ms)["dedup_id"]
                for sequence in range(230, 250)
            },
        )
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            249,
        )

    def test_trade_sequence_rehydrates_across_restart(self) -> None:
        collector, out = self.make_collector()
        received_ms = int(deribit.time.time() * 1000)
        collector.ingest_trade(self.trade_payload(200, received_ms=received_ms), received_ms)
        collector.sink.close()
        self.assertFalse(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).exists()
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        restarted.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        restarted.rehydrate_state()
        with self.assertRaises(deribit.TradeContinuityError):
            restarted.ingest_trade(
                self.trade_payload(202, received_ms=received_ms),
                received_ms + 1,
            )

    def test_trade_sequence_rehydrates_across_utc_rollover(self) -> None:
        collector, out = self.make_collector()
        first_ms = 1_800_000_000_000
        second_ms = first_ms + 86_400_000
        collector.ingest_trade(self.trade_payload(100, received_ms=first_ms), first_ms)
        collector.sink.close()

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        restarted.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        with mock.patch.object(
            deribit.time,
            "time",
            return_value=second_ms / 1000,
        ):
            restarted.rehydrate_state()
        with self.assertRaises(deribit.TradeContinuityError):
            restarted.ingest_trade(
                self.trade_payload(102, received_ms=second_ms),
                second_ms,
            )
        self.assertEqual(
            restarted._last_trade_seq[self.instrument_name],
            102,
        )

    def test_sequence_checkpoint_bounds_historical_startup_reads(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-seq-state-"))
        trade_paths = []
        for day_index in range(1, 13):
            day_dir = out / f"day202601{day_index:02d}"
            day_dir.mkdir(parents=True)
            row = self.trade_record(
                day_index,
                received_ms=day_index,
            )
            trade_path = day_dir / "trades.jsonl"
            trade_path.write_text(json.dumps(row) + "\n", encoding="utf-8")
            trade_paths.append(trade_path)
        latest_path = trade_paths[-1]
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 999},
            cursors={"trades": self.restart_cursor(out, latest_path)},
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        active_seconds = deribit.datetime(
            2026,
            2,
            1,
            tzinfo=deribit.timezone.utc,
        ).timestamp()
        original_open = deribit.open_regular_read
        trade_opens = []

        def tracking_open(path: Path, *args, **kwargs):
            if Path(path).name == "trades.jsonl":
                trade_opens.append(Path(path))
            return original_open(path, *args, **kwargs)

        with (
            mock.patch.object(
                deribit.time,
                "time",
                return_value=active_seconds,
            ),
            mock.patch.object(deribit, "open_regular_read", tracking_open),
        ):
            stats = collector.rehydrate_state()
            cached = collector.rehydrate_state()
        self.assertEqual(set(trade_opens), {latest_path})
        self.assertTrue(stats["sequence_checkpoint_loaded"])
        self.assertEqual(cached, stats)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            999,
        )

    def test_restart_cursor_replays_only_suffix_and_bounded_dedup_tail(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-cursor-"))
        day = "20260201"
        day_dir = out / f"day{day}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        ticker_path = day_dir / "ticker.jsonl"
        trade_prefix = json.dumps({"ignored": "x" * 131_072}) + "\n"
        ticker_prefix = json.dumps({"ignored": "y" * 131_072}) + "\n"
        trade_path.write_text(trade_prefix, encoding="utf-8")
        ticker_path.write_text(ticker_prefix, encoding="utf-8")
        trade_cursor = self.restart_cursor(out, trade_path)
        ticker_cursor = self.restart_cursor(out, ticker_path)
        trade_row = self.trade_record(
            101,
            received_ms=1_800_000_000_101,
        )
        ticker_timestamp = 1_800_000_000_000
        ticker_row = self.ticker_record(ticker_timestamp)
        with trade_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(trade_row) + "\n")
        with ticker_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(ticker_row) + "\n")
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 100},
            ticker_buckets={self.instrument_name: 1},
            cursors={
                "trades": trade_cursor,
                "ticker": ticker_cursor,
            },
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        active_seconds = deribit.datetime(
            2026,
            2,
            1,
            tzinfo=deribit.timezone.utc,
        ).timestamp()
        with mock.patch.object(
            deribit.time,
            "time",
            return_value=active_seconds,
        ):
            stats = collector.rehydrate_state(recent_trade_ids=1)
        self.assertEqual(collector._last_trade_seq[self.instrument_name], 101)
        self.assertEqual(
            collector.deduper.seen,
            {trade_row["dedup_id"]},
        )
        self.assertEqual(
            collector._ticker_last_bucket[self.instrument_name],
            ticker_timestamp // collector.ticker_sample_ms,
        )
        self.assertLess(stats["cursor_replay_bytes"], 10_000)
        self.assertLess(stats["tail_replay_bytes"], 100_000)

    def test_cursor_partial_tail_scan_is_anchored_at_eof(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-tail-cursor-"))
        day_dir = out / "day20260201"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        complete_prefix = b"x" * 4_096 + b"\n"
        trade_path.write_bytes(complete_prefix + b"partial")
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        with mock.patch.object(
            deribit,
            "RESTART_TAIL_MAX_BYTES",
            64,
        ):
            cursor = collector._cursor_for_path(trade_path)
        self.assertEqual(cursor["offset"], len(complete_prefix))
        self.assertGreater(cursor["offset"], 4_000)

    def test_cursor_suffix_over_byte_cap_is_terminal(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-suffix-cap-"))
        day = "20260201"
        day_dir = out / f"day{day}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        prefix_row = self.trade_record(
            100,
            received_ms=1_800_000_000_100,
        )
        trade_path.write_text(json.dumps(prefix_row) + "\n", encoding="utf-8")
        cursor = self.restart_cursor(out, trade_path)
        suffix_row = self.trade_record(
            101,
            received_ms=1_800_000_000_101,
        )
        with trade_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(suffix_row) + "\n")
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 100},
            cursors={"trades": cursor},
        )
        original_state = (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_bytes()
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        active_seconds = deribit.datetime(
            2026,
            2,
            1,
            tzinfo=deribit.timezone.utc,
        ).timestamp()
        with (
            mock.patch.object(
                deribit.time,
                "time",
                return_value=active_seconds,
            ),
            mock.patch.object(
                deribit,
                "RESTART_CURSOR_MAX_REPLAY_BYTES",
                32,
                create=True,
            ),
            self.assertRaisesRegex(OSError, "cursor suffix replay"),
        ):
            collector.rehydrate_state()
        self.assertIsNotNone(collector._sequence_state_error)
        self.assertTrue(collector._rehydration_in_progress)
        self.assertEqual(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_bytes(),
            original_state,
        )

    def test_fallback_rehydrate_honors_expired_deadline(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-deadline-"))
        now_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(now_ms)}"
        day_dir.mkdir(parents=True)
        row = self.trade_record(100, received_ms=now_ms)
        (day_dir / "trades.jsonl").write_text(
            json.dumps(row) + "\n",
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        collector._deadline_monotonic = deribit.time.monotonic() - 1
        with self.assertRaisesRegex(
            ConnectionError,
            "deadline.*restart replay",
        ):
            collector.rehydrate_state()
        self.assertEqual(collector._last_trade_seq, {})

    def test_restart_cursor_recovers_cross_day_suffix(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-rollover-cursor-"))
        prior_dir = out / "day20260131"
        current_dir = out / "day20260201"
        prior_dir.mkdir(parents=True)
        current_dir.mkdir(parents=True)
        prior_path = prior_dir / "trades.jsonl"
        current_path = current_dir / "trades.jsonl"
        prior_row = self.trade_record(
            100,
            received_ms=1_800_000_000_100,
        )
        suffix_row = self.trade_record(
            101,
            received_ms=1_800_000_000_101,
        )
        prior_path.write_text(json.dumps(prior_row) + "\n", encoding="utf-8")
        current_path.write_text(json.dumps(suffix_row) + "\n", encoding="utf-8")
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 100},
            cursors={"trades": self.restart_cursor(out, prior_path)},
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        active_seconds = deribit.datetime(
            2026,
            2,
            1,
            tzinfo=deribit.timezone.utc,
        ).timestamp()
        with mock.patch.object(
            deribit.time,
            "time",
            return_value=active_seconds,
        ):
            stats = collector.rehydrate_state(recent_trade_ids=1)
        self.assertEqual(collector._last_trade_seq[self.instrument_name], 101)
        self.assertEqual(
            collector.deduper.seen,
            {suffix_row["dedup_id"]},
        )
        self.assertGreater(stats["cursor_replay_bytes"], 0)
        saved = json.loads(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_text(encoding="utf-8")
        )
        self.assertEqual(
            saved["cursors"]["trades"]["path"],
            "day20260201/trades.jsonl",
        )

    def test_invalid_restart_cursor_never_enters_websocket_loop(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-cursor-error-"))
        now_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(now_ms)}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        trade_path.write_text(
            json.dumps(
                self.trade_record(100, received_ms=now_ms)
            ) + "\n",
            encoding="utf-8",
        )
        cursor = self.restart_cursor(out, trade_path)
        cursor["anchor_sha256"] = "0" * 64
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 100},
            cursors={"trades": cursor},
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        connect = mock.Mock(
            side_effect=AssertionError(
                "websocket entered after invalid restart cursor"
            )
        )
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                connect,
            ),
            mock.patch.object(deribit.sys, "stdout", io.StringIO()),
            mock.patch.object(deribit.sys, "stderr", io.StringIO()),
        ):
            rc = collector.run()
        self.assertEqual(rc, 5)
        connect.assert_not_called()

    def test_durable_write_cadence_survives_abrupt_restart(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        ticker = self.ticker_payload(received_ms)
        with mock.patch.object(
            deribit,
            "RESTART_CHECKPOINT_BYTE_INTERVAL",
            1,
            create=True,
        ):
            collector.ingest_trade(
                self.trade_payload(100, received_ms=received_ms),
                received_ms,
            )
            collector.ingest_ticker(
                self.instrument_name,
                ticker,
                received_ms,
            )
        state_path = out / deribit.TRADE_SEQUENCE_STATE_FILE
        checkpoint = json.loads(state_path.read_text(encoding="utf-8"))
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        self.assertEqual(
            checkpoint["sequences"][self.instrument_name],
            100,
        )
        self.assertEqual(
            checkpoint["ticker_buckets"][self.instrument_name],
            received_ms // collector.ticker_sample_ms,
        )
        self.assertEqual(
            checkpoint["cursors"]["trades"]["offset"],
            (day_dir / "trades.jsonl").stat().st_size,
        )
        self.assertEqual(
            checkpoint["cursors"]["ticker"]["offset"],
            (day_dir / "ticker.jsonl").stat().st_size,
        )

        collector.sink.close()
        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state(recent_trade_ids=1)
        self.assertEqual(
            restarted._last_trade_seq[self.instrument_name],
            100,
        )
        self.assertEqual(
            restarted._ticker_last_bucket[self.instrument_name],
            received_ms // restarted.ticker_sample_ms,
        )
        self.assertEqual(stats["cursor_replay_bytes"], 0)
        self.assertEqual(len(restarted.deduper.seen), 1)

    def test_sequence_checkpoint_write_failure_is_terminal(self) -> None:
        collector, _ = self.make_collector()
        collector._merge_trade_sequence(self.instrument_name, 100)
        collector._sequence_state_bytes_since_checkpoint = 123
        collector._sequence_state_checkpoint_due_monotonic = 17.0
        with mock.patch.object(
            deribit.os,
            "replace",
            side_effect=PermissionError("read-only checkpoint"),
        ):
            with self.assertRaises(OSError):
                collector._save_sequence_state()
        self.assertIsNotNone(collector._sequence_state_error)
        self.assertTrue(collector._sequence_state_dirty)
        self.assertEqual(
            collector._sequence_state_bytes_since_checkpoint,
            123,
        )
        self.assertEqual(
            collector._sequence_state_checkpoint_due_monotonic,
            17.0,
        )

        collector.max_minutes = 0
        stderr = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(deribit.sys, "stdout", io.StringIO()),
            mock.patch.object(deribit.sys, "stderr", stderr),
        ):
            self.assertEqual(collector.run(), 5)
        self.assertIn(
            "DERIBIT_CAPTURE_SEQUENCE_STATE_ERROR",
            stderr.getvalue(),
        )

    def test_startup_sequence_error_never_enters_websocket_loop(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-startup-"))
        (out / deribit.TRADE_SEQUENCE_STATE_FILE).write_text(
            "{malformed",
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        connect = mock.Mock(
            side_effect=AssertionError(
                "websocket entered after terminal startup error"
            )
        )
        stderr = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                connect,
            ),
            mock.patch.object(deribit.sys, "stdout", io.StringIO()),
            mock.patch.object(deribit.sys, "stderr", stderr),
        ):
            rc = collector.run()
        self.assertEqual(rc, 5)
        connect.assert_not_called()
        self.assertIn(
            "DERIBIT_CAPTURE_SEQUENCE_STATE_ERROR",
            stderr.getvalue(),
        )

    def test_deadline_bounds_ack_and_receive_waits(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.ws = ScriptedDeribitWs([])
        collector._deadline_monotonic = deribit.time.monotonic() + 0.01
        started = deribit.time.monotonic()
        self.assertFalse(
            collector._await_ack(
                1,
                {"deribit_price_index.btc_usd"},
                timeout_s=0.1,
            )
        )
        self.assertLess(deribit.time.monotonic() - started, 0.05)

        observed = {}

        class DeadlineWebSocket:
            def recv_message(self, timeout: float):
                observed["timeout"] = timeout
                collector.stop = True
                raise socket.timeout()

        collector.stop = False
        collector.ws = DeadlineWebSocket()
        collector._deadline_monotonic = deribit.time.monotonic() + 0.02
        collector.last_utc_day = deribit.utc_day(
            int(deribit.time.time() * 1000)
        )
        collector.last_chain_reconcile = deribit.time.monotonic()
        collector._loop()
        self.assertLessEqual(observed["timeout"], 0.02)

    def test_signal_interrupts_ack_and_persists_exact_stop_reason(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        observed = {"closed": 0}

        class SignalWebSocket:
            def recv_message(self, timeout: float):
                observed["timeout"] = timeout
                collector._handle_signal(signal.SIGTERM, None)
                raise socket.timeout()

            def close(self) -> None:
                observed["closed"] += 1

        collector.ws = SignalWebSocket()
        with mock.patch.object(deribit, "request_stop"):
            self.assertFalse(
                collector._await_ack(
                    1,
                    {"deribit_price_index.btc_usd"},
                )
            )
        self.assertLessEqual(
            observed["timeout"],
            deribit.INTERRUPT_POLL_SECONDS,
        )
        self.assertEqual(observed["closed"], 1)
        collector._emit_session_stop()
        collector._emit_session_stop()
        collector.sink.close()
        event_path = (
            out
            / f"day{deribit.utc_day(int(deribit.time.time() * 1000))}"
            / "events.jsonl"
        )
        stops = [
            row
            for row in (
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("event") == "session_stop"
        ]
        self.assertEqual(len(stops), 1)
        self.assertEqual(stops[0]["detail"]["reason"], "sigterm")

    def test_interrupt_polling_does_not_flood_idle_controls(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        now_ms = int(deribit.time.time() * 1000)
        now_monotonic = deribit.time.monotonic()
        collector.last_message_ms = now_ms
        collector.last_message_monotonic = now_monotonic
        collector.last_utc_day = deribit.utc_day(now_ms)
        collector.last_chain_reconcile = now_monotonic
        collector._next_chain_refresh_attempt = float("inf")
        collector._last_authoritative_chain_refresh_monotonic = now_monotonic

        class ImmediateTimeoutWebSocket:
            def __init__(self) -> None:
                self.calls = 0

            def recv_message(self, timeout: float):
                del timeout
                self.calls += 1
                if self.calls == 9:
                    collector.stop = True
                raise socket.timeout()

        collector.ws = ImmediateTimeoutWebSocket()
        collector._loop()
        event_path = (
            out
            / f"day{deribit.utc_day(now_ms)}"
            / "events.jsonl"
        )
        events = (
            [
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            ]
            if event_path.is_file()
            else []
        )
        self.assertEqual(
            len([
                row
                for row in events
                if row.get("event") == "connection_idle_ping"
            ]),
            0,
        )

    def test_websocket_handshake_failure_closes_tls_socket(self) -> None:
        class TlsSocket:
            def __init__(self):
                self.close_count = 0

            def settimeout(self, _timeout):
                return None

            def sendall(self, _request):
                return None

            def recv(self, _size):
                raise OSError("handshake receive failed")

            def close(self):
                self.close_count += 1

        tls_socket = TlsSocket()
        context = mock.Mock()
        context.wrap_socket.return_value = tls_socket
        with (
            mock.patch.object(
                deribit.socket,
                "create_connection",
                return_value=mock.Mock(),
            ),
            mock.patch.object(
                deribit.ssl,
                "create_default_context",
                return_value=context,
            ),
            self.assertRaisesRegex(OSError, "handshake receive failed"),
        ):
            deribit.WsClient(host="example.invalid")
        self.assertEqual(tls_socket.close_count, 1)

    def test_websocket_handshake_rejects_wrong_accept(self) -> None:
        class TlsSocket:
            def __init__(self):
                self.close_count = 0

            def settimeout(self, _timeout):
                return None

            def sendall(self, _request):
                return None

            def recv(self, _size):
                return (
                    b"HTTP/1.1 101 Switching Protocols\r\n"
                    b"Upgrade: websocket\r\n"
                    b"Connection: Upgrade\r\n"
                    b"Sec-WebSocket-Accept: wrong\r\n\r\n"
                )

            def close(self):
                self.close_count += 1

        tls_socket = TlsSocket()
        context = mock.Mock()
        context.wrap_socket.return_value = tls_socket
        with (
            mock.patch.object(
                deribit.socket,
                "create_connection",
                return_value=mock.Mock(),
            ),
            mock.patch.object(
                deribit.ssl,
                "create_default_context",
                return_value=context,
            ),
            self.assertRaisesRegex(ConnectionError, "accept"),
        ):
            deribit.WsClient(host="example.invalid")
        self.assertEqual(tls_socket.close_count, 1)

    def test_websocket_handshake_header_size_is_bounded(self) -> None:
        class TlsSocket:
            def __init__(self):
                self.recv_count = 0
                self.close_count = 0

            def settimeout(self, _timeout):
                return None

            def sendall(self, _request):
                return None

            def recv(self, _size):
                self.recv_count += 1
                return b"x" * 40

            def close(self):
                self.close_count += 1

        tls_socket = TlsSocket()
        context = mock.Mock()
        context.wrap_socket.return_value = tls_socket
        with (
            mock.patch.object(
                deribit.socket,
                "create_connection",
                return_value=mock.Mock(),
            ),
            mock.patch.object(
                deribit.ssl,
                "create_default_context",
                return_value=context,
            ),
            mock.patch.object(
                deribit,
                "WS_HANDSHAKE_MAX_HEADER_BYTES",
                64,
            ),
            self.assertRaisesRegex(ConnectionError, "headers exceed"),
        ):
            deribit.WsClient(host="example.invalid")
        self.assertEqual(tls_socket.recv_count, 2)
        self.assertEqual(tls_socket.close_count, 1)

    def test_websocket_valid_handshake_preserves_buffered_frame(self) -> None:
        nonce = b"k" * 16
        key = deribit.base64.b64encode(nonce).decode()
        accept = deribit.base64.b64encode(
            deribit.hashlib.sha1(
                (
                    key
                    + "258EAFA5-E914-47DA-95CA-C5AB0DC85B11"
                ).encode("ascii")
            ).digest()
        )
        first_frame = b"\x81\x02{}"

        class TlsSocket:
            def settimeout(self, _timeout):
                return None

            def sendall(self, _request):
                return None

            def recv(self, _size):
                return (
                    b"HTTP/1.1 101 Switching Protocols\r\n"
                    b"uPgRaDe: WebSocket\r\n"
                    b"Connection: keep-alive, Upgrade\r\n"
                    b"Sec-WebSocket-Accept: "
                    + accept
                    + b"\r\n\r\n"
                    + first_frame
                )

            def close(self):
                raise AssertionError("valid handshake socket was closed")

        tls_socket = TlsSocket()
        context = mock.Mock()
        context.wrap_socket.return_value = tls_socket
        with (
            mock.patch.object(
                deribit.socket,
                "create_connection",
                return_value=mock.Mock(),
            ),
            mock.patch.object(
                deribit.ssl,
                "create_default_context",
                return_value=context,
            ),
            mock.patch.object(
                deribit.os,
                "urandom",
                return_value=nonce,
            ),
        ):
            client = deribit.WsClient(host="example.invalid")
        self.assertIs(client.sock, tls_socket)
        self.assertEqual(client._buf, first_frame)

    def test_websocket_partial_frame_survives_timeout(self) -> None:
        class ScriptedSocket:
            def __init__(self, actions):
                self.actions = list(actions)

            def settimeout(self, _timeout):
                return None

            def recv(self, _size):
                action = self.actions.pop(0)
                if isinstance(action, BaseException):
                    raise action
                return action

        client = object.__new__(deribit.WsClient)
        client.sock = ScriptedSocket(
            [b"\x81\x05he", socket.timeout("partial"), b"llo"]
        )
        client._buf = b""
        client._frag = None
        with self.assertRaises(socket.timeout):
            client.recv_message(0.01)
        self.assertEqual(
            client.recv_message(0.01),
            (deribit.OP_TEXT, b"hello"),
        )

    def test_websocket_extended_header_eof_raises_promptly(self) -> None:
        class EofSocket:
            def __init__(self):
                self.calls = 0

            def settimeout(self, _timeout):
                return None

            def recv(self, _size):
                self.calls += 1
                if self.calls == 1:
                    return b"\x81\x7e"
                if self.calls == 2:
                    return b""
                raise AssertionError("EOF read loop did not stop")

        client = object.__new__(deribit.WsClient)
        client.sock = EofSocket()
        client._buf = b""
        client._frag = None
        with self.assertRaises(ConnectionError):
            client.recv_message(0.01)
        self.assertEqual(client.sock.calls, 2)

    def test_websocket_rejects_oversized_declared_payload_before_body(
        self,
    ) -> None:
        class OversizedHeaderSocket:
            def __init__(self) -> None:
                self.calls = 0

            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                self.calls += 1
                if self.calls == 1:
                    return (
                        b"\x81\x7f"
                        + deribit.struct.pack(
                            ">Q",
                            deribit.WS_MAX_MESSAGE_BYTES + 1,
                        )
                    )
                raise AssertionError(
                    "oversized frame body was requested"
                )

        client = object.__new__(deribit.WsClient)
        client.sock = OversizedHeaderSocket()
        client._buf = b""
        client._frag = None
        client._frag_size = 0

        with self.assertRaisesRegex(
            ConnectionError,
            "websocket message exceeds",
        ):
            client.recv_message(0.01)
        self.assertEqual(client.sock.calls, 1)

    def test_websocket_rejects_cumulative_fragments_before_next_body(
        self,
    ) -> None:
        class FragmentHeaderSocket:
            def __init__(self) -> None:
                self.actions = [
                    b"\x01\x03abc",
                    b"\x80\x03",
                ]
                self.calls = 0

            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                self.calls += 1
                if not self.actions:
                    raise AssertionError(
                        "overflowing continuation body was requested"
                    )
                return self.actions.pop(0)

        client = object.__new__(deribit.WsClient)
        client.sock = FragmentHeaderSocket()
        client._buf = b""
        client._frag = None
        client._frag_size = 0

        with (
            mock.patch.object(
                deribit,
                "WS_MAX_MESSAGE_BYTES",
                5,
                create=True,
            ),
            self.assertRaisesRegex(
                ConnectionError,
                "websocket message exceeds",
            ),
        ):
            client.recv_message(0.01)
        self.assertEqual(client.sock.calls, 2)

    def test_websocket_accepts_exact_message_size_boundary(self) -> None:
        class NoReadSocket:
            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                raise AssertionError("fully buffered frame read socket")

        with mock.patch.object(
            deribit,
            "WS_MAX_MESSAGE_BYTES",
            5,
            create=True,
        ):
            single = object.__new__(deribit.WsClient)
            single.sock = NoReadSocket()
            single._buf = b"\x81\x05hello"
            single._frag = None
            single._frag_size = 0
            self.assertEqual(
                single.recv_message(0.01),
                (deribit.OP_TEXT, b"hello"),
            )

            fragmented = object.__new__(deribit.WsClient)
            fragmented.sock = NoReadSocket()
            fragmented._buf = b"\x01\x02he\x80\x03llo"
            fragmented._frag = None
            fragmented._frag_size = 0
            self.assertEqual(
                fragmented.recv_message(0.01),
                (deribit.OP_TEXT, b"hello"),
            )

    def test_websocket_total_deadline_bounds_ping_flood(self) -> None:
        class PingFloodSocket:
            def __init__(self) -> None:
                self.recv_calls = 0
                self.sent = 0

            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                self.recv_calls += 1
                if self.recv_calls > 8:
                    raise AssertionError(
                        "ping flood escaped total receive deadline"
                    )
                return b"\x89\x00"

            def sendall(self, _payload: bytes) -> None:
                self.sent += 1

        client = object.__new__(deribit.WsClient)
        client.sock = PingFloodSocket()
        client._buf = b""
        client._frag = None
        client._frag_size = 0
        clock = {"now": 0.0}

        def advance_clock() -> float:
            clock["now"] += 0.04
            return clock["now"]

        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                side_effect=advance_clock,
            ),
            self.assertRaises(socket.timeout),
        ):
            client.recv_message(0.1)
        self.assertLessEqual(client.sock.recv_calls, 2)
        self.assertGreaterEqual(client.sock.sent, 1)

    def test_websocket_total_deadline_resumes_trickled_frame(
        self,
    ) -> None:
        class TrickleSocket:
            def __init__(self) -> None:
                self.actions = [
                    b"\x81\x05h",
                    b"e",
                    b"l",
                    b"l",
                    b"o",
                ]

            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                if not self.actions:
                    raise AssertionError("trickle frame over-read")
                return self.actions.pop(0)

        client = object.__new__(deribit.WsClient)
        client.sock = TrickleSocket()
        client._buf = b""
        client._frag = None
        client._frag_size = 0
        clock = {"now": 0.0}

        def advance_clock() -> float:
            clock["now"] += 0.04
            return clock["now"]

        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                side_effect=advance_clock,
            ),
            self.assertRaises(socket.timeout),
        ):
            client.recv_message(0.1)
        self.assertEqual(client._buf, b"\x81\x05he")

        with mock.patch.object(
            deribit.time,
            "monotonic",
            side_effect=lambda: clock["now"],
        ):
            self.assertEqual(
                client.recv_message(1.0),
                (deribit.OP_TEXT, b"hello"),
            )

    def test_pre_ready_oversized_frame_disconnects_without_market_state(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.max_minutes = 1.0

        class OversizedWebSocket:
            def __init__(self) -> None:
                self.calls = 0

            def sendall(self, _payload: bytes) -> None:
                return None

            def settimeout(self, _timeout: float) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                self.calls += 1
                if self.calls == 1:
                    collector.stop = True
                    return (
                        b"\x81\x7f"
                        + deribit.struct.pack(
                            ">Q",
                            deribit.WS_MAX_MESSAGE_BYTES + 1,
                        )
                    )
                return b""

            def close(self) -> None:
                return None

        socket_peer = OversizedWebSocket()
        websocket = object.__new__(deribit.WsClient)
        websocket.sock = socket_peer
        websocket._buf = b""
        websocket._frag = None
        websocket._frag_size = 0
        stdout = io.StringIO()
        with (
            mock.patch.object(collector, "seed_instruments"),
            mock.patch.object(
                deribit,
                "WsClient",
                return_value=websocket,
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
        ):
            rc = collector.run()

        records = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        lifecycle = [
            record["event"]
            for record in records
            if record.get("event")
            in {"connection_attempt", "connection_connect",
                "connection_disconnect", "connection_ready"}
        ]
        self.assertEqual(rc, 0)
        self.assertNotIn(deribit.READY_BANNER, stdout.getvalue())
        self.assertEqual(socket_peer.calls, 1)
        self.assertEqual(
            lifecycle,
            [
                "connection_attempt",
                "connection_connect",
                "connection_disconnect",
            ],
        )
        self.assertFalse(
            any(
                record.get("record_type")
                in {"option_trade", "option_ticker", "index_price"}
                for record in records
            )
        )

    def test_deribit_dedup_commits_only_after_durable_trade_write(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        trade = {
            **self.trade_payload(1),
            "timestamp": received_ms - 1,
        }
        dedup_id = deribit.sha256_id(
            f"trade|{self.instrument_name}|trade-1|1|{received_ms - 1}"
        )
        with (
            mock.patch.object(
                collector.sink,
                "write",
                side_effect=OSError("injected write failure"),
            ),
            self.assertRaises(OSError),
        ):
            collector.ingest_trade(trade, received_ms)
        self.assertNotIn(dedup_id, collector.deduper.seen)

        collector.ingest_trade(trade, received_ms)
        collector.ingest_trade(trade, received_ms)
        trade_path = (
            out / f"day{deribit.utc_day(received_ms)}" / "trades.jsonl"
        )
        self.assertEqual(
            len(trade_path.read_text(encoding="utf-8").splitlines()),
            1,
        )

    def test_catastrophic_chain_truncation_is_quarantined(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        baseline_rows = [
            self.instrument_payload(index)
            for index in range(100)
        ]
        collector.instruments = {}
        collector._persisted_instrument_fingerprints = {}
        collector._apply_instrument_snapshot(
            {"btc": baseline_rows},
            received_ms,
        )
        collector._last_trade_seq = {
            row["instrument_name"]: index
            for index, row in enumerate(baseline_rows, start=1)
        }
        collector.active_channels = {
            f"ticker.{row['instrument_name']}.agg2"
            for row in baseline_rows
        }
        collector._last_authoritative_chain_refresh_monotonic = 123.0
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        before_bytes = instrument_path.read_bytes()
        before_instruments = dict(collector.instruments)
        before_fingerprints = dict(
            collector._persisted_instrument_fingerprints
        )
        before_sequences = dict(collector._last_trade_seq)
        before_channels = set(collector.active_channels)

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=[baseline_rows[0]],
        ):
            collector._chain_refresh_worker(
                collector._chain_refresh_generation,
                None,
            )
        with mock.patch.object(
            deribit.time,
            "monotonic",
            return_value=500.0,
        ):
            collector._poll_chain_refresh(received_ms + 1)

        self.assertEqual(instrument_path.read_bytes(), before_bytes)
        self.assertEqual(collector.instruments, before_instruments)
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            before_fingerprints,
        )
        self.assertEqual(collector._last_trade_seq, before_sequences)
        self.assertEqual(collector.active_channels, before_channels)
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            123.0,
        )
        self.assertEqual(
            collector._next_chain_refresh_attempt,
            500.0 + deribit.CHAIN_REFRESH_RETRY_SECONDS,
        )
        events = [
            json.loads(line)
            for line in instrument_path.with_name(
                "events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        quarantines = [
            row
            for row in events
            if row.get("event") == "instrument_snapshot_quarantined"
        ]
        self.assertEqual(len(quarantines), 1)
        self.assertEqual(
            quarantines[0]["detail"]["unexpected_missing"],
            99,
        )

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=baseline_rows,
        ):
            collector._chain_refresh_worker(
                collector._chain_refresh_generation,
                None,
            )
        with mock.patch.object(
            deribit.time,
            "monotonic",
            return_value=501.0,
        ):
            collector._poll_chain_refresh(received_ms + 2)
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            501.0,
        )
        self.assertEqual(collector.instruments, before_instruments)

    def test_unexpected_delisting_requires_independent_confirmation(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        baseline_rows = [
            self.instrument_payload(index)
            for index in range(10)
        ]
        collector.instruments = {}
        collector._persisted_instrument_fingerprints = {}
        collector._apply_instrument_snapshot(
            {"btc": baseline_rows},
            received_ms,
        )
        candidate = baseline_rows[1:]

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=candidate,
        ):
            suspect = collector._fetch_instrument_snapshot(None)
        with self.assertRaises(
            deribit.InstrumentSnapshotError
        ) as first:
            collector._validate_instrument_snapshot_continuity(
                suspect,
                received_ms + 1,
            )
        self.assertEqual(
            first.exception.error_type,
            "UnexpectedRemovalConfirmationRequired",
        )
        self.assertEqual(
            set(collector.instruments),
            {row["instrument_name"] for row in baseline_rows},
        )

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=baseline_rows,
        ):
            full = collector._fetch_instrument_snapshot(None)
        collector._validate_instrument_snapshot_continuity(
            full,
            received_ms + 2,
        )
        self.assertEqual(full["btc"], baseline_rows)

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=candidate,
        ):
            suspect = collector._fetch_instrument_snapshot(None)
            with self.assertRaises(
                deribit.InstrumentSnapshotError
            ):
                collector._validate_instrument_snapshot_continuity(
                    suspect,
                    received_ms + 3,
                )
            confirmed = collector._fetch_instrument_snapshot(None)
        collector._validate_instrument_snapshot_continuity(
            confirmed,
            received_ms + 4,
        )
        collector._apply_instrument_snapshot(
            confirmed,
            received_ms + 1,
        )
        self.assertEqual(
            set(collector.instruments),
            {row["instrument_name"] for row in candidate},
        )
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        removals = [
            json.loads(line)
            for line in instrument_path.read_text(
                encoding="utf-8"
            ).splitlines()
            if json.loads(line).get("record_type")
            == "instrument_removed"
        ]
        self.assertEqual(len(removals), 1)
        self.assertEqual(
            removals[0]["instrument_name"],
            baseline_rows[0]["instrument_name"],
        )

    def test_explained_instrument_removals_apply_immediately(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        baseline_rows = [
            self.instrument_payload(index)
            for index in range(12)
        ]
        baseline_rows[0] = self.instrument_payload(
            0,
            expiration_timestamp=(
                received_ms
                - deribit.INSTRUMENT_EXPIRY_SETTLEMENT_GRACE_MS
                - 1
            ),
        )
        baseline_rows[1] = self.instrument_payload(
            1,
            state="closed",
            is_active=False,
        )
        collector.instruments = {}
        collector._persisted_instrument_fingerprints = {}
        collector._apply_instrument_snapshot(
            {"btc": baseline_rows},
            received_ms,
        )
        candidate = baseline_rows[2:]
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=candidate,
        ):
            accepted = collector._fetch_instrument_snapshot(None)
        collector._validate_instrument_snapshot_continuity(
            accepted,
            received_ms + 1,
        )
        collector._apply_instrument_snapshot(
            accepted,
            received_ms + 1,
        )
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        records = [
            json.loads(line)
            for line in instrument_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        removals = [
            row
            for row in records
            if row.get("record_type") == "instrument_removed"
        ]
        self.assertEqual(len(removals), 2)
        self.assertEqual(
            {row["instrument_name"] for row in removals},
            {
                baseline_rows[0]["instrument_name"],
                baseline_rows[1]["instrument_name"],
            },
        )

    def test_first_chain_seed_requires_plausible_confirmed_metadata(
        self,
    ) -> None:
        def new_collector():
            out = Path(tempfile.mkdtemp(
                prefix="collector-integrity-initial-chain-"
            ))
            instance = deribit.Collector(
                out,
                currencies=("btc",),
                day_cap=10_000_000,
                total_cap=100_000_000,
                heartbeat=False,
            )
            self.addCleanup(instance.sink.close)
            return instance

        valid_rows = [
            self.instrument_payload(index)
            for index in range(10)
        ]
        collector = new_collector()
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            side_effect=[valid_rows, [dict(row) for row in valid_rows]],
        ) as fetch:
            fetched = collector._fetch_instrument_snapshot(None)
        self.assertEqual(fetch.call_count, 2)
        self.assertEqual(fetched["btc"], valid_rows)
        mismatched_rows = [dict(row) for row in valid_rows]
        mismatched_rows[0]["tick_size"] *= 2
        mismatched = new_collector()
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                side_effect=[valid_rows, mismatched_rows],
            ),
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as confirmation,
        ):
            mismatched._fetch_instrument_snapshot(None)
        self.assertEqual(
            confirmation.exception.error_type,
            "InitialSnapshotConfirmationMismatch",
        )


        too_small = new_collector()
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=valid_rows[:9],
            ),
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as small,
        ):
            too_small._fetch_instrument_snapshot(None)
        self.assertEqual(
            small.exception.error_type,
            "InitialSnapshotImplausible",
        )

        incomplete_rows = [dict(row) for row in valid_rows]
        del incomplete_rows[0]["strike"]
        incomplete = new_collector()
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=incomplete_rows,
            ),
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as schema,
        ):
            incomplete._fetch_instrument_snapshot(None)
        self.assertEqual(schema.exception.error_type, "SchemaMismatch")

    def test_instrument_snapshot_requires_reconstruction_metadata(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        missing = object()
        numeric_fields = (
            "contract_size",
            "min_trade_amount",
            "lot_size",
            "tick_size",
        )
        numeric_invalid = (
            missing,
            None,
            False,
            0,
            -1,
            "1",
            float("nan"),
            float("inf"),
            float("-inf"),
        )
        string_invalid = {
            "settlement_currency": (
                missing,
                None,
                "",
                "btc",
                "USD",
                1,
            ),
            "price_index": (
                missing,
                None,
                "",
                "BTC_USD",
                "eth_usd",
                1,
            ),
        }

        cases = [
            (field, value)
            for field in numeric_fields
            for value in numeric_invalid
        ]
        cases.extend(
            (field, value)
            for field, values in string_invalid.items()
            for value in values
        )
        for field, value in cases:
            with self.subTest(field=field, value=repr(value)):
                row = self.instrument_payload()
                if value is missing:
                    del row[field]
                else:
                    row[field] = value
                with (
                    mock.patch.object(
                        deribit,
                        "fetch_instruments",
                        return_value=[row],
                    ),
                    self.assertRaises(
                        deribit.InstrumentSnapshotError
                    ) as raised,
                ):
                    collector._fetch_currency_instrument_rows(
                        "btc",
                        deadline=None,
                    )
                self.assertEqual(
                    raised.exception.error_type,
                    "SchemaMismatch",
                )

        for currency in ("btc", "eth"):
            with self.subTest(valid_currency=currency):
                row = self.instrument_payload(currency=currency)
                with mock.patch.object(
                    deribit,
                    "fetch_instruments",
                    return_value=[row],
                ):
                    self.assertEqual(
                        collector._fetch_currency_instrument_rows(
                            currency,
                            deadline=None,
                        ),
                        [row],
                    )

    def test_synchronous_seed_uses_post_fetch_completion_clock(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-seed-completion-clock-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        rows = [self.instrument_payload(index) for index in range(10)]
        completion = deribit.datetime(
            2026,
            2,
            2,
            0,
            0,
            0,
            100_000,
            tzinfo=deribit.timezone.utc,
        )
        completion_ms = int(completion.timestamp() * 1000)
        fetch_complete = False

        def complete_fetch():
            nonlocal fetch_complete
            fetch_complete = True
            return {"btc": rows}

        def completion_clock():
            self.assertTrue(
                fetch_complete,
                "wall time sampled before full snapshot completed",
            )
            return completion_ms / 1000.0

        with (
            mock.patch.object(
                collector,
                "_fetch_instrument_snapshot_interruptibly",
                side_effect=complete_fetch,
            ),
            mock.patch.object(
                deribit.time,
                "time",
                side_effect=completion_clock,
            ) as wall_clock,
        ):
            collector.seed_instruments()

        self.assertGreaterEqual(wall_clock.call_count, 1)
        self.assertEqual(
            collector._instrument_persist_day,
            "20260202",
        )
        prior_path = out / "day20260201" / "instruments.jsonl"
        self.assertFalse(prior_path.exists())
        instrument_path = out / "day20260202" / "instruments.jsonl"
        records = [
            deribit.strict_json_loads(line)
            for line in instrument_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(len(records), len(rows) + 1)
        self.assertEqual(
            records[-1]["record_type"],
            "instrument_snapshot_commit",
        )
        self.assertEqual(records[-1]["member_count"], len(rows))
        self.assertEqual(records[-1]["currencies"], ["btc"])
        self.assertEqual(
            {record["received_ts_ms"] for record in records},
            {completion_ms},
        )
        expected_snapshot_id = deribit.sha256_id(
            f"instrument-snapshot|20260202|{completion_ms}"
        )
        self.assertEqual(
            {record["snapshot_id"] for record in records},
            {expected_snapshot_id},
        )
        events = [
            deribit.strict_json_loads(line)
            for line in instrument_path.with_name(
                "events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        seed_done = [
            event
            for event in events
            if event.get("event") == "instrument_seed_done"
        ]
        self.assertEqual(len(seed_done), 1)
        self.assertEqual(
            seed_done[0]["received_ts_ms"],
            completion_ms,
        )
        self.assertEqual(
            deribit.utc_day(completion_ms),
            collector._instrument_persist_day,
        )

    def test_subscription_day_fence_persists_metadata_before_market_row(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        row = self.instrument_payload(name=self.instrument_name)
        before_midnight = deribit.datetime(
            2026,
            2,
            1,
            23,
            59,
            59,
            900_000,
            tzinfo=deribit.timezone.utc,
        )
        before_ms = int(before_midnight.timestamp() * 1000)
        after_ms = before_ms + 200
        after_day = deribit.utc_day(after_ms)
        collector._apply_instrument_snapshot({"btc": [row]}, before_ms)
        collector._last_authoritative_chain_refresh_monotonic = 77.0
        collector.last_utc_day = deribit.utc_day(before_ms)
        collector._next_chain_refresh_attempt = float("inf")
        collector.active_channels = {
            f"ticker.{self.instrument_name}.agg2"
        }
        instrument_path = out / f"day{after_day}" / "instruments.jsonl"
        metadata_counts_at_ticker_ingest = []
        ingest_ticker = collector.ingest_ticker

        def ingest_after_metadata(*args, **kwargs):
            metadata_counts_at_ticker_ingest.append(
                len(
                    instrument_path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                )
                if instrument_path.exists()
                else 0
            )
            return ingest_ticker(*args, **kwargs)

        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                return_value=123.0,
            ),
            mock.patch.object(
                collector,
                "ingest_ticker",
                side_effect=ingest_after_metadata,
            ),
        ):
            for offset in (0, 1):
                collector._handle_notification(
                    {
                        "method": "subscription",
                        "params": {
                            "channel": (
                                f"ticker.{self.instrument_name}.agg2"
                            ),
                            "data": self.ticker_payload(
                                timestamp=after_ms + offset
                            ),
                        },
                    },
                    after_ms + offset,
                )

        self.assertEqual(metadata_counts_at_ticker_ingest, [2, 2])
        instrument_rows = [
            deribit.strict_json_loads(line)
            for line in instrument_path.read_text(
                encoding="utf-8"
            ).splitlines()
        ]
        self.assertEqual(len(instrument_rows), 2)
        self.assertEqual(
            instrument_rows[-1]["record_type"],
            "instrument_snapshot_commit",
        )
        self.assertEqual(
            {record["received_ts_ms"] for record in instrument_rows},
            {after_ms},
        )
        self.assertEqual(collector._instrument_persist_day, after_day)
        self.assertEqual(collector.last_utc_day, after_day)
        self.assertLessEqual(collector._next_chain_refresh_attempt, 123.0)
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            77.0,
        )

    def test_subscription_day_fence_failure_precedes_market_state(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        row = self.instrument_payload(name=self.instrument_name)
        before_ms = int(
            deribit.datetime(
                2026,
                2,
                1,
                23,
                59,
                59,
                900_000,
                tzinfo=deribit.timezone.utc,
            ).timestamp()
            * 1000
        )
        after_ms = before_ms + 200
        collector._apply_instrument_snapshot({"btc": [row]}, before_ms)
        collector.active_channels = {"trades.option.BTC.100ms"}
        prior_day = collector._instrument_persist_day
        trade = self.trade_payload(101)
        original_write_batch = collector.sink.write_batch

        def reject_new_day(kind, records, now_ms=None):
            if kind == "instruments" and deribit.utc_day(now_ms) != prior_day:
                raise deribit.DayCapError("injected day-fence cap")
            return original_write_batch(kind, records, now_ms=now_ms)

        with (
            mock.patch.object(
                collector.sink,
                "write_batch",
                side_effect=reject_new_day,
            ),
            self.assertRaises(deribit.DayCapError),
        ):
            collector._handle_notification(
                {
                    "method": "subscription",
                    "params": {
                        "channel": "trades.option.BTC.100ms",
                        "data": [trade],
                    },
                },
                after_ms,
            )

        self.assertEqual(collector._instrument_persist_day, prior_day)
        self.assertEqual(collector.msg_seq, 0)
        self.assertEqual(collector.last_message_ms, 0)
        self.assertEqual(collector._last_trade_seq, {})
        self.assertEqual(collector.deduper.seen, set())
        self.assertEqual(collector.stats["trades_pushes"], 0)
        self.assertEqual(
            collector._required_last_received_monotonic,
            {},
        )
        trade_path = (
            out
            / f"day{deribit.utc_day(after_ms)}"
            / "trades.jsonl"
        )
        self.assertFalse(
            trade_path.exists() and trade_path.stat().st_size
        )

    def test_day_fence_requires_every_configured_currency(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-day-fence-currencies-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        before_ms = 1_769_990_399_900
        after_ms = before_ms + 200
        collector.instruments = {
            self.instrument_name: self.instrument_payload(
                name=self.instrument_name
            )
        }
        collector._instrument_persist_day = deribit.utc_day(before_ms)

        with self.assertRaises(
            deribit.InstrumentSnapshotError
        ) as raised:
            collector._ensure_instrument_day_fence(
                after_ms,
                received_monotonic=123.0,
            )

        self.assertEqual(raised.exception.currency, "eth")
        self.assertEqual(
            raised.exception.error_type,
            "CachedSnapshotIncomplete",
        )
        self.assertEqual(
            collector._instrument_persist_day,
            deribit.utc_day(before_ms),
        )
        self.assertFalse(
            (
                out
                / f"day{deribit.utc_day(after_ms)}"
                / "instruments.jsonl"
            ).exists()
        )

    def test_ready_day_fence_carries_metadata_and_keeps_refresh_due(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ready-day-fence-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
            max_minutes=1.0,
        )
        row = self.instrument_payload(name=self.instrument_name)
        before_ms = int(
            deribit.datetime(
                2026,
                2,
                1,
                23,
                59,
                59,
                900_000,
                tzinfo=deribit.timezone.utc,
            ).timestamp()
            * 1000
        )
        after_ms = before_ms + 200
        clock = {"ms": before_ms}
        observed = {}

        class ConnectedWebSocket:
            def close(self) -> None:
                return None

        def seed() -> None:
            collector._apply_instrument_snapshot(
                {"btc": [row]},
                before_ms,
            )
            collector._last_authoritative_chain_refresh_monotonic = 55.0

        def cross_midnight_during_coverage() -> bool:
            clock["ms"] = after_ms
            return True

        def inspect_ready_boundary() -> None:
            path = (
                out
                / f"day{deribit.utc_day(after_ms)}"
                / "instruments.jsonl"
            )
            observed["metadata_rows"] = (
                len(path.read_text(encoding="utf-8").splitlines())
                if path.exists()
                else 0
            )
            observed["next_refresh"] = (
                collector._next_chain_refresh_attempt
            )
            observed["last_message_monotonic"] = (
                collector.last_message_monotonic
            )
            observed["freshness"] = (
                collector._last_authoritative_chain_refresh_monotonic
            )
            collector.stop = True

        stdout = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                deribit.time,
                "time",
                side_effect=lambda: clock["ms"] / 1000.0,
            ),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                return_value=ConnectedWebSocket(),
            ),
            mock.patch.object(
                collector,
                "seed_instruments",
                side_effect=seed,
            ),
            mock.patch.object(
                collector,
                "subscribe_all",
                return_value=True,
            ),
            mock.patch.object(
                collector,
                "await_required_coverage",
                side_effect=cross_midnight_during_coverage,
            ),
            mock.patch.object(
                collector,
                "_loop",
                side_effect=inspect_ready_boundary,
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
        ):
            result = collector.run()

        self.assertEqual(result, 0)
        self.assertEqual(observed["metadata_rows"], 2)
        self.assertLessEqual(
            observed["next_refresh"],
            observed["last_message_monotonic"],
        )
        self.assertEqual(observed["freshness"], 55.0)
        self.assertIn(deribit.READY_BANNER, stdout.getvalue())
        ready_events = [
            deribit.strict_json_loads(line)
            for line in (
                out
                / f"day{deribit.utc_day(after_ms)}"
                / "events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
            if deribit.strict_json_loads(line).get("event")
            == "connection_ready"
        ]
        self.assertEqual(len(ready_events), 1)

    def test_ready_day_fence_cap_never_announces_ready(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ready-day-fence-cap-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
            max_minutes=1.0,
        )
        row = self.instrument_payload(name=self.instrument_name)
        before_ms = int(
            deribit.datetime(
                2026,
                2,
                1,
                23,
                59,
                59,
                900_000,
                tzinfo=deribit.timezone.utc,
            ).timestamp()
            * 1000
        )
        after_ms = before_ms + 200
        clock = {"ms": before_ms}
        original_write_batch = collector.sink.write_batch

        class ConnectedWebSocket:
            def close(self) -> None:
                return None

        def seed() -> None:
            collector._apply_instrument_snapshot(
                {"btc": [row]},
                before_ms,
            )

        def cross_midnight_during_coverage() -> bool:
            clock["ms"] = after_ms
            return True

        def reject_new_day(kind, records, now_ms=None):
            if (
                kind == "instruments"
                and deribit.utc_day(now_ms) != deribit.utc_day(before_ms)
            ):
                raise deribit.DayCapError("injected ready fence cap")
            return original_write_batch(kind, records, now_ms=now_ms)

        stdout = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                deribit.time,
                "time",
                side_effect=lambda: clock["ms"] / 1000.0,
            ),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                return_value=ConnectedWebSocket(),
            ),
            mock.patch.object(
                collector,
                "seed_instruments",
                side_effect=seed,
            ),
            mock.patch.object(
                collector,
                "subscribe_all",
                return_value=True,
            ),
            mock.patch.object(
                collector,
                "await_required_coverage",
                side_effect=cross_midnight_during_coverage,
            ),
            mock.patch.object(
                collector.sink,
                "write_batch",
                side_effect=reject_new_day,
            ),
            mock.patch.object(
                collector,
                "_loop",
                side_effect=AssertionError("READY loop reached"),
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
        ):
            result = collector.run()

        self.assertEqual(result, 3)
        self.assertNotIn(deribit.READY_BANNER, stdout.getvalue())
        self.assertFalse(collector._ready_printed)
        self.assertFalse(
            any(
                row.get("event") == "connection_ready"
                for path in out.glob("day*/events.jsonl")
                for row in (
                    deribit.strict_json_loads(line)
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                )
            )
        )

    def test_synchronous_seed_failure_uses_observation_clock(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.instruments = {}
        collector._persisted_instrument_fingerprints = {}
        failure = deribit.datetime(
            2026,
            2,
            2,
            0,
            0,
            1,
            tzinfo=deribit.timezone.utc,
        )
        failure_ms = int(failure.timestamp() * 1000)
        fetch_failed = False

        def fail_fetch():
            nonlocal fetch_failed
            fetch_failed = True
            raise deribit.InstrumentSnapshotError(
                "btc",
                "TimeoutError",
                "injected delayed failure",
            )

        def failure_clock():
            self.assertTrue(
                fetch_failed,
                "failure time sampled before fetch failed",
            )
            return failure_ms / 1000.0

        with (
            mock.patch.object(
                collector,
                "_fetch_instrument_snapshot_interruptibly",
                side_effect=fail_fetch,
            ),
            mock.patch.object(
                deribit.time,
                "time",
                side_effect=failure_clock,
            ),
            self.assertRaises(deribit.InstrumentSnapshotError),
        ):
            collector.seed_instruments()

        self.assertEqual(collector.instruments, {})
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            0.0,
        )
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        failures = [
            event
            for event in events
            if event.get("event") == "instrument_seed_error"
        ]
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0]["received_ts_ms"], failure_ms)
        self.assertEqual(
            failures[0]["detail"]["error_type"],
            "TimeoutError",
        )
        self.assertFalse(
            any(
                record.get("record_type") == "instrument"
                for path in out.glob("day*/instruments.jsonl")
                for record in (
                    deribit.strict_json_loads(line)
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                )
            )
        )

    def test_invalid_reconstruction_metadata_cannot_seed_or_ready(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-instrument-metadata-seed-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        rows = [self.instrument_payload(index) for index in range(10)]
        del rows[0]["contract_size"]

        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=rows,
            ),
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as raised,
        ):
            collector.seed_instruments()

        self.assertEqual(raised.exception.error_type, "SchemaMismatch")
        self.assertEqual(collector.instruments, {})
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(collector.active_channels, set())
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            0.0,
        )
        self.assertFalse(collector._ready_printed)
        instrument_records = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/instruments.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertFalse(
            any(
                row.get("record_type")
                in {"instrument", "instrument_removed"}
                for row in instrument_records
            ),
            "invalid seed persisted instrument metadata",
        )

    def test_invalid_reconstruction_metadata_refresh_preserves_chain(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_000
        baseline = [
            self.instrument_payload(index)
            for index in range(10)
        ]
        collector.instruments = {}
        collector._persisted_instrument_fingerprints = {}
        collector._apply_instrument_snapshot(
            {"btc": baseline},
            received_ms,
        )
        channels = {
            "trades.option.BTC.100ms",
            *(
                f"ticker.{row['instrument_name']}.agg2"
                for row in baseline
            ),
        }
        collector.active_channels = channels
        collector._required_last_received_monotonic = {
            channel: 100.0
            for channel in channels
            if channel.startswith("ticker.")
        }
        collector._required_last_source_timestamp_ms = {
            channel: received_ms
            for channel in channels
            if channel.startswith("ticker.")
        }
        collector._ticker_last_bucket = {
            row["instrument_name"]: 123
            for row in baseline
        }
        collector._last_authoritative_chain_refresh_monotonic = 101.0
        before_instruments = dict(collector.instruments)
        before_fingerprints = dict(
            collector._persisted_instrument_fingerprints
        )
        before_channels = set(collector.active_channels)
        before_received = dict(
            collector._required_last_received_monotonic
        )
        before_source = dict(
            collector._required_last_source_timestamp_ms
        )
        before_buckets = dict(collector._ticker_last_bucket)
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        before_bytes = instrument_path.read_bytes()
        invalid = [dict(row) for row in baseline]
        invalid[0]["tick_size"] = 0

        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=invalid,
        ):
            collector._chain_refresh_worker(
                collector._chain_refresh_generation,
                None,
            )
        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                return_value=500.0,
            ),
            mock.patch.object(
                collector,
                "reconcile_chain",
            ) as reconcile,
        ):
            collector._poll_chain_refresh(received_ms + 1)

        reconcile.assert_not_called()
        self.assertEqual(instrument_path.read_bytes(), before_bytes)
        self.assertEqual(collector.instruments, before_instruments)
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            before_fingerprints,
        )
        self.assertEqual(collector.active_channels, before_channels)
        self.assertEqual(
            collector._required_last_received_monotonic,
            before_received,
        )
        self.assertEqual(
            collector._required_last_source_timestamp_ms,
            before_source,
        )
        self.assertEqual(
            collector._ticker_last_bucket,
            before_buckets,
        )
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            101.0,
        )
        self.assertEqual(
            collector._next_chain_refresh_attempt,
            500.0 + deribit.CHAIN_REFRESH_RETRY_SECONDS,
        )
        events = [
            deribit.strict_json_loads(line)
            for line in instrument_path.with_name(
                "events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        quarantines = [
            event
            for event in events
            if event.get("event") == "instrument_snapshot_quarantined"
        ]
        self.assertEqual(len(quarantines), 1)
        self.assertEqual(
            quarantines[0]["detail"]["error_type"],
            "SchemaMismatch",
        )
        self.assertEqual(
            quarantines[0]["detail"]["candidate_count"],
            10,
        )
        self.assertNotIn("row", quarantines[0]["detail"])

    def test_instrument_persist_marker_follows_durable_write(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-seed-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        new_name = "BTC-30DEC30-90000-C"
        rows = [
            self.instrument_payload(index)
            for index in range(10)
        ]
        rows[0] = self.instrument_payload(0, name=new_name)
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=rows,
            ),
            mock.patch.object(
                deribit.time,
                "time",
                return_value=1_800_000_000.0,
            ),
        ):
            with (
                mock.patch.object(
                    collector.sink,
                    "write_batch",
                    side_effect=OSError("injected write failure"),
                ),
                self.assertRaises(OSError),
            ):
                collector.seed_instruments()
            self.assertNotIn(
                new_name,
                collector._persisted_instrument_fingerprints,
            )
            collector.seed_instruments()
        self.assertIn(
            new_name,
            collector._persisted_instrument_fingerprints,
        )
        instrument_path = (
            out
            / f"day{deribit.utc_day(1_800_000_000_000)}"
            / "instruments.jsonl"
        )
        self.assertEqual(
            len(instrument_path.read_text(encoding="utf-8").splitlines()),
            len(rows) + 1,
        )

    def test_max_minutes_bounds_startup_seed_retries(self) -> None:
        collector, _ = self.make_collector()
        collector.max_minutes = 0.001
        clock = {"now": 0.0}
        calls = {"fetch": 0}

        class RetryOverflow(BaseException):
            pass

        def fail_fetch(_currency: str, deadline=None):
            del deadline
            calls["fetch"] += 1
            if calls["fetch"] > 3:
                raise RetryOverflow
            raise ConnectionError("injected REST failure")

        def advance(seconds: float) -> None:
            clock["now"] += seconds

        with (
            mock.patch.object(deribit.time, "monotonic", side_effect=lambda: clock["now"]),
            mock.patch.object(deribit, "stop_aware_sleep", side_effect=advance),
            mock.patch.object(deribit, "fetch_instruments", side_effect=fail_fetch),
            mock.patch.object(deribit, "WsClient") as ws_client,
        ):
            try:
                rc = collector.run()
            except RetryOverflow:
                self.fail("startup retries ignored the max-minutes deadline")
        self.assertEqual(rc, 0)
        self.assertLessEqual(clock["now"], collector.max_minutes * 60)
        ws_client.assert_not_called()

    def test_rest_call_enforces_per_attempt_deadline_under_trickle(
        self,
    ) -> None:
        response = FakeHttpResponse([b"{"] + [b" "] * 120, delay=0.005)
        started = time.monotonic()
        with (
            mock.patch.object(
                deribit,
                "REST_ATTEMPT_TIMEOUT_SECONDS",
                0.05,
            ),
            mock.patch.object(deribit, "urlopen", return_value=response),
            self.assertRaises(ConnectionError),
        ):
            deribit.rest_call("/public/get_instruments", attempts=1)
        elapsed = time.monotonic() - started
        self.assertLess(elapsed, 0.4)
        self.assertTrue(response.chunks)

    def test_rest_call_bounds_response_size(self) -> None:
        chunk = b"x" * 4096
        response = FakeHttpResponse([chunk] * 8)
        with (
            mock.patch.object(deribit, "REST_MAX_RESPONSE_BYTES", 4096),
            mock.patch.object(deribit, "urlopen", return_value=response),
            self.assertRaises(ConnectionError),
        ):
            deribit.rest_call("/public/get_instruments", attempts=1)
        self.assertTrue(response.chunks)

    def test_unbounded_rest_call_fails_instead_of_hanging(self) -> None:
        responses = [
            FakeHttpResponse([b"{"] + [b" "] * 200, delay=0.005)
            for _ in range(3)
        ]
        started = time.monotonic()
        with (
            mock.patch.object(
                deribit,
                "REST_ATTEMPT_TIMEOUT_SECONDS",
                0.05,
            ),
            mock.patch.object(deribit, "stop_aware_sleep"),
            mock.patch.object(
                deribit,
                "urlopen",
                side_effect=list(responses),
            ),
            self.assertRaises(ConnectionError),
        ):
            deribit.fetch_instruments("btc", deadline=None)
        elapsed = time.monotonic() - started
        self.assertLess(elapsed, 0.9)
        self.assertTrue(all(response.chunks for response in responses))

    def test_chain_refresh_thread_is_replaceable_after_failure(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            side_effect=ConnectionError("injected refresh failure"),
        ):
            self.assertTrue(collector._start_chain_refresh())
            first = collector._chain_refresh_thread
            first.join(timeout=5.0)
            self.assertFalse(first.is_alive())
            self.assertFalse(collector._start_chain_refresh())
            collector._poll_chain_refresh(int(deribit.time.time() * 1000))
            self.assertEqual(
                collector.stats["instrument_refresh_errors"],
                1,
            )
            self.assertTrue(collector._start_chain_refresh())
            second = collector._chain_refresh_thread
            self.assertIsNot(second, first)
            second.join(timeout=5.0)
            self.assertFalse(second.is_alive())

    def trickling_handshake_socket(self, delay: float = 0.01):
        class TricklingHandshakeSocket:
            def __init__(self) -> None:
                self.reads = 0
                self.timeouts: list[float | None] = []
                self.closed = False

            def settimeout(self, timeout: float | None) -> None:
                self.timeouts.append(timeout)

            def gettimeout(self) -> float | None:
                return self.timeouts[-1] if self.timeouts else None

            def sendall(self, _payload: bytes) -> None:
                return None

            def recv(self, _size: int) -> bytes:
                self.reads += 1
                time.sleep(delay)
                return b"X-Pad: " + b"y" * 500 + b"\r\n"

            def close(self) -> None:
                self.closed = True

        return TricklingHandshakeSocket()

    def handshake_socket_patches(self, peer):
        class HandshakeContext:
            def wrap_socket(self, _raw, server_hostname=None):
                del server_hostname
                return peer

        return (
            mock.patch.object(
                deribit.socket,
                "create_connection",
                return_value=peer,
            ),
            mock.patch.object(
                deribit.ssl,
                "create_default_context",
                return_value=HandshakeContext(),
            ),
        )

    def test_websocket_handshake_enforces_total_deadline(self) -> None:
        peer = self.trickling_handshake_socket()
        create_connection, create_context = self.handshake_socket_patches(peer)
        started = time.monotonic()
        with (
            create_connection,
            create_context,
            self.assertRaises(ConnectionError),
        ):
            deribit.WsClient(timeout=0.1)
        elapsed = time.monotonic() - started
        self.assertLess(elapsed, 0.5)
        self.assertLess(
            peer.reads * 507,
            deribit.WS_HANDSHAKE_MAX_HEADER_BYTES,
        )
        self.assertTrue(peer.timeouts)
        self.assertLessEqual(max(peer.timeouts), 0.1)
        self.assertTrue(peer.closed)

    def test_unbounded_connect_surfaces_bounded_handshake_failure(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        self.assertIsNone(collector._deadline_monotonic)
        peer = self.trickling_handshake_socket()
        create_connection, create_context = self.handshake_socket_patches(peer)
        real_client = deribit.WsClient
        started = time.monotonic()
        with (
            create_connection,
            create_context,
            mock.patch.object(
                deribit,
                "WsClient",
                side_effect=lambda timeout=20.0: real_client(timeout=0.1),
            ),
            self.assertRaises(ConnectionError),
        ):
            collector._connect_websocket_interruptibly()
        self.assertLess(time.monotonic() - started, 1.5)

    def test_duplicate_json_member_is_rejected_at_every_depth(self) -> None:
        for payload in (
            '{"a":1,"a":2}',
            '{"o":{"x":1,"x":2}}',
            '[{"x":1,"x":2}]',
            '{"result":[],"result":["ticker.BTC-1JAN27-90000-C.agg2"]}',
        ):
            with self.subTest(payload=payload):
                with self.assertRaises(deribit.DuplicateJsonKeyError):
                    deribit.strict_json_loads(payload)
        self.assertEqual(
            deribit.strict_json_loads('{"a":1,"b":{"c":2}}'),
            {"a": 1, "b": {"c": 2}},
        )

    def test_duplicate_result_ack_never_activates_channels(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.max_minutes = 1.0
        channel = f"ticker.{self.instrument_name}.agg2"

        class RawFrameWebSocket:
            def __init__(self, frames: list[bytes]):
                self.frames = list(frames)

            def send_text(self, _text: str) -> None:
                return None

            def recv_message(self, timeout: float):
                del timeout
                if not self.frames:
                    raise socket.timeout()
                return deribit.OP_TEXT, self.frames.pop(0)

            def close(self) -> None:
                return None

        collector.ws = RawFrameWebSocket([
            (
                '{"jsonrpc":"2.0","id":7,"result":[],'
                f'"result":["{channel}"]}}'
            ).encode("utf-8")
        ])
        with self.assertRaises(deribit.DeribitProtocolFrameError):
            collector._await_ack(7, {channel}, timeout_s=1.0)
        self.assertEqual(collector.active_channels, set())
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            [
                event["detail"]["reason"]
                for event in events
                if event.get("event") == "invalid_jsonrpc_shape"
            ],
            ["duplicate_json_member"],
        )

    def test_duplicate_member_replay_line_is_malformed(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-duplicate-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "instruments.jsonl").write_text(
            '{"record_type":"instrument","record_type":"instrument"}\n',
            encoding="utf-8",
        )
        restarted, stats = self.restart_replay(out)
        self.assertEqual(stats["malformed"], 1)
        self.assertEqual(stats["metadata_replay_rejected"], 0)
        self.assertEqual(restarted.instruments, {})

    def test_duplicate_member_checkpoint_fails_closed(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-duplicate-checkpoint-"
        ))
        out.mkdir(parents=True, exist_ok=True)
        (out / deribit.TRADE_SEQUENCE_STATE_FILE).write_text(
            '{"schema_version":4,"schema_version":4}\n',
            encoding="utf-8",
        )
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        with self.assertRaises(OSError):
            collector.rehydrate_state()
        self.assertIsNotNone(collector._sequence_state_error)

    def test_instrument_expiry_requires_exact_settlement_instant(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        settlement_ms = int(
            deribit.datetime(
                2030,
                12,
                30,
                8,
                0,
                0,
                tzinfo=deribit.timezone.utc,
            ).timestamp()
            * 1000
        )
        self.assertEqual(
            deribit.OPTION_EXPIRY_SETTLEMENT_MS,
            8 * 3_600_000,
        )
        accepted = {
            **self.instrument_payload(0, currency="btc"),
            "expiration_timestamp": settlement_ms,
        }
        self.assertTrue(
            collector._instrument_row_matches_schema(accepted, "btc")
        )
        for label, expiration in (
            ("midnight", settlement_ms - 8 * 3_600_000),
            ("one_ms_early", settlement_ms - 1),
            ("one_ms_late", settlement_ms + 1),
            ("late_same_day", settlement_ms + 15 * 3_600_000 - 1),
        ):
            with self.subTest(offset=label):
                self.assertFalse(
                    collector._instrument_row_matches_schema(
                        {**accepted, "expiration_timestamp": expiration},
                        "btc",
                    )
                )

        baseline_instruments = dict(collector.instruments)
        wrong_time = {
            **accepted,
            "expiration_timestamp": settlement_ms - 1,
        }
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[wrong_time],
            ),
            self.assertRaises(deribit.InstrumentSnapshotError) as failure,
        ):
            collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(failure.exception.error_type, "SchemaMismatch")
        self.assertEqual(collector.instruments, baseline_instruments)

    def test_committed_wrong_expiry_instant_is_replay_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-expiry-instant-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        settlement_ms = int(
            deribit.datetime(
                2030,
                12,
                30,
                8,
                0,
                0,
                tzinfo=deribit.timezone.utc,
            ).timestamp()
            * 1000
        )
        row = {
            **self.instrument_payload(0, currency="btc"),
            "expiration_timestamp": settlement_ms - 1,
        }
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        fingerprint = builder._instrument_fingerprint(row)
        member = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(row, received_ms),
            "received_ts_ms": received_ms,
            "snapshot_id": "wrong-instant-snapshot",
            "instrument_fingerprint": fingerprint,
            "instrument": row,
        }
        commit = builder._instrument_snapshot_commit_record(
            "wrong-instant-snapshot",
            [member],
            received_ms,
        )
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "instruments.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(record, separators=(",", ":"))
                + "\n"
                for record in (member, commit)
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(stats["metadata_replay_rejected"], 2)
        self.assertEqual(stats["malformed"], 0)

    def test_option_contract_size_must_be_one_coin(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        self.assertEqual(deribit.OPTION_CONTRACT_SIZE_COINS, 1)
        base = self.instrument_payload(0, currency="btc")
        for accepted in (1, 1.0):
            with self.subTest(contract_size=repr(accepted)):
                self.assertTrue(
                    collector._instrument_row_matches_schema(
                        {**base, "contract_size": accepted},
                        "btc",
                    )
                )
        for rejected in (0.5, 2, 2.0, 1e-9, True, "1.0", None):
            with self.subTest(contract_size=repr(rejected)):
                self.assertFalse(
                    collector._instrument_row_matches_schema(
                        {**base, "contract_size": rejected},
                        "btc",
                    )
                )

        baseline_instruments = dict(collector.instruments)
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[{**base, "contract_size": 2.0}],
            ),
            self.assertRaises(deribit.InstrumentSnapshotError) as failure,
        ):
            collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(failure.exception.error_type, "SchemaMismatch")
        self.assertEqual(collector.instruments, baseline_instruments)

    def test_committed_wrong_contract_size_is_replay_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-contract-size-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        row = {
            **self.instrument_payload(0, currency="btc"),
            "contract_size": 2.0,
        }
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        member = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(row, received_ms),
            "received_ts_ms": received_ms,
            "snapshot_id": "wrong-size-snapshot",
            "instrument_fingerprint": builder._instrument_fingerprint(row),
            "instrument": row,
        }
        commit = builder._instrument_snapshot_commit_record(
            "wrong-size-snapshot",
            [member],
            received_ms,
        )
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "instruments.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(record, separators=(",", ":"))
                + "\n"
                for record in (member, commit)
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(stats["metadata_replay_rejected"], 2)
        self.assertEqual(stats["malformed"], 0)

    def test_unsolicited_channel_is_protocol_terminal(self) -> None:
        for channel, data in (
            (
                "deribit_price_index.eth_usd",
                {"index_name": "eth_usd", "price": 3000.0, "timestamp": 1},
            ),
            (
                "ticker.ETH-1JAN27-3000-C.agg2",
                self.ticker_payload(),
            ),
            (
                "instrument.creation.BTC-1JAN27-90000-C",
                self.instrument_payload(0, currency="btc"),
            ),
            (
                "instrument.state.BTC-1JAN27-90000-C",
                {"instrument_name": self.instrument_name, "state": "closed"},
            ),
        ):
            with self.subTest(channel=channel):
                collector, out = self.make_collector()
                self.addCleanup(collector.sink.close)
                collector.active_channels = {"trades.option.BTC.100ms"}
                before_instruments = dict(collector.instruments)
                received_ms = int(deribit.time.time() * 1000)
                with self.assertRaises(deribit.DeribitProtocolFrameError):
                    collector._handle_notification(
                        {
                            "method": "subscription",
                            "params": {"channel": channel, "data": data},
                        },
                        received_ms,
                    )
                self.assertEqual(collector.msg_seq, 0)
                self.assertEqual(collector.instruments, before_instruments)
                rows = [
                    deribit.strict_json_loads(line)
                    for path in out.glob("day*/*.jsonl")
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                self.assertFalse(
                    any(
                        row.get("record_type")
                        in {
                            "index_price",
                            "option_ticker",
                            "option_trade",
                            "instrument",
                            "instrument_state",
                        }
                        for row in rows
                    )
                )
                audits = [
                    row["detail"]
                    for row in rows
                    if row.get("event") == "invalid_jsonrpc_shape"
                ]
                self.assertEqual(
                    [audit["reason"] for audit in audits],
                    ["unsolicited_channel"],
                )
                self.assertEqual(audits[0]["channel"], channel)

    def test_pending_ack_channel_notification_is_durable(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.max_minutes = 1.0
        channel = "trades.option.BTC.100ms"
        received_ms = int(deribit.time.time() * 1000)
        trade = self.trade_payload(1, received_ms=received_ms)

        class InterleavedWs:
            def __init__(self):
                self.frames = [
                    json.dumps({
                        "method": "subscription",
                        "params": {"channel": channel, "data": [trade]},
                    }).encode("utf-8"),
                    json.dumps({
                        "jsonrpc": "2.0",
                        "id": 11,
                        "result": [channel],
                    }).encode("utf-8"),
                ]

            def send_text(self, _text: str) -> None:
                return None

            def recv_message(self, timeout: float):
                del timeout
                if not self.frames:
                    raise socket.timeout()
                return deribit.OP_TEXT, self.frames.pop(0)

            def close(self) -> None:
                return None

        collector.ws = InterleavedWs()
        self.assertTrue(
            collector._await_ack(11, {channel}, timeout_s=1.0)
        )
        self.assertEqual(collector.active_channels, {channel})
        self.assertEqual(collector._pending_channels, set())
        trades = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/trades.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(len(trades), 1)
        self.assertEqual(
            trades[0]["trade"]["trade_id"],
            trade["trade_id"],
        )

    def test_post_unsubscribe_frames_are_absorbed_once(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        collector.active_channels = {channel}
        clock = {"now": 1_000.0}

        class AckWs:
            def __init__(self):
                self.frames = [
                    json.dumps({
                        "jsonrpc": "2.0",
                        "id": 21,
                        "result": [channel],
                    }).encode("utf-8")
                ]

            def send_text(self, _text: str) -> None:
                return None

            def recv_message(self, timeout: float):
                del timeout
                if not self.frames:
                    raise socket.timeout()
                return deribit.OP_TEXT, self.frames.pop(0)

            def close(self) -> None:
                return None

        collector.ws = AckWs()
        with mock.patch.object(
            deribit.time,
            "monotonic",
            side_effect=lambda: clock["now"],
        ):
            self.assertTrue(
                collector._await_ack(
                    21,
                    {channel},
                    timeout_s=1.0,
                    method="public/unsubscribe",
                )
            )
            self.assertEqual(collector.active_channels, set())
            self.assertIn(channel, collector._recently_unsubscribed)
            received_ms = int(deribit.time.time() * 1000)
            for offset in (0, 1):
                collector._handle_notification(
                    {
                        "method": "subscription",
                        "params": {
                            "channel": channel,
                            "data": self.ticker_payload(
                                timestamp=received_ms + offset
                            ),
                        },
                    },
                    received_ms + offset,
                )
            self.assertEqual(collector.msg_seq, 0)
            self.assertEqual(collector.stats["protocol_frame_errors"], 0)
            self.assertEqual(
                collector.stats["unsubscribed_grace_frames"],
                2,
            )
            self.assertEqual(collector.last_message_ms, received_ms + 1)
            with self.assertRaises(deribit.DeribitProtocolFrameError):
                collector._handle_notification(
                    {
                        "method": "subscription",
                        "params": {
                            "channel": "deribit_price_index.eth_usd",
                            "data": {"index_name": "eth_usd"},
                        },
                    },
                    received_ms,
                )
            clock["now"] += deribit.UNSUBSCRIBE_GRACE_SECONDS + 1
            with self.assertRaises(deribit.DeribitProtocolFrameError):
                collector._handle_notification(
                    {
                        "method": "subscription",
                        "params": {
                            "channel": channel,
                            "data": self.ticker_payload(
                                timestamp=received_ms + 2
                            ),
                        },
                    },
                    received_ms + 2,
                )
        rows = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertFalse(
            any(row.get("record_type") == "option_ticker" for row in rows)
        )
        grace_audits = [
            row["detail"]
            for row in rows
            if row.get("event") == "unsubscribed_channel_frame"
        ]
        self.assertEqual(
            grace_audits,
            [{"channel": channel}],
        )
        self.assertEqual(
            [
                row["detail"]["reason"]
                for row in rows
                if row.get("event") == "invalid_jsonrpc_shape"
            ],
            ["unsolicited_channel", "unsolicited_channel"],
        )

    def test_subscribe_all_clears_unsubscribe_grace(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector._recently_unsubscribed = {
            "ticker.STALE.agg2": (time.monotonic() + 30.0, False)
        }

        def accept(_request_id: int, channels: set[str], **_kwargs) -> bool:
            collector.active_channels |= channels
            return True

        class SilentWs:
            def send_text(self, _text: str) -> None:
                return None

            def close(self) -> None:
                return None

        collector.ws = SilentWs()
        collector.heartbeat = False
        with mock.patch.object(collector, "_await_ack", side_effect=accept):
            self.assertTrue(collector.subscribe_all())
        self.assertEqual(collector._recently_unsubscribed, {})

    def test_stale_instruments_are_pruned_from_ticker_state(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-ticker-prune-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(3)
        ]
        names = [row["instrument_name"] for row in rows]
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        collector._apply_instrument_snapshot({"btc": rows}, received_ms)
        for offset, name in enumerate(names):
            self.assertTrue(
                collector.ingest_ticker(
                    name,
                    {
                        **self.ticker_payload(
                            timestamp=received_ms + offset
                        ),
                        "instrument_name": name,
                    },
                    received_ms + offset,
                )
            )
        self.assertTrue(
            collector.ingest_index(
                "btc_usd",
                {
                    "index_name": "btc_usd",
                    "price": 90_000.0,
                    "timestamp": received_ms,
                },
                received_ms,
            )
        )
        self.assertEqual(set(collector._ticker_last_bucket), set(names))

        collector._apply_instrument_snapshot(
            {"btc": rows[:2]},
            received_ms + 10,
        )
        self.assertEqual(
            set(collector._ticker_last_bucket),
            set(names[:2]),
        )
        self.assertEqual(
            set(collector._required_last_source_timestamp_ms),
            {f"ticker.{name}.agg2" for name in names[:2]}
            | {"deribit_price_index.btc_usd"},
        )
        collector._sequence_state_dirty = True
        self.assertTrue(collector._save_sequence_state())
        checkpoint = deribit.strict_json_loads(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            set(checkpoint["ticker_buckets"]),
            set(names[:2]),
        )
        self.assertNotIn(
            f"ticker.{names[2]}.agg2",
            checkpoint["required_source_timestamps"],
        )
        collector.sink.close()

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state()
        self.assertTrue(stats["sequence_checkpoint_loaded"])
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=rows[:2],
        ):
            restarted.seed_instruments()
        self.assertEqual(
            set(restarted._ticker_last_bucket),
            set(names[:2]),
        )
        self.assertNotIn(
            f"ticker.{names[2]}.agg2",
            restarted._required_last_source_timestamp_ms,
        )
        self.assertEqual(
            len(restarted._ticker_last_bucket),
            len(restarted.instruments),
        )

    def test_oversized_checkpoint_fails_closed_without_replacing_state(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-checkpoint-size-"
        ))
        self.write_restart_checkpoint(
            out,
            sequences={self.instrument_name: 5},
            ticker_buckets={self.instrument_name: 1},
        )
        state_path = out / deribit.TRADE_SEQUENCE_STATE_FILE
        original = state_path.read_bytes()
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        with mock.patch.object(
            deribit,
            "TRADE_SEQUENCE_STATE_MAX_BYTES",
            8,
        ):
            with self.assertRaises(OSError):
                collector.rehydrate_state()
        self.assertEqual(state_path.read_bytes(), original)

    def test_untrusted_error_payloads_are_bounded_in_audits(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        limit = deribit.AUDIT_TEXT_MAX_CHARS
        received_ms = int(deribit.time.time() * 1000)
        hostile = {
            "code": -32000,
            "message": "x" * 5_000_000,
            "data": {"nested": "y" * 1_000_000},
        }
        with self.assertRaises(ConnectionError) as raised:
            collector._handle_notification(
                {"error": hostile, "id": None},
                received_ms,
            )
        self.assertLessEqual(len(str(raised.exception)), 4 * limit)
        events = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        audits = [
            event["detail"]
            for event in events
            if event.get("event") == "server_error_notification"
        ]
        self.assertEqual(len(audits), 1)
        self.assertEqual(audits[0]["code"], -32000)
        self.assertEqual(len(audits[0]["msg"]), limit)
        self.assertNotIn("data", audits[0])

        collector._log_event(
            "connection_disconnect",
            {"reason": "z" * 5_000_000, "ready": True},
            now_ms=received_ms,
        )
        rows = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/events.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        disconnects = [
            row["detail"]
            for row in rows
            if row.get("event") == "connection_disconnect"
        ]
        self.assertEqual([len(d["reason"]) for d in disconnects], [limit])
        self.assertTrue(disconnects[0]["ready"])

    def test_hostile_error_cycles_do_not_exhaust_day_cap(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        hostile = {"code": -1, "message": "x" * 4_000_000}
        for cycle in range(20):
            with self.assertRaises(ConnectionError):
                collector._handle_notification(
                    {"error": hostile, "id": None},
                    received_ms + cycle,
                )
            collector._log_event(
                "reconnect_backoff",
                {
                    "error": f"ConnectionError: {'x' * 4_000_000}",
                    "backoff_s": 1.0,
                },
                now_ms=received_ms + cycle,
            )
        self.assertLess(collector.sink.day_bytes, 200_000)
        self.assertFalse(collector.stats.get("cap_caught", False))

    def test_pre_write_trade_rejection_preserves_batch_siblings(self) -> None:
        limit = deribit.REQUIRED_SOURCE_MAX_SKEW_MS
        for case in ("skew", "conflict", "chronology"):
            with self.subTest(rejection=case):
                out = Path(tempfile.mkdtemp(
                    prefix=f"collector-integrity-batch-{case}-"
                ))
                collector = deribit.Collector(
                    out,
                    currencies=("btc",),
                    day_cap=10_000_000,
                    total_cap=100_000_000,
                    heartbeat=False,
                )
                self.addCleanup(collector.sink.close)
                received_ms = int(deribit.time.time() * 1000)
                rows = [
                    self.instrument_payload(index, currency="btc")
                    for index in range(3)
                ]
                collector._apply_instrument_snapshot(
                    {"btc": rows},
                    received_ms,
                )
                names = sorted(row["instrument_name"] for row in rows)
                first, target, last = names

                def trade_for(
                    name: str,
                    sequence: int,
                    *,
                    trade_id: str,
                    timestamp: int,
                ) -> dict:
                    return {
                        **self.trade_payload(
                            sequence,
                            trade_id=trade_id,
                            received_ms=received_ms,
                        ),
                        "instrument_name": name,
                        "timestamp": timestamp,
                    }

                if case == "skew":
                    rejected = trade_for(
                        target,
                        100,
                        trade_id="rejected",
                        timestamp=received_ms + limit + 1,
                    )
                else:
                    seeded = trade_for(
                        target,
                        100,
                        trade_id="seeded",
                        timestamp=received_ms - 10,
                    )
                    collector.ingest_trade(seeded, received_ms)
                    rejected = trade_for(
                        target,
                        100 if case == "conflict" else 101,
                        trade_id="rejected",
                        timestamp=(
                            received_ms - 10
                            if case == "conflict"
                            else received_ms - 11
                        ),
                    )

                malformed = {
                    **trade_for(
                        last,
                        7,
                        trade_id="malformed",
                        timestamp=received_ms - 5,
                    ),
                    "amount": 0,
                }
                survivors = [
                    trade_for(
                        first,
                        50,
                        trade_id="survivor-first",
                        timestamp=received_ms - 5,
                    ),
                    trade_for(
                        last,
                        60,
                        trade_id="survivor-last",
                        timestamp=received_ms - 5,
                    ),
                ]
                with self.assertRaises(ConnectionError):
                    collector._route(
                        "trades.option.BTC.100ms",
                        [*survivors, rejected, malformed],
                        received_ms,
                    )
                day_dir = out / f"day{deribit.utc_day(received_ms)}"
                persisted = [
                    deribit.strict_json_loads(line)["trade"]["trade_id"]
                    for line in (day_dir / "trades.jsonl")
                    .read_text(encoding="utf-8")
                    .splitlines()
                ]
                self.assertNotIn("rejected", persisted)
                self.assertIn("survivor-first", persisted)
                self.assertIn("survivor-last", persisted)
                self.assertEqual(collector._last_trade_seq[first], 50)
                self.assertEqual(collector._last_trade_seq[last], 60)
                self.assertEqual(
                    collector._last_trade_seq.get(target),
                    None if case == "skew" else 100,
                )
                events = [
                    deribit.strict_json_loads(line)
                    for line in (day_dir / "events.jsonl")
                    .read_text(encoding="utf-8")
                    .splitlines()
                ]
                rejection_events = {
                    "skew": "trade_source_timestamp_error",
                    "conflict": "trade_sequence_conflict",
                    "chronology": "trade_source_timestamp_regression",
                }
                self.assertEqual(
                    sum(
                        event.get("event") == rejection_events[case]
                        for event in events
                    ),
                    1,
                )
                self.assertEqual(
                    sum(
                        event.get("event") == "trade_schema_error"
                        for event in events
                    ),
                    1,
                )
                self.assertEqual(collector.stats["trade_schema_errors"], 1)

    def test_ticker_schema_error_never_refreshes_freshness(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        channel = f"ticker.{self.instrument_name}.agg2"
        collector.active_channels = {channel}
        received_ms = int(deribit.time.time() * 1000)
        collector._route(
            channel,
            self.ticker_payload(timestamp=received_ms),
            received_ms,
            100.0,
        )
        malformed = self.ticker_payload(timestamp=received_ms + 1)
        del malformed["mark_iv"]
        collector._route(channel, malformed, received_ms + 1, 300.0)
        self.assertEqual(
            collector._required_last_received_monotonic[channel],
            100.0,
        )
        self.assertEqual(collector.stats["ticker_schema_errors"], 1)

    def test_deribit_unexpected_startup_failure_is_controlled_terminal(
        self,
    ) -> None:
        collector, out = self.make_collector()
        stdout = io.StringIO()
        stderr = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                collector,
                "seed_instruments",
                side_effect=ValueError("sensitive startup detail"),
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
            mock.patch.object(deribit.sys, "stderr", stderr),
        ):
            rc = collector.run()

        records = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        stops = [
            record
            for record in records
            if record.get("event") == "session_stop"
        ]
        self.assertEqual(rc, 5)
        self.assertNotIn(deribit.READY_BANNER, stdout.getvalue())
        self.assertIn(deribit.SHUTDOWN_BANNER, stdout.getvalue())
        self.assertIn(
            "DERIBIT_CAPTURE_INTERNAL_ERROR ValueError",
            stderr.getvalue(),
        )
        self.assertNotIn("sensitive startup detail", stderr.getvalue())
        self.assertFalse(
            any(
                str(record.get("event", "")).startswith("connection_")
                for record in records
            )
        )
        self.assertEqual(
            [record["detail"]["reason"] for record in stops],
            ["internal-error"],
        )

    def test_instrument_creation_must_precede_expiry(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        base = self.instrument_payload(0, currency="btc")
        expiration = base["expiration_timestamp"]
        self.assertTrue(
            collector._instrument_row_matches_schema(
                {**base, "creation_timestamp": expiration - 1},
                "btc",
            )
        )
        for label, creation in (
            ("equal", expiration),
            ("after_expiry", expiration + 1),
            ("far_after_expiry", expiration + 86_400_000),
        ):
            with self.subTest(creation=label):
                self.assertFalse(
                    collector._instrument_row_matches_schema(
                        {**base, "creation_timestamp": creation},
                        "btc",
                    )
                )

        baseline_instruments = dict(collector.instruments)
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[
                    {**base, "creation_timestamp": expiration},
                ],
            ),
            self.assertRaises(deribit.InstrumentSnapshotError) as failure,
        ):
            collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(failure.exception.error_type, "SchemaMismatch")
        self.assertEqual(collector.instruments, baseline_instruments)

    def test_committed_inverted_lifetime_is_replay_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-lifetime-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        base = self.instrument_payload(0, currency="btc")
        row = {
            **base,
            "creation_timestamp": base["expiration_timestamp"] + 1,
        }
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        member = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(row, received_ms),
            "received_ts_ms": received_ms,
            "snapshot_id": "inverted-lifetime-snapshot",
            "instrument_fingerprint": builder._instrument_fingerprint(row),
            "instrument": row,
        }
        commit = builder._instrument_snapshot_commit_record(
            "inverted-lifetime-snapshot",
            [member],
            received_ms,
        )
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "instruments.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(record, separators=(",", ":"))
                + "\n"
                for record in (member, commit)
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(stats["metadata_replay_rejected"], 2)
        self.assertEqual(stats["malformed"], 0)

    def test_snapshot_members_must_be_open_and_active(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        base = self.instrument_payload(0, currency="btc")
        self.assertTrue(
            collector._instrument_row_matches_schema(base, "btc")
        )
        for label, overrides in (
            ("closed_inactive", {"state": "closed", "is_active": False}),
            ("open_inactive", {"is_active": False}),
            ("closed_active", {"state": "closed"}),
            ("arbitrary_state", {"state": "halted"}),
            ("state_not_string", {"state": 1}),
            ("active_not_bool", {"is_active": "true"}),
        ):
            with self.subTest(row=label):
                self.assertFalse(
                    collector._instrument_row_matches_schema(
                        {**base, **overrides},
                        "btc",
                    )
                )

    def test_not_open_rows_are_filtered_not_fatal(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(4)
        ]
        closed = {
            **self.instrument_payload(4, currency="btc"),
            "state": "closed",
            "is_active": False,
        }
        inactive = {
            **self.instrument_payload(5, currency="btc"),
            "is_active": False,
        }
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            return_value=[*rows, closed, inactive],
        ):
            fetched = collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(
            [row["instrument_name"] for row in fetched["btc"]],
            [row["instrument_name"] for row in rows],
        )
        self.assertEqual(collector.stats["instrument_rows_not_open"], 2)

    def test_unusable_state_row_is_atomic_snapshot_failure(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        base = self.instrument_payload(0, currency="btc")
        baseline_instruments = dict(collector.instruments)
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[{**base, "state": "halted"}],
            ),
            self.assertRaises(deribit.InstrumentSnapshotError) as failure,
        ):
            collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(failure.exception.error_type, "SchemaMismatch")
        self.assertEqual(collector.instruments, baseline_instruments)
        self.assertEqual(collector.stats["instrument_rows_not_open"], 0)

    def test_committed_inactive_member_is_replay_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-inactive-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        row = {
            **self.instrument_payload(0, currency="btc"),
            "state": "closed",
            "is_active": False,
        }
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        member = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(row, received_ms),
            "received_ts_ms": received_ms,
            "snapshot_id": "inactive-snapshot",
            "instrument_fingerprint": builder._instrument_fingerprint(row),
            "instrument": row,
        }
        commit = builder._instrument_snapshot_commit_record(
            "inactive-snapshot",
            [member],
            received_ms,
        )
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "instruments.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(record, separators=(",", ":"))
                + "\n"
                for record in (member, commit)
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(stats["metadata_replay_rejected"], 2)

    def test_same_sequence_distinct_trade_is_connection_terminal(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = 1_800_000_000_500
        first = self.trade_payload(100, trade_id="A")
        collector.ingest_trade(first, received_ms)
        collector.ingest_trade(first, received_ms + 1)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            100,
        )
        conflict = self.trade_payload(100, trade_id="B")
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(conflict, received_ms + 2)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            100,
        )
        self.assertNotIn(
            collector._option_trade_dedup_id(conflict),
            collector.deduper.seen,
        )
        collector.ingest_trade(
            self.trade_payload(101, trade_id="C"),
            received_ms + 3,
        )
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            101,
        )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        trades = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "trades.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        self.assertEqual(
            [row["trade"]["trade_id"] for row in trades],
            ["A", "C"],
        )
        self.assertEqual(collector.stats["dup_trades"], 1)
        self.assertEqual(collector.stats["trade_sequence_conflicts"], 1)
        events = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        conflicts = [
            event["detail"]
            for event in events
            if event.get("event") == "trade_sequence_conflict"
        ]
        self.assertEqual(len(conflicts), 1)
        self.assertEqual(
            conflicts[0],
            {
                "instrument_name": self.instrument_name,
                "observed": 100,
                "last_persisted": 100,
                "reason": "repeat",
            },
        )
        regression = self.trade_payload(100, trade_id="D")
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(regression, received_ms + 4)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            101,
        )
        self.assertNotIn(
            collector._option_trade_dedup_id(regression),
            collector.deduper.seen,
        )
        self.assertEqual(collector.stats["trade_sequence_conflicts"], 2)
        regressions = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        self.assertEqual(
            [
                event["detail"]["reason"]
                for event in regressions
                if event.get("event") == "trade_sequence_conflict"
            ],
            ["repeat", "regression"],
        )
        self.assertEqual(
            [
                row["trade"]["trade_id"]
                for row in (
                    deribit.strict_json_loads(line)
                    for line in (day_dir / "trades.jsonl")
                    .read_text(encoding="utf-8")
                    .splitlines()
                )
            ],
            ["A", "C"],
        )

    def test_replayed_same_sequence_siblings_are_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-seq-conflict-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder, _ = self.make_collector()
        rows = []
        for trade_id in ("A", "B"):
            trade = {
                **self.trade_payload(100, trade_id=trade_id),
                "timestamp": received_ms - 1,
            }
            rows.append({
                "record_type": "option_trade",
                "dedup_id": builder._option_trade_dedup_id(trade),
                "received_ts_ms": received_ms,
                "trade": trade,
            })
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "trades.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(row, separators=(",", ":"))
                + "\n"
                for row in rows
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(
            restarted._last_trade_seq,
            {self.instrument_name: 100},
        )
        self.assertEqual(stats["trade_sequence_replay_rejected"], 1)
        self.assertEqual(stats["malformed"], 0)

    def test_replayed_sequence_regression_is_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-seq-regression-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder, _ = self.make_collector()
        rows = []
        for sequence, trade_id in ((100, "A"), (101, "B"), (100, "C")):
            trade = {
                **self.trade_payload(sequence, trade_id=trade_id),
                "timestamp": received_ms - 1,
            }
            rows.append({
                "record_type": "option_trade",
                "dedup_id": builder._option_trade_dedup_id(trade),
                "received_ts_ms": received_ms,
                "trade": trade,
            })
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "trades.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(row, separators=(",", ":"))
                + "\n"
                for row in rows
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(
            restarted._last_trade_seq,
            {self.instrument_name: 101},
        )
        self.assertEqual(stats["trade_sequence_replay_rejected"], 1)
        self.assertEqual(stats["malformed"], 0)

    def test_replay_never_downgrades_checkpoint_highwater(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        trade = {
            **self.trade_payload(100, trade_id="late"),
            "timestamp": received_ms - 1,
        }
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        path = day_dir / "trades.jsonl"
        path.write_text(
            deribit.strict_json_dumps(
                {
                    "record_type": "option_trade",
                    "dedup_id": collector._option_trade_dedup_id(trade),
                    "received_ts_ms": received_ms,
                    "trade": trade,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        collector._last_trade_seq = {self.instrument_name: 101}
        stats = {
            "malformed": 0,
            "cursor_replay_bytes": 0,
            "trade_sequence_replay_rejected": 0,
        }
        collector._read_trade_sequences(
            path,
            stats,
            0,
            "cursor_replay_bytes",
        )
        self.assertEqual(
            collector._last_trade_seq,
            {self.instrument_name: 101},
        )
        self.assertEqual(stats["trade_sequence_replay_rejected"], 1)

    def test_trade_source_timestamp_skew_is_enforced(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        limit = deribit.REQUIRED_SOURCE_MAX_SKEW_MS
        received_ms = int(deribit.time.time() * 1000)
        for sequence, offset in ((100, -limit), (101, limit)):
            collector.ingest_trade(
                {
                    **self.trade_payload(sequence),
                    "timestamp": received_ms + offset,
                },
                received_ms,
            )
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            101,
        )
        for label, offset in (
            ("future", limit + 1),
            ("stale", -limit - 1),
        ):
            with self.subTest(skew=label):
                rejected = {
                    **self.trade_payload(102, trade_id=f"bad-{label}"),
                    "timestamp": received_ms + offset,
                }
                with self.assertRaises(ConnectionError):
                    collector.ingest_trade(rejected, received_ms)
                self.assertNotIn(
                    collector._option_trade_dedup_id(rejected),
                    collector.deduper.seen,
                )
                self.assertEqual(
                    collector._last_trade_seq[self.instrument_name],
                    101,
                )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        self.assertEqual(
            [
                deribit.strict_json_loads(line)["trade"]["trade_seq"]
                for line in (day_dir / "trades.jsonl")
                .read_text(encoding="utf-8")
                .splitlines()
            ],
            [100, 101],
        )
        audits = [
            deribit.strict_json_loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        self.assertEqual(
            [
                event["detail"]["reason"]
                for event in audits
                if event.get("event") == "trade_source_timestamp_error"
            ],
            ["source_timestamp_future", "source_timestamp_stale"],
        )
        self.assertEqual(collector.stats["source_timestamp_errors"], 2)

    def test_trade_source_time_must_not_regress_with_sequence(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        base_ts = received_ms - 1_000
        collector.ingest_trade(
            {**self.trade_payload(100), "timestamp": base_ts},
            received_ms,
        )
        collector.ingest_trade(
            {**self.trade_payload(101), "timestamp": base_ts},
            received_ms + 1,
        )
        collector.ingest_trade(
            {**self.trade_payload(102), "timestamp": base_ts + 500},
            received_ms + 2,
        )
        self.assertEqual(
            collector._last_trade_source_ts[self.instrument_name],
            base_ts + 500,
        )
        inverted = {
            **self.trade_payload(103),
            "timestamp": base_ts + 499,
        }
        with self.assertRaises(ConnectionError):
            collector.ingest_trade(inverted, received_ms + 3)
        self.assertEqual(
            collector._last_trade_seq[self.instrument_name],
            102,
        )
        self.assertEqual(
            collector._last_trade_source_ts[self.instrument_name],
            base_ts + 500,
        )
        self.assertNotIn(
            collector._option_trade_dedup_id(inverted),
            collector.deduper.seen,
        )
        self.assertEqual(
            collector.stats["trade_chronology_violations"],
            1,
        )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        self.assertEqual(
            [
                deribit.strict_json_loads(line)["trade"]["trade_seq"]
                for line in (day_dir / "trades.jsonl")
                .read_text(encoding="utf-8")
                .splitlines()
            ],
            [100, 101, 102],
        )
        audits = [
            deribit.strict_json_loads(line)["detail"]
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
            if deribit.strict_json_loads(line).get("event")
            == "trade_source_timestamp_regression"
        ]
        self.assertEqual(
            audits,
            [
                {
                    "instrument_name": self.instrument_name,
                    "observed": base_ts + 499,
                    "last_persisted": base_ts + 500,
                    "sequence": 103,
                }
            ],
        )

    def test_chronology_survives_checkpoint_restart(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-chronology-restart-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        base_ts = received_ms - 1_000
        writer = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        writer.instruments = {
            self.instrument_name: {
                "instrument_name": self.instrument_name,
                "kind": "option",
            }
        }
        writer.ingest_trade(
            {**self.trade_payload(100), "timestamp": base_ts},
            received_ms,
        )
        writer._save_sequence_state()
        writer.sink.close()
        checkpoint = deribit.strict_json_loads(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            checkpoint["schema_version"],
            deribit.TRADE_SEQUENCE_STATE_SCHEMA_VERSION,
        )
        self.assertEqual(
            checkpoint["trade_source_timestamps"],
            {self.instrument_name: base_ts},
        )

        restarted, stats = self.restart_replay(out)
        self.assertTrue(stats["sequence_checkpoint_loaded"])
        self.assertEqual(
            restarted._last_trade_source_ts,
            {self.instrument_name: base_ts},
        )
        restarted.instruments = dict(writer.instruments)
        with self.assertRaises(ConnectionError):
            restarted.ingest_trade(
                {**self.trade_payload(101), "timestamp": base_ts - 1},
                received_ms + 5,
            )
        self.assertEqual(
            restarted._last_trade_seq,
            {self.instrument_name: 100},
        )

    def test_replayed_chronology_inversion_is_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-chronology-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        base_ts = received_ms - 1_000
        builder, _ = self.make_collector()
        rows = []
        for sequence, offset in ((100, 0), (101, -1)):
            trade = {
                **self.trade_payload(sequence),
                "timestamp": base_ts + offset,
            }
            rows.append({
                "record_type": "option_trade",
                "dedup_id": builder._option_trade_dedup_id(trade),
                "received_ts_ms": received_ms,
                "trade": trade,
            })
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "trades.jsonl").write_text(
            "".join(
                deribit.strict_json_dumps(row, separators=(",", ":"))
                + "\n"
                for row in rows
            ),
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(
            restarted._last_trade_seq,
            {self.instrument_name: 100},
        )
        self.assertEqual(
            restarted._last_trade_source_ts,
            {self.instrument_name: base_ts},
        )
        self.assertEqual(stats["trade_sequence_replay_rejected"], 1)

    def test_replayed_trade_skew_violation_is_rejected(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-trade-skew-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder, _ = self.make_collector()
        trade = {
            **self.trade_payload(100),
            "timestamp": received_ms
            + deribit.REQUIRED_SOURCE_MAX_SKEW_MS
            + 1,
        }
        builder.sink.close()
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        day_dir.mkdir(parents=True, exist_ok=True)
        (day_dir / "trades.jsonl").write_text(
            deribit.strict_json_dumps(
                {
                    "record_type": "option_trade",
                    "dedup_id": builder._option_trade_dedup_id(trade),
                    "received_ts_ms": received_ms,
                    "trade": trade,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted._last_trade_seq, {})
        self.assertEqual(stats["trade_sequence_replay_rejected"], 1)
        self.assertEqual(stats["trade_ids"], 0)

    def test_deribit_unexpected_post_ready_failure_is_controlled_terminal(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_minutes = 1.0
        stdout = io.StringIO()
        stderr = io.StringIO()

        class ReadyWebSocket:
            def close(self) -> None:
                return None

        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(collector, "seed_instruments"),
            mock.patch.object(
                collector,
                "_connect_websocket_interruptibly",
                return_value=ReadyWebSocket(),
            ),
            mock.patch.object(collector, "subscribe_all", return_value=True),
            mock.patch.object(
                collector,
                "await_required_coverage",
                return_value=True,
            ),
            mock.patch.object(
                collector,
                "_loop",
                side_effect=AssertionError("sensitive loop detail"),
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
            mock.patch.object(deribit.sys, "stderr", stderr),
        ):
            rc = collector.run()

        records = [
            deribit.strict_json_loads(line)
            for path in out.glob("day*/*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        disconnects = [
            record
            for record in records
            if record.get("event") == "connection_disconnect"
        ]
        stops = [
            record
            for record in records
            if record.get("event") == "session_stop"
        ]
        self.assertEqual(rc, 5)
        self.assertIn(deribit.READY_BANNER, stdout.getvalue())
        self.assertIn(deribit.SHUTDOWN_BANNER, stdout.getvalue())
        self.assertIn(
            "DERIBIT_CAPTURE_INTERNAL_ERROR AssertionError",
            stderr.getvalue(),
        )
        self.assertNotIn("sensitive loop detail", stderr.getvalue())
        self.assertEqual(collector.stats["reconnects"], 0)
        self.assertFalse(
            any(record.get("event") == "reconnect_backoff" for record in records)
        )
        self.assertEqual(
            [
                (
                    record["detail"]["reason"],
                    record["detail"]["ready"],
                )
                for record in disconnects
            ],
            [("internal-error:AssertionError", True)],
        )
        self.assertEqual(
            [record["detail"]["reason"] for record in stops],
            ["internal-error"],
        )

    def test_deribit_internal_stop_audit_failure_still_finalizes(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        collector.stop = True
        stdout = io.StringIO()
        stderr = io.StringIO()
        original_close = collector.sink.close
        closed = 0

        def fail_stop_audit(event: str, *_args, **_kwargs) -> None:
            self.assertEqual(event, "session_stop")
            raise TypeError("sensitive audit detail")

        def close_sink() -> None:
            nonlocal closed
            closed += 1
            original_close()

        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                collector,
                "_log_event",
                side_effect=fail_stop_audit,
            ),
            mock.patch.object(collector.sink, "close", side_effect=close_sink),
            mock.patch.object(deribit.sys, "stdout", stdout),
            mock.patch.object(deribit.sys, "stderr", stderr),
        ):
            rc = collector.run()

        self.assertEqual(rc, 5)
        self.assertEqual(closed, 1)
        self.assertIn(deribit.SHUTDOWN_BANNER, stdout.getvalue())
        self.assertIn(
            "DERIBIT_CAPTURE_INTERNAL_ERROR TypeError",
            stderr.getvalue(),
        )
        self.assertNotIn("sensitive audit detail", stderr.getvalue())

    def test_deribit_cap_before_readiness_exits_nonzero(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-cap-exit-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=1,
            total_cap=100_000_000,
            heartbeat=False,
            max_minutes=0.001,
        )
        output = io.StringIO()
        with mock.patch.object(deribit.sys, "stdout", output):
            rc = collector.run()
        self.assertEqual(rc, 3)
        self.assertIn("DERIBIT_CAPTURE_CAP_EXCEEDED", output.getvalue())
        self.assertIn(deribit.SHUTDOWN_BANNER, output.getvalue())
        self.assertNotIn(deribit.READY_BANNER, output.getvalue())

    def test_deribit_cap_first_discovered_at_session_stop_exits_nonzero(
        self,
    ) -> None:
        collector, out = self.make_collector()
        now_ms = int(deribit.time.time() * 1000)
        collector.sink.write(
            "events",
            {
                "record_type": "capture_event",
                "event": "preexisting",
                "received_ts": now_ms,
            },
            now_ms,
        )
        collector.sink.day_cap = collector.sink.day_bytes
        collector.stop = True
        output = io.StringIO()

        with mock.patch.object(deribit.sys, "stdout", output):
            rc = collector.run()

        records = [
            json.loads(line)
            for path in out.glob("day*/*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(rc, 3)
        self.assertTrue(collector.stats["cap_caught"])
        self.assertNotIn(deribit.READY_BANNER, output.getvalue())
        self.assertIn(
            "DERIBIT_CAPTURE_CAP_EXCEEDED during session stop",
            output.getvalue(),
        )
        self.assertIn(deribit.SHUTDOWN_BANNER, output.getvalue())
        self.assertFalse(
            any(record.get("event") == "session_stop" for record in records)
        )

    def test_deribit_session_stop_io_failure_remains_sink_error(self) -> None:
        collector, _ = self.make_collector()
        collector.stop = True
        output = io.StringIO()

        def fail_write(*_args, **_kwargs):
            error = OSError("injected session-stop write failure")
            collector.sink._latch_io_error("write", error)
            raise error

        with (
            mock.patch.object(deribit.sys, "stdout", output),
            mock.patch.object(collector.sink, "write", side_effect=fail_write),
        ):
            rc = collector.run()

        self.assertEqual(rc, 5)
        self.assertFalse(collector.stats.get("cap_caught", False))
        self.assertNotIn("DERIBIT_CAPTURE_CAP_EXCEEDED", output.getvalue())
        self.assertIn(deribit.SHUTDOWN_BANNER, output.getvalue())

    def test_deribit_cap_after_readiness_exits_nonzero(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        collector.max_minutes = 1.0
        output = io.StringIO()

        class ReadyWebSocket:
            def close(self) -> None:
                return None

        with (
            mock.patch.object(deribit.sys, "stdout", output),
            mock.patch.object(collector, "seed_instruments"),
            mock.patch.object(collector, "subscribe_all", return_value=True),
            mock.patch.object(
                collector,
                "await_required_coverage",
                return_value=True,
            ),
            mock.patch.object(
                collector,
                "_loop",
                side_effect=deribit.DayCapError("post-ready cap"),
            ),
            mock.patch.object(
                deribit,
                "WsClient",
                return_value=ReadyWebSocket(),
            ),
        ):
            rc = collector.run()
        self.assertEqual(rc, 3)
        self.assertIn(deribit.READY_BANNER, output.getvalue())
        self.assertIn("DERIBIT_CAPTURE_CAP_EXCEEDED", output.getvalue())
        self.assertIn(deribit.SHUTDOWN_BANNER, output.getvalue())

    def test_instrument_seed_is_atomic_across_required_currencies(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-seed-"))
        collector = deribit.Collector(
            out,
            currencies=("btc", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        prior = {
            "BTC-OLD": {"instrument_name": "BTC-OLD", "creation_timestamp": 1},
            "ETH-OLD": {"instrument_name": "ETH-OLD", "creation_timestamp": 1},
        }
        collector.instruments = dict(prior)
        collector._persisted_instrument_fingerprints = {
            name: collector._instrument_fingerprint(row)
            for name, row in prior.items()
        }
        collector._instrument_persist_day = "20260903"

        def fetch(currency: str, deadline=None):
            del deadline
            if currency == "eth":
                raise ConnectionError("simulated ETH REST failure")
            return [
                self.instrument_payload(
                    0,
                    name="BTC-NEW",
                )
            ]

        with (
            mock.patch.object(deribit, "fetch_instruments", side_effect=fetch),
            self.assertRaises(ConnectionError),
        ):
            collector.seed_instruments()
        self.assertEqual(collector.instruments, prior)
        self.assertEqual(
            set(collector._persisted_instrument_fingerprints),
            set(prior),
        )
        self.assertEqual(collector._instrument_persist_day, "20260903")

    def test_instrument_seed_commits_all_currencies_together(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-seed-"))
        collector = deribit.Collector(
            out,
            currencies=("btc", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)

        def fetch(currency: str, deadline=None):
            del deadline
            return [
                self.instrument_payload(
                    index,
                    currency=currency,
                )
                for index in range(10)
            ]

        with mock.patch.object(deribit, "fetch_instruments", side_effect=fetch):
            collector.seed_instruments()
        expected = {
            self.instrument_payload(
                index,
                currency=currency,
            )["instrument_name"]
            for currency in ("btc", "eth")
            for index in range(10)
        }
        self.assertEqual(set(collector.instruments), expected)
        self.assertEqual(
            set(collector._persisted_instrument_fingerprints),
            expected,
        )

    def test_instrument_snapshot_batch_uses_one_fsync(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-instrument-batch-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        collector.sink.write(
            "events",
            {
                "record_type": "capture_event",
                "dedup_id": "open-day",
                "event": "probe",
                "received_ts_ms": received_ms,
            },
            now_ms=received_ms,
        )
        day_dir = out / f"day{deribit.utc_day(received_ms)}"
        instrument_inode = (day_dir / "instruments.jsonl").stat().st_ino
        instrument_fsyncs = 0

        def count_fsync(fd: int) -> None:
            nonlocal instrument_fsyncs
            if deribit.os.fstat(fd).st_ino == instrument_inode:
                instrument_fsyncs += 1

        rows = [
            {
                "instrument_name": f"BTC-BATCH-{index:04d}",
                "creation_timestamp": index,
            }
            for index in range(2_000)
        ]
        with mock.patch.object(
            deribit.os,
            "fsync",
            side_effect=count_fsync,
        ):
            collector._apply_instrument_snapshot(
                {"btc": rows},
                received_ms,
            )
        self.assertEqual(instrument_fsyncs, 1)
        self.assertEqual(collector.stats["instrument_records"], 2_000)
        self.assertEqual(
            len(collector._persisted_instrument_fingerprints),
            2_000,
        )
        self.assertEqual(
            len(
                (day_dir / "instruments.jsonl")
                .read_text(encoding="utf-8")
                .splitlines()
            ),
            2_001,
        )

    def test_instrument_snapshot_batch_cap_failure_writes_zero_rows(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-instrument-cap-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=600,
            total_cap=100_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            {
                "instrument_name": f"BTC-CAP-{index:02d}",
                "creation_timestamp": index,
            }
            for index in range(20)
        ]
        with self.assertRaises(deribit.DayCapError):
            collector._apply_instrument_snapshot(
                {"btc": rows},
                received_ms,
            )
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        self.assertEqual(instrument_path.read_bytes(), b"")
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            {},
        )
        self.assertEqual(collector.stats["instrument_records"], 0)

    def test_instrument_snapshot_rejects_remote_row_flood(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        rows = [
            {
                "instrument_name": f"BTC-FLOOD-{index}",
                "creation_timestamp": index,
            }
            for index in range(3)
        ]
        with (
            mock.patch.object(
                deribit,
                "INSTRUMENT_SNAPSHOT_MAX_ROWS_PER_CURRENCY",
                2,
                create=True,
            ),
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=rows,
            ),
            self.assertRaises(deribit.InstrumentSnapshotError),
        ):
            collector._fetch_instrument_snapshot(deadline=None)

    def test_instrument_snapshot_versions_tombstones_and_replays(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-instrument-history-"))
        received_ms = int(deribit.time.time() * 1000)
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        first_a = self.instrument_payload(0, currency="btc")
        first_b = self.instrument_payload(1, currency="btc")
        second_a = {**first_a, "tick_size": 0.0002}
        collector._apply_instrument_snapshot(
            {"btc": [first_a, first_b]},
            received_ms,
        )
        collector._apply_instrument_snapshot(
            {"btc": [dict(first_a), dict(first_b)]},
            received_ms + 1,
        )
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        self.assertEqual(len(instrument_path.read_text().splitlines()), 3)

        collector._apply_instrument_snapshot(
            {"btc": [second_a]},
            received_ms + 2,
        )
        collector._apply_instrument_snapshot(
            {"btc": [dict(second_a)]},
            received_ms + 3,
        )
        records = [
            json.loads(line)
            for line in instrument_path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            [record["record_type"] for record in records],
            [
                "instrument",
                "instrument",
                "instrument_snapshot_commit",
                "instrument",
                "instrument_removed",
                "instrument_snapshot_commit",
            ],
        )
        self.assertEqual(records[3]["instrument"], second_a)
        self.assertEqual(
            records[4]["instrument_name"],
            first_b["instrument_name"],
        )
        self.assertEqual(
            collector.instruments,
            {first_a["instrument_name"]: second_a},
        )
        self.assertEqual(collector.stats["instrument_records"], 4)

        collector.sink.close()
        replay = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(replay.sink.close)
        replay.rehydrate_state()
        self.assertEqual(replay.instruments, collector.instruments)

    def test_semantic_metadata_poison_cannot_establish_chain_baseline(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-metadata-poison-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        complete = self.instrument_payload(0, currency="btc")
        incomplete = {
            "instrument_name": complete["instrument_name"],
            "kind": "option",
        }
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        fingerprint = builder._instrument_fingerprint(incomplete)
        poison = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(
                incomplete,
                received_ms,
            ),
            "received_ts_ms": received_ms,
            "snapshot_id": "forged-complete-snapshot",
            "instrument_fingerprint": fingerprint,
            "instrument": incomplete,
        }
        builder.sink.close()
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        instrument_path.parent.mkdir(parents=True, exist_ok=True)
        instrument_path.write_text(
            json.dumps(poison, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state()
        self.assertEqual(stats["metadata_replay_rejected"], 1)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[complete],
            ) as fetch,
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as failure,
        ):
            restarted._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(
            failure.exception.error_type,
            "InitialSnapshotImplausible",
        )
        fetch.assert_called_once_with("btc", deadline=None)

    def test_metadata_replay_requires_current_record_identity(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-metadata-identity-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )

        def upsert(
            row: dict,
            timestamp: int,
            snapshot_id: str = "snapshot",
        ) -> dict:
            fingerprint = builder._instrument_fingerprint(row)
            return {
                "record_type": "instrument",
                "dedup_id": builder._instrument_dedup_id(
                    row,
                    timestamp,
                ),
                "received_ts_ms": timestamp,
                "snapshot_id": snapshot_id,
                "instrument_fingerprint": fingerprint,
                "instrument": row,
            }

        def removal(
            row: dict,
            timestamp: int,
            *,
            fingerprint: str | None = None,
            snapshot_id: str = "snapshot-removal",
            reason: str = "absent_from_authoritative_snapshot",
        ) -> dict:
            prior_fingerprint = (
                builder._instrument_fingerprint(row)
                if fingerprint is None
                else fingerprint
            )
            instrument_name = row["instrument_name"]
            return {
                "record_type": "instrument_removed",
                "dedup_id": deribit.sha256_id(
                    f"instr|removed|{instrument_name}|"
                    f"{prior_fingerprint}|{snapshot_id}"
                ),
                "received_ts_ms": timestamp,
                "snapshot_id": snapshot_id,
                "instrument_name": instrument_name,
                "prior_instrument_fingerprint": prior_fingerprint,
                "reason": reason,
            }

        valid = self.instrument_payload(10, currency="btc")
        foreign = self.instrument_payload(11, currency="eth")
        incomplete = {
            "instrument_name": self.instrument_payload(
                12,
                currency="btc",
            )["instrument_name"],
            "kind": "option",
        }
        invalid_records = []
        invalid_records.append(upsert(incomplete, received_ms))
        bad_fingerprint = upsert(valid, received_ms + 1)
        bad_fingerprint["instrument_fingerprint"] = "0" * 64
        invalid_records.append(bad_fingerprint)
        bad_dedup = upsert(valid, received_ms + 2)
        bad_dedup["dedup_id"] = "1" * 64
        invalid_records.append(bad_dedup)
        invalid_records.append(upsert(valid, 0))
        invalid_records.append(
            upsert(valid, received_ms + 4, snapshot_id="")
        )
        invalid_records.append(
            removal(
                valid,
                received_ms + 5,
                reason="untrusted-reason",
            )
        )
        invalid_records.append(
            removal(
                valid,
                received_ms + 6,
                fingerprint="not-a-sha256",
            )
        )
        bad_removal_dedup = removal(valid, received_ms + 7)
        bad_removal_dedup["dedup_id"] = "2" * 64
        invalid_records.append(bad_removal_dedup)

        valid_upsert = upsert(valid, received_ms + 8)
        valid_removal = removal(valid, received_ms + 9)
        records = [
            *invalid_records,
            upsert(foreign, received_ms + 10),
            valid_upsert,
            valid_removal,
        ]
        builder.sink.close()
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        instrument_path.parent.mkdir(parents=True, exist_ok=True)
        instrument_path.write_text(
            "".join(
                json.dumps(record, separators=(",", ":")) + "\n"
                for record in records
            ),
            encoding="utf-8",
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state()
        self.assertEqual(
            stats["metadata_replay_rejected"],
            len(invalid_records),
        )
        self.assertEqual(stats["malformed"], 0)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )

    def make_snapshot_writer(self, out: Path, currencies=("btc",)):
        writer = deribit.Collector(
            out,
            currencies=currencies,
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        return writer

    def instrument_lines(self, out: Path, received_ms: int) -> list[str]:
        return (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        ).read_text(encoding="utf-8").splitlines()

    def rewrite_instrument_lines(
        self,
        out: Path,
        received_ms: int,
        lines: list[str],
    ) -> None:
        (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        ).write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )

    def restart_replay(self, out: Path, currencies=("btc",)):
        restarted = deribit.Collector(
            out,
            currencies=currencies,
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        return restarted, restarted.rehydrate_state()

    def test_torn_instrument_snapshot_prefix_never_establishes_baseline(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-torn-snapshot-"))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(4)
        ]
        writer = self.make_snapshot_writer(out)
        writer._apply_instrument_snapshot({"btc": rows}, received_ms)
        writer.sink.close()
        lines = self.instrument_lines(out, received_ms)
        self.assertEqual(
            json.loads(lines[-1])["record_type"],
            "instrument_snapshot_commit",
        )
        self.rewrite_instrument_lines(out, received_ms, lines[:1])

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(restarted._persisted_instrument_fingerprints, {})
        self.assertEqual(stats["instrument_ids"], 0)
        self.assertEqual(stats["metadata_replay_rejected"], 0)
        self.assertEqual(stats["metadata_uncommitted_snapshots"], 1)
        self.assertEqual(stats["metadata_uncommitted_records"], 1)
        self.assertEqual(stats["malformed"], 0)
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[rows[0]],
            ) as fetch,
            self.assertRaises(deribit.InstrumentSnapshotError) as failure,
        ):
            restarted._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(
            failure.exception.error_type,
            "InitialSnapshotImplausible",
        )
        fetch.assert_called_once_with("btc", deadline=None)

    def test_committed_instrument_snapshot_restores_exact_baseline(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-committed-"))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(4)
        ]
        writer = self.make_snapshot_writer(out)
        writer._apply_instrument_snapshot({"btc": rows}, received_ms)
        expected_instruments = dict(writer.instruments)
        expected_fingerprints = dict(
            writer._persisted_instrument_fingerprints
        )
        writer.sink.close()

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, expected_instruments)
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            expected_fingerprints,
        )
        self.assertEqual(stats["instrument_ids"], len(rows))
        self.assertEqual(stats["metadata_replay_rejected"], 0)
        self.assertEqual(stats["metadata_uncommitted_snapshots"], 0)
        self.assertEqual(stats["metadata_uncommitted_records"], 0)

    def test_torn_later_snapshot_cannot_mutate_committed_baseline(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-torn-later-"))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(4)
        ]
        writer = self.make_snapshot_writer(out)
        writer._apply_instrument_snapshot({"btc": rows}, received_ms)
        baseline_instruments = dict(writer.instruments)
        baseline_fingerprints = dict(
            writer._persisted_instrument_fingerprints
        )
        committed_lines = len(self.instrument_lines(out, received_ms))
        updated = {**rows[0], "tick_size": rows[0]["tick_size"] / 2}
        writer._apply_instrument_snapshot(
            {"btc": [updated, *rows[2:]]},
            received_ms + 1,
        )
        writer.sink.close()
        lines = self.instrument_lines(out, received_ms)
        self.assertEqual(
            json.loads(lines[-1])["record_type"],
            "instrument_snapshot_commit",
        )
        self.assertGreater(len(lines) - 1, committed_lines)
        self.rewrite_instrument_lines(out, received_ms, lines[:-1])

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, baseline_instruments)
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            baseline_fingerprints,
        )
        self.assertEqual(stats["metadata_replay_rejected"], 0)
        self.assertEqual(stats["metadata_uncommitted_snapshots"], 1)
        self.assertEqual(
            stats["metadata_uncommitted_records"],
            len(lines) - 1 - committed_lines,
        )

    def test_instrument_snapshot_commit_requires_exact_transaction(
        self,
    ) -> None:
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(3)
        ]

        def tamper_snapshot_id(commit: dict, members: list[dict]) -> list[dict]:
            del members
            return [{**commit, "snapshot_id": "forged-snapshot"}]

        def tamper_member_count(commit: dict, members: list[dict]) -> list[dict]:
            return [{**commit, "member_count": len(members) + 1}]

        def tamper_digest(commit: dict, members: list[dict]) -> list[dict]:
            del members
            return [{**commit, "member_digest": "0" * 64}]

        def tamper_dedup(commit: dict, members: list[dict]) -> list[dict]:
            del members
            return [{**commit, "dedup_id": "1" * 64}]

        def tamper_currencies(commit: dict, members: list[dict]) -> list[dict]:
            del members
            return [{**commit, "currencies": ["doge"]}]

        def tamper_received(commit: dict, members: list[dict]) -> list[dict]:
            del members
            return [{**commit, "received_ts_ms": 0}]

        def duplicate_member(commit: dict, members: list[dict]) -> list[dict]:
            return [members[0], {**commit, "member_count": len(members) + 1}]

        tampers = {
            "snapshot_id": tamper_snapshot_id,
            "member_count": tamper_member_count,
            "member_digest": tamper_digest,
            "dedup_id": tamper_dedup,
            "currencies": tamper_currencies,
            "received_ts_ms": tamper_received,
            "duplicate_member": duplicate_member,
        }
        for name, tamper in tampers.items():
            with self.subTest(tamper=name):
                out = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-commit-identity-"
                ))
                writer = self.make_snapshot_writer(out)
                writer._apply_instrument_snapshot({"btc": rows}, received_ms)
                writer.sink.close()
                lines = self.instrument_lines(out, received_ms)
                members = [json.loads(line) for line in lines[:-1]]
                commit = json.loads(lines[-1])
                self.rewrite_instrument_lines(
                    out,
                    received_ms,
                    [
                        *lines[:-1],
                        *(
                            json.dumps(record, separators=(",", ":"))
                            for record in tamper(commit, members)
                        ),
                    ],
                )
                restarted, stats = self.restart_replay(out)
                self.assertEqual(restarted.instruments, {})
                self.assertEqual(
                    restarted._persisted_instrument_fingerprints,
                    {},
                )
                self.assertEqual(stats["metadata_replay_rejected"], 1)
                self.assertEqual(stats["malformed"], 0)

    def test_commit_before_members_never_publishes_snapshot(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-commit-order-"))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(3)
        ]
        writer = self.make_snapshot_writer(out)
        writer._apply_instrument_snapshot({"btc": rows}, received_ms)
        writer.sink.close()
        lines = self.instrument_lines(out, received_ms)
        self.rewrite_instrument_lines(
            out,
            received_ms,
            [lines[-1], *lines[:-1]],
        )

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(restarted._persisted_instrument_fingerprints, {})
        self.assertEqual(stats["metadata_replay_rejected"], 1)
        self.assertEqual(stats["metadata_uncommitted_snapshots"], 2)
        self.assertEqual(stats["metadata_uncommitted_records"], len(rows))

    def test_failed_snapshot_chunk_leaves_no_trusted_baseline(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-chunk-fail-"))
        received_ms = int(deribit.time.time() * 1000)
        rows = [
            self.instrument_payload(index, currency="btc")
            for index in range(6)
        ]
        writer = self.make_snapshot_writer(out)
        writer._apply_instrument_snapshot({"btc": rows[:2]}, received_ms)
        baseline_instruments = dict(writer.instruments)
        baseline_fingerprints = dict(
            writer._persisted_instrument_fingerprints
        )
        handle = writer.sink._handles["instruments"]
        real_write = handle.write
        writes = {"n": 0}

        def fail_second_chunk(payload):
            writes["n"] += 1
            if writes["n"] > 1:
                raise OSError("injected chunk failure")
            return real_write(payload)

        with (
            mock.patch.object(deribit, "SINK_BATCH_CHUNK_BYTES", 256),
            mock.patch.object(handle, "write", side_effect=fail_second_chunk),
            self.assertRaises(OSError),
        ):
            writer._apply_instrument_snapshot(
                {"btc": rows},
                received_ms + 1,
            )
        self.assertEqual(writer.instruments, baseline_instruments)
        self.assertEqual(
            writer._persisted_instrument_fingerprints,
            baseline_fingerprints,
        )
        writer.sink.close()
        records = [
            json.loads(line)
            for line in self.instrument_lines(out, received_ms)
        ]
        commits = [
            record
            for record in records
            if record["record_type"] == "instrument_snapshot_commit"
        ]
        self.assertEqual(len(commits), 1)
        self.assertEqual(commits[0]["member_count"], 2)

        restarted, stats = self.restart_replay(out)
        self.assertEqual(restarted.instruments, baseline_instruments)
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            baseline_fingerprints,
        )
        self.assertEqual(stats["metadata_replay_rejected"], 0)
        self.assertGreaterEqual(stats["metadata_uncommitted_snapshots"], 1)

    def test_instrument_snapshot_rejects_name_row_identity_mismatch(
        self,
    ) -> None:
        mismatches = {
            "strike": (
                0,
                {"strike": 1.0},
            ),
            "call-name-put-row": (
                0,
                {"option_type": "put"},
            ),
            "put-name-call-row": (
                1,
                {"option_type": "call"},
            ),
        }
        for label, (poison_index, change) in mismatches.items():
            with self.subTest(label=label):
                out = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-name-identity-"
                ))
                collector = deribit.Collector(
                    out,
                    currencies=("btc",),
                    day_cap=10_000_000,
                    total_cap=100_000_000,
                    heartbeat=False,
                )
                self.addCleanup(collector.sink.close)
                rows = [
                    self.instrument_payload(index, currency="btc")
                    for index in range(10)
                ]
                rows[poison_index] = {
                    **rows[poison_index],
                    **change,
                }
                with (
                    mock.patch.object(
                        deribit,
                        "fetch_instruments",
                        return_value=rows,
                    ) as fetch,
                    self.assertRaises(
                        deribit.InstrumentSnapshotError
                    ) as failure,
                ):
                    collector._fetch_instrument_snapshot(
                        deadline=None
                    )
                self.assertEqual(
                    failure.exception.error_type,
                    "SchemaMismatch",
                )
                fetch.assert_called_once_with(
                    "btc",
                    deadline=None,
                )
                self.assertEqual(collector.instruments, {})
                self.assertFalse(
                    list(out.glob("day*/instruments.jsonl"))
                )
                self.assertFalse(
                    any(
                        channel.startswith("ticker.")
                        for channel in collector.default_channels()
                    )
                )

    def test_metadata_replay_rejects_name_row_identity_mismatch(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-name-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        strike_mismatch = {
            **self.instrument_payload(0, currency="btc"),
            "strike": 1.0,
        }
        type_mismatch = {
            **self.instrument_payload(1, currency="btc"),
            "option_type": "call",
        }
        records = []
        for index, row in enumerate(
            (strike_mismatch, type_mismatch)
        ):
            timestamp = received_ms + index
            fingerprint = builder._instrument_fingerprint(row)
            records.append({
                "record_type": "instrument",
                "dedup_id": builder._instrument_dedup_id(
                    row,
                    timestamp,
                ),
                "received_ts_ms": timestamp,
                "snapshot_id": f"snapshot-{index}",
                "instrument_fingerprint": fingerprint,
                "instrument": row,
            })
        builder.sink.close()
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        instrument_path.parent.mkdir(parents=True, exist_ok=True)
        instrument_path.write_text(
            "".join(
                deribit.strict_json_dumps(record) + "\n"
                for record in records
            ),
            encoding="utf-8",
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state()
        self.assertEqual(
            stats["metadata_replay_rejected"],
            len(records),
        )
        self.assertEqual(restarted.instruments, {})
        self.assertEqual(
            restarted._persisted_instrument_fingerprints,
            {},
        )

    def test_instrument_name_identity_accepts_equivalent_strikes(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        for strike in (90_000, 90_000.0, 9e4):
            with self.subTest(strike=repr(strike)):
                row = {
                    **self.instrument_payload(0, currency="btc"),
                    "strike": strike,
                }
                self.assertTrue(
                    collector._instrument_row_matches_schema(
                        row,
                        "btc",
                    )
                )
        self.assertTrue(
            collector._instrument_row_matches_schema(
                self.instrument_payload(1, currency="btc"),
                "btc",
            )
        )

    def test_instrument_snapshot_rejects_name_expiry_mismatch(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-expiry-identity-"
        ))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        mismatched = [
            {
                **self.instrument_payload(index, currency="btc"),
                "expiration_timestamp": 1_900_000_000_000,
            }
            for index in range(10)
        ]
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                side_effect=[mismatched, [dict(row) for row in mismatched]],
            ) as fetch,
            self.assertRaises(
                deribit.InstrumentSnapshotError
            ) as failure,
        ):
            collector._fetch_instrument_snapshot(deadline=None)
        self.assertEqual(
            failure.exception.error_type,
            "SchemaMismatch",
        )
        self.assertEqual(fetch.call_count, 1)
        self.assertEqual(collector.instruments, {})
        self.assertFalse(list(out.glob("day*/instruments.jsonl")))
        self.assertFalse(
            any(
                channel.startswith("ticker.")
                for channel in collector.default_channels()
            )
        )

    def test_metadata_replay_rejects_name_expiry_mismatch(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-expiry-replay-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        builder = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        row = {
            **self.instrument_payload(0, currency="btc"),
            "expiration_timestamp": 1_900_000_000_000,
        }
        fingerprint = builder._instrument_fingerprint(row)
        record = {
            "record_type": "instrument",
            "dedup_id": builder._instrument_dedup_id(
                row,
                received_ms,
            ),
            "received_ts_ms": received_ms,
            "snapshot_id": "mismatched-expiry",
            "instrument_fingerprint": fingerprint,
            "instrument": row,
        }
        builder.sink.close()
        instrument_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "instruments.jsonl"
        )
        instrument_path.parent.mkdir(parents=True, exist_ok=True)
        instrument_path.write_text(
            deribit.strict_json_dumps(record) + "\n",
            encoding="utf-8",
        )

        restarted = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(restarted.sink.close)
        stats = restarted.rehydrate_state()
        self.assertEqual(stats["metadata_replay_rejected"], 1)
        self.assertEqual(restarted.instruments, {})

    def test_restart_narrowing_currency_scope_excludes_prior_tickers(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-currency-scope-"
        ))
        received_ms = int(deribit.time.time() * 1000)
        btc = self.instrument_payload(0, currency="btc")
        eth = self.instrument_payload(1, currency="eth")
        writer = deribit.Collector(
            out,
            currencies=("btc", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        writer._apply_instrument_snapshot(
            {"btc": [btc], "eth": [eth]},
            received_ms,
        )
        writer.sink.close()
        day = deribit.utc_day(received_ms)
        historical_path = out / f"day{day}" / "instruments.jsonl"
        historical_bytes = historical_path.read_bytes()

        narrowed = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(narrowed.sink.close)
        stats = narrowed.rehydrate_state()
        self.assertEqual(stats["instrument_ids"], 1)
        self.assertEqual(set(narrowed.instruments), {btc["instrument_name"]})
        narrowed._apply_instrument_snapshot(
            {"btc": [dict(btc)]},
            received_ms + 1,
        )

        class CapturingWebSocket:
            def __init__(self) -> None:
                self.sent = []

            def send_text(self, text: str) -> None:
                self.sent.append(json.loads(text))

        websocket = CapturingWebSocket()
        narrowed.ws = websocket

        def accept(_request_id: int, channels: set[str]) -> bool:
            narrowed.active_channels |= channels
            return True

        with mock.patch.object(narrowed, "_await_ack", side_effect=accept):
            self.assertTrue(narrowed.subscribe_all())
        subscribed = {
            channel
            for request in websocket.sent
            for channel in request["params"]["channels"]
        }
        self.assertTrue(any(channel.startswith("ticker.BTC-") for channel in subscribed))
        self.assertFalse(any(channel.startswith("ticker.ETH-") for channel in subscribed))
        self.assertFalse(
            any(channel.startswith("ticker.ETH-") for channel in narrowed.active_channels)
        )

        next_day_ms = ((received_ms // 86_400_000) + 1) * 86_400_000 + 1_000
        narrowed._ensure_instrument_day_fence(next_day_ms, 1.0)
        self.assertFalse(
            any(channel.startswith("ticker.ETH-") for channel in narrowed._ticker_channels())
        )
        next_day_records = [
            json.loads(line)
            for line in (
                out
                / f"day{deribit.utc_day(next_day_ms)}"
                / "instruments.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        self.assertTrue(next_day_records)
        self.assertTrue(
            all(
                record["instrument"]["instrument_name"].startswith("BTC-")
                for record in next_day_records
                if record.get("record_type") == "instrument"
            )
        )
        self.assertEqual(historical_path.read_bytes(), historical_bytes)

    def test_ticker_channels_follow_normalized_currency_scope(self) -> None:
        btc = self.instrument_payload(0, currency="btc")
        eth = self.instrument_payload(1, currency="eth")
        eth_only = deribit.Collector(
            Path(tempfile.mkdtemp(prefix="collector-integrity-eth-scope-")),
            currencies=("ETH", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(eth_only.sink.close)
        eth_only.instruments = {
            btc["instrument_name"]: btc,
            eth["instrument_name"]: eth,
        }
        self.assertEqual(eth_only.currencies, ("eth",))
        self.assertEqual(
            eth_only._ticker_channels(),
            [f"ticker.{eth['instrument_name']}.agg2"],
        )

        dual = deribit.Collector(
            Path(tempfile.mkdtemp(prefix="collector-integrity-dual-scope-")),
            currencies=("btc", "eth"),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(dual.sink.close)
        dual.instruments = {
            btc["instrument_name"]: btc,
            eth["instrument_name"]: eth,
        }
        self.assertEqual(
            set(dual._ticker_channels()),
            {
                f"ticker.{btc['instrument_name']}.agg2",
                f"ticker.{eth['instrument_name']}.agg2",
            },
        )

    def test_instrument_snapshot_batch_failure_preserves_metadata_state(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        received_ms = int(deribit.time.time() * 1000)
        first = [
            {"instrument_name": "BTC-A", "creation_timestamp": 1},
            {"instrument_name": "BTC-B", "creation_timestamp": 2},
        ]
        collector._apply_instrument_snapshot(
            {"btc": first},
            received_ms,
        )
        prior_instruments = dict(collector.instruments)
        prior_fingerprints = dict(
            collector._persisted_instrument_fingerprints
        )
        with (
            mock.patch.object(
                collector.sink,
                "write_batch",
                side_effect=OSError("injected batch failure"),
            ),
            self.assertRaises(OSError),
        ):
            collector._apply_instrument_snapshot(
                {
                    "btc": [{
                        "instrument_name": "BTC-A",
                        "creation_timestamp": 1,
                        "tick_size": 0.2,
                    }],
                },
                received_ms + 1,
            )
        self.assertEqual(collector.instruments, prior_instruments)
        self.assertEqual(
            collector._persisted_instrument_fingerprints,
            prior_fingerprints,
        )

    def test_instrument_snapshot_rejects_duplicate_names(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        rows = [
            {"instrument_name": "BTC-DUP", "creation_timestamp": 1},
            {"instrument_name": "BTC-DUP", "creation_timestamp": 2},
        ]
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=rows,
            ),
            self.assertRaises(deribit.InstrumentSnapshotError),
        ):
            collector._fetch_instrument_snapshot(deadline=None)

    def test_daily_instrument_snapshots_have_distinct_state_ids(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-seed-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        rows = [
            self.instrument_payload(index)
            for index in range(10)
        ]
        first_ms = 1_800_000_000_000
        second_ms = first_ms + 86_400_000
        with (
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=rows,
            ),
            mock.patch.object(
                deribit.time,
                "time",
                side_effect=[
                    first_ms / 1000.0,
                    (first_ms + 1_000) / 1000.0,
                    second_ms / 1000.0,
                ],
            ),
        ):
            collector.seed_instruments()
            collector.seed_instruments()
            collector.seed_instruments()

        first_path = (
            out
            / f"day{deribit.utc_day(first_ms)}"
            / "instruments.jsonl"
        )
        second_path = (
            out
            / f"day{deribit.utc_day(second_ms)}"
            / "instruments.jsonl"
        )
        first_rows = first_path.read_text(encoding="utf-8").splitlines()
        second_rows = second_path.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(first_rows), len(rows) + 1)
        self.assertEqual(len(second_rows), len(rows) + 1)
        self.assertNotEqual(
            json.loads(first_rows[0])["dedup_id"],
            json.loads(second_rows[0])["dedup_id"],
        )

    def test_periodic_chain_refresh_does_not_block_websocket_frames(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-refresh-"))
        collector = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
        )
        self.addCleanup(collector.sink.close)
        fetch_started = deribit.threading.Event()
        fetch_release = deribit.threading.Event()
        fetch_completed = deribit.threading.Event()

        def blocked_fetch(_currency: str, deadline=None):
            del deadline
            fetch_started.set()
            fetch_release.wait(1.0)
            fetch_completed.set()
            return [
                self.instrument_payload(index)
                for index in range(10)
            ]

        class HeartbeatWebSocket:
            def __init__(self) -> None:
                self.sent = []

            def recv_message(self, timeout: float):
                del timeout
                collector.stop = True
                return (
                    deribit.OP_TEXT,
                    json.dumps({"method": "heartbeat"}).encode(),
                )

            def send_text(self, payload: str) -> None:
                self.sent.append(json.loads(payload))

        now_ms = int(deribit.time.time() * 1000)
        collector.ws = HeartbeatWebSocket()
        collector.last_message_ms = now_ms
        collector.last_message_monotonic = deribit.time.monotonic()
        collector.last_utc_day = deribit.utc_day(now_ms)
        collector.instruments = {
            self.instrument_name: self.instrument_payload(
                name=self.instrument_name
            )
        }
        collector._instrument_persist_day = deribit.utc_day(now_ms)
        collector.last_chain_reconcile = (
            deribit.time.monotonic() - deribit.CHAIN_REFRESH_SECONDS - 1
        )
        with mock.patch.object(
            deribit,
            "fetch_instruments",
            side_effect=blocked_fetch,
        ):
            try:
                collector._loop()
                self.assertTrue(fetch_started.wait(0.5))
                self.assertFalse(fetch_completed.is_set())
                self.assertEqual(
                    collector.ws.sent,
                    [{"jsonrpc": "2.0", "method": "public/test"}],
                )
            finally:
                fetch_release.set()
                refresh_thread = getattr(
                    collector,
                    "_chain_refresh_thread",
                    None,
                )
                if refresh_thread is not None:
                    refresh_thread.join(timeout=1.0)

    def test_stale_authoritative_chain_revokes_session_health(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        now_ms = int(deribit.time.time() * 1000)
        now_monotonic = 100.0 + 3 * deribit.CHAIN_REFRESH_SECONDS + 1
        collector.last_utc_day = deribit.utc_day(now_ms)
        collector.last_chain_reconcile = now_monotonic
        collector._next_chain_refresh_attempt = float("inf")
        collector._last_authoritative_chain_refresh_monotonic = 100.0
        ticker_channel = f"ticker.{self.instrument_name}.agg2"
        index_channel = "deribit_price_index.btc_usd"
        collector.active_channels = {ticker_channel, index_channel}
        collector._required_last_received_monotonic = {
            ticker_channel: now_monotonic,
            index_channel: now_monotonic,
        }
        with mock.patch.object(
            deribit.time,
            "monotonic",
            return_value=now_monotonic,
        ):
            with self.assertRaises(ConnectionError):
                collector._maybe_refresh_chain(now_ms)
            with self.assertRaises(ConnectionError):
                collector._maybe_refresh_chain(now_ms)
        event_path = (
            out
            / f"day{deribit.utc_day(now_ms)}"
            / "events.jsonl"
        )
        stale_events = [
            row
            for row in (
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("event") == "instrument_chain_stale"
        ]
        self.assertEqual(len(stale_events), 1)

    def test_authoritative_chain_recovers_after_transient_failure(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        now_ms = int(deribit.time.time() * 1000)
        now_monotonic = 1_000.0
        collector._last_authoritative_chain_refresh_monotonic = (
            now_monotonic - deribit.CHAIN_REFRESH_SECONDS
        )
        collector._chain_refresh_result = (
            collector._chain_refresh_generation,
            None,
            deribit.InstrumentSnapshotError(
                "btc",
                "TimeoutError",
                "transient refresh timeout",
            ),
        )
        with mock.patch.object(
            deribit.time,
            "monotonic",
            return_value=now_monotonic,
        ):
            collector._poll_chain_refresh(now_ms)
        self.assertFalse(
            getattr(collector, "_instrument_chain_stale_latched", False)
        )

        collector._instrument_chain_stale_latched = True
        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                return_value=now_monotonic + 1,
            ),
            mock.patch.object(
                deribit,
                "fetch_instruments",
                return_value=[
                    self.instrument_payload(
                        0,
                        name=self.instrument_name,
                    )
                ],
            ),
            mock.patch.object(
                deribit.time,
                "time",
                return_value=(now_ms + 1) / 1000.0,
            ),
        ):
            collector.seed_instruments()
        self.assertEqual(
            collector._last_authoritative_chain_refresh_monotonic,
            now_monotonic + 1,
        )
        self.assertFalse(collector._instrument_chain_stale_latched)
        event_path = (
            out
            / f"day{deribit.utc_day(now_ms)}"
            / "events.jsonl"
        )
        recovered = [
            row
            for row in (
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("event") == "instrument_chain_recovered"
        ]
        self.assertEqual(len(recovered), 1)

    def test_chain_refresh_retries_at_retry_deadline(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.sink.close)
        now_ms = int(deribit.time.time() * 1000)
        failure_monotonic = 1_000.0
        collector.last_utc_day = deribit.utc_day(now_ms)
        collector.last_chain_reconcile = failure_monotonic - 1
        collector._last_authoritative_chain_refresh_monotonic = (
            failure_monotonic - 1
        )
        collector._chain_refresh_result = (
            collector._chain_refresh_generation,
            None,
            deribit.InstrumentSnapshotError(
                "btc",
                "TimeoutError",
                "transient refresh timeout",
            ),
        )
        with mock.patch.object(
            deribit.time,
            "monotonic",
            return_value=failure_monotonic,
        ):
            collector._poll_chain_refresh(now_ms)
        start_refresh = mock.Mock()
        with (
            mock.patch.object(
                deribit.time,
                "monotonic",
                return_value=(
                    failure_monotonic
                    + deribit.CHAIN_REFRESH_RETRY_SECONDS
                ),
            ),
            mock.patch.object(
                collector,
                "_start_chain_refresh",
                start_refresh,
            ),
        ):
            collector._maybe_refresh_chain(now_ms)
        start_refresh.assert_called_once_with()

    def test_deribit_rows_have_run_and_connection_provenance(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-provenance-"))
        base_ms = int(deribit.time.time() * 1000)

        class FakeWebSocket:
            def close(self) -> None:
                return None

        def exercise_run(
            collector: object,
            sequence_start: int,
            connections: int,
        ) -> int:
            collector.instruments = {
                self.instrument_name: {
                    "instrument_name": self.instrument_name,
                    "kind": "option",
                }
            }
            calls = 0

            def loop() -> None:
                nonlocal calls
                calls += 1
                received_ms = base_ms + sequence_start + calls
                collector._route(
                    "trades.option.BTC.100ms",
                    [self.trade_payload(sequence_start + calls, received_ms=received_ms)],
                    received_ms,
                )
                collector._route(
                    f"ticker.{self.instrument_name}.agg2",
                    self.ticker_payload(timestamp=received_ms),
                    received_ms,
                    received_monotonic=100.0 + received_ms / 1000,
                )
                collector._route(
                    "deribit_price_index.btc_usd",
                    {
                        "timestamp": received_ms,
                        "price": 90_000.0 + calls,
                        "index_name": "btc_usd",
                    },
                    received_ms,
                    received_monotonic=100.0 + received_ms / 1000,
                )
                if calls < connections:
                    raise ConnectionError("injected reconnect")
                collector.stop = True

            with (
                mock.patch.object(collector, "seed_instruments"),
                mock.patch.object(
                    collector,
                    "_connect_websocket_interruptibly",
                    side_effect=lambda: FakeWebSocket(),
                ),
                mock.patch.object(
                    collector,
                    "subscribe_all",
                    return_value=True,
                ),
                mock.patch.object(
                    collector,
                    "await_required_coverage",
                    return_value=True,
                ),
                mock.patch.object(collector, "_loop", side_effect=loop),
                mock.patch.object(deribit, "stop_aware_sleep"),
            ):
                return collector.run()

        first = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
            ticker_sample_seconds=0.001,
        )
        self.addCleanup(first.sink.close)
        self.assertEqual(exercise_run(first, 0, 2), 0)
        first_run_id = first.run_id

        second = deribit.Collector(
            out,
            currencies=("btc",),
            day_cap=10_000_000,
            total_cap=100_000_000,
            heartbeat=False,
            ticker_sample_seconds=0.001,
        )
        self.addCleanup(second.sink.close)
        self.assertEqual(exercise_run(second, 2, 1), 0)
        self.assertNotEqual(second.run_id, first_run_id)

        day_dir = out / f"day{deribit.utc_day(base_ms)}"
        market_rows = []
        for filename in ("trades.jsonl", "ticker.jsonl", "index.jsonl"):
            market_rows.extend(
                json.loads(line)
                for line in (day_dir / filename)
                .read_text(encoding="utf-8")
                .splitlines()
            )
        self.assertEqual(len(market_rows), 9)
        connection_ids = {
            row["connection_id"]
            for row in market_rows
        }
        self.assertEqual(len(connection_ids), 3)
        self.assertEqual(
            {row["run_id"] for row in market_rows},
            {first_run_id, second.run_id},
        )
        for row in market_rows:
            self.assertTrue(
                row["connection_id"].startswith(f"{row['run_id']}:c")
            )

        events = [
            json.loads(line)
            for line in (day_dir / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        for connection_id in connection_ids:
            lifecycle = {
                row["event"]
                for row in events
                if row.get("connection_id") == connection_id
            }
            self.assertTrue({
                "connection_connect",
                "connection_ready",
                "connection_disconnect",
            }.issubset(lifecycle))
        process_events = {
            "day_rotation_open",
            "restart_state_rehydrated",
            "session_stop",
        }
        self.assertTrue(
            all(
                "connection_id" not in row
                for row in events
                if row.get("event") in process_events
            )
        )
        self.assertTrue(all("run_id" in row for row in events))

    def test_deribit_sink_owns_run_id_and_unique_rotation_identity(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-run-id-"))
        received_ms = 1_800_000_000_000
        for run_id in ("run-a", "run-b"):
            sink = deribit.Sink(
                out,
                10_000_000,
                100_000_000,
                run_id=run_id,
            )
            try:
                sink.write(
                    "events",
                    {
                        "record_type": "capture_event",
                        "dedup_id": f"probe-{run_id}",
                        "event": "probe",
                        "received_ts_ms": received_ms,
                        "run_id": "caller-must-not-own-this",
                    },
                    now_ms=received_ms,
                )
            finally:
                sink.close()
        events_path = (
            out
            / f"day{deribit.utc_day(received_ms)}"
            / "events.jsonl"
        )
        events = [
            json.loads(line)
            for line in events_path.read_text(encoding="utf-8").splitlines()
        ]
        rotations = [
            row
            for row in events
            if row.get("event") == "day_rotation_open"
        ]
        self.assertEqual(
            {row["run_id"] for row in rotations},
            {"run-a", "run-b"},
        )
        self.assertEqual(len({row["dedup_id"] for row in rotations}), 2)
        probes = [row for row in events if row.get("event") == "probe"]
        self.assertEqual(
            {row["run_id"] for row in probes},
            {"run-a", "run-b"},
        )

    def test_deribit_output_lease_rejects_real_process_contender(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-lease-"))
        day_dir = out / "day20260904"
        day_dir.mkdir(parents=True)
        capture_path = day_dir / "trades.jsonl"
        original = b"authoritative-capture-bytes\n"
        capture_path.write_bytes(original)
        lock_path = out / deribit.OUTPUT_LEASE_FILE
        holder = start_output_lock_holder(lock_path)
        self.addCleanup(stop_process, holder)
        result = subprocess.run(
            [
                sys.executable,
                str(HERE / "focused-deribit-gex" / "capture.py"),
                "--raw-dir",
                str(out),
                "--max-minutes",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=10,
            env={
                **deribit.os.environ,
                "PYTHONDONTWRITEBYTECODE": "1",
            },
        )
        self.assertEqual(result.returncode, 6)
        self.assertIn("DERIBIT_CAPTURE_LOCKED", result.stderr)
        self.assertNotIn(deribit.READY_BANNER, result.stdout)
        self.assertEqual(capture_path.read_bytes(), original)
        self.assertFalse(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).exists()
        )

        stop_process(holder)
        successor = deribit.Sink(
            out,
            day_cap_bytes=10_000_000,
            total_cap_bytes=100_000_000,
        )
        try:
            self.assertEqual(successor.prior_bytes, len(original))
            self.assertEqual(successor.day_bytes, 0)
        finally:
            successor.close()
        self.assertTrue(lock_path.is_file())

    def test_deribit_accounting_stat_failure_is_fail_closed(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-stat-"))
        day_dir = out / "day20260904"
        day_dir.mkdir(parents=True)
        capture_path = day_dir / "trades.jsonl"
        original = b"authoritative-capture-bytes\n"
        capture_path.write_bytes(original)
        original_lstat = deribit.os.lstat

        def guarded_lstat(target, *args, **kwargs):
            if Path(target) == capture_path:
                raise PermissionError("injected capture stat failure")
            return original_lstat(target, *args, **kwargs)

        with (
            mock.patch.object(deribit.os, "lstat", guarded_lstat),
            self.assertRaises(PermissionError),
        ):
            deribit.Sink(out, 10_000_000, 100_000_000)
        self.assertEqual(capture_path.read_bytes(), original)

    def test_deribit_special_capture_paths_are_rejected(self) -> None:
        received_ms = int(deribit.time.time() * 1000)
        day = deribit.utc_day(received_ms)
        outside = Path(tempfile.mkdtemp(prefix="collector-integrity-outside-"))
        external = outside / "external-target.jsonl"
        external_bytes = b'{"record_type":"external"}\n'

        for label in ("symlink", "fifo"):
            with self.subTest(capture=label):
                external.write_bytes(external_bytes)
                out = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-deribit-special-"
                ))
                day_dir = out / f"day{day}"
                day_dir.mkdir(parents=True)
                capture_path = day_dir / "trades.jsonl"
                if label == "symlink":
                    capture_path.symlink_to(external)
                else:
                    os.mkfifo(capture_path)
                with self.assertRaises(
                    deribit.NotARegularCaptureFileError
                ):
                    deribit.Sink(out, 10_000_000, 100_000_000)
                self.assertEqual(external.read_bytes(), external_bytes)

    def test_deribit_symlinked_day_directory_is_rejected(self) -> None:
        received_ms = int(deribit.time.time() * 1000)
        day = deribit.utc_day(received_ms)
        outside = Path(tempfile.mkdtemp(prefix="collector-integrity-outside-"))
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-deribit-day-link-"
        ))
        (out / f"day{day}").symlink_to(outside, target_is_directory=True)
        with self.assertRaises(deribit.NotARegularCaptureFileError):
            deribit.Sink(out, 10_000_000, 100_000_000)
        self.assertEqual(list(outside.iterdir()), [])

    def test_deribit_checkpoint_temp_symlink_is_replaced(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.sink.close)
        outside = Path(tempfile.mkdtemp(prefix="collector-integrity-outside-"))
        external = outside / "external-checkpoint.json"
        external_bytes = b'{"external":true}\n'
        external.write_bytes(external_bytes)
        temp_path = out / f".{deribit.TRADE_SEQUENCE_STATE_FILE}.tmp"
        temp_path.symlink_to(external)
        collector._sequence_state_dirty = True
        self.assertTrue(collector._save_sequence_state())
        self.assertEqual(external.read_bytes(), external_bytes)
        self.assertIsNone(collector._sequence_state_error)
        checkpoint = deribit.strict_json_loads(
            (out / deribit.TRADE_SEQUENCE_STATE_FILE).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            checkpoint["schema_version"],
            deribit.TRADE_SEQUENCE_STATE_SCHEMA_VERSION,
        )

    def test_new_output_roots_fsync_their_parent(self) -> None:
        for module in (deribit, liquidation):
            with self.subTest(module=module.__name__):
                parent = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-root-dir-"
                ))
                root = parent / "new-output"
                with mock.patch.object(
                    module,
                    "fsync_directory",
                ) as directory_fsync:
                    lease = module.OutputLease(root)
                try:
                    directory_fsync.assert_called_once_with(parent)
                finally:
                    lease.close()

    def test_nested_output_roots_fsync_every_created_ancestor(
        self,
    ) -> None:
        for module in (deribit, liquidation):
            with self.subTest(module=module.__name__):
                anchor = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-nested-root-"
                ))
                root = anchor / "a" / "b" / "out"
                operations = []
                original_open = module.os.open

                def track_open(path, *args, **kwargs):
                    operations.append(("lease_open", Path(path)))
                    return original_open(path, *args, **kwargs)

                with (
                    mock.patch.object(
                        module,
                        "fsync_directory",
                        side_effect=lambda path: operations.append(
                            ("parent_fsync", Path(path))
                        ),
                    ),
                    mock.patch.object(
                        module.os,
                        "open",
                        side_effect=track_open,
                    ),
                ):
                    lease = module.OutputLease(root)
                try:
                    self.assertEqual(
                        operations,
                        [
                            ("parent_fsync", anchor),
                            ("parent_fsync", anchor / "a"),
                            ("parent_fsync", anchor / "a" / "b"),
                            (
                                "lease_open",
                                root / module.OUTPUT_LEASE_FILE,
                            ),
                        ],
                    )
                finally:
                    lease.close()

    def test_output_root_parent_fsync_failure_precedes_lease_open(
        self,
    ) -> None:
        for module in (deribit, liquidation):
            for failure_position in (1, 2, 3):
                with self.subTest(
                    module=module.__name__,
                    failure_position=failure_position,
                ):
                    anchor = Path(tempfile.mkdtemp(
                        prefix="collector-integrity-root-fsync-failure-"
                    ))
                    root = anchor / "a" / "b" / "out"
                    fsync_calls = []
                    lease_open = mock.Mock(
                        side_effect=AssertionError(
                            "lease opened after ancestor fsync failure"
                        )
                    )

                    def fail_selected_parent(path):
                        fsync_calls.append(Path(path))
                        if len(fsync_calls) == failure_position:
                            raise OSError(
                                "injected ancestor fsync failure"
                            )

                    with (
                        mock.patch.object(
                            module,
                            "fsync_directory",
                            side_effect=fail_selected_parent,
                        ),
                        mock.patch.object(
                            module.os,
                            "open",
                            lease_open,
                        ),
                        self.assertRaisesRegex(
                            OSError,
                            "injected ancestor fsync failure",
                        ),
                    ):
                        module.OutputLease(root)
                    lease_open.assert_not_called()
                    self.assertFalse(
                        (root / module.OUTPUT_LEASE_FILE).exists()
                    )

    def test_existing_output_root_needs_no_creation_fsync(self) -> None:
        for module in (deribit, liquidation):
            with self.subTest(module=module.__name__):
                root = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-existing-root-"
                ))
                with mock.patch.object(
                    module,
                    "fsync_directory",
                ) as directory_fsync:
                    lease = module.OutputLease(root)
                try:
                    directory_fsync.assert_not_called()
                finally:
                    lease.close()

    def test_deribit_new_entries_fsync_directories_in_order(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-dir-"))
        sink = deribit.Sink(out, 10_000_000, 100_000_000)
        self.addCleanup(sink.close)
        first_ms = 1_800_000_000_000
        operations = []

        def file_fsync(fd: int) -> None:
            operations.append(("file", fd))

        def directory_fsync(path: Path) -> None:
            operations.append(("directory", Path(path)))

        record = {
            "record_type": "capture_event",
            "dedup_id": "directory-order",
            "event": "probe",
            "received_ts_ms": first_ms,
        }
        with (
            mock.patch.object(
                deribit.os,
                "fsync",
                side_effect=file_fsync,
            ),
            mock.patch.object(
                deribit,
                "fsync_directory",
                side_effect=directory_fsync,
            ),
        ):
            sink.write("trades", record, first_ms)
        first_day = out / f"day{deribit.utc_day(first_ms)}"
        self.assertEqual(
            operations,
            [
                ("file", sink._handles["events"].fileno()),
                ("directory", first_day),
                ("directory", out),
                ("file", sink._handles["trades"].fileno()),
            ],
        )

        operations.clear()
        with (
            mock.patch.object(
                deribit.os,
                "fsync",
                side_effect=file_fsync,
            ),
            mock.patch.object(
                deribit,
                "fsync_directory",
                side_effect=directory_fsync,
            ),
        ):
            sink.write(
                "trades",
                {**record, "dedup_id": "existing-append"},
                first_ms + 1,
            )
        self.assertEqual(
            operations,
            [("file", sink._handles["trades"].fileno())],
        )

        operations.clear()
        second_ms = first_ms + 86_400_000
        with (
            mock.patch.object(
                deribit.os,
                "fsync",
                side_effect=file_fsync,
            ),
            mock.patch.object(
                deribit,
                "fsync_directory",
                side_effect=directory_fsync,
            ),
        ):
            sink.write(
                "trades",
                {**record, "dedup_id": "next-day"},
                second_ms,
            )
        second_day = out / f"day{deribit.utc_day(second_ms)}"
        directories = [
            operation
            for operation in operations
            if operation[0] == "directory"
        ]
        self.assertEqual(
            directories,
            [
                ("directory", second_day),
                ("directory", out),
            ],
        )
        day_sync_index = operations.index(directories[0])
        self.assertEqual(
            operations[day_sync_index - 1],
            ("file", sink._handles["events"].fileno()),
        )
        self.assertEqual(
            operations[-1],
            ("file", sink._handles["trades"].fileno()),
        )

    def test_deribit_directory_fsync_failure_is_terminal_before_ready(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        collector.max_minutes = 1.0
        stdout = io.StringIO()
        with (
            mock.patch.object(deribit.signal, "signal"),
            mock.patch.object(
                deribit,
                "fsync_directory",
                side_effect=OSError("injected directory fsync failure"),
            ),
            mock.patch.object(deribit.sys, "stdout", stdout),
            mock.patch.object(deribit.sys, "stderr", io.StringIO()),
        ):
            rc = collector.run()
        self.assertEqual(rc, 5)
        self.assertNotIn(deribit.READY_BANNER, stdout.getvalue())
        self.assertEqual(collector._last_trade_seq, {})
        self.assertEqual(collector.deduper.seen, set())
        self.assertIsNotNone(collector.sink._io_error)

    def test_deribit_fsync_failure_latches_sink_terminal(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-io-"))
        sink = deribit.Sink(out, 10_000_000, 100_000_000)
        self.addCleanup(sink.close)
        now_ms = 1_800_000_000_000
        first = {
            "record_type": "capture_event",
            "dedup_id": "before-io-failure",
            "event": "probe",
            "received_ts_ms": now_ms,
        }
        sink.write("events", first, now_ms)
        with (
            mock.patch.object(
                deribit.os,
                "fsync",
                side_effect=OSError("injected fsync failure"),
            ),
            self.assertRaises(OSError),
        ):
            sink.write("events", {**first, "dedup_id": "failed-write"}, now_ms)
        frozen_size = sum(
            path.stat().st_size for path in out.glob("day*/*.jsonl")
        )
        with self.assertRaises(deribit.DayCapError):
            sink.write("events", {**first, "dedup_id": "must-not-write"}, now_ms)
        self.assertEqual(
            sum(path.stat().st_size for path in out.glob("day*/*.jsonl")),
            frozen_size,
        )

    def test_first_rotation_reconciles_against_target_day(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-cap-"))
        base_seconds = 1_800_000_000.0
        first_ms = int(base_seconds * 1000)
        second_ms = first_ms + 86_400_000
        first_day = deribit.utc_day(first_ms)
        prior_dir = out / f"day{first_day}"
        prior_dir.mkdir(parents=True)
        (prior_dir / "trades.jsonl").write_bytes(b"x" * 500)

        with mock.patch.object(deribit.time, "time", return_value=base_seconds):
            sink = deribit.Sink(out, 10_000_000, 650)
        self.addCleanup(sink.close)
        record = {
            "record_type": "capture_event",
            "dedup_id": "target-day-cap-probe",
            "event": "probe",
            "received_ts_ms": second_ms,
        }
        with self.assertRaises(deribit.DayCapError):
            sink.write("events", record, second_ms)
        self.assertEqual(sink.prior_bytes, 500)

    def test_deribit_rotation_uses_initiating_write_timestamp(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-clock-"))
        received_ms = 1_800_000_000_000
        later_seconds = received_ms / 1000 + 86_400
        sink = deribit.Sink(out, 10_000_000, 100_000_000)
        try:
            with mock.patch.object(
                deribit.time,
                "time",
                return_value=later_seconds,
            ):
                sink.write(
                    "events",
                    {
                        "record_type": "capture_event",
                        "dedup_id": "initiating-write",
                        "event": "probe",
                        "received_ts_ms": received_ms,
                    },
                    received_ms,
                )
        finally:
            sink.close()
        day = deribit.utc_day(received_ms)
        rows = [
            json.loads(line)
            for line in (
                out / f"day{day}" / "events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            [row["received_ts_ms"] for row in rows],
            [received_ms, received_ms],
        )
        self.assertEqual(rows[0]["day_utc"], day)

    def test_same_day_sink_reopens_have_distinct_rotation_ids(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-rotate-"))
        base_seconds = 1_800_000_000.0
        now_ms = int(base_seconds * 1000)
        with mock.patch.object(
            deribit.time,
            "time",
            side_effect=[
                base_seconds,
                base_seconds + 0.1,
                base_seconds + 1.0,
                base_seconds + 1.1,
            ],
        ):
            first = deribit.Sink(out, 10_000_000, 100_000_000)
            try:
                first.write(
                    "events",
                    {
                        "record_type": "capture_event",
                        "dedup_id": "probe-one",
                        "event": "probe",
                        "received_ts_ms": now_ms,
                    },
                    now_ms,
                )
            finally:
                first.close()
            second = deribit.Sink(out, 10_000_000, 100_000_000)
            try:
                second.write(
                    "events",
                    {
                        "record_type": "capture_event",
                        "dedup_id": "probe-two",
                        "event": "probe",
                        "received_ts_ms": now_ms + 1_000,
                    },
                    now_ms + 1_000,
                )
            finally:
                second.close()

        event_path = out / f"day{deribit.utc_day(now_ms)}" / "events.jsonl"
        rotation_rows = [
            row
            for row in (
                json.loads(line)
                for line in event_path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("event") == "day_rotation_open"
        ]
        self.assertEqual(len(rotation_rows), 2)
        self.assertNotEqual(
            rotation_rows[0]["received_ts_ms"],
            rotation_rows[1]["received_ts_ms"],
        )
        self.assertNotEqual(
            rotation_rows[0]["dedup_id"],
            rotation_rows[1]["dedup_id"],
        )

    def test_deribit_reopen_isolates_partial_jsonl_tail(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-deribit-tail-"))
        now_ms = int(deribit.time.time() * 1000)
        day = deribit.utc_day(now_ms)
        day_dir = out / f"day{day}"
        day_dir.mkdir(parents=True)
        trade_path = day_dir / "trades.jsonl"
        partial = b'{"record_type":"option_trade"'
        trade_path.write_bytes(partial)
        sink = deribit.Sink(out, 10_000_000, 100_000_000)
        try:
            sink.write(
                "trades",
                {
                    "record_type": "option_trade",
                    "dedup_id": "after-partial-tail",
                },
                now_ms,
            )
            disk_bytes = sum(
                path.stat().st_size
                for path in out.glob("day*/*.jsonl")
            )
            self.assertEqual(
                sink.prior_bytes + sink.day_bytes,
                disk_bytes,
            )
        finally:
            sink.close()
        raw = trade_path.read_bytes()
        self.assertTrue(raw.startswith(partial + b"\n"))
        lines = raw.splitlines()
        self.assertEqual(lines[0], partial)
        self.assertEqual(
            json.loads(lines[-1])["dedup_id"],
            "after-partial-tail",
        )
        event_path = day_dir / "events.jsonl"
        rotation = json.loads(
            event_path.read_text(encoding="utf-8").splitlines()[0]
        )
        self.assertEqual(
            rotation["recovered_partial_tail_kinds"],
            ["trades"],
        )


class LiquidationStateIntegrityTests(unittest.TestCase):
    def make_collector(self):
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-liquidation-"))

        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        collector.load_state()
        self.addCleanup(collector.release_lease)
        return collector, out

    @staticmethod
    def force_order_payload(
        event_ms: int = 1_800_000_000_000,
        symbol: str = "BTCUSDT",
    ) -> dict:
        return {
            "e": "forceOrder",
            "E": event_ms,
            "o": {
                "s": symbol,
                "S": "SELL",
                "o": "LIMIT",
                "f": "IOC",
                "q": "1.0",
                "p": "90000.0",
                "ap": "89999.0",
                "X": "FILLED",
                "l": "1.0",
                "z": "1.0",
                "T": event_ms - 1,
            },
        }
    @staticmethod
    def oi_payload(
        *,
        symbol: str = "BTCUSDT",
        timestamp: object = 1_800_000_000_000,
        sum_open_interest: object = "123.456",
        sum_open_interest_value: object = "11111111.25",
    ) -> dict:
        return {
            "symbol": symbol,
            "sumOpenInterest": sum_open_interest,
            "sumOpenInterestValue": sum_open_interest_value,
            "timestamp": timestamp,
        }

    def test_force_order_requires_positive_original_quantity(self) -> None:
        event_ms = 1_800_000_000_000
        for text in ("0", "0.0", "0e9", "-0"):
            with self.subTest(quantity=text):
                payload = self.force_order_payload(event_ms)
                payload["o"]["q"] = text
                payload["o"]["z"] = "0"
                payload["o"]["l"] = "0"
                self.assertEqual(
                    liquidation.force_order_schema_error(payload),
                    "order_quantity_nonpositive",
                )
                self.assertIsNone(
                    liquidation.canonical_force_order(
                        "!forceOrder@arr",
                        payload,
                    )
                )
        unfilled = self.force_order_payload(event_ms)
        unfilled["o"]["q"] = "2.5"
        unfilled["o"]["z"] = "0"
        unfilled["o"]["l"] = "0"
        self.assertIsNone(
            liquidation.force_order_schema_error(unfilled)
        )

        zero_quantity = self.force_order_payload(event_ms)
        zero_quantity["o"]["q"] = "0"
        zero_quantity["o"]["z"] = "0"
        zero_quantity["o"]["l"] = "0"
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        with (
            mock.patch.object(collector, "schedule_oi_snapshot") as oi,
            self.assertRaises(liquidation.ForceOrderSchemaError),
        ):
            collector.dispatch_market_message(
                zero_quantity,
                "c1",
                event_ms + 5,
            )
        oi.assert_not_called()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(any(row.get("kind") == "event" for row in rows))
        self.assertEqual(collector.seen_ids, set())
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            ["order_quantity_nonpositive"],
        )

        replay, replay_out = self.make_collector()
        self.addCleanup(replay.fsync_close)
        path = (
            replay_out
            / f"{liquidation.FILE_PREFIX}{liquidation.utc_date()}.jsonl"
        )
        path.write_text(
            json.dumps(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "event",
                    "ts_utc": "2027-01-15T08:00:00.000Z",
                    "recv_epoch_ms": event_ms + 5,
                    "conn": "c1",
                    "subscription_phase": "confirmed",
                    "stream": "!forceOrder@arr",
                    "event_id": liquidation.event_id_for(
                        "!forceOrder@arr",
                        self.force_order_payload(event_ms),
                    ),
                    "source": {
                        "venue": "binance-usds-m-futures",
                        "endpoint_kind": (
                            "public-unauthenticated-websocket"
                        ),
                        "url": replay.url,
                    },
                    "payload": zero_quantity,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        replay.seed_dedup_from_file(path)
        self.assertEqual(
            replay.stats["dedup_rehydrate_semantic_rejected"],
            1,
        )
        self.assertEqual(replay.seen_ids, set())

    def test_force_order_identity_delimiter_is_rejected(self) -> None:
        event_ms = 1_800_000_000_000
        ambiguous_a = self.force_order_payload(event_ms)
        ambiguous_a["o"]["s"] = "BTC|BUY"
        ambiguous_a["o"]["S"] = "SELL"
        ambiguous_a["o"]["o"] = "LIMIT"
        ambiguous_b = self.force_order_payload(event_ms)
        ambiguous_b["o"]["s"] = "BTC"
        ambiguous_b["o"]["S"] = "BUY"
        ambiguous_b["o"]["o"] = "SELL|LIMIT"
        self.assertEqual(
            liquidation.force_order_schema_error(ambiguous_a),
            "order_s_delimiter",
        )
        self.assertEqual(
            liquidation.force_order_schema_error(ambiguous_b),
            "order_o_delimiter",
        )
        for payload in (ambiguous_a, ambiguous_b):
            self.assertIsNone(
                liquidation.canonical_force_order(
                    "!forceOrder@arr",
                    payload,
                )
            )
            with self.assertRaises(liquidation.ForceOrderSchemaError):
                liquidation.event_id_for("!forceOrder@arr", payload)

        valid = self.force_order_payload(event_ms)
        order = valid["o"]
        expected = liquidation.hashlib.sha256(
            "|".join(
                [
                    "!forceOrder@arr",
                    str(valid["E"]),
                    *(
                        str(order[key])
                        for key in (
                            "s",
                            "S",
                            "o",
                            "f",
                            "q",
                            "p",
                            "ap",
                            "X",
                            "l",
                            "z",
                            "T",
                        )
                    ),
                ]
            ).encode("utf-8")
        ).hexdigest()[:32]
        self.assertEqual(
            liquidation.event_id_for("!forceOrder@arr", valid),
            expected,
        )

        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        with (
            mock.patch.object(collector, "schedule_oi_snapshot") as oi,
            self.assertRaises(liquidation.ForceOrderSchemaError),
        ):
            collector.dispatch_market_message(
                ambiguous_a,
                "c1",
                event_ms + 5,
            )
        oi.assert_not_called()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(any(row.get("kind") == "event" for row in rows))
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            ["order_s_delimiter"],
        )

        replay, replay_out = self.make_collector()
        self.addCleanup(replay.fsync_close)
        path = (
            replay_out
            / f"{liquidation.FILE_PREFIX}{liquidation.utc_date()}.jsonl"
        )
        path.write_text(
            json.dumps(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "event",
                    "ts_utc": "2027-01-15T08:00:00.000Z",
                    "recv_epoch_ms": event_ms + 5,
                    "conn": "c1",
                    "subscription_phase": "confirmed",
                    "stream": "!forceOrder@arr",
                    "event_id": expected,
                    "source": {
                        "venue": "binance-usds-m-futures",
                        "endpoint_kind": (
                            "public-unauthenticated-websocket"
                        ),
                        "url": replay.url,
                    },
                    "payload": ambiguous_a,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        replay.seed_dedup_from_file(path)
        self.assertEqual(
            replay.stats["dedup_rehydrate_semantic_rejected"],
            1,
        )
        self.assertEqual(replay.seen_ids, set())

        self.assertIsNotNone(
            liquidation.stream_config_error(["bt|c@forceOrder"])
        )
        with self.assertRaises(ValueError):
            liquidation.Collector(
                out_dir=Path(tempfile.mkdtemp(
                    prefix="collector-delimiter-stream-"
                )),
                url="wss://example.invalid/ws",
                streams=["bt|c@forceOrder"],
                max_seconds=None,
                max_bytes_per_day=10_000_000,
                max_total_bytes=100_000_000,
            )

    def test_malformed_order_object_returns_bounded_reason(self) -> None:
        event_ms = 1_800_000_000_000
        for label, order in (
            ("missing", None),
            ("string", "bad"),
            ("list", [{"s": "BTCUSDT"}]),
            ("empty", {}),
        ):
            with self.subTest(order=label):
                payload = {"e": "forceOrder", "E": event_ms}
                if order is not None:
                    payload["o"] = order
                self.assertEqual(
                    liquidation.force_order_schema_error(payload),
                    "order_object",
                )
                self.assertIsNone(
                    liquidation.canonical_force_order(
                        "!forceOrder@arr",
                        payload,
                    )
                )

    def test_latched_sink_terminates_instead_of_reconnecting(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = None

        class LatchingWebSocket:
            def __init__(self, owner):
                self.owner = owner

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                error = OSError("injected sink failure")
                self.owner._latch_sink_io("append", error)
                raise error

        class LatchingConnection:
            def __init__(self, owner):
                self.owner = owner

            async def __aenter__(self):
                return LatchingWebSocket(self.owner)

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        stderr = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    lambda _url, **_kwargs: LatchingConnection(collector),
                ),
                mock.patch.object(liquidation.sys, "stderr", stderr),
            ):
                rc = loop.run_until_complete(
                    asyncio.wait_for(collector.run(), timeout=5.0)
                )
        finally:
            loop.close()
        self.assertEqual(rc, 5)
        self.assertEqual(collector.stats["reconnects"], 1)
        self.assertIsNotNone(collector._sink_io_error)
        self.assertIn("COLLECTOR_INTERNAL_ERROR", stderr.getvalue())
        rows = [
            json.loads(line)
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertFalse(
            any(
                row.get("control", {}).get("type") == "session_stop"
                for row in rows
            )
        )

    def test_untrusted_ack_payloads_are_bounded_in_controls(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        limit = liquidation.AUDIT_TEXT_MAX_CHARS
        hostile = {"id": 1, "result": "x" * 900_000}

        class HostileAckWs:
            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                return json.dumps(hostile)

        loop = asyncio.new_event_loop()
        try:
            with self.assertRaises(RuntimeError) as raised:
                loop.run_until_complete(
                    collector.subscribe_and(
                        HostileAckWs(),
                        "c1",
                        ack_timeout=1.0,
                    )
                )
        finally:
            loop.close()
        self.assertLess(len(str(raised.exception)), 4 * limit)

        before_bytes = collector.actual_capture_bytes()
        for cycle in range(20):
            collector.emit_control(
                "disconnect",
                conn="c1",
                reason=liquidation.bounded_audit_text(
                    f"RuntimeError: {'x' * 900_000}"
                ),
                backoff_s=1.0,
                phase="pre-ack",
            )
        collector.fsync_close()
        rows = [
            json.loads(line)
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        reasons = [
            row["control"]["reason"]
            for row in rows
            if row.get("control", {}).get("type") == "disconnect"
        ]
        self.assertEqual(len(reasons), 20)
        self.assertEqual({len(reason) for reason in reasons}, {limit})
        self.assertLess(
            collector.actual_capture_bytes() - before_bytes,
            50_000,
        )

    def test_force_order_event_time_skew_is_enforced(self) -> None:
        limit = liquidation.FORCE_ORDER_MAX_SKEW_MS
        received_ms = 1_800_000_000_000
        for label, offset in (("future", limit), ("stale", -limit)):
            with self.subTest(boundary=label):
                collector, out = self.make_collector()
                self.addCleanup(collector.fsync_close)
                collector.open_for(liquidation.utc_date())
                frame = self.force_order_payload(received_ms + offset)
                with mock.patch.object(
                    collector,
                    "schedule_oi_snapshot",
                ) as oi:
                    collector.dispatch_market_message(
                        frame,
                        "c1",
                        received_ms,
                    )
                oi.assert_called_once_with("BTCUSDT")
                rows = [
                    json.loads(line)
                    for path in out.glob(
                        f"{liquidation.FILE_PREFIX}*.jsonl"
                    )
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                self.assertEqual(
                    sum(row.get("kind") == "event" for row in rows),
                    1,
                )
                collector.release_lease()

        for label, offset in (
            ("future", limit + 1),
            ("stale", -limit - 1),
        ):
            with self.subTest(violation=label):
                collector, out = self.make_collector()
                self.addCleanup(collector.fsync_close)
                collector.open_for(liquidation.utc_date())
                frame = self.force_order_payload(received_ms + offset)
                with (
                    mock.patch.object(
                        collector,
                        "schedule_oi_snapshot",
                    ) as oi,
                    self.assertRaises(liquidation.ForceOrderSchemaError),
                ):
                    collector.dispatch_market_message(
                        frame,
                        "c1",
                        received_ms,
                    )
                oi.assert_not_called()
                rows = [
                    json.loads(line)
                    for path in out.glob(
                        f"{liquidation.FILE_PREFIX}*.jsonl"
                    )
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                self.assertFalse(
                    any(row.get("kind") == "event" for row in rows)
                )
                self.assertEqual(collector.seen_ids, set())
                self.assertEqual(
                    [
                        row["control"]["reason"]
                        for row in rows
                        if row.get("control", {}).get("type")
                        == "force_order_schema_error"
                    ],
                    [f"source_timestamp_{label}"],
                )
                collector.release_lease()

    def test_force_order_publication_delay_is_bounded(self) -> None:
        limit = liquidation.FORCE_ORDER_MAX_PUBLICATION_DELAY_MS
        event_ms = 1_800_000_000_000
        for label, delay in (("immediate", 0), ("boundary", limit)):
            with self.subTest(accepted=label):
                payload = self.force_order_payload(event_ms)
                payload["o"]["T"] = event_ms - delay
                self.assertIsNone(
                    liquidation.force_order_schema_error(payload)
                )
        rejected = self.force_order_payload(event_ms)
        rejected["o"]["T"] = event_ms - limit - 1
        self.assertEqual(
            liquidation.force_order_schema_error(rejected),
            "order_publication_delay",
        )
        self.assertIsNone(
            liquidation.canonical_force_order("!forceOrder@arr", rejected)
        )

        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        with (
            mock.patch.object(collector, "schedule_oi_snapshot") as oi,
            self.assertRaises(liquidation.ForceOrderSchemaError),
        ):
            collector.dispatch_market_message(
                rejected,
                "c1",
                event_ms,
            )
        oi.assert_not_called()
        rows = [
            json.loads(line)
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertFalse(any(row.get("kind") == "event" for row in rows))
        self.assertEqual(collector.seen_ids, set())
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            ["order_publication_delay"],
        )

        replay, replay_out = self.make_collector()
        self.addCleanup(replay.fsync_close)
        replay_path = (
            replay_out
            / f"{liquidation.FILE_PREFIX}{liquidation.utc_date()}.jsonl"
        )
        replay_path.write_text(
            json.dumps(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "event",
                    "ts_utc": "2027-01-15T08:00:00.000Z",
                    "recv_epoch_ms": event_ms + 1,
                    "conn": "c1",
                    "subscription_phase": "confirmed",
                    "stream": "!forceOrder@arr",
                    "event_id": liquidation.event_id_for(
                        "!forceOrder@arr",
                        self.force_order_payload(event_ms),
                    ),
                    "source": {
                        "venue": "binance-usds-m-futures",
                        "endpoint_kind": (
                            "public-unauthenticated-websocket"
                        ),
                        "url": replay.url,
                    },
                    "payload": rejected,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        replay.seed_dedup_from_file(replay_path)
        self.assertEqual(
            replay.stats["dedup_rehydrate_semantic_rejected"],
            1,
        )
        self.assertEqual(replay.seen_ids, set())

    def test_replayed_force_order_receipt_is_validated(self) -> None:
        limit = liquidation.FORCE_ORDER_MAX_SKEW_MS
        received_ms = 1_800_000_000_000
        frame = self.force_order_payload(received_ms)
        event_id = liquidation.event_id_for("!forceOrder@arr", frame)
        cases = {
            "missing_receipt": None,
            "zero_receipt": 0,
            "float_receipt": float(received_ms),
            "skewed_receipt": received_ms + limit + 1,
        }
        for label, receipt in cases.items():
            with self.subTest(receipt=label):
                replay, out = self.make_collector()
                self.addCleanup(replay.fsync_close)
                record = {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "event",
                    "ts_utc": "2027-01-15T08:00:00.000Z",
                    "conn": "c1",
                    "subscription_phase": "confirmed",
                    "stream": "!forceOrder@arr",
                    "event_id": event_id,
                    "source": {
                        "venue": "binance-usds-m-futures",
                        "endpoint_kind": (
                            "public-unauthenticated-websocket"
                        ),
                        "url": replay.url,
                    },
                    "payload": frame,
                }
                if receipt is not None:
                    record["recv_epoch_ms"] = receipt
                path = (
                    out
                    / f"{liquidation.FILE_PREFIX}"
                    f"{liquidation.utc_date()}.jsonl"
                )
                path.write_text(
                    json.dumps(record, separators=(",", ":")) + "\n",
                    encoding="utf-8",
                )
                replay.seed_dedup_from_file(path)
                self.assertEqual(
                    replay.stats["dedup_rehydrate_semantic_rejected"],
                    1,
                )
                self.assertEqual(replay.seen_ids, set())
                replay.release_lease()

    def test_duplicate_json_member_is_rejected_at_every_depth(self) -> None:
        for payload in (
            '{"id":999,"id":1,"result":null}',
            '{"o":{"s":"BTCUSDT","s":"ETHUSDT"}}',
            '[{"x":1,"x":2}]',
        ):
            with self.subTest(payload=payload):
                with self.assertRaises(liquidation.DuplicateJsonKeyError):
                    liquidation.strict_json_loads(payload)
        self.assertEqual(
            liquidation.strict_json_loads('{"id":1,"result":null}'),
            {"id": 1, "result": None},
        )

    def test_duplicate_id_ack_never_reaches_ready(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())

        class DuplicateAckWebSocket:
            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                return '{"id":999,"id":1,"result":null}'

        loop = asyncio.new_event_loop()
        try:
            with self.assertRaises(liquidation.ForceOrderSchemaError):
                loop.run_until_complete(
                    collector.subscribe_and(
                        DuplicateAckWebSocket(),
                        "c1",
                        ack_timeout=1.0,
                    )
                )
        finally:
            loop.close()
        self.assertEqual(collector.stats["subscription_acks"], 0)
        controls = [
            row["control"]
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
            if row.get("kind") == "control"
        ]
        self.assertEqual(
            [
                control["reason"]
                for control in controls
                if control["type"] == "duplicate_json_frame"
            ],
            ["duplicate_json_member"],
        )

    def test_duplicate_member_market_frame_forces_reconnect(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        valid_frame = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )

        class DuplicateWebSocket:
            def __init__(self, attempt: int):
                self.attempt = attempt
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                if self.attempt == 1:
                    return (
                        '{"e":"forceOrder","E":1800000000000,'
                        '"o":{"s":"BTCUSDT","s":"ETHUSDT"}}'
                    )
                collector.request_stop()
                return json.dumps(valid_frame)

        class DuplicateConnection:
            def __init__(self, attempt: int):
                self.attempt = attempt

            async def __aenter__(self):
                return DuplicateWebSocket(self.attempt)

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        attempts = iter((1, 2))
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    lambda _url, **_kwargs: DuplicateConnection(
                        next(attempts)
                    ),
                ),
                mock.patch.object(
                    liquidation.random,
                    "uniform",
                    return_value=0,
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        controls = [
            row["control"] for row in rows if row.get("kind") == "control"
        ]
        self.assertEqual(rc, 0)
        self.assertEqual(collector.stats["reconnects"], 1)
        self.assertEqual(
            [
                control["reason"]
                for control in controls
                if control["type"] == "duplicate_json_frame"
            ],
            ["duplicate_json_member"],
        )
        self.assertEqual(
            [
                control.get("phase")
                for control in controls
                if control["type"] == "disconnect"
            ],
            ["post-ready", None],
        )
        self.assertEqual(
            sum(row.get("kind") == "event" for row in rows),
            1,
        )

    def test_duplicate_member_state_file_fails_closed(self) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-duplicate-state-"
        ))
        (out / liquidation.STATE_FILE).write_text(
            '{"schema_version":1,"schema_version":1,'
            '"total_bytes":0,"day_bytes":{}}',
            encoding="utf-8",
        )
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        self.addCleanup(collector.release_lease)
        with mock.patch.object(liquidation.sys, "stderr", io.StringIO()):
            collector.load_state()
        self.assertGreater(
            collector.state_total_bytes,
            collector.max_total_bytes,
        )
        self.assertEqual(collector.day_bytes, {})

    @staticmethod
    def oi_http_response(row: dict) -> FakeHttpResponse:
        return FakeHttpResponse([json.dumps([row]).encode("utf-8")])

    def test_force_order_rejects_transaction_after_event_time(self) -> None:
        event_ms = 1_800_000_000_000
        future = self.force_order_payload(event_ms)
        future["o"]["T"] = event_ms + 1
        self.assertEqual(
            liquidation.force_order_schema_error(future),
            "order_timestamp_after_event",
        )
        self.assertIsNone(
            liquidation.canonical_force_order("!forceOrder@arr", future)
        )
        simultaneous = self.force_order_payload(event_ms)
        simultaneous["o"]["T"] = event_ms
        self.assertIsNone(
            liquidation.force_order_schema_error(simultaneous)
        )
        self.assertIsNone(
            liquidation.force_order_schema_error(
                self.force_order_payload(event_ms)
            )
        )

        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        with (
            mock.patch.object(collector, "schedule_oi_snapshot") as oi,
            self.assertRaises(liquidation.ForceOrderSchemaError),
        ):
            collector.dispatch_market_message(future, "c1", event_ms + 5)
        oi.assert_not_called()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(any(row.get("kind") == "event" for row in rows))
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            ["order_timestamp_after_event"],
        )

        replay, replay_out = self.make_collector()
        self.addCleanup(replay.fsync_close)
        path = (
            replay_out
            / f"{liquidation.FILE_PREFIX}{liquidation.utc_date()}.jsonl"
        )
        path.write_text(
            json.dumps(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "event",
                    "ts_utc": "2027-01-15T08:00:00.000Z",
                    "recv_epoch_ms": event_ms + 5,
                    "conn": "c1",
                    "subscription_phase": "confirmed",
                    "stream": "!forceOrder@arr",
                    "event_id": liquidation.event_id_for(
                        "!forceOrder@arr",
                        self.force_order_payload(event_ms),
                    ),
                    "source": {
                        "venue": "binance-usds-m-futures",
                        "endpoint_kind": (
                            "public-unauthenticated-websocket"
                        ),
                        "url": replay.url,
                    },
                    "payload": future,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        replay.seed_dedup_from_file(path)
        self.assertEqual(
            replay.stats["dedup_rehydrate_semantic_rejected"],
            1,
        )
        self.assertEqual(replay.seen_ids, set())

    def test_oi_fetch_enforces_total_deadline_under_trickle(self) -> None:
        response = FakeHttpResponse(
            [b"["] + [b" "] * 120,
            delay=0.005,
        )
        started = time.monotonic()
        with mock.patch(
            "urllib.request.urlopen",
            return_value=response,
        ):
            snapshot, failure = liquidation.fetch_open_interest_snapshot(
                "BTCUSDT",
                timeout=0.1,
            )
        elapsed = time.monotonic() - started
        self.assertIsNone(snapshot)
        self.assertEqual(
            failure,
            {"reason": "timeout", "scope": "network"},
        )
        self.assertLess(elapsed, 0.4)
        self.assertTrue(response.chunks)

    def test_oi_fetch_bounds_response_size(self) -> None:
        chunk = b"x" * (64 * 1024)
        chunks = [chunk] * (
            liquidation.OI_MAX_RESPONSE_BYTES // len(chunk) + 2
        )
        response = FakeHttpResponse(list(chunks))
        with mock.patch(
            "urllib.request.urlopen",
            return_value=response,
        ):
            snapshot, failure = liquidation.fetch_open_interest_snapshot(
                "BTCUSDT",
            )
        self.assertIsNone(snapshot)
        self.assertEqual(
            failure,
            {"reason": "response-too-large", "scope": "endpoint"},
        )
        self.assertTrue(response.chunks)

    def test_inflight_oi_fetch_is_quiesced_before_run_returns(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 0.1
        liquidation_frame = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )
        active = {"running": 0, "started": 0}

        class OIWebSocket:
            def __init__(self):
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                if self.calls == 2:
                    return json.dumps(liquidation_frame)
                await asyncio.Future()
                raise AssertionError("unreachable")

        class OIConnection:
            async def __aenter__(self):
                return OIWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        def slow_fetch(_symbol: str):
            active["running"] += 1
            active["started"] += 1
            try:
                time.sleep(0.3)
                return {"symbol": "BTCUSDT", "row": {}}, None
            finally:
                active["running"] -= 1

        started = time.monotonic()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=OIConnection(),
                ),
                mock.patch.object(
                    collector,
                    "_oi_fetch_sync",
                    side_effect=slow_fetch,
                ),
            ):
                rc = loop.run_until_complete(
                    asyncio.wait_for(collector.run(), timeout=3.0)
                )
        finally:
            loop.close()
        elapsed = time.monotonic() - started
        self.assertEqual(rc, 0)
        self.assertEqual(active["started"], 1)
        self.assertEqual(active["running"], 0)
        self.assertIsNone(collector._oi_executor)
        self.assertLess(elapsed, 2.0)
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(
            any(row.get("kind") == "oi_snapshot" for row in rows)
        )
        stops = [
            row
            for row in rows
            if row.get("kind") == "control"
            and row.get("control", {}).get("type") == "session_stop"
        ]
        self.assertEqual(len(stops), 1)


    def test_semantically_corrupt_liquidation_state_fails_closed(self) -> None:
        bad_states = (
            [],
            {
                "schema_version": liquidation.SCHEMA_VERSION,
                "total_bytes": None,
                "day_bytes": {},
            },
            {
                "schema_version": liquidation.SCHEMA_VERSION,
                "total_bytes": 0,
                "day_bytes": [],
            },
            {
                "schema_version": liquidation.SCHEMA_VERSION + 1,
                "total_bytes": 0,
                "day_bytes": {},
            },
            {
                "schema_version": liquidation.SCHEMA_VERSION,
                "total_bytes": 0,
                "day_bytes": {"20260904": 1.5},
            },
        )
        for index, bad_state in enumerate(bad_states):
            with self.subTest(index=index):
                out = Path(tempfile.mkdtemp(prefix="collector-state-schema-"))
                (out / liquidation.STATE_FILE).write_text(
                    json.dumps(bad_state),
                    encoding="utf-8",
                )
                collector = liquidation.Collector(
                    out_dir=out,
                    url="wss://example.invalid/ws",
                    streams=["!forceOrder@arr"],
                    max_seconds=None,
                    max_bytes_per_day=10_000_000,
                    max_total_bytes=100_000_000,
                )
                with mock.patch.object(
                    liquidation.sys,
                    "stderr",
                    io.StringIO(),
                ):
                    collector.load_state()
                self.assertGreater(
                    collector.state_total_bytes,
                    collector.max_total_bytes,
                )
                self.assertEqual(collector.day_bytes, {})
                self.assertIsNone(collector.file)

    def test_liquidation_stream_config_requires_force_order(self) -> None:
        for streams in (
            [],
            ["btcusdt@aggTrade"],
            ["!forceOrder@arr", "btcusdt@forceOrder"],
        ):
            with self.subTest(streams=streams), self.assertRaises(ValueError):
                liquidation.Collector(
                    out_dir=Path(tempfile.mkdtemp(prefix="collector-stream-config-")),
                    url="wss://example.invalid/ws",
                    streams=streams,
                    max_seconds=None,
                    max_bytes_per_day=10_000_000,
                    max_total_bytes=100_000_000,
                )
        collector = liquidation.Collector(
            out_dir=Path(tempfile.mkdtemp(prefix="collector-stream-config-")),
            url="wss://example.invalid/ws",
            streams=[
                "btcusdt@aggTrade",
                "btcusdt@forceOrder",
                "ethusdt@forceOrder",
            ],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        self.addCleanup(collector.fsync_close)
        self.assertEqual(
            collector.streams,
            [
                "btcusdt@aggTrade",
                "btcusdt@forceOrder",
                "ethusdt@forceOrder",
            ],
        )

    def test_liquidation_stream_config_rejects_colliding_identities(
        self,
    ) -> None:
        for streams in (
            ["btcusdt@forceOrder", "BTCUSDT@forceOrder"],
            ["btcusdt@forceOrder", "btcusdt@forceOrder"],
            [
                "BTCUSDT@forceOrder",
                "ethusdt@forceOrder",
                "btcusdt@forceOrder",
            ],
        ):
            with self.subTest(streams=streams):
                self.assertIsNotNone(
                    liquidation.stream_config_error(streams)
                )
                out = Path(tempfile.mkdtemp(
                    prefix="collector-stream-collision-"
                ))
                with self.assertRaises(ValueError):
                    liquidation.Collector(
                        out_dir=out,
                        url="wss://example.invalid/ws",
                        streams=streams,
                        max_seconds=None,
                        max_bytes_per_day=10_000_000,
                        max_total_bytes=100_000_000,
                    )
                self.assertEqual(list(out.iterdir()), [])
        allowed = liquidation.Collector(
            out_dir=Path(tempfile.mkdtemp(
                prefix="collector-stream-distinct-"
            )),
            url="wss://example.invalid/ws",
            streams=["BTCUSDT@forceOrder", "ethusdt@forceOrder"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        self.addCleanup(allowed.release_lease)
        self.assertEqual(
            allowed.streams,
            ["BTCUSDT@forceOrder", "ethusdt@forceOrder"],
        )

    def test_liquidation_backward_day_change_rotates_and_charges_day(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(
            prefix="collector-integrity-liquidation-backward-day-"
        ))
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        collector.load_state()
        self.addCleanup(collector.release_lease)
        self.addCleanup(collector.fsync_close)
        later = ("2026-09-04T00:00:01.000Z", 1_788_480_001.0)
        earlier = ("2026-09-03T23:59:59.000Z", 1_788_479_999.0)
        for index, clock in enumerate((later, earlier, later)):
            collector.write_line(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "control",
                    "control": {"type": "probe", "seq": index},
                },
                write_clock=clock,
            )
        collector.fsync_close()

        for day in ("20260903", "20260904"):
            path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
            rows = [
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            ]
            self.assertTrue(rows)
            self.assertEqual(
                {row["ts_utc"][:10].replace("-", "") for row in rows},
                {day},
            )
            self.assertEqual(collector.day_bytes[day], path.stat().st_size)
        forward_rows = [
            json.loads(line)
            for line in (
                out / f"{liquidation.FILE_PREFIX}20260904.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            [
                row["control"]["seq"]
                for row in forward_rows
                if row["control"]["type"] == "probe"
            ],
            [0, 2],
        )
        self.assertEqual(
            [
                row["control"]["to_day"]
                for row in forward_rows + [
                    json.loads(line)
                    for line in (
                        out / f"{liquidation.FILE_PREFIX}20260903.jsonl"
                    ).read_text(encoding="utf-8").splitlines()
                ]
                if row["control"]["type"] == "rotate"
            ],
            ["20260904", "20260903"],
        )

        collector.max_bytes_per_day = collector.day_bytes["20260903"] + 1
        with self.assertRaises(liquidation.CapExceeded):
            collector.write_line(
                {
                    "v": liquidation.SCHEMA_VERSION,
                    "kind": "control",
                    "control": {"type": "probe", "seq": 3},
                },
                write_clock=earlier,
            )

    def test_liquidation_reopen_isolates_partial_jsonl_tail(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-liquidation-tail-"))
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        partial = b'{"kind":"event"'
        path.write_bytes(partial)
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        collector.load_state()
        self.addCleanup(collector.fsync_close)
        collector.open_for(day)
        collector.emit_control("after_partial_tail")
        collector.fsync_close()
        collector.save_state()

        raw = path.read_bytes()
        self.assertTrue(raw.startswith(partial + b"\n"))
        lines = raw.splitlines()
        self.assertEqual(lines[0], partial)
        reopened = json.loads(lines[1])
        self.assertEqual(reopened["control"]["type"], "file_reopened")
        self.assertTrue(
            reopened["control"]["partial_tail_isolated"]
        )
        self.assertEqual(
            json.loads(lines[-1])["control"]["type"],
            "after_partial_tail",
        )
        state = json.loads(
            (out / liquidation.STATE_FILE).read_text(encoding="utf-8")
        )
        self.assertEqual(state["total_bytes"], len(raw))
        self.assertEqual(state["day_bytes"][day], len(raw))

    def test_liquidation_dedup_rehydrate_reads_only_bounded_tail(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.dedup_max_ids = 2
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        prefix = json.dumps({
            "kind": "control",
            "padding": "x" * 131_072,
        }) + "\n"
        tail_payloads = [
            self.force_order_payload(
                int(liquidation.time.time() * 1000) - 5_000 + index
            )
            for index in range(3)
        ]
        tail_ids = [
            liquidation.event_id_for("!forceOrder@arr", payload)
            for payload in tail_payloads
        ]
        tail = [
            json.dumps({
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "recv_epoch_ms": payload["E"] + 1,
                "stream": "!forceOrder@arr",
                "event_id": event_id,
                "payload": payload,
            })
            for event_id, payload in zip(tail_ids, tail_payloads)
        ]
        path.write_text(
            prefix + "\n".join(tail) + "\n",
            encoding="utf-8",
        )
        collector.seed_dedup_from_file(path)
        self.assertEqual(collector.seen_ids, set(tail_ids[-2:]))
        self.assertEqual(
            list(collector._seen_id_order),
            tail_ids[-2:],
        )
        self.assertEqual(collector.stats["dedup_rehydrate_ids"], 2)
        self.assertLess(collector.stats["dedup_rehydrate_bytes"], 100_000)

    def test_liquidation_rehydrate_recomputes_canonical_event_id(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        stream = "!forceOrder@arr"

        corrected_schema_payload = self.force_order_payload(
            1_800_000_000_100
        )
        corrected_schema_id = liquidation.event_id_for(
            stream,
            corrected_schema_payload,
        )
        invalid_payload = self.force_order_payload(
            1_800_000_000_100
        )
        del invalid_payload["o"]["q"]

        stored_payload = self.force_order_payload(
            1_800_000_000_200
        )
        corrected_identity_payload = self.force_order_payload(
            1_800_000_000_201
        )
        corrected_identity_id = liquidation.event_id_for(
            stream,
            corrected_identity_payload,
        )
        unconfigured_payload = self.force_order_payload(
            1_800_000_000_300
        )
        unconfigured_stream = "btcusdt@forceOrder"
        unconfigured_id = liquidation.event_id_for(
            unconfigured_stream,
            unconfigured_payload,
        )
        wrong_version_payload = self.force_order_payload(
            1_800_000_000_400
        )
        wrong_version_id = liquidation.event_id_for(
            stream,
            wrong_version_payload,
        )
        poisoned_records = [
            {
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "stream": stream,
                "event_id": corrected_schema_id,
                "payload": invalid_payload,
            },
            {
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "stream": stream,
                "event_id": corrected_identity_id,
                "payload": stored_payload,
            },
            {
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "stream": unconfigured_stream,
                "event_id": unconfigured_id,
                "payload": unconfigured_payload,
            },
            {
                "v": liquidation.SCHEMA_VERSION + 1,
                "kind": "event",
                "stream": stream,
                "event_id": wrong_version_id,
                "payload": wrong_version_payload,
            },
        ]
        path.write_text(
            "".join(
                liquidation.strict_json_dumps(record) + "\n"
                for record in poisoned_records
            ),
            encoding="utf-8",
        )

        collector.open_for(day)

        self.assertEqual(collector.seen_ids, set())
        self.assertEqual(
            collector.stats["dedup_rehydrate_semantic_rejected"],
            len(poisoned_records),
        )
        self.assertEqual(
            collector.stats["dedup_rehydrate_malformed"],
            0,
        )
        collector.emit_event(
            stream,
            corrected_schema_payload,
            received_ms=1_800_000_000_500,
            connection_id="c1",
        )
        collector.emit_event(
            stream,
            corrected_identity_payload,
            received_ms=1_800_000_000_501,
            connection_id="c1",
        )
        self.assertEqual(
            collector.seen_ids,
            {corrected_schema_id, corrected_identity_id},
        )
        self.assertEqual(collector.stats["events_written"], 2)
        self.assertEqual(collector.stats["events_deduped"], 0)

        collector.emit_event(stream, corrected_schema_payload)
        collector.emit_event(stream, corrected_identity_payload)
        self.assertEqual(collector.stats["events_written"], 2)
        self.assertEqual(collector.stats["events_deduped"], 2)

    def test_liquidation_live_dedup_uses_same_fifo_bound(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.dedup_max_ids = 2
        event_ids = []
        for index in range(3):
            payload = self.force_order_payload(
                1_800_000_000_000 + index
            )
            event_ids.append(
                liquidation.event_id_for("!forceOrder@arr", payload)
            )
            collector.emit_event("!forceOrder@arr", payload)
        self.assertEqual(collector.seen_ids, set(event_ids[-2:]))
        self.assertEqual(list(collector._seen_id_order), event_ids[-2:])
        self.assertEqual(collector.stats["dedup_ids_evicted"], 1)

    def test_expired_deadline_prevents_liquidation_session_start(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-liquidation-deadline-"))
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        path.write_text(
            json.dumps({
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "event_id": "existing-event",
            }) + "\n",
            encoding="utf-8",
        )
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=0,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        connect = mock.Mock(
            side_effect=AssertionError(
                "websocket entered after startup deadline"
            )
        )
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(liquidation.websockets, "connect", connect),
                mock.patch.object(liquidation.sys, "stdout", io.StringIO()),
                mock.patch.object(liquidation.sys, "stderr", io.StringIO()),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        rows = [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(rc, 0)
        connect.assert_not_called()
        self.assertFalse(any(
            row.get("control", {}).get("type") == "session_start"
            for row in rows
        ))

    def test_liquidation_output_lease_rejects_real_process_contender(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-liquidation-lease-"))
        day = liquidation.utc_date()
        capture_path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        original = b'{"kind":"control"}\n'
        capture_path.write_bytes(original)
        lock_path = out / liquidation.OUTPUT_LEASE_FILE
        holder = start_output_lock_holder(lock_path)
        self.addCleanup(stop_process, holder)
        result = subprocess.run(
            [
                sys.executable,
                str(HERE / "focused-liquidation-overlay" / "capture.py"),
                "--out",
                str(out),
                "--max-seconds",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=10,
            env={
                **liquidation.os.environ,
                "PYTHONDONTWRITEBYTECODE": "1",
            },
        )
        self.assertEqual(result.returncode, 6)
        self.assertIn("LIQUIDATION_CAPTURE_LOCKED", result.stderr)
        self.assertNotIn(liquidation.READY_BANNER, result.stdout)
        self.assertEqual(capture_path.read_bytes(), original)
        self.assertFalse((out / liquidation.STATE_FILE).exists())

        stop_process(holder)
        successor = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=0,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        try:
            successor.load_state()
            self.assertEqual(
                successor.state_total_bytes,
                len(original),
            )
            self.assertEqual(
                successor.day_bytes[day],
                len(original),
            )
        finally:
            successor.release_lease()
        self.assertTrue(lock_path.is_file())

    def test_liquidation_accounting_stat_failure_is_fail_closed(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-liquidation-stat-"))
        day = liquidation.utc_date()
        capture_path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        original = b'{"kind":"control"}\n'
        capture_path.write_bytes(original)
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        original_lstat = liquidation.os.lstat

        def guarded_lstat(target, *args, **kwargs):
            if Path(target) == capture_path:
                raise PermissionError("injected capture stat failure")
            return original_lstat(target, *args, **kwargs)

        with (
            mock.patch.object(liquidation.os, "lstat", guarded_lstat),
            self.assertRaises(PermissionError),
        ):
            collector.load_state()
        self.assertIsNone(collector.file)
        self.assertEqual(capture_path.read_bytes(), original)

    def test_special_capture_paths_are_rejected(self) -> None:
        day = liquidation.utc_date()
        outside = Path(tempfile.mkdtemp(prefix="collector-integrity-outside-"))
        external = outside / "external-target.jsonl"
        external_bytes = b'{"kind":"external"}\n'

        for label in ("symlink", "fifo"):
            with self.subTest(capture=label):
                external.write_bytes(external_bytes)
                out = Path(tempfile.mkdtemp(
                    prefix="collector-integrity-special-capture-"
                ))
                capture_path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
                if label == "symlink":
                    capture_path.symlink_to(external)
                else:
                    os.mkfifo(capture_path)
                collector = liquidation.Collector(
                    out_dir=out,
                    url="wss://example.invalid/ws",
                    streams=["!forceOrder@arr"],
                    max_seconds=None,
                    max_bytes_per_day=10_000_000,
                    max_total_bytes=100_000_000,
                )
                self.addCleanup(collector.release_lease)
                with self.assertRaises(liquidation.SinkIOError):
                    collector.open_for(day)
                self.assertIsNone(collector.file)
                self.assertIsNotNone(collector._sink_io_error)
                self.assertEqual(external.read_bytes(), external_bytes)

    def test_state_temp_symlink_is_replaced_not_followed(self) -> None:
        outside = Path(tempfile.mkdtemp(prefix="collector-integrity-outside-"))
        external = outside / "external-state.json"
        external_bytes = b'{"external":true}\n'
        external.write_bytes(external_bytes)
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        temp_path = (out / liquidation.STATE_FILE).with_suffix(".json.tmp")
        temp_path.symlink_to(external)
        collector.save_state()
        self.assertEqual(external.read_bytes(), external_bytes)
        state = json.loads(
            (out / liquidation.STATE_FILE).read_text(encoding="utf-8")
        )
        self.assertEqual(state["schema_version"], liquidation.SCHEMA_VERSION)

    def test_liquidation_durable_event_rehydrates_without_state(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector._state_save_due_monotonic = float("inf")
        payload = self.force_order_payload()
        event_id = liquidation.event_id_for(
            "!forceOrder@arr",
            payload,
        )
        with mock.patch.object(collector, "maybe_save_state"):
            collector.emit_event(
                "!forceOrder@arr",
                payload,
                received_ms=1_800_000_000_500,
                connection_id="c1",
            )
        self.assertFalse((out / liquidation.STATE_FILE).exists())
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        collector.fsync_close()
        collector.release_lease()

        successor = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        self.addCleanup(successor.release_lease)
        successor.load_state()
        successor.seed_dedup_from_file(path)
        self.assertIn(event_id, successor.seen_ids)

    def test_liquidation_new_file_and_state_fsync_directory_in_order(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector._state_save_due_monotonic = float("inf")
        operations = []

        def file_fsync(fd: int) -> None:
            operations.append(("file", fd))

        def directory_fsync(path: Path) -> None:
            operations.append(("directory", Path(path)))

        day = liquidation.utc_date()
        with (
            mock.patch.object(
                liquidation.os,
                "fsync",
                side_effect=file_fsync,
            ),
            mock.patch.object(
                liquidation,
                "fsync_directory",
                side_effect=directory_fsync,
            ),
        ):
            collector.open_for(day)
        self.assertEqual(
            operations,
            [
                ("file", collector.file.fileno()),
                ("directory", out),
            ],
        )

        operations.clear()
        with (
            mock.patch.object(
                liquidation.os,
                "fsync",
                side_effect=file_fsync,
            ),
            mock.patch.object(
                liquidation,
                "fsync_directory",
                side_effect=directory_fsync,
            ),
        ):
            collector.emit_control("existing-file-append")
        self.assertEqual(
            operations,
            [("file", collector.file.fileno())],
        )

        state_operations = []
        real_replace = liquidation.os.replace

        def replace(source, destination) -> None:
            state_operations.append("replace")
            real_replace(source, destination)

        with (
            mock.patch.object(
                liquidation.os,
                "replace",
                side_effect=replace,
            ),
            mock.patch.object(
                liquidation,
                "fsync_directory",
                side_effect=lambda path: state_operations.append(
                    ("directory", Path(path))
                ),
            ),
        ):
            collector.save_state()
        self.assertEqual(
            state_operations,
            ["replace", ("directory", out)],
        )

    def test_liquidation_directory_fsync_failure_precedes_dedup(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector._state_save_due_monotonic = float("inf")
        payload = self.force_order_payload()
        event_id = liquidation.event_id_for(
            "!forceOrder@arr",
            payload,
        )
        with (
            mock.patch.object(
                liquidation,
                "fsync_directory",
                side_effect=OSError("injected directory fsync failure"),
            ),
            self.assertRaises(liquidation.SinkIOError),
        ):
            collector.dispatch_market_message(
                payload,
                "c1",
                1_800_000_000_500,
            )
        self.assertNotIn(event_id, collector.seen_ids)
        self.assertIsNotNone(collector._sink_io_error)

    def test_liquidation_directory_fsync_failure_never_announces_ready(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        collector.max_seconds = 1.0
        collector._state_save_due_monotonic = float("inf")
        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation,
                    "fsync_directory",
                    side_effect=OSError(
                        "injected directory fsync failure"
                    ),
                ),
                mock.patch.object(
                    liquidation.sys,
                    "stdout",
                    stdout,
                ),
                mock.patch.object(
                    liquidation.sys,
                    "stderr",
                    io.StringIO(),
                ),
            ):
                rc = loop.run_until_complete(
                    asyncio.wait_for(
                        collector.run(),
                        timeout=1.0,
                    )
                )
        finally:
            loop.close()
        self.assertEqual(rc, 5)
        self.assertNotIn(
            liquidation.READY_BANNER,
            stdout.getvalue(),
        )
        self.assertEqual(collector.seen_ids, set())
        self.assertIsNotNone(collector._sink_io_error)

    def test_liquidation_fsync_failure_latches_sink_terminal(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        with (
            mock.patch.object(
                liquidation.os,
                "fsync",
                side_effect=OSError("injected fsync failure"),
            ),
            self.assertRaises(OSError),
        ):
            collector.emit_control("failed-write")
        frozen_size = collector.actual_capture_bytes()
        with self.assertRaises(liquidation.SinkIOError):
            collector.emit_control("must-not-write")
        self.assertEqual(collector.actual_capture_bytes(), frozen_size)
        collector.fsync_close()
        collector.save_state()
        state = json.loads(
            (out / liquidation.STATE_FILE).read_text(encoding="utf-8")
        )
        self.assertEqual(state["total_bytes"], frozen_size)

    def test_liquidation_dedup_commits_only_after_durable_write(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        payload = self.force_order_payload(1_800_000_000_000)
        event_id = liquidation.event_id_for("!forceOrder@arr", payload)
        with (
            mock.patch.object(
                collector,
                "write_line",
                side_effect=OSError("injected write failure"),
            ),
            self.assertRaises(OSError),
        ):
            collector.emit_event("!forceOrder@arr", payload)
        self.assertNotIn(event_id, collector.seen_ids)

        collector.emit_event("!forceOrder@arr", payload)
        collector.emit_event("!forceOrder@arr", payload)
        events = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            events.extend(
                row
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "event"
            )
        self.assertEqual(len(events), 1)

    def test_first_cap_trip_latches_against_smaller_later_writes(self) -> None:
        collector, _ = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        collector.max_bytes_per_day = collector.actual_capture_bytes() + 50

        with self.assertRaises(liquidation.CapExceeded) as first:
            collector.emit_control("oversized")
        frozen_bytes = collector.actual_capture_bytes()
        with self.assertRaises(liquidation.CapExceeded) as second:
            collector.emit_control("would-otherwise-fit")
        self.assertEqual(str(second.exception), str(first.exception))
        self.assertEqual(collector.actual_capture_bytes(), frozen_bytes)

    def test_active_writes_refresh_exact_state_when_due(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        collector._state_save_due_monotonic = 0.0
        collector.emit_control("periodic-state-regression")

        state_path = out / liquidation.STATE_FILE
        self.assertTrue(state_path.is_file())
        state = json.loads(state_path.read_text(encoding="utf-8"))
        files = list(out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"))
        on_disk_total = sum(path.stat().st_size for path in files)
        self.assertEqual(state["total_bytes"], on_disk_total)
        self.assertEqual(
            state["day_bytes"][liquidation.utc_date()],
            files[0].stat().st_size,
        )

    def test_byte_threshold_refreshes_state_before_time_deadline(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        collector.save_state()
        before = json.loads((out / liquidation.STATE_FILE).read_text(encoding="utf-8"))
        collector._state_save_due_monotonic = float("inf")
        collector._bytes_at_last_state_save = collector.total_bytes
        with mock.patch.object(
            liquidation,
            "STATE_SAVE_BYTE_INTERVAL",
            8,
        ):
            collector.emit_control("state-byte-threshold")
        after = json.loads((out / liquidation.STATE_FILE).read_text(encoding="utf-8"))
        self.assertGreater(after["total_bytes"], before["total_bytes"])
        self.assertEqual(after["total_bytes"], collector.actual_capture_bytes())

    def test_unexpected_child_task_failure_is_nonzero_and_auditable(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1e-9

        async def broken_worker() -> None:
            raise RuntimeError("injected worker failure")

        collector.oi_worker = broken_worker
        loop = asyncio.new_event_loop()
        try:
            rc = loop.run_until_complete(
                asyncio.wait_for(collector.run(), timeout=1.0)
            )
        finally:
            loop.close()
        self.assertEqual(rc, 5)
        reasons = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            reasons.extend(
                row["control"]["reason"]
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "control"
                and row.get("control", {}).get("type") == "session_stop"
            )
        self.assertEqual(reasons, ["internal-error"])



    def test_max_seconds_bounds_subscription_ack_wait(self) -> None:
        collector, _ = self.make_collector()
        collector.max_seconds = 0.05

        class NoAckWebSocket:
            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                await asyncio.Future()
                raise AssertionError("unreachable")

        class NoAckConnection:
            async def __aenter__(self):
                return NoAckWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        started = liquidation.time.monotonic()
        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(
                liquidation.websockets,
                "connect",
                return_value=NoAckConnection(),
            ):
                rc = loop.run_until_complete(
                    asyncio.wait_for(collector.run(), timeout=0.5)
                )
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        self.assertLess(liquidation.time.monotonic() - started, 0.5)
        self.assertEqual(collector.stats["subscription_acks"], 0)

    def test_bounded_run_records_max_seconds_stop_reason(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1e-9
        loop = asyncio.new_event_loop()
        try:
            rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        stop_records = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            stop_records.extend(
                row
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "control"
                and row.get("control", {}).get("type") == "session_stop"
            )
        self.assertEqual(
            [row["control"]["reason"] for row in stop_records],
            ["max-seconds"],
        )

    def test_signal_stop_records_exact_signal_reason(self) -> None:
        for signum, expected in (
            (signal.SIGTERM, "sigterm"),
            (signal.SIGINT, "sigint"),
        ):
            with self.subTest(signal=expected):
                collector, out = self.make_collector()
                collector.request_stop(signum)
                loop = asyncio.new_event_loop()
                try:
                    rc = loop.run_until_complete(collector.run())
                finally:
                    loop.close()
                self.assertEqual(rc, 0)
                reasons = []
                for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
                    reasons.extend(
                        row["control"]["reason"]
                        for row in (
                            json.loads(line)
                            for line in path.read_text(
                                encoding="utf-8"
                            ).splitlines()
                        )
                        if row.get("kind") == "control"
                        and row.get("control", {}).get("type")
                        == "session_stop"
                    )
                self.assertEqual(reasons, [expected])

    def test_websocket_keepalive_has_bounded_half_open_detection(self) -> None:
        collector, _ = self.make_collector()
        collector.max_seconds = 0.1
        observed: dict[str, object] = {}

        class FakeWebSocket:
            def __init__(self):
                self.receive_count = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.receive_count += 1
                if self.receive_count == 1:
                    return json.dumps({"id": 1, "result": None})
                await asyncio.sleep(1)
                raise AssertionError("bounded receive timeout did not fire")

        class FakeConnection:
            async def __aenter__(self):
                return FakeWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        def connect(url: str, **kwargs):
            observed["url"] = url
            observed.update(kwargs)
            return FakeConnection()

        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(liquidation.websockets, "connect", connect):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        self.assertEqual(collector.stats["subscription_acks"], 1)
        self.assertEqual(observed["ping_interval"], 20)
        self.assertEqual(observed["ping_timeout"], 20)

    def test_symbol_specific_raw_frame_keeps_stream_provenance(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.streams = ["btcusdt@forceOrder"]
        frame = self.force_order_payload(1_800_000_000_000)
        collector.dispatch_market_message(frame, "c1", 1_800_000_000_123)
        events = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            events.extend(
                row
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "event"
            )
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["stream"], "btcusdt@forceOrder")
        self.assertEqual(
            events[0]["event_id"],
            liquidation.event_id_for("btcusdt@forceOrder", frame),
        )

    def test_force_order_stream_rejects_non_force_event_type(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.streams = ["btcusdt@forceOrder"]
        collector.dispatch_market_message(
            {
                "stream": "btcusdt@forceOrder",
                "data": {
                    "e": "aggTrade",
                    "E": 1_800_000_000_000,
                    "s": "BTCUSDT",
                    "a": 1,
                },
            },
            "c1",
            1_800_000_000_001,
        )
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertEqual(
            len([row for row in rows if row.get("kind") == "event"]),
            0,
        )
        self.assertEqual(collector.stats["ignored_market_frames"], 1)
        self.assertEqual(
            len([
                row
                for row in rows
                if row.get("control", {}).get("type")
                == "ignored_non_force_order_stream"
            ]),
            1,
        )

    def test_invalid_force_orders_are_audited_without_events(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        event_ms = 1_800_000_000_000
        missing_identity = self.force_order_payload(event_ms)
        del missing_identity["o"]["q"]
        invalid_frames = [
            {"e": "forceOrder", "E": event_ms},
            {"e": "forceOrder", "E": event_ms, "o": {}},
            {"e": "forceOrder", "E": event_ms, "o": "bad"},
            missing_identity,
            {
                **self.force_order_payload(event_ms),
                "E": 0,
            },
            {
                **self.force_order_payload(event_ms),
                "o": {
                    **self.force_order_payload(event_ms)["o"],
                    "q": "NaN",
                },
            },
        ]
        for frame in invalid_frames:
            with self.assertRaises(liquidation.ForceOrderSchemaError):
                collector.dispatch_market_message(
                    frame,
                    "c1",
                    event_ms + 1,
                )
        with self.assertRaises(liquidation.ForceOrderSchemaError):
            collector.dispatch_market_message(
                {
                    "stream": "btcusdt@forceOrder",
                    "data": self.force_order_payload(event_ms),
                },
                "c1",
                event_ms + 1,
            )
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(
            any(row.get("kind") == "event" for row in rows)
        )
        schema_errors = [
            row
            for row in rows
            if row.get("control", {}).get("type")
            == "force_order_schema_error"
        ]
        self.assertEqual(
            len(schema_errors),
            len(invalid_frames) + 1,
        )
        self.assertTrue(
            all(
                len(json.dumps(row)) < 1_000
                for row in schema_errors
            )
        )

    def test_force_order_quantity_relationships_gate_event_and_oi(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        event_ms = 1_800_000_000_000
        invalid = [
            (
                "order_last_fill_exceeds_accumulated",
                {"q": "2.0", "z": "1e0", "l": "1.0001"},
            ),
            (
                "order_accumulated_fill_exceeds_quantity",
                {"q": "1e0", "z": "1.0001", "l": "0.5"},
            ),
        ]
        for index, (reason, quantities) in enumerate(invalid):
            frame = self.force_order_payload(event_ms + index)
            frame["o"].update(quantities)
            with self.assertRaises(
                liquidation.ForceOrderSchemaError
            ) as failure:
                collector.dispatch_market_message(
                    frame,
                    "c1",
                    event_ms + index + 100,
                )
            self.assertEqual(failure.exception.reason, reason)

        self.assertEqual(collector.seen_ids, set())
        self.assertEqual(collector._oi_queue.qsize(), 0)
        rows = [
            liquidation.strict_json_loads(line)
            for path in out.glob(
                f"{liquidation.FILE_PREFIX}*.jsonl"
            )
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertFalse(
            any(row.get("kind") == "event" for row in rows)
        )
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("kind") == "control"
                and row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            [reason for reason, _quantities in invalid],
        )

        accepted = [
            ("BTCUSDT", {"q": "1", "z": "1.0", "l": "1e0"}),
            ("ETHUSDT", {"q": "2.5", "z": "0", "l": "0.0"}),
            ("SOLUSDT", {"q": "3e0", "z": "2.0", "l": "0.25"}),
        ]
        for index, (symbol, quantities) in enumerate(accepted, start=10):
            frame = self.force_order_payload(
                event_ms + index,
                symbol=symbol,
            )
            frame["o"].update(quantities)
            collector.dispatch_market_message(
                frame,
                "c1",
                event_ms + index + 100,
            )
        self.assertEqual(len(collector.seen_ids), len(accepted))
        self.assertEqual(collector._oi_queue.qsize(), len(accepted))

    def test_impossible_force_order_is_not_trusted_on_rehydrate(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        stream = "!forceOrder@arr"
        payload = self.force_order_payload()
        payload["o"].update({"q": "1", "z": "2", "l": "1"})
        order = payload["o"]
        parts = [stream, str(payload["E"])]
        parts.extend(
            str(order[key])
            for key in (
                "s",
                "S",
                "o",
                "f",
                "q",
                "p",
                "ap",
                "X",
                "l",
                "z",
                "T",
            )
        )
        legacy_event_id = liquidation.hashlib.sha256(
            "|".join(parts).encode("utf-8")
        ).hexdigest()[:32]
        path.write_text(
            liquidation.strict_json_dumps({
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "stream": stream,
                "event_id": legacy_event_id,
                "payload": payload,
            })
            + "\n",
            encoding="utf-8",
        )

        collector.seed_dedup_from_file(path)

        self.assertNotIn(legacy_event_id, collector.seen_ids)
        self.assertEqual(
            collector.stats["dedup_rehydrate_semantic_rejected"],
            1,
        )
        self.assertEqual(collector.stats["dedup_rehydrate_ids"], 0)

    def test_impossible_force_order_disconnects_confirmed_session(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        impossible = self.force_order_payload()
        impossible["o"].update({"q": "1", "z": "2", "l": "1"})

        class ImpossibleOrderWebSocket:
            def __init__(self) -> None:
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                collector.request_stop()
                return json.dumps(impossible)

        class ImpossibleOrderConnection:
            async def __aenter__(self):
                return ImpossibleOrderWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=ImpossibleOrderConnection(),
                ),
                mock.patch.object(sys, "stdout", stdout),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()

        rows = [
            liquidation.strict_json_loads(line)
            for path in out.glob(
                f"{liquidation.FILE_PREFIX}*.jsonl"
            )
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        lifecycle = [
            row["control"]["type"]
            for row in rows
            if row.get("kind") == "control"
            and row.get("control", {}).get("type")
            in {"connection_attempt", "connect", "disconnect"}
        ]
        self.assertEqual(rc, 0)
        self.assertIn(liquidation.READY_BANNER, stdout.getvalue())
        self.assertEqual(
            lifecycle,
            ["connection_attempt", "connect", "disconnect"],
        )
        self.assertEqual(
            [
                row["control"]["reason"]
                for row in rows
                if row.get("kind") == "control"
                and row.get("control", {}).get("type")
                == "force_order_schema_error"
            ],
            ["order_accumulated_fill_exceeds_quantity"],
        )
        self.assertFalse(
            any(row.get("kind") == "event" for row in rows)
        )
        self.assertFalse(
            any(row.get("kind") == "oi_snapshot" for row in rows)
        )

    def test_valid_raw_and_wrapped_force_orders_share_identity(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        payload = self.force_order_payload()
        collector.dispatch_market_message(
            payload,
            "c1",
            payload["E"] + 1,
        )
        collector.dispatch_market_message(
            {
                "stream": "!forceOrder@arr",
                "data": payload,
            },
            "c1",
            payload["E"] + 2,
        )
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        events = [row for row in rows if row.get("kind") == "event"]
        self.assertEqual(len(events), 1)
        self.assertEqual(
            events[0]["event_id"],
            liquidation.event_id_for("!forceOrder@arr", payload),
        )

    def test_pre_ack_auxiliary_flood_does_not_buffer(self) -> None:
        collector, _ = self.make_collector()
        collector.streams = [
            "!forceOrder@arr",
            "btcusdt@aggTrade",
        ]
        frame_count = 4_097

        class AuxiliaryFloodWebSocket:
            def __init__(self):
                self.index = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                if self.index < frame_count:
                    self.index += 1
                    return json.dumps({
                        "stream": "btcusdt@aggTrade",
                        "data": {
                            "e": "aggTrade",
                            "E": 1_800_000_000_000 + self.index,
                            "s": "BTCUSDT",
                            "a": self.index,
                        },
                    })
                return json.dumps({"id": 1, "result": None})

        written = []
        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(
                collector,
                "write_line",
                side_effect=lambda record, **_kwargs: written.append(record),
            ):
                result = loop.run_until_complete(
                    collector.subscribe_and(
                        AuxiliaryFloodWebSocket(),
                        "c1",
                    )
                )
        finally:
            loop.close()
            collector.fsync_close()
        self.assertIsNone(result)
        self.assertEqual(
            collector.stats["ignored_market_frames"],
            frame_count,
        )
        self.assertFalse(
            any(row.get("kind") == "event" for row in written)
        )
        self.assertEqual(
            sum(
                row.get("control", {}).get("type")
                == "ignored_non_force_order_stream"
                for row in written
            ),
            1,
        )

    def test_pre_ack_force_order_flood_is_persisted_without_buffer(
        self,
    ) -> None:
        collector, _ = self.make_collector()
        frame_count = 4_097
        base_payload = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )

        class ForceOrderFloodWebSocket:
            def __init__(self):
                self.index = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                if self.index < frame_count:
                    event_ms = base_payload["E"] + self.index
                    self.index += 1
                    return json.dumps({
                        **base_payload,
                        "E": event_ms,
                        "o": {
                            **base_payload["o"],
                            "T": event_ms - 1,
                        },
                    })
                return json.dumps({"id": 1, "result": None})

        written = []
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    collector,
                    "write_line",
                    side_effect=(
                        lambda record, **_kwargs: written.append(record)
                    ),
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
            ):
                result = loop.run_until_complete(
                    collector.subscribe_and(
                        ForceOrderFloodWebSocket(),
                        "c1",
                    )
                )
        finally:
            loop.close()
            collector.fsync_close()
        events = [
            row for row in written if row.get("kind") == "event"
        ]
        self.assertIsNone(result)
        self.assertEqual(len(events), frame_count)
        self.assertEqual(
            len({row["event_id"] for row in events}),
            frame_count,
        )
        self.assertTrue(
            all(
                row["subscription_phase"] == "pre_ack"
                for row in events
            )
        )

    def test_ack_error_preserves_unconfirmed_event_and_attempt_identity(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        pre_ack_event = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )

        class AckErrorWebSocket:
            def __init__(self):
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps(pre_ack_event)
                collector.request_stop()
                return json.dumps({
                    "id": 1,
                    "error": {"code": 400},
                })

        class AckErrorConnection:
            async def __aenter__(self):
                return AckErrorWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=AckErrorConnection(),
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
                mock.patch.object(sys, "stdout", stdout),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        self.assertNotIn(liquidation.READY_BANNER, stdout.getvalue())
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        events = [row for row in rows if row.get("kind") == "event"]
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["subscription_phase"], "pre_ack")
        controls = [
            row
            for row in rows
            if row.get("kind") == "control"
        ]
        self.assertFalse(
            any(
                row["control"]["type"] == "connect"
                for row in controls
            )
        )
        attempts = [
            row
            for row in controls
            if row["control"]["type"] == "connection_attempt"
        ]
        disconnects = [
            row
            for row in controls
            if row["control"]["type"] == "disconnect"
        ]
        self.assertEqual(len(attempts), 1)
        self.assertEqual(len(disconnects), 1)
        self.assertEqual(
            {
                events[0]["conn"],
                attempts[0]["conn"],
                disconnects[0]["conn"],
            },
            {events[0]["conn"]},
        )


    def test_non_exact_subscription_ack_never_reaches_ready(self) -> None:
        hybrid = {
            "id": 1,
            "result": None,
            "stream": "!forceOrder@arr",
            "data": self.force_order_payload(),
        }
        invalid_acks = {
            "boolean-id": {"id": True, "result": None},
            "float-id": {"id": 1.0, "result": None},
            "missing-result": {"id": 1},
            "extra-key": {"id": 1, "result": None, "extra": "ignored"},
            "hybrid-market-frame": hybrid,
        }

        for label, ack in invalid_acks.items():
            with self.subTest(label=label):
                collector, out = self.make_collector()
                collector.max_seconds = 1.0

                class InvalidAckWebSocket:
                    async def send(self, _message: str) -> None:
                        return None

                    async def recv(self) -> str:
                        collector.request_stop()
                        return json.dumps(ack)

                class InvalidAckConnection:
                    async def __aenter__(self):
                        return InvalidAckWebSocket()

                    async def __aexit__(self, _exc_type, _exc, _tb):
                        return False

                stdout = io.StringIO()
                loop = asyncio.new_event_loop()
                try:
                    with (
                        mock.patch.object(
                            liquidation.websockets,
                            "connect",
                            return_value=InvalidAckConnection(),
                        ),
                        mock.patch.object(sys, "stdout", stdout),
                    ):
                        rc = loop.run_until_complete(collector.run())
                finally:
                    loop.close()

                rows = [
                    liquidation.strict_json_loads(line)
                    for path in out.glob(
                        f"{liquidation.FILE_PREFIX}*.jsonl"
                    )
                    for line in path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                ]
                lifecycle = [
                    row
                    for row in rows
                    if row.get("kind") == "control"
                    and row.get("control", {}).get("type")
                    in {"connection_attempt", "connect", "disconnect"}
                ]
                self.assertEqual(rc, 0)
                self.assertEqual(collector.stats["subscription_acks"], 0)
                self.assertNotIn(liquidation.READY_BANNER, stdout.getvalue())
                self.assertEqual(
                    [
                        row["control"]["type"]
                        for row in lifecycle
                    ],
                    ["connection_attempt", "disconnect"],
                )
                self.assertEqual(
                    len({row["conn"] for row in lifecycle}),
                    1,
                )
                self.assertEqual(
                    lifecycle[-1]["control"]["phase"],
                    "pre-ack",
                )
                self.assertFalse(
                    any(row.get("kind") == "event" for row in rows)
                )
                self.assertFalse(
                    any(row.get("kind") == "oi_snapshot" for row in rows)
                )
    def test_sync_connect_factory_failure_has_attempt_provenance(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0

        def fail_connect(*_args, **_kwargs):
            collector.request_stop()
            raise OSError("injected synchronous connect failure")

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    side_effect=fail_connect,
                ),
                mock.patch.object(sys, "stdout", stdout),
            ):
                result = loop.run_until_complete(collector.run())
        finally:
            loop.close()

        self.assertEqual(result, 0)
        self.assertNotIn(liquidation.READY_BANNER, stdout.getvalue())
        controls = [
            row
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                liquidation.strict_json_loads(line)
                for line in path.read_text(
                    encoding="utf-8"
                ).splitlines()
            )
            if row.get("kind") == "control"
            and row.get("control", {}).get("type")
            in {"connection_attempt", "connect", "disconnect"}
        ]
        self.assertEqual(
            [row["control"]["type"] for row in controls],
            ["connection_attempt", "disconnect"],
        )
        self.assertEqual(len({row["conn"] for row in controls}), 1)
        self.assertTrue(controls[0]["conn"].startswith(collector.run_id))
        self.assertEqual(
            controls[-1]["control"]["phase"],
            "pre-ack",
        )

    def test_async_connect_enter_failure_has_attempt_provenance(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        entered = False

        class EnterFailure:
            async def __aenter__(self):
                nonlocal entered
                entered = True
                collector.request_stop()
                raise OSError("injected async enter failure")

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=EnterFailure(),
                ),
                mock.patch.object(sys, "stdout", stdout),
            ):
                result = loop.run_until_complete(collector.run())
        finally:
            loop.close()

        self.assertEqual(result, 0)
        self.assertTrue(entered)
        self.assertNotIn(liquidation.READY_BANNER, stdout.getvalue())
        controls = [
            row
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                liquidation.strict_json_loads(line)
                for line in path.read_text(
                    encoding="utf-8"
                ).splitlines()
            )
            if row.get("kind") == "control"
            and row.get("control", {}).get("type")
            in {"connection_attempt", "connect", "disconnect"}
        ]
        self.assertEqual(
            [row["control"]["type"] for row in controls],
            ["connection_attempt", "disconnect"],
        )
        self.assertEqual(len({row["conn"] for row in controls}), 1)
        self.assertEqual(
            controls[-1]["control"]["phase"],
            "pre-ack",
        )

    def test_successful_connection_control_order_is_attempt_connect_disconnect(
        self,
    ) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0

        class SuccessfulWebSocket:
            def __init__(self):
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                collector.request_stop()
                raise OSError("injected post-ready stop")

        class SuccessfulConnection:
            async def __aenter__(self):
                return SuccessfulWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=SuccessfulConnection(),
                ),
                mock.patch.object(sys, "stdout", stdout),
            ):
                result = loop.run_until_complete(collector.run())
        finally:
            loop.close()

        self.assertEqual(result, 0)
        self.assertIn(liquidation.READY_BANNER, stdout.getvalue())
        controls = [
            row
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                liquidation.strict_json_loads(line)
                for line in path.read_text(
                    encoding="utf-8"
                ).splitlines()
            )
            if row.get("kind") == "control"
            and row.get("control", {}).get("type")
            in {"connection_attempt", "connect", "disconnect"}
        ]
        self.assertEqual(
            [row["control"]["type"] for row in controls],
            ["connection_attempt", "connect", "disconnect"],
        )
        self.assertEqual(len({row["conn"] for row in controls}), 1)
        self.assertEqual(
            controls[-1]["control"]["phase"],
            "post-ready",
        )

    def test_attempt_cap_failure_prevents_connector_and_ready(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        connect = mock.Mock()
        emit_control = collector.emit_control

        def reject_attempt(control_type: str, **fields):
            if control_type == "connection_attempt":
                raise liquidation.CapExceeded(
                    "injected attempt control cap"
                )
            return emit_control(control_type, **fields)

        stdout = io.StringIO()
        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    connect,
                ),
                mock.patch.object(
                    collector,
                    "emit_control",
                    side_effect=reject_attempt,
                ),
                mock.patch.object(sys, "stdout", stdout),
            ):
                result = loop.run_until_complete(collector.run())
        finally:
            loop.close()

        self.assertEqual(result, 3)
        connect.assert_not_called()
        self.assertNotIn(liquidation.READY_BANNER, stdout.getvalue())
        controls = [
            row["control"]["type"]
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                liquidation.strict_json_loads(line)
                for line in path.read_text(
                    encoding="utf-8"
                ).splitlines()
            )
            if row.get("kind") == "control"
        ]
        self.assertNotIn("connection_attempt", controls)
        self.assertNotIn("connect", controls)
        self.assertNotIn("disconnect", controls)

    def test_pre_ack_cap_failure_is_terminal(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        collector.max_bytes_per_day = 4_000
        oversized_event = {
                    **self.force_order_payload(
                        int(liquidation.time.time() * 1000)
                    ),
            "padding": "x" * 10_000,
        }

        class CapWebSocket:
            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                return json.dumps(oversized_event)

        class CapConnection:
            async def __aenter__(self):
                return CapWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=CapConnection(),
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
                mock.patch.object(sys, "stderr", io.StringIO()),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 3)
        self.assertEqual(collector.stats["reconnects"], 0)
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(
            any(row.get("kind") == "event" for row in rows)
        )

    def test_pre_ack_sink_failure_is_terminal(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        event = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )
        connection_entries = 0

        class SinkWebSocket:
            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                return json.dumps(event)

        class SinkConnection:
            async def __aenter__(self):
                nonlocal connection_entries
                connection_entries += 1
                return SinkWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        original_write_line = collector.write_line

        def fail_event_write(record: dict, **kwargs) -> None:
            if record.get("kind") == "event":
                raise collector._latch_sink_io(
                    "append",
                    OSError("injected pre-ack failure"),
                )
            original_write_line(record, **kwargs)

        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    return_value=SinkConnection(),
                ),
                mock.patch.object(
                    collector,
                    "write_line",
                    side_effect=fail_event_write,
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
                mock.patch.object(sys, "stderr", io.StringIO()),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 5)
        self.assertEqual(connection_entries, 1)
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(
            any(row.get("kind") == "event" for row in rows)
        )

    def test_invalid_force_order_causes_one_controlled_reconnect(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        valid_frame = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )
        invalid_frame = {
            "e": "forceOrder",
            "E": valid_frame["E"],
            "o": "bad",
        }

        class SchemaWebSocket:
            def __init__(self, attempt: int):
                self.attempt = attempt
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                if self.attempt == 1:
                    return json.dumps(invalid_frame)
                collector.request_stop()
                return json.dumps(valid_frame)

        class SchemaConnection:
            def __init__(self, attempt: int):
                self.attempt = attempt

            async def __aenter__(self):
                return SchemaWebSocket(self.attempt)

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        attempts = iter((1, 2))

        def connect(_url: str, **_kwargs):
            return SchemaConnection(next(attempts))

        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    connect,
                ),
                mock.patch.object(
                    liquidation.random,
                    "uniform",
                    return_value=0,
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        self.assertEqual(collector.stats["reconnects"], 1)
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertEqual(
            sum(
                row.get("control", {}).get("type")
                == "force_order_schema_error"
                for row in rows
            ),
            1,
        )
        self.assertEqual(
            sum(row.get("kind") == "event" for row in rows),
            1,
        )

    def test_post_ready_error_frame_forces_audited_reconnect(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0
        valid_frame = self.force_order_payload(
            int(liquidation.time.time() * 1000)
        )
        error_frame = {
            "error": {
                "code": -1121,
                "msg": "Invalid symbol." + "x" * 600,
                "extra": {"nested": [1, 2, 3]},
            }
        }

        class ErrorWebSocket:
            def __init__(self, attempt: int):
                self.attempt = attempt
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                if self.attempt == 1:
                    return json.dumps(error_frame)
                collector.request_stop()
                return json.dumps(valid_frame)

        class ErrorConnection:
            def __init__(self, attempt: int):
                self.attempt = attempt

            async def __aenter__(self):
                return ErrorWebSocket(self.attempt)

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        attempts = iter((1, 2))
        stdout = io.StringIO()

        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(
                    liquidation.websockets,
                    "connect",
                    lambda _url, **_kwargs: ErrorConnection(next(attempts)),
                ),
                mock.patch.object(
                    liquidation.random,
                    "uniform",
                    return_value=0,
                ),
                mock.patch.object(collector, "schedule_oi_snapshot"),
                mock.patch.object(liquidation.sys, "stdout", stdout),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        controls = [
            row["control"]
            for row in rows
            if row.get("kind") == "control"
        ]
        errors = [
            control
            for control in controls
            if control["type"] == "error_frame"
        ]
        disconnects = [
            control
            for control in controls
            if control["type"] == "disconnect"
        ]
        self.assertEqual(rc, 0)
        self.assertEqual(collector.stats["reconnects"], 1)
        self.assertEqual(
            stdout.getvalue().count(liquidation.READY_BANNER),
            1,
        )
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0]["code"], -1121)
        self.assertEqual(len(errors[0]["msg"]), 256)
        self.assertEqual(errors[0]["phase"], "confirmed")
        self.assertEqual(
            [control.get("phase") for control in disconnects],
            ["post-ready", None],
        )
        self.assertEqual(disconnects[1]["reason"], "planned")
        self.assertIn("server error frame", disconnects[0]["reason"])
        self.assertLess(len(disconnects[0]["reason"]), 512)
        self.assertEqual(
            sum(row.get("kind") == "event" for row in rows),
            1,
        )

    def test_connection_identity_is_unique_across_process_restarts(self) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-connection-provenance-"))
        for run_index in range(2):
            collector = liquidation.Collector(
                out_dir=out,
                url="wss://example.invalid/ws",
                streams=["!forceOrder@arr"],
                max_seconds=None,
                max_bytes_per_day=10_000_000,
                max_total_bytes=100_000_000,
            )
            collector.load_state()
            collector.conn_seq = 1
            connection_id = collector._connection_id()
            collector.emit_control(
                "connect",
                conn=connection_id,
            )
            event_ms = 1_800_000_000_000 + run_index
            collector.dispatch_market_message(
                {
                    "e": "forceOrder",
                    "E": event_ms,
                    "o": {
                        "s": "BTCUSDT",
                        "S": "SELL",
                        "o": "LIMIT",
                        "f": "IOC",
                        "q": "1",
                        "p": "90000",
                        "ap": "89999",
                        "X": "FILLED",
                        "l": "1",
                        "z": "1",
                        "T": event_ms,
                    },
                },
                connection_id,
                event_ms,
            )
            collector.emit_control(
                "disconnect",
                conn=connection_id,
                reason="planned",
            )
            collector.fsync_close()
            collector.release_lease()

        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        connection_rows = [
            row
            for row in rows
            if row.get("kind") == "event"
            or row.get("control", {}).get("type")
            in {"connect", "disconnect"}
        ]
        run_ids = {row.get("run_id") for row in connection_rows}
        connection_ids = {row.get("conn") for row in connection_rows}
        self.assertEqual(len(run_ids), 2)
        self.assertNotIn(None, run_ids)
        self.assertEqual(len(connection_ids), 2)
        for run_id in run_ids:
            run_rows = [
                row for row in connection_rows
                if row.get("run_id") == run_id
            ]
            self.assertEqual(len(run_rows), 3)
            self.assertEqual(
                len({row["conn"] for row in run_rows}),
                1,
            )
            for row in run_rows:
                if row.get("kind") == "control":
                    self.assertEqual(
                        row["control"]["conn"],
                        row["conn"],
                    )

    def test_auxiliary_stream_flood_writes_one_control_per_connection(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.streams = [
            "btcusdt@aggTrade",
            "btcusdt@forceOrder",
        ]
        for trade_id in range(100):
            collector.dispatch_market_message(
                {
                    "stream": "btcusdt@aggTrade",
                    "data": {
                        "e": "aggTrade",
                        "E": 1_800_000_000_000 + trade_id,
                        "s": "BTCUSDT",
                        "a": trade_id,
                    },
                },
                "c1",
                1_800_000_000_000 + trade_id,
            )
        collector.dispatch_market_message(
            self.force_order_payload(1_800_000_001_000),
            "c1",
            1_800_000_001_001,
        )
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        ignored_controls = [
            row
            for row in rows
            if (
                row.get("control", {}).get("type")
                == "ignored_non_force_order_stream"
            )
        ]
        self.assertEqual(len(ignored_controls), 1)
        self.assertEqual(collector.stats["ignored_market_frames"], 100)
        self.assertEqual(
            len([row for row in rows if row.get("kind") == "event"]),
            1,
        )

    def test_disconnect_phase_is_per_connection_attempt(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 1.0

        class AttemptWebSocket:
            def __init__(self, attempt: int):
                self.attempt = attempt
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.attempt == 1 and self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                if self.attempt == 2:
                    collector.request_stop()
                raise OSError(f"attempt-{self.attempt}-drop")

        class AttemptConnection:
            def __init__(self, attempt: int):
                self.attempt = attempt

            async def __aenter__(self):
                return AttemptWebSocket(self.attempt)

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        attempts = iter((1, 2))

        def connect(_url: str, **_kwargs):
            return AttemptConnection(next(attempts))

        loop = asyncio.new_event_loop()
        try:
            with (
                mock.patch.object(liquidation.websockets, "connect", connect),
                mock.patch.object(liquidation.random, "uniform", return_value=0),
            ):
                rc = loop.run_until_complete(collector.run())
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        phases = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            phases.extend(
                row["control"]["phase"]
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "control"
                and row.get("control", {}).get("type") == "disconnect"
                and "phase" in row["control"]
            )
        self.assertEqual(phases, ["post-ready", "pre-ack"])

    def test_signal_shutdown_persists_disconnect_before_session_stop(self) -> None:
        collector, out = self.make_collector()
        collector.max_seconds = 0.3
        receive_started = asyncio.Event()

        class SignalWebSocket:
            def __init__(self):
                self.calls = 0

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                self.calls += 1
                if self.calls == 1:
                    return json.dumps({"id": 1, "result": None})
                receive_started.set()
                await asyncio.Future()
                raise AssertionError("unreachable")

        class SignalConnection:
            async def __aenter__(self):
                return SignalWebSocket()

            async def __aexit__(self, _exc_type, _exc, _tb):
                return False

        async def scenario() -> int:
            async def signal_after_receive() -> None:
                await receive_started.wait()
                collector.request_stop(signal.SIGTERM)

            signal_task = asyncio.create_task(signal_after_receive())
            rc = await collector.run()
            await signal_task
            return rc

        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(
                liquidation.websockets,
                "connect",
                return_value=SignalConnection(),
            ):
                rc = loop.run_until_complete(
                    asyncio.wait_for(scenario(), timeout=1.0)
                )
        finally:
            loop.close()
        self.assertEqual(rc, 0)
        controls = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            controls.extend(
                row["control"]
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "control"
            )
        types = [control["type"] for control in controls]
        self.assertLess(types.index("connect"), types.index("disconnect"))
        self.assertLess(types.index("disconnect"), types.index("session_stop"))
        self.assertEqual(controls[-1]["reason"], "sigterm")

    def test_pre_ack_liquidation_is_persisted_before_ack(self) -> None:
        collector, out = self.make_collector()
        received_ms = 1_800_000_000_123
        liquidation_frame = self.force_order_payload(received_ms - 10)

        class PreAckWebSocket:
            def __init__(self):
                self.messages = [
                    json.dumps(liquidation_frame),
                    json.dumps({"id": 1, "result": None}),
                ]

            async def send(self, _message: str) -> None:
                return None

            async def recv(self) -> str:
                return self.messages.pop(0)

        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(
                liquidation.time,
                "time",
                return_value=received_ms / 1000,
            ):
                result = loop.run_until_complete(
                    collector.subscribe_and(PreAckWebSocket(), "c1")
                )
                self.assertIsNone(result)
                collector.emit_control("connect", conn="c1")
                collector.dispatch_market_message(
                    liquidation_frame,
                    "c1",
                    received_ms + 1,
                )
        finally:
            loop.close()
            collector.fsync_close()

        events = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            events.extend(
                row
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "event"
            )
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["recv_epoch_ms"], received_ms)
        self.assertEqual(events[0]["subscription_phase"], "pre_ack")
        self.assertEqual(collector._oi_queue.qsize(), 1)


    def test_midnight_write_uses_one_clock_sample_for_day_and_timestamp(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector._state_save_due_monotonic = float("inf")
        timestamps = iter(
            [
                ("2026-09-03T23:59:59.999Z", 1_788_479_999.999),
                ("2026-09-04T00:00:00.001Z", 1_788_480_000.001),
            ]
        )
        with (
            mock.patch.object(
                liquidation,
                "utc_date",
                return_value="20260903",
            ),
            mock.patch.object(
                liquidation,
                "utc_now",
                side_effect=lambda: next(timestamps),
            ),
        ):
            collector.emit_event(
                "!forceOrder@arr",
                self.force_order_payload(2),
            )
        path = out / f"{liquidation.FILE_PREFIX}20260903.jsonl"
        rows = [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        self.assertTrue(rows)
        self.assertTrue(
            all(
                row["ts_utc"].startswith("2026-09-03")
                for row in rows
            )
        )
        self.assertEqual(rows[-1]["kind"], "event")
        self.assertFalse(
            (out / f"{liquidation.FILE_PREFIX}20260904.jsonl").exists()
        )

    def test_liquidation_rejects_invalid_utf8_without_replacement(
        self,
    ) -> None:
        class OneFrameWebSocket:
            async def send(self, _payload: str) -> None:
                return None

            async def recv(self):
                return frame

        payload = self.force_order_payload()
        payload["extra"] = "INVALID_UTF8"
        frame = json.dumps(payload).encode("utf-8").replace(
            b"INVALID_UTF8",
            b"\xff",
        )
        collector, out = self.make_collector()
        loop = asyncio.new_event_loop()
        try:
            with self.assertRaises(liquidation.ForceOrderSchemaError):
                loop.run_until_complete(
                    collector.subscribe_and(
                        OneFrameWebSocket(),
                        "c1",
                    )
                )
        finally:
            loop.close()
            collector.fsync_close()
        rows = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            rows.extend(
                json.loads(line)
                for line in path.read_text(encoding="utf-8").splitlines()
            )
        self.assertFalse(any(row.get("kind") == "event" for row in rows))
        audits = [
            row
            for row in rows
            if row.get("kind") == "control"
            and row.get("control", {}).get("type")
            == "invalid_utf8_frame"
        ]
        self.assertEqual(len(audits), 1)

    def test_oi_fetch_rejects_invalid_utf8_without_replacement(
        self,
    ) -> None:
        row = self.oi_payload()
        row["extra"] = "INVALID_UTF8"
        response = FakeHttpResponse([
            json.dumps([row]).encode("utf-8").replace(
                b"INVALID_UTF8",
                b"\xff",
            )
        ])
        with (
            mock.patch(
                "urllib.request.urlopen",
                return_value=response,
            ),
            mock.patch.object(
                liquidation,
                "utc_now",
                return_value=(
                    "2027-01-15T08:00:00.000Z",
                    1_800_000_000.0,
                ),
            ),
        ):
            snapshot, failure = (
                liquidation.fetch_open_interest_snapshot("BTCUSDT")
            )
        self.assertIsNone(snapshot)
        self.assertEqual(
            failure,
            {"reason": "invalid-json", "scope": "endpoint"},
        )

    def test_liquidation_rejects_non_rfc_force_order_constants(
        self,
    ) -> None:
        class OneFrameWebSocket:
            def __init__(self, frame: str):
                self.frame = frame
                self.received = False

            async def send(self, _payload: str) -> None:
                return None

            async def recv(self):
                if self.received:
                    raise AssertionError(
                        "non-RFC frame did not trigger reconnect"
                    )
                self.received = True
                return self.frame

        for constant in (
            float("nan"),
            float("inf"),
            float("-inf"),
        ):
            with self.subTest(constant=constant):
                collector, out = self.make_collector()
                payload = self.force_order_payload()
                payload["extra"] = constant
                loop = asyncio.new_event_loop()
                try:
                    with self.assertRaises(
                        liquidation.ForceOrderSchemaError
                    ):
                        loop.run_until_complete(
                            collector.subscribe_and(
                                OneFrameWebSocket(json.dumps(payload)),
                                "c1",
                            )
                        )
                finally:
                    loop.close()
                    collector.fsync_close()
                rows = []
                for path in out.glob(
                    f"{liquidation.FILE_PREFIX}*.jsonl"
                ):
                    rows.extend(
                        json.loads(line)
                        for line in path.read_text(
                            encoding="utf-8"
                        ).splitlines()
                    )
                self.assertFalse(
                    any(row.get("kind") == "event" for row in rows)
                )
                audits = [
                    row
                    for row in rows
                    if row.get("kind") == "control"
                    and row.get("control", {}).get("type")
                    == "non_finite_json_frame"
                ]
                self.assertEqual(len(audits), 1)
                self.assertEqual(
                    audits[0]["control"].get("reason"),
                    "non_finite_json_constant",
                )

    def test_oi_fetch_rejects_non_rfc_extra_constants(self) -> None:
        for constant in (
            float("nan"),
            float("inf"),
            float("-inf"),
        ):
            with self.subTest(constant=constant):
                row = self.oi_payload()
                row["extra"] = constant
                with mock.patch(
                    "urllib.request.urlopen",
                    return_value=self.oi_http_response(row),
                ):
                    snapshot, failure = (
                        liquidation.fetch_open_interest_snapshot("BTCUSDT")
                    )
                self.assertIsNone(snapshot)
                self.assertEqual(
                    failure,
                    {
                        "reason": "invalid-json",
                        "scope": "endpoint",
                        "error_type": "non-finite-constant",
                    },
                )

    def test_liquidation_sink_rejects_non_finite_before_append(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        with self.assertRaises(
            liquidation.JsonSerializationIntegrityError
        ):
            collector.write_line({
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "payload": {"extra": float("nan")},
            })
        self.assertFalse(path.exists())
        self.assertEqual(collector.total_bytes, 0)
        self.assertEqual(collector.seen_ids, set())

        collector.emit_control("after_non_finite_rejection")
        for line in path.read_text(encoding="utf-8").splitlines():
            liquidation.strict_json_loads(line)

    def test_liquidation_non_rfc_history_never_seeds_dedup(
        self,
    ) -> None:
        out = Path(tempfile.mkdtemp(prefix="collector-integrity-non-rfc-"))
        day = liquidation.utc_date()
        path = out / f"{liquidation.FILE_PREFIX}{day}.jsonl"
        path.write_text(
            json.dumps({
                "v": liquidation.SCHEMA_VERSION,
                "kind": "event",
                "event_id": "poisoned-event-id",
                "payload": {"extra": float("nan")},
            })
            + "\n",
            encoding="utf-8",
        )
        collector = liquidation.Collector(
            out_dir=out,
            url="wss://example.invalid/ws",
            streams=["!forceOrder@arr"],
            max_seconds=None,
            max_bytes_per_day=10_000_000,
            max_total_bytes=100_000_000,
        )
        self.addCleanup(collector.release_lease)
        collector.load_state()
        collector.seed_dedup_from_file(path)
        self.assertNotIn("poisoned-event-id", collector.seen_ids)
        self.assertEqual(
            collector.stats["dedup_rehydrate_malformed"],
            1,
        )

    def test_oi_fetch_rejects_schema_and_numeric_poison(self) -> None:
        missing = object()
        cases = (
            ("symbol-missing", "symbol", missing, "symbol"),
            ("symbol-wrong", "symbol", "ETHUSDT", "symbol"),
            ("timestamp-missing", "timestamp", missing, "timestamp"),
            ("timestamp-bool", "timestamp", True, "timestamp"),
            ("timestamp-zero", "timestamp", 0, "timestamp"),
            ("oi-missing", "sumOpenInterest", missing, "sumOpenInterest"),
            ("oi-number", "sumOpenInterest", 1.0, "sumOpenInterest"),
            ("oi-empty", "sumOpenInterest", "", "sumOpenInterest"),
            ("oi-nan", "sumOpenInterest", "NaN", "sumOpenInterest"),
            ("oi-infinite", "sumOpenInterest", "Infinity", "sumOpenInterest"),
            ("oi-negative", "sumOpenInterest", "-1", "sumOpenInterest"),
            (
                "notional-missing",
                "sumOpenInterestValue",
                missing,
                "sumOpenInterestValue",
            ),
            (
                "notional-number",
                "sumOpenInterestValue",
                1.0,
                "sumOpenInterestValue",
            ),
            (
                "notional-nan",
                "sumOpenInterestValue",
                "NaN",
                "sumOpenInterestValue",
            ),
            (
                "notional-infinite",
                "sumOpenInterestValue",
                "Infinity",
                "sumOpenInterestValue",
            ),
            (
                "notional-negative",
                "sumOpenInterestValue",
                "-1",
                "sumOpenInterestValue",
            ),
        )
        for name, field, value, expected_field in cases:
            with self.subTest(name=name):
                row = self.oi_payload()
                if value is missing:
                    del row[field]
                else:
                    row[field] = value
                with (
                    mock.patch(
                        "urllib.request.urlopen",
                        return_value=self.oi_http_response(row),
                    ),
                    mock.patch.object(
                        liquidation,
                        "utc_now",
                        return_value=(
                            "2027-01-15T08:00:00.000Z",
                            1_800_000_000.0,
                        ),
                    ),
                ):
                    snapshot, failure = (
                        liquidation.fetch_open_interest_snapshot("BTCUSDT")
                    )
                self.assertIsNone(snapshot)
                self.assertEqual(failure["reason"], "schema-mismatch")
                self.assertEqual(failure["scope"], "endpoint")
                self.assertEqual(failure["field"], expected_field)

    def test_oi_fetch_rejects_stale_and_future_source_time(self) -> None:
        now_s = 1_800_000_000.0
        cases = (
            (
                "stale",
                int(now_s * 1000)
                - liquidation.OI_SOURCE_MAX_AGE_MS
                - 1,
            ),
            (
                "future",
                int(now_s * 1000)
                + liquidation.OI_SOURCE_MAX_FUTURE_SKEW_MS
                + 1,
            ),
        )
        for violation, timestamp in cases:
            with self.subTest(violation=violation):
                with (
                    mock.patch(
                        "urllib.request.urlopen",
                        return_value=self.oi_http_response(
                            self.oi_payload(timestamp=timestamp)
                        ),
                    ),
                    mock.patch.object(
                        liquidation,
                        "utc_now",
                        return_value=(
                            "2027-01-15T08:00:00.000Z",
                            now_s,
                        ),
                    ),
                ):
                    snapshot, failure = (
                        liquidation.fetch_open_interest_snapshot("BTCUSDT")
                    )
                self.assertIsNone(snapshot)
                self.assertEqual(
                    failure["reason"],
                    "source-timestamp-integrity",
                )
                self.assertEqual(failure["scope"], "endpoint")
                self.assertEqual(failure["violation"], violation)

    def test_oi_fetch_preserves_valid_row_with_one_clock_sample(self) -> None:
        row = self.oi_payload()
        response = self.oi_http_response(row)
        with (
            mock.patch(
                "urllib.request.urlopen",
                return_value=response,
            ),
            mock.patch.object(
                liquidation,
                "utc_now",
                return_value=(
                    "2027-01-15T08:00:00.123Z",
                    1_800_000_000.123,
                ),
            ) as clock,
        ):
            snapshot, failure = liquidation.fetch_open_interest_snapshot(
                "BTCUSDT"
            )
        self.assertIsNone(failure)
        self.assertEqual(snapshot["symbol"], "BTCUSDT")
        self.assertEqual(snapshot["row"], row)
        self.assertEqual(
            snapshot["retrieved_utc"],
            "2027-01-15T08:00:00.123Z",
        )
        self.assertEqual(snapshot["received_epoch_ms"], 1_800_000_000_123)
        self.assertEqual(
            snapshot["source_timestamp_ms"],
            row["timestamp"],
        )
        clock.assert_called_once_with()

    def test_oi_endpoint_failure_does_not_block_healthy_sibling(
        self,
    ) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        calls = []

        def fetch(symbol: str, timeout: float = 5.0):
            del timeout
            calls.append(symbol)
            if symbol == "BADUSDT":
                return None, {
                    "reason": "schema-mismatch",
                    "scope": "endpoint",
                    "field": "symbol",
                }
            return {
                "symbol": symbol,
                "row": self.oi_payload(symbol=symbol),
                "retrieved_utc": "2027-01-15T08:00:00.000Z",
                "received_epoch_ms": 1_800_000_000_000,
                "source_timestamp_ms": 1_800_000_000_000,
                "url": "https://example.invalid",
            }, None

        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(
                liquidation,
                "fetch_open_interest_snapshot",
                side_effect=fetch,
            ):
                collector.schedule_oi_snapshot("BADUSDT")
                collector.schedule_oi_snapshot("BTCUSDT")
                self.assertEqual(collector._oi_queue.qsize(), 2)
                collector.stop.set()
                loop.run_until_complete(collector.oi_worker())
        finally:
            loop.close()

        self.assertEqual(calls, ["BADUSDT", "BTCUSDT"])
        self.assertEqual(collector._oi_health, "healthy")
        self.assertEqual(collector.stats["oi_failed"], 1)
        self.assertEqual(collector.stats["oi_ok"], 1)
        self.assertEqual(collector.stats["oi_circuit_opened"], 0)
        self.assertEqual(collector.stats["oi_circuit_recovered"], 0)
        self.assertEqual(collector.stats["oi_symbol_quarantined"], 1)
        self.assertIn(
            "BADUSDT",
            collector._oi_symbol_retry_after_monotonic,
        )
        self.assertNotIn(
            "BTCUSDT",
            collector._oi_symbol_retry_after_monotonic,
        )
        snapshots = [
            row
            for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl")
            for row in (
                liquidation.strict_json_loads(line)
                for line in path.read_text(
                    encoding="utf-8"
                ).splitlines()
            )
            if row.get("kind") == "oi_snapshot"
        ]
        self.assertEqual(
            [
                (row["symbol"], row["status"])
                for row in snapshots
            ],
            [
                ("BADUSDT", "failed-closed"),
                ("BTCUSDT", "ok"),
            ],
        )
        self.assertEqual(
            snapshots[0]["failure"],
            {
                "reason": "schema-mismatch",
                "scope": "endpoint",
                "field": "symbol",
            },
        )

        collector.stop.clear()
        btc_last_attempt = collector._oi_last_attempt["BTCUSDT"]
        with mock.patch.object(
            liquidation.time,
            "monotonic",
            return_value=(
                btc_last_attempt + collector.OI_MIN_INTERVAL_S + 1
            ),
        ):
            collector.schedule_oi_snapshot("BTCUSDT")
        self.assertEqual(collector._oi_queue.qsize(), 1)

    def test_oi_failure_scope_classification(self) -> None:
        global_failures = (
            {"reason": "timeout", "scope": "network"},
            {"reason": "http-error", "scope": "egress"},
            {"reason": "unexpected-error", "scope": "client"},
            {"reason": "collector-deadline", "scope": "collector"},
            {"reason": "service-error", "scope": "service"},
            {
                "reason": "http-error",
                "scope": "endpoint",
                "http_status": 429,
            },
            {
                "reason": "http-error",
                "scope": "endpoint",
                "http_status": 500,
            },
            {
                "reason": "http-error",
                "scope": "endpoint",
                "http_status": 503,
            },
        )
        symbol_failures = (
            {
                "reason": "schema-mismatch",
                "scope": "endpoint",
                "field": "symbol",
            },
            {
                "reason": "source-timestamp-integrity",
                "scope": "endpoint",
                "violation": "stale",
            },
            {"reason": "invalid-json", "scope": "endpoint"},
            {
                "reason": "http-error",
                "scope": "endpoint",
                "http_status": 400,
            },
            {
                "reason": "http-error",
                "scope": "endpoint",
                "http_status": 404,
            },
        )
        for failure in global_failures:
            with self.subTest(global_failure=failure):
                self.assertTrue(
                    liquidation.oi_failure_is_global(failure)
                )
        for failure in symbol_failures:
            with self.subTest(symbol_failure=failure):
                self.assertFalse(
                    liquidation.oi_failure_is_global(failure)
                )

    def test_oi_fetch_reports_http_451_without_response_body(self) -> None:
        error = urllib.error.HTTPError(
            "https://fapi.binance.com/example",
            451,
            "Unavailable For Legal Reasons",
            {},
            None,
        )
        with mock.patch("urllib.request.urlopen", side_effect=error):
            snapshot, failure = liquidation.fetch_open_interest_snapshot("BTCUSDT")
        self.assertIsNone(snapshot)
        self.assertEqual(
            failure,
            {"reason": "http-error", "http_status": 451, "scope": "egress"},
        )

    def test_oi_failure_circuit_suppresses_flood_and_recovers(self) -> None:
        collector, out = self.make_collector()
        self.addCleanup(collector.fsync_close)
        collector.open_for(liquidation.utc_date())
        calls: list[str] = []

        def blocked(symbol: str, timeout: float = 5.0):
            del timeout
            calls.append(symbol)
            return None, {
                "reason": "http-error",
                "http_status": 451,
                "scope": "egress",
            }

        loop = asyncio.new_event_loop()
        try:
            with mock.patch.object(liquidation, "fetch_open_interest_snapshot", blocked):
                collector.schedule_oi_snapshot("BTCUSDT")
                collector.schedule_oi_snapshot("ETHUSDT")
                self.assertEqual(collector._oi_queue.qsize(), 2)
                collector.stop.set()
                loop.run_until_complete(collector.oi_worker())

            self.assertEqual(calls, ["BTCUSDT"])
            self.assertEqual(collector._oi_health, "blocked")
            self.assertGreater(collector._oi_retry_after_monotonic, 0)
            self.assertGreaterEqual(collector.stats["oi_circuit_skipped"], 1)

            skipped_before = collector.stats["oi_circuit_skipped"]
            collector.stop.clear()
            collector.schedule_oi_snapshot("BNBUSDT")
            collector.schedule_oi_snapshot("SOLUSDT")
            self.assertTrue(collector._oi_queue.empty())
            self.assertEqual(calls, ["BTCUSDT"])
            self.assertEqual(
                collector.stats["oi_circuit_skipped"],
                skipped_before + 2,
            )

            collector.stop.clear()
            collector._oi_retry_after_monotonic = 0.0

            def recovered(symbol: str, timeout: float = 5.0):
                del timeout
                calls.append(symbol)
                return {
                    "symbol": symbol,
                    "row": {"sumOpenInterest": "1"},
                    "retrieved_utc": "2026-09-04T00:00:00.000Z",
                    "url": "https://example.invalid",
                }, None

            with mock.patch.object(liquidation, "fetch_open_interest_snapshot", recovered):
                collector.schedule_oi_snapshot("ETHUSDT")
                collector.stop.set()
                loop.run_until_complete(collector.oi_worker())
        finally:
            loop.close()

        self.assertEqual(collector._oi_health, "healthy")
        snapshots = []
        for path in out.glob(f"{liquidation.FILE_PREFIX}*.jsonl"):
            snapshots.extend(
                row
                for row in (
                    json.loads(line)
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                if row.get("kind") == "oi_snapshot"
            )
        self.assertEqual(snapshots[0]["retry_after_seconds"], 3600)
        self.assertEqual([row["status"] for row in snapshots], ["failed-closed", "ok"])
        self.assertEqual(snapshots[0]["failure"]["http_status"], 451)


if __name__ == "__main__":
    unittest.main()
