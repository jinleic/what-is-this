#!/usr/bin/env python3
"""Append-only, unauthenticated public liquidation capture collector.

Contract (focused-campaign rank-3 `focused-liquidation-overlay`):
  - Public, UNAUTHENTICATED websocket capture only. No orders, no credentials,
    no private endpoints, no messages, no infrastructure mutation.
  - Append-only JSON Lines, UTC everywhere, deterministic event ids,
    daily rotation by UTC date. Re-delivery across reconnects/restarts is
    possible: the raw corpus is AT-LEAST-ONCE and replay-side deduplication
    (run.py, deterministic event ids) is authoritative.
  - Reconnect with capped exponential backoff + jitter; the append sink stays
    open across reconnects and is reopened/rotated by write_line on demand.
    Readiness requires the EXACT Binance ack {"id":1,"result":null} (no error
    key) with a bounded, stop-aware wait.
  - Fail-closed byte caps: default 256 MiB per UTC day and 10 GiB cumulative
    (persisted across restarts, floored by real on-disk bytes). Once a cap
    trips, no further writes go through the sink; shutdown is clean. Never
    auto-delete prior capture.
  - Bounded smoke mode via --max-seconds. Stop-aware SIGTERM/SIGINT shutdown.
  - Readiness: prints and flushes the exact ASCII line `LIQUIDATION_CAPTURE_READY`
    on stdout only after connection + exact subscription ack.

Line schema (one JSON object per physical line; never modified after write):
  event:   {"v":1,"kind":"event","ts_utc":...,"recv_epoch_ms":...,"conn":...,
            "stream":...,"event_id":...,"source":{...},"payload":{...}}
  control: {"v":1,"kind":"control","ts_utc":...,"conn":...,
            "control":{"type":"...","reason"/fields...}}
"""

from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import os
import random
import signal
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    import websockets  # type: ignore
except ImportError:  # pragma: no cover
    print("capture.py requires the 'websockets' package", file=sys.stderr)
    raise SystemExit(2)

SCHEMA_VERSION = 1
MI = 1024 * 1024
GI = 1024 * MI
DEFAULT_MAX_BYTES_PER_DAY = 256 * MI
DEFAULT_MAX_TOTAL_BYTES = 10 * GI
DEFAULT_URL = "wss://fstream.binance.com/market/ws"
DEFAULT_STREAMS = ["!forceOrder@arr"]
READY_BANNER = "LIQUIDATION_CAPTURE_READY"

FILE_PREFIX = "liquidation-capture-"
STATE_FILE = "collector-state.json"


def utc_now() -> tuple[str, float]:
    now = datetime.now(timezone.utc)
    return now.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z", now.timestamp()


def utc_date() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d")


def canonical_force_order(stream: str, payload: dict) -> tuple[str, dict] | None:
    """Normalize a liquidation frame to (canonical_stream, body) regardless of
    transport shape. Returns None for frames that are not recognizable
    liquidation payloads (the caller decides how to fail closed).

    Shapes accepted (Binance emits BOTH depending on endpoint):
      wrapped: {"stream":"!forceOrder@arr","data":{"e":"forceOrder","E":..,"o":{..}}}
      raw:     {"e":"forceOrder","E":..,"o":{..}}   (dynamic SUBSCRIBE on /ws)
    """
    if not isinstance(payload, dict):
        return None
    body = payload.get("data") if isinstance(payload.get("data"), dict) else payload
    stream = str(payload.get("stream") or stream or "")
    if stream.endswith("@forceOrder"):
        return stream, body
    if body.get("e") == "forceOrder":
        # Raw (/ws) frames carry no stream name; map to the canonical
        # all-market stream for our single-stream default so wrapped and raw
        # forms produce the SAME semantic stream and event id.
        return "!forceOrder@arr", body
    return None


def event_id_for(stream: str, payload: dict) -> str:
    """Deterministic id from exchange-semantic identity, not from receipt
    time. Wrapped and raw frames of the SAME liquidation canonicalize to the
    identical id (canonical stream + E + o.* fields), so dedup holds across
    transport shapes."""
    canon = canonical_force_order(stream, payload)
    if canon is not None:
        canon_stream, body = canon
        o = body.get("o", body)
        parts = [canon_stream, str(body.get("E", ""))]
        for key in ("s", "S", "o", "f", "q", "p", "ap", "X", "l", "z", "T"):
            parts.append(str(o.get(key, "")))
        return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]
    body = payload.get("data", payload) if isinstance(payload, dict) else {}
    stream = str(payload.get("stream") or stream or "?")
    if isinstance(body, dict) and body.get("e") == "aggTrade":
        parts = [stream, str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("a", ""))]
        return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]
    if isinstance(body, dict) and body.get("e") == "markPriceUpdate":
        parts = [stream, str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("T", ""))]
        return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]
    parts = [stream, json.dumps(body, sort_keys=True, separators=(",", ":"))]
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


class CapExceeded(Exception):
    pass


OI_REST_BASE = "https://fapi.binance.com/futures/data/openInterestHist"
OI_SNAPSHOT_PREFIX = "oi-snapshot-"


def fetch_open_interest_snapshot(symbol: str, timeout: float = 5.0) -> dict | None:
    """Bounded, unauthenticated OI snapshot. FAILS CLOSED: returns None on any
    error/region-block/non-200; never raises, never blocks the event loop
    caller beyond the socket timeout. Persisted by run.py augmentation only
    when a real JSON payload with sumOpenInterest fields is present."""
    try:
        import urllib.parse
        import urllib.request
        url = (f"{OI_REST_BASE}?" + urllib.parse.urlencode(
            {"symbol": symbol, "period": "5m", "limit": 1}))
        req = urllib.request.Request(url, headers={"User-Agent": "research-capture/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if resp.status != 200:
                return None
            data = json.loads(resp.read().decode("utf-8", "replace"))
        if (isinstance(data, list) and data and isinstance(data[0], dict)
                and "sumOpenInterest" in data[0]):
            return {"symbol": symbol, "row": data[0], "retrieved_utc": utc_now()[0],
                    "url": url}
        return None
    except Exception:
        # Restricted region, DNS failure, timeout, schema drift: fail closed.
        return None


class Collector:
    def __init__(self, out_dir: Path, url: str, streams: list[str], max_seconds: float | None,
                 max_bytes_per_day: int, max_total_bytes: int):
        self.out_dir = out_dir
        self.url = url.rstrip("/")
        self.streams = list(streams)
        self.max_seconds = max_seconds
        self.max_bytes_per_day = max_bytes_per_day
        self.max_total_bytes = max_total_bytes

        self.stop = asyncio.Event()
        self.file = None
        self.file_date = ""
        self.total_bytes = 0      # bytes written by THIS process
        self.state_total_bytes = 0  # bytes persisted by all previous runs
        self.day_bytes: dict[str, int] = {}
        self.seen_ids: set[str] = set()
        self.conn_seq = 0
        self._oi_queue: asyncio.Queue = asyncio.Queue(maxsize=64)
        self._oi_last_attempt: dict[str, float] = {}
        self.stats = {
            "events_written": 0,
            "events_deduped": 0,
            "controls_written": 0,
            "reconnects": 0,
            "conns": 0,
            "bytes_written": 0,
            "subscription_acks": 0,
            "oi_ok": 0,
            "oi_failed": 0,
            "oi_rate_limited": 0,
            "oi_dropped_full_queue": 0,
        }

    # ---------- state (cumulative cap persistence) ----------
    def actual_capture_bytes(self) -> int:
        """Total bytes on disk across all capture files (fail-closed floor)."""
        total = 0
        if self.out_dir.exists():
            for path in self.out_dir.glob(f"{FILE_PREFIX}*.jsonl"):
                try:
                    total += path.stat().st_size
                except OSError:
                    continue
        return total

    def load_state(self) -> None:
        path = self.out_dir / STATE_FILE
        if path.exists():
            try:
                st = json.loads(path.read_text(encoding="utf-8"))
                self.state_total_bytes = int(st.get("total_bytes", 0))
                for day, nb in (st.get("day_bytes") or {}).items():
                    self.day_bytes[str(day)] = max(nb, int(self.day_bytes.get(day, 0)))
            except (json.JSONDecodeError, ValueError, OSError):
                # Fail closed: unknown prior usage must be treated as maxed out.
                print("collector-state.json unreadable; failing closed on caps",
                      file=sys.stderr)
                self.state_total_bytes = self.max_total_bytes + 1
        # Reconcile with real capture files in ALL paths (only-raise semantics):
        # 1) per-day accounting uses max(persisted, on-disk) so day caps cannot
        #    be dodged by a missing/rolled-back day_bytes;
        # 2) the cumulative prior total is floored by the sum of on-disk
        #    capture-file sizes, so a crash that loses collector-state.json
        #    cannot undercount prior bytes against the cumulative cap.
        self.reconcile_day_bytes()
        self.state_total_bytes = max(self.state_total_bytes, self.actual_capture_bytes())

    def save_state(self) -> None:
        path = self.out_dir / STATE_FILE
        tmp = path.with_suffix(".json.tmp")
        payload = {
            "schema_version": SCHEMA_VERSION,
            "total_bytes": self.state_total_bytes + self.total_bytes,
            "day_bytes": dict(sorted(self.day_bytes.items())),
            "updated_utc": utc_now()[0],
        }
        with tmp.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)

    # ---------- byte accounting / caps ----------
    # `day_bytes` in collector-state.json is the SINGLE owner of per-UTC-day
    # byte accounting. On startup it is reconciled against the actual on-disk
    # capture-file sizes with a per-day max() so either source can only raise,
    # never lower, the counted usage (fail-closed reconciliation). The current
    # process adds to it via commit(); total = persisted prior total + current
    # process total, counted exactly once.
    _file_day = ""

    def reconcile_day_bytes(self) -> None:
        """Merge real capture-file sizes into day_bytes (never decrease)."""
        if not self.out_dir.exists():
            return
        for path in self.out_dir.glob(f"{FILE_PREFIX}*.jsonl"):
            day = path.stem[len(FILE_PREFIX):]
            if not day.isdigit():
                continue
            try:
                size = path.stat().st_size
            except OSError:
                continue
            self.day_bytes[day] = max(int(self.day_bytes.get(day, 0)), size)

    def _day_used(self, day: str) -> int:
        return int(self.day_bytes.get(day, 0))

    def check_caps(self, incoming: int, day: str | None = None) -> None:
        """Caps are checked against the TARGET day (the open file's day after
        rotation), never the wall clock, so a midnight control write cannot
        exceed the old day's budget with new-day bytes or vice versa."""
        target_day = day or self._file_day or utc_date()
        proj_day = self._day_used(target_day) + incoming
        if proj_day > self.max_bytes_per_day:
            raise CapExceeded(
                f"per-day cap: {proj_day} > {self.max_bytes_per_day} bytes on {target_day}")
        proj_total = self.state_total_bytes + self.total_bytes + incoming
        if proj_total > self.max_total_bytes:
            raise CapExceeded(f"cumulative cap: {proj_total} > {self.max_total_bytes} bytes")

    def commit(self, raw_bytes: int, day: str) -> None:
        self.total_bytes += raw_bytes
        self.day_bytes[day] = int(self.day_bytes.get(day, 0)) + raw_bytes

    _in_open_for = False

    def open_for(self, day: str) -> None:
        if self._in_open_for:
            return  # guard: control writes from inside open_for re-enter here
        self._in_open_for = True
        try:
            prev_day = self._file_day if self.file is not None else None
            if self.file is not None:
                self.fsync_close()
            path = self.out_dir / f"{FILE_PREFIX}{day}.jsonl"
            existed = path.exists()
            # Open the NEW day file BEFORE emitting any control so the rotate
            # record itself is written into the new file, never a closed handle.
            self.file = path.open("a", encoding="utf-8", buffering=1)
            self.file_date = day
            self._file_day = day
            if not existed:
                self.emit_control("file_created", path=path.name)
            else:
                self.seed_dedup_from_file(path)
                # Reconcile this day's usage with the real size (never lower).
                try:
                    self.day_bytes[day] = max(int(self.day_bytes.get(day, 0)),
                                              path.stat().st_size)
                except OSError:
                    pass
                self.emit_control("file_reopened", path=path.name,
                                  byte_count=path.stat().st_size if path.exists() else 0)
            if prev_day is not None and prev_day != day:
                self.emit_control("rotate", reason="utc midnight day change",
                                  from_day=prev_day, to_day=day)
        finally:
            self._in_open_for = False

    def seed_dedup_from_file(self, path: Path) -> None:
        try:
            with path.open(encoding="utf-8") as fh:
                for line in fh:
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if rec.get("kind") == "event" and rec.get("event_id"):
                        self.seen_ids.add(str(rec["event_id"]))
        except OSError:
            pass

    def fsync_close(self) -> None:
        if self.file is not None:
            try:
                self.file.flush()
                os.fsync(self.file.fileno())
            finally:
                self.file.close()
                self.file = None

    def write_line(self, record: dict) -> None:
        # The sink day is authoritative after an explicit open_for(); it only
        # auto-advances at the UTC boundary (target == wall-clock day) or on
        # the very first write of the process (sink never opened). This
        # prevents a reopened-day file from being silently pulled back to the
        # wall-clock day, which would misattribute bytes and dodge that day's
        # cap.
        wall_day = utc_date()
        if self.file is None or (
                self._file_day != wall_day and wall_day > self._file_day):
            self.open_for(wall_day)
        line = json.dumps(record, separators=(",", ":"), ensure_ascii=False)
        nbytes = len(line.encode("utf-8")) + 1
        # Cap check against the file day the bytes will be attributed to.
        self.check_caps(nbytes, self._file_day)
        assert self.file is not None
        self.file.write(line + "\n")
        self.file.flush()
        os.fsync(self.file.fileno())
        self.commit(nbytes, self._file_day)
        if record.get("kind") == "event":
            self.stats["events_written"] += 1

    def emit_event(self, stream: str, raw_payload: dict) -> None:
        ts, _ = utc_now()
        eid = event_id_for(stream, raw_payload)
        if eid in self.seen_ids:
            self.stats["events_deduped"] += 1
            return
        if len(self.seen_ids) > 2_000_000:
            # Bounded memory; clear only in a documented, replay-auditable way.
            self.seen_ids.clear()
            self.emit_control("dedup_cache_cleared", reason="2_000_000 id bound")
        self.seen_ids.add(eid)
        self.write_line({
            "v": SCHEMA_VERSION,
            "kind": "event",
            "ts_utc": ts,
            "recv_epoch_ms": int(time.time() * 1000),
            "conn": f"c{self.conn_seq}",
            "stream": stream,
            "event_id": eid,
            "source": {
                "venue": "binance-usds-m-futures",
                "endpoint_kind": "public-unauthenticated-websocket",
                "url": self.url,
            },
            "payload": raw_payload,
        })

    def emit_control(self, ctype: str, **fields) -> None:
        ts, _ = utc_now()
        self.write_line({
            "v": SCHEMA_VERSION,
            "kind": "control",
            "ts_utc": ts,
            "conn": f"c{self.conn_seq}",
            "control": {"type": ctype, **fields},
        })
        self.stats["controls_written"] += 1

    # ---------- websocket session ----------
    async def subscribe_and(self, ws, ack_timeout: float = 30.0) -> None:
        """Send SUBSCRIBE and await the EXACT authoritative ack:
        {"id":1,"result":null} with no error key. Anything else (non-null
        result, error key, different id, timeout, or stop) never satisfies
        readiness. Wait is bounded and stop-aware."""
        sub = {"method": "SUBSCRIBE", "id": 1, "params": self.streams}
        await ws.send(json.dumps(sub))
        deadline = time.monotonic() + ack_timeout
        while not self.stop.is_set():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"subscription ack not received within {ack_timeout}s")
            try:
                frame = await asyncio.wait_for(ws.recv(), timeout=min(remaining, 2.0))
            except asyncio.TimeoutError:
                continue
            text = frame if isinstance(frame, str) else frame.decode("utf-8", "replace")
            try:
                msg = json.loads(text)
            except json.JSONDecodeError:
                continue  # not the ack; keep waiting within the bound
            if not isinstance(msg, dict):
                continue
            if "error" in msg:
                raise RuntimeError(f"subscription error response: {msg}")
            if msg.get("id") != 1:
                if "result" in msg:
                    # A result-bearing response for a DIFFERENT id while we
                    # await OUR ack is a session anomaly: reject, never ready.
                    raise RuntimeError(
                        f"response for unexpected id while awaiting ack: {msg}")
                continue  # unrelated market frame before our ack
            # EXACT ack: id==1, no error key, result is precisely null.
            if "result" in msg and msg["result"] is None:
                self.stats["subscription_acks"] += 1
                return
            raise RuntimeError(f"non-authoritative subscription ack: {msg}")
        raise RuntimeError("stopped while awaiting subscription ack")

    def request_stop(self, *_args) -> None:
        self.stop.set()
    # ---------- bounded, non-blocking OI enrichment ----------
    OI_MIN_INTERVAL_S = 60.0      # rate limit per symbol (public REST courtesy)
    OI_QUEUE_MAX = 64             # bounded queue: stream receipt is never delayed

    def schedule_oi_snapshot(self, symbol: str) -> None:
        """Queue a bounded OI snapshot fetch for this symbol. Deduped by rate
        limit; counted-drop when the bounded queue is full so stream receipt
        is never delayed or blocked by enrichment. Idempotent if the queue
        task has not started yet: the worker simply awaits items."""
        if not symbol:
            return
        now = time.monotonic()
        last = self._oi_last_attempt.get(symbol)
        if last is not None and (now - last) < self.OI_MIN_INTERVAL_S:
            self.stats["oi_rate_limited"] = self.stats.get("oi_rate_limited", 0) + 1
            return
        if self._oi_queue.qsize() >= self.OI_QUEUE_MAX:
            self.stats["oi_dropped_full_queue"] = self.stats.get("oi_dropped_full_queue", 0) + 1
            return
        self._oi_last_attempt[symbol] = now
        # True waiting semantics: put_nowait on a bounded asyncio.Queue;
        # the worker consumes with a plain await (no deque.popleft on an
        # empty idle queue — the rereview's IndexError case).
        try:
            self._oi_queue.put_nowait(symbol)
        except asyncio.QueueFull:
            self.stats["oi_dropped_full_queue"] = self.stats.get("oi_dropped_full_queue", 0) + 1

    def _oi_fetch_sync(self, symbol: str) -> dict | None:
        """Executor target: runs the blocking fetch in a worker thread."""
        try:
            return fetch_open_interest_snapshot(symbol, 5.0)
        except TypeError:
            # Test doubles / legacy signatures may take only the symbol.
            return fetch_open_interest_snapshot(symbol)

    async def oi_worker(self) -> None:
        """Single background worker: awaits symbols from the bounded queue
        (properly idle-blocked when empty), fetches bounded snapshots, and
        appends success OR explicit fail-closed evidence through the same
        cap-checked sink. Never blocks the event receive loop; exits when
        stopped and drained. CapExceeded propagates to the run() envelope."""
        while True:
            if self.stop.is_set() and self._oi_queue.empty():
                return
            get_task = asyncio.ensure_future(self._oi_queue.get())
            stop_task = asyncio.ensure_future(self.stop.wait())
            done, pending = await asyncio.wait({get_task, stop_task},
                                               return_when=asyncio.FIRST_COMPLETED)
            for t in pending:
                t.cancel()
            if stop_task in done and not (get_task in done and not get_task.exception()):
                # Stop requested and no item became available: drained exit.
                get_task.cancel()
                return
            symbol = get_task.result()
            # Blocking fetch runs in the executor thread; the event receive
            # loop is never delayed beyond normal scheduling.
            snapshot = await asyncio.get_event_loop().run_in_executor(
                None, self._oi_fetch_sync, symbol)
            record = {
                "v": SCHEMA_VERSION,
                "kind": "oi_snapshot",
                "ts_utc": utc_now()[0],
                "symbol": symbol,
            }
            if snapshot is not None:
                record["status"] = "ok"
                record["data"] = snapshot
                self.stats["oi_ok"] = self.stats.get("oi_ok", 0) + 1
            else:
                record["status"] = "failed-closed"
                record["note"] = ("open interest snapshot unavailable from this "
                                  "egress (region block/timeout/schema drift)")
                self.stats["oi_failed"] = self.stats.get("oi_failed", 0) + 1
            self.write_line(record)  # cap trips propagate to the run envelope

    async def run(self) -> int:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.load_state()
        if self.state_total_bytes >= self.max_total_bytes:
            print("total cap already reached; failing closed", file=sys.stderr)
            return 4
        day = utc_date()
        if self._day_used(day) >= self.max_bytes_per_day:
            print("per-day cap already reached; failing closed", file=sys.stderr)
            return 4

        deadline = (time.monotonic() + self.max_seconds) if self.max_seconds else None

        def left() -> float:
            return (deadline - time.monotonic()) if deadline else float("inf")

        async def stop_aware_sleep(seconds: float) -> None:
            # Sleeps in small slices so SIGTERM interrupts promptly.
            end = time.monotonic() + seconds
            while not self.stop.is_set():
                slice_s = min(0.5, max(0.0, end - time.monotonic()))
                if slice_s <= 0:
                    return
                await asyncio.sleep(slice_s)

        attempt = 0
        ready_announced = False
        nonlocal_cap_reason: list[str | None] = [None]
        self.exit_code = 0

        async def regulated_connect_loop() -> None:
            """The whole session lifecycle lives inside ONE CapExceeded
            envelope, so every control write — startup file_created/
            file_reopened/session_start, each connect/disconnect, planned or
            error disconnects — passes through the same cap-checked sink path.
            Nothing escapes the boundary."""
            nonlocal attempt, ready_announced
            try:
                self.open_for(day)
                self.emit_control("session_start", url=self.url,
                                  streams=self.streams, pid=os.getpid(),
                                  max_bytes_per_day=self.max_bytes_per_day,
                                  max_total_bytes=self.max_total_bytes)
                self.save_state()
                while not self.stop.is_set() and left() > 0:
                    try:
                        self.conn_seq += 1
                        cid = f"c{self.conn_seq}"
                        async with websockets.connect(
                            self.url, ping_interval=180, ping_timeout=600,
                            open_timeout=20, close_timeout=10,
                        ) as ws:
                            await self.subscribe_and(ws)
                            attempt = 0
                            self.emit_control("connect", conn=cid, url=self.url,
                                              streams=self.streams)
                            self.save_state()
                            # Durable-gate the banner: readiness is claimed
                            # only after the connect control AND state save
                            # have hit the cap-checked sink successfully.
                            if not ready_announced:
                                sys.stdout.write(READY_BANNER + "\n")
                                sys.stdout.flush()
                                ready_announced = True
                            else:
                                print("resubscribed", file=sys.stderr)
                            # event loop for this connection
                            while not self.stop.is_set():
                                remaining = left()
                                if remaining <= 0:
                                    break
                                try:
                                    frame = await asyncio.wait_for(
                                        ws.recv(),
                                        timeout=min(2.0, max(0.05, remaining)))
                                except asyncio.TimeoutError:
                                    continue
                                text = frame if isinstance(frame, str) else frame.decode("utf-8", "replace")
                                try:
                                    msg = json.loads(text)
                                except json.JSONDecodeError:
                                    self.emit_control("bad_frame", conn=cid,
                                                      snippet=text[:200])
                                    continue
                                if isinstance(msg, dict) and "error" in msg:
                                    self.emit_control("error_frame", conn=cid,
                                                      error=msg["error"])
                                    continue
                                if isinstance(msg, dict) and msg.get("id") == 1:
                                    continue  # ack echoes are not events
                                # Dispatch BOTH transport shapes:
                                #   wrapped {"stream":..,"data":{..}}
                                #   raw {"e":"forceOrder",..} (dynamic
                                #   SUBSCRIBE on /market/ws emits this shape)
                                canon = canonical_force_order("", msg)
                                if canon is not None:
                                    canon_stream, canon_body = canon
                                    self.emit_event(canon_stream, msg)
                                    # OI enrichment: symbol from EITHER shape.
                                    o_obj = canon_body.get("o")
                                    ev_symbol = (o_obj or {}).get("s") \
                                        if isinstance(o_obj, dict) else None
                                    if ev_symbol:
                                        self.schedule_oi_snapshot(str(ev_symbol))
                                    continue
                                # Fail closed on unknown/ambiguous RAW event
                                # types: never guess a stream for an
                                # unrecognized payload.
                                body_probe = msg.get("data") \
                                    if isinstance(msg.get("data"), dict) else msg
                                if isinstance(body_probe, dict) and body_probe.get("e"):
                                    self.emit_control("unknown_raw_event_type",
                                                      conn=cid,
                                                      event_type=str(body_probe.get("e")))
                                    continue
                                if msg.get("stream"):
                                    # Wrapped non-forceOrder frame (multi-
                                    # stream subscription): record for audit,
                                    # never misclassify as a liquidation.
                                    self.emit_control("non_force_order_stream",
                                                      conn=cid,
                                                      stream=str(msg.get("stream")))
                                    continue
                                self.emit_control("unwrapped_frame", conn=cid)
                                continue
                            self.emit_control("disconnect", conn=cid,
                                              reason="planned")
                            self.save_state()
                            if self.stop.is_set():
                                break
                            continue
                    except (websockets.exceptions.ConnectionClosed, OSError,
                            socket.gaierror, asyncio.TimeoutError, TimeoutError,
                            RuntimeError) as exc:
                        # Failure-before-ack and failure-after-ready both land
                        # here; the sink remains open through the control
                        # write below (all inside the cap envelope).
                        self.stats["reconnects"] += 1
                        attempt += 1
                        backoff = min(60.0, (2 ** min(attempt, 6)) * random.uniform(0.8, 1.2))
                        self.emit_control(
                            "disconnect",
                            conn=f"c{self.conn_seq}",
                            reason=f"{type(exc).__name__}: {exc}",
                            backoff_s=round(backoff, 2),
                            phase=("pre-ack" if self.stats["subscription_acks"] == 0
                                   else "post-ready"))
                        self.save_state()
                        if self.stop.is_set() or left() <= 0:
                            break
                        await stop_aware_sleep(min(backoff, max(left(), 0.1)))
            except CapExceeded as exc:
                # Terminal from ANY control write inside the envelope: freeze
                # accounting, close, save, signal out-of-band via return code.
                nonlocal_cap_reason[0] = str(exc)

        oi_task: asyncio.Task | None = None
        try:
            oi_task = asyncio.ensure_future(self.oi_worker())
            connect_task = asyncio.ensure_future(regulated_connect_loop())
            # Coordinate both tasks: whichever finishes FIRST with an
            # exception (typically an OI-side CapExceeded) becomes the shared
            # terminal cap reason; we then stop the loop immediately rather
            # than swallowing the cap at an outer gather.
            # Coordinate both tasks: whichever finishes FIRST with an
            # exception (typically an OI-side CapExceeded) becomes the shared
            # terminal cap reason. A normal connect-loop exit (max-seconds,
            # SIGTERM) also terminates: the idle worker is then cancelled —
            # never awaited to completion, since a stopped-and-drained worker
            # returns but an in-flight fetch may hold the queue briefly.
            while True:
                done, pending = await asyncio.wait(
                    {connect_task, oi_task}, return_when=asyncio.FIRST_COMPLETED)
                terminal = None
                for t in done:
                    if t.cancelled():
                        continue
                    exc = t.exception()
                    if isinstance(exc, CapExceeded):
                        terminal = str(exc)
                if terminal:
                    nonlocal_cap_reason[0] = terminal
                    self.stop.set()
                    break
                if connect_task in done:
                    # Normal connect-loop end (stop/max-seconds). Grace-drain
                    # the OI worker briefly, then cancel it.
                    self.stop.set()
                    try:
                        await asyncio.wait_for(asyncio.shield(oi_task), timeout=5.0)
                    except (asyncio.TimeoutError, asyncio.CancelledError,
                            CapExceeded, Exception):
                        pass
                    break
                if oi_task in done:
                    # Worker exited on its own (queued drained post-stop or
                    # unexpected non-cap crash): keep running the loop.
                    if self.stop.is_set():
                        break
                    continue
            for t in (connect_task, oi_task):
                if not t.done():
                    t.cancel()
            for t in (connect_task, oi_task):
                try:
                    await t
                except (asyncio.CancelledError, CapExceeded, Exception):
                    # Terminal states captured above; cancel hygiene here.
                    pass
        finally:
            # Stop the OI worker (if still alive) and give it a short bounded
            # grace to flush the current snapshot; pending symbols may be
            # dropped at shutdown (counted in stats).
            self.stop.set()
            for t in (locals().get("oi_task"), locals().get("connect_task")):
                if t is not None and not t.done():
                    t.cancel()
                    try:
                        await t
                    except (asyncio.CancelledError, CapExceeded, Exception):
                        pass
            cap_reason = nonlocal_cap_reason[0]
            reason = ("cap" if cap_reason
                      else "sigterm" if self.stop.is_set()
                      else "max-seconds" if self.max_seconds else "manual")
            if cap_reason:
                # Cap terminal: sink is already closed; NO further capped
                # writes through it (cap report goes to stderr + exit code,
                # not through the sink).
                self.fsync_close()
                self.save_state()
                print(f"CAP_EXCEEDED: {cap_reason}; fail-closed shutdown "
                      "(no data deleted)", file=sys.stderr)
                self.exit_code = 3
            else:
                # Normal stop: the session_stop control still goes through the
                # cap-checked sink (write_line reopens/rotates if needed). If
                # THAT write trips a newly-exhausted cap (near-cap SIGTERM),
                # fail closed identically: close + save + exit 3 out-of-band.
                try:
                    self.emit_control("session_stop", reason=reason)
                    self.fsync_close()
                    self.save_state()
                    self.exit_code = 0
                except CapExceeded as exc:
                    self.fsync_close()
                    self.save_state()
                    print(f"CAP_EXCEEDED during shutdown control: {exc}; "
                          "fail-closed (no data deleted)", file=sys.stderr)
                    self.exit_code = 3
        return self.exit_code


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, required=True,
                   help="capture output directory (raw, git-ignored)")
    p.add_argument("--url", default=DEFAULT_URL,
                   help="routed public websocket endpoint (/market path)")
    p.add_argument("--streams", default=", ".join(DEFAULT_STREAMS),
                   help="comma-separated stream names")
    p.add_argument("--max-seconds", type=float, default=None,
                   help="bounded smoke duration; omit for persistent capture")
    p.add_argument("--max-bytes-per-day", type=int, default=DEFAULT_MAX_BYTES_PER_DAY)
    p.add_argument("--max-total-bytes", type=int, default=DEFAULT_MAX_TOTAL_BYTES)
    args = p.parse_args()
    args.streams = [s.strip() for s in args.streams.split(",") if s.strip()]
    return args


def main() -> int:
    args = parse_args()
    collector = Collector(args.out, args.url, args.streams, args.max_seconds,
                          args.max_bytes_per_day, args.max_total_bytes)
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    for sig in (signal.SIGTERM, signal.SIGINT):
        try:
            loop.add_signal_handler(sig, collector.request_stop, sig)
        except NotImplementedError:  # pragma: no cover
            signal.signal(sig, collector.request_stop)
    try:
        # run() handles every terminal condition internally: caps return 3,
        # pre-existing exhaustion returns 4, SIGTERM/max-seconds return 0.
        return loop.run_until_complete(collector.run())
    finally:
        loop.close()


if __name__ == "__main__":
    raise SystemExit(main())
