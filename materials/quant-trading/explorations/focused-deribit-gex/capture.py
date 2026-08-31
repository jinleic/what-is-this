#!/usr/bin/env python3
"""Deribit public options capture collector (BTC + ETH).

Research-only; uses ONLY public unauthenticated Deribit endpoints:
  WS   wss://www.deribit.com/ws/api/v2
  REST https://www.deribit.com/api/v2

Semantics source-verified on 2026-08-30 (see explorations/focused-deribit-gex/
source_verification.md, section 8 for the corrective round):
  - trades.option.BTC.100ms / trades.option.ETH.100ms (canonical dot form; the
    dash form `trades.option-btc_usd.100ms` acks with an EMPTY result array and
    subscribes nothing — never used): per-trade pushes with `direction` = taker
    direction (`buy`/`sell`), `amount` in coins (options), `price` in BTC/ETH
    (inverse option premium), `index_price` USD, `iv` %, `trade_id` (unique per
    currency), `trade_seq` (per instrument). `.raw` interval is NOT available
    unauthenticated (error 13778); 100ms is the finest public aggregation.
  - ticker.{instrument}.agg2 : full ticker incl. options `greeks`
    {delta,gamma,theta,vega,rho} (standard Black-Scholes, gamma per 1 USD
    underlying move), `index_price`, `underlying_price`, `mark_iv`, `open_interest`
    in coins, `timestamp` ms. CHAIN-WIDE TICKERS ARE OPT-IN ONLY
    (--capture-tickers) and are NOT part of the sustainable default.
  - public/get_instruments seeds the option chain (strike, option_type,
    expiration_timestamp, contract_size(=1.0 coin), lot_size, tick_size,
    settlement_currency, price_index). Day-scoped persistence: every new UTC
    capture day re-persists the chain so BS reconstruction always has metadata;
    reconnects inside a day never rewrite unchanged rows.
  - deribit_price_index.{index}: pushed index price also recorded as its own record.

Lifecycle (review-hardened):
  - append-only JSONL under data/raw/focused/live/deribit/dayYYYYMMDD/
    (at-least-once raw writes; REPLAY dedup-by-deterministic-id is authoritative)
  - deterministic dedup ids (sha256 of channel+ids) recorded with every record
  - sink rotates BEFORE every write whenever the capture UTC day changes; the
    rotation-open record itself is checked against both caps
  - subscribe acks REQUIRE exact returned==requested set equality; empty,
    partial, superset, or errored acks fail readiness and trigger a reconnect
    with capped, stop-aware backoff. DERIBIT_CAPTURE_READY is printed only
    after ALL exact acks plus a successful set_heartbeat ack.
  - reconnect backoff resets ONLY after a stable ready session, not on TCP setup
  - every wait (reconnect sleep, REST retry sleep, ack waits) is SIGTERM-bounded
    (<= 0.1 s wake); 256 MiB per-UTC-day + 10 GiB cumulative caps fail closed;
    prior capture is NEVER auto-deleted
  - SIGTERM/SIGINT: bounded clean stop; files flushed; exit 0.

Deterministic dedup id examples:
  trade:      sha256("trade|"+instrument_name+"|"+trade_id+"|"+trade_seq+"|"+timestamp)
  ticker:     sha256("ticker|"+instrument_name+"|"+timestamp)
  index:      sha256("index|"+index_name+"|"+timestamp)
  instrument: sha256("instr|"+instrument_name+"|"+creation_timestamp+"|"+state_ts)

No orders, no auth, no private endpoints. This file must remain stdlib-only.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import threading
import signal
import socket
import ssl
import struct
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

WS_HOST = "www.deribit.com"
WS_PATH = "/ws/api/v2/"
REST_BASE = "https://www.deribit.com/api/v2"

CURRENCIES = ("btc", "eth")  # lowercase for channel names

READY_BANNER = "DERIBIT_CAPTURE_READY"
SHUTDOWN_BANNER = "DERIBIT_CAPTURE_STOPPED"

OP_CONT, OP_TEXT, OP_CLOSE, OP_PING, OP_PONG = 0x0, 0x1, 0x8, 0x9, 0xA


def utc_day(ts_ms: int) -> str:
    return datetime.fromtimestamp(ts_ms / 1000.0, tz=timezone.utc).strftime("%Y%m%d")


def utc_iso(ts_ms: int) -> str:
    return datetime.fromtimestamp(ts_ms / 1000.0, tz=timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S.%f"
    )[:-3] + "Z"


def sha256_id(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class DayCapError(RuntimeError):
    """Raised when the per-UTC-day (or cumulative) byte cap is exhausted."""


class WsClient:
    """Minimal RFC6455 client (stdlib-only) for one WS connection."""

    def __init__(self, host: str = WS_HOST, path: str = WS_PATH, timeout: float = 20.0):
        self.host = host
        self.path = path
        self.timeout = timeout
        self._connect()

    def _connect(self):
        raw = socket.create_connection((self.host, 443), timeout=self.timeout)
        try:
            sock = ssl.create_default_context().wrap_socket(raw, server_hostname=self.host)
        except Exception:
            raw.close()
            raise
        self.sock = sock
        key = base64.b64encode(os.urandom(16)).decode()
        req = (
            f"GET {self.path} HTTP/1.1\r\n"
            f"Host: {self.host}\r\n"
            "Upgrade: websocket\r\n"
            "Connection: Upgrade\r\n"
            f"Sec-WebSocket-Key: {key}\r\n"
            "Sec-WebSocket-Version: 13\r\n"
            "\r\n"
        )
        sock.sendall(req.encode())
        resp = b""
        while b"\r\n\r\n" not in resp:
            chunk = sock.recv(4096)
            if not chunk:
                raise ConnectionError("handshake closed")
            resp += chunk
        status_line = resp.split(b"\r\n", 1)[0]
        head = resp.split(b"\r\n\r\n", 1)[0]
        if b"101" not in status_line:
            sock.close()
            raise ConnectionError(f"ws handshake failed: {status_line[:120]!r}")
        self._leftover = resp[len(head) + 4:]

    def send_text(self, text: str) -> None:
        self._send_frame(OP_TEXT, text.encode("utf-8"))

    def _send_frame(self, opcode: int, payload: bytes) -> None:
        mask = os.urandom(4)
        hdr = bytearray([0x80 | opcode])
        n = len(payload)
        if n < 126:
            hdr.append(0x80 | n)
        elif n < 65536:
            hdr.append(0x80 | 126)
            hdr += struct.pack(">H", n)
        else:
            hdr.append(0x80 | 127)
            hdr += struct.pack(">Q", n)
        hdr += mask
        masked = bytes(b ^ mask[i % 4] for i, b in enumerate(payload))
        self.sock.sendall(hdr + masked)


    def recv_message(self, timeout: float) -> tuple[int, bytes]:
        """Read one complete WS message; leftovers persist across calls."""
        # Persistent parse buffer: bytes belonging to the next frame must survive
        # every return path, else the stream desynchronizes.
        data = getattr(self, "_buf", b"")
        frag = getattr(self, "_frag", None)
        self.sock.settimeout(timeout)
        while True:
            while len(data) < 2:
                chunk = self.sock.recv(65536)
                if not chunk:
                    raise ConnectionError("closed")
                data += chunk
            op = data[0] & 0x0F
            fin = bool(data[0] & 0x80)
            masked = bool(data[1] & 0x80)
            ln = data[1] & 0x7F
            pos = 2
            if ln == 126:
                while len(data) < 4:
                    data += self.sock.recv(65536)
                ln = struct.unpack(">H", data[2:4])[0]
                pos = 4
            elif ln == 127:
                while len(data) < 10:
                    data += self.sock.recv(65536)
                ln = struct.unpack(">Q", data[2:10])[0]
                pos = 10
            need = pos + (4 if masked else 0) + ln
            while len(data) < need:
                chunk = self.sock.recv(65536)
                if not chunk:
                    raise ConnectionError("closed mid-frame")
                data += chunk
            if masked:
                mask = data[pos:pos + 4]
                body = bytes(
                    b ^ mask[i % 4] for i, b in enumerate(data[pos + 4:pos + 4 + ln])
                )
            else:
                body = data[pos:pos + ln]
            data = data[need:]
            self._buf = data
            if op == OP_PING:
                self._send_frame(OP_PONG, body)
                continue
            if op == OP_PONG:
                continue
            if op == OP_CLOSE:
                self._frag = None
                return OP_CLOSE, body
            if op == OP_TEXT and frag is None:
                if fin:
                    self._frag = None
                    self._buf = data
                    return OP_TEXT, body
                frag = [body]
                self._frag = frag
                continue
            if op == OP_CONT and frag is not None:
                frag.append(body)
                if fin:
                    merged = b"".join(frag)
                    self._frag = None
                    self._buf = data
                    return OP_TEXT, merged
                continue
            self._frag = None
            self._buf = data
            # binary or unknown opcode: skip payload and continue

    def close(self) -> None:
        try:
            self.sock.close()
        except Exception:
            pass


_STOP_REQUESTED = threading.Event()


def request_stop() -> None:
    """Mark the process stop-requested; makes every wait stop-aware so SIGTERM
    latency stays bounded (no multi-minute blocking sleeps)."""
    _STOP_REQUESTED.set()


def stop_requested() -> bool:
    return _STOP_REQUESTED.is_set()


def stop_aware_sleep(seconds: float) -> None:
    """Sleep in <=0.1 s slices so SIGTERM lands within 0.1 s of being set."""
    deadline = time.monotonic() + max(seconds, 0.0)
    while time.monotonic() < deadline and not stop_requested():
        time.sleep(min(0.1, max(0.0, deadline - time.monotonic())))


def rest_call(path: str, params: dict | None = None, attempts: int = 3) -> dict:
    query = ""
    if params:
        query = "?" + "&".join(f"{k}={v}" for k, v in params.items())
    url = f"{REST_BASE}{path}{query}"
    last: Exception | None = None
    for attempt in range(attempts):
        if stop_requested():
            raise ConnectionError("rest_call aborted: stop requested")
        try:
            req = Request(url, headers={"User-Agent": "quant-research-gex-capture/1.0"})
            with urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001 — collector must self-heal
            last = exc
            stop_aware_sleep(1.5 * (attempt + 1))
    raise ConnectionError(f"rest_call failed {path}: {last}")


def fetch_instruments(currency: str) -> list[dict]:
    payload = rest_call(
        "/public/get_instruments", {"currency": currency.upper(), "kind": "option", "expired": "false"}
    )
    if "result" not in payload:
        raise ConnectionError(f"get_instruments error: {payload.get('error')}")
    return payload["result"]


class Sink:
    """Append-only JSONL sink with startup-reconciled byte accounting,
    daily rotation, and fail-closed caps. Prior capture is never deleted."""

    def __init__(self, root: Path, day_cap_bytes: int, total_cap_bytes: int):
        self.root = root
        self.day_cap = day_cap_bytes
        self.total_cap = total_cap_bytes
        self.day_key: str | None = None
        self.day_bytes = 0          # exact bytes on disk for the ACTIVE day
        self.prior_bytes = 0        # exact bytes of every other (finished) day on disk
        self._handles: dict[str, object] = {}
        self.reconcile_existing()

    def _paths(self, day_key: str) -> dict[str, Path]:
        base = self.root / f"day{day_key}"
        base.mkdir(parents=True, exist_ok=True)
        return {
            "trades": base / "trades.jsonl",
            "ticker": base / "ticker.jsonl",
            "index": base / "index.jsonl",
            "instruments": base / "instruments.jsonl",
            "events": base / "events.jsonl",
        }

    def reconcile_existing(self) -> None:
        """startup reconciliation: prior_bytes = every day dir EXCEPT today's;
        today's existing bytes are loaded on first rotate via _rotate."""
        self.prior_bytes = 0
        self.day_key = None
        self.day_bytes = 0
        today = utc_day(int(time.time() * 1000))
        if self.root.is_dir():
            for entry in sorted(self.root.iterdir()):
                if not (entry.is_dir() and entry.name.startswith("day")):
                    continue
                if entry.name == f"day{today}":
                    continue  # handled by first _rotate
                for f in entry.iterdir():
                    try:
                        self.prior_bytes += f.stat().st_size
                    except OSError:
                        pass

    def write(self, kind: str, record: dict, now_ms: int) -> None:
        if kind not in ("trades", "ticker", "index", "instruments", "events"):
            raise ValueError(kind)
        target_day = utc_day(now_ms)
        # rotate BEFORE every write when the capture UTC day changes (midnight etc.)
        if self.day_key is None or self.day_key != target_day:
            self._rotate(target_day)
        payload = (json.dumps(record, separators=(",", ":")) + "\n").encode("utf-8")
        n = len(payload)
        # cumulative cap: all bytes on disk (finished days + active day + this record)
        if self.total_cap is not None:
            if self.prior_bytes + self.day_bytes + n > self.total_cap:
                raise DayCapError(
                    f"cumulative cap {self.total_cap} bytes exceeded: "
                    f"disk={self.prior_bytes} active={self.day_bytes} new={n}"
                )
        if self.day_bytes + n > self.day_cap:
            raise DayCapError(
                f"daily cap {self.day_cap} bytes exceeded for UTC day {self.day_key}"
            )
        h = self._handles[kind]
        h.write(payload)
        h.flush()
        os.fsync(h.fileno())
        self.day_bytes += n

    def _rotate(self, day_key: str) -> None:
        # fold current active-day bytes into prior pool (files remain on disk)
        self.prior_bytes += self.day_bytes
        for h in self._handles.values():
            try:
                h.flush()
                os.fsync(h.fileno())
            except Exception:
                pass
            try:
                h.close()
            except Exception:
                pass
        self._handles = {}
        self.day_key = day_key
        self.day_bytes = 0
        paths = self._paths(day_key)
        for kind, path in paths.items():
            # count existing bytes of the new active day on disk (restart-safe)
            try:
                self.day_bytes += path.stat().st_size
            except OSError:
                pass
            self._handles[kind] = open(path, "ab")
        rec = {
            "record_type": "capture_event",
            "dedup_id": sha256_id(f"event|rotate|{day_key}|{paths['trades'].name}"),
            "event": "day_rotation_open",
            "day_utc": day_key,
            "received_ts_ms": int(time.time() * 1000),
        }
        ev_payload = (json.dumps(rec, separators=(",", ":")) + "\n").encode("utf-8")
        # the rotation event itself must respect BOTH caps (restart-safety)
        if self.day_bytes + len(ev_payload) > self.day_cap:
            raise DayCapError(
                f"daily cap {self.day_cap} bytes exceeded for UTC day {self.day_key} "
                f"(pre-existing {self.day_bytes} bytes on restart)"
            )
        if self.total_cap is not None and (
            self.prior_bytes + self.day_bytes + len(ev_payload) > self.total_cap
        ):
            raise DayCapError(
                f"cumulative cap {self.total_cap} bytes exceeded "
                f"(pre-existing disk={self.prior_bytes} active={self.day_bytes})"
            )
        h = self._handles["events"]
        h.write(ev_payload)
        h.flush()
        os.fsync(h.fileno())
        self.day_bytes += len(ev_payload)

    def close(self) -> None:
        for h in self._handles.values():
            try:
                h.flush()
                os.fsync(h.fileno())
            except Exception:
                pass
            try:
                h.close()
            except Exception:
                pass
        self._handles = {}


class Deduper:
    """Deterministic dedup: set of ids persisted per process run; the REPLAY
    (run.py) is where cross-session dedup is authoritative via persisted id files.
    The collector keeps a bounded in-memory recent-ids LRU to drop obvious dupes
    from quick reconnect replays without unbounded memory."""

    def __init__(self, max_ids: int = 500_000):
        self.seen: set[str] = set()
        self.order: list[str] = []
        self.max_ids = max_ids
        self.duplicates = 0

    def check(self, dedup_id: str) -> bool:
        if dedup_id in self.seen:
            self.duplicates += 1
            return False
        self.seen.add(dedup_id)
        self.order.append(dedup_id)
        if len(self.order) > self.max_ids // 10:
            for _ in range(self.max_ids // 20):
                old = self.order.pop(0)
                self.seen.discard(old)
        return True


class Collector:
    def __init__(self, raw_dir: Path, currencies=("btc", "eth"), day_cap=256 * 1024 * 1024,
                 total_cap=10 * 1024 * 1024 * 1024, heartbeat=True, max_minutes=None,
                 capture_tickers=False):
        self.raw_dir = raw_dir
        self.currencies = currencies
        self.day_cap = day_cap
        self.total_cap = total_cap
        self.instruments: dict[str, dict] = {}
        self._persisted_instruments: set[str] = set()
        self.heartbeat = heartbeat
        self.max_minutes = max_minutes
        self.capture_tickers = capture_tickers
        self.stop = False
        self.sink = Sink(raw_dir, day_cap, total_cap)
        self.deduper = Deduper()
        self.active_channels: set[str] = set()
        self.ws: WsClient | None = None
        self.req_id = 0
        self.msg_seq = 0
        self._ready_printed = False
        self.last_message_ms = 0
        self._started = time.time()
        self.last_utc_day = None
        self._stable_session = False
        self.stats = {
            "trades_pushes": 0, "trades_records": 0, "dup_trades": 0,
            "ticker_pushes": 0, "ticker_records": 0, "index_records": 0,
            "instrument_records": 0, "reconnects": 0, "events": 0,
            "primary_subscribe_ok": False,
        }

    # ---------- restart-state rehydration ----------
    def rehydrate_state(self, recent_trade_ids: int = 20_000) -> dict:
        """Restore cross-process state from the ACTIVE UTC day's append-only files
        BEFORE the first REST seed or subscription:

        - `_persisted_instruments` from instruments.jsonl so a restart re-writes
          zero unchanged instrument rows.
        - the trade Deduper from the tail of trades.jsonl so a restart cannot
          duplicate a trade already captured today (bounded to the newest
          `recent_trade_ids` entries to keep memory bounded).

        Read-only; never mutates or deletes prior capture.
        """
        stats = {"instrument_ids": 0, "trade_ids": 0, "day": None, "malformed": 0}
        day = utc_day(int(time.time() * 1000))
        stats["day"] = day
        day_dir = self.raw_dir / f"day{day}"
        instr_path = day_dir / "instruments.jsonl"
        trades_path = day_dir / "trades.jsonl"
        if instr_path.is_file():
            with instr_path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        stats["malformed"] += 1
                        continue
                    if rec.get("record_type") != "instrument":
                        continue
                    iname = (rec.get("instrument") or {}).get("instrument_name")
                    if iname:
                        self._persisted_instruments.add(iname)
            stats["instrument_ids"] = len(self._persisted_instruments)
        if trades_path.is_file():
            tail: list[str] = []
            with trades_path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        stats["malformed"] += 1
                        continue
                    did = rec.get("dedup_id")
                    if did:
                        tail.append(did)
                        if len(tail) > recent_trade_ids:
                            del tail[: len(tail) - recent_trade_ids]
            for did in tail:
                self.deduper.check(did)
            stats["trade_ids"] = len(tail)
        # a same-day restart must NOT re-treat the just-loaded ids as stale:
        # mark today's instrument book as already persisted for this capture day
        self._instrument_persist_day = day
        self._rehydrated = stats
        return stats

    # ---------- signal handling ----------
    def _handle_signal(self, signum, _frame):
        self.stop = True
        request_stop()

    # ---------- subscriptions ----------
    def _next_id(self) -> int:
        self.req_id += 1
        return self.req_id

    def default_channels(self) -> list[str]:
        """Sustainable default: per-currency option trades + price index only.
        Official doc syntax: trades.option.BTC.100ms (deribit docs), index push
        deribit_price_index.{cur}. Chain-wide tickers are opt-in only."""
        channels = []
        for cur in self.currencies:
            C = cur.upper()
            channels.append(f"trades.option.{C}.100ms")
            channels.append(f"deribit_price_index.{cur}_usd")
        return channels

    def _ticker_channels(self) -> list[str]:
        return [f"ticker.{iname}.agg2" for iname in sorted(self.instruments)
                if iname.startswith(("BTC-", "ETH-"))]

    def _await_ack(self, rid: int, chunk: set[str], timeout_s: float = 20.0) -> bool:
        """Wait for the JSON-RPC ack of one request. Succeeds ONLY on exact set
        equality: returned == requested (empty/partial/superset all fail).
        Stop-aware: SIGTERM aborts the wait within <=0.1 s."""
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            if stop_requested():
                self._log_event("subscribe_aborted_stop")
                return False
            try:
                op, data = self.ws.recv_message(timeout=min(2.0, max(0.1, deadline - time.time())))
            except socket.timeout:
                continue
            if op == OP_CLOSE:
                raise ConnectionError("closed during subscribe")
            try:
                obj = json.loads(data.decode())
            except Exception:
                continue
            if obj.get("id") != rid:
                continue
            if "error" in obj:
                self._log_event("subscribe_error", {"error": obj["error"], "n": len(chunk)})
                return False
            returned = set(obj.get("result") or [])
            if returned != chunk:
                self._log_event("subscribe_set_mismatch", {
                    "missing": sorted(chunk - returned)[:20],
                    "extra": sorted(returned - chunk)[:20],
                    "result_empty": "result_empty" if not returned else None,
                })
                return False
            self.active_channels |= chunk
            return True
        self._log_event("subscribe_timeout", {"n": len(chunk)})
        return False

    def subscribe_all(self) -> bool:
        """Subscribe in deterministic batches; readiness requires ALL batches to
        return exact channel-set equality AND a successful heartbeat ack. Any
        failure or stop request -> False (caller closes socket and retries with
        capped backoff)."""
        channels = self.default_channels()
        if self.capture_tickers:
            channels += self._ticker_channels()
        CHUNK = 220
        self.active_channels = set()
        ok_all = True
        for i in range(0, len(channels), CHUNK):
            if stop_requested():
                return False
            chunk = channels[i:i + CHUNK]
            rid = self._next_id()
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0", "id": rid, "method": "public/subscribe",
                "params": {"channels": chunk},
            }))
            if not self._await_ack(rid, set(chunk)):
                ok_all = False
                break
        if ok_all and self.heartbeat:
            if stop_requested():
                return False
            hb_id = self._next_id()
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0", "id": hb_id, "method": "public/set_heartbeat",
                "params": {"interval": 10},
            }))
            deadline = time.time() + 15
            hb_ok = False
            while time.time() < deadline and not hb_ok:
                if stop_requested():
                    return False
                try:
                    op, data = self.ws.recv_message(timeout=min(2.0, max(0.1, deadline - time.time())))
                except socket.timeout:
                    continue
                if op == OP_CLOSE:
                    raise ConnectionError("closed during heartbeat setup")
                try:
                    obj = json.loads(data.decode())
                except Exception:
                    continue
                if obj.get("id") != hb_id:
                    continue
                if "error" in obj:
                    self._log_event("heartbeat_setup_error", {"error": obj["error"]})
                    break
                if obj.get("result") == "ok":
                    hb_ok = True
            if not hb_ok:
                self._log_event("heartbeat_setup_timeout")
                ok_all = False
        # final exact verification of the whole held set
        if ok_all and self.active_channels != set(channels):
            self._log_event("subscribe_set_mismatch_final", {
                "missing": sorted(set(channels) - self.active_channels)[:20],
                "extra": sorted(self.active_channels - set(channels))[:20]})
            ok_all = False
        return ok_all

    def seed_instruments(self, received_ms: int) -> None:
        """Fetch chain per currency, update the in-memory map, and write rows not
        yet persisted for the CURRENT capture day. Reconnects within a day write
        nothing when nothing changed. When the capture UTC day rolls over, the
        sink rotates before the write and `_persisted_instruments` is reset so
        the new day's file always carries the strike/expiry/contract_size fields
        needed for BS reconstruction in run.py (day-scoped metadata)."""
        if stop_requested():
            raise ConnectionError("stop requested before instrument seed")
        target_day = utc_day(received_ms)
        if target_day != getattr(self, "_instrument_persist_day", None):
            # new capture day: forget every instrument as persisted
            self._persisted_instruments = set()
            self._instrument_persist_day = target_day
        for cur in self.currencies:
            if stop_requested():
                self._log_event("instrument_seed_aborted_stop", {"currency": cur})
                return
            try:
                rows = fetch_instruments(cur)
            except ConnectionError as exc:
                if stop_requested():
                    self._log_event("instrument_seed_aborted_stop", {"currency": cur})
                    return
                self._log_event("instrument_seed_error", {"currency": cur, "error": str(exc)})
                continue
            except Exception as exc:  # noqa: BLE001
                self._log_event("instrument_seed_error", {"currency": cur, "error": str(exc)})
                continue
            new_rows = []
            for row in rows:
                iname = row["instrument_name"]
                self.instruments[iname] = row
                if iname not in self._persisted_instruments:
                    dedup = self._instrument_dedup_id(row)
                    new_rows.append((dedup, row))
                    self._persisted_instruments.add(iname)
            for dedup, row in new_rows:
                rec = {
                    "record_type": "instrument",
                    "dedup_id": dedup,
                    "received_ts_ms": received_ms,
                    "instrument": row,  # strike, option_type, expiration_timestamp,
                                        # contract_size, lot_size, tick_size, ...
                }
                try:
                    self.sink.write("instruments", rec, now_ms=received_ms)
                except DayCapError:
                    self._persisted_instruments.discard(row["instrument_name"])
                    raise
                self.stats["instrument_records"] += 1
            self.stats["instrument_seed_new"] = self.stats.get("instrument_seed_new", 0) + len(new_rows)
        self._log_event("instrument_seed_done", {"batches": list(self.currencies),
                                                 "n": len(self.instruments)})

    def _log_event(self, event: str, detail: dict | None = None) -> None:
        now_ms = int(time.time() * 1000)
        rec = {
            "record_type": "capture_event",
            "dedup_id": sha256_id(f"event|{event}|{now_ms}|{self.msg_seq}"),
            "event": event,
            "received_ts_ms": now_ms,
        }
        if detail:
            rec["detail"] = detail
        try:
            self.sink.write("events", rec, now_ms=now_ms)
        except DayCapError:
            raise
        self.stats["events"] += 1

    def ingest_trade(self, trade: dict, received_ms: int) -> None:
        iname = trade.get("instrument_name", "")
        dedup = sha256_id(
            f"trade|{iname}|{trade.get('trade_id')}|{trade.get('trade_seq')}|{trade.get('timestamp')}"
        )
        if not self._dedupe_and_count("trades", dedup):
            return
        rec = {
            "record_type": "option_trade",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "trade": trade,  # verbatim: direction(taker), amount(coins), price(BTC/ETH premium),
                             # index_price(USD), iv(%), contract info joined at replay from instruments
        }
        self.sink.write("trades", rec, now_ms=received_ms)
        self.stats["trades_records"] += 1

    def ingest_ticker(self, iname: str, data: dict, received_ms: int) -> None:
        ts = data.get("timestamp")
        if ts is None:
            return
        dedup = sha256_id(f"ticker|{iname}|{ts}")
        if not self._dedupe_and_count("ticker", dedup):
            return
        rec = {
            "record_type": "option_ticker",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "instrument_name": iname,
            "ticker": data,  # verbatim: greeks{gamma,...} etc.
        }
        self.sink.write("ticker", rec, now_ms=received_ms)
        self.stats["ticker_records"] += 1

    def ingest_index(self, index_name: str, data: dict, received_ms: int) -> None:
        ts = data.get("timestamp")
        if ts is None:
            return
        dedup = sha256_id(f"index|{index_name}|{ts}")
        if not self._dedupe_and_count("index", dedup):
            return
        rec = {
            "record_type": "index_price",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "index_name": index_name,
            "data": data,
        }
        self.sink.write("index", rec, now_ms=received_ms)
        self.stats["index_records"] += 1

    def _dedupe_and_count(self, kind: str, dedup: str) -> bool:
        ok = self.deduper.check(dedup)
        if not ok:
            self._log_event("duplicate_suppressed", {"dedup_id": dedup, "kind": kind})
        return ok

    def _instrument_dedup_id(self, row: dict) -> str:
        return sha256_id(f"instr|{row['instrument_name']}|{row.get('creation_timestamp')}|seed")


    def reconcile_chain(self, received_ms: int) -> None:
        """Resubscribe tickers for new instruments ONLY in ticker-capture mode."""
        if not self.capture_tickers:
            return
        want = set()
        for iname in self.instruments:
            if iname.startswith(("BTC-", "ETH-")):
                want.add(f"ticker.{iname}.agg2")
        have = {c for c in self.active_channels if c.startswith("ticker.")}
        add = want - have
        remove = have - want
        if add:
            self._modify(list(add), "public/subscribe")
            self._log_event("tickers_added", {"n": len(add)})
        if remove:
            self._modify(list(remove), "public/unsubscribe")
            self._log_event("tickers_removed", {"n": len(remove)})

    def run(self) -> int:
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)
        backoff = 1.0
        started = time.time()
        try:
            # ---- startup lifecycle (rehydrate + seed) shares the same
            # stop/cap/error handling as the session loop: a SIGTERM arriving
            # during rehydrate/REST seed aborts here with no traceback, the
            # sink is closed exactly once, and DERIBIT_CAPTURE_STOPPED prints.
            while not self.stop:
                if stop_requested():
                    break
                try:
                    rehydrated = self.rehydrate_state()
                    self._log_event("restart_state_rehydrated", rehydrated)
                    self.seed_instruments(int(time.time() * 1000))
                    break  # startup complete; fall through to the session loop
                except DayCapError as exc:
                    self.stats["cap_caught"] = True
                    print(f"DERIBIT_CAPTURE_CAP_EXCEEDED {exc}", flush=True)
                    return 0
                except (ConnectionError, socket.timeout, ssl.SSLError, OSError) as exc:
                    if stop_requested():
                        break
                    self.stats["reconnects"] += 1
                    self._log_event("startup_retry_backoff",
                                    {"error": f"{type(exc).__name__}: {exc}",
                                     "backoff_s": backoff})
                    stop_aware_sleep(backoff)
                    backoff = min(backoff * 2.0, 60.0)

            # ---- session loop
            while not self.stop:
                if stop_requested():
                    break
                try:
                    self.ws = WsClient()
                    self.active_channels = set()
                    self.seed_instruments(int(time.time() * 1000))
                    ok = self.subscribe_all()
                    self._log_event("subscribe_ack", {"ok": ok,
                                                      "channels": sorted(self.active_channels),
                                                      "n_channels": len(self.active_channels)})
                    if not ok:
                        # failed/incomplete subscribe: never enter the receive
                        # loop; reconnect with capped, stop-aware backoff
                        raise ConnectionError(
                            "subscription ack not exact or heartbeat setup failed")
                    self._stable_session = True
                    if not self._ready_printed:
                        print(READY_BANNER, flush=True)
                        self._ready_printed = True
                    backoff = 1.0  # reset only after a stable ready session
                    self.last_chain_reconcile = time.time()
                    self._loop()
                except DayCapError as exc:
                    self.stats["cap_caught"] = True
                    print(f"DERIBIT_CAPTURE_CAP_EXCEEDED {exc}", flush=True)
                    break
                except (ConnectionError, socket.timeout, ssl.SSLError, OSError) as exc:
                    self._stable_session = False
                    self.stats["reconnects"] += 1
                    try:
                        self._log_event("reconnect_backoff",
                                        {"error": f"{type(exc).__name__}: {exc}",
                                         "backoff_s": backoff})
                    except DayCapError:
                        self.stats["cap_caught"] = True
                        print("DERIBIT_CAPTURE_CAP_EXCEEDED during reconnect logging",
                              flush=True)
                        break
                    if stop_requested():
                        break
                    stop_aware_sleep(backoff)
                    backoff = min(backoff * 2.0, 60.0)
                finally:
                    if self.ws is not None:
                        self.ws.close()
                        self.ws = None
                if self.max_minutes is not None and time.time() - started > self.max_minutes * 60:
                    self.stop = True
        finally:
            # single sink close for EVERY exit path (startup, loop, cap, stop)
            self.sink.close()
            print(SHUTDOWN_BANNER, flush=True)
        return 0

    def _modify(self, channels: list[str], method: str) -> None:
        CHUNK = 220
        for i in range(0, len(channels), CHUNK):
            if stop_requested():
                raise ConnectionError("stop requested during modify")
            chunk = channels[i:i + CHUNK]
            rid = self._next_id()
            self.ws.send_text(json.dumps({"jsonrpc": "2.0", "id": rid, "method": method,
                                          "params": {"channels": chunk}}))
            deadline = time.time() + 15
            done = False
            while time.time() < deadline and not done:
                if stop_requested():
                    raise ConnectionError("stop requested during modify")
                try:
                    op, data = self.ws.recv_message(timeout=min(2.0, max(0.1, deadline - time.time())))
                except socket.timeout:
                    continue
                if op == OP_CLOSE:
                    raise ConnectionError("closed during modify")
                try:
                    obj = json.loads(data.decode())
                except Exception:
                    continue
                if obj.get("id") == rid:
                    if "error" in obj:
                        self._log_event("modify_error", {"method": method,
                                                         "error": obj["error"]})
                        break
                    if method == "public/subscribe":
                        if set(obj.get("result") or []) == set(chunk):
                            self.active_channels |= set(chunk)
                        else:
                            self._log_event("modify_set_mismatch", {"expected": len(chunk),
                                                                    "returned": len(obj.get("result") or [])})
                    elif method == "public/unsubscribe":
                        for ch in chunk:
                            self.active_channels.discard(ch)
                    done = True

    def _loop(self) -> None:
        idle_rounds = 0
        while not self.stop:
            received_ms = int(time.time() * 1000)
            # daily rotation + midnight chain refresh
            self._maybe_refresh_chain(received_ms)
            try:
                op, data = self.ws.recv_message(timeout=2.0)
            except socket.timeout:
                idle_rounds += 1
                if idle_rounds >= 8:
                    elapsed = int(time.time() * 1000) - self.last_message_ms
                    self._log_event("connection_idle_ping", {"idle_ms": min(elapsed, 2**31)})
                    idle_rounds = 0
                continue
            except (ConnectionError, OSError) as exc:
                raise ConnectionError(f"recv: {exc}")
            idle_rounds = 0
            if op == OP_CLOSE:
                raise ConnectionError("server closed")
            if op != OP_TEXT:
                continue
            try:
                obj = json.loads(data.decode())
            except json.JSONDecodeError:
                self._log_event("malformed_json", {"frame_bytes": len(data)})
                continue
            self.msg_seq += 1
            self.last_message_ms = int(time.time() * 1000)
            method = obj.get("method")
            if method == "subscription":
                params = obj.get("params", {})
                channel = params.get("channel", "")
                body = params.get("data")
                self._route(channel, body, self.last_message_ms)
            elif method == "heartbeat":
                self.ws.send_text(json.dumps({"jsonrpc": "2.0", "method": "public/test"}))
                self._log_event("heartbeat_ok", None)
            elif "error" in obj:
                self._log_event("server_error_notification", {"error": obj["error"]})
            # cap-safety: enforce budget in terms of schedule too
            if self.max_minutes is not None and time.time() - self._started > self.max_minutes * 60:
                self.stop = True

    def _maybe_refresh_chain(self, now_ms: int) -> None:
        day = utc_day(now_ms)
        if day != self.last_utc_day:
            self.last_utc_day = day
            self.seed_instruments(now_ms)
            if self.ws is not None and self.active_channels:
                self.reconcile_chain(now_ms)

    def _route(self, channel: str, body, received_ms: int) -> None:
        if channel.startswith("trades.option.") or channel.startswith("trades.option-"):
            self.stats["trades_pushes"] += 1
            if isinstance(body, list):
                for trade in body:
                    self.ingest_trade(trade, received_ms)
        elif channel.startswith("deribit_price_index."):
            idx_name = channel[len("deribit_price_index."):]
            if isinstance(body, dict):
                self.ingest_index(idx_name, body, received_ms)
            elif isinstance(body, list):
                for item in body:
                    self.ingest_index(idx_name, item, received_ms)
        elif channel.startswith("instrument.creation."):
            # same payload shape as get_instruments rows; dedup-before-write
            if isinstance(body, dict):
                row = body
                iname = row.get("instrument_name")
                if iname and iname not in self._persisted_instruments:
                    self.instruments[iname] = row
                    dedup = self._instrument_dedup_id(row)
                    rec = {"record_type": "instrument", "dedup_id": dedup,
                           "received_ts_ms": received_ms, "instrument": row}
                    try:
                        self.sink.write("instruments", rec, now_ms=received_ms)
                    except DayCapError:
                        raise
                    self._persisted_instruments.add(iname)
                    self.stats["instrument_records"] += 1
                    self.reconcile_chain(received_ms)
        elif channel.startswith("instrument.state."):
            if isinstance(body, dict):
                iname = body.get("instrument_name")
                ts = body.get("timestamp", received_ms)
                dedup = sha256_id(f"instr|state|{iname}|{ts}|{body.get('state')}")
                rec = {"record_type": "instrument_state", "dedup_id": dedup,
                       "received_ts_ms": received_ms, "state": body}
                self.sink.write("instruments", rec, now_ms=received_ms)
                self.stats["instrument_records"] += 1


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Deribit public options capture collector")
    p.add_argument("--raw-dir", type=Path,
                   default=Path("quant-trading/data/raw/focused/live/deribit"))
    p.add_argument("--currencies", default="BTC,ETH")
    p.add_argument("--max-day-bytes", type=int, default=256 * 1024 * 1024,
                   help="fail-closed per-UTC-day cap (default 256 MiB)")
    p.add_argument("--max-total-bytes", type=int, default=10 * 1024 * 1024 * 1024,
                   help="fail-closed cumulative cap (default 10 GiB); prior capture never deleted")
    p.add_argument("--max-minutes", type=float, default=None,
                   help="bounded smoke duration; unset = persistent")
    p.add_argument("--capture-tickers", action="store_true",
                   help="opt-in: subscribe ticker.{instrument}.agg2 for the full open "
                        "BTC+ETH option chain (heavy bandwidth; NOT default)")
    p.add_argument("--no-heartbeat", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    currencies = tuple(c.strip().lower() for c in args.currencies.split(",") if c.strip())
    for c in currencies:
        if c not in ("btc", "eth"):
            print(f"unsupported currency {c} (research scope: btc, eth)", file=sys.stderr)
            return 2
    coll = Collector(args.raw_dir, currencies, day_cap=args.max_day_bytes,
                     total_cap=args.max_total_bytes, heartbeat=not args.no_heartbeat,
                     max_minutes=args.max_minutes, capture_tickers=args.capture_tickers)
    return coll.run()


if __name__ == "__main__":
    raise SystemExit(main())


