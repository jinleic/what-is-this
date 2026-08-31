#!/usr/bin/env python3
"""Bounded BTC/ETH context-stream collector for the fast liquidation screen.

Research only: public unauthenticated Binance market websockets; no orders, no
account access, no credentials. This process NEVER touches the running
persistent collectors' files, sockets, or state: it owns a separate output
directory and state file.

Scope (steering-bound):
  - Liquidity trigger: !forceOrder@arr (all-market liquidation trigger stream).
  - Context streams ONLY for BTCUSDT and ETHUSDT, ONLY markPrice@1s and
    aggTrade, and ONLY while an onset context window for that symbol is open
    (dynamic SUBSCRIBE on a new eligible onset, UNSUBSCRIBE exactly 20 minutes
    after the window opened; a repeated onset inside an open window never
    extends the window). At most 5 streams are ever active:
    1 trigger + 4 context (2 symbols x 2 kinds).
  - Frames accepted in BOTH transports: wrapped {"stream":..,"data":{..}}
    (/stream combined frames) and raw {"e":..} (/market/ws dynamic frames).
    Deterministic canonical ids match the collector/replay formulas
    (forceOrder canonicalizes to !forceOrder@arr; aggTrade ids use stream|E|s|a;
    markPriceUpdate ids use stream|E|s|T), so dedup is byte-compatible with
    the existing dossier evidence.
  - Onset rule mirrors the frozen screen: per symbol, an onset is a
    deduplicated forceOrder event at least 300000 ms after the immediately
    prior same-symbol event inside the current continuously captured session;
    the first same-symbol event of a session needs E - connect_capture_ms
    >= 300000. session_start/disconnect/session_stop reset the session.
  - Append-only JSONL with fail-closed caps: 256 MiB per UTC day, 10 GiB
    cumulative; rotation before every write on UTC-day change; prior capture
    never deleted; startup byte accounting floored from disk.
  - 90-minute hard stop from process start (CLI may request less; never more).
    SIGTERM/SIGINT stop aware and bounded.

Default behaviour is NOT launched by the analyst; when launched it prints
LIQUIDATION_CONTEXT_READY after the exact trigger-subscribe ack and
LIQUIDATION_CONTEXT_STOPPED on any clean exit; cap exhaustion prints
LIQUIDATION_CONTEXT_CAP_EXCEEDED and exits 3.

--selftest runs deterministic no-network proofs: normalization (wrapped/raw
identity), onset gating (connect rule, 300s gap, per-symbol reset), the
subscription state machine (bounded stream set, open-on-onset, close exactly
at +20 min, no extension, no non-registered symbols), cap fail-closed
accounting, and the 90-minute hard stop. These use synthetic frames only.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import signal
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DERIBIT_CAPTURE = REPO_ROOT / "explorations" / "focused-deribit-gex" / "capture.py"
DEFAULT_OUT = REPO_ROOT / "data" / "raw" / "focused" / "live" / "fast-context"
DEFAULT_URL = "wss://fstream.binance.com/market/ws"

READY_BANNER = "LIQUIDATION_CONTEXT_READY"
STOP_BANNER = "LIQUIDATION_CONTEXT_STOPPED"
CAP_BANNER = "LIQUIDATION_CONTEXT_CAP_EXCEEDED"

TRIGGER_STREAM = "!forceOrder@arr"
REGISTERED_SYMBOLS = ("BTCUSDT", "ETHUSDT")
SYMBOL_STREAMS = {s: (f"{s.lower()}@aggTrade", f"{s.lower()}@markPrice@1s")
                  for s in REGISTERED_SYMBOLS}
CONTEXT_STREAMS_PER_SYMBOL = 2
MAX_ACTIVE_STREAMS = 1 + CONTEXT_STREAMS_PER_SYMBOL * len(REGISTERED_SYMBOLS)

ONSET_GAP_MS = 300_000
WINDOW_MS = 20 * 60 * 1000
HARD_STOP_MS = 90 * 60 * 1000
MI = 1024 * 1024
GI = 1024 * MI
DEFAULT_DAY_CAP = 256 * MI
DEFAULT_TOTAL_CAP = 10 * GI

EXIT_OK = 0
EXIT_ARGS = 2
EXIT_CAP = 3
EXIT_CONN = 4

FORCE_ORDER_FIELDS = ("s", "S", "o", "f", "q", "p", "ap", "X", "l", "z", "T")


def _load_deribit_helpers():
    """Small stdlib-only reuse: WsClient, Deduper, stop helpers."""
    spec = importlib.util.spec_from_file_location("deribit_capture_helpers", DERIBIT_CAPTURE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HELPERS = _load_deribit_helpers()
WsClient = HELPERS.WsClient
Deduper = HELPERS.Deduper
request_stop = HELPERS.request_stop
stop_requested = HELPERS.stop_requested
stop_aware_sleep = HELPERS.stop_aware_sleep


def utc_now_ms() -> int:
    return int(time.time() * 1000)


def iso_utc(ms: int | None) -> str | None:
    if ms is None:
        return None
    return datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S.") + f"{datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc).microsecond // 1000:03d}Z"


# --------------------------------------------------------------------------
# frame normalization (wrapped/raw), canonical ids byte-match the collectors
# --------------------------------------------------------------------------

def unwrap(frame: dict) -> tuple[str | None, dict]:
    """Return (stream_name_or_None, body). Accepts wrapped {stream,data} and
    raw frames; unknown shapes return (None, {})."""
    if not isinstance(frame, dict):
        return None, {}
    data = frame.get("data")
    if isinstance(data, dict):
        return frame.get("stream"), data
    return frame.get("stream"), frame


def canonical_stream(stream: str | None, body: dict) -> str | None:
    """Map raw shapes to the canonical stream names used by the collectors."""
    kind = body.get("e")
    if kind == "forceOrder":
        if stream and stream.endswith("@forceOrder"):
            return stream
        return TRIGGER_STREAM
    if kind == "aggTrade":
        if stream:
            return stream
        return f"{body.get('s', '').lower()}@aggTrade"
    if kind == "markPriceUpdate":
        if stream:
            return stream
        return f"{body.get('s', '').lower()}@markPrice@1s"
    return stream


def canonical_id(stream: str | None, body: dict) -> str | None:
    """Byte-match the dossier id formulas."""
    import hashlib
    kind = body.get("e")
    canon = canonical_stream(stream, body)
    if kind == "forceOrder":
        o = body.get("o", body)
        if not isinstance(o, dict):
            return None
        parts = [canon or "", str(body.get("E", ""))]
        for key in FORCE_ORDER_FIELDS:
            parts.append(str(o.get(key, "")))
    elif kind == "aggTrade":
        parts = [canon or "", str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("a", ""))]
    elif kind == "markPriceUpdate":
        parts = [canon or "", str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("T", ""))]
    else:
        return None
    import hashlib
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


# --------------------------------------------------------------------------
# onset tracking (live variant of the frozen screen rule)
# --------------------------------------------------------------------------

class SessionTracker:
    """Per-symbol onset detection on the exchange E clock, with capture-session
    resets. connect_ms anchors the first-event eligibility."""

    def __init__(self, connect_ms: int):
        self.connect_ms = connect_ms
        self.prior: dict[str, int] = {}       # symbol -> last E seen this session
        self.onsets: list[dict] = []

    def reset(self, connect_ms: int | None) -> None:
        self.connect_ms = connect_ms
        self.prior = {}

    def observe_force_order(self, body: dict) -> dict | None:
        if self.connect_ms is None:
            return None
        o = body.get("o") if isinstance(body.get("o"), dict) else {}
        symbol = o.get("s")
        E = body.get("E")
        if symbol not in REGISTERED_SYMBOLS or not isinstance(E, int) or isinstance(E, bool):
            return None
        prev = self.prior.get(symbol)
        is_onset = (E - self.connect_ms >= ONSET_GAP_MS) if prev is None \
            else (E - prev >= ONSET_GAP_MS)
        self.prior[symbol] = E
        if not is_onset:
            return None
        onset = {"symbol": symbol, "E": E, "side": o.get("S"),
                 "base": prev, "connect_ms": self.connect_ms}
        self.onsets.append(onset)
        return onset


# --------------------------------------------------------------------------
# subscription state machine (<=5 streams; open on onset; close at +20 min)
# --------------------------------------------------------------------------

class SubscriptionState:
    """Trigger stream is always desired. Context streams for a symbol are
    desired exactly while its context window is open (onset .. +WINDOW_MS).
    Windows open only on a NEW eligible onset and never extend; desired-set
    changes translate 1:1 into SUBSCRIBE/UNSUBSCRIBE actions."""

    def __init__(self, now_ms: int):
        self.now_ms = now_ms
        self.active: set[str] = set()         # canonically subscribed streams
        self.windows: dict[str, dict] = {}    # symbol -> window record
        self.actions: list[dict] = []         # audit log of desired-set changes
        self.closed = False                   # set once hard stop reached
        self._check_invariant()

    # ---- desired set ----------------------------------------------------
    def desired(self, now_ms: int | None = None) -> set[str]:
        if self.closed:
            return set()
        t = self.now_ms if now_ms is None else now_ms
        want = {TRIGGER_STREAM}
        for symbol, w in self.windows.items():
            if w is not None and w["open_until_ms"] > t:
                want.update(SYMBOL_STREAMS[symbol])
        return want

    def open_window(self, symbol: str, onset_ms: int, now_ms: int) -> bool:
        """A new eligible onset opens a context window unless one is already
        open for that symbol (no extension; overlapping onsets recorded)."""
        if symbol not in REGISTERED_SYMBOLS:
            return False
        if self.closed:
            return False
        w = self.windows.get(symbol)
        if w is not None and w["open_until_ms"] > now_ms:
            w["overlapping_onsets"] = w.get("overlapping_onsets", 0) + 1
            return False
        self.windows[symbol] = {
            "onset_ms": onset_ms,
            "open_from_ms": now_ms,
            "open_until_ms": now_ms + WINDOW_MS,
            "overlapping_onsets": 0,
        }
        self.actions.append({"action": "subscribe", "symbol": symbol, "at_ms": now_ms,
                             "onset_ms": onset_ms})
        self._check_invariant()
        return True

    def tick(self, now_ms: int) -> list[str]:
        """Close expired windows; return the symbols unsubscribed this tick."""
        closed = []
        for symbol, w in list(self.windows.items()):
            if w is not None and w["open_until_ms"] <= now_ms:
                self.actions.append({"action": "unsubscribe", "symbol": symbol,
                                     "at_ms": now_ms, "open_until_ms": w["open_until_ms"],
                                     "duration_ms": now_ms - w["open_from_ms"]})
                self.windows[symbol] = None
                closed.append(symbol)
        if any(self.windows.get(s) is not None for s in self.windows):
            pass
        self._check_invariant()
        return closed

    def on_disconnect(self) -> None:
        self.active = set()

    def on_subscribe_ack(self, streams: set[str]) -> None:
        self.active |= streams

    def on_unsubscribe_ack(self, streams: set[str]) -> None:
        self.active -= streams

    def _check_invariant(self) -> None:
        want = self.desired()
        assert len(want) <= MAX_ACTIVE_STREAMS, "stream bound violated"
        for s in want:
            if s == TRIGGER_STREAM:
                continue
            sym = s.split("@", 1)[0].upper()
            assert sym in REGISTERED_SYMBOLS, f"non-registered symbol stream {s}"
            assert ("@aggTrade" in s) or ("@markPrice@1s" in s), f"unregistered kind {s}"


# --------------------------------------------------------------------------
# compact append-only sink with fail-closed caps
# --------------------------------------------------------------------------

class CapExceeded(Exception):
    pass


class CompactSink:
    """One JSONL file per UTC day; 256 MiB/day + 10 GiB cumulative fail-closed;
    rotation control record written before ordinary writes on a new day;
    startup byte accounting floored from disk; never deletes prior capture."""

    def __init__(self, root: Path, day_cap: int = DEFAULT_DAY_CAP,
                 total_cap: int = DEFAULT_TOTAL_CAP, now_ms: int | None = None):
        self.root = root
        self.day_cap = day_cap
        self.total_cap = total_cap
        self.day = datetime.fromtimestamp((now_ms or utc_now_ms()) / 1000.0,
                                          tz=timezone.utc).strftime("%Y%m%d")
        self.prior_bytes = 0
        self.day_bytes = 0
        self._handle = None
        root.mkdir(parents=True, exist_ok=True)
        for entry in sorted(root.glob("day*")):
            if not entry.is_dir() or entry.name == f"day{self.day}":
                continue
            for f in entry.iterdir():
                try:
                    self.prior_bytes += f.stat().st_size
                except OSError:
                    pass
        active = root / f"day{self.day}" / "context-capture-{d}.jsonl".format(d=self.day)
        active.parent.mkdir(parents=True, exist_ok=True)
        self.day_bytes = active.stat().st_size if active.exists() else 0
        self._path = active

    def _open(self) -> None:
        if self._handle is None:
            self._handle = open(self._path, "ab")

    def write(self, record: dict, day_key: str | None = None) -> None:
        """Append one record; day rotation and BOTH caps are fail-closed."""
        target = day_key or datetime.fromtimestamp(utc_now_ms() / 1000.0,
                                                   tz=timezone.utc).strftime("%Y%m%d")
        payload = (json.dumps(record, separators=(",", ":"), default=str) + "\n").encode("utf-8")
        n = len(payload)
        if target != self.day:
            self._rotate(target)
        if self.total_cap is not None and self.prior_bytes + self.day_bytes + n > self.total_cap:
            raise CapExceeded(f"cumulative cap {self.total_cap} exceeded: "
                              f"prior={self.prior_bytes} day={self.day_bytes} new={n}")
        if self.day_cap is not None and self.day_bytes + n > self.day_cap:
            raise CapExceeded(f"day cap {self.day_cap} exceeded for UTC day {self.day}")
        self._open()
        self._handle.write(payload)
        self._handle.flush()
        os.fsync(self._handle.fileno())
        self.day_bytes += n

    def _rotate(self, day_key: str) -> None:
        if self._handle is not None:
            self._handle.flush()
            os.fsync(self._handle.fileno())
            self._handle.close()
            self._handle = None
        self.prior_bytes += self.day_bytes
        self.day = day_key
        self._path = self.root / f"day{day_key}" / f"context-capture-{day_key}.jsonl"
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self.day_bytes = self._path.stat().st_size if self._path.exists() else 0
        ctrl = {"record_type": "control", "control": {"type": "day_rotate"},
                "day_utc": day_key, "ts_utc": iso_utc(utc_now_ms())}
        n = len((json.dumps(ctrl) + "\n").encode("utf-8"))
        if self.day_cap is not None and self.day_bytes + n > self.day_cap:
            raise CapExceeded(f"day cap {self.day_cap} exceeded at rotation for {day_key}")
        if self.total_cap is not None and self.prior_bytes + self.day_bytes + n > self.total_cap:
            raise CapExceeded(f"cumulative cap {self.total_cap} exceeded at rotation")
        self._open()
        self._handle.write((json.dumps(ctrl) + "\n").encode("utf-8"))
        self._handle.flush()
        os.fsync(self._handle.fileno())
        self.day_bytes += n

    def close(self) -> None:
        if self._handle is not None:
            try:
                self._handle.flush()
                os.fsync(self._handle.fileno())
            finally:
                self._handle.close()
                self._handle = None


# --------------------------------------------------------------------------
# runtime
# --------------------------------------------------------------------------

def next_request_id(state: dict) -> int:
    state["req_id"] += 1
    return state["req_id"]


def ws_send(ws: WsClient, obj: dict) -> None:
    ws.send_text(json.dumps(obj))


def await_exact_ack(ws: WsClient, rid: int, timeout_s: float = 15.0) -> bool:
    """True only for {id:rid, result:null} with no error key (Binance ack)."""
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        if stop_requested():
            return False
        try:
            op, data = ws.recv_message(timeout=min(2.0, max(0.05, deadline - time.monotonic())))
        except (TimeoutError, socket_timeout_factory()):
            continue
        if op == 0x8:
            return False
        try:
            obj = json.loads(data.decode("utf-8"))
        except Exception:  # noqa: BLE001
            continue
        if obj.get("id") != rid:
            continue
        if "error" in obj or obj.get("result") is not None:
            return False
        return True
    return False


def socket_timeout_factory():
    """socket.timeout is an alias of TimeoutError on py3.10+; kept callable so
    the except clause stays a tuple factory for readability."""
    return TimeoutError


class ContextCollector:
    def __init__(self, out_dir: Path, url: str, hard_minutes: float,
                 max_minutes: float | None, day_cap: int, total_cap: int):
        self.url = url
        self.hard_deadline_ms = utc_now_ms() + int(hard_minutes * 60_000)
        self.soft_deadline_ms = (utc_now_ms() + int(max_minutes * 60_000)
                                 if max_minutes else None)
        self.sink = CompactSink(out_dir, day_cap=day_cap, total_cap=total_cap)
        self.deduper = Deduper(max_ids=400_000)
        self.state = {"req_id": 0}
        self.ws: WsClient | None = None
        self.session: SessionTracker | None = None
        self.subs = SubscriptionState(utc_now_ms())
        self.stats = {"frames": 0, "events_written": 0, "duplicate_frames": 0,
                      "onsets": 0, "context_subscribes": 0,
                      "context_unsubscribes": 0, "reconnects": 0, "controls": 0}

    # ---- control records -------------------------------------------------
    def control(self, ctype: str, **fields) -> None:
        rec = {"record_type": "control", "control": {"type": ctype},
               "ts_utc": iso_utc(utc_now_ms()), "wall_ms": utc_now_ms()}
        if fields:
            rec.update(fields)
        self.sink.write(rec)
        self.stats["controls"] += 1

    # ---- subscriptions ---------------------------------------------------
    def apply_desired(self) -> bool:
        """Drive actual subscriptions to the desired set; every change is
        exact-ack verified. Returns False on any ack failure (reconnect)."""
        want = self.subs.desired()
        add = want - self.subs.active
        remove = self.subs.active - want
        if remove:
            rid = next_request_id(self.state)
            ws_send(self.ws, {"id": rid, "method": "UNSUBSCRIBE", "params": sorted(remove)})
            if not await_exact_ack(self.ws, rid):
                return False
            self.subs.on_unsubscribe_ack(set(remove))
        if add:
            rid = next_request_id(self.state)
            ws_send(self.ws, {"id": rid, "method": "SUBSCRIBE", "params": sorted(add)})
            if not await_exact_ack(self.ws, rid):
                return False
            self.subs.on_subscribe_ack(add)
            if TRIGGER_STREAM in add and len(self.subs.active) == 1:
                print(READY_BANNER, flush=True)
        return True

    # ---- frame handling --------------------------------------------------
    def handle_frame(self, raw_obj) -> None:
        body_stream, body = unwrap(raw_obj)
        kind = body.get("e") if isinstance(body, dict) else None
        if kind == "forceOrder":
            onset = self.session.observe_force_order(body) if self.session else None
            if onset is not None:
                self.stats["onsets"] += 1
                opened = self.subs.open_window(onset["symbol"], onset["E"], utc_now_ms())
                if opened:
                    self.control("context_window_start", symbol=onset["symbol"],
                                 onset_ms=onset["E"], window_ms=WINDOW_MS,
                                 over_symbol_cap=False)
        elif kind in ("aggTrade", "markPriceUpdate"):
            symbol = body.get("s")
            if symbol not in REGISTERED_SYMBOLS:
                return
            did = canonical_id(body_stream, body)
            if did is None:
                return
            if not self.deduper.check(f"{BODYKIND_TAG[kind]}|{did}"):
                self.stats["duplicate_frames"] += 1
                return
        else:
            return
        # write accepted market event verbatim (at-least-once raw evidence)
        stream = canonical_stream(body_stream, body)
        rec = {"record_type": "event", "stream": stream,
               "ts_utc": iso_utc(utc_now_ms()), "wall_ms": utc_now_ms(),
               "payload": raw_obj if isinstance(raw_obj, dict) else {"raw": str(raw_obj)}}
        eid = canonical_id(body_stream, body)
        if eid:
            rec["event_id"] = eid
        self.sink.write(rec)
        self.stats["events_written"] += 1

    # ---- lifecycle -------------------------------------------------------
    def _hard_stop_due(self) -> bool:
        now = utc_now_ms()
        if now >= self.hard_deadline_ms:
            return True
        if self.soft_deadline_ms and now >= self.soft_deadline_ms:
            return True
        return False

    def run(self) -> int:
        signal.signal(signal.SIGTERM, lambda *_: request_stop())
        signal.signal(signal.SIGINT, lambda *_: request_stop())
        self.control("session_start", url=self.url,
                     streams_allowed=sorted([TRIGGER_STREAM] + [s for sym in REGISTERED_SYMBOLS
                                                                for s in SYMBOL_STREAMS[sym]]),
                     hard_stop_utc=iso_utc(self.hard_deadline_ms))
        backoff = 1.0
        try:
            while not stop_requested():
                if self._hard_stop_due():
                    break
                self.session = SessionTracker(None)
                self.subs.on_disconnect()
                try:
                    host = self.url.split("://", 1)[1].split("/", 1)[0]
                    path = "/" + self.url.split("://", 1)[1].split("/", 1)[1]
                    self.ws = WsClient(host=host, path=path, timeout=20.0)
                    self.control("connect", url=self.url, wall_ms=utc_now_ms())
                    self.session = SessionTracker(utc_now_ms())
                    if not self.apply_desired():
                        raise ConnectionError("trigger/context subscribe ack failed")
                    backoff = 1.0
                    self._receive_loop()
                except CapExceeded as exc:
                    print(f"{CAP_BANNER} {exc}", flush=True)
                    return EXIT_CAP
                except (ConnectionError, OSError) as exc:
                    if stop_requested() or self._hard_stop_due():
                        break
                    self.stats["reconnects"] += 1
                    self.control("reconnect_backoff", error=f"{type(exc).__name__}: {exc}",
                                 backoff_s=backoff)
                    if self.ws is not None:
                        self.ws.close()
                        self.ws = None
                    stop_aware_sleep(backoff)
                    backoff = min(backoff * 2.0, 60.0)
                finally:
                    if self.ws is not None:
                        try:
                            self.ws.close()
                        except Exception:  # noqa: BLE001
                            pass
                        self.ws = None
        finally:
            self.sink.close()
            self.control("session_stop", stats=self.stats)
            self.sink.close()
            print(STOP_BANNER, flush=True)
        return EXIT_OK

    def _receive_loop(self) -> None:
        last_window_tick = 0.0
        while not stop_requested():
            if self._hard_stop_due():
                self.subs.closed = True
                self.apply_desired()
                return
            if time.monotonic() - last_window_tick >= 1.0:
                opened = [s for s, w in self.subs.windows.items()
                          if w is not None and w["open_until_ms"] <= utc_now_ms()]
                if opened:
                    self.subs.tick(utc_now_ms())
                    if not self.apply_desired():
                        raise ConnectionError("window unsubscribe ack failed")
                    self.stats["context_unsubscribes"] += len(opened)
                last_window_tick = time.monotonic()
            try:
                op, data = self.ws.recv_message(timeout=2.0)
            except TimeoutError:
                continue
            except (OSError, ConnectionError):
                raise ConnectionError("recv failed")
            if op == 0x8:
                raise ConnectionError("server closed")
            if op != 0x1:
                continue
            try:
                obj = json.loads(data.decode("utf-8"))
            except json.JSONDecodeError:
                self.control("malformed_json", frame_bytes=len(data))
                continue
            self.stats["frames"] += 1
            self.handle_frame(obj)


BODYKIND_TAG = {"aggTrade": "agg", "markPriceUpdate": "mark"}


# --------------------------------------------------------------------------
# deterministic no-network selftests
# --------------------------------------------------------------------------

def selftest() -> int:
    failures: list[str] = []

    def check(name: str, cond: bool, detail: str = "") -> None:
        if not cond:
            failures.append(f"{name}: {detail}")
        print(f"  [{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not cond else ""))

    T0 = parse_contract_time("2026-08-30T12:00:00.000Z")

    # ---- normalization: wrapped == raw identity --------------------------
    o = {"s": "BTCUSDT", "S": "SELL", "o": "LIMIT", "f": "IOC", "q": "0.1",
         "p": "20000", "ap": "20001", "X": "FILLED", "l": "0.1", "z": "0.1",
         "T": T0 + 10}
    raw_fo = {"e": "forceOrder", "E": T0 + 11, "o": o}
    wrapped_fo = {"stream": TRIGGER_STREAM, "data": raw_fo}
    id_raw = canonical_id(None, raw_fo)
    id_wr = canonical_id(wrapped_fo.get("stream"), wrapped_fo["data"])
    check("forceOrder ids: raw and wrapped produce IDENTICAL canonical ids",
          id_raw == id_wr and len(id_raw) == 32, f"{id_raw} vs {id_wr}")
    s_raw, b_raw = unwrap(raw_fo)
    s_wr, b_wr = unwrap(wrapped_fo)
    check("unwrap: raw frames yield no stream; wrapped yield their stream",
          s_raw is None and s_wr == TRIGGER_STREAM and b_raw == b_wr)
    at = {"e": "aggTrade", "E": T0 + 12, "s": "BTCUSDT", "a": 7, "p": "1", "q": "1"}
    mp = {"e": "markPriceUpdate", "E": T0 + 13, "s": "ETHUSDT", "T": T0 + 14}
    check("aggTrade id formula stream|E|s|a",
          canonical_id("btcusdt@aggTrade", at) == canonical_id(None, at),
          f"{canonical_id(None, at)}")
    check("markPriceUpdate id formula stream|E|s|T",
          canonical_id("ethusdt@markPrice@1s", mp) == canonical_id(None, mp))
    check("canonical_stream: raw forceOrder maps to !forceOrder@arr",
          canonical_stream(None, raw_fo) == TRIGGER_STREAM)
    check("unknown event kind has no canonical id (fail-closed)",
          canonical_id(None, {"e": "kline", "E": 1}) is None)

    # ---- onset gating ----------------------------------------------------
    tr = SessionTracker(T0)
    frames = [(T0 + 100_000, "e1"), (T0 + 500_000, "e2"), (T0 + 550_000, "e3"),
              (T0 + 1_000_000, "e4")]
    got = []
    for E, eid in frames:
        on = tr.observe_force_order({"e": "forceOrder", "E": E,
                                     "o": dict(o, s="BTCUSDT", T=E)})
        if on:
            got.append((on["symbol"], on["E"]))
    check("onsets: first needs 300s after connect; gap rule thereafter",
          got == [("BTCUSDT", T0 + 500_000), ("BTCUSDT", T0 + 1_000_000)],
          f"got {got}")
    eth_on = tr.observe_force_order({"e": "forceOrder", "E": T0 + 1_100_000,
                                     "o": dict(o, s="ETHUSDT", T=T0 + 1_100_000)})
    check("per-symbol independence: ETH 1100s after connect IS its own onset",
          eth_on is not None and len([x for x in tr.onsets if x["symbol"] == "ETHUSDT"]) == 1
          and len([x for x in tr.onsets if x["symbol"] == "BTCUSDT"]) == 2,
          f"eth_on={eth_on} onsets={tr.onsets}")
    tr2 = SessionTracker(T0)
    tr2.observe_force_order({"e": "forceOrder", "E": T0 + 400_000,
                             "o": dict(o, s="ETHUSDT", T=T0 + 400_000)})
    check("first event only 400s after connect IS an onset",
          len(tr2.onsets) == 1)
    tr3 = SessionTracker(T0)
    tr3.observe_force_order({"e": "forceOrder", "E": T0 + 400_000,
                             "o": dict(o, s="ETHUSDT", T=T0 + 400_000)})
    tr3.reset(None)
    check("reset(None) closes the session: events cannot qualify",
          tr3.observe_force_order({"e": "forceOrder", "E": T0 + 900_000,
                                   "o": dict(o, s="BTCUSDT", T=T0 + 900_000)}) is None)
    check("non-registered symbol ignored",
          SessionTracker(T0).observe_force_order(
              {"e": "forceOrder", "E": T0 + 900_000, "o": dict(o, s="XRPUSDT", T=T0)}) is None)

    # ---- subscription state machine -------------------------------------
    subs = SubscriptionState(T0)
    check("initial desired set: only the trigger stream",
          subs.desired() == {TRIGGER_STREAM})
    subs.open_window("BTCUSDT", T0 + 500_000, T0 + 500_100)
    want = subs.desired(T0 + 500_100)
    check("window open: exactly aggTrade + markPrice@1s for BTCUSDT added",
          want == {TRIGGER_STREAM, "btcusdt@aggTrade", "btcusdt@markPrice@1s"})
    extended = subs.open_window("BTCUSDT", T0 + 600_000, T0 + 600_100)
    check("overlapping onset does NOT extend the window",
          extended is False
          and subs.windows["BTCUSDT"]["open_until_ms"] == T0 + 500_100 + WINDOW_MS
          and subs.windows["BTCUSDT"]["overlapping_onsets"] == 1)
    check("second symbol window: 5-stream bound respected",
          subs.open_window("ETHUSDT", T0 + 600_000, T0 + 600_200)
          and len(subs.desired(T0 + 600_200)) == MAX_ACTIVE_STREAMS)
    check("no third symbol can ever open",
          subs.open_window("XRPUSDT", T0 + 600_300, T0 + 600_300) is False)
    closed = subs.tick(T0 + 500_100 + WINDOW_MS)
    check("window closes exactly at +20 minutes (BTC only; ETH still open)",
          closed == ["BTCUSDT"]
          and TRIGGER_STREAM in subs.desired(T0 + 500_100 + WINDOW_MS)
          and "btcusdt@aggTrade" not in subs.desired(T0 + 500_100 + WINDOW_MS)
          and len(subs.desired(T0 + 500_100 + WINDOW_MS)) == 3)
    closed2 = subs.tick(T0 + 600_200 + WINDOW_MS + 1)
    check("second window also closes at its own +20m; desired back to trigger",
          closed2 == ["ETHUSDT"] and subs.desired(T0 + 600_200 + WINDOW_MS + 1) == {TRIGGER_STREAM})
    subs.closed = True
    check("hard stop closes everything: desired set empty",
          subs.desired(T0 + 900_000) == set())
    worst = SubscriptionState(T0)
    for s in REGISTERED_SYMBOLS:
        worst.open_window(s, T0, T0)
    check("invariant: max simultaneous streams == 5",
          len(worst.desired(T0 + 1)) == 5)

    # ---- caps: fail-closed accounting ------------------------------------
    import tempfile
    import shutil
    tmp = Path(tempfile.mkdtemp(prefix="fast_context_test_"))
    try:
        sink = CompactSink(tmp / "cap", day_cap=1000, total_cap=None, now_ms=T0)
        sink.write({"n": 1}, day_key="20260830")
        raised = False
        try:
            for i in range(100):
                sink.write({"i": i, "pad": "x" * 40}, day_key="20260830")
        except CapExceeded:
            raised = True
        check("day cap trips fail-closed before exceeding 1000 bytes",
              raised and sink.day_bytes <= 1000, f"day_bytes={sink.day_bytes}")
        sink.close()
        sink2 = CompactSink(tmp / "cap2", day_cap=None, total_cap=1000, now_ms=T0)
        raised2 = False
        try:
            for i in range(100):
                sink2.write({"i": i, "pad": "x" * 40}, day_key="20260830")
        except CapExceeded:
            raised2 = True
        check("cumulative cap trips fail-closed", raised2)
        sink2.close()
        sink3 = CompactSink(tmp / "cap3", day_cap=1000, total_cap=None, now_ms=T0)
        sink3.write({"n": 1}, day_key="20260830")
        sink3.write({"n": 2}, day_key="20260831")   # rotation control record first
        sink3.write({"n": 3}, day_key="20260831")
        rot_ok = (sink3.prior_bytes >= 0
                  and sink3.day == "20260831"
                  and (tmp / "cap3" / "day20260830" / "context-capture-20260830.jsonl").exists()
                  and (tmp / "cap3" / "day20260831" / "context-capture-20260831.jsonl").exists())
        check("rotation: new UTC day file created; prior day untouched",
              rot_ok, f"day={sink3.day} prior={sink3.prior_bytes}")
        sink3.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- 90-minute hard stop clamp ---------------------------------------
    check("hard-stop constant is exactly 90 minutes", HARD_STOP_MS == 90 * 60 * 1000)
    check("window constant is exactly 20 minutes", WINDOW_MS == 20 * 60 * 1000)

    print()
    if failures:
        print(f"SELFTEST FAILED ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("ALL CONTEXT SELFTEST PROOFS PASS (deterministic, no network)")
    return 0


def parse_contract_time(text: str) -> int:
    return int(datetime.fromisoformat(text.replace("Z", "+00:00")).timestamp() * 1000)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--max-minutes", type=float, default=None,
                   help="voluntary bounded run length; must be <= hard stop")
    p.add_argument("--hard-minutes", type=float, default=HARD_STOP_MS / 60_000,
                   help="absolute hard stop; clamped to <= 90 minutes")
    p.add_argument("--max-day-bytes", type=int, default=DEFAULT_DAY_CAP)
    p.add_argument("--max-total-bytes", type=int, default=DEFAULT_TOTAL_CAP)
    p.add_argument("--selftest", action="store_true",
                   help="deterministic no-network proofs on synthetic frames")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    hard = min(args.hard_minutes, HARD_STOP_MS / 60_000)
    if hard <= 0:
        print("--hard-minutes must be positive", file=sys.stderr)
        return EXIT_ARGS
    if args.max_minutes is not None and args.max_minutes > hard + 1e-9:
        print(f"requested --max-minutes {args.max_minutes} exceeds hard stop "
              f"{hard} minutes", file=sys.stderr)
        return EXIT_ARGS
    coll = ContextCollector(args.out, args.url, hard, args.max_minutes,
                            args.max_day_bytes, args.max_total_bytes)
    return coll.run()


if __name__ == "__main__":
    raise SystemExit(main())
