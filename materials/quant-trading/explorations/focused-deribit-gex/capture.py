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
    underlying move), `index_price`, `underlying_price`, `mark_iv`, and
    `open_interest` in coins. Full-chain subscription is required; persistence
    keeps one compact, source-timestamped sample per instrument per 15-minute
    bucket. Ticker bytes have isolated daily and cumulative caps; exhaustion
    pauses ticker persistence without stopping trade, index, metadata, or
    lifecycle capture.
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
    with capped, stop-aware backoff. Full-chain ticker readiness additionally
    requires one schema-valid payload from every subscribed ticker channel.
    DERIBIT_CAPTURE_READY prints only after exact acks, ticker coverage, and a
    successful set_heartbeat ack.
  - reconnect backoff resets ONLY after a stable ready session, not on TCP setup
  - every wait (reconnect sleep, REST retry sleep, ack waits) is SIGTERM-bounded
    (<= 0.1 s wake); critical and ticker storage each have independent
    256 MiB per-UTC-day + 10 GiB cumulative defaults. Critical cap exhaustion
    fails closed; ticker cap exhaustion keeps critical streams running. Prior
    capture is NEVER auto-deleted.
  - oversized restart checkpoints fail closed without replacement: durable cap
    pauses cannot be reconstructed from captured ticker rows.
  - SIGTERM/SIGINT: bounded clean stop; files flushed; exit 0.

Deterministic dedup id examples:
  trade:      sha256("trade|"+instrument_name+"|"+trade_id+"|"+trade_seq+"|"+timestamp)
  ticker:     sha256("ticker|"+instrument_name+"|"+timestamp)
  index:      sha256("index|"+index_name+"|"+timestamp)
  instrument: sha256("instr|"+instrument_name+"|"+creation_timestamp+"|"+state_ts)

No orders, no auth, no private endpoints. This file must remain stdlib-only.
"""

from __future__ import annotations

from collections import deque
import argparse
import base64
import fcntl
import hashlib
import json
import math
import os
import threading
import signal
import socket
import ssl
import stat
import struct
import sys
import time
import uuid
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import NamedTuple
from urllib.request import Request, urlopen

WS_HOST = "www.deribit.com"
WS_PATH = "/ws/api/v2/"
REST_BASE = "https://www.deribit.com/api/v2"
WS_HANDSHAKE_MAX_HEADER_BYTES = 64 * 1024
WS_MAX_MESSAGE_BYTES = 16 * 1024 * 1024
WS_ACCEPT_GUID = "258EAFA5-E914-47DA-95CA-C5AB0DC85B11"

CURRENCIES = ("btc", "eth")  # lowercase for channel names
OPTION_EXPIRY_SETTLEMENT_MS = 8 * 3_600_000  # Deribit settles options 08:00Z
OPTION_CONTRACT_SIZE_COINS = 1  # Deribit options are one-coin contracts
OPTION_OPEN_STATE = "open"  # get_instruments(expired=false) open chain
INSTRUMENT_TERMINAL_STATES = frozenset(
    {"closed", "settled", "expired", "inactive"}
)
OPTION_EXPIRY_MONTHS = {
    "JAN": 1,
    "FEB": 2,
    "MAR": 3,
    "APR": 4,
    "MAY": 5,
    "JUN": 6,
    "JUL": 7,
    "AUG": 8,
    "SEP": 9,
    "OCT": 10,
    "NOV": 11,
    "DEC": 12,
}

class OptionInstrumentIdentity(NamedTuple):
    currency: str
    expiry: datetime
    strike: Decimal
    option_type: str

TICKER_SAMPLE_SECONDS = 15 * 60
DEFAULT_DAY_CAP_BYTES = 256 * 1024 * 1024
DEFAULT_TOTAL_CAP_BYTES = 10 * 1024 * 1024 * 1024
DEFAULT_TICKER_DAY_CAP_BYTES = DEFAULT_DAY_CAP_BYTES
DEFAULT_TICKER_TOTAL_CAP_BYTES = DEFAULT_TOTAL_CAP_BYTES
TICKER_COVERAGE_TIMEOUT_SECONDS = 45.0
TICKER_STALE_SECONDS = 60.0
TICKER_HEALTH_CHECK_SECONDS = 5.0
REQUIRED_SOURCE_MAX_SKEW_MS = int(TICKER_STALE_SECONDS * 1000)
CHAIN_REFRESH_SECONDS = 5 * 60
CHAIN_REFRESH_RETRY_SECONDS = 30.0
CHAIN_MAX_STALE_SECONDS = 3 * CHAIN_REFRESH_SECONDS
TRADE_SEQUENCE_STATE_FILE = "trade-sequence-state.json"
TRADE_SEQUENCE_STATE_SCHEMA_VERSION = 6
OUTPUT_LEASE_FILE = ".collector.lock"
SESSION_IDLE_TIMEOUT_SECONDS = 30.0
CONNECTION_IDLE_LOG_SECONDS = 16.0
INTERRUPT_POLL_SECONDS = 0.1
RESTART_CURSOR_ANCHOR_BYTES = 64
RESTART_CURSOR_MAX_REPLAY_FILES = 2
RESTART_CURSOR_MAX_REPLAY_BYTES = 32 * 1024 * 1024
RESTART_TAIL_BLOCK_BYTES = 64 * 1024
RESTART_TAIL_MAX_BYTES = 32 * 1024 * 1024
RESTART_CHECKPOINT_BYTE_INTERVAL = 4 * 1024 * 1024
UNSUBSCRIBE_GRACE_SECONDS = 30.0
TRADE_SEQUENCE_STATE_MAX_BYTES = 32 * 1024 * 1024
AUDIT_TEXT_MAX_CHARS = 256
RESTART_CHECKPOINT_TIME_INTERVAL = 60.0
SINK_BATCH_CHUNK_BYTES = 1024 * 1024
INSTRUMENT_SNAPSHOT_MAX_ROWS_PER_CURRENCY = 20_000
INSTRUMENT_INITIAL_MIN_ROWS_PER_CURRENCY = 10
INSTRUMENT_EXPIRY_SETTLEMENT_GRACE_MS = 5 * 60 * 1000
INSTRUMENT_MAX_UNEXPECTED_REMOVAL_FRACTION = 0.20
INSTRUMENT_SNAPSHOT_QUARANTINE_ERROR_TYPES = frozenset({
    "CatastrophicSnapshotTruncation",
    "UnexpectedRemovalConfirmationRequired",
    "SchemaMismatch",
})
TICKER_FIELDS = (
    "timestamp",
    "state",
    "instrument_name",
    "index_price",
    "underlying_price",
    "underlying_index",
    "estimated_delivery_price",
    "open_interest",
    "mark_price",
    "mark_iv",
    "greeks",
)

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

def is_finite_number(value) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
    ) or (
        isinstance(value, float)
        and math.isfinite(value)
    )


class DayCapError(RuntimeError):
    """Raised when the per-UTC-day (or cumulative) byte cap is exhausted."""

class TickerCapError(RuntimeError):
    """Raised when only the isolated ticker byte budget is exhausted."""

    def __init__(self, scope: str, message: str):
        super().__init__(message)
        self.scope = scope

class CollectorLeaseError(RuntimeError):
    """Raised when another process owns the output directory."""

class TradeContinuityError(ConnectionError):
    """Raised after the gap boundary and observed trade are durable."""


class TradeRejectionError(ConnectionError):
    """One trade was refused before any durable write of that trade."""


class UnresolvedInstrumentTradeError(ConnectionError):
    """Raised after a trade without current metadata is durable."""

class SourceTimestampIntegrityError(ConnectionError):
    """Raised before an invalid required-channel source time mutates state."""

class DeribitProtocolFrameError(ConnectionError):
    """A JSON-RPC frame has an unsafe structural shape."""
class NonFiniteJsonConstantError(ValueError):
    """An untrusted JSON payload used a non-RFC numeric constant."""


class JsonSerializationIntegrityError(ConnectionError):
    """A record cannot be represented as strict RFC JSON."""


class DuplicateJsonKeyError(ValueError):
    """An untrusted JSON object repeated a member name."""


def _reject_non_finite_json_constant(_constant: str) -> None:
    raise NonFiniteJsonConstantError("non-finite JSON constant")


def _reject_duplicate_json_keys(pairs):
    seen: set[str] = set()
    for key, _value in pairs:
        if key in seen:
            raise DuplicateJsonKeyError(
                f"duplicate JSON member: {key[:64]}"
            )
        seen.add(key)
    return dict(pairs)


def strict_json_loads(payload):
    return json.loads(
        payload,
        parse_constant=_reject_non_finite_json_constant,
        object_pairs_hook=_reject_duplicate_json_keys,
    )


def strict_json_dumps(value, **kwargs) -> str:
    kwargs["allow_nan"] = False
    try:
        return json.dumps(value, **kwargs)
    except ValueError as exc:
        raise JsonSerializationIntegrityError(
            "record contains a non-RFC JSON value"
        ) from exc


class NotARegularCaptureFileError(OSError):
    """A capture, state, or temp path is not a regular file."""


def _nofollow_flags(base: int) -> int:
    flags = base
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    flags |= getattr(os, "O_NONBLOCK", 0)
    return flags


def _regular_fdopen(fd: int, path: Path, mode: str, **kwargs):
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            raise NotARegularCaptureFileError(
                f"capture path is not a regular file: {path}"
            )
    except BaseException:
        os.close(fd)
        raise
    return os.fdopen(fd, mode, **kwargs)


def open_regular_append(path: Path):
    """Append to a regular file, never a symlink, FIFO, or device."""
    return _regular_fdopen(
        os.open(
            path,
            _nofollow_flags(os.O_RDWR | os.O_CREAT | os.O_APPEND),
            0o600,
        ),
        path,
        "ab",
    )


def open_regular_read(path: Path):
    """Read a regular capture file, never a symlink, FIFO, or device."""
    return _regular_fdopen(
        os.open(path, _nofollow_flags(os.O_RDONLY)),
        path,
        "rb",
    )


def create_regular_temp(path: Path):
    """Create an exclusive regular temp file for atomic replacement."""
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    return _regular_fdopen(
        os.open(
            path,
            _nofollow_flags(os.O_WRONLY | os.O_CREAT | os.O_EXCL),
            0o600,
        ),
        path,
        "w",
        encoding="utf-8",
    )


def require_capture_directory(path: Path) -> None:
    """Fail closed when a capture directory entry is a link or special."""
    try:
        mode = os.lstat(path).st_mode
    except FileNotFoundError:
        return
    if not stat.S_ISDIR(mode):
        raise NotARegularCaptureFileError(
            f"capture directory is not a directory: {path}"
        )



def bounded_error_provenance(error: object) -> dict:
    """Bounded, typed provenance for one untrusted server error payload."""
    if isinstance(error, dict):
        code = error.get("code")
        message = error.get("message", error.get("msg"))
        return {
            "code": (
                code
                if isinstance(code, int) and not isinstance(code, bool)
                else None
            ),
            "msg": (
                None
                if message is None
                else str(message)[:AUDIT_TEXT_MAX_CHARS]
            ),
        }
    return {"code": None, "msg": str(error)[:AUDIT_TEXT_MAX_CHARS]}


def bounded_audit_text(text: object) -> str:
    """Clamp one untrusted string before it reaches a durable audit row."""
    return str(text)[:AUDIT_TEXT_MAX_CHARS]


def bounded_audit_detail(detail, depth: int = 0):
    """Clamp every string inside one capture-event detail payload."""
    if depth > 4:
        return None
    if isinstance(detail, str):
        return detail[:AUDIT_TEXT_MAX_CHARS]
    if isinstance(detail, dict):
        return {
            str(key)[:AUDIT_TEXT_MAX_CHARS]: bounded_audit_detail(
                value,
                depth + 1,
            )
            for key, value in detail.items()
        }
    if isinstance(detail, (list, tuple)):
        return [
            bounded_audit_detail(value, depth + 1)
            for value in detail
        ]
    return detail


def fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_DIRECTORY", 0)
    flags |= getattr(os, "O_CLOEXEC", 0)
    directory_fd = os.open(path, flags)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
def create_durable_directory(path: Path) -> None:
    """Create each missing path component and commit its parent entry."""
    missing = []
    current = path
    while True:
        try:
            metadata = current.stat()
        except FileNotFoundError:
            if current.parent == current:
                raise
            missing.append(current)
            current = current.parent
            continue
        if not stat.S_ISDIR(metadata.st_mode):
            raise NotADirectoryError(str(current))
        break
    for directory in reversed(missing):
        try:
            directory.mkdir()
        except FileExistsError:
            if not directory.is_dir():
                raise
        fsync_directory(directory.parent)





class InstrumentSnapshotError(ConnectionError):
    """Required-currency REST snapshot fetch or validation failure."""

    def __init__(
        self,
        currency: str,
        error_type: str,
        message: str,
        detail: dict | None = None,
    ):
        super().__init__(message)
        self.currency = currency
        self.error_type = error_type
        self.detail = dict(detail or {})




class WsClient:
    """Minimal RFC6455 client (stdlib-only) for one WS connection."""

    def __init__(self, host: str = WS_HOST, path: str = WS_PATH, timeout: float = 20.0):
        self.host = host
        self.path = path
        self.timeout = timeout
        self._connect()

    def _connect(self):
        deadline = time.monotonic() + max(0.01, self.timeout)
        raw = socket.create_connection(
            (self.host, 443),
            timeout=max(0.01, deadline - time.monotonic()),
        )
        raw.settimeout(max(0.01, deadline - time.monotonic()))
        try:
            sock = ssl.create_default_context().wrap_socket(
                raw,
                server_hostname=self.host,
            )
        except BaseException:
            raw.close()
            raise
        try:
            key = base64.b64encode(os.urandom(16)).decode("ascii")
            request = (
                f"GET {self.path} HTTP/1.1\r\n"
                f"Host: {self.host}\r\n"
                "Upgrade: websocket\r\n"
                "Connection: Upgrade\r\n"
                f"Sec-WebSocket-Key: {key}\r\n"
                "Sec-WebSocket-Version: 13\r\n"
                "\r\n"
            )
            sock.sendall(request.encode("ascii"))
            response = bytearray()
            delimiter = -1
            while delimiter < 0:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise ConnectionError(
                        "websocket handshake deadline reached"
                    )
                sock.settimeout(remaining)
                chunk = sock.recv(4096)
                if not chunk:
                    raise ConnectionError("websocket handshake closed")
                response.extend(chunk)
                delimiter = response.find(b"\r\n\r\n")
                if (
                    delimiter < 0
                    and len(response) > WS_HANDSHAKE_MAX_HEADER_BYTES
                ):
                    raise ConnectionError(
                        "websocket handshake headers exceed limit"
                    )
            if delimiter > WS_HANDSHAKE_MAX_HEADER_BYTES:
                raise ConnectionError(
                    "websocket handshake headers exceed limit"
                )
            head = bytes(response[:delimiter])
            lines = head.split(b"\r\n")
            status = lines[0].split(b" ", 2)
            if (
                len(status) < 2
                or status[0] != b"HTTP/1.1"
                or status[1] != b"101"
            ):
                raise ConnectionError(
                    f"websocket handshake status invalid: "
                    f"{lines[0][:120]!r}"
                )
            headers: dict[bytes, list[bytes]] = {}
            for line in lines[1:]:
                if b":" not in line:
                    raise ConnectionError(
                        "websocket handshake header is malformed"
                    )
                name, value = line.split(b":", 1)
                headers.setdefault(name.strip().lower(), []).append(
                    value.strip()
                )
            upgrade_tokens = {
                token.strip().lower()
                for value in headers.get(b"upgrade", [])
                for token in value.split(b",")
            }
            connection_tokens = {
                token.strip().lower()
                for value in headers.get(b"connection", [])
                for token in value.split(b",")
            }
            if b"websocket" not in upgrade_tokens:
                raise ConnectionError(
                    "websocket handshake lacks Upgrade: websocket"
                )
            if b"upgrade" not in connection_tokens:
                raise ConnectionError(
                    "websocket handshake lacks Connection: upgrade"
                )
            accept_values = headers.get(b"sec-websocket-accept", [])
            expected_accept = base64.b64encode(
                hashlib.sha1(
                    (key + WS_ACCEPT_GUID).encode("ascii")
                ).digest()
            )
            if (
                len(accept_values) != 1
                or accept_values[0] != expected_accept
            ):
                raise ConnectionError(
                    "websocket handshake accept mismatch"
                )
        except BaseException:
            try:
                sock.close()
            except Exception:
                pass
            raise
        self.sock = sock
        self._buf = bytes(response[delimiter + 4:])
        self._frag = None
        self._frag_size = 0

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


    def _recv_at_least(
        self,
        data: bytes,
        required: int,
        deadline: float | None = None,
    ) -> bytes:
        if len(data) >= required:
            return data
        buffered = bytearray(data)
        while len(buffered) < required:
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    self._buf = bytes(buffered)
                    raise socket.timeout(
                        "websocket receive deadline"
                    )
                self.sock.settimeout(remaining)
            try:
                chunk = self.sock.recv(
                    min(65536, required - len(buffered))
                )
            except OSError:
                self._buf = bytes(buffered)
                raise
            if not chunk:
                self._buf = bytes(buffered)
                raise ConnectionError("closed mid-frame")
            buffered.extend(chunk)
        return bytes(buffered)

    def recv_message(self, timeout: float) -> tuple[int, bytes]:
        """Read one complete message without discarding partial socket data."""
        data = self._buf
        frag = self._frag
        frag_size = getattr(
            self,
            "_frag_size",
            sum(len(part) for part in frag) if frag is not None else 0,
        )
        deadline = time.monotonic() + max(0.0, timeout)

        def check_deadline() -> None:
            if time.monotonic() >= deadline:
                raise socket.timeout("websocket receive deadline")

        while True:
            data = self._recv_at_least(data, 2, deadline)
            op = data[0] & 0x0F
            fin = bool(data[0] & 0x80)
            masked = bool(data[1] & 0x80)
            length = data[1] & 0x7F
            pos = 2
            if length == 126:
                data = self._recv_at_least(data, 4, deadline)
                length = struct.unpack(">H", data[2:4])[0]
                pos = 4
            elif length == 127:
                data = self._recv_at_least(data, 10, deadline)
                length = struct.unpack(">Q", data[2:10])[0]
                pos = 10
            is_control = bool(op & 0x08)
            if is_control:
                if not fin or length > 125:
                    self._frag = None
                    self._frag_size = 0
                    raise ConnectionError(
                        "invalid websocket control frame"
                    )
            else:
                message_size = (
                    frag_size + length
                    if op == OP_CONT and frag is not None
                    else length
                )
                if message_size > WS_MAX_MESSAGE_BYTES:
                    self._frag = None
                    self._frag_size = 0
                    raise ConnectionError(
                        "websocket message exceeds "
                        f"{WS_MAX_MESSAGE_BYTES} bytes"
                    )
            need = pos + (4 if masked else 0) + length
            data = self._recv_at_least(data, need, deadline)
            if masked:
                mask = data[pos:pos + 4]
                body = bytes(
                    byte ^ mask[index % 4]
                    for index, byte in enumerate(
                        data[pos + 4:pos + 4 + length]
                    )
                )
            else:
                body = data[pos:pos + length]
            data = data[need:]
            self._buf = data
            if op == OP_PING:
                self._send_frame(OP_PONG, body)
                check_deadline()
                continue
            if op == OP_PONG:
                check_deadline()
                continue
            if op == OP_CLOSE:
                self._frag = None
                self._frag_size = 0
                return OP_CLOSE, body
            if op == OP_TEXT and frag is None:
                if fin:
                    self._frag = None
                    self._frag_size = 0
                    return OP_TEXT, body
                frag = [body]
                frag_size = len(body)
                self._frag = frag
                self._frag_size = frag_size
                check_deadline()
                continue
            if op == OP_CONT and frag is not None:
                frag.append(body)
                frag_size += len(body)
                self._frag_size = frag_size
                if fin:
                    merged = b"".join(frag)
                    self._frag = None
                    self._frag_size = 0
                    return OP_TEXT, merged
                check_deadline()
                continue
            frag = None
            frag_size = 0
            self._frag = None
            self._frag_size = 0
            check_deadline()

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


REST_ATTEMPT_TIMEOUT_SECONDS = 20.0
REST_MAX_RESPONSE_BYTES = 64 * 1024 * 1024
REST_READ_BLOCK_BYTES = 256 * 1024


def bound_response_read_deadline(response: object, remaining: float) -> None:
    """Best-effort per-read socket timeout honoring the total deadline."""
    sock = getattr(getattr(response, "fp", None), "raw", None)
    sock = getattr(sock, "_sock", None)
    if sock is None:
        return
    try:
        sock.settimeout(max(0.001, remaining))
    except OSError:
        pass


def read_bounded_response(response, deadline: float, path: str) -> bytes:
    """Read one response under a hard total deadline and byte ceiling."""
    body = bytearray()
    while True:
        if stop_requested():
            raise ConnectionError(f"rest_call aborted: stop requested {path}")
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise TimeoutError(f"rest_call read deadline reached for {path}")
        bound_response_read_deadline(response, remaining)
        chunk = response.read(REST_READ_BLOCK_BYTES)
        if not chunk:
            return bytes(body)
        body += chunk
        if len(body) > REST_MAX_RESPONSE_BYTES:
            raise ValueError(
                f"rest_call response exceeds {REST_MAX_RESPONSE_BYTES} "
                f"bytes for {path}"
            )



def rest_call(
    path: str,
    params: dict | None = None,
    attempts: int = 3,
    deadline: float | None = None,
) -> dict:
    query = ""
    if params:
        query = "?" + "&".join(f"{key}={value}" for key, value in params.items())
    url = f"{REST_BASE}{path}{query}"
    last: Exception | None = None
    for attempt in range(attempts):
        if stop_requested():
            raise ConnectionError("rest_call aborted: stop requested")
        remaining = (
            deadline - time.monotonic()
            if deadline is not None
            else float("inf")
        )
        if remaining <= 0:
            raise TimeoutError(f"rest_call deadline reached for {path}")
        attempt_deadline = time.monotonic() + min(
            REST_ATTEMPT_TIMEOUT_SECONDS,
            max(0.01, remaining),
        )
        try:
            req = Request(
                url,
                headers={"User-Agent": "quant-research-gex-capture/1.0"},
            )
            with urlopen(
                req,
                timeout=max(0.01, attempt_deadline - time.monotonic()),
            ) as resp:
                return strict_json_loads(
                    read_bounded_response(
                        resp,
                        attempt_deadline,
                        path,
                    ).decode()
                )
        except Exception as exc:  # collector retries within the same deadline
            last = exc
            remaining = (
                deadline - time.monotonic()
                if deadline is not None
                else float("inf")
            )
            if remaining <= 0 or attempt == attempts - 1:
                break
            stop_aware_sleep(
                min(1.5 * (attempt + 1), max(0.0, remaining))
            )
    raise ConnectionError(f"rest_call failed {path}: {last}")


def fetch_instruments(
    currency: str,
    deadline: float | None = None,
) -> list[dict]:
    payload = rest_call(
        "/public/get_instruments",
        {
            "currency": currency.upper(),
            "kind": "option",
            "expired": "false",
        },
        deadline=deadline,
    )
    if "result" not in payload:
        raise ConnectionError(f"get_instruments error: {payload.get('error')}")
    return payload["result"]


class OutputLease:
    """Kernel-owned exclusive lease; the stable file is never unlinked."""

    def __init__(self, root: Path):
        create_durable_directory(root)
        lock_path = root / OUTPUT_LEASE_FILE
        flags = os.O_RDWR | os.O_CREAT
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        fd = os.open(lock_path, flags, 0o600)
        try:
            if not stat.S_ISREG(os.fstat(fd).st_mode):
                raise CollectorLeaseError(
                    f"output lease is not a regular file: {lock_path}"
                )
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise CollectorLeaseError(
                    f"output directory already has a writer: {root}"
                ) from exc
            self._handle = os.fdopen(fd, "a+")
        except BaseException:
            os.close(fd)
            raise

    def close(self) -> None:
        handle = getattr(self, "_handle", None)
        if handle is None:
            return
        self._handle = None
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()
    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass


class Sink:
    """Append-only JSONL sink with isolated critical/ticker byte accounting,
    daily rotation, and fail-closed writes. Prior capture is never deleted."""

    def __init__(
        self,
        root: Path,
        day_cap_bytes: int,
        total_cap_bytes: int,
        run_id: str | None = None,
        ticker_day_cap_bytes: int | None = None,
        ticker_total_cap_bytes: int | None = None,
    ):
        self.root = root
        self.day_cap = day_cap_bytes
        self.total_cap = total_cap_bytes
        self.ticker_day_cap = (
            day_cap_bytes
            if ticker_day_cap_bytes is None
            else ticker_day_cap_bytes
        )
        self.ticker_total_cap = (
            total_cap_bytes
            if ticker_total_cap_bytes is None
            else ticker_total_cap_bytes
        )
        self.run_id = run_id or uuid.uuid4().hex
        self.day_key: str | None = None
        self.day_bytes = 0
        self.prior_bytes = 0
        self.ticker_day_bytes = 0
        self.ticker_prior_bytes = 0
        self._handles: dict[str, object] = {}
        self._io_error: str | None = None
        self._ticker_tail_cap_error: TickerCapError | None = None
        self._lease = OutputLease(root)
        try:
            self.reconcile_existing()
        except BaseException:
            self._lease.close()
            raise

    def _paths(self, day_key: str) -> dict[str, Path]:
        base = self.root / f"day{day_key}"
        return {
            "trades": base / "trades.jsonl",
            "ticker": base / "ticker.jsonl",
            "index": base / "index.jsonl",
            "instruments": base / "instruments.jsonl",
            "events": base / "events.jsonl",
        }

    def reconcile_existing(self, active_day: str | None = None) -> None:
        """Recount critical and ticker bytes into their isolated budgets."""
        self.prior_bytes = 0
        self.day_bytes = 0
        self.ticker_prior_bytes = 0
        self.ticker_day_bytes = 0
        try:
            root_mode = self.root.stat().st_mode
        except FileNotFoundError:
            return
        if not stat.S_ISDIR(root_mode):
            raise NotADirectoryError(str(self.root))
        for entry in sorted(self.root.iterdir()):
            if not entry.name.startswith("day"):
                continue
            entry_mode = os.lstat(entry).st_mode
            if stat.S_ISLNK(entry_mode):
                raise NotARegularCaptureFileError(
                    f"capture day directory is a symlink: {entry}"
                )
            if not stat.S_ISDIR(entry_mode):
                continue
            critical_size = 0
            ticker_size = 0
            for path in entry.iterdir():
                entry_stat = os.lstat(path)
                if not stat.S_ISREG(entry_stat.st_mode):
                    raise NotARegularCaptureFileError(
                        f"capture path is not a regular file: {path}"
                    )
                if path.name == "ticker.jsonl":
                    ticker_size += entry_stat.st_size
                else:
                    critical_size += entry_stat.st_size
            if (
                active_day is not None
                and entry.name == f"day{active_day}"
            ):
                self.day_bytes += critical_size
                self.ticker_day_bytes += ticker_size
            else:
                self.prior_bytes += critical_size
                self.ticker_prior_bytes += ticker_size
    def _latch_io_error(self, operation: str, exc: OSError) -> None:
        if self._io_error is None:
            self._io_error = f"{operation}:{type(exc).__name__}"

    def _check_ticker_cap(self, new_bytes: int) -> None:
        if (
            self.ticker_total_cap is not None
            and self.ticker_prior_bytes
            + self.ticker_day_bytes
            + new_bytes
            > self.ticker_total_cap
        ):
            raise TickerCapError(
                "cumulative",
                f"ticker cumulative cap {self.ticker_total_cap} bytes "
                f"exceeded: disk={self.ticker_prior_bytes} "
                f"active={self.ticker_day_bytes} new={new_bytes}",
            )
        if self.ticker_day_bytes + new_bytes > self.ticker_day_cap:
            raise TickerCapError(
                "daily",
                f"ticker daily cap {self.ticker_day_cap} bytes exceeded "
                f"for UTC day {self.day_key}",
            )


    def write(self, kind: str, record: dict, now_ms: int) -> int:
        return self.write_batch(kind, [record], now_ms)

    def write_batch(
        self,
        kind: str,
        records: list[dict],
        now_ms: int,
    ) -> int:
        if self._io_error is not None:
            raise DayCapError(f"sink I/O terminal: {self._io_error}")
        if kind not in ("trades", "ticker", "index", "instruments", "events"):
            raise ValueError(kind)
        if not records:
            return 0
        payloads = [
            (
                strict_json_dumps(
                    {**record, "run_id": self.run_id},
                    separators=(",", ":"),
                )
                + "\n"
            ).encode("utf-8")
            for record in records
        ]
        total_bytes = sum(len(payload) for payload in payloads)
        target_day = utc_day(now_ms)
        if self.day_key is None or self.day_key != target_day:
            try:
                self._rotate(target_day, now_ms)
            except OSError as exc:
                self._latch_io_error("rotate", exc)
                raise
        if kind == "ticker":
            if self._ticker_tail_cap_error is not None:
                raise self._ticker_tail_cap_error
            self._check_ticker_cap(total_bytes)
        else:
            if (
                self.total_cap is not None
                and self.prior_bytes + self.day_bytes + total_bytes
                > self.total_cap
            ):
                raise DayCapError(
                    f"cumulative cap {self.total_cap} bytes exceeded: "
                    f"disk={self.prior_bytes} active={self.day_bytes} "
                    f"new={total_bytes}"
                )
            if self.day_bytes + total_bytes > self.day_cap:
                raise DayCapError(
                    f"daily cap {self.day_cap} bytes exceeded for UTC day "
                    f"{self.day_key}"
                )
        handle = self._handles[kind]
        try:
            chunk = bytearray()
            for payload in payloads:
                if (
                    chunk
                    and len(chunk) + len(payload) > SINK_BATCH_CHUNK_BYTES
                ):
                    written = handle.write(chunk)
                    if written != len(chunk):
                        raise OSError(
                            f"short batch append: {written} of "
                            f"{len(chunk)} bytes"
                        )
                    chunk.clear()
                if len(payload) > SINK_BATCH_CHUNK_BYTES:
                    written = handle.write(payload)
                    if written != len(payload):
                        raise OSError(
                            f"short batch append: {written} of "
                            f"{len(payload)} bytes"
                        )
                else:
                    chunk.extend(payload)
            if chunk:
                written = handle.write(chunk)
                if written != len(chunk):
                    raise OSError(
                        f"short batch append: {written} of "
                        f"{len(chunk)} bytes"
                    )
            handle.flush()
            os.fsync(handle.fileno())
        except OSError as exc:
            self._latch_io_error("append", exc)
            raise
        if kind == "ticker":
            self.ticker_day_bytes += total_bytes
        else:
            self.day_bytes += total_bytes
        return total_bytes
    def _isolate_partial_tails(
        self,
        paths: dict[str, Path],
    ) -> list[str]:
        partial_kinds = []
        for kind, path in paths.items():
            try:
                size = os.lstat(path).st_size
            except FileNotFoundError:
                continue
            except OSError as exc:
                self._latch_io_error("tail-inspect", exc)
                raise
            if size == 0:
                continue
            try:
                with open_regular_read(path) as handle:
                    handle.seek(-1, os.SEEK_END)
                    final_byte = handle.read(1)
            except OSError as exc:
                self._latch_io_error("tail-inspect", exc)
                raise
            if final_byte != b"\n":
                partial_kinds.append(kind)

        critical_recovery_bytes = sum(
            kind != "ticker" for kind in partial_kinds
        )
        if self.day_bytes + critical_recovery_bytes > self.day_cap:
            raise DayCapError(
                f"daily cap {self.day_cap} bytes exceeded while "
                f"isolating partial tails for UTC day {self.day_key}"
            )
        if self.total_cap is not None and (
            self.prior_bytes
            + self.day_bytes
            + critical_recovery_bytes
            > self.total_cap
        ):
            raise DayCapError(
                f"cumulative cap {self.total_cap} bytes exceeded while "
                "isolating partial tails"
            )

        self._ticker_tail_cap_error = None
        recoverable_kinds = list(partial_kinds)
        if "ticker" in partial_kinds:
            try:
                self._check_ticker_cap(1)
            except TickerCapError as exc:
                self._ticker_tail_cap_error = exc
                recoverable_kinds.remove("ticker")

        for kind in recoverable_kinds:
            try:
                with open_regular_append(paths[kind]) as handle:
                    written = handle.write(b"\n")
                    if written != 1:
                        raise OSError(
                            f"short tail separator append: {written} of 1 byte"
                        )
                    handle.flush()
                    os.fsync(handle.fileno())
            except OSError as exc:
                self._latch_io_error("tail-recovery", exc)
                raise
            if kind == "ticker":
                self.ticker_day_bytes += 1
            else:
                self.day_bytes += 1
        return sorted(recoverable_kinds)


    def _rotate(self, day_key: str, received_ms: int) -> None:
        close_error = None
        for handle in self._handles.values():
            try:
                handle.flush()
                os.fsync(handle.fileno())
            except OSError as exc:
                close_error = close_error or exc
            try:
                handle.close()
            except OSError as exc:
                close_error = close_error or exc
        if close_error is not None:
            self._latch_io_error("rotate-close", close_error)
            raise close_error
        self._handles = {}
        self.day_key = day_key
        base = self.root / f"day{day_key}"
        require_capture_directory(base)
        day_created = not base.exists()
        base.mkdir(parents=True, exist_ok=True)
        paths = self._paths(day_key)
        child_entries_created = any(
            not path.exists()
            for path in paths.values()
        )
        self.reconcile_existing(active_day=day_key)
        recovered_partial_tail_kinds = self._isolate_partial_tails(paths)
        for kind, path in paths.items():
            if kind == "ticker" and self._ticker_tail_cap_error is not None:
                continue
            self._handles[kind] = open_regular_append(path)
        rec = {
            "record_type": "capture_event",
            "dedup_id": sha256_id(
                f"event|day_rotation_open|{day_key}|{received_ms}|"
                f"{self.run_id}"
            ),
            "event": "day_rotation_open",
            "day_utc": day_key,
            "received_ts_ms": received_ms,
            "recovered_partial_tail_kinds": recovered_partial_tail_kinds,
            "run_id": self.run_id,
        }
        if self._ticker_tail_cap_error is not None:
            rec["ticker_tail_recovery_blocked"] = True
        ev_payload = (
            strict_json_dumps(rec, separators=(",", ":")) + "\n"
        ).encode("utf-8")
        # The critical rotation event itself must respect both critical caps.
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
        written = h.write(ev_payload)
        if written != len(ev_payload):
            raise OSError(
                f"short rotation append: {written} of {len(ev_payload)} bytes"
            )
        h.flush()
        os.fsync(h.fileno())
        if child_entries_created:
            fsync_directory(base)
        if day_created:
            fsync_directory(self.root)
        self.day_bytes += len(ev_payload)

    def close(self) -> None:
        close_error = None
        for handle in self._handles.values():
            try:
                handle.flush()
                os.fsync(handle.fileno())
            except OSError as exc:
                close_error = close_error or exc
            try:
                handle.close()
            except OSError as exc:
                close_error = close_error or exc
        self._handles = {}
        try:
            self._lease.close()
        except OSError as exc:
            close_error = close_error or exc
        if close_error is not None:
            self._latch_io_error("close", close_error)
            raise close_error


class Deduper:
    """Bounded recent-id FIFO for obvious reconnect duplicates.

    Cross-session dedup remains authoritative in replay.
    """

    def __init__(self, max_ids: int = 500_000):
        if max_ids <= 0:
            raise ValueError("max_ids must be positive")
        self.seen: set[str] = set()
        self.order: deque[str] = deque()
        self.max_ids = max_ids
        self.duplicates = 0

    def is_duplicate(self, dedup_id: str) -> bool:
        if dedup_id in self.seen:
            self.duplicates += 1
            return True
        return False

    def commit(self, dedup_id: str) -> None:
        if dedup_id in self.seen:
            return
        self.seen.add(dedup_id)
        self.order.append(dedup_id)
        if len(self.order) > self.max_ids:
            old = self.order.popleft()
            self.seen.discard(old)


class Collector:
    def __init__(
        self,
        raw_dir: Path,
        currencies=("btc", "eth"),
        day_cap=DEFAULT_DAY_CAP_BYTES,
        total_cap=DEFAULT_TOTAL_CAP_BYTES,
        heartbeat=True,
        max_minutes=None,
        capture_tickers=True,
        ticker_sample_seconds=TICKER_SAMPLE_SECONDS,
        ticker_day_cap=DEFAULT_TICKER_DAY_CAP_BYTES,
        ticker_total_cap=DEFAULT_TICKER_TOTAL_CAP_BYTES,
    ):
        self.raw_dir = raw_dir
        self.run_id = uuid.uuid4().hex
        self._connection_seq = 0
        self._current_connection_id: str | None = None
        normalized_currencies = []
        for currency in currencies:
            if not isinstance(currency, str):
                raise ValueError("currencies must contain strings")
            normalized = currency.strip().lower()
            if normalized not in CURRENCIES:
                raise ValueError(
                    f"unsupported currency {currency!r} "
                    f"(research scope: {', '.join(CURRENCIES)})"
                )
            if normalized not in normalized_currencies:
                normalized_currencies.append(normalized)
        if not normalized_currencies:
            raise ValueError("at least one currency is required")
        self.currencies = tuple(normalized_currencies)
        self._instrument_prefixes = tuple(
            f"{currency.upper()}-" for currency in self.currencies
        )
        self.day_cap = day_cap
        self.total_cap = total_cap
        self.ticker_day_cap = ticker_day_cap
        self.ticker_total_cap = ticker_total_cap
        self.instruments: dict[str, dict] = {}
        self._persisted_instrument_fingerprints: dict[str, str] = {}
        self.heartbeat = heartbeat
        self.max_minutes = max_minutes
        self.capture_tickers = capture_tickers
        self.ticker_sample_ms = int(ticker_sample_seconds * 1000)
        if self.ticker_sample_ms <= 0:
            raise ValueError(
                "ticker_sample_seconds must be at least one millisecond"
            )
        self.stop = False
        self.sink = Sink(
            raw_dir,
            day_cap,
            total_cap,
            run_id=self.run_id,
            ticker_day_cap_bytes=ticker_day_cap,
            ticker_total_cap_bytes=ticker_total_cap,
        )
        self.deduper = Deduper()
        self.active_channels: set[str] = set()
        self._pending_channels: set[str] = set()
        self.ws: WsClient | None = None
        self.req_id = 0
        self.msg_seq = 0
        self._capture_event_seq = 0
        self._ready_printed = False
        self._requested_stop_reason: str | None = None
        self._session_stop_emitted = False
        self.last_message_ms = 0
        self.last_message_monotonic = 0.0
        self._deadline_monotonic: float | None = None
        self.last_utc_day = None
        self.last_chain_reconcile = 0.0
        self._next_chain_refresh_attempt = 0.0
        self._last_authoritative_chain_refresh_monotonic = 0.0
        self._instrument_chain_stale_latched = False
        self._chain_refresh_lock = threading.Lock()
        self._chain_refresh_thread: threading.Thread | None = None
        self._chain_refresh_result: tuple[
            int,
            dict[str, list[dict]] | None,
            Exception | None,
        ] | None = None
        self._chain_refresh_generation = 0
        self._pending_unexpected_instrument_removals: dict[
            str,
            frozenset[str],
        ] = {}
        self._recently_unsubscribed: dict[str, tuple[float, bool]] = {}
        self._stable_session = False
        self._ticker_last_bucket: dict[str, int] = {}
        self._ticker_cap_day: str | None = None
        self._ticker_cap_last_dropped_bucket: dict[str, int] = {}
        self._ticker_total_cap_exhausted = False
        self._required_last_received_monotonic: dict[str, float] = {}
        self._required_last_source_timestamp_ms: dict[
            str,
            int | float,
        ] = {}
        self._checkpoint_source_timestamp_rejected = 0
        self._next_required_health_check = 0.0
        self._last_trade_seq: dict[str, int] = {}
        self._last_trade_source_ts: dict[str, int] = {}
        self._sequence_state_dirty = False
        self._sequence_state_bytes_since_checkpoint = 0
        self._sequence_state_checkpoint_due_monotonic = (
            time.monotonic() + RESTART_CHECKPOINT_TIME_INTERVAL
        )
        self._sequence_state_error: str | None = None
        self._sequence_state_cursors: dict[str, tuple[Path, int]] = {}
        self._sequence_state_rebuild_reason: str | None = None
        self._rehydration_in_progress = False
        self._rehydrated: dict | None = None
        self.stats = {
            "trades_pushes": 0, "trades_records": 0, "dup_trades": 0,
            "trades_unresolved": 0, "trade_sequence_gaps": 0,
            "trade_schema_errors": 0,
            "ticker_pushes": 0, "ticker_records": 0,
            "ticker_samples_suppressed": 0, "ticker_schema_errors": 0,
            "ticker_cap_exhaustions": 0,
            "ticker_records_dropped_cap": 0,
            "source_timestamp_errors": 0,
            "protocol_frame_errors": 0,
            "source_timestamp_repeats": 0,
            "source_timestamp_regressions": 0,
            "index_records": 0, "index_schema_errors": 0,
            "instrument_records": 0, "reconnects": 0,
            "instrument_chain_stale": 0,
            "instrument_snapshot_quarantined": 0,
            "trade_sequence_conflicts": 0,
            "instrument_rows_not_open": 0,
            "unsubscribed_grace_frames": 0,
            "trade_chronology_violations": 0,
            "events": 0, "primary_subscribe_ok": False,
        }

    # ---------- restart-state rehydration ----------
    def _sequence_state_failure(
        self,
        operation: str,
        exc: BaseException,
    ) -> OSError:
        error = f"{operation}: {type(exc).__name__}: {exc}"
        if self._sequence_state_error is None:
            self._sequence_state_error = error
        return OSError(error)

    def _merge_trade_sequence(
        self,
        instrument_name: str,
        sequence: int,
        source_ts_ms: int | None = None,
    ) -> None:
        previous = self._last_trade_seq.get(instrument_name)
        if previous is None or sequence > previous:
            self._last_trade_seq[instrument_name] = sequence
            if source_ts_ms is not None:
                self._last_trade_source_ts[instrument_name] = source_ts_ms
            self._sequence_state_dirty = True

    def _record_restart_source_write(self, written_bytes: int) -> None:
        self._sequence_state_bytes_since_checkpoint += written_bytes
        if (
            self._sequence_state_bytes_since_checkpoint
            >= RESTART_CHECKPOINT_BYTE_INTERVAL
            or time.monotonic()
            >= self._sequence_state_checkpoint_due_monotonic
        ):
            self._save_sequence_state()

    def _capture_paths(
        self,
        kind: str,
        through_day: str | None = None,
    ) -> list[Path]:
        filename = f"{kind}.jsonl"
        paths = []
        for path in self.raw_dir.glob(f"day*/{filename}"):
            day_name = path.parent.name
            day = day_name[3:] if day_name.startswith("day") else ""
            if (
                len(day) == 8
                and day.isdigit()
                and (through_day is None or day <= through_day)
            ):
                paths.append(path)
        return sorted(paths)

    def _cursor_for_path(self, path: Path) -> dict:
        metadata = path.lstat()
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"restart cursor path is not a regular file: {path}")
        offset = metadata.st_size
        with open_regular_read(path) as fh:
            if offset:
                fh.seek(offset - 1)
                if fh.read(1) != b"\n":
                    lower_bound = max(
                        0,
                        offset - RESTART_TAIL_MAX_BYTES,
                    )
                    position = offset
                    found = False
                    while position > lower_bound:
                        size = min(
                            RESTART_TAIL_BLOCK_BYTES,
                            position - lower_bound,
                        )
                        position -= size
                        fh.seek(position)
                        block = fh.read(size)
                        newline = block.rfind(b"\n")
                        if newline >= 0:
                            offset = position + newline + 1
                            found = True
                            break
                    if not found:
                        if metadata.st_size <= RESTART_TAIL_MAX_BYTES:
                            offset = 0
                        else:
                            raise ValueError(
                                f"partial restart record exceeds "
                                f"{RESTART_TAIL_MAX_BYTES} bytes: {path}"
                            )
            anchor_start = max(0, offset - RESTART_CURSOR_ANCHOR_BYTES)
            fh.seek(anchor_start)
            anchor = fh.read(offset - anchor_start)
        return {
            "path": path.relative_to(self.raw_dir).as_posix(),
            "offset": offset,
            "anchor_sha256": hashlib.sha256(anchor).hexdigest(),
        }

    def _validate_cursor(
        self,
        kind: str,
        cursor: object,
    ) -> tuple[Path, int]:
        if not isinstance(cursor, dict) or set(cursor) != {
            "path",
            "offset",
            "anchor_sha256",
        }:
            raise ValueError(f"invalid {kind} restart cursor schema")
        relative = cursor["path"]
        offset = cursor["offset"]
        anchor_sha256 = cursor["anchor_sha256"]
        if (
            not isinstance(relative, str)
            or not isinstance(offset, int)
            or isinstance(offset, bool)
            or offset < 0
            or not isinstance(anchor_sha256, str)
            or len(anchor_sha256) != 64
        ):
            raise ValueError(f"invalid {kind} restart cursor value")
        relative_path = Path(relative)
        parts = relative_path.parts
        expected_name = f"{kind}.jsonl"
        if (
            relative_path.is_absolute()
            or len(parts) != 2
            or parts[1] != expected_name
            or not parts[0].startswith("day")
            or len(parts[0]) != 11
            or not parts[0][3:].isdigit()
        ):
            raise ValueError(f"invalid {kind} restart cursor path")
        path = self.raw_dir / relative_path
        metadata = path.lstat()
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"{kind} restart cursor is not a regular file")
        if offset > metadata.st_size:
            raise ValueError(f"{kind} restart cursor exceeds file size")
        anchor_start = max(0, offset - RESTART_CURSOR_ANCHOR_BYTES)
        with open_regular_read(path) as fh:
            fh.seek(anchor_start)
            anchor = fh.read(offset - anchor_start)
        if offset and not anchor.endswith(b"\n"):
            raise ValueError(f"{kind} restart cursor is not at a record boundary")
        if hashlib.sha256(anchor).hexdigest() != anchor_sha256:
            raise ValueError(f"{kind} restart cursor anchor mismatch")
        return path, offset

    def _rebuild_untrusted_sequence_state(self, reason: str) -> None:
        """Rebuild replayable state without clearing durable storage-cap pauses."""
        self._last_trade_seq = {}
        self._last_trade_source_ts = {}
        self._ticker_last_bucket = {}
        self._required_last_source_timestamp_ms = {}
        self._sequence_state_cursors = {}
        self._checkpoint_source_timestamp_rejected = 0
        self._sequence_state_rebuild_reason = reason
        self._sequence_state_dirty = True

    def _load_sequence_state(self, day: str) -> bool:
        path = self.raw_dir / TRADE_SEQUENCE_STATE_FILE
        try:
            if os.lstat(path).st_size > TRADE_SEQUENCE_STATE_MAX_BYTES:
                raise OSError(
                    "sequence checkpoint exceeds size limit; "
                    "cannot restore durable capture state"
                )
            with open_regular_read(path) as fh:
                payload = strict_json_loads(fh.read().decode("utf-8"))
        except FileNotFoundError:
            return False
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
            NonFiniteJsonConstantError,
            DuplicateJsonKeyError,
        ) as exc:
            raise self._sequence_state_failure("load sequence checkpoint", exc)
        try:
            if not isinstance(payload, dict):
                raise ValueError("sequence checkpoint root is not an object")
            schema_version = payload.get("schema_version")
            if type(schema_version) is not int:
                raise ValueError("invalid sequence checkpoint schema version")
            if schema_version in {2, 3, 4, 5}:
                self._rebuild_untrusted_sequence_state(
                    f"legacy_schema_v{schema_version}"
                )
                return False
            if schema_version != TRADE_SEQUENCE_STATE_SCHEMA_VERSION:
                raise ValueError("unsupported sequence checkpoint schema")
            checkpoint_sample_ms = payload.get("ticker_sample_ms")
            if (
                type(checkpoint_sample_ms) is not int
                or checkpoint_sample_ms <= 0
            ):
                raise ValueError("invalid checkpoint ticker sample interval")
            if (
                not isinstance(payload.get("sequences"), dict)
                or not isinstance(payload.get("ticker_buckets"), dict)
                or not isinstance(payload.get("required_source_timestamps"), dict)
                or not isinstance(payload.get("trade_source_timestamps"), dict)
                or not isinstance(payload.get("cursors"), dict)
                or not isinstance(payload.get("ticker_cap_state"), dict)
            ):
                raise ValueError("unsupported sequence checkpoint schema")
            ticker_cap_state = payload["ticker_cap_state"]
            if set(ticker_cap_state) != {"day_utc", "total_exhausted"}:
                raise ValueError("invalid ticker cap checkpoint schema")
            ticker_cap_day = ticker_cap_state["day_utc"]
            if ticker_cap_day is not None:
                if (
                    not isinstance(ticker_cap_day, str)
                    or len(ticker_cap_day) != 8
                    or not ticker_cap_day.isdigit()
                ):
                    raise ValueError("invalid ticker cap checkpoint day")
                datetime.strptime(ticker_cap_day, "%Y%m%d")
            ticker_total_cap_exhausted = ticker_cap_state[
                "total_exhausted"
            ]
            if type(ticker_total_cap_exhausted) is not bool:
                raise ValueError("invalid ticker total cap checkpoint state")
            self._ticker_cap_day = ticker_cap_day
            self._ticker_total_cap_exhausted = ticker_total_cap_exhausted
            if checkpoint_sample_ms != self.ticker_sample_ms:
                self._rebuild_untrusted_sequence_state(
                    "ticker_sample_ms_mismatch"
                )
                return False
            sequences = {}
            for instrument_name, sequence in payload["sequences"].items():
                if (
                    self._option_instrument_currency(instrument_name) is None
                    or not isinstance(sequence, int)
                    or isinstance(sequence, bool)
                    or sequence <= 0
                ):
                    raise ValueError("invalid sequence checkpoint entry")
                sequences[instrument_name] = sequence
            trade_source_timestamps = {}
            for (
                instrument_name,
                source_ts_ms,
            ) in payload["trade_source_timestamps"].items():
                if (
                    instrument_name not in sequences
                    or not isinstance(source_ts_ms, int)
                    or isinstance(source_ts_ms, bool)
                    or source_ts_ms <= 0
                ):
                    raise ValueError(
                        "invalid trade source timestamp entry"
                    )
                trade_source_timestamps[instrument_name] = source_ts_ms
            if set(trade_source_timestamps) != set(sequences):
                raise ValueError(
                    "trade source timestamps do not cover every sequence"
                )
            ticker_buckets = {}
            checkpoint_rejected = 0
            future_limit_ms = (
                int(time.time() * 1000)
                + REQUIRED_SOURCE_MAX_SKEW_MS
            )
            max_ticker_bucket = (
                future_limit_ms // self.ticker_sample_ms
            )
            for instrument_name, bucket in payload["ticker_buckets"].items():
                if (
                    self._option_instrument_currency(instrument_name) is None
                    or not isinstance(bucket, int)
                    or isinstance(bucket, bool)
                    or bucket < 0
                ):
                    raise ValueError("invalid ticker checkpoint entry")
                if bucket > max_ticker_bucket:
                    checkpoint_rejected += 1
                    continue
                ticker_buckets[instrument_name] = bucket
            required_source_timestamps = {}
            for (
                channel,
                source_ts_ms,
            ) in payload["required_source_timestamps"].items():
                valid_channel = (
                    isinstance(channel, str)
                    and (
                        (
                            channel.startswith("ticker.")
                            and channel.endswith(".agg2")
                            and self._option_instrument_currency(
                                channel[len("ticker."):-len(".agg2")]
                            )
                            is not None
                        )
                        or channel.startswith("deribit_price_index.")
                    )
                )
                if (
                    not valid_channel
                    or not is_finite_number(source_ts_ms)
                    or source_ts_ms <= 0
                ):
                    raise ValueError(
                        "invalid required source timestamp entry"
                    )
                if source_ts_ms > future_limit_ms:
                    checkpoint_rejected += 1
                    continue
                required_source_timestamps[channel] = source_ts_ms
            for instrument_name, bucket in ticker_buckets.items():
                channel = f"ticker.{instrument_name}.agg2"
                required_source_timestamps.setdefault(
                    channel,
                    bucket * self.ticker_sample_ms,
                )
            cursor_payloads = payload["cursors"]
            if not set(cursor_payloads).issubset({"trades", "ticker"}):
                raise ValueError("invalid restart cursor kind")
            cursors = {
                kind: self._validate_cursor(kind, cursor)
                for kind, cursor in cursor_payloads.items()
            }
            for kind in ("trades", "ticker"):
                if self._capture_paths(kind, day) and kind not in cursors:
                    raise ValueError(
                        f"checkpoint lacks cursor for existing {kind} files"
                    )
        except (OSError, TypeError, ValueError) as exc:
            raise self._sequence_state_failure(
                "validate sequence checkpoint",
                exc,
            )
        self._last_trade_seq = sequences
        self._last_trade_source_ts = trade_source_timestamps
        self._ticker_last_bucket = ticker_buckets
        self._required_last_source_timestamp_ms = (
            required_source_timestamps
        )
        self._checkpoint_source_timestamp_rejected = checkpoint_rejected
        self._sequence_state_cursors = cursors
        self._sequence_state_dirty = checkpoint_rejected > 0
        return True

    def _save_sequence_state(self) -> bool:
        if (
            not self._sequence_state_dirty
            or self._rehydration_in_progress
        ):
            return False
        path = self.raw_dir / TRADE_SEQUENCE_STATE_FILE
        temporary_path = self.raw_dir / f".{TRADE_SEQUENCE_STATE_FILE}.tmp"
        try:
            day = utc_day(int(time.time() * 1000))
            cursor_payloads = {}
            for kind in ("trades", "ticker"):
                capture_paths = self._capture_paths(kind, day)
                if capture_paths:
                    cursor_payloads[kind] = self._cursor_for_path(
                        capture_paths[-1]
                    )
            payload = {
                "schema_version": TRADE_SEQUENCE_STATE_SCHEMA_VERSION,
                "ticker_sample_ms": self.ticker_sample_ms,
                "sequences": self._last_trade_seq,
                "trade_source_timestamps": self._last_trade_source_ts,
                "ticker_buckets": self._ticker_last_bucket,
                "ticker_cap_state": {
                    "day_utc": self._ticker_cap_day,
                    "total_exhausted": self._ticker_total_cap_exhausted,
                },
                "required_source_timestamps": (
                    self._required_last_source_timestamp_ms
                ),
                "cursors": cursor_payloads,
                "updated_utc": utc_iso(int(time.time() * 1000)),
            }
            self.raw_dir.mkdir(parents=True, exist_ok=True)
            with create_regular_temp(temporary_path) as fh:
                json.dump(
                    payload,
                    fh,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
                fh.write("\n")
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(temporary_path, path)
            fsync_directory(self.raw_dir)
        except (OSError, TypeError, ValueError) as exc:
            raise self._sequence_state_failure("save sequence checkpoint", exc)
        self._sequence_state_cursors = {
            kind: (
                self.raw_dir / cursor["path"],
                cursor["offset"],
            )
            for kind, cursor in cursor_payloads.items()
        }
        self._sequence_state_dirty = False
        self._sequence_state_bytes_since_checkpoint = 0
        self._sequence_state_checkpoint_due_monotonic = (
            time.monotonic() + RESTART_CHECKPOINT_TIME_INTERVAL
        )
        return True

    def _check_rehydrate_abort(self) -> None:
        if self.stop or stop_requested():
            raise ConnectionError("stop requested during restart replay")
        if self._deadline_expired():
            raise ConnectionError("deadline reached during restart replay")

    def _replay_json_records(
        self,
        path: Path,
        offset: int,
        stats: dict,
        byte_counter: str,
    ):
        self._check_rehydrate_abort()
        try:
            fh = open_regular_read(path)
        except FileNotFoundError:
            return
        with fh:
            fh.seek(offset)
            fragments = []
            record_bytes = 0
            discard_record = False
            while True:
                self._check_rehydrate_abort()
                fragment = fh.readline(RESTART_TAIL_BLOCK_BYTES)
                if not fragment:
                    break
                if (
                    byte_counter == "cursor_replay_bytes"
                    and stats[byte_counter] + len(fragment)
                    > RESTART_CURSOR_MAX_REPLAY_BYTES
                ):
                    raise self._sequence_state_failure(
                        "cursor suffix replay",
                        ValueError(
                            f"exceeds {RESTART_CURSOR_MAX_REPLAY_BYTES} bytes"
                        ),
                    )
                stats[byte_counter] += len(fragment)
                if byte_counter != "metadata_replay_bytes":
                    self._sequence_state_dirty = True
                if discard_record:
                    if fragment.endswith(b"\n"):
                        discard_record = False
                    continue
                record_bytes += len(fragment)
                if record_bytes > RESTART_TAIL_BLOCK_BYTES * 16:
                    stats["malformed"] += 1
                    fragments.clear()
                    record_bytes = 0
                    discard_record = not fragment.endswith(b"\n")
                    continue
                fragments.append(fragment)
                if not fragment.endswith(b"\n"):
                    continue
                raw_record = b"".join(fragments).strip()
                fragments.clear()
                record_bytes = 0
                if not raw_record:
                    continue
                try:
                    record = strict_json_loads(raw_record)
                except (
                    UnicodeError,
                    json.JSONDecodeError,
                    NonFiniteJsonConstantError,
                    DuplicateJsonKeyError,
                ):
                    stats["malformed"] += 1
                    continue
                if not isinstance(record, dict):
                    stats["malformed"] += 1
                    continue
                yield record
            if fragments:
                raw_record = b"".join(fragments).strip()
                if raw_record:
                    try:
                        record = strict_json_loads(raw_record)
                    except (
                        UnicodeError,
                        json.JSONDecodeError,
                        NonFiniteJsonConstantError,
                        DuplicateJsonKeyError,
                    ):
                        stats["malformed"] += 1
                    else:
                        if isinstance(record, dict):
                            yield record
                        else:
                            stats["malformed"] += 1

    def _persisted_option_trade_identity(
        self,
        record: object,
    ) -> tuple[str, int, str] | None:
        if (
            not isinstance(record, dict)
            or record.get("record_type") != "option_trade"
        ):
            return None
        trade = record.get("trade")
        metadata_status = record.get("instrument_metadata_status")
        source_channel = None
        if metadata_status is not None:
            if metadata_status != "unresolved_at_receipt":
                return None
            source_channel = record.get("source_channel")
            if not isinstance(source_channel, str):
                return None
        elif "source_channel" in record:
            return None
        if (
            self._option_trade_payload_error_field(
                trade,
                source_channel,
            )
            is not None
        ):
            return None
        assert isinstance(trade, dict)
        dedup_id = record.get("dedup_id")
        if dedup_id != self._option_trade_dedup_id(trade):
            return None
        received_ms = record.get("received_ts_ms")
        if (
            not isinstance(received_ms, int)
            or isinstance(received_ms, bool)
            or received_ms <= 0
            or abs(trade["timestamp"] - received_ms)
            > REQUIRED_SOURCE_MAX_SKEW_MS
        ):
            return None
        return (
            trade["instrument_name"],
            trade["trade_seq"],
            dedup_id,
        )

    def _persisted_option_ticker_identity(
        self,
        record: object,
    ) -> tuple[str, int | float] | None:
        if (
            not isinstance(record, dict)
            or record.get("record_type") != "option_ticker"
        ):
            return None
        instrument_name = record.get("instrument_name")
        ticker = record.get("ticker")
        if (
            not isinstance(instrument_name, str)
            or self._option_instrument_currency(instrument_name) is None
            or self._option_ticker_schema_error_field(
                instrument_name,
                ticker,
            )
            is not None
        ):
            return None
        assert isinstance(ticker, dict)
        timestamp = ticker["timestamp"]
        received_ms = record.get("received_ts_ms")
        sample_interval_ms = record.get("sample_interval_ms")
        if (
            not is_finite_number(received_ms)
            or received_ms <= 0
            or abs(timestamp - received_ms)
            > REQUIRED_SOURCE_MAX_SKEW_MS
            or not isinstance(sample_interval_ms, int)
            or isinstance(sample_interval_ms, bool)
            or sample_interval_ms != self.ticker_sample_ms
            or record.get("dedup_id")
            != self._option_ticker_dedup_id(
                instrument_name,
                timestamp,
            )
        ):
            return None
        return instrument_name, timestamp

    def _read_trade_sequences(
        self,
        path: Path,
        stats: dict,
        offset: int,
        byte_counter: str,
        seen_last_dedup: dict[str, str] | None = None,
    ) -> None:
        seen_last_dedup = (
            {} if seen_last_dedup is None else seen_last_dedup
        )
        for rec in self._replay_json_records(
            path,
            offset,
            stats,
            byte_counter,
        ):
            identity = self._persisted_option_trade_identity(rec)
            if identity is None:
                stats["trade_sequence_replay_rejected"] += 1
                continue
            instrument_name, sequence, dedup_id = identity
            source_ts_ms = rec["trade"]["timestamp"]
            prior_sequence = self._last_trade_seq.get(instrument_name)
            prior_dedup = seen_last_dedup.get(instrument_name)
            prior_source_ts = self._last_trade_source_ts.get(
                instrument_name
            )
            if prior_sequence is not None and (
                sequence < prior_sequence
                or (
                    sequence == prior_sequence
                    and prior_dedup is not None
                    and prior_dedup != dedup_id
                )
                or (
                    sequence > prior_sequence
                    and prior_source_ts is not None
                    and source_ts_ms < prior_source_ts
                )
            ):
                stats["trade_sequence_replay_rejected"] += 1
                continue
            seen_last_dedup[instrument_name] = dedup_id
            self._merge_trade_sequence(
                instrument_name,
                sequence,
                source_ts_ms,
            )

    def _read_ticker_buckets(
        self,
        path: Path,
        stats: dict,
        offset: int,
        byte_counter: str,
    ) -> None:
        for rec in self._replay_json_records(
            path,
            offset,
            stats,
            byte_counter,
        ):
            identity = self._persisted_option_ticker_identity(rec)
            if identity is None:
                stats["ticker_replay_rejected"] += 1
                if isinstance(rec, dict):
                    ticker = rec.get("ticker")
                    timestamp = (
                        ticker.get("timestamp")
                        if isinstance(ticker, dict)
                        else None
                    )
                    received_ms = rec.get("received_ts_ms")
                    if (
                        not is_finite_number(timestamp)
                        or timestamp <= 0
                        or not is_finite_number(received_ms)
                        or received_ms <= 0
                        or abs(timestamp - received_ms)
                        > REQUIRED_SOURCE_MAX_SKEW_MS
                    ):
                        stats["source_timestamp_rejected"] += 1
                continue
            instrument_name, timestamp = identity
            channel = f"ticker.{instrument_name}.agg2"
            prior = self._required_last_source_timestamp_ms.get(channel)
            bucket = int(timestamp) // self.ticker_sample_ms
            if bucket > self._ticker_last_bucket.get(instrument_name, -1):
                self._ticker_last_bucket[instrument_name] = bucket
                self._sequence_state_dirty = True
            if prior is None or timestamp > prior:
                self._required_last_source_timestamp_ms[channel] = timestamp
                self._sequence_state_dirty = True

    def _persisted_instrument_member(
        self,
        rec: dict,
    ) -> tuple[str, str, str] | None:
        """Return (snapshot_id, instrument_name, currency) for a valid member."""
        snapshot_id = rec.get("snapshot_id")
        received_ms = rec.get("received_ts_ms")
        if (
            not isinstance(snapshot_id, str)
            or not snapshot_id
            or not isinstance(received_ms, int)
            or isinstance(received_ms, bool)
            or received_ms <= 0
        ):
            return None
        record_type = rec.get("record_type")
        if record_type == "instrument":
            instrument = rec.get("instrument")
            if not isinstance(instrument, dict):
                return None
            instrument_name = instrument.get("instrument_name")
            currency = self._option_instrument_currency(instrument_name)
            if currency is None:
                return None
            if (
                not self._instrument_row_matches_schema(
                    instrument,
                    currency.lower(),
                )
                or rec.get("instrument_fingerprint")
                != self._instrument_fingerprint(instrument)
                or rec.get("dedup_id")
                != self._instrument_dedup_id(instrument, received_ms)
            ):
                return None
            return snapshot_id, instrument_name, currency.lower()
        if record_type == "instrument_removed":
            instrument_name = rec.get("instrument_name")
            currency = self._option_instrument_currency(instrument_name)
            prior_fingerprint = rec.get("prior_instrument_fingerprint")
            if (
                currency is None
                or not isinstance(prior_fingerprint, str)
                or len(prior_fingerprint) != 64
                or any(
                    character not in "0123456789abcdef"
                    for character in prior_fingerprint
                )
                or rec.get("reason")
                != "absent_from_authoritative_snapshot"
                or rec.get("dedup_id")
                != self._instrument_removal_dedup_id(
                    instrument_name,
                    prior_fingerprint,
                    snapshot_id,
                )
            ):
                return None
            return snapshot_id, instrument_name, currency.lower()
        return None

    def _instrument_snapshot_commit_publishes(
        self,
        rec: dict,
        snapshot_id: str,
        staged: list[tuple[dict, str, str]],
    ) -> bool:
        """Whether one commit record proves a complete durable snapshot."""
        received_ms = rec.get("received_ts_ms")
        currencies = rec.get("currencies")
        member_count = rec.get("member_count")
        if (
            not staged
            or rec.get("snapshot_id") != snapshot_id
            or not isinstance(received_ms, int)
            or isinstance(received_ms, bool)
            or received_ms <= 0
            or rec.get("dedup_id")
            != sha256_id(f"instr|commit|{snapshot_id}")
            or not isinstance(currencies, list)
            or not currencies
            or currencies != sorted(set(currencies))
            or any(currency not in CURRENCIES for currency in currencies)
            or not isinstance(member_count, int)
            or isinstance(member_count, bool)
            or member_count != len(staged)
        ):
            return False
        if len({name for _, name, _ in staged}) != len(staged):
            return False
        if any(currency not in currencies for _, _, currency in staged):
            return False
        return rec.get("member_digest") == (
            self._instrument_snapshot_member_digest(
                [record for record, _, _ in staged]
            )
        )

    def _publish_replayed_instrument_snapshot(
        self,
        staged: list[tuple[dict, str, str]],
    ) -> None:
        for record, instrument_name, currency in staged:
            if currency not in self.currencies:
                continue
            if record["record_type"] == "instrument":
                self.instruments[instrument_name] = record["instrument"]
                self._persisted_instrument_fingerprints[
                    instrument_name
                ] = record["instrument_fingerprint"]
                continue
            self.instruments.pop(instrument_name, None)
            self._persisted_instrument_fingerprints.pop(
                instrument_name,
                None,
            )

    def _read_instrument_state(self, path: Path, stats: dict) -> None:
        """Restore only fully committed durable instrument snapshots."""
        open_snapshot: str | None = None
        staged: list[tuple[dict, str, str]] = []
        poisoned = False
        for rec in self._replay_json_records(
            path,
            0,
            stats,
            "metadata_replay_bytes",
        ):
            record_type = rec.get("record_type")
            if record_type == "instrument_snapshot_commit":
                if (
                    open_snapshot is not None
                    and not poisoned
                    and self._instrument_snapshot_commit_publishes(
                        rec,
                        open_snapshot,
                        staged,
                    )
                ):
                    self._publish_replayed_instrument_snapshot(staged)
                else:
                    stats["metadata_replay_rejected"] += 1
                    stats["metadata_uncommitted_snapshots"] += 1
                    stats["metadata_uncommitted_records"] += len(staged)
                open_snapshot = None
                staged = []
                poisoned = False
                continue
            if record_type not in ("instrument", "instrument_removed"):
                continue
            member = self._persisted_instrument_member(rec)
            snapshot_id = (
                member[0] if member is not None else rec.get("snapshot_id")
            )
            if not isinstance(snapshot_id, str) or not snapshot_id:
                stats["metadata_replay_rejected"] += 1
                continue
            if snapshot_id != open_snapshot:
                if open_snapshot is not None:
                    stats["metadata_uncommitted_snapshots"] += 1
                    stats["metadata_uncommitted_records"] += len(staged)
                open_snapshot = snapshot_id
                staged = []
                poisoned = False
            if member is None:
                stats["metadata_replay_rejected"] += 1
                poisoned = True
                continue
            staged.append((rec, member[1], member[2]))
        if open_snapshot is not None:
            stats["metadata_uncommitted_snapshots"] += 1
            stats["metadata_uncommitted_records"] += len(staged)

    def _read_recent_trade_ids(
        self,
        path: Path,
        limit: int,
        stats: dict,
    ) -> list[str]:
        if limit <= 0:
            return []
        self._check_rehydrate_abort()
        try:
            fh = open_regular_read(path)
        except FileNotFoundError:
            return []
        with fh:
            fh.seek(0, os.SEEK_END)
            position = fh.tell()
            chunks: deque[bytes] = deque()
            bytes_read = 0
            ids = []
            newline_count = 0
            while position > 0 and bytes_read < RESTART_TAIL_MAX_BYTES:
                self._check_rehydrate_abort()
                size = min(
                    RESTART_TAIL_BLOCK_BYTES,
                    position,
                    RESTART_TAIL_MAX_BYTES - bytes_read,
                )
                position -= size
                fh.seek(position)
                block = fh.read(size)
                chunks.appendleft(block)
                bytes_read += len(block)
                stats["tail_replay_bytes"] += len(block)
                newline_count += block.count(b"\n")
                data = b"".join(chunks)
                if newline_count < limit and position > 0:
                    continue
                lines = data.splitlines()
                if position > 0 and lines:
                    lines = lines[1:]
                ids = []
                for raw_line in lines:
                    try:
                        rec = strict_json_loads(raw_line)
                    except (
                        UnicodeError,
                        json.JSONDecodeError,
                        NonFiniteJsonConstantError,
                        DuplicateJsonKeyError,
                    ):
                        stats["malformed"] += 1
                        continue
                    identity = self._persisted_option_trade_identity(rec)
                    if identity is None:
                        stats["trade_tail_replay_rejected"] += 1
                        continue
                    ids.append(identity[2])
                if len(ids) >= limit:
                    break
            return ids[-limit:]

    def _cursor_replay_paths(
        self,
        kind: str,
        day: str,
    ) -> list[tuple[Path, int]]:
        paths = self._capture_paths(kind, day)
        cursor_path, offset = self._sequence_state_cursors[kind]
        try:
            cursor_index = paths.index(cursor_path)
        except ValueError as exc:
            raise ValueError(
                f"{kind} restart cursor is outside the active history"
            ) from exc
        replay_paths = paths[cursor_index:]
        if len(replay_paths) > RESTART_CURSOR_MAX_REPLAY_FILES:
            raise ValueError(
                f"{kind} restart cursor spans {len(replay_paths)} files"
            )
        return [
            (path, offset if index == 0 else 0)
            for index, path in enumerate(replay_paths)
        ]

    def rehydrate_state(self, recent_trade_ids: int = 20_000) -> dict:
        """Restore restart state before REST seed or subscriptions.

        A validated checkpoint replays only append suffixes. Missing state gets
        one interruptible all-history trade bootstrap; invalid state is
        terminal rather than silently widening startup replay.
        """
        if self._rehydrated is not None:
            return dict(self._rehydrated)
        if recent_trade_ids < 0:
            raise ValueError("recent_trade_ids must be nonnegative")
        self._rehydration_in_progress = True
        day = utc_day(int(time.time() * 1000))
        stats = {
            "instrument_ids": 0,
            "trade_ids": 0,
            "ticker_instruments": 0,
            "trade_sequences": 0,
            "sequence_checkpoint_loaded": False,
            "sequence_checkpoint_rebuilt": False,
            "sequence_checkpoint_rebuild_reason": None,
            "sequence_checkpoint_written": False,
            "sequence_history_files_read": 0,
            "cursor_replay_bytes": 0,
            "fallback_replay_bytes": 0,
            "tail_replay_bytes": 0,
            "metadata_replay_bytes": 0,
            "metadata_replay_rejected": 0,
            "metadata_uncommitted_snapshots": 0,
            "metadata_uncommitted_records": 0,
            "source_timestamp_rejected": 0,
            "trade_sequence_replay_rejected": 0,
            "trade_tail_replay_rejected": 0,
            "ticker_replay_rejected": 0,
            "day": day,
            "malformed": 0,
        }
        checkpoint_loaded = self._load_sequence_state(day)
        stats["sequence_checkpoint_loaded"] = checkpoint_loaded
        stats["sequence_checkpoint_rebuilt"] = (
            self._sequence_state_rebuild_reason is not None
        )
        stats["sequence_checkpoint_rebuild_reason"] = (
            self._sequence_state_rebuild_reason
        )
        stats["source_timestamp_rejected"] += (
            self._checkpoint_source_timestamp_rejected
        )
        day_dir = self.raw_dir / f"day{day}"
        self._read_instrument_state(day_dir / "instruments.jsonl", stats)
        stats["instrument_ids"] = len(
            self._persisted_instrument_fingerprints
        )

        trade_paths = self._capture_paths("trades", day)
        seen_last_trade_dedup: dict[str, str] = {}
        if checkpoint_loaded and trade_paths:
            try:
                replay_trade_paths = self._cursor_replay_paths("trades", day)
            except (KeyError, ValueError) as exc:
                raise self._sequence_state_failure(
                    "validate trade restart cursor",
                    exc,
                )
            for path, offset in replay_trade_paths:
                self._read_trade_sequences(
                    path,
                    stats,
                    offset,
                    "cursor_replay_bytes",
                    seen_last_trade_dedup,
                )
                if path.parent.name[3:] < day:
                    stats["sequence_history_files_read"] += 1
        else:
            for path in trade_paths:
                self._read_trade_sequences(
                    path,
                    stats,
                    0,
                    "fallback_replay_bytes",
                    seen_last_trade_dedup,
                )
                if path.parent.name[3:] < day:
                    stats["sequence_history_files_read"] += 1
        if trade_paths:
            tail = self._read_recent_trade_ids(
                trade_paths[-1],
                recent_trade_ids,
                stats,
            )
            for dedup_id in tail:
                self.deduper.commit(dedup_id)
            stats["trade_ids"] = len(tail)

        ticker_paths = self._capture_paths("ticker", day)
        if checkpoint_loaded and ticker_paths:
            try:
                replay_ticker_paths = self._cursor_replay_paths("ticker", day)
            except (KeyError, ValueError) as exc:
                raise self._sequence_state_failure(
                    "validate ticker restart cursor",
                    exc,
                )
            for path, offset in replay_ticker_paths:
                self._read_ticker_buckets(
                    path,
                    stats,
                    offset,
                    "cursor_replay_bytes",
                )
        else:
            active_ticker_path = day_dir / "ticker.jsonl"
            if active_ticker_path in ticker_paths:
                self._read_ticker_buckets(
                    active_ticker_path,
                    stats,
                    0,
                    "fallback_replay_bytes",
                )

        stats["trade_sequences"] = len(self._last_trade_seq)
        stats["ticker_instruments"] = len(self._ticker_last_bucket)
        self._rehydration_in_progress = False
        stats["sequence_checkpoint_written"] = self._save_sequence_state()
        self._instrument_persist_day = day
        self._rehydrated = dict(stats)
        return stats

    # ---------- signal handling ----------
    def _handle_signal(self, signum, _frame):
        try:
            self._requested_stop_reason = signal.Signals(signum).name.lower()
        except (TypeError, ValueError):
            self._requested_stop_reason = "external"
        self.stop = True
        request_stop()
        if self.ws is not None:
            self.ws.close()

    def _terminal_stop_reason(
        self,
        cap_terminal: bool = False,
        internal_terminal: bool = False,
    ) -> str:
        if internal_terminal:
            return "internal-error"
        if self._requested_stop_reason is not None:
            return self._requested_stop_reason
        if self._deadline_expired():
            return "deadline"
        if cap_terminal:
            return "cap"
        if self.sink._io_error is not None:
            return "sink_io_error"
        if self._sequence_state_error is not None:
            return "sequence_state_error"
        if stop_requested():
            return "external"
        if self.stop:
            return "manual"
        return "completed"

    def _emit_session_stop(
        self,
        cap_terminal: bool = False,
        internal_terminal: bool = False,
    ) -> None:
        if self._session_stop_emitted:
            return
        self._log_event(
            "session_stop",
            {
                "reason": self._terminal_stop_reason(
                    cap_terminal,
                    internal_terminal,
                )
            },
        )
        self._session_stop_emitted = True

    def _remaining_seconds(self) -> float:
        if self._deadline_monotonic is None:
            return float("inf")
        return self._deadline_monotonic - time.monotonic()

    def _deadline_expired(self) -> bool:
        return self._remaining_seconds() <= 0

    def _bounded_timeout(self, requested: float) -> float:
        remaining = self._remaining_seconds()
        if remaining <= 0:
            raise TimeoutError("collector deadline reached")
        return min(requested, remaining)

    # ---------- subscriptions ----------
    def _next_id(self) -> int:
        self.req_id += 1
        return self.req_id

    def default_channels(self) -> list[str]:
        """Required public channels for one complete collector session."""
        channels = []
        for cur in self.currencies:
            C = cur.upper()
            channels.append(f"trades.option.{C}.100ms")
            channels.append(f"deribit_price_index.{cur}_usd")
        if self.capture_tickers:
            channels.extend(self._ticker_channels())
        return channels

    def _ticker_channels(self) -> list[str]:
        return [
            f"ticker.{instrument_name}.agg2"
            for instrument_name in sorted(self.instruments)
            if instrument_name.startswith(self._instrument_prefixes)
        ]

    def _raise_invalid_jsonrpc_shape(
        self,
        phase: str,
        reason: str,
        value,
        received_ms: int,
        frame_bytes: int | None = None,
        channel: str | None = None,
    ) -> None:
        detail = {
            "phase": phase,
            "reason": reason,
            "value_type": type(value).__name__,
        }
        if channel is not None:
            detail["channel"] = str(channel)[:128]
        if frame_bytes is not None:
            detail["frame_bytes"] = frame_bytes
        self.stats["protocol_frame_errors"] += 1
        self._log_event(
            "invalid_jsonrpc_shape",
            detail,
            now_ms=received_ms,
        )
        raise DeribitProtocolFrameError(
            f"{phase}: {reason} ({type(value).__name__})"
        )

    def _decode_jsonrpc_object(
        self,
        data: bytes,
        phase: str,
        received_ms: int,
        *,
        audit_malformed: bool = False,
    ) -> dict | None:
        try:
            obj = strict_json_loads(data.decode("utf-8"))
        except NonFiniteJsonConstantError:
            self._raise_invalid_jsonrpc_shape(
                phase,
                "non_finite_json_constant",
                None,
                received_ms,
                frame_bytes=len(data),
            )
        except DuplicateJsonKeyError:
            self._raise_invalid_jsonrpc_shape(
                phase,
                "duplicate_json_member",
                None,
                received_ms,
                frame_bytes=len(data),
            )
        except (UnicodeDecodeError, json.JSONDecodeError):
            if audit_malformed:
                self._log_event(
                    "malformed_json",
                    {"frame_bytes": len(data)},
                    now_ms=received_ms,
                )
            return None
        if not isinstance(obj, dict):
            self._raise_invalid_jsonrpc_shape(
                phase,
                "top_level_not_object",
                obj,
                received_ms,
                frame_bytes=len(data),
            )
        return obj

    def _absorb_unsubscribed_channel(
        self,
        channel: str,
        received_ms: int,
    ) -> bool:
        """Absorb one in-flight frame for a just-unsubscribed channel."""
        now = time.monotonic()
        entry = self._recently_unsubscribed.get(channel)
        for name in [
            name
            for name, (deadline, _audited)
            in self._recently_unsubscribed.items()
            if deadline <= now
        ]:
            del self._recently_unsubscribed[name]
        if entry is None or entry[0] <= now:
            return False
        deadline, audited = entry
        self.stats["unsubscribed_grace_frames"] += 1
        if not audited:
            self._log_event(
                "unsubscribed_channel_frame",
                {"channel": channel[:128]},
                now_ms=received_ms,
            )
            self._recently_unsubscribed[channel] = (deadline, True)
        self.last_message_ms = received_ms
        self.last_message_monotonic = now
        return True

    def _handle_notification(
        self,
        obj: dict,
        received_ms: int,
        *,
        phase: str = "notification",
    ) -> bool:
        """Validate and route one asynchronous JSON-RPC notification."""
        if not isinstance(obj, dict):
            self._raise_invalid_jsonrpc_shape(
                phase,
                "top_level_not_object",
                obj,
                received_ms,
            )
        method = obj.get("method")
        if method == "subscription":
            params = obj.get("params")
            if not isinstance(params, dict):
                self._raise_invalid_jsonrpc_shape(
                    phase,
                    "subscription_params_not_object",
                    params,
                    received_ms,
                )
            channel = params.get("channel")
            if not isinstance(channel, str) or not channel:
                self._raise_invalid_jsonrpc_shape(
                    phase,
                    "subscription_channel_not_string",
                    channel,
                    received_ms,
                )
            if channel not in (
                self.active_channels | self._pending_channels
            ):
                if not self._absorb_unsubscribed_channel(
                    channel,
                    received_ms,
                ):
                    self._raise_invalid_jsonrpc_shape(
                        phase,
                        "unsolicited_channel",
                        channel,
                        received_ms,
                        channel=channel,
                    )
                return True
            received_monotonic = time.monotonic()
            self._ensure_instrument_day_fence(
                received_ms,
                received_monotonic,
            )
            self.msg_seq += 1
            self.last_message_ms = received_ms
            self.last_message_monotonic = received_monotonic
            self._route(
                channel,
                params.get("data"),
                received_ms,
                received_monotonic,
            )
            return True
        if method == "heartbeat":
            self.msg_seq += 1
            self.last_message_ms = received_ms
            self.last_message_monotonic = time.monotonic()
            assert self.ws is not None
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0",
                "method": "public/test",
            }))
            self._log_event("heartbeat_ok", None)
            return True
        if "error" in obj and obj.get("id") is None:
            error = obj["error"]
            if not isinstance(error, dict):
                self._raise_invalid_jsonrpc_shape(
                    phase,
                    "error_not_object",
                    error,
                    received_ms,
                )
            provenance = bounded_error_provenance(error)
            self._log_event(
                "server_error_notification",
                provenance,
            )
            raise ConnectionError(
                "server error notification: "
                f"code={provenance['code']} msg={provenance['msg']}"
            )
        return False

    def _await_ack(
        self,
        rid: int,
        chunk: set[str],
        timeout_s: float = 20.0,
        method: str = "public/subscribe",
    ) -> bool:
        """Require an exact channel-set ack while preserving notifications."""
        self._pending_channels = set(chunk)
        try:
            return self._await_exact_ack(rid, chunk, timeout_s, method)
        finally:
            self._pending_channels = set()

    def _await_exact_ack(
        self,
        rid: int,
        chunk: set[str],
        timeout_s: float,
        method: str,
    ) -> bool:
        """Consume frames until the exact ack, a failure, or the deadline."""
        event_prefix = (
            "subscribe"
            if method == "public/subscribe"
            else "unsubscribe"
        )
        try:
            deadline = time.monotonic() + self._bounded_timeout(timeout_s)
        except TimeoutError:
            self._log_event(f"{event_prefix}_timeout", {"n": len(chunk)})
            return False
        while time.monotonic() < deadline and not self._deadline_expired():
            if stop_requested() or self.stop:
                self._log_event(f"{event_prefix}_aborted_stop")
                return False
            remaining = min(
                deadline - time.monotonic(),
                self._remaining_seconds(),
            )
            if remaining <= 0:
                break
            try:
                assert self.ws is not None
                op, data = self.ws.recv_message(
                    timeout=min(INTERRUPT_POLL_SECONDS, remaining)
                )
            except socket.timeout:
                continue
            if op == OP_CLOSE:
                raise ConnectionError(f"closed during {event_prefix}")
            received_ms = int(time.time() * 1000)
            obj = self._decode_jsonrpc_object(
                data,
                f"{event_prefix}_ack",
                received_ms,
            )
            if obj is None:
                continue
            if obj.get("id") != rid:
                self._handle_notification(
                    obj,
                    received_ms,
                    phase=f"{event_prefix}_ack",
                )
                continue
            if "error" in obj:
                error = obj["error"]
                if not isinstance(error, dict):
                    self._raise_invalid_jsonrpc_shape(
                        f"{event_prefix}_ack",
                        "error_not_object",
                        error,
                        received_ms,
                        frame_bytes=len(data),
                    )
                self._log_event(
                    f"{event_prefix}_error",
                    {"error": error, "n": len(chunk)},
                )
                return False
            result = obj.get("result")
            if (
                not isinstance(result, list)
                or any(
                    not isinstance(channel, str) or not channel
                    for channel in result
                )
            ):
                self._raise_invalid_jsonrpc_shape(
                    f"{event_prefix}_ack",
                    "result_not_channel_list",
                    result,
                    received_ms,
                    frame_bytes=len(data),
                )
            returned = set(result)
            if returned != chunk:
                self._log_event(
                    f"{event_prefix}_set_mismatch",
                    {
                        "missing": sorted(chunk - returned)[:20],
                        "extra": sorted(returned - chunk)[:20],
                        "result_empty": (
                            "result_empty"
                            if not returned
                            else None
                        ),
                    },
                )
                return False
            if method == "public/subscribe":
                self.active_channels |= chunk
            else:
                self.active_channels -= chunk
                grace_deadline = (
                    time.monotonic() + UNSUBSCRIBE_GRACE_SECONDS
                )
                for unsubscribed in chunk:
                    self._recently_unsubscribed[unsubscribed] = (
                        grace_deadline,
                        False,
                    )
            return True
        self._log_event(f"{event_prefix}_timeout", {"n": len(chunk)})
        return False

    def subscribe_all(self) -> bool:
        """Require exact channel acks and heartbeat setup for the full session."""
        channels = self.default_channels()
        chunk_size = 220
        self.active_channels = set()
        self._recently_unsubscribed = {}
        self._required_last_received_monotonic = {}
        ok_all = True
        for i in range(0, len(channels), chunk_size):
            if stop_requested():
                return False
            chunk = channels[i:i + chunk_size]
            rid = self._next_id()
            assert self.ws is not None
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0",
                "id": rid,
                "method": "public/subscribe",
                "params": {"channels": chunk},
            }))
            if not self._await_ack(rid, set(chunk)):
                ok_all = False
                break
        if ok_all and self.heartbeat:
            if stop_requested():
                return False
            hb_id = self._next_id()
            assert self.ws is not None
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0",
                "id": hb_id,
                "method": "public/set_heartbeat",
                "params": {"interval": 10},
            }))
            try:
                deadline = (
                    time.monotonic()
                    + self._bounded_timeout(15.0)
                )
            except TimeoutError:
                self._log_event("heartbeat_setup_timeout")
                return False
            hb_ok = False
            while (
                time.monotonic() < deadline
                and not self._deadline_expired()
                and not hb_ok
            ):
                if stop_requested() or self.stop:
                    return False
                remaining = min(
                    deadline - time.monotonic(),
                    self._remaining_seconds(),
                )
                if remaining <= 0:
                    break
                try:
                    op, data = self.ws.recv_message(
                        timeout=min(INTERRUPT_POLL_SECONDS, remaining)
                    )
                except socket.timeout:
                    continue
                if op == OP_CLOSE:
                    raise ConnectionError("closed during heartbeat setup")
                received_ms = int(time.time() * 1000)
                obj = self._decode_jsonrpc_object(
                    data,
                    "heartbeat_ack",
                    received_ms,
                )
                if obj is None:
                    continue
                if obj.get("id") != hb_id:
                    self._handle_notification(
                        obj,
                        received_ms,
                        phase="heartbeat_ack",
                    )
                    continue
                if "error" in obj:
                    error = obj["error"]
                    if not isinstance(error, dict):
                        self._raise_invalid_jsonrpc_shape(
                            "heartbeat_ack",
                            "error_not_object",
                            error,
                            received_ms,
                            frame_bytes=len(data),
                        )
                    self._log_event(
                        "heartbeat_setup_error",
                        {"error": error},
                    )
                    break
                result = obj.get("result")
                if not isinstance(result, str):
                    self._raise_invalid_jsonrpc_shape(
                        "heartbeat_ack",
                        "result_not_string",
                        result,
                        received_ms,
                        frame_bytes=len(data),
                    )
                if result == "ok":
                    hb_ok = True
            if not hb_ok:
                self._log_event("heartbeat_setup_timeout")
                ok_all = False
        if ok_all and self.active_channels != set(channels):
            self._log_event(
                "subscribe_set_mismatch_final",
                {
                    "missing": sorted(set(channels) - self.active_channels)[:20],
                    "extra": sorted(self.active_channels - set(channels))[:20],
                },
            )
            ok_all = False
        return ok_all

    def _active_ticker_channels(self) -> set[str]:
        return {
            channel
            for channel in self.active_channels
            if channel.startswith("ticker.")
        }

    def _active_required_channels(self) -> set[str]:
        return {
            channel
            for channel in self.active_channels
            if channel.startswith("ticker.")
            or channel.startswith("deribit_price_index.")
        }

    def _fresh_required_channels(
        self,
        expected: set[str],
        now_monotonic: float,
    ) -> set[str]:
        cutoff = now_monotonic - TICKER_STALE_SECONDS
        return {
            channel
            for channel in expected
            if self._required_last_received_monotonic.get(
                channel,
                float("-inf"),
            )
            >= cutoff
        }

    def await_required_coverage(
        self,
        channels: set[str] | None = None,
        timeout_s: float = TICKER_COVERAGE_TIMEOUT_SECONDS,
    ) -> bool:
        """Require a fresh usable payload from every high-rate input."""
        expected = (
            self._active_required_channels()
            if channels is None
            else set(channels)
        )
        if not expected:
            self._log_event("required_coverage_empty")
            return False
        try:
            deadline = (
                time.monotonic()
                + self._bounded_timeout(timeout_s)
            )
        except TimeoutError:
            self._log_event(
                "required_coverage_timeout",
                {
                    "expected": len(expected),
                    "missing": len(expected),
                    "sample": sorted(expected)[:20],
                },
            )
            return False
        while True:
            now_monotonic = time.monotonic()
            fresh = self._fresh_required_channels(
                expected,
                now_monotonic,
            )
            if expected.issubset(fresh):
                break
            if stop_requested() or self.stop:
                self._log_event("required_coverage_aborted_stop")
                return False
            remaining = min(
                deadline - now_monotonic,
                self._remaining_seconds(),
            )
            if remaining <= 0:
                missing = expected - fresh
                self._log_event(
                    "required_coverage_timeout",
                    {
                        "expected": len(expected),
                        "missing": len(missing),
                        "sample": sorted(missing)[:20],
                    },
                )
                return False
            try:
                assert self.ws is not None
                op, data = self.ws.recv_message(
                    timeout=min(INTERRUPT_POLL_SECONDS, remaining)
                )
            except socket.timeout:
                continue
            if op == OP_CLOSE:
                raise ConnectionError(
                    "closed during required channel coverage"
                )
            if op != OP_TEXT:
                continue
            received_ms = int(time.time() * 1000)
            obj = self._decode_jsonrpc_object(
                data,
                "required_coverage",
                received_ms,
            )
            if obj is None:
                continue
            self._handle_notification(
                obj,
                received_ms,
                phase="required_coverage",
            )
        self._log_event(
            "required_coverage_ready",
            {"n_channels": len(expected)},
        )
        self._next_required_health_check = (
            time.monotonic() + TICKER_HEALTH_CHECK_SECONDS
        )
        return True

    def _check_required_continuity(
        self,
        now_monotonic: float | None = None,
    ) -> None:
        now_monotonic = (
            time.monotonic()
            if now_monotonic is None
            else now_monotonic
        )
        if now_monotonic < self._next_required_health_check:
            return
        self._next_required_health_check = (
            now_monotonic + TICKER_HEALTH_CHECK_SECONDS
        )
        expected = self._active_required_channels()
        fresh = self._fresh_required_channels(expected, now_monotonic)
        stale = sorted(expected - fresh)
        if stale:
            self._log_event(
                "required_channel_continuity_gap",
                {
                    "expected": len(expected),
                    "stale": len(stale),
                    "stale_after_s": TICKER_STALE_SECONDS,
                    "sample": sorted(stale)[:20],
                },
            )
            raise ConnectionError(
                f"{len(stale)} required high-rate channels stale"
            )

    def _instrument_row_matches_schema(
        self,
        row: object,
        currency: str,
    ) -> bool:
        if not isinstance(row, dict):
            return False
        identity = self._option_instrument_identity(
            row.get("instrument_name")
        )
        strike = row.get("strike")
        expiration_timestamp = row.get("expiration_timestamp")
        if (
            identity is None
            or not is_finite_number(strike)
            or strike <= 0
            or not isinstance(expiration_timestamp, int)
            or isinstance(expiration_timestamp, bool)
            or expiration_timestamp <= 0
        ):
            return False
        try:
            row_strike = Decimal(str(strike))
            settlement_ms = int(
                identity.expiry.timestamp() * 1000
            ) + OPTION_EXPIRY_SETTLEMENT_MS
        except (
            InvalidOperation,
            OSError,
            OverflowError,
            ValueError,
        ):
            return False
        upper = currency.upper()
        return (
            identity.currency == upper
            and identity.strike == row_strike
            and identity.option_type == row.get("option_type")
            and expiration_timestamp == settlement_ms
            and row.get("kind") == "option"
            and isinstance(row.get("instrument_id"), int)
            and not isinstance(row.get("instrument_id"), bool)
            and row["instrument_id"] > 0
            and isinstance(row.get("creation_timestamp"), int)
            and not isinstance(row.get("creation_timestamp"), bool)
            and row["creation_timestamp"] > 0
            and row["creation_timestamp"] < expiration_timestamp
            and row.get("base_currency") == upper
            and row.get("state") == OPTION_OPEN_STATE
            and row.get("is_active") is True
            and is_finite_number(row.get("contract_size"))
            and row["contract_size"] == OPTION_CONTRACT_SIZE_COINS
            and is_finite_number(row.get("min_trade_amount"))
            and row["min_trade_amount"] > 0
            and is_finite_number(row.get("lot_size"))
            and row["lot_size"] > 0
            and is_finite_number(row.get("tick_size"))
            and row["tick_size"] > 0
            and row.get("settlement_currency") == upper
            and row.get("price_index") == f"{currency.lower()}_usd"
        )

    def _fetch_currency_instrument_rows(
        self,
        currency: str,
        deadline: float | None,
    ) -> list[dict]:
        if stop_requested():
            raise InstrumentSnapshotError(
                currency,
                "StopRequested",
                "stop requested during instrument snapshot",
            )
        try:
            rows = self._filter_not_open_rows(
                fetch_instruments(currency, deadline=deadline),
            )
        except Exception as exc:
            raise InstrumentSnapshotError(
                currency,
                type(exc).__name__,
                f"required {currency.upper()} instrument seed failed",
            ) from exc
        if (
            isinstance(rows, list)
            and len(rows) > INSTRUMENT_SNAPSHOT_MAX_ROWS_PER_CURRENCY
        ):
            raise InstrumentSnapshotError(
                currency,
                "RowLimitExceeded",
                f"required {currency.upper()} instrument seed exceeds "
                f"{INSTRUMENT_SNAPSHOT_MAX_ROWS_PER_CURRENCY} rows",
            )
        if (
            not isinstance(rows, list)
            or not rows
            or any(
                not self._instrument_row_matches_schema(row, currency)
                for row in rows
            )
        ):
            raise InstrumentSnapshotError(
                currency,
                "SchemaMismatch",
                f"required {currency.upper()} instrument seed schema mismatch",
                detail={
                    "candidate_count": (
                        len(rows) if isinstance(rows, list) else None
                    ),
                },
            )
        instrument_names = [
            row["instrument_name"]
            for row in rows
        ]
        if len(instrument_names) != len(set(instrument_names)):
            raise InstrumentSnapshotError(
                currency,
                "DuplicateInstrumentName",
                f"required {currency.upper()} instrument seed "
                "contains duplicate names",
            )
        return rows

    @staticmethod
    def _instrument_row_is_not_open(row: object) -> bool:
        """Whether the venue explicitly marks one typed row as not open."""
        if not isinstance(row, dict):
            return False
        state = row.get("state")
        is_active = row.get("is_active")
        if not isinstance(state, str) or not state:
            return False
        if not isinstance(is_active, bool):
            return False
        return (
            state.lower() in INSTRUMENT_TERMINAL_STATES
            or not is_active
        )

    def _filter_not_open_rows(self, rows: object) -> object:
        """Drop rows the venue marks closed or inactive from the chain."""
        if not isinstance(rows, list):
            return rows
        kept = [
            row
            for row in rows
            if not self._instrument_row_is_not_open(row)
        ]
        self.stats["instrument_rows_not_open"] += len(rows) - len(kept)
        return kept

    def _instrument_snapshot_signature(
        self,
        rows: list[dict],
    ) -> str:
        canonical = [
            (
                row["instrument_name"],
                self._instrument_fingerprint(row),
            )
            for row in sorted(
                rows,
                key=lambda item: item["instrument_name"],
            )
        ]
        return sha256_id(
            strict_json_dumps(canonical, separators=(",", ":"))
        )

    @staticmethod
    def _instrument_snapshot_member_identity(
        record: dict,
    ) -> tuple[str, str, str] | None:
        record_type = record.get("record_type")
        if record_type == "instrument":
            instrument = record.get("instrument")
            instrument_name = (
                instrument.get("instrument_name")
                if isinstance(instrument, dict)
                else None
            )
            fingerprint = record.get("instrument_fingerprint")
        elif record_type == "instrument_removed":
            instrument_name = record.get("instrument_name")
            fingerprint = record.get("prior_instrument_fingerprint")
        else:
            return None
        if (
            not isinstance(instrument_name, str)
            or not instrument_name
            or not isinstance(fingerprint, str)
            or not fingerprint
        ):
            return None
        return record_type, instrument_name, fingerprint

    @classmethod
    def _instrument_snapshot_member_digest(
        cls,
        records: list[dict],
    ) -> str | None:
        identities = []
        for record in records:
            identity = cls._instrument_snapshot_member_identity(record)
            if identity is None:
                return None
            identities.append(list(identity))
        return sha256_id(
            strict_json_dumps(identities, separators=(",", ":"))
        )

    def _instrument_snapshot_commit_record(
        self,
        snapshot_id: str,
        members: list[dict],
        received_ms: int,
    ) -> dict:
        digest = self._instrument_snapshot_member_digest(members)
        if digest is None:
            raise InstrumentSnapshotError(
                "*",
                "SnapshotCommitIdentityMissing",
                "instrument snapshot member identity is incomplete",
            )
        return {
            "record_type": "instrument_snapshot_commit",
            "dedup_id": sha256_id(f"instr|commit|{snapshot_id}"),
            "received_ts_ms": received_ms,
            "snapshot_id": snapshot_id,
            "currencies": sorted(self.currencies),
            "member_count": len(members),
            "member_digest": digest,
        }

    def _fetch_instrument_snapshot(
        self,
        deadline: float | None,
    ) -> dict[str, list[dict]]:
        """Fetch and validate every currency without mutating chain state."""
        if stop_requested():
            raise InstrumentSnapshotError(
                "*",
                "StopRequested",
                "stop requested before instrument snapshot",
            )
        fetched: dict[str, list[dict]] = {}
        for currency in self.currencies:
            rows = self._fetch_currency_instrument_rows(
                currency,
                deadline,
            )
            prefix = f"{currency.upper()}-"
            has_baseline = any(
                name.startswith(prefix)
                for name in self.instruments
            )
            if not has_baseline:
                if len(rows) < INSTRUMENT_INITIAL_MIN_ROWS_PER_CURRENCY:
                    raise InstrumentSnapshotError(
                        currency,
                        "InitialSnapshotImplausible",
                        f"required {currency.upper()} initial instrument "
                        "snapshot is implausibly small",
                        detail={
                            "candidate_count": len(rows),
                            "minimum_count": (
                                INSTRUMENT_INITIAL_MIN_ROWS_PER_CURRENCY
                            ),
                        },
                    )
                confirmation = self._fetch_currency_instrument_rows(
                    currency,
                    deadline,
                )
                if (
                    self._instrument_snapshot_signature(rows)
                    != self._instrument_snapshot_signature(confirmation)
                ):
                    raise InstrumentSnapshotError(
                        currency,
                        "InitialSnapshotConfirmationMismatch",
                        f"required {currency.upper()} initial instrument "
                        "snapshots disagree",
                        detail={
                            "candidate_count": len(rows),
                            "confirmation_count": len(confirmation),
                        },
                    )
                rows = confirmation
            fetched[currency] = rows
        return fetched

    @staticmethod
    def _instrument_removal_is_explained(
        row: dict,
        received_ms: int,
    ) -> bool:
        state = row.get("state")
        if (
            isinstance(state, str)
            and state.lower()
            in INSTRUMENT_TERMINAL_STATES
        ):
            return True
        if row.get("is_active") is False:
            return True
        expiration = row.get("expiration_timestamp")
        return (
            isinstance(expiration, int)
            and not isinstance(expiration, bool)
            and expiration > 0
            and expiration + INSTRUMENT_EXPIRY_SETTLEMENT_GRACE_MS
            <= received_ms
        )

    def _validate_instrument_snapshot_continuity(
        self,
        fetched: dict[str, list[dict]],
        received_ms: int,
    ) -> None:
        candidates: dict[str, frozenset[str]] = {}
        details: dict[str, dict] = {}
        next_pending = dict(
            self._pending_unexpected_instrument_removals
        )
        for currency, rows in fetched.items():
            prefix = f"{currency.upper()}-"
            prior = {
                name: row
                for name, row in self.instruments.items()
                if name.startswith(prefix)
            }
            if not prior:
                next_pending.pop(currency, None)
                continue
            names = {row["instrument_name"] for row in rows}
            missing = set(prior) - names
            unexpected = frozenset(
                name
                for name in missing
                if not self._instrument_removal_is_explained(
                    prior[name],
                    received_ms,
                )
            )
            if not unexpected:
                next_pending.pop(currency, None)
                continue
            detail = {
                "prior_count": len(prior),
                "candidate_count": len(rows),
                "missing_count": len(missing),
                "unexpected_missing": len(unexpected),
                "unexpected_fraction": round(
                    len(unexpected) / len(prior),
                    6,
                ),
                "sample": sorted(unexpected)[:20],
            }
            if (
                len(unexpected) / len(prior)
                > INSTRUMENT_MAX_UNEXPECTED_REMOVAL_FRACTION
            ):
                next_pending.pop(currency, None)
                self._pending_unexpected_instrument_removals = (
                    next_pending
                )
                raise InstrumentSnapshotError(
                    currency,
                    "CatastrophicSnapshotTruncation",
                    f"required {currency.upper()} instrument snapshot "
                    "has catastrophic unexplained omissions",
                    detail=detail,
                )
            candidates[currency] = unexpected
            details[currency] = detail

        unconfirmed = [
            currency
            for currency, candidate in candidates.items()
            if self._pending_unexpected_instrument_removals.get(
                currency
            )
            != candidate
        ]
        for currency in fetched:
            if currency in candidates:
                next_pending[currency] = candidates[currency]
            else:
                next_pending.pop(currency, None)
        if unconfirmed:
            self._pending_unexpected_instrument_removals = next_pending
            currency = unconfirmed[0]
            raise InstrumentSnapshotError(
                currency,
                "UnexpectedRemovalConfirmationRequired",
                f"required {currency.upper()} instrument removals need "
                "an independent confirming snapshot",
                detail=details[currency],
            )
        for currency in candidates:
            next_pending.pop(currency, None)
        self._pending_unexpected_instrument_removals = next_pending

    def _apply_instrument_snapshot(
        self,
        fetched: dict[str, list[dict]],
        received_ms: int,
    ) -> None:
        target_day = utc_day(received_ms)
        unexpected_currencies = set(fetched) - set(self.currencies)
        if unexpected_currencies:
            currency = sorted(unexpected_currencies)[0]
            raise InstrumentSnapshotError(
                currency,
                "UnexpectedCurrency",
                f"instrument snapshot contains unconfigured currency "
                f"{currency.upper()}",
            )
        day_changed = (
            target_day
            != getattr(self, "_instrument_persist_day", None)
        )
        next_fingerprints = (
            {}
            if day_changed
            else {
                name: fingerprint
                for name, fingerprint
                in self._persisted_instrument_fingerprints.items()
                if name.startswith(self._instrument_prefixes)
            }
        )
        next_instruments = {
            name: row
            for name, row in self.instruments.items()
            if name.startswith(self._instrument_prefixes)
        }
        snapshot_id = sha256_id(
            f"instrument-snapshot|{target_day}|{received_ms}"
        )
        records = []
        upsert_count = 0
        removed_count = 0
        for currency, rows in fetched.items():
            if self.stop or stop_requested() or self._deadline_expired():
                raise ConnectionError(
                    "stop or deadline during instrument snapshot apply"
                )
            names = {row["instrument_name"] for row in rows}
            if len(names) != len(rows):
                raise InstrumentSnapshotError(
                    currency,
                    "DuplicateInstrumentName",
                    f"required {currency.upper()} instrument seed "
                    "contains duplicate names",
                )
            prefix = f"{currency.upper()}-"
            stale = {
                instrument_name
                for instrument_name in next_instruments
                if instrument_name.startswith(prefix)
                and instrument_name not in names
            }
            for row in rows:
                if self.stop or stop_requested() or self._deadline_expired():
                    raise ConnectionError(
                        "stop or deadline during instrument snapshot apply"
                    )
                instrument_name = row["instrument_name"]
                fingerprint = self._instrument_fingerprint(row)
                next_instruments[instrument_name] = row
                if next_fingerprints.get(instrument_name) == fingerprint:
                    continue
                next_fingerprints[instrument_name] = fingerprint
                records.append({
                    "record_type": "instrument",
                    "dedup_id": self._instrument_dedup_id(
                        row,
                        received_ms,
                    ),
                    "received_ts_ms": received_ms,
                    "snapshot_id": snapshot_id,
                    "instrument_fingerprint": fingerprint,
                    "instrument": row,
                })
                upsert_count += 1
            for instrument_name in sorted(stale):
                prior_row = next_instruments.pop(instrument_name)
                prior_fingerprint = next_fingerprints.pop(
                    instrument_name,
                    None,
                )
                if prior_fingerprint is None:
                    prior_fingerprint = self._instrument_fingerprint(
                        prior_row
                    )
                records.append({
                    "record_type": "instrument_removed",
                    "dedup_id": self._instrument_removal_dedup_id(
                        instrument_name,
                        prior_fingerprint,
                        snapshot_id,
                    ),
                    "received_ts_ms": received_ms,
                    "snapshot_id": snapshot_id,
                    "instrument_name": instrument_name,
                    "prior_instrument_fingerprint": prior_fingerprint,
                    "reason": "absent_from_authoritative_snapshot",
                })
                removed_count += 1
        if self.stop or stop_requested() or self._deadline_expired():
            raise ConnectionError(
                "stop or deadline before instrument snapshot commit"
            )
        member_count = len(records)
        if member_count:
            records.append(
                self._instrument_snapshot_commit_record(
                    snapshot_id,
                    records,
                    received_ms,
                )
            )
        self.sink.write_batch(
            "instruments",
            records,
            now_ms=received_ms,
        )
        self.instruments = next_instruments
        self._instrument_persist_day = target_day
        self._persisted_instrument_fingerprints = next_fingerprints
        self.stats["instrument_records"] += member_count
        self.stats["instrument_seed_new"] = (
            self.stats.get("instrument_seed_new", 0)
            + upsert_count
        )
        self.stats["instrument_seed_removed"] = (
            self.stats.get("instrument_seed_removed", 0)
            + removed_count
        )
        stale_names = (
            set(self._last_trade_seq)
            | set(self._ticker_last_bucket)
        ) - set(self.instruments)
        if stale_names:
            for instrument_name in stale_names:
                self._last_trade_seq.pop(instrument_name, None)
                self._last_trade_source_ts.pop(instrument_name, None)
                self._ticker_last_bucket.pop(instrument_name, None)
                self._required_last_source_timestamp_ms.pop(
                    f"ticker.{instrument_name}.agg2",
                    None,
                )
            self._sequence_state_dirty = True
        self._save_sequence_state()
        self._log_event(
            "instrument_seed_done",
            {"batches": list(self.currencies), "n": len(self.instruments)},
            now_ms=received_ms,
            include_connection=False,
        )
    def _ensure_instrument_day_fence(
        self,
        received_ms: int,
        received_monotonic: float | None = None,
    ) -> bool:
        """Persist cached reconstruction metadata before a new-day boundary."""
        target_day = utc_day(received_ms)
        if (
            getattr(self, "_instrument_persist_day", None)
            == target_day
        ):
            return False
        cached_snapshot = {}
        for currency in self.currencies:
            prefix = f"{currency.upper()}-"
            rows = [
                self.instruments[instrument_name]
                for instrument_name in sorted(self.instruments)
                if instrument_name.startswith(prefix)
                and isinstance(self.instruments[instrument_name], dict)
                and self.instruments[instrument_name].get("kind") == "option"
            ]
            if not rows:
                raise InstrumentSnapshotError(
                    currency,
                    "CachedSnapshotIncomplete",
                    f"no cached {currency.upper()} option metadata "
                    f"for UTC-day fence {target_day}",
                )
            cached_snapshot[currency] = rows
        self._apply_instrument_snapshot(cached_snapshot, received_ms)
        now = (
            time.monotonic()
            if received_monotonic is None
            else received_monotonic
        )
        self.last_utc_day = target_day
        self._next_chain_refresh_attempt = now
        self.stats["instrument_day_carries"] = (
            self.stats.get("instrument_day_carries", 0) + 1
        )
        return True


    def _record_instrument_snapshot_error(
        self,
        event: str,
        exc: Exception,
        received_ms: int,
    ) -> None:
        detail = {"error_type": type(exc).__name__}
        if isinstance(exc, InstrumentSnapshotError):
            detail = {
                "currency": exc.currency,
                "error_type": exc.error_type,
                **exc.detail,
            }
        self._log_event(
            event,
            detail,
            now_ms=received_ms,
            include_connection=False,
        )

    def _fetch_instrument_snapshot_interruptibly(
        self,
    ) -> dict[str, list[dict]]:
        result: dict[str, object] = {}
        done = threading.Event()

        def fetch_worker() -> None:
            try:
                result["fetched"] = self._fetch_instrument_snapshot(
                    self._deadline_monotonic,
                )
            except BaseException as exc:
                result["error"] = exc
            finally:
                done.set()

        threading.Thread(
            target=fetch_worker,
            name="deribit-startup-seed",
            daemon=True,
        ).start()
        while not done.wait(INTERRUPT_POLL_SECONDS):
            if self.stop or stop_requested():
                raise InstrumentSnapshotError(
                    "*",
                    "StopRequested",
                    "stop requested during instrument seed",
                )
            if self._deadline_expired():
                raise TimeoutError(
                    "collector deadline reached during instrument seed"
                )
        if self.stop or stop_requested():
            raise InstrumentSnapshotError(
                "*",
                "StopRequested",
                "stop requested after instrument seed",
            )
        error = result.get("error")
        if error is not None:
            raise error
        fetched = result.get("fetched")
        assert isinstance(fetched, dict)
        return fetched

    def _mark_authoritative_chain_fresh(
        self,
        received_ms: int,
        now_monotonic: float | None = None,
    ) -> None:
        now_monotonic = (
            time.monotonic()
            if now_monotonic is None
            else now_monotonic
        )
        recovered = self._instrument_chain_stale_latched
        prior_success = self._last_authoritative_chain_refresh_monotonic
        self._last_authoritative_chain_refresh_monotonic = now_monotonic
        self._instrument_chain_stale_latched = False
        if recovered:
            self._log_event(
                "instrument_chain_recovered",
                {
                    "stale_age_s": round(
                        max(0.0, now_monotonic - prior_success),
                        3,
                    ),
                },
                now_ms=received_ms,
                include_connection=False,
            )

    def _require_authoritative_chain_fresh(
        self,
        received_ms: int,
        now_monotonic: float,
    ) -> None:
        last_success = self._last_authoritative_chain_refresh_monotonic
        if not last_success:
            return
        age = max(0.0, now_monotonic - last_success)
        if age <= CHAIN_MAX_STALE_SECONDS:
            return
        if not self._instrument_chain_stale_latched:
            self._log_event(
                "instrument_chain_stale",
                {
                    "age_s": round(age, 3),
                    "max_age_s": CHAIN_MAX_STALE_SECONDS,
                },
                now_ms=received_ms,
                include_connection=False,
            )
            self.stats["instrument_chain_stale"] += 1
            self._instrument_chain_stale_latched = True
        raise ConnectionError(
            "authoritative instrument chain stale for "
            f"{age:.1f}s"
        )

    def seed_instruments(self) -> None:
        """Synchronously seed every required currency before readiness."""
        try:
            fetched = self._fetch_instrument_snapshot_interruptibly()
            received_ms = int(time.time() * 1000)
            self._validate_instrument_snapshot_continuity(
                fetched,
                received_ms,
            )
            self._apply_instrument_snapshot(fetched, received_ms)
        except Exception as exc:
            failure_ms = int(time.time() * 1000)
            quarantined = (
                isinstance(exc, InstrumentSnapshotError)
                and exc.error_type
                in INSTRUMENT_SNAPSHOT_QUARANTINE_ERROR_TYPES
            )
            if quarantined:
                self.stats["instrument_snapshot_quarantined"] += 1
            self._record_instrument_snapshot_error(
                (
                    "instrument_snapshot_quarantined"
                    if quarantined
                    else "instrument_seed_error"
                ),
                exc,
                failure_ms,
            )
            raise
        self._mark_authoritative_chain_fresh(received_ms)

    def _chain_refresh_worker(
        self,
        generation: int,
        deadline: float | None,
    ) -> None:
        fetched = None
        error = None
        try:
            fetched = self._fetch_instrument_snapshot(deadline)
        except Exception as exc:
            error = exc
        with self._chain_refresh_lock:
            self._chain_refresh_result = (generation, fetched, error)

    def _start_chain_refresh(self) -> bool:
        with self._chain_refresh_lock:
            thread = self._chain_refresh_thread
            if thread is not None and thread.is_alive():
                return False
            if self._chain_refresh_result is not None:
                return False
            generation = self._chain_refresh_generation
            thread = threading.Thread(
                target=self._chain_refresh_worker,
                args=(generation, self._deadline_monotonic),
                name="deribit-chain-refresh",
                daemon=True,
            )
            self._chain_refresh_thread = thread
            thread.start()
            return True

    def _poll_chain_refresh(self, received_ms: int) -> None:
        with self._chain_refresh_lock:
            result = self._chain_refresh_result
            if result is None:
                return
            self._chain_refresh_result = None
            self._chain_refresh_thread = None
        generation, fetched, error = result
        if generation != self._chain_refresh_generation:
            return
        now = time.monotonic()
        if error is None:
            assert fetched is not None
            try:
                self._validate_instrument_snapshot_continuity(
                    fetched,
                    received_ms,
                )
            except InstrumentSnapshotError as exc:
                error = exc
        if error is not None:
            self.stats["instrument_refresh_errors"] = (
                self.stats.get("instrument_refresh_errors", 0) + 1
            )
            quarantined = (
                isinstance(error, InstrumentSnapshotError)
                and error.error_type
                in INSTRUMENT_SNAPSHOT_QUARANTINE_ERROR_TYPES
            )
            if quarantined:
                self.stats["instrument_snapshot_quarantined"] += 1
            self._record_instrument_snapshot_error(
                (
                    "instrument_snapshot_quarantined"
                    if quarantined
                    else "instrument_refresh_error"
                ),
                error,
                received_ms,
            )
            self._next_chain_refresh_attempt = (
                now + CHAIN_REFRESH_RETRY_SECONDS
            )
            return
        assert fetched is not None
        self._apply_instrument_snapshot(fetched, received_ms)
        self._mark_authoritative_chain_fresh(received_ms, now)
        if self.ws is not None and self.active_channels:
            self.reconcile_chain(received_ms)
        self.last_utc_day = utc_day(received_ms)
        self.last_chain_reconcile = now
        self._next_chain_refresh_attempt = now + CHAIN_REFRESH_SECONDS

    def _invalidate_chain_refresh(self) -> None:
        with self._chain_refresh_lock:
            self._chain_refresh_generation += 1
            self._chain_refresh_result = None
            thread = self._chain_refresh_thread
            if thread is not None and not thread.is_alive():
                self._chain_refresh_thread = None

    def _log_event(
        self,
        event: str,
        detail: dict | None = None,
        now_ms: int | None = None,
        include_connection: bool = True,
    ) -> None:
        now_ms = now_ms if now_ms is not None else int(time.time() * 1000)
        self._capture_event_seq += 1
        capture_event_seq = self._capture_event_seq
        rec = {
            "record_type": "capture_event",
            "dedup_id": sha256_id(
                f"event|{self.run_id}|{event}|{now_ms}|{self.msg_seq}|"
                f"{capture_event_seq}"
            ),
            "event": event,
            "received_ts_ms": now_ms,
            "capture_event_seq": capture_event_seq,
        }
        if include_connection and self._current_connection_id is not None:
            rec["connection_id"] = self._current_connection_id
        if detail:
            rec["detail"] = bounded_audit_detail(detail)
        try:
            self.sink.write("events", rec, now_ms=now_ms)
        except DayCapError:
            raise
        self.stats["events"] += 1

    def _attach_connection_provenance(self, record: dict) -> None:
        if self._current_connection_id is not None:
            record["connection_id"] = self._current_connection_id

    @staticmethod
    def _option_instrument_identity(
        instrument_name: object,
    ) -> OptionInstrumentIdentity | None:
        if not isinstance(instrument_name, str):
            return None
        parts = instrument_name.split("-")
        if len(parts) != 4:
            return None
        currency, expiry, strike_text, option_code = parts
        if currency.lower() not in CURRENCIES or currency != currency.upper():
            return None
        if option_code not in {"C", "P"} or len(expiry) not in {6, 7}:
            return None
        day_text = expiry[:-5]
        month_text = expiry[-5:-2]
        year_text = expiry[-2:]
        if (
            not day_text.isdigit()
            or not year_text.isdigit()
            or month_text not in OPTION_EXPIRY_MONTHS
        ):
            return None
        try:
            expiry_date = datetime(
                2000 + int(year_text),
                OPTION_EXPIRY_MONTHS[month_text],
                int(day_text),
                tzinfo=timezone.utc,
            )
            strike = Decimal(strike_text)
        except (
            InvalidOperation,
            OverflowError,
            ValueError,
        ):
            return None
        if not strike.is_finite() or strike <= 0:
            return None
        return OptionInstrumentIdentity(
            currency=currency,
            expiry=expiry_date,
            strike=strike,
            option_type=(
                "call" if option_code == "C" else "put"
            ),
        )

    @classmethod
    def _option_instrument_currency(
        cls,
        instrument_name: object,
    ) -> str | None:
        identity = cls._option_instrument_identity(
            instrument_name
        )
        return identity.currency if identity is not None else None

    @staticmethod
    def _option_trade_channel_currency(channel: object) -> str | None:
        prefix = "trades.option."
        suffix = ".100ms"
        if (
            not isinstance(channel, str)
            or not channel.startswith(prefix)
            or not channel.endswith(suffix)
        ):
            return None
        currency = channel[len(prefix):-len(suffix)]
        if currency.lower() not in CURRENCIES or currency != currency.upper():
            return None
        return currency

    def _option_trade_payload_error_field(
        self,
        trade: object,
        source_channel: str | None,
    ) -> str | None:
        """Return the first invalid canonical option-trade payload field."""
        if not isinstance(trade, dict):
            return "body"
        instrument_name = trade.get("instrument_name")
        instrument_currency = self._option_instrument_currency(
            instrument_name
        )
        if instrument_currency is None:
            return "instrument_name"
        if source_channel is not None:
            channel_currency = self._option_trade_channel_currency(
                source_channel
            )
            if channel_currency is None:
                return "source_channel"
            if channel_currency != instrument_currency:
                return "instrument_name"
        trade_id = trade.get("trade_id")
        if not isinstance(trade_id, str) or not trade_id:
            return "trade_id"
        for field in ("trade_seq", "timestamp"):
            value = trade.get(field)
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value <= 0
            ):
                return field
        if trade.get("direction") not in {"buy", "sell"}:
            return "direction"
        for field in ("amount", "price", "index_price"):
            value = trade.get(field)
            if not is_finite_number(value) or value <= 0:
                return field
        iv = trade.get("iv")
        if not is_finite_number(iv) or iv < 0:
            return "iv"
        return None

    def _option_trade_schema_error_field(self, trade: object) -> str | None:
        """Validate a direct trade call, including current metadata."""
        invalid_field = self._option_trade_payload_error_field(trade, None)
        if invalid_field is not None:
            return invalid_field
        assert isinstance(trade, dict)
        instrument = self.instruments.get(trade["instrument_name"])
        if (
            not isinstance(instrument, dict)
            or instrument.get("kind") != "option"
        ):
            return "instrument_name"
        return None
    def _option_ticker_schema_error_field(
        self,
        instrument_name: str,
        ticker: object,
    ) -> str | None:
        """Return the first invalid canonical option-ticker field."""
        if not isinstance(ticker, dict):
            return "body"
        if ticker.get("instrument_name") != instrument_name:
            return "instrument_name"
        timestamp = ticker.get("timestamp")
        if not is_finite_number(timestamp) or timestamp <= 0:
            return "timestamp"
        for field in ("index_price", "underlying_price"):
            value = ticker.get(field)
            if not is_finite_number(value) or value <= 0:
                return field
        for field in ("open_interest", "mark_iv"):
            value = ticker.get(field)
            if not is_finite_number(value) or value < 0:
                return field
        greeks = ticker.get("greeks")
        if not isinstance(greeks, dict):
            return "greeks"
        for field in ("delta", "gamma", "theta", "vega", "rho"):
            value = greeks.get(field)
            if not is_finite_number(value):
                return f"greeks.{field}"
        if greeks["gamma"] < 0:
            return "greeks.gamma"
        return None


    @staticmethod
    def _option_trade_dedup_id(trade: dict) -> str:
        return sha256_id(
            f"trade|{trade['instrument_name']}|{trade['trade_id']}|"
            f"{trade['trade_seq']}|{trade['timestamp']}"
        )

    @staticmethod
    def _option_ticker_dedup_id(
        instrument_name: str,
        timestamp: int | float,
    ) -> str:
        return sha256_id(f"ticker|{instrument_name}|{timestamp}")

    def ingest_trade(
        self,
        trade: dict,
        received_ms: int,
        source_channel: str | None = None,
    ) -> None:
        invalid_field = self._option_trade_payload_error_field(
            trade,
            source_channel,
        )
        if invalid_field is None and source_channel is None:
            invalid_field = self._option_trade_schema_error_field(trade)
        if invalid_field is not None:
            self.stats["trade_schema_errors"] += 1
            detail = {
                "field": invalid_field,
                "body_type": type(trade).__name__,
            }
            if source_channel is not None:
                detail["channel"] = source_channel
            self._log_event(
                "trade_schema_error",
                detail,
                now_ms=received_ms,
            )
            raise ConnectionError(
                f"option trade has invalid {invalid_field}"
            )
        iname = trade["instrument_name"]
        seq = trade["trade_seq"]
        dedup = self._option_trade_dedup_id(trade)
        if not self._dedupe_and_count("trades", dedup):
            return
        instrument = self.instruments.get(iname)
        metadata_resolved = (
            isinstance(instrument, dict)
            and instrument.get("kind") == "option"
        )
        source_ts_ms = trade["timestamp"]
        skew_ms = source_ts_ms - received_ms
        if abs(skew_ms) > REQUIRED_SOURCE_MAX_SKEW_MS:
            self.stats["source_timestamp_errors"] += 1
            self._log_event(
                "trade_source_timestamp_error",
                {
                    "instrument_name": iname,
                    "reason": (
                        "source_timestamp_future"
                        if skew_ms > 0
                        else "source_timestamp_stale"
                    ),
                    "source_ts_ms": source_ts_ms,
                    "received_ts_ms": received_ms,
                },
                now_ms=received_ms,
            )
            raise TradeRejectionError(
                f"option trade source timestamp skew for {iname}: "
                f"source={source_ts_ms} receipt={received_ms}"
            )
        last_seq = self._last_trade_seq.get(iname)
        gap_expected = (
            last_seq + 1
            if last_seq is not None and seq > last_seq + 1
            else None
        )
        if last_seq is not None and seq <= last_seq:
            self.stats["trade_sequence_conflicts"] += 1
            self._log_event(
                "trade_sequence_conflict",
                {
                    "instrument_name": iname,
                    "observed": seq,
                    "last_persisted": last_seq,
                    "reason": (
                        "repeat" if seq == last_seq else "regression"
                    ),
                },
                now_ms=received_ms,
            )
            raise TradeRejectionError(
                f"trade sequence conflict for {iname}: observed {seq} "
                f"at or below persisted {last_seq} with another trade id"
            )
        last_source_ts = self._last_trade_source_ts.get(iname)
        if last_source_ts is not None and source_ts_ms < last_source_ts:
            self.stats["trade_chronology_violations"] += 1
            self._log_event(
                "trade_source_timestamp_regression",
                {
                    "instrument_name": iname,
                    "observed": source_ts_ms,
                    "last_persisted": last_source_ts,
                    "sequence": seq,
                },
                now_ms=received_ms,
            )
            raise TradeRejectionError(
                f"trade source timestamp regression for {iname}: "
                f"sequence {seq} at {source_ts_ms} predates "
                f"persisted {last_source_ts}"
            )
        if gap_expected is not None:
            self.stats["trade_sequence_gaps"] += 1
            self._log_event(
                "trade_continuity_gap",
                {
                    "instrument_name": iname,
                    "expected": gap_expected,
                    "observed": seq,
                    "last_persisted": last_seq,
                },
                now_ms=received_ms,
            )
        rec = {
            "record_type": "option_trade",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "trade": trade,
        }
        if not metadata_resolved:
            rec["instrument_metadata_status"] = "unresolved_at_receipt"
            rec["source_channel"] = source_channel
        self._attach_connection_provenance(rec)
        written_bytes = self.sink.write(
            "trades",
            rec,
            now_ms=received_ms,
        )
        self._sequence_state_dirty = True
        self.deduper.commit(dedup)
        self._merge_trade_sequence(iname, seq, source_ts_ms)
        self.stats["trades_records"] += 1
        if not metadata_resolved:
            self.stats["trades_unresolved"] += 1
        self._record_restart_source_write(written_bytes)
        if not metadata_resolved:
            self._log_event(
                "trade_instrument_unresolved",
                {
                    "instrument_name": iname,
                    "source_channel": source_channel,
                    "dedup_id": dedup,
                },
                now_ms=received_ms,
            )
            raise UnresolvedInstrumentTradeError(
                f"authoritative metadata absent for option trade {iname}"
            )
        if gap_expected is not None:
            raise TradeContinuityError(
                f"trade sequence gap for {iname}: "
                f"expected {gap_expected}, observed {seq}"
            )

    def _validate_required_source_timestamp(
        self,
        channel: str,
        source_ts_ms: int | float,
        received_ms: int,
    ) -> tuple[bool, bool]:
        """Return (new sample usable, channel proven live) for one payload."""
        skew_ms = source_ts_ms - received_ms
        prior = self._required_last_source_timestamp_ms.get(channel)
        reason = None
        if skew_ms > REQUIRED_SOURCE_MAX_SKEW_MS:
            reason = "source_timestamp_future"
        elif skew_ms < -REQUIRED_SOURCE_MAX_SKEW_MS:
            reason = "source_timestamp_stale"
        if reason is not None:
            self.stats["source_timestamp_errors"] += 1
            self._log_event(
                "source_timestamp_integrity_error",
                {
                    "channel": channel,
                    "reason": reason,
                    "source_ts_ms": source_ts_ms,
                    "received_ts_ms": received_ms,
                    "prior_source_ts_ms": prior,
                },
                now_ms=received_ms,
            )
            raise SourceTimestampIntegrityError(
                f"{channel} {reason}: source={source_ts_ms} "
                f"receipt={received_ms} prior={prior}"
            )
        if prior is not None and source_ts_ms <= prior:
            repeated = source_ts_ms == prior
            self.stats[
                "source_timestamp_repeats"
                if repeated
                else "source_timestamp_regressions"
            ] += 1
            return False, repeated
        return True, True

    def _record_ticker_cap_drop(
        self,
        iname: str,
        channel: str,
        source_ts_ms: int | float,
        bucket: int,
    ) -> bool:
        self._required_last_source_timestamp_ms[channel] = source_ts_ms
        self._sequence_state_dirty = True
        if (
            bucket
            > self._ticker_cap_last_dropped_bucket.get(iname, -1)
        ):
            self._ticker_cap_last_dropped_bucket[iname] = bucket
            self.stats["ticker_records_dropped_cap"] += 1
        self._record_restart_source_write(0)
        return True

    def _handle_ticker_cap(
        self,
        exc: TickerCapError,
        iname: str,
        channel: str,
        source_ts_ms: int | float,
        received_ms: int,
        bucket: int,
    ) -> bool:
        day = utc_day(received_ms)
        if exc.scope == "cumulative":
            first_exhaustion = not self._ticker_total_cap_exhausted
            self._ticker_total_cap_exhausted = True
        else:
            first_exhaustion = self._ticker_cap_day != day
            self._ticker_cap_day = day
        self._record_ticker_cap_drop(
            iname,
            channel,
            source_ts_ms,
            bucket,
        )
        self._save_sequence_state()
        if first_exhaustion:
            self.stats["ticker_cap_exhaustions"] += 1
            self._log_event(
                "ticker_cap_exceeded",
                {
                    "scope": exc.scope,
                    "day_utc": day,
                    "reason": bounded_audit_text(exc),
                    "day_cap_bytes": self.sink.ticker_day_cap,
                    "total_cap_bytes": self.sink.ticker_total_cap,
                    "day_bytes": self.sink.ticker_day_bytes,
                    "prior_bytes": self.sink.ticker_prior_bytes,
                },
                now_ms=received_ms,
            )
        return True

    def ingest_ticker(
        self,
        iname: str,
        data: dict,
        received_ms: int,
    ) -> bool:
        ts = data["timestamp"]
        channel = f"ticker.{iname}.agg2"
        usable, live = self._validate_required_source_timestamp(
            channel,
            ts,
            received_ms,
        )
        if not usable:
            return live
        bucket = int(ts) // self.ticker_sample_ms
        if self._ticker_last_bucket.get(iname, -1) >= bucket:
            self._required_last_source_timestamp_ms[channel] = ts
            self._sequence_state_dirty = True
            self.stats["ticker_samples_suppressed"] += 1
            self._record_restart_source_write(0)
            return True
        target_day = utc_day(received_ms)
        if (
            self._ticker_total_cap_exhausted
            or self._ticker_cap_day == target_day
        ):
            return self._record_ticker_cap_drop(
                iname,
                channel,
                ts,
                bucket,
            )
        paused_day = self._ticker_cap_day
        compact = {
            field: data[field]
            for field in TICKER_FIELDS
            if field in data
        }
        dedup = self._option_ticker_dedup_id(iname, ts)
        rec = {
            "record_type": "option_ticker",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "sample_interval_ms": self.ticker_sample_ms,
            "instrument_name": iname,
            "ticker": compact,
        }
        self._attach_connection_provenance(rec)
        try:
            written_bytes = self.sink.write(
                "ticker",
                rec,
                now_ms=received_ms,
            )
        except TickerCapError as exc:
            return self._handle_ticker_cap(
                exc,
                iname,
                channel,
                ts,
                received_ms,
                bucket,
            )
        self._required_last_source_timestamp_ms[channel] = ts
        self._sequence_state_dirty = True
        self.deduper.commit(dedup)
        self._ticker_last_bucket[iname] = bucket
        self.stats["ticker_records"] += 1
        self._record_restart_source_write(written_bytes)
        if paused_day is not None and paused_day != target_day:
            self._ticker_cap_day = None
            self._ticker_cap_last_dropped_bucket.clear()
            self._sequence_state_dirty = True
            self._save_sequence_state()
            self._log_event(
                "ticker_cap_recovered",
                {
                    "prior_day_utc": paused_day,
                    "day_utc": target_day,
                },
                now_ms=received_ms,
            )
        return True

    def ingest_index(
        self,
        index_name: str,
        data: dict,
        received_ms: int,
    ) -> bool:
        ts = data.get("timestamp") if isinstance(data, dict) else None
        price = data.get("price") if isinstance(data, dict) else None
        payload_name = (
            data.get("index_name")
            if isinstance(data, dict)
            else None
        )
        if (
            not is_finite_number(ts)
            or ts <= 0
            or not is_finite_number(price)
            or price <= 0
            or payload_name not in (None, index_name)
        ):
            self.stats["index_schema_errors"] += 1
            self._log_event(
                "index_schema_error",
                {
                    "index_name": index_name,
                    "body_type": type(data).__name__,
                },
                now_ms=received_ms,
            )
            return False
        channel = f"deribit_price_index.{index_name}"
        usable, live = self._validate_required_source_timestamp(
            channel,
            ts,
            received_ms,
        )
        if not usable:
            return live
        dedup = sha256_id(f"index|{index_name}|{ts}")
        if not self._dedupe_and_count("index", dedup):
            return False
        rec = {
            "record_type": "index_price",
            "dedup_id": dedup,
            "received_ts_ms": received_ms,
            "index_name": index_name,
            "data": data,
        }
        self._attach_connection_provenance(rec)
        written_bytes = self.sink.write(
            "index",
            rec,
            now_ms=received_ms,
        )
        self._required_last_source_timestamp_ms[channel] = ts
        self._sequence_state_dirty = True
        self.deduper.commit(dedup)
        self.stats["index_records"] += 1
        self._record_restart_source_write(written_bytes)
        return True

    def _dedupe_and_count(self, kind: str, dedup: str) -> bool:
        if self.deduper.is_duplicate(dedup):
            if kind == "trades":
                self.stats["dup_trades"] += 1
            self._log_event(
                "duplicate_suppressed",
                {"dedup_id": dedup, "kind": kind},
            )
            return False
        return True

    def _instrument_fingerprint(self, row: dict) -> str:
        canonical = strict_json_dumps(
            row,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        return sha256_id(canonical)

    def _instrument_dedup_id(
        self,
        row: dict,
        state_ts: str | int,
    ) -> str:
        fingerprint = self._instrument_fingerprint(row)
        return sha256_id(
            f"instr|upsert|{row['instrument_name']}|"
            f"{fingerprint}|{state_ts}"
        )

    @staticmethod
    def _instrument_removal_dedup_id(
        instrument_name: str,
        prior_fingerprint: str,
        snapshot_id: str,
    ) -> str:
        return sha256_id(
            f"instr|removed|{instrument_name}|"
            f"{prior_fingerprint}|{snapshot_id}"
        )


    def reconcile_chain(self, received_ms: int) -> None:
        """Make the held ticker set exactly match the current open option chain."""
        if not self.capture_tickers:
            return
        want = set(self._ticker_channels())
        have = self._active_ticker_channels()
        add = want - have
        remove = have - want
        if remove:
            self._modify(sorted(remove), "public/unsubscribe")
            for channel in remove:
                self._required_last_received_monotonic.pop(channel, None)
            self._log_event("tickers_removed", {"n": len(remove)})
        if add:
            self._modify(sorted(add), "public/subscribe")
            if not self.await_required_coverage(add):
                raise ConnectionError(
                    f"ticker coverage incomplete after adding {len(add)} channels"
                )
            self._log_event("tickers_added", {"n": len(add)})

    def _connect_websocket_interruptibly(self) -> WsClient:
        result: dict[str, object] = {}
        done = threading.Event()
        connect_timeout = self._bounded_timeout(20.0)

        def connect_worker() -> None:
            try:
                client = WsClient(timeout=connect_timeout)
                if (
                    self.stop
                    or stop_requested()
                    or self._deadline_expired()
                ):
                    client.close()
                    result["error"] = ConnectionError(
                        "stop or deadline reached during websocket connect"
                    )
                else:
                    result["client"] = client
            except BaseException as exc:
                result["error"] = exc
            finally:
                done.set()

        threading.Thread(
            target=connect_worker,
            name="deribit-websocket-connect",
            daemon=True,
        ).start()
        while not done.wait(INTERRUPT_POLL_SECONDS):
            if self.stop or stop_requested():
                raise ConnectionError(
                    "stop requested during websocket connect"
                )
            if self._deadline_expired():
                raise TimeoutError(
                    "collector deadline reached during websocket connect"
                )
        if self.stop or stop_requested():
            client = result.get("client")
            if client is not None:
                client.close()
            raise ConnectionError(
                "stop requested after websocket connect"
            )
        error = result.get("error")
        if error is not None:
            raise error
        client = result.get("client")
        assert client is not None
        return client

    def run(self) -> int:
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)
        self._deadline_monotonic = (
            time.monotonic() + self.max_minutes * 60
            if self.max_minutes is not None
            else None
        )
        backoff = 1.0
        cap_terminal = False
        internal_terminal_type: str | None = None
        configuration_logged = False
        try:
            while (
                not self.stop
                and not stop_requested()
                and not self._deadline_expired()
            ):
                try:
                    if not configuration_logged:
                        self._log_event(
                            "collector_configuration",
                            {
                                "currencies": list(self.currencies),
                                "capture_tickers": self.capture_tickers,
                                "ticker_sample_ms": self.ticker_sample_ms,
                                "critical_day_cap_bytes": self.day_cap,
                                "critical_total_cap_bytes": self.total_cap,
                                "ticker_day_cap_bytes": self.ticker_day_cap,
                                "ticker_total_cap_bytes": self.ticker_total_cap,
                                "heartbeat": self.heartbeat,
                                "max_minutes": self.max_minutes,
                            },
                            include_connection=False,
                        )
                        configuration_logged = True
                    rehydrated = self.rehydrate_state()
                    self._log_event("restart_state_rehydrated", rehydrated)
                    self.seed_instruments()
                    break
                except DayCapError as exc:
                    self.stats["cap_caught"] = True
                    cap_terminal = True
                    print(f"DERIBIT_CAPTURE_CAP_EXCEEDED {exc}", flush=True)
                    self.stop = True
                except (
                    ConnectionError,
                    socket.timeout,
                    ssl.SSLError,
                    OSError,
                ) as exc:
                    if (
                        stop_requested()
                        or self._deadline_expired()
                        or self.sink._io_error is not None
                        or self._sequence_state_error is not None
                    ):
                        break
                    self.stats["reconnects"] += 1
                    try:
                        self._log_event(
                            "startup_retry_backoff",
                            {
                                "error": f"{type(exc).__name__}: {exc}",
                                "backoff_s": backoff,
                            },
                        )
                    except DayCapError:
                        self.stats["cap_caught"] = True
                        cap_terminal = True
                        print(
                            "DERIBIT_CAPTURE_CAP_EXCEEDED during "
                            "startup retry logging",
                            flush=True,
                        )
                        break
                    stop_aware_sleep(
                        min(backoff, max(0.0, self._remaining_seconds()))
                    )
                    backoff = min(backoff * 2.0, 60.0)
                except Exception as exc:
                    internal_terminal_type = type(exc).__name__[:128]
                    break

            while (
                not self.stop
                and not stop_requested()
                and not self._deadline_expired()
                and not cap_terminal
                and self.sink._io_error is None
                and self._sequence_state_error is None
                and internal_terminal_type is None
            ):
                self._connection_seq += 1
                self._current_connection_id = (
                    f"{self.run_id}:c{self._connection_seq}"
                )
                disconnect_reason = "session-loop-ended"
                disconnect_emitted = False
                connection_ready = False
                try:
                    self._invalidate_chain_refresh()
                    self._log_event(
                        "connection_attempt",
                        {"attempt": self._connection_seq},
                    )
                    self.ws = self._connect_websocket_interruptibly()
                    self._log_event(
                        "connection_connect",
                        {"attempt": self._connection_seq},
                    )
                    self.active_channels = set()
                    self._recently_unsubscribed = {}
                    self.seed_instruments()
                    ok = self.subscribe_all()
                    self._log_event(
                        "subscribe_ack",
                        {
                            "ok": ok,
                            "channels": sorted(self.active_channels),
                            "n_channels": len(self.active_channels),
                        },
                    )
                    if not ok:
                        raise ConnectionError(
                            "subscription ack not exact or heartbeat setup failed"
                        )
                    if not self.await_required_coverage():
                        raise ConnectionError(
                            "subscription acknowledged without complete "
                            "high-rate channel coverage"
                        )
                    if self._deadline_expired():
                        raise TimeoutError(
                            "collector deadline reached before readiness"
                        )
                    ready_ms = int(time.time() * 1000)
                    ready_monotonic = time.monotonic()
                    day_carried = self._ensure_instrument_day_fence(
                        ready_ms,
                        ready_monotonic,
                    )
                    self._stable_session = True
                    self._log_event(
                        "connection_ready",
                        {"attempt": self._connection_seq},
                        now_ms=ready_ms,
                    )
                    connection_ready = True
                    if not self._ready_printed:
                        print(READY_BANNER, flush=True)
                        self._ready_printed = True
                    backoff = 1.0
                    self.last_message_ms = ready_ms
                    self.last_message_monotonic = ready_monotonic
                    self.last_utc_day = utc_day(ready_ms)
                    self.last_chain_reconcile = ready_monotonic
                    if not day_carried:
                        self._next_chain_refresh_attempt = (
                            ready_monotonic + CHAIN_REFRESH_SECONDS
                        )
                    self._loop()
                    disconnect_reason = (
                        self._requested_stop_reason
                        or (
                            "deadline"
                            if self._deadline_expired()
                            else "session-loop-ended"
                        )
                    )
                except DayCapError as exc:
                    disconnect_reason = "cap"
                    self.stats["cap_caught"] = True
                    cap_terminal = True
                    print(
                        f"DERIBIT_CAPTURE_CAP_EXCEEDED {exc}",
                        flush=True,
                    )
                    break
                except (
                    ConnectionError,
                    socket.timeout,
                    ssl.SSLError,
                    OSError,
                ) as exc:
                    self._stable_session = False
                    disconnect_reason = f"{type(exc).__name__}: {exc}"
                    try:
                        self._log_event(
                            "connection_disconnect",
                            {
                                "reason": disconnect_reason,
                                "ready": connection_ready,
                            },
                        )
                        disconnect_emitted = True
                    except DayCapError:
                        self.stats["cap_caught"] = True
                        cap_terminal = True
                    except OSError:
                        pass
                    if self.ws is not None:
                        self.ws.close()
                        self.ws = None
                    if (
                        cap_terminal
                        or stop_requested()
                        or self._deadline_expired()
                        or self.sink._io_error is not None
                        or self._sequence_state_error is not None
                    ):
                        break
                    self.stats["reconnects"] += 1
                    try:
                        self._log_event(
                            "reconnect_backoff",
                            {
                                "error": disconnect_reason,
                                "backoff_s": backoff,
                            },
                        )
                    except DayCapError:
                        self.stats["cap_caught"] = True
                        cap_terminal = True
                        print(
                            "DERIBIT_CAPTURE_CAP_EXCEEDED during "
                            "reconnect logging",
                            flush=True,
                        )
                        break
                    stop_aware_sleep(
                        min(
                            backoff,
                            max(0.0, self._remaining_seconds()),
                        )
                    )
                    backoff = min(backoff * 2.0, 60.0)
                except Exception as exc:
                    internal_terminal_type = type(exc).__name__[:128]
                    disconnect_reason = (
                        f"internal-error:{internal_terminal_type}"
                    )
                    break
                finally:
                    self._invalidate_chain_refresh()
                    if (
                        not disconnect_emitted
                        and not cap_terminal
                        and self.sink._io_error is None
                    ):
                        try:
                            self._log_event(
                                "connection_disconnect",
                                {
                                    "reason": disconnect_reason,
                                    "ready": connection_ready,
                                },
                            )
                        except DayCapError:
                            self.stats["cap_caught"] = True
                            cap_terminal = True
                        except OSError:
                            pass
                        except Exception as exc:
                            if internal_terminal_type is None:
                                internal_terminal_type = (
                                    type(exc).__name__[:128]
                                )
                    if self.ws is not None:
                        self.ws.close()
                        self.ws = None
                    self._stable_session = False
                    self._current_connection_id = None
        except Exception as exc:
            if internal_terminal_type is None:
                internal_terminal_type = type(exc).__name__[:128]
        finally:
            try:
                self._emit_session_stop(
                    cap_terminal,
                    internal_terminal_type is not None,
                )
            except DayCapError:
                self.stats["cap_caught"] = True
                cap_terminal = True
                print(
                    "DERIBIT_CAPTURE_CAP_EXCEEDED during session stop",
                    flush=True,
                )
            except OSError:
                pass
            except Exception as exc:
                if internal_terminal_type is None:
                    internal_terminal_type = type(exc).__name__[:128]
            try:
                self._save_sequence_state()
            except OSError:
                pass
            try:
                self.sink.close()
            except OSError:
                pass
            print(SHUTDOWN_BANNER, flush=True)
        if self.sink._io_error is not None:
            print(
                f"DERIBIT_CAPTURE_IO_ERROR {self.sink._io_error}",
                file=sys.stderr,
                flush=True,
            )
        if self._sequence_state_error is not None:
            print(
                "DERIBIT_CAPTURE_SEQUENCE_STATE_ERROR "
                f"{self._sequence_state_error}",
                file=sys.stderr,
                flush=True,
            )
        if internal_terminal_type is not None:
            print(
                "DERIBIT_CAPTURE_INTERNAL_ERROR "
                f"{internal_terminal_type}",
                file=sys.stderr,
                flush=True,
            )
        if (
            self.sink._io_error is not None
            or self._sequence_state_error is not None
        ):
            return 5
        if cap_terminal:
            return 3
        if internal_terminal_type is not None:
            return 5
        return 0

    def _modify(self, channels: list[str], method: str) -> None:
        chunk_size = 220
        for i in range(0, len(channels), chunk_size):
            if stop_requested() or self.stop or self._deadline_expired():
                raise ConnectionError("stop or deadline during modify")
            chunk = channels[i:i + chunk_size]
            rid = self._next_id()
            assert self.ws is not None
            self.ws.send_text(json.dumps({
                "jsonrpc": "2.0",
                "id": rid,
                "method": method,
                "params": {"channels": chunk},
            }))
            if not self._await_ack(rid, set(chunk), timeout_s=15.0, method=method):
                raise ConnectionError(
                    f"{method} acknowledgement incomplete for {len(chunk)} channels"
                )

    def _check_connection_liveness(
        self,
        now_monotonic: float | None = None,
    ) -> None:
        if not self.last_message_monotonic:
            return
        now_monotonic = (
            time.monotonic()
            if now_monotonic is None
            else now_monotonic
        )
        idle_ms = int(
            max(0.0, now_monotonic - self.last_message_monotonic) * 1000
        )
        if idle_ms > int(SESSION_IDLE_TIMEOUT_SECONDS * 1000):
            self._log_event(
                "connection_stale",
                {"idle_ms": idle_ms, "limit_s": SESSION_IDLE_TIMEOUT_SECONDS},
            )
            raise ConnectionError(f"session silent for {idle_ms} ms")

    def _loop(self) -> None:
        idle_diagnostic_emitted = False
        while (
            not self.stop
            and not stop_requested()
            and not self._deadline_expired()
        ):
            received_ms = int(time.time() * 1000)
            received_monotonic = time.monotonic()
            self._maybe_refresh_chain(received_ms)
            self._check_connection_liveness(received_monotonic)
            self._check_required_continuity(received_monotonic)
            timeout = self._bounded_timeout(INTERRUPT_POLL_SECONDS)
            try:
                assert self.ws is not None
                op, data = self.ws.recv_message(timeout=timeout)
            except socket.timeout:
                if self._deadline_expired():
                    break
                received_ms = int(time.time() * 1000)
                received_monotonic = time.monotonic()
                self._check_connection_liveness(received_monotonic)
                self._check_required_continuity(received_monotonic)
                elapsed_seconds = max(
                    0.0,
                    received_monotonic - self.last_message_monotonic,
                )
                if (
                    not idle_diagnostic_emitted
                    and elapsed_seconds >= CONNECTION_IDLE_LOG_SECONDS
                ):
                    self._log_event(
                        "connection_idle_ping",
                        {
                            "idle_ms": min(
                                int(elapsed_seconds * 1000),
                                2**31,
                            )
                        },
                    )
                    idle_diagnostic_emitted = True
                continue
            except (ConnectionError, OSError) as exc:
                raise ConnectionError(f"recv: {exc}")
            idle_diagnostic_emitted = False
            if op == OP_CLOSE:
                raise ConnectionError("server closed")
            if op != OP_TEXT:
                continue
            received_ms = int(time.time() * 1000)
            obj = self._decode_jsonrpc_object(
                data,
                "steady_state",
                received_ms,
                audit_malformed=True,
            )
            if obj is None:
                continue
            self._handle_notification(
                obj,
                received_ms,
                phase="steady_state",
            )

    def _maybe_refresh_chain(self, now_ms: int) -> None:
        self._poll_chain_refresh(now_ms)
        if self.stop or stop_requested() or self._deadline_expired():
            return
        now = time.monotonic()
        self._ensure_instrument_day_fence(now_ms, now)
        self._require_authoritative_chain_fresh(now_ms, now)
        if now < self._next_chain_refresh_attempt:
            return
        if self._start_chain_refresh():
            self._next_chain_refresh_attempt = (
                now + CHAIN_REFRESH_SECONDS
            )

    def _route(
        self,
        channel: str,
        body,
        received_ms: int,
        received_monotonic: float | None = None,
    ) -> None:
        received_monotonic = (
            time.monotonic()
            if received_monotonic is None
            else received_monotonic
        )
        if channel.startswith("trades.option.") or channel.startswith("trades.option-"):
            self.stats["trades_pushes"] += 1
            if not isinstance(body, list):
                self.stats["trade_schema_errors"] += 1
                self._log_event(
                    "trade_schema_error",
                    {
                        "field": "batch",
                        "body_type": type(body).__name__,
                    },
                    now_ms=received_ms,
                )
                raise ConnectionError("option trade notification is not a list")
            valid_trades = []
            invalid_trades = []
            for trade in body:
                valid = (
                    self._option_trade_payload_error_field(
                        trade,
                        channel,
                    )
                    is None
                )
                (valid_trades if valid else invalid_trades).append(trade)
            valid_trades.sort(
                key=lambda trade: (
                    trade["instrument_name"],
                    trade["trade_seq"],
                )
            )
            batch_error = None
            for trade in valid_trades:
                try:
                    self.ingest_trade(
                        trade,
                        received_ms,
                        source_channel=channel,
                    )
                except (
                    TradeContinuityError,
                    TradeRejectionError,
                    UnresolvedInstrumentTradeError,
                ) as exc:
                    batch_error = batch_error or exc
            for trade in invalid_trades:
                try:
                    self.ingest_trade(
                        trade,
                        received_ms,
                        source_channel=channel,
                    )
                except ConnectionError as exc:
                    batch_error = batch_error or exc
            if batch_error is not None:
                raise batch_error
        elif channel.startswith("ticker.") and channel.endswith(".agg2"):
            self.stats["ticker_pushes"] += 1
            iname = channel[len("ticker."):-len(".agg2")]
            invalid_field = self._option_ticker_schema_error_field(
                iname,
                body,
            )
            if invalid_field is not None:
                self.stats["ticker_schema_errors"] += 1
                self._log_event(
                    "ticker_schema_error",
                    {
                        "field": invalid_field,
                        "channel": channel,
                        "body_type": type(body).__name__,
                    },
                    now_ms=received_ms,
                )
                return
            if self.ingest_ticker(iname, body, received_ms):
                self._required_last_received_monotonic[channel] = (
                    received_monotonic
                )
        elif channel.startswith("deribit_price_index."):
            idx_name = channel[len("deribit_price_index."):]
            items = body if isinstance(body, list) else [body]
            usable = False
            for item in items:
                usable = (
                    self.ingest_index(idx_name, item, received_ms)
                    or usable
                )
            if usable:
                self._required_last_received_monotonic[channel] = (
                    received_monotonic
                )


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Deribit public options capture collector")
    p.add_argument("--raw-dir", type=Path,
                   default=Path("quant-trading/data/raw/focused/live/deribit"))
    p.add_argument("--currencies", default="BTC,ETH")
    p.add_argument(
        "--max-day-bytes",
        type=int,
        default=DEFAULT_DAY_CAP_BYTES,
        help="critical-stream per-UTC-day cap (default 256 MiB)",
    )
    p.add_argument(
        "--max-total-bytes",
        type=int,
        default=DEFAULT_TOTAL_CAP_BYTES,
        help=(
            "critical-stream cumulative cap (default 10 GiB); "
            "prior capture never deleted"
        ),
    )
    p.add_argument(
        "--max-ticker-day-bytes",
        type=int,
        default=DEFAULT_TICKER_DAY_CAP_BYTES,
        help=(
            "isolated ticker per-UTC-day cap (default 256 MiB); "
            "exhaustion does not stop critical streams"
        ),
    )
    p.add_argument(
        "--max-ticker-total-bytes",
        type=int,
        default=DEFAULT_TICKER_TOTAL_CAP_BYTES,
        help=(
            "isolated ticker cumulative cap (default 10 GiB); "
            "exhaustion does not stop critical streams"
        ),
    )
    p.add_argument("--max-minutes", type=float, default=None,
                   help="bounded smoke duration; unset = persistent")
    p.add_argument(
        "--no-capture-tickers",
        action="store_false",
        dest="capture_tickers",
        help="explicit degraded mode: omit full-chain ticker/GEX inputs",
    )
    p.add_argument(
        "--ticker-sample-seconds",
        type=float,
        default=TICKER_SAMPLE_SECONDS,
        help="persist at most one compact ticker per instrument per interval",
    )
    p.add_argument("--no-heartbeat", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    currencies = tuple(c.strip().lower() for c in args.currencies.split(",") if c.strip())
    if not currencies:
        print("at least one currency is required", file=sys.stderr)
        return 2
    for c in currencies:
        if c not in ("btc", "eth"):
            print(f"unsupported currency {c} (research scope: btc, eth)", file=sys.stderr)
            return 2
    try:
        coll = Collector(
            args.raw_dir,
            currencies,
            day_cap=args.max_day_bytes,
            total_cap=args.max_total_bytes,
            heartbeat=not args.no_heartbeat,
            max_minutes=args.max_minutes,
            capture_tickers=args.capture_tickers,
            ticker_sample_seconds=args.ticker_sample_seconds,
            ticker_day_cap=args.max_ticker_day_bytes,
            ticker_total_cap=args.max_ticker_total_bytes,
        )
    except CollectorLeaseError as exc:
        print(f"DERIBIT_CAPTURE_LOCKED {exc}", file=sys.stderr)
        return 6
    return coll.run()


if __name__ == "__main__":
    raise SystemExit(main())


