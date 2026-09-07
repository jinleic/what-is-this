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

from collections import deque
from concurrent.futures import ThreadPoolExecutor
import argparse
import asyncio
import fcntl
import hashlib
import json
import os
import random
import signal
import socket
import stat
import sys
import time
import uuid
from decimal import Decimal, InvalidOperation
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
WS_PING_INTERVAL_S = 20
WS_PING_TIMEOUT_S = 20

FILE_PREFIX = "liquidation-capture-"
STATE_FILE = "collector-state.json"
OUTPUT_LEASE_FILE = ".collector.lock"
STATE_SAVE_INTERVAL_S = 5.0
STATE_SAVE_BYTE_INTERVAL = MI
DEDUP_MAX_IDS = 2_000_000
DEDUP_REHYDRATE_BLOCK_BYTES = 64 * 1024
DEDUP_REHYDRATE_MAX_BYTES = 32 * MI


def is_force_order_stream(stream: object) -> bool:
    if not isinstance(stream, str):
        return False
    return stream == "!forceOrder@arr" or (
        stream.endswith("@forceOrder")
        and len(stream) > len("@forceOrder")
    )
def stream_config_error(streams: list[str]) -> str | None:
    force_order_streams = [
        stream
        for stream in streams
        if is_force_order_stream(stream)
    ]
    if not force_order_streams:
        return (
            "streams must include !forceOrder@arr or a "
            "symbol-specific *@forceOrder stream"
        )
    if (
        "!forceOrder@arr" in force_order_streams
        and len(force_order_streams) > 1
    ):
        return (
            "!forceOrder@arr cannot overlap a symbol-specific "
            "*@forceOrder stream"
        )
    identities = [stream.lower() for stream in force_order_streams]
    if len(set(identities)) != len(identities):
        return (
            "force-order streams must have distinct case-insensitive "
            "identities"
        )
    if any(
        FORCE_ORDER_IDENTITY_DELIMITER in stream
        for stream in force_order_streams
    ):
        return (
            "force-order streams must not contain the identity "
            f"delimiter {FORCE_ORDER_IDENTITY_DELIMITER!r}"
        )
    return None


def bounded_error_provenance(error: object) -> dict:
    """Bounded, typed provenance for one untrusted server error payload."""
    if isinstance(error, dict):
        code = error.get("code")
        message = error.get("msg")
        return {
            "code": (
                code
                if isinstance(code, int) and not isinstance(code, bool)
                else None
            ),
            "msg": None if message is None else str(message)[:256],
        }
    return {"code": None, "msg": str(error)[:256]}


def utc_now() -> tuple[str, float]:
    now = datetime.now(timezone.utc)
    return now.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z", now.timestamp()


def utc_date() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d")


class ForceOrderSchemaError(ConnectionError):
    """A force-order frame cannot be assigned a safe semantic identity."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason

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


def open_regular_append_text(path: Path):
    """Append UTF-8 lines to a regular file, never a link or special."""
    return _regular_fdopen(
        os.open(
            path,
            _nofollow_flags(os.O_RDWR | os.O_CREAT | os.O_APPEND),
            0o600,
        ),
        path,
        "a",
        encoding="utf-8",
        buffering=1,
    )


def open_regular_read(path: Path):
    """Read a regular capture file, never a link or special."""
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


def regular_capture_size(path: Path) -> int:
    """Size of one regular capture file; fail closed on links/specials."""
    metadata = os.lstat(path)
    if not stat.S_ISREG(metadata.st_mode):
        raise NotARegularCaptureFileError(
            f"capture path is not a regular file: {path}"
        )
    return metadata.st_size


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




FORCE_ORDER_TEXT_FIELDS = ("s", "S", "o", "f", "X")
FORCE_ORDER_NUMERIC_FIELDS = ("q", "p", "ap", "l", "z")
FORCE_ORDER_IDENTITY_DELIMITER = "|"
FORCE_ORDER_MAX_SKEW_MS = 60_000
FORCE_ORDER_MAX_PUBLICATION_DELAY_MS = 60_000
AUDIT_TEXT_MAX_CHARS = 256


def bounded_audit_text(text: object) -> str:
    """Clamp one untrusted string before it reaches a durable audit row."""
    return str(text)[:AUDIT_TEXT_MAX_CHARS]


def bounded_audit_detail(detail, depth: int = 0):
    """Clamp every string inside one control-row field payload."""
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


def force_order_schema_error(payload: object) -> str | None:
    if not isinstance(payload, dict):
        return "frame_not_object"
    body = (
        payload.get("data")
        if isinstance(payload.get("data"), dict)
        else payload
    )
    if body.get("e") != "forceOrder":
        return "event_type"
    event_time = body.get("E")
    if (
        not isinstance(event_time, int)
        or isinstance(event_time, bool)
        or event_time <= 0
    ):
        return "event_timestamp"
    order = body.get("o")
    if not isinstance(order, dict) or not order:
        return "order_object"
    order_time = order.get("T")
    if (
        not isinstance(order_time, int)
        or isinstance(order_time, bool)
        or order_time <= 0
    ):
        return "order_timestamp"
    if order_time > event_time:
        return "order_timestamp_after_event"
    if event_time - order_time > FORCE_ORDER_MAX_PUBLICATION_DELAY_MS:
        return "order_publication_delay"
    for field in FORCE_ORDER_TEXT_FIELDS:
        value = order.get(field)
        if not isinstance(value, str) or not value:
            return f"order_{field}_text"
        if FORCE_ORDER_IDENTITY_DELIMITER in value:
            return f"order_{field}_delimiter"
    if order["S"] not in {"BUY", "SELL"}:
        return "order_side"
    numbers: dict[str, Decimal] = {}
    for field in FORCE_ORDER_NUMERIC_FIELDS:
        value = order.get(field)
        if not isinstance(value, str) or not value:
            return f"order_{field}_numeric_text"
        try:
            number = Decimal(value)
        except InvalidOperation:
            return f"order_{field}_numeric_text"
        if not number.is_finite() or number < 0:
            return f"order_{field}_numeric_value"
        numbers[field] = number
    if numbers["q"] <= 0:
        return "order_quantity_nonpositive"
    if numbers["l"] > numbers["z"]:
        return "order_last_fill_exceeds_accumulated"
    if numbers["z"] > numbers["q"]:
        return "order_accumulated_fill_exceeds_quantity"
    return None


def canonical_force_order(
    stream: str,
    payload: dict,
) -> tuple[str, dict] | None:
    """Return a validated canonical stream/body pair, or None."""
    if force_order_schema_error(payload) is not None:
        return None
    body = (
        payload.get("data")
        if isinstance(payload.get("data"), dict)
        else payload
    )
    canonical_stream = str(stream or payload.get("stream") or "")
    if not is_force_order_stream(canonical_stream):
        return None
    order = body["o"]
    if (
        canonical_stream != "!forceOrder@arr"
        and canonical_stream.removesuffix("@forceOrder").lower()
        != order["s"].lower()
    ):
        return None
    return canonical_stream, body


def event_id_for(stream: str, payload: dict) -> str:
    """Hash one validated exchange-semantic identity, never receipt time."""
    canon = canonical_force_order(stream, payload)
    body = (
        payload.get("data", payload)
        if isinstance(payload, dict)
        else {}
    )
    if canon is not None:
        canon_stream, body = canon
        order = body["o"]
        parts = [canon_stream, str(body["E"])]
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
        ):
            parts.append(str(order[key]))
        return hashlib.sha256(
            "|".join(parts).encode("utf-8")
        ).hexdigest()[:32]
    if isinstance(body, dict) and body.get("e") == "forceOrder":
        raise ForceOrderSchemaError(
            force_order_schema_error(payload)
            or "force_order_stream_identity"
        )
    stream = str(payload.get("stream") or stream or "?")
    if isinstance(body, dict) and body.get("e") == "aggTrade":
        parts = [
            stream,
            str(body.get("E", "")),
            str(body.get("s", "")),
            str(body.get("a", "")),
        ]
        return hashlib.sha256(
            "|".join(parts).encode("utf-8")
        ).hexdigest()[:32]
    if isinstance(body, dict) and body.get("e") == "markPriceUpdate":
        parts = [
            stream,
            str(body.get("E", "")),
            str(body.get("s", "")),
            str(body.get("T", "")),
        ]
        return hashlib.sha256(
            "|".join(parts).encode("utf-8")
        ).hexdigest()[:32]
    parts = [
        stream,
        strict_json_dumps(body, sort_keys=True, separators=(",", ":")),
    ]
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


class CapExceeded(Exception):
    pass

class SinkIOError(OSError):
    pass

class RestartReplayAborted(Exception):
    pass

class CollectorLeaseError(RuntimeError):
    """Raised when another process owns the output directory."""


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


OI_REST_BASE = "https://fapi.binance.com/futures/data/openInterestHist"
OI_SNAPSHOT_PREFIX = "oi-snapshot-"
OI_EGRESS_RETRY_S = 60 * 60
OI_TRANSIENT_RETRY_BASE_S = 60
OI_MAX_RETRY_S = 60 * 60
OI_EGRESS_HTTP_STATUSES = frozenset({403, 418, 451})
OI_GLOBAL_FAILURE_SCOPES = frozenset({
    "network",
    "egress",
    "client",
    "collector",
    "service",
})


def oi_failure_is_global(failure: object) -> bool:
    """Whether one OI failure plausibly affects every symbol."""
    if not isinstance(failure, dict):
        return True
    status = failure.get("http_status")
    global_status = (
        isinstance(status, int)
        and not isinstance(status, bool)
        and (
            status in OI_EGRESS_HTTP_STATUSES
            or status == 429
            or status >= 500
        )
    )
    return (
        failure.get("scope") in OI_GLOBAL_FAILURE_SCOPES
        or global_status
    )


def oi_failure_retry_seconds(
    failure: dict,
    failure_count: int,
) -> int:
    status = failure.get("http_status")
    if status in OI_EGRESS_HTTP_STATUSES:
        return OI_EGRESS_RETRY_S
    retry_s = min(
        OI_MAX_RETRY_S,
        OI_TRANSIENT_RETRY_BASE_S
        * (2 ** min(max(0, failure_count - 1), 6)),
    )
    if status == 429:
        retry_s = max(retry_s, 5 * 60)
    return retry_s


OI_SOURCE_MAX_AGE_MS = 15 * 60 * 1000
OI_SOURCE_MAX_FUTURE_SKEW_MS = 60 * 1000
OI_MAX_RESPONSE_BYTES = 256 * 1024
OI_READ_BLOCK_BYTES = 32 * 1024


def validate_open_interest_row(
    data: object,
    symbol: str,
    received_epoch_ms: int,
) -> tuple[dict | None, dict | None]:
    """Return one unchanged Binance OI row or a bounded integrity failure."""
    if (
        not isinstance(data, list)
        or len(data) != 1
        or not isinstance(data[0], dict)
    ):
        return None, {
            "reason": "schema-mismatch",
            "scope": "endpoint",
            "field": "response",
        }
    row = data[0]
    if row.get("symbol") != symbol:
        return None, {
            "reason": "schema-mismatch",
            "scope": "endpoint",
            "field": "symbol",
        }
    source_timestamp_ms = row.get("timestamp")
    if (
        not isinstance(source_timestamp_ms, int)
        or isinstance(source_timestamp_ms, bool)
        or source_timestamp_ms <= 0
    ):
        return None, {
            "reason": "schema-mismatch",
            "scope": "endpoint",
            "field": "timestamp",
        }
    for field in ("sumOpenInterest", "sumOpenInterestValue"):
        value = row.get(field)
        if not isinstance(value, str) or not value or value.strip() != value:
            return None, {
                "reason": "schema-mismatch",
                "scope": "endpoint",
                "field": field,
            }
        try:
            parsed = Decimal(value)
        except InvalidOperation:
            parsed = None
        if parsed is None or not parsed.is_finite() or parsed < 0:
            return None, {
                "reason": "schema-mismatch",
                "scope": "endpoint",
                "field": field,
            }
    age_ms = received_epoch_ms - source_timestamp_ms
    if age_ms > OI_SOURCE_MAX_AGE_MS:
        return None, {
            "reason": "source-timestamp-integrity",
            "scope": "endpoint",
            "violation": "stale",
        }
    if age_ms < -OI_SOURCE_MAX_FUTURE_SKEW_MS:
        return None, {
            "reason": "source-timestamp-integrity",
            "scope": "endpoint",
            "violation": "future",
        }
    return row, None



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



def fetch_open_interest_snapshot(
    symbol: str,
    timeout: float = 5.0,
) -> tuple[dict | None, dict | None]:
    """Return either a public OI row or structured fail-closed provenance."""
    import urllib.error
    import urllib.parse
    import urllib.request

    url = f"{OI_REST_BASE}?" + urllib.parse.urlencode(
        {"symbol": symbol, "period": "5m", "limit": 1}
    )
    req = urllib.request.Request(url, headers={"User-Agent": "research-capture/1.0"})
    deadline = time.monotonic() + max(0.0, timeout)
    try:
        with urllib.request.urlopen(
            req,
            timeout=max(0.001, deadline - time.monotonic()),
        ) as resp:
            if resp.status != 200:
                return None, {
                    "reason": "http-error",
                    "http_status": int(resp.status),
                    "scope": "endpoint",
                }
            body = bytearray()
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    return None, {"reason": "timeout", "scope": "network"}
                bound_response_read_deadline(resp, remaining)
                chunk = resp.read(OI_READ_BLOCK_BYTES)
                if not chunk:
                    break
                body += chunk
                if len(body) > OI_MAX_RESPONSE_BYTES:
                    return None, {
                        "reason": "response-too-large",
                        "scope": "endpoint",
                    }
            data = strict_json_loads(bytes(body))
    except urllib.error.HTTPError as exc:
        return None, {
            "reason": "http-error",
            "http_status": int(exc.code),
            "scope": (
                "egress"
                if exc.code in OI_EGRESS_HTTP_STATUSES
                else "endpoint"
            ),
        }
    except (TimeoutError, socket.timeout):
        return None, {"reason": "timeout", "scope": "network"}
    except urllib.error.URLError as exc:
        return None, {
            "reason": "network-error",
            "error_type": type(exc.reason).__name__,
            "scope": "network",
        }
    except DuplicateJsonKeyError:
        return None, {
            "reason": "invalid-json",
            "scope": "endpoint",
            "error_type": "duplicate-member",
        }
    except NonFiniteJsonConstantError:
        return None, {
            "reason": "invalid-json",
            "scope": "endpoint",
            "error_type": "non-finite-constant",
        }
    except (json.JSONDecodeError, UnicodeDecodeError):
        return None, {"reason": "invalid-json", "scope": "endpoint"}
    except Exception as exc:  # fail closed without persisting sensitive text
        return None, {
            "reason": "unexpected-error",
            "error_type": type(exc).__name__,
            "scope": "client",
        }
    retrieved_utc, received_epoch_s = utc_now()
    received_epoch_ms = int(received_epoch_s * 1000)
    row, failure = validate_open_interest_row(
        data,
        symbol,
        received_epoch_ms,
    )
    if failure is not None:
        return None, failure
    assert row is not None
    return {
        "symbol": symbol,
        "row": row,
        "retrieved_utc": retrieved_utc,
        "received_epoch_ms": received_epoch_ms,
        "source_timestamp_ms": row["timestamp"],
        "url": url,
    }, None


class Collector:
    def __init__(self, out_dir: Path, url: str, streams: list[str], max_seconds: float | None,
                 max_bytes_per_day: int, max_total_bytes: int):
        self.out_dir = out_dir
        self.url = url.rstrip("/")
        self.streams = list(streams)
        config_error = stream_config_error(self.streams)
        if config_error is not None:
            raise ValueError(config_error)
        self._lease = OutputLease(out_dir)
        self.max_seconds = max_seconds
        self.max_bytes_per_day = max_bytes_per_day
        self.max_total_bytes = max_total_bytes

        self.stop = asyncio.Event()
        self._requested_stop_reason: str | None = None
        self._run_deadline_monotonic: float | None = None
        self._cap_reason: str | None = None
        self._sink_io_error: str | None = None
        self.file = None
        self.file_date = ""
        self.total_bytes = 0      # bytes written by THIS process
        self.state_total_bytes = 0  # bytes persisted by all previous runs
        self.day_bytes: dict[str, int] = {}
        self.seen_ids: set[str] = set()
        self.dedup_max_ids = DEDUP_MAX_IDS
        self._seen_id_order: deque[str] = deque()
        self.conn_seq = 0
        self.run_id = uuid.uuid4().hex
        self._oi_queue: asyncio.Queue = asyncio.Queue(maxsize=64)
        self._oi_executor: ThreadPoolExecutor | None = None
        self._oi_last_attempt: dict[str, float] = {}
        self._oi_health = "healthy"
        self._oi_probe_pending = False
        self._oi_failure_count = 0
        self._oi_retry_after_monotonic = 0.0
        self._oi_symbol_failure_count: dict[str, int] = {}
        self._oi_symbol_retry_after_monotonic: dict[str, float] = {}
        self._state_save_due_monotonic = 0.0
        self._bytes_at_last_state_save = 0
        self._ignored_market_control_conn: str | None = None
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
            "oi_circuit_skipped": 0,
            "oi_circuit_opened": 0,
            "oi_circuit_recovered": 0,
            "oi_symbol_circuit_skipped": 0,
            "oi_symbol_quarantined": 0,
            "oi_symbol_recovered": 0,
            "state_saves": 0,
            "state_save_errors": 0,
            "ignored_market_frames": 0,
            "dedup_ids_evicted": 0,
            "dedup_rehydrate_bytes": 0,
            "dedup_rehydrate_ids": 0,
            "dedup_rehydrate_truncated": 0,
            "dedup_rehydrate_malformed": 0,
            "dedup_rehydrate_semantic_rejected": 0,
        }
    def _connection_id(self) -> str:
        return f"{self.run_id}:c{self.conn_seq}"

    def _qualified_connection_id(self, connection_id: str | None) -> str:
        prefix = f"{self.run_id}:c"
        if (
            isinstance(connection_id, str)
            and connection_id.startswith(prefix)
        ):
            return connection_id
        return self._connection_id()


    def release_lease(self) -> None:
        self._lease.close()


    # ---------- state (cumulative cap persistence) ----------
    def _capture_paths(self) -> list[Path]:
        try:
            root_mode = self.out_dir.stat().st_mode
        except FileNotFoundError:
            return []
        if not stat.S_ISDIR(root_mode):
            raise NotADirectoryError(str(self.out_dir))
        return list(self.out_dir.glob(f"{FILE_PREFIX}*.jsonl"))

    def actual_capture_bytes(self) -> int:
        """Total bytes on disk across all capture files (fail-closed floor)."""
        return sum(
            regular_capture_size(path)
            for path in self._capture_paths()
        )

    @staticmethod
    def _validate_state_payload(payload) -> tuple[int, dict[str, int]]:
        if not isinstance(payload, dict):
            raise ValueError("state root must be an object")
        if payload.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("unsupported state schema")
        total_bytes = payload.get("total_bytes")
        day_bytes = payload.get("day_bytes")
        if type(total_bytes) is not int or total_bytes < 0:
            raise ValueError("total_bytes must be a nonnegative integer")
        if not isinstance(day_bytes, dict):
            raise ValueError("day_bytes must be an object")
        validated_days = {}
        for day, byte_count in day_bytes.items():
            if (
                not isinstance(day, str)
                or len(day) != 8
                or not day.isdigit()
            ):
                raise ValueError("day_bytes key must be YYYYMMDD")
            datetime.strptime(day, "%Y%m%d")
            if type(byte_count) is not int or byte_count < 0:
                raise ValueError(
                    "day_bytes values must be nonnegative integers"
                )
            validated_days[day] = byte_count
        return total_bytes, validated_days

    def load_state(self) -> None:
        path = self.out_dir / STATE_FILE
        try:
            payload = strict_json_loads(
                path.read_text(encoding="utf-8")
            )
            state_total_bytes, state_day_bytes = (
                self._validate_state_payload(payload)
            )
            self.state_total_bytes = state_total_bytes
            for day, byte_count in state_day_bytes.items():
                self.day_bytes[day] = max(
                    byte_count,
                    int(self.day_bytes.get(day, 0)),
                )
        except FileNotFoundError:
            pass
        except (
            json.JSONDecodeError,
            OSError,
            OverflowError,
            TypeError,
            UnicodeError,
            ValueError,
        ):
            print(
                "collector-state.json unreadable; failing closed on caps",
                file=sys.stderr,
            )
            self.day_bytes = {}
            self.state_total_bytes = self.max_total_bytes + 1
        self.reconcile_day_bytes()
        self.state_total_bytes = max(
            self.state_total_bytes,
            self.actual_capture_bytes(),
        )

    def save_state(self) -> None:
        self.reconcile_day_bytes()
        path = self.out_dir / STATE_FILE
        tmp = path.with_suffix(".json.tmp")
        payload = {
            "schema_version": SCHEMA_VERSION,
            "total_bytes": max(
                self.state_total_bytes + self.total_bytes,
                self.actual_capture_bytes(),
            ),
            "day_bytes": dict(sorted(self.day_bytes.items())),
            "updated_utc": utc_now()[0],
        }
        with create_regular_temp(tmp) as fh:
            json.dump(
                payload,
                fh,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        self._fsync_output_directory("state-directory-fsync")
        self._bytes_at_last_state_save = self.total_bytes
        self._state_save_due_monotonic = time.monotonic() + STATE_SAVE_INTERVAL_S
        self.stats["state_saves"] += 1

    def maybe_save_state(self, force: bool = False) -> None:
        """Refresh the derived sidecar with bounded time and byte lag."""
        now = time.monotonic()
        byte_due = (
            self.total_bytes - self._bytes_at_last_state_save
            >= STATE_SAVE_BYTE_INTERVAL
        )
        if not force and now < self._state_save_due_monotonic and not byte_due:
            return
        try:
            self.save_state()
        except SinkIOError:
            self.stats["state_save_errors"] += 1
            raise
        except OSError as exc:
            self.stats["state_save_errors"] += 1
            self._state_save_due_monotonic = now + STATE_SAVE_INTERVAL_S
            print(f"STATE_SAVE_FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)

    # ---------- byte accounting / caps ----------
    # Append-only capture files are authoritative. `day_bytes` is a monotonic
    # in-memory cap ledger reconciled against those files; collector-state.json
    # is an atomically replaced derived cache. Either persisted or on-disk
    # accounting may only raise cap usage, never lower it.
    _file_day = ""

    def reconcile_day_bytes(self) -> None:
        """Merge real capture-file sizes into day_bytes (never decrease)."""
        for path in self._capture_paths():
            day = path.stem[len(FILE_PREFIX):]
            if not day.isdigit():
                continue
            size = regular_capture_size(path)
            self.day_bytes[day] = max(
                int(self.day_bytes.get(day, 0)),
                size,
            )

    def _day_used(self, day: str) -> int:
        return int(self.day_bytes.get(day, 0))

    def check_caps(self, incoming: int, day: str | None = None) -> None:
        """Latch the first cap breach; no later append may fit around it."""
        if self._cap_reason is not None:
            raise CapExceeded(self._cap_reason)
        target_day = day or self._file_day or utc_date()
        proj_day = self._day_used(target_day) + incoming
        if proj_day > self.max_bytes_per_day:
            self._cap_reason = (
                f"per-day cap: {proj_day} > {self.max_bytes_per_day} "
                f"bytes on {target_day}"
            )
            raise CapExceeded(self._cap_reason)
        proj_total = self.state_total_bytes + self.total_bytes + incoming
        if proj_total > self.max_total_bytes:
            self._cap_reason = (
                f"cumulative cap: {proj_total} > {self.max_total_bytes} bytes"
            )
            raise CapExceeded(self._cap_reason)

    def commit(self, raw_bytes: int, day: str) -> None:
        self.total_bytes += raw_bytes
        self.day_bytes[day] = int(self.day_bytes.get(day, 0)) + raw_bytes
    def _isolate_partial_tail(self, path: Path, day: str) -> bool:
        try:
            size = regular_capture_size(path)
        except FileNotFoundError:
            return False
        except OSError as exc:
            raise self._latch_sink_io("tail-inspect", exc) from exc
        if size == 0:
            return False
        try:
            with open_regular_read(path) as handle:
                handle.seek(-1, os.SEEK_END)
                final_byte = handle.read(1)
        except OSError as exc:
            raise self._latch_sink_io("tail-inspect", exc) from exc
        if final_byte == b"\n":
            return False
        self.day_bytes[day] = max(
            int(self.day_bytes.get(day, 0)),
            size,
        )
        authoritative_prior = max(
            0,
            self.actual_capture_bytes() - self.total_bytes,
        )
        self.state_total_bytes = max(
            self.state_total_bytes,
            authoritative_prior,
        )
        self.check_caps(1, day)
        try:
            with open_regular_append_text(path) as handle:
                written = handle.write("\n")
                if written != 1:
                    raise OSError(
                        f"short tail separator append: {written} of 1 byte"
                    )
                handle.flush()
                os.fsync(handle.fileno())
        except OSError as exc:
            raise self._latch_sink_io("tail-recovery", exc) from exc
        self.commit(1, day)
        return True


    _in_open_for = False

    def open_for(
        self,
        day: str,
        write_clock: tuple[str, float] | None = None,
    ) -> None:
        if self._in_open_for:
            return  # guard: control writes from inside open_for re-enter here
        self._in_open_for = True
        try:
            write_clock = write_clock or utc_now()
            prev_day = self._file_day if self.file is not None else None
            if self.file is not None:
                self.fsync_close()
            path = self.out_dir / f"{FILE_PREFIX}{day}.jsonl"
            try:
                os.lstat(path)
                existed = True
            except FileNotFoundError:
                existed = False
            except OSError as exc:
                raise self._latch_sink_io("open-stat", exc) from exc
            partial_tail_isolated = self._isolate_partial_tail(path, day)
            # Open the NEW day file BEFORE emitting any control so the rotate
            # record itself is written into the new file, never a closed handle.
            try:
                self.file = open_regular_append_text(path)
            except OSError as exc:
                raise self._latch_sink_io("open", exc) from exc
            self.file_date = day
            self._file_day = day
            if not existed:
                self.emit_control(
                    "file_created",
                    write_clock=write_clock,
                    path=path.name,
                )
                self._fsync_output_directory(
                    "capture-directory-fsync"
                )
            else:
                self.seed_dedup_from_file(path)
                # Reconcile this day's usage with the real size (never lower).
                try:
                    reopened_size = regular_capture_size(path)
                except OSError as exc:
                    raise self._latch_sink_io(
                        "reopen-stat",
                        exc,
                    ) from exc
                self.day_bytes[day] = max(
                    int(self.day_bytes.get(day, 0)),
                    reopened_size,
                )
                self.emit_control(
                    "file_reopened",
                    write_clock=write_clock,
                    path=path.name,
                    byte_count=reopened_size,
                    partial_tail_isolated=partial_tail_isolated,
                )
            if prev_day is not None and prev_day != day:
                self.emit_control(
                    "rotate",
                    write_clock=write_clock,
                    reason="utc capture day change",
                    from_day=prev_day,
                    to_day=day,
                )
        finally:
            self._in_open_for = False

    def _remember_event_id(self, event_id: str) -> bool:
        if event_id in self.seen_ids:
            return False
        if self.dedup_max_ids <= 0:
            raise ValueError("dedup_max_ids must be positive")
        while len(self._seen_id_order) >= self.dedup_max_ids:
            evicted = self._seen_id_order.popleft()
            self.seen_ids.discard(evicted)
            self.stats["dedup_ids_evicted"] += 1
        self.seen_ids.add(event_id)
        self._seen_id_order.append(event_id)
        return True

    def _check_dedup_rehydrate_abort(self) -> None:
        if self.stop.is_set():
            raise RestartReplayAborted(
                "stop requested during dedup restart replay"
            )
        if (
            self._run_deadline_monotonic is not None
            and time.monotonic() >= self._run_deadline_monotonic
        ):
            self.stop.set()
            raise RestartReplayAborted(
                "deadline reached during dedup restart replay"
            )

    def seed_dedup_from_file(self, path: Path) -> None:
        self._check_dedup_rehydrate_abort()
        newest_ids: list[str] = []
        found_ids: set[str] = set()

        def collect(raw_line: bytes) -> None:
            if not raw_line or len(newest_ids) >= self.dedup_max_ids:
                return
            try:
                rec = strict_json_loads(raw_line)
            except (
                UnicodeError,
                json.JSONDecodeError,
                NonFiniteJsonConstantError,
                DuplicateJsonKeyError,
            ):
                self.stats["dedup_rehydrate_malformed"] += 1
                return
            if not isinstance(rec, dict) or rec.get("kind") != "event":
                return
            event_id = rec.get("event_id")
            stream = rec.get("stream")
            payload = rec.get("payload")
            valid_stream = (
                isinstance(stream, str)
                and stream in self.streams
                and is_force_order_stream(stream)
            )
            canonical = (
                canonical_force_order(stream, payload)
                if valid_stream and isinstance(payload, dict)
                else None
            )
            expected_event_id = (
                event_id_for(stream, payload)
                if canonical is not None
                else None
            )
            received_ms = rec.get("recv_epoch_ms")
            trusted_receipt = (
                canonical is not None
                and isinstance(received_ms, int)
                and not isinstance(received_ms, bool)
                and received_ms > 0
                and abs(canonical[1]["E"] - received_ms)
                <= FORCE_ORDER_MAX_SKEW_MS
            )
            if (
                rec.get("v") != SCHEMA_VERSION
                or not isinstance(event_id, str)
                or not event_id
                or expected_event_id != event_id
                or not trusted_receipt
            ):
                self.stats[
                    "dedup_rehydrate_semantic_rejected"
                ] += 1
                return
            if (
                event_id not in found_ids
                and event_id not in self.seen_ids
            ):
                found_ids.add(event_id)
                newest_ids.append(event_id)

        try:
            with open_regular_read(path) as fh:
                fh.seek(0, os.SEEK_END)
                position = fh.tell()
                carry = b""
                bytes_read = 0
                while (
                    position > 0
                    and bytes_read < DEDUP_REHYDRATE_MAX_BYTES
                    and len(newest_ids) < self.dedup_max_ids
                ):
                    self._check_dedup_rehydrate_abort()
                    size = min(
                        DEDUP_REHYDRATE_BLOCK_BYTES,
                        position,
                        DEDUP_REHYDRATE_MAX_BYTES - bytes_read,
                    )
                    position -= size
                    fh.seek(position)
                    block = fh.read(size)
                    bytes_read += len(block)
                    self.stats["dedup_rehydrate_bytes"] += len(block)
                    pieces = (block + carry).split(b"\n")
                    carry = pieces[0]
                    for raw_line in reversed(pieces[1:]):
                        collect(raw_line)
                        if len(newest_ids) >= self.dedup_max_ids:
                            break
                if (
                    position == 0
                    and carry
                    and len(newest_ids) < self.dedup_max_ids
                ):
                    collect(carry)
                if position > 0:
                    self.stats["dedup_rehydrate_truncated"] += 1
        except OSError as exc:
            raise self._latch_sink_io("dedup-rehydrate", exc) from exc
        for event_id in reversed(newest_ids):
            self._remember_event_id(event_id)
        self.stats["dedup_rehydrate_ids"] += len(newest_ids)

    def _latch_sink_io(self, operation: str, exc: OSError) -> SinkIOError:
        if self._sink_io_error is None:
            self._sink_io_error = f"{operation}:{type(exc).__name__}"
        return SinkIOError(f"sink I/O terminal: {self._sink_io_error}")
    def _fsync_output_directory(self, operation: str) -> None:
        try:
            fsync_directory(self.out_dir)
        except OSError as exc:
            raise self._latch_sink_io(operation, exc) from exc


    def fsync_close(self) -> None:
        if self.file is None:
            return
        error = None
        try:
            self.file.flush()
            os.fsync(self.file.fileno())
        except OSError as exc:
            error = exc
        try:
            self.file.close()
        except OSError as exc:
            error = error or exc
        finally:
            self.file = None
        if error is not None:
            raise self._latch_sink_io("close", error) from error

    def write_line(
        self,
        record: dict,
        write_clock: tuple[str, float] | None = None,
    ) -> None:
        if self._cap_reason is not None:
            raise CapExceeded(self._cap_reason)
        if self._sink_io_error is not None:
            raise SinkIOError(f"sink I/O terminal: {self._sink_io_error}")
        write_clock = write_clock or utc_now()
        record = {**record, "run_id": self.run_id}
        if "ts_utc" not in record:
            record["ts_utc"] = write_clock[0]
        line = strict_json_dumps(
            record,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        wall_day = write_clock[0][:10].replace("-", "")
        rotated = self.file is None or self._file_day != wall_day
        if rotated:
            try:
                self.open_for(wall_day, write_clock=write_clock)
            except SinkIOError:
                raise
            except OSError as exc:
                raise self._latch_sink_io("open", exc) from exc
        nbytes = len(line.encode("utf-8")) + 1
        self.check_caps(nbytes, self._file_day)
        assert self.file is not None
        try:
            written = self.file.write(line + "\n")
            if written != len(line) + 1:
                raise OSError(
                    f"short append: {written} of {len(line) + 1} characters"
                )
            self.file.flush()
            os.fsync(self.file.fileno())
        except OSError as exc:
            raise self._latch_sink_io("append", exc) from exc
        self.commit(nbytes, self._file_day)
        self.maybe_save_state(force=rotated)
        if record.get("kind") == "event":
            self.stats["events_written"] += 1

    def emit_event(
        self,
        stream: str,
        raw_payload: dict,
        received_ms: int | None = None,
        connection_id: str | None = None,
        subscription_phase: str = "confirmed",
    ) -> None:
        if subscription_phase not in {"pre_ack", "confirmed"}:
            raise ValueError(
                f"invalid subscription phase: {subscription_phase}"
            )
        received_ms = (
            int(time.time() * 1000) if received_ms is None else received_ms
        )
        eid = event_id_for(stream, raw_payload)
        if eid in self.seen_ids:
            self.stats["events_deduped"] += 1
            return
        self.write_line({
            "v": SCHEMA_VERSION,
            "kind": "event",
            "recv_epoch_ms": received_ms,
            "conn": self._qualified_connection_id(connection_id),
            "subscription_phase": subscription_phase,
            "stream": stream,
            "event_id": eid,
            "source": {
                "venue": "binance-usds-m-futures",
                "endpoint_kind": "public-unauthenticated-websocket",
                "url": self.url,
            },
            "payload": raw_payload,
        })
        self._remember_event_id(eid)

    def emit_control(
        self,
        ctype: str,
        *,
        write_clock: tuple[str, float] | None = None,
        **fields,
    ) -> None:
        fields = bounded_audit_detail(dict(fields))
        if "conn" in fields:
            fields["conn"] = self._qualified_connection_id(fields["conn"])
        self.write_line(
            {
                "v": SCHEMA_VERSION,
                "kind": "control",
                "conn": self._connection_id(),
                "control": {"type": ctype, **fields},
            },
            write_clock=write_clock,
        )
        self.stats["controls_written"] += 1

    # ---------- websocket session ----------
    async def subscribe_and(
        self,
        ws,
        connection_id: str,
        ack_timeout: float = 30.0,
    ) -> None:
        """Await the exact ack while durably handling interleaved frames."""
        sub = {"method": "SUBSCRIBE", "id": 1, "params": self.streams}
        await ws.send(json.dumps(sub))
        deadline = time.monotonic() + ack_timeout
        while not self.stop.is_set():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(
                    f"subscription ack not received within {ack_timeout}s"
                )
            try:
                frame = await asyncio.wait_for(
                    ws.recv(),
                    timeout=min(remaining, 2.0),
                )
            except asyncio.TimeoutError:
                continue
            received_ms = int(time.time() * 1000)
            try:
                msg = strict_json_loads(frame)
            except UnicodeDecodeError:
                self.emit_control(
                    "invalid_utf8_frame",
                    conn=connection_id,
                    reason="invalid_utf8",
                )
                raise ForceOrderSchemaError("invalid_utf8")
            except NonFiniteJsonConstantError:
                self.emit_control(
                    "non_finite_json_frame",
                    conn=connection_id,
                    reason="non_finite_json_constant",
                )
                raise ForceOrderSchemaError(
                    "non_finite_json_constant"
                )
            except DuplicateJsonKeyError:
                self.emit_control(
                    "duplicate_json_frame",
                    conn=connection_id,
                    reason="duplicate_json_member",
                )
                raise ForceOrderSchemaError("duplicate_json_member")
            except json.JSONDecodeError:
                continue
            if not isinstance(msg, dict):
                continue
            if "error" in msg:
                raise RuntimeError(
                    "subscription error response: "
                    f"{bounded_error_provenance(msg['error'])}"
                )
            if (
                set(msg) == {"id", "result"}
                and type(msg["id"]) is int
                and msg["id"] == 1
                and msg["result"] is None
            ):
                self.stats["subscription_acks"] += 1
                return
            if "id" in msg or "result" in msg:
                raise RuntimeError(
                    "non-authoritative subscription ack: "
                    f"id={bounded_audit_text(msg.get('id'))} "
                    f"result={bounded_audit_text(msg.get('result'))}"
                )
            self.dispatch_market_message(
                msg,
                connection_id,
                received_ms,
                subscription_phase="pre_ack",
            )
        raise RuntimeError("stopped while awaiting subscription ack")

    def _configured_force_order_stream(
        self,
        msg: dict,
        body: dict,
    ) -> str:
        configured = {
            stream
            for stream in self.streams
            if is_force_order_stream(stream)
        }
        wrapped_stream = msg.get("stream")
        if wrapped_stream is not None:
            if not isinstance(wrapped_stream, str):
                raise ForceOrderSchemaError("force_order_stream_type")
            candidates = {
                stream
                for stream in configured
                if stream.lower() == wrapped_stream.lower()
            }
            if len(candidates) != 1:
                raise ForceOrderSchemaError(
                    "unconfigured_force_order_stream"
                )
            return next(iter(candidates))
        aggregate = {
            stream
            for stream in configured
            if stream.lower() == "!forceorder@arr"
        }
        if aggregate:
            if len(aggregate) != 1:
                raise ForceOrderSchemaError(
                    "ambiguous_force_order_stream"
                )
            return next(iter(aggregate))
        symbol_stream = f"{body['o']['s'].lower()}@forceorder"
        candidates = {
            stream
            for stream in configured
            if stream.lower() == symbol_stream
        }
        if len(candidates) != 1:
            raise ForceOrderSchemaError(
                "unconfigured_force_order_stream"
            )
        return next(iter(candidates))

    def _reject_force_order(
        self,
        reason: str,
        msg: dict,
        body: dict,
        conn_id: str,
    ) -> None:
        self.emit_control(
            "force_order_schema_error",
            conn=conn_id,
            reason=reason,
            stream=str(msg.get("stream") or "")[:128],
            event_type=str(body.get("e") or "")[:64],
        )
        raise ForceOrderSchemaError(reason)

    def _ignore_non_force_order_frame(
        self,
        msg: dict,
        body: dict,
        conn_id: str,
    ) -> None:
        self.stats["ignored_market_frames"] += 1
        if self._ignored_market_control_conn == conn_id:
            return
        self._ignored_market_control_conn = conn_id
        self.emit_control(
            "ignored_non_force_order_stream",
            conn=conn_id,
            stream=str(msg.get("stream") or ""),
            event_type=str(body.get("e") or ""),
        )

    def dispatch_market_message(
        self,
        msg,
        conn_id: str,
        received_ms: int,
        subscription_phase: str = "confirmed",
    ) -> None:
        """Persist one confirmed or pre-ack public market frame."""
        if not isinstance(msg, dict):
            self.emit_control(
                "invalid_frame_shape",
                conn=conn_id,
                frame_type=type(msg).__name__,
            )
            return
        if "error" in msg:
            provenance = bounded_error_provenance(msg["error"])
            self.emit_control(
                "error_frame",
                conn=conn_id,
                phase=subscription_phase,
                **provenance,
            )
            raise RuntimeError(
                "server error frame: "
                f"code={provenance['code']} msg={provenance['msg']}"
            )
        if msg.get("id") == 1:
            return
        body = (
            msg.get("data")
            if isinstance(msg.get("data"), dict)
            else msg
        )
        if body.get("e") != "forceOrder":
            if body.get("e") or msg.get("stream"):
                self._ignore_non_force_order_frame(msg, body, conn_id)
                return
            self.emit_control("unwrapped_frame", conn=conn_id)
            return
        schema_error = force_order_schema_error(msg)
        if schema_error is not None:
            self._reject_force_order(
                schema_error,
                msg,
                body,
                conn_id,
            )
        event_skew_ms = body["E"] - received_ms
        if abs(event_skew_ms) > FORCE_ORDER_MAX_SKEW_MS:
            self._reject_force_order(
                "source_timestamp_future"
                if event_skew_ms > 0
                else "source_timestamp_stale",
                msg,
                body,
                conn_id,
            )
        try:
            canonical_stream = self._configured_force_order_stream(
                msg,
                body,
            )
        except ForceOrderSchemaError as exc:
            self._reject_force_order(
                exc.reason,
                msg,
                body,
                conn_id,
            )
        canon = canonical_force_order(canonical_stream, msg)
        if canon is None:
            self._reject_force_order(
                "force_order_stream_identity",
                msg,
                body,
                conn_id,
            )
        canon_stream, canon_body = canon
        self.emit_event(
            canon_stream,
            msg,
            received_ms=received_ms,
            connection_id=conn_id,
            subscription_phase=subscription_phase,
        )
        self.schedule_oi_snapshot(canon_body["o"]["s"])

    def request_stop(self, *args) -> None:
        reason = "manual"
        if args:
            try:
                reason = signal.Signals(args[0]).name.lower()
            except (TypeError, ValueError):
                reason = "external"
        self._requested_stop_reason = reason
        self.stop.set()
    # ---------- bounded, non-blocking OI enrichment ----------
    OI_MIN_INTERVAL_S = 60.0      # rate limit per symbol (public REST courtesy)
    OI_QUEUE_MAX = 64             # bounded queue: stream receipt is never delayed

    def schedule_oi_snapshot(self, symbol: str) -> None:
        """Queue OI without spreading one symbol's quarantine."""
        if not symbol:
            return
        now = time.monotonic()
        last = self._oi_last_attempt.get(symbol)
        if last is not None and (now - last) < self.OI_MIN_INTERVAL_S:
            self.stats["oi_rate_limited"] += 1
            return
        if (
            now < self._oi_retry_after_monotonic
            or (
                self._oi_health == "blocked"
                and self._oi_probe_pending
            )
        ):
            self.stats["oi_circuit_skipped"] += 1
            return
        if now < self._oi_symbol_retry_after_monotonic.get(symbol, 0.0):
            self.stats["oi_circuit_skipped"] += 1
            self.stats["oi_symbol_circuit_skipped"] += 1
            return
        if self._oi_queue.qsize() >= self.OI_QUEUE_MAX:
            self.stats["oi_dropped_full_queue"] += 1
            return
        try:
            self._oi_queue.put_nowait(symbol)
        except asyncio.QueueFull:
            self.stats["oi_dropped_full_queue"] += 1
            return
        self._oi_last_attempt[symbol] = now
        if self._oi_health == "blocked":
            self._oi_probe_pending = True

    def _oi_fetch_executor(self) -> ThreadPoolExecutor:
        """Own the OI fetch thread so shutdown can quiesce it."""
        if self._oi_executor is None:
            self._oi_executor = ThreadPoolExecutor(
                max_workers=1,
                thread_name_prefix="oi-fetch",
            )
        return self._oi_executor

    def _quiesce_oi_executor(self) -> None:
        """Join the deadline-bounded in-flight fetch before returning."""
        executor = self._oi_executor
        if executor is None:
            return
        self._oi_executor = None
        executor.shutdown(wait=True, cancel_futures=True)

    def _oi_fetch_sync(self, symbol: str) -> tuple[dict | None, dict | None]:
        """Executor target with a total deadline bounded by the run deadline."""
        timeout = 5.0
        if self._run_deadline_monotonic is not None:
            remaining = self._run_deadline_monotonic - time.monotonic()
            if remaining <= 0:
                return None, {
                    "reason": "collector-deadline",
                    "scope": "collector",
                }
            timeout = min(timeout, remaining)
        return fetch_open_interest_snapshot(symbol, timeout)

    async def oi_worker(self) -> None:
        """Persist real OI or one structured probe failure per retry window."""
        while True:
            if self.stop.is_set() and self._oi_queue.empty():
                return
            get_task = asyncio.ensure_future(self._oi_queue.get())
            stop_task = asyncio.ensure_future(self.stop.wait())
            child_tasks = {get_task, stop_task}
            try:
                done, pending = await asyncio.wait(
                    child_tasks,
                    return_when=asyncio.FIRST_COMPLETED,
                )
            except BaseException:
                for task in child_tasks:
                    if not task.done():
                        task.cancel()
                await asyncio.gather(
                    *child_tasks,
                    return_exceptions=True,
                )
                raise
            for task in pending:
                task.cancel()
            if pending:
                await asyncio.gather(
                    *pending,
                    return_exceptions=True,
                )
            if stop_task in done and not (
                get_task in done and not get_task.exception()
            ):
                get_task.cancel()
                return
            symbol = get_task.result()
            now = time.monotonic()
            if now < self._oi_retry_after_monotonic:
                self.stats["oi_circuit_skipped"] += 1
                self._oi_last_attempt.pop(symbol, None)
                continue
            if now < self._oi_symbol_retry_after_monotonic.get(
                symbol,
                0.0,
            ):
                self.stats["oi_circuit_skipped"] += 1
                self.stats["oi_symbol_circuit_skipped"] += 1
                self._oi_last_attempt.pop(symbol, None)
                continue
            global_probe = self._oi_health == "blocked"
            if (
                self._run_deadline_monotonic is not None
                and time.monotonic() >= self._run_deadline_monotonic
            ):
                await self.stop.wait()
                return
            snapshot, failure = await asyncio.get_event_loop().run_in_executor(
                self._oi_fetch_executor(),
                self._oi_fetch_sync,
                symbol,
            )
            if (
                self._run_deadline_monotonic is not None
                and time.monotonic() >= self._run_deadline_monotonic
            ):
                await self.stop.wait()
                return
            if global_probe:
                self._oi_probe_pending = False
            record = {
                "v": SCHEMA_VERSION,
                "kind": "oi_snapshot",
                "symbol": symbol,
            }
            if snapshot is not None:
                symbol_was_quarantined = (
                    symbol in self._oi_symbol_failure_count
                    or symbol in self._oi_symbol_retry_after_monotonic
                )
                self._oi_symbol_failure_count.pop(symbol, None)
                self._oi_symbol_retry_after_monotonic.pop(symbol, None)
                if symbol_was_quarantined:
                    self.stats["oi_symbol_recovered"] += 1
                if global_probe:
                    self._oi_health = "healthy"
                    self._oi_failure_count = 0
                    self._oi_retry_after_monotonic = 0.0
                    self.stats["oi_circuit_recovered"] += 1
                record["status"] = "ok"
                record["data"] = snapshot
                self.stats["oi_ok"] += 1
            else:
                failure = failure or {
                    "reason": "missing-result",
                    "scope": "client",
                }
                if oi_failure_is_global(failure):
                    self._oi_health = "blocked"
                    self._oi_failure_count += 1
                    retry_s = oi_failure_retry_seconds(
                        failure,
                        self._oi_failure_count,
                    )
                    self._oi_retry_after_monotonic = (
                        time.monotonic() + retry_s
                    )
                    self.stats["oi_circuit_opened"] += 1
                else:
                    failure_count = (
                        self._oi_symbol_failure_count.get(symbol, 0)
                        + 1
                    )
                    self._oi_symbol_failure_count[symbol] = failure_count
                    retry_s = oi_failure_retry_seconds(
                        failure,
                        failure_count,
                    )
                    self._oi_symbol_retry_after_monotonic[symbol] = (
                        time.monotonic() + retry_s
                    )
                    self.stats["oi_symbol_quarantined"] += 1
                record["status"] = "failed-closed"
                record["failure"] = failure
                record["retry_after_seconds"] = retry_s
                self.stats["oi_failed"] += 1
            self.write_line(record)

    async def run(self) -> int:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.load_state()
        if self.state_total_bytes >= self.max_total_bytes:
            print("total cap already reached; failing closed", file=sys.stderr)
            self.release_lease()
            return 4
        run_clock = utc_now()
        day = run_clock[0][:10].replace("-", "")
        if self._day_used(day) >= self.max_bytes_per_day:
            print("per-day cap already reached; failing closed", file=sys.stderr)
            self.release_lease()
            return 4

        deadline = (
            time.monotonic() + self.max_seconds
            if self.max_seconds is not None
            else None
        )
        self._run_deadline_monotonic = deadline

        def left() -> float:
            return (
                deadline - time.monotonic()
                if deadline is not None
                else float("inf")
            )

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
        nonlocal_internal_error: list[str | None] = [None]
        self.exit_code = 0

        async def regulated_connect_loop() -> None:
            """The whole session lifecycle lives inside ONE CapExceeded
            envelope, so every control write — startup file_created/
            file_reopened/session_start, each connect/disconnect, planned or
            error disconnects — passes through the same cap-checked sink path.
            Nothing escapes the boundary."""
            nonlocal attempt, ready_announced
            try:
                self.open_for(day, write_clock=run_clock)
                if self.stop.is_set() or left() <= 0:
                    return
                self.emit_control(
                    "session_start",
                    write_clock=run_clock,
                    url=self.url,
                    streams=self.streams,
                    pid=os.getpid(),
                    max_bytes_per_day=self.max_bytes_per_day,
                    max_total_bytes=self.max_total_bytes,
                )
                self.save_state()
                while not self.stop.is_set() and left() > 0:
                    self.conn_seq += 1
                    cid = self._connection_id()
                    self.emit_control(
                        "connection_attempt",
                        conn=cid,
                        url=self.url,
                        streams=self.streams,
                    )
                    attempt_ready = False
                    try:
                        remaining = max(0.001, left())
                        async with websockets.connect(
                            self.url,
                            ping_interval=WS_PING_INTERVAL_S,
                            ping_timeout=WS_PING_TIMEOUT_S,
                            open_timeout=min(20.0, remaining),
                            close_timeout=min(10.0, remaining),
                        ) as ws:
                            await self.subscribe_and(
                                ws,
                                cid,
                                ack_timeout=min(
                                    30.0,
                                    max(0.001, left()),
                                ),
                            )
                            attempt = 0
                            self.emit_control(
                                "connect",
                                conn=cid,
                                url=self.url,
                                streams=self.streams,
                            )
                            self.save_state()
                            attempt_ready = True
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
                                recv_task = asyncio.ensure_future(ws.recv())
                                stop_wait_task = asyncio.ensure_future(
                                    self.stop.wait()
                                )
                                done, pending = await asyncio.wait(
                                    {recv_task, stop_wait_task},
                                    timeout=min(
                                        2.0,
                                        max(0.001, remaining),
                                    ),
                                    return_when=asyncio.FIRST_COMPLETED,
                                )
                                for pending_task in pending:
                                    pending_task.cancel()
                                if pending:
                                    await asyncio.gather(
                                        *pending,
                                        return_exceptions=True,
                                    )
                                if recv_task not in done:
                                    if stop_wait_task in done:
                                        break
                                    continue
                                frame = recv_task.result()
                                received_ms = int(time.time() * 1000)
                                try:
                                    msg = strict_json_loads(frame)
                                except UnicodeDecodeError:
                                    self.emit_control(
                                        "invalid_utf8_frame",
                                        conn=cid,
                                        reason="invalid_utf8",
                                    )
                                    raise ForceOrderSchemaError(
                                        "invalid_utf8"
                                    )
                                except NonFiniteJsonConstantError:
                                    self.emit_control(
                                        "non_finite_json_frame",
                                        conn=cid,
                                        reason="non_finite_json_constant",
                                    )
                                    raise ForceOrderSchemaError(
                                        "non_finite_json_constant"
                                    )
                                except DuplicateJsonKeyError:
                                    self.emit_control(
                                        "duplicate_json_frame",
                                        conn=cid,
                                        reason="duplicate_json_member",
                                    )
                                    raise ForceOrderSchemaError(
                                        "duplicate_json_member"
                                    )
                                except json.JSONDecodeError:
                                    self.emit_control(
                                        "bad_frame",
                                        conn=cid,
                                        frame_bytes=(
                                            len(frame.encode("utf-8"))
                                            if isinstance(frame, str)
                                            else len(frame)
                                        ),
                                    )
                                    continue
                                self.dispatch_market_message(
                                    msg,
                                    cid,
                                    received_ms,
                                )
                            self.emit_control("disconnect", conn=cid,
                                              reason="planned")
                            self.save_state()
                            if self.stop.is_set():
                                break
                            continue
                    except (
                        ForceOrderSchemaError,
                        websockets.exceptions.ConnectionClosed,
                        OSError,
                        socket.gaierror,
                        asyncio.TimeoutError,
                        TimeoutError,
                        RuntimeError,
                    ) as exc:
                        # Failure-before-ack and failure-after-ready both land
                        # here; the sink remains open through the control
                        # write below (all inside the cap envelope).
                        self.stats["reconnects"] += 1
                        attempt += 1
                        backoff = min(60.0, (2 ** min(attempt, 6)) * random.uniform(0.8, 1.2))
                        self.emit_control(
                            "disconnect",
                            conn=cid,
                            reason=bounded_audit_text(
                                f"{type(exc).__name__}: {exc}"
                            ),
                            backoff_s=round(backoff, 2),
                            phase=("post-ready" if attempt_ready else "pre-ack"))
                        self.save_state()
                        if self.stop.is_set() or left() <= 0:
                            break
                        await stop_aware_sleep(
                            min(backoff, max(0.0, left()))
                        )
            except RestartReplayAborted:
                self.stop.set()
            except CapExceeded as exc:
                # Terminal from ANY control write inside the envelope: freeze
                # accounting, close, save, signal out-of-band via return code.
                nonlocal_cap_reason[0] = str(exc)
                self.stop.set()

        oi_task: asyncio.Task | None = None
        try:
            oi_task = asyncio.ensure_future(self.oi_worker())
            connect_task = asyncio.ensure_future(regulated_connect_loop())
            while True:
                done, _pending = await asyncio.wait(
                    {connect_task, oi_task},
                    return_when=asyncio.FIRST_COMPLETED,
                )
                cap_terminal = None
                internal_terminal = None
                for task in done:
                    if task.cancelled():
                        continue
                    exc = task.exception()
                    if isinstance(exc, CapExceeded):
                        cap_terminal = str(exc)
                    elif exc is not None:
                        internal_terminal = type(exc).__name__
                if cap_terminal:
                    nonlocal_cap_reason[0] = cap_terminal
                    self.stop.set()
                    break
                if internal_terminal:
                    nonlocal_internal_error[0] = internal_terminal
                    self.stop.set()
                    break
                if nonlocal_cap_reason[0]:
                    self.stop.set()
                    break
                if connect_task in done:
                    self.stop.set()
                    oi_grace = 5.0
                    if deadline is not None and left() <= 0:
                        oi_grace = 0.0
                    if oi_grace > 0 and not oi_task.done():
                        try:
                            await asyncio.wait_for(
                                asyncio.shield(oi_task),
                                timeout=oi_grace,
                            )
                        except CapExceeded as exc:
                            nonlocal_cap_reason[0] = str(exc)
                        except (
                            asyncio.TimeoutError,
                            asyncio.CancelledError,
                        ):
                            pass
                        except Exception as exc:
                            nonlocal_internal_error[0] = type(exc).__name__
                    break
                if oi_task in done:
                    if self.stop.is_set():
                        if (
                            self._requested_stop_reason
                            and not connect_task.done()
                        ):
                            try:
                                await asyncio.wait_for(
                                    asyncio.shield(connect_task),
                                    timeout=min(
                                        2.5,
                                        max(0.001, left()),
                                    ),
                                )
                            except CapExceeded as exc:
                                nonlocal_cap_reason[0] = str(exc)
                            except (
                                asyncio.TimeoutError,
                                asyncio.CancelledError,
                            ):
                                pass
                            except Exception as exc:
                                nonlocal_internal_error[0] = type(exc).__name__
                        break
                    nonlocal_internal_error[0] = "UnexpectedOIWorkerExit"
                    self.stop.set()
                    break
            for task in (connect_task, oi_task):
                if not task.done():
                    task.cancel()
            for task in (connect_task, oi_task):
                try:
                    await task
                except (asyncio.CancelledError, CapExceeded, Exception):
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
            internal_error = nonlocal_internal_error[0]
            reason = (
                "cap"
                if cap_reason
                else "internal-error"
                if internal_error
                else self._requested_stop_reason
                or ("max-seconds" if self.max_seconds is not None else "manual")
            )
            if cap_reason:
                try:
                    self.fsync_close()
                except SinkIOError:
                    pass
                try:
                    self.save_state()
                except OSError as exc:
                    print(
                        f"STATE_SAVE_FAILED: {type(exc).__name__}",
                        file=sys.stderr,
                    )
                print(
                    f"CAP_EXCEEDED: {cap_reason}; fail-closed shutdown "
                    "(no data deleted)",
                    file=sys.stderr,
                )
                self.exit_code = 3
            elif self._sink_io_error is not None:
                try:
                    self.fsync_close()
                except SinkIOError:
                    pass
                try:
                    self.save_state()
                except OSError as exc:
                    print(
                        f"STATE_SAVE_FAILED: {type(exc).__name__}",
                        file=sys.stderr,
                    )
                print(
                    f"COLLECTOR_INTERNAL_ERROR: {self._sink_io_error}; "
                    "sink frozen",
                    file=sys.stderr,
                )
                self.exit_code = 5
            else:
                try:
                    self.emit_control("session_stop", reason=reason)
                    self.fsync_close()
                    self.save_state()
                    if internal_error:
                        print(
                            f"COLLECTOR_INTERNAL_ERROR: {internal_error}; "
                            "fail-closed shutdown",
                            file=sys.stderr,
                        )
                        self.exit_code = 5
                    else:
                        self.exit_code = 0
                except SinkIOError:
                    try:
                        self.save_state()
                    except OSError as exc:
                        print(
                            f"STATE_SAVE_FAILED: {type(exc).__name__}",
                            file=sys.stderr,
                        )
                    print(
                        f"COLLECTOR_INTERNAL_ERROR: {self._sink_io_error}; "
                        "sink frozen",
                        file=sys.stderr,
                    )
                    self.exit_code = 5
                except CapExceeded as exc:
                    try:
                        self.fsync_close()
                    except SinkIOError:
                        pass
                    self.save_state()
                    print(
                        f"CAP_EXCEEDED during shutdown control: {exc}; "
                        "fail-closed (no data deleted)",
                        file=sys.stderr,
                    )
                    self.exit_code = 3
            self._quiesce_oi_executor()
        self.release_lease()
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
    config_error = stream_config_error(args.streams)
    if config_error is not None:
        p.error(f"--streams {config_error}")
    return args


def main() -> int:
    args = parse_args()
    try:
        collector = Collector(
            args.out,
            args.url,
            args.streams,
            args.max_seconds,
            args.max_bytes_per_day,
            args.max_total_bytes,
        )
    except CollectorLeaseError as exc:
        print(f"LIQUIDATION_CAPTURE_LOCKED {exc}", file=sys.stderr)
        return 6
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, collector.request_stop)
    try:
        # run() handles every terminal condition internally: caps return 3,
        # pre-existing exhaustion returns 4, SIGTERM/max-seconds return 0.
        return loop.run_until_complete(collector.run())
    finally:
        collector.release_lease()
        loop.close()


if __name__ == "__main__":
    raise SystemExit(main())
