#!/usr/bin/env python3
"""Locking and durable-write primitives shared by the r55 certificate campaigns.

Campaign singleton locks use exclusive nonblocking ``flock``. Host-heavy
critical sections may request a blocking lock. Lock files are never unlinked,
so exclusion is kernel-enforced and cannot be defeated by stale-lock cleanup.
"""

from __future__ import annotations

import errno
import fcntl
import json
import os
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
SUPERVISION_DIR = WORKSPACE / "scratch/r55-involution-f5"
DRIVER_LOCK = SUPERVISION_DIR / ".certificate-driver.lock"
SCHEDULER_LOCK = SUPERVISION_DIR / ".certificate-scheduler.lock"
HOST_HEAVY_LOCK = WORKSPACE / "scratch/.host-heavy-job.lock"


class LockBusy(RuntimeError):
    """Another live process already holds the requested campaign lock."""


class FileLock:
    """Exclusive advisory lock on a permanent lock file."""

    def __init__(self, path: Path, *, blocking: bool = False) -> None:
        self.path = Path(path)
        self.blocking = blocking
        self._stream = None

    def acquire(self) -> None:
        if self._stream is not None:
            raise RuntimeError(f"{self.path}: already held by this instance")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        stream = self.path.open("a")
        try:
            operation = fcntl.LOCK_EX
            if not self.blocking:
                operation |= fcntl.LOCK_NB
            fcntl.flock(stream.fileno(), operation)
        except OSError as error:
            stream.close()
            if error.errno in {
                    errno.EACCES, errno.EAGAIN, errno.EWOULDBLOCK}:
                raise LockBusy(
                    f"{self.path} is held by another process") from error
            raise
        self._stream = stream

    def fileno(self) -> int:
        if self._stream is None:
            raise RuntimeError(f"{self.path}: lock is not held")
        return self._stream.fileno()

    def release(self) -> None:
        stream, self._stream = self._stream, None
        if stream is None:
            return
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        finally:
            stream.close()

    def __enter__(self) -> "FileLock":
        self.acquire()
        return self

    def __exit__(self, *_exception) -> None:
        self.release()


def atomic_write_json(path: Path, document: Any) -> None:
    """Durably replace *path* with canonical JSON via a sibling temp file.

    The temp file is fsynced, atomically renamed over the destination, and the
    containing directory is fsynced so the replacement survives power loss.
    """
    path = Path(path)
    temp = path.with_name(f"{path.name}.tmp")
    try:
        payload = json.dumps(
            document, allow_nan=False, indent=2, sort_keys=True) + "\n"
        with temp.open("w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    except BaseException:
        temp.unlink(missing_ok=True)
        raise
    directory = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def load_json_object(path: Path, error: type[Exception]) -> dict:
    """Strictly load *path* as a JSON object or raise *error* with context."""
    path = Path(path)

    def unique_object(pairs):
        document = {}
        for key, value in pairs:
            if key in document:
                raise error(f"duplicate JSON member {key!r} in {path}")
            document[key] = value
        return document

    def reject_constant(value):
        raise error(f"nonstandard JSON constant {value} in {path}")

    try:
        document = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as detail:
        raise error(f"cannot load {path}: {detail}") from detail
    if type(document) is not dict:
        raise error(f"{path} root is not an object")
    return document
