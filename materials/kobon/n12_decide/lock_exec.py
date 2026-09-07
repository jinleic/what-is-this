#!/usr/bin/env python3
"""Run one command while holding a permanent kernel-enforced file lock."""

from __future__ import annotations

import argparse
import errno
import fcntl
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

BUSY_EXIT = 75
TERM_GRACE_SECONDS = 3.0


def _terminate_group(
        child: subprocess.Popen,
        grace_seconds: float = TERM_GRACE_SECONDS) -> None:
    """Terminate the whole child group, escalate, and reap its leader."""
    group = child.pid
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        pass
    deadline = time.monotonic() + grace_seconds
    while True:
        try:
            os.killpg(group, 0)
        except ProcessLookupError:
            break
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            try:
                os.killpg(group, signal.SIGKILL)
            except ProcessLookupError:
                pass
            break
        child.poll()
        time.sleep(min(0.05, remaining))
    if child.poll() is None:
        child.wait()


def run(lock_path: Path, command: list[str], *, nonblocking: bool) -> int:
    if not command:
        raise ValueError("a command is required")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o600)
    operation = fcntl.LOCK_EX | (fcntl.LOCK_NB if nonblocking else 0)
    try:
        try:
            fcntl.flock(descriptor, operation)
        except OSError as error:
            if error.errno in {
                errno.EACCES, errno.EAGAIN, errno.EWOULDBLOCK
            }:
                print(f"REFUSED: lock held: {lock_path}", file=sys.stderr)
                return BUSY_EXIT
            raise

        watched_signals = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
        watched_set = set(watched_signals)
        previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, watched_set)
        previous_handlers = {
            signum: signal.getsignal(signum) for signum in watched_signals
        }
        pending_signal = None
        child = None

        def defer_interrupt(signum, _frame) -> None:
            nonlocal pending_signal
            if pending_signal is None:
                pending_signal = signum

        try:
            for signum in watched_signals:
                signal.signal(signum, defer_interrupt)
            signal.pthread_sigmask(signal.SIG_UNBLOCK, watched_set)
            if pending_signal is not None:
                return 128 + pending_signal
            child = subprocess.Popen(
                command, start_new_session=True, pass_fds=(descriptor,))
            while True:
                if pending_signal is not None:
                    _terminate_group(child)
                    return 128 + pending_signal
                try:
                    return child.wait(timeout=0.1)
                except subprocess.TimeoutExpired:
                    continue
        except BaseException:
            if child is not None and child.poll() is None:
                _terminate_group(child)
            raise
        finally:
            signal.pthread_sigmask(signal.SIG_BLOCK, watched_set)
            try:
                for signum, handler in previous_handlers.items():
                    signal.signal(signum, handler)
            finally:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nonblocking", action="store_true")
    parser.add_argument("lock", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    try:
        return run(args.lock, command, nonblocking=args.nonblocking)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())
