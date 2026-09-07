#!/usr/bin/env python3
"""Finish both fixed-five certificate batches with one CakePB process at a time.

The two batch drivers are independently resumable at source-orbit boundaries, but
running them together drove this 96 GiB host to 72 GiB of swap and consumed about
30 GiB of disk per hour.  This scheduler preserves the measured safe invariant:
only one driver (therefore one CakePB process) runs at a time.  It gives the larger
W-negative partition three source orbits per cycle and the hierarchy one.

Safety envelope: one scheduler instance at a time holds a permanent
scheduler-instance lock, and both batch drivers share a separate permanent
driver lock, so a supervised child can never deadlock against its parent and a
direct driver can never overlap a supervised one.  While a child runs, the
free-disk floor is polled and a breach terminates and reaps the child's whole
process group.  Timeout retries deduct the records the failed attempt already
persisted from the per-cycle allowance, any terminal failure is journaled, and
a later restart stops with success at that journal entry until a human passes
--resume-after-failure.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import campaign_runtime

ROOT = Path(__file__).resolve().parents[1]
MATH = ROOT.parent
DATA = ROOT / "data"
W_OUTPUT = DATA / "involution_f5_w_dfs_certificates.json"
H_OUTPUT = DATA / "involution_f5_hierarchical_certificates.json"
W_DRIVER = ROOT / "src" / "certify_involution_f5_w_dfs_batch.py"
H_DRIVER = ROOT / "src" / "certify_involution_f5_hierarchical_batch.py"
EVENTS = DATA / "involution_f5_certificate_scheduler.jsonl"
MIN_FREE_GIB = 20
POLL_SECONDS = 30.0
TERM_GRACE_SECONDS = 10.0
GROUP_POLL_SECONDS = 0.1


class SchedulerFailure(RuntimeError):
    """A batch rejected a proof or the host crossed a safety boundary."""


def _coverage(path: Path, key: str, total: int) -> tuple[int, int]:
    if not path.is_file():
        return 0, total
    document = campaign_runtime.load_json_object(path, SchedulerFailure)
    coverage = document.get("coverage")
    records = document.get("records")
    if type(coverage) is not dict or type(records) is not list:
        raise SchedulerFailure(f"invalid coverage metadata in {path}")
    done = coverage.get("certified_support_orbits")
    target = coverage.get(key)
    remaining = coverage.get("remaining_support_orbits")
    complete = coverage.get("campaign_complete")
    if (type(done) is not int or type(target) is not int
            or type(remaining) is not int
            or not 0 <= done <= total
            or target != total
            or remaining != total - done
            or complete is not (done == total)
            or len(records) != done):
        raise SchedulerFailure(f"invalid coverage metadata in {path}")
    seen: set[int] = set()
    for record in records:
        source = record.get("source_index") if type(record) is dict else None
        if type(source) is not int or source in seen:
            raise SchedulerFailure(f"invalid coverage records in {path}")
        seen.add(source)
    return done, total


def _snapshot() -> dict[str, int]:
    w_done, w_total = _coverage(
        W_OUTPUT, "target_w_negative_support_orbits", 550)
    h_done, h_total = _coverage(
        H_OUTPUT, "target_deeper_support_orbits", 155)
    return {
        "w_done": w_done,
        "w_remaining": w_total - w_done,
        "h_done": h_done,
        "h_remaining": h_total - h_done,
    }


def _emit(event: str, **fields) -> None:
    record = {
        "time_unix": int(time.time()),
        "event": event,
        **fields,
    }
    line = json.dumps(record, sort_keys=True)
    with EVENTS.open("a") as stream:
        stream.write(line + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    print(line, flush=True)


def _guard_host() -> None:
    active = subprocess.run(
        ["pgrep", "-x", "cake_pb"], check=False,
        capture_output=True, text=True)
    if active.returncode == 0 and active.stdout.strip():
        raise SchedulerFailure(
            f"refusing concurrent CakePB processes: {active.stdout.strip()}")
    free_gib = shutil.disk_usage(MATH).free / (1 << 30)
    if free_gib < MIN_FREE_GIB:
        raise SchedulerFailure(
            f"only {free_gib:.1f} GiB free; require {MIN_FREE_GIB} GiB")


def _last_event() -> str | None:
    """Return the event name of the latest parseable scheduler log record."""
    if not EVENTS.is_file():
        return None
    event = None
    with EVENTS.open() as stream:
        for line in stream:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if type(record) is dict and type(record.get("event")) is str:
                event = record["event"]
    return event


def _terminate_group(child: subprocess.Popen) -> None:
    """Gracefully stop the whole child group, then reap its direct leader."""
    group = child.pid
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        pass
    deadline = time.monotonic() + TERM_GRACE_SECONDS
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
        time.sleep(min(GROUP_POLL_SECONDS, remaining))
    if child.poll() is None:
        child.wait()


def _raise_interrupt(signum, _frame) -> None:
    raise KeyboardInterrupt(f"received {signal.Signals(signum).name}")


def _supervised_run(
        command: list[str],
        environment: dict[str, str],
) -> tuple[subprocess.CompletedProcess, str | None]:
    """Run a driver child while policing disk and every exceptional exit."""
    watched_signals = (signal.SIGINT, signal.SIGTERM)
    watched_set = set(watched_signals)
    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, watched_set)
    previous_handlers = {
        signal_number: signal.getsignal(signal_number)
        for signal_number in watched_signals
    }
    pending_signal = None
    child = None

    def defer_interrupt(signum, _frame) -> None:
        nonlocal pending_signal
        if pending_signal is None:
            pending_signal = signum

    try:
        for signal_number in watched_signals:
            signal.signal(signal_number, defer_interrupt)
        signal.pthread_sigmask(signal.SIG_UNBLOCK, watched_set)
        if pending_signal is not None:
            _raise_interrupt(pending_signal, None)
        try:
            child = subprocess.Popen(
                command, cwd=MATH, env=environment,
                start_new_session=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True)
        except OSError as error:
            raise SchedulerFailure(
                f"cannot start batch driver: {error}") from error
        for signal_number in watched_signals:
            signal.signal(signal_number, _raise_interrupt)
        if pending_signal is not None:
            _raise_interrupt(pending_signal, None)
        while True:
            try:
                stdout, stderr = child.communicate(timeout=POLL_SECONDS)
            except subprocess.TimeoutExpired:
                free_gib = shutil.disk_usage(MATH).free / (1 << 30)
                if free_gib >= MIN_FREE_GIB:
                    continue
                _terminate_group(child)
                stdout, stderr = child.communicate()
                breach = (
                    f"only {free_gib:.1f} GiB free; require "
                    f"{MIN_FREE_GIB} GiB; batch process group terminated")
                return (subprocess.CompletedProcess(
                            command, child.returncode, stdout, stderr),
                        breach)
            completed = subprocess.CompletedProcess(
                command, child.returncode, stdout, stderr)
            _terminate_group(child)
            return completed, None
    except BaseException:
        if child is not None:
            _terminate_group(child)
            stdout, stderr = child.communicate()
            if stdout:
                print(stdout, end="", flush=True)
            if stderr:
                print(stderr, end="", file=sys.stderr, flush=True)
        raise
    finally:
        signal.pthread_sigmask(signal.SIG_BLOCK, watched_set)
        try:
            for signal_number, handler in previous_handlers.items():
                signal.signal(signal_number, handler)
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


def _run_batch(kind: str, limit: int, timeout: float,
               max_timeout: float) -> None:
    driver = W_DRIVER if kind == "w" else H_DRIVER
    done_key = f"{kind}_done"
    before = _snapshot()
    current_timeout = timeout
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(MATH), str(ROOT / "src")])

    while True:
        _guard_host()
        command = [
            sys.executable, "-B", str(driver),
            "--timeout-per-step", str(current_timeout),
            "--max-new-records", str(limit),
        ]
        _emit("batch_start", kind=kind, limit=limit,
              timeout_per_step=current_timeout, **before)
        completed, breach = _supervised_run(command, environment)
        if completed.stdout:
            print(completed.stdout, end="", flush=True)
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr, flush=True)
        after = _snapshot()
        _emit("batch_exit", kind=kind, returncode=completed.returncode,
              timeout_per_step=current_timeout,
              terminated=breach is not None, **after)
        if breach is not None:
            raise SchedulerFailure(breach)
        if completed.returncode == 0:
            return
        output = completed.stdout + completed.stderr
        if "TimeoutExpired" not in output:
            raise SchedulerFailure(
                f"{kind} batch rejected or failed; not retrying as timeout")
        persisted = after[done_key] - before[done_key]
        if not 0 <= persisted <= limit:
            raise SchedulerFailure(
                f"{kind} batch persisted impossible progress {persisted} "
                f"against allowance {limit}")
        limit -= persisted
        if limit == 0:
            return
        if current_timeout >= max_timeout:
            raise SchedulerFailure(
                f"{kind} batch exceeded maximum timeout {max_timeout}s")
        current_timeout = min(current_timeout * 2, max_timeout)
        _emit("timeout_retry", kind=kind, remaining_allowance=limit,
              timeout_per_step=current_timeout, **after)
        before = after


def _validate_complete_ledgers() -> None:
    import certify_involution_f5_hierarchical_batch as h_batch
    import certify_involution_f5_w_dfs_batch as w_batch

    try:
        with campaign_runtime.FileLock(campaign_runtime.DRIVER_LOCK):
            w_document = w_batch._validate_output_unlocked(
                W_OUTPUT, require_complete=True)
            h_document = h_batch._validate_output_unlocked(
                H_OUTPUT, require_complete=True)
    except (campaign_runtime.LockBusy, h_batch.HierarchyViolation,
            w_batch.BatchViolation, OSError) as error:
        raise SchedulerFailure(
            f"complete ledger validation failed: {error}") from error
    w_sources = {record["source_index"] for record in w_document["records"]}
    h_sources = {record["source_index"] for record in h_document["records"]}
    if (len(w_sources) != 550 or len(h_sources) != 155
            or w_sources & h_sources
            or len(w_sources | h_sources) != 705):
        raise SchedulerFailure(
            "complete certificate ledgers do not partition 705 sources")


def run(w_per_cycle: int, h_per_cycle: int, timeout: float,
        max_timeout: float, one_cycle: bool) -> dict[str, int]:
    """Run the fair 3:1 campaign; the caller holds the scheduler-instance lock."""
    try:
        return _campaign(w_per_cycle, h_per_cycle, timeout, max_timeout,
                         one_cycle)
    except SchedulerFailure as error:
        try:
            state = _snapshot()
        except SchedulerFailure:
            state = {}
        _emit("terminal_failure", detail=str(error), **state)
        raise


def _campaign(w_per_cycle: int, h_per_cycle: int, timeout: float,
              max_timeout: float, one_cycle: bool) -> dict[str, int]:
    while True:
        state = _snapshot()
        if state["w_remaining"] == 0 and state["h_remaining"] == 0:
            _validate_complete_ledgers()
            _emit("campaign_complete", **state)
            return state
        if state["w_remaining"]:
            _run_batch(
                "w", min(w_per_cycle, state["w_remaining"]),
                timeout, max_timeout)
        state = _snapshot()
        if state["h_remaining"]:
            _run_batch(
                "h", min(h_per_cycle, state["h_remaining"]),
                timeout, max_timeout)
        if one_cycle:
            state = _snapshot()
            _emit("cycle_complete", **state)
            return state


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--w-per-cycle", type=int, default=3)
    parser.add_argument("--h-per-cycle", type=int, default=1)
    parser.add_argument("--timeout-per-step", type=float, default=21600.0)
    parser.add_argument("--max-timeout-per-step", type=float, default=86400.0)
    parser.add_argument("--one-cycle", action="store_true")
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--resume-after-failure", action="store_true")
    args = parser.parse_args(argv)
    if args.status:
        print(json.dumps(_snapshot(), sort_keys=True))
        return 0
    if (args.w_per_cycle <= 0 or args.h_per_cycle <= 0
            or args.timeout_per_step <= 0
            or args.max_timeout_per_step < args.timeout_per_step):
        parser.error("positive batch sizes and ordered positive timeouts required")
    with campaign_runtime.FileLock(campaign_runtime.SCHEDULER_LOCK):
        if _last_event() == "terminal_failure":
            if not args.resume_after_failure:
                print(json.dumps(
                    {"terminal_failure_suppressed": True}, sort_keys=True))
                return 0
            _emit("failure_acknowledged", resuming=True)
        print(json.dumps(run(
            args.w_per_cycle, args.h_per_cycle, args.timeout_per_step,
            args.max_timeout_per_step, args.one_cycle), sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
