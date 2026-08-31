#!/usr/bin/env python3
"""Detached scheduler for qendpoint-lift slices 02..16 and their lineage closes.

Reproduces the committed slice-01 recipe per slice K in 2..16 (per job):

  1. lift:   python -I -B uc/liu9_qendpoint_lift.py
                 --source  uc/verification/results/liu9-complement-residual-rho1-1728-100k.json
                 --budget 5000001 --seconds 86400
                 --slice-index K --slice-count 16
                 --checkpoint uc/verification/results/queue_state/liu9-qendpoint-slice-<KK>-checkpoint.json
                 --output     uc/verification/results/liu9-qendpoint-slice-<KK>-of-16.json

  2. close:  python -I -B uc/verification/liu9_qendpoint_slice_lineage_close.py
                 --source uc/verification/results/liu9-qendpoint-slice-<KK>-of-16.json
                 --output uc/verification/results/liu9-qendpoint-slice<KK>-lineage-close.json

Scheduling. A fixed pool of 5 single-core workers owns the fixed 15-job queue
(slices 2..16, ascending).  Worker assignment is FIFO over the static job
list: job K is handed to the first worker whose slot is free.  No wall time,
load average, or random jitter enters any decision, so the schedule is a
deterministic function of job completion order.  Each slice job additionally
holds an exclusive lock file (queue_state/locks/slice<KK>.lock) so a slice can
never run twice concurrently.

Resumability. The lift restarts from scratch whenever the same job-level
invocation is re-run (its run() is not checkpoint-resumable; the checkpoint
file is an observation log).  This scheduler therefore treats an existing
final report as the resumable boundary: re-invoking the queue skips finished
stages and only rebuilds missing artifacts.  Slice-level single-flight is
enforced by the lock file.

State. queue_state/queue.status is a JSON object rewritten atomically
(tmp+rename) on each deterministic transition.  queue_state/queue.done gets
one line per finished job: "<KK> <exit_code>".

Exit of the top-level queue: 0 iff every slice's lineage close exists with
compute_unresolved_lineages == 0; exit 2 otherwise.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
UC = HERE.parent
MATH = UC.parent
VENV_PY = MATH / ".venv/bin/python"
RESULTS = UC / "verification/results"
QUEUE_DIR = RESULTS / "queue_state"
LOCKS = QUEUE_DIR / "locks"
LOGS = RESULTS / "queue_logs"

LIFT = UC / "liu9_qendpoint_lift.py"
CLOSER = HERE / "liu9_qendpoint_slice_lineage_close.py"
SOURCE = RESULTS / "liu9-complement-residual-rho1-1728-100k.json"
BUDGET = 5_000_001
SECONDS_CAP = 86_400.0
SLICES = list(range(2, 17))
MAX_CONCURRENT = 5

STATE_PATH = QUEUE_DIR / "queue.status"
DONE_PATH = QUEUE_DIR / "queue.done"


def slice_label(slice_index: int) -> str:
    """2..16 -> 02..16."""
    return "%02d" % slice_index


def report_name(slice_index: int) -> str:
    return "liu9-qendpoint-slice-%s-of-16.json" % slice_label(slice_index)


def close_name(slice_index: int) -> str:
    return "liu9-qendpoint-slice%s-lineage-close.json" % slice_label(slice_index)


def checkpoint_name(slice_index: int) -> str:
    return "liu9-qendpoint-slice-%s-checkpoint.json" % slice_label(slice_index)


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    os.replace(tmp, path)


def read_state() -> dict:
    try:
        return json.loads(STATE_PATH.read_text())
    except (OSError, ValueError):
        return {}


def log(line: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    print("%s %s" % (stamp, line), flush=True)


def append_done(slice_index: int, exit_code: int) -> None:
    with DONE_PATH.open("a", encoding="utf-8") as handle:
        handle.write("%s %d\n" % (slice_label(slice_index), exit_code))


def acquire_lock(slice_index: int) -> bool:
    LOCKS.mkdir(parents=True, exist_ok=True)
    path = LOCKS / ("slice%s.lock" % slice_label(slice_index))
    try:
        handle = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    os.write(handle, b"%d\n" % os.getpid())
    os.close(handle)
    return True


def release_lock(slice_index: int) -> None:
    path = LOCKS / ("slice%s.lock" % slice_label(slice_index))
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def launch(cmd: list[str], log_path: Path, banner: str) -> int:
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n===== %s =====\n" % banner)
        handle.flush()
        process = subprocess.Popen(
            cmd, stdout=handle, stderr=subprocess.STDOUT,
            env=child_env(), cwd=str(MATH),
        )
        return process.wait()


def child_env() -> dict[str, str]:
    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "PYTHONHASHSEED": "0",
    })
    return env


def run_slice(slice_index: int) -> int:
    """Lift then lineage-close one slice. Returns 0 iff fully closed."""
    label = slice_label(slice_index)
    report_path = RESULTS / report_name(slice_index)
    checkpoint_path = QUEUE_DIR / checkpoint_name(slice_index)
    close_path = RESULTS / close_name(slice_index)
    slice_log = LOGS / ("slice%s.log" % label)

    if not acquire_lock(slice_index):
        log("slice %s: another run holds the lock" % label)
        return 3
    try:
        if not report_path.exists():
            code = launch(
                [
                    str(VENV_PY), "-I", "-B", str(LIFT),
                    "--source", str(SOURCE),
                    "--budget", str(BUDGET),
                    "--seconds", str(SECONDS_CAP),
                    "--slice-index", str(slice_index),
                    "--slice-count", "16",
                    "--checkpoint", str(checkpoint_path),
                    "--output", str(report_path),
                ],
                slice_log, "lift slice %s" % label)
            if code != 0 or not report_path.exists():
                log("slice %s: lift exit %d" % (label, code))
                return 1
        if not close_path.exists():
            code = launch(
                [
                    str(VENV_PY), "-I", "-B", str(CLOSER),
                    "--source", str(report_path),
                    "--output", str(close_path),
                ],
                slice_log, "lineage close slice %s" % label)
            if code != 0 or not close_path.exists():
                log("slice %s: lineage close exit %d" % (label, code))
                return 2
        return 0
    finally:
        release_lock(slice_index)


def mark(slice_index: int, patch: dict) -> dict:
    state = read_state()
    slices = state.setdefault("slices", {})
    entry = slices.setdefault(str(slice_index), {})
    entry.update(patch)
    _atomic_json(STATE_PATH, state)
    return state


def worker_main(slice_index: int) -> int:
    code = run_slice(slice_index)
    log("slice %s: job exit %d" % (slice_label(slice_index), code))
    return code


def queue_main() -> int:
    QUEUE_DIR.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    if not SOURCE.exists():
        raise SystemExit("missing source %s" % SOURCE)
    if not VENV_PY.exists():
        raise SystemExit("missing interpreter %s" % VENV_PY)
    state = read_state()
    state.setdefault("started_utc", time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    state["budget"] = BUDGET
    state["slice_count"] = 16
    state["max_concurrent"] = MAX_CONCURRENT
    state["slices"] = state.get("slices", {})
    _atomic_json(STATE_PATH, state)
    log("queue start: slices %d..%d, %d workers, budget %d"
        % (SLICES[0], SLICES[-1], MAX_CONCURRENT, BUDGET))

    running: dict[int, subprocess.Popen] = {}
    codes: dict[int, int] = {}
    pending = [s for s in SLICES if str(s) not in state["slices"]
               or state["slices"][str(s)].get("stage") != "done"]
    log("pending after resume scan: %s" % pending)
    while pending or running:
        while pending and len(running) < MAX_CONCURRENT:
            slice_index = pending.pop(0)
            if not acquire_lock(slice_index):
                log("slice %s: lock held; deferring to queue tail"
                    % slice_label(slice_index))
                pending.append(slice_index)
                time.sleep(20)
                continue
            release_lock(slice_index)  # owned by the worker, not the queue
            mark(slice_index, {"stage": "running"})
            handle = subprocess.Popen(
                [str(VENV_PY), "-I", "-B", "-u", os.path.abspath(__file__),
                 "worker", str(slice_index)],
                cwd=str(MATH),
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                env=child_env(),
            )
            running[slice_index] = handle
            log("slice %s: worker pid %d launched"
                % (slice_label(slice_index), handle.pid))
        for slice_index, process in list(running.items()):
            if process.poll() is not None:
                del running[slice_index]
                codes[slice_index] = process.returncode or 0
                mark(slice_index, {"stage": "done",
                                   "exit_code": codes[slice_index]})
                append_done(slice_index, codes[slice_index])
        time.sleep(15)

    composite = ("CLOSED" if all(codes.get(s) == 0 for s in SLICES)
                 else "INCOMPLETE")
    state = read_state()
    state["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    state["composite"] = composite
    _atomic_json(STATE_PATH, state)
    log("queue done: %s" % composite)
    return 0 if composite == "CLOSED" else 2


def main(argv: list[str]) -> int:
    """`worker K` runs one slice; no-args runs the queue; anything else errors."""
    if not argv:
        return queue_main()
    if len(argv) == 2 and argv[0] == "worker":
        return worker_main(int(argv[1]))
    raise SystemExit(
        "usage: run_qendpoint_slice_queue.py [worker SLICE_INDEX] "
        "(no flags; the queue takes no options)")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
