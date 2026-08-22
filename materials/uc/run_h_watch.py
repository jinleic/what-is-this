"""Watch campaign H: replay each committed slice, then run the frozen collector.

Runs every cert3_replay invocation as
  python -B <campaign>/snapshot/cert3_replay.py
so the pinned snapshot is never polluted by bytecode caches.  Marker files
.replay_slice<i>_OK and collect_output.txt land in the campaign directory,
NOT inside snapshot/, so the pristine-snapshot listing is untouched.  The
collector itself still re-replays all eight traces before declaring the
certificate; the per-slice replays here are an early-warning lane plus
independent confirmation.
"""

import glob
import json
import os
import subprocess
import sys
import time

CAMPAIGN = (sys.argv[1] if len(sys.argv) > 1 else
            "/Users/jinleic/jinleic-workspace/math/uc/campaigns/cert3_20260817T212001Z_4e6268eb7ce44aa3883ef15ea401c7c7_2f23a58ebdb8")  # noqa: E501
SNAP = os.path.join(CAMPAIGN, "snapshot")
MATH = "/Users/jinleic/jinleic-workspace/math"
VERIFY_PYTHON = os.path.join(MATH, ".venv", "bin", "python")
NSLICES = 8
POLL_S = 60


def log(msg):
    print("[watch %s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def replay_slice(s):
    res = glob.glob(os.path.join(CAMPAIGN, "result_slice%d_*.json" % s))
    trc = glob.glob(os.path.join(CAMPAIGN, "trace_slice%d_*.bin" % s))
    marker = os.path.join(CAMPAIGN, ".replay_slice%d_OK" % s)
    if os.path.exists(marker):
        # Markers are only ever written after a verified pass (rc==0 and a
        # replayed tally JSON, validated here or manually against the worker
        # record).  The collector gate requires "pass" on all eight slices.
        return "pass"
    if not res or not trc:
        return "absent"
    cmd = [VERIFY_PYTHON, "-B", os.path.join(SNAP, "cert3_replay.py"),
           os.path.join(CAMPAIGN, "launch.json"), res[0], trc[0]]
    out = os.path.join(CAMPAIGN, "replay_slice%d.txt" % s)
    log("replay slice %d: %s" % (s, os.path.basename(res[0])))
    with open(out, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=MATH)
    text = open(out).read()
    # cert3_replay prints the replayed-tally JSON and returns 0 on success;
    # it validates the worker's tallies itself and raises on any mismatch.
    ok = False
    if proc.returncode == 0:
        try:
            replayed = json.loads(text.strip().splitlines()[-1])
            ok = replayed.get("processed", 0) > 0 and replayed["residual"] == 0
        except (json.JSONDecodeError, KeyError, IndexError):
            ok = False
    if ok:
        with open(marker, "w") as fh:
            fh.write(text)
        log("replay slice %d PASS" % s)
        return "pass"
    log("replay slice %d FAIL (rc=%d); last lines: %s"
        % (s, proc.returncode, text.strip().splitlines()[-1]
           if text.strip() else "<empty>"))
    return "fail"


def run_collector():
    out_path = os.path.join(CAMPAIGN, "collect_output.txt")
    if os.path.exists(out_path):
        # A stored file alone is no certificate: the first run may have been
        # rejected before all slices committed.  Only its content counts.
        stored = open(out_path).read()
        if "COMPOSITE CERTIFICATE" in stored:
            return True
    cmd = [VERIFY_PYTHON, "-B", os.path.join(SNAP, "cert3_collect.py"),
           os.path.join(CAMPAIGN, "launch.json")]
    log("running frozen collector")
    with open(out_path, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=MATH)
    text = open(out_path).read()
    ok = proc.returncode == 0 and "COMPOSITE CERTIFICATE" in text
    log("collector: rc=%d certificate=%s" % (proc.returncode, ok))
    return ok


def main():
    states = {}
    while True:
        for s in range(NSLICES):
            if states.get(s) in ("pass", "fail"):
                continue
            states[s] = replay_slice(s)
        if all(states.get(s) == "pass" for s in range(NSLICES)):
            if run_collector():
                log("CAMPAIGN H CERTIFIED")
                return 0
            log("collector did not certify; continuing to poll")
        time.sleep(POLL_S)


if __name__ == "__main__":
    sys.exit(main())
