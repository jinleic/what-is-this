#!/usr/bin/env python3
"""Two-rung replay launcher (v11 corrected interval core, static freeze).

COMPUTATIONAL corrected-core DIAGNOSTIC only: measures the old-protocol
outputs faithfully and compares them; historical `certified` field names
are SOURCE SCHEMA ONLY; formal two-rung/omega claims remain SUSPENDED
pending (1) outward end-to-end aggregate rewrite and (2) mathematical
justification of the feasibility absorption.

Hard-requires:
  - five thread env vars set to 1
  - __debug__ mode (asserts active)
  - current process niceness >= 10
  - sequential only (rung1 -> rung2 -> rung2 without)
  - 900s alarm/cap per rung

Usage:
  launcher.py <rung1_out> <rung2_out> <rung2_without_out> <compare_out>
"""
import hashlib
import json
import os
import sys
import subprocess
import time

REQUIRED_ENVS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                 "NUMEXPR_NUM_THREADS")
HERE = os.path.dirname(os.path.abspath(__file__))
BASE_TARGET = 2.371866
RELEASE_ENV = "OMEGA_TWO_RUNG_RELEASED"
EXPECTED_OUTPUT_NAMES = {
    "rung1_output": "rung1_replay.json",
    "rung2_with_output": "rung2_with_replay.json",
    "rung2_without_output": "rung2_without_replay.json",
    "comparison_summary": "comparison_summary.json",
    "rung1_stderr_log": "rung1_stderr.log",
    "rung1_stdout_log": "rung1_stdout.json",
    "rung2_with_stderr_log": "rung2_stderr.log",
    "rung2_with_stdout_log": "rung2_stdout.json",
    "rung2_without_stderr_log": "rung2_without_stderr.log",
    "rung2_without_stdout_log": "rung2_without_stdout.json",
}

# Hard-coded expected dependency hashes (MAIN directive: never mutate).
EXPECTED_INPUT_SHA256 = {
    "interval_core.py": "ce70f95922f3d5098e892e1b7047570ebb5267826b51f2776c91b56d9934adb5",
    "alman25_float.py": "02910b117c1c689ea29ee1fb343fa272bc3ff4c29ca1c353f7a9b3f76da1c448",
    "vxxz24_float.py": "dd8edaa6a6be50d9bd9709c4fc0733cbd74436db577164ae3431e7d0b5bab1b1",
    "stage_b_rung1_certified.py": "10f0396d13a9bc5d0b2f68b750cb7ed45af26f3ee7a199d6ca49e90b59a818fa",
    "stage_b_rung2.py": "442f6b8008fd8d2cc44708f8505aae24a4e85d76ca44392636ecd3b5bde3f9d4",
    "K100_2.37155181.mat": "f15684a43d0193fdd2c1ceba9e1001192f7b2bad36f0ca715636dcecf9e66eba",
    "W1.00_2.371339.mat": "783353fda82acb3fb93c247dcad857b2db5f61944f5d0e91ae5f9e5a6c7feec3",
}

PINNED_ENV = {k: "1" for k in REQUIRED_ENVS}


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


EXPECTED_PARAMS_BYTES_SHA256 = {
    "W1.00_2.371339.mat":
        "f23369136314cc51497d0c468c00d69a02b74d2b665073b9f934240ab428b46f",
    "K100_2.37155181.mat":
        "df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0722c3e25871da93d",
}


def verify_and_lock_hashes():
    """Hard-verify file hashes AND params-byte hashes, return all."""
    actual = {}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        path = os.path.join(HERE, name)
        if not os.path.isfile(path):
            raise RuntimeError(f"missing frozen dep {name}")
        actual[name] = _sha(path)
        if actual[name] != expected:
            raise RuntimeError(
                f"{name}: sha256 {actual[name]} != expected {expected}; "
                "refusing to launch")
    # Verify params-bytes hashes (Main directive: distinct from file sha
    # because the MAT container adds metadata; historical docs conflated
    # these.  Do not inherit the mislabeling.)
    from scipy.io import loadmat
    import numpy as np
    for mat_file, expected_bytes in EXPECTED_PARAMS_BYTES_SHA256.items():
        mat_path = os.path.join(HERE, mat_file)
        d = loadmat(mat_path)
        params = np.asarray(d['params']).flatten()
        bytes_sha = hashlib.sha256(params.tobytes()).hexdigest()
        if bytes_sha != expected_bytes:
            raise RuntimeError(
                f"{mat_file}: params bytes SHA {bytes_sha} != expected "
                f"{expected_bytes}; refusing to launch")
        actual[f"{mat_file}.params_bytes"] = bytes_sha
    return actual


def _nice_value():
    try:
        current_nice = os.nice(0)
    except Exception:
        current_nice = 0
    return current_nice


def _hard_requirements():
    if not __debug__:
        raise RuntimeError("__debug__ required (python -O refused)")
    if os.environ.get(RELEASE_ENV) != "1":
        raise RuntimeError(
            f"{RELEASE_ENV}=1 required for Main's released diagnostic")
    for env in REQUIRED_ENVS:
        if os.environ.get(env) != "1":
            raise RuntimeError(f"{env} must be '1' (got {os.environ.get(env)!r})")
    nice = _nice_value()
    if nice < 10:
        raise RuntimeError(
            "current process nice must be >= 10 (got "
            f"{nice}); use nice -n 10 in the launch command")
    print(f"perimeter OK: 5 envs=1, nice={nice}, __debug__=True")


def _resolve_local(path, expected_name):
    """Require one exact direct campaign output and refuse prior evidence."""
    path = os.path.abspath(path)
    expected = os.path.join(os.path.abspath(HERE), expected_name)
    if path != expected:
        raise RuntimeError(
            f"output path must be exactly {expected!r}, got {path!r}")
    for reserved in (path, path + ".tmp"):
        if os.path.lexists(reserved):
            raise RuntimeError(
                f"output path {reserved} already exists (write-once)")
    return path


def _atomic_write(path, text):
    """Create `path` atomically without overwriting prior evidence."""
    path = os.path.abspath(path)
    if os.path.lexists(path):
        raise RuntimeError(f"output path {path} already exists (write-once)")
    tmp = path + ".tmp"
    if os.path.lexists(tmp):
        raise RuntimeError(f"temporary path {tmp} already exists")
    descriptor = os.open(
        tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    if os.path.lexists(path):
        raise RuntimeError(
            f"target appeared during write: {path}; preserving temp")
    os.replace(tmp, path)
    directory_fd = os.open(os.path.dirname(path), os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)

def _parse_trailing_object(text):
    """Parse the final indented JSON object from child stdout."""
    decoder = json.JSONDecoder()
    starts = (i for i, char in enumerate(text) if char == "{")
    for start in reversed(list(starts)):
        try:
            candidate, end = decoder.raw_decode(text, start)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, dict) and not text[end:].strip():
            return candidate
    raise RuntimeError("no trailing JSON object on stdout")


def _run_child(script, extra_argv, reserved_output, stderr_log, stdout_log,
               alarm=900):
    """Run one frozen child and preserve its raw and parsed output."""
    env = dict(os.environ)
    env.update(PINNED_ENV)
    cmd = [sys.executable, os.path.join(HERE, script), *extra_argv]
    t0 = time.time()
    timed_out = False
    with open(stderr_log, "x") as stderr_fh:
        process = subprocess.Popen(
            cmd, cwd=HERE, env=env, stdout=subprocess.PIPE,
            stderr=stderr_fh)
        try:
            raw_stdout, _ = process.communicate(timeout=alarm)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.terminate()
            try:
                raw_stdout, _ = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                raw_stdout, _ = process.communicate()
        stderr_fh.flush()
        os.fsync(stderr_fh.fileno())
    directory_fd = os.open(os.path.dirname(stderr_log), os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    elapsed = time.time() - t0
    text = raw_stdout.decode("utf-8", errors="strict")
    _atomic_write(stdout_log, text)
    if timed_out:
        raise RuntimeError(f"{script} exceeded {alarm}s; TERMINATED")
    if process.returncode != 0:
        raise RuntimeError(
            f"{script} exited {process.returncode}; see {stderr_log}")

    parsed = _parse_trailing_object(text)
    _atomic_write(
        reserved_output,
        json.dumps(parsed, indent=1, sort_keys=True) + "\n")
    return parsed, elapsed, stderr_log, stdout_log


def _compare(out1, out2_with, base_target=BASE_TARGET):
    """Compare rung1 vs rung2_with after both completed."""
    key1 = "omega_cert_upper"      # rung1 key (alman25)
    key2 = "omega_with_absorbed_defects"  # rung2 with-Lemma key
    o1 = out1[key1]
    o2 = out2_with[key2]
    rung1_less = o1 < o2
    both_below = o1 < base_target and o2 < base_target
    return {"rung1_key": key1,
            "rung2_with_key": key2,
            "rung1_value": o1,
            "rung2_with_value": o2,
            "rung1_less_than_rung2": rung1_less,
            "both_below_target": both_below}


def main():
    _hard_requirements()
    if len(sys.argv) != 5:
        print("usage: launcher.py <rung1_out> <rung2_out> "
              "<rung2_without_out> <compare_out>")
        sys.exit(2)
    out1, out2_with, out2_without, out_compare = sys.argv[1:5]


    frozen = verify_and_lock_hashes()
    print(f"frozen dep hashes: "
          f"{json.dumps({k: v[:8] for k, v in frozen.items()}, indent=1)}")

    # Preflight: resolve ALL outputs + stderr + stdout logs + comparison
    # before starting rung1 (a pre-existing later path cannot leave an
    # avoidable partial run).
    resolved_paths = {}
    for key, path in [
            ("rung1_output", out1),
            ("rung2_with_output", out2_with),
            ("rung2_without_output", out2_without),
            ("comparison_summary", out_compare),
            ("rung1_stderr_log", os.path.join(HERE, "rung1_stderr.log")),
            ("rung1_stdout_log", os.path.join(HERE, "rung1_stdout.json")),
            ("rung2_with_stderr_log", os.path.join(HERE, "rung2_stderr.log")),
            ("rung2_with_stdout_log", os.path.join(HERE, "rung2_stdout.json")),
            ("rung2_without_stderr_log",
             os.path.join(HERE, "rung2_without_stderr.log")),
            ("rung2_without_stdout_log",
             os.path.join(HERE, "rung2_without_stdout.json"))]:
        resolved_paths[key] = _resolve_local(
            path, EXPECTED_OUTPUT_NAMES[key])

    r1, t1, err1, raw1 = _run_child(
        "stage_b_rung1_certified.py", [],
        resolved_paths["rung1_output"],
        resolved_paths["rung1_stderr_log"],
        resolved_paths["rung1_stdout_log"], alarm=900)
    r2_with, t2, err2, raw2 = _run_child(
        "stage_b_rung2.py", [],
        resolved_paths["rung2_with_output"],
        resolved_paths["rung2_with_stderr_log"],
        resolved_paths["rung2_with_stdout_log"], alarm=900)
    r2_without, t3, err3, raw3 = _run_child(
        "stage_b_rung2.py", ["without"],
        resolved_paths["rung2_without_output"],
        resolved_paths["rung2_without_stderr_log"],
        resolved_paths["rung2_without_stdout_log"], alarm=900)

    comparison = _compare(r1, r2_with)
    summary = {
        "rung1_elapsed": t1, "rung2_with_elapsed": t2,
        "rung2_without_elapsed": t3,
        "rung1_stderr_log": err1, "rung1_stdout_log": raw1,
        "rung2_stderr_log": err2, "rung2_stdout_log": raw2,
        "rung2_without_stderr_log": err3,
        "rung2_without_stdout_log": raw3,
        "rung2_diagnostic_off_recorded": r2_without,
        "comparison": comparison,
        "frozen_dep_hashes": frozen,
    }
    _atomic_write(
        resolved_paths["comparison_summary"],
        json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
