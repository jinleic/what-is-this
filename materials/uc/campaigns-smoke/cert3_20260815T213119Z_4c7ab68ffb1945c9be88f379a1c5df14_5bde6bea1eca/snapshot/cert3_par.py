"""Create and run an immutable, partitioned cert3 certification campaign.

The campaign protocol separates mutable development from proof execution:

1. ``create`` makes a fresh UUID-named directory, copies every executable and
   reviewed proof/checker source into ``snapshot/``, and atomically writes a
   launch manifest containing their SHA-256 hashes, exact parameters, eight
   exact dyadic slice intervals, and one unique run UUID per slice.
2. ``run`` must be invoked from the copied ``snapshot/cert3_par.py``.  A parent
   process supervises one worker, records its real exit code, and commits a
   run-unique result with fsync + atomic rename only after exit code 0.
3. ``cert3_collect.py`` accepts only the eight exact result filenames and run
   UUIDs named by that fresh launch manifest.

No file is overwritten or deleted.  A killed run can leave a hidden worker
artifact, but can never leave a committed result.  Make a new campaign rather
than reusing a run UUID.

Usage:
  cert3_par.py create OFFSET SECONDS [CAMPAIGN_ROOT]
  cert3_par.py run LAUNCH_JSON SLICE

Example:
  python uc/cert3_par.py create 0.0001 28800 uc/campaigns
  python uc/campaigns/<id>/snapshot/cert3_par.py run \
      uc/campaigns/<id>/launch.json 0
"""

import datetime
import hashlib
import json
import os
import math
import platform
import subprocess
import sys
import uuid

# Pinned-source execution: never materialize bytecode caches inside a hashed
# snapshot.  Import order (extension modules, then cached bytecode) could
# otherwise execute bytes the manifest never hashed.
sys.dont_write_bytecode = True
from fractions import Fraction

import flint
import mpmath
import numpy
import scipy
import sympy
from flint import ctx
from mpmath import mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cert2
import cert3
from arbcore import get_rh_gmax
from entropy import PSI

NSLICES = 8
MIN_WIDTH = 1e-3
WORK_PREC = 80
COLLAR_MIN_WIDTH = 1.25e-4
FACE_MIN_WIDTH = 1e-3
FACE_BOX_BUDGET = 100000
MAX_BOXES = 2_000_000_000
PROGRESS_EVERY = 4_000_000

# Every local Python module whose bytes can change a numerical verdict.
EXECUTABLE_FILES = (
    "arbcore.py",
    "bound_kkt.py",
    "cert2.py",
    "cert3.py",
    "cert3_par.py",
    "diag_exhaust.py",
    "entropy.py",
)

# These do not execute inside the certifier.  They are copied and hashed
# separately so a reviewer can pin the independent checker and every local
# artifact in the analytic theorem chain, including Theorem A/A' dependencies
# and the exact bridge to Cambie's cited implication.
REVIEW_FILES = (
    "bridge_uc.py",
    "cert3_collect.py",
    "cert3_replay.py",
    "decomposition.py",
    "lemma_rh_proof.py",
    "margin_lemma.py",
    "reduction.py",
    "thmB3_proof.py",
)


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("ascii")


def object_sha256(value, omitted_key=None):
    if omitted_key is not None:
        value = {k: v for k, v in value.items() if k != omitted_key}
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def hash_files(directory, names):
    return {name: file_sha256(os.path.join(directory, name))
            for name in sorted(names)}


def combined_source_hash(per_file):
    digest = hashlib.sha256()
    for name in sorted(per_file):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(per_file[name].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def environment_record():
    return {
        "flint": getattr(flint, "__version__", "unknown"),
        "mpmath": mpmath.__version__,
        "numpy": numpy.__version__,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "scipy": scipy.__version__,
        "sympy": sympy.__version__,
    }


def fraction_text(value):
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def expected_interval(index):
    lo = Fraction(1, 2) + Fraction(index, 2 * NSLICES)
    hi = Fraction(1, 2) + Fraction(index + 1, 2 * NSLICES)
    return fraction_text(lo), fraction_text(hi)


def valid_run_id(value):
    if type(value) is not str or len(value) != 32:
        return False
    try:
        parsed = uuid.UUID(hex=value)
    except (ValueError, AttributeError):
        return False
    return parsed.hex == value and parsed.version == 4


def campaign_artifact_path(campaign_dir, filename, label):
    """Return a direct child of ``campaign_dir`` or reject it."""
    if type(filename) is not str or os.path.basename(filename) != filename:
        raise RuntimeError("%s must be a basename" % label)
    campaign_dir = os.path.abspath(campaign_dir)
    path = os.path.abspath(os.path.join(campaign_dir, filename))
    try:
        confined = os.path.commonpath((campaign_dir, path)) == campaign_dir
    except ValueError:
        confined = False
    if not confined or os.path.dirname(path) != campaign_dir:
        raise RuntimeError("%s escapes campaign directory" % label)
    return path


def validate_campaign_path(campaign_dir, path, label):
    """Validate a supervisor-provided direct-child path before worker use."""
    campaign_dir = os.path.abspath(campaign_dir)
    path = os.path.abspath(path)
    try:
        confined = os.path.commonpath((campaign_dir, path)) == campaign_dir
    except ValueError:
        confined = False
    if not confined or os.path.dirname(path) != campaign_dir:
        raise RuntimeError("%s escapes campaign directory" % label)
    return path


def target_relation_holds(t_decimal, offset):
    """Exactly prove t_decimal >= (3-sqrt(5))/2 + offset.

    For r=t-offset in [0,3/2], r is at least the smaller root psi iff
    r^2-3r+1 <= 0.  All quantities here are exact rationals.
    """
    try:
        delta = Fraction(offset)
        r = Fraction(t_decimal) - delta
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return (delta > 0 and 0 <= r <= Fraction(3, 2)
            and r * r - 3 * r + 1 <= 0)


def certified_t_decimal(offset):
    """Choose a short decimal at or above the exact algebraic target."""
    delta = Fraction(offset)
    if delta <= 0:
        raise ValueError("OFFSET must be a positive rational decimal")
    candidate = float(mp.mpf(PSI) + mp.mpf(offset))
    for _ in range(16):
        text = repr(candidate)
        if target_relation_holds(text, offset):
            return text
        candidate = math.nextafter(candidate, math.inf)
    raise RuntimeError("failed to bracket exact target from above")


def runtime_parameters(offset):
    t_decimal = certified_t_decimal(offset)
    t = mp.mpf(t_decimal)
    old_prec = ctx.prec
    try:
        ctx.prec = WORK_PREC
        lambdas, _ = cert3.lambda_family(t)
    finally:
        ctx.prec = old_prec
    return {
        "alpha": "0.0356069",
        "center_gate": cert3.CENTER_GATE,
        "centered5_max_width": cert3.CENTERED5_MAX_WIDTH,
        "collar_edge": cert3.COLLAR_EDGE,
        "collar_min_width": COLLAR_MIN_WIDTH,
        "face_min_width": FACE_MIN_WIDTH,
        "face_box_budget": FACE_BOX_BUDGET,
        "lambda_family": [repr(float(value)) for value in lambdas],
        "t_decimal": t_decimal,
        "level_offset": offset,
        "max_boxes": MAX_BOXES,
        "min_width": MIN_WIDTH,
        "nslices": NSLICES,
        "progress_every": PROGRESS_EVERY,
        "q2pin": cert3.Q2PIN,
        "root_w_hi": "1",
        "root_w_lo": "1/2",
        "work_prec_bits": WORK_PREC,
    }


def fsync_directory(path):
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    fd = os.open(path, flags)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def install_no_replace(temp, path):
    """Atomically install ``temp`` at ``path`` without replacing a peer."""
    import ctypes
    libc = ctypes.CDLL(None, use_errno=True)
    rename_exclusive = libc.renamex_np
    rename_exclusive.argtypes = (ctypes.c_char_p, ctypes.c_char_p,
                                 ctypes.c_uint)
    rename_exclusive.restype = ctypes.c_int
    # Darwin <stdio.h>: RENAME_EXCL = 0x00000004.
    if rename_exclusive(os.fsencode(temp), os.fsencode(path), 0x00000004) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), path)


def write_new_bytes(path, data):
    """Create ``path`` through a durable, atomically no-replace rename.

    A check followed by ``os.rename`` is racy because POSIX rename replaces an
    existing destination.  ``RENAME_EXCL`` makes the no-overwrite invariant a
    single kernel operation even under concurrent same-slice supervisors.
    """
    path = os.path.abspath(path)
    parent = os.path.dirname(path)
    os.makedirs(parent, exist_ok=True)
    temp = os.path.join(parent, ".tmp-%s-%s" %
                        (os.path.basename(path), uuid.uuid4().hex))
    with open(temp, "xb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    try:
        install_no_replace(temp, path)
    except FileExistsError:
        # Preserve the run-unique hidden temp as evidence of the rejected
        # duplicate attempt; no file is deleted or overwritten.
        raise RuntimeError("refusing to overwrite existing artifact: %s; "
                           "rejected temp retained at %s" % (path, temp))
    fsync_directory(parent)


def write_new_json(path, value):
    write_new_bytes(path, canonical_bytes(value) + b"\n")


def validate_launch_structure(launch):
    problems = []
    if type(launch) is not dict:
        return ["launch document must be a JSON object"]
    if launch.get("schema") != "cert3-campaign-v2":
        problems.append("unsupported launch schema")
    if type(launch.get("nslices")) is not int \
            or launch.get("nslices") != NSLICES:
        problems.append("nslices must equal integer 8")
    if launch.get("launch_sha256") != object_sha256(launch, "launch_sha256"):
        problems.append("launch_sha256 mismatch")
    t_decimal = launch.get("t_decimal")
    offset = launch.get("level_offset", "")
    if not target_relation_holds(t_decimal, offset):
        problems.append("t_decimal does not dominate exact psi+offset")
    params = launch.get("parameters")
    if type(params) is not dict:
        problems.append("parameters must be a JSON object")
    else:
        if params.get("t_decimal") != t_decimal:
            problems.append("parameter/launch t_decimal mismatch")
        for key, expected in (
                ("nslices", NSLICES),
                ("max_boxes", MAX_BOXES),
                ("work_prec_bits", WORK_PREC),
                ("face_box_budget", FACE_BOX_BUDGET)):
            if type(params.get(key)) is not int or params.get(key) != expected:
                problems.append(
                    "parameters.%s must equal integer %d" % (key, expected))
    runs = launch.get("runs")
    if type(runs) is not list or len(runs) != NSLICES:
        problems.append("launch must contain exactly 8 runs")
        return problems
    if any(type(run) is not dict for run in runs):
        problems.append("every run must be a JSON object")
        return problems

    indices = [run.get("slice") for run in runs]
    if any(type(idx) is not int for idx in indices):
        problems.append("every run slice must be an integer")
    elif sorted(indices) != list(range(NSLICES)):
        problems.append("run indices must be unique and exactly 0..7")

    run_ids = [run.get("run_id") for run in runs]
    if any(not valid_run_id(run_id) for run_id in run_ids):
        problems.append("run IDs must be canonical UUID4 lowercase hex")
    elif len(set(run_ids)) != NSLICES:
        problems.append("run IDs must be eight unique values")

    result_names = [run.get("result_filename") for run in runs]
    if any(type(name) is not str or os.path.basename(name) != name
           for name in result_names):
        problems.append("result filenames must be basenames")
    elif len(set(result_names)) != NSLICES:
        problems.append("result filenames must be unique")
    trace_names = [run.get("trace_filename") for run in runs]
    if any(type(name) is not str or os.path.basename(name) != name
           for name in trace_names):
        problems.append("trace filenames must be basenames")
    elif len(set(trace_names)) != NSLICES:
        problems.append("trace filenames must be unique")

    for position, run in enumerate(runs):
        idx = run.get("slice")
        run_id = run.get("run_id")
        if type(idx) is int and 0 <= idx < NSLICES:
            lo, hi = expected_interval(idx)
            if run.get("w_lo") != lo or run.get("w_hi") != hi:
                problems.append("slice %d has wrong exact interval" % idx)
            if valid_run_id(run_id):
                expected_result = "result_slice%d_%s.json" % (idx, run_id)
                if run.get("result_filename") != expected_result:
                    problems.append(
                        "slice %d has wrong result filename" % idx)
                expected_trace = "trace_slice%d_%s.bin" % (idx, run_id)
                if run.get("trace_filename") != expected_trace:
                    problems.append(
                        "slice %d has wrong trace filename" % idx)
        elif type(idx) is int:
            problems.append("run %d slice is outside 0..7" % position)
        if type(run.get("time_budget_s")) is not int \
                or run.get("time_budget_s") <= 0:
            problems.append("slice %r has invalid time budget" % idx)
    return problems


def verify_launch(launch_path, require_current_sources=True):
    with open(launch_path) as fh:
        launch = json.load(fh)
    problems = validate_launch_structure(launch)
    if type(launch) is not dict:
        raise RuntimeError("invalid launch manifest:\n  - " +
                           "\n  - ".join(problems))
    snapshot = os.path.join(os.path.dirname(os.path.abspath(launch_path)),
                            "snapshot")
    if os.path.realpath(HERE) != os.path.realpath(snapshot):
        problems.append("driver is not executing from this campaign snapshot")
    expected_entries = sorted(set(EXECUTABLE_FILES) | set(REVIEW_FILES))
    try:
        entries = sorted(os.listdir(snapshot))
    except OSError as exc:
        problems.append("cannot list snapshot directory: %s" % exc)
        entries = expected_entries
    if entries != expected_entries:
        problems.append(
            "snapshot must contain exactly the pinned sources; extra=%s "
            "missing=%s" % (sorted(set(entries) - set(expected_entries)),
                            sorted(set(expected_entries) - set(entries))))
    if require_current_sources:
        executable = hash_files(snapshot, EXECUTABLE_FILES)
        reviewed = hash_files(snapshot, REVIEW_FILES)
        if executable != launch.get("executable_files"):
            problems.append("snapshot executable hashes differ from launch")
        if combined_source_hash(executable) != launch.get("code_sha256"):
            problems.append("snapshot code_sha256 differs from launch")
        if reviewed != launch.get("review_files"):
            problems.append("snapshot review-file hashes differ from launch")
        if environment_record() != launch.get("environment"):
            problems.append("runtime environment differs from launch")
        try:
            expected_params = runtime_parameters(
                launch.get("level_offset", ""))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
            problems.append("cannot derive runtime parameters: %s" % exc)
        else:
            if expected_params != launch.get("parameters"):
                problems.append("runtime parameters differ from launch")
    if problems:
        raise RuntimeError("invalid launch manifest:\n  - " +
                           "\n  - ".join(problems))
    return launch


def create_campaign(offset, seconds, campaign_root):
    if seconds <= 0 or int(seconds) != seconds:
        raise ValueError("SECONDS must be a positive integer")
    seconds = int(seconds)
    executable = hash_files(HERE, EXECUTABLE_FILES)
    reviewed = hash_files(HERE, REVIEW_FILES)
    code_hash = combined_source_hash(executable)
    now = datetime.datetime.now(datetime.timezone.utc)
    campaign_id = "%s_%s" % (now.strftime("%Y%m%dT%H%M%SZ"), uuid.uuid4().hex)
    campaign_dir = os.path.abspath(os.path.join(
        campaign_root, "cert3_%s_%s" % (campaign_id, code_hash[:12])))
    snapshot = os.path.join(campaign_dir, "snapshot")
    os.makedirs(snapshot, exist_ok=False)
    for name in sorted(set(EXECUTABLE_FILES + REVIEW_FILES)):
        source = os.path.join(HERE, name)
        target = os.path.join(snapshot, name)
        with open(source, "rb") as fh:
            write_new_bytes(target, fh.read())
    copied_exec = hash_files(snapshot, EXECUTABLE_FILES)
    copied_review = hash_files(snapshot, REVIEW_FILES)
    if copied_exec != executable or copied_review != reviewed:
        raise RuntimeError("snapshot copy hash mismatch")
    runs = []
    for idx in range(NSLICES):
        run_id = uuid.uuid4().hex
        lo, hi = expected_interval(idx)
        runs.append({
            "result_filename": "result_slice%d_%s.json" % (idx, run_id),
            "trace_filename": "trace_slice%d_%s.bin" % (idx, run_id),
            "run_id": run_id,
            "slice": idx,
            "time_budget_s": seconds,
            "w_hi": hi,
            "w_lo": lo,
        })
    launch = {
        "campaign_id": campaign_id,
        "code_sha256": code_hash,
        "created_utc": now.isoformat(),
        "environment": environment_record(),
        "executable_files": executable,
        "level_offset": offset,
        "t_decimal": runtime_parameters(offset)["t_decimal"],
        "nslices": NSLICES,
        "parameters": runtime_parameters(offset),
        "review_files": reviewed,
        "runs": runs,
        "schema": "cert3-campaign-v2",
    }
    launch["launch_sha256"] = object_sha256(launch)
    launch_path = os.path.join(campaign_dir, "launch.json")
    write_new_json(launch_path, launch)
    fsync_directory(campaign_dir)
    print("CAMPAIGN_PATH %s" % campaign_dir, flush=True)
    print("LAUNCH_PATH %s" % launch_path, flush=True)
    print("CODE_SHA256 %s" % code_hash, flush=True)
    print("LAUNCH_SHA256 %s" % launch["launch_sha256"], flush=True)
    print("COLLECTOR_SHA256 %s" % reviewed["cert3_collect.py"], flush=True)
    driver = os.path.join(snapshot, "cert3_par.py")
    for idx in range(NSLICES):
        print("RUN %d %s %s run %s %d" %
              (idx, sys.executable, driver, launch_path, idx), flush=True)
    return launch_path


def run_record_common(launch, run):
    return {
        "campaign_id": launch["campaign_id"],
        "code_sha256": launch["code_sha256"],
        "launch_sha256": launch["launch_sha256"],
        "level_offset": launch["level_offset"],
        "nslices": launch["nslices"],
        "result_filename": run["result_filename"],
        "run_id": run["run_id"],
        "slice": run["slice"],
        "t_decimal": launch["t_decimal"],
        "time_budget_s": run["time_budget_s"],
        "trace_filename": run["trace_filename"],
        "w_hi": run["w_hi"],
        "w_lo": run["w_lo"],
    }


def worker_main(launch_path, slice_index, worker_path, trace_temp_path):
    launch = verify_launch(launch_path)
    run = launch["runs"][slice_index]
    if type(run["slice"]) is not int or run["slice"] != slice_index:
        raise RuntimeError("slice/run mismatch")
    campaign_dir = os.path.dirname(os.path.abspath(launch_path))
    worker_path = validate_campaign_path(
        campaign_dir, worker_path, "worker record path")
    trace_temp_path = validate_campaign_path(
        campaign_dir, trace_temp_path, "trace temp path")
    t_text = launch["t_decimal"]
    t = mp.mpf(t_text)
    # Preserve the validated order used by the certifier: compute the global
    # covers at their modules' native precision, then set working precision.
    gmax = get_rh_gmax()
    rho_gmax = cert2.get_rho_gmax()
    ctx.prec = WORK_PREC
    lambdas, _ = cert3.lambda_family(t)
    if [repr(float(value)) for value in lambdas] \
            != launch["parameters"]["lambda_family"]:
        raise RuntimeError("lambda family changed after launch")
    lo = float(Fraction(run["w_lo"]))
    hi = float(Fraction(run["w_hi"]))
    root = ((0.0, 1.0),) * 4 + ((lo, hi),)
    residual_name = "residual_slice%d_%s_%s.npy" % (
        slice_index, run["run_id"], uuid.uuid4().hex)
    residual_dump = campaign_artifact_path(
        campaign_dir, residual_name, "residual dump path")
    print("WORKER_START slice=%d/8 run_id=%s w=[%s,%s] t=psi+%s "
          "budget=%ds code_sha256=%s" %
          (slice_index, run["run_id"], run["w_lo"], run["w_hi"],
           launch["level_offset"], run["time_budget_s"],
           launch["code_sha256"]), flush=True)
    with open(trace_temp_path, "xb") as trace_fh:
        try:
            complete, stats = cert3.certify(
                t_text,
                root,
                max_boxes=MAX_BOXES,
                time_budget=float(run["time_budget_s"]),
                min_width=MIN_WIDTH,
                face_min_width=FACE_MIN_WIDTH,
                face_box_budget=FACE_BOX_BUDGET,
                collar_min_width=COLLAR_MIN_WIDTH,
                gmax=gmax,
                rho_gmax=rho_gmax,
                lambdas=lambdas,
                progress_every=PROGRESS_EVERY,
                residual_dump=residual_dump,
                trace=trace_fh,
            )
        finally:
            trace_fh.flush()
            os.fsync(trace_fh.fileno())
    fsync_directory(campaign_dir)
    # Detect any mid-run edit to the immutable snapshot before recording.
    verify_launch(launch_path)
    tallies = dict(stats)
    # Counter returns zero for an absent key without inserting it.  Materialize
    # every completion-critical tally so the independent checker can require
    # explicit integer zeros rather than treating absence as success.
    for key in ("stack", "residual", "budget_boxes", "budget_time"):
        tallies[key] = int(stats[key])
    processed = tallies.get("processed")
    trace_size = os.path.getsize(trace_temp_path)
    if type(processed) is not int or trace_size != processed:
        raise RuntimeError(
            "trace size %d does not equal integer processed tally %r" %
            (trace_size, processed))
    record = run_record_common(launch, run)
    record.update({
        "complete": bool(complete),
        "returned_normally": True,
        "t_decimal": t_text,
        "tallies": tallies,
        "trace_schema": "cert3-trace-v1",
        "trace_sha256": file_sha256(trace_temp_path),
        "trace_size": trace_size,
        "verdict": "COMPLETE" if complete else "INCOMPLETE",
        "worker_schema": "cert3-worker-v2",
    })
    record["worker_record_sha256"] = object_sha256(
        record, "worker_record_sha256")
    write_new_json(worker_path, record)
    print("WORKER_RECORD %s" % worker_path, flush=True)
    return 0


def validate_worker_record(record, launch, run):
    problems = []
    if type(record) is not dict:
        raise RuntimeError("invalid worker record:\n  - record is not an object")
    if record.get("worker_schema") != "cert3-worker-v2":
        problems.append("wrong worker schema")
    if record.get("worker_record_sha256") != object_sha256(
            record, "worker_record_sha256"):
        problems.append("worker record SHA-256 mismatch")
    for key, value in run_record_common(launch, run).items():
        if record.get(key) != value:
            problems.append("worker field %s differs from launch" % key)
    if record.get("returned_normally") is not True:
        problems.append("worker did not record normal return")
    if record.get("trace_schema") != "cert3-trace-v1":
        problems.append("wrong trace schema")
    trace_hash = record.get("trace_sha256")
    if type(trace_hash) is not str or len(trace_hash) != 64 \
            or any(char not in "0123456789abcdef" for char in trace_hash):
        problems.append("trace SHA-256 is not lowercase hexadecimal")
    trace_size = record.get("trace_size")
    if type(trace_size) is not int or trace_size < 0:
        problems.append("trace size is not a nonnegative integer")
    tallies = record.get("tallies")
    if type(tallies) is not dict:
        problems.append("worker tallies are missing/non-object")
    else:
        for key, value in tallies.items():
            if key != "residual_geometry" and type(value) is not int:
                problems.append("tally %s is not an integer" % key)
        if type(tallies.get("processed")) is not int:
            problems.append("processed tally is not an integer")
        elif type(trace_size) is int and trace_size != tallies["processed"]:
            problems.append("trace size differs from processed tally")
    if problems:
        raise RuntimeError("invalid worker record:\n  - " +
                           "\n  - ".join(problems))


def supervisor_main(launch_path, slice_index):
    launch_path = os.path.abspath(launch_path)
    launch = verify_launch(launch_path)
    if type(slice_index) is not int or not 0 <= slice_index < NSLICES:
        raise ValueError("SLICE must be an integer in 0..7")
    run = launch["runs"][slice_index]
    campaign_dir = os.path.dirname(launch_path)
    result_path = campaign_artifact_path(
        campaign_dir, run["result_filename"], "result path")
    trace_path = campaign_artifact_path(
        campaign_dir, run["trace_filename"], "trace path")
    if os.path.exists(result_path):
        raise RuntimeError("refusing to overwrite committed result: %s" %
                           result_path)
    if os.path.exists(trace_path):
        raise RuntimeError("refusing to overwrite committed trace: %s" %
                           trace_path)
    worker_name = ".worker_slice%d_%s_%s.json" % (
        slice_index, run["run_id"], uuid.uuid4().hex)
    trace_temp_name = ".trace_slice%d_%s_%s.bin" % (
        slice_index, run["run_id"], uuid.uuid4().hex)
    worker_path = campaign_artifact_path(
        campaign_dir, worker_name, "worker record path")
    trace_temp_path = campaign_artifact_path(
        campaign_dir, trace_temp_name, "trace temp path")
    command = [sys.executable, "-B", os.path.abspath(__file__), "_worker",
               launch_path, str(slice_index), worker_path, trace_temp_path]
    print("SUPERVISOR_START " + " ".join(command), flush=True)
    process = subprocess.Popen(command)
    returncode = process.wait()
    print("WORKER_EXIT_CODE %d" % returncode, flush=True)
    if returncode != 0:
        print("NO_RESULT_COMMITTED worker failed", flush=True)
        return returncode if returncode > 0 else 128 - returncode
    # Verify the environment and frozen bytes again after child exit.
    launch = verify_launch(launch_path)
    with open(worker_path) as fh:
        worker_record = json.load(fh)
    validate_worker_record(worker_record, launch, run)
    trace_size = os.path.getsize(trace_temp_path)
    if trace_size != worker_record["trace_size"]:
        raise RuntimeError("trace temp size differs from worker record")
    if trace_size != worker_record["tallies"]["processed"]:
        raise RuntimeError("trace temp size differs from processed tally")
    if file_sha256(trace_temp_path) != worker_record["trace_sha256"]:
        raise RuntimeError("trace temp SHA-256 differs from worker record")
    try:
        install_no_replace(trace_temp_path, trace_path)
    except FileExistsError:
        raise RuntimeError("refusing to overwrite committed trace: %s; "
                           "rejected temp retained at %s" %
                           (trace_path, trace_temp_path))
    fsync_directory(campaign_dir)
    result = dict(worker_record)
    result.update({
        "result_schema": "cert3-result-v2",
        "supervisor_committed": True,
        "worker_exit_code": returncode,
    })
    result["record_sha256"] = object_sha256(result, "record_sha256")
    write_new_json(result_path, result)
    print("TRACE_PATH %s" % trace_path, flush=True)
    print("RESULT_PATH %s" % result_path, flush=True)
    print("RESULT_JSON " + canonical_bytes(result).decode("ascii"), flush=True)
    return 0


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    mode = argv[0]
    if mode == "create":
        if len(argv) not in (3, 4):
            raise SystemExit(__doc__)
        root = argv[3] if len(argv) == 4 else os.path.join(HERE, "campaigns")
        create_campaign(argv[1], float(argv[2]), root)
        return 0
    if mode == "run":
        if len(argv) != 3:
            raise SystemExit(__doc__)
        return supervisor_main(argv[1], int(argv[2]))
    if mode == "_worker":
        if len(argv) != 5:
            raise SystemExit("internal worker arguments invalid")
        return worker_main(argv[1], int(argv[2]), argv[3], argv[4])
    raise SystemExit("unknown mode %r\n%s" % (mode, __doc__))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
