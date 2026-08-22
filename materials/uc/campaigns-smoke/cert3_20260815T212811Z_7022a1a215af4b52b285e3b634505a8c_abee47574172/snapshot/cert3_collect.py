"""Independently audit one immutable cert3 campaign.

The checker accepts a campaign only by re-proving it.  JSON provenance fields
(hashes, exit codes, verdicts, tallies) gate admission, but the mathematical
acceptance criterion is the replay of all eight committed proof traces: every
trace is re-executed from its exact dyadic root with the frozen snapshot
modules, each claimed discharge rule is recomputed in interval arithmetic, and
the DFS must end with an empty stack, no residual event, and tallies equal to
the committed record.  A hand-authored record without a mathematically valid
trace is rejected; a valid trace without matching provenance is also rejected.

Usage:
  python <campaign>/snapshot/cert3_collect.py <campaign>/launch.json

It claims a composite certificate only when the fresh launch manifest names
exactly eight unique runs over the exact dyadic partition of [1/2,1], every
run has exactly one atomic, hash-valid result with a real worker exit code 0,
normal-return markers, COMPLETE verdict, positive work,
stack=residual=budget_boxes=budget_time=0, one committed proof trace whose
bytes re-hash to the recorded digest and whose size equals the processed
tally, and every one of the eight traces REPLAYS: the frozen
``cert3_replay.replay_trace`` re-derives each box and re-proves each clearing
step.  Missing, duplicate, unexpected, stale, malformed, mismatched,
incomplete, or non-replaying artifacts all block the claim.
"""

import hashlib
import json
import os
import sys
import time
import uuid
from fractions import Fraction

NSLICES = 8
ZERO_TALLIES = ("stack", "residual", "budget_boxes", "budget_time")

# The checker's own authoritative inventory of load-bearing artifacts.  The
# audited manifest must name exactly these files; taking the name sets from
# the manifest itself would prove only "whatever the launch listed is
# unmodified", never that the inventory is complete.
EXPECTED_EXECUTABLE_FILES = frozenset((
    "arbcore.py",
    "bound_kkt.py",
    "cert2.py",
    "cert3.py",
    "cert3_par.py",
    "diag_exhaust.py",
    "entropy.py",
))
EXPECTED_REVIEW_FILES = frozenset((
    "bridge_uc.py",
    "cert3_collect.py",
    "cert3_replay.py",
    "decomposition.py",
    "lemma_rh_proof.py",
    "margin_lemma.py",
    "reduction.py",
    "thmB3_proof.py",
))
REPLAY_TALLY_KEYS = (
    "processed", "infeasible", "corner", "ratio", "center", "center_mixed",
    "center_mixed_swap", "center_w", "face", "residual", "split")


def valid_run_id(value):
    """Exactly the producer's rule: canonical lowercase-hex UUID4."""
    if type(value) is not str or len(value) != 32:
        return False
    try:
        parsed = uuid.UUID(hex=value)
    except (ValueError, AttributeError):
        return False
    return parsed.hex == value and parsed.version == 4


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


def fraction_text(value):
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def expected_interval(index):
    return (fraction_text(Fraction(1, 2) + Fraction(index, 2 * NSLICES)),
            fraction_text(Fraction(1, 2) +
                          Fraction(index + 1, 2 * NSLICES)))


def target_relation_holds(t_decimal, offset):
    """Exact rational check: t_decimal >= (3-sqrt(5))/2 + offset."""
    try:
        delta = Fraction(offset)
        r = Fraction(t_decimal) - delta
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return (delta > 0 and 0 <= r <= Fraction(3, 2)
            and r * r - 3 * r + 1 <= 0)


def launch_problems(launch, launch_path):
    problems = []
    if launch.get("schema") != "cert3-campaign-v2":
        problems.append("unsupported launch schema")
    if launch.get("launch_sha256") != object_sha256(launch, "launch_sha256"):
        problems.append("launch_sha256 mismatch")
    if type(launch.get("nslices")) is not int or launch.get("nslices") != NSLICES:
        problems.append("nslices is %r, expected integer 8" % launch.get("nslices"))
    params = launch.get("parameters")
    if type(params) is not dict:
        problems.append("parameters must be a JSON object")
        params = {}
    if type(params.get("nslices")) is not int or params.get("nslices") != NSLICES:
        problems.append("parameters.nslices is not integer 8")
    if type(params.get("work_prec_bits")) is not int \
            or params.get("work_prec_bits") != 80:
        problems.append("parameters.work_prec_bits is not integer 80")
    if type(params.get("face_box_budget")) is not int \
            or params.get("face_box_budget") <= 0:
        problems.append("parameters.face_box_budget is not a positive integer")
    if params.get("level_offset") != launch.get("level_offset"):
        problems.append("parameter/launch level offsets differ")
    t_decimal = launch.get("t_decimal")
    if not target_relation_holds(t_decimal, launch.get("level_offset", "")):
        problems.append("t_decimal does not dominate exact psi+offset")
    if params.get("t_decimal") != t_decimal:
        problems.append("parameter/launch t_decimal mismatch")
    runs = launch.get("runs")
    if type(runs) is not list or len(runs) != NSLICES:
        problems.append("launch must contain exactly 8 runs")
        return problems
    if any(type(run) is not dict for run in runs):
        problems.append("every run must be a JSON object")
        return problems
    indices = [run.get("slice") for run in runs]
    indices_ok = all(type(idx) is int for idx in indices)
    if not indices_ok:
        problems.append("every run slice must be an integer")
    elif sorted(indices) != list(range(NSLICES)):
        indices_ok = False
        problems.append("slice indices must be unique and exactly 0..7")
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
        if type(idx) is int and 0 <= idx < NSLICES:
            lo, hi = expected_interval(idx)
            if run.get("w_lo") != lo or run.get("w_hi") != hi:
                problems.append("slice %d interval [%r,%r] != [%s,%s]" %
                                (idx, run.get("w_lo"), run.get("w_hi"),
                                 lo, hi))
            if valid_run_id(run.get("run_id")):
                expected_name = "result_slice%d_%s.json" % (idx, run["run_id"])
                if run.get("result_filename") != expected_name:
                    problems.append(
                        "slice %d result filename is not run-unique" % idx)
                expected_trace = "trace_slice%d_%s.bin" % (idx, run["run_id"])
                if run.get("trace_filename") != expected_trace:
                    problems.append(
                        "slice %d trace filename is not run-unique" % idx)
        elif type(idx) is int:
            problems.append("run %d slice is outside 0..7" % position)
        if type(run.get("time_budget_s")) is not int \
                or run.get("time_budget_s") <= 0:
            problems.append("slice %r time budget is not a positive integer" % idx)
    # Exact Fraction arithmetic proves no gap or overlap in the launch tiling.
    if indices_ok:
        ordered = sorted(runs, key=lambda run: run["slice"])
        try:
            intervals = [(Fraction(run["w_lo"]), Fraction(run["w_hi"]))
                         for run in ordered]
        except (TypeError, ValueError, ZeroDivisionError, KeyError) as exc:
            problems.append("slice intervals are not exact rationals: %s" % exc)
            intervals = None
        if intervals is not None:
            if intervals[0][0] != Fraction(1, 2) or intervals[-1][1] != 1:
                problems.append("slice intervals do not span exactly [1/2,1]")
            for left, right in zip(intervals, intervals[1:]):
                if left[1] != right[0]:
                    problems.append("slice gap/overlap between %s and %s" %
                                    (left, right))
    campaign_dir = os.path.dirname(os.path.abspath(launch_path))
    snapshot = os.path.join(campaign_dir, "snapshot")
    executable_expected = launch.get("executable_files")
    review_expected = launch.get("review_files")
    if not isinstance(executable_expected, dict):
        problems.append("missing executable file hashes")
    elif set(executable_expected) != EXPECTED_EXECUTABLE_FILES:
        problems.append(
            "executable inventory %s differs from the checker's pinned set %s"
            % (sorted(executable_expected), sorted(EXPECTED_EXECUTABLE_FILES)))
    else:
        try:
            executable_actual = hash_files(snapshot, executable_expected)
        except OSError as exc:
            problems.append("cannot hash executable snapshot: %s" % exc)
        else:
            if executable_actual != executable_expected:
                problems.append("executable snapshot hashes differ from launch")
            if combined_source_hash(executable_actual) != launch.get("code_sha256"):
                problems.append("snapshot code_sha256 differs from launch")
    if not isinstance(review_expected, dict):
        problems.append("missing separately pinned review file hashes")
    elif set(review_expected) != EXPECTED_REVIEW_FILES:
        problems.append(
            "review inventory %s differs from the checker's pinned set %s"
            % (sorted(review_expected), sorted(EXPECTED_REVIEW_FILES)))
    else:
        try:
            review_actual = hash_files(snapshot, review_expected)
        except OSError as exc:
            problems.append("cannot hash review snapshot: %s" % exc)
        else:
            if review_actual != review_expected:
                problems.append("review/proof snapshot hashes differ from launch")
            own_name = os.path.basename(os.path.abspath(__file__))
            own_expected = review_expected.get(own_name)
            if own_expected != file_sha256(os.path.abspath(__file__)):
                problems.append("running checker bytes differ from pinned checker")
    return problems


def result_common(launch, run):
    return {
        "campaign_id": launch["campaign_id"],
        "code_sha256": launch["code_sha256"],
        "launch_sha256": launch["launch_sha256"],
        "level_offset": launch["level_offset"],
        "t_decimal": launch["t_decimal"],
        "nslices": launch["nslices"],
        "run_id": run["run_id"],
        "slice": run["slice"],
        "time_budget_s": run["time_budget_s"],
        "w_hi": run["w_hi"],
        "w_lo": run["w_lo"],
        "result_filename": run["result_filename"],
        "trace_filename": run["trace_filename"],
    }


def result_problems(record, launch, run, filename):
    tag = "slice %d (%s)" % (run["slice"], filename)
    problems = []
    if record.get("result_schema") != "cert3-result-v2":
        problems.append("%s: wrong result schema" % tag)
    if record.get("worker_schema") != "cert3-worker-v2":
        problems.append("%s: wrong worker schema" % tag)
    if record.get("record_sha256") != object_sha256(record, "record_sha256"):
        problems.append("%s: result SHA-256 mismatch" % tag)
    worker_view = {k: v for k, v in record.items()
                   if k not in ("record_sha256", "result_schema",
                                "supervisor_committed", "worker_exit_code")}
    if worker_view.get("worker_record_sha256") != object_sha256(
            worker_view, "worker_record_sha256"):
        problems.append("%s: embedded worker SHA-256 mismatch" % tag)
    for key, value in result_common(launch, run).items():
        if record.get(key) != value:
            problems.append("%s: field %s differs from launch" % (tag, key))
    if type(record.get("worker_exit_code")) is not int \
            or record.get("worker_exit_code") != 0:
        problems.append("%s: worker exit code is not integer 0" % tag)
    if record.get("returned_normally") is not True:
        problems.append("%s: worker did not return normally" % tag)
    if record.get("supervisor_committed") is not True:
        problems.append("%s: supervisor commit marker absent" % tag)
    if record.get("complete") is not True or record.get("verdict") != "COMPLETE":
        problems.append("%s: verdict is not unambiguously COMPLETE" % tag)
    tallies = record.get("tallies")
    if not isinstance(tallies, dict):
        problems.append("%s: tallies missing/non-object" % tag)
        return problems
    for key in ZERO_TALLIES:
        if key not in tallies:
            problems.append("%s: tallies[%s] is missing" % (tag, key))
        elif type(tallies[key]) is not int or tallies[key] != 0:
            problems.append("%s: tallies[%s]=%r, expected integer 0" %
                            (tag, key, tallies[key]))
    if type(tallies.get("processed")) is not int or tallies.get("processed", 0) <= 0:
        problems.append("%s: processed must be a positive integer" % tag)
    if tallies.get("residual_geometry") != {}:
        problems.append("%s: residual_geometry must be empty" % tag)
    if record.get("trace_schema") != "cert3-trace-v1":
        problems.append("%s: wrong trace schema" % tag)
    trace_hash = record.get("trace_sha256")
    if type(trace_hash) is not str or len(trace_hash) != 64 \
            or any(char not in "0123456789abcdef" for char in trace_hash):
        problems.append("%s: trace SHA-256 is not lowercase hexadecimal" % tag)
    trace_size = record.get("trace_size")
    if type(trace_size) is not int or trace_size <= 0:
        problems.append("%s: trace size is not a positive integer" % tag)
    elif trace_size != tallies.get("processed"):
        problems.append("%s: trace size differs from processed tally" % tag)
    return problems


def audit(launch_path):
    with open(launch_path) as fh:
        launch = json.load(fh)
    problems = launch_problems(launch, launch_path)
    runs = launch.get("runs") if isinstance(launch.get("runs"), list) else []
    campaign_dir = os.path.dirname(os.path.abspath(launch_path))
    expected_names = {run.get("result_filename") for run in runs
                      if isinstance(run, dict)}
    actual_names = {name for name in os.listdir(campaign_dir)
                    if name.startswith("result_") and name.endswith(".json")}
    unexpected = sorted(actual_names - expected_names)
    missing = sorted(expected_names - actual_names)
    if unexpected:
        problems.append("unexpected/duplicate result artifacts: %s" % unexpected)
    if missing:
        problems.append("missing result artifacts: %s" % missing)
    if len(actual_names) != NSLICES:
        problems.append("found %d committed results, expected exactly 8" %
                        len(actual_names))
    expected_traces = {run.get("trace_filename") for run in runs
                       if isinstance(run, dict)}
    actual_traces = {name for name in os.listdir(campaign_dir)
                     if name.startswith("trace_") and name.endswith(".bin")}
    unexpected_traces = sorted(actual_traces - expected_traces)
    missing_traces = sorted(expected_traces - actual_traces)
    if unexpected_traces:
        problems.append("unexpected/duplicate trace artifacts: %s" %
                        unexpected_traces)
    if missing_traces:
        problems.append("missing trace artifacts: %s" % missing_traces)
    if len(actual_traces) != NSLICES:
        problems.append("found %d committed traces, expected exactly 8" %
                        len(actual_traces))
    records = {}
    if not problems:
        for run in sorted(runs, key=lambda item: item["slice"]):
            name = run["result_filename"]
            path = os.path.join(campaign_dir, name)
            try:
                with open(path) as fh:
                    record = json.load(fh)
            except (OSError, ValueError) as exc:
                problems.append("%s cannot be parsed: %s" % (name, exc))
                continue
            problems.extend(result_problems(record, launch, run, name))
            trace_name = run["trace_filename"]
            trace_path = os.path.join(campaign_dir, trace_name)
            try:
                actual_trace_size = os.path.getsize(trace_path)
            except OSError as exc:
                problems.append("%s cannot be sized: %s" % (trace_name, exc))
                continue
            if actual_trace_size != record.get("trace_size"):
                problems.append(
                    "%s: file size %d differs from recorded trace size %r" %
                    (trace_name, actual_trace_size, record.get("trace_size")))
            if file_sha256(trace_path) != record.get("trace_sha256"):
                problems.append(
                    "%s: file SHA-256 differs from recorded trace digest" %
                    trace_name)
            records[run["slice"]] = record
    if not problems:
        problems.extend(replay_all_traces(launch, campaign_dir, runs, records))
    return launch, problems


def replay_all_traces(launch, campaign_dir, runs, records):
    """Mathematically re-prove all eight traces with the frozen snapshot.

    This runs only after every hash in the snapshot has been validated, so the
    imported modules are exactly the pinned bytes.  Any import or replay
    failure is a rejection, not a crash.
    """
    problems = []
    snapshot = os.path.realpath(os.path.join(campaign_dir, "snapshot"))
    sys.path.insert(0, snapshot)
    try:
        import cert3_replay
    except Exception as exc:
        return ["cannot import frozen replay module: %r" % exc]
    loaded = os.path.dirname(os.path.realpath(cert3_replay.__file__))
    if loaded != snapshot:
        return ["frozen replay module resolved outside the snapshot: %s"
                % loaded]
    for run in sorted(runs, key=lambda item: item["slice"]):
        idx = run["slice"]
        record = records[idx]
        trace_path = os.path.join(campaign_dir, run["trace_filename"])
        root = ((0.0, 1.0),) * 4 + (
            (float(Fraction(run["w_lo"])), float(Fraction(run["w_hi"]))),)
        started = time.monotonic()
        try:
            replayed = cert3_replay.replay_trace(
                launch["t_decimal"], root, trace_path,
                parameters=launch["parameters"],
                expected_tallies=record["tallies"])
        except Exception as exc:
            problems.append("slice %d trace replay FAILED: %s" % (idx, exc))
            continue
        print("REPLAY slice=%d PASS processed=%d split=%d face=%d "
              "infeasible=%d elapsed=%.1fs" %
              (idx, replayed["processed"], replayed["split"],
               replayed["face"], replayed["infeasible"],
               time.monotonic() - started), flush=True)
    return problems
    return launch, problems


def main(argv):
    if len(argv) != 1:
        raise SystemExit(__doc__)
    launch_path = os.path.abspath(argv[0])
    own_hash = file_sha256(os.path.abspath(__file__))
    print("CHECKER_SHA256 %s" % own_hash)
    launch, problems = audit(launch_path)
    print("CAMPAIGN_ID %s" % launch.get("campaign_id"))
    print("CODE_SHA256 %s" % launch.get("code_sha256"))
    print("LAUNCH_SHA256 %s" % launch.get("launch_sha256"))
    if problems:
        print("NO COMPOSITE CERTIFICATE -- %d problem(s):" % len(problems))
        for problem in problems:
            print("  - %s" % problem)
        return 1
    print("ALL MACHINE CHECKS PASS:")
    print("  exactly one fresh, run-unique result for each slice 0..7")
    print("  exact non-overlapping dyadic partition of w in [1/2,1]")
    print("  frozen executable/review hashes and launch/result hashes match")
    print("  every worker exited 0 and returned normally")
    print("  every verdict COMPLETE with processed>0 and no residual/budget")
    print("  every committed trace re-hashes to its record with size equal")
    print("  to the processed tally")
    print("  every trace REPLAYED: boxes re-derived from the exact dyadic")
    print("  root and every claimed discharge re-proved in Arb arithmetic")
    print("COMPOSITE CERTIFICATE: Phi >= 0 on the feasible 5-parameter family")
    print("at certified rational t = %s >= exact psi + %s over w in [1/2,1]."
          % (launch["t_decimal"], launch["level_offset"]))
    print("The proved orbit-swap symmetry covers w in [0,1].")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
