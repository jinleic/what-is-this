#!/usr/bin/env python
"""Released vectorized-RREF census runner — 2026-08-31T16-12-46Z.

STRICT RELEASE GATE: this runner refuses to execute unless all three hold:
  1. env TOY_FAST_CENSUS_RELEASED == "1" (Main's explicit release),
  2. env PARENT_CENSUS_TERMINATED == "1", and
  3. CLI flag --i-have-main-release.
Static-only checks (argument parse, refusal path) are
allowed without release; ANY arithmetic (controls, census, P1-P4) requires
the gate.

Order under release:
  1. compute the five frozen source hashes FIRST (exact dict for
     run_census/run_controls headers),
  2. run_controls() ONCE; require ALL_OK; write checksummed
     controls artifact atomically (write-once semantics),
  3. run_census() per cell with per-cell ledger headers carrying BOTH
     source and control hashes,
  4. machine summary over all 24 ledgers: separate
     P1_U3/P1_U4/P2/P3/P4 pass/fail/inconclusive counts; genericity_pass
     is True ONLY if every required predicate passed in every cell; a
     completed-but-failing census reports status FAILED, never generic OK.
  5. final manifest maps every ledger FILE NAME to ITS OWN sha256 hash
     (not merely a name list); manifest/checksums NEVER overwrite
     existing evidence (write-once or byte-identical).
Durability: atomic writes fsync the file AND the parent directory;
if either fsync fails the write is refused and the claim weakened —
never silently accepted.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)      # FROZEN module only — no live mceliece/src

RELEASE_KEY = "TOY_FAST_CENSUS_RELEASED"
PARENT_RELEASE_KEY = "PARENT_CENSUS_TERMINATED"
STATIC_LEDGER = os.path.join(HERE, "sha256s_static.txt")
EXPECTED_LEDGER_DIR = os.path.abspath(os.path.join(HERE, "..", "state"))
WALL_CAP_SECONDS = 2 * 60 * 60


SOURCES = ["toy_census.py", "toy_static_guards.py", "run_toy_census.py",
           "fastfield_frozen.py", "gfield_frozen.py"]


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def source_hashes() -> dict:
    """exact frozen source hash map, computed FIRST so the census and
    control headers both carry it."""
    return {src: sha256_of(os.path.join(HERE, src)) for src in SOURCES}

def verify_static_ledger() -> dict:
    expected = {}
    with open(STATIC_LEDGER, encoding="utf-8") as stream:
        for line in stream:
            digest, name = line.rstrip("\n").split("  ", 1)
            expected[name] = digest
    if set(expected) != set(SOURCES):
        raise RuntimeError(
            f"static source ledger names mismatch: {sorted(expected)}")
    actual = source_hashes()
    if actual != expected:
        raise RuntimeError(
            "static source hash mismatch before frozen-module import")
    return actual


def assert_release_environment() -> None:
    if not __debug__:
        raise RuntimeError("optimized Python is forbidden for census release")
    if os.getpriority(os.PRIO_PROCESS, 0) < 10:
        raise RuntimeError("process niceness must be at least 10")
    required = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    for name in required:
        if os.environ.get(name) != "1":
            raise RuntimeError(f"{name} must be exactly 1")

def wall_cap_handler(_signum, _frame) -> None:
    raise TimeoutError(
        f"McEliece census exceeded {WALL_CAP_SECONDS} seconds; "
        "preserving every completed or staged artifact")

def fsync_directory(path: str) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def prepare_ledger_dir(path: str) -> str:
    path = os.path.abspath(path)
    if path != EXPECTED_LEDGER_DIR:
        raise RuntimeError(
            f"ledger must be exactly {EXPECTED_LEDGER_DIR!r}")
    if os.path.lexists(path):
        if os.path.islink(path) or not os.path.isdir(path):
            raise RuntimeError(
                f"campaign state path is non-directory evidence: {path}")
    else:
        os.mkdir(path, 0o700)
        fsync_directory(os.path.dirname(path))
    allowed = {
        "control_plants", "controls_result.json", "manifest.json",
        "checksums.sha256",
        *(f"cell_{beta}.json" for beta in range(8, 32)),
    }
    unexpected = sorted(set(os.listdir(path)) - allowed)
    if unexpected:
        raise RuntimeError(
            f"unexpected preserved campaign-state entries: {unexpected}")
    for name in os.listdir(path):
        entry = os.path.join(path, name)
        if name == "control_plants":
            if os.path.islink(entry) or not os.path.isdir(entry):
                raise RuntimeError(
                    f"control_plants is non-directory evidence: {entry}")
        elif os.path.islink(entry) or not os.path.isfile(entry):
            raise RuntimeError(
                f"state entry is non-regular evidence: {entry}")
    return path




def write_once_atomic(path: str, payload_bytes: bytes) -> None:
    """Atomic write with byte-idempotent, non-overwrite evidence semantics."""
    if os.path.lexists(path):
        if os.path.islink(path) or not os.path.isfile(path):
            raise RuntimeError(
                f"NONREGULAR EVIDENCE REFUSED: {path}; preserving it")
        with open(path, "rb") as stream:
            existing = stream.read()
        if existing == payload_bytes:
            return
        raise RuntimeError(
            f"MANIFEST EVIDENCE OVERWRITE REFUSED: {path} already exists "
            f"with different bytes; preserving existing evidence")
    tmp = path + ".tmp"
    if os.path.lexists(tmp):
        raise RuntimeError(
            f"STAGED WRITE REFUSED: {tmp} already exists; refusing to "
            f"silently truncate; preserve and inspect it first")
    descriptor = os.open(
        tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(payload_bytes)
        stream.flush()
        os.fsync(stream.fileno())
    dirfd = os.open(os.path.dirname(path) or ".", os.O_RDONLY)
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)
    if os.path.lexists(path):
        raise RuntimeError(
            f"TARGET APPEARED DURING WRITE: {path}; preserving staged temp")
    os.replace(tmp, path)
    dirfd = os.open(os.path.dirname(path) or ".", os.O_RDONLY)
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)


def summarize(cell_recs: list) -> dict:
    """Machine summary over all validated cell ledgers.

    OWNER AUDIT SEMANTICS:
      - P1 is aggregated at CELL level: a cell's P1_U3 passes iff ALL 5
        of its branch verdicts are P1_PASS; likewise P1_U4; the cell P1
        passes iff BOTH shapes pass.  Precedence: fail > inconclusive >
        pass.
      - Branch-level totals are INFORMATIONAL (P1_U3_branches /
        P1_U4_branches), never used for genericity_pass.
      - A P2 set inequality is a definitive FAIL.
      - P3 pass requires the run_cell boolean (exactly five mapped
        branches, zero unmapped, all singletons)."""
    summary = {"cells": len(cell_recs),
               "P1_U3": {"pass": 0, "fail": 0, "inconclusive": 0},
               "P1_U4": {"pass": 0, "fail": 0, "inconclusive": 0},
               "P2": {"pass": 0, "fail": 0, "inconclusive": 0},
               "P3": {"pass": 0, "fail": 0, "inconclusive": 0},
               "P4": {"pass": 0, "fail": 0, "inconclusive": 0},
               "P1_U3_branches": {"pass": 0, "fail": 0,
                                  "inconclusive": 0},
               "P1_U4_branches": {"pass": 0, "fail": 0,
                                  "inconclusive": 0}}

    def bump(d, verdict):
        v = {"P1_PASS": "pass", "P1_FAIL": "fail",
             "P1_INCONCLUSIVE": "inconclusive"}.get(verdict)
        if v is None:
            v = "fail"
        d[v] += 1

    def shape_verdict(shape):
        verdicts = [b["verdict"]["verdict"] for b in shape["branches"]]
        if any(v not in ("P1_PASS", "P1_INCONCLUSIVE") for v in verdicts):
            return "P1_FAIL"
        if all(v == "P1_PASS" for v in verdicts):
            return "P1_PASS"
        return "P1_INCONCLUSIVE"

    for rec in cell_recs:
        for shape_key, branch_key in (("P1_U3", "P1_U3_branches"),
                                      ("P1_U4", "P1_U4_branches")):
            sv = shape_verdict(rec[shape_key])
            bump(summary[shape_key], sv)
            for b in rec[shape_key]["branches"]:
                bump(summary[branch_key], b["verdict"]["verdict"])
        bump(summary["P2"], "P1_PASS" if rec["P2"]["passes"] else "P1_FAIL")
        bump(summary["P3"], "P1_PASS" if rec["P3"]["passes"] else "P1_FAIL")
        bump(summary["P4"], "P1_PASS" if rec["P4"]["passes"] else "P1_FAIL")
    required = ("P1_U3", "P1_U4", "P2", "P3", "P4")
    summary["genericity_pass"] = (
        summary["cells"] == 24
        and all(summary[k]["pass"] == 24 and summary[k]["fail"] == 0
                and summary[k]["inconclusive"] == 0 for k in required))
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--i-have-main-release", action="store_true")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--ledger", default=EXPECTED_LEDGER_DIR)
    args = ap.parse_args()
    if (os.environ.get(RELEASE_KEY, "") != "1"
            or os.environ.get(PARENT_RELEASE_KEY, "") != "1"
            or not args.i_have_main_release):
        print(json.dumps({
            "status": "REFUSED",
            "reason": "need TOY_FAST_CENSUS_RELEASED=1, "
                      "PARENT_CENSUS_TERMINATED=1, and "
                      "--i-have-main-release",
            "static_only_next_step": "await parent termination and Main release",
        }, indent=1))
        return 2
    assert_release_environment()
    requested_ledger_dir = os.path.abspath(args.ledger)
    if requested_ledger_dir != EXPECTED_LEDGER_DIR:
        raise RuntimeError(
            f"ledger must be exactly the campaign state directory "
            f"{EXPECTED_LEDGER_DIR!r}")
    sh = verify_static_ledger()
    tc = importlib.import_module("toy_census")
    if (tc.RELEASE_KEY != RELEASE_KEY
            or tc.PARENT_RELEASE_KEY != PARENT_RELEASE_KEY):
        raise RuntimeError("frozen module release-key mismatch")
    if args.check_only:
        print(json.dumps({
            "status": "STATIC_RELEASE_CHECK_PASS",
            "sha256_sources": sh,
            "state_created": False,
            "controls_run": False,
            "census_run": False,
        }, indent=1))
        return 0
    signal.signal(signal.SIGALRM, wall_cap_handler)
    signal.alarm(WALL_CAP_SECONDS)
    ledger_dir = prepare_ledger_dir(requested_ledger_dir)
    control_plant_root = os.path.join(ledger_dir, "control_plants")
    t0 = time.time()
    # ---- 1. controls once; a complete exact artifact is reusable on resume ----
    controls_path = os.path.join(ledger_dir, "controls_result.json")
    if os.path.lexists(controls_path):
        if os.path.islink(controls_path) or not os.path.isfile(controls_path):
            raise RuntimeError(
                "controls_result.json is non-regular evidence; preserving it")
        with open(controls_path, encoding="utf-8") as stream:
            controls_res = json.load(stream)
        if not controls_res.get("ALL_OK"):
            raise RuntimeError(
                "existing controls artifact is not ALL_OK; preserving it")
        if controls_res.get("source_sha256") != sh:
            raise RuntimeError(
                "existing controls artifact source hashes mismatch")
        observed_plants = {
            os.path.relpath(path, control_plant_root): sha256_of(path)
            for path in tc.control_plant_paths(control_plant_root)
        }
        if controls_res.get("control_plant_sha256") != observed_plants:
            raise RuntimeError(
                "existing controls artifact plant hashes mismatch")
    else:
        controls_res = tc.run_controls(control_plant_root)
        if not controls_res.get("ALL_OK"):
            print(json.dumps({"status": "CONTROLS_FAILED",
                              "controls": controls_res}, indent=1))
            return 3
        controls_res["source_sha256"] = sh
        controls_payload = json.dumps(
            controls_res, sort_keys=True, default=str).encode("utf-8")
        write_once_atomic(controls_path, controls_payload)
    controls_artifact_sha = sha256_of(controls_path)

    # ---- 2. census with control+source hashes in each ledger header ----
    res = tc.run_census(ledger_dir, sh,
                        extra_header={"controls_artifact_sha256":
                                      controls_artifact_sha})
    # Final collection uses the records VALIDATED by run_census through
    # read_ledger_state (owner audit: never a raw json.load that bypasses
    # the semantic checker).
    validated = res.get("validated", {})
    missing = [b for b in range(8, 32) if b not in validated]
    if missing:
        raise RuntimeError(
            f"VALIDATED-RECORD GAP: cells {missing} have no validated "
            f"record; refusing to summarize from unvalidated bytes")
    cell_recs = [validated[b] for b in range(8, 32)]
    summary = summarize(cell_recs)

    # ---- 3. final manifest (owner-audit rules) ----
    # The ledger map covers the controls artifact, five persistent control
    # plants, and 24 cell files — never manifest/checksums themselves.
    manifest_path = os.path.join(ledger_dir, "manifest.json")
    checksums_path = os.path.join(ledger_dir, "checksums.sha256")
    plant_hashes = {}
    for path in tc.control_plant_paths(control_plant_root):
        relative = os.path.relpath(path, control_plant_root)
        plant_hashes[relative] = sha256_of(path)
    if controls_res.get("control_plant_sha256") != plant_hashes:
        raise RuntimeError(
            "control result plant hashes do not match persistent artifacts")
    ledger_hashes = {"controls_result.json": controls_artifact_sha}
    ledger_hashes.update({
        os.path.join("control_plants", name): digest
        for name, digest in plant_hashes.items()
    })
    for beta in range(8, 32):
        name = os.path.basename(tc.ledger_path(ledger_dir, beta))
        ledger_hashes[name] = sha256_of(os.path.join(ledger_dir, name))
    if os.path.lexists(manifest_path):
        if os.path.islink(manifest_path) or not os.path.isfile(manifest_path):
            raise RuntimeError(
                "manifest path is non-regular evidence; preserving it")
        # RESUME-AFTER-MANIFEST: the existing manifest is immutable
        # evidence.  Validate it (cell/controls hashes must match what we
        # just observed) and finish only the missing checksum file; never
        # regenerate a timestamp-different manifest over existing bytes.
        with open(manifest_path) as fh:
            prior = json.load(fh)
        if prior.get("ledgers_sha256") != ledger_hashes:
            raise RuntimeError(
                "MANIFEST MISMATCH: existing manifest.json records "
                "different ledger hashes than the current state; "
                "preserving existing evidence and refusing to rewrite")
        if prior.get("sha256_sources") != sh:
            raise RuntimeError(
                "MANIFEST MISMATCH: existing manifest.json records "
                "different source hashes; preserving evidence")
        if prior.get("controls_artifact_sha256") != controls_artifact_sha:
            raise RuntimeError(
                "MANIFEST MISMATCH: existing manifest.json records a "
                "different controls_artifact_sha256; preserving evidence")
        if prior.get("summary") != summary:
            raise RuntimeError(
                "MANIFEST MISMATCH: existing manifest.json records a "
                "different summary than the validated ledgers imply; "
                "preserving evidence")
        manifest = prior
        manifest_sha = sha256_of(manifest_path)
    else:
        manifest = {
            "campaign": "2026-08-31T16-12-46Z_vectorized-rref-replay",
            "released_at": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                         time.gmtime()),
            "wall_seconds": time.time() - t0,
            "cells_completed": res["cells_completed"],
            "n_cells": res["n"],
            "of": res["of"],
            "sha256_sources": sh,
            "controls_artifact_sha256": controls_artifact_sha,
            "summary": summary,
            "ledgers_sha256": ledger_hashes,
            "engine": "exact NumPy table-indexed RREF with per-cell E/kernel cache",
        }
        manifest_bytes = json.dumps(manifest, indent=1, sort_keys=True,
                                    default=str).encode("utf-8")
        write_once_atomic(manifest_path, manifest_bytes)
        manifest_sha = sha256_of(manifest_path)
    # checksums cover sources + controls + the 24 cells + the manifest
    cs_lines = []
    for src in SOURCES:
        cs_lines.append(f"{sh[src]}  code/{src}\n")
    for name in sorted(ledger_hashes):
        cs_lines.append(f"{ledger_hashes[name]}  state/{name}\n")
    cs_lines.append(f"{manifest_sha}  state/manifest.json\n")
    write_once_atomic(checksums_path, "".join(cs_lines).encode("utf-8"))

    # ---- 4. status: FAILED census never reports generic OK ----
    status = "OK" if summary["genericity_pass"] else "CENSUS_FAILED"
    print(json.dumps({"status": status,
                      "summary": summary,
                      "manifest_path": os.path.join(ledger_dir,
                                                    "manifest.json")},
                     indent=1))
    return 0 if summary["genericity_pass"] else 4


if __name__ == "__main__":
    sys.exit(main())
