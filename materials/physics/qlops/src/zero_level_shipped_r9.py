"""Revision-9 four-cell shipped experiment; execute only on the assigned mini.

The R8 adapter remains an unchanged library of pinned source, smoke, sampling,
and binomial-statistic primitives. This runner never calls its R8 production
recovery or analyzer and never overwrites an author's shipped circuit.
"""
from __future__ import annotations

import argparse
import io
import json
import math
import re
import socket
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pymatching
import stim
import zero_level_repo_repro as Z

GATE = "optional-follow-up-new-cycle-would-need-funding-"
OBJECTIVE = "zero-level-shipped-fz4-four-cell-r9"
PREREG = "pre_statement_optional-follow-up-new-cycle-would-need-funding-.md"
PREREG_SHA = "a2cec34d09471d2664cef152d0322106404808fe4db019b3293b80bca26a791d"
CELLS = ("0.0008", "0.0006", "0.0004", "0.0001")
REFERENCE_RUN = "20260908T195555Z_67d0f894_cfb28f9ee00f"
SOURCE_MEMBERS_SHA = "1a270099609a3376e607386b5d1b0cc7732e330e04bade81a9dd0b7cd06301b7"


def require(ok, message):
    if not ok:
        raise Z.Refusal(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def retain(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == data, f"refusing changed existing artifact: {path}")
    else:
        Z.atomic_write_bytes(path, data)


def relative_source_paths(value, source_root):
    if isinstance(value, dict):
        return {k: relative_source_paths(v, source_root) for k, v in value.items()}
    if isinstance(value, list):
        return [relative_source_paths(v, source_root) for v in value]
    if isinstance(value, str) and value.startswith(str(source_root) + "/"):
        return "source/" + value[len(str(source_root)) + 1:]
    return value


class Progress:
    def __init__(self, path):
        self.path = path

    def emit(self, event, **fields):
        row = {"utc": now(), "event": event, **fields}
        with self.path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
            stream.flush()
        print(json.dumps(row, sort_keys=True), flush=True)

    def sampler_factory(self, cell, phase):
        progress = self

        class ObservedSampler:
            def __init__(self, circuit, seed):
                self.inner = circuit.compile_detector_sampler(seed=seed)
                self.shots = 0
                self.chunks = 0

            def sample(self, n, **kwargs):
                arrays = self.inner.sample(n, **kwargs)
                self.shots += n
                self.chunks += 1
                if self.chunks == 1 or self.chunks % 10 == 0:
                    progress.emit("sampled", phase=phase, cell=cell,
                                  shots=self.shots, chunks=self.chunks)
                return arrays

        return ObservedSampler


def capture(run, tar_path):
    require(socket.gethostname() == "mini-pro.local", "assigned host must be mini-pro.local")
    require(run.parent == Path(__file__).resolve().parents[1] / "campaigns",
            "run must be in this target's campaigns directory")
    manifest_bytes = (run / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    prereg_bytes = (run.parent.parent / PREREG).read_bytes()
    require(Z.sha256_bytes(prereg_bytes) == PREREG_SHA, "preregistered bytes changed")
    require(set(manifest) == Z.MANIFEST_KEYS | {"prereg_file"},
            "noncanonical current campaign manifest schema")
    require(manifest["prereg_file"] == "prereg.md" and
            (run / "prereg.md").read_bytes() == prereg_bytes,
            "campaign-owned prereg snapshot mismatch")
    require(manifest["run_id"] == run.name and manifest["gate"] == GATE and
            manifest["target"] == Z.TARGET and manifest["status"] == "RUNNING" and
            manifest["prereg_sha256"] == PREREG_SHA, "campaign/prereg binding mismatch")
    match = Z.RUN_ID_RE.fullmatch(run.name)
    require(match is not None, "invalid minted run name")
    stamp = datetime.strptime(manifest["created_utc"], Z.TS_Z).strftime("%Y%m%dT%H%M%SZ")
    binder = "\x1f".join((GATE, manifest["agent"], PREREG_SHA, stamp, match["uuid8"]))
    require(match["stamp"] == stamp and
            match["hash12"] == Z.sha256_bytes(binder.encode())[:12], "mint binder mismatch")
    Z.assert_run_live(run, manifest_bytes)
    environment = Z.validate_env()
    tar_bytes = tar_path.read_bytes()
    require(len(tar_bytes) == Z.TAR_SIZE and Z.sha256_bytes(tar_bytes) == Z.TAR_SHA256,
            "canonical source tar changed")
    code = {name: (Path(__file__).parent / name).read_bytes() for name in
            ("zero_level_shipped_r9.py", "zero_level_repo_repro.py")}
    identity = {"run": run.name, "gate": GATE, "objective": OBJECTIVE,
                "prereg_sha256": PREREG_SHA, "adapter_version": "r9-shipped-direct",
                "commit": Z.COMMIT, "tree_sha": Z.TREE_SHA, "tar_sha256": Z.TAR_SHA256,
                "source_members_sha256": SOURCE_MEMBERS_SHA, "env": environment,
                "host": socket.gethostname(), "stream_note": Z.STREAM_NOTE,
                "code_sha256": {name: Z.sha256_bytes(data) for name, data in code.items()}}
    return {"manifest": manifest_bytes, "prereg": prereg_bytes,
            "tar": tar_bytes, "code": code, "identity": identity}


def live(run, ctx):
    Z.assert_run_live(run, ctx["manifest"])
    Z.assert_prereg_unchanged(run.parent.parent / PREREG, ctx["prereg"])
    require((run / "prereg.md").read_bytes() == ctx["prereg"],
            "campaign-owned prereg snapshot changed")
    for name, data in ctx["code"].items():
        require((Path(__file__).parent / name).read_bytes() == data, "running code changed")


def extraction(tar_bytes, destination):
    result = Z.extract_pinned_members(io.BytesIO(tar_bytes), destination, Z.MEMBER_PINS)
    require(Z.source_members_digest(result) == SOURCE_MEMBERS_SHA, "source member digest changed")
    return result


def recover_shipped(source, p):
    require(p in CELLS, "R9 recovery is restricted to four funded F-Z4 cells")
    vdir = source / Z.REPO_ROOT_NAME / "stim" / "grown"
    driver_path, driver_blob = Z.MEMBER_PINS["driver_grown"]
    driver_bytes = (source / Z.REPO_ROOT_NAME / driver_path).read_bytes()
    require(Z.git_blob_sha1(driver_bytes) == driver_blob, "driver blob mismatch")
    builder, built, prefix_sha = Z._exec_driver_prefix(
        Z.split_driver_prefix(driver_bytes.decode()), vdir, float(p))
    require(isinstance(built, list) and len(built) == 1, "unexpected builder result")
    rebuilt = built[0]
    flags = builder.postselct_numbers()
    mask = [i for i, flag in enumerate(flags) if flag]
    shipped_rel, shipped_blob = Z.MEMBER_PINS[f"shipped_grown_{p}"]
    shipped_bytes = (source / Z.REPO_ROOT_NAME / shipped_rel).read_bytes()
    require(Z.git_blob_sha1(shipped_bytes) == shipped_blob, "shipped blob mismatch")
    circuit = stim.Circuit(shipped_bytes.decode())
    require(circuit.num_detectors == rebuilt.num_detectors == len(flags) == 841,
            "detector index/mask length mismatch")
    require(circuit.num_observables == rebuilt.num_observables == 3, "observable count mismatch")
    defs = lambda c, name: [str(inst) for inst in c.flattened() if inst.name == name]
    definitions = {}
    for name in ("DETECTOR", "OBSERVABLE_INCLUDE"):
        shipped_defs, rebuilt_defs = defs(circuit, name), defs(rebuilt, name)
        require(shipped_defs == rebuilt_defs, f"ordered {name} definitions differ")
        definitions[name] = {"count": len(shipped_defs), "ordered_definitions": shipped_defs,
                             "sha256": Z.sha256_bytes(canonical(shipped_defs)), "equal": True}
    equal = circuit.flattened() == rebuilt.flattened()
    require(not equal, "expected F-Z4 source-generation split absent")
    dem = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    require(matcher.num_fault_ids == 3, "shipped DEM matcher must expose three fault IDs")
    provenance = {"recovery_mode": "shipped_direct_r9", "p_label": p,
                  "sampled_circuit_source": "pinned_shipped_bytes",
                  "dem_source": "same_pinned_shipped_circuit_object",
                  "shipped_oracle_used": True, "shipped_rebuilt_flattened_equal": equal,
                  "shipped_path": f"source/{Z.REPO_ROOT_NAME}/{shipped_rel}",
                  "shipped_blob_sha1": shipped_blob,
                  "shipped_bytes_sha256": Z.sha256_bytes(shipped_bytes),
                  "sampled_circuit_sha256": Z.sha256_bytes(str(circuit).encode()),
                  "sampled_flattened_sha256": Z.sha256_bytes(str(circuit.flattened()).encode()),
                  "rebuilt_flattened_sha256": Z.sha256_bytes(str(rebuilt.flattened()).encode()),
                  "driver_blob_sha1": driver_blob, "prefix_exec_sha256": prefix_sha,
                  "mask": mask, "mask_sha256": Z.sha256_bytes(canonical(mask)),
                  "mask_source": "pinned driver prefix postselct_numbers; upstream reader index semantics",
                  "num_detectors": circuit.num_detectors, "num_observables": 3,
                  "num_fault_ids": matcher.num_fault_ids,
                  "dem_sha256": Z.sha256_bytes(str(dem).encode()),
                  "definition_checks": definitions,
                  "no_generation_claim": "output methods not executed; historical generator not recovered"}
    return circuit, mask, matcher, provenance


def smoke(source, progress):
    progress.emit("inherited_smoke_start")
    gates, inherited = Z.run_smoke_gates(source)
    direct = []
    for p in CELLS:
        circuit, mask, matcher, recovery = recover_shipped(source, p)
        progress.emit("direct_smoke_start", cell=p)
        a = Z.sample_counts(circuit, mask, matcher, 20000, 5000, 3302,
                            crosscheck_accepted_rows=True)
        b = Z.sample_counts(circuit, mask, matcher, 20000, 5000, 3302)
        require(a["crosscheck"]["rows"] > 0 and a["crosscheck"]["mismatches"] == 0,
                f"{p} scalar/batch smoke failed")
        require(a["chunk_records"] == b["chunk_records"], f"{p} direct smoke replay failed")
        direct.append({"p_label": p, "recovery": recovery, "scalar_batch": a,
                       "repeat": b, "pass": True})
    return relative_source_paths({"inherited_gates": gates, "inherited_evidence": inherited,
                                  "direct_split_cells": direct, "pass": True}, source)


def spec(identity, p):
    shots = Z._author_row("grown", p)[1]
    return dict(identity, variant="grown", p_label=p, shots=shots,
                seed=Z.FULL_SEEDS["grown"][p], chunk=10000,
                schedule=Z.chunk_schedule(shots, 10000))


def validate_row(row, expected, recovery):
    require(row["identity"] == expected, "point identity mismatch")
    require(row["recovery"] == recovery and row["oracle"] == "shipped_direct_r9",
            "point recovery is not bound to freshly loaded shipped artifact")
    require(row["status"] == "ok" and row["reason"] is None and row["crosscheck"] is None,
            "non-production or failed point row")
    for key in ("shots", "seed", "schedule"):
        require(row[key] == expected[key], f"point {key} mismatch")
    records, schedule = row["chunk_records"], expected["schedule"]
    require(len(records) == schedule["num_chunks"], "incomplete chunk ledger")
    for index, record in enumerate(records):
        n = schedule["final_chunk"] if index == len(records) - 1 else schedule["chunk"]
        require(set(record) == {"n", "accepted", "errors", "syndrome_obs_sha256"},
                "unexpected chunk schema")
        require(all(type(record[key]) is int for key in ("n", "accepted", "errors")) and
                record["n"] == n and 0 <= record["errors"] <= record["accepted"] <= n,
                "chunk count invariant failure")
        require(re.fullmatch(r"[0-9a-f]{64}", record["syndrome_obs_sha256"]) is not None,
                "invalid chunk digest")
    for key in ("accepted", "errors"):
        require(type(row[key]) is int and row[key] == sum(r[key] for r in records),
                f"chunk/top-level {key} mismatch")
    require(row["discarded"] == row["shots"] - row["accepted"], "discard count mismatch")
    require(row["acceptance"] == row["accepted"] / row["shots"], "acceptance mismatch")
    require(row["ler"] == (row["errors"] / row["accepted"] if row["accepted"] else None),
            "LER mismatch")


def analyze(rows, run):
    require([row["identity"]["p_label"] for row in rows] == list(CELLS), "four-cell set mismatch")
    comparisons = []
    for row in rows:
        p = row["identity"]["p_label"]
        _, shots, author_ler, author_acc = Z._author_row("grown", p)
        for kind, author in (("ler", author_ler), ("acceptance", author_acc)):
            value = row[kind]
            author_sigma = (Z.ler_sigma(shots, author_ler, author_acc) if kind == "ler"
                            else Z.acc_sigma(shots, author_acc))
            sigma = (Z.ler_sigma(row["shots"], value, row["acceptance"])
                     if kind == "ler" and value is not None else
                     Z.acc_sigma(row["shots"], value) if kind == "acceptance" else None)
            combined = math.hypot(author_sigma, sigma) if sigma is not None else None
            z = Z.zscore(value, author, combined) if combined is not None else None
            comparisons.append({"p_label": p, "quantity": kind, "author": author,
                                "repro": value, "author_shots": shots, "repro_shots": row["shots"],
                                "sigma_author": author_sigma, "sigma_repro": sigma,
                                "sigma_combined": combined, "z": z,
                                "pass": z is not None and abs(z) <= 3.53})
    ref_path = run / "reference_r8" / "grown_0.0006.json"
    ref = json.loads(ref_path.read_bytes())
    reference_ledger = dict(line.split("  ", 1)[::-1] for line in
                            (run / "reference_r8" / "sha256s.txt").read_text().splitlines())
    require(Z.sha256_file(ref_path) ==
            reference_ledger["results/points/grown_0.0006.json"],
            "historical reference bytes differ from frozen R8 checksum ledger")
    require(ref["identity"]["run_id"] == REFERENCE_RUN and ref["oracle"] == "rebuilt" and
            ref["p_label"] == "0.0006", "wrong historical signature reference")
    require(ref["shots"] == Z._author_row("grown", "0.0006")[1] and
            ref["accepted"] == sum(chunk["accepted"] for chunk in ref["chunk_records"]) and
            ref["acceptance"] == ref["accepted"] / ref["shots"],
            "historical reference acceptance/count mismatch")
    new = next(c for c in comparisons if c["p_label"] == "0.0006" and c["quantity"] == "acceptance")
    old_sigma = math.hypot(new["sigma_author"], Z.acc_sigma(ref["shots"], ref["acceptance"]))
    defined = [c["z"] for c in comparisons if c["z"] is not None]
    return {"run": run.name, "objective": OBJECTIVE, "cells": list(CELLS),
            "primary_shots": sum(r["shots"] for r in rows), "comparisons": comparisons,
            "pass_count": sum(c["pass"] for c in comparisons),
            "undefined_count": len(comparisons) - len(defined),
            "beyond_threshold": sum(abs(z) > 3.53 for z in defined),
            "decisive_count": sum(abs(z) >= 5 for z in defined),
            "threshold": 3.53, "retained_family": 24,
            "r8_signature": {"path": "reference_r8/grown_0.0006.json",
                             "sha256": Z.sha256_file(ref_path), "reference_run": REFERENCE_RUN,
                             "old_acceptance": ref["acceptance"],
                             "old_z": Z.zscore(ref["acceptance"], new["author"], old_sigma),
                             "shipped_acceptance": new["repro"], "shipped_z": new["z"],
                             "signed_difference_shipped_minus_rebuilt": new["repro"] - ref["acceptance"],
                             "interpretation": "descriptive cross-generation comparison; no historic-generator or causal attribution"},
            "limitations": ["four-cell finite check only", "normal-approximate inherited binomial gate",
                            "no sequential combined familywise claim", "same-seed replay is not new independent evidence",
                            "Gate A/B unchanged", "paper d7 and physical c~300 NOT-REPRODUCED"],
            "verdict_deferred_until_independent_review": True}


def primary(run, ctx, progress):
    source = run / "source"
    ext = extraction(ctx["tar"], source)
    retain(source / "source_tar.gz", ctx["tar"])
    for name, data in ctx["code"].items():
        retain(run / "provenance" / "snapshots" / name, data)
    retain(run / "provenance" / "snapshots" / "prereg.md", ctx["prereg"])
    retain(run / "provenance" / "identity.json", canonical(ctx["identity"]) + b"\n")
    retain(run / "provenance" / "source_pins.json", canonical(ext["pins"]) + b"\n")
    readiness = smoke(source, progress)
    live(run, ctx)
    retain(run / "results" / "smoke.json", canonical(readiness) + b"\n")
    recovered = {p: recover_shipped(source, p) for p in CELLS}
    rows = {}
    for p in CELLS:
        path = run / "results" / "points" / f"grown_{p}.json"
        if path.exists():
            row = json.loads(path.read_bytes())
            validate_row(row, spec(ctx["identity"], p), recovered[p][3])
            rows[p] = row
    for p in CELLS:
        live(run, ctx)
        if p in rows:
            progress.emit("validated_existing_primary", cell=p)
            continue
        circuit, mask, matcher, recovery = recovered[p]
        expected = spec(ctx["identity"], p)
        progress.emit("primary_cell_start", cell=p, shots=expected["shots"])
        counts = Z.sample_counts(circuit, mask, matcher, expected["shots"], 10000,
                                 expected["seed"], _sampler_factory=progress.sampler_factory(p, "primary"))
        row = {"identity": expected, "recovery": recovery, "oracle": "shipped_direct_r9", **counts}
        validate_row(row, expected, recovery)
        live(run, ctx)
        retain(run / "results" / "points" / f"grown_{p}.json", canonical(row) + b"\n")
        rows[p] = row
        progress.emit("primary_cell_complete", cell=p, shots=row["shots"],
                      accepted=row["accepted"], errors=row["errors"])
    analysis = analyze([rows[p] for p in CELLS], run)
    live(run, ctx)
    retain(run / "results" / "analysis.json", canonical(analysis) + b"\n")
    progress.emit("primary_complete", primary_shots=analysis["primary_shots"],
                  comparisons=len(analysis["comparisons"]), verdict="deferred to independent review")


def verification_context(run, ctx):
    producer = json.loads((run / "provenance" / "identity.json").read_bytes())
    require(producer == dict(ctx["identity"], code_sha256=producer["code_sha256"]),
            "producer/current environment or input identity differs")
    for name, digest in producer["code_sha256"].items():
        require(Z.sha256_file(run / "provenance" / "snapshots" / name) == digest,
                "retained producer code does not match primary identity")
    require(producer["code_sha256"]["zero_level_repo_repro.py"] ==
            ctx["identity"]["code_sha256"]["zero_level_repo_repro.py"],
            "shared scientific primitives changed after primary sampling")
    for name, data in ctx["code"].items():
        retain(run / "provenance" / "verification_snapshots" / name, data)
    retain(run / "provenance" / "verification_identity.json",
           canonical(ctx["identity"]) + b"\n")
    return producer


def analyze_existing(run, ctx, progress):
    producer = verification_context(run, ctx)
    rows = []
    with tempfile.TemporaryDirectory(prefix="r9-analysis-", dir=run.parent.parent) as tmp:
        fresh = Path(tmp)
        extraction(ctx["tar"], fresh)
        for p in CELLS:
            recovery = recover_shipped(fresh, p)[3]
            row = json.loads((run / "results" / "points" / f"grown_{p}.json").read_bytes())
            validate_row(row, spec(producer, p), recovery)
            rows.append(row)
    analysis = analyze(rows, run)
    live(run, ctx)
    retain(run / "results" / "analysis.json", canonical(analysis) + b"\n")
    progress.emit("analysis_complete_without_primary_resampling",
                  primary_shots=analysis["primary_shots"], comparisons=analysis["comparisons"],
                  producer_code_sha256=producer["code_sha256"],
                  analysis_code_sha256=ctx["identity"]["code_sha256"])


def audit(run, ctx, progress):
    producer = verification_context(run, ctx)
    require((run / "source" / "source_tar.gz").read_bytes() == ctx["tar"], "retained tar differs")
    require((run / "provenance" / "snapshots" / "prereg.md").read_bytes() == ctx["prereg"],
            "prereg snapshot differs")
    rows, replay = [], []
    with tempfile.TemporaryDirectory(prefix="r9-fresh-audit-", dir=run.parent.parent) as tmp:
        fresh = Path(tmp)
        extraction(ctx["tar"], fresh)
        readiness = smoke(fresh, progress)
        stored_smoke = json.loads((run / "results" / "smoke.json").read_bytes())
        require(readiness == stored_smoke, "fresh complete smoke evidence differs")
        for p in CELLS:
            live(run, ctx)
            circuit, mask, matcher, recovery = recover_shipped(fresh, p)
            row = json.loads((run / "results" / "points" / f"grown_{p}.json").read_bytes())
            expected = spec(producer, p)
            validate_row(row, expected, recovery)
            progress.emit("full_replay_start", cell=p, shots=expected["shots"])
            fresh_counts = Z.sample_counts(circuit, mask, matcher, expected["shots"], 10000,
                                           expected["seed"], _sampler_factory=progress.sampler_factory(p, "audit"))
            mismatches = [key for key, value in fresh_counts.items() if row.get(key) != value]
            require(not mismatches, f"{p} full replay mismatch: {mismatches}")
            replay.append({"p_label": p, "shots": row["shots"], "accepted": row["accepted"],
                           "errors": row["errors"], "chunks": len(row["chunk_records"]),
                           "chunk_ledger_sha256": Z.sha256_bytes(canonical(row["chunk_records"])),
                           "every_count_and_chunk_digest_equal": True})
            rows.append(row)
            progress.emit("full_replay_complete", **replay[-1])
    analysis = analyze(rows, run)
    require(analysis == json.loads((run / "results" / "analysis.json").read_bytes()),
            "fresh eight-comparison analysis differs")
    live(run, ctx)
    report = {"run": run.name, "producer_identity": producer,
              "verification_identity": ctx["identity"], "finished_utc": now(),
              "fresh_extraction_and_pins": True, "fresh_smoke_evidence_equal": True,
              "point_identity_recovery_counts_validated": True,
              "replay": replay, "all_eight_comparisons_recomputed_equal": True,
              "analysis_sha256": Z.sha256_file(run / "results" / "analysis.json"),
              "pass": True, "statistical_verdict_not_selected_by_audit": True}
    retain(run / "results" / "audit.json", canonical(report) + b"\n")
    progress.emit("fresh_audit_complete", replayed_cells=len(replay), pass_all=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("primary", "analyze", "audit"))
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--source-tar", required=True, type=Path)
    ns = parser.parse_args()
    run = ns.run_dir.resolve()
    progress = Progress(run / f"{ns.mode}-progress.jsonl")
    try:
        ctx = capture(run, ns.source_tar.resolve())
        progress.emit("start", mode=ns.mode, identity=ctx["identity"])
        {"primary": primary, "analyze": analyze_existing, "audit": audit}[ns.mode](run, ctx, progress)
        return 0
    except Z.Refusal as exc:
        progress.emit("administrative_refusal", reason=str(exc))
        print(f"REJECTED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
