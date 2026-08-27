#!/usr/bin/env python3
"""Fail-closed aggregation for eight independent UC arithmetic reports."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
from fractions import Fraction
import json
from pathlib import Path
import sys

SCHEMA = "uc-independent-arithmetic-composite-v1"
REPORT_SCHEMA = "uc-independent-arithmetic-report-v1"
EXPECTED_LOCK_SHA256 = "d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48"
EXPECTED_CAMPAIGN_ID = "20260818T212601Z_425f109c15b64a6198785c6cebbbdaab"
EXPECTED_VERIFIER_SHA256 = "e67058df1347545c39b6c2d27c53e758d68f4ff338e42225acff51836aa91a9a"
EXPECTED_RUNTIME_SEAL_RAW_SHA256 = "58c49334a92b2d9e8b02b672d9e67da01b772a626b0ec29e88347d1cc7f92a0c"
EXPECTED_RUNTIME_SEAL_REPORT_SHA256 = "2bf686b6d0ba496383c2dd2afe51e4bffb34ee7b9e18f3bd49e7f68f4bb5184c"
EXPECTED_TARGET = "0.3820660112501052"
EXPECTED_TARGET_UPPER = "115472366449356389981023/302231454903657293676544"
TALLY_NAMES = (
    "processed",
    "infeasible",
    "corner",
    "ratio",
    "center",
    "center_mixed",
    "center_mixed_swap",
    "center_w",
    "face",
    "residual",
    "split",
)
EXPECTED_TOTAL_PROCESSED = 488_465_854
EXPECTED_TOTAL_SPLITS = 244_232_923
EXPECTED_TOTAL_TERMINALS = 244_232_931


class VerificationError(RuntimeError):
    """A fail-closed aggregate rejection."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError("duplicate JSON key %r" % key)
        result[key] = value
    return result


def reject_constant(token):
    raise VerificationError("non-finite JSON constant %s" % token)


def decode(raw, label):
    try:
        return json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError("%s is malformed JSON" % label) from exc


def canonical_bytes(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def object_hash(value):
    reduced = {key: item for key, item in value.items() if key != "report_sha256"}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def read_regular(path, label):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), "%s must be a regular file" % label)
    raw = path.read_bytes()
    return raw, decode(raw, label)

def load_report_lock(path, report_directory):
    raw, record = read_regular(path, "independent report lock")
    require(
        isinstance(record, dict)
        and record.get("schema") == "uc-independent-report-lock-v1",
        "independent report lock schema differs",
    )
    require(record.get("report_sha256") == object_hash(record), "report lock self-hash differs")
    require(record.get("campaign_id") == EXPECTED_CAMPAIGN_ID, "report lock campaign differs")
    require(record.get("campaign_lock_sha256") == EXPECTED_LOCK_SHA256, "report lock campaign digest differs")
    require(record.get("verifier_sha256") == EXPECTED_VERIFIER_SHA256, "report lock verifier differs")
    require(
        record.get("runtime_seal_raw_sha256") == EXPECTED_RUNTIME_SEAL_RAW_SHA256,
        "report lock runtime seal differs",
    )
    reports = record.get("reports")
    require(isinstance(reports, list) and len(reports) == 8, "report lock needs eight reports")
    by_slice = {}
    for item in reports:
        require(isinstance(item, dict), "report lock item is malformed")
        slice_index = item.get("slice")
        require(
            isinstance(slice_index, int) and not isinstance(slice_index, bool)
            and 0 <= slice_index < 8 and slice_index not in by_slice,
            "report lock slice index is invalid",
        )
        expected_path = report_directory / ("slice-%d.json" % slice_index)
        require(item.get("path") == expected_path.name, "report lock path differs")
        require(
            isinstance(item.get("size"), int) and item["size"] > 0,
            "report lock size is invalid",
        )
        require(
            isinstance(item.get("raw_sha256"), str)
            and len(item["raw_sha256"]) == 64,
            "report lock raw hash is invalid",
        )
        by_slice[slice_index] = item
    launch = record.get("launch_manifest")
    require(isinstance(launch, dict), "report lock launch manifest is missing")
    launch_path = report_directory / launch.get("path", "")
    launch_raw, launch_report = read_regular(launch_path, "secure launch manifest")
    require(len(launch_raw) == launch.get("size"), "launch manifest size differs")
    require(
        hashlib.sha256(launch_raw).hexdigest() == launch.get("raw_sha256"),
        "launch manifest raw hash differs",
    )
    expected_arguments = [
        "-I",
        "-B",
        "verification/independent_arithmetic_replay_secure.py",
        "--runtime-seal",
        "verification/results/runtime-environment-seal-secure.json",
        "--lock",
        "verification/campaign-lock.json",
        "--campaign",
        "campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8",
        "--slice",
        "<0..7>",
        "--output",
        "verification/results/independent-arithmetic/secure-full/slice-<N>.json",
        "--checkpoint",
        "verification/results/independent-arithmetic/secure-full/checkpoint-<N>.json",
    ]
    require(
        isinstance(launch_report, dict)
        and launch_report.get("schema")
        == "uc-secure-independent-arithmetic-launch-v1"
        and launch_report.get("campaign_id") == EXPECTED_CAMPAIGN_ID
        and launch_report.get("campaign_lock_sha256") == EXPECTED_LOCK_SHA256
        and launch_report.get("verifier_sha256") == EXPECTED_VERIFIER_SHA256
        and launch_report.get("runtime_seal_raw_sha256")
        == EXPECTED_RUNTIME_SEAL_RAW_SHA256
        and launch_report.get("resume_used") is False
        and launch_report.get("initial_trace_offset") == 0
        and launch_report.get("arguments") == expected_arguments,
        "secure launch manifest semantics differ",
    )
    live = record.get("live_runtime_attestation")
    require(isinstance(live, dict), "report lock live attestation is missing")
    live_path = report_directory / live.get("path", "")
    live_raw, live_report = read_regular(live_path, "live runtime attestation")
    require(len(live_raw) == live.get("size"), "live attestation size differs")
    require(
        hashlib.sha256(live_raw).hexdigest() == live.get("raw_sha256"),
        "live attestation raw hash differs",
    )
    require(
        isinstance(live_report, dict)
        and live_report.get("outcome") == "PASS"
        and live_report.get("report_sha256") == live.get("report_sha256")
        and live_report.get("report_sha256") == object_hash(live_report),
        "live runtime attestation envelope differs",
    )
    require(
        live_report.get("schema") == "uc-live-replay-runtime-attestation-v1"
        and live_report.get("claim_status") == "MACHINE-VERIFIED"
        and live_report.get("process_count") == 8
        and live_report.get("runtime_seal_raw_sha256")
        == EXPECTED_RUNTIME_SEAL_RAW_SHA256
        and live_report.get("runtime_seal_report_sha256")
        == EXPECTED_RUNTIME_SEAL_REPORT_SHA256
        and live_report.get("launch_manifest_raw_sha256")
        == hashlib.sha256(launch_raw).hexdigest(),
        "live runtime attestation pins differ",
    )
    seal_path = (
        report_directory.parent.parent / "runtime-environment-seal-secure.json"
    )
    seal_raw, seal_report = read_regular(seal_path, "secure runtime seal")
    require(
        hashlib.sha256(seal_raw).hexdigest() == EXPECTED_RUNTIME_SEAL_RAW_SHA256
        and seal_report.get("report_sha256") == EXPECTED_RUNTIME_SEAL_REPORT_SHA256
        and seal_report.get("report_sha256") == object_hash(seal_report),
        "secure runtime seal envelope differs",
    )
    sealed_images = {
        item["path"]: item for item in seal_report.get("files", [])
        if isinstance(item, dict) and isinstance(item.get("path"), str)
    }
    allowed_additional_images = {
        (
            "/opt/homebrew/Cellar/python@3.14/3.14.3_1/Frameworks/Python.framework/Versions/3.14/Resources/Python.app/Contents/MacOS/Python",
            51392,
            "7ecc1ecbf9daa9303c4bf502ff62ffdd9010ed5c08729d470ae9380c10ce1211",
            False,
        ),
        (
            "/opt/homebrew/Cellar/python@3.14/3.14.3_1/Frameworks/Python.framework/Versions/3.14/lib/python3.14/lib-dynload/pyexpat.cpython-314-darwin.so",
            99456,
            "28284cd7ac83d5dafda621ee75401c3c144b75f15e410480f2f5037ce441287d",
            False,
        ),
    }
    processes = live_report.get("processes")
    require(isinstance(processes, list) and len(processes) == 8, "live process inventory differs")
    seen_slices = set()
    common_images = None
    for process in processes:
        require(isinstance(process, dict), "live process record is malformed")
        slice_index = process.get("slice")
        require(
            isinstance(slice_index, int) and not isinstance(slice_index, bool)
            and 0 <= slice_index < 8 and slice_index not in seen_slices,
            "live process slice is invalid",
        )
        seen_slices.add(slice_index)
        command = process.get("command")
        expected_command = (
            "/opt/homebrew/Cellar/python@3.14/3.14.3_1/Frameworks/"
            "Python.framework/Versions/3.14/Resources/Python.app/Contents/"
            "MacOS/Python -I -B verification/independent_arithmetic_replay_secure.py "
            "--runtime-seal verification/results/runtime-environment-seal-secure.json "
            "--lock verification/campaign-lock.json --campaign "
            "campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8 "
            "--slice %d --output verification/results/independent-arithmetic/"
            "secure-full/slice-%d.json --checkpoint verification/results/"
            "independent-arithmetic/secure-full/checkpoint-%d.json"
            % (slice_index, slice_index, slice_index)
        )
        require(command == expected_command, "live process command differs")
        images = process.get("loaded_non_system")
        require(isinstance(images, list) and len(images) == 64, "live image inventory differs")
        image_set = set()
        for image in images:
            require(
                isinstance(image, dict)
                and isinstance(image.get("path"), str)
                and isinstance(image.get("size"), int)
                and image["size"] > 0
                and isinstance(image.get("sha256"), str)
                and len(image["sha256"]) == 64
                and isinstance(image.get("present_in_prelaunch_seal"), bool),
                "live image record is malformed",
            )
            image_set.add((
                image["path"],
                image["size"],
                image["sha256"],
                image["present_in_prelaunch_seal"],
            ))
        require(len(image_set) == 64, "live image inventory contains duplicates")
        for path, size, digest, present in image_set:
            if present:
                sealed = sealed_images.get(path)
                require(
                    sealed is not None
                    and sealed.get("size") == size
                    and sealed.get("sha256") == digest,
                    "live sealed image differs from runtime seal",
                )
            else:
                require(
                    (path, size, digest, present) in allowed_additional_images,
                    "live unsealed image is not explicitly pinned",
                )
        require(
            {item for item in image_set if item[3] is False}
            == allowed_additional_images,
            "live additional image set differs",
        )
        if common_images is None:
            common_images = image_set
        require(image_set == common_images, "live worker image sets differ")
        require(
            process.get("system_unsealed") == ["/usr/lib/dyld"],
            "live system image boundary differs",
        )
    return raw, record, by_slice




def validate_report(
    path,
    expected_report_record,
    expected_slice,
    expected_run_id,
    expected,
    expected_trace_size,
    expected_campaign_id,
    expected_trace_sha,
    expected_launch_sha,
    expected_result_sha,
    expected_verifier_sha,
    lock_sha,
    common,
):
    raw, report = read_regular(path, "slice %d report" % expected_slice)
    require(len(raw) == expected_report_record["size"], "slice raw report size differs")
    require(
        hashlib.sha256(raw).hexdigest() == expected_report_record["raw_sha256"],
        "slice raw report digest differs from external report lock",
    )
    require(isinstance(report, dict), "slice report is not an object")
    require(report.get("schema") == REPORT_SCHEMA, "slice report schema differs")
    require(report.get("outcome") == "PASS", "slice report is not PASS")
    require(
        report.get("claim_status") == "MACHINE-VERIFIED",
        "slice report claim status differs",
    )
    require(report.get("slice") == expected_slice, "slice number differs")
    require(report.get("run_id") == expected_run_id, "slice run id differs")
    require(report.get("campaign_id") == expected_campaign_id, "slice campaign id differs")
    require(report.get("finite_event_limit") is None, "slice used a finite event limit")
    require(report.get("arithmetic_complete") is True, "slice arithmetic is incomplete")
    require(report.get("trace_eof") is True, "slice did not reach trace EOF")
    require(report.get("pending_stack") == 0, "slice retains DFS work")
    require(report.get("initial_trace_offset") == 0, "slice did not start at trace byte zero")
    require(report.get("resume_used") is False, "slice used checkpoint resume")
    require(
        report.get("expected_tallies_authenticated") is True,
        "slice tallies were not authenticated",
    )
    require(report.get("report_sha256") == object_hash(report), "slice self-hash differs")

    tallies = report.get("tallies")
    require(isinstance(tallies, dict), "slice tallies are missing")
    require(set(tallies) == set(TALLY_NAMES), "slice tally fields differ")
    require(tallies == expected, "slice tallies differ from campaign lock")
    require(report.get("processed") == expected["processed"], "processed total differs")
    require(report.get("trace_offset") == expected_trace_size, "slice trace offset differs")
    require(tallies["processed"] == expected_trace_size, "slice trace size differs")
    terminal_total = sum(
        tallies[name] for name in (
            "infeasible",
            "corner",
            "ratio",
            "center",
            "center_mixed",
            "center_mixed_swap",
            "center_w",
            "face",
            "residual",
        )
    )
    require(
        terminal_total == tallies["processed"] - tallies["split"],
        "slice terminal accounting differs",
    )
    require(
        report.get("target") == {
            "decimal": EXPECTED_TARGET,
            "upper_fraction": EXPECTED_TARGET_UPPER,
        },
        "slice target differs",
    )
    require(
        report.get("precision") == {
            "global_cover_bits": 160,
            "replay_bits": 80,
            "cover_cells": 20000,
        },
        "slice precision differs",
    )
    require(expected["residual"] == 0, "locked residual is nonzero")
    require(
        expected["processed"] == 2 * expected["split"] + 1,
        "slice tree identity differs",
    )

    inputs = report.get("input_sha256")
    require(isinstance(inputs, dict), "slice input hashes are missing")
    require(inputs.get("lock") == lock_sha, "slice lock hash differs")
    require(inputs.get("trace") == expected_trace_sha, "slice trace hash differs")
    require(inputs.get("launch") == expected_launch_sha, "slice launch hash differs")
    require(inputs.get("result") == expected_result_sha, "slice result hash differs")
    require(inputs.get("verifier") == expected_verifier_sha, "slice verifier hash differs")
    implementation = report.get("implementation")
    require(isinstance(implementation, dict), "slice implementation record is missing")
    require(
        implementation.get("imports_frozen_arithmetic") is False,
        "slice imported frozen arithmetic",
    )
    require(
        implementation.get("gradient_engine") == "interval automatic differentiation",
        "slice gradient implementation differs",
    )
    require(
        implementation.get("rediscover_face_proofs") is True,
        "slice did not rediscover face proofs",
    )
    require(
        implementation.get("shares_python_flint") is True,
        "slice primitive trust declaration differs",
    )
    require(
        implementation.get("byte_zero_required_for_pass") is True,
        "slice byte-zero PASS policy differs",
    )
    runtime_seal = report.get("runtime_seal")
    require(isinstance(runtime_seal, dict), "slice runtime seal is missing")
    require(
        runtime_seal.get("raw_sha256") == EXPECTED_RUNTIME_SEAL_RAW_SHA256,
        "slice runtime seal raw hash differs",
    )
    require(
        runtime_seal.get("report_sha256") == EXPECTED_RUNTIME_SEAL_REPORT_SHA256,
        "slice runtime seal report hash differs",
    )

    for field in ("campaign_id", "environment", "runtime_seal"):
        if field not in common:
            common[field] = report.get(field)
        require(report.get(field) == common[field], "slice %s differs" % field)
    verifier_hash = inputs.get("verifier")
    if "verifier_hash" not in common:
        common["verifier_hash"] = verifier_hash
    require(verifier_hash == common["verifier_hash"], "slice verifier hash differs")

    return {
        "slice": expected_slice,
        "report_path": str(Path(path)),
        "report_raw_sha256": hashlib.sha256(raw).hexdigest(),
        "report_sha256": report["report_sha256"],
        "processed": expected["processed"],
        "split": expected["split"],
        "terminal": expected["processed"] - expected["split"],
        "elapsed_seconds": report.get("elapsed_seconds"),
        "trace_sha256": inputs["trace"],
    }


def aggregate_reports(
    report_directory, lock_path, report_lock_path, expected_report_lock_sha256
):
    report_directory = Path(report_directory)
    report_lock_raw, report_lock, report_records = load_report_lock(
        report_lock_path, report_directory
    )
    report_lock_raw_sha256 = hashlib.sha256(report_lock_raw).hexdigest()
    require(
        isinstance(expected_report_lock_sha256, str)
        and len(expected_report_lock_sha256) == 64
        and report_lock_raw_sha256 == expected_report_lock_sha256,
        "report lock raw digest differs from the external trusted input",
    )
    lock_raw, lock = read_regular(lock_path, "campaign lock")
    lock_sha = hashlib.sha256(lock_raw).hexdigest()
    require(
        isinstance(lock, dict) and lock.get("schema") == "uc-campaign-lock-v1",
        "campaign lock schema differs",
    )
    require(lock_sha == EXPECTED_LOCK_SHA256, "campaign lock digest is not pinned")
    require(lock.get("campaign_id") == EXPECTED_CAMPAIGN_ID, "campaign id is not pinned")
    slices = lock.get("slices")
    require(isinstance(slices, list) and len(slices) == 8, "lock does not contain eight slices")
    launch_record = lock.get("launch")
    require(isinstance(launch_record, dict), "lock launch record is missing")
    expected_launch_sha = launch_record.get("raw_sha256")
    verifier_path = Path(__file__).with_name("independent_arithmetic_replay_secure.py")
    require(verifier_path.is_file(), "secure independent verifier source is missing")
    expected_verifier_sha = file_hash(verifier_path)
    require(
        expected_verifier_sha == EXPECTED_VERIFIER_SHA256,
        "secure independent verifier digest is not pinned",
    )

    common = {}
    summaries = []
    expected_w_lower = Fraction(1, 2)
    for index, lock_slice in enumerate(slices):
        require(
            isinstance(lock_slice, dict) and lock_slice.get("slice") == index,
            "lock slice order differs",
        )
        try:
            w_lower = Fraction(lock_slice["w_lo"])
            w_upper = Fraction(lock_slice["w_hi"])
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            raise VerificationError("lock weight slice is malformed") from exc
        require(w_lower == expected_w_lower, "weight slices have a gap or overlap")
        require(w_lower < w_upper <= 1, "weight slice is empty or out of range")
        expected_w_lower = w_upper
        expected = dict(lock_slice.get("expected_replay_tallies", {}))
        require(set(expected) == set(TALLY_NAMES), "locked tally fields differ")
        expected_run_id = lock_slice.get("run_id")
        expected_result_sha = lock_slice.get("result", {}).get("raw_sha256")
        expected_trace_sha = lock_slice.get("trace", {}).get("sha256")
        expected_trace_size = lock_slice.get("trace", {}).get("size")
        summaries.append(
            validate_report(
                report_directory / ("slice-%d.json" % index),
                report_records[index],
                index,
                expected_run_id,
                expected,
                expected_trace_size,
                EXPECTED_CAMPAIGN_ID,
                expected_trace_sha,
                expected_launch_sha,
                expected_result_sha,
                expected_verifier_sha,
                lock_sha,
                common,
            )
        )
    require(expected_w_lower == 1, "weight slices do not end at one")

    total_processed = sum(item["processed"] for item in summaries)
    total_splits = sum(item["split"] for item in summaries)
    total_terminals = sum(item["terminal"] for item in summaries)
    require(total_processed == EXPECTED_TOTAL_PROCESSED, "composite processed total differs")
    require(total_splits == EXPECTED_TOTAL_SPLITS, "composite split total differs")
    require(total_terminals == EXPECTED_TOTAL_TERMINALS, "composite terminal total differs")

    source = Path(__file__).resolve()
    result = {
        "schema": SCHEMA,
        "report_type": "independent_arithmetic_composite",
        "claim_status": "MACHINE-VERIFIED",
        "outcome": "PASS",
        "campaign_id": common["campaign_id"],
        "finished_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "lock": {
            "sha256": lock_sha,
            "size": len(lock_raw),
        },
        "report_lock": {
            "raw_sha256": hashlib.sha256(report_lock_raw).hexdigest(),
            "report_sha256": report_lock["report_sha256"],
            "size": len(report_lock_raw),
        },
        "independent_verifier_sha256": common["verifier_hash"],
        "environment": common["environment"],
        "slices": summaries,
        "total_processed": total_processed,
        "total_splits": total_splits,
        "total_terminals": total_terminals,
        "total_residual": 0,
        "aggregator": {
            "sha256": file_hash(source),
            "size": source.stat().st_size,
        },
        "limitations": [
            "The independent verifier shares the authenticated trace partition.",
            "The producer and independent verifier both trust python-flint/Arb primitives.",
            "The aggregate is a hash-bound public report, not authenticated execution attestation.",
        ],
    }
    result["report_sha256"] = object_hash(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Aggregate eight independent UC replay reports")
    parser.add_argument("--reports", required=True)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--report-lock", required=True)
    parser.add_argument("--report-lock-sha256", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        result = aggregate_reports(
            args.reports,
            args.lock,
            args.report_lock,
            args.report_lock_sha256,
        )
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("xb") as stream:
            stream.write(canonical_bytes(result) + b"\n")
    except (VerificationError, OSError, ValueError, TypeError) as exc:
        print("independent aggregate FAIL: %s" % exc, file=sys.stderr)
        return 1
    print(
        canonical_bytes(
            {
                "claim_status": result["claim_status"],
                "outcome": result["outcome"],
                "report": str(output),
                "report_sha256": result["report_sha256"],
                "total_processed": result["total_processed"],
            }
        ).decode("ascii")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
