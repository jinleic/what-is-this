#!/usr/bin/env python3
"""Check exact coverage and immutable hashes of the full fixed-five proof bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
W_CENSUS = ROOT / "data" / "involution_f5_w_square_census.json"
SIGNED_CENSUS = ROOT / "data" / "involution_f5_signed_square_census.json"
RAMSEY_CENSUS = ROOT / "data" / "involution_f5_ramsey_square_census.json"
DEFAULT_W_CERTIFICATES = ROOT / "data" / "involution_f5_w_dfs_certificates.json"
DEFAULT_HIERARCHY = ROOT / "data" / "involution_f5_hierarchical_certificates.json"
EXPECTED = {
    W_CENSUS: "e99963ed0494b90a7cd3cb117113b317132b711e0592867d27353e76e8a11a6a",
    SIGNED_CENSUS: "01a805f5a783d572966785cedb5c31c89b69ff986ed6bdc77bf3a5dc4c3b8150",
    RAMSEY_CENSUS: "6069740d3c5339dc1d734e10f70d07be53f7046ddcbf07c37974fab322f359c1",
}


class CoverageViolation(RuntimeError):
    """A manifest, target partition, proof archive, or pin is invalid."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise CoverageViolation(f"duplicate JSON member {key!r}")
        result[key] = value
    return result


def _load(path: Path) -> dict:
    try:
        document = json.loads(
            path.read_text(), object_pairs_hook=_unique_object,
            parse_constant=lambda value: (_ for _ in ()).throw(
                CoverageViolation(f"nonstandard JSON constant {value}")))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CoverageViolation(f"cannot load {path}: {error}") from error
    if type(document) is not dict:
        raise CoverageViolation(f"{path} root is not an object")
    return document


def _archive(relative: str, record: dict, nested: bool = False) -> None:
    if type(relative) is not str or not relative.startswith("r55/data/"):
        raise CoverageViolation("invalid proof archive path")
    path = ROOT / relative.removeprefix("r55/")
    expected = record.get("gzip") if nested else record
    if type(expected) is not dict or not path.is_file():
        raise CoverageViolation(f"missing proof archive {relative}")
    if (path.stat().st_size != expected.get("bytes")
            or _sha256(path) != expected.get("sha256")):
        raise CoverageViolation(f"proof archive hash mismatch: {relative}")


def _check_w_record(record: dict) -> int:
    if (record.get("status") != "VERIPB_CAKEPB_VERIFIED"
            or record.get("source_index") is None):
        raise CoverageViolation("invalid W certificate record")
    if "seed_manifest" in record:
        for key in ("veripb_proof", "cakepb_kernel_proof"):
            artifact = record[key]
            _archive(artifact["relative_path"], artifact)
    else:
        if record.get("unchecked_deletion_present") is not False:
            raise CoverageViolation("W proof permits unchecked deletion")
        proof = record["veripb_proof"]
        kernel = record["cakepb_kernel_proof"]
        _archive(proof["relative_gzip_path"], proof, nested=True)
        _archive(kernel["relative_gzip_path"], kernel, nested=True)
    return record["source_index"]


def _check_hierarchy_record(record: dict, ramsey_sources: set[int],
                            expected_w: dict[int, int],
                            expected_signed: dict[int, int]) -> int:
    source = record.get("source_index")
    if (record.get("status") != "VERIPB_CAKEPB_VERIFIED"
            or type(source) is not int
            or record.get("unchecked_deletion_present") is not False):
        raise CoverageViolation("invalid hierarchical record")
    t_records = record.get("t_certificates")
    if type(t_records) is not list or len(t_records) != expected_w[source]:
        raise CoverageViolation("hierarchical W-terminal coverage mismatch")
    formula_ids = set()
    signed_terminals = 0
    for expected_index, t_record in enumerate(t_records):
        if t_record.get("w_terminal_index") != expected_index:
            raise CoverageViolation("noncontiguous W terminal indices")
        formula_id = t_record.get("w_formula_id")
        if type(formula_id) is not int or formula_id in formula_ids:
            raise CoverageViolation("invalid W terminal formula IDs")
        formula_ids.add(formula_id)
        certificate = t_record.get("certificate")
        if type(certificate) is not dict:
            raise CoverageViolation("missing T certificate")
        _archive(certificate["proof_relative_path"],
                 certificate["proof_gzip"])
        _archive(certificate["kernel_relative_path"],
                 certificate["kernel_gzip"])
        terminals = t_record.get("ramsey_terminals")
        if source in ramsey_sources:
            if type(terminals) is not list:
                raise CoverageViolation("missing Ramsey terminal witnesses")
            for terminal in terminals:
                if (terminal.get("kind") not in ("K", "I")
                        or type(terminal.get("vertices")) is not list
                        or len(terminal["vertices"]) != 5):
                    raise CoverageViolation("invalid Ramsey terminal witness")
            signed_terminals += len(terminals)
        elif terminals != []:
            raise CoverageViolation("signed-negative T proof has Ramsey terminals")
    if source in ramsey_sources and signed_terminals != expected_signed[source]:
        raise CoverageViolation("Ramsey terminal aggregate mismatch")
    w_certificate = record.get("w_certificate")
    if type(w_certificate) is not dict:
        raise CoverageViolation("missing final W certificate")
    _archive(w_certificate["proof_relative_path"],
             w_certificate["proof_gzip"])
    _archive(w_certificate["kernel_relative_path"],
             w_certificate["kernel_gzip"])
    sidecar = record.get("w_terminal_sidecar")
    _archive(sidecar["relative_path"], sidecar)
    return source


def verify(w_path: Path, hierarchy_path: Path) -> dict:
    for path, digest in EXPECTED.items():
        if not path.is_file() or _sha256(path) != digest:
            raise CoverageViolation(f"input artifact pin mismatch: {path}")
    w_census = _load(W_CENSUS)
    signed_census = _load(SIGNED_CENSUS)
    ramsey_census = _load(RAMSEY_CENSUS)
    w_targets = {
        record["source_index"]
        for record in w_census["census"]["candidates"]
        if not record["w_square_satisfiable"]
    }
    deeper_targets = {
        record["source_index"]
        for record in w_census["census"]["candidates"]
        if record["w_square_satisfiable"]
    }
    ramsey_sources = {
        record["source_index"] for record in ramsey_census["census"]["candidates"]
    }
    expected_w = {
        record["source_index"]: record["w_completions_tested"]
        for record in signed_census["census"]["candidates"]
    }
    expected_w.update({
        record["source_index"]: record["w_completions_tested"]
        for record in ramsey_census["census"]["candidates"]
    })
    expected_signed = {
        record["source_index"]: record["signed_completions_tested"]
        for record in ramsey_census["census"]["candidates"]
    }
    w_document = _load(w_path)
    hierarchy = _load(hierarchy_path)
    if (w_document.get("coverage", {}).get("campaign_complete") is not True
            or hierarchy.get("coverage", {}).get("campaign_complete") is not True):
        raise CoverageViolation("certificate campaigns are not complete")
    w_sources = [_check_w_record(record) for record in w_document["records"]]
    deeper_sources = [
        _check_hierarchy_record(
            record, ramsey_sources, expected_w, expected_signed)
        for record in hierarchy["records"]
    ]
    if set(w_sources) != w_targets or len(w_sources) != len(w_targets):
        raise CoverageViolation("W certificate target partition mismatch")
    if (set(deeper_sources) != deeper_targets
            or len(deeper_sources) != len(deeper_targets)):
        raise CoverageViolation("hierarchical target partition mismatch")
    if w_targets & deeper_targets or len(w_targets | deeper_targets) != 705:
        raise CoverageViolation("certificate partitions do not cover 705 supports")
    return {
        "w_negative_certificates": len(w_sources),
        "hierarchical_certificates": len(deeper_sources),
        "support_orbits_certified": len(w_sources) + len(deeper_sources),
        "fixed_five_branch_certificate_complete": True,
        "general_ramsey_bound_claimed": False,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--w-certificates", type=Path,
                        default=DEFAULT_W_CERTIFICATES)
    parser.add_argument("--hierarchy", type=Path, default=DEFAULT_HIERARCHY)
    args = parser.parse_args(argv)
    print(json.dumps(
        verify(args.w_certificates, args.hierarchy), sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
