#!/usr/bin/env python3
"""Enumerate all exact linear W row domains for the 705 fixed-five survivors."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
from pathlib import Path

import involution_f5_signed_completion as completion

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_w_domain_census"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_w_domain_census.json"




def _signature_bytes(signature: tuple[tuple[int, ...], tuple[int, ...]]) -> bytes:
    return bytes(signature[0]) + bytes(signature[1])


def build_analysis() -> dict:
    records = []
    all_counts = []
    candidate_digest = hashlib.sha256()
    labelled_digest = hashlib.sha256()
    residual_histogram = collections.Counter()
    orbit_histogram = collections.Counter()

    for instance in completion.load_all_instances():
        row_counts = []
        domain_digest = hashlib.sha256()
        for row in range(20):
            domains = completion.enumerate_w_row_domains(instance, row)
            row_counts.append(len(domains))
            all_counts.append(len(domains))
            domain_digest.update(row.to_bytes(1, "big"))
            domain_digest.update(len(domains).to_bytes(4, "big"))
            for domain in domains:
                domain_digest.update(
                    completion.w_domain_code(domain).to_bytes(4, "big"))
        d5_images = completion.explicit_d5_images(instance)
        for signature in d5_images:
            labelled_digest.update(_signature_bytes(signature))
        residual_group_order = 1 << len(completion.residual_generators(instance))
        residual_histogram[residual_group_order] += 1
        orbit_histogram[instance.orbit_size] += 1
        record = {
            "source_index": instance.source_index,
            "orbit_size": instance.orbit_size,
            "row_domain_counts": row_counts,
            "row_domains_sha256": domain_digest.hexdigest(),
            "d5_image_count": len(d5_images),
            "residual_identical_column_group_order": residual_group_order,
        }
        encoded = json.dumps(
            record, sort_keys=True, separators=(",", ":")).encode()
        candidate_digest.update(encoded)
        candidate_digest.update(b"\n")
        records.append(record)

    if len(records) != 705 or sum(record["orbit_size"] for record in records) != 6627:
        raise completion.CompletionViolation("unexpected exact survivor census")
    if len(all_counts) != 14_100 or not all(all_counts):
        raise completion.CompletionViolation("a survivor has an empty linear W row domain")

    total_domains = sum(all_counts)
    disposition = (
        "F5_W_LINEAR_ROW_DOMAINS_EXACT_"
        f"{total_domains}_ACROSS_705_D5_SUPPORT_ORBITS"
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": disposition,
        "input_dependencies": {
            "ramsey_frontier": {
                "relative_path": "r55/data/involution_f5_ramsey_filter.json",
                "bytes": completion.FRONTIER_BYTES,
                "sha256": completion.FRONTIER_SHA256,
            },
            "support_census": {
                "relative_path": "r55/data/involution_f5_support_census.json",
                "bytes": completion.SUPPORT_BYTES,
                "sha256": completion.SUPPORT_SHA256,
            },
        },
        "domain_model": {
            "matrix": "symmetric W with shared unordered-pair relation states",
            "entry_domain": "W_ij/2 in {-1,0,+1}",
            "constraints": [
                "Q_F R + R W = -J",
                "exact W=+2 and W=-2 row degrees",
                "recorded four-vertex R(3,3) relation restrictions",
            ],
            "w_square_enforced_here": False,
            "t_square_enforced_here": False,
        },
        "census": {
            "support_orbits": len(records),
            "labelled_supports_after_explicit_d5_expansion": 6627,
            "rows": len(all_counts),
            "total_row_domains": total_domains,
            "minimum_row_domains": min(all_counts),
            "maximum_row_domains": max(all_counts),
            "orbit_size_histogram": {
                str(key): orbit_histogram[key] for key in (1, 5, 10)
            },
            "residual_group_order_histogram": {
                str(key): residual_histogram[key] for key in (1, 4, 16, 64)
            },
            "explicit_d5_images_sha256": labelled_digest.hexdigest(),
            "candidate_records_sha256": candidate_digest.hexdigest(),
            "candidates": records,
        },
        "next_gate": {
            "w_square": "OPEN_EXACT_W_SQUARE_COMPLETION",
            "t_square": "OPEN_EXACT_T_SQUARE_COMPLETION",
            "ramsey_graph_check": "OPEN",
            "general_ramsey_bound_claimed": False,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    document = build_analysis()
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "disposition": document["disposition"],
        "support_orbits": document["census"]["support_orbits"],
        "total_row_domains": document["census"]["total_row_domains"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
