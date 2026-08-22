#!/usr/bin/env python3
"""Standalone checks for the 3x3 Onsager-quotient corollary.

Clean-room discipline: this test never imports the producer experiment.  It
re-verifies every antecedent from disk (source hashes recomputed, consumed
artifacts' own checks re-evaluated, factor arithmetic redone from the raw
rank/dimension integers) and then re-validates the emitted artifact envelope.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "e100_oa_3x3.py"
RESULT = ROOT / "results" / "algebra_structure" / "oa_3x3.json"
STRUCTURE = ROOT / "results" / "algebra_structure" / "char0_3x3.json"
OA_AUDIT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"
MANIFEST = ROOT / "sources" / "manifest.yaml"

EXPECTED_SOURCES = {
    "date_roan2000_onsager": (
        "sources/fulltext/date_roan2000_onsager.pdf",
        "c1334b58ed30c77b1592e628d8bf26c40d52ef63b9226ecad1ddec0e96771b06",
    ),
    "roan1991_onsager": (
        "sources/fulltext/roan1991_onsager.pdf",
        "82988421997ae8fb388c79c4f46381def00ee15ecaaed04922bd7ab94c352313",
    ),
}

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"{'PASS' if passed else 'FAIL'}: {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    # 1. sources: recompute hashes, confirm manifest pins
    manifest_text = MANIFEST.read_text(encoding="utf-8")
    for key, (relpath, expected) in EXPECTED_SOURCES.items():
        actual = digest(ROOT / relpath)
        check(
            f"source_hash_{key}",
            actual == expected and expected in manifest_text,
            actual[:16],
        )

    # 2. antecedent artifacts: their own embedded checks all pass
    structure = json.loads(STRUCTURE.read_text(encoding="utf-8"))
    check(
        "structure_artifact_self_checks",
        bool(structure["checks"]) and all(c["passed"] for c in structure["checks"]),
        f"{len(structure['checks'])} checks",
    )
    oa = json.loads(OA_AUDIT.read_text(encoding="utf-8"))
    check(
        "oa_audit_artifact_self_checks",
        bool(oa["checks"]) and all(c["passed"] for c in oa["checks"]),
        f"{len(oa['checks'])} checks",
    )
    stmt = oa["data"]["levi_corollary"]["statement"]
    check(
        "external_levi_corollary_statement",
        "Levi factor" in stmt and "direct sum of sl2 copies" in stmt
        and oa["data"]["levi_corollary"]["claim_tag"] == "[THEOREM]",
    )

    # 3. factor arithmetic redone from raw integers: dim = (rank+1)^2 - 1, rank >= 2
    levi = structure["data"]["Levi_type"]
    ranks = list(levi["factor_ranks"])
    dims = list(levi["factor_dimensions"])
    check("six_factors", len(ranks) == 6 and len(dims) == 6, f"{len(ranks)} factors")
    check(
        "factor_dimension_identity",
        all(dim == (rank + 1) ** 2 - 1 for rank, dim in zip(ranks, dims)),
        f"ranks {ranks}",
    )
    check(
        "no_factor_is_sl2",
        all(rank >= 2 for rank in ranks),
        f"min rank {min(ranks)} (sl2 would be rank 1)",
    )
    check(
        "radical_is_centre_dim1",
        structure["data"]["solvable_radical"]["equals_center"] is True
        and structure["data"]["solvable_radical"]["dimension_Q"] == 1
        and structure["data"]["derived_algebra"]["semisimple"] is True,
    )
    # total-dimension consistency: 1 + sum of factor dims = certified dimension
    check(
        "dimension_accounting",
        1 + sum(dims) == structure["data"]["dimension_Q"]["value"],
        f"1 + {sum(dims)} = {structure['data']['dimension_Q']['value']}",
    )

    # 4. emitted artifact: envelope, chain fields, embedded checks
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "result_envelope",
        set(result) == {"provenance", "data", "checks"}
        and bool(result["checks"])
        and all(c["passed"] for c in result["checks"]),
        f"{len(result['checks'])} embedded checks",
    )
    thm = result["data"]["onsager_quotient_no_go_3x3"]
    check(
        "theorem_fields",
        thm["claim_tag"] == "[THEOREM]"
        and "not isomorphic to any" in thm["statement"]
        and "Dolan--Grady" in thm["statement"]
        and len(thm["proof_chain"]) == 5,
    )
    check(
        "artifact_factor_table_matches_structure",
        [row["dimension"] for row in result["data"]["levi_factors"]] == dims
        and not any(row["is_sl2"] for row in result["data"]["levi_factors"]),
    )
    check(
        "artifact_source_hashes_match_recomputation",
        all(
            result["data"]["external_sources"][key]["sha256"] == expected
            for key, (_, expected) in EXPECTED_SOURCES.items()
        ),
    )

    # 5. producer replay: exits 0 and leaves an identical mathematical payload
    replay = subprocess.run(
        [sys.executable, str(SCRIPT)],
        capture_output=True,
        text=True,
        timeout=600,
        cwd=ROOT,
    )
    check("producer_replay_exit0", replay.returncode == 0, replay.stdout.strip().splitlines()[-1] if replay.stdout else "")
    replayed = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "replay_payload_stable",
        replayed["data"] == result["data"]
        and [c["passed"] for c in replayed["checks"]] == [c["passed"] for c in result["checks"]],
        "data section identical; timestamps excluded by construction",
    )

    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {FAILURES}")
        return 1
    print("PASS test_oa_3x3: all standalone checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
