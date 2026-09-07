#!/usr/bin/env python3
"""Regenerate the n=12 target-39 family and verify byte identity and cover."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KOBON = HERE.parent
WORKSPACE = KOBON.parents[1]
MANIFEST = HERE / "manifest.json"
RESULT = HERE / "provenance.json"

sys.path[:0] = [str(HERE), str(KOBON)]
import engine  # noqa: E402
import make_manifest  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, document: dict) -> None:
    text = json.dumps(document, indent=2, sort_keys=True) + "\n"
    temporary = path.with_name(f"{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def orbit(triple: tuple[int, int, int], n: int) -> set[tuple[int, ...]]:
    return {
        tuple(sorted(((-value if reflected else value) + shift) % n
                     for value in triple))
        for shift in range(n)
        for reflected in (False, True)
    }


def verify_cover(n: int) -> dict:
    representatives = engine.dihedral_orbits_3(n)
    seen: set[tuple[int, ...]] = set()
    disjoint = True
    for representative in representatives:
        current = orbit(representative, n)
        if seen & current:
            disjoint = False
        seen |= current
    expected = set(itertools.combinations(range(n), 3))
    return {
        "representatives": [list(item) for item in representatives],
        "representative_count": len(representatives),
        "covered_triples": len(seen),
        "expected_triples": len(expected),
        "pairwise_disjoint_orbits": disjoint,
        "exact_cover": disjoint and seen == expected,
    }


def verify(output_dir: Path, *, reuse_generated: bool = False) -> dict:
    manifest_bytes = MANIFEST.read_bytes()
    current_manifest_bytes = (
        json.dumps(make_manifest.build(), indent=2, sort_keys=True) + "\n"
    ).encode()
    if manifest_bytes != current_manifest_bytes:
        raise RuntimeError(
            "manifest.json is stale for the current engine, CNFs, runner, "
            "supervisor, or resolved Kissat binary")
    manifest = json.loads(manifest_bytes)
    manifest_digest = hashlib.sha256(manifest_bytes).hexdigest()
    expected_instances = {
        record["cube"]: WORKSPACE / record["path"]
        for record in manifest["instances"]
    }
    if set(expected_instances) != set(manifest["case_order"]):
        raise RuntimeError("manifest instance set disagrees with case_order")
    if reuse_generated:
        prefix = output_dir / "kobon_n12_t39_proof"
        generated = {
            "base": str(prefix) + ".cnf",
            **{cube: f"{prefix}_{cube}.cnf"
               for cube in manifest["case_order"]},
        }
    else:
        if output_dir.exists() and any(output_dir.iterdir()):
            raise RuntimeError(
                f"refusing to overwrite nonempty output directory: {output_dir}")
        output_dir.mkdir(parents=True, exist_ok=True)
        generated = engine.dump_instances(12, 39, str(output_dir), cubes=True)
    comparisons = []
    for cube in manifest["case_order"]:
        actual = Path(generated[cube])
        expected = expected_instances[cube]
        actual_hash = sha256(actual)
        expected_hash = sha256(expected)
        comparisons.append({
            "cube": cube,
            "generated_path": str(actual.relative_to(WORKSPACE)),
            "generated_bytes": actual.stat().st_size,
            "generated_sha256": actual_hash,
            "expected_path": str(expected.relative_to(WORKSPACE)),
            "expected_bytes": expected.stat().st_size,
            "expected_sha256": expected_hash,
            "byte_identical": (actual.stat().st_size == expected.stat().st_size
                               and actual_hash == expected_hash),
        })
    base = Path(generated["base"])
    canonical_base = WORKSPACE / "scratch/kobon/n12/kobon_n12_t39_proof.cnf"
    base_record = {
        "generated_path": str(base.relative_to(WORKSPACE)),
        "generated_bytes": base.stat().st_size,
        "generated_sha256": sha256(base),
        "canonical_path": str(canonical_base.relative_to(WORKSPACE)),
        "canonical_bytes": canonical_base.stat().st_size,
        "canonical_sha256": sha256(canonical_base),
    }
    base_record["byte_identical"] = (
        base_record["generated_bytes"] == base_record["canonical_bytes"]
        and base_record["generated_sha256"] == base_record["canonical_sha256"])
    cover = verify_cover(12)
    passed = (all(item["byte_identical"] for item in comparisons)
              and base_record["byte_identical"]
              and cover["exact_cover"]
              and len(comparisons) == 15)
    verifier_path = Path(__file__).resolve()
    return {
        "schema": "kobon-n12-provenance/2",
        "status": "PASS" if passed else "FAIL",
        "manifest_sha256": manifest_digest,
        "manifest_prevalidated_against_current_inputs": True,
        "generation_mode": (
            "REUSED_PRIOR_FRESH_GENERATION"
            if reuse_generated else "FRESH_REGENERATION"),
        "verifier": {
            "path": str(verifier_path.relative_to(WORKSPACE)),
            "bytes": verifier_path.stat().st_size,
            "sha256": sha256(verifier_path),
        },
        "engine": manifest["engine"],
        "claim_scope": "byte regeneration of the existing CNF family plus exact D12 triple-orbit cover",
        "base": base_record,
        "cases": comparisons,
        "cover": cover,
        "limits": [
            "This check proves byte identity and the finite D12 triple-orbit partition.",
            "It does not replace the clause-soundness argument for engine.dump_instances.",
            "Solver UNSAT remains discovery-only until complete DRAT proofs verify.",
        ],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument(
        "--reuse-generated", action="store_true",
        help="recheck a prior fresh output directory without overwriting it")
    args = parser.parse_args(argv)
    manifest_digest = hashlib.sha256(MANIFEST.read_bytes()).hexdigest()
    output_dir = (args.output_dir or
                  WORKSPACE / "scratch/kobon/n12" /
                  f"regenerated_{manifest_digest[:16]}")
    document = verify(output_dir, reuse_generated=args.reuse_generated)
    atomic_json(args.result, document)
    print(json.dumps({
        "case_count": len(document["cases"]),
        "cover_exact": document["cover"]["exact_cover"],
        "manifest_sha256": document["manifest_sha256"],
        "status": document["status"],
    }, sort_keys=True))
    return 0 if document["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
