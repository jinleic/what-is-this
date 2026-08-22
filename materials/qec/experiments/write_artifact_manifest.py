"""Generate the canonical research-artifact integrity manifest.

Only durable, claim-bearing outputs are included.  Partial and quarantined runs
are deliberately excluded and remain discoverable under their named directories.
"""
from __future__ import annotations

import argparse
import hashlib
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "manifest.yaml"

CANONICAL_ROOTS = (
    Path("results/certificates"),
    Path("results/processed"),
    Path("results/raw"),
    Path("proofs"),
    Path("reports"),
    Path("sources"),
)
CANONICAL_FILES = (
    Path("README.md"),
    Path("LICENSE"),
    Path("CITATION.cff"),
    Path("pyproject.toml"),
    Path("uv.lock"),
    Path("third_party/manifest.yaml"),
    Path("notes/open_status.md"),
    Path("notes/novelty_matrix.md"),
    Path("notes/literature_matrix.md"),
    Path("notes/hypotheses.csv"),
    Path("notes/failed_routes.md"),
    Path("checkpoints/verified_results.json"),
)
EXCLUDED_TREES = (
    "results/partial_runs/",
    "results/quarantine/",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_class(relative: Path) -> str:
    if relative.parts[:2] == ("results", "certificates"):
        return "certificate"
    if relative.parts[:2] == ("results", "raw"):
        return "canonical_raw_data"
    if relative.parts[:2] == ("results", "processed"):
        return "processed_result"
    if relative.parts and relative.parts[0] == "proofs":
        return "proof"
    if relative.parts and relative.parts[0] == "reports":
        return "report"
    if relative.parts and relative.parts[0] == "sources":
        return "source_provenance"
    if relative == Path("third_party/manifest.yaml"):
        return "software_provenance"
    if relative.parts and relative.parts[0] == "notes":
        return "research_ledger"
    if relative == Path("checkpoints/verified_results.json"):
        return "verified_claim_ledger"
    return "repository_metadata"


def canonical_paths() -> list[Path]:
    paths: set[Path] = set()
    for relative in CANONICAL_FILES:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"required canonical artifact is missing: {relative}")
        paths.add(path)
    for relative_root in CANONICAL_ROOTS:
        root = ROOT / relative_root
        if not root.is_dir():
            raise FileNotFoundError(f"required canonical directory is missing: {relative_root}")
        paths.update(path for path in root.rglob("*") if path.is_file())
    empty = [path.relative_to(ROOT) for path in paths if path.stat().st_size == 0]
    if empty:
        rendered = ", ".join(path.as_posix() for path in sorted(empty))
        raise ValueError(f"canonical artifacts must be nonempty: {rendered}")
    return sorted(paths, key=lambda path: path.relative_to(ROOT).as_posix())


def build_manifest(checked: str) -> dict:
    entries = []
    for path in canonical_paths():
        relative = path.relative_to(ROOT)
        entries.append(
            {
                "path": relative.as_posix(),
                "artifact_class": artifact_class(relative),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return {
        "schema_version": 1,
        "date_checked": checked,
        "hash_algorithm": "SHA-256",
        "scope": "canonical claim-bearing repository artifacts",
        "canonical_only": True,
        "excluded_trees": list(EXCLUDED_TREES),
        "entries": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default=date.today().isoformat())
    args = parser.parse_args()
    payload = build_manifest(args.date)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".yaml.tmp")
    temporary.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    temporary.replace(OUTPUT)
    print(f"wrote {OUTPUT.relative_to(ROOT)} with {len(payload['entries'])} entries")


if __name__ == "__main__":
    main()
