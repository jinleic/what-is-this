"""Integrity checks for the canonical artifact manifest.

Two admissible states:
- complete manifest: every entry hashes and sizes exactly, no partial/quarantine
  paths, no zero-byte artifacts;
- truthful blocked placeholder: ``complete: false`` naming the zero-byte
  canonical files that prevent certification, each of which must actually be
  empty on disk (a placeholder that hides ready artifacts is a defect too).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "artifacts" / "manifest.yaml"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_manifest_is_complete_and_valid_or_a_truthful_blocked_placeholder() -> None:
    payload = yaml.safe_load(MANIFEST.read_text())
    assert payload["schema_version"] == 1
    assert payload["canonical_only"] is True

    if payload.get("complete", True):
        assert payload["hash_algorithm"] == "SHA-256"
        assert set(payload["excluded_trees"]) == {
            "results/partial_runs/",
            "results/quarantine/",
        }
        entries = payload["entries"]
        paths = [entry["path"] for entry in entries]
        assert paths == sorted(paths)
        assert len(paths) == len(set(paths))
        assert paths
        assert not any(
            "partial_runs" in path or "quarantine" in path for path in paths
        )
        for entry in entries:
            path = ROOT / entry["path"]
            assert path.is_file(), entry["path"]
            assert entry["bytes"] > 0, f"zero-byte canonical artifact: {entry['path']}"
            assert path.stat().st_size == entry["bytes"], entry["path"]
            assert sha256(path) == entry["sha256"], entry["path"]
    else:
        assert payload["status"] == "BLOCKED_PENDING_CANONICAL_OUTPUTS"
        blocking = payload["blocking_empty_files"]
        assert blocking, "blocked placeholder must name at least one blocking file"
        for rel in blocking:
            path = ROOT / rel
            assert not path.is_file() or path.stat().st_size == 0, (
                f"{rel} is nonempty; placeholder is stale - regenerate the manifest"
            )
        assert "entries" not in payload, (
            "placeholder must not carry hash entries; regenerate instead"
        )
