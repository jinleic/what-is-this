"""Campaign runner: immutable snapshots under campaigns/<UTC>_<uuid8>_<hash12>/.

For each declared (arms x points) run the runner:
  1. names the campaign dir <ISO-UTC>_<uuid8>_<configsha12>,
  2. freezes a manifest.json BEFORE decoding (config, seeds, versions,
     circuit shas, DEM derivation),
  3. writes results incrementally as gz-JSON,
  4. closes with summary.json (LERs, CIs, timing percentiles) and an
     inventory of artifacts; the dir is never edited after close.

Naming contract shared with the fss-bb target (2026-08-29): same shape,
per-target directory.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import stim

TARGET_ROOT = Path(__file__).resolve().parent.parent.parent  # qldpc-dec/
CAMPAIGNS_DIR = TARGET_ROOT / "campaigns"


def _nice_bounded(cmd: list[str]) -> list[str]:
    return ["nice", "-n", "10", *cmd]


def config_hash(config: dict) -> str:
    blob = json.dumps(
        config, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(blob).hexdigest()[:12]


def _json_text(value: object, *, indent: int) -> str:
    return json.dumps(value, indent=indent, allow_nan=False) + "\n"


class Campaign:
    def __init__(self, config: dict, root: Path | None = None, dry_name: str | None = None):
        self.config = config
        self.hash12 = config_hash(config)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        uuid8 = uuid.uuid4().hex[:8]
        base = root or CAMPAIGNS_ROOT()
        self.dir = base / f"{stamp}_{uuid8}_{self.hash12}"
        if dry_name:
            self.dir = base / dry_name
        self.dir.mkdir(parents=True, exist_ok=False)
        self.results_path = self.dir / "results.json.gz"
        self._results: list[dict] = []
        self._manifest_written = False
        self._closed = False

    def _ensure_open(self) -> None:
        if self._closed:
            raise RuntimeError(f"campaign is already closed: {self.dir}")


    def write_manifest(self, extra: dict | None = None):
        self._ensure_open()
        if self._manifest_written:
            raise RuntimeError(f"campaign manifest is already frozen: {self.dir}")
        man = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "platform": platform.platform(),
            "python": sys.version,
            "config": self.config,
            "config_hash12": self.hash12,
        }
        if extra:
            man.update(extra)
        (self.dir / "manifest.json").write_text(_json_text(man, indent=2))
        self._manifest_written = True
        return man

    def append_result(self, result: dict):
        self._ensure_open()
        if not self._manifest_written:
            raise RuntimeError("write the campaign manifest before appending results")
        updated = [*self._results, result]
        payload = _json_text(updated, indent=1)
        compressed = gzip.compress(payload.encode(), mtime=0)
        self.results_path.write_bytes(compressed)
        self._results = updated

    def close(self, summary: dict):
        self._ensure_open()
        if not self._manifest_written:
            raise RuntimeError("write the campaign manifest before closing")
        summary_payload = _json_text(summary, indent=2)
        files = {p.name for p in self.dir.iterdir() if p.name != "inventory.json"}
        files.update(("summary.json", "inventory.json"))
        inv = {
            "files": sorted(files),
            "closed_utc": datetime.now(timezone.utc).isoformat(),
            "results_count": len(self._results),
            "policy": "immutable-after-close",
        }
        inventory_payload = _json_text(inv, indent=2)
        (self.dir / "summary.json").write_text(summary_payload)
        (self.dir / "inventory.json").write_text(inventory_payload)
        self._closed = True


def CAMPAIGNS_ROOT() -> Path:
    return CAMPAIGNS_DIR


if __name__ == "__main__":
    demo = Campaign({"demo": True}, dry_name="dryrun_demo")
    demo.write_manifest({"note": "smoke"})
    demo.append_result({"x": 1})
    demo.close({"ok": True})
    print("campaign dir:", demo.dir)
