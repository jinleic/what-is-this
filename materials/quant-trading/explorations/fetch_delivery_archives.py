#!/usr/bin/env python3
"""Fetch Binance Data Vision dated delivery-futures 1h kline archives named in a plan.

Keyless public GETs only. Every archive is stored next to its upstream .CHECKSUM
sidecar under quant-trading/data/raw/binance/delivery-1h/<margin>/<SYMBOL>/ and
verified (sha256 of the zip == sidecar digest, size == listing size). Existing
archives that already verify are reused, never re-downloaded. Nothing inside any
archive is parsed here; this script acquires bytes, it computes no outcome.

Every path comes from the plan record, so the plan is the single source of truth
for what is fetched and where it lands.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import time
import urllib.request
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
BASE_URL = "https://data.binance.vision/"
UA = {"User-Agent": "quant-frontier-inventory/1.0"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str, dest: Path, attempts: int = 5) -> int:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as response:
                data = response.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return len(data)
        except Exception as exc:  # noqa: BLE001 - retried, then reported
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"{url}: {last}")


def sidecar_digest(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"empty checksum sidecar {path}")
    return text.split()[0].lower()


def process(entry: dict) -> dict:
    size = entry["size_bytes"]
    zip_path = WORKSPACE_ROOT / entry["path"]
    sidecar = WORKSPACE_ROOT / entry["checksum_sidecar"]
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    record = {"margin": entry["margin"], "symbol": entry["symbol"], "month": entry["month"],
              "path": entry["path"], "checksum_path": entry["checksum_sidecar"], "listing_size": size}
    downloaded = False
    try:
        if not sidecar.exists():
            fetch(BASE_URL + entry["key"] + ".CHECKSUM", sidecar)
            downloaded = True
        expected = sidecar_digest(sidecar)
        ok = zip_path.exists() and zip_path.stat().st_size == size and sha256_file(zip_path) == expected
        if not ok:
            fetch(BASE_URL + entry["key"], zip_path)
            downloaded = True
            observed = sha256_file(zip_path)
            if observed != expected:
                raise ValueError(f"sha256 mismatch {entry['path']}: {observed} != sidecar {expected}")
            if zip_path.stat().st_size != size:
                raise ValueError(f"size mismatch {entry['path']}: {zip_path.stat().st_size} != listing {size}")
        record.update({"status": "verified", "downloaded": downloaded, "sha256": expected,
                       "checksum_file_sha256": sha256_file(sidecar), "size_bytes": zip_path.stat().st_size})
    except Exception as exc:  # noqa: BLE001 - recorded, never hidden
        record.update({"status": "failed", "downloaded": downloaded, "error": str(exc)})
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args(argv)
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    entries = plan["records"]
    started = time.time()
    records: list[dict] = []
    with concurrent.futures.ThreadPoolExecutor(args.threads) as pool:
        for index, record in enumerate(pool.map(process, entries), 1):
            records.append(record)
            if index % 200 == 0 or index == len(entries):
                failed = sum(1 for r in records if r["status"] != "verified")
                print(f"{index}/{len(entries)} verified={index - failed} failed={failed} "
                      f"{time.time() - started:.0f}s", flush=True)
    failed = [r for r in records if r["status"] != "verified"]
    log = {"plan": str(args.plan), "plan_id": plan.get("plan_id"), "files": len(records),
           "verified": len(records) - len(failed), "failed": len(failed),
           "downloaded": sum(1 for r in records if r.get("downloaded")),
           "bytes_verified": sum(r.get("size_bytes", 0) for r in records if r["status"] == "verified"),
           "seconds": round(time.time() - started, 1), "records": records}
    args.log.write_text(json.dumps(log, indent=1) + "\n", encoding="utf-8")
    print(f"done: verified {log['verified']} failed {log['failed']} bytes {log['bytes_verified']} -> {args.log}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
