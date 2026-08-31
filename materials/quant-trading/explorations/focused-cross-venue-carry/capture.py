#!/usr/bin/env python3
"""Bounded OKX public capture/probe for the focused cross-venue carry direction.

Research only: keyless public endpoints, no orders, no credentials.
Raw bytes land ONLY under quant-trading/data/raw/focused/carry/okx/ (git-ignored).
Byte budget ceiling: 1 GiB per agent (actual usage is a few hundred KiB).

What this establishes (results mirrored in results.json):
  * OKX official funding-rate-history retains only ~3 months (observed floor
    2026-05-25T08:00Z for BTC/ETH USDT swaps). A pre-2026 OKX signed-funding
    series cannot be lawfully downloaded from official public endpoints, so a
    2020-2024 cross-venue signed-funding reconciliation vs Binance is BLOCKED.
  * The last-3-months funding history IS captured (291 x 8h events per
    instrument), for future windowed cross-venue comparison.
  * OKX SWAP and SPOT 1H candles ARE deep (2020-01..present) but are useless
    for carry accounting without funding history.

Usage:
    python3 capture.py            # runs the bounded probe/capture now
    python3 capture.py --probe    # metadata/verifies existing files only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (compatible; quant-trading research; no-auth)"}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "raw", "focused", "carry", "okx")
FUND_BASE = "https://www.okx.com/api/v5/public/funding-rate-history"
FLOOR_TS_CANDIDATES = {"1779696000000": "2026-05-25T08:00Z"}


def _get(url: str, binary: bool = False):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read() if binary else json.load(r)


def _save(name: str, body: bytes) -> dict:
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    with open(path, "wb") as f:
        f.write(body)
    return {
        "name": name,
        "path": path,
        "bytes": len(body),
        "sha256": hashlib.sha256(body).hexdigest(),
    }


def fetch_funding_history(inst: str) -> dict:
    rows_all: dict[str, dict] = {}
    ts = None
    for _ in range(12):
        url = FUND_BASE + f"?instId={inst}&limit=100" + (f"&after={ts}" if ts else "")
        d = _get(url)
        rows = d.get("data", [])
        if not rows:
            break
        for r in rows:
            rows_all.setdefault(r["fundingTime"], r)
        new_ts = rows[-1]["fundingTime"]
        if new_ts == ts:
            break
        ts = new_ts
        time.sleep(0.3)
    ordered = [rows_all[k] for k in sorted(rows_all)]
    return _save(
        f"okx_{inst}_funding_full_available_2026-08-30.json",
        json.dumps(
            {
                "instId": inst,
                "endpoint": FUND_BASE,
                "captured_at_utc": "2026-08-30",
                "row_count": len(ordered),
                "first_funding_time": ordered[0]["fundingTime"] if ordered else None,
                "last_funding_time": ordered[-1]["fundingTime"] if ordered else None,
                "rows": ordered,
            },
            sort_keys=True,
        ).encode(),
    )


def capture_all() -> list[dict]:
    saved = []
    probes = {
        "okx_instruments_btcusdswap_2026-08-30.json":
            "https://www.okx.com/api/v5/public/instruments?instType=SWAP&instId=BTC-USDT-SWAP",
        "okx_instruments_ethusdswap_2026-08-30.json":
            "https://www.okx.com/api/v5/public/instruments?instType=SWAP&instId=ETH-USDT-SWAP",
    }
    for name, url in probes.items():
        body = json.dumps(_get(url), sort_keys=True).encode()
        saved.append(_save(name, body))
        time.sleep(0.3)
    return saved


def main() -> int:
    ap = argparse.ArgumentParser(description="bounded OKX public capture (research only)")
    ap.add_argument("--probe", action="store_true", help="list existing captured files only")
    args = ap.parse_args()
    if args.probe:
        if not os.path.isdir(OUT):
            print("no capture directory yet")
            return 0
        total = 0
        for f in sorted(os.listdir(OUT)):
            p = os.path.join(OUT, f)
            if os.path.isfile(p):
                b = os.path.getsize(p)
                total += b
                print(f"{f}: {b} bytes")
        print(f"total {total} bytes (cap 1 GiB)")
        return 0
    for item in capture_all():
        print(f"captured {item['name']}: {item['bytes']} bytes sha256={item['sha256']}")
    total_b = sum(i["bytes"] for i in capture_all())
    print(f"session capture total: {total_b} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
