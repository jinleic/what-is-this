#!/usr/bin/env python3
"""Bounded local inspection of event-semantics partial artifacts. No network.

Reads the 2 GDELT export slices already on disk (one canonical under data/raw,
one unregistered under explorations/event-semantics/data/gdelt), verifies the
canonical file's upstream MD5 (GDELT publishes upstream MD5 as ETag), and emits
the timing/coverage facts the audit needs. Read-only; prints JSON to stdout.
"""
import csv
import io
import json
import os
import zipfile
import zlib

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FILES = [
    ("exploration_unregistered",
     os.path.join(ROOT, "explorations", "event-semantics", "data", "gdelt",
                  "20220809030000.export.CSV.zip"),
     {"upstream_md5": None, "upstream_bytes": None,
      "manifest_path": None, "url": None, "period": "20220809030000"}),
    ("canonical_manifest_sample",
     os.path.join(ROOT, "data", "raw", "gdelt", "20260829231500.export.CSV.zip"),
     {"upstream_md5": "56a02e3b605f81fbf2e440b2cc676496", "upstream_bytes": 29549,
      "manifest_path": "data/raw/gdelt/20260829231500.export.CSV.zip",
      "url": "https://data.gdeltproject.org/gdeltv2/20260829231500.export.CSV.zip",
      "period": "20260829231500"}),
]

# GDELT 2.0 Event Export v2.0 column map (61 cols) for the fields the audit uses.
COLS = {
    "GLOBALEVENTID": 0,
    "SQLDATE": 1,            # event date (day, UTC)
    "Actor1Code": 6,
    "Actor1Name": 7,
    "Actor2Code": 16,
    "IsRootEvent": 26,
    "EventCode": 27,         # CAMEO with subclass decor if quad-class used
    "QuadClass": 29,
    "GoldsteinScale": 30,
    "NumMentions": 31,
    "NumSources": 32,
    "NumArticles": 33,
    "AvgTone": 34,
    "Actor1Type1Code": 35,
    "Actor1Geo_CountryCode": 37,
    "Actor2Geo_CountryCode": 55,
    "DATEADDED": 59,         # ingest minute timestamp, UTC
    "SOURCEURL": 60,
}


def md5_fp(chunks):
    try:
        import hashlib
        return hashlib.md5()
    except ImportError:
        return None


def inspect(tag, path, expected):
    out = {"tag": tag, "path": os.path.relpath(path, ROOT)}
    out["exists"] = os.path.exists(path)
    if not out["exists"]:
        return out
    size = os.path.getsize(path)
    out["zip_bytes"] = size
    out["expected_upstream_bytes"] = expected["upstream_bytes"]
    out["bytes_match"] = (expected["upstream_bytes"] is None
                          or size == expected["upstream_bytes"])
    md5 = md5_fp(None)
    with zipfile.ZipFile(path) as zf:
        members = zf.namelist()
        payloads = {}
        for m in members:
            with zf.open(m) as fh:
                raw = fh.read()
                if md5 is not None:
                    md5.update(raw)
                payloads[m] = raw.decode("utf-8", errors="replace")
    digest_hex = md5.hexdigest() if md5 is not None else None
    out["member"] = members[0] if len(members) == 1 else members
    out["computed_md5_of_csv"] = digest_hex
    out["expected_upstream_md5"] = expected["upstream_md5"]
    out["md5_match"] = (expected["upstream_md5"] is None
                        or digest_hex == expected["upstream_md5"])
    rows = [r for r in csv.reader(io.StringIO("".join(payloads.values())),
                                  delimiter="\t") if r]
    out["rows"] = len(rows)
    out["cols_per_row"] = sorted({len(r) for r in rows})
    # DATEADDED: minute of ingest
    da = sorted({r[COLS["DATEADDED"]] for r in rows})
    out["dateadded"] = {"min": da[0], "max": da[-1],
                        "distinct": len(da), "values": da}
    # SQLDATE: event-day stated in the record
    days = {}
    for r in rows:
        days[r[COLS["SQLDATE"]]] = days.get(r[COLS["SQLDATE"]], 0) + 1
    out["day_record"] = dict(sorted(days.items()))
    # Backward event-day spread relative to ingest minute
    def day_of(ts):
        return ts[:8]
    ingest_day = day_of(da[0])
    lag = {}
    dd_max = 0
    missing_dates = 0
    for r in rows:
        d = r[COLS["SQLDATE"]]
        if not d:
            missing_dates += 1
            continue
        delta = (int(ingest_day) - int(d)) // 1
        dd = (int(ingest_day[:4]) * 10000 + int(ingest_day[4:6]) * 100
              + int(ingest_day[6:8])) - int(d)
        dd_days = dd / 1  # coarse day-digit comparison handles same-month drift only
        lag_key = "same_day" if d == ingest_day else "earlier_day"
        lag[lag_key] = lag.get(lag_key, 0) + 1
        if dd > dd_max:
            dd_max = dd
    out["event_day_relative_to_ingest"] = {
        **{k: v for k, v in lag.items()},
        "max_day_digit_delta_in_month": dd_max,
        "missing_sqldate_rows": missing_dates,
    }
    # Actor observability
    a1 = sum(1 for r in rows if r[COLS["Actor1Code"]].strip())
    a2 = sum(1 for r in rows if r[COLS["Actor2Code"]].strip())
    out["actors"] = {"actor1_present": a1, "actor2_present": a2,
                     "both_present": sum(1 for r in rows
                                         if r[COLS["Actor1Code"]].strip()
                                         and r[COLS["Actor2Code"]].strip())}
    # Timestamp plausibility: DATEADDED inside canonical 15-min bucket
    minute_prefix = expected["period"][:12] if expected["period"] else None
    out["dateadded_in_period_minute"] = (
        all(v[:12] == minute_prefix for v in da) if minute_prefix else None)
    # Dedup-relevant identity collisions
    keys = {}
    for r in rows:
        k = (r[COLS["DATEADDED"]], r[COLS["SOURCEURL"]].strip())
        keys[k] = keys.get(k, 0) + 1
    out["dedup_key_dateadded_sourceurl"] = {
        "rows": len(rows), "distinct_keys": len(keys),
        "collided_keys": sum(1 for v in keys.values() if v > 1)}
    ids = {}
    for r in rows:
        ids[r[COLS["GLOBALEVENTID"]]] = ids.get(r[COLS["GLOBALEVENTID"]], 0) + 1
    out["globaleventid"] = {"distinct": len(ids),
                            "duplicates": sum(1 for v in ids.values() if v > 1)}
    # Event-code spread (structure kind for the audit, no returns)
    codes = {}
    for r in rows:
        c = r[COLS["EventCode"]]
        codes[c] = codes.get(c, 0) + 1
    out["eventcode_top"] = dict(sorted(codes.items(),
                                       key=lambda kv: -kv[1])[:8])
    return out


for tag, path, expected in FILES:
    print(json.dumps(inspect(tag, path, expected), indent=2, sort_keys=True))
