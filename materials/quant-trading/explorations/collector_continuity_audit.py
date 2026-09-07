#!/usr/bin/env python3
"""Outcome-blind continuity audit of the two append-only collectors.

Reads ONLY timestamps, record kinds, byte lengths and hashes. Never parses a
payload field that could carry an outcome (price, side, quantity, iv, ...).

Deterministic: every figure is byte-pinned (file size read first, then exactly
that many bytes are scanned), so the record is reproducible with
`head -c <bytes_pinned>` as long as the same --now instant is supplied.

    nice -n 10 python3 collector_continuity_audit.py \
        --now 2026-09-06T02:05:09Z \
        --out collector-continuity-audit-20260906.json

Resource policy: single thread, RLIMIT_CPU hard cap, no retry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent
LIQ_DIR = WORKSPACE / "data/raw/focused/live/liquidations-persistent"
DERIBIT_DIR = WORKSPACE / "data/raw/focused/live/deribit"
REGISTRATION = HERE / "prospective-forward-test-registration.json"
PROGRAM_NEGATIVE = HERE / "program_negative.json"
STATE = HERE / "state.json"
CAMPAIGNS = HERE / "campaigns"

TS_UTC_RE = re.compile(rb'"ts_utc":"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3})Z"')
KIND_RE = re.compile(rb'"kind":"([a-z_]+)"')
RECV_MS_RE = re.compile(rb'"received_ts_ms":(\d+)')
EVENT_RE = re.compile(rb'"event":"([a-z_]+)"')
CHANNELS_RE = re.compile(rb'"channels":\[([^\]]*)\]')


def parse_iso(s: str) -> datetime:
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def iso(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{dt.microsecond // 1000:03d}Z"


def ts_ms_to_dt(ms: int) -> datetime:
    return datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc)


def parse_ts_utc(b: bytes) -> datetime:
    return datetime.strptime(b.decode(), "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)


def sha256_pinned(path: Path, nbytes: int) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        remaining = nbytes
        while remaining > 0:
            chunk = fh.read(min(1 << 20, remaining))
            if not chunk:
                break
            h.update(chunk)
            remaining -= len(chunk)
    return h.hexdigest()


def read_pinned_lines(path: Path, nbytes: int):
    """Yield complete lines inside the first nbytes; report trailing partial."""
    with path.open("rb") as fh:
        data = fh.read(nbytes)
    trailing_partial = len(data) > 0 and not data.endswith(b"\n")
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    return lines, trailing_partial


def scan_liquidation_file(path: Path, gap_limit_s: float) -> dict:
    nbytes = path.stat().st_size
    lines, trailing_partial = read_pinned_lines(path, nbytes)
    kinds: dict[str, int] = {}
    data_ts: list[datetime] = []
    control_ts: list[datetime] = []
    unparsed = 0
    for line in lines:
        km = KIND_RE.search(line)
        tm = TS_UTC_RE.search(line)
        if not km or not tm:
            unparsed += 1
            continue
        kind = km.group(1).decode()
        kinds[kind] = kinds.get(kind, 0) + 1
        ts = parse_ts_utc(tm.group(1))
        (control_ts if kind == "control" else data_ts).append(ts)
    violations = 0
    max_gap_s = 0.0
    max_gap_at = None
    gaps_over_limit = []
    for prev, cur in zip(data_ts, data_ts[1:]):
        d = (cur - prev).total_seconds()
        if d < 0:
            violations += 1
        if d > max_gap_s:
            max_gap_s, max_gap_at = d, (iso(prev), iso(cur))
        if d > gap_limit_s:
            gaps_over_limit.append({"from": iso(prev), "to": iso(cur), "seconds": round(d, 3)})
    return {
        "file": path.name,
        "bytes_pinned": nbytes,
        "sha256_pinned": sha256_pinned(path, nbytes),
        "records_all": len(lines),
        "records_data": len(data_ts),
        "records_control": len(control_ts),
        "records_unparsed": unparsed,
        "kinds": dict(sorted(kinds.items())),
        "first_data_ts_utc": iso(data_ts[0]) if data_ts else None,
        "last_data_ts_utc": iso(data_ts[-1]) if data_ts else None,
        "last_any_ts_utc": iso(max(data_ts[-1:] + control_ts[-1:])) if (data_ts or control_ts) else None,
        "data_record_monotonicity_violations": violations,
        "max_data_gap_seconds": round(max_gap_s, 3),
        "max_data_gap_at": max_gap_at,
        "data_gaps_over_limit": gaps_over_limit,
        "trailing_partial_line": trailing_partial,
        "_first": data_ts[0] if data_ts else None,
        "_last": data_ts[-1] if data_ts else None,
    }


def scan_deribit_stream(path: Path, gap_limit_s: float, want_events: bool) -> dict:
    nbytes = path.stat().st_size
    out = {"bytes_pinned": nbytes, "sha256_pinned": sha256_pinned(path, nbytes) if nbytes else None}
    if nbytes == 0:
        out.update({"records": 0, "first_received_utc": None, "last_received_utc": None})
        return out
    lines, trailing_partial = read_pinned_lines(path, nbytes)
    ts: list[int] = []
    events: dict[str, int] = {}
    subscribe_acks: list[dict] = []
    unparsed = 0
    for line in lines:
        m = RECV_MS_RE.search(line)
        if not m:
            unparsed += 1
            continue
        ts.append(int(m.group(1)))
        if want_events:
            em = EVENT_RE.search(line)
            if em:
                ev = em.group(1).decode()
                events[ev] = events.get(ev, 0) + 1
                if ev == "subscribe_ack":
                    cm = CHANNELS_RE.search(line)
                    chans = [c.strip().strip('"') for c in cm.group(1).decode().split(",")] if cm else []
                    subscribe_acks.append({"received_utc": iso(ts_ms_to_dt(ts[-1])), "channels": chans})
    violations = 0
    max_gap_ms = 0
    max_gap_at = None
    gaps_over_limit = []
    for prev, cur in zip(ts, ts[1:]):
        d = cur - prev
        if d < 0:
            violations += 1
        if d > max_gap_ms:
            max_gap_ms, max_gap_at = d, (iso(ts_ms_to_dt(prev)), iso(ts_ms_to_dt(cur)))
        if d > gap_limit_s * 1000:
            gaps_over_limit.append({"from": iso(ts_ms_to_dt(prev)), "to": iso(ts_ms_to_dt(cur)), "seconds": d / 1000.0})
    out.update({
        "records": len(lines),
        "records_unparsed": unparsed,
        "first_received_utc": iso(ts_ms_to_dt(ts[0])) if ts else None,
        "last_received_utc": iso(ts_ms_to_dt(ts[-1])) if ts else None,
        "received_monotonicity_violations": violations,
        "max_gap_seconds": max_gap_ms / 1000.0,
        "max_gap_at": max_gap_at,
        "gaps_over_limit": gaps_over_limit,
        "trailing_partial_line": trailing_partial,
        "_first_ms": ts[0] if ts else None,
        "_last_ms": ts[-1] if ts else None,
    })
    if want_events:
        out["event_census"] = dict(sorted(events.items()))
        out["subscribe_acks"] = subscribe_acks
    return out


def strip_private(obj):
    if isinstance(obj, dict):
        return {k: strip_private(v) for k, v in obj.items() if not k.startswith("_")}
    if isinstance(obj, list):
        return [strip_private(v) for v in obj]
    return obj


def verify_campaign_checksums() -> dict:
    dirs = sorted(p for p in CAMPAIGNS.iterdir() if p.is_dir())
    total = ok = bad = missing = 0
    per_dir = {}
    for d in dirs:
        sums = d / "sha256s.txt"
        if not sums.exists():
            per_dir[d.name] = "NO sha256s.txt"
            continue
        n_ok = n_bad = n_missing = 0
        for raw in sums.read_text().splitlines():
            if not raw.strip():
                continue
            digest, _, name = raw.partition("  ")
            target = d / name.strip()
            if not target.exists():
                n_missing += 1
                continue
            actual = sha256_pinned(target, target.stat().st_size)
            if actual == digest.strip():
                n_ok += 1
            else:
                n_bad += 1
        per_dir[d.name] = {"ok": n_ok, "mismatch": n_bad, "missing": n_missing}
        total += n_ok + n_bad + n_missing
        ok += n_ok
        bad += n_bad
        missing += n_missing
    return {"run_dirs": len(dirs), "files": total, "ok": ok, "mismatch": bad, "missing": missing, "per_dir": per_dir}


def alpha_ledger() -> dict:
    budget = 0.05
    per_k = {k: budget / (k * (k + 1)) for k in range(1, 9)}
    spent_7 = sum(per_k[k] for k in range(1, 8))
    pn = json.loads(PROGRAM_NEGATIVE.read_text())
    st = json.loads(STATE.read_text())
    reg = json.loads(REGISTRATION.read_text())
    recorded = {
        "program_negative.counts": pn.get("counts"),
        "state.alpha_spent": st.get("alpha_spent"),
        "registration.status": reg.get("status"),
        "registration.ordinal_clause": reg["binding_constraints_if_it_is_ever_run"]["ordinal"],
    }
    checks = {
        "sum_k1_7_equals_0.04375": abs(spent_7 - 0.04375) < 1e-15,
        "share_is_0.875": abs(spent_7 / budget - 0.875) < 1e-15,
        "alpha_8_equals_registration": abs(per_k[8] - 0.0006944444444444445) < 1e-18,
        "state_alpha_spent_matches": st.get("alpha_spent", {}).get("spent_through_cycle_7") == 0.04375
        and st.get("alpha_spent", {}).get("alpha_8_if_selected") == per_k[8],
        "registration_unspent": reg.get("status", "").startswith("UNSPENT"),
    }
    return {"familywise_budget": budget, "alpha_k": {str(k): v for k, v in per_k.items()},
            "spent_through_k7": spent_7, "remaining_alpha_8": per_k[8], "recorded": recorded, "checks": checks}


def boot_time_utc() -> str | None:
    try:
        raw = subprocess.run(["sysctl", "-n", "kern.boottime"], capture_output=True, text=True, timeout=5).stdout
        m = re.search(r"sec = (\d+)", raw)
        return iso(ts_ms_to_dt(int(m.group(1)) * 1000)) if m else None
    except Exception:  # noqa: BLE001 - evidence capture only
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--now", help="audit instant, ISO8601 Z (default: current UTC)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpu-seconds", type=int, default=30)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (args.cpu_seconds, args.cpu_seconds))
    now = parse_iso(args.now) if args.now else datetime.now(timezone.utc).replace(microsecond=0)

    reg = json.loads(REGISTRATION.read_text())
    reg_sha = sha256_pinned(REGISTRATION, REGISTRATION.stat().st_size)
    horizon = reg["pre_registered_collection_horizon"]
    elig = reg["eligibility_conditions_for_the_horizon_to_be_readable"]
    h_start = parse_iso(horizon["start_utc_inclusive"])
    h_cap = parse_iso(horizon["hard_cap_utc"])
    gap_limit_s = elig["max_tolerated_single_gap_minutes"] * 60.0
    uptime_min = elig["collector_uptime_fraction_min"]

    # ---- liquidation collector -------------------------------------------------
    liq_files = sorted(LIQ_DIR.glob("liquidation-capture-*.jsonl"))
    liq = [scan_liquidation_file(p, gap_limit_s) for p in liq_files]
    seams = []
    for a, b in zip(liq, liq[1:]):
        if a["_last"] and b["_first"]:
            seams.append({"from_file": a["file"], "to_file": b["file"],
                          "seconds": round((b["_first"] - a["_last"]).total_seconds(), 3),
                          "monotonic": b["_first"] >= a["_last"]})
    liq_last = max(x["_last"] for x in liq if x["_last"])
    liq_first = min(x["_first"] for x in liq if x["_first"])
    liq_days_present = sorted(p.name[len("liquidation-capture-"):-len(".jsonl")] for p in liq_files)

    # ---- deribit collector -----------------------------------------------------
    deribit_days = {}
    deribit_last_ms = 0
    for d in sorted(DERIBIT_DIR.glob("day*")):
        entry = {}
        for stream in ("index", "trades", "events", "instruments", "ticker"):
            p = d / f"{stream}.jsonl"
            if not p.exists():
                entry[stream] = None
                continue
            if stream in ("index", "trades", "events"):
                entry[stream] = scan_deribit_stream(p, gap_limit_s, want_events=(stream == "events"))
                if entry[stream].get("_last_ms"):
                    deribit_last_ms = max(deribit_last_ms, entry[stream]["_last_ms"])
            else:
                n = p.stat().st_size
                entry[stream] = {"bytes_pinned": n, "sha256_pinned": sha256_pinned(p, n) if n else None}
        deribit_days[d.name] = entry
    deribit_last = ts_ms_to_dt(deribit_last_ms)
    ticker_census = {day: e["ticker"]["bytes_pinned"] for day, e in deribit_days.items() if e.get("ticker")}
    ticker_channel_ever_subscribed = any(
        any("ticker" in c for c in ack["channels"])
        for e in deribit_days.values() if e.get("events")
        for ack in e["events"].get("subscribe_acks", [])
    )

    # ---- horizon coverage -------------------------------------------------------
    def bytes_records_at_or_after(files: list[dict], start: datetime) -> tuple[int, int]:
        # Files are UTC-day rotated; a record at/after start can only live in a day file dated >= start.
        n_rec = n_bytes = 0
        for f in files:
            day = f["file"][len("liquidation-capture-"):-len(".jsonl")]
            if day >= start.strftime("%Y%m%d"):
                n_bytes += f["bytes_pinned"]
                n_rec += f["records_data"]
        return n_bytes, n_rec

    liq_h_bytes, liq_h_records = bytes_records_at_or_after(liq, h_start)
    deribit_h_days = [d for d in deribit_days if d[len("day"):] >= h_start.strftime("%Y%m%d")]
    last_any = max(liq_last, deribit_last)
    dead_for = now - min(liq_last, deribit_last)
    horizon_elapsed = now - h_start
    horizon_days_elapsed = [(h_start + timedelta(days=i)).strftime("%Y%m%d")
                            for i in range(int(horizon_elapsed.total_seconds() // 86400) + 1)]
    horizon_days_with_data = [d for d in horizon_days_elapsed if d in liq_days_present]
    # Eligibility arithmetic: uptime >= u over a horizon of L days needs downtime D <= (1-u) L.
    max_downtime_at_cap = (1 - uptime_min) * (h_cap - h_start)
    restart_deadline_for_cap = h_start + max_downtime_at_cap
    downtime_so_far = now - h_start  # zero bytes inside the horizon so far
    min_horizon_if_restarted_now = downtime_so_far / (1 - uptime_min)
    earliest_close_if_restarted_now = h_start + timedelta(
        days=max(horizon["minimum_duration_days"],
                 int(-(-min_horizon_if_restarted_now.total_seconds() // 86400))))

    record = {
        "schema_version": 1,
        "id": "collector-continuity-audit",
        "audit_instant_utc": iso(now),
        "purpose": "Outcome-blind continuity, coverage and eligibility audit of both append-only collectors against the pre-registered prospective horizon. No payload field read; no signal computed; no outcome formed.",
        "command": f"nice -n 10 python3 explorations/collector_continuity_audit.py --now {now.strftime('%Y-%m-%dT%H:%M:%SZ')} --out <this file>",
        "inputs": {
            "registration": str(REGISTRATION.relative_to(WORKSPACE)),
            "registration_sha256": reg_sha,
            "liquidation_dir": str(LIQ_DIR.relative_to(WORKSPACE)),
            "deribit_dir": str(DERIBIT_DIR.relative_to(WORKSPACE)),
            "horizon_start_utc": horizon["start_utc_inclusive"],
            "horizon_hard_cap_utc": horizon["hard_cap_utc"],
            "eligibility": elig,
        },
        "host_evidence": {
            "kern_boottime_utc": boot_time_utc(),
            "collector_processes_running_at_audit": "none: `pgrep -fl 'liquidation|deribit|collector'` empty and neither quant-liquidation-capture nor quant-deribit-gex-capture appears in hub ps",
        },
        "liquidation_collector": {
            "files": strip_private(liq),
            "day_seams": seams,
            "first_data_ts_utc": iso(liq_first),
            "last_data_ts_utc": iso(liq_last),
            "days_present": liq_days_present,
            "total_bytes": sum(f["bytes_pinned"] for f in liq),
            "total_data_records": sum(f["records_data"] for f in liq),
            "total_data_record_monotonicity_violations": sum(f["data_record_monotonicity_violations"] for f in liq),
            "data_gaps_over_limit_total": sum(len(f["data_gaps_over_limit"]) for f in liq),
            "any_trailing_partial_line": any(f["trailing_partial_line"] for f in liq),
        },
        "deribit_collector": {
            "days": strip_private(deribit_days),
            "last_received_utc_any_stream": iso(deribit_last),
            "ticker_bytes_by_day": ticker_census,
            "ticker_channel_ever_in_subscribe_ack": ticker_channel_ever_subscribed,
        },
        "horizon_coverage": {
            "horizon_start_utc": iso(h_start),
            "elapsed_in_horizon_seconds": int(horizon_elapsed.total_seconds()),
            "liquidation_bytes_at_or_after_start": liq_h_bytes,
            "liquidation_data_records_at_or_after_start": liq_h_records,
            "deribit_day_dirs_at_or_after_start": deribit_h_days,
            "horizon_utc_days_elapsed": horizon_days_elapsed,
            "horizon_utc_days_with_any_liquidation_data": horizon_days_with_data,
            "last_record_utc_both_collectors": iso(last_any),
            "missing_interval_utc": {"from": iso(min(liq_last, deribit_last)), "to": iso(now),
                                     "seconds": int(dead_for.total_seconds()),
                                     "hours": round(dead_for.total_seconds() / 3600, 2)},
            "uptime_fraction_inside_horizon_so_far": 0.0 if liq_h_records == 0 else None,
        },
        "eligibility_arithmetic": {
            "rule": "uptime >= u over a horizon of length L requires cumulative downtime D <= (1-u)*L; with the 90-day hard cap, D <= (1-u)*90 d",
            "max_cumulative_downtime_under_hard_cap_seconds": int(max_downtime_at_cap.total_seconds()),
            "restart_deadline_for_uptime_condition_to_remain_satisfiable_utc": iso(restart_deadline_for_cap),
            "downtime_inside_horizon_so_far_seconds": int(downtime_so_far.total_seconds()),
            "min_horizon_length_if_restarted_at_audit_instant_days": round(min_horizon_if_restarted_now.total_seconds() / 86400, 3),
            "earliest_utc_midnight_close_if_restarted_at_audit_instant_with_zero_further_downtime": iso(earliest_close_if_restarted_now),
            "void_horizon_days_under_15_minute_rule": horizon_days_elapsed,
            "denominator_caveat": "The registration does not pin the denominator of collector_uptime_fraction_min (horizon length vs. calendar since collection start). This arithmetic uses the horizon length, the only reading under which the fraction is a property of the horizon.",
        },
        "alpha_ledger": alpha_ledger(),
        "frozen_campaigns": verify_campaign_checksums(),
    }
    Path(args.out).write_text(json.dumps(record, indent=2) + "\n")
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print(json.dumps({"out": args.out, "cpu_seconds": round(ru.ru_utime + ru.ru_stime, 2)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
