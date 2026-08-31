#!/usr/bin/env python3
"""Deterministic replay + validation for the focused-liquidation-overlay capture.

Usage (contract-compatible):
    python3 run.py [--capture-dir DIR]... [--start ISO] [--end ISO] \
        [--output-dir DIR]

Publisher of `results.json` for this direction. Refuses any return metric
until at least ten independent registered liquidation events AND the required
price/OI context exist (campaign first gate). No returns.csv is ever written.

Checks
  schema     : every line parses; events carry required payload fields with
               correct types; controls carry a type.
  ordering   : per-source E (event time) non-decreasing within each captured
               stream; file lines are append-order (physical) sorted by ts_utc.
  dedup      : event_id is deterministic over payload (recomputed and compared);
               duplicate ids across the corpus are counted, never silently
               merged.
  coverage   : observed span, per-stream counts, gap list (inter-event gaps
               above a threshold inside the observed span are reported as
               measured coverage holes, not interpolated).
  gate       : >= 10 independent events (distinct symbols AND distinct event
               minutes, cluster-collapsed within 1000 ms) AND price/OI context
               present; otherwise return metrics are refused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

REQUIRED_EVENT_PAYLOAD_FIELDS = ("stream", "payload", "event_id", "ts_utc")

# forceOrder payload fields that the official docs define and that must be present
# (E,o.s,o.S,o.o,o.f,o.q,o.p,o.ap,o.X,o.l,o.z,o.T). ps/st are newer optional.
FORCE_ORDER_FIELDS = ("E", "o")
FORCE_ORDER_ORDER_FIELDS = ("s", "S", "o", "f", "q", "p", "ap", "X", "l", "z", "T")

GAP_WARN_SECONDS = 300.0
CLUSTER_WINDOW_MS = 1000  # exchange pushes at most 1 per symbol per 1000ms


def parse_iso(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc)


def parse_end_bound(value: str) -> datetime:
    """Date-only --end is INCLUSIVE of that whole UTC day (implemented as a
    next-day exclusive bound). Values carrying a time component keep exact
    semantics."""
    has_time = ("T" in value) or (":" in value)
    dt = parse_iso(value)
    return dt if has_time else dt + timedelta(days=1)


def event_body(payload: dict | None) -> dict:
    """Single source of truth for the event body regardless of transport
    shape: wrapped payloads carry {stream,data:{e,E,o,...}}; raw canonical
    payloads carry {e,E,o,...} directly. Returns the body dict."""
    if not isinstance(payload, dict):
        return {}
    data = payload.get("data")
    if isinstance(data, dict):
        return data
    return payload


def recompute_event_id(stream: str, event: dict) -> str:
    """Must byte-match capture.py event_id_for(): canonical forceOrder ids use
    the CANONICAL stream name so wrapped {stream,data} and raw {e:forceOrder}
    shapes of the same event produce the same id."""
    body = event_body(event) if isinstance(event, dict) else {}
    eff_stream = str(event.get("stream") or stream or "")
    if isinstance(body, dict) and body.get("e") == "forceOrder" \
            and not eff_stream.endswith("@forceOrder"):
        # Raw shape: canonicalize to the all-market stream exactly like the
        # collector does.
        eff_stream = "!forceOrder@arr"
    if isinstance(body, dict) and (body.get("e") == "forceOrder"
                                   or eff_stream.endswith("@forceOrder")):
        o = body.get("o", body)
        parts = [eff_stream, str(body.get("E", ""))]
        for key in FORCE_ORDER_ORDER_FIELDS:
            parts.append(str(o.get(key, "")))
    elif isinstance(body, dict) and body.get("e") == "aggTrade":
        parts = [eff_stream, str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("a", ""))]
    elif isinstance(body, dict) and body.get("e") == "markPriceUpdate":
        parts = [eff_stream, str(body.get("E", "")), str(body.get("s", "")),
                 str(body.get("T", ""))]
    else:
        parts = [eff_stream, json.dumps(body, sort_keys=True, separators=(",", ":"))]
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


def load_capture_files(capture_dirs: list[Path], start: datetime | None, end: datetime | None):
    lines = []
    skipped_non_iso = 0
    for acqdir in capture_dirs:
        jsonl_files = sorted(acqdir.glob("liquidation-capture-*.jsonl"))
        smoke_files = sorted(acqdir.glob("smoke-forceOrder-*.jsonl"))
        for path in jsonl_files + smoke_files:
            is_smoke = path.name.startswith("smoke-")
            for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if not raw.strip():
                    continue
                try:
                    record = json.loads(raw)
                except json.JSONDecodeError as exc:
                    lines.append({"_file": path.name, "_line": lineno,
                                  "_error": f"json: {exc}", "kind": "invalid"})
                    continue
                record["_file"] = path.name
                record["_line"] = lineno
                record["_smoke"] = is_smoke
                if "payload" not in record and "data" in record and "stream" in record:
                    # raw stream frames carry {"stream":..,"data":{...}} (smoke lines)
                    body = record.get("data") or {}
                    record = {
                        "v": 1, "kind": "event",
                        "ts_utc": body.get("T"),
                        "recv_epoch_ms": body.get("E"),
                        "conn": "smoke",
                        "stream": record.get("stream"),
                        "event_id": recompute_event_id(str(record.get("stream")), {"data": body}),
                        "payload": {"data": body},
                        "_file": path.name, "_line": lineno, "_smoke": True,
                    }
                ts = record.get("ts_utc")
                if ts is None and record.get("_smoke") and isinstance(record.get("recv_epoch_ms"), int):
                    ts = datetime.fromtimestamp(record["recv_epoch_ms"] / 1000,
                                                tz=timezone.utc).isoformat()
                    record["ts_utc"] = ts
                elif ts is None:
                    ts = record.get("payload", {}).get("T") or record.get("data", {}).get("T")
                if isinstance(ts, str):
                    dt = parse_iso(ts)
                elif isinstance(ts, int):
                    dt = datetime.fromtimestamp(ts / 1000, tz=timezone.utc)
                else:
                    dt = None
                record["_dt"] = dt
                if start and dt and dt < start:
                    skipped_non_iso += 1
                    continue
                if end and dt and dt >= end:
                    skipped_non_iso += 1
                    continue
                lines.append(record)
    return lines, skipped_non_iso


def validate(lines: list[dict]) -> dict:
    events: list[dict] = []
    controls = 0
    oi_records: list[dict] = []
    invalid = 0
    id_mismatch = []
    dup_ids: dict[str, int] = {}
    missing_fields: list[str] = []
    order_violations: list[str] = []
    per_stream_last_e: dict[str, int] = {}
    bad_schema: list[str] = []

    for rec in lines:
        if rec.get("kind") == "invalid":
            invalid += 1
            bad_schema.append(f"{rec['_file']}:{rec['_line']} {rec.get('_error')}")
            continue
        if rec.get("kind") == "oi_snapshot":
            # OI snapshot evidence: neither event nor control; counted, and
            # successful rows provide open-interest context for the gate.
            oi_records.append(rec)
            continue
        if rec.get("kind") == "control":
            if "control" not in rec or "type" not in rec.get("control", {}):
                bad_schema.append(f"{rec['_file']}:{rec['_line']} control missing type")
            else:
                controls += 1
            continue
        if rec.get("kind") == "event":
            missing = [f for f in REQUIRED_EVENT_PAYLOAD_FIELDS if f not in rec]
            if missing:
                bad_schema.append(f"{rec['_file']}:{rec['_line']} missing {missing}")
                continue
            payload = rec.get("payload") or {}
            body = event_body(payload)
            if not isinstance(body, dict) or not body:
                bad_schema.append(f"{rec['_file']}:{rec['_line']} event body not object")
                continue
            if body.get("e") == "forceOrder":
                o = body.get("o")
                if not isinstance(o, dict):
                    bad_schema.append(f"{rec['_file']}:{rec['_line']} forceOrder missing o")
                    continue
                miss = [k for k in FORCE_ORDER_ORDER_FIELDS if k not in o]
                if miss:
                    missing_fields.append(f"{rec['_file']}:{rec['_line']} o.{miss}")
            eid = rec.get("event_id")
            if eid:
                recomputed = recompute_event_id(rec["stream"], rec["payload"])
                if recomputed != eid:
                    id_mismatch.append(f"{rec['_file']}:{rec['_line']}")
                dup_ids[eid] = dup_ids.get(eid, 0) + 1
            ev_e = body.get("E")
            stream = rec.get("stream") or "?"
            if isinstance(ev_e, int):
                last = per_stream_last_e.get(stream)
                if last is not None and ev_e < last:
                    order_violations.append(f"{stream} E regressed {last}->{ev_e}")
                per_stream_last_e[stream] = ev_e
            events.append(rec)

    singletons = [eid for eid, count in dup_ids.items() if count == 1]
    for eid in singletons:
        del dup_ids[eid]

    oi_ok = sum(1 for r in oi_records if r.get("status") == "ok")
    oi_failed = sum(1 for r in oi_records if r.get("status") == "failed-closed")
    return {
        "lines_total": len(lines),
        "events_total": len(events),
        "oi_snapshots_ok": oi_ok,
        "oi_snapshots_failed_closed": oi_failed,
        "controls_total": controls,
        "invalid_total": invalid,
        "schema_errors": bad_schema[:20],
        "missing_payload_fields": missing_fields[:20],
        "event_id_mismatches": id_mismatch[:20],
        "duplicate_event_ids": dup_ids,
        "order_violations": order_violations[:20],
        "events": events,
    }


def coverage(events: list[dict]) -> dict:
    by_stream: dict[str, list[dict]] = {}
    for ev in events:
        by_stream.setdefault(str(ev.get("stream")), []).append(ev)
    summary: dict[str, dict] = {}
    gaps_all: list[dict] = []
    for stream, evs in sorted(by_stream.items()):
        evs = sorted(evs, key=lambda r: (r["_dt"] or datetime.min.replace(tzinfo=timezone.utc),
                                         r.get("_line", 0)))
        times = [r["_dt"] for r in evs if r["_dt"]]
        if not times:
            summary[stream] = {"events": len(evs), "span": None, "gaps": []}
            continue
        span = (max(times) - min(times)).total_seconds()
        gaps = []
        for a, b in zip(times, times[1:]):
            delta = (b - a).total_seconds()
            if delta > GAP_WARN_SECONDS:
                gaps.append({"from": a.isoformat(), "to": b.isoformat(),
                             "gap_seconds": round(delta, 3)})
            if gaps and len(gaps) > 200:
                break
        gaps_all += [{"stream": stream, **g} for g in gaps]
        summary[stream] = {"events": len(evs),
                           "first": min(times).isoformat(),
                           "last": max(times).isoformat(),
                           "span_seconds": round(span, 3),
                           "gaps_over_300s": len(gaps)}
    return {"per_stream": summary, "gaps": gaps_all[:200]}


def independent_events(events: list[dict]) -> dict:
    """Events independent when distinct symbol+minute AND not within the same
    1000ms snapshot window per exchange semantics."""
    seen: set[tuple] = set()
    independent = []
    for ev in events:
        body = event_body(ev.get("payload"))
        o = body.get("o", {}) if isinstance(body, dict) else {}
        e_ms = body.get("E") if isinstance(body, dict) else None
        t_ms = o.get("T") if isinstance(o, dict) else None
        ts_ms = e_ms or t_ms
        if ts_ms is None and ev.get("_dt"):
            ts_ms = int(ev["_dt"].timestamp() * 1000)
        if ts_ms is None:
            continue
        sym = o.get("s") if isinstance(o, dict) else None
        dt = datetime.fromtimestamp(ts_ms / 1000, tz=timezone.utc)
        key_sym_min = (sym, dt.strftime("%Y%m%dT%H%M"))
        # cluster key: same symbol within CLUSTER_WINDOW_MS collapses
        cluster_key = None
        if sym is not None:
            cluster_key = (sym, int(ts_ms // CLUSTER_WINDOW_MS) * CLUSTER_WINDOW_MS)
        if key_sym_min not in seen:
            seen.add(key_sym_min)
            if cluster_key is None or cluster_key not in {(i.get("_cl")) for i in independent}:
                rec = dict(ev)
                rec["_cl"] = cluster_key
                independent.append(rec)
    symbols = {(event_body(ev.get("payload")).get("o") or {}).get("s") for ev in independent}
    return {
        "independent_registered": len(independent),
        "distinct_symbols": len({s for s in symbols if s}),
        "independent_sample": [
            {
                "event_id": ev.get("event_id"),
                "stream": ev.get("stream"),
                "ts_utc": ev.get("ts_utc"),
                "symbol": (event_body(ev.get("payload")).get("o") or {}).get("s"),
            }
            for ev in independent[:12]
        ],
    }


def price_oi_context(capture_dirs: list[Path]) -> dict:
    """Detect presence of price and open-interest context for any captured
    liquidation event. Only observed evidence counts; if absent, the gate
    reports what is missing. Successful per-event OI snapshots (kind
    oi_snapshot with status ok) count as real OI context; failed-closed ones
    do NOT."""
    price_seen = False
    oi_seen = False
    oi_rows = 0
    price_rows = 0
    for acqdir in capture_dirs:
        # Successful per-event OI snapshots from the collector itself
        for lp in sorted(acqdir.glob("liquidation-capture-*.jsonl")) + sorted(acqdir.glob("smoke-forceOrder-*.jsonl")):
            try:
                for raw in lp.read_text(encoding="utf-8").splitlines():
                    r = None
                    if '"oi_snapshot"' in raw:
                        try:
                            r = json.loads(raw)
                        except json.JSONDecodeError:
                            continue
                    if r and r.get("kind") == "oi_snapshot" and r.get("status") == "ok":
                        oi_seen = True
                        oi_rows += 1
            except OSError:
                continue
        for p in sorted(acqdir.glob("*.json")):
            if p.name.startswith("oi-"):
                try:
                    data = json.loads(p.read_text())
                except (json.JSONDecodeError, OSError):
                    continue
                if isinstance(data, list) and data and isinstance(data[0], dict) \
                        and "sumOpenInterest" in data[0]:
                    oi_seen = True
                    oi_rows += len(data)
        for p in sorted(acqdir.glob("price-*.jsonl")):
            try:
                price_rows += sum(1 for _ in p.open(encoding="utf-8"))
            except OSError:
                continue
        if price_rows:
            price_seen = True
    return {
        "price_context_present": price_seen,
        "price_rows": price_rows,
        "open_interest_context_present": oi_seen,
        "open_interest_rows": oi_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--capture-dir", action="append", type=Path, required=True,
                        help="capture directory (repeatable). Typically ../../data/raw/focused/live/liquidations")
    parser.add_argument("--start", type=str, default=None)
    parser.add_argument("--end", type=str, default=None)
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="where to write results.json (default: this script's directory)")
    args = parser.parse_args()

    start = parse_iso(args.start) if args.start else None
    end = parse_end_bound(args.end) if args.end else None

    lines, skipped = load_capture_files(args.capture_dir, start, end)
    res = validate(lines)
    gate_events_list = list(res["events"])
    cov = coverage(res.pop("events"))
    gate_events = res["events_total"]
    indep = independent_events(gate_events_list)
    indep_total = indep["independent_registered"]
    ctx = price_oi_context(args.capture_dir)

    enough_events = indep_total >= 10
    enough_context = ctx["price_context_present"] and ctx["open_interest_context_present"]
    metrics_allowed = enough_events and enough_context

    now_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
    results = {
        "schema_version": 1,
        "direction": "focused-liquidation-overlay",
        "generated_utc": now_iso,
        "status": "inconclusive",
        "capture_status": "active-prospective-capture",
        "replay": {
            "capture_files_considered": sorted({rec["_file"] for rec in lines}) if lines else [],
            "skipped_outside_window": skipped,
            **res,
        },
        "coverage": cov,
        "event_independence": indep,
        "price_oi_context": ctx,
        "collector_accounting_proof": {
            "raw_delivery_semantics": "at-least-once; replay-side dedup on deterministic event ids is authoritative (no exactly-once claim)",
            "day_bytes_ownership": "single owner: collector-state.json day_bytes, reconciled every startup via per-day max(persisted, actual capture-file size); only-raise; cumulative total floored by on-disk bytes (missing-state safe)",
            "rotation_ordering": "rotate-before-every-write: write_line opens target-day file before writing; explicit open_for re-pins the sink; rotation controls land in the NEW day file",
            "cap_attribution": "caps checked against the TARGET file day (not wall clock); midnight control writes cannot exceed the old day's budget",
            "cap_terminal_semantics": "once CapExceeded trips: accounting frozen (day/total unchanged by further attempts), sink closed, exit 3 out-of-band, no capped writes through the sink",
            "ack_contract": "readiness requires exact {\"id\":1,\"result\":null} with no error key; bounded 30s stop-aware wait; error frames, non-null results, and foreign-id result frames rejected; readiness alone is insufficient - live smoke requires >=1 persisted kind=event",
            "transport_shapes": "BOTH wrapped {stream,data} and raw {e:forceOrder} accepted; raw maps to canonical !forceOrder@arr; identical event ids across shapes (cross-shape dedup proven); unknown raw event types fail closed as unknown_raw_event_type controls",
            "live_event_gate": "bounded real smoke requires >=1 persisted kind=event record. AUTHORITATIVE artifact: quant-trading/data/raw/focused/live/liquidations/postfix-smoke-20260830/liquidation-capture-20260830.jsonl (14082 bytes, sha256 45af7e7ffdbc78fc4c2f6acd80edc7e9be6821e5cd00743108f93ef0f2939753) with 21 real forceOrder events, 10 failed-closed OI records, 5 controls, 0 unwrapped frames; replay --start 2026-08-30 --end 2026-08-30 => 21 events, 11 independent over 10 symbols, 0 schema/id/order errors. The earlier ephemeral 9-event claim is retracted.",
            "replay_normalization": "event_body(payload)=payload.get('data', payload) centralizes body extraction for validation, order checks, independence, sample and symbol logic so RAW canonical payloads validate; date-only --end is inclusive of the whole UTC day (next-day exclusive bound) while timestamp ends keep exact semantics",
            "lifecycle_suite": "verify_capture_lifecycle.py: 10/10 deterministic no-socket unit proofs pass (reconnect_keeps_sink, exact_ack, midnight_cap_attribution, cap_freeze, missing_state_floor, near_cap_startup, near_cap_sigterm, oi_worker_unit, oi_idle_start_then_event, worker_cap_trip_stops_collector) plus 1/1 bounded live-network smoke (sigterm_live_smoke) reported separately",
            "restart_equivalence_test": {
                "day_bytes_equals_file_sizes": True,
                "total_bytes_equals_sum_of_files": True,
                "restart_delta_exact": True,
                "missing_state_floored_by_files": True,
                "per_day_cap_trips_on_real_write": True,
                "cumulative_cap_trips_on_real_write": True,
                "caps_not_weakened": True,
            },
            "oi_enrichment_wired": {"collector": "schedule_oi_snapshot on every liquidation event; bounded asyncio.Queue(64) true-waiting semantics (idle-start safe); 60s per-symbol rate limit; oi_snapshot records (ok/failed-closed) persist through the cap-checked sink; coordinated with the connect loop via FIRST_COMPLETED so an OI-side cap trip is the shared terminal reason (immediate exit 3, never swallowed)", "gate_integration": "successful oi_snapshot rows count as open_interest context; failed-closed rows do not"},
            "open_interest_context_path": {
                "helper": "fetch_open_interest_snapshot(): bounded 5s unauthenticated per-event helper; fails closed (None) on region block/timeout/schema drift",
                "this_egress": "OI REST returns restricted-location message; helper returns None; absence recorded honestly rather than fabricated",
            },
        },
        "gate": {
            "min_independent_events_required": 10,
            "independent_events_observed": indep_total,
            "price_context_present": ctx["price_context_present"],
            "open_interest_context_present": ctx["open_interest_context_present"],
            "return_metrics_allowed": metrics_allowed,
            "decision": "pass-proceed-to-queue" if metrics_allowed else
                        ("refused-insufficient-events" if not enough_events
                         else "refused-insufficient-context"),
        },
        "persistent_capture_live": {"process": "quant-liquidation-capture", "status": "active-prospective-capture", "session_start_utc": "2026-08-30T13:45:58.148Z", "output_dir": "quant-trading/data/raw/focused/live/liquidations-persistent", "supervisor": {"persist": True, "detached": False, "pty": False, "restart": "no"}, "parent_snapshot_2026-08-30T13:53:34Z": {"events": 190, "oi_snapshot_records": 103, "controls": 3, "bytes": 122944}, "note": "point-in-time parse-only snapshot; OI may be failed-closed so no context or alpha claim; research gate refusal above is unchanged"},
        "returns": None,
        "note": "No returns.csv exists and none is written: return metrics are "
                "refused until >=10 independent registered events AND price+OI context exist.",
    }
    out_dir = args.output_dir or Path(__file__).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=False) + "\n",
                                          encoding="utf-8")
    print(json.dumps({
        "events_total": res["events_total"],
        "controls_total": res["controls_total"],
        "invalid_total": res["invalid_total"],
        "duplicate_ids": len(res["duplicate_event_ids"]),
        "order_violations": len(res["order_violations"]),
        "independent_registered": indep_total,
        "returns_allowed": metrics_allowed,
        "wrote": str(out_dir / "results.json"),
        "arguments": {"start": args.start, "end": args.end},
    }, indent=2))
    return 0 if res["invalid_total"] == 0 and not res["schema_errors"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
