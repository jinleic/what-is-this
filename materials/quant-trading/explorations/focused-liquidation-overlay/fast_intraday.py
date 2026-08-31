#!/usr/bin/env python3
"""Frozen fast screen: liquidation_btc_eth_post_onset.

Implements, exactly, the screen registered in
quant-trading/explorations/accelerated-intraday-contract.json (frozen before
outcome analysis; status frozen-amended-before-outcome-analysis), including
semantic_amendment (2026-08-30T17:10:49Z) and final_semantic_amendment
(2026-08-30T17:14:49Z, schema clarified 2026-08-30T17:15:56Z):

  - Inputs: the persistent Binance !forceOrder@arr capture and the concurrent
    public Deribit btc_usd / eth_usd index records. Research only: no orders,
    no account access, no credentials, no deployment claim.
  - Snapshot: for every live JSONL input the byte size is statted FIRST, then
    exactly that byte prefix is hashed (SHA-256) and parsed; one trailing
    partial line is ignored and its byte count recorded.
  - Clocks: Binance forceOrder body.E and Deribit index data.timestamp are the
    only exchange clocks used in math; receipt timestamps are provenance only.
  - Dedup: existing deterministic event ids (collector-verified formulas).
  - Onsets: per BTCUSDT/ETHUSDT symbol, the first deduplicated forceOrder
    event after >= 300000 ms without a prior deduplicated same-symbol event
    inside the same continuously captured segment. Segments begin ONLY at
    top-level kind=control control.type=connect records (top-level ts_utc);
    disconnect / session_stop / a later connect end the segment; session_start
    resets but does not begin. First same-symbol event in a segment requires
    body.E - connect_capture_epoch_ms >= 300000.
  - Alignment: initial proxy P_0 = latest same-underlying Deribit index
    observation with timestamp <= body.E and no more than 5000 ms old. Horizon
    P_h = first observation with timestamp >= body.E + h and <= body.E + h +
    5000 ms. No interpolation.
  - Outcome: net_h = direction * (P_h / P_0 - 1) - 0.0002, direction +1 for
    o.S == BUY and -1 for o.S == SELL; unrounded throughout.
  - Cohort: single gated cohort = eligible onsets with valid P_0, P_5m and
    P_15m alignments by the candidate stop. All sample counts and every
    5m/15m gate use this same cohort; the 1m outcome is descriptive only.
  - Earliest stop: the earliest +15m outcome-availability time (all
    sample conditions true first) among observations at or before
    final_semantic_amendment.deadline_utc; otherwise the full common prefix
    and inconclusive.
  - Truth: sample met + every gate observed true => fast-screen-pass;
    sample met + any observed false gate => falsified; any sample condition
    unmet/unknown or required alignment unavailable => inconclusive. Undefined
    values serialize as null ("UNKNOWN"); never a non-standard JSON Infinity.

This program never mutates raw capture, never talks to the network, and emits
one deterministic report. Deterministic reruns against the same frozen
prefixes produce byte-identical math blocks.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = REPO_ROOT / "explorations" / "accelerated-intraday-contract.json"
EVENTS_DIR = REPO_ROOT / "data" / "raw" / "focused" / "live" / "liquidations-persistent"
INDEX_ROOT = REPO_ROOT / "data" / "raw" / "focused" / "live" / "deribit"
DEFAULT_OUT = Path(__file__).resolve().parent / "fast_intraday_results.json"

SCREEN_ID = "liquidation_btc_eth_post_onset"
TARGET_SYMBOLS = ("BTCUSDT", "ETHUSDT")
SYMBOL_UNDERLYING = {"BTCUSDT": "btc_usd", "ETHUSDT": "eth_usd"}
ONSET_GAP_MS = 300_000
MAX_INDEX_AGE_MS = 5_000
HORIZON_MS = {"1m": 60_000, "5m": 300_000, "15m": 900_000}
BLOCK_MS = 900_000
HAIRCUT = 0.0002
FORCE_ORDER_ORDER_FIELDS = ("s", "S", "o", "f", "q", "p", "ap", "X", "l", "z", "T")
CANON_STREAM = "!forceOrder@arr"

EXIT_OK = 0            # deterministic report written (any truth status)
EXIT_CONTRACT = 2      # contract missing/unusable
EXIT_INPUTS = 3        # input absence or prefix integrity failure
EXIT_SELFTEST = 5      # selftest failed


# --------------------------------------------------------------------------
# time + money helpers
# --------------------------------------------------------------------------

def parse_iso_ms(text: str) -> int | None:
    if not isinstance(text, str) or not text:
        return None
    try:
        dt = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return int(round(dt.timestamp() * 1000))


def iso_utc(ms: int | None) -> str | None:
    if ms is None:
        return None
    dt = datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc)
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{dt.microsecond // 1000:03d}Z"


def median(values: list[float]) -> float | None:
    """Average of the two middle values for even n (deterministic, unrounded)."""
    if not values:
        return None
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def finite(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and x == x and x not in (
        float("inf"), float("-inf"))


# --------------------------------------------------------------------------
# snapshot: stat-first byte-prefix freeze of live inputs
# --------------------------------------------------------------------------

def snapshot_file(path: Path) -> dict:
    """Stat the byte size first, then hash and parse exactly that prefix."""
    prefix_bytes = path.stat().st_size          # authoritative size, taken FIRST
    with path.open("rb") as fh:
        raw = fh.read(prefix_bytes)             # exactly prefix_bytes, no more
    if len(raw) != prefix_bytes:                # file shrank mid-read (defensive)
        raise RuntimeError(f"{path}: size changed during snapshot read")
    try:
        display_path = str(path.relative_to(REPO_ROOT))
    except ValueError:
        display_path = str(path)
    report = {
        "path": display_path,
        "prefix_bytes": prefix_bytes,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "trailing_partial_line_bytes": 0,
        "lines_total": 0,
        "json_objects": 0,
        "malformed_lines": 0,
        "first_exchange_ts_ms": None,
        "latest_exchange_ts_ms": None,
    }
    lines = raw.split(b"\n")
    if lines and lines[-1] != b"":
        report["trailing_partial_line_bytes"] = len(lines[-1])
        lines = lines[:-1]                      # ignore ONE trailing partial line
    elif lines:
        lines = lines[:-1]                      # drop empty element after final \n
    parsed = []
    for line in lines:
        report["lines_total"] += 1
        s = line.strip()
        if not s:
            continue
        try:
            rec = json.loads(s.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            report["malformed_lines"] += 1
            continue
        report["json_objects"] += 1
        parsed.append(rec)
    return report, parsed


def canonical_event_id(payload: dict) -> str | None:
    """Byte-matches capture.py event_id_for() for forceOrder payloads."""
    if not isinstance(payload, dict):
        return None
    body = payload.get("data") if isinstance(payload.get("data"), dict) else payload
    if not isinstance(body, dict) or body.get("e") != "forceOrder":
        return None
    o = body.get("o", body)
    if not isinstance(o, dict):
        o = {}
    parts = [CANON_STREAM, str(body.get("E", ""))]
    for key in FORCE_ORDER_ORDER_FIELDS:
        parts.append(str(o.get(key, "")))
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


def index_dedup_id(rec: dict) -> str | None:
    """Byte-matches focused-deribit-gex capture.py: sha256('index|name|ts')."""
    data = rec.get("data")
    if not isinstance(data, dict):
        return None
    name = data.get("index_name") or rec.get("index_name")
    ts = data.get("timestamp")
    if name is None or ts is None:
        return None
    return hashlib.sha256(f"index|{name}|{ts}".encode("utf-8")).hexdigest()


# --------------------------------------------------------------------------
# contract loading (thresholds are read from the frozen file, never hardcoded)
# --------------------------------------------------------------------------

def load_contract(path: Path) -> dict:
    raw = path.read_bytes()
    contract = json.loads(raw.decode("utf-8"))
    screen = contract["screens"][SCREEN_ID]
    for key in ("event_definition", "price_alignment", "outcome", "sample_stop",
                "gate", "truth_table"):
        if key not in screen:
            raise KeyError(f"contract screen missing {key}")
    fsa = contract["final_semantic_amendment"]
    for key in ("deadline_utc", "deadline_rule"):
        if not fsa.get(key):
            raise KeyError(f"final_semantic_amendment missing {key}")
    for key in ("return_formula", "single_cohort", "session_boundary_clock_exception"):
        if key not in fsa.get("liquidation", {}):
            raise KeyError(f"final_semantic_amendment.liquidation missing {key}")
    return {
        "raw": contract,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "screen": screen,
        "sample_stop": screen["sample_stop"],
        "gate": screen["gate"],
        "deadline_ms": parse_iso_ms(fsa["deadline_utc"]),
    }


# --------------------------------------------------------------------------
# liquidation corpus -> segments -> onsets
# --------------------------------------------------------------------------

class Segment:
    __slots__ = ("connect_ms", "connect_iso", "start_line", "events")

    def __init__(self, connect_ms: int, connect_iso: str | None, start_line: int):
        self.connect_ms = connect_ms
        self.connect_iso = connect_iso
        self.start_line = start_line
        self.events: list[dict] = []            # capture-order forceOrder events


def build_onsets(event_records: list[dict], cutoff_ms: int) -> dict:
    """Dedup, segment walk, and onset qualification exactly per amendment."""
    stats = {
        "parsed_events": 0,
        "dedup_duplicate_ids": 0,
        "id_mismatch_records": 0,
        "non_forceOrder_records": 0,
        "missing_exchange_E": 0,
        "non_target_symbol_events": 0,
        "after_cutoff_events": 0,
        "orphan_events_outside_segment": 0,
        "control_connect": 0,
        "control_session_start": 0,
        "control_disconnect": 0,
        "control_session_stop": 0,
        "control_other": 0,
        "segments": [],
        "segment_events_total": 0,
        "eligible_onsets": [],
        "onsets_by_symbol": {s: 0 for s in TARGET_SYMBOLS},
        "onsets_by_side": {"BUY": 0, "SELL": 0, "other": 0},
        "incomplete_alignments": {"initial": 0, "horizon_1m": 0, "horizon_5m": 0,
                                  "horizon_15m": 0},
        "incomplete_onsets": [],
        "side_missing_or_unknown": 0,
    }

    # ---- dedup by deterministic event_id (first occurrence wins, stable order)
    seen_ids: set[str] = set()
    events: list[dict] = []
    for rec in event_records:
        kind = rec.get("kind")
        if kind not in ("event", "control"):
            continue                    # oi_snapshot etc. handled where relevant
        stats["parsed_events"] += 1
        if kind == "control":
            events.append({"rec": rec, "event_id": None})   # flow through to segmentation
            continue
        eid = rec.get("event_id")
        recomputed = canonical_event_id(rec.get("payload"))
        if recomputed is not None and eid is not None and eid != recomputed:
            stats["id_mismatch_records"] += 1
        if recomputed is not None and eid is None:
            eid = recomputed
        if eid is None:
            stats["missing_exchange_E"] += 1
            continue
        if eid in seen_ids:
            stats["dedup_duplicate_ids"] += 1
            continue
        seen_ids.add(eid)
        events.append({"rec": rec, "event_id": eid})

    # ---- capture segments (connect begins; disconnect/session_stop/later
    # connect end; session_start resets without beginning)
    segments: list[Segment] = []
    current: Segment | None = None
    for ev in events:
        rec = ev["rec"]
        kind = rec.get("kind")
        if kind == "control":
            if ev.get("event_id") is not None:
                continue                # hybrid shape; treated as event elsewhere
            ctype = (rec.get("control") or {}).get("type")
            if ctype == "connect":
                stats["control_connect"] += 1
                connect_ms = parse_iso_ms(rec.get("ts_utc") or "")
                current = Segment(connect_ms, rec.get("ts_utc"), None)
                segments.append(current)
            elif ctype == "session_start":
                stats["control_session_start"] += 1
                current = None                      # resets, does not begin
            elif ctype == "disconnect":
                stats["control_disconnect"] += 1
                current = None
            elif ctype == "session_stop":
                stats["control_session_stop"] += 1
                current = None
            else:
                stats["control_other"] += 1
            continue
        if kind != "event" or ev.get("event_id") is None:
            continue                                # oi_snapshot and others ignored
        body = rec.get("payload")
        body = body.get("data") if isinstance(body, dict) and isinstance(
            body.get("data"), dict) else body
        o = body.get("o") if isinstance(body, dict) else None
        if not isinstance(body, dict) or body.get("e") != "forceOrder":
            stats["non_forceOrder_records"] += 1
            continue
        if not isinstance(o, dict):
            o = {}
        E = body.get("E")
        if not isinstance(E, int) or isinstance(E, bool):
            stats["missing_exchange_E"] += 1
            continue
        symbol = o.get("s")
        if symbol not in TARGET_SYMBOLS:
            stats["non_target_symbol_events"] += 1
            continue
        if E > cutoff_ms:
            stats["after_cutoff_events"] += 1
            continue
        if current is None:
            stats["orphan_events_outside_segment"] += 1
            continue
        current.events.append({
            "event_id": ev["event_id"],
            "symbol": symbol,
            "side": o.get("S"),
            "E": int(E),
            "seq": len(segments) * 10_000_000 + len(current.events),
        })
        stats["segment_events_total"] += 1

    # ---- onset qualification per symbol inside each segment
    def add_onset(seg: Segment, ev: dict, prior_e: int | None, reason: str) -> None:
        side = ev["side"]
        if side not in ("BUY", "SELL"):
            stats["side_missing_or_unknown"] += 1
            return
        onset = {
            "event_id": ev["event_id"],
            "symbol": ev["symbol"],
            "side": side,
            "onset_ms": ev["E"],
            "segment_connect_ms": seg.connect_ms,
            "eligibility": reason,
            "prior_same_symbol_E_ms": prior_e,
        }
        stats["eligible_onsets"].append(onset)
        stats["onsets_by_symbol"][ev["symbol"]] += 1
        stats["onsets_by_side"][side] += 1

    for seg in segments:
        by_symbol: dict[str, list[dict]] = {}
        for ev in seg.events:
            by_symbol.setdefault(ev["symbol"], []).append(ev)
        for symbol, evs in by_symbol.items():
            evs.sort(key=lambda e: (e["E"], e["seq"]))      # exchange-clock order
            prev_e: int | None = None
            for ev in evs:
                if prev_e is None:
                    if seg.connect_ms is not None and ev["E"] - seg.connect_ms >= ONSET_GAP_MS:
                        add_onset(seg, ev, None, f"first_in_segment_silence_from_connect_{ONSET_GAP_MS}ms")
                else:
                    if ev["E"] - prev_e >= ONSET_GAP_MS:
                        add_onset(seg, ev, prev_e, "gap_since_prior_same_symbol")
                prev_e = ev["E"]

    stats["eligible_onsets"].sort(key=lambda x: (x["onset_ms"], x["event_id"]))
    stats["segments"] = [{
        "connect_ts_utc": seg.connect_iso,
        "connect_capture_epoch_ms": seg.connect_ms,
        "forceOrder_events": len(seg.events),
    } for seg in segments]
    return stats


# --------------------------------------------------------------------------
# Deribit index series
# --------------------------------------------------------------------------

def load_index_series(index_records: list[dict], cutoff_ms: int) -> dict:
    """Per-underlying deduplicated series sorted by exchange timestamp.
    Duplicate timestamps keep the LAST occurrence in stable input order."""
    stats = {
        "records": 0,
        "dedup_duplicate_ids": 0,
        "dedup_id_missing": 0,
        "dedup_id_mismatch": 0,
        "after_cutoff_observations": 0,
        "nonpositive_or_bad_price": 0,
        "per_underlying": {},
        "latest_observation_ms": None,
        "first_observation_ms": None,
    }
    seen_ids: set[str] = set()
    per_key: dict[str, list[tuple[int, int, float]]] = {}   # key -> (ts, seq, price)
    for seq, rec in enumerate(index_records):
        stats["records"] += 1
        data = rec.get("data")
        if not isinstance(data, dict):
            continue
        name = data.get("index_name") or rec.get("index_name") or ""
        if "btc_usd" in name:
            key = "btc_usd"
        elif "eth_usd" in name:
            key = "eth_usd"
        else:
            continue
        ts = data.get("timestamp")
        price = data.get("price")
        if not isinstance(ts, int) or isinstance(ts, bool) or not finite(price) or price <= 0:
            stats["nonpositive_or_bad_price"] += 1
            continue
        did = rec.get("dedup_id")
        expected = index_dedup_id(rec)
        if did is None:
            stats["dedup_id_missing"] += 1
            did = expected
        if expected is not None and did != expected:
            stats["dedup_id_mismatch"] += 1
        if did in seen_ids:
            stats["dedup_duplicate_ids"] += 1
            continue
        seen_ids.add(did)
        if ts > cutoff_ms:
            stats["after_cutoff_observations"] += 1
            continue
        per_key.setdefault(key, []).append((int(ts), seq, float(price)))
        if stats["first_observation_ms"] is None or ts < stats["first_observation_ms"]:
            stats["first_observation_ms"] = ts
        if stats["latest_observation_ms"] is None or ts > stats["latest_observation_ms"]:
            stats["latest_observation_ms"] = ts

    series: dict[str, dict] = {}
    for key, rows in per_key.items():
        rows.sort(key=lambda r: (r[0], r[1]))
        # duplicate timestamps keep the LAST stable occurrence
        deduped: list[tuple[int, float]] = []
        for ts, _seq, price in rows:
            if deduped and deduped[-1][0] == ts:
                deduped[-1] = (ts, price)
            else:
                deduped.append((ts, price))
        times = [ts for ts, _p in deduped]
        series[key] = {
            "times_ms": times,
            "prices": [p for _ts, p in deduped],
            "observations": len(deduped),
            "duplicate_timestamp_collapses": len(rows) - len(deduped),
        }
        stats["per_underlying"][key] = series[key]["observations"]
    return {"series": series, "stats": stats}


def align_initial(times: list[int], prices: list[float], target_ms: int) -> tuple[int, float] | None:
    """Latest observation <= target and no more than 5000 ms old."""
    i = bisect.bisect_right(times, target_ms) - 1
    if i < 0:
        return None
    ts = times[i]
    if target_ms - ts > MAX_INDEX_AGE_MS:
        return None
    return ts, prices[i]


def align_horizon(times: list[int], prices: list[float], target_ms: int) -> tuple[int, float] | None:
    """First observation >= target and at most 5000 ms late."""
    i = bisect.bisect_left(times, target_ms)
    if i >= len(times):
        return None
    ts = times[i]
    if ts - target_ms > MAX_INDEX_AGE_MS:
        return None
    return ts, prices[i]


def compute_onset_outcomes(onsets: list[dict], index: dict, cutoff_ms: int) -> list[dict]:
    """Attach P_0 / P_1m / P_5m / P_15m alignments and net returns. Unrounded."""
    out = []
    for onset in onsets:
        symbol = onset["symbol"]
        E = onset["onset_ms"]
        key = SYMBOL_UNDERLYING[symbol]
        s = index["series"].get(key, {"times_ms": [], "prices": []})
        direction = 1 if onset["side"] == "BUY" else -1
        mark = {"descriptive_only_horizons": ["1m"]}

        p0 = align_initial(s["times_ms"], s["prices"], E)

        def outcome(hname: str) -> dict:
            tgt = E + HORIZON_MS[hname]
            if tgt > cutoff_ms:                     # horizon beyond common cutoff
                return {"target_ms": tgt, "observed": None, "missing": "beyond_common_cutoff"}
            ph = align_horizon(s["times_ms"], s["prices"], tgt)
            if ph is None:
                return {"target_ms": tgt, "observed": None, "missing": "no_index_observation_in_5s_window"}
            o = {"target_ms": tgt, "observed_ts_ms": ph[0], "price": ph[1]}
            if p0 is not None:
                o["signed_return"] = direction * (ph[1] / p0[1] - 1.0)
                o["net_return"] = o["signed_return"] - HAIRCUT
            return o

        for hname in ("1m", "5m", "15m"):
            mark[f"P_{hname}"] = outcome(hname)
        if p0 is None:
            mark["P_0"] = {"observed": None, "missing": "no_index_observation_within_5000ms_before"}
            for hname in ("1m", "5m", "15m"):
                mark[f"P_{hname}"].pop("signed_return", None)
                mark[f"P_{hname}"].pop("net_return", None)
        else:
            mark["P_0"] = {"observed_ts_ms": p0[0], "price": p0[1]}

        complete = (p0 is not None
                    and mark["P_5m"].get("price") is not None
                    and mark["P_15m"].get("price") is not None)
        mark["cohort_complete"] = bool(complete)
        if complete:
            mark["availability_ms"] = mark["P_15m"]["observed_ts_ms"]
        row = dict(onset)
        row.update(mark)
        out.append(row)
    return out


# --------------------------------------------------------------------------
# earliest stop + sample conditions + gates on the single cohort
# --------------------------------------------------------------------------

def sample_conditions(cohort: list[dict], stop_cfg: dict) -> dict:
    by_symbol = {s: 0 for s in TARGET_SYMBOLS}
    by_side = {"BUY": 0, "SELL": 0}
    blocks = set()
    for o in cohort:
        by_symbol[o["symbol"]] += 1
        if o["side"] in by_side:
            by_side[o["side"]] += 1
        blocks.add(o["onset_ms"] // BLOCK_MS)
    n = len(cohort)
    return {
        "by_symbol": by_symbol,
        "by_side": by_side,
        "nonempty_15m_blocks": len(blocks),
        "conditions": [
            {"name": "minimum_complete_onsets_at_15m", "required": stop_cfg["minimum_complete_onsets_at_15m"],
             "observed": n, "met": n >= stop_cfg["minimum_complete_onsets_at_15m"]},
            {"name": "minimum_per_asset", "required": stop_cfg["minimum_per_asset"],
             "observed": {"BTCUSDT": by_symbol["BTCUSDT"], "ETHUSDT": by_symbol["ETHUSDT"]},
             "met": by_symbol["BTCUSDT"] >= stop_cfg["minimum_per_asset"]
                    and by_symbol["ETHUSDT"] >= stop_cfg["minimum_per_asset"]},
            {"name": "minimum_per_side", "required": stop_cfg["minimum_per_side"],
             "observed": by_side, "met": by_side["BUY"] >= stop_cfg["minimum_per_side"]
                    and by_side["SELL"] >= stop_cfg["minimum_per_side"]},
            {"name": "minimum_nonempty_15m_onset_blocks", "required": stop_cfg["minimum_nonempty_15m_onset_blocks"],
             "observed": len(blocks),
             "met": len(blocks) >= stop_cfg["minimum_nonempty_15m_onset_blocks"]},
        ],
    }


def find_earliest_stop(onsets: list[dict], stop_cfg: dict, deadline_ms: int) -> dict:
    """Cohort availability time = timestamp of the onset's 15m outcome
    observation; search ascending candidate stops meeting every sample
    condition, restricted to availability at or before the deadline."""
    complete = sorted((o for o in onsets if o["cohort_complete"]),
                      key=lambda o: (o["availability_ms"], o["event_id"]))
    tried = []
    for k in range(1, len(complete) + 1):
        T = complete[k - 1]["availability_ms"]
        if T > deadline_ms:
            break
        if k > 1 and complete[k - 2]["availability_ms"] == T:
            continue                                # same stop instant already tried
        cohort = [o for o in complete if o["availability_ms"] <= T]
        cond = sample_conditions(cohort, stop_cfg)
        all_met = all(c["met"] for c in cond["conditions"])
        tried.append({"candidate_ms": T, "cohort_size": len(cohort), "all_met": all_met})
        if all_met:
            return {"found": True, "stop_ms": T, "cohort": cohort, "tried": tried,
                    "conditions": cond}
    cohort = complete                                # full common prefix fallback
    cond = sample_conditions(cohort, stop_cfg)
    return {"found": False, "stop_ms": None, "cohort": cohort, "tried": tried,
            "conditions": cond}


def evaluate_gates(cohort: list[dict], gate_cfg: dict, sample_cond: dict) -> dict:
    """All 5m/15m gates on the single cohort. observed=None means UNKNOWN."""
    def gate(name: str, threshold, observed, passed: bool | None) -> dict:
        return {"gate": name, "threshold": threshold, "observed": observed, "result": passed}

    n = len(cohort)
    net5 = [o["P_5m"]["net_return"] for o in cohort]
    net15 = [o["P_15m"]["net_return"] for o in cohort]
    net1 = [o["P_1m"]["net_return"] for o in cohort if isinstance(o["P_1m"].get("net_return"), (int, float))]

    med5 = median(net5) if net5 else None
    med15 = median(net15) if net15 else None
    if n and med5 is not None:
        med5_gate = med5 > 0.0
    else:
        med5_gate = None
    med15_gate = (med15 > 0.0) if (n and med15 is not None) else None

    hits5 = sum(1 for v in net5 if v > 0.0)
    hits15 = sum(1 for v in net15 if v > 0.0)
    hr5 = hits5 / n if n else None
    hr15 = hits15 / n if n else None
    hr5_gate = (hr5 >= gate_cfg["hit_rate_net_positive_5m_gte"]) if hr5 is not None else None
    hr15_gate = (hr15 >= gate_cfg["hit_rate_net_positive_15m_gte"]) if hr15 is not None else None

    min_asset = effective_required("minimum_per_asset", sample_cond)
    asset_medians = {}
    both_ok = True
    both_defined = True
    for sym in TARGET_SYMBOLS:
        vals = [o["P_15m"]["net_return"] for o in cohort if o["symbol"] == sym]
        m = median(vals) if vals else None
        defined = len(vals) >= min_asset and m is not None
        asset_medians[sym] = {"count": len(vals), "median_net_15m": m, "defined": defined}
        if not defined:
            both_defined = False
        elif m <= 0.0:
            both_ok = False
    both_gate = (both_ok and both_defined) if both_defined else None

    block_med: dict[int, float] = {}
    for o in cohort:
        b = o["onset_ms"] // BLOCK_MS
        block_med.setdefault(b, []).append(o["P_15m"]["net_return"])
    block_positive = 0
    for b, vals in sorted(block_med.items()):
        if median(vals) > 0.0:
            block_positive += 1
    frac = block_positive / len(block_med) if block_med else None
    frac_gate = (frac >= gate_cfg["positive_block_median_fraction_15m_gte"]) if frac is not None else None

    gates = [
        gate("median_net_signed_return_5m_gt_zero", True, med5, med5_gate),
        gate("median_net_signed_return_15m_gt_zero", True, med15, med15_gate),
        gate("hit_rate_net_positive_5m_gte", gate_cfg["hit_rate_net_positive_5m_gte"], hr5, hr5_gate),
        gate("hit_rate_net_positive_15m_gte", gate_cfg["hit_rate_net_positive_15m_gte"], hr15, hr15_gate),
        gate("both_assets_median_net_15m_gt_zero", True,
             {"BTCUSDT": asset_medians["BTCUSDT"]["median_net_15m"],
              "ETHUSDT": asset_medians["ETHUSDT"]["median_net_15m"]}, both_gate),
        gate("positive_block_median_fraction_15m_gte", gate_cfg["positive_block_median_fraction_15m_gte"],
             frac, frac_gate),
    ]
    descriptive = {
        "median_net_signed_return_1m_descriptive_only": median(net1) if net1 else None,
        "hit_rate_net_positive_1m_descriptive_only": (
            sum(1 for v in net1 if v > 0.0) / len(net1)) if net1 else None,
    }
    return {
        "cohort_size": n,
        "gates": gates,
        "descriptive_1m": descriptive,
        "asset_detail_15m": asset_medians,
        "cohort_blocks": [
            {"block_start_ms": b * BLOCK_MS, "onsets": len(v), "median_net_15m": median(v)}
            for b, v in sorted(block_med.items())
        ],
        "any_unknown": any(g["result"] is None for g in gates),
        "any_false": any(g["result"] is False for g in gates),
    }


def effective_required(name: str, sample_cond: dict) -> int:
    for c in sample_cond["conditions"]:
        if c["name"] == name:
            return c["required"]
    raise KeyError(name)


def resolve_truth(sample_all_met: bool, gates: dict) -> dict:
    if not sample_all_met:
        status = "inconclusive"
        reason = "registered sample condition unmet or unknown"
    elif gates["any_false"]:
        status = "falsified"
        reason = "sample conditions met and at least one gate observed false"
    elif gates["any_unknown"]:
        status = "inconclusive"
        reason = "sample conditions met but at least one gate value UNKNOWN"
    else:
        status = "fast-screen-pass"
        reason = "sample conditions met and every gate observed true"
    return {
        "status": status,
        "rule": ("fast-screen-pass iff every sample condition and gate is observed "
                 "true; falsified if sample conditions pass and any gate is false; "
                 "inconclusive if a sample condition is unmet or required alignment "
                 "is unavailable"),
        "reason": reason,
        "no_claim": ("fast screens cannot establish durable alpha, replace "
                     "prospective confirmation, or authorize capital deployment; "
                     "no returns or deployment promotion is implied"),
    }


# --------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------

def collect_input_files() -> tuple[list[Path], list[Path]]:
    event_files = sorted(EVENTS_DIR.glob("liquidation-capture-*.jsonl"))
    index_files = sorted(INDEX_ROOT.glob("day*/index.jsonl"))
    return event_files, index_files


def run(contract_path: Path, out_path: Path) -> tuple[dict, int]:
    try:
        contract = load_contract(contract_path)
    except (OSError, KeyError, ValueError, json.JSONDecodeError) as exc:
        return {"error": f"contract unusable: {exc}"}, EXIT_CONTRACT

    event_files, index_files = collect_input_files()
    if not event_files or not index_files:
        return {"error": f"inputs missing: events={len(event_files)} index={len(index_files)} "
                         f"(events_dir={EVENTS_DIR} index_root={INDEX_ROOT})"}, EXIT_INPUTS

    invoked_wall = datetime.now(timezone.utc)
    snapshot_wall_iso = invoked_wall.strftime("%Y-%m-%dT%H:%M:%S.") + \
        f"{invoked_wall.microsecond // 1000:03d}Z"

    events_report = []
    event_records: list[dict] = []
    index_report = []
    index_records: list[dict] = []
    integrity_errors = 0
    for path in event_files:
        rep, parsed = snapshot_file(path)
        integrity_errors += rep["malformed_lines"]
        first = None
        latest = None
        for rec in parsed:
            body = rec.get("payload")
            body = body.get("data") if isinstance(body, dict) and isinstance(body.get("data"), dict) else body
            E = body.get("E") if isinstance(body, dict) else None
            if isinstance(E, int) and not isinstance(E, bool):
                if first is None or E < first:
                    first = E
                latest = E if latest is None else max(latest, E)
        rep["kind"] = "binance_forceOrder_capture"
        rep["first_exchange_E_ms"] = first
        rep["first_exchange_E_utc"] = iso_utc(first)
        rep["latest_exchange_E_ms"] = latest
        rep["latest_exchange_E_utc"] = iso_utc(latest)
        events_report.append(rep)
        event_records.extend(parsed)
    for path in index_files:
        rep, parsed = snapshot_file(path)
        integrity_errors += rep["malformed_lines"]
        first = None
        latest = None
        per_asset_latest: dict[str, int] = {}
        per_asset_first: dict[str, int] = {}
        for rec in parsed:
            data = rec.get("data")
            if not isinstance(data, dict):
                continue
            ts = data.get("timestamp")
            name = data.get("index_name") or rec.get("index_name") or ""
            if not isinstance(ts, int) or isinstance(ts, bool):
                continue
            asset = "btc_usd" if "btc_usd" in name else ("eth_usd" if "eth_usd" in name else None)
            if asset is None:
                continue
            if first is None or ts < first:
                first = ts
            latest = ts if latest is None else max(latest, ts)
            if asset not in per_asset_latest or ts > per_asset_latest[asset]:
                per_asset_latest[asset] = ts
            if asset not in per_asset_first or ts < per_asset_first[asset]:
                per_asset_first[asset] = ts
        rep["kind"] = "deribit_index_price_records"
        rep["first_data_timestamp_ms"] = first
        rep["first_data_timestamp_utc"] = iso_utc(first)
        rep["latest_data_timestamp_ms"] = latest
        rep["latest_data_timestamp_utc"] = iso_utc(latest)
        rep["latest_per_asset_ms"] = {a: per_asset_latest[a] for a in sorted(per_asset_latest)}
        rep["latest_per_asset_utc"] = {a: iso_utc(t) for a, t in sorted(per_asset_latest.items())}
        rep["first_per_asset_ms"] = {a: per_asset_first[a] for a in sorted(per_asset_first)}
        rep["first_per_asset_utc"] = {a: iso_utc(t) for a, t in sorted(per_asset_first.items())}
        index_report.append(rep)
        index_records.extend(parsed)
    if integrity_errors:
        return {"error": f"{integrity_errors} malformed JSON lines inside frozen prefixes; "
                         "failing closed (no partial-corpus math)"}, EXIT_INPUTS

    latest_event_E = max((r["latest_exchange_E_ms"] for r in events_report
                          if r.get("latest_exchange_E_ms") is not None), default=None)
    index_latest_ms = max((r.get("latest_data_timestamp_ms") for r in index_report
                           if r.get("latest_data_timestamp_ms") is not None), default=None)
    asset_latest: dict[str, int] = {}
    for r in index_report:
        for asset, ts in (r.get("latest_per_asset_ms") or {}).items():
            if asset not in asset_latest or ts > asset_latest[asset]:
                asset_latest[asset] = ts
    if latest_event_E is None or index_latest_ms is None:
        return {"error": "no exchange timestamps found in frozen prefixes"}, EXIT_INPUTS
    required_asset_latest = [asset_latest.get(a) for a in ("btc_usd", "eth_usd")]
    if any(t is None for t in required_asset_latest):
        return {"error": f"per-asset index latests unavailable: {asset_latest}"}, EXIT_INPUTS
    cutoff_ms = min([latest_event_E, index_latest_ms] + required_asset_latest)
    per_asset_cutoff_ms = {a: min(latest_event_E, t) for a, t in asset_latest.items()}
    deadline_ms = contract["deadline_ms"]
    search_deadline_ms = min(cutoff_ms, deadline_ms) if deadline_ms else cutoff_ms
    onsets_info = build_onsets(event_records, cutoff_ms)
    index = load_index_series(index_records, cutoff_ms)
    outcomes = compute_onset_outcomes(onsets_info["eligible_onsets"], index, cutoff_ms)
    for o in outcomes:
        if not o["cohort_complete"]:
            onsets_info["incomplete_onsets"].append({
                "event_id": o["event_id"], "symbol": o["symbol"], "side": o["side"],
                "onset_ms": o["onset_ms"],
                "missing": {h: o[f"P_{h}"].get("observed") is None for h in ("1m", "5m", "15m")}
                           | ({"initial": o["P_0"].get("observed") is None}),
            })
            for h in ("1m", "5m", "15m"):
                if o[f"P_{h}"].get("observed") is None:
                    onsets_info["incomplete_alignments"][f"horizon_{h}"] += 1
            if o["P_0"].get("observed") is None:
                onsets_info["incomplete_alignments"]["initial"] += 1

    stop_cfg = contract["sample_stop"]
    stop = find_earliest_stop(outcomes, stop_cfg, search_deadline_ms)
    cohort = stop["cohort"]
    conditions = stop["conditions"]
    sample_all_met = stop["found"]
    gates = evaluate_gates(cohort, contract["gate"], conditions)
    truth = resolve_truth(sample_all_met, gates)

    report = {
        "schema_version": 1,
        "campaign": contract["raw"].get("campaign"),
        "screen": SCREEN_ID,
        "owner_dir": "quant-trading/explorations/focused-liquidation-overlay",
        "contract": {
            "path": str(contract_path.relative_to(REPO_ROOT)),
            "sha256": contract["sha256"],
            "status": contract["raw"].get("status"),
            "created_at_utc": contract["raw"].get("created_at_utc"),
            "semantic_amendment_created_at_utc": (contract["raw"].get("semantic_amendment") or {}).get("created_at_utc"),
            "final_semantic_amendment_created_at_utc": (contract["raw"].get("final_semantic_amendment") or {}).get("created_at_utc"),
            "final_semantic_amendment_schema_clarified_at_utc": (contract["raw"].get("final_semantic_amendment") or {}).get("schema_clarified_at_utc"),
            "frozen_before_outcome_analysis": True,
            "thresholds_loaded_from_contract": {
                "sample_stop": stop_cfg,
                "gate": contract["gate"],
            },
        },
        "run": {
            "snapshot_wall_clock_utc": snapshot_wall_iso,
            "exit_code_semantics": "0 when the deterministic report is written (any truth "
                                   "status, including inconclusive); nonzero only for "
                                   "contract/input integrity failures",
            "research_only": True,
            "network_calls": 0,
            "raw_files_mutated": 0,
            "returns_or_deployment_claim": None,
        },
        "inputs": {
            "snapshot_rule": "stat byte size first, then hash and parse exactly that "
                             "byte prefix; one trailing partial line ignored and counted",
            "events": {"dir": str(EVENTS_DIR.relative_to(REPO_ROOT)), "files": events_report},
            "deribit_index": {"root": str(INDEX_ROOT.relative_to(REPO_ROOT)), "files": index_report},
            "common_cutoff": {
                "rule": "minimum latest complete exchange timestamp across the event "
                        "input and EACH required price input (btc_usd and eth_usd "
                        "index latests); never a max over the pooled index union",
                "latest_event_E_ms": latest_event_E,
                "index_latest_pooled_ms": index_latest_ms,
                "latest_per_asset_ms": {a: asset_latest[a] for a in sorted(asset_latest)},
                "latest_per_asset_utc": {a: iso_utc(t) for a, t in sorted(asset_latest.items())},
                "cutoff_ms": cutoff_ms,
                "cutoff_utc": iso_utc(cutoff_ms),
                "search_deadline_ms": search_deadline_ms,
                "per_asset_cutoff_ms": {a: per_asset_cutoff_ms[a] for a in sorted(per_asset_cutoff_ms)},
                "hard_deadline_utc": iso_utc(deadline_ms) if deadline_ms else None,
                "deadline_rule": contract["raw"]["final_semantic_amendment"]["deadline_rule"],
            },
        },
        "liquidation_corpus": {
            "parsed_events": onsets_info["parsed_events"],
            "dedup_duplicate_ids": onsets_info["dedup_duplicate_ids"],
            "event_id_mismatch_records": onsets_info["id_mismatch_records"],
            "excluded": {
                "non_forceOrder_records": onsets_info["non_forceOrder_records"],
                "missing_exchange_E": onsets_info["missing_exchange_E"],
                "non_target_symbol": onsets_info["non_target_symbol_events"],
                "after_common_cutoff": onsets_info["after_cutoff_events"],
                "orphan_outside_connect_segment": onsets_info["orphan_events_outside_segment"],
                "side_missing_or_unknown": onsets_info["side_missing_or_unknown"],
            },
            "controls": {k: onsets_info[k] for k in
                         ("control_connect", "control_session_start", "control_disconnect",
                          "control_session_stop", "control_other")},
            "segments": onsets_info["segments"],
            "segment_forceOrder_events": onsets_info["segment_events_total"],
        },
        "deribit_index": {
            "series": {k: {"observations": v["observations"],
                           "duplicate_timestamp_collapses": v["duplicate_timestamp_collapses"]}
                       for k, v in index["series"].items()},
            "dedup_duplicate_ids": index["stats"]["dedup_duplicate_ids"],
            "dedup_id_missing": index["stats"]["dedup_id_missing"],
            "dedup_id_mismatch": index["stats"]["dedup_id_mismatch"],
            "after_common_cutoff": index["stats"]["after_cutoff_observations"],
            "nonpositive_or_bad_price": index["stats"]["nonpositive_or_bad_price"],
        },
        "onsets": {
            "total_eligible": len(onsets_info["eligible_onsets"]),
            "by_symbol": onsets_info["onsets_by_symbol"],
            "by_side": onsets_info["onsets_by_side"],
            "incomplete_alignments": onsets_info["incomplete_alignments"],
            "incomplete_onsets": onsets_info["incomplete_onsets"],
            "eligible_detail": [
                {"event_id": o["event_id"], "symbol": o["symbol"], "side": o["side"],
                 "onset_ms": o["onset_ms"], "onset_utc": iso_utc(o["onset_ms"]),
                 "eligibility": o["eligibility"], "cohort_complete": o["cohort_complete"],
                 "P_0": o["P_0"], "P_1m": o["P_1m"], "P_5m": o["P_5m"], "P_15m": o["P_15m"]}
                for o in outcomes
            ],
        },
        "earliest_stop": {
            "rule": contract["raw"]["semantic_amendment"]["immutable_input_prefix"]["earliest_stop"],
            "search_deadline_utc": iso_utc(search_deadline_ms),
            "found": stop["found"],
            "candidate_stop_ms": stop["stop_ms"],
            "candidate_stop_utc": iso_utc(stop["stop_ms"]),
            "candidate_attempts": stop["tried"],
        },
        "cohort": {
            "definition": contract["raw"]["final_semantic_amendment"]["liquidation"]["single_cohort"],
            "size": len(cohort),
            "by_symbol": conditions["by_symbol"],
            "by_side": conditions["by_side"],
            "nonempty_15m_blocks": conditions["nonempty_15m_blocks"],
        },
        "sample_conditions": {
            "all_met": sample_all_met,
            "conditions": conditions["conditions"],
        },
        "gates": gates,
        "truth": truth,
        "boundaries": {
            "no_orders_or_account_access": True,
            "no_credentials": True,
            "raw_capture_untouched": True,
            "persistent_collectors_not_restarted": True,
            "oi_not_imputed": True,
            "prediction_of_this_run_used_for_sampling_extension": False,
        },
    }

    out_path.write_text(json.dumps(report, indent=2, sort_keys=False,
                                   allow_nan=False) + "\n", encoding="utf-8")
    return report, EXIT_OK


# --------------------------------------------------------------------------
# deterministic no-network selftest (synthetic fixtures only)
# --------------------------------------------------------------------------

def _write_jsonl(path: Path, rows: list[dict]) -> None:
    path.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")


def selftest() -> int:
    import shutil
    import tempfile

    global EVENTS_DIR, INDEX_ROOT, REPO_ROOT
    failures: list[str] = []

    def check(name: str, cond: bool, detail: str = "") -> None:
        if not cond:
            failures.append(f"{name}: {detail}")
        print(f"  [{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not cond else ""))

    tmp = Path(tempfile.mkdtemp(prefix="fast_intraday_selftest_"))
    try:
        # ---------------- segment / onset / dedup semantics ----------------
        ev_dir = tmp / "liquidations-persistent"
        ev_dir.mkdir()
        T0 = parse_iso_ms("2026-08-30T13:00:00.000Z")

        def ctrl(t_ms: int, ctype: str) -> dict:
            return {"v": 1, "kind": "control", "ts_utc": iso_utc(t_ms), "conn": "c0",
                    "control": {"type": ctype}}

        def ev(t_ms: int, symbol: str, side: str, e_ms: int | None = None) -> dict:
            body = {"e": "forceOrder", "E": e_ms if e_ms is not None else t_ms,
                    "o": {"s": symbol, "S": side, "o": "LIMIT", "f": "IOC", "q": "1",
                          "p": "1", "ap": "1", "X": "FILLED", "l": "0", "z": "1",
                          "T": e_ms if e_ms is not None else t_ms, "ps": symbol, "st": 1}}
            payload = {"e": "forceOrder", "E": body["E"], "o": body["o"]}
            return {"v": 1, "kind": "event", "ts_utc": iso_utc(t_ms),
                    "recv_epoch_ms": t_ms, "conn": "c1", "stream": CANON_STREAM,
                    "event_id": canonical_event_id(payload), "payload": payload}

        rows = [
            ctrl(T0 - 1000, "file_created"),
            ctrl(T0 - 900, "session_start"),
            ctrl(T0, "connect"),                     # segment 1 begins
            ev(T0 + 100_000, "BTCUSDT", "BUY"),      # <300s from connect: not onset
            ev(T0 + 500_000, "BTCUSDT", "BUY"),      # gap 400s: ONSET
            ev(T0 + 501_500, "BTCUSDT", "BUY", e_ms=T0 + 500_000),  # redelivery, same E: deduped
            ev(T0 + 590_000, "BTCUSDT", "SELL"),     # gap 90s: not onset
            ev(T0 + 1_000_000, "ETHUSDT", "SELL"),   # first, 1000s after connect: ONSET
            ev(T0 + 1_200_000, "ETHUSDT", "SELL"),   # gap 200s: not onset
            ctrl(T0 + 1_300_000, "disconnect"),      # segment 1 ends
            ev(T0 + 1_350_000, "BTCUSDT", "BUY"),    # outside segment: orphan
            ctrl(T0 + 1_400_000, "connect"),         # segment 2 begins
            ev(T0 + 1_500_000, "BTCUSDT", "BUY"),    # 100s after connect2: not onset
                                                     # (pre-disconnect silence NOT used)
            ev(T0 + 1_800_000, "BTCUSDT", "SELL"),   # 400s after connect2: ONSET
            ctrl(T0 + 1_900_000, "session_start"),   # resets, does not begin
            ev(T0 + 2_300_000, "BTCUSDT", "BUY"),    # no open segment: orphan
            ctrl(T0 + 2_400_000, "connect"),         # segment 3 begins
            ev(T0 + 2_500_000, "BTCUSDT", "BUY"),    # 100s: not onset
            ev(T0 + 2_900_000, "XRPUSDT", "BUY"),    # non-target symbol
        ]
        _write_jsonl(ev_dir / "liquidation-capture-20260830.jsonl", rows)

        on = build_onsets(rows, cutoff_ms=10**13)
        got = [(o["symbol"], o["side"], o["onset_ms"]) for o in on["eligible_onsets"]]
        want = [("BTCUSDT", "BUY", T0 + 500_000),
                ("ETHUSDT", "SELL", T0 + 1_000_000),
                ("BTCUSDT", "SELL", T0 + 1_800_000)]
        check("onsets: connect-gated first events and >=300s same-symbol gaps",
              got == want, f"got {got} want {want}")
        check("dedup: duplicate delivery collapsed to first occurrence",
              on["dedup_duplicate_ids"] == 1 and on["parsed_events"] == len(rows),
              f"dupes={on['dedup_duplicate_ids']} parsed={on['parsed_events']} want={len(rows)}")
        check("segments: session_start resets without beginning; orphans counted",
              on["orphan_events_outside_segment"] == 2 and on["control_session_start"] == 2,
              f"orphans={on['orphan_events_outside_segment']} starts={on['control_session_start']}")
        check("non-target symbols excluded by count",
              on["non_target_symbol_events"] == 1,
              f"got {on['non_target_symbol_events']}")

        # ---------------- alignment windows + formula ----------------
        idx0 = T0 + 500_000                        # BTC onset E
        idx0_eth = T0 + 1_000_000                  # ETH onset E
        idx_rows = []
        for ts, price, key in [
            (idx0 - 5_000, 100.0, "btc_usd"),      # exactly 5000ms old: OK
            (idx0 + 300_000, 101.0, "btc_usd"),    # horizon target exactly: OK
            (idx0 + 900_000 + 5_000, 102.0, "btc_usd"),  # exactly 5000ms late: OK
            (idx0_eth - 5_001, 200.0, "eth_usd"),  # 5001ms old: too old
            (idx0_eth + 605_005, 201.0, "eth_usd"),  # 5ms late: OK
            (idx0_eth + 1_050_005, 205.0, "eth_usd"),  # 5ms late: OK
        ]:
            did = hashlib.sha256(f"index|{key}|{ts}".encode()).hexdigest()
            idx_rows.append({"record_type": "index_price", "dedup_id": did,
                             "received_ts_ms": ts + 3, "index_name": key,
                             "data": {"timestamp": ts, "price": price,
                                      "index_name": key}})
        # duplicate timestamp keeps LAST stable occurrence
        did = hashlib.sha256(f"index|btc_usd|{idx0 - 5_000}".encode()).hexdigest()
        idx_rows.append({"record_type": "index_price", "dedup_id": did + "x",
                         "received_ts_ms": idx0 - 5_000 + 1, "index_name": "btc_usd",
                         "data": {"timestamp": idx0 - 5_000, "price": 100.5,
                                  "index_name": "btc_usd"}})

        ix = load_index_series(idx_rows, cutoff_ms=10**13)
        check("index: deterministic dedup ids verified and duplicates collapsed",
              ix["stats"]["dedup_id_mismatch"] == 1
              and ix["series"]["btc_usd"]["duplicate_timestamp_collapses"] == 1,
              f"mismatch={ix['stats']['dedup_id_mismatch']}")

        buy = {"event_id": "buy1", "symbol": "BTCUSDT", "side": "BUY",
               "onset_ms": idx0, "segment_connect_ms": None, "eligibility": "t",
               "prior_same_symbol_E_ms": None}
        sell = {"event_id": "sell1", "symbol": "ETHUSDT", "side": "SELL",
                "onset_ms": idx0_eth, "segment_connect_ms": None, "eligibility": "t",
                "prior_same_symbol_E_ms": None}
        outs = compute_onset_outcomes([buy, sell], ix, 10**13)
        b, s = outs
        check("P_0: exactly-5000ms-old accepted; duplicate ts keeps last price",
              b["P_0"]["observed_ts_ms"] == idx0 - 5_000 and b["P_0"]["price"] == 100.5,
              f"got {b['P_0']}")
        check("P_5m: horizon target boundary accepted",
              b["P_5m"]["observed_ts_ms"] == idx0 + 300_000 and b["P_5m"]["price"] == 101.0,
              f"got {b['P_5m']}")
        check("P_15m: exactly-5000ms-late accepted",
              b["P_15m"]["observed_ts_ms"] == idx0 + 905_000, f"got {b['P_15m']}")
        net5 = 1 * (101.0 / 100.5 - 1) - 0.0002
        net15 = 1 * (102.0 / 100.5 - 1) - 0.0002
        check("net formula: direction*(P/P0-1)-0.0002 unrounded",
              abs(b["P_5m"]["net_return"] - net5) < 1e-15
              and abs(b["P_15m"]["net_return"] - net15) < 1e-15, f"got {b['P_5m']}")
        check("SELL direction sign: net = -(P/P0-1)-0.0002",
              s["P_0"]["observed"] is None and s["cohort_complete"] is False,
              f"got P_0={s['P_0']} complete={s['cohort_complete']}")
        late = {"event_id": "late1", "symbol": "BTCUSDT", "side": "BUY",
                "onset_ms": T0 + 900_000, "segment_connect_ms": None,
                "eligibility": "t", "prior_same_symbol_E_ms": None}
        outs2 = compute_onset_outcomes([late], ix, T0 + 1_200_000)   # cutoff before 15m target
        check("horizon beyond common cutoff is unavailable",
              outs2[0]["P_15m"].get("missing") == "beyond_common_cutoff",
              f"got {outs2[0]['P_15m']}")

        # alignment strictness: 5001ms-old initial rejected
        old = {"event_id": "old1", "symbol": "BTCUSDT", "side": "BUY",
               "onset_ms": idx0 + 1, "segment_connect_ms": None, "eligibility": "t",
               "prior_same_symbol_E_ms": None}
        outs3 = compute_onset_outcomes([old], ix,
                                       cutoff_ms=idx0 + 2_000_000)
        check("P_0: observation 5001ms old rejected",
              outs3[0]["P_0"].get("observed") is None and not outs3[0]["cohort_complete"],
              f"got {outs3[0]['P_0']}")

        # ---------------- truth resolution + earliest stop ----------------
        coh = []
        for i in range(6):
            o = dict(b)
            o["event_id"] = f"e{i}"
            o = json.loads(json.dumps(o))
            o["P_15m"]["net_return"] = 0.001 if i % 2 == 0 else -0.001
            o["P_5m"]["net_return"] = 0.001 if i % 2 == 0 else -0.001
            coh.append(o)
        g = evaluate_gates(coh, {"hit_rate_net_positive_5m_gte": 0.6,
                                 "hit_rate_net_positive_15m_gte": 0.6,
                                 "positive_block_median_fraction_15m_gte": 0.6},
                           sample_conditions(coh, {"minimum_complete_onsets_at_15m": 20,
                                                   "minimum_per_asset": 5, "minimum_per_side": 5,
                                                   "minimum_nonempty_15m_onset_blocks": 6}))
        med15 = g["gates"][1]["observed"]
        check("median: mean of two middles; zero median not positive (5m and 15m)",
              med15 == 0.0 and g["gates"][1]["result"] is False
              and g["gates"][0]["observed"] == 0.0 and g["gates"][0]["result"] is False,
              f"med15={med15} med5={g['gates'][0]['observed']}")
        check("hit fractions: 3/6 = 0.5 < 0.6 false; zero return is not a hit",
              g["gates"][2]["observed"] == 0.5 and g["gates"][2]["result"] is False
              and g["gates"][3]["observed"] == 0.5 and g["gates"][3]["result"] is False,
              f"{g['gates'][2]} {g['gates'][3]}")
        even_pos = json.loads(json.dumps(coh[:2]))
        even_pos[0]["P_5m"]["net_return"] = 0.001
        even_pos[1]["P_5m"]["net_return"] = 0.002
        g_even = evaluate_gates(even_pos, {"hit_rate_net_positive_5m_gte": 0.6,
                                           "hit_rate_net_positive_15m_gte": 0.6,
                                           "positive_block_median_fraction_15m_gte": 0.6},
                                sample_conditions(even_pos, {"minimum_complete_onsets_at_15m": 20,
                                                             "minimum_per_asset": 5,
                                                             "minimum_per_side": 5,
                                                             "minimum_nonempty_15m_onset_blocks": 6}))
        check("median: even-n mean of middles (0.0015) > 0 -> gate true",
              g_even["gates"][0]["observed"] == 0.0015 and g_even["gates"][0]["result"] is True,
              f"{g_even['gates'][0]}")
        blk = json.loads(json.dumps(coh[:4]))
        blk[0]["onset_ms"] = 0; blk[1]["onset_ms"] = 1            # block 0
        blk[2]["onset_ms"] = 900_000; blk[3]["onset_ms"] = 900_001  # block 1
        blk[0]["P_15m"]["net_return"] = 0.001; blk[1]["P_15m"]["net_return"] = -0.001
        blk[2]["P_15m"]["net_return"] = 0.001; blk[3]["P_15m"]["net_return"] = 0.002
        g_blk = evaluate_gates(blk, {"hit_rate_net_positive_5m_gte": 0.6,
                                     "hit_rate_net_positive_15m_gte": 0.6,
                                     "positive_block_median_fraction_15m_gte": 0.6},
                               sample_conditions(blk, {"minimum_complete_onsets_at_15m": 20,
                                                       "minimum_per_asset": 5,
                                                       "minimum_per_side": 5,
                                                       "minimum_nonempty_15m_onset_blocks": 6}))
        check("block fraction: only positive-median blocks count (1/2 = 0.5 < 0.6)",
              g_blk["gates"][5]["observed"] == 0.5 and g_blk["gates"][5]["result"] is False,
              f"{g_blk['gates'][5]}")
        check("any_false -> falsified", resolve_truth(True, g)["status"] == "falsified")
        check("sample unmet -> inconclusive (gates not shown as false)",
              resolve_truth(False, g)["status"] == "inconclusive")
        stop_cfg_min = {"minimum_complete_onsets_at_15m": 20, "minimum_per_asset": 5,
                        "minimum_per_side": 5, "minimum_nonempty_15m_onset_blocks": 6}
        g2 = evaluate_gates(coh[:1], {"hit_rate_net_positive_5m_gte": 0.6,
                                      "hit_rate_net_positive_15m_gte": 0.6,
                                      "positive_block_median_fraction_15m_gte": 0.6},
                            sample_conditions(coh[:1], stop_cfg_min))
        check("empty/under-minimum cohort => UNKNOWN serialized as null",
              g2["gates"][4]["result"] is None and g2["any_unknown"] is True
              and json.dumps(g2).count("null") >= 2, f"{g2['gates'][4]}")

        # earliest stop: first availability instant meeting every condition
        stop_cfg = {"minimum_complete_onsets_at_15m": 3, "minimum_per_asset": 3,
                    "minimum_per_side": 3, "minimum_nonempty_15m_onset_blocks": 3}
        onsets = []
        for i in range(6):
            row = json.loads(json.dumps(b))
            row["event_id"] = f"s{i}"
            row["symbol"] = "BTCUSDT" if i < 3 else "ETHUSDT"
            row["side"] = "BUY" if i < 3 else "SELL"
            row["onset_ms"] = idx0 + 960_000 * i          # distinct 15m blocks
            row["cohort_complete"] = True
            row["availability_ms"] = idx0 + 900_000 + 60_000 * i
            onsets.append(row)
        st = find_earliest_stop(onsets, stop_cfg, deadline_ms=idx0 + 2_000_000)
        check("earliest stop found at first all-met availability instant (n=6)",
              st["found"] and st["stop_ms"] == idx0 + 900_000 + 300_000
              and len(st["cohort"]) == 6,
              f"found={st['found']} stop={st['stop_ms']} n={len(st['cohort'])}")
        st2 = find_earliest_stop(onsets, stop_cfg, deadline_ms=idx0 + 900_000)
        check("deadline before availability => full-prefix inconclusive",
              st2["found"] is False and st2["stop_ms"] is None,
              f"found={st2['found']} n={len(st2['cohort'])}")

        # common cutoff: min(latest event E, latest btc_usd, latest eth_usd) with
        # per-asset asset splitting exactly as run() computes it
        probe_rows = []
        probe_event_latest = idx0 + 1_200_000
        for key, latest_ts in [("btc_usd", idx0 + 2_000_000), ("eth_usd", idx0 + 1_500_000)]:
            did = hashlib.sha256(f"index|{key}|{latest_ts}".encode()).hexdigest()
            probe_rows.append({"record_type": "index_price", "dedup_id": did,
                               "received_ts_ms": latest_ts + 1, "index_name": key,
                               "data": {"timestamp": latest_ts, "price": 1.0,
                                        "index_name": key}})
        rep_probe = {"latest_exchange_E_ms": probe_event_latest,
                     "latest_per_asset_ms": {"btc_usd": idx0 + 2_000_000,
                                             "eth_usd": idx0 + 1_500_000}}
        asset_probe: dict[str, int] = {}
        for rr in [rep_probe]:
            for asset, ts in (rr.get("latest_per_asset_ms") or {}).items():
                if asset not in asset_probe or ts > asset_probe[asset]:
                    asset_probe[asset] = ts
        probe_cutoff = min([rep_probe["latest_exchange_E_ms"]]
                           + [asset_probe[a] for a in ("btc_usd", "eth_usd")])
        pooled = max(rep_probe["latest_exchange_E_ms"], asset_probe["btc_usd"],
                     asset_probe["eth_usd"])
        check("common cutoff: min(event E, btc_usd latest, eth_usd latest) "
              "NOT the pooled max",
              probe_cutoff == idx0 + 1_200_000 and pooled == idx0 + 2_000_000
              and probe_cutoff < pooled,
              f"cutoff={probe_cutoff} pooled={pooled}")

        # ---------------- byte-prefix snapshot semantics ----------------
        p = ev_dir / "liquidation-capture-20260830.jsonl"
        snap2_rep, snap2_rows = snapshot_file(p)    # clean file: no partial line
        full = p.read_bytes()
        check("snapshot: clean file parsed whole; sha256 over exact bytes",
              snap2_rep["prefix_bytes"] == len(full)
              and snap2_rep["trailing_partial_line_bytes"] == 0
              and snap2_rep["sha256"] == hashlib.sha256(full).hexdigest(),
              f"prefix={snap2_rep['prefix_bytes']} partial={snap2_rep['trailing_partial_line_bytes']}")
        with p.open("ab") as fh:
            fh.write(b'{"v":1,"kind":"event","trunc')      # dangling partial line
        snap_rep, snap_rows = snapshot_file(p)
        with_partial = p.read_bytes()
        check("snapshot: prefix bytes hashed exactly; one trailing partial ignored",
              snap_rep["prefix_bytes"] == len(with_partial)
              and snap_rep["trailing_partial_line_bytes"] == len(b'{"v":1,"kind":"event","trunc')
              and snap_rep["sha256"] == hashlib.sha256(with_partial).hexdigest()
              and len(snap_rows) == snap2_rep["json_objects"],
              f"prefix={snap_rep['prefix_bytes']} partial={snap_rep['trailing_partial_line_bytes']}")
        p.write_bytes(full)                         # restore for deterministic state

        print()
        if failures:
            print(f"SELFTEST FAILED ({len(failures)}):")
            for f in failures:
                print(f"  - {f}")
            return 1
        print("ALL SELFTEST PROOFS PASS (deterministic, no network)")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    global EVENTS_DIR, INDEX_ROOT
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--events-dir", type=Path, default=None)
    ap.add_argument("--index-root", type=Path, default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="deterministic no-network proofs on synthetic fixtures only")
    args = ap.parse_args()
    if args.events_dir:
        EVENTS_DIR = args.events_dir
    if args.index_root:
        INDEX_ROOT = args.index_root
    if args.selftest:
        return selftest()
    report, code = run(args.contract, args.out)
    if code == EXIT_OK:
        truth = report.get("truth", {})
        print(f"FROZEN_SCREEN_REPORT {args.out}")
        print(f" screen={SCREEN_ID} status={truth.get('status')} "
              f"cohort={report.get('cohort', {}).get('size')} "
              f"sample_all_met={report.get('sample_conditions', {}).get('all_met')} "
              f"stop_found={report.get('earliest_stop', {}).get('found')}")
        return EXIT_OK
    print(json.dumps(report), file=sys.stderr)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
