#!/usr/bin/env python3
"""Frozen two-window test of calendar-spread convergence on Binance coin-margined quarterly delivery futures.

Mechanism (cycle k=7). For three consecutive quarterly contracts (prev, near, far) on one base:
enter at the scheduled contract hand-over (first complete UTC day after prev delivery with both legs
tradable at the 23:00 close), take the position AGAINST the sign of c = ln(P_far / P_near), and exit at
the near contract's last complete UTC day before its delivery. Inverse-contract arithmetic in base
units, four taker legs charged at their own transaction prices.

All frozen hashes and upstream archive checksums are verified before any price is parsed. Results use
exclusive-create semantics. A pass is replicated; a fail closes the family. Net-of-cost expected-return
evidence on frozen historical calendars, not deployable alpha.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import zipfile
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
EXPLORATIONS = WORKSPACE_ROOT / "quant-trading" / "explorations"
CONTRACT_PATH = EXPLORATIONS / "pre_statement_delivery_term_structure.json"
INVENTORY_PATH = EXPLORATIONS / "delivery-term-structure-inputs.json"
AUDIT_PATH = EXPLORATIONS / "delivery-term-structure-prior-window-audit.json"
POLICY_PATH = EXPLORATIONS / "frontier-loop-policy.json"

HEADER = ("open_time,open,high,low,close,volume,close_time,quote_volume,count,"
          "taker_buy_volume,taker_buy_quote_volume,ignore")
EXPECTED_FIELDS = 12
HOUR_MS = 3_600_000
DAY_MS = 86_400_000

TAKER_FEE = 0.0005
HOURS_PER_DAY = 24
CLOSE_HOUR = 23
MAX_DELIVERY_GAP_DAYS = 100
MAGNITUDE_THRESHOLD = 0.0010
QUARTER_STABILITY_NUMERATOR = 3
QUARTER_STABILITY_DENOMINATOR = 5
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20_260_903
BOOTSTRAP_BLOCK_QUARTERS = 4
MIN_POSITIONS = 70
MIN_QUARTERS = 8
MIN_BASES = 6
MIN_SEGMENT_POSITIONS = 8
MAJORS = ("BTCUSD", "ETHUSD")
LOOP_ALPHA = 0.05 / (7 * 8)

FROZEN_PARAMETERS = {
    "TAKER_FEE": TAKER_FEE, "HOURS_PER_DAY": HOURS_PER_DAY, "CLOSE_HOUR": CLOSE_HOUR,
    "MAX_DELIVERY_GAP_DAYS": MAX_DELIVERY_GAP_DAYS, "MAGNITUDE_THRESHOLD": MAGNITUDE_THRESHOLD,
    "QUARTER_STABILITY_NUMERATOR": QUARTER_STABILITY_NUMERATOR,
    "QUARTER_STABILITY_DENOMINATOR": QUARTER_STABILITY_DENOMINATOR,
    "BOOTSTRAP_REPLICATES": BOOTSTRAP_REPLICATES, "BOOTSTRAP_SEED": BOOTSTRAP_SEED,
    "BOOTSTRAP_BLOCK_QUARTERS": BOOTSTRAP_BLOCK_QUARTERS, "MIN_POSITIONS": MIN_POSITIONS,
    "MIN_QUARTERS": MIN_QUARTERS, "MIN_BASES": MIN_BASES,
    "MIN_SEGMENT_POSITIONS": MIN_SEGMENT_POSITIONS, "loop_alpha": LOOP_ALPHA,
}
PIN_FIELDS = ("input_inventory_sha256", "prior_window_audit_sha256", "policy_sha256", "runner_sha256")
REJECTION_KEYS = ("no-complete-entry-day", "no-complete-exit-day", "late-handover-entry",
                  "entry-not-before-exit", "zero-spread-at-entry", "nonpositive-or-nonfinite-close")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _inside_root(path: Path, root: Path) -> Path:
    resolved = Path(path).resolve()
    if not resolved.is_relative_to(Path(root).resolve()):
        raise ValueError(f"input escapes workspace root: {path}")
    return resolved


def check_contract_parameters(contract: dict) -> None:
    """The contract's parameter block must equal this module's constants, exactly."""
    declared = contract["frozen_rule"]["parameters"]
    if declared != FROZEN_PARAMETERS:
        differing = {k: (declared.get(k), FROZEN_PARAMETERS.get(k))
                     for k in set(declared) | set(FROZEN_PARAMETERS)
                     if declared.get(k) != FROZEN_PARAMETERS.get(k)}
        raise ValueError(f"contract/runner parameter mismatch: {differing}")
    if contract["loop_multiplicity"]["alpha_k"] != LOOP_ALPHA:
        raise ValueError("contract alpha_k does not match the runner loop alpha")
    if "USD_" not in contract["frozen_rule"]["universe"]:
        raise ValueError("contract universe clause does not describe coin-margined dated contracts")


def delivery_date(symbol: str) -> date:
    stamp = symbol.split("_", 1)[1]
    if len(stamp) != 6 or not stamp.isdigit() or not stamp.isascii():
        raise ValueError(f"malformed delivery stamp: {symbol!r}")
    year, month, day = 2000 + int(stamp[:2]), int(stamp[2:4]), int(stamp[4:6])
    if month not in (3, 6, 9, 12):
        raise ValueError(f"non-quarterly delivery month: {symbol!r}")
    return date(year, month, day)


def quarter_label(day: date) -> str:
    return f"{day.year}Q{(day.month - 1) // 3 + 1}"


def derive_pair_universe(symbols: list[str]) -> list[dict]:
    """Rebuild the (prev, near, far) candidate set from contract names alone."""
    by_base: dict[str, list[str]] = {}
    for entry in symbols:
        margin, _, symbol = entry.partition(":")
        if margin != "cm":
            continue
        by_base.setdefault(symbol.split("_", 1)[0], []).append(symbol)
    out: list[dict] = []
    for base in sorted(by_base):
        ordered = sorted(by_base[base], key=delivery_date)
        for index in range(1, len(ordered) - 1):
            prev, near, far = ordered[index - 1], ordered[index], ordered[index + 1]
            if (delivery_date(near) - delivery_date(prev)).days > MAX_DELIVERY_GAP_DAYS:
                continue
            if (delivery_date(far) - delivery_date(near)).days > MAX_DELIVERY_GAP_DAYS:
                continue
            out.append({"margin": "cm", "base": base, "near": near, "far": far,
                        "prev_delivery": delivery_date(prev).isoformat(),
                        "near_delivery": delivery_date(near).isoformat(),
                        "quarter": quarter_label(delivery_date(near))})
    return out


def verify_archive(entry: dict, root: Path = WORKSPACE_ROOT) -> str:
    archive_path = _inside_root(root / entry["path"], root)
    sidecar_path = _inside_root(root / entry["checksum_sidecar"], root)
    if not archive_path.is_file() or not sidecar_path.is_file():
        raise ValueError(f"archive or sidecar absent: {entry['path']}")
    if archive_path.stat().st_size != entry["size_bytes"]:
        raise ValueError(f"archive byte size drift: {entry['path']}")
    if sha256_file(archive_path) != entry["sha256"]:
        raise ValueError(f"archive hash drift: {entry['path']}")
    if sha256_file(sidecar_path) != entry["checksum_sidecar_sha256"]:
        raise ValueError(f"sidecar file hash drift: {entry['checksum_sidecar']}")
    fields = sidecar_path.read_text(encoding="utf-8").split()
    if len(fields) < 2:
        raise ValueError(f"malformed upstream sidecar: {entry['checksum_sidecar']}")
    if fields[0].lower() != entry["sha256"]:
        raise ValueError(f"upstream sidecar hash mismatch: {entry['path']}")
    if Path(fields[1]).name != archive_path.name:
        raise ValueError(f"upstream sidecar filename mismatch: {entry['path']}")
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.namelist()
        if members != [entry["zip_member"]]:
            raise ValueError(f"zip member mismatch: {entry['path']}")
        if archive.testzip() is not None:
            raise ValueError(f"zip CRC failure: {entry['path']}")
    return entry["sha256"]


def verify_static_inputs(contract: dict) -> tuple[dict, int]:
    """Verify every pin and every archive before any price is parsed."""
    check_contract_parameters(contract)
    expected = {
        "input_inventory_sha256": contract["assets_and_data"]["input_inventory_sha256"],
        "prior_window_audit_sha256": contract["novelty_and_consumption"]["prior_window_audit_sha256"],
        "policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
    }
    observed = {"input_inventory_sha256": sha256_file(INVENTORY_PATH),
                "prior_window_audit_sha256": sha256_file(AUDIT_PATH),
                "policy_sha256": sha256_file(POLICY_PATH)}
    for field, value in expected.items():
        if observed[field] != value:
            raise ValueError(f"{field} mismatch: expected {value}, observed {observed[field]}")
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("candidate_outcomes_accessed") is not False:
        raise ValueError("input inventory does not certify pre-outcome construction")
    if inventory["archive_count"] != len(inventory["files"]):
        raise ValueError("inventory archive count mismatch")
    if len({f["path"] for f in inventory["files"]}) != len(inventory["files"]):
        raise ValueError("inventory contains duplicate paths")
    derived = derive_pair_universe(inventory["symbols"])
    if derived != inventory["pair_universe"]:
        raise ValueError("recomputed pair universe does not match the pinned pair universe")
    if len(derived) != inventory["pair_universe_count"]:
        raise ValueError("pair universe count mismatch")
    for entry in inventory["files"]:
        verify_archive(entry)
    return inventory, len(inventory["files"])


def _chain_first(first: list[str], reader):
    yield first
    yield from reader


def load_contract_days(inventory: dict, root: Path = WORKSPACE_ROOT) -> dict[tuple[str, str], dict]:
    """Per contract, per UTC day: hours seen, the 23:00 close, and the 23:00 trade count.

    Reads open_time, close and count only. Duplicate timestamps, off-hour-grid rows, malformed
    field counts and schema drift all abort; nothing is imputed or repaired.
    """
    days: dict[tuple[str, str], dict] = {}
    for entry in inventory["files"]:
        key = (entry["margin"], entry["symbol"])
        table = days.setdefault(key, {})
        seen: set[int] = set()
        archive_path = _inside_root(root / entry["path"], root)
        with zipfile.ZipFile(archive_path) as archive:
            with archive.open(entry["zip_member"]) as raw:
                reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8"))
                first = next(reader, None)
                if first is None:
                    raise ValueError(f"empty archive: {entry['path']}")
                header_present = ",".join(first) == HEADER
                if header_present != entry["header_present"]:
                    raise ValueError(f"header form drift: {entry['path']}")
                rows = reader if header_present else _chain_first(first, reader)
                count = 0
                for line_number, row in enumerate(rows, 2):
                    if len(row) != EXPECTED_FIELDS:
                        raise ValueError(f"schema column count at {entry['path']}:{line_number}")
                    stamp = int(row[0])
                    if stamp % HOUR_MS:
                        raise ValueError(f"off-hour timestamp at {entry['path']}:{line_number}")
                    if stamp in seen:
                        raise ValueError(f"duplicate timestamp at {entry['path']}:{line_number}")
                    seen.add(stamp)
                    count += 1
                    day = datetime.fromtimestamp(stamp / 1000, timezone.utc).date().isoformat()
                    hour = (stamp % DAY_MS) // HOUR_MS
                    slot = table.setdefault(day, {"hours": set(), "close": None, "trades": None})
                    slot["hours"].add(hour)
                    if hour == CLOSE_HOUR:
                        slot["close"] = float(row[4])
                        slot["trades"] = int(row[8])
                if count != entry["rows"]:
                    raise ValueError(f"row count drift: {entry['path']}")
    return days


def tradable_days(table: dict) -> set[str]:
    """A day is tradable iff all 24 hourly bars exist and the 23:00 bar has a real trade.

    Deliberately structural: bar presence and trade count only. A price VALUE never selects a day.
    The `close is not None` test asserts the 23:00 bar exists, which 24 present hours already imply;
    an invalid close value rejects the whole candidate in position_return (nonpositive-or-nonfinite-close)
    rather than silently shifting entry or exit to a neighbouring day.
    """
    out = set()
    for day, slot in table.items():
        if len(slot["hours"]) != HOURS_PER_DAY:
            continue
        if slot["trades"] is None or slot["trades"] <= 0:
            continue
        if slot["close"] is None:
            continue
        out.add(day)
    return out


def locate_position(candidate: dict, days: dict) -> dict:
    near = days.get((candidate["margin"], candidate["near"]), {})
    far = days.get((candidate["margin"], candidate["far"]), {})
    both = tradable_days(near) & tradable_days(far)
    prev_delivery = date.fromisoformat(candidate["prev_delivery"])
    near_delivery = date.fromisoformat(candidate["near_delivery"])
    after = sorted(d for d in both if date.fromisoformat(d) > prev_delivery)
    before = sorted(d for d in both if date.fromisoformat(d) < near_delivery)
    if not after:
        return {"rejected": "no-complete-entry-day"}
    if not before:
        return {"rejected": "no-complete-exit-day"}
    entry_day, exit_day = after[0], before[-1]
    midpoint = prev_delivery + (near_delivery - prev_delivery) / 2
    if date.fromisoformat(entry_day) > midpoint:
        return {"rejected": "late-handover-entry"}
    if not entry_day < exit_day:
        return {"rejected": "entry-not-before-exit"}
    return {"entry_day": entry_day, "exit_day": exit_day,
            "entry_lag_days": (date.fromisoformat(entry_day) - prev_delivery).days,
            "hold_days": (date.fromisoformat(exit_day) - date.fromisoformat(entry_day)).days,
            "near_entry": near[entry_day]["close"], "near_exit": near[exit_day]["close"],
            "far_entry": far[entry_day]["close"], "far_exit": far[exit_day]["close"]}


def position_return(located: dict) -> dict:
    """Exact inverse-contract calendar spread, equal contract counts, multiplier cancelled."""
    pn_e, pn_x = located["near_entry"], located["near_exit"]
    pf_e, pf_x = located["far_entry"], located["far_exit"]
    for value in (pn_e, pn_x, pf_e, pf_x):
        if not math.isfinite(value) or value <= 0.0:
            return {"rejected": "nonpositive-or-nonfinite-close"}
    spread_entry = math.log(pf_e / pn_e)
    if spread_entry == 0.0:
        return {"rejected": "zero-spread-at-entry"}
    if spread_entry > 0.0:                      # far is rich: long near, short far
        long_entry, long_exit, short_entry, short_exit, side = pn_e, pn_x, pf_e, pf_x, "long-near"
    else:                                       # near is rich: long far, short near
        long_entry, long_exit, short_entry, short_exit, side = pf_e, pf_x, pn_e, pn_x, "long-far"
    pnl = (1.0 / long_entry - 1.0 / long_exit) - (1.0 / short_entry - 1.0 / short_exit)
    fees = TAKER_FEE * (1.0 / long_entry + 1.0 / short_entry + 1.0 / long_exit + 1.0 / short_exit)
    gross_exposure = 1.0 / long_entry + 1.0 / short_entry
    return {"side": side, "spread_entry": spread_entry, "spread_exit": math.log(pf_x / pn_x),
            "gross_return": pnl / gross_exposure, "net_return": (pnl - fees) / gross_exposure,
            "fee_drag": fees / gross_exposure}


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def block_bootstrap_p(by_quarter: dict[str, list[float]], replicates: int = BOOTSTRAP_REPLICATES,
                      seed: int = BOOTSTRAP_SEED, block: int = BOOTSTRAP_BLOCK_QUARTERS) -> float | None:
    """One-sided moving-block bootstrap over ordered delivery quarters; positions pooled."""
    quarters = sorted(by_quarter)
    n = len(quarters)
    if n < block:
        return None
    rng = random.Random(seed)
    blocks_per_draw = -(-n // block)
    at_or_below_zero = 0
    for _ in range(replicates):
        pooled: list[float] = []
        picked = 0
        for _ in range(blocks_per_draw):
            start = rng.randrange(n)
            for offset in range(block):
                if picked >= n:
                    break
                pooled.extend(by_quarter[quarters[(start + offset) % n]])
                picked += 1
        if pooled and _mean(pooled) <= 0.0:
            at_or_below_zero += 1
    return (1 + at_or_below_zero) / (replicates + 1)


def evaluate_gates(rows: list[dict]) -> dict:
    nets = [r["net_return"] for r in rows]
    mean_net = _mean(nets)
    by_quarter: dict[str, list[float]] = {}
    for row in rows:
        by_quarter.setdefault(row["quarter"], []).append(row["net_return"])
    positive_quarters = sum(1 for values in by_quarter.values() if _mean(values) > 0.0)
    segments = {"major": [r["net_return"] for r in rows if r["base"] in MAJORS],
                "alt": [r["net_return"] for r in rows if r["base"] not in MAJORS]}
    seg_means = {k: (_mean(v) if v else None) for k, v in segments.items()}
    seg_pass = all(seg_means[k] is not None and seg_means[k] > 0.0
                   for k, v in segments.items() if len(v) >= MIN_SEGMENT_POSITIONS)
    evaluated = [k for k, v in segments.items() if len(v) >= MIN_SEGMENT_POSITIONS]
    p_value = block_bootstrap_p(by_quarter)
    return {
        "G1_net_magnitude": {"pass": mean_net >= MAGNITUDE_THRESHOLD, "mean_net_return": mean_net,
                             "threshold": MAGNITUDE_THRESHOLD},
        "G2_segment_consistency": {"pass": bool(seg_pass) and len(evaluated) == 2,
                                   "segment_mean_net": seg_means,
                                   "segment_positions": {k: len(v) for k, v in segments.items()},
                                   "segments_evaluated": evaluated},
        "G3_quarter_stability": {
            "pass": positive_quarters * QUARTER_STABILITY_DENOMINATOR
                    >= len(by_quarter) * QUARTER_STABILITY_NUMERATOR,
            "positive_quarters": positive_quarters, "included_quarters": len(by_quarter),
            "rule": "positive_quarters * 5 >= included_quarters * 3"},
        "G4_dependence_robustness": {
            "pass": p_value is not None and p_value < LOOP_ALPHA, "one_sided_block_bootstrap_p": p_value,
            "alpha": LOOP_ALPHA, "replicates": BOOTSTRAP_REPLICATES, "seed": BOOTSTRAP_SEED,
            "block_quarters": BOOTSTRAP_BLOCK_QUARTERS},
    }


def run_window(contract: dict, label: str, inventory: dict, days: dict) -> dict:
    window = contract["windows_by_near_delivery"][label]
    low, high = window["start_inclusive"], window["end_exclusive"]
    rejected = {key: 0 for key in REJECTION_KEYS}
    rows: list[dict] = []
    for candidate in inventory["pair_universe"]:
        if not low <= candidate["near_delivery"] < high:
            continue
        located = locate_position(candidate, days)
        if "rejected" in located:
            rejected[located["rejected"]] += 1
            continue
        priced = position_return(located)
        if "rejected" in priced:
            rejected[priced["rejected"]] += 1
            continue
        rows.append({**candidate, **{k: located[k] for k in
                                     ("entry_day", "exit_day", "entry_lag_days", "hold_days")}, **priced})
    bases = sorted({r["base"] for r in rows})
    quarters = sorted({r["quarter"] for r in rows})
    majors = sum(1 for r in rows if r["base"] in MAJORS)
    holds = sorted(r["hold_days"] for r in rows)
    sample = {"positions": len(rows), "included_quarters": len(quarters), "bases": len(bases),
              "base_list": bases, "quarter_list": quarters, "major_positions": majors,
              "alt_positions": len(rows) - majors,
              "hold_days_min_median_max": [holds[0], holds[len(holds) // 2], holds[-1]] if holds else [],
              "sides": {side: sum(1 for r in rows if r["side"] == side) for side in ("long-near", "long-far")}}
    minimum_met = (len(rows) >= MIN_POSITIONS and len(quarters) >= MIN_QUARTERS and len(bases) >= MIN_BASES
                   and majors >= MIN_SEGMENT_POSITIONS and len(rows) - majors >= MIN_SEGMENT_POSITIONS)
    if not minimum_met:
        return {"window": label, "status": "inconclusive-insufficient-sample", "sample_minimum_met": False,
                "sample": sample, "rejected": rejected, "candidates_in_window":
                    sum(1 for c in inventory["pair_universe"] if low <= c["near_delivery"] < high)}
    gates = evaluate_gates(rows)
    by_quarter = {}
    for row in rows:
        by_quarter.setdefault(row["quarter"], []).append(row["net_return"])
    diagnostics = {
        "pooled_bps": {
            "mean_net_return": 1e4 * _mean([r["net_return"] for r in rows]),
            "mean_gross_return": 1e4 * _mean([r["gross_return"] for r in rows]),
            "mean_fee_drag": 1e4 * _mean([r["fee_drag"] for r in rows]),
            "net_return_sd": 1e4 * (math.sqrt(sum((r["net_return"] - _mean([x["net_return"] for x in rows])) ** 2
                                                  for r in rows) / (len(rows) - 1)) if len(rows) > 1 else 0.0)},
        "mean_abs_spread_entry_bps": 1e4 * _mean([abs(r["spread_entry"]) for r in rows]),
        "mean_abs_spread_exit_bps": 1e4 * _mean([abs(r["spread_exit"]) for r in rows]),
        "mean_abs_spread_change_bps": 1e4 * _mean([abs(r["spread_exit"]) - abs(r["spread_entry"]) for r in rows]),
        "quarter_mean_net_bps": {q: 1e4 * _mean(v) for q, v in sorted(by_quarter.items())},
        "by_side_mean_net_bps": {side: (1e4 * _mean([r["net_return"] for r in rows if r["side"] == side])
                                        if any(r["side"] == side for r in rows) else None)
                                 for side in ("long-near", "long-far")},
        "by_base_mean_net_bps": {b: 1e4 * _mean([r["net_return"] for r in rows if r["base"] == b]) for b in bases},
    }
    passed = all(gate["pass"] for gate in gates.values())
    return {"window": label, "status": "pass-on-this-sample" if passed else "falsified",
            "sample_minimum_met": True, "sample": sample, "rejected": rejected,
            "gate_results": gates, "diagnostics": diagnostics,
            "positions": [{k: r[k] for k in ("base", "near", "far", "quarter", "entry_day", "exit_day",
                                             "hold_days", "side", "spread_entry", "spread_exit",
                                             "gross_return", "net_return")} for r in rows]}


def write_json_new(path: Path, payload: dict) -> None:
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1, sort_keys=True)
        handle.write("\n")


def validate_primary_artifact(primary: dict, contract_sha: str, pins: dict) -> None:
    """Bind a replication run to the exact primary artifact this runner writes.

    Reads pins from `static_inputs`, which is where main() nests them; checking the top level
    instead would abort every genuine hand-off.
    """
    if primary.get("role") != "primary":
        raise ValueError(f"replication artifact is not role primary: {primary.get('role')!r}")
    if primary.get("contract_sha256") != contract_sha:
        raise ValueError(f"contract_sha256 mismatch: expected {contract_sha}, "
                         f"observed {primary.get('contract_sha256')}")
    observed = primary.get("static_inputs")
    if not isinstance(observed, dict):
        raise ValueError("primary artifact carries no static_inputs block")
    for field in PIN_FIELDS:
        if observed.get(field) != pins[field]:
            raise ValueError(f"{field} mismatch: expected {pins[field]}, observed {observed.get(field)}")
    status = primary.get("result", {}).get("status")
    if status != "pass-on-this-sample":
        raise ValueError(f"replication requires a committed primary pass; observed {status!r}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--window", choices=("primary", "replication"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-inputs-only", action="store_true")
    parser.add_argument("--primary-artifact", type=Path)
    parser.add_argument("--primary-artifact-sha256")
    args = parser.parse_args(argv)

    observed_contract = sha256_file(CONTRACT_PATH)
    if observed_contract != args.contract_sha256:
        raise SystemExit(f"contract hash mismatch: expected {args.contract_sha256}, observed {observed_contract}")
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite existing result: {args.output}")
    _inside_root(args.output, WORKSPACE_ROOT)

    inventory, verified = verify_static_inputs(contract)
    runner_sha = sha256_file(Path(__file__).resolve())
    pins = {"input_inventory_sha256": contract["assets_and_data"]["input_inventory_sha256"],
            "prior_window_audit_sha256": contract["novelty_and_consumption"]["prior_window_audit_sha256"],
            "policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
            "runner_sha256": runner_sha}
    if args.verify_inputs_only:
        print(json.dumps({**{f: pins[f] for f in PIN_FIELDS}, "contract_sha256": observed_contract,
                          "archives_verified": verified}, indent=2, sort_keys=True))
        return 0

    primary_sha = None
    if args.window == "replication":
        if args.primary_artifact is None or args.primary_artifact_sha256 is None:
            raise SystemExit("replication requires --primary-artifact and --primary-artifact-sha256")
        primary_sha = sha256_file(args.primary_artifact)
        if primary_sha != args.primary_artifact_sha256:
            raise SystemExit("primary artifact hash mismatch")
        primary = json.loads(Path(args.primary_artifact).read_text(encoding="utf-8"))
        validate_primary_artifact(primary, observed_contract, pins)

    days = load_contract_days(inventory)
    result = run_window(contract, args.window, inventory, days)
    payload = {"schema_version": 1, "contract_id": contract["id"], "contract_sha256": observed_contract,
               "role": args.window, "static_inputs": {**pins, "archives_verified": verified},
               "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
               "result": result}
    if primary_sha is not None:
        payload["primary_artifact_sha256"] = primary_sha
    write_json_new(args.output, payload)
    print(json.dumps({f"{args.window}_status": result["status"],
                      "positions": result["sample"]["positions"],
                      "included_quarters": result["sample"]["included_quarters"],
                      "output": str(args.output)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
