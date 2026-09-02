#!/usr/bin/env python3
"""Frozen two-window test of non-hour quarter-boundary market activity.

The runner verifies every preregistration and raw-input hash before parsing a
candidate value, applies one fixed UTC clock contrast, and writes each result
with exclusive-create semantics. It is descriptive microstructure research,
not a return strategy.
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
from datetime import datetime, timedelta, timezone
from pathlib import Path
from statistics import median

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
EXPLORATIONS = WORKSPACE_ROOT / "quant-trading/explorations"
CONTRACT_PATH = EXPLORATIONS / "pre_statement_quarter_hour.json"
INVENTORY_PATH = EXPLORATIONS / "quarter-hour-inputs.json"
AUDIT_PATH = EXPLORATIONS / "quarter-hour-prior-window-audit.json"
POLICY_PATH = EXPLORATIONS / "frontier-loop-policy.json"

ASSETS = ("BTCUSDT", "ETHUSDT")
QUARTER_MINUTES = {15, 30, 45}
CONTROL_MINUTES = {5, 10, 20, 25, 35, 40, 50, 55}
MINUTE_MS = 60_000
DAY_MS = 86_400_000
MICROSECOND_THRESHOLD = 100_000_000_000_000
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 20_260_901
BOOTSTRAP_BLOCK_DATES = 7
EXPECTED_HEADER = (
    "open_time", "open", "high", "low", "close", "volume", "close_time",
    "quote_volume", "count", "taker_buy_volume", "taker_buy_quote_volume",
    "ignore",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_iso_ms(value: str) -> int:
    return int(datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
               * 1000)


def iso_ms(ms: int) -> str:
    return datetime.fromtimestamp(ms / 1000, timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S") + "Z"


def normalize_epoch(raw: str) -> int:
    value = int(raw)
    return value // 1000 if value > MICROSECOND_THRESHOLD else value


def write_json_new(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def required_offsets() -> set[int]:
    return {
        (hour * 60 + minute) * MINUTE_MS
        for hour in range(24)
        for minute in QUARTER_MINUTES | CONTROL_MINUTES
    }


def _inside_root(path: Path, root: Path) -> Path:
    resolved_root = root.resolve()
    resolved = path.resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"input escapes workspace root: {path}") from exc
    return resolved


def verify_archive(entry: dict, root: Path = WORKSPACE_ROOT) -> str:
    archive_path = _inside_root(root / entry["path"], root)
    sidecar_path = _inside_root(root / entry["checksum_sidecar"], root)
    if not archive_path.is_file() or not sidecar_path.is_file():
        raise ValueError(f"archive or sidecar absent: {entry['path']}")
    if archive_path.stat().st_size != entry["size_bytes"]:
        raise ValueError(f"archive byte size drift: {entry['path']}")
    if sha256_file(sidecar_path) != entry["checksum_sidecar_sha256"]:
        raise ValueError(f"sidecar file hash drift: {entry['checksum_sidecar']}")
    sidecar_fields = sidecar_path.read_text(encoding="utf-8").split()
    if len(sidecar_fields) < 2:
        raise ValueError(f"malformed upstream sidecar: {entry['checksum_sidecar']}")
    if sidecar_fields[0].lower() != entry["sha256"]:
        raise ValueError(f"upstream sidecar hash mismatch: {entry['path']}")
    if Path(sidecar_fields[-1]).name != archive_path.name:
        raise ValueError(f"upstream sidecar filename mismatch: {entry['path']}")
    if sha256_file(archive_path) != entry["sha256"]:
        raise ValueError(f"archive hash drift: {entry['path']}")

    with zipfile.ZipFile(archive_path) as archive:
        members = archive.namelist()
        if members != [entry["zip_member"]]:
            raise ValueError(f"zip member mismatch: {entry['path']}")
        corrupt_member = archive.testzip()
        if corrupt_member is not None:
            raise ValueError(f"zip CRC failure: {entry['path']}:{corrupt_member}")
        with archive.open(members[0]) as raw:
            first = next(csv.reader(io.TextIOWrapper(raw, encoding="utf-8")))
    has_header = tuple(first) == EXPECTED_HEADER
    if has_header != entry["header_present"] or not has_header:
        raise ValueError(f"zip header form mismatch: {entry['path']}")
    return members[0]


def verify_static_inputs(contract: dict) -> tuple[dict, dict]:
    expected = {
        "input_inventory_sha256": contract["assets_and_data"]
        ["input_inventory_sha256"],
        "prior_window_audit_sha256": contract["novelty_and_consumption"]
        ["prior_window_audit_sha256"],
        "loop_policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
    }
    observed = {
        "input_inventory_sha256": sha256_file(INVENTORY_PATH),
        "prior_window_audit_sha256": sha256_file(AUDIT_PATH),
        "loop_policy_sha256": sha256_file(POLICY_PATH),
    }
    for field, value in expected.items():
        if observed[field] != value:
            raise ValueError(
                f"{field} mismatch: expected {value}, observed {observed[field]}")

    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("candidate_outcomes_accessed") is not False:
        raise ValueError("input inventory does not certify pre-outcome construction")
    if tuple(inventory.get("assets", ())) != ASSETS:
        raise ValueError("input inventory asset order mismatch")
    entries = inventory.get("files", [])
    if inventory.get("archive_count") != len(entries) or not entries:
        raise ValueError("input inventory archive count mismatch")
    paths = [entry["path"] for entry in entries]
    if len(paths) != len(set(paths)):
        raise ValueError("input inventory contains duplicate archive paths")
    for entry in entries:
        verify_archive(entry)
    return inventory, {**observed, "archives_verified": len(entries)}


def window_months(start_ms: int, end_ms: int) -> list[str]:
    cursor = datetime.fromtimestamp(start_ms / 1000, timezone.utc).replace(
        day=1, hour=0, minute=0, second=0, microsecond=0)
    out: list[str] = []
    while int(cursor.timestamp() * 1000) < end_ms:
        out.append(cursor.strftime("%Y-%m"))
        cursor = (cursor.replace(day=28) + timedelta(days=4)).replace(day=1)
    return out


def _entry_asset_month(entry: dict) -> tuple[str, str]:
    name = Path(entry["path"]).name
    try:
        asset, suffix = name.split("-1m-", 1)
    except ValueError as exc:
        raise ValueError(f"unexpected inventory filename: {name}") from exc
    if not suffix.endswith(".zip"):
        raise ValueError(f"unexpected inventory filename: {name}")
    return asset, suffix[:-4]


def _parse_archive(entry: dict, start_ms: int, end_ms: int,
                   seen: set[int], days: dict[int, dict[int, dict]]) -> None:
    path = WORKSPACE_ROOT / entry["path"]
    with zipfile.ZipFile(path) as archive:
        with archive.open(entry["zip_member"]) as raw:
            reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8"))
            header = tuple(next(reader))
            if header != EXPECTED_HEADER:
                raise ValueError(f"schema header drift: {entry['path']}")
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(EXPECTED_HEADER):
                    raise ValueError(
                        f"schema column count at {entry['path']}:{line_number}")
                try:
                    open_ms = normalize_epoch(row[0])
                    open_price = float(row[1])
                    close_price = float(row[4])
                    trades = int(row[8])
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"schema value at {entry['path']}:{line_number}") from exc
                if not start_ms <= open_ms < end_ms:
                    continue
                if open_ms % MINUTE_MS != 0:
                    raise ValueError(
                        f"non-minute timestamp at {entry['path']}:{line_number}")
                if open_ms in seen:
                    raise ValueError(f"duplicate timestamp {open_ms} for input asset")
                seen.add(open_ms)
                if (not math.isfinite(open_price) or open_price <= 0.0
                        or not math.isfinite(close_price) or close_price <= 0.0):
                    raise ValueError(
                        f"nonpositive/nonfinite price at {entry['path']}:{line_number}")
                if trades < 0:
                    raise ValueError(
                        f"negative trade count at {entry['path']}:{line_number}")
                minute = (open_ms // MINUTE_MS) % 60
                if minute not in QUARTER_MINUTES | CONTROL_MINUTES:
                    continue
                day_ms = open_ms - open_ms % DAY_MS
                days.setdefault(day_ms, {})[open_ms] = {
                    "open": open_price,
                    "close": close_price,
                    "number_of_trades": trades,
                }


def load_asset_days(inventory: dict, asset: str, start_ms: int,
                    end_ms: int) -> dict[int, dict[int, dict]]:
    wanted_months = window_months(start_ms, end_ms)
    entries_by_month: dict[str, dict] = {}
    for entry in inventory["files"]:
        entry_asset, month = _entry_asset_month(entry)
        if entry_asset == asset and month in wanted_months:
            if month in entries_by_month:
                raise ValueError(f"duplicate inventory month for {asset} {month}")
            entries_by_month[month] = entry
    if sorted(entries_by_month) != sorted(wanted_months):
        raise ValueError(
            f"inventory months for {asset}: expected {wanted_months}, "
            f"observed {sorted(entries_by_month)}")

    days: dict[int, dict[int, dict]] = {}
    seen: set[int] = set()
    for month in wanted_months:
        _parse_archive(entries_by_month[month], start_ms, end_ms, seen, days)
    return days


def daily_observation(asset: str, day_ms: int,
                      bars: dict[int, dict]) -> dict:
    offsets = required_offsets()
    required = {day_ms + offset for offset in offsets}
    if not required.issubset(bars):
        missing = len(required - set(bars))
        raise ValueError(f"required bars missing ({missing})")

    quarter = [bars[day_ms + offset] for offset in sorted(offsets)
               if ((offset // MINUTE_MS) % 60) in QUARTER_MINUTES]
    control = [bars[day_ms + offset] for offset in sorted(offsets)
               if ((offset // MINUTE_MS) % 60) in CONTROL_MINUTES]
    if len(quarter) != 72 or len(control) != 192:
        raise ValueError("required bars have wrong clock-group cardinality")
    if any(bar["number_of_trades"] == 0 for bar in quarter + control):
        raise ValueError("required bar has zero trades")

    q_trades = sum(bar["number_of_trades"] for bar in quarter) / len(quarter)
    c_trades = sum(bar["number_of_trades"] for bar in control) / len(control)
    q_absret = sum(abs(math.log(bar["close"] / bar["open"]))
                   for bar in quarter) / len(quarter)
    c_absret = sum(abs(math.log(bar["close"] / bar["open"]))
                   for bar in control) / len(control)
    if min(q_trades, c_trades, q_absret, c_absret) <= 0.0:
        raise ValueError("nonpositive daily ratio component")

    return {
        "asset": asset,
        "date": datetime.fromtimestamp(day_ms / 1000, timezone.utc).date().isoformat(),
        "quarter_bars": len(quarter),
        "control_bars": len(control),
        "x_count": math.log(q_trades / c_trades),
        "x_absret": math.log(q_absret / c_absret),
    }


def _draw_with_rng(dates: list[str], rng: random.Random,
                   block_dates: int) -> list[str]:
    picked: list[str] = []
    while len(picked) < len(dates):
        start = rng.randrange(len(dates))
        for offset in range(block_dates):
            picked.append(dates[(start + offset) % len(dates)])
    return picked[:len(dates)]


def draw_block_dates(dates: list[str], seed: int = BOOTSTRAP_SEED,
                     block_dates: int = BOOTSTRAP_BLOCK_DATES) -> list[str]:
    if not dates:
        raise ValueError("cannot bootstrap zero dates")
    return _draw_with_rng(list(dates), random.Random(seed), block_dates)


def block_bootstrap_p(rows: list[dict], field: str,
                      replicates: int = BOOTSTRAP_REPLICATES,
                      seed: int = BOOTSTRAP_SEED,
                      block_dates: int = BOOTSTRAP_BLOCK_DATES) -> float:
    by_date: dict[str, list[float]] = {}
    for row in rows:
        by_date.setdefault(row["date"], []).append(row[field])
    dates = sorted(by_date)
    if len(dates) < block_dates * 3:
        raise ValueError("too few dates for frozen block bootstrap")
    if any(len(by_date[date]) != len(ASSETS) for date in dates):
        raise ValueError("bootstrap rows are not asset-paired by date")

    rng = random.Random(seed)
    at_or_below_zero = 0
    for _ in range(replicates):
        picked = _draw_with_rng(dates, rng, block_dates)
        values = [value for date in picked for value in by_date[date]]
        if median(values) <= 0.0:
            at_or_below_zero += 1
    return (1 + at_or_below_zero) / (replicates + 1)


def _month_stability(rows: list[dict], field: str) -> dict:
    by_month: dict[str, list[float]] = {}
    for row in rows:
        by_month.setdefault(row["date"][:7], []).append(row[field])
    medians = {month: median(values) for month, values in sorted(by_month.items())}
    positive = sum(value > 0.0 for value in medians.values())
    return {
        "included_months": len(medians),
        "positive_months": positive,
        "positive_fraction": positive / len(medians) if medians else None,
        "monthly_pooled_medians": medians,
    }


def _rejection_key(exc: ValueError) -> str:
    message = str(exc)
    if message.startswith("required bars missing"):
        return "required-bars-missing"
    if message == "required bar has zero trades":
        return "zero-trade-required-bar"
    if message == "nonpositive daily ratio component":
        return "nonpositive-daily-ratio-component"
    raise exc


def run_window(contract: dict, label: str, inventory: dict) -> dict:
    window = contract["windows"][label]
    minimum = contract["minimum_sample"][label]
    start_ms = parse_iso_ms(window["start_utc_inclusive"])
    end_ms = parse_iso_ms(window["end_utc_exclusive"])
    expected_days = list(range(start_ms, end_ms, DAY_MS))

    complete: dict[str, dict[int, dict]] = {}
    rejected: dict[str, dict[str, int]] = {asset: {} for asset in ASSETS}
    for asset in ASSETS:
        raw_days = load_asset_days(inventory, asset, start_ms, end_ms)
        asset_complete: dict[int, dict] = {}
        for day_ms in expected_days:
            try:
                asset_complete[day_ms] = daily_observation(
                    asset, day_ms, raw_days.get(day_ms, {}))
            except ValueError as exc:
                key = _rejection_key(exc)
                rejected[asset][key] = rejected[asset].get(key, 0) + 1
        complete[asset] = asset_complete

    paired_days = sorted(set(complete[ASSETS[0]]) & set(complete[ASSETS[1]]))
    rows = [complete[asset][day_ms]
            for day_ms in paired_days for asset in ASSETS]
    per_asset_n = {
        asset: sum(row["asset"] == asset for row in rows) for asset in ASSETS
    }
    included_months = len({row["date"][:7] for row in rows})
    sample = {
        "expected_dates": len(expected_days),
        "paired_dates": len(paired_days),
        "pooled_observations": len(rows),
        "per_asset_observations": per_asset_n,
        "included_calendar_months": included_months,
        "rejected_asset_dates": rejected,
        "unpaired_complete_dates": {
            asset: len(set(complete[asset]) - set(paired_days)) for asset in ASSETS
        },
        "minimum": minimum,
    }
    sample_ok = (
        len(paired_days) >= minimum["paired_dates"]
        and len(rows) >= minimum["pooled"]
        and included_months >= minimum["included_calendar_months"]
        and all(per_asset_n[asset] >= minimum["per_asset"] for asset in ASSETS)
    )
    if not sample_ok:
        return {
            "window": window,
            "status": "inconclusive",
            "reason": "frozen minimum sample not met",
            "sample": sample,
            "gate_results": None,
        }

    pooled = {
        "median_x_count": median(row["x_count"] for row in rows),
        "median_x_absret": median(row["x_absret"] for row in rows),
    }
    per_asset = {}
    for asset in ASSETS:
        subset = [row for row in rows if row["asset"] == asset]
        per_asset[asset] = {
            "observations": len(subset),
            "median_x_count": median(row["x_count"] for row in subset),
            "median_x_absret": median(row["x_absret"] for row in subset),
        }
    count_stability = _month_stability(rows, "x_count")
    absret_stability = _month_stability(rows, "x_absret")
    # Each call starts the same frozen RNG stream. Date count and block length
    # are identical, so both endpoints receive exactly the same date draws.
    p_count = block_bootstrap_p(rows, "x_count")
    p_absret = block_bootstrap_p(rows, "x_absret")
    alpha = contract["loop_multiplicity"]["alpha_k"]

    gate_results = {
        "G1_count_magnitude": {
            "pooled_median_x_count": pooled["median_x_count"],
            "threshold": math.log(1.10),
            "pass": pooled["median_x_count"] >= math.log(1.10),
        },
        "G2_count_asset_consistency": {
            "per_asset_median_x_count": {
                asset: per_asset[asset]["median_x_count"] for asset in ASSETS
            },
            "threshold": math.log(1.08),
            "pass": all(per_asset[asset]["median_x_count"] >= math.log(1.08)
                        for asset in ASSETS),
        },
        "G3_absret_magnitude": {
            "pooled_median_x_absret": pooled["median_x_absret"],
            "threshold": math.log(1.05),
            "pass": pooled["median_x_absret"] >= math.log(1.05),
        },
        "G4_absret_asset_consistency": {
            "per_asset_median_x_absret": {
                asset: per_asset[asset]["median_x_absret"] for asset in ASSETS
            },
            "threshold": 0.0,
            "pass": all(per_asset[asset]["median_x_absret"] > 0.0
                        for asset in ASSETS),
        },
        "G5_month_stability": {
            "x_count": count_stability,
            "x_absret": absret_stability,
            "threshold": 0.75,
            "pass": (count_stability["positive_fraction"] >= 0.75
                     and absret_stability["positive_fraction"] >= 0.75),
        },
        "G6_dependence_robustness": {
            "one_sided_block_bootstrap_p_count": p_count,
            "one_sided_block_bootstrap_p_absret": p_absret,
            "alpha": alpha,
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "block_dates": BOOTSTRAP_BLOCK_DATES,
            "pass": p_count < alpha and p_absret < alpha,
        },
    }
    all_pass = all(gate["pass"] for gate in gate_results.values())
    return {
        "window": window,
        "status": "pass-on-this-sample" if all_pass else "falsified",
        "reason": None if all_pass else "at least one frozen gate is false",
        "sample": sample,
        "pooled": pooled,
        "per_asset": per_asset,
        "gate_results": gate_results,
        "gate_descriptions": contract["gates"],
    }


def validate_primary_artifact(primary: dict, expected: dict) -> None:
    if primary.get("role") != "primary":
        raise ValueError("replication artifact is not role primary")
    for field, value in expected.items():
        if primary.get(field) != value:
            raise ValueError(
                f"replication primary {field} mismatch: expected {value}, "
                f"observed {primary.get(field)}")
    status = primary.get("result", {}).get("status")
    if status != "pass-on-this-sample":
        raise ValueError(
            "replication requires a committed primary pass; "
            f"observed {status!r}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT_PATH)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--window", choices=("primary", "replication"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-artifact", type=Path)
    parser.add_argument("--primary-artifact-sha256")
    parser.add_argument("--verify-inputs-only", action="store_true")
    args = parser.parse_args(argv)

    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    observed_contract = sha256_file(args.contract.resolve())
    if observed_contract != args.contract_sha256:
        raise SystemExit(
            f"contract hash mismatch: expected {args.contract_sha256}, "
            f"observed {observed_contract}")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    try:
        inventory, integrity = verify_static_inputs(contract)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"static integrity failure: {exc}") from exc

    runner_sha = sha256_file(Path(__file__).resolve())
    expected_primary = {
        "contract_sha256": observed_contract,
        "input_inventory_sha256": integrity["input_inventory_sha256"],
        "prior_window_audit_sha256": integrity["prior_window_audit_sha256"],
        "loop_policy_sha256": integrity["loop_policy_sha256"],
        "runner_sha256": runner_sha,
    }
    if args.verify_inputs_only:
        print(json.dumps({**expected_primary,
                          "archives_verified": integrity["archives_verified"]},
                         indent=2, sort_keys=True))
        return 0

    primary_sha = None
    if args.window == "replication":
        if args.primary_artifact is None or args.primary_artifact_sha256 is None:
            raise SystemExit(
                "replication requires --primary-artifact and "
                "--primary-artifact-sha256")
        if not args.primary_artifact.is_file():
            raise SystemExit(
                f"replication primary artifact absent: {args.primary_artifact}")
        primary_sha = sha256_file(args.primary_artifact)
        if primary_sha != args.primary_artifact_sha256:
            raise SystemExit(
                "replication primary artifact hash does not match committed hash")
        primary_doc = json.loads(
            args.primary_artifact.read_text(encoding="utf-8"))
        try:
            validate_primary_artifact(primary_doc, expected_primary)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

    try:
        result = run_window(contract, args.window, inventory)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"outcome input/schema failure: {exc}") from exc
    now_ms = int(datetime.now(timezone.utc).timestamp() * 1000)
    payload = {
        "schema_version": 1,
        "id": contract["id"],
        "created_at_utc": iso_ms(now_ms),
        **expected_primary,
        "archives_verified": integrity["archives_verified"],
        "role": args.window,
        "result": result,
        "claim_boundary": contract["claim_boundary"],
    }
    if args.window == "replication":
        payload["primary_artifact"] = str(
            args.primary_artifact.resolve().relative_to(WORKSPACE_ROOT))
        payload["primary_artifact_sha256"] = primary_sha
    write_json_new(args.output, payload)
    print(json.dumps({
        f"{args.window}_status": result["status"],
        "paired_dates": result["sample"]["paired_dates"],
        "pooled_observations": result["sample"]["pooled_observations"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
