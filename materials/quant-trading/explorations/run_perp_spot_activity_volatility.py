#!/usr/bin/env python3
"""Frozen test: perpetual-versus-spot activity and next-day volatility surprise.

Executes the contract at quant-trading/explorations/pre_statement.json (hash
pinned by the campaign manifest). Sample-first: the primary window resolves the
verdict; the preregistered replication window runs only on a primary pass.

Integrity rules enforced here, not by convention:
  * contract hash, input-inventory hash, and every archive hash are re-verified
    against upstream .CHECKSUM sidecars before any outcome value is computed;
  * hourly keys must be hour-aligned, unique, and present for every hour the
    frozen rule requires, else the day is rejected with a counted reason;
  * outputs are write-once; an existing path aborts the run.
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
CONTRACT_PATH = EXPLORATIONS / "pre_statement.json"
INVENTORY_PATH = EXPLORATIONS / "perp-spot-activity-inputs.json"
ASSETS = ("BTCUSDT", "ETHUSDT")
HOUR_MS = 3_600_000
SIGNAL_HOURS = 24
OUTCOME_HOURS = 24
QUINTILE_FRACTION = 0.20
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20260901
BOOTSTRAP_BLOCK_DATES = 7
RV_FLOOR = 1e-12


# ---------------------------------------------------------------- utilities ---

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iso_ms(ms: int) -> str:
    return datetime.fromtimestamp(ms / 1000, timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S"
    ) + "Z"


def parse_iso_ms(value: str) -> int:
    return int(
        datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp() * 1000
    )


def normalize_epoch(raw: str) -> int:
    """Binance archives are ms before 2025 and microseconds from 2025."""
    value = int(raw)
    return value // 1000 if value > 10_000_000_000_000 else value


def write_json_new(path: Path, payload: dict) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite frozen artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


# ------------------------------------------------------------ input loading ---

def month_keys(start_ms: int, end_ms: int) -> list[str]:
    """Calendar months intersecting [start-2h, end+2h), context included."""
    first = datetime.fromtimestamp((start_ms - 2 * HOUR_MS) / 1000, timezone.utc)
    last = datetime.fromtimestamp((end_ms + 2 * HOUR_MS) / 1000, timezone.utc)
    keys: list[str] = []
    year, month = first.year, first.month
    while (year, month) <= (last.year, last.month):
        keys.append(f"{year:04d}-{month:02d}")
        year, month = (year + 1, 1) if month == 12 else (year, month + 1)
    return keys


def verify_archive(rel: str) -> str:
    path = WORKSPACE_ROOT / rel
    sidecar = path.with_name(path.name + ".CHECKSUM")
    expected = sidecar.read_text(encoding="utf-8").split()[0].lower()
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"checksum mismatch for {rel}: {observed} != {expected}")
    return observed


def load_hours(asset: str, kind: str, months: list[str]) -> dict[int, dict]:
    """open_time -> {close, quote_volume}; duplicates and misalignment fail."""
    out: dict[int, dict] = {}
    for month in months:
        rel = (
            f"quant-trading/data/raw/binance/{kind}/{asset}/"
            f"{asset}-1h-{month}.zip"
        )
        path = WORKSPACE_ROOT / rel
        if not path.is_file():
            raise FileNotFoundError(f"required archive absent: {rel}")
        verify_archive(rel)
        with zipfile.ZipFile(path) as archive:
            members = archive.namelist()
            if len(members) != 1:
                raise ValueError(f"{rel}: expected exactly one member, got {members}")
            name = members[0]
            with archive.open(name) as handle:
                text = io.TextIOWrapper(handle, encoding="utf-8").read()
        for row in csv.reader(io.StringIO(text)):
            if not row or not row[0] or not row[0].lstrip("-").isdigit():
                continue  # header row in 2025+ archives
            open_ms = normalize_epoch(row[0])
            if open_ms % HOUR_MS != 0:
                raise ValueError(f"{rel}: non-hour-aligned open_time {open_ms}")
            if open_ms in out:
                raise ValueError(f"{rel}: duplicate open_time {open_ms}")
            close = float(row[4])
            quote_volume = float(row[7])
            if not math.isfinite(close) or close <= 0.0:
                raise ValueError(f"{rel}: nonpositive close at {open_ms}")
            if not math.isfinite(quote_volume) or quote_volume < 0.0:
                raise ValueError(f"{rel}: negative quote_volume at {open_ms}")
            out[open_ms] = {"close": close, "quote_volume": quote_volume}
    return out


# -------------------------------------------------------------- frozen rule ---

def daily_observation(asset: str, decision_ms: int, spot: dict[int, dict],
                      perp: dict[int, dict]) -> dict:
    """One frozen observation. Raises ValueError on any missing required hour."""
    signal_hours = [decision_ms - (k + 1) * HOUR_MS for k in range(SIGNAL_HOURS)]
    baseline_hours = [decision_ms - (k + 1) * HOUR_MS for k in range(SIGNAL_HOURS + 1)]
    outcome_hours = [decision_ms + k * HOUR_MS for k in range(OUTCOME_HOURS)]

    for hour in sorted(set(baseline_hours + outcome_hours)):
        if hour not in spot:
            raise ValueError(f"missing required spot hour {hour} for {asset}")
    for hour in sorted(set(signal_hours)):
        if hour not in perp:
            raise ValueError(f"missing required perp hour {hour} for {asset}")

    spot_volume = sum(spot[h]["quote_volume"] for h in signal_hours)
    perp_volume = sum(perp[h]["quote_volume"] for h in signal_hours)
    if spot_volume <= 0.0 or perp_volume <= 0.0:
        raise ValueError(f"nonpositive signal volume for {asset} at {decision_ms}")

    baseline_closes = [spot[h]["close"] for h in sorted(baseline_hours)]
    trailing_returns = [
        math.log(baseline_closes[i + 1] / baseline_closes[i])
        for i in range(len(baseline_closes) - 1)
    ]
    future_closes = [spot[decision_ms - HOUR_MS]["close"]] + [
        spot[h]["close"] for h in outcome_hours
    ]
    future_returns = [
        math.log(future_closes[i + 1] / future_closes[i])
        for i in range(len(future_closes) - 1)
    ]

    trailing_rv = math.sqrt(sum(r * r for r in trailing_returns))
    future_rv = math.sqrt(sum(r * r for r in future_returns))
    return {
        "asset": asset,
        "date": iso_ms(decision_ms)[:10],
        "decision_ms": decision_ms,
        "x": math.log(perp_volume / spot_volume),
        "y": math.log(max(future_rv, RV_FLOOR) / max(trailing_rv, RV_FLOOR)),
        "perp_quote_volume": perp_volume,
        "spot_quote_volume": spot_volume,
        "trailing_rv": trailing_rv,
        "future_rv": future_rv,
        "trailing_return_count": len(trailing_returns),
        "future_return_count": len(future_returns),
        "signal_first_ms": min(signal_hours),
        "signal_last_ms": max(signal_hours),
        "outcome_first_return_ms": min(outcome_hours),
        "outcome_last_return_ms": max(outcome_hours),
    }


# --------------------------------------------------------------- statistics ---

def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    dx = math.sqrt(sum((v - mx) ** 2 for v in xs))
    dy = math.sqrt(sum((v - my) ** 2 for v in ys))
    if dx == 0.0 or dy == 0.0:
        return 0.0
    return num / (dx * dy)


def spearman(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 3:
        raise ValueError("spearman needs paired inputs of length >= 3")
    return pearson(average_ranks(list(xs)), average_ranks(list(ys)))


def quintile_gap(rows: list[dict]) -> dict:
    ranked = sorted(rows, key=lambda r: r["x"])
    size = int(QUINTILE_FRACTION * len(ranked))
    if size < 2:
        raise ValueError("quintile tails too small")
    lower = [r["y"] for r in ranked[:size]]
    upper = [r["y"] for r in ranked[-size:]]
    return {
        "tail_size": size,
        "lower_median_y": median(lower),
        "upper_median_y": median(upper),
        "gap": median(upper) - median(lower),
    }


def block_bootstrap_p(rows: list[dict], replicates: int = BOOTSTRAP_REPLICATES,
                      seed: int = BOOTSTRAP_SEED,
                      block_dates: int = BOOTSTRAP_BLOCK_DATES) -> float:
    """One-sided circular moving-block bootstrap over dates (assets paired)."""
    by_date: dict[str, list[dict]] = {}
    for row in rows:
        by_date.setdefault(row["date"], []).append(row)
    dates = sorted(by_date)
    n_dates = len(dates)
    if n_dates < block_dates * 3:
        raise ValueError("too few dates for the frozen block bootstrap")
    rng = random.Random(seed)
    at_or_below_zero = 0
    for _ in range(replicates):
        picked: list[str] = []
        while len(picked) < n_dates:
            start = rng.randrange(n_dates)
            for offset in range(block_dates):
                picked.append(dates[(start + offset) % n_dates])
        picked = picked[:n_dates]  # exact original length; no n inflation
        xs: list[float] = []
        ys: list[float] = []
        for date in picked:
            for row in by_date[date]:
                xs.append(row["x"])
                ys.append(row["y"])
        if pearson(average_ranks(xs), average_ranks(ys)) <= 0.0:
            at_or_below_zero += 1
    return (1 + at_or_below_zero) / (replicates + 1)


# ------------------------------------------------------------ window driver ---

def decision_days(start_ms: int, end_ms: int) -> list[int]:
    out = []
    cursor = datetime.fromtimestamp(start_ms / 1000, timezone.utc)
    while int(cursor.timestamp() * 1000) < end_ms:
        out.append(int(cursor.timestamp() * 1000))
        cursor += timedelta(days=1)
    return out


def run_window(contract: dict, label: str) -> dict:
    window = contract["windows"][label]
    start_ms = parse_iso_ms(window["start_utc_inclusive"])
    end_ms = parse_iso_ms(window["end_utc_exclusive"])
    months = month_keys(start_ms, end_ms + OUTCOME_HOURS * HOUR_MS)
    gates = contract["gates"]
    minimum = contract["minimum_sample"][label]

    rows: list[dict] = []
    rejects: dict[str, dict[str, int]] = {a: {} for a in ASSETS}
    for asset in ASSETS:
        spot = load_hours(asset, "spot-1h", months)
        perp = load_hours(asset, "um-perp-1h", months)
        for decision_ms in decision_days(start_ms, end_ms):
            try:
                rows.append(daily_observation(asset, decision_ms, spot, perp))
            except ValueError as exc:
                key = str(exc).split(" for ")[0].split(" hour ")[0]
                rejects[asset][key] = rejects[asset].get(key, 0) + 1

    per_asset = {}
    for asset in ASSETS:
        subset = [r for r in rows if r["asset"] == asset]
        per_asset[asset] = {
            "observations": len(subset),
            "spearman_rho": spearman([r["x"] for r in subset],
                                     [r["y"] for r in subset]) if len(subset) >= 3 else None,
            "median_x": median([r["x"] for r in subset]) if subset else None,
            "median_y": median([r["y"] for r in subset]) if subset else None,
            "rejected": rejects[asset],
        }

    pooled_n = len(rows)
    sample_met = (
        pooled_n >= minimum["pooled"]
        and all(per_asset[a]["observations"] >= minimum["per_asset"] for a in ASSETS)
    )
    if not sample_met:
        return {
            "window": window,
            "status": "inconclusive",
            "reason": "frozen minimum sample not met",
            "pooled_observations": pooled_n,
            "minimum_sample": minimum,
            "per_asset": per_asset,
            "gate_results": None,
        }

    pooled_rho = spearman([r["x"] for r in rows], [r["y"] for r in rows])
    tails = quintile_gap(rows)
    p_value = block_bootstrap_p(rows)
    threshold_g1 = 0.08
    threshold_g3 = math.log(1.10)
    gate_results = {
        "G1_pooled_effect": {"rho": pooled_rho, "threshold": threshold_g1,
                             "pass": pooled_rho >= threshold_g1},
        "G2_asset_consistency": {
            "per_asset_rho": {a: per_asset[a]["spearman_rho"] for a in ASSETS},
            "pass": all((per_asset[a]["spearman_rho"] or 0.0) > 0.0 for a in ASSETS),
        },
        "G3_economic_magnitude": {**tails, "threshold": threshold_g3,
                                  "pass": tails["gap"] >= threshold_g3},
        "G4_dependence_robustness": {
            "one_sided_block_bootstrap_p": p_value,
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "block_dates": BOOTSTRAP_BLOCK_DATES,
            "pass": p_value < 0.05,
        },
    }
    all_pass = all(g["pass"] for g in gate_results.values())
    return {
        "window": window,
        "status": "pass-on-this-sample" if all_pass else "falsified",
        "reason": None if all_pass else "at least one frozen gate is false",
        "pooled_observations": pooled_n,
        "minimum_sample": minimum,
        "pooled_spearman_rho": pooled_rho,
        "per_asset": per_asset,
        "gate_results": gate_results,
        "gate_descriptions": gates,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT_PATH)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--window", choices=("primary", "replication"),
                        required=True,
                        help="replication requires a committed primary pass")
    parser.add_argument("--output", type=Path, required=True,
                        help="write-once artifact for the selected window")
    parser.add_argument("--primary-artifact", type=Path,
                        help="committed primary artifact; required for replication")
    parser.add_argument("--verify-inputs-only", action="store_true")
    args = parser.parse_args(argv)

    contract_path = args.contract.resolve()
    observed_contract = sha256_file(contract_path)
    if observed_contract != args.contract_sha256:
        raise SystemExit(
            f"contract hash mismatch: expected {args.contract_sha256}, "
            f"observed {observed_contract}"
        )
    contract = json.loads(contract_path.read_text(encoding="utf-8"))

    inventory_sha = sha256_file(INVENTORY_PATH)
    if inventory_sha != contract["assets_and_data"]["input_inventory_sha256"]:
        raise SystemExit("input inventory hash mismatch against the contract")
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    verified = []
    for entry in inventory["files"]:
        observed = verify_archive(entry["path"])
        if observed != entry["sha256"]:
            raise SystemExit(f"archive drifted from inventory: {entry['path']}")
        verified.append(entry["path"])

    audit_sha = sha256_file(EXPLORATIONS / "prior-window-consumption-audit.json")
    if audit_sha != contract["novelty_and_consumption"]["prior_window_audit_sha256"]:
        raise SystemExit("prior-window audit hash mismatch against the contract")

    if args.verify_inputs_only:
        print(json.dumps({
            "contract_sha256": observed_contract,
            "input_inventory_sha256": inventory_sha,
            "prior_window_audit_sha256": audit_sha,
            "archives_verified": len(verified),
        }, indent=2))
        return 0

    primary_status = None
    if args.window == "replication":
        if args.primary_artifact is None:
            raise SystemExit("--primary-artifact is required for replication")
        if not args.primary_artifact.is_file():
            raise SystemExit(
                "replication refused: committed primary artifact absent at "
                f"{args.primary_artifact}")
        primary_doc = json.loads(args.primary_artifact.read_text(encoding="utf-8"))
        if primary_doc.get("role") != "primary":
            raise SystemExit("replication refused: artifact is not the primary role")
        if primary_doc.get("contract_sha256") != observed_contract:
            raise SystemExit("replication refused: primary ran a different contract")
        primary_status = primary_doc.get("result", {}).get("status")
        if primary_status != "pass-on-this-sample":
            raise SystemExit(
                "replication refused: the contract permits replication only after "
                f"a committed primary pass; observed {primary_status!r}")

    result = run_window(contract, args.window)
    payload = {
        "schema_version": 1,
        "id": contract["id"],
        "created_at_utc": iso_ms(int(datetime.now(timezone.utc).timestamp() * 1000)),
        "contract_sha256": observed_contract,
        "input_inventory_sha256": inventory_sha,
        "prior_window_audit_sha256": audit_sha,
        "archives_verified": len(verified),
        "role": args.window,
        "result": result,
        "claim_boundary": contract["claim_boundary"],
    }
    if args.window == "replication":
        payload["primary_status"] = primary_status
        payload["primary_artifact"] = str(
            args.primary_artifact.resolve().relative_to(WORKSPACE_ROOT))
    write_json_new(args.output, payload)
    print(json.dumps({f"{args.window}_status": result["status"],
                      "pooled_observations": result["pooled_observations"],
                      "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
