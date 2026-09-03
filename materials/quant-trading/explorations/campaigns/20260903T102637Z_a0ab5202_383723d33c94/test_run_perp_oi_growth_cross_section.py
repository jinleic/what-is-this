#!/usr/bin/env python3
"""Behavioral tests for the frozen cross-sectional open-interest-growth runner."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_perp_cross_section_daily_reversal as parent
import run_perp_oi_growth_cross_section as runner

DAY_MS = runner.DAY_MS
HOUR_MS = runner.HOUR_MS
FEE = runner.TAKER_FEE
FULL = parent.FULL_DAY_MASK
START_2024 = runner.parse_iso_ms("2024-01-01T00:00:00Z")
END_2025 = runner.parse_iso_ms("2025-01-01T00:00:00Z")
START_2025 = END_2025
END_2026 = runner.parse_iso_ms("2026-01-01T00:00:00Z")


def day_record(open0: float, close: float, quote: float = 1e7, open1: float | None = None,
               complete: bool = True) -> dict:
    return {"hours": FULL if complete else FULL - 1, "duplicate": False, "bad_price": False,
            "open0": open0, "open1": open1 if open1 is not None else open0, "close23": close,
            "last_hour": 23, "last_close": close, "quote_volume": quote, "trades": 100}


def flat_records(symbols: list[str], first_day: int, days: int, price: float = 100.0) -> dict:
    records = {}
    for symbol in symbols:
        for offset in range(-runner.ELIGIBILITY_DAYS, days):
            records[(symbol, first_day + offset * DAY_MS)] = day_record(price, price)
    return records


def flat_open_interest(symbols: list[str], first_day: int, days: int) -> dict:
    return {(symbol, first_day + offset * DAY_MS): 1000.0
            for symbol in symbols for offset in range(-2, days + 1)}


def synthetic_contract(**overrides) -> dict:
    contract = {
        "id": "perp-oi-growth-cross-section-test",
        "frozen_rule": {"parameters": dict(runner.FROZEN_PARAMETERS)},
        "assets_and_data": {"input_inventory_sha256": "i" * 64, "kline_inventory_sha256": "k" * 64},
        "novelty_and_consumption": {"prior_window_audit_sha256": "a" * 64},
        "loop_multiplicity": {"policy_sha256": "p" * 64},
        "governance": {"parent_runner_sha256": "r" * 64},
        "windows": {"primary": {"start_utc_inclusive": "2024-01-01T00:00:00Z", "end_utc_exclusive": "2025-01-01T00:00:00Z"},
                    "replication": {"start_utc_inclusive": "2025-01-01T00:00:00Z", "end_utc_exclusive": "2026-01-01T00:00:00Z"}},
        "claim_boundary": "test boundary",
    }
    contract.update(overrides)
    return contract


def synthetic_inventory(universe_primary: list[str], universe_replication: list[str], files: list[dict]) -> dict:
    return {"candidate_outcomes_accessed": False, "kline_inventory_sha256": "k" * 64,
            "windows": {"primary": {"start_utc_inclusive": "2024-01-01T00:00:00Z", "end_utc_exclusive": "2025-01-01T00:00:00Z",
                                    "universe": universe_primary},
                        "replication": {"start_utc_inclusive": "2025-01-01T00:00:00Z", "end_utc_exclusive": "2026-01-01T00:00:00Z",
                                        "universe": universe_replication}},
            "archive_count": len(files), "files": files}


def write_metrics_zip(root: Path, symbol: str, date: str, rows: list[tuple[str, str, str]],
                      header: tuple[str, ...] = runner.METRICS_HEADER) -> dict:
    """rows: (create_time, symbol, sum_open_interest); other columns filled with 1."""
    body = ",".join(header) + "\n" + "".join(
        f"{stamp},{sym},{oi},1,1,1,1,1\n" for stamp, sym, oi in rows)
    member = f"{symbol}-metrics-{date}.csv"
    rel = Path("quant-trading/data/raw/binance/um-perp-metrics-5m") / symbol / f"{symbol}-metrics-{date}.zip"
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(member, body)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    sidecar = path.with_name(path.name + ".CHECKSUM")
    sidecar.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return {"symbol": symbol, "date": date, "path": str(rel), "checksum_sidecar": str(rel) + ".CHECKSUM",
            "size_bytes": path.stat().st_size, "sha256": digest,
            "checksum_sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(), "zip_member": member}


class ContractBindingTests(unittest.TestCase):
    def test_contract_parameters_must_match_runner_constants(self) -> None:
        contract = synthetic_contract()
        runner.check_contract_parameters(contract)
        contract["frozen_rule"]["parameters"]["universe_size"] = 73
        with self.assertRaisesRegex(ValueError, "differ from runner constants"):
            runner.check_contract_parameters(contract)

    def test_frozen_contract_binds_runner_and_windows(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        runner.check_contract_parameters(contract)
        self.assertEqual(contract["loop_multiplicity"]["alpha_k"], runner.LOOP_ALPHA)
        for label in ("primary", "replication"):
            start = runner.parse_iso_ms(contract["windows"][label]["start_utc_inclusive"])
            runner.universe_lookback(start)
        self.assertEqual(contract["governance"]["parent_runner_sha256"], runner.sha256_file(runner.PARENT_PATH))

    def test_static_input_hash_drift_aborts_before_any_archive_check(self) -> None:
        contract = synthetic_contract()
        good = {runner.INVENTORY_PATH: "i" * 64, runner.KLINE_INVENTORY_PATH: "k" * 64,
                runner.AUDIT_PATH: "a" * 64, runner.POLICY_PATH: "p" * 64, runner.PARENT_PATH: "r" * 64}
        for path in good:
            drifted = dict(good)
            drifted[path] = "0" * 64
            with mock.patch.object(runner, "sha256_file", side_effect=lambda p, d=drifted: d[Path(p)]), \
                    mock.patch.object(runner, "verify_metrics_archive") as verify:
                with self.assertRaisesRegex(ValueError, "mismatch"):
                    runner.verify_static_inputs(contract)
            verify.assert_not_called()

    def test_inventory_must_certify_pre_outcome_state_and_universes(self) -> None:
        contract = synthetic_contract()
        primary = [f"S{i:02d}USDT" for i in range(runner.UNIVERSE_SIZE)]
        replication = [f"R{i:02d}USDT" for i in range(runner.UNIVERSE_SIZE)]
        files = [{"symbol": "S00USDT", "date": "2024-01-01", "path": "x/S00USDT-metrics-2024-01-01.zip"}]
        pins = {runner.INVENTORY_PATH: "i" * 64, runner.KLINE_INVENTORY_PATH: "k" * 64,
                runner.AUDIT_PATH: "a" * 64, runner.POLICY_PATH: "p" * 64, runner.PARENT_PATH: "r" * 64}
        mutations = [
            (lambda inv: inv.update(candidate_outcomes_accessed=True), "pre-outcome"),
            (lambda inv: inv.update(kline_inventory_sha256="0" * 64), "bind the frozen kline"),
            (lambda inv: inv["windows"]["primary"].update(universe=primary[:-1]), "distinct symbols"),
            (lambda inv: inv["windows"].pop("replication"), "exactly the primary and replication"),
            (lambda inv: inv.update(archive_count=5), "archive count"),
            (lambda inv: inv["files"].append(dict(files[0])), "duplicate paths"),
            (lambda inv: inv["files"].append({"symbol": "ZZZUSDT", "date": "2024-01-01", "path": "x/ZZZUSDT-metrics-2024-01-01.zip"}), "outside both universes"),
        ]
        for mutate, message in mutations:
            inventory = synthetic_inventory(primary, replication, [dict(files[0])])
            mutate(inventory)
            if "archive count" not in message:
                inventory["archive_count"] = len(inventory["files"])
            with mock.patch.object(runner, "sha256_file", side_effect=lambda p: pins[Path(p)]), \
                    mock.patch.object(runner.INVENTORY_PATH.__class__, "read_text", lambda self, encoding=None: json.dumps(inventory)), \
                    mock.patch.object(runner, "verify_metrics_archive"), \
                    mock.patch.object(parent, "verify_archive"):
                with self.assertRaisesRegex(ValueError, message):
                    runner.verify_static_inputs(contract)


class MetricsParsingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = Path(tempfile.mkdtemp(prefix=".tmp-oi5-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)

    def parse(self, rows, symbol="ABCUSDT", date="2024-03-05", **kw):
        entry = write_metrics_zip(self.root, symbol, date, rows, **kw)
        return runner._parse_metrics_archive(entry, root=self.root), entry

    def test_last_row_at_or_after_cutoff_is_end_of_day(self) -> None:
        value, _ = self.parse([("2024-03-05 00:00:00", "ABCUSDT", "5"), ("2024-03-05 23:45:00", "ABCUSDT", "7"),
                               ("2024-03-05 23:55:00", "ABCUSDT", "9"), ("2024-03-05 23:50:00", "ABCUSDT", "8")])
        self.assertEqual(value, 9.0)

    def test_archive_without_row_at_cutoff_is_missing_not_imputed(self) -> None:
        value, _ = self.parse([("2024-03-05 00:00:00", "ABCUSDT", "5"), ("2024-03-05 23:40:00", "ABCUSDT", "7")])
        self.assertIsNone(value)

    def test_nonpositive_end_of_day_value_is_missing(self) -> None:
        value, _ = self.parse([("2024-03-05 23:55:00", "ABCUSDT", "0")])
        self.assertIsNone(value)

    def test_schema_drift_aborts(self) -> None:
        with self.assertRaisesRegex(ValueError, "header drift"):
            self.parse([("2024-03-05 23:55:00", "ABCUSDT", "1")], header=runner.METRICS_HEADER[:-1] + ("other",))
        with self.assertRaisesRegex(ValueError, "outside archive date"):
            self.parse([("2024-03-06 23:55:00", "ABCUSDT", "1")])
        with self.assertRaisesRegex(ValueError, "symbol column mismatch"):
            self.parse([("2024-03-05 23:55:00", "XYZUSDT", "1")])
        with self.assertRaisesRegex(ValueError, "schema value"):
            self.parse([("2024-03-05 23:55:00", "ABCUSDT", "abc")])

    def test_verify_metrics_archive_detects_every_pinned_drift(self) -> None:
        _, entry = self.parse([("2024-03-05 23:55:00", "ABCUSDT", "1")])
        self.assertEqual(runner.verify_metrics_archive(entry, root=self.root), entry["zip_member"])
        for field, value, message in (("size_bytes", 1, "byte size drift"), ("sha256", "0" * 64, "sidecar hash mismatch"),
                                      ("checksum_sidecar_sha256", "0" * 64, "sidecar file hash drift"),
                                      ("zip_member", "other.csv", "zip member mismatch")):
            with self.assertRaisesRegex(ValueError, message):
                runner.verify_metrics_archive({**entry, field: value}, root=self.root)
        with self.assertRaisesRegex(ValueError, "symbol/date disagree"):
            runner.verify_metrics_archive({**entry, "date": "2024-03-06"}, root=self.root)

    def test_loader_counts_missing_end_of_day_rows_and_rejects_duplicates(self) -> None:
        good = write_metrics_zip(self.root, "ABCUSDT", "2024-03-05", [("2024-03-05 23:55:00", "ABCUSDT", "3")])
        empty = write_metrics_zip(self.root, "ABCUSDT", "2024-03-06", [("2024-03-06 12:00:00", "ABCUSDT", "3")])
        inventory = {"files": [good, empty]}
        counters = {"archive-without-end-of-day-row": 0}
        original = runner._parse_metrics_archive
        with mock.patch.object(runner, "WORKSPACE_ROOT", self.root), \
                mock.patch.object(runner, "_parse_metrics_archive", lambda entry: original(entry, root=self.root)):
            values = runner.load_end_of_day_open_interest(inventory, ["ABCUSDT"], 0, 4 * 10 ** 13, counters)
        self.assertEqual(values, {("ABCUSDT", runner.parse_iso_ms("2024-03-05T00:00:00Z")): 3.0})
        self.assertEqual(counters["archive-without-end-of-day-row"], 1)
        with self.assertRaisesRegex(ValueError, "duplicate inventory date"):
            with mock.patch.object(runner, "_parse_metrics_archive", lambda entry: 1.0):
                runner.load_end_of_day_open_interest({"files": [good, dict(good)]}, ["ABCUSDT"], 0, 4 * 10 ** 13, counters)


class UniverseRuleTests(unittest.TestCase):
    def test_window_start_must_be_first_of_january(self) -> None:
        with self.assertRaisesRegex(ValueError, "1 January"):
            runner.universe_lookback(runner.parse_iso_ms("2024-02-01T00:00:00Z"))

    def test_top_symbols_by_preceding_year_volume_with_completeness_floor(self) -> None:
        lookback_start = runner.parse_iso_ms("2023-01-01T00:00:00Z")
        records = {}
        for index in range(runner.UNIVERSE_SIZE + 3):
            symbol = f"S{index:03d}USDT"
            days = runner.MIN_COMPLETE_LOOKBACK_DAYS if index != 0 else runner.MIN_COMPLETE_LOOKBACK_DAYS - 1
            for offset in range(days):
                records[(symbol, lookback_start + offset * DAY_MS)] = day_record(1.0, 1.0, quote=1e6 * (1000 - index))
        records[("S001USDT", lookback_start + 364 * DAY_MS)] = day_record(1.0, 1.0, quote=1e6, complete=False)
        records[("S002USDT", lookback_start + 364 * DAY_MS)] = day_record(1.0, 1.0, quote=1e6 * 997)
        with mock.patch.object(parent, "load_day_records", return_value=records) as load:
            universe, stats = runner.window_universe({"files": []}, START_2024)
        load.assert_called_once_with({"files": []}, lookback_start, START_2024)
        self.assertNotIn("S000USDT", universe)
        self.assertEqual(len(universe), runner.UNIVERSE_SIZE)
        self.assertEqual(universe[:3], ["S002USDT", "S001USDT", "S003USDT"])
        self.assertEqual(stats["S002USDT"]["quote_volume_usd"], 1e6 * 998 * runner.MIN_COMPLETE_LOOKBACK_DAYS + 1e6 * 997)
        self.assertEqual(stats["S001USDT"]["complete_days"], runner.MIN_COMPLETE_LOOKBACK_DAYS)


class DateConstructionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.symbols = [f"S{i:02d}USDT" for i in range(40)]
        self.day = START_2024
        self.records = flat_records(self.symbols, self.day, 2)
        self.open_interest = flat_open_interest(self.symbols, self.day, 2)
        for index, symbol in enumerate(self.symbols):
            self.open_interest[(symbol, self.day - DAY_MS)] = 1000.0 * math.exp(0.01 * index)
        self.counters = {key: 0 for key in ("kline-history-incomplete", "below-liquidity-floor",
                                             "missing-open-interest", "archive-without-end-of-day-row")}
        self.rejected = {key: 0 for key in runner.REJECTION_KEYS}

    def build(self):
        return runner.build_date(self.records, self.open_interest, self.symbols, self.day, self.counters, self.rejected)

    def test_highest_growth_is_long_lowest_is_short_in_quintiles(self) -> None:
        row = self.build()
        self.assertEqual(row["leg_size"], 8)
        longs = set(row["symbols"][:8])
        shorts = set(row["symbols"][8:])
        self.assertEqual(longs, set(self.symbols[-8:]))
        self.assertEqual(shorts, set(self.symbols[:8]))

    def test_ties_break_by_symbol_deterministically(self) -> None:
        for symbol in self.symbols:
            self.open_interest[(symbol, self.day - DAY_MS)] = 1000.0
        row = self.build()
        self.assertEqual(row["symbols"][:8], self.symbols[-8:])
        self.assertEqual(row["symbols"][8:], self.symbols[:8])

    def test_insufficient_cross_section_excludes_date(self) -> None:
        self.assertIsNone(runner.build_date(self.records, self.open_interest, self.symbols[:29], self.day,
                                            self.counters, self.rejected))
        self.assertEqual(self.rejected["insufficient-cross-section"], 1)

    def test_eligibility_exclusions_are_counted_by_cause(self) -> None:
        self.records[("S39USDT", self.day - 5 * DAY_MS)] = day_record(100.0, 100.0, complete=False)
        for back in range(1, runner.ELIGIBILITY_DAYS + 1):
            self.records[("S38USDT", self.day - back * DAY_MS)] = day_record(100.0, 100.0, quote=4.9e6)
        del self.open_interest[("S37USDT", self.day - 2 * DAY_MS)]
        signals = runner.eligible_symbols(self.records, self.open_interest, self.symbols, self.day, self.counters)
        self.assertEqual(len(signals), 37)
        self.assertEqual(self.counters, {"kline-history-incomplete": 1, "below-liquidity-floor": 1,
                                         "missing-open-interest": 1, "archive-without-end-of-day-row": 0})

    def test_fee_algebra_and_excess_definitions(self) -> None:
        for symbol in self.symbols:
            self.records[(symbol, self.day + DAY_MS)] = day_record(101.0, 101.0)
        row = self.build()
        gross = 0.01
        self.assertAlmostEqual(row["long_net"], gross * (1 - FEE) - 2 * FEE)
        self.assertAlmostEqual(row["short_net"], -gross * (1 + FEE) - 2 * FEE)
        self.assertAlmostEqual(row["market_gross"], gross)
        self.assertAlmostEqual(row["long_excess"], row["long_net"] - gross)
        self.assertAlmostEqual(row["short_excess"], row["short_net"] + gross)
        self.assertAlmostEqual(row["spread_net"], (row["long_net"] + row["short_net"]) / 2)
        self.assertAlmostEqual(row["spread_gross"], 0.0)

    def test_missing_entry_bar_drops_position_and_exit_fallback_is_counted(self) -> None:
        self.records[("S39USDT", self.day)]["open0"] = None
        del self.records[("S00USDT", self.day + DAY_MS)]
        row = self.build()
        self.assertEqual(row["long_positions"], 7)
        self.assertEqual(row["short_positions"], 8)
        self.assertEqual(self.rejected["missing-entry-bar"], 1)
        self.assertEqual(self.rejected["early-exit-last-bar"], 1)

    def test_lagged_diagnostic_uses_hour_one_entry(self) -> None:
        for symbol in self.symbols:
            self.records[(symbol, self.day)]["open1"] = 100.0
            self.records[(symbol, self.day + DAY_MS)] = day_record(100.0, 100.0, open1=102.0)
        row = self.build()
        self.assertAlmostEqual(row["spread_net"], -2 * FEE)
        lag_long = 0.02 * (1 - FEE) - 2 * FEE
        lag_short = -0.02 * (1 + FEE) - 2 * FEE
        self.assertAlmostEqual(row["lag_spread_net"], (lag_long + lag_short) / 2)


def synthetic_rows(spreads: list[float], long_excess: float = 0.001, short_excess: float = 0.001) -> list[dict]:
    rows = []
    for index, spread in enumerate(spreads):
        day = START_2024 + index * DAY_MS
        date = runner.date_label(day)
        rows.append({"date": date, "quarter": runner.quarter_label(date), "spread_net": spread, "spread_gross": spread + 2 * FEE,
                     "long_excess": long_excess, "short_excess": short_excess, "long_net": spread, "short_net": spread,
                     "long_gross": spread, "short_gross": spread, "market_gross": 0.0, "lag_spread_net": spread,
                     "eligible": 40, "leg_size": 8, "positions": 16, "long_positions": 8, "short_positions": 8, "symbols": []})
    return rows


class GateTests(unittest.TestCase):
    def test_all_gates_pass_on_strong_stable_positive_spread(self) -> None:
        gates = runner.evaluate_gates(synthetic_rows([0.01] * 366))
        self.assertTrue(all(gate["pass"] for gate in gates.values()))
        self.assertEqual(gates["G4_dependence_robustness"]["alpha"], 0.05 / 30)

    def test_each_gate_falsifies_alone(self) -> None:
        base = [0.01] * 366
        below = runner.evaluate_gates(synthetic_rows([0.0004] * 366))
        self.assertFalse(below["G1_net_magnitude"]["pass"])
        self.assertTrue(below["G4_dependence_robustness"]["pass"])
        legs = runner.evaluate_gates(synthetic_rows(base, short_excess=-0.0001))
        self.assertFalse(legs["G2_leg_contribution"]["pass"])
        self.assertTrue(legs["G1_net_magnitude"]["pass"])
        spreads = [0.05] * 182 + [-0.0001] * 184
        unstable = runner.evaluate_gates(synthetic_rows(spreads))
        self.assertFalse(unstable["G3_quarter_stability"]["pass"])
        self.assertTrue(unstable["G1_net_magnitude"]["pass"])
        noisy = runner.evaluate_gates(synthetic_rows([0.02 if index % 2 else -0.018 for index in range(366)]))
        self.assertTrue(noisy["G1_net_magnitude"]["pass"])
        self.assertFalse(noisy["G4_dependence_robustness"]["pass"])

    def test_quarter_stability_uses_exact_three_fifths_rule(self) -> None:
        three_of_four = runner.evaluate_gates(synthetic_rows([0.01] * 274 + [-0.01] * 92))
        self.assertEqual(three_of_four["G3_quarter_stability"]["positive_quarters"], 3)
        self.assertTrue(three_of_four["G3_quarter_stability"]["pass"])
        two_of_four = runner.evaluate_gates(synthetic_rows([0.01] * 182 + [-0.01] * 184))
        self.assertEqual(two_of_four["G3_quarter_stability"]["positive_quarters"], 2)
        self.assertFalse(two_of_four["G3_quarter_stability"]["pass"])

    def test_bootstrap_is_deterministic_and_one_sided(self) -> None:
        rows = synthetic_rows([0.02 if index % 3 else -0.01 for index in range(366)])
        first = runner.evaluate_gates(rows)["G4_dependence_robustness"]["one_sided_block_bootstrap_p_spread_net"]
        second = runner.evaluate_gates(rows)["G4_dependence_robustness"]["one_sided_block_bootstrap_p_spread_net"]
        self.assertEqual(first, second)
        negative = runner.evaluate_gates(synthetic_rows([-0.01] * 366))["G4_dependence_robustness"]
        self.assertGreater(negative["one_sided_block_bootstrap_p_spread_net"], 0.99)


class WindowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.universe = [f"S{i:02d}USDT" for i in range(runner.UNIVERSE_SIZE)]
        self.contract = synthetic_contract()
        self.inventory = synthetic_inventory(self.universe, self.universe, [])

    def run_primary(self, records, open_interest, universe=None):
        with mock.patch.object(runner, "window_universe", return_value=(universe or self.universe, {})), \
                mock.patch.object(parent, "load_day_records", return_value=records), \
                mock.patch.object(runner, "load_end_of_day_open_interest", return_value=open_interest):
            return runner.run_window(self.contract, "primary", self.inventory, {"files": []})

    def test_universe_or_window_mismatch_aborts(self) -> None:
        with self.assertRaisesRegex(ValueError, "universe recomputed"):
            self.run_primary({}, {}, universe=list(reversed(self.universe)))
        self.inventory["windows"]["primary"]["end_utc_exclusive"] = "2024-12-31T00:00:00Z"
        with self.assertRaisesRegex(ValueError, "window differs"):
            self.run_primary({}, {})

    def test_small_sample_is_inconclusive_without_gates(self) -> None:
        records = flat_records(self.universe, START_2024, 40)
        open_interest = flat_open_interest(self.universe, START_2024, 40)
        result = self.run_primary(records, open_interest)
        self.assertEqual(result["status"], "inconclusive")
        self.assertFalse(result["sample_minimum_met"])
        self.assertEqual(result["gate_results"], {})
        self.assertEqual(result["sample"]["included_dates"], 40)
        self.assertEqual(result["sample"]["positions"], 40 * 2 * (runner.UNIVERSE_SIZE // 5))

    def test_full_year_negative_spread_is_falsified_and_pass_label_is_frozen(self) -> None:
        days = 366
        records = flat_records(self.universe, START_2024, days)
        open_interest = flat_open_interest(self.universe, START_2024, days)
        result = self.run_primary(records, open_interest)
        self.assertEqual(result["status"], "falsified")
        self.assertTrue(result["sample_minimum_met"])
        self.assertAlmostEqual(result["pooled"]["mean_spread_net"], -2 * FEE)
        self.assertEqual(result["sample"]["included_dates"], days)
        self.assertEqual(result["diagnostics"]["funding_timestamp_crossings"]["subtracted"], False)
        for index, symbol in enumerate(self.universe):
            for offset in range(-2, days + 1):
                open_interest[(symbol, START_2024 + offset * DAY_MS)] = 1000.0 * math.exp(0.001 * index * (offset % 2))
        for symbol in self.universe[-14:]:
            for offset in range(days + 1):
                records[(symbol, START_2024 + offset * DAY_MS)] = day_record(100.0, 100.0)
                records[(symbol, START_2024 + (offset + 1) * DAY_MS)] = day_record(105.0, 105.0)
        winner = self.run_primary(records, open_interest)
        self.assertIn(winner["status"], {"pass-on-this-sample", "falsified"})


class LifecycleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-oi5-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)
        self.contract_path = self.directory / "contract.json"
        self.contract_path.write_text(json.dumps(synthetic_contract()), encoding="utf-8")
        self.contract_sha = runner.sha256_file(self.contract_path)
        self.verified = {"input_inventory_sha256": "i" * 64, "kline_inventory_sha256": "k" * 64,
                         "prior_window_audit_sha256": "a" * 64, "loop_policy_sha256": "p" * 64,
                         "parent_runner_sha256": "r" * 64, "metrics_archives_verified": 51224,
                         "kline_archives_verified": 14535}

    def fake_result(self, label, status="falsified"):
        return {"role": label, "window": ["a", "b"], "status": status, "sample_minimum_met": True,
                "sample": {"included_dates": 366, "positions": 10000}, "gate_results": {}, "pooled": {}, "diagnostics": {},
                "claim_boundary": "test", "dates": []}

    def run_main(self, args, status="falsified"):
        with mock.patch.object(runner, "CONTRACT_PATH", self.contract_path), \
                mock.patch.object(runner, "verify_static_inputs", return_value=({}, {}, self.verified)), \
                mock.patch.object(runner, "run_window", side_effect=lambda c, label, i, k: self.fake_result(label, status)):
            return runner.main(args)

    def test_contract_hash_mismatch_and_existing_output_refuse(self) -> None:
        output = self.directory / "primary.json"
        with self.assertRaisesRegex(SystemExit, "contract hash mismatch"):
            self.run_main(["--contract-sha256", "0" * 64, "--window", "primary", "--output", str(output)])
        output.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary", "--output", str(output)])
        self.assertEqual(output.read_text(encoding="utf-8"), "{}")

    def test_verify_inputs_only_writes_nothing(self) -> None:
        output = self.directory / "never.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                                        "--output", str(output), "--verify-inputs-only"]), 0)
        self.assertFalse(output.exists())

    def test_replication_requires_a_committed_primary_pass(self) -> None:
        primary = self.directory / "primary.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary", "--output", str(primary)]), 0)
        payload = json.loads(primary.read_text(encoding="utf-8"))
        self.assertEqual(payload["role"], "primary")
        self.assertEqual(payload["result"]["status"], "falsified")
        replication = self.directory / "replication.json"
        with self.assertRaisesRegex(SystemExit, "requires --primary-artifact"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication)])
        with self.assertRaisesRegex(ValueError, "committed primary pass"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication),
                           "--primary-artifact", str(primary), "--primary-artifact-sha256", runner.sha256_file(primary)])
        with self.assertRaisesRegex(SystemExit, "primary artifact hash mismatch"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication),
                           "--primary-artifact", str(primary), "--primary-artifact-sha256", "0" * 64])
        self.assertFalse(replication.exists())

    def test_primary_pass_then_replication_end_to_end(self) -> None:
        primary = self.directory / "primary.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary", "--output", str(primary)],
                                       status="pass-on-this-sample"), 0)
        replication = self.directory / "replication.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication),
                                        "--primary-artifact", str(primary), "--primary-artifact-sha256", runner.sha256_file(primary)]), 0)
        payload = json.loads(replication.read_text(encoding="utf-8"))
        self.assertEqual(payload["role"], "replication")
        self.assertEqual(payload["primary_artifact_sha256"], runner.sha256_file(primary))
        self.assertEqual(payload["contract_sha256"], self.contract_sha)
        self.assertEqual(payload["runner_sha256"], runner.sha256_file(Path(runner.__file__).resolve()))

    def test_replication_refuses_a_primary_from_a_different_runner(self) -> None:
        primary = self.directory / "primary.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary", "--output", str(primary)],
                                       status="pass-on-this-sample"), 0)
        payload = json.loads(primary.read_text(encoding="utf-8"))
        payload["runner_sha256"] = "9" * 64
        edited = self.directory / "primary-edited.json"
        edited.write_text(json.dumps(payload), encoding="utf-8")
        replication = self.directory / "replication.json"
        with self.assertRaisesRegex(ValueError, "runner_sha256 mismatch"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication),
                           "--primary-artifact", str(edited), "--primary-artifact-sha256", runner.sha256_file(edited)])
        self.assertFalse(replication.exists())

    def test_replication_refuses_drifted_static_input_pins(self) -> None:
        primary = self.directory / "primary.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary", "--output", str(primary)],
                                       status="pass-on-this-sample"), 0)
        payload = json.loads(primary.read_text(encoding="utf-8"))
        payload["static_inputs"]["input_inventory_sha256"] = "8" * 64
        edited = self.directory / "primary-drifted.json"
        edited.write_text(json.dumps(payload), encoding="utf-8")
        replication = self.directory / "replication.json"
        with self.assertRaisesRegex(SystemExit, "static-input pins differ"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication", "--output", str(replication),
                           "--primary-artifact", str(edited), "--primary-artifact-sha256", runner.sha256_file(edited)])
        self.assertFalse(replication.exists())

    def test_verify_inputs_only_prints_every_pin(self) -> None:
        stream = io.StringIO()
        output = self.directory / "never2.json"
        with contextlib.redirect_stdout(stream):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                           "--output", str(output), "--verify-inputs-only"])
        printed = json.loads(stream.getvalue())
        for field in runner.PIN_FIELDS:
            self.assertEqual(printed[field], self.verified[field])
        self.assertEqual(printed["runner_sha256"], runner.sha256_file(Path(runner.__file__).resolve()))


if __name__ == "__main__":
    unittest.main()
