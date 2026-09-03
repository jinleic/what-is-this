#!/usr/bin/env python3
"""Behavioral tests for the frozen cross-sectional daily-reversal runner."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_perp_cross_section_daily_reversal as runner


DAY_MS = 86_400_000
HOUR_MS = 3_600_000
FEE = 0.0005
FULL = runner.FULL_DAY_MASK
DAY0 = runner.parse_iso_ms("2022-01-01T00:00:00Z")


def rec(open0: float = 100.0, open1: float = 100.0, close23: float = 100.0,
        last_close: float | None = None, qvol: float = 1e7,
        hours: int = FULL, duplicate: bool = False, bad: bool = False) -> dict:
    return {"hours": hours, "duplicate": duplicate, "bad_price": bad,
            "open0": open0, "open1": open1, "close23": close23,
            "last_hour": 23 if hours == FULL else 0,
            "last_close": close23 if last_close is None else last_close,
            "quote_volume": qvol * 24.0 / 24.0 if hours == FULL else qvol,
            "trades": 24}


def make_universe(day_ms: int, signals: list[float], outcomes: list[float] | None = None,
                  qvol: float = 1e7) -> tuple[dict, list[str]]:
    """Complete 20-day histories, a signal day, an entry day, and an exit day."""
    outcomes = outcomes or [0.0] * len(signals)
    records: dict = {}
    symbols = [f"S{index:03d}USDT" for index in range(len(signals))]
    for symbol, signal, outcome in zip(symbols, signals, outcomes):
        for back in range(2, 22):
            records[(symbol, day_ms - back * DAY_MS)] = rec(qvol=qvol)
        records[(symbol, day_ms - DAY_MS)] = rec(open0=100.0, close23=100.0 * (1 + signal), qvol=qvol)
        records[(symbol, day_ms)] = rec(open0=100.0, open1=100.0, close23=100.0 * (1 + outcome),
                                        last_close=100.0 * (1 + outcome))
        records[(symbol, day_ms + DAY_MS)] = rec(open0=100.0 * (1 + outcome),
                                                 open1=100.0 * (1 + outcome))
    return records, symbols


def frozen_test_contract(days: int = 730) -> dict:
    contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
    end = DAY0 + days * DAY_MS
    contract["windows"]["primary"] = {
        "start_utc_inclusive": runner.iso_ms(DAY0),
        "end_utc_exclusive": runner.iso_ms(end),
        "expected_dates": days,
    }
    return contract


def synthetic_day_rows(days: int, spread_fn, long_excess: float | None = None,
                       short_excess: float | None = None, positions: int = 20) -> list[dict]:
    rows = []
    for index in range(days):
        day_ms = DAY0 + index * DAY_MS
        date = runner.date_label(day_ms)
        spread = spread_fn(index)
        le = spread if long_excess is None else long_excess
        se = spread if short_excess is None else short_excess
        held = [{"date": date, "symbol": f"S{i:03d}USDT", "side": "long" if i % 2 == 0 else "short",
                 "signal": 0.0, "gross": spread, "net": spread, "lag_net": spread}
                for i in range(positions)]
        rows.append({"date": date, "quarter": runner.quarter_label(date), "eligible": 100,
                     "leg_size": 10, "long_positions": 10, "short_positions": 10,
                     "market_gross": 0.0, "long_net": le, "short_net": se,
                     "spread_net": spread, "spread_gross": spread + 2 * FEE,
                     "long_gross": spread + FEE, "short_gross": spread + FEE,
                     "long_excess": le, "short_excess": se,
                     "lag_spread_net": spread, "positions": held})
    return rows


def run_with_rows(rows: list[dict], days: int = 730) -> dict:
    contract = frozen_test_contract(days)
    by_day = {runner.parse_iso_ms(row["date"] + "T00:00:00Z"): row for row in rows}

    def fake_build(records, symbols, day_ms, counters):
        return by_day.get(day_ms)

    with mock.patch.object(runner, "load_day_records", return_value={}), \
            mock.patch.object(runner, "build_date", side_effect=fake_build):
        return runner.run_window(contract, "primary", {"symbols": ["AUSDT"], "files": []})


class DateGridTests(unittest.TestCase):
    def test_expected_dates_cover_frozen_windows_exactly(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        for label, count in (("primary", 730), ("replication", 731)):
            window = contract["windows"][label]
            dates = runner.expected_dates(
                runner.parse_iso_ms(window["start_utc_inclusive"]),
                runner.parse_iso_ms(window["end_utc_exclusive"]))
            self.assertEqual(len(dates), count)
            self.assertEqual(window["expected_dates"], count)
            self.assertTrue(all(day % DAY_MS == 0 for day in dates))
        primary = contract["windows"]["primary"]
        replication = contract["windows"]["replication"]
        self.assertEqual(primary["end_utc_exclusive"], replication["start_utc_inclusive"])
        with self.assertRaisesRegex(ValueError, "midnight"):
            runner.expected_dates(DAY0 + HOUR_MS, DAY0 + DAY_MS)

    def test_quarter_labels_follow_calendar_quarters(self) -> None:
        self.assertEqual(runner.quarter_label("2022-01-01"), "2022Q1")
        self.assertEqual(runner.quarter_label("2022-03-31"), "2022Q1")
        self.assertEqual(runner.quarter_label("2022-04-01"), "2022Q2")
        self.assertEqual(runner.quarter_label("2025-12-31"), "2025Q4")

    def test_contract_read_boundary_months_are_the_frozen_span(self) -> None:
        months = runner.window_months(DAY0 - 20 * DAY_MS,
                                      runner.parse_iso_ms("2024-01-01T00:00:00Z") + 2 * HOUR_MS)
        self.assertEqual(months[0], "2021-12")
        self.assertEqual(months[-1], "2024-01")
        self.assertEqual(len(months), 26)


class CostTests(unittest.TestCase):
    def test_net_return_charges_both_taker_fees_on_traded_notional(self) -> None:
        self.assertAlmostEqual(runner.net_return(0.0, "long"), -2 * FEE, places=15)
        self.assertAlmostEqual(runner.net_return(0.0, "short"), -2 * FEE, places=15)
        self.assertAlmostEqual(runner.net_return(0.01, "long"),
                               1.01 * (1 - FEE) - (1 + FEE), places=15)
        self.assertAlmostEqual(runner.net_return(0.01, "short"),
                               -0.01 - FEE - FEE * 1.01, places=15)
        # A short of a faller gains less than the fall by the fee on both notionals.
        self.assertLess(runner.net_return(-0.02, "short"), 0.02)
        with self.assertRaises(ValueError):
            runner.net_return(0.0, "flat")


class ParseTests(unittest.TestCase):
    def _archive(self, rows: list[list[str]], header: bool) -> tuple[Path, dict]:
        directory = Path(tempfile.mkdtemp(prefix=".tmp-xs4-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, directory, ignore_errors=True)
        name = "ZZZUSDT-1h-2022-01.zip"
        member = "ZZZUSDT-1h-2022-01.csv"
        buffer = io.StringIO()
        if header:
            buffer.write(",".join(runner.EXPECTED_HEADER) + "\n")
        for row in rows:
            buffer.write(",".join(row) + "\n")
        path = directory / name
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(member, buffer.getvalue())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        sidecar = directory / f"{name}.CHECKSUM"
        sidecar.write_text(f"{digest}  {name}\n", encoding="utf-8")
        entry = {"path": str(path.relative_to(runner.WORKSPACE_ROOT)),
                 "checksum_sidecar": str(sidecar.relative_to(runner.WORKSPACE_ROOT)),
                 "checksum_sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
                 "sha256": digest, "size_bytes": path.stat().st_size, "zip_member": member,
                 "header": ",".join(runner.EXPECTED_HEADER), "header_present": header}
        return path, entry

    @staticmethod
    def _bar(open_ms: int, open_price: float = 100.0, close: float = 101.0,
             qvol: float = 5e5, micro: bool = False) -> list[str]:
        stamp = str(open_ms * 1000 if micro else open_ms)
        return [stamp, f"{open_price}", "102", "99", f"{close}", "10",
                str(open_ms + HOUR_MS - 1), f"{qvol}", "7", "5", f"{qvol / 2}", "0"]

    def test_parse_archive_reads_both_header_forms_and_normalizes_microseconds(self) -> None:
        for header, micro in ((True, False), (False, True)):
            rows = [self._bar(DAY0 + h * HOUR_MS, close=101.0 if h == 23 else 100.0,
                              micro=micro) for h in range(24)]
            _, entry = self._archive(rows, header)
            records: dict = {}
            runner._parse_archive(entry, DAY0 - DAY_MS, DAY0 + DAY_MS, records)
            record = records[("ZZZUSDT", DAY0)]
            self.assertTrue(runner.is_complete(record))
            self.assertEqual(record["open0"], 100.0)
            self.assertEqual(record["close23"], 101.0)
            self.assertEqual(record["trades"], 24 * 7)
            self.assertAlmostEqual(record["quote_volume"], 24 * 5e5)

    def test_duplicate_hour_keeps_first_occurrence_and_flags_day(self) -> None:
        rows = [self._bar(DAY0 + h * HOUR_MS) for h in range(24)]
        rows.insert(5, self._bar(DAY0 + 4 * HOUR_MS, open_price=999.0))
        _, entry = self._archive(rows, True)
        records: dict = {}
        runner._parse_archive(entry, DAY0, DAY0 + DAY_MS, records)
        record = records[("ZZZUSDT", DAY0)]
        self.assertTrue(record["duplicate"])
        self.assertFalse(runner.is_complete(record))
        self.assertEqual(record["open0"], 100.0)
        self.assertEqual(record["hours"], FULL)

    def test_nonpositive_price_flags_day_but_never_aborts(self) -> None:
        rows = [self._bar(DAY0 + h * HOUR_MS, open_price=0.0 if h == 3 else 100.0)
                for h in range(24)]
        _, entry = self._archive(rows, True)
        records: dict = {}
        runner._parse_archive(entry, DAY0, DAY0 + DAY_MS, records)
        self.assertTrue(records[("ZZZUSDT", DAY0)]["bad_price"])
        self.assertFalse(runner.is_complete(records[("ZZZUSDT", DAY0)]))

    def test_bars_outside_read_boundary_are_ignored(self) -> None:
        rows = [self._bar(DAY0 + h * HOUR_MS) for h in range(24)]
        rows += [self._bar(DAY0 + DAY_MS + h * HOUR_MS) for h in range(24)]
        _, entry = self._archive(rows, True)
        records: dict = {}
        runner._parse_archive(entry, DAY0, DAY0 + DAY_MS + 2 * HOUR_MS, records)
        self.assertEqual(records[("ZZZUSDT", DAY0 + DAY_MS)]["hours"], 0b11)
        self.assertEqual(sorted(day for _, day in records), [DAY0, DAY0 + DAY_MS])

    def test_schema_drift_aborts(self) -> None:
        rows = [self._bar(DAY0 + h * HOUR_MS) for h in range(24)]
        rows[7] = rows[7][:-1]
        _, entry = self._archive(rows, True)
        with self.assertRaisesRegex(ValueError, "column count"):
            runner._parse_archive(entry, DAY0, DAY0 + DAY_MS, {})
        rows = [self._bar(DAY0 + h * HOUR_MS) for h in range(24)]
        _, entry = self._archive(rows, True)
        entry["header_present"] = False
        with self.assertRaisesRegex(ValueError, "first-row drift"):
            runner._parse_archive(entry, DAY0, DAY0 + DAY_MS, {})

    def test_verify_archive_rejects_every_drift(self) -> None:
        rows = [self._bar(DAY0 + h * HOUR_MS) for h in range(24)]
        path, entry = self._archive(rows, True)
        self.assertEqual(runner.verify_archive(entry), entry["zip_member"])
        for field, value, message in (
                ("size_bytes", entry["size_bytes"] + 1, "byte size drift"),
                ("sha256", "0" * 64, "sidecar hash mismatch"),
                ("checksum_sidecar_sha256", "0" * 64, "sidecar file hash drift"),
                ("zip_member", "other.csv", "zip member mismatch"),
                ("header_present", False, "first-row form mismatch"),
                ("header", "open_time,open", "inventory header mismatch")):
            broken = {**entry, field: value}
            with self.assertRaisesRegex(ValueError, message):
                runner.verify_archive(broken)
        path.write_bytes(path.read_bytes() + b"\0")
        with self.assertRaisesRegex(ValueError, "byte size drift"):
            runner.verify_archive(entry)


class EligibilityTests(unittest.TestCase):
    def test_eligibility_needs_twenty_complete_days_and_liquidity_floor(self) -> None:
        records, symbols = make_universe(DAY0, [0.01, -0.01])
        eligible = runner.eligible_symbols(records, symbols, DAY0)
        self.assertEqual(set(eligible), set(symbols))
        self.assertAlmostEqual(eligible["S000USDT"], 0.01, places=12)
        records[("S000USDT", DAY0 - 20 * DAY_MS)] = rec(hours=FULL - 1)
        self.assertNotIn("S000USDT", runner.eligible_symbols(records, symbols, DAY0))
        low, symbols_low = make_universe(DAY0, [0.0], qvol=runner.LIQUIDITY_FLOOR_USD - 1.0)
        self.assertEqual(runner.eligible_symbols(low, symbols_low, DAY0), {})
        exact, symbols_exact = make_universe(DAY0, [0.0], qvol=runner.LIQUIDITY_FLOOR_USD)
        self.assertEqual(list(runner.eligible_symbols(exact, symbols_exact, DAY0)), symbols_exact)

    def test_eligibility_and_signal_use_only_pre_midnight_bars(self) -> None:
        records, symbols = make_universe(DAY0, [0.02])
        before = runner.eligible_symbols(records, symbols, DAY0)
        records[("S000USDT", DAY0)] = rec(open0=0.0, hours=0, bad=True)
        del records[("S000USDT", DAY0 + DAY_MS)]
        self.assertEqual(runner.eligible_symbols(records, symbols, DAY0), before)

    def test_duplicate_or_bad_history_day_blocks_eligibility(self) -> None:
        for flag in ("duplicate", "bad"):
            records, symbols = make_universe(DAY0, [0.0])
            records[("S000USDT", DAY0 - 7 * DAY_MS)] = rec(**{flag: True}) if flag == "duplicate" \
                else rec(bad=True)
            self.assertEqual(runner.eligible_symbols(records, symbols, DAY0), {})


class RankingAndExecutionTests(unittest.TestCase):
    def test_losers_go_long_winners_go_short_in_floor_deciles(self) -> None:
        signals = [(index - 30) / 1000.0 for index in range(63)]
        records, symbols = make_universe(DAY0, signals)
        counters = {key: 0 for key in ("insufficient-cross-section", "empty-leg",
                                       "missing-entry-bar", "unvaluable-exit",
                                       "early-exit-last-bar", "lag-missing-entry-bar",
                                       "lag-unvaluable-exit", "lag-exit-at-next-midnight",
                                       "lag-early-exit-last-bar")}
        row = runner.build_date(records, symbols, DAY0, counters)
        self.assertEqual(row["eligible"], 63)
        self.assertEqual(row["leg_size"], 6)
        longs = sorted(p["symbol"] for p in row["positions"] if p["side"] == "long")
        shorts = sorted(p["symbol"] for p in row["positions"] if p["side"] == "short")
        self.assertEqual(longs, [f"S{i:03d}USDT" for i in range(6)])
        self.assertEqual(shorts, [f"S{i:03d}USDT" for i in range(57, 63)])
        self.assertEqual(sum(counters.values()), 0)

    def test_ties_break_by_symbol_name(self) -> None:
        records, symbols = make_universe(DAY0, [0.0] * 60)
        counters = {key: 0 for key in ("insufficient-cross-section", "empty-leg",
                                       "missing-entry-bar", "unvaluable-exit",
                                       "early-exit-last-bar", "lag-missing-entry-bar",
                                       "lag-unvaluable-exit", "lag-exit-at-next-midnight",
                                       "lag-early-exit-last-bar")}
        row = runner.build_date(records, symbols, DAY0, counters)
        longs = [p["symbol"] for p in row["positions"] if p["side"] == "long"]
        self.assertEqual(longs, symbols[:6])

    def test_insufficient_cross_section_excludes_date(self) -> None:
        records, symbols = make_universe(DAY0, [0.0] * 49)
        counters = {"insufficient-cross-section": 0}
        self.assertIsNone(runner.build_date(records, symbols, DAY0, counters))
        self.assertEqual(counters["insufficient-cross-section"], 1)

    def test_positions_earn_reversal_net_of_fees_and_market_excess_decomposes(self) -> None:
        signals = [(index - 30) / 1000.0 for index in range(60)]
        outcomes = [-signal for signal in signals]  # perfect reversal
        records, symbols = make_universe(DAY0, signals, outcomes)
        counters = {key: 0 for key in ("insufficient-cross-section", "empty-leg",
                                       "missing-entry-bar", "unvaluable-exit",
                                       "early-exit-last-bar", "lag-missing-entry-bar",
                                       "lag-unvaluable-exit", "lag-exit-at-next-midnight",
                                       "lag-early-exit-last-bar")}
        row = runner.build_date(records, symbols, DAY0, counters)
        long_gross = sum(-signals[i] for i in range(6)) / 6
        short_gross = sum(signals[i] for i in range(54, 60)) / 6
        self.assertAlmostEqual(row["long_gross"], long_gross, places=12)
        self.assertAlmostEqual(row["short_gross"], short_gross, places=12)
        self.assertAlmostEqual(row["long_net"], long_gross * (1 - FEE) - 2 * FEE, places=12)
        self.assertAlmostEqual(row["short_net"], -(-short_gross) * (1 + FEE) - 2 * FEE, places=12)
        self.assertAlmostEqual(row["spread_net"], (row["long_net"] + row["short_net"]) / 2, places=15)
        self.assertAlmostEqual(row["long_excess"] + row["short_excess"], 2 * row["spread_net"], places=15)
        self.assertAlmostEqual(row["market_gross"], sum(outcomes) / 60, places=12)
        self.assertAlmostEqual(row["lag_spread_net"], row["spread_net"], places=12)

    def test_missing_entry_bar_skips_position_and_counts(self) -> None:
        records, symbols = make_universe(DAY0, [(index - 30) / 1000.0 for index in range(60)])
        records[("S000USDT", DAY0)] = rec(open0=0.0, hours=FULL - 1)
        counters = {key: 0 for key in ("insufficient-cross-section", "empty-leg",
                                       "missing-entry-bar", "unvaluable-exit",
                                       "early-exit-last-bar", "lag-missing-entry-bar",
                                       "lag-unvaluable-exit", "lag-exit-at-next-midnight",
                                       "lag-early-exit-last-bar")}
        row = runner.build_date(records, symbols, DAY0, counters)
        self.assertEqual(row["long_positions"], 5)
        self.assertEqual(counters["missing-entry-bar"], 1)
        self.assertNotIn("S000USDT", [p["symbol"] for p in row["positions"]])

    def test_exit_chain_prefers_next_midnight_then_last_bar_then_skips(self) -> None:
        records = {("AUSDT", DAY0): rec(open0=100.0, last_close=97.0),
                   ("AUSDT", DAY0 + DAY_MS): rec(open0=103.0, open1=104.0)}
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0), (103.0, None))
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0, lagged=True), (104.0, None))
        records[("AUSDT", DAY0 + DAY_MS)]["open1"] = None
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0, lagged=True),
                         (103.0, "lag-exit-at-next-midnight"))
        del records[("AUSDT", DAY0 + DAY_MS)]
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0), (97.0, "early-exit-last-bar"))
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0, lagged=True),
                         (97.0, "lag-early-exit-last-bar"))
        records[("AUSDT", DAY0)]["last_close"] = 0.0
        self.assertEqual(runner.exit_price(records, "AUSDT", DAY0), (None, None))

    def test_delisted_symbol_closes_at_last_bar_and_is_counted(self) -> None:
        signals = [(index - 30) / 1000.0 for index in range(60)]
        records, symbols = make_universe(DAY0, signals)
        del records[("S000USDT", DAY0 + DAY_MS)]
        records[("S000USDT", DAY0)]["last_close"] = 90.0
        counters = {key: 0 for key in ("insufficient-cross-section", "empty-leg",
                                       "missing-entry-bar", "unvaluable-exit",
                                       "early-exit-last-bar", "lag-missing-entry-bar",
                                       "lag-unvaluable-exit", "lag-exit-at-next-midnight",
                                       "lag-early-exit-last-bar")}
        row = runner.build_date(records, symbols, DAY0, counters)
        held = next(p for p in row["positions"] if p["symbol"] == "S000USDT")
        self.assertAlmostEqual(held["gross"], -0.10, places=12)
        self.assertEqual(counters["early-exit-last-bar"], 1)
        self.assertEqual(counters["lag-early-exit-last-bar"], 1)
        self.assertEqual(row["long_positions"], 6)


class GateTests(unittest.TestCase):
    def test_all_four_gates_must_pass_for_promotion(self) -> None:
        result = run_with_rows(synthetic_day_rows(730, lambda i: 0.0020))
        self.assertEqual(result["status"], "pass-on-this-sample")
        self.assertTrue(all(gate["pass"] for gate in result["gate_results"].values()))
        self.assertEqual(set(result["gate_results"]),
                         {"G1_net_magnitude", "G2_leg_contribution",
                          "G3_quarter_stability", "G4_dependence_robustness"})
        self.assertEqual(result["sample"]["positions"], 730 * 20)

    def test_fee_sized_spread_fails_magnitude_gate_only(self) -> None:
        result = run_with_rows(synthetic_day_rows(730, lambda i: 0.0004))
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items() if not gate["pass"]}
        self.assertEqual(failed, {"G1_net_magnitude"})

    def test_one_negative_leg_excess_falsifies_window(self) -> None:
        rows = synthetic_day_rows(730, lambda i: 0.0020, long_excess=0.0050,
                                  short_excess=-0.0010)
        result = run_with_rows(rows)
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items() if not gate["pass"]}
        self.assertEqual(failed, {"G2_leg_contribution"})

    def test_quarter_instability_alone_falsifies_window(self) -> None:
        def spread(index: int) -> float:
            quarter = runner.quarter_label(runner.date_label(DAY0 + index * DAY_MS))
            return 0.0030 if quarter in ("2022Q1", "2022Q3", "2023Q1", "2023Q3") else -0.0005
        result = run_with_rows(synthetic_day_rows(730, spread))
        self.assertEqual(result["status"], "falsified")
        failed = {name for name, gate in result["gate_results"].items() if not gate["pass"]}
        self.assertEqual(failed, {"G3_quarter_stability"})
        self.assertEqual(result["gate_results"]["G3_quarter_stability"]["spread_net"]["positive_quarters"], 4)

    def test_quarter_stability_uses_exact_three_fifths_rule(self) -> None:
        def spread_for(positive: set[str]):
            return lambda i: 0.0030 if runner.quarter_label(
                runner.date_label(DAY0 + i * DAY_MS)) in positive else -0.0005
        five = {"2022Q1", "2022Q2", "2022Q3", "2022Q4", "2023Q1"}
        self.assertTrue(run_with_rows(synthetic_day_rows(730, spread_for(five)))
                        ["gate_results"]["G3_quarter_stability"]["pass"])
        four = {"2022Q1", "2022Q2", "2022Q3", "2022Q4"}
        self.assertFalse(run_with_rows(synthetic_day_rows(730, spread_for(four)))
                         ["gate_results"]["G3_quarter_stability"]["pass"])

    def test_insufficient_sample_is_inconclusive_without_gates(self) -> None:
        result = run_with_rows(synthetic_day_rows(249, lambda i: 0.0020), days=249)
        self.assertEqual(result["status"], "inconclusive")
        self.assertIsNone(result["gate_results"])
        thin = synthetic_day_rows(260, lambda i: 0.0020, positions=4)
        result = run_with_rows(thin, days=260)
        self.assertEqual(result["status"], "inconclusive")

    def test_diagnostics_never_rescue_a_false_gate(self) -> None:
        rows = synthetic_day_rows(730, lambda i: -0.0003)
        for row in rows:
            row["spread_gross"] = 0.0100
            row["lag_spread_net"] = 0.0100
        result = run_with_rows(rows)
        self.assertEqual(result["status"], "falsified")
        self.assertGreater(result["diagnostics"]["gross_spread"]["mean"], 0.0)
        self.assertFalse(result["cost_model"]["funding_subtracted"])
        self.assertEqual(result["diagnostics"]["funding_timestamp_crossings"]["total_upper_bound"],
                         3 * result["sample"]["positions"])

    def test_block_bootstrap_is_deterministic_and_one_sided(self) -> None:
        dates = {runner.date_label(DAY0 + i * DAY_MS): 0.001 for i in range(100)}
        p_positive = runner.block_bootstrap_p(dates, replicates=500)
        self.assertEqual(p_positive, 1 / 501)
        self.assertEqual(runner.block_bootstrap_p(dates, replicates=500), p_positive)
        negative = {date: -0.001 for date in dates}
        self.assertEqual(runner.block_bootstrap_p(negative, replicates=500), 1.0)
        mixed = {date: (0.02 if i % 2 else -0.019) for i, date in enumerate(dates)}
        p_mixed = runner.block_bootstrap_p(mixed, replicates=2000)
        self.assertGreater(p_mixed, 0.05)
        with self.assertRaisesRegex(ValueError, "too few dates"):
            runner.block_bootstrap_p({d: 0.0 for d in list(dates)[:20]}, replicates=10)


class ContractBindingTests(unittest.TestCase):
    def test_contract_parameters_must_match_runner_constants(self) -> None:
        contract = frozen_test_contract()
        runner.check_contract_parameters(contract)
        contract["frozen_rule"]["parameters"]["liquidity_floor_mean_quote_volume_usd"] = 1.0
        with self.assertRaisesRegex(ValueError, "differ from runner constants"):
            runner.check_contract_parameters(contract)

    def test_static_input_hash_drift_aborts_before_any_archive_check(self) -> None:
        contract = frozen_test_contract()
        contract["assets_and_data"]["input_inventory_sha256"] = "0" * 64
        with mock.patch.object(runner, "verify_archive") as verify:
            with self.assertRaisesRegex(ValueError, "input_inventory_sha256 mismatch"):
                runner.verify_static_inputs(contract)
        verify.assert_not_called()

    def test_inventory_must_certify_pre_outcome_state_and_universe(self) -> None:
        contract = frozen_test_contract()
        good = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        directory = Path(tempfile.mkdtemp(prefix=".tmp-xs4-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, directory, ignore_errors=True)
        fake_inventory = directory / "inventory.json"
        pins = {fake_inventory: contract["assets_and_data"]["input_inventory_sha256"],
                runner.AUDIT_PATH: contract["novelty_and_consumption"]["prior_window_audit_sha256"],
                runner.POLICY_PATH: contract["loop_multiplicity"]["policy_sha256"]}
        for mutate, message in (
                (lambda inv: inv.update(candidate_outcomes_accessed=True), "pre-outcome"),
                (lambda inv: inv.update(symbols=list(reversed(inv["symbols"]))), "sorted unique"),
                (lambda inv: inv.update(archive_count=inv["archive_count"] - 1), "archive count"),
                (lambda inv: inv["files"].append(dict(inv["files"][0])), "duplicate paths"),
                (lambda inv: inv["files"].append({**inv["files"][0], "path": inv["files"][0]["path"].replace(
                    inv["files"][0]["path"].rsplit("/", 1)[1], "NOPEUSDT-1h-2022-01.zip")}), "outside universe")):
            inventory = json.loads(json.dumps(good))
            mutate(inventory)
            if message != "archive count":
                inventory["archive_count"] = len(inventory["files"])
            fake_inventory.write_text(json.dumps(inventory), encoding="utf-8")
            with mock.patch.object(runner, "INVENTORY_PATH", fake_inventory), \
                    mock.patch.object(runner, "sha256_file", side_effect=lambda p: pins[p]), \
                    mock.patch.object(runner, "verify_archive"):
                with self.assertRaisesRegex(ValueError, message):
                    runner.verify_static_inputs(contract)


class LifecycleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-xs4-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)
        self.contract_path = self.directory / "contract.json"
        contract = frozen_test_contract()
        contract["windows"]["replication"] = {
            "start_utc_inclusive": runner.iso_ms(DAY0 + 730 * DAY_MS),
            "end_utc_exclusive": runner.iso_ms(DAY0 + 1460 * DAY_MS),
            "expected_dates": 730}
        self.contract_path.write_text(json.dumps(contract), encoding="utf-8")
        self.contract_sha = hashlib.sha256(self.contract_path.read_bytes()).hexdigest()
        self.integrity = {"input_inventory_sha256": "i" * 64, "prior_window_audit_sha256": "a" * 64,
                          "loop_policy_sha256": "p" * 64, "archives_verified": 1}

    def _patched(self, spread: float):
        rows = synthetic_day_rows(1460, lambda i: spread)
        by_day = {runner.parse_iso_ms(row["date"] + "T00:00:00Z"): row for row in rows}
        return (mock.patch.object(runner, "verify_static_inputs",
                                  return_value=({"symbols": ["AUSDT"], "files": []}, self.integrity)),
                mock.patch.object(runner, "load_day_records", return_value={}),
                mock.patch.object(runner, "build_date",
                                  side_effect=lambda r, s, d, c: by_day.get(d)))

    def test_main_refuses_existing_output_and_contract_hash_mismatch(self) -> None:
        output = self.directory / "out.json"
        output.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "already exists"):
            runner.main(["--contract", str(self.contract_path), "--contract-sha256",
                         self.contract_sha, "--window", "primary", "--output", str(output)])
        with self.assertRaisesRegex(SystemExit, "contract hash mismatch"):
            runner.main(["--contract", str(self.contract_path), "--contract-sha256", "0" * 64,
                         "--window", "primary", "--output", str(self.directory / "new.json")])

    def test_main_primary_then_replication_end_to_end(self) -> None:
        primary = self.directory / "primary.json"
        replication = self.directory / "replication.json"
        patches = self._patched(0.0020)
        with patches[0], patches[1], patches[2]:
            self.assertEqual(runner.main([
                "--contract", str(self.contract_path), "--contract-sha256", self.contract_sha,
                "--window", "primary", "--output", str(primary)]), 0)
            artifact = json.loads(primary.read_text(encoding="utf-8"))
            self.assertEqual(artifact["role"], "primary")
            self.assertEqual(artifact["result"]["status"], "pass-on-this-sample")
            self.assertEqual(artifact["contract_sha256"], self.contract_sha)
            primary_sha = hashlib.sha256(primary.read_bytes()).hexdigest()
            with self.assertRaisesRegex(SystemExit, "requires --primary-artifact"):
                runner.main(["--contract", str(self.contract_path), "--contract-sha256",
                             self.contract_sha, "--window", "replication",
                             "--output", str(replication)])
            with self.assertRaisesRegex(SystemExit, "does not match committed hash"):
                runner.main(["--contract", str(self.contract_path), "--contract-sha256",
                             self.contract_sha, "--window", "replication",
                             "--output", str(replication), "--primary-artifact", str(primary),
                             "--primary-artifact-sha256", "0" * 64])
            self.assertEqual(runner.main([
                "--contract", str(self.contract_path), "--contract-sha256", self.contract_sha,
                "--window", "replication", "--output", str(replication),
                "--primary-artifact", str(primary), "--primary-artifact-sha256", primary_sha]), 0)
            replicated = json.loads(replication.read_text(encoding="utf-8"))
            self.assertEqual(replicated["role"], "replication")
            self.assertEqual(replicated["primary_artifact_sha256"], primary_sha)
            self.assertEqual(replicated["result"]["window"]["start_utc_inclusive"],
                             runner.iso_ms(DAY0 + 730 * DAY_MS))

    def test_main_replication_refuses_falsified_primary(self) -> None:
        primary = self.directory / "primary.json"
        patches = self._patched(-0.0010)
        with patches[0], patches[1], patches[2]:
            runner.main(["--contract", str(self.contract_path), "--contract-sha256",
                         self.contract_sha, "--window", "primary", "--output", str(primary)])
            self.assertEqual(json.loads(primary.read_text())["result"]["status"], "falsified")
            primary_sha = hashlib.sha256(primary.read_bytes()).hexdigest()
            with self.assertRaisesRegex(SystemExit, "requires a committed primary pass"):
                runner.main(["--contract", str(self.contract_path), "--contract-sha256",
                             self.contract_sha, "--window", "replication",
                             "--output", str(self.directory / "replication.json"),
                             "--primary-artifact", str(primary),
                             "--primary-artifact-sha256", primary_sha])

    def test_verify_inputs_only_writes_nothing(self) -> None:
        output = self.directory / "never.json"
        patches = self._patched(0.0)
        with patches[0], patches[1], patches[2]:
            self.assertEqual(runner.main([
                "--contract", str(self.contract_path), "--contract-sha256", self.contract_sha,
                "--window", "primary", "--output", str(output), "--verify-inputs-only"]), 0)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
