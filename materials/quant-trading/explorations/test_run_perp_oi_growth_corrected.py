#!/usr/bin/env python3
"""Behavioral tests for the corrected cross-sectional open-interest-growth runner (cycle 6).

The boundary-row tests are the load-bearing new ones: one drives a real 289-row archive from
the frozen inventory, the rest pin the drop-and-count semantics, the never-selectable rule, and
the surviving fail-closed checks.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import run_perp_cross_section_daily_reversal as parent
import run_perp_oi_growth_cross_section as oi5
import run_perp_oi_growth_corrected as runner

DAY_MS = runner.DAY_MS
HOUR_MS = runner.HOUR_MS
FULL = parent.FULL_DAY_MASK
CUTOFF = runner.END_OF_DAY_CUTOFF_MS
REAL_BOUNDARY_ARCHIVE = ("1000SHIBUSDT", "2024-04-03")


def fresh_census() -> dict[str, int]:
    return {key: 0 for key in ("archives-read", "archives-with-off-date-row",
                               "off-date-boundary-row", "off-date-row-at-or-after-cutoff")}


def day_record(open0: float, close: float, quote: float = 1e7, open1: float | None = None,
               complete: bool = True) -> dict:
    return {"hours": FULL if complete else FULL - 1, "duplicate": False, "bad_price": False,
            "open0": open0, "open1": open1 if open1 is not None else open0, "close23": close,
            "last_hour": 23, "last_close": close, "quote_volume": quote, "trades": 100}


def write_metrics_zip(directory: Path, symbol: str, date: str, rows: list[tuple[str, str, str]],
                      header: tuple[str, ...] = runner.METRICS_HEADER) -> dict:
    """One metrics archive plus its inventory entry; rows are (create_time, symbol, open_interest)."""
    member = f"{symbol}-metrics-{date}.csv"
    path = directory / f"{symbol}-metrics-{date}.zip"
    body = io.StringIO()
    writer = csv.writer(body, lineterminator="\n")
    writer.writerow(header)
    for stamp, row_symbol, value in rows:
        writer.writerow([stamp, row_symbol, value] + ["0"] * (len(header) - 3))
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(member, body.getvalue())
    return {"symbol": symbol, "date": date, "path": str(path), "zip_member": member}


def full_day(symbol: str, date: str, value: str = "1000", count: int = 288) -> list[tuple[str, str, str]]:
    rows = []
    for index in range(count):
        minute = index * 5
        rows.append((f"{date} {minute // 60:02d}:{minute % 60:02d}:00", symbol, value))
    return rows


class BoundaryRowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-oi6-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)

    def parse(self, rows, symbol="AAAUSDT", date="2024-03-05", **kwargs):
        entry = write_metrics_zip(self.directory, symbol, date, rows, **kwargs)
        census = fresh_census()
        return runner._parse_metrics_archive(entry, census, root=Path("/")), census

    def test_real_289_row_archive_parses_and_counts_one_boundary_row(self) -> None:
        symbol, date = REAL_BOUNDARY_ARCHIVE
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        entry = next(e for e in inventory["files"]
                     if (e["symbol"], e["date"]) == REAL_BOUNDARY_ARCHIVE)
        with zipfile.ZipFile(runner.WORKSPACE_ROOT / entry["path"]) as archive:
            with archive.open(entry["zip_member"]) as raw:
                rows = list(csv.reader(io.TextIOWrapper(raw, encoding="utf-8")))
        self.assertEqual(len(rows), 289, "fixture drift: expected a header plus 288 data rows")
        off_date = [row for row in rows[1:] if row[0][:10] != date]
        self.assertEqual(len(off_date), 1)
        self.assertEqual(off_date[0][0], "2024-04-04 00:00:00")
        self.assertEqual(len(rows) - 1 - len(off_date), 287, "287 in-date rows survive the drop")

        census = fresh_census()
        value = runner._parse_metrics_archive(entry, census)
        self.assertEqual(census["off-date-boundary-row"], 1)
        self.assertEqual(census["off-date-offset-days:+1"], 1)
        self.assertEqual(census["archives-with-off-date-row"], 1)
        self.assertEqual(census["off-date-row-at-or-after-cutoff"], 0)

        in_date = [row for row in rows[1:] if row[0][:10] == date]
        expected = float(in_date[-1][2])
        self.assertEqual(value, expected)
        self.assertGreater(value, 0.0)

    def test_cycle_five_runner_still_aborts_on_the_same_real_archive(self) -> None:
        """The correction is a real behavioral change, not a no-op."""
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        entry = next(e for e in inventory["files"]
                     if (e["symbol"], e["date"]) == REAL_BOUNDARY_ARCHIVE)
        with self.assertRaisesRegex(ValueError, "create_time outside archive date"):
            oi5._parse_metrics_archive(entry)

    def test_off_date_row_is_dropped_and_never_selected(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05", value="900")
        rows.append(("2024-03-06 00:00:00", "AAAUSDT", "5000"))
        value, census = self.parse(rows)
        self.assertEqual(value, 900.0)
        self.assertEqual(census["off-date-boundary-row"], 1)
        self.assertEqual(census["off-date-offset-days:+1"], 1)
        self.assertEqual(census["off-date-row-at-or-after-cutoff"], 0)

    def test_off_date_row_after_cutoff_is_counted_and_still_not_selected(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05", value="900")
        rows.append(("2024-03-06 23:50:00", "AAAUSDT", "5000"))
        value, census = self.parse(rows)
        self.assertEqual(value, 900.0)
        self.assertEqual(census["off-date-row-at-or-after-cutoff"], 1)
        self.assertEqual(census["off-date-boundary-row"], 1)

    def test_previous_day_boundary_row_is_counted_with_a_negative_offset(self) -> None:
        rows = [("2024-03-04 23:55:00", "AAAUSDT", "7777")] + full_day("AAAUSDT", "2024-03-05", value="900")
        value, census = self.parse(rows)
        self.assertEqual(value, 900.0)
        self.assertEqual(census["off-date-offset-days:-1"], 1)
        self.assertEqual(census["off-date-row-at-or-after-cutoff"], 1)

    def test_archive_of_entirely_off_date_rows_yields_no_observation(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-06")
        entry = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-05", rows)
        census = fresh_census()
        self.assertIsNone(runner._parse_metrics_archive(entry, census, root=Path("/")))
        self.assertEqual(census["off-date-boundary-row"], 288)
        self.assertEqual(census["archives-with-off-date-row"], 1)

    def test_archive_without_a_row_after_the_cutoff_is_unchanged_behavior(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05", count=200)
        value, census = self.parse(rows)
        self.assertIsNone(value)
        self.assertEqual(census["off-date-boundary-row"], 0)

    def test_nonpositive_end_of_day_value_is_missing_not_zero(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05")
        rows[-1] = ("2024-03-05 23:55:00", "AAAUSDT", "0")
        value, _ = self.parse(rows)
        self.assertIsNone(value)

    def test_last_in_date_row_after_cutoff_wins_over_earlier_rows(self) -> None:
        rows = [("2024-03-05 23:45:00", "AAAUSDT", "111"), ("2024-03-05 23:50:00", "AAAUSDT", "222"),
                ("2024-03-05 23:55:00", "AAAUSDT", "333"), ("2024-03-06 00:00:00", "AAAUSDT", "444")]
        value, census = self.parse(rows)
        self.assertEqual(value, 333.0)
        self.assertEqual(census["off-date-boundary-row"], 1)

    def test_malformed_stamp_still_aborts(self) -> None:
        for stamp in ("2024-03-05T23:55:00", "2024-03-05 23:55", "not-a-time", "2024-3-05 23:55:00"):
            rows = full_day("AAAUSDT", "2024-03-05") + [(stamp, "AAAUSDT", "1")]
            with self.assertRaisesRegex(ValueError, "schema value"):
                self.parse(rows)

    def test_out_of_range_stamp_fields_still_abort(self) -> None:
        for stamp in ("2024-03-05 24:00:00", "2024-03-05 23:60:00", "2024-03-05 23:55:60"):
            rows = full_day("AAAUSDT", "2024-03-05") + [(stamp, "AAAUSDT", "1")]
            with self.assertRaisesRegex(ValueError, "schema value"):
                self.parse(rows)

    def test_impossible_calendar_date_aborts_rather_than_counting_as_off_date(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05") + [("2024-02-30 00:00:00", "AAAUSDT", "1")]
        with self.assertRaisesRegex(ValueError, "schema value"):
            self.parse(rows)

    def test_symbol_column_mismatch_still_aborts(self) -> None:
        rows = full_day("AAAUSDT", "2024-03-05") + [("2024-03-05 23:59:00", "BBBUSDT", "1")]
        with self.assertRaisesRegex(ValueError, "symbol column mismatch"):
            self.parse(rows)

    def test_column_count_drift_still_aborts(self) -> None:
        entry = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-05", full_day("AAAUSDT", "2024-03-05"),
                                  header=runner.METRICS_HEADER + ("extra",))
        with self.assertRaisesRegex(ValueError, "metrics header drift"):
            runner._parse_metrics_archive(entry, fresh_census(), root=Path("/"))

    def test_header_drift_still_aborts(self) -> None:
        bad = ("created_at",) + runner.METRICS_HEADER[1:]
        entry = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-05",
                                  full_day("AAAUSDT", "2024-03-05"), header=bad)
        with self.assertRaisesRegex(ValueError, "metrics header drift"):
            runner._parse_metrics_archive(entry, fresh_census(), root=Path("/"))

    def test_parse_stamp_returns_date_and_time_of_day(self) -> None:
        self.assertEqual(runner._parse_stamp("2024-03-05 23:45:00"), ("2024-03-05", CUTOFF))
        self.assertEqual(runner._parse_stamp("2024-03-05 00:00:00"), ("2024-03-05", 0))


class LoaderCensusTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-oi6-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)

    def test_loader_aggregates_the_census_and_counts_missing_days(self) -> None:
        good = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-05",
                                 full_day("AAAUSDT", "2024-03-05", value="900")
                                 + [("2024-03-06 00:00:00", "AAAUSDT", "1")])
        short = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-06",
                                  full_day("AAAUSDT", "2024-03-06", count=10))
        inventory = {"files": [good, short]}
        counters = {key: 0 for key in runner.ELIGIBILITY_KEYS}
        census = fresh_census()
        with mock.patch.object(runner, "WORKSPACE_ROOT", Path("/")):
            values = runner.load_end_of_day_open_interest(
                inventory, ["AAAUSDT"], runner.parse_iso_ms("2024-03-01T00:00:00Z"),
                runner.parse_iso_ms("2024-03-10T00:00:00Z"), counters, census)
        self.assertEqual(list(values.values()), [900.0])
        self.assertEqual(census["archives-read"], 2)
        self.assertEqual(census["off-date-boundary-row"], 1)
        self.assertEqual(census["archives-with-off-date-row"], 1)
        self.assertEqual(counters["archive-without-end-of-day-row"], 1)

    def test_loader_rejects_a_duplicate_symbol_date(self) -> None:
        entry = write_metrics_zip(self.directory, "AAAUSDT", "2024-03-05", full_day("AAAUSDT", "2024-03-05"))
        inventory = {"files": [entry, dict(entry)]}
        with mock.patch.object(runner, "WORKSPACE_ROOT", Path("/")):
            with self.assertRaisesRegex(ValueError, "duplicate inventory date"):
                runner.load_end_of_day_open_interest(
                    inventory, ["AAAUSDT"], runner.parse_iso_ms("2024-03-01T00:00:00Z"),
                    runner.parse_iso_ms("2024-03-10T00:00:00Z"),
                    {key: 0 for key in runner.ELIGIBILITY_KEYS}, fresh_census())

    def test_loader_skips_symbols_outside_the_universe_and_dates_outside_the_span(self) -> None:
        other = write_metrics_zip(self.directory, "BBBUSDT", "2024-03-05", full_day("BBBUSDT", "2024-03-05"))
        late = write_metrics_zip(self.directory, "AAAUSDT", "2024-04-05", full_day("AAAUSDT", "2024-04-05"))
        census = fresh_census()
        with mock.patch.object(runner, "WORKSPACE_ROOT", Path("/")):
            values = runner.load_end_of_day_open_interest(
                {"files": [other, late]}, ["AAAUSDT"], runner.parse_iso_ms("2024-03-01T00:00:00Z"),
                runner.parse_iso_ms("2024-03-10T00:00:00Z"),
                {key: 0 for key in runner.ELIGIBILITY_KEYS}, census)
        self.assertEqual(values, {})
        self.assertEqual(census["archives-read"], 0)


class ContractBindingTests(unittest.TestCase):
    def test_frozen_contract_binds_runner_parameters(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["frozen_rule"]["parameters"], runner.FROZEN_PARAMETERS)
        runner.check_contract_parameters(contract)

    def test_only_loop_alpha_differs_from_the_cycle_five_parameter_set(self) -> None:
        drifted = {**runner.FROZEN_PARAMETERS, "universe_size": 71}
        with mock.patch.object(runner, "FROZEN_PARAMETERS", drifted):
            with self.assertRaisesRegex(ValueError, "only loop_alpha may differ"):
                runner.check_contract_parameters({"frozen_rule": {"parameters": drifted}})

    def test_parameter_drift_aborts(self) -> None:
        with self.assertRaisesRegex(ValueError, "contract parameters differ"):
            runner.check_contract_parameters(
                {"frozen_rule": {"parameters": {**runner.FROZEN_PARAMETERS, "taker_fee_per_side": 0.0}}})

    def test_loop_alpha_is_the_sixth_ordinal_and_stricter_than_cycle_five(self) -> None:
        self.assertEqual(runner.LOOP_ALPHA, 0.05 / (6 * 7))
        self.assertLess(runner.LOOP_ALPHA, oi5.LOOP_ALPHA)
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["loop_multiplicity"]["alpha_k"], runner.LOOP_ALPHA)
        self.assertEqual(contract["loop_multiplicity"]["cycle_ordinal_k"], 6)

    def test_reused_rule_sections_are_byte_identical_to_cycle_five(self) -> None:
        c6 = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        c5 = json.loads(oi5.CONTRACT_PATH.read_text(encoding="utf-8"))
        for section in ("universe", "daily_bars", "signal", "eligibility", "portfolio",
                        "execution_and_costs"):
            self.assertEqual(c6["frozen_rule"][section], c5["frozen_rule"][section], section)
        self.assertEqual(c6["windows"], c5["windows"])
        self.assertEqual(c6["minimum_sample"], c5["minimum_sample"])
        self.assertNotEqual(c6["frozen_rule"]["open_interest"], c5["frozen_rule"]["open_interest"])

    def test_contract_pins_the_frozen_cycle_five_runner_and_parent(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["governance"]["oi5_runner_sha256"],
                         runner.sha256_file(runner.OI5_PATH))
        self.assertEqual(contract["governance"]["parent_runner_sha256"],
                         runner.sha256_file(runner.PARENT_PATH))

    def test_static_input_hash_drift_aborts_before_any_archive_check(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        contract["assets_and_data"]["input_inventory_sha256"] = "0" * 64
        with mock.patch.object(runner, "verify_metrics_archive") as verify:
            with self.assertRaisesRegex(ValueError, "input_inventory_sha256 mismatch"):
                runner.verify_static_inputs(contract)
        verify.assert_not_called()

    def test_inventory_must_certify_pre_outcome_state_and_universes(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        good = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        directory = Path(tempfile.mkdtemp(prefix=".tmp-oi6-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, directory, ignore_errors=True)
        pins = {runner.INVENTORY_PATH: contract["assets_and_data"]["input_inventory_sha256"],
                runner.KLINE_INVENTORY_PATH: contract["assets_and_data"]["kline_inventory_sha256"],
                runner.AUDIT_PATH: contract["novelty_and_consumption"]["prior_window_audit_sha256"],
                runner.POLICY_PATH: contract["loop_multiplicity"]["policy_sha256"],
                runner.PARENT_PATH: contract["governance"]["parent_runner_sha256"],
                runner.OI5_PATH: contract["governance"]["oi5_runner_sha256"]}
        for mutate, message in (
                (lambda inv: inv.update(candidate_outcomes_accessed=True), "pre-outcome"),
                (lambda inv: inv.update(symbols=list(reversed(inv["symbols"]))), "archive count mismatch"),
                (lambda inv: inv["windows"].pop("replication"), "primary and replication universes"),
                (lambda inv: inv["windows"]["primary"]["universe"].pop(), "72 distinct symbols")):
            inventory = json.loads(json.dumps(good))
            mutate(inventory)
            if "archive count mismatch" not in message and "duplicate" not in message:
                inventory["archive_count"] = len(inventory["files"])
            fake = directory / "inventory.json"
            fake.write_text(json.dumps(inventory), encoding="utf-8")
            with mock.patch.object(runner, "INVENTORY_PATH", fake), \
                 mock.patch.object(runner, "sha256_file",
                                   side_effect=lambda p, pin=pins[runner.INVENTORY_PATH]: pins.get(p, pin)), \
                 mock.patch.object(runner, "verify_metrics_archive"), \
                 mock.patch.object(parent, "verify_archive"):
                with self.assertRaisesRegex(ValueError, message):
                    runner.verify_static_inputs(contract)


class GateTests(unittest.TestCase):
    @staticmethod
    def rows(spread: float, dates: int = 40, excess: float = 0.001) -> list[dict]:
        out = []
        for index in range(dates):
            day = f"2024-{1 + index // 10:02d}-{1 + index % 10:02d}"
            out.append({"date": day, "quarter": runner.oi5.quarter_label(day), "spread_net": spread,
                        "spread_gross": spread + 0.001, "long_excess": excess, "short_excess": excess})
        return out

    def test_gates_use_the_sixth_ordinal_alpha(self) -> None:
        gates = runner.evaluate_gates(self.rows(0.01))
        self.assertEqual(gates["G4_dependence_robustness"]["alpha"], runner.LOOP_ALPHA)
        self.assertTrue(all(gate["pass"] for gate in gates.values()))

    def test_gate_four_is_stricter_than_cycle_five(self) -> None:
        rows = self.rows(0.01)
        p = runner.evaluate_gates(rows)["G4_dependence_robustness"]["one_sided_block_bootstrap_p_spread_net"]
        self.assertEqual(p, oi5.evaluate_gates(rows)["G4_dependence_robustness"]
                         ["one_sided_block_bootstrap_p_spread_net"])
        with mock.patch.object(oi5, "block_bootstrap_p", return_value=0.0015):
            self.assertTrue(oi5.evaluate_gates(rows)["G4_dependence_robustness"]["pass"])
            self.assertFalse(runner.evaluate_gates(rows)["G4_dependence_robustness"]["pass"])

    def test_negative_spread_fails_magnitude_and_leg_gates(self) -> None:
        gates = runner.evaluate_gates(self.rows(-0.01, excess=-0.001))
        self.assertFalse(gates["G1_net_magnitude"]["pass"])
        self.assertFalse(gates["G2_leg_contribution"]["pass"])

    def test_one_negative_leg_fails_the_leg_gate(self) -> None:
        rows = self.rows(0.01)
        for row in rows:
            row["short_excess"] = -0.001
        self.assertFalse(runner.evaluate_gates(rows)["G2_leg_contribution"]["pass"])


class WindowTests(unittest.TestCase):
    def contract(self) -> dict:
        return {"id": "test", "frozen_rule": {"parameters": dict(runner.FROZEN_PARAMETERS)},
                "windows": {"primary": {"start_utc_inclusive": "2024-01-01T00:00:00Z",
                                        "end_utc_exclusive": "2025-01-01T00:00:00Z"}},
                "claim_boundary": "test boundary"}

    def test_universe_drift_from_the_frozen_inventory_aborts(self) -> None:
        contract = self.contract()
        inventory = {"windows": {"primary": {"start_utc_inclusive": "2024-01-01T00:00:00Z",
                                             "end_utc_exclusive": "2025-01-01T00:00:00Z",
                                             "universe": ["ZZZUSDT"]}}}
        with mock.patch.object(runner, "window_universe", return_value=(["AAAUSDT"], {})):
            with self.assertRaisesRegex(ValueError, "universe recomputed from frozen klines differs"):
                runner.run_window(contract, "primary", inventory, {})

    def test_window_drift_between_contract_and_inventory_aborts(self) -> None:
        inventory = {"windows": {"primary": {"start_utc_inclusive": "2023-01-01T00:00:00Z",
                                             "end_utc_exclusive": "2024-01-01T00:00:00Z",
                                             "universe": []}}}
        with self.assertRaisesRegex(ValueError, "window differs between contract and input inventory"):
            runner.run_window(self.contract(), "primary", inventory, {})

    def test_small_sample_is_inconclusive_without_gates_and_reports_the_census(self) -> None:
        contract = self.contract()
        symbols = [f"S{index:03d}USDT" for index in range(40)]
        inventory = {"windows": {"primary": {"start_utc_inclusive": "2024-01-01T00:00:00Z",
                                             "end_utc_exclusive": "2025-01-01T00:00:00Z",
                                             "universe": symbols}}}
        start = runner.parse_iso_ms("2024-01-01T00:00:00Z")
        records = {(symbol, start + offset * DAY_MS): day_record(100.0, 101.0)
                   for symbol in symbols for offset in range(-runner.ELIGIBILITY_DAYS, 40)}
        open_interest = {(symbol, start + offset * DAY_MS): 1000.0 + index
                         for index, symbol in enumerate(symbols) for offset in range(-2, 40)}
        with mock.patch.object(runner, "window_universe", return_value=(symbols, {})), \
             mock.patch.object(parent, "load_day_records", return_value=records), \
             mock.patch.object(runner, "load_end_of_day_open_interest", return_value=open_interest):
            result = runner.run_window(contract, "primary", inventory, {})
        self.assertEqual(result["status"], "inconclusive")
        self.assertFalse(result["sample_minimum_met"])
        self.assertEqual(result["gate_results"], {})
        self.assertIn("archive_census", result["diagnostics"])
        self.assertIn("off-date-boundary-row", result["diagnostics"]["archive_census"])


class LifecycleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-oi6-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)
        self.contract_path = self.directory / "contract.json"
        self.contract_path.write_text(json.dumps({"id": "test"}), encoding="utf-8")
        self.contract_sha = runner.sha256_file(self.contract_path)
        self.verified = {"input_inventory_sha256": "i" * 64, "kline_inventory_sha256": "k" * 64,
                         "prior_window_audit_sha256": "a" * 64, "loop_policy_sha256": "p" * 64,
                         "parent_runner_sha256": "r" * 64, "oi5_runner_sha256": "o" * 64,
                         "metrics_archives_verified": 51224, "kline_archives_verified": 14535}

    def fake_result(self, label: str, status: str = "falsified") -> dict:
        return {"role": label, "status": status, "claim_boundary": "test",
                "sample": {"included_dates": 0, "positions": 0}, "dates": []}

    def run_main(self, args, status: str = "falsified"):
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

    def test_output_outside_the_workspace_refuses(self) -> None:
        with self.assertRaisesRegex(ValueError, "escapes workspace root"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                           "--output", "/tmp/oi6-never-written.json"])

    def test_verify_inputs_only_prints_every_pin_and_writes_nothing(self) -> None:
        output = self.directory / "unused.json"
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                                            "--output", str(output), "--verify-inputs-only"]), 0)
        payload = json.loads(printed.getvalue())
        for field in runner.PIN_FIELDS:
            self.assertEqual(payload[field], self.verified[field])
        self.assertEqual(payload["runner_sha256"], runner.sha256_file(Path(runner.__file__).resolve()))
        self.assertFalse(output.exists())

    def test_primary_then_replication_end_to_end(self) -> None:
        primary = self.directory / "primary.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                                        "--output", str(primary)], status="pass-on-this-sample"), 0)
        replication = self.directory / "replication.json"
        self.assertEqual(self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                                        "--output", str(replication), "--primary-artifact", str(primary),
                                        "--primary-artifact-sha256", runner.sha256_file(primary)]), 0)
        payload = json.loads(replication.read_text(encoding="utf-8"))
        self.assertEqual(payload["role"], "replication")
        self.assertEqual(payload["primary_artifact_sha256"], runner.sha256_file(primary))
        self.assertEqual(payload["runner_sha256"], runner.sha256_file(Path(runner.__file__).resolve()))

    def test_replication_refuses_a_primary_from_a_different_runner(self) -> None:
        primary = self.directory / "primary.json"
        self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                       "--output", str(primary)], status="pass-on-this-sample")
        payload = json.loads(primary.read_text(encoding="utf-8"))
        payload["runner_sha256"] = "9" * 64
        edited = self.directory / "edited.json"
        edited.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "runner_sha256 mismatch"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                           "--output", str(self.directory / "r2.json"), "--primary-artifact", str(edited),
                           "--primary-artifact-sha256", runner.sha256_file(edited)])

    def test_replication_refuses_drifted_static_input_pins(self) -> None:
        primary = self.directory / "primary.json"
        self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                       "--output", str(primary)], status="pass-on-this-sample")
        self.verified = {**self.verified, "oi5_runner_sha256": "1" * 64}
        with self.assertRaisesRegex(SystemExit, "static-input pins differ"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                           "--output", str(self.directory / "r3.json"), "--primary-artifact", str(primary),
                           "--primary-artifact-sha256", runner.sha256_file(primary)])

    def test_replication_requires_the_primary_artifact_arguments(self) -> None:
        with self.assertRaisesRegex(SystemExit, "replication requires"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                           "--output", str(self.directory / "r4.json")])


if __name__ == "__main__":
    unittest.main()
