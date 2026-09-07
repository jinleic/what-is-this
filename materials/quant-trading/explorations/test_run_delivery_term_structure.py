#!/usr/bin/env python3
"""Behavioral tests for the frozen delivery-futures calendar-spread runner (cycle 7).

The load-bearing tests are: the pair universe rebuilt from contract names alone, the calendar-derived
entry/exit search, the exact inverse-contract arithmetic and its sign convention, the four gates, and
the write-once / replication-binding refusals. One test binds the real frozen inventory.
"""
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
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from unittest import mock

import run_delivery_term_structure as runner

HOUR_MS = runner.HOUR_MS
DAY_MS = runner.DAY_MS
FEE = runner.TAKER_FEE


def day_ms(day: str) -> int:
    return int(datetime.fromisoformat(day + "T00:00:00+00:00").timestamp() * 1000)


def full_day(day: str, close: float = 100.0, trades: int = 5, hours: int = 24) -> list[list[str]]:
    """A complete 24-bar UTC day; the 23:00 bar carries `close` and `trades`."""
    base = day_ms(day)
    rows = []
    for hour in range(hours):
        is_close = hour == runner.CLOSE_HOUR
        rows.append([str(base + hour * HOUR_MS), "1", "1", "1",
                     f"{close}" if is_close else "1", "1", str(base + hour * HOUR_MS + HOUR_MS - 1),
                     "1", str(trades) if is_close else "1", "1", "1", "0"])
    return rows


def write_archive(directory: Path, symbol: str, month: str, rows: list[list[str]],
                  header: bool = True) -> dict:
    name = f"{symbol}-1h-{month}.csv"
    body = ([runner.HEADER] if header else []) + [",".join(r) for r in rows]
    payload = ("\n".join(body) + "\n").encode()
    zip_path = directory / f"{name}.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(name, payload)
    sidecar = directory / f"{name}.zip.CHECKSUM"
    digest = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    sidecar.write_text(f"{digest}  {zip_path.name}\n", encoding="utf-8")
    return {"margin": "cm", "symbol": symbol, "month": month,
            "path": str(zip_path.relative_to(runner.WORKSPACE_ROOT)),
            "zip_member": name, "sha256": digest, "size_bytes": zip_path.stat().st_size,
            "checksum_sidecar": str(sidecar.relative_to(runner.WORKSPACE_ROOT)),
            "checksum_sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            "header_present": header, "rows": len(rows)}


class PairUniverseTests(unittest.TestCase):
    def test_universe_is_rebuilt_from_contract_names_alone(self) -> None:
        symbols = ["cm:AAAUSD_240329", "cm:AAAUSD_240628", "cm:AAAUSD_240927", "cm:AAAUSD_241227"]
        derived = runner.derive_pair_universe(symbols)
        self.assertEqual([(p["near"], p["far"], p["prev_delivery"]) for p in derived],
                         [("AAAUSD_240628", "AAAUSD_240927", "2024-03-29"),
                          ("AAAUSD_240927", "AAAUSD_241227", "2024-06-28")])

    def test_first_and_last_contract_are_never_a_near_leg(self) -> None:
        derived = runner.derive_pair_universe(["cm:AAAUSD_240329", "cm:AAAUSD_240628", "cm:AAAUSD_240927"])
        self.assertEqual([p["near"] for p in derived], ["AAAUSD_240628"])

    def test_a_missing_quarter_breaks_the_triple(self) -> None:
        """A >100-day delivery gap is not a hand-over, so no candidate spans it."""
        symbols = ["cm:AAAUSD_240329", "cm:AAAUSD_240628", "cm:AAAUSD_241227", "cm:AAAUSD_250328"]
        derived = runner.derive_pair_universe(symbols)
        self.assertEqual([p["near"] for p in derived], [])

    def test_usd_margined_contracts_are_excluded_by_rule(self) -> None:
        symbols = ["um:BTCUSDT_240329", "um:BTCUSDT_240628", "um:BTCUSDT_240927"]
        self.assertEqual(runner.derive_pair_universe(symbols), [])

    def test_non_quarterly_and_malformed_stamps_abort(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-quarterly"):
            runner.delivery_date("BTCBUSD_210129")
        for bad in ("BTCUSD_2103", "BTCUSD_21032a", "BTCUSD_2103２6"):
            with self.assertRaisesRegex(ValueError, "malformed delivery stamp"):
                runner.delivery_date(bad)

    def test_real_inventory_pair_universe_matches_the_pinned_list(self) -> None:
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        self.assertEqual(runner.derive_pair_universe(inventory["symbols"]), inventory["pair_universe"])
        self.assertEqual(len(inventory["pair_universe"]), inventory["pair_universe_count"])


class TradableDayTests(unittest.TestCase):
    def table(self, **days) -> dict:
        return {d: {"hours": set(range(h)), "close": c, "trades": t} for d, (h, c, t) in days.items()}

    def test_day_selection_is_structural_and_never_looks_at_the_price_value(self) -> None:
        """Bar presence, a real closing trade, and an existing 23:00 bar select days.

        A zero or non-finite close does NOT remove a day, so a price value can never shift entry
        or exit; position_return rejects such a candidate outright instead.
        """
        table = self.table(**{"2024-01-01": (24, 100.0, 3), "2024-01-02": (23, 100.0, 3),
                              "2024-01-03": (24, 100.0, 0), "2024-01-04": (24, 0.0, 3),
                              "2024-01-05": (24, float("nan"), 3), "2024-01-06": (24, None, 3)})
        self.assertEqual(runner.tradable_days(table),
                         {"2024-01-01", "2024-01-04", "2024-01-05"})
        self.assertEqual(runner.position_return({"near_entry": 100.0, "near_exit": 0.0,
                                                 "far_entry": 101.0, "far_exit": 105.0}),
                         {"rejected": "nonpositive-or-nonfinite-close"})


class LocatePositionTests(unittest.TestCase):
    CANDIDATE = {"margin": "cm", "base": "AAAUSD", "near": "AAAUSD_240628", "far": "AAAUSD_240927",
                 "prev_delivery": "2024-03-29", "near_delivery": "2024-06-28", "quarter": "2024Q2"}

    def days_for(self, near_days, far_days) -> dict:
        def table(spec):
            return {d: {"hours": set(range(24)), "close": c, "trades": 5} for d, c in spec.items()}
        return {("cm", "AAAUSD_240628"): table(near_days), ("cm", "AAAUSD_240927"): table(far_days)}

    def test_entry_is_the_first_common_tradable_day_after_the_prior_delivery(self) -> None:
        near = {"2024-03-29": 100.0, "2024-03-30": 100.0, "2024-06-27": 110.0}
        far = {"2024-03-30": 101.0, "2024-06-27": 110.5}
        located = runner.locate_position(self.CANDIDATE, self.days_for(near, far))
        self.assertEqual((located["entry_day"], located["exit_day"]), ("2024-03-30", "2024-06-27"))
        self.assertEqual(located["entry_lag_days"], 1)

    def test_prior_delivery_day_itself_is_never_the_entry(self) -> None:
        near = {"2024-03-29": 100.0, "2024-04-05": 100.0, "2024-06-27": 110.0}
        far = {"2024-03-29": 101.0, "2024-04-05": 101.0, "2024-06-27": 110.5}
        located = runner.locate_position(self.CANDIDATE, self.days_for(near, far))
        self.assertEqual(located["entry_day"], "2024-04-05")

    def test_delivery_day_itself_is_never_the_exit(self) -> None:
        near = {"2024-03-30": 100.0, "2024-06-27": 110.0, "2024-06-28": 111.0}
        far = {"2024-03-30": 101.0, "2024-06-27": 110.5, "2024-06-28": 111.5}
        located = runner.locate_position(self.CANDIDATE, self.days_for(near, far))
        self.assertEqual(located["exit_day"], "2024-06-27")

    def test_entry_after_the_interval_midpoint_is_rejected_as_late_handover(self) -> None:
        """Guards the degenerate short hold: a late-listing far leg is not an instance."""
        near = {"2024-06-01": 100.0, "2024-06-27": 110.0}
        far = {"2024-06-01": 101.0, "2024-06-27": 110.5}
        self.assertEqual(runner.locate_position(self.CANDIDATE, self.days_for(near, far)),
                         {"rejected": "late-handover-entry"})

    def test_midpoint_boundary_is_inclusive(self) -> None:
        midpoint = date(2024, 3, 29) + (date(2024, 6, 28) - date(2024, 3, 29)) / 2
        near = {midpoint.isoformat(): 100.0, "2024-06-27": 110.0}
        far = {midpoint.isoformat(): 101.0, "2024-06-27": 110.5}
        located = runner.locate_position(self.CANDIDATE, self.days_for(near, far))
        self.assertEqual(located["entry_day"], midpoint.isoformat())

    def test_a_leg_that_is_never_jointly_tradable_is_rejected(self) -> None:
        near = {"2024-03-30": 100.0, "2024-06-27": 110.0}
        self.assertEqual(runner.locate_position(self.CANDIDATE, self.days_for(near, {})),
                         {"rejected": "no-complete-entry-day"})

    def test_single_common_day_cannot_be_both_entry_and_exit(self) -> None:
        near = {"2024-03-30": 100.0}
        far = {"2024-03-30": 101.0}
        self.assertEqual(runner.locate_position(self.CANDIDATE, self.days_for(near, far)),
                         {"rejected": "entry-not-before-exit"})


class ArithmeticTests(unittest.TestCase):
    def test_sign_convention_shorts_the_rich_leg(self) -> None:
        rich_far = runner.position_return({"near_entry": 100.0, "near_exit": 105.0,
                                           "far_entry": 101.0, "far_exit": 105.0})
        self.assertEqual(rich_far["side"], "long-near")
        rich_near = runner.position_return({"near_entry": 101.0, "near_exit": 105.0,
                                            "far_entry": 100.0, "far_exit": 105.0})
        self.assertEqual(rich_near["side"], "long-far")

    def test_convergence_pays_and_divergence_loses(self) -> None:
        converge = runner.position_return({"near_entry": 100.0, "near_exit": 105.0,
                                           "far_entry": 102.0, "far_exit": 105.0})
        diverge = runner.position_return({"near_entry": 100.0, "near_exit": 105.0,
                                          "far_entry": 102.0, "far_exit": 110.0})
        self.assertGreater(converge["gross_return"], 0.0)
        self.assertLess(diverge["gross_return"], 0.0)

    def test_exact_inverse_arithmetic_matches_a_hand_computation(self) -> None:
        pn_e, pn_x, pf_e, pf_x = 100.0, 105.0, 102.0, 105.0
        out = runner.position_return({"near_entry": pn_e, "near_exit": pn_x,
                                      "far_entry": pf_e, "far_exit": pf_x})
        pnl = (1 / pn_e - 1 / pn_x) - (1 / pf_e - 1 / pf_x)
        fees = FEE * (1 / pn_e + 1 / pf_e + 1 / pn_x + 1 / pf_x)
        gross = 1 / pn_e + 1 / pf_e
        self.assertAlmostEqual(out["gross_return"], pnl / gross, places=15)
        self.assertAlmostEqual(out["net_return"], (pnl - fees) / gross, places=15)
        self.assertAlmostEqual(out["fee_drag"], fees / gross, places=15)

    def test_fee_drag_is_about_two_taker_legs_of_gross(self) -> None:
        out = runner.position_return({"near_entry": 100.0, "near_exit": 100.0,
                                      "far_entry": 101.0, "far_exit": 100.5})
        self.assertAlmostEqual(out["fee_drag"], 2 * FEE, delta=2e-5)

    def test_a_flat_spread_still_pays_the_fee(self) -> None:
        out = runner.position_return({"near_entry": 100.0, "near_exit": 110.0,
                                      "far_entry": 101.0, "far_exit": 111.1})
        self.assertLess(out["net_return"], 0.0)
        self.assertAlmostEqual(out["net_return"], out["gross_return"] - out["fee_drag"], places=15)

    def test_zero_spread_and_bad_prices_are_rejected(self) -> None:
        self.assertEqual(runner.position_return({"near_entry": 100.0, "near_exit": 105.0,
                                                 "far_entry": 100.0, "far_exit": 105.0}),
                         {"rejected": "zero-spread-at-entry"})
        for bad in ({"near_entry": 0.0}, {"near_exit": -1.0}, {"far_exit": float("inf")},
                    {"far_entry": float("nan")}):
            prices = {"near_entry": 100.0, "near_exit": 105.0, "far_entry": 101.0, "far_exit": 105.0}
            prices.update(bad)
            self.assertEqual(runner.position_return(prices),
                             {"rejected": "nonpositive-or-nonfinite-close"})

    def test_multiplier_cancels_so_the_return_is_scale_free(self) -> None:
        """Doubling both legs' price level must not change the spread return."""
        a = runner.position_return({"near_entry": 100.0, "near_exit": 105.0,
                                    "far_entry": 101.0, "far_exit": 105.0})
        b = runner.position_return({"near_entry": 200.0, "near_exit": 210.0,
                                    "far_entry": 202.0, "far_exit": 210.0})
        self.assertAlmostEqual(a["net_return"], b["net_return"], places=12)


class BootstrapTests(unittest.TestCase):
    def test_is_deterministic_and_seed_bound(self) -> None:
        """Heterogeneous quarters: identical quarters would make every resample equal and hide the seed."""
        data = {f"20{20 + i // 4}Q{i % 4 + 1}": [0.01 * (1 if i % 3 else -2.5)] for i in range(13)}
        first = runner.block_bootstrap_p(data)
        self.assertEqual(first, runner.block_bootstrap_p(data))
        self.assertNotEqual(first, runner.block_bootstrap_p(data, seed=1))
        self.assertGreater(first, 0.0)

    def test_too_few_quarters_for_one_block_returns_none(self) -> None:
        self.assertIsNone(runner.block_bootstrap_p({"2024Q1": [0.01], "2024Q2": [0.01], "2024Q3": [0.01]}))

    def test_all_positive_quarters_give_a_small_p_and_all_negative_a_large_one(self) -> None:
        good = {f"2024Q{i}": [0.05] for i in range(1, 5)} | {f"2025Q{i}": [0.05] for i in range(1, 5)}
        bad = {q: [-v[0]] for q, v in good.items()}
        self.assertLess(runner.block_bootstrap_p(good, replicates=2000), runner.LOOP_ALPHA)
        self.assertGreater(runner.block_bootstrap_p(bad, replicates=2000), 0.9)

    def test_p_is_never_exactly_zero(self) -> None:
        good = {f"2024Q{i}": [1.0] for i in range(1, 5)} | {f"2025Q{i}": [1.0] for i in range(1, 5)}
        self.assertGreater(runner.block_bootstrap_p(good, replicates=100), 0.0)


class GateTests(unittest.TestCase):
    def rows(self, net: float, quarters: int = 13, bases=("BTCUSD", "ETHUSD", "ADAUSD", "XRPUSD")) -> list[dict]:
        out = []
        for index in range(quarters):
            for base in bases:
                out.append({"base": base, "quarter": f"20{20 + index // 4}Q{index % 4 + 1}",
                            "net_return": net, "gross_return": net + 2 * FEE, "fee_drag": 2 * FEE,
                            "side": "long-near", "spread_entry": 0.01, "spread_exit": 0.0})
        return out

    def test_a_strong_consistent_effect_passes_every_gate(self) -> None:
        gates = runner.evaluate_gates(self.rows(0.004))
        self.assertTrue(all(gate["pass"] for gate in gates.values()), gates)

    def test_magnitude_below_the_threshold_fails_g1(self) -> None:
        gates = runner.evaluate_gates(self.rows(runner.MAGNITUDE_THRESHOLD / 2))
        self.assertFalse(gates["G1_net_magnitude"]["pass"])

    def test_threshold_boundary_is_inclusive(self) -> None:
        gates = runner.evaluate_gates(self.rows(runner.MAGNITUDE_THRESHOLD))
        self.assertTrue(gates["G1_net_magnitude"]["pass"])

    def test_one_negative_segment_fails_g2_even_when_the_pool_is_positive(self) -> None:
        rows = self.rows(0.004)
        for row in rows:
            if row["base"] in runner.MAJORS:
                row["net_return"] = -0.001
        gates = runner.evaluate_gates(rows)
        self.assertFalse(gates["G2_segment_consistency"]["pass"])
        self.assertGreater(gates["G1_net_magnitude"]["mean_net_return"], 0.0)

    def test_g2_requires_both_segments_to_be_evaluable(self) -> None:
        rows = self.rows(0.004, bases=("ADAUSD", "XRPUSD"))
        gates = runner.evaluate_gates(rows)
        self.assertFalse(gates["G2_segment_consistency"]["pass"])
        self.assertEqual(gates["G2_segment_consistency"]["segments_evaluated"], ["alt"])

    def test_g3_uses_exact_three_fifths_arithmetic(self) -> None:
        rows = self.rows(0.004, quarters=10)
        negatives = {"2020Q1", "2020Q2", "2020Q3", "2020Q4"}      # 6 of 10 positive -> 30 >= 30
        for row in rows:
            if row["quarter"] in negatives:
                row["net_return"] = -0.01
        gates = runner.evaluate_gates(rows)
        self.assertEqual(gates["G3_quarter_stability"]["positive_quarters"], 6)
        self.assertTrue(gates["G3_quarter_stability"]["pass"])
        for row in rows:
            if row["quarter"] == "2021Q1":
                row["net_return"] = -0.01                          # 5 of 10 -> 25 < 30
        self.assertFalse(runner.evaluate_gates(rows)["G3_quarter_stability"]["pass"])

    def test_g4_compares_against_the_cycle_seven_alpha(self) -> None:
        gates = runner.evaluate_gates(self.rows(0.004))
        self.assertEqual(gates["G4_dependence_robustness"]["alpha"], 0.05 / 56)
        self.assertLess(gates["G4_dependence_robustness"]["one_sided_block_bootstrap_p"], 0.05 / 56)


class ContractBindingTests(unittest.TestCase):
    def test_frozen_contract_binds_the_runner_parameters(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["frozen_rule"]["parameters"], runner.FROZEN_PARAMETERS)
        runner.check_contract_parameters(contract)

    def test_a_single_parameter_edit_desynchronises_and_aborts(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        contract["frozen_rule"]["parameters"]["TAKER_FEE"] = 0.0001
        with self.assertRaisesRegex(ValueError, "parameter mismatch"):
            runner.check_contract_parameters(contract)

    def test_static_input_hash_drift_aborts_before_any_archive_check(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        contract["assets_and_data"]["input_inventory_sha256"] = "0" * 64
        with mock.patch.object(runner, "verify_archive") as verify:
            with self.assertRaisesRegex(ValueError, "input_inventory_sha256 mismatch"):
                runner.verify_static_inputs(contract)
        verify.assert_not_called()

    def test_inventory_mutations_are_each_refused(self) -> None:
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        good = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        directory = Path(tempfile.mkdtemp(prefix=".tmp-dts-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, directory, ignore_errors=True)
        fake = directory / "inventory.json"
        pins = {runner.INVENTORY_PATH: contract["assets_and_data"]["input_inventory_sha256"],
                runner.AUDIT_PATH: contract["novelty_and_consumption"]["prior_window_audit_sha256"],
                runner.POLICY_PATH: contract["loop_multiplicity"]["policy_sha256"]}
        inventory_pin = pins[runner.INVENTORY_PATH]
        for mutate, message, fix_count in (
                (lambda inv: inv.update(candidate_outcomes_accessed=True), "pre-outcome", True),
                (lambda inv: inv["files"].pop(), "archive count mismatch", False),
                (lambda inv: inv["files"].append(dict(inv["files"][0])), "duplicate paths", True),
                (lambda inv: inv["pair_universe"].pop(), "pinned pair universe", True),
                (lambda inv: inv["symbols"].remove("cm:BTCUSD_210625"), "pinned pair universe", True)):
            inventory = json.loads(json.dumps(good))
            mutate(inventory)
            if fix_count:
                inventory["archive_count"] = len(inventory["files"])
            fake.write_text(json.dumps(inventory), encoding="utf-8")
            with mock.patch.object(runner, "INVENTORY_PATH", fake), \
                 mock.patch.object(runner, "sha256_file",
                                   side_effect=lambda p: pins.get(p, inventory_pin)), \
                 mock.patch.object(runner, "verify_archive"):
                with self.assertRaisesRegex(ValueError, message):
                    runner.verify_static_inputs(contract)


class ParsingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-dts-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)

    def load(self, csv_rows, header=True, **overrides):
        entry = write_archive(self.directory, "AAAUSD_240628", "2024-04", csv_rows, header=header)
        entry.update(overrides)
        return runner.load_contract_days({"files": [entry]})

    def test_a_complete_day_yields_its_closing_bar(self) -> None:
        table = self.load(full_day("2024-04-01", close=123.5, trades=7))
        slot = table[("cm", "AAAUSD_240628")]["2024-04-01"]
        self.assertEqual((slot["close"], slot["trades"], len(slot["hours"])), (123.5, 7, 24))

    def test_both_header_forms_parse_identically(self) -> None:
        rows = full_day("2024-04-01", close=99.0)
        with_header = self.load(rows, header=True)[("cm", "AAAUSD_240628")]["2024-04-01"]
        headerless = self.load(rows, header=False)[("cm", "AAAUSD_240628")]["2024-04-01"]
        self.assertEqual(with_header["close"], headerless["close"])
        self.assertEqual(len(with_header["hours"]), len(headerless["hours"]))

    def test_header_form_drift_aborts(self) -> None:
        with self.assertRaisesRegex(ValueError, "header form drift"):
            self.load(full_day("2024-04-01"), header=True, header_present=False)

    def test_duplicate_timestamp_aborts(self) -> None:
        rows = full_day("2024-04-01")
        with self.assertRaisesRegex(ValueError, "duplicate timestamp"):
            self.load(rows + [list(rows[0])])

    def test_off_hour_timestamp_aborts(self) -> None:
        rows = full_day("2024-04-01")
        rows.append([str(day_ms("2024-04-02") + 61_000)] + rows[0][1:])
        with self.assertRaisesRegex(ValueError, "off-hour timestamp"):
            self.load(rows)

    def test_short_row_aborts(self) -> None:
        rows = full_day("2024-04-01")
        rows.append(rows[0][:-1])
        with self.assertRaisesRegex(ValueError, "schema column count"):
            self.load(rows)

    def test_row_count_drift_aborts(self) -> None:
        with self.assertRaisesRegex(ValueError, "row count drift"):
            self.load(full_day("2024-04-01"), rows=999)

    def test_a_row_is_assigned_to_its_own_utc_day_not_the_archive_month(self) -> None:
        """Off-month rows are defined, not silently dropped; the census found none."""
        rows = full_day("2024-04-01") + [[str(day_ms("2024-05-01") + 23 * HOUR_MS)] + ["1"] * 3
                                         + ["55.0", "1", "1", "1", "9", "1", "1", "0"]]
        table = self.load(rows)[("cm", "AAAUSD_240628")]
        self.assertEqual(table["2024-05-01"]["close"], 55.0)
        self.assertEqual(len(table["2024-05-01"]["hours"]), 1)

    def test_a_partial_day_is_parsed_but_never_tradable(self) -> None:
        table = self.load(full_day("2024-04-01", hours=20))[("cm", "AAAUSD_240628")]
        self.assertEqual(len(table["2024-04-01"]["hours"]), 20)
        self.assertEqual(runner.tradable_days(table), set())


class RealArchiveTests(unittest.TestCase):
    """One test bound to the real frozen corpus, so a silent parser change cannot pass."""

    def test_a_real_archive_parses_to_its_inventoried_row_count_and_complete_days(self) -> None:
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        entry = next(f for f in inventory["files"]
                     if f["symbol"] == "BTCUSD_210326" and f["month"] == "2021-01")
        table = runner.load_contract_days({"files": [entry]})[("cm", "BTCUSD_210326")]
        self.assertEqual(sum(len(v["hours"]) for v in table.values()), entry["rows"])
        self.assertEqual(entry["rows"], 744)                       # 31 days x 24 hours, January 2021
        self.assertEqual(len([d for d in table if len(table[d]["hours"]) == 24]), 31)
        self.assertEqual(len(runner.tradable_days(table)), 31)
        self.assertTrue(all(v["close"] > 0 for v in table.values()))

    def test_the_real_frozen_archives_verify_against_their_upstream_sidecars(self) -> None:
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        for entry in inventory["files"][:12]:
            self.assertEqual(runner.verify_archive(entry), entry["sha256"])

    def test_a_flipped_inventory_hash_is_caught(self) -> None:
        inventory = json.loads(runner.INVENTORY_PATH.read_text(encoding="utf-8"))
        entry = dict(inventory["files"][0])
        entry["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "archive hash drift"):
            runner.verify_archive(entry)


class WindowTests(unittest.TestCase):
    """run_window is glue; these pin window selection, rejection counting and the sample floor.

    The builder derives days FROM the candidates, so both legs are always jointly tradable on the
    hand-over day and on the near contract's last day. An earlier version keyed days off each
    contract's own life and produced zero positions, which made every assertion here vacuous.
    """

    def synthetic(self, quarters: int, bases: tuple[str, ...], convergence: float = 1.0):
        start = date(2024, 3, 29)
        universe: list[dict] = []
        for index in range(1, quarters + 1):
            deliv = start + timedelta(days=91 * index)
            for base in bases:
                universe.append({"margin": "cm", "base": base,
                                 "near": f"{base}_{deliv:%y%m%d}",
                                 "far": f"{base}_{deliv + timedelta(days=91):%y%m%d}",
                                 "prev_delivery": (deliv - timedelta(days=91)).isoformat(),
                                 "near_delivery": deliv.isoformat(),
                                 "quarter": runner.quarter_label(deliv)})
        days: dict = {}

        def put(symbol: str, day: str, close: float) -> None:
            table = days.setdefault(("cm", symbol), {})
            table[day] = {"hours": set(range(24)), "trades": 5, "close": close}

        for candidate in universe:
            entry = (date.fromisoformat(candidate["prev_delivery"]) + timedelta(days=1)).isoformat()
            exit_ = (date.fromisoformat(candidate["near_delivery"]) - timedelta(days=1)).isoformat()
            put(candidate["near"], entry, 100.0)
            put(candidate["near"], exit_, 100.0)
            put(candidate["far"], entry, 102.0)                     # contango at the hand-over
            put(candidate["far"], exit_, 102.0 - 2.0 * convergence)  # 1.0 -> full convergence
        return universe, days

    WINDOW = {"windows_by_near_delivery": {"primary": {"start_inclusive": "2000-01-01",
                                                       "end_exclusive": "2100-01-01"}}}
    # 13 quarters x 6 bases = 78 positions, clearing MIN_POSITIONS=70 with 26 majors / 52 alts.
    SIX = ("BTCUSD", "ETHUSD", "ADAUSD", "XRPUSD", "LTCUSD", "BNBUSD")

    def test_the_builder_actually_produces_positions(self) -> None:
        universe, days = self.synthetic(13, ("BTCUSD", "ETHUSD", "ADAUSD", "XRPUSD"))
        result = runner.run_window(self.WINDOW, "primary", {"pair_universe": universe}, days)
        self.assertEqual(result["sample"]["positions"], len(universe))
        self.assertEqual(sum(result["rejected"].values()), 0)
        self.assertEqual(result["sample"]["sides"], {"long-near": len(universe), "long-far": 0})

    def test_a_converging_synthetic_panel_passes_and_a_static_one_is_falsified(self) -> None:
        universe, days = self.synthetic(13, self.SIX)
        good = runner.run_window(self.WINDOW, "primary", {"pair_universe": universe}, days)
        self.assertEqual(good["status"], "pass-on-this-sample")
        universe, days = self.synthetic(13, self.SIX, convergence=0.0)
        flat = runner.run_window(self.WINDOW, "primary", {"pair_universe": universe}, days)
        self.assertEqual(flat["status"], "falsified")
        self.assertLess(flat["diagnostics"]["pooled_bps"]["mean_net_return"], 0.0)

    def test_an_insufficient_sample_is_inconclusive_not_falsified(self) -> None:
        universe, days = self.synthetic(3, ("BTCUSD", "ETHUSD"))
        result = runner.run_window(self.WINDOW, "primary", {"pair_universe": universe}, days)
        self.assertEqual(result["status"], "inconclusive-insufficient-sample")
        self.assertFalse(result["sample_minimum_met"])
        self.assertNotIn("gate_results", result)

    def test_window_bounds_select_by_near_delivery_date(self) -> None:
        universe, days = self.synthetic(13, ("BTCUSD", "ETHUSD", "ADAUSD", "XRPUSD"))
        contract = {"windows_by_near_delivery": {"primary": {"start_inclusive": "2000-01-01",
                                                             "end_exclusive": "2025-01-01"}}}
        result = runner.run_window(contract, "primary", {"pair_universe": universe}, days)
        self.assertTrue(result["sample"]["positions"] > 0)
        self.assertTrue(all(c["near_delivery"] < "2025-01-01"
                            for c in universe if c["quarter"] in result["sample"]["quarter_list"]))
        self.assertLess(result["sample"]["positions"], len(universe))

    def test_rejections_are_counted_never_silently_dropped(self) -> None:
        """The dropped contract is the far leg of one candidate and the near leg of the next,
        so exactly two candidates lose a tradable day and both are counted."""
        universe, days = self.synthetic(13, self.SIX)
        days.pop(("cm", universe[0]["far"]))
        result = runner.run_window(self.WINDOW, "primary", {"pair_universe": universe}, days)
        self.assertEqual(result["rejected"]["no-complete-entry-day"], 2)
        self.assertEqual(sum(result["rejected"].values()), 2)
        self.assertEqual(result["sample"]["positions"], len(universe) - 2)


class LifecycleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.directory = Path(tempfile.mkdtemp(prefix=".tmp-dts-", dir=runner.EXPLORATIONS))
        self.addCleanup(shutil.rmtree, self.directory, ignore_errors=True)
        self.contract_sha = runner.sha256_file(runner.CONTRACT_PATH)
        contract = json.loads(runner.CONTRACT_PATH.read_text(encoding="utf-8"))
        # main() builds pins from the contract itself, so a real artifact carries these values.
        self.pins = {"input_inventory_sha256": contract["assets_and_data"]["input_inventory_sha256"],
                     "prior_window_audit_sha256":
                         contract["novelty_and_consumption"]["prior_window_audit_sha256"],
                     "policy_sha256": contract["loop_multiplicity"]["policy_sha256"],
                     "runner_sha256": runner.sha256_file(Path(runner.__file__).resolve())}

    def fake_result(self, status: str) -> dict:
        return {"window": "primary", "status": status, "sample": {"positions": 112,
                "included_quarters": 13}}

    def run_main(self, args, status="pass-on-this-sample"):
        with mock.patch.object(runner, "verify_static_inputs", return_value=({}, 1662)), \
             mock.patch.object(runner, "sha256_file", side_effect=self.hashes), \
             mock.patch.object(runner, "load_contract_days", return_value={}), \
             mock.patch.object(runner, "run_window", return_value=self.fake_result(status)):
            return runner.main(args)

    def hashes(self, path):
        path = Path(path)
        if path == runner.CONTRACT_PATH:
            return self.contract_sha
        if path.name.startswith("run_delivery_term_structure"):
            return self.pins["runner_sha256"]
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_verify_inputs_only_writes_nothing_and_prints_every_pin(self) -> None:
        out = io.StringIO()
        never = self.directory / "never.json"
        with contextlib.redirect_stdout(out), \
             mock.patch.object(runner, "verify_static_inputs", return_value=({}, 1662)), \
             mock.patch.object(runner, "sha256_file", side_effect=self.hashes), \
             mock.patch.object(runner, "load_contract_days") as loader:
            runner.main(["--contract-sha256", self.contract_sha, "--window", "primary",
                         "--output", str(never), "--verify-inputs-only"])
        loader.assert_not_called()
        self.assertFalse(never.exists())
        printed = json.loads(out.getvalue())
        for field in runner.PIN_FIELDS:
            self.assertIn(field, printed)

    def test_contract_hash_mismatch_aborts(self) -> None:
        with self.assertRaisesRegex(SystemExit, "contract hash mismatch"):
            self.run_main(["--contract-sha256", "0" * 64, "--window", "primary",
                           "--output", str(self.directory / "x.json")])

    def test_existing_output_is_never_overwritten(self) -> None:
        target = self.directory / "taken.json"
        target.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                           "--output", str(target)])

    def test_output_outside_the_workspace_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "escapes workspace root"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                           "--output", "/tmp/dts-never-written.json"])

    def test_primary_writes_a_pinned_artifact(self) -> None:
        target = self.directory / "primary.json"
        self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                       "--output", str(target)])
        payload = json.loads(target.read_text(encoding="utf-8"))
        self.assertEqual(payload["role"], "primary")
        self.assertEqual(payload["contract_sha256"], self.contract_sha)
        self.assertEqual(payload["static_inputs"]["runner_sha256"], self.pins["runner_sha256"])

    def write_real_primary(self, name: str = "primary.json", status: str = "pass-on-this-sample") -> Path:
        """Produce a primary artifact through main(), so downstream tests see the real schema."""
        target = self.directory / name
        self.run_main(["--contract-sha256", self.contract_sha, "--window", "primary",
                       "--output", str(target)], status=status)
        return target

    def test_a_real_primary_artifact_is_accepted_by_the_replication_validator(self) -> None:
        """Round trip: the writer's schema must satisfy the validator that reads it.

        Regression for a pins-nesting bug: main() nests pins under static_inputs while the
        validator originally read them at the top level, so every genuine hand-off aborted
        while both suites passed against hand-built fixtures.
        """
        primary = json.loads(self.write_real_primary().read_text(encoding="utf-8"))
        runner.validate_primary_artifact(primary, self.contract_sha, self.pins)

    def test_replication_runs_end_to_end_against_its_own_primary(self) -> None:
        primary_path = self.write_real_primary()
        primary_sha = hashlib.sha256(primary_path.read_bytes()).hexdigest()
        target = self.directory / "replication.json"
        with mock.patch.object(runner, "verify_static_inputs", return_value=({}, 1662)), \
             mock.patch.object(runner, "sha256_file", side_effect=self.hashes), \
             mock.patch.object(runner, "load_contract_days", return_value={}), \
             mock.patch.object(runner, "run_window",
                               return_value={"window": "replication", "status": "pass-on-this-sample",
                                             "sample": {"positions": 77, "included_quarters": 10}}):
            runner.main(["--contract-sha256", self.contract_sha, "--window", "replication",
                         "--output", str(target), "--primary-artifact", str(primary_path),
                         "--primary-artifact-sha256", primary_sha])
        payload = json.loads(target.read_text(encoding="utf-8"))
        self.assertEqual(payload["role"], "replication")
        self.assertEqual(payload["primary_artifact_sha256"], primary_sha)

    def test_replication_requires_the_primary_artifact_and_its_hash(self) -> None:
        with self.assertRaisesRegex(SystemExit, "requires --primary-artifact"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                           "--output", str(self.directory / "r.json")])

    def test_replication_refuses_a_drifted_primary_file(self) -> None:
        primary_path = self.write_real_primary()
        with self.assertRaisesRegex(SystemExit, "primary artifact hash mismatch"):
            self.run_main(["--contract-sha256", self.contract_sha, "--window", "replication",
                           "--output", str(self.directory / "r.json"),
                           "--primary-artifact", str(primary_path),
                           "--primary-artifact-sha256", "0" * 64])

    def test_replication_refuses_a_falsified_primary(self) -> None:
        primary = json.loads(self.write_real_primary(status="falsified").read_text(encoding="utf-8"))
        with self.assertRaisesRegex(ValueError, "requires a committed primary pass"):
            runner.validate_primary_artifact(primary, self.contract_sha, self.pins)

    def test_replication_refuses_a_primary_from_a_different_runner(self) -> None:
        primary = json.loads(self.write_real_primary().read_text(encoding="utf-8"))
        primary["static_inputs"]["runner_sha256"] = "9" * 64
        with self.assertRaisesRegex(ValueError, "runner_sha256 mismatch"):
            runner.validate_primary_artifact(primary, self.contract_sha, self.pins)

    def test_replication_refuses_a_primary_under_a_different_contract(self) -> None:
        primary = json.loads(self.write_real_primary().read_text(encoding="utf-8"))
        with self.assertRaisesRegex(ValueError, "contract_sha256 mismatch"):
            runner.validate_primary_artifact(primary, "0" * 64, self.pins)

    def test_replication_refuses_an_artifact_with_no_static_inputs(self) -> None:
        primary = json.loads(self.write_real_primary().read_text(encoding="utf-8"))
        primary.pop("static_inputs")
        with self.assertRaisesRegex(ValueError, "no static_inputs block"):
            runner.validate_primary_artifact(primary, self.contract_sha, self.pins)

    def test_replication_refuses_a_non_primary_role(self) -> None:
        with self.assertRaisesRegex(ValueError, "not role primary"):
            runner.validate_primary_artifact({"role": "replication"}, self.contract_sha, self.pins)


if __name__ == "__main__":
    unittest.main()
