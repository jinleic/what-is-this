"""Contract tests for the memory-safe fixed-five certificate scheduler."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import campaign_runtime  # noqa: E402
import certify_involution_f5_hierarchical_batch as h_batch  # noqa: E402
import certify_involution_f5_w_dfs_batch as w_batch  # noqa: E402
import run_involution_f5_certificate_campaign as scheduler  # noqa: E402

_STATE = {
    "w_done": 24,
    "w_remaining": 526,
    "h_done": 1,
    "h_remaining": 154,
}


class TestCoverageMetadata(unittest.TestCase):
    def test_valid_coverage_is_read(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.json"
            path.write_text(json.dumps({
                "coverage": {
                    "certified_support_orbits": 7,
                    "target_w_negative_support_orbits": 10,
                    "remaining_support_orbits": 3,
                    "campaign_complete": False,
                },
                "records": [
                    {"source_index": index} for index in range(7)],
            }))
            self.assertEqual(
                scheduler._coverage(
                    path, "target_w_negative_support_orbits", 10),
                (7, 10))

    def test_wrong_target_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.json"
            path.write_text(json.dumps({
                "coverage": {
                    "certified_support_orbits": 7,
                    "target_w_negative_support_orbits": 9,
                    "remaining_support_orbits": 3,
                    "campaign_complete": False,
                },
                "records": [
                    {"source_index": index} for index in range(7)],
            }))
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "invalid coverage metadata"):
                scheduler._coverage(
                    path, "target_w_negative_support_orbits", 10)

    def test_malformed_json_is_normalized_to_scheduler_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.json"
            path.write_text("{not json")
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "cannot load"):
                scheduler._coverage(
                    path, "target_w_negative_support_orbits", 10)


class TestFairScheduling(unittest.TestCase):
    def test_one_cycle_runs_w_then_h_with_bounded_counts(self):
        states = iter([
            {"w_done": 24, "w_remaining": 526,
             "h_done": 1, "h_remaining": 154},
            {"w_done": 27, "w_remaining": 523,
             "h_done": 1, "h_remaining": 154},
            {"w_done": 27, "w_remaining": 523,
             "h_done": 2, "h_remaining": 153},
        ])
        calls = []

        def run_batch(kind, limit, timeout, max_timeout):
            calls.append((kind, limit, timeout, max_timeout))

        with (mock.patch.object(
                scheduler, "_snapshot", side_effect=lambda: next(states)),
              mock.patch.object(
                  scheduler, "_run_batch", side_effect=run_batch),
              mock.patch.object(scheduler, "_emit")):
            result = scheduler.run(3, 1, 21600.0, 86400.0, True)

        self.assertEqual(calls, [
            ("w", 3, 21600.0, 86400.0),
            ("h", 1, 21600.0, 86400.0),
        ])
        self.assertEqual(result["w_done"], 27)
        self.assertEqual(result["h_done"], 2)


class TestCompletionValidation(unittest.TestCase):
    def test_full_validation_precedes_campaign_complete_event(self):
        state = {
            "w_done": 550,
            "w_remaining": 0,
            "h_done": 155,
            "h_remaining": 0,
        }
        order = []

        def validate():
            order.append("validate")

        def emit(event, **_fields):
            order.append(event)

        with (mock.patch.object(scheduler, "_snapshot", return_value=state),
              mock.patch.object(
                  scheduler, "_validate_complete_ledgers",
                  side_effect=validate),
              mock.patch.object(scheduler, "_emit", side_effect=emit)):
            result = scheduler.run(3, 1, 100.0, 400.0, False)

        self.assertEqual(result, state)
        self.assertEqual(order, ["validate", "campaign_complete"])

    def test_stale_complete_ledger_is_journaled_as_terminal_failure(self):
        state = {
            "w_done": 550,
            "w_remaining": 0,
            "h_done": 155,
            "h_remaining": 0,
        }
        with (mock.patch.object(scheduler, "_snapshot", return_value=state),
              mock.patch.object(
                  scheduler, "_validate_complete_ledgers",
                  side_effect=scheduler.SchedulerFailure(
                      "stale complete ledger")),
              mock.patch.object(scheduler, "_emit") as emit):
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "stale complete ledger"):
                scheduler.run(3, 1, 100.0, 400.0, False)

        events = [call.args[0] for call in emit.call_args_list]
        self.assertNotIn("campaign_complete", events)
        self.assertEqual(events[-1], "terminal_failure")


    def test_complete_gate_calls_both_strict_driver_validators(self):
        w_document = {
            "records": [
                {"source_index": source} for source in range(550)]}
        h_document = {
            "records": [
                {"source_index": source} for source in range(550, 705)]}
        with tempfile.TemporaryDirectory() as directory:
            with (mock.patch.object(
                    campaign_runtime, "DRIVER_LOCK",
                    Path(directory) / "driver.lock"),
                  mock.patch.object(
                      w_batch, "_validate_output_unlocked",
                      return_value=w_document) as validate_w,
                  mock.patch.object(
                      h_batch, "_validate_output_unlocked",
                      return_value=h_document) as validate_h):
                scheduler._validate_complete_ledgers()
        validate_w.assert_called_once_with(
            scheduler.W_OUTPUT, require_complete=True)
        validate_h.assert_called_once_with(
            scheduler.H_OUTPUT, require_complete=True)

    def test_complete_gate_rejects_cross_ledger_overlap(self):
        w_document = {
            "records": [
                {"source_index": source} for source in range(550)]}
        h_document = {
            "records": [
                {"source_index": source} for source in range(549, 704)]}
        with tempfile.TemporaryDirectory() as directory:
            with (mock.patch.object(
                    campaign_runtime, "DRIVER_LOCK",
                    Path(directory) / "driver.lock"),
                  mock.patch.object(
                      w_batch, "_validate_output_unlocked",
                      return_value=w_document),
                  mock.patch.object(
                      h_batch, "_validate_output_unlocked",
                      return_value=h_document)):
                with self.assertRaisesRegex(
                        scheduler.SchedulerFailure, "partition 705"):
                    scheduler._validate_complete_ledgers()


class TestRetryAllowance(unittest.TestCase):
    def test_timeout_retry_deducts_persisted_records(self):
        snapshots = iter([
            _STATE,
            {**_STATE, "w_done": 25, "w_remaining": 525},
            {**_STATE, "w_done": 27, "w_remaining": 523},
        ])
        commands = []

        def supervised(command, _environment):
            commands.append(command)
            if len(commands) == 1:
                return (subprocess.CompletedProcess(
                            command, 1, "Traceback: TimeoutExpired\n", ""),
                        None)
            return subprocess.CompletedProcess(command, 0, "ok\n", ""), None

        with (mock.patch.object(
                scheduler, "_snapshot", side_effect=lambda: next(snapshots)),
              mock.patch.object(scheduler, "_guard_host"),
              mock.patch.object(
                  scheduler, "_supervised_run", side_effect=supervised),
              mock.patch.object(scheduler, "_emit") as emit,
              mock.patch("builtins.print")):
            scheduler._run_batch("w", 3, 100.0, 400.0)

        self.assertEqual(len(commands), 2)
        self.assertEqual(
            commands[0][commands[0].index("--max-new-records") + 1], "3")
        self.assertEqual(
            commands[1][commands[1].index("--max-new-records") + 1], "2")
        self.assertEqual(
            commands[1][commands[1].index("--timeout-per-step") + 1],
            "200.0")
        retry = next(
            call for call in emit.call_args_list
            if call.args[0] == "timeout_retry")
        self.assertEqual(retry.kwargs["remaining_allowance"], 2)

    def test_spent_allowance_ends_batch_without_retry(self):
        snapshots = iter([
            _STATE,
            {**_STATE, "w_done": 27, "w_remaining": 523},
        ])
        completed = subprocess.CompletedProcess(
            ["driver"], 1, "TimeoutExpired", "")
        supervised = mock.Mock(return_value=(completed, None))
        with (mock.patch.object(
                scheduler, "_snapshot", side_effect=lambda: next(snapshots)),
              mock.patch.object(scheduler, "_guard_host"),
              mock.patch.object(
                  scheduler, "_supervised_run", supervised),
              mock.patch.object(scheduler, "_emit"),
              mock.patch("builtins.print")):
            scheduler._run_batch("w", 3, 100.0, 400.0)
        supervised.assert_called_once()

    def test_non_timeout_failure_is_not_retried(self):
        completed = subprocess.CompletedProcess(
            ["driver"], 1, "VeriPB rejected proof", "")
        supervised = mock.Mock(return_value=(completed, None))
        with (mock.patch.object(scheduler, "_snapshot", return_value=_STATE),
              mock.patch.object(scheduler, "_guard_host"),
              mock.patch.object(
                  scheduler, "_supervised_run", supervised),
              mock.patch.object(scheduler, "_emit"),
              mock.patch("builtins.print")):
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "not retrying"):
                scheduler._run_batch("w", 3, 100.0, 400.0)
        supervised.assert_called_once()

    def test_maximum_timeout_failure_is_not_retried(self):
        completed = subprocess.CompletedProcess(
            ["driver"], 1, "TimeoutExpired", "")
        supervised = mock.Mock(return_value=(completed, None))
        with (mock.patch.object(scheduler, "_snapshot", return_value=_STATE),
              mock.patch.object(scheduler, "_guard_host"),
              mock.patch.object(
                  scheduler, "_supervised_run", supervised),
              mock.patch.object(scheduler, "_emit"),
              mock.patch("builtins.print")):
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "maximum timeout"):
                scheduler._run_batch("w", 3, 400.0, 400.0)
        supervised.assert_called_once()


class TestDiskFloorSupervision(unittest.TestCase):
    def test_supervised_run_terminates_on_polled_floor_and_keeps_output(self):
        child = mock.Mock(pid=4321, returncode=-15)
        child.communicate.side_effect = [
            subprocess.TimeoutExpired(
                ["driver"], scheduler.POLL_SECONDS),
            ("retained stdout", "retained stderr"),
        ]
        disk = mock.Mock(free=1 << 30)
        with (mock.patch.object(
                scheduler.subprocess, "Popen", return_value=child) as popen,
              mock.patch.object(
                  scheduler.shutil, "disk_usage", return_value=disk),
              mock.patch.object(
                  scheduler, "_terminate_group") as terminate):
            completed, breach = scheduler._supervised_run(
                ["driver"], {"ONLY": "TEST"})

        terminate.assert_called_once_with(child)
        self.assertEqual(completed.returncode, -15)
        self.assertEqual(completed.stdout, "retained stdout")
        self.assertEqual(completed.stderr, "retained stderr")
        self.assertIn("batch process group terminated", breach)
        self.assertTrue(popen.call_args.kwargs["start_new_session"])

    def test_supervised_run_reaps_group_after_normal_leader_exit(self):
        child = mock.Mock(pid=4321, returncode=0)
        child.communicate.return_value = ("stdout", "stderr")
        with (
            mock.patch.object(
                scheduler.subprocess, "Popen", return_value=child),
            mock.patch.object(
                scheduler, "_terminate_group") as terminate,
        ):
            completed, breach = scheduler._supervised_run(
                ["driver"], {"ONLY": "TEST"})
        terminate.assert_called_once_with(child)
        self.assertEqual(completed.returncode, 0)
        self.assertIsNone(breach)

    def test_keyboard_interrupt_cleans_child_spawned_with_signals_unblocked(self):
        child = mock.Mock(pid=4321)
        child.communicate.side_effect = [
            KeyboardInterrupt(),
            ("drained stdout", "drained stderr"),
        ]
        initial_mask = {scheduler.signal.SIGTERM}
        current_mask = set(initial_mask)
        spawn_masks = []

        def update_mask(operation, signals):
            previous = set(current_mask)
            if operation == scheduler.signal.SIG_BLOCK:
                current_mask.update(signals)
            elif operation == scheduler.signal.SIG_UNBLOCK:
                current_mask.difference_update(signals)
            elif operation == scheduler.signal.SIG_SETMASK:
                current_mask.clear()
                current_mask.update(signals)
            else:
                self.fail(f"unexpected signal-mask operation {operation}")
            return previous

        def spawn(*_args, **_kwargs):
            spawn_masks.append(set(current_mask))
            return child

        with (mock.patch.object(
                scheduler.subprocess, "Popen", side_effect=spawn),
              mock.patch.object(
                  scheduler, "_terminate_group") as terminate,
              mock.patch.object(
                  scheduler.signal, "getsignal",
                  return_value=scheduler.signal.SIG_DFL),
              mock.patch.object(scheduler.signal, "signal") as install,
              mock.patch.object(
                  scheduler.signal, "pthread_sigmask",
                  side_effect=update_mask),
              mock.patch("builtins.print") as printed):
            with self.assertRaises(KeyboardInterrupt):
                scheduler._supervised_run(["driver"], {"ONLY": "TEST"})

        terminate.assert_called_once_with(child)
        self.assertEqual(child.communicate.call_count, 2)
        self.assertEqual(spawn_masks, [set()])
        self.assertEqual(current_mask, initial_mask)
        install.assert_any_call(
            scheduler.signal.SIGINT, scheduler._raise_interrupt)
        install.assert_any_call(
            scheduler.signal.SIGINT, scheduler.signal.SIG_DFL)
        install.assert_any_call(
            scheduler.signal.SIGTERM, scheduler._raise_interrupt)
        install.assert_any_call(
            scheduler.signal.SIGTERM, scheduler.signal.SIG_DFL)
        printed.assert_any_call("drained stdout", end="", flush=True)
        printed.assert_any_call(
            "drained stderr", end="", file=scheduler.sys.stderr, flush=True)

    def test_signal_during_spawn_is_deferred_until_child_handle_exists(self):
        child = mock.Mock(pid=4321)
        child.communicate.return_value = (
            "spawn-interrupt stdout", "spawn-interrupt stderr")
        installed_handlers = {}

        def install(signal_number, handler):
            installed_handlers[signal_number] = handler

        def spawn(*_args, **_kwargs):
            installed_handlers[scheduler.signal.SIGTERM](
                scheduler.signal.SIGTERM, None)
            return child

        with (mock.patch.object(
                scheduler.subprocess, "Popen", side_effect=spawn),
              mock.patch.object(
                  scheduler, "_terminate_group") as terminate,
              mock.patch.object(
                  scheduler.signal, "getsignal",
                  return_value=scheduler.signal.SIG_DFL),
              mock.patch.object(
                  scheduler.signal, "signal", side_effect=install),
              mock.patch.object(
                  scheduler.signal, "pthread_sigmask",
                  return_value=set()),
              mock.patch("builtins.print") as printed):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "received SIGTERM"):
                scheduler._supervised_run(["driver"], {"ONLY": "TEST"})

        terminate.assert_called_once_with(child)
        child.communicate.assert_called_once_with()
        printed.assert_any_call(
            "spawn-interrupt stdout", end="", flush=True)
        printed.assert_any_call(
            "spawn-interrupt stderr", end="", file=scheduler.sys.stderr,
            flush=True)


    def test_group_kill_reaches_descendants_after_full_term_grace(self):
        child = mock.Mock(pid=4321)
        child.poll.return_value = 0
        with (mock.patch.object(scheduler.os, "killpg") as killpg,
              mock.patch.object(
                  scheduler.time, "monotonic",
                  side_effect=[100.0, 111.0]),
              mock.patch.object(scheduler.time, "sleep") as asleep):
            scheduler._terminate_group(child)
        self.assertEqual(killpg.call_args_list, [
            mock.call(4321, scheduler.signal.SIGTERM),
            mock.call(4321, 0),
            mock.call(4321, scheduler.signal.SIGKILL),
        ])
        asleep.assert_not_called()
        child.wait.assert_not_called()

    def test_descendants_receive_term_grace_after_leader_exits(self):
        child = mock.Mock(pid=4321)
        child.poll.return_value = 0
        killpg = mock.Mock(side_effect=[
            None,
            None,
            ProcessLookupError(),
        ])
        with (mock.patch.object(scheduler.os, "killpg", killpg),
              mock.patch.object(
                  scheduler.time, "monotonic",
                  side_effect=[100.0, 100.0]),
              mock.patch.object(scheduler.time, "sleep") as asleep):
            scheduler._terminate_group(child)
        self.assertEqual(killpg.call_args_list, [
            mock.call(4321, scheduler.signal.SIGTERM),
            mock.call(4321, 0),
            mock.call(4321, 0),
        ])
        asleep.assert_called_once_with(scheduler.GROUP_POLL_SECONDS)
        child.wait.assert_not_called()

    def test_floor_breach_is_retained_and_raised_without_retry(self):
        completed = subprocess.CompletedProcess(
            ["driver"], -15, "partial child output\n", "partial child error\n")
        breach = (
            "only 3.0 GiB free; require 20 GiB; "
            "batch process group terminated")
        supervised = mock.Mock(return_value=(completed, breach))
        with (mock.patch.object(scheduler, "_snapshot", return_value=_STATE),
              mock.patch.object(scheduler, "_guard_host"),
              mock.patch.object(
                  scheduler, "_supervised_run", supervised),
              mock.patch.object(scheduler, "_emit") as emit,
              mock.patch("builtins.print") as printed):
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "process group terminated"):
                scheduler._run_batch("w", 3, 100.0, 400.0)

        supervised.assert_called_once()
        stdout_calls = [
            call for call in printed.call_args_list
            if call.args and call.args[0] == "partial child output\n"]
        stderr_calls = [
            call for call in printed.call_args_list
            if call.args and call.args[0] == "partial child error\n"]
        self.assertEqual(len(stdout_calls), 1)
        self.assertEqual(len(stderr_calls), 1)
        batch_exit = next(
            call for call in emit.call_args_list
            if call.args[0] == "batch_exit")
        self.assertTrue(batch_exit.kwargs["terminated"])


class TestTerminalFailureJournal(unittest.TestCase):
    def test_run_journals_terminal_failure_before_reraising(self):
        with (mock.patch.object(scheduler, "_snapshot", return_value=_STATE),
              mock.patch.object(
                  scheduler, "_run_batch",
                  side_effect=scheduler.SchedulerFailure(
                      "w batch rejected")),
              mock.patch.object(scheduler, "_emit") as emit):
            with self.assertRaisesRegex(
                    scheduler.SchedulerFailure, "w batch rejected"):
                scheduler.run(3, 1, 100.0, 400.0, True)

        failures = [
            call for call in emit.call_args_list
            if call.args[0] == "terminal_failure"]
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0].kwargs["detail"], "w batch rejected")


class TestLatestEvent(unittest.TestCase):
    def test_missing_event_log_has_no_latest_event(self):
        with tempfile.TemporaryDirectory() as directory:
            events = Path(directory) / "events.jsonl"
            with mock.patch.object(scheduler, "EVENTS", events):
                self.assertIsNone(scheduler._last_event())

    def test_latest_parseable_event_survives_a_torn_trailing_line(self):
        with tempfile.TemporaryDirectory() as directory:
            events = Path(directory) / "events.jsonl"
            events.write_text(
                '{"event": "batch_start"}\n'
                "not json\n"
                '{"event": "terminal_failure"}\n'
                '{"event": "batch_start"')
            with mock.patch.object(scheduler, "EVENTS", events):
                self.assertEqual(
                    scheduler._last_event(), "terminal_failure")


class TestTerminalFailureRestart(unittest.TestCase):
    def test_restart_after_terminal_failure_succeeds_without_rerunning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            events = root / "events.jsonl"
            events.write_text(
                json.dumps({"event": "batch_start"}) + "\n"
                + json.dumps({"event": "terminal_failure"}) + "\n")
            with (mock.patch.object(scheduler, "EVENTS", events),
                  mock.patch.object(
                      scheduler, "_snapshot",
                      side_effect=scheduler.SchedulerFailure(
                          "unreadable coverage")) as snapshot,
                  mock.patch.object(scheduler, "run") as run_campaign,
                  mock.patch.object(
                      campaign_runtime, "SCHEDULER_LOCK",
                      root / "scheduler.lock"),
                  mock.patch("builtins.print")):
                self.assertEqual(scheduler.main([]), 0)
            run_campaign.assert_not_called()
            snapshot.assert_not_called()

    def test_resume_after_failure_acknowledges_before_proceeding(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            events = root / "events.jsonl"
            events.write_text(
                json.dumps({"event": "terminal_failure"}) + "\n")
            order = []

            def emit(event, **_fields):
                order.append(event)

            def run_campaign(*_args):
                order.append("run")
                return dict(_STATE)

            with (mock.patch.object(scheduler, "EVENTS", events),
                  mock.patch.object(
                      scheduler, "_snapshot", return_value=_STATE),
                  mock.patch.object(
                      scheduler, "_emit", side_effect=emit),
                  mock.patch.object(
                      scheduler, "run", side_effect=run_campaign),
                  mock.patch.object(
                      campaign_runtime, "SCHEDULER_LOCK",
                      root / "scheduler.lock"),
                  mock.patch("builtins.print")):
                self.assertEqual(scheduler.main(
                    ["--one-cycle", "--resume-after-failure"]), 0)
            self.assertEqual(order, ["failure_acknowledged", "run"])

    def test_acknowledged_failure_does_not_block_later_restart(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            events = root / "events.jsonl"
            events.write_text(
                json.dumps({"event": "terminal_failure"}) + "\n"
                + json.dumps({"event": "failure_acknowledged"}) + "\n")
            with (mock.patch.object(scheduler, "EVENTS", events),
                  mock.patch.object(
                      scheduler, "_snapshot", return_value=_STATE),
                  mock.patch.object(
                      scheduler, "run", return_value=dict(_STATE))
                  as run_campaign,
                  mock.patch.object(
                      campaign_runtime, "SCHEDULER_LOCK",
                      root / "scheduler.lock"),
                  mock.patch("builtins.print")):
                self.assertEqual(scheduler.main([]), 0)
            run_campaign.assert_called_once()


class TestSchedulerInstanceLock(unittest.TestCase):
    def test_second_scheduler_instance_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / "scheduler.lock"
            with campaign_runtime.FileLock(lock):
                with (mock.patch.object(
                        campaign_runtime, "SCHEDULER_LOCK", lock),
                      mock.patch.object(
                          scheduler, "EVENTS", root / "events.jsonl"),
                      mock.patch.object(
                          scheduler, "_snapshot", return_value=_STATE),
                      mock.patch.object(scheduler, "run") as run_campaign):
                    with self.assertRaises(campaign_runtime.LockBusy):
                        scheduler.main([])
            run_campaign.assert_not_called()


if __name__ == "__main__":
    unittest.main()
