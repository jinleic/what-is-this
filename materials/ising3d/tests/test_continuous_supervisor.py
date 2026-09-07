#!/usr/bin/env python3
"""Isolated acceptance regressions; no model calls or numerical jobs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('supervisor_test_subject', ROOT / 'tools' / 'continuous_supervisor.py')
coordinator = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = coordinator
spec.loader.exec_module(coordinator)


class AcceptanceBoundaries(unittest.TestCase):
    def setUp(self):
        self.scratch = Path(tempfile.mkdtemp(prefix='ising3d-supervisor-regression-'))
        self.root = self.scratch / 'workspace'
        self.target = self.root / 'math' / 'ising3d'
        for subdir in ('experiments', 'tests', 'results'):
            (self.target / subdir).mkdir(parents=True)
        self.paths = patch.multiple(coordinator, ROOT=self.root, DOMAIN_DIR=self.root / 'math', TARGET_DIR=self.target)
        self.paths.start()
        self.addCleanup(self.paths.stop)
        self.sup = coordinator.Supervisor(coordinator.build_parser().parse_args(['check']))
        (self.target / 'experiments' / 'candidate.py').write_text('BOUND = 7\n')
        (self.target / 'tests' / 'verify.py').write_text('raise SystemExit(0)\n')
        self.plan = {
            'decision': 'run', 'front': 'norm', 'hypothesis_id': 'H681',
            'next_prompt': 'A bounded exact trial with explicitly pinned artifact data and no wider claim.',
            'artifact': 'results/canary.json',
            'producer': {'script': 'experiments/candidate.py', 'arguments': ['--output', 'results/canary.json']},
            'verifier': {'script': 'tests/verify.py', 'arguments': [], 'pass_marker': 'CANARY_VERIFIED PASS'},
            'verifier_requirements': 'Reconstruct the exact bounded rank from independent inputs.',
            'artifact_expect': {'data.rank': 9}, 'checkpoint_evidence': ['checkpoint preregistration'],
            'source_changes': [],
        }

    def artifact(self, *, rank=9, version=1):
        data = {'rank': rank}
        source = self.target / 'experiments' / 'candidate.py'
        doc = {'meta': {
            'schema_version': version,
            'source_sha256': {'math/ising3d/experiments/candidate.py': hashlib.sha256(source.read_bytes()).hexdigest()},
            'data_sha256': hashlib.sha256(coordinator.canonical(data).encode()).hexdigest(),
        }, 'data': data}
        (self.target / self.plan['artifact']).write_text(coordinator.canonical(doc) + '\n')

    def test_failed_trial_cannot_be_retried_by_renaming_outputs_or_adding_unrelated_code(self):
        fingerprint = self.sup.fingerprint(self.plan)
        self.sup.status['history'] = [{'fingerprint': fingerprint, 'outcome': 'PRODUCER_RESOURCE_BLOCKED', 'task_id': 'prior'}]
        (self.target / 'experiments' / 'renamed.py').write_bytes((self.target / 'experiments' / 'candidate.py').read_bytes())
        (self.target / 'tests' / 'unrelated.py').write_text('OTHER_BOUND = 31\n')
        self.plan['producer']['script'] = 'experiments/renamed.py'
        self.plan['artifact'] = 'results/renamed.json'
        self.plan['producer']['arguments'] = ['--output', self.plan['artifact']]
        self.assertIsNotNone(self.sup.history_conflict(self.sup.fingerprint(self.plan)))

    def test_a_changed_mathematical_input_is_a_distinct_trial(self):
        self.plan['producer']['arguments'] += ['--size', '2']
        fingerprint = self.sup.fingerprint(self.plan)
        self.sup.status['history'] = [{'fingerprint': fingerprint, 'outcome': 'PRODUCER_RESOURCE_BLOCKED', 'task_id': 'prior'}]
        self.plan['producer']['arguments'][-1] = '3'
        self.assertIsNone(self.sup.history_conflict(self.sup.fingerprint(self.plan)))

    def test_verifier_only_seed_rejects_resealed_wrong_pinned_data_before_launch(self):
        self.artifact(rank=7)
        self.plan['producer'] = None
        errors = self.sup.validate_plan(self.plan, seed=True, floor=0, ids={'H681'})
        self.assertTrue(errors, 'existing artifact bypassed its pinned scalar preflight')

    def test_boolean_schema_version_is_not_integer_version_one(self):
        self.artifact(version=True)
        report = self.sup.validate_artifact(self.plan['artifact'], 0, self.plan['artifact_expect'])
        self.assertFalse(report['ok'])

    def test_guard_accounting_failures_never_land_in_the_resource_blocked_bucket(self):
        # RESOURCE_BLOCKED writes the trial into history and retires it forever,
        # so an accounting or teardown anomaly must classify separately.
        unmeasurable = {'guard_exit_code': 1, 'child_exit_code': -15,
                        'reason': 'unmeasurable descendant outlived the 30 s accounting window'}
        uncleaned = {'guard_exit_code': 1, 'child_exit_code': 0,
                     'reason': 'owned group survived its bounded cleanup'}
        exhausted = {'guard_exit_code': 1, 'child_exit_code': -9,
                     'reason': '2 GiB memory ceiling reached'}
        self.assertEqual(coordinator.Supervisor.classify(unmeasurable, 1), 'UNACCOUNTABLE')
        self.assertEqual(coordinator.Supervisor.classify(uncleaned, 1), 'UNACCOUNTABLE')
        self.assertEqual(coordinator.Supervisor.classify(exhausted, 1), 'RESOURCE_BLOCKED')

    def test_rejected_quota_window_is_read_from_the_cli_stream(self):
        rejected = json.dumps([
            {'type': 'rate_limit_event',
             'rate_limit_info': {'status': 'rejected', 'resetsAt': 1788592800}},
            {'type': 'result', 'subtype': 'success'},
        ]).encode()
        allowed = json.dumps([
            {'type': 'rate_limit_event',
             'rate_limit_info': {'status': 'allowed', 'resetsAt': 1788592800}},
        ]).encode()
        self.assertEqual(coordinator.quota_reset_unix(rejected), 1788592800.0)
        self.assertIsNone(coordinator.quota_reset_unix(allowed))
        self.assertIsNone(coordinator.quota_reset_unix(b'not json at all'))

    def test_quota_wait_is_bounded_and_costs_no_failure_budget(self):
        self.sup.state_dir = self.scratch / 'state'
        with patch.object(coordinator, 'QUOTA_WAIT_MIN_SECONDS', 0.0):
            self.sup.wait_for_quota(time.time() - 3600)
        self.assertEqual(self.sup.status['phase'], 'WAITING_QUOTA')
        self.assertIsNone(self.sup.status['quota_resets_utc'])
        self.assertEqual(self.sup.status['agent']['consecutive_failures'], 0)
        events = (self.sup.state_dir / coordinator.TRANSITIONS_NAME).read_text()
        self.assertIn('AGENT_RATE_LIMITED', events)

    def test_quota_wait_never_polls_the_window_back_to_back(self):
        self.sup.state_dir = self.scratch / 'floor'
        with patch.object(coordinator.Supervisor, 'sleep', autospec=True) as slept:
            self.sup.wait_for_quota(time.time() - 86400)
        self.assertGreaterEqual(slept.call_args[0][1], coordinator.QUOTA_WAIT_MIN_SECONDS)

    def test_unreadable_recovery_snapshot_halts_instead_of_raising(self):
        self.sup.state_dir = self.scratch / 'recovery'
        turn_dir = self.sup.state_dir / 'turns' / '0001_plan'
        turn_dir.mkdir(parents=True)
        self.sup.status['active_agent'] = {
            'pid': 2 ** 22, 'turn_dir': str(turn_dir), 'kind': 'plan',
            'result_path': str(turn_dir / 'guard.json'), 'deadline_unix': time.time() - 1,
        }
        self.sup.recover_agent_turn()
        self.assertEqual(self.sup.status['phase'], 'HALTED')
        self.assertEqual(self.sup.status['halt_reason'], 'RECOVERY_SNAPSHOT_INVALID')

    def test_claim_identity_separates_a_real_supervisor_from_a_reused_pid(self):
        launched = ("nice -n 19 ../.venv/bin/python -I -B tools/continuous_supervisor.py run "
                    "--seed-plan tools/seed.json")
        absolute = "/m/.venv/bin/python -I -B /m/ising3d/tools/continuous_supervisor.py run"
        self.assertTrue(coordinator.command_runs_this_script(launched))
        self.assertTrue(coordinator.command_runs_this_script(absolute))
        self.assertFalse(coordinator.command_runs_this_script(
            "/m/.venv/bin/python -I -B tests/test_continuous_supervisor.py"))
        self.assertFalse(coordinator.command_runs_this_script("/usr/bin/vim notes/plan.md"))
        # This live test process must never be mistaken for the supervisor.
        self.assertTrue(coordinator.pid_alive(os.getpid()))
        self.assertFalse(coordinator.supervisor_alive(os.getpid()))

    def test_checkpoint_notes_land_where_the_planner_actually_reads(self):
        head = '# Current state\n\n**Snapshot:** old.\n'
        body = '\n'.join(f'- filler line {index}' for index in range(2000))
        text = f'{head}\n## Open questions\n\n{body}\n'
        note = '\n### Supervisor outcome 0007_H700 — 2026-09-05\n\n- CONFIRMED.\n'
        updated = coordinator.insert_note(text, note)
        self.assertIn(note.strip(), updated[:6000])
        self.assertTrue(updated.startswith(head))
        self.assertIn(body, updated)
        self.assertEqual(coordinator.insert_note(updated, note), updated)


if __name__ == '__main__':
    unittest.main()
