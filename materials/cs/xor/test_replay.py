"""A deadline must never leave a successful or absent terminal receipt."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('xor_replay', Path(__file__).with_name('replay.py'))
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class ReceiptTests(unittest.TestCase):
    def run_interrupted(self, code, wall_seconds):
        # Keep interruption evidence for diagnosis; tests delete no files.
        root = Path(tempfile.mkdtemp(prefix='xor-receipt-regression-'))
        receipt = root / 'result.json'
        stage = replay.supervise([sys.executable, '-I', '-B', '-c', code], receipt, wall_seconds)
        payload = json.loads(receipt.read_text())
        self.assertEqual(stage['guard_status'], 'INCONCLUSIVE')
        self.assertEqual(payload['status'], 'INCONCLUSIVE')
        self.assertFalse(payload['primary_completed'])
        return stage

    def test_sigalrm_leaves_inconclusive_receipt(self):
        stage = self.run_interrupted('import signal,time; signal.alarm(1); time.sleep(4)', 5)
        self.assertEqual(stage['reason'], 'SIGNAL: SIGALRM')

    def test_cpu_limit_leaves_inconclusive_receipt(self):
        stage = self.run_interrupted(
            'import resource\nresource.setrlimit(resource.RLIMIT_CPU,(1,2))\nwhile True: pass', 5)
        self.assertEqual(stage['reason'], 'SIGNAL: SIGXCPU')

    def test_outer_deadline_leaves_inconclusive_receipt(self):
        stage = self.run_interrupted(
            'import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(5)', 0.2)
        self.assertEqual(stage['reason'], 'WALL-TIMEOUT')


if __name__ == '__main__':
    unittest.main()
