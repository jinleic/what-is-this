"""Focused contracts for the manifest-bound Kobon UNSAT certifier."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import lzma
import os
import signal
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

N12_DIR = Path(__file__).resolve().parents[1] / "n12_decide"
if str(N12_DIR) not in sys.path:
    sys.path.insert(0, str(N12_DIR))

import certify_unsat  # noqa: E402


def _write_executable(path: Path, source: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{source}", encoding="utf-8")
    path.chmod(0o755)


def _file_record(path: Path, workspace: Path) -> dict:
    payload = path.read_bytes()
    try:
        name = str(path.resolve().relative_to(workspace.resolve()))
    except ValueError:
        name = str(path.resolve())
    return {
        "bytes": len(payload),
        "path": name,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


class CertificateFixture:
    def __init__(self, root: Path):
        self.workspace = root
        self.here = root / "math/kobon/n12_decide"
        self.here.mkdir(parents=True)
        self.bin = root / "tools"
        self.counter = root / "solver-count.txt"
        self.pid_file = root / "slow-pids.txt"
        self.cake_marker = root / "cake-active.txt"
        self.solver = self.bin / "kissat"
        self.xz = self.bin / "xz"
        self.checker = self.bin / "drat-trim"
        self.pgrep = self.bin / "pgrep"
        self.ps = self.bin / "ps"
        self.driver_lock = self.here / ".certify-unsat.lock"
        self.host_lock = root / "scratch/.host-heavy-job.lock"
        self.cases: dict[str, Path] = {}
        self._install_tools()
        self.add_case("proof", "PROOF")
        self.rebuild_manifest()

    def _install_tools(self) -> None:
        solver_source = f'''import lzma
import os
from pathlib import Path
import subprocess
import sys
import time
if sys.argv[1:] == ["--version"]:
    print("fake-kissat 1.0")
    raise SystemExit(0)
with open({str(self.counter)!r}, "a", encoding="ascii") as stream:
    stream.write("launch\\n")
    stream.flush()
    os.fsync(stream.fileno())
cnf = Path(sys.argv[1]).read_text(encoding="ascii")
proof = Path(sys.argv[2])
if "MODE SAT" in cnf:
    print("s SATISFIABLE", flush=True)
    print("v 1 0", flush=True)
    raise SystemExit(10)
if "MODE SLOW" in cnf:
    proof.write_bytes(b"partial-xz")
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    Path({str(self.pid_file)!r}).write_text(
        f"{{os.getpid()}} {{child.pid}}\\n", encoding="ascii"
    )
    time.sleep(120)
payload = b"GOOD"
if "MODE BADSTREAM" in cnf:
    payload = b"BADSTREAM"
elif "MODE CHECKER_EXIT" in cnf:
    payload = b"CHECKER_EXIT"
with lzma.open(proof, "wb", format=lzma.FORMAT_XZ) as stream:
    stream.write(payload)
print("s UNSATISFIABLE", flush=True)
raise SystemExit(20)
'''
        xz_source = '''import lzma
import os
from pathlib import Path
import sys
expected = {"LANG", "LC_ALL", "PATH", "TZ"}
observed = set(os.environ)
# macOS injects this key when launching a Python executable even from an
# explicit four-key env; it is not inherited from the certifier.
observed.discard("__CF_USER_TEXT_ENCODING")
if observed != expected:
    print("unsafe environment: " + repr(sorted(os.environ)), file=sys.stderr)
    raise SystemExit(7)
if sys.argv[1:] == ["--version"]:
    print("xz fake 1.0")
    raise SystemExit(0)
if "--format=xz" not in sys.argv:
    raise SystemExit(8)
path = Path(sys.argv[-1])
try:
    with lzma.open(path, "rb", format=lzma.FORMAT_XZ) as stream:
        payload = stream.read()
except Exception as error:
    print(f"invalid xz: {error}", file=sys.stderr)
    raise SystemExit(1)
if "--test" in sys.argv:
    raise SystemExit(0)
if "--decompress" not in sys.argv or "--stdout" not in sys.argv:
    raise SystemExit(2)
sys.stdout.buffer.write(payload)
sys.stdout.buffer.flush()
raise SystemExit(9 if payload == b"BADSTREAM" else 0)
'''
        checker_source = '''from pathlib import Path
import sys
if "-i" not in sys.argv or "-w" not in sys.argv:
    print("s NOT VERIFIED")
    raise SystemExit(2)
cnf = Path(sys.argv[1]).read_text(encoding="ascii")
proof = sys.stdin.buffer.read()
if "MODE UP" in cnf:
    print("c input unit propagation conflict")
    print("s VERIFIED")
    raise SystemExit(0)
if proof == b"GOOD":
    print("c consumed binary proof")
    print("s VERIFIED")
    raise SystemExit(0)
if proof == b"CHECKER_EXIT":
    print("s VERIFIED")
    raise SystemExit(9)
print("s NOT VERIFIED")
raise SystemExit(1)
'''
        pgrep_source = f'''from pathlib import Path
marker = Path({str(self.cake_marker)!r})
if marker.exists():
    print(marker.read_text(encoding="ascii").strip() or "99999")
    raise SystemExit(0)
raise SystemExit(1)
'''
        ps_source = '''import os
print(f"{os.getpid()} {os.getpgrp()} 1")
'''
        _write_executable(self.solver, solver_source)
        _write_executable(self.xz, xz_source)
        _write_executable(self.checker, checker_source)
        _write_executable(self.pgrep, pgrep_source)
        _write_executable(self.ps, ps_source)

    def add_case(self, name: str, mode: str) -> None:
        cnf_dir = self.workspace / "cnf"
        cnf_dir.mkdir(exist_ok=True)
        path = cnf_dir / f"{name}.cnf"
        if mode == "UP":
            body = f"c MODE {mode}\np cnf 1 2\n1 0\n-1 0\n"
        else:
            body = f"c MODE {mode}\np cnf 1 1\n1 0\n"
        path.write_text(body, encoding="ascii")
        self.cases[name] = path

    def rebuild_manifest(self) -> None:
        self.document = self._manifest_document()
        self.manifest_text = json.dumps(
            self.document, indent=2, sort_keys=True
        ) + "\n"
        self.manifest_path = self.here / "manifest.json"
        self.manifest_path.write_text(self.manifest_text, encoding="utf-8")
        self.digest = hashlib.sha256(self.manifest_text.encode()).hexdigest()

    def _manifest_document(self) -> dict:
        instances = []
        for name, path in self.cases.items():
            record = _file_record(path, self.workspace)
            record.update({
                "clauses": 2 if "MODE UP" in path.read_text() else 1,
                "cube": name,
                "variables": 1,
            })
            instances.append(record)
        solver = _file_record(self.solver, self.workspace)
        solver.update({
            "path": str(self.solver.resolve()),
            "version": "fake-kissat 1.0",
        })
        return {
            "case_count": len(instances),
            "case_order": [item["cube"] for item in instances],
            "instances": instances,
            "schema": certify_unsat.DECISION_SCHEMA,
            "solver": solver,
        }

    def certifier(self, **overrides) -> certify_unsat.Certifier:
        config_values = {
            "workspace": self.workspace,
            "here": self.here,
            "manifest_path": self.manifest_path,
            "checker_path": self.checker,
            "driver_lock_path": self.driver_lock,
            "host_lock_path": self.host_lock,
            "expected_decision_sha256": self.digest,
            "expected_cases": tuple(self.cases),
            "xz_path": self.xz,
            "pgrep_path": self.pgrep,
            "ps_path": self.ps,
            "disk_floor_bytes": 1,
            "disk_poll_seconds": 0.02,
            "process_poll_seconds": 0.01,
            "memory_limit_bytes": None,
            "resource_wait_seconds": 0.2,
            "kill_grace_seconds": 0.2,
        }
        constructor = {
            "manifest_builder": lambda: copy.deepcopy(self.document),
        }
        for key, value in overrides.items():
            if key in certify_unsat.Config.__dataclass_fields__:
                config_values[key] = value
            else:
                constructor[key] = value
        return certify_unsat.Certifier(
            certify_unsat.Config(**config_values), **constructor
        )

    def launch_count(self) -> int:
        if not self.counter.exists():
            return 0
        return len(self.counter.read_text(encoding="ascii").splitlines())


class CertifyUnsatTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.fixture = CertificateFixture(Path(self.temporary.name))

    def _run_success(self, case: str = "proof"):
        certifier = self.fixture.certifier()
        self.assertEqual(certifier.run([case]), certify_unsat.EXIT_OK)
        ledger = json.loads(certifier.ledger_path.read_bytes())
        self.assertEqual(len(ledger["rows"]), 1)
        return certifier, ledger["rows"][0]

    def test_stale_manifest_and_nonmanifest_cases_refuse_before_solver(self):
        self.fixture.manifest_path.write_text(
            self.fixture.manifest_text + " ", encoding="utf-8"
        )
        with self.assertRaisesRegex(certify_unsat.Refused, "stale"):
            self.fixture.certifier().run(["proof"])
        self.fixture.manifest_path.write_text(
            self.fixture.manifest_text, encoding="utf-8"
        )
        with self.assertRaisesRegex(certify_unsat.Refused, "case_order"):
            self.fixture.certifier().run(["../cnf/proof.cnf"])
        with self.assertRaisesRegex(certify_unsat.Refused, "duplicate"):
            self.fixture.certifier().run(["proof", "proof"])
        self.assertEqual(self.fixture.launch_count(), 0)

    def test_corrupt_foreign_and_tampered_ledgers_are_refused(self):
        certifier, row = self._run_success()
        ledger_path = certifier.ledger_path
        ledger_path.write_bytes(b"{not-json\n")
        with self.assertRaisesRegex(certify_unsat.Refused, "corrupt"):
            self.fixture.certifier().run([], status=True)
        ledger = certifier._new_ledger()
        foreign = copy.deepcopy(row)
        foreign["case"] = "not-in-manifest"
        ledger["rows"] = [foreign]
        ledger_path.write_bytes(certify_unsat.canonical_json(ledger))
        with self.assertRaisesRegex(certify_unsat.Refused, "foreign case"):
            self.fixture.certifier().run([], status=True)
        ledger["rows"] = [row]
        ledger_path.write_bytes(certify_unsat.canonical_json(ledger))
        proof = certifier._artifact_path(row["proof"])
        proof.write_bytes(proof.read_bytes() + b"tamper")
        with self.assertRaisesRegex(certify_unsat.Refused, "artifact bytes"):
            self.fixture.certifier().run([], status=True)

    def test_duplicate_attempt_rows_are_refused(self):
        certifier, row = self._run_success()
        ledger = certifier._new_ledger()
        ledger["rows"] = [row, copy.deepcopy(row)]
        certifier.ledger_path.write_bytes(certify_unsat.canonical_json(ledger))
        with self.assertRaisesRegex(certify_unsat.Refused, "duplicate attempt"):
            self.fixture.certifier().run([], status=True)

    def test_resume_and_status_freshly_recheck_without_kissat(self):
        first, row = self._run_success()
        self.assertEqual(row["outcome"], "CERTIFIED_UNSAT")
        second = self.fixture.certifier()
        self.assertEqual(second.run(["proof"]), certify_unsat.EXIT_OK)
        output = io.StringIO()
        status = self.fixture.certifier(out=output)
        self.assertEqual(status.run([], status=True), certify_unsat.EXIT_OK)
        self.assertIn("FRESHLY_REVERIFIED_TERMINAL", output.getvalue())
        self.assertEqual(self.fixture.launch_count(), 1)
        persisted = json.loads(first.ledger_path.read_bytes())
        self.assertEqual(len(persisted["rows"]), 1)
        self.assertGreaterEqual(len(list(first.recheck_dir.glob("*.json"))), 2)

    def test_status_rejects_self_consistent_but_forged_logged_claim(self):
        certifier, row = self._run_success()
        proof = certifier._artifact_path(row["proof"])
        proof.write_bytes(lzma.compress(b"FORGED", format=lzma.FORMAT_XZ))
        row["proof"] = certify_unsat.file_record(
            proof, self.fixture.workspace
        )
        ledger = certifier._new_ledger()
        ledger["rows"] = [row]
        certifier.ledger_path.write_bytes(certify_unsat.canonical_json(ledger))
        with self.assertRaisesRegex(certify_unsat.Refused, "changed evidence"):
            self.fixture.certifier().run([], status=True)
        self.assertEqual(self.fixture.launch_count(), 1)

    def test_binary_xz_pipeline_requires_both_process_exits(self):
        self.fixture.add_case("bad_xz", "BADSTREAM")
        self.fixture.add_case("bad_checker", "CHECKER_EXIT")
        self.fixture.rebuild_manifest()
        certifier = self.fixture.certifier()
        self.assertEqual(
            certifier.run(["bad_xz"]), certify_unsat.EXIT_INCOMPLETE
        )
        row = json.loads(certifier.ledger_path.read_bytes())["rows"][-1]
        self.assertEqual(row["outcome"], "XZ_STREAM_FAILED")
        self.assertNotEqual(row["proof_check"]["xz_exit_code"], 0)
        self.assertEqual(
            certifier._artifact_path(row["proof"]).read_bytes()[:6],
            bytes.fromhex("fd377a585a00"),
        )
        self.assertIn("--format=xz", row["proof_check"]["xz_argv"])
        certifier = self.fixture.certifier()
        self.assertEqual(
            certifier.run(["bad_checker"]), certify_unsat.EXIT_INCOMPLETE
        )
        row = json.loads(certifier.ledger_path.read_bytes())["rows"][-1]
        self.assertEqual(row["outcome"], "CHECKER_EXIT_FAILED")
        self.assertEqual(row["proof_check"]["xz_exit_code"], 0)
        self.assertNotEqual(row["proof_check"]["checker_exit_code"], 0)
        self.assertEqual(row["proof_check"]["checker_argv"][-2:], ["-i", "-w"])

    def test_input_up_has_a_distinct_evidence_class(self):
        self.fixture.add_case("input_up", "UP")
        self.fixture.rebuild_manifest()
        certifier = self.fixture.certifier()
        self.assertEqual(
            certifier.run(["proof", "input_up"]), certify_unsat.EXIT_OK
        )
        rows = {
            row["case"]: row
            for row in json.loads(certifier.ledger_path.read_bytes())["rows"]
        }
        self.assertEqual(rows["proof"]["outcome"], "CERTIFIED_UNSAT")
        self.assertEqual(
            rows["proof"]["evidence_class"],
            "PROOF_CONSUMING_BINARY_DRAT_XZ",
        )
        self.assertEqual(
            rows["input_up"]["outcome"], "VERIFIED_INPUT_UP_UNSAT"
        )
        self.assertEqual(
            rows["input_up"]["evidence_class"],
            "INPUT_UNIT_PROPAGATION_UNSAT",
        )

    def test_environment_is_minimal_hashed_and_immune_to_xz_poisoning(self):
        with mock.patch.dict(
            os.environ,
            {
                "XZ_DEFAULTS": "--format=raw",
                "XZ_OPT": "--format=raw",
                "DYLD_INSERT_LIBRARIES": "/bad",
                "LD_PRELOAD": "/bad",
            },
            clear=False,
        ):
            _, row = self._run_success()
        expected = {
            "LANG": "C",
            "LC_ALL": "C",
            "PATH": str(self.fixture.xz.parent.resolve()),
            "TZ": "UTC",
        }
        self.assertEqual(row["child_environment"], expected)
        self.assertEqual(
            row["child_environment_sha256"],
            hashlib.sha256(certify_unsat.canonical_json(expected)).hexdigest(),
        )

    def test_nonterminal_requires_explicit_retry(self):
        self.fixture.add_case("sat", "SAT")
        self.fixture.rebuild_manifest()
        first = self.fixture.certifier()
        self.assertEqual(first.run(["sat"]), certify_unsat.EXIT_INCOMPLETE)
        launches = self.fixture.launch_count()
        with self.assertRaisesRegex(certify_unsat.Refused, "retry-nonterminal"):
            self.fixture.certifier().run(["sat"])
        self.assertEqual(self.fixture.launch_count(), launches)
        self.assertEqual(
            self.fixture.certifier().run(
                ["sat"], retry_nonterminal=True
            ),
            certify_unsat.EXIT_INCOMPLETE,
        )
        self.assertEqual(self.fixture.launch_count(), launches + 1)

    def test_abandoned_start_proof_and_log_are_ledgered_as_orphan(self):
        abandoned = self.fixture.certifier()
        abandoned.prepare()
        case, attempt = "proof", "abandoned.12345678"
        paths = abandoned._paths(case, attempt)
        started = time.time_ns()
        argv = [
            abandoned.tools["kissat"]["path"],
            abandoned.bindings[case]["resolved_path"],
            str(paths["proof"]),
        ]
        abandoned._write_start(
            case, attempt, started, argv, paths, abandoned._options()
        )
        paths["proof"].write_bytes(b"partial immutable proof")
        paths["solver"].write_bytes(b"partial immutable log")
        output = io.StringIO()
        recovered = self.fixture.certifier(out=output)
        self.assertEqual(recovered.run([], status=True), certify_unsat.EXIT_OK)
        rows = json.loads(recovered.ledger_path.read_bytes())["rows"]
        self.assertEqual(rows[-1]["outcome"], "ORPHANED_ATTEMPT")
        self.assertEqual(len(rows[-1]["orphan_artifacts"]), 3)
        self.assertTrue(all(
            certify_unsat._SHA_RE.fullmatch(record["sha256"])
            for record in rows[-1]["orphan_artifacts"]
        ))
        with self.assertRaisesRegex(certify_unsat.Refused, "retry-nonterminal"):
            self.fixture.certifier().run(["proof"])

    def test_driver_and_host_lock_fds_are_inherited(self):
        certifier = self.fixture.certifier()
        certifier.driver_lock.acquire()
        self.addCleanup(certifier.driver_lock.release)
        host = certify_unsat.FileLock(self.fixture.host_lock)
        host.acquire()
        certifier.host_lock = host
        self.addCleanup(host.release)
        script = (
            "import os,sys; "
            "[os.fstat(int(value)) for value in sys.argv[1:]]; print('held')"
        )
        process = certifier._spawn(
            [
                sys.executable, "-c", script,
                str(certifier.driver_lock.fd), str(host.fd),
            ],
            stdout=certify_unsat.subprocess.PIPE,
            stderr=certify_unsat.subprocess.PIPE,
        )
        stdout, _ = process.communicate(timeout=2)
        certifier._complete(process)
        self.assertEqual(process.returncode, 0)
        self.assertEqual(stdout.strip(), b"held")

    def test_driver_lock_excludes_a_second_certifier(self):
        holder = certify_unsat.FileLock(self.fixture.driver_lock)
        holder.acquire()
        self.addCleanup(holder.release)
        with self.assertRaisesRegex(certify_unsat.LockBusy, "lock held"):
            self.fixture.certifier().run(["proof"])
        self.assertEqual(self.fixture.launch_count(), 0)
        self.assertTrue(self.fixture.driver_lock.exists())

    def test_host_lock_rechecks_cake_before_verification(self):
        certifier = self.fixture.certifier(resource_wait_seconds=0.03)
        certifier.prepare()
        observed = iter(([12345], []))
        certifier._cake_processes = lambda: list(next(observed))
        with certifier._host_exclusion():
            competitor = certify_unsat.FileLock(self.fixture.host_lock)
            with self.assertRaises(certify_unsat.LockBusy):
                competitor.acquire()
        self.assertEqual(self.fixture.launch_count(), 0)

    def test_host_probe_failure_reaps_before_unlock(self):
        certifier = self.fixture.certifier()
        certifier.prepare()

        def assert_still_locked() -> None:
            self.assertIsNotNone(certifier.host_lock)
            competitor = certify_unsat.FileLock(self.fixture.host_lock)
            with self.assertRaises(certify_unsat.LockBusy):
                competitor.acquire()

        with (
            mock.patch.object(
                certifier, "_cake_processes",
                side_effect=certify_unsat.Refused("probe failed")),
            mock.patch.object(
                certifier, "_terminate_all",
                side_effect=assert_still_locked) as terminate,
        ):
            with self.assertRaisesRegex(certify_unsat.Refused, "probe failed"):
                with certifier._host_exclusion():
                    self.fail("probe failure must prevent verifier entry")
        terminate.assert_called_once_with()
        with certify_unsat.FileLock(self.fixture.host_lock):
            pass

    def test_disk_floor_kills_and_reaps_solver_process_group(self):
        self.fixture.add_case("slow", "SLOW")
        self.fixture.rebuild_manifest()

        def free_space(_path: Path) -> int:
            return 0 if self.fixture.pid_file.exists() else 1 << 40

        certifier = self.fixture.certifier(
            disk_free=free_space,
            disk_floor_bytes=1 << 30,
            disk_poll_seconds=0.02,
        )
        self.assertEqual(certifier.run(["slow"]), certify_unsat.EXIT_INCOMPLETE)
        row = json.loads(certifier.ledger_path.read_bytes())["rows"][-1]
        self.assertEqual(row["outcome"], "DISK_FLOOR_ABORT")
        self.assertEqual(row["solver"]["stop_reason"], "DISK_FLOOR")
        self.assertEqual(certifier.active, {})
        pids = [int(value) for value in self.fixture.pid_file.read_text().split()]
        deadline = time.monotonic() + 3
        for pid in pids:
            while time.monotonic() < deadline:
                try:
                    os.kill(pid, 0)
                except ProcessLookupError:
                    break
                time.sleep(0.01)
            else:
                self.fail(f"process-group member {pid} survived disk-floor abort")

    def test_default_resource_contract_is_eight_gib(self):
        self.assertEqual(certify_unsat.DEFAULT_MEMORY_LIMIT_GIB, 8.0)
        self.assertEqual(
            certify_unsat.default_config().memory_limit_bytes, 8 << 30
        )

    def test_default_disk_floor_launches_with_current_headroom(self):
        self.assertEqual(certify_unsat.DEFAULT_DISK_FLOOR_GIB, 64.0)
        current_headroom = 123 << 30
        certifier = self.fixture.certifier(
            disk_floor_bytes=certify_unsat.default_config().disk_floor_bytes,
            disk_free=lambda _path: current_headroom,
        )
        self.assertEqual(certifier.run(["proof"]), certify_unsat.EXIT_OK)
        row = json.loads(certifier.ledger_path.read_bytes())["rows"][-1]
        self.assertEqual(row["outcome"], "CERTIFIED_UNSAT")
        self.assertGreater(current_headroom, row["options"]["disk_floor_bytes"])

    def test_eight_gib_limit_launches_a_real_child(self):
        certifier = self.fixture.certifier(memory_limit_bytes=8 << 30)
        self.assertEqual(certifier.run(["proof"]), certify_unsat.EXIT_OK)
        row = json.loads(certifier.ledger_path.read_bytes())["rows"][-1]
        self.assertEqual(row["outcome"], "CERTIFIED_UNSAT")
        self.assertEqual(row["solver"]["exit_code"], 20)
        self.assertIsNone(row["solver"]["launch_error"])

    def test_darwin_uses_aggregate_rss_without_broken_preexec(self):
        with mock.patch.object(certify_unsat.sys, "platform", "darwin"):
            self.assertIsNone(certify_unsat._memory_limit(8 << 30))


if __name__ == "__main__":
    unittest.main()
