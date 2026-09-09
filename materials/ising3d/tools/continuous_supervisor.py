#!/usr/bin/env python3
"""continuous_supervisor.py -- persistent, checkpointed exact-research coordinator for math/ising3d.

One long-lived process owns one campaign claim and drives, until a safe stop,

    plan turn -> fresh independent review/verifier-author turn
      -> producer job -> verifier job -> checkpoint + control-plane update -> plan turn ...
    Every model and numerical subprocess uses tools/resource_guard.py, sequentially.

Invocation (Main launches it detached through hub; use absolute paths):

    math/.venv/bin/python -I -B math/ising3d/tools/continuous_supervisor.py run \
        --seed-plan /abs/seed.json [--max-cycles N] [--max-turns N] [--plan-only] [--resume-halted]
    ... status | check [--seed-plan F] | validate-plan FILE [--seed] | render-prompt --kind plan|review

Seed plan file: the plan schema (tools/supervisor_plan_schema.json) plus a "verifier" object in
the review-schema shape. producer null + existing artifact + independent verifier runs only the
verifier (H681: e254 artifact + tests/test_trace_nine_bounded_canary.py). A seed is consumed once.

Safeguards, each enforced in code rather than by prompt text:
  1. one advisory flock (tools/.continuous_supervisor.lock); a second instance exits 0 untouched
  2. campaign lifecycle through scripts/campaign.py core APIs in-process (claim.pid is this
     process); adopts the single live campaign carrying its own agent label, never inits over a
     live claim or run dir, never closes; heartbeats every --heartbeat-seconds; the prereg bytes
     are snapshotted into <run>/prereg_snapshot.txt right after init
  3. every numeric process is `python -I -B <script> <args>` from math/ising3d under the guard;
     acceptance = guard reason 'completed' AND guard/child exit 0 AND artifact/verifier checks;
     a printed PASS never overrides a failed guard verdict
  4. never restarts a healthy worker: guard busy (exit 75, single worker slot) waits; SIGTERM to
     the supervisor detaches (the job keeps its own session) and the next start re-adopts the pid;
     a job is never signalled unless --job-kill-after-seconds is set explicitly
  5. completed trial identities pin source bytes and mathematical input arguments; output
     renaming and unrelated file changes cannot authorize retries. Resource-blocked norm and
     W-law tasks must advance to W-law and Callen respectively.
  6. agent turns run `claude -p --restricted --safe-mode --strict-mcp-config` with
     Read/Glob/Grep/Edit/Write only, cwd math/ising3d, every frozen path denied, and are
     validated afterwards (frozen hashes, append-only ledgers, no new Markdown, schema, paths)
  7. 50 GiB free reserve on the workspace and /tmp, a cumulative write allowance for supervisor
     files, bounded prompt/response/log/transition files; a set-user-ID descendant that cannot be
     measured is reported and allowed at most 30 s of life, and teardown is bounded so an
     unkillable member can never hang the guard; either anomaly HALTS as UNACCOUNTABLE without
     writing the trial into history; any breach is HALTED with a reason and no eviction,
     deletion, or loosened limit
  8. agent/auth failures back off exponentially (bounded) and are recorded as unavailable
     planning; a rejected model quota window is waited out instead of counted as a failure; a CLI
     incompatibility halts immediately; no mathematical claim is fabricated
  9. an interrupted agent turn is audited on restart and halts with evidence when its
     pre-snapshot is unreadable; a claim is only respected when its pid is really a live
     supervisor, so a reused pid can never strand the campaign
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import hashlib
import io
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path, PurePosixPath

sys.dont_write_bytecode = True

TOOLS_DIR = Path(__file__).resolve().parent
TARGET_DIR = TOOLS_DIR.parent
DOMAIN_DIR = TARGET_DIR.parent
ROOT = DOMAIN_DIR.parent
DOMAIN, TARGET = DOMAIN_DIR.name, TARGET_DIR.name
IDENT = f"{DOMAIN}/{TARGET}"
sys.path.insert(0, str(ROOT / "scripts"))
import campaign  # noqa: E402  canonical control plane (stdlib only)

GUARD_PATH = TOOLS_DIR / "resource_guard.py"
PLAN_SCHEMA_PATH = TOOLS_DIR / "supervisor_plan_schema.json"
REVIEW_SCHEMA_PATH = TOOLS_DIR / "supervisor_review_schema.json"
LOCK_PATH = TOOLS_DIR / ".continuous_supervisor.lock"
INTERPRETER = TARGET_DIR / ".venv" / "bin" / "python"
CLAUDE_DEFAULT = Path.home() / ".local" / "bin" / "claude"
STATUS_SCHEMA = "ising3d-continuous-supervisor-status/1"
SUPERVISOR_DIRNAME = "supervisor"
STATUS_NAME = "supervisor_status.json"
TRANSITIONS_NAME = "transitions.jsonl"

GIB = 1024 ** 3
MIB = 1024 ** 2
FREE_RESERVE_DEFAULT = 50 * GIB
ARTIFACT_MAX_BYTES = 16 * MIB
WRITE_ALLOWANCE_DEFAULT = 64 * MIB
TRANSITIONS_MAX_BYTES = 8 * MIB
PROMPT_MAX_BYTES = 96 * 1024
RESPONSE_MAX_BYTES = 12 * MIB  # must stay under the guard's 16 MiB per-file cap
STDERR_MAX_BYTES = 64 * 1024
LOG_TAIL_BYTES = 4096
GUARD_BUSY_MARK = "single worker slot"
GUARD_UNMEASURABLE_MARK = "unmeasurable descendant"
GUARD_CLEANUP_MARK = "survived its bounded cleanup"
READY_MARK = "RESOURCE_GUARD READY"
HYP_RE = re.compile(r"^H[0-9]{3,4}$")
FATAL_CLI_MARKS = ("unknown option", "too many arguments", "not a valid JSON Schema",
                   "error: option", "strict mode")
QUOTA_WAIT_MIN_SECONDS = 60.0
QUOTA_WAIT_MAX_SECONDS = 8 * 3600.0

# Files the agent may edit; everything else that pre-exists is denied and hash-checked.
LEDGER_APPEND_ONLY = ("notes/hypotheses.csv", "checkpoints/failed_routes.md", "research_log.md")
EDITABLE_CHECKPOINTS = ("checkpoints/current_state.md", "checkpoints/next_actions.md",
                        "checkpoints/failed_routes.md")
NEW_FILE_ROOTS = ("experiments/", "tests/")
FROZEN_DIRS = ("tools", "results", "src", "deliverables", "sources", "campaigns", "reports", "proofs")
# Explicitly named frozen inputs (also covered by the enumerated deny list).
FROZEN_EXPLICIT = (
    "experiments/e248_replica_trace_nine.py", "experiments/e251_trace_nine_projection.py",
    "experiments/e252_trace_nine_norm_envelope.py",
    "experiments/e253_trace_nine_structured_determinant.py",
    "tests/test_replica_trace_nine.py", "tests/test_trace_nine_projection.py",
    "tests/test_trace_nine_norm_envelope.py", "tests/test_trace_nine_structured_determinant.py",
    "state.json", "README.md", "pyproject.toml", "uv.lock",
)
SNAPSHOT_SKIP_DIRS = {"__pycache__", ".git", ".venv"}
SNAPSHOT_SKIP_FILES = {".DS_Store"}

COMPLETED_OUTCOMES = {
    "CONFIRMED", "PRODUCER_RESOURCE_BLOCKED", "PRODUCER_FAILED", "ARTIFACT_INVALID",
    "REVIEW_REJECTED", "VERIFIER_RESOURCE_BLOCKED", "VERIFIER_FAILED", "JOB_INTEGRITY",
}
MAX_RETRYABLE_REPEATS = 2

LIMITS_TEXT = ("nice 19, one numerical thread, 35% of one CPU core by stop/continue pacing, "
               "monitored 2 GiB memory stop, 1800 CPU-seconds, 16 MiB per file, 64 MiB total "
               "writes, 50 GiB free-space reserve on the workspace and /tmp filesystems, and a "
               "30-second limit on any descendant the guard is not permitted to measure")

SYSTEM_APPEND = (
    "You are one bounded, fresh turn of a persistent exact-mathematics supervisor for the 3D "
    "Ising exact program in this directory. You cannot run code; a separate resource-guarded "
    "process executes independently reviewed scripts. Never delete anything, create Markdown "
    "files, commit, push, access unrelated personal data, change global configuration, or "
    "mutate infrastructure. Authored scripts must not spawn processes, access the network, "
    "delete files, modify existing files, or write outside their one declared artifact. "
    "Keep all numerical work single-threaded. Make only surgical file edits. Never edit "
    "protected paths. Report only what you have read; never invent computed results, digests, "
    "degrees, timings, or memory figures. Return exactly one schema-conforming JSON object."
)

# Exec replaces this tiny redirector: Claude remains inside the guard's owned
# process group, while its JSON cannot interleave with resource-control stdout.
AGENT_REDIRECT_CODE = """import os, resource, sys
limit, out, err, *command = sys.argv[1:]
resource.setrlimit(resource.RLIMIT_FSIZE, (int(limit), int(limit)))
for descriptor, path in ((1, out), (2, err)):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.dup2(fd, descriptor)
    os.close(fd)
os.execv(command[0], command)
"""


class SupervisorError(Exception):
    """Configuration or precondition error; CLI maps to exit code 3."""


# ----------------------------------------------------------------- utilities

def utc_now() -> str:
    return campaign.fmt_z(campaign.now_utc())


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return campaign.sha256_file(path)


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def tail_text(path: Path, limit: int = LOG_TAIL_BYTES) -> str:
    try:
        size = path.stat().st_size
        with path.open("rb") as fh:
            if size > limit:
                fh.seek(size - limit)
            data = fh.read(limit)
    except OSError:
        return ""
    text = data.decode("utf-8", "replace")
    return ("[...tail...]\n" + text) if size > limit else text


def excerpt(text: str, limit: int) -> str:
    data = text.encode("utf-8")
    if len(data) <= limit:
        return text
    head = data[:limit].decode("utf-8", "ignore")
    return f"{head}\n[... truncated {len(data) - limit} bytes ...]"


def md_section(text: str, keyword: str) -> str:
    """Body of the first '## ' heading containing keyword (case-insensitive), heading included."""
    lines = text.splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.startswith("## ") and keyword.lower() in line.lower():
            start = index
            break
    if start is None:
        return ""
    end = len(lines)
    for index in range(start + 1, len(lines)):
        if lines[index].startswith("## "):
            end = index
            break
    return "\n".join(lines[start:end]).strip()


def free_bytes() -> int:
    return min(shutil.disk_usage(TARGET_DIR).free, shutil.disk_usage("/tmp").free)


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def pid_command(pid: int) -> str:
    try:
        out = subprocess.run(["ps", "-o", "command=", "-p", str(pid)], capture_output=True,
                             text=True, timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return ""
    return out.stdout.strip()


def command_runs_this_script(command: str) -> bool:
    """Whole-path-token match, so a test module, editor, or grep is not a supervisor."""
    script = Path(__file__).name
    return any(PurePosixPath(token).name == script for token in command.split())


def supervisor_alive(pid: int) -> bool:
    """True only when that pid is really another supervisor, not a reused number.

    A false negative would let a second supervisor adopt a live run, so the match
    stays permissive; the exclusive flock taken before this check is the backstop.
    """
    return pid_alive(pid) and command_runs_this_script(pid_command(pid))


def quota_reset_unix(stdout: bytes) -> float | None:
    """Unix time a rejected model quota window reopens, else None.

    The CLI streams a rate_limit_event before failing; a rejected window is a
    wait, not an agent fault, so it must never count against the failure budget.
    """
    try:
        events = json.loads(stdout or b"null")
    except ValueError:
        return None
    resets: list[float] = []
    for event in events if isinstance(events, list) else [events]:
        if not isinstance(event, dict) or event.get("type") != "rate_limit_event":
            continue
        info = event.get("rate_limit_info")
        if not isinstance(info, dict):
            continue
        if "rejected" not in (str(info.get("status")), str(info.get("overageStatus"))):
            continue
        with contextlib.suppress(TypeError, ValueError):
            resets.append(float(info["resetsAt"]))
    return max(resets) if resets else None


def insert_note(text: str, note: str) -> str:
    """Place a checkpoint note where newest-first readers actually look.

    Appending at EOF hid every outcome from the planner: context_blocks() reads
    only the head of current_state.md, and PROGRESS.md is newest-first.
    """
    if note in text:
        return text
    marks = [index for index in (text.find("\n## "), text.find("\n### ")) if index > 0]
    if not marks:
        return text.rstrip("\n") + "\n" + note
    boundary = min(marks) + 1
    return text[:boundary] + note.lstrip("\n") + "\n" + text[boundary:]


def get_dotted(doc, dotted: str):
    node = doc
    for part in dotted.split("."):
        if isinstance(node, list):
            if not part.isdigit() or int(part) >= len(node):
                raise KeyError(dotted)
            node = node[int(part)]
        elif isinstance(node, dict):
            if part not in node:
                raise KeyError(dotted)
            node = node[part]
        else:
            raise KeyError(dotted)
    return node


def scalar_equal(actual, expected) -> bool:
    if isinstance(actual, bool) or isinstance(expected, bool):
        return actual is expected
    if isinstance(actual, (int, float)) and isinstance(expected, (int, float)):
        return actual == expected
    return type(actual) is type(expected) and actual == expected


def hypothesis_number(hyp: str) -> int:
    return int(hyp[1:])


def read_lock_holder() -> dict | None:
    if not LOCK_PATH.exists():
        return None
    with contextlib.suppress(OSError, ValueError):
        holder = json.loads(LOCK_PATH.read_text(encoding="utf-8") or "null")
        return holder if isinstance(holder, dict) else None
    return None


# ------------------------------------------------------ minimal JSON schema

_TYPE_CHECK = {
    "object": lambda v: isinstance(v, dict),
    "array": lambda v: isinstance(v, list),
    "string": lambda v: isinstance(v, str),
    "integer": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "number": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "boolean": lambda v: isinstance(v, bool),
    "null": lambda v: v is None,
}


def schema_errors(value, schema: dict, where: str = "$") -> list[str]:
    """Validate the subset of JSON Schema used by the supervisor schemas (no $ref)."""
    errors: list[str] = []
    types = schema.get("type")
    if types is not None:
        allowed = [types] if isinstance(types, str) else list(types)
        if not any(_TYPE_CHECK[t](value) for t in allowed):
            return [f"{where}: expected type {'/'.join(allowed)}"]
    if "anyOf" in schema and all(schema_errors(value, option, where)
                                 for option in schema["anyOf"]):
        # The CLI validator enforces anyOf; ignoring it here would silently let
        # a locally built plan carry a shape the model could never return.
        return [f"{where}: matches none of the allowed forms"]
    if value is None:
        return errors
    if "enum" in schema and value not in schema["enum"]:
        errors.append(f"{where}: not one of {schema['enum']}")
    if isinstance(value, str):
        if "minLength" in schema and len(value) < schema["minLength"]:
            errors.append(f"{where}: shorter than {schema['minLength']}")
        if "maxLength" in schema and len(value) > schema["maxLength"]:
            errors.append(f"{where}: longer than {schema['maxLength']}")
        if "pattern" in schema and not re.search(schema["pattern"], value):
            errors.append(f"{where}: does not match {schema['pattern']}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            errors.append(f"{where}: below {schema['minimum']}")
        if "maximum" in schema and value > schema["maximum"]:
            errors.append(f"{where}: above {schema['maximum']}")
    if isinstance(value, list):
        if "minItems" in schema and len(value) < schema["minItems"]:
            errors.append(f"{where}: fewer than {schema['minItems']} items")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            errors.append(f"{where}: more than {schema['maxItems']} items")
        if "items" in schema:
            for index, item in enumerate(value):
                errors += schema_errors(item, schema["items"], f"{where}[{index}]")
    if isinstance(value, dict):
        props = schema.get("properties", {})
        for key in schema.get("required", []):
            if key not in value:
                errors.append(f"{where}: missing required '{key}'")
        if "maxProperties" in schema and len(value) > schema["maxProperties"]:
            errors.append(f"{where}: more than {schema['maxProperties']} properties")
        extra = schema.get("additionalProperties", True)
        for key, item in value.items():
            if key in props:
                errors += schema_errors(item, props[key], f"{where}.{key}")
            elif extra is False:
                errors.append(f"{where}: unexpected property '{key}'")
            elif isinstance(extra, dict):
                errors += schema_errors(item, extra, f"{where}.{key}")
    return errors


# ---------------------------------------------------------------- supervisor

class Supervisor:
    def __init__(self, cfg: argparse.Namespace) -> None:
        self.cfg = cfg
        self.lock_fd: int | None = None
        self.run_dir: Path | None = None
        self.state_dir: Path | None = None
        self.status: dict = self.fresh_status()
        self.stop_requested = False
        self.exit_code = 0
        self.last_heartbeat = 0.0
        self.child: subprocess.Popen | None = None
        self.plan_schema = read_json(PLAN_SCHEMA_PATH)
        self.review_schema = read_json(REVIEW_SCHEMA_PATH)
        self.editable_extra = set(cfg.editable or [])
        self._halting = False
        # While a turn or job window is open, state.json writes are deferred so that
        # integrity diffs attribute every change to the agent or the job alone.
        self.quiet_depth = 0
        self.deferred_state_note: str | None = None

    # ---------------------------------------------------------- status IO
    def fresh_status(self) -> dict:
        return {
            "schema": STATUS_SCHEMA,
            "phase": "INIT",
            "halt_reason": None,
            "halt_detail": None,
            "supervisor_pid": os.getpid(),
            "started_utc": utc_now(),
            "updated_utc": utc_now(),
            "campaign_run_id": None,
            "cycle": 0,
            "turns": 0,
            "jobs": 0,
            "bytes_written": 0,
            "write_allowance_bytes": self.cfg.write_allowance_bytes,
            "seed": None,
            "current_task": None,
            "active_job": None,
            "last_outcome": None,
            "history": [],
            "agent_created_files": [],
            "quota_resets_utc": None,
            "agent": {"consecutive_failures": 0, "next_attempt_utc": None, "last_error": None},
        }

    def status_path(self) -> Path:
        assert self.state_dir is not None
        return self.state_dir / STATUS_NAME

    def save_status(self) -> None:
        if self.state_dir is None:
            return
        self.status["updated_utc"] = utc_now()
        self.status["supervisor_pid"] = os.getpid()
        self.write_json(self.status_path(), self.status, "status")

    def load_status(self) -> None:
        path = self.status_path()
        if path.is_file():
            loaded = read_json(path)
            if loaded.get("schema") != STATUS_SCHEMA:
                raise SupervisorError(f"{path}: unexpected schema {loaded.get('schema')!r}")
            base = self.fresh_status()
            base.update(loaded)
            base["write_allowance_bytes"] = self.cfg.write_allowance_bytes
            self.status = base
        else:
            self.status = self.fresh_status()

    # ------------------------------------------------- bounded evidence IO
    def charge(self, nbytes: int, what: str) -> None:
        self.status["bytes_written"] = int(self.status.get("bytes_written", 0)) + int(nbytes)
        if self.status["bytes_written"] > self.cfg.write_allowance_bytes and not self._halting:
            self.halt("WRITE_ALLOWANCE",
                      f"{self.status['bytes_written']} bytes written by the supervisor exceed "
                      f"the {self.cfg.write_allowance_bytes}-byte allowance (last: {what})")

    def write_text(self, path: Path, text: str, what: str) -> Path:
        data = text.encode("utf-8")
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.parent / f".{path.name}.{os.getpid()}.tmp"
        with tmp.open("xb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        self.charge(len(data), what)
        return path

    def write_json(self, path: Path, payload, what: str) -> Path:
        return self.write_text(path, json.dumps(payload, indent=1, sort_keys=True) + "\n", what)

    def event(self, event_type: str, **detail) -> None:
        record = {"utc": utc_now(), "cycle": self.status.get("cycle"), "event": event_type}
        record.update(detail)
        line = canonical(record) + "\n"
        print(f"SUPERVISOR {event_type} {canonical(detail)}", flush=True)
        if self.state_dir is None:
            return
        path = self.state_dir / TRANSITIONS_NAME
        size = path.stat().st_size if path.exists() else 0
        if size + len(line) > TRANSITIONS_MAX_BYTES:
            if not self._halting:
                self.halt("TRANSITIONS_BOUND",
                          f"{path} would exceed {TRANSITIONS_MAX_BYTES} bytes")
            return
        with path.open("ab") as fh:
            fh.write(line.encode("utf-8"))
        self.charge(len(line), "transitions")

    def halt(self, reason: str, detail: str) -> None:
        if self._halting:
            return
        self._halting = True
        self.status["phase"] = "HALTED"
        self.status["halt_reason"] = reason
        self.status["halt_detail"] = detail[:2000]
        self.exit_code = 2
        self.event("HALTED", reason=reason, detail=detail[:500])
        self.save_status()
        self.update_target_state(f"continuous_supervisor HALTED ({reason}): {detail}"[:600])
        self._halting = False

    # ------------------------------------------------------------ locking
    def lock(self) -> bool:
        LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            holder = ""
            with contextlib.suppress(OSError):
                holder = os.read(fd, 512).decode("utf-8", "replace").strip()
            os.close(fd)
            print(f"SUPERVISOR YIELD reason=lock_held holder={holder or '?'}", flush=True)
            return False
        os.ftruncate(fd, 0)
        os.lseek(fd, 0, os.SEEK_SET)
        os.write(fd, canonical({"pid": os.getpid(), "started_utc": utc_now()}).encode())
        self.lock_fd = fd
        return True

    # ----------------------------------------------------------- campaign
    @staticmethod
    def campaign_call(func, *args, **kwargs):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            result = func(*args, **kwargs)
        return result, buffer.getvalue().strip()

    def attach_campaign(self) -> bool:
        claim = campaign.claim_path(ROOT, DOMAIN, TARGET)
        lives = campaign.live_run_dirs(TARGET_DIR)
        mine = [run for run in lives
                if campaign.run_manifest(run).get("agent") == self.cfg.agent]
        if self.cfg.campaign_run:
            explicit = campaign.runs_dir(TARGET_DIR) / self.cfg.campaign_run
            if explicit not in lives:
                raise SupervisorError(f"--campaign-run {self.cfg.campaign_run} is not a live run")
            mine = [explicit]
        if claim.exists():
            holder = read_json(claim)
            if not isinstance(holder, dict) or holder.get("agent") != self.cfg.agent:
                who = holder if not isinstance(holder, dict) else \
                    f"agent={holder.get('agent')} pid={holder.get('pid')}"
                print(f"SUPERVISOR YIELD reason=claim_held {who}", flush=True)
                return False
            try:
                holder_pid = int(holder.get("pid") or 0)
            except (TypeError, ValueError):
                holder_pid = 0
            if holder_pid and holder_pid != os.getpid() and supervisor_alive(holder_pid):
                print(f"SUPERVISOR YIELD reason=healthy_supervisor pid={holder_pid}", flush=True)
                return False
            if len(mine) != 1:
                raise SupervisorError(
                    f"claim held by agent {self.cfg.agent} but {len(mine)} live run dirs carry "
                    f"that label; pass --campaign-run")
            self.run_dir = mine[0]
            holder["pid"] = os.getpid()
            holder["heartbeat_utc"] = utc_now()
            campaign.atomic_write_json(claim, holder)
            adopted = True
        else:
            if lives:
                names = ", ".join(run.name for run in lives)
                raise SupervisorError(
                    f"live run dir(s) without a claim [{names}]; refusing to init a new campaign")
            prereg = Path(self.cfg.prereg)
            if not prereg.is_absolute():
                prereg = TARGET_DIR / prereg
            if not prereg.is_file():
                raise SupervisorError(f"prereg file absent: {prereg}")
            prereg_bytes = prereg.read_bytes()
            run_id, _ = self.campaign_call(
                campaign.core_init, ROOT, DOMAIN, TARGET, gate=self.cfg.gate,
                prereg=str(prereg), agent=self.cfg.agent)
            self.run_dir = campaign.campaign_dir(ROOT, DOMAIN, TARGET) / run_id
            snapshot = self.run_dir / "prereg_snapshot.txt"
            if not snapshot.exists():
                snapshot.write_bytes(prereg_bytes)
            adopted = False
        self.state_dir = self.run_dir / SUPERVISOR_DIRNAME
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.load_status()
        if not adopted:
            self.charge(len(prereg_bytes), "prereg_snapshot")
        self.status["campaign_run_id"] = self.run_dir.name
        self.save_status()
        self.event("CAMPAIGN_ATTACHED", run_id=self.run_dir.name,
                   mode="adopted" if adopted else "initialized", agent=self.cfg.agent)
        self.last_heartbeat = time.monotonic()
        return True

    def heartbeat_if_due(self, force: bool = False) -> None:
        if self.run_dir is None:
            return
        if not force and time.monotonic() - self.last_heartbeat < self.cfg.heartbeat_seconds:
            return
        self.last_heartbeat = time.monotonic()
        try:
            self.campaign_call(campaign.core_heartbeat, ROOT, DOMAIN, TARGET)
        except campaign.CampaignError as exc:
            self.halt("CLAIM_LOST", f"heartbeat refused: {exc}")
            return
        if free_bytes() < self.cfg.free_reserve_bytes:
            self.halt("DISK_RESERVE", f"free space {free_bytes()} below reserve "
                                      f"{self.cfg.free_reserve_bytes}")

    def update_target_state(self, next_action: str) -> None:
        """Refresh latest_campaign via campaign.py; with --state-updates full also import the
        supervisor note and next_action (the original next_action is kept in the note)."""
        if self.state_dir is None:
            return
        if self.quiet_depth:
            self.deferred_state_note = next_action
            return
        try:
            if self.cfg.state_updates == "full":
                state = read_json(TARGET_DIR / "state.json")
                if not isinstance(state, dict):
                    raise campaign.CampaignError("state.json is not an object")
                note = state.get("supervisor") if isinstance(state.get("supervisor"), dict) \
                    else {}
                if "original_next_action" not in note:
                    note["original_next_action"] = str(state.get("next_action") or "")[:1200]
                note.update({
                    "campaign_run_id": self.status.get("campaign_run_id"),
                    "phase": self.status.get("phase"),
                    "cycle": self.status.get("cycle"),
                    "last_outcome": (self.status.get("last_outcome") or {}).get("outcome"),
                    "halt_reason": self.status.get("halt_reason"),
                    "updated_utc": utc_now(),
                })
                state["supervisor"] = note
                state["next_action"] = next_action[:600]
                state["updated_by"] = self.cfg.agent
                state["updated_utc"] = utc_now()
                if state.get("latest_campaign_verdict") not in campaign.TERMINAL:
                    state["latest_campaign_verdict"] = None
                import_path = self.state_dir / "state_import.json"
                self.write_json(import_path, {"targets": {IDENT: state}}, "state_import")
                self.campaign_call(campaign.core_state_import, ROOT, str(import_path))
            self.campaign_call(campaign.core_state_refresh, ROOT,
                               campaign.registry_domains(ROOT), IDENT, self.cfg.agent)
        except (campaign.CampaignError, OSError, ValueError) as exc:
            self.event("STATE_UPDATE_FAILED", error=str(exc)[:300])

    def flush_deferred_state(self) -> None:
        if self.deferred_state_note is not None and not self.quiet_depth:
            note, self.deferred_state_note = self.deferred_state_note, None
            self.update_target_state(note)

    # ----------------------------------------------------------- waiting
    def sleep(self, seconds: float) -> None:
        deadline = time.monotonic() + seconds
        while not self.stop_requested and time.monotonic() < deadline:
            time.sleep(min(5.0, max(0.0, deadline - time.monotonic())))
            self.heartbeat_if_due()
            if self.status.get("phase") == "HALTED":
                return

    def request_stop(self, signum, frame) -> None:  # noqa: ARG002
        self.stop_requested = True
        print(f"SUPERVISOR STOP_REQUESTED signal={signum}", flush=True)

    # ------------------------------------------------------- tree policy
    def snapshot_tree(self) -> dict[str, list]:
        own = str(self.state_dir) if self.state_dir else None
        out: dict[str, list] = {}
        for dirpath, dirnames, filenames in os.walk(TARGET_DIR):
            dirnames[:] = sorted(d for d in dirnames if d not in SNAPSHOT_SKIP_DIRS
                                 and not (own and os.path.join(dirpath, d) == own))
            for name in filenames:
                if name in SNAPSHOT_SKIP_FILES:
                    continue
                full = os.path.join(dirpath, name)
                if full == str(LOCK_PATH):
                    continue
                try:
                    st = os.stat(full)
                except OSError:
                    continue
                if not os.path.isfile(full):
                    continue
                rel = os.path.relpath(full, TARGET_DIR)
                out[rel] = [st.st_size, st.st_mtime_ns, sha256_file(Path(full))]
        return out

    def is_editable(self, rel: str) -> bool:
        if rel in LEDGER_APPEND_ONLY or rel in EDITABLE_CHECKPOINTS or rel in self.editable_extra:
            return True
        if rel.startswith(NEW_FILE_ROOTS):
            return rel in set(self.status.get("agent_created_files", []))
        return False

    def deny_rules(self, snapshot: dict) -> list[str]:
        base = "//" + str(TARGET_DIR).lstrip("/")
        rules = [f"{tool}({base}/{d}/**)" for tool in ("Edit", "Write") for d in FROZEN_DIRS]
        for rel in sorted(snapshot):
            top = rel.split("/", 1)[0]
            if top in FROZEN_DIRS or self.is_editable(rel):
                continue
            rules.extend(f"{tool}({base}/{rel})" for tool in ("Edit", "Write"))
        return rules

    @staticmethod
    def explicit_deny_rules() -> list[str]:
        base = "//" + str(TARGET_DIR).lstrip("/")
        return [f"{tool}({base}/{d}/**)" for tool in ("Edit", "Write") for d in FROZEN_DIRS] + \
               [f"{tool}({base}/{rel})" for tool in ("Edit", "Write") for rel in FROZEN_EXPLICIT]

    def settings_file(self, snapshot: dict, kind: str) -> Path:
        assert self.state_dir is not None
        rules = self.deny_rules(snapshot)
        if kind == "review":
            base = "//" + str(TARGET_DIR).lstrip("/")
            rules += [f"{tool}({base}/experiments/**)" for tool in ("Edit", "Write")]
        # verbose:false keeps --output-format json to a single result object. With
        # the config default the CLI streams every tool result into stdout, and a
        # real planning turn (which reads the large checkpoints) then blows the
        # response bound and dies on SIGXFSZ before returning its plan.
        payload = {"permissions": {"deny": rules}, "verbose": False}
        text = json.dumps(payload, indent=1) + "\n"
        path = self.state_dir / "agent_settings.json"
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            self.write_text(path, text, "agent_settings")
        return path

    @staticmethod
    def ledger_bytes() -> dict[str, bytes]:
        out = {}
        for rel in LEDGER_APPEND_ONLY:
            path = TARGET_DIR / rel
            out[rel] = path.read_bytes() if path.is_file() else b""
        return out

    def check_turn_changes(self, before: dict, ledgers: dict[str, bytes], after: dict) -> dict:
        violations, created, changed = [], [], []
        for rel in sorted(set(after) - set(before)):
            created.append(rel)
            if rel.endswith(".md"):
                violations.append(f"new Markdown file {rel}")
            elif not rel.startswith(NEW_FILE_ROOTS):
                violations.append(f"new file outside experiments/tests: {rel}")
        for rel in sorted(set(before) - set(after)):
            violations.append(f"deleted file {rel}")
        for rel in sorted(set(before) & set(after)):
            if before[rel] != after[rel]:
                changed.append(rel)
                if not self.is_editable(rel):
                    violations.append(f"protected path modified: {rel}")
        for rel, old in ledgers.items():
            path = TARGET_DIR / rel
            new = path.read_bytes() if path.is_file() else b""
            if new[:len(old)] != old:
                violations.append(f"ledger not append-only: {rel}")
        return {"ok": not violations, "violations": violations, "created": created,
                "changed": changed}

    @staticmethod
    def check_job_changes(before: dict, after: dict, artifact_rel: str | None,
                          allow_results: bool) -> dict:
        violations, created = [], []
        for rel in sorted(set(after) - set(before)):
            created.append(rel)
            if allow_results and rel == artifact_rel:
                continue
            violations.append(f"job created an unexpected file: {rel}")
        for rel in sorted(set(before) - set(after)):
            violations.append(f"job deleted file {rel}")
        for rel in sorted(set(before) & set(after)):
            if before[rel] != after[rel]:
                violations.append(f"job modified pre-existing file {rel}")
        return {"ok": not violations, "violations": violations, "created": created}

    # ------------------------------------------------------ research data
    @staticmethod
    def hypothesis_ids() -> tuple[set[str], int]:
        ids: set[str] = set()
        highest = 0
        path = TARGET_DIR / "notes" / "hypotheses.csv"
        if not path.is_file():
            return ids, highest
        with path.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.reader(fh):
                if row and HYP_RE.match(row[0]):
                    ids.add(row[0])
                    highest = max(highest, hypothesis_number(row[0]))
        return ids, highest


    # ------------------------------------------------------- validation
    @staticmethod
    def resolve_rel(value: str, roots: tuple[str, ...], suffix: str, label: str,
                    errors: list[str]) -> str | None:
        pure = PurePosixPath(value)
        if pure.is_absolute() or ".." in pure.parts or pure.as_posix() != value \
                or "\\" in value:
            errors.append(f"{label}: must be a normalized target-relative path: {value!r}")
            return None
        if not (TARGET_DIR / value).resolve().is_relative_to(TARGET_DIR.resolve()):
            errors.append(f"{label}: path resolves outside the target directory")
            return None
        if roots and not value.startswith(roots):
            errors.append(f"{label}: must be under {'/'.join(roots)}: {value}")
            return None
        if suffix and not value.endswith(suffix):
            errors.append(f"{label}: must end with {suffix}: {value}")
            return None
        return value

    @staticmethod
    def check_arguments(args: list[str], label: str, errors: list[str]) -> None:
        for arg in args:
            if "\n" in arg or "\0" in arg:
                errors.append(f"{label}: argument contains a control character")
            if arg.startswith("/"):
                try:
                    Path(arg).resolve().relative_to(TARGET_DIR)
                except ValueError:
                    errors.append(f"{label}: absolute argument outside the target dir: {arg}")
            elif ".." in PurePosixPath(arg).parts:
                errors.append(f"{label}: argument climbs out of the target dir: {arg}")

    def validate_verifier(self, verifier, producer_script: str | None,
                          errors: list[str]) -> None:
        sub = self.review_schema["properties"]["verifier"]
        errors += schema_errors(verifier, sub, "$.verifier")
        if not isinstance(verifier, dict):
            errors.append("verifier object required")
            return
        rel = self.resolve_rel(str(verifier.get("script", "")), ("tests/",),
                               ".py", "verifier.script", errors)
        if rel is None:
            return
        path = TARGET_DIR / rel
        if not path.is_file():
            errors.append(f"verifier.script does not exist: {rel}")
            return
        if not str(verifier.get("pass_marker") or "").strip():
            errors.append("verifier.pass_marker required")
        self.check_arguments(list(verifier.get("arguments", [])), "verifier", errors)
        if producer_script:
            if rel == producer_script:
                errors.append("verifier.script must differ from producer.script")
            stem = PurePosixPath(producer_script).stem
            text = path.read_text(encoding="utf-8", errors="replace")
            if re.search(rf"^\s*(from|import)\s+(\S+\.)?{re.escape(stem)}\b", text, re.M) \
                    or re.search(rf"import_module\(\s*['\"]([\w.]+\.)?{re.escape(stem)}['\"]", text):
                errors.append(f"verifier imports the producer module {stem}; not independent")

    def fingerprint(self, plan: dict) -> str:
        producer = plan.get("producer")
        if producer is not None:
            # Renaming an output or adding an unrelated verifier does not change
            # the mathematical trial. Input arguments and source bytes do.
            arguments = []
            for argument in producer["arguments"]:
                if argument == plan["artifact"]:
                    arguments.append("<artifact-output>")
                elif "=" in argument and argument.split("=", 1)[1] == plan["artifact"]:
                    arguments.append(argument.split("=", 1)[0] + "=<artifact-output>")
                else:
                    arguments.append(argument)
            material = {"producer_sha256": sha256_file(TARGET_DIR / producer["script"]),
                        "arguments": arguments}
        else:
            verifier = plan.get("verifier")
            material = {
                "artifact_sha256": sha256_file(TARGET_DIR / plan["artifact"]),
                "verifier_sha256": None if verifier is None else
                    sha256_file(TARGET_DIR / verifier["script"]),
                "arguments": [] if verifier is None else verifier["arguments"],
            }
        return sha256_bytes(canonical(material).encode())

    def history_conflict(self, fp: str) -> str | None:
        for record in self.status.get("history", []):
            if record.get("fingerprint") == fp:
                return (f"identical retry of executed task {record.get('task_id')} "
                        f"(outcome {record.get('outcome')}); change the actual trial or "
                        "advance the ranked queue")
        return None

    def validate_plan(self, plan, *, seed: bool, floor: int, ids: set[str]) -> list[str]:
        errors: list[str] = []
        if not isinstance(plan, dict):
            return ["plan is not an object"]
        verifier = plan.get("verifier")
        body = {k: v for k, v in plan.items() if k != "verifier"}
        errors += schema_errors(body, self.plan_schema)
        if errors:
            return errors
        if plan["decision"] == "stop":
            if not str(plan.get("stop_reason") or "").strip():
                errors.append("stop_reason required with decision stop")
            return errors
        if plan["front"] == "other" and (not plan.get("front_rank") or not plan.get("front_label")):
            errors.append("front 'other' requires front_rank and front_label")
        hyp = plan["hypothesis_id"]
        if seed:
            if hyp not in ids:
                self.event("SEED_HYPOTHESIS_ROW_MISSING", hypothesis_id=hyp)
        else:
            if hypothesis_number(hyp) <= floor:
                errors.append(f"hypothesis_id must be a new row above H{floor:03d}")
            if hyp not in ids:
                errors.append(f"{hyp} has no row in notes/hypotheses.csv (append it)")
        artifact = self.resolve_rel(plan["artifact"], ("results/",), ".json", "artifact", errors)
        producer = plan["producer"]
        producer_rel = None
        if producer is not None:
            producer_rel = self.resolve_rel(producer["script"], ("experiments/",), ".py",
                                            "producer.script", errors)
            if producer_rel and not (TARGET_DIR / producer_rel).is_file():
                errors.append(f"producer.script does not exist: {producer_rel}")
            self.check_arguments(list(producer["arguments"]), "producer", errors)
            if artifact and (TARGET_DIR / artifact).exists():
                errors.append(f"artifact already exists; a producer must create a new path: "
                              f"{artifact}")
        elif artifact and not (TARGET_DIR / artifact).is_file():
            errors.append(f"producer is null but artifact does not exist: {artifact}")
        if seed:
            if verifier is None:
                errors.append("seed plan requires a verifier object")
            else:
                self.validate_verifier(verifier, producer_rel, errors)
        elif verifier is not None:
            errors.append("verifier must come from the independent review turn, not the plan")
        for key in (plan.get("artifact_expect") or {}):
            if not re.fullmatch(r"[A-Za-z0-9_]+(\.[A-Za-z0-9_]+)*", key):
                errors.append(f"artifact_expect key is not a dotted path: {key!r}")
        if errors:
            return errors
        if producer is None:
            report = self.validate_artifact(artifact, 0, plan.get("artifact_expect") or {})
            errors.extend(report["errors"])
            if errors:
                return errors
        previous = self.status.get("last_outcome") or {}
        if previous.get("outcome") in {"PRODUCER_RESOURCE_BLOCKED", "VERIFIER_RESOURCE_BLOCKED"}:
            next_front = {"norm": "wlaw", "wlaw": "callen"}.get(previous.get("front"))
            if next_front and plan["front"] != next_front:
                errors.append(f"resource-blocked {previous['front']} must advance to {next_front}")
        conflict = self.history_conflict(self.fingerprint(plan))
        if conflict:
            errors.append(conflict)
        return errors

    def make_task(self, plan: dict, source: str) -> dict:
        number = int(self.status.get("cycle", 0)) + 1
        task_id = f"{number:04d}_{plan['hypothesis_id']}"
        task = {
            "task_id": task_id,
            "source": source,
            "created_utc": utc_now(),
            "front": plan["front"],
            "front_rank": plan.get("front_rank"),
            "front_label": plan.get("front_label"),
            "hypothesis_id": plan["hypothesis_id"],
            "next_prompt": plan["next_prompt"],
            "artifact": plan["artifact"],
            "producer": plan["producer"],
            "verifier": plan.get("verifier"),
            "verifier_requirements": plan["verifier_requirements"],
            "artifact_expect": plan.get("artifact_expect") or {},
            "checkpoint_evidence": plan["checkpoint_evidence"],
            "source_changes": plan.get("source_changes") or [],
            "fingerprint": self.fingerprint(plan),
            "stage": "PRODUCER" if plan.get("verifier") else "REVIEW",
            "attempts": {},
            "producer_result": None,
            "artifact_validation": None,
            "review": None,
            "verifier_result": None,
            "outcome": None,
        }
        self.write_json(self.state_dir / "tasks" / task_id / "plan.json", plan, "task_plan")
        return task

    # --------------------------------------------------------- artifacts
    def validate_artifact(self, rel: str, started_unix: float, expect: dict) -> dict:
        path = TARGET_DIR / rel
        report = {"ok": False, "errors": [], "path": rel, "size": None,
                  "sha256": None, "data_sha256": None, "expect": {}}
        if not path.is_file():
            report["errors"].append("artifact missing")
            return report
        st = path.stat()
        report["size"] = st.st_size
        if st.st_size > ARTIFACT_MAX_BYTES:
            report["errors"].append(f"artifact larger than {ARTIFACT_MAX_BYTES} bytes")
            return report
        if started_unix and st.st_mtime < started_unix - 2:
            report["errors"].append("artifact predates the producer run")
        report["sha256"] = sha256_file(path)
        try:
            doc = read_json(path)
        except (OSError, ValueError) as exc:
            report["errors"].append(f"artifact is not JSON: {exc}")
            return report
        if not isinstance(doc, dict) or not isinstance(doc.get("meta"), dict) \
                or "data" not in doc:
            report["errors"].append("artifact lacks meta/data objects")
            return report
        meta, data = doc["meta"], doc["data"]
        if type(meta.get("schema_version")) is not int or meta["schema_version"] != 1:
            report["errors"].append(f"meta.schema_version is {meta.get('schema_version')!r}")
        sources = meta.get("source_sha256")
        if not isinstance(sources, dict) or not sources:
            report["errors"].append("meta.source_sha256 missing or empty")
        else:
            for src_rel, digest in sources.items():
                pure = PurePosixPath(src_rel)
                if pure.is_absolute() or ".." in pure.parts or pure.as_posix() != src_rel \
                        or not (ROOT / src_rel).resolve().is_relative_to(DOMAIN_DIR.resolve()):
                    report["errors"].append(f"source path outside the research domain: {src_rel}")
                    continue
                src = ROOT / src_rel
                if not src.is_file():
                    report["errors"].append(f"meta.source_sha256 path missing: {src_rel}")
                elif sha256_file(src) != digest:
                    report["errors"].append(f"meta.source_sha256 mismatch: {src_rel}")
        report["data_sha256"] = sha256_bytes(canonical(data).encode())
        if meta.get("data_sha256") != report["data_sha256"]:
            report["errors"].append("meta.data_sha256 does not match canonical data")
        for key, expected in expect.items():
            try:
                actual = get_dotted(doc, key)
            except KeyError:
                report["expect"][key] = "MISSING"
                report["errors"].append(f"artifact_expect {key}: missing")
                continue
            ok = scalar_equal(actual, expected)
            report["expect"][key] = "OK" if ok else f"MISMATCH:{str(actual)[:80]}"
            if not ok:
                report["errors"].append(f"artifact_expect {key}: got {str(actual)[:80]!r}")
        report["ok"] = not report["errors"]
        return report

    # ------------------------------------------------------ guard jobs
    @staticmethod
    def job_command(script_rel: str, arguments: list[str], result_path: Path) -> list[str]:
        return [str(INTERPRETER), "-I", "-B", str(GUARD_PATH), "--cwd", str(TARGET_DIR),
                "--result", str(result_path), "--", str(INTERPRETER), "-I", "-B", script_rel,
                *arguments]

    def launch_job(self, role: str, task: dict) -> dict:
        spec = task["producer"] if role == "producer" else task["verifier"]
        attempts = task.setdefault("attempts", {})
        attempt = int(attempts.get(role, 0)) + 1
        attempts[role] = attempt
        task_dir = self.state_dir / "tasks" / task["task_id"]
        task_dir.mkdir(parents=True, exist_ok=True)
        result_path = task_dir / f"{role}.attempt{attempt}.guard.json"
        log_path = task_dir / f"{role}.log"
        pre_path = task_dir / f"{role}.pre.json"
        command = self.job_command(spec["script"], list(spec["arguments"]), result_path)
        pre_text = canonical(self.snapshot_tree()) + "\n"
        if not pre_path.exists() or pre_path.read_text(encoding="utf-8") != pre_text:
            self.write_text(pre_path, pre_text, "job_pre_snapshot")
        log_size_before = log_path.stat().st_size if log_path.exists() else 0
        env = {k: v for k, v in os.environ.items() if not k.startswith("CLAUDECODE")}
        env.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONUNBUFFERED": "1"})
        with log_path.open("ab") as log:
            self.child = subprocess.Popen(command, cwd=TARGET_DIR, stdin=subprocess.DEVNULL,
                                          stdout=log, stderr=subprocess.STDOUT, env=env,
                                          start_new_session=True)
        job = {
            "task_id": task["task_id"], "role": role, "attempt": attempt,
            "pid": self.child.pid, "started_utc": utc_now(), "started_unix": time.time(),
            "result_path": str(result_path), "log_path": str(log_path),
            "pre_path": str(pre_path), "log_size_before": log_size_before,
            "command": command, "ready": False,
        }
        self.status["active_job"] = job
        self.status["jobs"] = int(self.status.get("jobs", 0)) + 1
        self.status["phase"] = f"EXECUTING_{role.upper()}"
        self.save_status()
        self.event("JOB_LAUNCHED", task_id=task["task_id"], role=role, attempt=attempt,
                   pid=job["pid"], script=spec["script"])
        return job

    def adopt_job(self, job: dict) -> bool:
        """Re-attach to a guard launched by an earlier supervisor process; never relaunch.

        Returns True when the recorded pid is still that guard and must be waited for.
        """
        pid = int(job["pid"])
        if Path(job["result_path"]).exists() or not pid_alive(pid):
            return False
        if job["result_path"] not in pid_command(pid):
            self.event("JOB_PID_REUSED", pid=pid, task_id=job.get("task_id"))
            return False
        self.event("JOB_ADOPTED", pid=pid, task_id=job.get("task_id"), role=job.get("role"))
        return True

    def wait_job(self, job: dict) -> str:
        """Block until the guard exits ('EXITED') or a stop is requested ('DETACHED')."""
        pid = int(job["pid"])
        own = self.child is not None and self.child.pid == pid
        last_status = time.monotonic()
        kill_after = float(self.cfg.job_kill_after_seconds or 0)
        overrun_logged = False
        while True:
            if own:
                if self.child.poll() is not None:
                    return "EXITED"
            elif not pid_alive(pid):
                return "EXITED"
            if self.stop_requested:
                self.status["phase"] = "DETACHED"
                self.save_status()
                self.event("JOB_DETACHED", pid=pid, task_id=job.get("task_id"))
                return "DETACHED"
            self.heartbeat_if_due()
            if not job.get("ready"):
                log = Path(job["log_path"])
                if log.exists() and READY_MARK in tail_text(log, 64 * 1024):
                    job["ready"] = True
                    job["ready_utc"] = utc_now()
                    self.save_status()
            elapsed = time.time() - float(job["started_unix"])
            if elapsed > self.cfg.job_wall_ceiling_seconds and not overrun_logged:
                overrun_logged = True
                self.event("JOB_OVERRUN", pid=pid, elapsed_seconds=int(elapsed),
                           action="kill" if kill_after else "wait")
            if kill_after and elapsed > kill_after:
                self.event("JOB_KILLED_BY_POLICY", pid=pid, elapsed_seconds=int(elapsed))
                with contextlib.suppress(OSError):
                    os.kill(pid, signal.SIGTERM)
                kill_after = 0.0
            if time.monotonic() - last_status > 60:
                self.save_status()
                last_status = time.monotonic()
            time.sleep(10)

    @staticmethod
    def read_result(job: dict) -> dict | None:
        path = Path(job["result_path"])
        for _ in range(30):
            if path.exists():
                with contextlib.suppress(OSError, ValueError):
                    loaded = read_json(path)
                    if isinstance(loaded, dict):
                        return loaded
            time.sleep(1)
        return None

    @staticmethod
    def classify(result: dict | None, exit_code: int | None) -> str:
        if result is None:
            return "BUSY" if exit_code == 75 else "NO_RESULT"
        guard = result.get("guard_exit_code")
        child = result.get("child_exit_code")
        reason = str(result.get("reason") or "")
        if guard == 75 and child is None:
            if GUARD_BUSY_MARK in reason:
                return "BUSY"
            if "free disk" in reason:
                return "DISK"
            return "GUARD_SETUP_FAILED"
        if reason == "cancelled":
            return "INTERRUPTED"
        if guard == 0 and child == 0 and reason == "completed":
            return "ACCEPTED" if exit_code in (None, 0) else "GUARD_SETUP_FAILED"
        if GUARD_UNMEASURABLE_MARK in reason or GUARD_CLEANUP_MARK in reason:
            # The guard could not account for or clean up the group. That is an
            # operator-level anomaly, not an exhausted resource: keep it out of
            # the RESOURCE_BLOCKED bucket so it never retires a valid trial.
            return "UNACCOUNTABLE"
        if reason != "completed":
            return "RESOURCE_BLOCKED"
        return "CHILD_FAILED"

    def run_guarded(self, role: str, task: dict) -> dict:
        busy_since: float | None = None
        while True:
            if self.status.get("phase") == "HALTED":
                return {"kind": "HALTED"}
            job = self.status.get("active_job")
            if job is not None and job.get("task_id") != task["task_id"]:
                self.halt("ACTIVE_JOB_MISMATCH", f"active job belongs to {job.get('task_id')}")
                return {"kind": "HALTED"}
            self.quiet_depth += 1
            try:
                if job is None:
                    if free_bytes() < self.cfg.free_reserve_bytes:
                        self.halt("DISK_RESERVE", "free space below reserve before launch")
                        return {"kind": "HALTED"}
                    job = self.launch_job(role, task)
                    must_wait = True
                else:
                    must_wait = self.adopt_job(job)
                if must_wait and self.wait_job(job) == "DETACHED":
                    return {"kind": "DETACHED"}
                exit_code = self.child.returncode \
                    if self.child is not None and self.child.pid == job["pid"] else None
                self.child = None
                result = self.read_result(job)
                kind = self.classify(result, exit_code)
                log_path = Path(job["log_path"])
                log_size = log_path.stat().st_size if log_path.exists() else 0
                self.charge(max(0, log_size - int(job.get("log_size_before", 0))),
                            f"{role}_log")
                result_file = Path(job["result_path"])
                if result_file.exists():
                    self.charge(result_file.stat().st_size, "guard_result")
                self.save_status()
                self.event("JOB_EXITED", task_id=task["task_id"], role=role, kind=kind,
                           guard_exit_code=None if result is None else result.get("guard_exit_code"),
                           child_exit_code=None if result is None else result.get("child_exit_code"),
                           reason=None if result is None else str(result.get("reason"))[:200])
                integrity = {"ok": False, "violations": ["pre-snapshot unavailable"],
                             "created": []}
                pre_path = Path(job.get("pre_path") or "")
                if kind != "BUSY" and pre_path.is_file():
                    integrity = self.check_job_changes(
                        json.loads(pre_path.read_text(encoding="utf-8")), self.snapshot_tree(),
                        task["artifact"] if role == "producer" else None,
                        allow_results=(role == "producer"))
            finally:
                self.quiet_depth -= 1
                self.flush_deferred_state()
            if kind == "BUSY":
                self.status["active_job"] = None
                busy_since = busy_since or time.time()
                waited = time.time() - busy_since
                cap = float(self.cfg.guard_busy_max_seconds or 0)
                if cap and waited > cap:
                    self.halt("GUARD_BUSY_TOO_LONG", f"worker slot busy for {int(waited)} s")
                    return {"kind": "HALTED"}
                self.status["phase"] = f"WAITING_SLOT_{role.upper()}"
                self.save_status()
                self.sleep(self.cfg.guard_busy_retry_seconds)
                if self.stop_requested:
                    return {"kind": "DETACHED"}
                continue
            if kind == "DISK":
                self.halt("DISK_RESERVE", str(result.get("reason")))
                return {"kind": "HALTED"}
            if kind == "UNACCOUNTABLE":
                # Never finalize: an unaccounted or uncleaned group must not
                # write this trial into history and retire it forever.
                self.halt("GUARD_UNACCOUNTABLE", str(result.get("reason")))
                return {"kind": "HALTED"}
            return {"kind": kind, "result": result, "attempt": job["attempt"],
                    "started_unix": job["started_unix"], "log_tail": tail_text(log_path),
                    "log_path": str(log_path), "integrity": integrity,
                    "ready": bool(job.get("ready"))}

    # ------------------------------------------------------ agent turns
    def agent_command(self, schema: dict, settings_path: Path) -> list[str]:
        cmd = [str(self.cfg.claude), "-p", "--restricted", "--safe-mode", "--strict-mcp-config",
               "--no-chrome", "--prompt-suggestions", "false",
               "--tools", "Read,Glob,Grep,Edit,Write",
               "--disallowedTools", ",".join(self.explicit_deny_rules()),
               "--settings", str(settings_path),
               "--permission-mode", "acceptEdits", "--permission-prompts", "none",
               "--no-session-persistence", "--output-format", "json",
               "--json-schema", canonical(schema),
               "--append-system-prompt", SYSTEM_APPEND]
        if self.cfg.model:
            cmd += ["--model", self.cfg.model]
        if self.cfg.effort:
            cmd += ["--effort", self.cfg.effort]
        if self.cfg.max_budget_usd:
            cmd += ["--max-budget-usd", str(self.cfg.max_budget_usd)]
        return cmd

    def agent_backoff_active(self) -> bool:
        next_utc = self.status["agent"].get("next_attempt_utc")
        if not next_utc:
            return False
        return campaign.parse_z(next_utc) > campaign.now_utc()

    def record_agent_failure(self, error: str, fatal: bool = False) -> None:
        if fatal:
            self.halt("AGENT_CLI_INCOMPATIBLE", error[:600])
            return
        agent = self.status["agent"]
        agent["consecutive_failures"] = int(agent.get("consecutive_failures", 0)) + 1
        agent["last_error"] = error[:500]
        delay = min(self.cfg.agent_backoff_seconds * 2 ** (agent["consecutive_failures"] - 1),
                    self.cfg.agent_backoff_max_seconds)
        agent["next_attempt_utc"] = campaign.fmt_z(
            campaign.now_utc() + dt.timedelta(seconds=delay))
        self.status["phase"] = "BACKOFF"
        self.save_status()
        self.event("AGENT_FAILURE", failures=agent["consecutive_failures"],
                   retry_after_seconds=int(delay), error=error[:300])
        self.update_target_state(
            f"continuous_supervisor: planning unavailable ({agent['consecutive_failures']} "
            f"consecutive agent failures; next attempt {agent['next_attempt_utc']}); no "
            f"mathematical work performed by the supervisor itself. Last error: {error[:200]}")
        if agent["consecutive_failures"] >= self.cfg.agent_max_consecutive_failures:
            self.halt("AGENT_UNAVAILABLE", f"{agent['consecutive_failures']} consecutive "
                                           f"agent failures; last: {error[:300]}")

    def wait_for_quota(self, resets_unix: float) -> None:
        """Sleep out a rejected model quota window; not a failure, so no backoff."""
        # A rejected window that already looks expired still gets a floor: each
        # retry costs a real agent turn, so never poll the quota back-to-back.
        seconds = min(max(resets_unix + 30.0 - time.time(), QUOTA_WAIT_MIN_SECONDS),
                      QUOTA_WAIT_MAX_SECONDS)
        resets_utc = campaign.fmt_z(dt.datetime.fromtimestamp(resets_unix, dt.UTC))
        self.status["phase"] = "WAITING_QUOTA"
        self.status["quota_resets_utc"] = resets_utc
        self.save_status()
        self.event("AGENT_RATE_LIMITED", resets_utc=resets_utc, wait_seconds=int(seconds))
        self.update_target_state(
            f"continuous_supervisor: model quota window rejected; waiting until {resets_utc} "
            f"({int(seconds)} s) before the next planning turn. No mathematical work is "
            f"performed by the supervisor itself and no limit is relaxed.")
        self.sleep(seconds)
        self.status["quota_resets_utc"] = None
        self.save_status()

    def run_agent_turn(self, kind: str, prompt: str, schema: dict) -> dict:
        """Run a fresh restricted model turn inside the same exclusive resource guard."""
        assert self.state_dir is not None
        if free_bytes() < self.cfg.free_reserve_bytes:
            self.halt("DISK_RESERVE", "free space below reserve before agent turn")
            return {"ok": False, "error": "disk reserve", "fatal": False, "turn_dir": None}
        self.status["turns"] = int(self.status.get("turns", 0)) + 1
        turn_no = self.status["turns"]
        turn_dir = self.state_dir / "turns" / f"{turn_no:04d}_{kind}"
        prompt = excerpt(prompt, PROMPT_MAX_BYTES)
        self.write_text(turn_dir / "prompt.txt", prompt, "prompt")
        snapshot = self.snapshot_tree()
        self.write_json(turn_dir / "pre.json", {
            "tree": snapshot,
            "ledgers": {rel: {"size": len(data), "sha256": sha256_bytes(data)}
                        for rel, data in self.ledger_bytes().items()},
        }, "agent_pre_snapshot")
        settings_path = self.settings_file(snapshot, kind)
        agent_command = self.agent_command(schema, settings_path)
        stdout_path, stderr_path = turn_dir / "stdout.json", turn_dir / "stderr.txt"
        result_path = turn_dir / "guard.json"
        command = [str(INTERPRETER), "-I", "-B", str(GUARD_PATH),
                   "--cwd", str(TARGET_DIR), "--result", str(result_path), "--",
                   str(INTERPRETER), "-I", "-B", "-c", AGENT_REDIRECT_CODE,
                   str(RESPONSE_MAX_BYTES), str(stdout_path), str(stderr_path), *agent_command]
        self.write_json(turn_dir / "command.json", command, "agent_command")
        self.status["phase"] = f"PLANNING_{kind.upper()}"
        self.save_status()
        self.event("AGENT_TURN_STARTED", turn=turn_no, kind=kind, model=self.cfg.model)
        env = {k: v for k, v in os.environ.items() if not k.startswith("CLAUDECODE")}
        env["NO_COLOR"] = "1"
        started = time.time()
        deadline = started + self.cfg.agent_turn_timeout_seconds
        guard_log = turn_dir / "guard.log"
        try:
            with guard_log.open("xb") as log:
                proc = subprocess.Popen(command, cwd=TARGET_DIR, stdin=subprocess.PIPE,
                                        stdout=log, stderr=subprocess.STDOUT, env=env,
                                        start_new_session=True)
        except OSError as exc:
            return {"ok": False, "error": f"cannot start agent guard: {exc}", "fatal": True,
                    "turn_dir": turn_dir}
        self.status["active_agent"] = {
            "pid": proc.pid, "turn_dir": str(turn_dir), "kind": kind,
            "result_path": str(result_path), "deadline_unix": deadline,
        }
        self.save_status()
        pending_input = prompt.encode("utf-8")
        timed_out = False
        signalled = False
        while True:
            try:
                proc.communicate(input=pending_input, timeout=10)
                break
            except subprocess.TimeoutExpired:
                pending_input = None
                self.heartbeat_if_due()
                timed_out = timed_out or time.time() >= deadline
                if (timed_out or self.stop_requested or self.status["phase"] == "HALTED") \
                        and not signalled:
                    # Signal only our guard. It owns, stops, and reaps its whole group.
                    proc.send_signal(signal.SIGTERM)
                    signalled = True
        guard_result = self.read_result({"result_path": str(result_path)})
        self.status["active_agent"] = None
        self.save_status()
        stdout = stdout_path.read_bytes() if stdout_path.is_file() else b""
        stderr_text = tail_text(stderr_path, STDERR_MAX_BYTES)
        logical_bytes = sum(path.stat().st_size for path in
                            (stdout_path, stderr_path, result_path, guard_log) if path.is_file())
        self.charge(max(logical_bytes, (guard_result or {}).get("written_bytes", 0)),
                    "guarded_agent_writes")
        meta = {"exit_code": proc.returncode, "elapsed_seconds": round(time.time() - started, 1),
                "timed_out": timed_out, "stdout_bytes": len(stdout), "guard": guard_result}
        classification = self.classify(guard_result, proc.returncode)
        if classification != "ACCEPTED":
            if classification == "DISK":
                self.halt("DISK_RESERVE", str((guard_result or {}).get("reason")))
            if classification == "UNACCOUNTABLE":
                self.halt("GUARD_UNACCOUNTABLE", str((guard_result or {}).get("reason")))
            quota = quota_reset_unix(stdout)
            fatal = any(mark in stderr_text for mark in FATAL_CLI_MARKS)
            return {"ok": False, "busy": classification == "BUSY", "fatal": fatal,
                    "rate_limited_until": quota,
                    "error": f"agent {classification}: {(guard_result or {}).get('reason')}; "
                             f"{'rate limited; ' if quota else ''}{stderr_text[-300:]}",
                    "meta": meta, "turn_dir": turn_dir}
        if timed_out or len(stdout) >= RESPONSE_MAX_BYTES:
            return {"ok": False, "fatal": False, "error": "agent time or output bound reached",
                    "meta": meta, "turn_dir": turn_dir}
        try:
            envelope = json.loads(stdout)
            if isinstance(envelope, list):
                results = [event for event in envelope
                           if isinstance(event, dict) and event.get("type") == "result"]
                if len(results) != 1:
                    raise ValueError("expected exactly one final result event")
                envelope = results[0]
            if not isinstance(envelope, dict) or envelope.get("type") != "result":
                raise ValueError("expected a final result envelope")
            if envelope.get("is_error"):
                raise ValueError(f"agent error: {str(envelope.get('result'))[:300]}")
            output = envelope.get("structured_output")
            errors = schema_errors(output, schema)
            if errors:
                raise ValueError("; ".join(errors)[:500])
        except (ValueError, TypeError) as exc:
            return {"ok": False, "error": str(exc), "fatal": False, "meta": meta,
                    "turn_dir": turn_dir}
        meta.update({key: envelope.get(key) for key in
                     ("session_id", "subtype", "num_turns", "duration_ms", "total_cost_usd")})
        meta["permission_denials"] = len(envelope.get("permission_denials") or [])
        self.write_json(turn_dir / "output.json", output, "agent_output")
        return {"ok": True, "output": output, "fatal": False, "meta": meta, "turn_dir": turn_dir}

    def recover_agent_turn(self) -> None:
        """Wait for an already-owned turn; never overlap it or trust an interrupted reply."""
        active = self.status.get("active_agent")
        if not active:
            return
        pid = int(active["pid"])
        result_path = Path(active["result_path"])
        self.event("AGENT_RECOVERING", pid=pid, kind=active["kind"])
        signalled = False
        while not result_path.exists() and pid_alive(pid):
            if str(result_path) not in pid_command(pid):
                self.halt("AGENT_PID_REUSED", "recorded pid is not the owned guard")
                return
            if self.stop_requested:
                return
            if time.time() >= active["deadline_unix"] and not signalled:
                os.kill(pid, signal.SIGTERM)
                signalled = True
            self.heartbeat_if_due()
            time.sleep(10)
        turn_dir = Path(active["turn_dir"])
        try:
            before = read_json(turn_dir / "pre.json")
            integrity = self.check_turn_changes(before["tree"], {}, self.snapshot_tree())
        except (OSError, ValueError, KeyError, TypeError) as exc:
            # A missing or truncated pre-snapshot means the interrupted turn
            # cannot be cleared; halt with evidence instead of escaping run().
            self.status["active_agent"] = None
            self.save_status()
            self.halt("RECOVERY_SNAPSHOT_INVALID",
                      f"cannot audit interrupted {active['kind']} turn: {type(exc).__name__}: {exc}")
            return
        for rel, expected in before["ledgers"].items():
            path = TARGET_DIR / rel
            data = path.read_bytes() if path.is_file() else b""
            if sha256_bytes(data[:expected["size"]]) != expected["sha256"]:
                integrity["violations"].append(f"ledger not append-only: {rel}")
        self.status["active_agent"] = None
        self.write_json(turn_dir / "recovery.json", integrity, "agent_recovery")
        if integrity["violations"]:
            self.halt("AGENT_INTEGRITY", "; ".join(integrity["violations"])[:800])
            return
        files = set(self.status.get("agent_created_files", []))
        files.update(integrity["created"])
        self.status["agent_created_files"] = sorted(files)
        self.save_status()
        self.event("AGENT_INTERRUPTED_CHECKPOINTED", action="fresh bounded turn; no reply accepted")

    # ------------------------------------------------------------ prompts
    @staticmethod
    def context_blocks() -> dict[str, str]:
        def read(rel: str) -> str:
            path = TARGET_DIR / rel
            return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""
        next_actions = read("checkpoints/next_actions.md")
        current_state = read("checkpoints/current_state.md")
        failed_routes = read("checkpoints/failed_routes.md")
        hyp_rows = read("notes/hypotheses.csv").splitlines()
        return {
            "state_json": excerpt(read("state.json"), 3000),
            "ranked_queue": excerpt(md_section(next_actions, "ranked queue"), 9000),
            "prepared_prompt": excerpt(md_section(next_actions, "prepared next exact prompt"), 7000),
            "environment": excerpt(md_section(next_actions, "environment state"), 4000),
            "current_state_head": excerpt(current_state, 6000),
            "failed_routes_tail": excerpt("\n".join(failed_routes.splitlines()[-40:]), 3500),
            "hypotheses_header": hyp_rows[1] if len(hyp_rows) > 1 else "",
            "hypotheses_tail": excerpt("\n".join(hyp_rows[-4:]), 6000),
        }

    def history_block(self) -> str:
        rows = []
        for record in self.status.get("history", [])[-12:]:
            rows.append(f"- {record.get('task_id')} front={record.get('front')} "
                        f"outcome={record.get('outcome')} "
                        f"fingerprint={str(record.get('fingerprint', ''))[:16]} "
                        f"producer={record.get('producer')} verifier={record.get('verifier')} "
                        f"artifact={record.get('artifact')} at {record.get('utc')}")
        return "\n".join(rows) or "- none yet"

    def last_outcome_block(self) -> str:
        last = self.status.get("last_outcome")
        if not last:
            return "No task has completed under this supervisor yet."
        return excerpt(json.dumps(last, indent=1, sort_keys=True), 9000)

    @staticmethod
    def protected_block() -> str:
        return ("Protected (denied and hash-checked): " + ", ".join(f"{d}/" for d in FROZEN_DIRS)
                + ", every pre-existing file under experiments/ and tests/ (new files there are "
                  "allowed), state.json, README.md, pyproject.toml, uv.lock, all of notes/ except "
                  "notes/hypotheses.csv, and checkpoints/ except current_state.md, "
                  "next_actions.md, failed_routes.md. Append-only: "
                  "notes/hypotheses.csv, checkpoints/failed_routes.md, research_log.md.")

    def build_plan_prompt(self, floor: int, rejections: list[dict]) -> str:
        blocks = self.context_blocks()
        rejected = ""
        if rejections:
            rejected = "\n# Previous plan attempts in this cycle were rejected\n" + "\n".join(
                f"- attempt {r['attempt']}: " + "; ".join(r["errors"][:8]) for r in rejections)
        return f"""# Role
You are the planner and producer author for one bounded exact-research task (supervisor cycle {self.status.get('cycle', 0) + 1}, turn {self.status.get('turns', 0) + 1}, {utc_now()}), working in {TARGET_DIR}.

# Hard rules
1. You cannot run code. The supervisor runs `python -I -B <script> <arguments>` from this directory under tools/resource_guard.py: {LIMITS_TEXT}. A guard stop is a blocker even when the script prints PASS; a guard-refused worker slot is waited for by the supervisor, never retried by you.
2. The producer is a script under experiments/ that writes exactly one new JSON artifact under results/ at a path that does not exist yet. Artifact layout: meta {{schema_version: 1, source_sha256: {{workspace-relative path: sha256}}, data_sha256: sha256 of json.dumps(data, sort_keys=True, separators=(',',':'))}}, data {{...}}; the file has a trailing newline only. Live memory and timings go to stdout, never into data. Deterministic data only.
3. Surgical edits only. Never create Markdown files. {self.protected_block()} results/ is written only by guarded runs.
4. Append one new row to notes/hypotheses.csv (columns: {blocks['hypotheses_header']}) with a hypothesis id above H{floor:03d}; when the latest observed result is a failure or blocker, append a dated entry to checkpoints/failed_routes.md; update checkpoints/current_state.md and next_actions.md surgically if scope changed. Never rewrite or delete existing rows or entries.
5. Accurate scope: label executed measurements OBSERVED with the guard figures and unlaunched estimates PREFLIGHT. Never infer characteristic-zero reconstruction, aggregate height, real roots, W10 closure, or thermodynamic behavior without a separate exact proof. A finite-field degree is not a characteristic-zero bound; polynomial cancellation is not height-monotone; a canary is not a worst-case envelope.
6. No identical retries: the source bytes and mathematical input arguments identify an executed task, regardless of output names or unrelated source changes. A resource-blocked norm task must advance to W-law (rank 2); a resource-blocked W-law task must advance to Callen (rank 3). Otherwise justify a concrete source reduction or advance the ranked queue. Forbidden launches: blind full W10 closure, exponent fits from two sizes, dense CRT, full test suite, bulk exports.
7. Before any new producer runs, a separate fresh reviewer examines its safety and authors the independent verifier. Your verifier_requirements must describe independently checkable invariants and exact scope. Pin known falsifiable scalars in artifact_expect; never invent an unknown computed digest or degree. Include the producer and every actual dependency in meta.source_sha256, but do not name a verifier that does not yet exist.
8. Return exactly one JSON object matching the schema. Use decision "stop" with stop_reason only when no honest bounded task exists.

# Current target state (state.json)
{blocks['state_json']}

# Ranked queue (checkpoints/next_actions.md)
{blocks['ranked_queue']}

# Prepared exact prompt (checkpoints/next_actions.md)
{blocks['prepared_prompt']}

# Environment state (checkpoints/next_actions.md)
{blocks['environment']}

# Current state head (checkpoints/current_state.md)
{blocks['current_state_head']}

# Latest observed result under this supervisor (guard verdicts are authoritative)
{self.last_outcome_block()}

# Executed task fingerprints under this supervisor (may not be repeated)
{self.history_block()}

# Hypothesis ledger tail (notes/hypotheses.csv)
{blocks['hypotheses_tail']}

# Failed routes tail (checkpoints/failed_routes.md)
{blocks['failed_routes_tail']}
{rejected}

# Required output
One JSON object with: decision, front (norm|wlaw|callen|other with front_rank 1-11 and front_label), hypothesis_id (new row above H{floor:03d}), next_prompt (context-rich task statement with predictions and scope limits), artifact (new results/... .json path, or an existing artifact when producer is null), producer {{script under experiments/, arguments}} or null, verifier_requirements, artifact_expect (dotted path -> exact scalar), checkpoint_evidence (brief citations), source_changes (files you created or edited).
"""

    def build_review_prompt(self, task: dict, rejections: list[dict]) -> str:
        blocks = self.context_blocks()
        artifact_rel = task["artifact"]
        artifact_path = TARGET_DIR / artifact_rel
        summary = "Not produced yet: this is the pre-execution safety and verifier-authoring gate."
        if artifact_path.is_file():
            try:
                doc = read_json(artifact_path)
                meta = doc.get("meta", {}) if isinstance(doc, dict) else {}
                data = doc.get("data", {}) if isinstance(doc, dict) else {}
                scalars = {k: v for k, v in data.items()
                           if isinstance(v, (str, int, float, bool)) or v is None} \
                    if isinstance(data, dict) else {}
                summary = json.dumps({"meta": meta, "data_scalar_fields": scalars,
                                      "data_keys": sorted(data) if isinstance(data, dict) else None,
                                      "size_bytes": artifact_path.stat().st_size},
                                     indent=1, sort_keys=True)
            except (OSError, ValueError) as exc:
                summary = f"artifact unreadable: {exc}"
        producer = task.get("producer")
        producer_line = "none (verify an existing artifact)" if producer is None else \
            f"{producer['script']} {' '.join(producer['arguments'])}"
        producer_result = json.dumps(task.get("producer_result") or {}, sort_keys=True)
        validation = json.dumps(task.get("artifact_validation") or {}, sort_keys=True)
        existing = sorted(p.name for p in (TARGET_DIR / "tests").glob("test_trace_nine*.py"))
        rejected = ""
        if rejections:
            rejected = "\n# Previous review attempts for this task were rejected\n" + "\n".join(
                f"- attempt {r['attempt']}: " + "; ".join(r["errors"][:8]) for r in rejections)
        return f"""# Role
You are a fresh independent code reviewer and verifier author ({utc_now()}) in {TARGET_DIR}. You did not write the producer. New producers have NOT RUN yet: inspect their code before execution and author an independent verifier against the frozen inputs. Do not invent observations.

# Hard rules
1. Inspect the complete producer path for bounded resource use, one thread, no subprocesses or network, no deletion or existing-file mutation, and writes only to its declared new artifact. Reject unsafe code; never edit anything under experiments/. Write the independent verifier under tests/ without importing the producer or copying its reduction code. Recompute claims through an independent route, validate artifact provenance and every pinned scalar, exit 0 only on full PASS, print the declared marker only then, and never write a file.
2. The supervisor runs it as `python -I -B <script> <arguments>` from this directory under tools/resource_guard.py: {LIMITS_TEXT}. Keep the verifier inside those limits; a guard stop is a failed verification.
3. Surgical edits only. Never create Markdown files. {self.protected_block()}
4. State accepted_scope honestly: exactly what a successful independent run would certify. Until that run passes the guard this is a proposed scope, not an observed result. Do not infer wider characteristic-zero, height, root, closure, or thermodynamic claims without separate exact proof.
5. Use decision "reject" with reject_reason when the artifact or plan cannot be verified independently or its claims exceed their evidence; then append a dated entry to checkpoints/failed_routes.md. Never edit existing notes/hypotheses.csv rows; you may append a row only to record a rejection.
6. Return exactly one JSON object matching the schema.

# Task under review
- task_id: {task['task_id']}  hypothesis: {task['hypothesis_id']}  front: {task['front']} {task.get('front_label') or ''}
- producer: {producer_line}
- artifact: {artifact_rel}
- planner's task statement: {excerpt(task['next_prompt'], 5000)}
- verifier requirements from the planner: {excerpt(task['verifier_requirements'], 3000)}
- artifact_expect (predictions; checked only after production): {json.dumps(task.get('artifact_expect') or {}, sort_keys=True)}
- checkpoint evidence cited: {json.dumps(task.get('checkpoint_evidence'))}
- source changes declared by the planner: {json.dumps(task.get('source_changes'))}

# Producer guard result (authoritative)
{producer_result}

# Supervisor artifact validation
{validation}

# Artifact summary
{excerpt(summary, 6000)}

# Existing trace-nine verifiers (conventions only; frozen ones must not be edited)
{', '.join(existing) or 'none'}

# Current target state (state.json)
{blocks['state_json']}

# Current state head (checkpoints/current_state.md)
{excerpt(blocks['current_state_head'], 3500)}
{rejected}

# Required output
One JSON object with: decision (verify|reject), verifier {{script under tests/, arguments, pass_marker}} or null when rejecting, independence_statement, checked_claims, accepted_scope, concerns.
"""

    # ------------------------------------------------------------ planning
    def turn_with_integrity(self, kind: str, prompt: str, schema: dict) -> dict:
        self.quiet_depth += 1
        try:
            before = self.snapshot_tree()
            ledgers = self.ledger_bytes()
            result = self.run_agent_turn(kind, prompt, schema)
            after = self.snapshot_tree()
            integrity = self.check_turn_changes(before, ledgers, after)
            if kind == "review":
                integrity["violations"] += [
                    f"reviewer modified producer area: {rel}"
                    for rel in integrity["created"] + integrity["changed"]
                    if rel.startswith("experiments/")
                ]
                integrity["ok"] = not integrity["violations"]
        finally:
            self.quiet_depth -= 1
            self.flush_deferred_state()
        created = [rel for rel in integrity["created"] if rel.startswith(NEW_FILE_ROOTS)]
        if created:
            files = set(self.status.get("agent_created_files", []))
            files.update(created)
            self.status["agent_created_files"] = sorted(files)
        if result.get("turn_dir"):
            self.write_json(Path(result["turn_dir"]) / "integrity.json",
                            {"integrity": integrity, "meta": result.get("meta")}, "turn_integrity")
        result["integrity"] = integrity
        if not integrity["ok"]:
            self.halt("AGENT_INTEGRITY", "; ".join(integrity["violations"])[:800])
            result["ok"] = False
        self.event("AGENT_TURN_FINISHED", kind=kind, ok=result["ok"],
                   integrity_ok=integrity["ok"], violations=integrity["violations"][:6],
                   created=created[:8], meta=result.get("meta"))
        return result

    def safe_stop(self, reason: str, detail: str = "") -> None:
        self.status["phase"] = "STOPPED"
        self.status["halt_reason"] = reason
        self.status["halt_detail"] = detail[:2000] or None
        self.save_status()
        self.event("SAFE_STOP", reason=reason, detail=detail[:300])
        self.stop_requested = True

    def plan_next_task(self) -> dict | None:
        seed = self.status.get("seed")
        if seed and not seed.get("consumed"):
            try:
                plan = read_json(Path(seed["path"]))
            except (OSError, ValueError) as exc:
                self.halt("SEED_INVALID", f"cannot read seed plan: {exc}")
                return None
            ids, _ = self.hypothesis_ids()
            errors = self.validate_plan(plan, seed=True, floor=0, ids=ids)
            if errors:
                self.halt("SEED_INVALID", "; ".join(errors)[:1500])
                return None
            seed["consumed"] = True
            seed["consumed_utc"] = utc_now()
            self.save_status()
            self.event("SEED_ACCEPTED", hypothesis_id=plan["hypothesis_id"],
                       artifact=plan["artifact"])
            return self.make_task(plan, source="seed")
        if self.agent_backoff_active():
            self.status["phase"] = "BACKOFF"
            self.save_status()
            self.sleep(min(60.0, self.cfg.agent_backoff_seconds))
            return None
        if self.cfg.max_turns and self.status.get("turns", 0) >= self.cfg.max_turns:
            self.safe_stop("MAX_TURNS")
            return None
        ids, floor = self.hypothesis_ids()
        rejections: list[dict] = []
        for attempt in range(1, self.cfg.max_plan_rejections + 2):
            prompt = self.build_plan_prompt(floor, rejections)
            result = self.turn_with_integrity("plan", prompt, self.plan_schema)
            if not result["ok"]:
                if self.status.get("phase") == "HALTED":
                    return None
                if result.get("busy"):
                    self.sleep(self.cfg.guard_busy_retry_seconds)
                    return None
                if result.get("rate_limited_until"):
                    self.wait_for_quota(result["rate_limited_until"])
                    return None
                self.record_agent_failure(result.get("error", "unknown agent failure"),
                                          fatal=bool(result.get("fatal")))
                return None
            self.status["agent"].update({"consecutive_failures": 0, "next_attempt_utc": None})
            plan = result["output"]
            ids, _ = self.hypothesis_ids()
            errors = list(result["integrity"]["violations"]) + \
                self.validate_plan(plan, seed=False, floor=floor, ids=ids)
            if errors:
                rejections.append({"attempt": attempt, "errors": errors})
                self.event("PLAN_REJECTED", attempt=attempt, errors=errors[:8])
                self.save_status()
                continue
            if plan["decision"] == "stop":
                self.safe_stop("AGENT_STOP", str(plan.get("stop_reason")))
                self.update_target_state("continuous_supervisor stopped on planner decision: "
                                         + str(plan.get("stop_reason"))[:500])
                return None
            task = self.make_task(plan, source=f"turn-{self.status['turns']:04d}")
            self.event("TASK_PLANNED", task_id=task["task_id"], front=task["front"],
                       hypothesis_id=task["hypothesis_id"],
                       producer=(plan.get("producer") or {}).get("script"),
                       artifact=task["artifact"])
            return task
        self.halt("PLANNING_STALLED", "every plan attempt in this cycle was rejected: "
                  + "; ".join(rejections[-1]["errors"][:5]))
        return None

    def review_task(self, task: dict) -> str:
        """Fresh reviewer turn; returns 'VERIFY', 'REJECTED', or 'UNAVAILABLE'."""
        if self.agent_backoff_active():
            self.status["phase"] = "BACKOFF"
            self.save_status()
            self.sleep(min(60.0, self.cfg.agent_backoff_seconds))
            return "UNAVAILABLE"
        if self.cfg.max_turns and self.status.get("turns", 0) >= self.cfg.max_turns:
            self.safe_stop("MAX_TURNS")
            return "UNAVAILABLE"
        rejections: list[dict] = []
        producer_rel = (task.get("producer") or {}).get("script")
        for attempt in range(1, self.cfg.max_plan_rejections + 2):
            prompt = self.build_review_prompt(task, rejections)
            result = self.turn_with_integrity("review", prompt, self.review_schema)
            if not result["ok"]:
                if self.status.get("phase") == "HALTED":
                    return "UNAVAILABLE"
                if result.get("busy"):
                    self.sleep(self.cfg.guard_busy_retry_seconds)
                    return "UNAVAILABLE"
                if result.get("rate_limited_until"):
                    self.wait_for_quota(result["rate_limited_until"])
                    return "UNAVAILABLE"
                self.record_agent_failure(result.get("error", "unknown agent failure"),
                                          fatal=bool(result.get("fatal")))
                return "UNAVAILABLE"
            self.status["agent"].update({"consecutive_failures": 0, "next_attempt_utc": None})
            review = result["output"]
            errors = list(result["integrity"]["violations"]) + \
                schema_errors(review, self.review_schema)
            if not errors:
                if review["decision"] == "reject":
                    if not str(review.get("reject_reason") or "").strip():
                        errors.append("reject_reason required with decision reject")
                elif review.get("verifier") is None:
                    errors.append("verifier object required with decision verify")
                else:
                    self.validate_verifier(review["verifier"], producer_rel, errors)
            if errors:
                rejections.append({"attempt": attempt, "errors": errors})
                self.event("REVIEW_REJECTED_BY_SUPERVISOR", attempt=attempt, errors=errors[:8])
                self.save_status()
                continue
            task["review"] = {k: review.get(k) for k in ("decision", "reject_reason",
                                                          "independence_statement",
                                                          "checked_claims", "accepted_scope",
                                                          "concerns")}
            task["review"]["turn"] = self.status["turns"]
            if review["decision"] == "reject":
                return "REJECTED"
            task["verifier"] = review["verifier"]
            self.save_status()
            self.event("REVIEW_ACCEPTED", task_id=task["task_id"],
                       verifier=review["verifier"]["script"])
            return "VERIFY"
        self.halt("REVIEW_STALLED", "every review attempt for this task was rejected: "
                  + "; ".join(rejections[-1]["errors"][:5]))
        return "UNAVAILABLE"

    # ---------------------------------------------------------- execution
    def execute_task(self, task: dict) -> None:
        while not self.stop_requested and self.status.get("phase") != "HALTED":
            stage = task["stage"]
            if stage == "CHECKPOINT":
                self.finalize_task(task, task["outcome"])
                return
            if stage == "PRODUCER":
                if task["producer"] is None:
                    report = self.validate_artifact(task["artifact"], 0, task["artifact_expect"])
                    task["artifact_validation"] = report
                    if not report["ok"]:
                        self.finalize_task(task, "ARTIFACT_INVALID")
                        return
                    task["stage"] = "VERIFIER"
                    self.save_status()
                    continue
                run = self.run_guarded("producer", task)
                if run["kind"] in ("DETACHED", "HALTED"):
                    return
                task["producer_result"] = {
                    "guard": run["result"], "kind": run["kind"], "attempt": run["attempt"],
                    "log_tail": run["log_tail"], "log_path": run["log_path"],
                    "ready_seen": run["ready"], "integrity": run["integrity"],
                }
                self.status["active_job"] = None
                if not run["integrity"]["ok"]:
                    self.finalize_task(task, "JOB_INTEGRITY")
                    self.halt("JOB_INTEGRITY", "; ".join(run["integrity"]["violations"])[:800])
                    return
                if run["kind"] == "ACCEPTED":
                    report = self.validate_artifact(task["artifact"], run["started_unix"],
                                                    task.get("artifact_expect") or {})
                    task["artifact_validation"] = report
                    if not report["ok"]:
                        self.finalize_task(task, "ARTIFACT_INVALID")
                        return
                    task["stage"] = "VERIFIER"
                    self.save_status()
                    continue
                self.finalize_task(task, {
                    "RESOURCE_BLOCKED": "PRODUCER_RESOURCE_BLOCKED",
                    "CHILD_FAILED": "PRODUCER_FAILED",
                    "GUARD_SETUP_FAILED": "GUARD_SETUP_FAILED",
                }.get(run["kind"], "INTERRUPTED"))
                return
            if stage == "REVIEW":
                verdict = self.review_task(task)
                if verdict == "UNAVAILABLE":
                    return
                if verdict == "REJECTED":
                    self.finalize_task(task, "REVIEW_REJECTED")
                    return
                task["stage"] = "PRODUCER"
                self.save_status()
                continue
            if stage == "VERIFIER":
                artifact_path = TARGET_DIR / task["artifact"]
                before_sha = sha256_file(artifact_path) if artifact_path.is_file() else None
                run = self.run_guarded("verifier", task)
                if run["kind"] in ("DETACHED", "HALTED"):
                    return
                after_sha = sha256_file(artifact_path) if artifact_path.is_file() else None
                marker = (task["verifier"].get("pass_marker") or "").strip()
                log_file = Path(run["log_path"])
                log_text = log_file.read_text(encoding="utf-8", errors="replace") \
                    if log_file.exists() else ""
                marker_pattern = re.compile(re.escape(marker) + r"(?: [A-Za-z_][A-Za-z0-9_]*=\S+)*")
                marker_seen = bool(marker) and any(
                    marker_pattern.fullmatch(line.strip()) for line in log_text.splitlines())
                task["verifier_result"] = {
                    "guard": run["result"], "kind": run["kind"], "attempt": run["attempt"],
                    "log_tail": run["log_tail"], "log_path": run["log_path"],
                    "pass_marker": marker, "pass_marker_seen": marker_seen,
                    "artifact_sha256_before": before_sha, "artifact_sha256_after": after_sha,
                    "artifact_unchanged": before_sha is not None and before_sha == after_sha,
                    "integrity": run["integrity"],
                }
                self.status["active_job"] = None
                if not run["integrity"]["ok"]:
                    self.finalize_task(task, "JOB_INTEGRITY")
                    self.halt("JOB_INTEGRITY", "; ".join(run["integrity"]["violations"])[:800])
                    return
                if run["kind"] == "ACCEPTED":
                    ok = marker_seen and task["verifier_result"]["artifact_unchanged"]
                    self.finalize_task(task, "CONFIRMED" if ok else "VERIFIER_FAILED")
                    return
                self.finalize_task(task, {
                    "RESOURCE_BLOCKED": "VERIFIER_RESOURCE_BLOCKED",
                    "CHILD_FAILED": "VERIFIER_FAILED",
                    "GUARD_SETUP_FAILED": "GUARD_SETUP_FAILED",
                }.get(run["kind"], "INTERRUPTED"))
                return
            self.halt("UNKNOWN_STAGE", f"task stage {stage!r}")
            return

    def checkpoint_task(self, task: dict, summary: dict) -> None:
        """Persist the outcome and the next prompt before another model turn can start."""
        outcome = summary["outcome"]
        previous_front = task["front"]
        next_front = {"norm": "wlaw", "wlaw": "callen"}.get(previous_front)
        if "RESOURCE_BLOCKED" in outcome and next_front:
            action = f"Advance to the ranked {next_front} front; do not repeat {previous_front}."
        else:
            action = "Choose the highest admissible bounded exact test in the ranked queue."
        continuation = (
            f"Continue from {task['hypothesis_id']} ({previous_front}): {outcome}. {action} "
            "Read the linked evaluation, current checkpoint, failed routes and hypothesis ledger. "
            "State falsifiable predictions and scope; obtain independent pre-execution code review, "
            "run one bounded producer and its independent verifier, then checkpoint and continue. "
            "Do not repeat an unchanged failed input, rerun completed canaries, or interpret an inner "
            "PASS as a successful resource exit. Preserve nice 19, one thread, 35% of one core, "
            "the monitored 2 GiB stop, 50 GiB free reserve, bounded writes, and no deletions."
        )
        evaluation = self.state_dir / "tasks" / task["task_id"] / "evaluation.json"
        evidence = str(evaluation.relative_to(ROOT))
        guards = canonical({"producer": summary["producer_guard"],
                            "verifier": summary["verifier_guard"]})
        note = (f"\n### Supervisor outcome {task['task_id']} — {summary['utc']}\n\n"
                f"- {task['hypothesis_id']}: **{outcome}**; bounded scope only.\n"
                f"- Evidence: `{evidence}`.\n"
                f"- Artifact: `{task['artifact']}`.\n"
                f"- Guard verdicts: `{guards}`.\n"
                f"- Next prompt: {continuation}\n")
        current = TARGET_DIR / "checkpoints" / "current_state.md"
        current_text = current.read_text()
        current_updated = insert_note(current_text, note)
        if current_updated != current_text:
            self.write_text(current, current_updated, "current_checkpoint")
        failed = TARGET_DIR / "checkpoints" / "failed_routes.md"
        if outcome != "CONFIRMED":
            failed_text = failed.read_text()
            if note not in failed_text:
                self.write_text(failed, failed_text + note, "failed_route")
        next_actions = TARGET_DIR / "checkpoints" / "next_actions.md"
        text = next_actions.read_text()
        old = md_section(text, "prepared next exact prompt")
        if not old:
            raise SupervisorError("next_actions.md lacks its prepared next exact prompt section")
        replacement = old.splitlines()[0] + "\n\n" + continuation + f"\n\nEvidence: `{evidence}`."
        replaced = text.replace(old, replacement, 1)
        if replaced != text:
            self.write_text(next_actions, replaced, "next_prompt_checkpoint")
        ids, highest = self.hypothesis_ids()
        observation_id = task["hypothesis_id"] if task["hypothesis_id"] not in ids else f"H{highest + 1:03d}"
        row = io.StringIO()
        csv.writer(row, lineterminator="\n").writerow([
            observation_id, summary["utc"][:10],
            f"[OBSERVED][SUPERVISOR] Bounded task {task['hypothesis_id']} ended {outcome}.",
            "Independent evidence and the outer resource exit control acceptance.",
            task["verifier_requirements"][:1000],
            "Sequential bounded producer and independent verifier.",
            guards if outcome == "CONFIRMED" else "No accepted new result.",
            "No wider mathematical claim." if outcome == "CONFIRMED" else guards,
            "confirmed" if outcome == "CONFIRMED" else
                ("resource_blocked" if "RESOURCE_BLOCKED" in outcome else "observed_failure"),
            f"{task['artifact']};{evidence}", continuation,
        ])
        ledger = TARGET_DIR / "notes" / "hypotheses.csv"
        ledger_text = ledger.read_text()
        recorded = any(len(row) == 11 and evidence in row[9].split(";")
                       for row in csv.reader(io.StringIO(ledger_text)))
        if not recorded:
            self.write_text(ledger, ledger_text + row.getvalue(), "hypothesis_observation")
        progress = DOMAIN_DIR / "PROGRESS.md"
        progress_text = progress.read_text()
        progress_updated = insert_note(progress_text, note)
        if progress_updated != progress_text:
            self.write_text(progress, progress_updated, "domain_progress")
        self.refresh_allocations()
        self.status["continuation_prompt"] = continuation

    def refresh_allocations(self) -> None:
        _, highest = self.hypothesis_ids()
        numbers = [int(match.group(1)) for path in (TARGET_DIR / "experiments").glob("e*.py")
                   if (match := re.match(r"e(\d+)_", path.name))]
        path = TARGET_DIR / "checkpoints" / "next_actions.md"
        old = path.read_text()
        new, ids = re.subn(r"(LEDGER next free ID: `)H\d+(`)",
                          lambda match: f"{match[1]}H{highest + 1:03d}{match[2]}", old)
        new, scripts = re.subn(r"(Next unused experiment number: `)e\d+(`)",
                              lambda match: f"{match[1]}e{max(numbers, default=0) + 1}{match[2]}", new)
        if ids != 1 or scripts != 1:
            raise SupervisorError("checkpoint allocation fields are missing or ambiguous")
        if new != old:
            self.write_text(path, new, "allocation_checkpoint")

    def finalize_task(self, task: dict, outcome: str) -> None:
        task["outcome"] = outcome
        task["stage"] = "CHECKPOINT"
        if not task.get("finished_utc"):
            task["finished_utc"] = utc_now()
        if "executed_source_sha256" not in task:
            task["executed_source_sha256"] = {
                spec["script"]: sha256_file(TARGET_DIR / spec["script"])
                for spec in (task.get("producer"), task.get("verifier"))
                if spec and (TARGET_DIR / spec["script"]).is_file()
            }
        self.status["active_job"] = None
        self.save_status()
        task_dir = self.state_dir / "tasks" / task["task_id"]
        self.write_json(task_dir / "evaluation.json", {**task, "stage": "DONE"}, "evaluation")
        producer = task.get("producer") or {}
        verifier = task.get("verifier") or {}
        guard_p = (task.get("producer_result") or {}).get("guard") or {}
        guard_v = (task.get("verifier_result") or {}).get("guard") or {}
        summary = {
            "task_id": task["task_id"], "source": task["source"], "outcome": outcome,
            "front": task["front"], "front_label": task.get("front_label"),
            "hypothesis_id": task["hypothesis_id"], "artifact": task["artifact"],
            "producer": producer.get("script"), "producer_arguments": producer.get("arguments"),
            "verifier": verifier.get("script"), "verifier_arguments": verifier.get("arguments"),
            "producer_guard": guard_p, "verifier_guard": guard_v,
            "artifact_validation": task.get("artifact_validation"),
            "verifier_pass_marker_seen": (task.get("verifier_result") or {}).get("pass_marker_seen"),
            "artifact_unchanged_by_verifier": (task.get("verifier_result") or {}).get(
                "artifact_unchanged"),
            "review": task.get("review"),
            "producer_log_tail": (task.get("producer_result") or {}).get("log_tail"),
            "verifier_log_tail": (task.get("verifier_result") or {}).get("log_tail"),
            "utc": task["finished_utc"],
            "scope_note": "CONFIRMED means the independent verifier completed under its guard "
                          "with exit 0, artifact provenance and pinned data passed, its exact "
                          "PASS marker appeared, and artifact bytes stayed unchanged. A producer "
                          "run is asserted only when this task launched one; existing artifacts "
                          "require separately recorded producer evidence. No wider claim follows.",
        }
        self.status["last_outcome"] = summary
        history = self.status.setdefault("history", [])
        if not any(record["task_id"] == task["task_id"] for record in history):
            history.append({"task_id": task["task_id"], "fingerprint": task["fingerprint"],
                            "hypothesis_id": task["hypothesis_id"], "front": task["front"],
                            "outcome": outcome, "utc": task["finished_utc"],
                            "producer": producer.get("script"), "verifier": verifier.get("script"),
                            "artifact": task["artifact"]})
        if outcome == "CONFIRMED":
            frozen = {producer.get("script"), verifier.get("script")}
            self.status["agent_created_files"] = [
                path for path in self.status.get("agent_created_files", []) if path not in frozen
            ]
        if outcome != "JOB_INTEGRITY":
            self.checkpoint_task(task, summary)
        task["stage"] = "DONE"
        self.status["current_task"] = None
        self.status["phase"] = "PLANNING"
        self.save_status()
        self.event("TASK_FINISHED", task_id=task["task_id"], outcome=outcome,
                   hypothesis_id=task["hypothesis_id"], front=task["front"],
                   producer_reason=guard_p.get("reason"), verifier_reason=guard_v.get("reason"),
                   peak_memory_mib=round(max(guard_p.get("peak_memory_bytes") or 0,
                                             guard_v.get("peak_memory_bytes") or 0) / MIB, 1))
        self.update_target_state(
            f"continuous_supervisor cycle {self.status.get('cycle')}: {task['hypothesis_id']} "
            f"({task['front']}) {outcome}; producer={producer.get('script') or 'none'} "
            f"verifier={verifier.get('script') or 'none'} artifact={task['artifact']}; "
            f"guard reasons producer={guard_p.get('reason')} verifier={guard_v.get('reason')}; "
            f"bounded observed scope only; next: fresh plan turn")

    # --------------------------------------------------------------- main
    def preflight_seed(self) -> list[str]:
        """Validate the seed file against a fresh history before touching the campaign."""
        if not self.cfg.seed_plan:
            return []
        path = Path(self.cfg.seed_plan).resolve()
        if not path.is_file():
            return [f"seed plan absent: {path}"]
        try:
            plan = read_json(path)
        except (OSError, ValueError) as exc:
            return [f"seed plan unreadable: {exc}"]
        ids, _ = self.hypothesis_ids()
        return self.validate_plan(plan, seed=True, floor=0, ids=ids)

    def attach_seed(self) -> None:
        if not self.cfg.seed_plan:
            return
        path = Path(self.cfg.seed_plan).resolve()
        digest = sha256_file(path)
        seed = self.status.get("seed")
        if seed and seed.get("sha256") == digest:
            return
        if seed and not seed.get("consumed"):
            self.event("SEED_REPLACED", previous=seed.get("sha256"))
        self.status["seed"] = {"path": str(path), "sha256": digest, "consumed": False}
        self.save_status()
        self.event("SEED_REGISTERED", path=str(path), sha256=digest)

    def run(self) -> int:
        if not self.lock():
            return 0
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(sig, self.request_stop)
        seed_errors = self.preflight_seed()
        if seed_errors:
            print("SUPERVISOR ERROR seed plan invalid: " + "; ".join(seed_errors), flush=True)
            return 3
        try:
            if not self.attach_campaign():
                return 0
        except (campaign.CampaignError, SupervisorError) as exc:
            print(f"SUPERVISOR ERROR {exc}", flush=True)
            return 3
        if self.status.get("phase") == "HALTED" and not self.cfg.resume_halted:
            print(f"SUPERVISOR YIELD reason=halted halt_reason={self.status.get('halt_reason')} "
                  f"detail={self.status.get('halt_detail')}", flush=True)
            return 0
        if self.status.get("phase") in ("HALTED", "STOPPED"):
            self.event("RESUMED", previous=self.status.get("halt_reason"))
            self.status["halt_reason"] = None
            self.status["halt_detail"] = None
        self.status["phase"] = "READY"
        self.attach_seed()
        self.save_status()
        try:
            self.recover_agent_turn()
        except Exception:  # noqa: BLE001  recovery must halt with evidence, never traceback out
            self.halt("RECOVERY_ERROR", traceback.format_exc()[-1800:])
        if self.status.get("phase") == "HALTED" or self.stop_requested:
            return self.exit_code
        print(f"SUPERVISOR READY run={self.run_dir.name} pid={os.getpid()} "
              f"state_dir={self.state_dir}", flush=True)
        startup_note = (
            f"continuous_supervisor running (campaign {self.run_dir.name}, pid {os.getpid()}); "
            f"phase READY; cycle {self.status.get('cycle')}")
        if self.status.get("active_job"):
            self.deferred_state_note = startup_note
        else:
            self.update_target_state(startup_note)
        cycles_started = 0
        resumed_counted = False
        while not self.stop_requested and self.status.get("phase") != "HALTED":
            try:
                if free_bytes() < self.cfg.free_reserve_bytes:
                    self.halt("DISK_RESERVE", f"free space {free_bytes()} below reserve")
                    break
                task = self.status.get("current_task")
                if task is None:
                    if self.cfg.max_cycles and cycles_started >= self.cfg.max_cycles:
                        self.safe_stop("MAX_CYCLES")
                        break
                    task = self.plan_next_task()
                    if task is None:
                        continue
                    cycles_started += 1
                    self.status["cycle"] = int(self.status.get("cycle", 0)) + 1
                    self.status["current_task"] = task
                    self.save_status()
                    if self.cfg.plan_only:
                        self.safe_stop("PLAN_ONLY",
                                       f"task {task['task_id']} planned, not launched")
                        break
                elif not resumed_counted:
                    resumed_counted = True
                    cycles_started += 1
                    self.event("TASK_RESUMED", task_id=task.get("task_id"),
                               stage=task.get("stage"))
                self.execute_task(task)
            except Exception:  # noqa: BLE001  a daemon must halt with evidence, not vanish
                self.quiet_depth = 0
                self.halt("INTERNAL_ERROR", traceback.format_exc()[-1800:])
                break
        if self.stop_requested and self.status.get("phase") not in ("HALTED", "STOPPED",
                                                                    "DETACHED"):
            self.safe_stop("SIGNAL")
        self.flush_deferred_state()
        self.heartbeat_if_due(force=True)
        return self.exit_code


# ------------------------------------------------------------------------ CLI

def locate_state_dir(cfg: argparse.Namespace) -> Path | None:
    if cfg.campaign_run:
        run = campaign.runs_dir(TARGET_DIR) / cfg.campaign_run
        return run / SUPERVISOR_DIRNAME if run.is_dir() else None
    mine = [run for run in campaign.live_run_dirs(TARGET_DIR)
            if campaign.run_manifest(run).get("agent") == cfg.agent]
    return mine[0] / SUPERVISOR_DIRNAME if len(mine) == 1 else None


def offline_supervisor(cfg: argparse.Namespace) -> Supervisor:
    sup = Supervisor(cfg)
    state_dir = locate_state_dir(cfg)
    if state_dir and (state_dir / STATUS_NAME).is_file():
        sup.state_dir = state_dir
        sup.load_status()
        sup.state_dir = None  # read-only helpers must never write evidence
    return sup


def cmd_status(cfg: argparse.Namespace) -> int:
    state_dir = locate_state_dir(cfg)
    if state_dir is None or not (state_dir / STATUS_NAME).is_file():
        print(json.dumps({"supervisor": None, "reason": "no live campaign with supervisor state"}))
        return 0
    status = read_json(state_dir / STATUS_NAME)
    status["active_job_alive"] = bool(status.get("active_job")) and \
        pid_alive(int(status["active_job"]["pid"]))
    holder = read_lock_holder()
    status["lock_holder"] = holder
    status["lock_holder_alive"] = bool(holder) and pid_alive(int(holder.get("pid") or 0))
    print(json.dumps(status, indent=1, sort_keys=True))
    return 0


def cmd_check(cfg: argparse.Namespace) -> int:
    report = {
        "interpreter": str(INTERPRETER), "interpreter_ok": INTERPRETER.is_file(),
        "guard": str(GUARD_PATH), "guard_ok": GUARD_PATH.is_file(),
        "claude": str(cfg.claude), "claude_ok": Path(cfg.claude).is_file(),
        "plan_schema_ok": PLAN_SCHEMA_PATH.is_file(),
        "review_schema_ok": REVIEW_SCHEMA_PATH.is_file(),
        "free_bytes": free_bytes(), "free_reserve_bytes": cfg.free_reserve_bytes,
        "free_ok": free_bytes() >= cfg.free_reserve_bytes,
        "target_dir": str(TARGET_DIR), "ident": IDENT,
    }
    try:
        report["registry_ok"] = DOMAIN in campaign.registry_domains(ROOT)
    except campaign.CampaignError as exc:
        report["registry_ok"] = False
        report["registry_error"] = str(exc)
    claim = campaign.claim_path(ROOT, DOMAIN, TARGET)
    report["claim"] = read_json(claim) if claim.exists() else None
    report["live_runs"] = [run.name for run in campaign.live_run_dirs(TARGET_DIR)]
    prereg = Path(cfg.prereg)
    if not prereg.is_absolute():
        prereg = TARGET_DIR / prereg
    report["prereg"] = str(prereg)
    report["prereg_ok"] = prereg.is_file()
    holder = read_lock_holder()
    report["lock_holder"] = holder
    report["lock_holder_alive"] = bool(holder) and pid_alive(int(holder.get("pid") or 0))
    if cfg.seed_plan:
        report["seed_errors"] = offline_supervisor(cfg).preflight_seed()
    report["ok"] = all(report.get(k) for k in ("interpreter_ok", "guard_ok", "claude_ok",
                                               "plan_schema_ok", "review_schema_ok", "free_ok",
                                               "registry_ok", "prereg_ok")) \
        and not report.get("seed_errors")
    print(json.dumps(report, indent=1, sort_keys=True))
    return 0 if report["ok"] else 3


def cmd_validate_plan(cfg: argparse.Namespace) -> int:
    """Offline validation; the cycle floor (id above the pre-turn maximum) is a live-loop check."""
    sup = offline_supervisor(cfg)
    plan = read_json(Path(cfg.file))
    ids, highest = sup.hypothesis_ids()
    errors = sup.validate_plan(plan, seed=cfg.seed, floor=0, ids=ids)
    print(json.dumps({"file": cfg.file, "seed": cfg.seed, "errors": errors,
                      "hypothesis_highest": highest,
                      "fingerprint": None if errors else sup.fingerprint(plan)}, indent=1))
    return 0 if not errors else 3


def cmd_render_prompt(cfg: argparse.Namespace) -> int:
    sup = offline_supervisor(cfg)
    _, floor = sup.hypothesis_ids()
    if cfg.kind == "plan":
        print(sup.build_plan_prompt(floor, []))
        return 0
    task = sup.status.get("current_task")
    if not task:
        print("no current task to review", file=sys.stderr)
        return 3
    print(sup.build_review_prompt(task, []))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="continuous_supervisor.py", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--agent", default="continuous_supervisor",
                        help="campaign agent label (claim/manifest owner)")
    common.add_argument("--gate", default="wave29-exact-continuation")
    common.add_argument("--prereg", default="checkpoints/next_actions.md",
                        help="existing preregistration file (hash-pinned by campaign init)")
    common.add_argument("--campaign-run", default=None, help="explicit live run id to adopt")
    common.add_argument("--claude", default=str(CLAUDE_DEFAULT))
    common.add_argument("--model", default="fable")
    common.add_argument("--effort", choices=("low", "medium", "high", "xhigh", "max"), default=None)
    common.add_argument("--max-budget-usd", type=float, default=None)
    common.add_argument("--seed-plan", default=None, help="seed plan JSON (plan schema + verifier)")
    common.add_argument("--editable", action="append", default=[],
                        help="extra target-relative file the agent may edit (repeatable)")
    common.add_argument("--free-reserve-bytes", type=int, default=FREE_RESERVE_DEFAULT)
    common.add_argument("--write-allowance-bytes", type=int, default=WRITE_ALLOWANCE_DEFAULT)
    common.add_argument("--heartbeat-seconds", type=float, default=300.0)
    common.add_argument("--state-updates", choices=("full", "refresh-only"), default="full",
                        help="full = campaign.py state import (next_action + supervisor note) "
                             "then refresh; refresh-only = never rewrite next_action")
    common.add_argument("--job-wall-ceiling-seconds", type=float, default=6 * 3600.0,
                        help="log JOB_OVERRUN past this; the job is never signalled unless "
                             "--job-kill-after-seconds is set")
    common.add_argument("--job-kill-after-seconds", type=float, default=0.0)
    common.add_argument("--guard-busy-retry-seconds", type=float, default=600.0)
    common.add_argument("--guard-busy-max-seconds", type=float, default=48 * 3600.0,
                        help="0 = wait for the worker slot indefinitely")
    common.add_argument("--agent-turn-timeout-seconds", type=float, default=2700.0)
    common.add_argument("--agent-backoff-seconds", type=float, default=60.0)
    common.add_argument("--agent-backoff-max-seconds", type=float, default=3600.0)
    common.add_argument("--agent-max-consecutive-failures", type=int, default=8)
    common.add_argument("--max-plan-rejections", type=int, default=2,
                        help="re-prompts per cycle after a rejected plan/review")
    sub = parser.add_subparsers(dest="command", required=True)
    p_run = sub.add_parser("run", parents=[common], help="persistent loop")
    p_run.add_argument("--max-cycles", type=int, default=0, help="0 = unbounded")
    p_run.add_argument("--max-turns", type=int, default=0, help="0 = unbounded agent turns")
    p_run.add_argument("--plan-only", action="store_true",
                       help="plan (or accept the seed) and stop before any job launch")
    p_run.add_argument("--resume-halted", action="store_true",
                       help="continue after a HALTED status instead of yielding")
    p_run.set_defaults(func=lambda cfg: Supervisor(cfg).run())
    p_status = sub.add_parser("status", parents=[common], help="print machine status JSON")
    p_status.set_defaults(func=cmd_status)
    p_check = sub.add_parser("check", parents=[common], help="preflight report, no launches")
    p_check.set_defaults(func=cmd_check)
    p_val = sub.add_parser("validate-plan", parents=[common], help="validate a plan file")
    p_val.add_argument("file")
    p_val.add_argument("--seed", action="store_true", help="validate as a seed plan")
    p_val.set_defaults(func=cmd_validate_plan)
    p_prompt = sub.add_parser("render-prompt", parents=[common],
                              help="print the prompt a turn would receive")
    p_prompt.add_argument("--kind", choices=("plan", "review"), default="plan")
    p_prompt.set_defaults(func=cmd_render_prompt)
    return parser


def main(argv: list[str] | None = None) -> int:
    cfg = build_parser().parse_args(argv)
    try:
        return int(cfg.func(cfg) or 0)
    except SupervisorError as exc:
        print(f"SUPERVISOR ERROR {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())
