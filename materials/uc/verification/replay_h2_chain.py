#!/usr/bin/env python3
"""Replay driver for the Liu Hypothesis 2 / unconditional-constant certificate chain.

The chain is eight machine-checked artifacts, each produced by one module:

    reduction        uc/liu9_h2_reduction.py          -> liu9-h2-reduction.json
    psi-reduction    uc/liu9_psi_reduction_audit.py   -> liu9-psi-reduction.json
    twovar           uc/liu9_h2_twovar_lemmas.py      -> liu9-h2-twovar.json
    boundary         uc/liu9_h2_boundary_layer.py     -> liu9-h2-boundary.json
    phi-audit        uc/liu9_h2_phi_audit.py          -> liu9-h2-phi-audit.json
    general-lift     uc/liu9_h2_general_lift.py       -> liu9-h2-general-lift.json
    mixture-theorem  uc/liu9_h2_mixture_theorem.py    -> liu9-h2-mixture-theorem.json
    four-fifths-ab   uc/liu9_cprime_four_fifths_ab.py -> liu9-cprime-four-fifths-ab.json

The first seven prove Liu's Hypothesis 2 and the unconditional constant
c' = 1 - m*.  The eighth proves A >= 0 and B >= 0 for the scaled Example-5
protocol f(x) = (4/5)x(1-x) and hence the better constant c'' = 1 - m_{16/25}.
Four of its six dependency pins are chain links and are cross-checked in [B]
below (twovar_module, boundary_module, mixture_module against those links'
module hashes; mixture_artifact against the mixture-theorem file hash).  The
remaining two, frontier_module (uc/liu9_cprime_frontier.py) and
frontier_artifact (uc/verification/results/liu9-cprime-frontier.json), are
outside the chain: the certificate module asserts them itself at run time and
refuses to certify on a mismatch, so [C] covers them indirectly.

For every link this driver

  A. checks that each pinned SHA-256 below (artifact file, internal report
     digest, module source) occurs verbatim in PROGRESS.md, so the ledger and
     this table cannot drift apart silently;
  B. recomputes the artifact file hash, the internal digest (from the artifact
     body, in the module's own canonicalisation) and the module hash from disk,
     and asserts them against the pins and against the hashes the artifact
     itself embeds (`report_sha256`, `tool_sha256`, dependency pins);
  C. re-runs the module the requested number of times (default 2) from the
     repository root with the invoking interpreter under `-I -B`, and requires
     every run to reproduce the on-disk artifact byte for byte.

Runs never overwrite a pinned artifact: modules with `--output` write to a
scratch directory, and the one module that writes a fixed path
(`liu9_psi_reduction_audit.py`) is snapshotted to the scratch directory before
each run and restored from that snapshot whenever the file on disk differs
afterwards, also when the run is interrupted.  On success the scratch
directory is deleted; on any failure it is kept and its path printed, so a
divergent output remains available as the witness.  Exit status is 0 only when
every check passes; any mismatch is listed and the status is nonzero.

Standard library only.  Invoke from the repository root as

    nice -n 19 ./.venv/bin/python -I -B uc/verification/replay_h2_chain.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

HERE = Path(__file__).resolve().parent          # uc/verification
ROOT = HERE.parent.parent                       # math/
LEDGER = ROOT / "PROGRESS.md"


@dataclass(frozen=True)
class Link:
    name: str
    module: str                 # relative to ROOT
    artifact: str               # relative to ROOT
    claim_status: str
    output_flag: bool           # module accepts --output; else it writes its fixed path
    digest_scope: str           # "omit" / "omit+nl" / "empty": how the module digests its body
    file_sha256: str
    report_sha256: str
    tool_sha256: str


# Pinned values.  Every one must also appear verbatim in PROGRESS.md (check A).
CHAIN: Tuple[Link, ...] = (
    Link("reduction", "uc/liu9_h2_reduction.py",
         "uc/verification/results/liu9-h2-reduction.json",
         "REDUCTION-CERTIFIED", True, "omit",
         "c63ae94e4fbd2397e19ea9e2a4e2c155f1e02265852d0aacad5c97d70f078e3a",
         "c40d27be152217bd0d862bf54c4a957ab7bcad22e6b2f7747cf7462e9bc7a397",
         "99840c239eb60c3db3a1093827f3c9849c5d010debec70eb13661a87e0418fe0"),
    Link("psi-reduction", "uc/liu9_psi_reduction_audit.py",
         "uc/verification/results/liu9-psi-reduction.json",
         "PROVED", False, "empty",
         "9fb5766f1c3110cce9b231d2f2f1b7c1eec67bc3d6ac7602be50fc72b00c07a3",
         "405eea478f2146fde762f749dd4a8f6986be79f0be47f56cb2cc24daf64b2af0",
         "18602dc7c8599ab6752583535ed5809df571deb8d8851e4accb4d07362e927bb"),
    Link("twovar", "uc/liu9_h2_twovar_lemmas.py",
         "uc/verification/results/liu9-h2-twovar.json",
         "CERTIFIED_PROVED", True, "omit+nl",
         "a028ff17ae8492211a2ef9de5721de7f3a0afd8bd45bf7d8d9b77b7c90d2e6eb",
         "4ddbf0f38ca881dde204f086076424617edcdbe6a2826dfa945ef963ebdf70c6",
         "2afe8e242aa5f506c6d54967a133fe246a7e5cb8bad4961b9860e93b86f67163"),
    Link("boundary", "uc/liu9_h2_boundary_layer.py",
         "uc/verification/results/liu9-h2-boundary.json",
         "PROVED", True, "omit+nl",
         "e4d6d96df3ae435800cef5739fb1ae45b98febf539ab418df116ceb03bd481c9",
         "389804b249cd8b16a47f3c6ab64a70ad4a9afe78097b4a8c6739759804c9eaaa",
         "fcb3ed55e8bd9ad6d2f180f5ffa79f192dacc940a5c7506eca2a85aa612905e6"),
    Link("phi-audit", "uc/liu9_h2_phi_audit.py",
         "uc/verification/results/liu9-h2-phi-audit.json",
         "AUDIT-PASSED", True, "omit+nl",
         "f2359fa3b3b08bab29965253d777dcf32e4b1aee1e0bacbfb0924b86ac163996",
         "a515102c1807849b606bf4a681bd0a67173445b1a22610d397bb779440222cab",
         "b1d9ed3d9da9e7c99c026843b683447809fb38d28dba4c49a087729bf66cb50c"),
    Link("general-lift", "uc/liu9_h2_general_lift.py",
         "uc/verification/results/liu9-h2-general-lift.json",
         "PROVED", True, "omit+nl",
         "eb4874d39904593bddc3f2ae638d1bad6756a19d4ebf6e06a60427b6bb7ac1ce",
         "3ff1e797cd5dbcc4273d9f002b9e540b58829d237b2d344f266e0bea46c02bc8",
         "b265e726d110098a09c00de331c4433454973ba741966b3018b77afcc225b5c3"),
    Link("mixture-theorem", "uc/liu9_h2_mixture_theorem.py",
         "uc/verification/results/liu9-h2-mixture-theorem.json",
         "PROVED", True, "omit+nl",
         "329f7e2d71af8cd78d1a921c72b9ae4d05b71113134eb3341d69932f19ea24b3",
         "b108b781224dff601ffcf3f562703d3b620e7da6d89cf6d000b6ccbd06a4773b",
         "3098a1ca30a0f16582df02e1fb7dfd8fa27970c031167125e696e58ac43dfbab"),
    Link("four-fifths-ab", "uc/liu9_cprime_four_fifths_ab.py",
         "uc/verification/results/liu9-cprime-four-fifths-ab.json",
         "PROVED", True, "omit+nl",
         "9812fa9f64f99b8117da7bd6d31a3dfd5458deeaaddf764355012d04645939d4",
         "9a69041ae55a7d1e8e5d2fb3622753f323627a206cb18e101cd3eeabebad073e",
         "dbbd83a59e16ac44d48de67faef403f1318dc7fa6fb5824df5238e032691c6be"),
)
BY_NAME: Dict[str, Link] = {link.name: link for link in CHAIN}
BY_ARTIFACT_STEM: Dict[str, Link] = {Path(link.artifact).stem: link for link in CHAIN}
# A dependency block may name an artifact (hyphenated stem) or a module source
# (underscored file name); index both so cross-pin checks never silently skip.
BY_MODULE_NAME: Dict[str, Link] = {Path(link.module).name: link for link in CHAIN}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def internal_digest(report: dict, scope: str) -> str:
    body = dict(report)
    if scope == "empty":
        body["report_sha256"] = ""
    else:
        body.pop("report_sha256", None)
    blob = json.dumps(body, sort_keys=True, separators=(",", ":"))
    if scope == "omit+nl":
        blob += "\n"
    return sha256_bytes(blob.encode("utf-8"))


class Checker:
    def __init__(self) -> None:
        self.failures: List[str] = []
        self.passes = 0

    def expect(self, condition: bool, message: str) -> bool:
        if condition:
            self.passes += 1
        else:
            self.failures.append(message)
        return condition


def check_ledger_pins(chk: Checker) -> None:
    ledger = LEDGER.read_text(encoding="utf-8")
    for link in CHAIN:
        for kind, value in (("file", link.file_sha256), ("internal", link.report_sha256),
                            ("tool", link.tool_sha256)):
            chk.expect(value in ledger,
                       f"[A] {link.name}: pinned {kind} sha256 {value[:16]}... is not in PROGRESS.md")


def check_on_disk(chk: Checker) -> Dict[str, dict]:
    reports: Dict[str, dict] = {}
    for link in CHAIN:
        artifact = ROOT / link.artifact
        module = ROOT / link.module
        if not chk.expect(artifact.is_file(), f"[B] {link.name}: artifact missing {link.artifact}"):
            continue
        if not chk.expect(module.is_file(), f"[B] {link.name}: module missing {link.module}"):
            continue
        file_hash = sha256_path(artifact)
        chk.expect(file_hash == link.file_sha256,
                   f"[B] {link.name}: artifact sha256 on disk {file_hash[:16]}... != pinned {link.file_sha256[:16]}...")
        tool_hash = sha256_path(module)
        chk.expect(tool_hash == link.tool_sha256,
                   f"[B] {link.name}: module sha256 on disk {tool_hash[:16]}... != pinned {link.tool_sha256[:16]}...")
        report = json.loads(artifact.read_text(encoding="utf-8"))
        reports[link.name] = report
        chk.expect(report.get("claim_status") == link.claim_status,
                   f"[B] {link.name}: claim_status {report.get('claim_status')!r} != expected {link.claim_status!r}")
        stated = report.get("report_sha256")
        recomputed = internal_digest(report, link.digest_scope)
        chk.expect(stated == link.report_sha256,
                   f"[B] {link.name}: stated report_sha256 {str(stated)[:16]}... != pinned {link.report_sha256[:16]}...")
        chk.expect(recomputed == link.report_sha256,
                   f"[B] {link.name}: recomputed internal digest {recomputed[:16]}... != pinned {link.report_sha256[:16]}...")
        embedded_tool = report.get("tool_sha256")
        if embedded_tool is not None:
            chk.expect(embedded_tool == link.tool_sha256,
                       f"[B] {link.name}: embedded tool_sha256 {embedded_tool[:16]}... != module on disk {link.tool_sha256[:16]}...")
    # Cross-artifact pins: the audit pins the boundary certificate; the lift and the
    # mixture theorem pin every artifact they depend on.
    audit = reports.get("phi-audit")
    boundary = BY_NAME["boundary"]
    if audit is not None:
        expected = audit.get("module_hashes_expected", {})
        chk.expect(expected.get("file_sha256") == boundary.file_sha256
                   and expected.get("report_sha256") == boundary.report_sha256
                   and expected.get("tool_sha256") == boundary.tool_sha256,
                   "[B] phi-audit: module_hashes_expected does not pin the boundary certificate")
        chk.expect(audit.get("audit_outcome") == "AUDIT-PASSED",
                   f"[B] phi-audit: audit_outcome {audit.get('audit_outcome')!r}")
    resolved = 0
    for name in ("general-lift", "mixture-theorem", "four-fifths-ab"):
        report = reports.get(name)
        if report is None:
            continue
        for dep, block in report.get("dependencies", {}).items():
            chk.expect(block.get("match") is True, f"[B] {name}: dependency {dep} recorded match={block.get('match')!r}")
            # A block may be keyed by artifact stem, by a role name carrying the path, or
            # by a module file name; and it may spell the pins with or without the
            # `_expected` suffix.  Resolve every shape, so no cross-pin is skipped.
            path = block.get("path")
            stem = Path(str(path)).stem if isinstance(path, str) else None
            base = Path(str(path)).name if isinstance(path, str) else None
            pinned_sha = block.get("sha256_expected", block.get("sha256"))
            link = BY_ARTIFACT_STEM.get(dep) or (BY_ARTIFACT_STEM.get(stem) if stem else None)
            if link is not None:
                chk.expect(pinned_sha == link.file_sha256,
                           f"[B] {name}: dependency {dep} pins artifact {str(pinned_sha)[:16]}... != chain pin {link.file_sha256[:16]}...")
                pinned_status = block.get("claim_status_expected", block.get("claim_status"))
                chk.expect(pinned_status == link.claim_status,
                           f"[B] {name}: dependency {dep} expects claim_status {pinned_status!r}")
                resolved += 1
                continue
            link = BY_MODULE_NAME.get(dep) or (BY_MODULE_NAME.get(base) if base else None)
            if link is not None:
                chk.expect(pinned_sha == link.tool_sha256,
                           f"[B] {name}: dependency {dep} pins module {str(pinned_sha)[:16]}... != chain pin {link.tool_sha256[:16]}...")
                resolved += 1
    # Every one of these three reports pins at least one chain link; a resolver that
    # silently matches nothing is itself a failure.
    chk.expect(resolved >= 8,
               f"[B] cross-artifact pins: only {resolved} dependency blocks resolved to a chain link")
    return reports


def run_module(link: Link, run_index: int, scratch: Path, log: List[str]) -> Tuple[Optional[bytes], float, int]:
    """Run one module once; return (produced bytes, wall seconds, returncode)."""
    module = ROOT / link.module
    artifact = ROOT / link.artifact
    cmd = [sys.executable, "-I", "-B", str(module)]
    if link.output_flag:
        target = scratch / f"{link.name}.run{run_index}.json"
        cmd += ["--output", str(target)]
        snapshot = None
    else:
        # Fixed-path module: keep the pinned bytes on disk in scratch before the run,
        # so they survive even if this driver itself is killed mid-run.
        target = artifact
        snapshot = scratch / f"{link.name}.pinned.json"
        snapshot.write_bytes(artifact.read_bytes())
    env = dict(os.environ)
    env.setdefault("PYTHONHASHSEED", "0")
    produced: Optional[bytes] = None
    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, cwd=str(ROOT), env=env, capture_output=True, text=True)
        elapsed = time.monotonic() - start
        produced = target.read_bytes() if target.is_file() else None
    finally:
        if snapshot is not None:
            # Never leave anything but the pinned bytes in place of the artifact.
            original = snapshot.read_bytes()
            current = artifact.read_bytes() if artifact.is_file() else None
            if current != original:
                (scratch / f"{link.name}.run{run_index}.json").write_bytes(current or b"")
                artifact.write_bytes(original)
    tail = (proc.stdout.strip().splitlines() or [""])[-1]
    log.append(f"    run {run_index}: rc={proc.returncode} {elapsed:8.2f}s  {tail[:100]}")
    if proc.returncode != 0:
        log.append("    stderr: " + proc.stderr.strip()[-600:].replace("\n", "\n            "))
    return produced, elapsed, proc.returncode


def check_replay(chk: Checker, runs: int, scratch: Path, summary: Dict[str, dict]) -> None:
    for link in CHAIN:
        artifact = ROOT / link.artifact
        on_disk = artifact.read_bytes() if artifact.is_file() else None
        log: List[str] = []
        times: List[float] = []
        hashes: List[str] = []
        for k in range(1, runs + 1):
            produced, elapsed, rc = run_module(link, k, scratch, log)
            times.append(elapsed)
            hashes.append(sha256_bytes(produced) if produced is not None else "missing")
            chk.expect(rc == 0, f"[C] {link.name}: run {k} exited {rc}")
            chk.expect(produced is not None and produced == on_disk,
                       f"[C] {link.name}: run {k} output ({hashes[-1][:16]}...) is not byte-identical to the on-disk artifact")
        print(f"  {link.name:<16} pinned file {link.file_sha256[:16]}...  runs {' '.join(h[:12] for h in hashes)}")
        for line in log:
            print(line)
        summary[link.name] = {"artifact": link.artifact, "module": link.module,
                              "file_sha256": link.file_sha256, "report_sha256": link.report_sha256,
                              "tool_sha256": link.tool_sha256, "run_sha256": hashes,
                              "byte_identical": all(h == link.file_sha256 for h in hashes),
                              "seconds": [round(t, 2) for t in times]}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--runs", type=int, default=2, help="replays per module (default 2)")
    parser.add_argument("--report", type=Path, default=None,
                        help="optional JSON log of the replay (wall times included; not a certificate)")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.runs < 1:
        raise SystemExit("--runs must be >= 1")
    chk = Checker()
    print(f"REPLAY_H2_CHAIN root={ROOT} interpreter={sys.executable} runs={args.runs}")
    print("[A] pins present in PROGRESS.md")
    check_ledger_pins(chk)
    print("[B] on-disk artifacts, digests, module hashes, cross-artifact pins")
    check_on_disk(chk)
    print("[C] replays")
    summary: Dict[str, dict] = {}
    scratch = Path(tempfile.mkdtemp(prefix="replay_h2_chain_"))
    completed = False
    try:
        check_replay(chk, args.runs, scratch, summary)
        completed = True
    finally:
        if completed and not chk.failures:
            shutil.rmtree(scratch, ignore_errors=True)
        else:
            print(f"scratch kept (pinned snapshots and any divergent outputs): {scratch}")
    for message in chk.failures:
        print("FAIL", message)
    verdict = "PASS" if not chk.failures else f"FAIL ({len(chk.failures)} failures)"
    print(f"checks passed: {chk.passes}  failed: {len(chk.failures)}")
    print(f"REPLAY_H2_CHAIN {verdict}")
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps({"verdict": verdict, "runs": args.runs,
                                           "checks_passed": chk.passes, "failures": chk.failures,
                                           "links": summary}, indent=1, sort_keys=True) + "\n",
                               encoding="utf-8")
    return 0 if not chk.failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
