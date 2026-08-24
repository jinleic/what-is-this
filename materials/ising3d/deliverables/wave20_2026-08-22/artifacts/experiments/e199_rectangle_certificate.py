#!/usr/bin/env python3
"""Assemble the exact boundary-row route-closure certificate.

Only the two owned producers e197 and e198 are imported.  In particular, no
numbered sibling-front module is used.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from e197_rectangle_words import run_word_audit
from e198_rectangle_induction import run_induction_audit

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "rectangle_twogen.json"
CPU_BUDGET_SECONDS = 120.0
RSS_CAP_BYTES = 2_000_000_000

SOURCE_PATHS = (
    ROOT / "experiments" / "e197_rectangle_words.py",
    ROOT / "experiments" / "e198_rectangle_induction.py",
    Path(__file__),
    ROOT / "tests" / "test_rectangle_twogen.py",
    ROOT / "proofs" / "rectangle_twogen.md",
    ROOT / "proofs" / "twogen_allsize.md",
    ROOT / "proofs" / "ladder_w8.md",
    ROOT / "proofs" / "dla_3xl.md",
)


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append(
        {
            "name": name,
            "passed": bool(passed),
            "detail": f"[COMPUTATION] {detail}",
        }
    )
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def run_all() -> tuple[dict[str, object], list[dict[str, object]], float]:
    started = time.process_time()
    word_data, word_checks = run_word_audit()
    induction_data, induction_checks = run_induction_audit()
    checks = [*word_checks, *induction_checks]

    word_shapes = {
        tuple(int(value) for value in anchor["core_shape"])
        for anchor in word_data["finite_anchors"]
    }
    restriction_widths = {
        int(anchor["width"]) for anchor in induction_data["finite_anchors"]
    }
    record(
        checks,
        "word_anchor_set",
        word_shapes == {(2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3)},
        f"bounded exact core shapes are {sorted(word_shapes)}",
    )
    record(
        checks,
        "restriction_anchor_set",
        restriction_widths == set(range(1, 7)),
        f"bounded exact ancillary widths are {sorted(restriction_widths)}",
    )
    record(
        checks,
        "child_checks_pass",
        all(bool(item["passed"]) for item in [*word_checks, *induction_checks]),
        f"all {len(word_checks) + len(induction_checks)} independently computed child checks pass",
    )
    record(
        checks,
        "source_files_present",
        all(path.is_file() for path in SOURCE_PATHS),
        f"all {len(SOURCE_PATHS)} owned and cited source files exist",
    )

    used = time.process_time() - started
    record(
        checks,
        "process_time_budget",
        used < CPU_BUDGET_SECONDS,
        f"process time {used:.6f}s is below {CPU_BUDGET_SECONDS:.0f}s",
    )
    observed_rss = max_rss_bytes()
    record(
        checks,
        "rss_budget",
        observed_rss < RSS_CAP_BYTES,
        f"peak RSS {observed_rss} bytes is below {RSS_CAP_BYTES} bytes",
    )

    data = {
        "setting": {
            "tag": "[EXTERNAL]",
            "generators": "A=sum_v X_v and B=sum_(uv in E) Z_u Z_v on the open r by L grid",
            "lie_algebra": "Lie_Q<iA,iB> on the full physical spin space",
            "normalized_bracket": "{S,T}=[S,T]/(2i)",
        },
        "rank_one_restriction": induction_data,
        "word_compression": word_data,
        "depth_three_identity": {
            "tag": "[THEOREM]",
            "quantifiers": "core height r>=1 and width L>=1",
            "statement": (
                "Phi({B_full,{A_full,B_full}})={B_core,{A_core,B_core}}+"
                "S_boundary+(3L-2)I for final-row |+> compression"
            ),
            "ordinary_commutator_form": (
                "Phi([B_full,[A_full,B_full]])=[B_core,[A_core,B_core]]-"
                "4 S_boundary-(12L-8)I"
            ),
        },
        "reflection_exclusion": {
            "tag": "[THEOREM]",
            "quantifiers": "core height r>=2 and width L>=1",
            "statement": (
                "The compressed depth-three word is not in g_(r,L)+Q I because every core "
                "Lie word is fixed by row reflection whereas S_boundary is not."
            ),
        },
        "conclusion": {
            "tag": "[THEOREM]",
            "headline": (
                "The natural ladder-to-grid row restriction fails for every rank-one frozen "
                "row, and |+>-row compression fails as a Lie map at an explicit depth-three word."
            ),
            "route_closed": (
                "inherit a ladder certificate by freezing one added row or by applying the same "
                "Lie words upstairs and taking their |+> partial matrix elements"
            ),
            "counterexample_scope": "all core heights r>=2 and all widths L>=1",
        },
        "dimension_status": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "No all-3xL or all-rectangle bound dim g_(r,L)>rL(2rL-1) is proved or disproved."
            ),
            "remaining_routes": (
                "direct stable-signature Lie words, higher-rank ancillary codes, smaller invariant "
                "submodules, and nonlocal representation restrictions"
            ),
            "finite_scope_guard": (
                "the Pauli anchors stop at full shape 4x3 and the ancillary enumeration at width 6; "
                "the all-size route closure comes from the symbolic identities in the proof note"
            ),
        },
        "method_separation": {
            "tag": "[THEOREM]",
            "statement": (
                "The obstruction uses row leakage and nonmultiplicativity of partial matrix "
                "elements, not ad_A Clifford-grade separation."
            ),
        },
        "external_context": {
            "tag": "[EXTERNAL]",
            "ladder": "proofs/ladder_w8.md proves the dimension obstruction for all open 2xL, L>=3",
            "finite_3xl": "proofs/dla_3xl.md gives finite 3xL anchors but no all-L induction",
            "grade_wall": "proofs/twogen_allsize.md closes ad_A-only grade separation",
        },
        "ledger_ids": ["H518", "H519", "H520", "H521", "H522", "H523"],
    }
    return data, checks, started


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]], started: float
) -> dict[str, object]:
    source_hashes = {
        str(path.relative_to(ROOT)): sha256_file(path) for path in SOURCE_PATHS
    }
    return {
        "meta": {
            "experiment": "e199_rectangle_certificate",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": "nice -n 10 .venv/bin/python experiments/e199_rectangle_certificate.py",
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": "exact integer Hermitian-Pauli brackets and exact rational elimination",
            "modular_arithmetic": "not used",
            "finite_anchor_policy": "bounded validation only; no finite row is promoted to an all-size dimension theorem",
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "observed_process_time_seconds": round(time.process_time() - started, 6),
                "observed_peak_rss_bytes": max_rss_bytes(),
            },
            "source_sha256": source_hashes,
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks, started = run_all()
    artifact = make_artifact(data, checks, started)
    record(
        artifact["checks"],
        "artifact_shape",
        set(artifact) == {"meta", "data", "checks"},
        "top-level keys are exactly meta, data, checks",
    )
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if all(bool(item["passed"]) for item in artifact["checks"]):
        print("PASS e199 rectangle two-generator certificate", flush=True)
        return 0
    print("FAIL e199 rectangle two-generator certificate", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
