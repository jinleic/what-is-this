"""GB8 - post-hoc re-targeting of the frozen GB7 paired decisions.

pre_statement.md Revision GB8 (2026-09-02). No sampling, no decoding. Reads the
FROZEN GB7 summary (campaigns/20260901T145247Z_b7ea9ac4_ac3f6689e03b, terminal
verdict FROZEN-NEGATIVE) and re-applies GB7's own three-condition rule
(McNemar width effect AND paired interval excludes 1 AND paired interval
excludes the published reciprocal) against the published factor for the
configuration that was actually run.

Why: the original Gate B section of pre_statement.md (frozen 2026-08-29,
"Pinned published quantities" table) pinned beam64_640iters -> "17x lower",
and GB5a's band [11, 22] and GB7's PUBLISHED_RECIPROCALS inherited it.
arXiv:2512.07057 Table 1 lists FOUR configurations; Section III assigns 5.6x
to beam32_340iters, 7.0x to beam64_640iters (num_results=1) and 17x to
beam64_32res_640iters (num_results=32). Our beam64 rung is (64, 40, 30, 20,
num_results=1) = beam64_640iters, so its published factor is 7.0x, not 17x.
The 17x configuration was never run. Identical text in v1 (2025-12-08) and v2
(2025-12-17): a prereg transcription error, not a version drift.

The GB7 pre-registered decision and every frozen artifact are untouched; this
module writes src/evidence/gb8_retarget.json, and verify_gateb_gb5.py
re-derives every number in it from the frozen summary.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import MappingProxyType

TARGET = Path(__file__).resolve().parents[1]
GB7_CAMPAIGN_ID = "20260901T145247Z_b7ea9ac4_ac3f6689e03b"
GB7_CAMPAIGN = TARGET / "campaigns" / GB7_CAMPAIGN_ID
ARTIFACT = TARGET / "src" / "evidence" / "gb8_retarget.json"

PAPER = MappingProxyType({
    "arxiv": "2512.07057",
    "versions_read": ["v1 2025-12-08", "v2 2025-12-17"],
    "table1_rows": {
        "beam8_230iters": {"max_rounds": 10, "beam_width": 8, "initial_iters": 30,
                           "iters_per_round": 20, "num_results": 1},
        "beam32_340iters": {"max_rounds": 10, "beam_width": 32, "initial_iters": 40,
                            "iters_per_round": 30, "num_results": 1},
        "beam64_640iters": {"max_rounds": 20, "beam_width": 64, "initial_iters": 40,
                            "iters_per_round": 30, "num_results": 1},
        "beam64_32res_640iters": {"max_rounds": 20, "beam_width": 64, "initial_iters": 40,
                                  "iters_per_round": 30, "num_results": 32},
    },
    # Section III, paragraph 2 (Fig. 2 discussion), p=1e-3, vs bp30+osd.
    "section3_factors_vs_bposd": {
        "beam8_230iters": 1.3, "beam32_340iters": 5.6,
        "beam64_640iters": 7.0, "beam64_32res_640iters": 17.0,
    },
    # Abstract: "With a beam width of 8, we reach the same logical error rate
    # as BP-OSD" -- conflicts with Section III's 1.3x. Recorded, not resolved.
    "abstract_beam8_factor_vs_bposd": 1.0,
    "published_uncertainty": "none: no shots, failure counts, or error bars anywhere",
})

# Harness rung name -> paper configuration name (parameters pinned in the
# frozen GB5a/GB7 manifests: cpp_config_line ... num_results=1).
RUNG_TO_PAPER = MappingProxyType({"beam32": "beam32_340iters", "beam64": "beam64_640iters"})
# Targets as pinned on 2026-08-29 and applied by GB7 (frozen in its summary).
GB7_TARGETS = MappingProxyType({"beam32": 5.6, "beam64": 17.0})


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_hash(campaign: Path, name: str) -> str:
    for line in (campaign / "sha256s.txt").read_text().splitlines():
        digest, _, fname = line.strip().partition("  ")
        if fname == name:
            return digest
    raise FileNotFoundError(f"{name} not listed in {campaign / 'sha256s.txt'}")


def apply_rule(decision: dict, published_factor: float) -> dict:
    """GB7's three-condition rule against an arbitrary published factor.

    Uses only frozen quantities: McNemar width-effect flag, the paired ratio
    interval, and its exclusion of 1. Mirrors gateb_gb7_paired.rung_decision.
    """
    lo, hi = decision["paired_ratio_rung_over_beam8"]["ci95"]
    width_effect = bool(decision["mcnemar"]["width_effect_confirmed"])
    excludes_one = bool(decision["interval_excludes_one"])
    reciprocal = 1.0 / published_factor
    published_excluded = bool(reciprocal < lo or reciprocal > hi)
    satisfied = [label for label, ok in (
        ("WIDTH_EFFECT_CONFIRMED", width_effect),
        ("INTERVAL_EXCLUDES_ONE", excludes_one),
        ("PUBLISHED_FACTOR_EXCLUDED", published_excluded),
    ) if ok]
    if width_effect and excludes_one and published_excluded:
        outcome = "GAP_CONFIRMED"
    elif satisfied:
        outcome = "NO_GAP_" + "_AND_".join(satisfied)
    else:
        outcome = "UNRESOLVED_AT_N"
    measured = float(decision["measured_factor_beam8_over_rung"])
    return {
        "published_factor": float(published_factor),
        "published_reciprocal": reciprocal,
        "published_factor_excluded": published_excluded,
        "reciprocal_inside_ci95": not published_excluded,
        "measured_factor_minus_published": measured - float(published_factor),
        "outcome": outcome,
    }


def retarget_all(summary: dict) -> dict:
    rungs: dict[str, dict] = {}
    for rung, paper_name in RUNG_TO_PAPER.items():
        d = summary["decisions"][rung]
        corrected = float(PAPER["section3_factors_vs_bposd"][paper_name])
        as_run = float(GB7_TARGETS[rung])
        frozen = {
            "paired_cells": dict(d["paired_cells"]),
            "ratio_rung_over_beam8": d["paired_ratio_rung_over_beam8"]["ratio"],
            "ci95": list(d["paired_ratio_rung_over_beam8"]["ci95"]),
            "mcnemar_exact_two_sided_p": d["mcnemar"]["exact_two_sided_p"],
            "width_effect_confirmed": d["mcnemar"]["width_effect_confirmed"],
            "interval_excludes_one": d["interval_excludes_one"],
            "measured_factor_beam8_over_rung": d["measured_factor_beam8_over_rung"],
            "measured_factor_ci95": list(d["measured_factor_ci95"]),
            "gb7_outcome": d["outcome"],
        }
        as_run_check = apply_rule(d, as_run)
        if as_run_check["outcome"] != d["outcome"]:
            raise AssertionError(
                f"{rung}: rule replay under the target as run gives "
                f"{as_run_check['outcome']} but the frozen GB7 outcome is {d['outcome']}")
        # Sensitivity: the paper's own beam8 denominator conflict. Under the
        # Section III reading (beam8 = 1.3x better than bp30+osd) the
        # paper-implied rung/beam8 factor is F/1.3. Derived from two curve
        # read-outs with no stated uncertainty: a recorded conflict, not a
        # verdict.
        implied = corrected / float(PAPER["section3_factors_vs_bposd"]["beam8_230iters"])
        sens = apply_rule(d, implied)
        sens["direction"] = (
            "measured factor exceeds implied published factor"
            if sens["measured_factor_minus_published"] > 0
            else "measured factor below implied published factor")
        rungs[rung] = {
            "paper_configuration": paper_name,
            "paper_parameters": dict(PAPER["table1_rows"][paper_name]),
            "frozen_gb7": frozen,
            "gb7_target_as_run": {
                "note": "target pinned 2026-08-29 and applied by GB7 (frozen)",
                **as_run_check},
            "corrected_target": {
                "note": "published factor for the configuration actually run "
                        "(Table 1 + Section III), beam8 == bp30+osd denominator "
                        "per the 2026-08-29 pre-statement and the abstract",
                **apply_rule(d, corrected)},
            "sensitivity_beam8_1p3x_denominator": {
                "note": "Section III beam8 = 1.3x reading; implied rung/beam8 factor = F/1.3",
                "implied_factor_vs_beam8": implied,
                **sens},
        }
    return {
        "revision": "GB8",
        "kind": "post-hoc re-targeting of frozen decisions; no sampling; no decoding",
        "gb7_campaign": GB7_CAMPAIGN_ID,
        "paper": dict(PAPER),
        "rungs": rungs,
        "not_tested": {
            "beam64_32res_640iters": {
                **dict(PAPER["table1_rows"]["beam64_32res_640iters"]),
                "published_factor_vs_bposd": 17.0,
                "status": "NOT_TESTED_configuration_never_run",
                "harness_gap": "beam8_cpp pins num_results=1 with no CLI flag",
            }
        },
        "what_changes": (
            "interpretation only: the GB7 beam64 GAP_CONFIRMED refutes the "
            "mis-attributed 17x target, not the paper's claim for beam64_640iters"
        ),
        "what_does_not_change": (
            "every frozen artifact, the GB7 pre-registered rule and its "
            "FROZEN-NEGATIVE verdict, all failure counts, cells and intervals"
        ),
    }


def main() -> None:
    campaign = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else GB7_CAMPAIGN
    status = json.loads((campaign / "status.json").read_text())
    if not str(status.get("verdict", "")).startswith("FROZEN-"):
        raise SystemExit(f"refusing to re-target a non-frozen campaign: {status}")
    summary_path = campaign / "summary.json"
    listed = frozen_hash(campaign, "summary.json")
    actual = sha256_file(summary_path)
    if listed != actual:
        raise SystemExit(f"summary.json hash drift: listed {listed} actual {actual}")
    summary = json.loads(summary_path.read_text())
    payload = retarget_all(summary)
    payload["source"] = {
        "campaign": campaign.name,
        "terminal_verdict": status["verdict"],
        "summary_sha256": actual,
    }
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT}")
    print(f"sha256 {sha256_file(ARTIFACT)}")
    for rung, r in payload["rungs"].items():
        c = r["corrected_target"]
        print(f"{rung}: {r['paper_configuration']} published {c['published_factor']}x "
              f"reciprocal {c['published_reciprocal']:.4f} in CI95 {r['frozen_gb7']['ci95']} "
              f"-> {c['outcome']}")


if __name__ == "__main__":
    main()
