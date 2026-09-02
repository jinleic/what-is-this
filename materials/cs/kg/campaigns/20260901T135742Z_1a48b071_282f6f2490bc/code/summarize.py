#!/usr/bin/env python3
"""Summarize frozen per-tile artifacts after all run phases complete."""
from __future__ import annotations
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
LOGS = RUN / "logs"


def load(name):
    p = LOGS / name
    return json.loads(p.read_text()) if p.exists() else None


def compact_tiles(progress):
    rows = []
    if not progress:
        return rows
    for r in progress["tiles_completed"]:
        rows.append({
            "tile": r["tile_index"], "c_lo": r["c_lo"], "c_hi": r["c_hi"],
            "passed": r["passed"], "worst_margin": r.get("certified_margin_min") if r["passed"] else r.get("open_margin"),
            "worst_margin_radius": r.get("certified_margin_min_radius") if r["passed"] else r.get("open_margin_radius"),
            "worst_cell": r.get("worst_certified_cell") if r["passed"] else r.get("open_cell"),
            "max_depth": r.get("max_depth_evaluated"),
            "max_or_open_panels": r.get("worst_certified_panels") if r["passed"] else r.get("open_panels"),
            "fallback_runs": r.get("fallback_runs", 0),
            "panel_histogram": r.get("panel_histogram", {}),
            "envelope_evals": r.get("envelope_evals"),
            "process_seconds": r.get("process_seconds"),
        })
    return rows


def main():
    primary = load("band_1p75_3p5_result.json")
    if primary is None:
        raise SystemExit("primary result absent")
    pprog = load("band_1p75_3p5_progress.json")
    close1 = load("closeout_1p0_1p3_result.json")
    close2 = load("closeout_1p45_1p75_result.json")
    reason = str((primary.get("frontier") or {}).get("open_reason", ""))
    if primary["all_passed"]:
        verdict = "FROZEN-CERTIFIED"
    elif "budget" in reason or "ceiling" in reason:
        verdict = "FROZEN-INCONCLUSIVE"
    else:
        verdict = "FROZEN-NEGATIVE"
    out = {
        "verdict": verdict,
        "gate": "band-1p75-3p5-corr-cap2p17",
        "instrument_sha256": "94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30",
        "controls": {
            "startup": load("gate_battery.json"),
            "reconciliation": load("gate_reconciliation.json"),
            "near_miss": load("gate_near_miss_flip.json"),
            "planted_ratio_1p05": load("gate_planted_1p05_tile.json"),
        },
        "primary": primary,
        "primary_tiles": compact_tiles(pprog),
        "closeout_1p0_1p3": close1,
        "closeout_1p0_1p3_tiles": compact_tiles(load("closeout_1p0_1p3_progress.json")),
        "closeout_1p45_1p75": close2,
        "closeout_1p45_1p75_tiles": compact_tiles(load("closeout_1p45_1p75_progress.json")),
        "campaign_budget": load("campaign_budget.json"),
        "refutation_gate_fired": False,
        "claim_4p083_6p0": "NONE — Campaign C unchanged and out of scope",
    }
    (RUN / "result.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"verdict": verdict,
                      "primary_pass": primary["all_passed"],
                      "prefix": primary["certified_prefix_right"],
                      "tiles": primary["tiles_processed"],
                      "cpu_seconds": primary["cpu_seconds"]}, indent=2))

if __name__ == "__main__":
    main()
