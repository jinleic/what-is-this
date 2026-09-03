#!/usr/bin/env python3
"""Checkpointed Campaign-A envelope driver for KgBand175.

Prereg: cs/kg/prereg/20260902_closeout_1p45_1p75/pre_statement.md
(commit 82935f7; sha256 946d8b89... in provenance.json).
Only A's panel cap/c-range knobs are used. The envelope and all geometry
helpers are imported byte-identically from Campaign A's gate_corr_env.py.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util as ilu
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# Must precede every frozen-campaign import.
sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
LOGS = RUN / "logs"
LOGS.mkdir(parents=True, exist_ok=True)
MANIFEST = json.loads((RUN / "manifest.json").read_text())
CREATED = dt.datetime.strptime(MANIFEST["created_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)

EXPECTED_SHA = {
    "gate_corr_env.py": "94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30",
    "startup_controls.py": "cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3",
    "decompose_cell.py": "291eddc1d76b5a60338855012629077a7d5799db9c7434e07b18ae19df6204ff",
    "gate_subdiv.py": "6f1b8ee0f4bd11e3a3dec1c0a584dfe72975c9b91ab53ff2b26366cb2e0b5df5",
}
STARTUP_JSON_SHA = "211d8e9afd6491af62d7bb89282e586f300f5a32535320caa31aa342fc63beac"
PRIMARY_CAP = 131072
FALLBACK_N = 524288
BAND_EVAL_BUDGET = 260000
CPU_BUDGET_SECONDS = 60 * 3600
WALL_CEILING_SECONDS = 72 * 3600
A_TILE1_MARGIN = "1.04499016763434518123845756638e-6"
NEAR = {
    "cL": "3.16988277218086405264471296000",
    "cR": "3.23328042762448133369760721920",
    "a0lo": "-0.849527994791666666666667",
    "a0hi": "-0.849446614583333333333333",
    "a2lo": "0.527669270833333333333333",
    "a2hi": "0.527750651041666666666667",
}
PAD = "1e-20"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, obj: Any) -> None:
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    os.replace(tmp, path)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def instrument_identity() -> dict:
    rows = {}
    for name, expected in EXPECTED_SHA.items():
        actual = sha(HERE / name)
        rows[name] = {"expected": expected, "actual": actual, "ok": actual == expected}
    ok = all(row["ok"] for row in rows.values())
    rec = {"ok": ok, "files": rows, "checked_before_import": True}
    atomic_json(LOGS / "instrument_identity.json", rec)
    if not ok:
        raise SystemExit("frozen instrument sha mismatch; see instrument_identity.json")
    return rec


def load_instrument():
    # A's module registers/executes its frozen d4core dependency by the same
    # __file__-relative construction it used in Campaign A.
    spec = ilu.spec_from_file_location("gce", HERE / "gate_corr_env.py")
    gce = ilu.module_from_spec(spec)
    sys.modules["gce"] = gce
    spec.loader.exec_module(gce)
    d4core = sys.modules["d4core"]
    if gce.PANEL_CAP != PRIMARY_CAP:
        raise SystemExit(f"PANEL_CAP changed before import: {gce.PANEL_CAP}")
    return gce, d4core


def ball(v) -> dict:
    return {"ball": v.str(30), "radius": v.rad().str(12)}


def box_json(gce, b) -> list[str]:
    return [x.str(30) for x in b]


def stack_json(gce, stack) -> list[dict]:
    return [{"box": gce.box_key(b), "depth": depth,
             "coordinates": box_json(gce, b)} for b, depth in stack]


def budget_record() -> dict:
    path = LOGS / "campaign_budget.json"
    if path.exists():
        return read_json(path)
    return {
        "cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "wall_ceiling_seconds": WALL_CEILING_SECONDS,
        "controls_cpu_seconds": 0.0,
        "tile_cpu_seconds_by_tag": {},
        "cpu_used_seconds": 0.0,
    }



def save_budget(rec: dict) -> None:
    rec["cpu_used_seconds"] = (rec.get("controls_cpu_seconds", 0.0) +
                               rec.get("excluded_launch_cpu_seconds", 0.0) +
                               sum(rec.get("tile_cpu_seconds_by_tag", {}).values()))
    rec["cpu_remaining_seconds"] = CPU_BUDGET_SECONDS - rec["cpu_used_seconds"]
    rec["wall_used_seconds"] = (dt.datetime.now(dt.timezone.utc) - CREATED).total_seconds()
    rec["wall_remaining_seconds"] = WALL_CEILING_SECONDS - rec["wall_used_seconds"]
    atomic_json(LOGS / "campaign_budget.json", rec)


class BudgetGuard:
    def __init__(self, base_cpu: float):
        self.base_cpu = base_cpu
        self.session_start = time.process_time()

    def total_cpu(self) -> float:
        return self.base_cpu + (time.process_time() - self.session_start)

    def wall_used(self) -> float:
        return (dt.datetime.now(dt.timezone.utc) - CREATED).total_seconds()

    def reason(self) -> str | None:
        if self.total_cpu() >= CPU_BUDGET_SECONDS:
            return "cpu budget exhausted"
        if self.wall_used() >= WALL_CEILING_SECONDS:
            return "wall ceiling exhausted"
        return None


# --------------------------------------------------------------- controls

def run_startup_battery() -> dict:
    before_child = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.time()
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, str(HERE / "startup_controls.py")],
        cwd=RUN, env=env, capture_output=True, text=True)
    after_child = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = ((after_child.ru_utime - before_child.ru_utime) +
           (after_child.ru_stime - before_child.ru_stime))
    out_path = LOGS / "startup_controls.json"
    if not out_path.exists():
        raise SystemExit("startup_controls.py produced no JSON")
    data = read_json(out_path)
    actual_sha = sha(out_path)
    ok = (proc.returncode == 0 and data.get("aborted") is None and
          len(data.get("controls", {})) == 16 and actual_sha == STARTUP_JSON_SHA)
    rec = {
        "ok": ok, "controls_count": len(data.get("controls", {})),
        "aborted": data.get("aborted"), "sha256": actual_sha,
        "expected_sha256": STARTUP_JSON_SHA, "process_seconds": cpu,
        "wall_seconds": time.time() - started, "stdout": proc.stdout,
        "stderr": proc.stderr,
    }
    atomic_json(LOGS / "gate_battery.json", rec)
    if not ok:
        raise SystemExit("startup battery failed or was not byte-identical")
    return rec


def reconciliation(gce, d4core) -> dict:
    from flint import arb, fmpq
    started = time.process_time()
    with d4core.Prec(gce.PREC):
        cL, cR = arb("4.083"), arb("4.115664")
        box = (arb(fmpq(-13, 16)), arb(fmpq(-77, 96)),
               arb(fmpq(55, 96)), arb(fmpq(7, 12)))
        upper = gce.cell_envelope(cL, box, PRIMARY_CAP, gce.PREC)
        margin = gce.strict_margin_lower(cR, upper, gce.PREC)
        frozen = arb(A_TILE1_MARGIN)
        diff = margin - frozen
        threshold = arb(2) ** -40
        ok = bool(margin.lower() > 0) and bool(diff.rad() < threshold)
    rec = {
        "ok": ok, "c_pair": [cL.str(30), cR.str(30)],
        "cell": box_json(gce, box), "panels": PRIMARY_CAP,
        "upper": ball(upper), "margin": ball(margin),
        "frozen_margin": A_TILE1_MARGIN, "difference": ball(diff),
        "threshold_2^-40": threshold.str(20),
        "process_seconds": time.process_time() - started,
    }
    atomic_json(LOGS / "gate_reconciliation.json", rec)
    if not ok:
        raise SystemExit("instrument-identity reconciliation failed")
    return rec


def near_box(gce, d4core):
    from flint import arb
    with d4core.Prec(gce.PREC):
        eps = arb(PAD)
        return (arb(NEAR["a0lo"]) - eps, arb(NEAR["a0hi"]) + eps,
                arb(NEAR["a2lo"]) - eps, arb(NEAR["a2hi"]) + eps)


def near_miss(gce, d4core) -> dict:
    from flint import arb
    started = time.process_time()
    with d4core.Prec(gce.PREC):
        cL, cR = arb(NEAR["cL"]), arb(NEAR["cR"])
        b = near_box(gce, d4core)
        rows = {}
        for n, need in ((32768, "negative"), (PRIMARY_CAP, "positive")):
            t0 = time.process_time()
            upper = gce.cell_envelope(cL, b, n, gce.PREC)
            margin = gce.strict_margin_lower(cR, upper, gce.PREC)
            sign_ok = bool(margin.upper() < 0) if need == "negative" else bool(margin.lower() > 0)
            rows[str(n)] = {"need": need, "ok": sign_ok,
                            "upper": ball(upper), "margin": ball(margin),
                            "process_seconds": time.process_time() - t0}
    ok = rows["32768"]["ok"] and rows[str(PRIMARY_CAP)]["ok"]
    rec = {"ok": ok, "c_pair": [cL.str(30), cR.str(30)],
           "superset_box": box_json(gce, b), "pad": PAD,
           "panels": rows, "process_seconds": time.process_time() - started}
    atomic_json(LOGS / "gate_near_miss_flip.json", rec)
    if not ok:
        raise SystemExit("near-miss flip failed; run adapted D1-D6 obstruction battery")
    return rec


def planted(gce, d4core) -> dict:
    from flint import arb
    started = time.process_time()
    with d4core.Prec(gce.PREC):
        cL, cR = arb(NEAR["cL"]), arb(NEAR["cL"]) * arb("1.05")
        b = near_box(gce, d4core)
        upper = gce.cell_envelope(cL, b, PRIMARY_CAP, gce.PREC)
        margin = gce.strict_margin_lower(cR, upper, gce.PREC)
        ok = bool(margin.upper() < 0)
    rec = {"ok": ok, "ratio": "1.05", "c_pair": [cL.str(30), cR.str(30)],
           "superset_box": box_json(gce, b), "panels": PRIMARY_CAP,
           "upper": ball(upper), "margin": ball(margin),
           "required": "margin upper endpoint < 0",
           "process_seconds": time.process_time() - started}
    atomic_json(LOGS / "gate_planted_1p05_tile.json", rec)
    if not ok:
        raise SystemExit("ratio-1.05 planted-failure control failed")
    return rec


def run_controls(gce, d4core) -> None:
    # Exact prereg order: startup -> reconciliation -> near-miss -> planted.
    started = time.process_time()
    battery = run_startup_battery()
    r = reconciliation(gce, d4core)
    n = near_miss(gce, d4core)
    p = planted(gce, d4core)
    parent_cpu = time.process_time() - started
    child_cpu = battery["process_seconds"]
    control_cpu = parent_cpu + child_cpu
    result = {"ok": all(x["ok"] for x in (battery, r, n, p)),
              "order": ["startup", "reconciliation", "near_miss", "ratio_1p05"],
              "process_seconds_parent": parent_cpu,
              "process_seconds_child": child_cpu,
              "process_seconds_total": control_cpu}
    atomic_json(LOGS / "controls_result.json", result)
    budget = budget_record()
    budget["controls_cpu_seconds"] = control_cpu
    save_budget(budget)
    print(f"CONTROLS PASS cpu={control_cpu:.3f}s", flush=True)


# ------------------------------------------------------------ DFS / budget

def evaluate_primary(gce, c_left, c_right, box, stats: dict, guard: BudgetGuard) -> dict:
    last_upper = None
    last_margin = None
    last_n = 0
    for n in gce.panel_ladder(c_left):
        if stats["envelope_evals"] >= stats["budget"]:
            return {"passed": False, "reason": "band budget reached",
                    "margin": last_margin, "upper": last_upper, "n_panels": last_n}
        reason = guard.reason()
        if reason is not None and last_margin is not None:
            return {"passed": False, "reason": reason,
                    "margin": last_margin, "upper": last_upper, "n_panels": last_n}
        upper = gce.cell_envelope(c_left, box, n, gce.PREC)
        margin = gce.strict_margin_lower(c_right, upper, gce.PREC)
        stats["envelope_evals"] += 1
        stats["panel_histogram"][str(n)] = stats["panel_histogram"].get(str(n), 0) + 1
        last_upper, last_margin, last_n = upper, margin, n
        if margin > 0:
            return {"passed": True, "reason": "strict endpoint separation",
                    "margin": margin, "upper": upper, "n_panels": n}
    return {"passed": False, "reason": "panel cap reached",
            "margin": last_margin, "upper": last_upper, "n_panels": last_n}


def certify_tile(gce, d4core, c_left, c_right, stats: dict, guard: BudgetGuard) -> dict:
    roots = gce.root_boxes(gce.N_SIDE)
    assert len(roots) == 132
    stack = [(b, 0) for b in reversed(roots)]
    ts = {"root_boxes": 132, "boxes_evaluated": 0, "certified_leaves": 0,
          "margin_min_ball": None, "margin_max_ball": None,
          "worst_certified_cell": None, "worst_certified_panels": None,
          "fallback_runs": 0, "max_depth_evaluated": 0}
    while stack:
        box, depth = stack.pop()
        ts["max_depth_evaluated"] = max(ts["max_depth_evaluated"], depth)
        result = evaluate_primary(gce, c_left, c_right, box, stats, guard)
        ts["boxes_evaluated"] += 1
        if result["passed"]:
            ts["certified_leaves"] += 1
            gce._update_leaf_extrema(ts, box, result)
            continue
        if result["reason"] == "panel cap reached" and gce.can_split(box, depth):
            children = gce.split_box(box)
            for child in reversed(children):
                if gce.box_intersects_disk(child):
                    stack.append((child, depth + 1))
            continue
        if result["reason"] == "panel cap reached":
            # Exactly the preregistered direct n=2^19 fallback rung.
            if stats["envelope_evals"] >= stats["budget"]:
                result["reason"] = "band budget reached before fallback"
            elif guard.reason() is not None:
                result["reason"] = guard.reason()
            else:
                with d4core.Prec(gce.PREC):
                    upper = gce.cell_envelope(c_left, box, FALLBACK_N, gce.PREC)
                    margin = gce.strict_margin_lower(c_right, upper, gce.PREC)
                stats["envelope_evals"] += 1
                stats["panel_histogram"][str(FALLBACK_N)] = stats["panel_histogram"].get(str(FALLBACK_N), 0) + 1
                ts["fallback_runs"] += 1
                result = {"passed": bool(margin > 0), "reason": "registered 2^19 fallback",
                          "margin": margin, "upper": upper, "n_panels": FALLBACK_N}
                if result["passed"]:
                    ts["certified_leaves"] += 1
                    gce._update_leaf_extrema(ts, box, result)
                    continue
        margin = result.get("margin")
        return {
            **ts, "passed": False, "open_reason": result["reason"],
            "open_cell": gce.box_key(box), "open_cell_coordinates": box_json(gce, box),
            "open_depth": depth, "open_panels": result.get("n_panels"),
            "open_margin": margin.str(30) if margin is not None else None,
            "open_margin_radius": margin.rad().str(12) if margin is not None else None,
            "open_margin_ball": ball(margin) if margin is not None else None,
            "remaining_stack_boxes": len(stack), "remaining_stack": stack_json(gce, stack),
            "budget_cpu_at_frontier": guard.total_cpu(),
            "budget_wall_at_frontier": guard.wall_used(),
        }
    return {
        **ts, "passed": True, "open_reason": None, "open_cell": None,
        "open_depth": None, "open_panels": None, "open_margin": None,
        "open_margin_radius": None, "open_margin_ball": None,
        "remaining_stack_boxes": 0, "remaining_stack": [],
    }


def hist_delta(after: dict[str, int], before: dict[str, int]) -> dict[str, int]:
    return {k: after.get(k, 0) - before.get(k, 0)
            for k in sorted(set(after) | set(before)) if after.get(k, 0) != before.get(k, 0)}


def band_spec(tag: str):
    specs = {
        "band_1p75_3p5": ("1.75", "3.5", 87),
        "closeout_1p0_1p3": ("1.0", "1.3", 33),
        "closeout_1p45_1p75": ("1.45", "1.75", 24),
    }
    return specs[tag]


def run_band(gce, d4core, tag: str, start_tile: int | None,
             end_tile: int | None = None) -> dict:
    c_lo, c_hi, expected = band_spec(tag)
    tiles = gce.c_tiles(c_lo, c_hi)
    gce.assert_tile_cover(c_lo, c_hi, tiles)
    assert len(tiles) == expected, (tag, len(tiles), expected)
    assert len(gce.root_boxes(12)) == 132
    progress_path = LOGS / f"{tag}_progress.json"
    if progress_path.exists():
        progress = read_json(progress_path)
        rows = progress["tiles_completed"]
        evals = int(progress["envelope_evals"])
        hist = {str(k): int(v) for k, v in progress["panel_histogram"].items()}
    else:
        rows, evals, hist = [], 0, {}
    if start_tile is None:
        start_tile = len(rows) + 1
    if start_tile != len(rows) + 1:
        raise SystemExit(f"resume mismatch: --start-tile {start_tile}, expected {len(rows)+1}")
    budget = budget_record()
    base_cpu = float(budget["cpu_used_seconds"])
    guard = BudgetGuard(base_cpu)
    stats = {"budget": BAND_EVAL_BUDGET, "envelope_evals": evals, "panel_histogram": hist}
    frontier = None
    stop = min(len(tiles), end_tile if end_tile is not None else len(tiles))
    for index in range(start_tile - 1, stop):
        if guard.reason() is not None:
            frontier = {"open_reason": guard.reason(), "tile_index": index + 1,
                        "c_pair": [tiles[index][0].str(30), tiles[index][1].str(30)],
                        "at_tile_boundary": True, "certified_prefix_tiles": len(rows)}
            atomic_json(LOGS / f"{tag}_frontier.json", frontier)
            break
        c_left, c_right = tiles[index]
        if tag == "band_1p75_3p5":
            assert gce.panel_ladder(c_left) == [1024, 4096, 16384, 65536, 131072]
        eval_before = stats["envelope_evals"]
        hist_before = dict(stats["panel_histogram"])
        t0 = time.process_time()
        tile = certify_tile(gce, d4core, c_left, c_right, stats, guard)
        tile_cpu = time.process_time() - t0
        row = gce.serialise_tile(c_left, c_right, tile)
        row.update({
            "tile_index": index + 1,
            "process_seconds": tile_cpu,
            "envelope_evals": stats["envelope_evals"] - eval_before,
            "panel_histogram": hist_delta(stats["panel_histogram"], hist_before),
            "fallback_runs": tile["fallback_runs"],
            "max_depth_evaluated": tile["max_depth_evaluated"],
            "remaining_stack": tile["remaining_stack"],
            "open_cell_coordinates": tile.get("open_cell_coordinates"),
            "open_margin_ball": tile.get("open_margin_ball"),
        })
        rows.append(row)
        atomic_json(LOGS / f"{tag}_tile_{index+1:03d}.json", row)
        budget = budget_record()
        budget.setdefault("tile_cpu_seconds_by_tag", {})[tag] = sum(
            float(r["process_seconds"]) for r in rows)
        save_budget(budget)
        prefix_right = c_right.str(30) if tile["passed"] else c_left.str(30)
        progress = {
            "tag": tag, "c_lo": c_lo, "c_hi": c_hi,
            "tiles_total": expected, "tiles_completed": rows,
            "tiles_completed_count": len(rows), "next_tile": index + 2,
            "certified_prefix_right": prefix_right,
            "envelope_evals": stats["envelope_evals"],
            "panel_histogram": stats["panel_histogram"],
            "campaign_budget": budget_record(),
        }
        atomic_json(progress_path, progress)
        verdict = "PASS" if tile["passed"] else "OPEN"
        margin = row.get("certified_margin_min") if tile["passed"] else row.get("open_margin")
        print(f"{tag} tile {index+1}/{expected} {verdict} boxes={row['boxes_evaluated']} "
              f"evals={row['envelope_evals']} cpu={tile_cpu:.3f}s margin={margin}", flush=True)
        if not tile["passed"]:
            frontier = row
            atomic_json(LOGS / f"{tag}_frontier.json", frontier)
            break
    all_passed = len(rows) == expected and all(bool(r["passed"]) for r in rows)
    prefix = c_hi if all_passed else (rows[-1]["c_hi"] if rows and rows[-1]["passed"] else
                                     (rows[-1]["c_lo"] if rows else c_lo))
    result = {
        "tag": tag, "c_lo": c_lo, "c_hi": c_hi,
        "tiles_total": expected, "tiles_processed": len(rows),
        "tiles_certified": sum(bool(r["passed"]) for r in rows),
        "all_passed": all_passed, "certified_prefix_right": prefix,
        "frontier": frontier, "envelope_evals": stats["envelope_evals"],
        "panel_histogram": stats["panel_histogram"],
        "cpu_seconds": sum(float(r["process_seconds"]) for r in rows),
        "campaign_budget": budget_record(),
    }
    atomic_json(LOGS / f"{tag}_result.json", result)
    return result


def controls_valid() -> bool:
    path = LOGS / "controls_result.json"
    return path.exists() and bool(read_json(path).get("ok"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("controls", "band"), required=True)
    parser.add_argument("--tag", choices=("band_1p75_3p5", "closeout_1p0_1p3", "closeout_1p45_1p75"))
    parser.add_argument("--start-tile", type=int)
    parser.add_argument("--end-tile", type=int)
    args = parser.parse_args()
    instrument_identity()
    gce, d4core = load_instrument()
    if args.phase == "controls":
        run_controls(gce, d4core)
        return
    if not controls_valid():
        raise SystemExit("controls_result.json absent/non-green; refusing band compute")
    if args.tag is None:
        raise SystemExit("--tag required for --phase band")
    if args.tag.startswith("closeout_"):
        primary = LOGS / "band_1p75_3p5_result.json"
        if not primary.exists() or not read_json(primary).get("all_passed"):
            raise SystemExit("close-out forbidden before primary [1.75,3.5] PASS")
    result = run_band(gce, d4core, args.tag, args.start_tile, args.end_tile)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
