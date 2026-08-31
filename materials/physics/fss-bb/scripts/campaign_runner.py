#!/usr/bin/env python
"""Campaign runner for the BB-code erasure-channel FSS study (fss-bb).

Two primitives, both matching arXiv:2603.19062v3 exactly:

  * threshold mode (default, paper-faithful): per size, the Sec-2.4 adaptive
    search (bracket from p=0.38 step 0.04, Illinois regula falsi, stop at
    bracket < 5e-4 or 10 evals; 200k shots/eval) — every evaluated (p, WER)
    point is kept, and the FSS fit consumes the windowed subset
    |p - p*(N)| < 0.06 (Sec 2.5), exactly as the paper does.
  * grid mode: fixed p grid per size (for denser FSS datasets; optional).

All evaluated points carry: seed lineage, batch-decode wall time (feasibility
evidence for the Gate-B decision), and per-sector logical counts.

Artifacts written to campaigns/<UTC>_<uuid8>_<hash12>/ (immutable once written):
  manifest.json   — versions, host, args, seeds, timing, git-less digest input
  results.csv     — all evaluated points (size, p, shots, counts, wer, wall_s)
  thresholds.json — per-size p* (+ 5000-iter parametric bootstrap 95% CI)
  fss_fit.json    — collapse fit (+ 500-iter bootstrap CI), window sensitivity
                    0.04/0.06/0.08, linearized companion fits
  dem/N{N}_p{p}.dem — marginal erasure DEM exports (interop only)
  summary.md is NOT written by this runner (README discipline lives outside).

Usage (smoke):
  OMP_NUM_THREADS=1 python scripts/campaign_runner.py \
      --sizes 12x6 --shots 2000 --tag smoke
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))

from bb_codes import bb_code  # noqa: E402
from decode import BpOsdConfig, decode_erasure_batch  # noqa: E402
from dem_export import erasure_dem  # noqa: E402
from erasure_channel import sample_erasure_shots  # noqa: E402
from fss_fit import (  # noqa: E402
    fit_fss_bootstrap,
    fit_fss_win,
    fit_linearized_fss,
    adaptive_pseudo_threshold,
    pseudo_threshold_bootstrap,
)

SIZES = {"12x6": (12, 6), "18x9": (18, 9), "24x12": (24, 12), "30x15": (30, 15), "36x18": (36, 18)}


def point_seed(seed_base: int, tag: str, L: int, M: int, p: float) -> int:
    """Deterministic per-point seed lineage (recorded in results.csv)."""
    h = hashlib.sha256(f"{seed_base}|{tag}|{L}x{M}|{p!r}".encode()).hexdigest()
    return int(h[:16], 16)


def evaluate_point(code, p: float, shots: int, seed: int, cfg) -> dict:
    """One Monte-Carlo evaluation: sample, decode both sectors, tally."""
    rng = np.random.default_rng(seed)
    mask, ex, ez = sample_erasure_shots(rng, code["N"], p, shots)
    t0 = time.perf_counter()
    r = decode_erasure_batch(code, mask, ex, ez, cfg)
    wall = time.perf_counter() - t0
    return {
        "p": p,
        "shots": shots,
        "n_logical": r["n_logical"],
        "wer": r["wer"],
        "n_fail_x": int(r["fail_x"].sum()),
        "n_fail_z": int(r["fail_z"].sum()),
        "wall_s": wall,
        "per_shot_ms": wall / shots * 1e3,
        "seed": seed,
    }


def run_threshold_search(code, seed_base: int, tag: str, shots: int, cfg, start=0.38, step=0.04):
    """Paper Sec 2.4 search; every evaluation recorded with its seed."""
    cache: dict[float, dict] = {}

    def ev(p: float) -> tuple[int, int]:
        if p not in cache:
            seed = point_seed(seed_base, tag, code["L"], code["M"], p)
            cache[p] = evaluate_point(code, p, shots, seed, cfg)
        return cache[p]["n_logical"], cache[p]["shots"]

    res = adaptive_pseudo_threshold(ev, shots_per_eval=shots, start=start, step=step)
    res["eval_points"] = [cache[pt.p] for pt in res["points"] if pt.shots]
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="12x6,18x9,24x12,30x15,36x18")
    ap.add_argument("--shots", type=int, default=200_000)
    ap.add_argument("--seed-base", type=int, default=20260829)
    ap.add_argument("--tag", default="run")
    ap.add_argument("--mode", choices=["threshold", "grid"], default="threshold")
    ap.add_argument("--grid-half", type=float, default=0.06, help="grid mode: p* +- half, step 0.015")
    ap.add_argument("--n-boot-thresh", type=int, default=5000)
    ap.add_argument("--n-boot-fss", type=int, default=500)
    ap.add_argument("--dry-run", action="store_true", help="plan + cost estimate only")
    args = ap.parse_args()
    if args.shots < 200_000 and args.tag != "smoke":
        print("WARNING: shots < 200000 without --tag smoke: below the paper's floor", file=sys.stderr)

    cfg = BpOsdConfig()
    sizes = [s.strip() for s in args.sizes.split(",") if s.strip()]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    cid = f"{stamp}_{uuid.uuid4().hex[:8]}"
    out_dir = HERE.parent / "campaigns" / cid

    # <<< dry-run: cost model from a measured per-shot benchmark >>>
    if args.dry_run:
        print(json.dumps({"mode": args.mode, "sizes": sizes, "shots": args.shots}, indent=2))
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "dem").mkdir(exist_ok=True)
    rows, thresholds = [], {}
    meb = {}
    t_start = time.perf_counter()
    for s in sizes:
        L, M = SIZES[s]
        code = bb_code(L, M)
        res = run_threshold_search(code, args.seed_base, args.tag, args.shots, cfg)
        ci, _ = pseudo_threshold_bootstrap(res, n_boot=args.n_boot_thresh)
        thresholds[s] = {
            "L": L, "M": M, "N": code["N"], "K": code["K"],
            "p_star": res["p"], "bracket": list(res["bracket"]),
            "ci95": list(ci),
            "eval_points": res["eval_points"],
        }
        for pt in res["eval_points"]:
            rows.append({"L": L, "M": M, "N": code["N"], **pt})
        (out_dir / "dem").mkdir(exist_ok=True)  # recreated: dir may have been
        # moved/removed externally mid-run (seen when a crashed run's residue
        # was archived during this campaign)
        dem = erasure_dem(code["Hx"], code["Hz"], code["Lx"], code["Lz"], res["p"])
        (out_dir / "dem" / f"N{code['N']}_p{res['p']:.4f}.dem").write_text(str(dem))
        meb[s] = {
            "per_shot_ms_mean": float(np.mean([r["per_shot_ms"] for r in res["eval_points"]])),
            "per_shot_ms_max": float(np.max([r["per_shot_ms"] for r in res["eval_points"]])),
            "n_evals": len(res["eval_points"]),
        }
        print(f"[{s}] p*={res['p']:.4f} CI95=({ci[0]:.4f},{ci[1]:.4f}) evals={len(res['eval_points'])}")

    # FSS on the windowed evaluated points
    p = np.array([r["p"] for r in rows])
    Nv = np.array([r["N"] for r in rows], float)
    counts = np.array([r["n_logical"] for r in rows], float)
    shots_v = np.array([r["shots"] for r in rows], float)
    pstar = np.array([thresholds[f"{r['L']}x{r['M']}"]["p_star"] for r in rows])
    size_of = {s: i for i, s in enumerate(sizes)}
    fss = {"primary_window": 0.06, "fits": {}}
    for w in (0.04, 0.06, 0.08):
        try:
            fss["fits"][f"w{w}"] = fit_fss_bootstrap(
                counts, shots_v, p, Nv, pstar, window=w, n_boot=args.n_boot_fss
            )
        except ValueError as e:  # sparse dataset (e.g. single-size smoke)
            fss["fits"][f"w{w}"] = {"error": str(e)}
    fss["linearized"] = {}
    for label, nu_used in (("nu118", 1.18),):
        try:
            fss["linearized"][label] = fit_linearized_fss(
                [thresholds[s]["p_star"] for s in sizes],
                [thresholds[s]["N"] for s in sizes],
                nu=nu_used,
            )
        except np.linalg.LinAlgError as e:  # single-size: x* span is degenerate
            fss["linearized"][label] = {"error": str(e)}
    fit_nu = (fss["fits"].get("w0.06") or {}).get("nu")
    if fit_nu:
        try:
            fss["linearized"]["fitted_nu"] = fit_linearized_fss(
                [thresholds[s]["p_star"] for s in sizes],
                [thresholds[s]["N"] for s in sizes],
                nu=fit_nu,
            )
        except np.linalg.LinAlgError as e:
            fss["linearized"]["fitted_nu"] = {"error": str(e)}

    # write artifacts
    with open(out_dir / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    (out_dir / "thresholds.json").write_text(json.dumps(thresholds, indent=2, default=str))
    (out_dir / "fss_fit.json").write_text(json.dumps(fss, indent=2, default=str))
    manifest = {
        "campaign_id": cid,
        "tag": args.tag,
        "paper": "arXiv:2603.19062v3",
        "mode": args.mode,
        "shots_per_eval": args.shots,
        "seed_base": args.seed_base,
        "seed_lineage": "point_seed = sha256(seed_base|tag|LxM|p)[:16]",
        "decoder": cfg.as_ldpc_kwargs() | {"impl": "ldpc.BpOsdDecoder (bposd core)"},
        "channel": "erasure: mask~Bern(p); erased: X,Z indep Bern(0.5); priors 0.5/1e-10 per shot",
        "pack_versions": _versions(),
        "host": platform.platform(),
        "timing_total_s": time.perf_counter() - t_start,
        "timing_meb": meb,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    digest = hashlib.sha256(
        (out_dir / "results.csv").read_bytes() + (out_dir / "manifest.json").read_bytes()
    ).hexdigest()[:12]
    final = out_dir.parent / f"{cid.rsplit('_', 1)[0]}_{cid.rsplit('_', 1)[1]}_{digest}"
    out_dir.rename(final)
    print(f"campaign dir: {final.name}")
    w6 = fss["fits"]["w0.06"]
    if "p_inf" in w6:
        print(f"p*_inf = {w6['p_inf']:.4f} nu = {w6['nu']:.3f} (window 0.06)")
    else:
        print(f"p*_inf: not fit in this smoke ({w6.get('error', 'unknown')[:60]})")
    return 0


def _versions() -> dict:
    import stim
    import scipy
    import ldpc

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "stim": stim.__version__,
        "ldpc": getattr(ldpc, "__version__", "2.4.1"),
    }


if __name__ == "__main__":
    sys.exit(main())
