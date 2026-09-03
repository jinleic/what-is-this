"""GB9 pre-registration sizing (pre_statement.md Revision GB9).

Computed BEFORE any GB9 sample from the frozen GB7 rates only. Two
hypotheses for the beam64_32res rung, both expressed as per-shot paired-cell
probabilities against beam8 (reference) on the shared p=1e-3 Z stream:

  H_PAPER  beam64_32res is 17x better than bp30+osd (== beam8 per the
           abstract and Revision GB1): rung/beam8 ratio R = 1/17, with the
           frozen GB7 beam64 discordance structure (c/f_rung = 6/27).
  H_NULL   beam64_32res is no better than beam64_640iters: the frozen GB7
           beam64 cells verbatim (both 21, b 202, c 6 per 1e8).

For each candidate N we draw the four-cell table, apply the SAME paired
interval code GB7 used (Wilson when c=0, four-cell bootstrap otherwise; the
bootstrap draw count is reduced for Monte-Carlo speed and stated), and
report the probability of each pre-registered outcome. Wall time uses the
measured beam_nr_cpp cost at num_results=32 (7.926 ms/shot on 24 threads,
188 ms/shot/thread, p=1e-3, 2000-shot gate stream).

Writes src/evidence/gb9_sizing.json.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

TARGET = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TARGET / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qldpc_dec.bootstrap import wilson_interval  # noqa: E402

GB7_CAMPAIGN = TARGET / "campaigns" / "20260901T145247Z_b7ea9ac4_ac3f6689e03b"
ARTIFACT = TARGET / "src" / "evidence" / "gb9_sizing.json"
PUBLISHED_FACTOR = 17.0
BAND = (PUBLISHED_FACTOR / 1.3, PUBLISHED_FACTOR * 1.3)   # +-30% pre-registered band
MS_PER_SHOT_24T = 7.926
K1_RATIO = 27.0 / 223.0     # frozen GB7 beam64_640iters / beam8 paired ratio
BEAM8_MS_PER_SHOT_24T = 0.35   # GB7: 1e8 beam8 shots in ~5.5 h on 26 threads => 0.2 ms; margin
REPLICATES = 2000
BOOT_DRAWS = 4000
CANDIDATES = (30_000_000, 50_000_000, 70_000_000, 100_000_000)


def interval(cells: np.ndarray, shots: int, rng: np.random.Generator) -> tuple[float, float, float]:
    both, b, c, _ = (int(x) for x in cells)
    f8 = both + b
    fr = both + c
    if f8 == 0:
        return float("nan"), float("nan"), float("nan")
    ratio = fr / f8
    if c == 0:
        lo, hi = wilson_interval(both, f8)
        return ratio, float(lo), float(hi)
    samples = rng.multinomial(shots, cells / shots, size=BOOT_DRAWS)
    den = samples[:, 0] + samples[:, 1]
    num = samples[:, 0] + samples[:, 2]
    valid = den > 0
    boot = num[valid] / den[valid]
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return ratio, float(lo), float(hi)


def mcnemar_exact(b: int, c: int) -> float:
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    from math import comb
    tail = sum(comb(n, i) for i in range(0, k + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


def outcomes(cells: np.ndarray, shots: int, rng: np.random.Generator) -> dict[str, bool]:
    both, b, c, _ = (int(x) for x in cells)
    ratio, lo, hi = interval(cells, shots, rng)
    if not np.isfinite(ratio):
        return {"defined": False}
    width = mcnemar_exact(b, c) < 0.05 and b > c
    excl_one = hi < 1.0 or lo > 1.0
    recip = 1.0 / PUBLISHED_FACTOR
    published_excluded = recip < lo or recip > hi
    band_lo, band_hi = 1.0 / BAND[1], 1.0 / BAND[0]       # on the ratio scale
    inside = band_lo <= lo and hi <= band_hi
    outside = hi < band_lo or lo > band_hi
    return {
        "defined": True,
        "gap_confirmed": bool(width and excl_one and published_excluded),
        "published_not_excluded": bool(not published_excluded),
        "k1_ratio_excluded": bool(K1_RATIO < lo or K1_RATIO > hi),
        "band_reproduced": bool(inside),
        "band_not_reproduced": bool(outside),
        "band_inconclusive": bool(not inside and not outside),
    }


def main() -> None:
    summary = json.loads((GB7_CAMPAIGN / "summary.json").read_text())
    d64 = summary["decisions"]["beam64"]
    n8 = int(d64["paired_ratio_rung_over_beam8"]["beam8_failures"])       # 223
    cells64 = d64["paired_cells"]
    shots7 = int(summary["shots_shared"])
    p8 = n8 / shots7
    # H_NULL: frozen beam64 cells verbatim
    null = np.array([cells64["both_fail"], cells64["reference_fails_candidate_ok"],
                     cells64["reference_ok_candidate_fails"], 0], dtype=float) / shots7
    # H_PAPER: R = 1/17 with the frozen c/f_rung structure
    c_frac = cells64["reference_ok_candidate_fails"] / (
        cells64["both_fail"] + cells64["reference_ok_candidate_fails"])
    fr = p8 / PUBLISHED_FACTOR
    paper = np.array([fr * (1 - c_frac), p8 - fr * (1 - c_frac), fr * c_frac, 0.0])
    for h in (null, paper):
        h[3] = 1.0 - h[:3].sum()

    rng = np.random.default_rng(20260902)
    table = {}
    for shots in CANDIDATES:
        row = {}
        for name, probs in (("H_PAPER", paper), ("H_NULL", null)):
            acc: dict[str, int] = {}
            defined = 0
            for _ in range(REPLICATES):
                cells = rng.multinomial(shots, probs).astype(float)
                out = outcomes(cells, shots, rng)
                if not out["defined"]:
                    continue
                defined += 1
                for k, v in out.items():
                    if k != "defined":
                        acc[k] = acc.get(k, 0) + int(v)
            row[name] = {k: v / defined for k, v in acc.items()}
            row[name]["expected_beam8_failures"] = p8 * shots
            row[name]["expected_rung_failures"] = float(probs[0] + probs[2]) * shots
        row["wall_hours_32res_24threads"] = shots * MS_PER_SHOT_24T / 3.6e6
        row["wall_hours_beam8_24threads"] = shots * BEAM8_MS_PER_SHOT_24T / 3.6e6
        row["rss_gb_per_decode"] = shots * (117 + 8 + 24) / 1e9
        table[str(shots)] = row

    payload = {
        "revision": "GB9",
        "computed_from": {
            "gb7_campaign": GB7_CAMPAIGN.name,
            "gb7_summary_sha256": hashlib.sha256((GB7_CAMPAIGN / "summary.json").read_bytes()).hexdigest(),
            "beam8_failures_per_1e8": n8,
            "beam64_cells_per_1e8": cells64,
        },
        "hypotheses": {
            "H_PAPER": {"ratio_rung_over_beam8": 1.0 / PUBLISHED_FACTOR,
                        "cell_probabilities": paper.tolist(),
                        "c_fraction_of_rung_failures": c_frac},
            "H_NULL": {"ratio_rung_over_beam8": cells64["both_fail"] + cells64["reference_ok_candidate_fails"],
                       "cell_probabilities": null.tolist()},
        },
        "rule": {
            "published_factor": PUBLISHED_FACTOR,
            "band_factor": list(BAND),
            "gap": "McNemar p<0.05 and b>c AND CI excludes 1 AND CI excludes 1/17",
            "band": "three-way: CI inside [1/22.1, 1/13.1] REPRODUCED; wholly outside NOT_REPRODUCED; else INCONCLUSIVE",
        },
        "monte_carlo": {"replicates": REPLICATES, "bootstrap_draws": BOOT_DRAWS,
                        "seed": 20260902,
                        "note": "campaign uses 100000 draws; reduced here for speed"},
        "cost_model": {"ms_per_shot_32res_24threads": MS_PER_SHOT_24T,
                       "ms_per_shot_beam8_24threads": BEAM8_MS_PER_SHOT_24T},
        "candidates": table,
    }
    ARTIFACT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT}")
    print("sha256", hashlib.sha256(ARTIFACT.read_bytes()).hexdigest())
    for shots, row in table.items():
        hp, hn = row["H_PAPER"], row["H_NULL"]
        print(f"N={int(shots):>11,d}  wall={row['wall_hours_32res_24threads']:6.1f}h  "
              f"E[f8]={hp['expected_beam8_failures']:6.1f}  "
              f"H_PAPER: E[fr]={hp['expected_rung_failures']:5.1f} "
              f"P(17x not excluded)={hp['published_not_excluded']:.3f} "
              f"P(band REPRODUCED)={hp['band_reproduced']:.3f} | "
              f"H_NULL: E[fr]={hn['expected_rung_failures']:5.1f} "
              f"P(GAP_CONFIRMED)={hn['gap_confirmed']:.3f}")


if __name__ == "__main__":
    main()
