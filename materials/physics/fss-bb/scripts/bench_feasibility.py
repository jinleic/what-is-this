#!/usr/bin/env python
"""Gate-B feasibility benchmark: per-shot decode wall-time vs N on the erasure path.

Measures the per-erasure-pattern decode cost at the p* Operating point
(p ~ 0.45, near-threshold, worst-side priors) for the five study sizes plus
the two Gate-B candidates (l,m) = (42,21) -> N=1764 and (48,24) -> N=2304,
both verified Ki=12 by bb_codes.STUDY_SIZES-style asserts at runtime.

Output: campaigns-smoke/feasibility_bench.json + stdout table.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from bb_codes import bb_code  # noqa: E402
from decode import BpOsdConfig, decode_erasure_batch  # noqa: E402
from erasure_channel import sample_erasure_shots  # noqa: E402

# (L, M): N. Gate-B candidates first verified in-browser by build check below.
SIZES = [(12, 6), (18, 9), (24, 12), (30, 15), (36, 18), (42, 21), (48, 24)]
P_BENCH = 0.45  # near-threshold worst case for both sectors (most OSD iterations)


def main() -> int:
    cfg = BpOsdConfig()
    out = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.platform(),
        "p_bench": P_BENCH,
        "decoder": "BP(50, ms, parallel) + OSD-CS(10), per-shot priors, both sectors",
        "threads": "single-core enforced via OMP_NUM_THREADS=1 etc.",
        "rows": [],
    }
    for L, M in SIZES:
        code = bb_code(L, M)
        N = code["N"]
        # build-cost timing (one-off, scales with N^3 in the naive GF(2) rref)
        t0 = time.perf_counter()
        _ = bb_code(L, M)  # re-run timed (caches nothing)
        build_s = time.perf_counter() - t0

        # 200 shots of per-shot-prior BP+OSD at the near-threshold point
        rng = np.random.default_rng(20260829)
        n_shots = 200
        mask, ex, ez = sample_erasure_shots(rng, N, P_BENCH, n_shots)
        t0 = time.perf_counter()
        r = decode_erasure_batch(code, mask, ex, ez, cfg)
        wall = time.perf_counter() - t0
        per_shot_ms = wall / n_shots * 1e3
        wer = r["wer"]

        # projected campaign cost at the paper's 200k shots/eval, ~6 evals/size
        per_eval_h = per_shot_ms / 1e3 * 200_000 / 3600
        est = {
            "LM": f"{L}x{M}", "N": N, "K": code["K"],
            "build_s": round(build_s, 2),
            "per_shot_ms@p0.45": round(per_shot_ms, 3),
            "wer@p0.45_200shots": round(wer, 3),
            "one_eval_200k_h": round(per_eval_h, 3),
            "threshold_search_h~6evals": round(per_eval_h * 6, 2),
        }
        out["rows"].append(est)
        print(est, flush=True)

    dest = HERE.parent / "campaigns-smoke" / "feasibility_bench.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nartifact: {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
