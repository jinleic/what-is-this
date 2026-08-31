"""Diagnose the 6.4x mean-timing gap: distribution + protocol forensics.

Gate A timing closed with p99.9 REPRODUCED (1.13x) but mean MISSED (6.4x)
at p=1e-3. This script isolates WHY by decomposing the single-decode time
distribution into its regime structure, then testing the four named
hypotheses from PROGRESS.md:

  H1 (call floor): per-decode call has a fixed Python/binding floor; the
     mean is floor-dominated, the tail is solver-dominated.
  H2 (solution path): failures (nontrivial OSD work) cost differently
     from successes; heavy failures inflate the tail only.
  H3 (OMP/threading): omp_thread_count / BLAS threads inflate the mean
     via spin-waits; tail unaffected (single dominant solve).
  H4 (allocator/cache): first-touch pages + cold caches inflate a
     transient initial (warmup) segment, dragging the mean.

Protocol mirrors run_gate.run_timing_point (same seed derivation, same
DEM, same decoder kwargs) so numbers are comparable to frozen campaign
20260830T071750Z_00327d68_a05388ffcf42.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qldpc_dec.bp_osd import make_bp_osd_decoder  # noqa: E402
from qldpc_dec.circuits import circuit_sha, dem_sha, load_circuit, load_dem  # noqa: E402
from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.seeds import sampling_seed  # noqa: E402


def timing_distribution(p: float, basis: str, num: int = 4000) -> dict:
    """Time single decodes and split the distribution by regime."""
    circuit = load_circuit(p, basis)
    dem = load_dem(p, basis)
    seed = sampling_seed(20260829, p, basis, "timing")
    det, obs = circuit.compile_detector_sampler(seed=seed).sample(
        num, separate_observables=True
    )

    # same prediction mapping as the batch adapter: correction vector ->
    # observable frame via the merged-convention A matrix
    _H, A, _p = dem_to_matrices(dem, merge=True)
    A = A.tocsr()

    dec = make_bp_osd_decoder(dem, max_iter=30)
    times = np.empty(num)
    preds = np.empty((num, A.shape[0]), dtype=np.int64)
    for i in range(num):
        t0 = time.perf_counter()
        corr = dec.decode(det[i].astype(np.uint8))
        times[i] = time.perf_counter() - t0
        preds[i] = (A @ (np.asarray(corr).astype(np.int64) % 2)) % 2

    fails = (preds != obs).any(axis=1)
    ms = times * 1e3

    stats: dict = {}
    qs = [1, 5, 10, 25, 50, 75, 90, 95, 99, 99.9, 100]
    stats["percentiles_ms"] = {q: float(np.percentile(ms, q)) for q in qs}
    stats["mean_ms"] = float(ms.mean())
    stats["median_ms"] = float(np.median(ms))
    stats["std_ms"] = float(ms.std())
    stats["n_failures"] = int(fails.sum())

    # bimodality: largest gap in sorted times
    s = np.sort(ms)
    gaps = np.diff(s)
    gi = int(np.argmax(gaps))
    stats["bimodal_gap_after_ms"] = float(s[gi])
    stats["bimodal_gap_width_ms"] = float(gaps[gi])
    stats["frac_below_gap"] = float((gi + 1) / num)
    stats["mean_ms_below_gap"] = float(s[: gi + 1].mean())
    stats["mean_ms_above_gap"] = float(s[gi + 1 :].mean())

    # success vs failure shot timing (H2)
    stats["mean_ms_success_shots"] = float(ms[~fails].mean())
    stats["mean_ms_failed_shots"] = (
        float(ms[fails].mean()) if fails.any() else None
    )
    stats["p999_ms_success_only"] = float(np.percentile(ms[~fails], 99.9))

    # warmup / autocorrelation on successes (H4)
    sol = ms[~fails]
    if len(sol) > 1100:
        stats["lag1_autocorr_successes"] = float(
            np.corrcoef(sol[:-1], sol[1:])[0, 1]
        )
        stats["first100_mean_ms"] = float(sol[:100].mean())
        stats["last1000_mean_ms"] = float(sol[-1000:].mean())
    return {"p": p, "basis": basis, "num": num, **stats}


def omp_probe(p: float, basis: str = "Z", num: int = 600) -> dict:
    """Record ldpc's own thread knob as configured."""
    circuit = load_circuit(p, basis)
    dem = load_dem(p, basis)
    dec = make_bp_osd_decoder(dem, max_iter=30)
    return {
        "ldpc_omp_thread_count_attr": int(getattr(dec, "omp_thread_count", -1)),
        "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS", "?"),
    }


if __name__ == "__main__":
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
    basis = sys.argv[2] if len(sys.argv) > 2 else "Z"
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    res = timing_distribution(p, basis, n)
    res.update(omp_probe(p, basis))
    res["circuit_sha256_12"] = circuit_sha(p, basis)[:12]
    res["dem_sha256_12"] = dem_sha(p, basis)[:12]
    print(json.dumps(res, indent=1))
