"""Inspect the corrected-prior BP+OSD single-decode timing distribution.

The original 2026-08-30 use of this script diagnosed a 6.4x mean gap under
the scalar-mean channel later invalidated by pre-statement Revision GB3.
Its inference that the paper used OSD0 semantics is retracted: with the full
heterogeneous channel, the published OSD-CS10 configuration reproduces both
mean and p99.9, while a focused corrected-prior OSD0 probe does not.

The live function uses the current full-vector decoder and asserts that
channel before timing. It does not reconstruct the invalidated historical
campaign; those measurements remain documented in ``physics/PROGRESS.md``.
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
from qldpc_dec.run_gate import IONQ_BPOSD_DECODER_CONFIG  # noqa: E402
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
    _H, A, priors = dem_to_matrices(dem, merge=True)
    A = A.tocsr()

    dec = make_bp_osd_decoder(dem, **IONQ_BPOSD_DECODER_CONFIG)
    configured = np.asarray(dec.channel_probs, dtype=np.float64)
    priors = np.asarray(priors, dtype=np.float64)
    np.testing.assert_array_equal(configured, priors)
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
    return {
        "p": p,
        "basis": basis,
        "num": num,
        "decoder": dict(IONQ_BPOSD_DECODER_CONFIG),
        "bposd_prior_check": {
            "columns": int(priors.size),
            "distinct_probabilities": int(np.unique(priors).size),
            "configured_channel_exact": True,
        },
        **stats,
    }


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
