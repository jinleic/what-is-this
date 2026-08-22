"""Parallel circuit-level Monte Carlo with a persistent per-worker decoder.

sinter rebuilds decoder state per batch, which dominates runtime for BP+OSD on
a dense hypergraph DEM (measured 2.5 core-seconds per shot versus 0.079 s for
a warm decoder).  This harness constructs the decoder exactly once per worker
process and then streams shots through it.

All timing numbers reported by :func:`collect` are wall-clock on a fixed
worker count so that decoder comparisons are normalised.
"""

from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, asdict

import numpy as np
import stim

__all__ = ["MonteCarloResult", "collect"]

_STATE: dict = {}


def _init(circuit_text: str, dec_cfg: dict) -> None:
    import numpy as np
    import scipy.sparse as sp
    import stim
    from ldpc import BpOsdDecoder

    from .bposd_dem import dem_to_matrices

    circ = stim.Circuit(circuit_text)
    dem = circ.detector_error_model(decompose_errors=False)
    M = dem_to_matrices(dem)
    pr = np.clip(M.priors, 1e-12, 1 - 1e-12)
    _STATE["circ"] = circ
    _STATE["L"] = M.L
    _STATE["dec"] = BpOsdDecoder(
        sp.csr_matrix(M.H), error_channel=list(pr),
        max_iter=dec_cfg["max_iter"], bp_method=dec_cfg["bp_method"],
        ms_scaling_factor=dec_cfg["ms_scaling_factor"],
        osd_method=dec_cfg["osd_method"], osd_order=dec_cfg["osd_order"],
    )


def _work(args) -> tuple[int, int, float, list[float]]:
    shots, seed = args
    circ = _STATE["circ"]
    dec = _STATE["dec"]
    L = _STATE["L"]
    sampler = circ.compile_detector_sampler(seed=seed)
    det, obs = sampler.sample(shots, separate_observables=True)
    det = det.astype(np.uint8)
    obs = obs.astype(np.uint8)
    fails = 0
    lat: list[float] = []
    t0 = time.perf_counter()
    for i in range(shots):
        s = time.perf_counter()
        corr = dec.decode(det[i])
        lat.append(time.perf_counter() - s)
        if np.any(((L @ corr) % 2).astype(np.uint8) != obs[i]):
            fails += 1
    return shots, fails, time.perf_counter() - t0, lat


@dataclass
class MonteCarloResult:
    shots: int
    failures: int
    ler_per_shot: float
    ci95_lo: float
    ci95_hi: float
    wall_s: float
    workers: int
    decode_s_per_shot_1core: float
    latency_p50_ms: float
    latency_p95_ms: float
    latency_p99_ms: float
    dem_errors: int
    dem_detectors: int

    def to_dict(self) -> dict:
        return asdict(self)


def collect(circuit: stim.Circuit, *, max_shots: int, max_errors: int = 10**9,
            workers: int = 12, chunk: int = 250, seed: int = 20260811,
            max_iter: int = 30, osd_order: int = 0, osd_method: str = "osd0",
            bp_method: str = "ms", ms_scaling_factor: float = 0.625,
            ) -> MonteCarloResult:
    from .bposd_dem import clopper_pearson

    dem = circuit.detector_error_model(decompose_errors=False)
    cfg = dict(max_iter=max_iter, bp_method=bp_method,
               ms_scaling_factor=ms_scaling_factor,
               osd_method=osd_method, osd_order=osd_order)
    txt = str(circuit)
    total = fails = 0
    lat_all: list[float] = []
    core_s = 0.0
    t0 = time.time()
    rng = np.random.default_rng(seed)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init,
                             initargs=(txt, cfg)) as ex:
        while total < max_shots and fails < max_errors:
            batch = min(workers * 2, max(1, (max_shots - total) // chunk + 1))
            jobs = [(min(chunk, max_shots - total - i * chunk),
                     int(rng.integers(1, 2**31)))
                    for i in range(batch)]
            jobs = [j for j in jobs if j[0] > 0]
            if not jobs:
                break
            for s, f, cs, lat in ex.map(_work, jobs):
                total += s
                fails += f
                core_s += cs
                lat_all.extend(lat)
    wall = time.time() - t0
    lo, hi = clopper_pearson(fails, total)
    a = np.asarray(lat_all) * 1e3 if lat_all else np.zeros(1)
    return MonteCarloResult(
        shots=total, failures=fails, ler_per_shot=fails / total if total else float("nan"),
        ci95_lo=lo, ci95_hi=hi, wall_s=round(wall, 1), workers=workers,
        decode_s_per_shot_1core=core_s / total if total else float("nan"),
        latency_p50_ms=float(np.percentile(a, 50)),
        latency_p95_ms=float(np.percentile(a, 95)),
        latency_p99_ms=float(np.percentile(a, 99)),
        dem_errors=dem.num_errors, dem_detectors=dem.num_detectors,
    )
