"""Erasure-channel DEM export for campaign interop (agreed with qldpc-dec 2026-08-29).

The study's decode path NEVER goes through a DEM: the erasure mask differs per
shot and the paper's decoder updates posteriors per shot (erased -> 0.5,
non-erased -> 1e-10), which a static DetectorErrorModel cannot express. DEMs
are exported ONLY as interop/metadata artifacts under
campaigns/<UTC>_<uuid>_<hash>/dem/, named N{N}_p{p}.dem.

Faithful marginal flattening of the paper's channel (arXiv:2603.19062v3
Sec 2.2: erased w.p. p, then independent X and Z errors w.p. 0.5 each):
conditional on erasure, P(X-only) = P(Z-only) = P(Y) = 1/4, so each qubit
carries three mechanisms at p/4: X -> flips Z-checks (Hz columns) and Z-logical
observables; Z -> flips X-checks (Hx) and X-logical observables; Y -> both.
(X,Z joint-erasure correlations are not representable in a DEM; documented
here, not hidden.)

Detector rows: X-check rows first (indices 0..n_x-1 = rows of Hx, detecting
Z-errors), then Z-check rows (n_x..n_x+n_z-1 = rows of Hz, detecting
X-errors). Observable rows: Lx rows (X-logicals, flipped by Z-errors) then Lz
rows (Z-logicals, flipped by X-errors) — the conventional CSS ordering used
by stim circuits generated for CSS codes.
"""

from __future__ import annotations

import numpy as np


def erasure_dem(Hx: np.ndarray, Hz: np.ndarray, Lx: np.ndarray, Lz: np.ndarray, p: float):
    """Build a stim.DetectorErrorModel for the marginal erasure channel."""
    import stim

    n = Hx.shape[1]
    n_x, n_z = Hx.shape[0], Hz.shape[0]
    n_obs = Lx.shape[0] + Lz.shape[0]
    dem = stim.DetectorErrorModel()

    def targets_for(which: str, q: int) -> list:
        t = []
        if which in ("X", "Y"):  # flips Z-checks (rows of Hz) and Z-logicals (Lz)
            for r in np.nonzero(Hz[:, q])[0]:
                t.append(stim.target_relative_detector_id(n_x + int(r)))
            for r in np.nonzero(Lz[:, q])[0]:
                t.append(stim.target_logical_observable_id(int(r) + Lx.shape[0]))
        if which in ("Z", "Y"):  # flips X-checks (rows of Hx) and X-logicals (Lx)
            for r in np.nonzero(Hx[:, q])[0]:
                t.append(stim.target_relative_detector_id(int(r)))
            for r in np.nonzero(Lx[:, q])[0]:
                t.append(stim.target_logical_observable_id(int(r)))
        return t

    qid = 0
    for q in range(n):
        qid += 1
        for which in ("X", "Z", "Y"):
            t = targets_for(which, q)
            if t:
                dem.append("error", p / 4.0, t)
            else:  # gauge-like qubit touching no detector: record as zero-weight
                dem.append("error", 0.0, [])
    # shift_info None; detector/observable counts declared implicitly. Assert:
    assert dem.num_detectors == n_x + n_z, (dem.num_detectors, n_x + n_z)
    assert dem.num_observables == n_obs, (dem.num_observables, n_obs)
    return dem
