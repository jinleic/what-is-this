"""DEMatrices: (H, A, p) triples extracted from stim DEMs.

This stim build (1.16.0) does not expose DEM.to_checker_matrix(), so the
canonical decoding objects are built via ldpc's
detector_error_model_to_check_matrices
(ldpc.ckt_noise.dem_matrices), which returns:

    check_matrix           -> H (detectors x error mechanisms, hyperedges)
    observables_matrix     -> A (observables x error mechanisms)
    priors                 -> p (probability of each error mechanism)

IMPORTANT subtlety handled here: that helper MERGES identical hyperedge
columns and accumulates their probabilities (P(at least one)). The
beam-search paper's (H, A, p) convention keeps stim's raw error mechanisms
(8784 columns for the gross Z-memory DEM, one per circuit fault location
after decomposition). Both conventions are exposed:

    dem_to_matrices(dem, merge=False) -- raw stim mechanisms (default)
    dem_to_matrices(dem, merge=True)  -- ldpc merged convention

The BP+OSD module uses merge=True (the ldpc decoder consumes a plain
check matrix); the beam-search and NMS modules use the raw convention so
that per-mechanism priors match the DEM exactly.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import stim


def _raw_dem_arrays(dem: stim.DetectorErrorModel):
    """H, A, lam directly from stim error instructions (no column merging).

    Mirrors the (H, A, p) definition in arXiv:2512.07057 Appendix A:
    H[:, j] has 1s where error mechanism j flips detectors; A[:, j] where
    j flips logical observables; p_j is the mechanism probability.
    """
    num_dets = dem.num_detectors
    num_obs = dem.num_observables
    rows, cols = [], []
    obs_rows, obs_cols = [], []
    pvec = []
    j = 0
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        args = inst.args_copy()
        p = float(args[0]) if args else 0.0
        pvec.append(p)
        dets, obs = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                obs.append(t.val)
        for d in dets:
            rows.append(d)
            cols.append(j)
        for o in obs:
            obs_rows.append(o)
            obs_cols.append(j)
        j += 1
    H = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, cols)),
        shape=(num_dets, j),
    )
    A = sp.csc_matrix(
        (np.ones(len(obs_rows), dtype=np.uint8), (obs_rows, obs_cols)),
        shape=(num_obs, j),
    )
    lam = np.log((1.0 - np.asarray(pvec)) / np.clip(np.asarray(pvec), 1e-12, None))
    H = H.sorted_indices()
    A = A.sorted_indices()
    return H.tocsr(), A.tocsr(), lam


def _merged_dem_matrices(dem: stim.DetectorErrorModel):
    """H, A, p via ldpc's detector_error_model_to_check_matrices.

    Uses the edge decomposition (hyperedge_to_edge) conversion: the check
    matrix given to BP+OSD must be graphlike for OSD-CS; the
    allow_undecomposed_hyperedges=True mode keeps undecomposed Y-type
    hyperedges in H (they stay as weighted hyperedges -- the standard
    ldpc-stack convention for BB DEMs).
    """
    from ldpc.ckt_noise.dem_matrices import detector_error_model_to_check_matrices

    mats = detector_error_model_to_check_matrices(
        dem, allow_undecomposed_hyperedges=True
    )
    H = mats.check_matrix.tocsr()
    A = mats.observables_matrix.tocsr()
    p = np.asarray(mats.priors)
    return H, A, p


def dem_to_matrices(dem: stim.DetectorErrorModel, merge: bool = False):
    """Return (H, A, p-or-lam) from a stim DEM.

    merge=False -> raw stim mechanisms; priors returned as LLRs lam.
    merge=True  -> ldpc merged mechanisms; priors returned as probabilities p.
    """
    if merge:
        return _merged_dem_matrices(dem)
    return _raw_dem_arrays(dem)


if __name__ == "__main__":
    dem = stim.Circuit(
        "X_ERROR(0.1) 0 1 \n M 0 1 \n DETECTOR rec[-2] rec[-1] \n"
        " DETECTOR rec[-1] \n OBSERVABLE_INCLUDE(0) rec[-1]"
    ).detector_error_model()
    H, A, lam = dem_to_matrices(dem)
    print("raw H shape", H.shape, "A shape", A.shape, "lam", lam)
    Hd, Ad, pd = dem_to_matrices(dem, merge=True)
    print("merged H shape", Hd.shape, "A shape", Ad.shape, "p", pd)
