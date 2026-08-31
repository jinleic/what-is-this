"""BP+OSD decoding paths for the BB erasure-channel study.

Decoder configuration matches arXiv:2603.19062v3 Sec 2.3 verbatim:
min-sum BP, OSD-CS post-processing, OSD order 10, max 50 BP iterations,
implemented with the same decoder core the paper uses (the bposd package is a
thin wrapper around ldpc's BpOsdDecoder; `ldpc` 2.4.1 is installed and
`bposd` is not, so we drive the core directly, with per-shot channel-prior
updates via `update_channel_probs`).

Two decode paths:

1. `decode_erasure_batch` — primary study path. Direct erasure sampling, two
   sector decoders (Hx on Z-errors, Hz on X-errors), per-shot priors
   (erased -> 0.5, non-erased -> 1e-10). Wall-clock timing instrumented for
   the Gate-B feasibility analysis.

2. `decode_dem_shots` — adapter matching the shared qldpc-dec harness shape
   contract (agreed 2026-08-29, recorded in ../README.md):
       dem    : stim.DetectorErrorModel
       shots  : np.ndarray[bool, shape=(num_shots, num_dem_det)]
       returns: np.ndarray[int, shape=(num_shots, num_obs)]
   so a later merge onto BBCodeHarness is a drop-in swap.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from ldpc import BpOsdDecoder

try:  # package import when shipped as fss-bb package
    from .erasure_channel import ERASED_PRIOR, NON_ERASED_PRIOR, sector_failures
except ImportError:  # flat sys.path usage inside src/
    from erasure_channel import ERASED_PRIOR, NON_ERASED_PRIOR, sector_failures


@dataclass(frozen=True)
class BpOsdConfig:
    """Paper-verbatim decoder settings (arXiv:2603.19062v3 Sec 2.3)."""

    bp_method: str = "ms"  # min-sum BP
    osd_method: str = "osd_cs"  # OSD-CS post-processing
    osd_order: int = 10  # OSD order 10 (Sec 2.3); mixed-channel used 5 (Sec 2.6)
    max_iter: int = 50  # maximum 50 BP iterations
    ms_scaling_factor: float = 1.0
    schedule: str = "parallel"

    def as_ldpc_kwargs(self) -> dict:
        return dict(
            bp_method=self.bp_method,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
            max_iter=self.max_iter,
            ms_scaling_factor=self.ms_scaling_factor,
            schedule=self.schedule,
        )


def _make_decoder(H: np.ndarray, cfg: BpOsdConfig, base_probs: np.ndarray) -> BpOsdDecoder:
    return BpOsdDecoder(H, channel_probs=base_probs, **cfg.as_ldpc_kwargs())


def _decode_sector(H: np.ndarray, errors: np.ndarray, mask: np.ndarray, cfg: BpOsdConfig):
    """Decode one sector. `errors` (num_shots, N) is the error vector for this
    sector (Z-errors for the Hx/X-sector decoder; X-errors for Z-sector).
    Priors per shot: erased -> 0.5, else 1e-10 (paper Sec 2.2 step 3).
    Returns (corrections uint8 (num_shots, N), wall_seconds, converged_flags).
    """
    num_shots, n = errors.shape
    dec = _make_decoder(H, cfg, np.full(n, ERASED_PRIOR))
    # syndrome of this sector's decoder = errors @ H^T (X-checks detect Z-errors)
    syn = (errors @ H.T % 2).astype(np.uint8)
    corr = np.zeros((num_shots, n), dtype=np.uint8)
    converged = np.zeros(num_shots, dtype=bool)
    t0 = time.perf_counter()
    probs = np.empty(n, dtype=np.float64)
    for s in range(num_shots):
        if mask is not None:
            np.copyto(probs, NON_ERASED_PRIOR)
            probs[mask[s]] = ERASED_PRIOR
            dec.update_channel_probs(probs)
        corr[s] = dec.decode(syn[s])
        try:
            converged[s] = dec.converge
        except Exception:
            converged[s] = True
    wall = time.perf_counter() - t0
    return corr, wall, converged


def decode_erasure_batch(code: dict, mask, ex, ez, cfg: BpOsdConfig | None = None) -> dict:
    """Full two-sector erasure decode of a batch of shots.

    code: dict from bb_codes.bb_code (Hx, Hz, Lz, Lx, ...).
    Returns dict with:
        wer          float  — fraction of shots failing either sector
        fail_x/fail_z (num_shots,) bool
        corrections  (2, num_shots, N) uint8 — [x-sector output (Z-type), z-sector output]
        wall_s       (2,) float — per-sector decode wall time
        per_shot_s   (2,) float — wall/num_shots per sector
        converged    (2, num_shots)
    """
    cfg = cfg or BpOsdConfig()
    Hx, Hz = code["Hx"], code["Hz"]
    # X sector: checks Hx, decodes Z-errors ez
    cz, wall_x, conv_x = _decode_sector(Hx, ez, mask, cfg)
    # Z sector: checks Hz, decodes X-errors ex
    cx, wall_z, conv_z = _decode_sector(Hz, ex, mask, cfg)
    fail_x, fail_z = sector_failures(code["Lz"], code["Lx"], Hx, Hz, ex, ez, cx, cz)
    ns = mask.shape[0]
    return {
        "wer": float((fail_x | fail_z).mean()),
        "fail_x": fail_x,
        "fail_z": fail_z,
        "n_logical": int((fail_x | fail_z).sum()),
        "corrections": np.stack([cz, cx]),
        "wall_s": (wall_x, wall_z),
        "per_shot_s": (wall_x / max(ns, 1), wall_z / max(ns, 1)),
        "converged": np.stack([conv_x, conv_z]),
    }


def sector_failures(Lz, Lx, Hx, Hz, ex, ez, cx, cz):
    """Evaluate sector failures given corrections (num_shots, N each).

    Residual wx = ex + cx (X-type): success iff Hz @ wx == 0 AND it is a
    trivial stabilizer coset. The coset test pairs the X-type residual against
    Z-type logical rows Lz (X·Z anticommutation); pairing against X-type rows
    (Lx) is meaningless. Failure iff Lz @ wx != 0.
    Residual wz = ez + cz (Z-type): failure iff Hx @ wz != 0 or Lx @ wz != 0.
    """
    wz = ez ^ cz
    wx = ex ^ cx
    fail_z = ((wz @ Hx.T % 2).any(axis=1)) | ((wz @ Lx.T % 2).any(axis=1))
    fail_x = ((wx @ Hz.T % 2).any(axis=1)) | ((wx @ Lz.T % 2).any(axis=1))
    return fail_x, fail_z


def _dem_matrices(dem) -> tuple[np.ndarray, np.ndarray]:
    """Dense parity-check matrix (num_det x n_mech) and observable-flip matrix
    (num_obs x n_mech) of a stim DetectorErrorModel."""
    num_det = dem.num_detectors
    num_obs = dem.num_observables
    H_cols = []
    obs_cols = []
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        col = np.zeros(num_det, dtype=np.uint8)
        oc = np.zeros(num_obs, dtype=np.uint8)
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                col[t.val] = 1
            elif t.is_logical_observable_id():
                oc[t.val] = 1
        H_cols.append(col)
        obs_cols.append(oc)
    H = np.stack(H_cols, axis=1) if H_cols else np.zeros((num_det, 0), np.uint8)
    O = np.stack(obs_cols, axis=1) if obs_cols else np.zeros((num_obs, 0), np.uint8)
    return H, O


def decode_dem_shots(
    dem,
    shots: np.ndarray,
    cfg: BpOsdConfig | None = None,
    channel_probs: np.ndarray | None = None,
) -> np.ndarray:
    """qldpc-dec BBCodeHarness shape contract:

        dem   : stim.DetectorErrorModel
        shots : bool (num_shots, dem.num_detectors)
        return: int (num_shots, dem.num_observables) in {0,1}

    BP+OSD on the DEM check matrix; observable prediction = correction applied
    to the DEM's observable-flip table.
    """
    cfg = cfg or BpOsdConfig()
    shots = np.asarray(shots, dtype=bool)
    H, O = _dem_matrices(dem)
    dec = _make_decoder(H.T, cfg, channel_probs if channel_probs is not None else np.full(H.shape[1], 1e-3))
    num_shots = shots.shape[0]
    preds = np.zeros((num_shots, O.shape[0]), dtype=int)
    for s in range(num_shots):
        corr = dec.decode(shots[s].astype(np.uint8))
        preds[s] = (O @ corr) % 2
    return preds
