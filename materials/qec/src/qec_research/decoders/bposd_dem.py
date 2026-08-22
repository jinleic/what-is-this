"""BP+OSD decoding directly from a Stim detector error model.

The DEM is converted to
    * check matrix  H  (num_detectors x num_errors)
    * observable matrix  L  (num_observables x num_errors)
    * prior probability vector  p
No graphlike decomposition is performed, so hyperedges (which BB and PBB
circuits produce in abundance) are preserved exactly and X/Z/Y correlations
are kept -- this is the correlation-aware path required for non-CSS codes.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import stim

__all__ = ["DemMatrices", "dem_to_matrices", "BposdDemDecoder", "clopper_pearson"]


@dataclass
class DemMatrices:
    H: np.ndarray            # (num_detectors, num_errors) uint8
    L: np.ndarray            # (num_observables, num_errors) uint8
    priors: np.ndarray       # (num_errors,) float
    num_detectors: int
    num_observables: int


def dem_to_matrices(dem: stim.DetectorErrorModel) -> DemMatrices:
    nd, no = dem.num_detectors, dem.num_observables
    dets: list[list[int]] = []
    obs: list[list[int]] = []
    prob: list[float] = []

    def handle(inst: stim.DemInstruction, shift: int) -> None:
        if inst.type != "error":
            return
        p = inst.args_copy()[0]
        ds, os_ = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                ds.append(t.val + shift)
            elif t.is_logical_observable_id():
                os_.append(t.val)
        dets.append(ds)
        obs.append(os_)
        prob.append(p)

    def walk(block: stim.DetectorErrorModel, shift: int) -> None:
        for inst in block:
            if isinstance(inst, stim.DemRepeatBlock):
                body = inst.body_copy()
                d_per = body.num_detectors
                for rep in range(inst.repeat_count):
                    walk(body, shift + rep * d_per)
            elif inst.type == "shift_detectors":
                shift += inst.targets_copy()[0] if inst.targets_copy() else 0
            else:
                handle(inst, shift)

    flat = dem.flattened()
    for inst in flat:
        handle(inst, 0)

    ne = len(prob)
    H = np.zeros((nd, ne), dtype=np.uint8)
    L = np.zeros((no, ne), dtype=np.uint8)
    for e, (ds, os_) in enumerate(zip(dets, obs)):
        for d in ds:
            H[d, e] ^= 1
        for o in os_:
            L[o, e] ^= 1
    return DemMatrices(H, L, np.asarray(prob, dtype=float), nd, no)


class BposdDemDecoder:
    """BP + ordered-statistics post-processing on the full hypergraph DEM."""

    def __init__(self, dem: stim.DetectorErrorModel, *, max_iter: int = 30,
                 osd_order: int = 10, osd_method: str = "osd_cs",
                 bp_method: str = "ms", ms_scaling_factor: float = 0.625):
        from ldpc import BpOsdDecoder

        self.mats = dem_to_matrices(dem)
        m = self.mats
        # guard against p=0 or p=1 priors which BP cannot use
        pr = np.clip(m.priors, 1e-12, 1 - 1e-12)
        self.dec = BpOsdDecoder(
            m.H.astype(np.uint8), error_channel=list(pr),
            max_iter=max_iter, bp_method=bp_method,
            ms_scaling_factor=ms_scaling_factor,
            osd_method=osd_method, osd_order=osd_order,
        )

    def decode_batch(self, detection_events: np.ndarray) -> np.ndarray:
        """Return predicted observable flips, shape (shots, num_observables)."""
        L = self.mats.L
        out = np.zeros((detection_events.shape[0], self.mats.num_observables),
                       dtype=np.uint8)
        for i, syn in enumerate(detection_events):
            corr = self.dec.decode(syn.astype(np.uint8))
            out[i] = (L @ corr) % 2
        return out


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact binomial confidence interval.  k = failures out of n trials."""
    from scipy.stats import beta

    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi
