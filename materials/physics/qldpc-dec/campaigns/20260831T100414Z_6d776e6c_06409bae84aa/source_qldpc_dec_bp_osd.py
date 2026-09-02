"""BP+OSD baseline decoder (ldpc package), gate A.

Configuration pinned in pre_statement.md (2026-08-29), matching the IonQ
beam-search paper Sec. III baseline ("30 min-sum BP iterations followed by
order-10 combination-sweep OSD") and their bp30+osd entry in
simulation_functions.py:

    BPOSD(dem, max_bp_iters=30, bp_method="min_sum", osd_order=10,
          osd_method="osd_cs")

The ldpc v2 BpOsdDecoder equivalent uses the merged DEM's heterogeneous
per-column priors:
    BpOsdDecoder(csrs, error_channel=p_vec, max_iter=30, bp_method="ms",
                 schedule="parallel", osd_method="osd_cs", osd_order=10)

Second arm recorded in the pre-statement (Relay-BP paper baseline):
max_iter=10000 (their "10,000 BP iterations with CS order 10"), runnable
via this module with an explicit parameter override.
"""
from __future__ import annotations

import numpy as np
import stim
from ldpc import BpOsdDecoder

from .dem_matrices import dem_to_matrices
from .harness import BBCodeHarness


def make_bp_osd_decoder(
    dem: stim.DetectorErrorModel,
    max_iter: int = 30,
    bp_method: str = "ms",
    schedule: str = "parallel",
    osd_method: str = "osd_cs",
    osd_order: int = 10,
    error_rate: float | None = None,
    **kwargs,
) -> BpOsdDecoder:
    """Build the pinned-reference BP+OSD decoder from a stim DEM.

    ``stimbposd.BPOSD`` passes the merged DEM's full prior vector through
    ``error_channel``. A scalar ``error_rate`` remains available only as an
    explicit experimental override.
    """
    H, _A, p_vec = dem_to_matrices(dem, merge=True)
    prior: dict[str, object]
    if error_rate is None:
        prior = {"error_channel": np.asarray(p_vec, dtype=float).tolist()}
    else:
        prior = {"error_rate": float(error_rate)}
    return BpOsdDecoder(
        H,
        max_iter=max_iter,
        bp_method=bp_method,
        schedule=schedule,
        osd_method=osd_method,
        osd_order=osd_order,
        **prior,
        **kwargs,
    )


class BpOsdBatchDecoder:
    """Adapter of ldpc's per-shot decode to the harness batch interface."""

    def __init__(self, dem: stim.DetectorErrorModel, **decoder_kwargs):
        self.dem = dem
        self.decoder = make_bp_osd_decoder(dem, **decoder_kwargs)
        # merged (ldpc) convention: A over the SAME columns the decoder saw
        _H, A, _p = dem_to_matrices(dem, merge=True)
        self.A = A.tocsr()

    def decode_batch(self, shots: np.ndarray) -> np.ndarray:
        n = len(shots)
        out = np.empty((n, self.A.shape[0]), dtype=np.int64)
        for i in range(n):
            corr = self.decoder.decode(shots[i].astype(np.uint8))
            out[i] = (self.A @ (corr.astype(np.int64) % 2)) % 2
        return out



def bp_osd_harness(error_rate: float | None = None, **decoder_kwargs) -> BBCodeHarness:
    """Harness bound to the pinned BP+OSD decoder."""
    if error_rate is not None:
        decoder_kwargs.setdefault("error_rate", error_rate)
    return BBCodeHarness(lambda dem: BpOsdBatchDecoder(dem, **decoder_kwargs))


if __name__ == "__main__":
    heterogeneous = stim.DetectorErrorModel(
        "error(0.1) D0\n"
        "error(0.2) D1"
    )
    prior_probe = make_bp_osd_decoder(heterogeneous, osd_order=0)
    np.testing.assert_array_equal(prior_probe.channel_probs, [0.1, 0.2])
    # smoke: tiny repetition-code DEM
    dem = stim.Circuit("""
        X_ERROR(0.001) 0
        M 0
        DETECTOR rec[-1]
        OBSERVABLE_INCLUDE(0) rec[-1]
    """).detector_error_model()
    h = bp_osd_harness()
    shots = np.zeros((4, dem.num_detectors), dtype=np.bool_)
    shots[1, 0] = True  # single flipped detector
    print(h.decode(dem, shots))
