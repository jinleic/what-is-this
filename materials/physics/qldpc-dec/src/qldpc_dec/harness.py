"""BBCodeHarness: shared batch decode interface.

Contract (agreed with fss-bb sibling target, 2026-08-29):
    decode(dem: stim.DetectorErrorModel,
           shots: np.ndarray[bool, shape=(num_shots, num_dets)])
        -> np.ndarray[int, shape=(num_shots, num_obs)]

Returned predictions are 0/1 ints. All harness decoders (BP+OSD, beam
search, ensemble-NMS) implement this interface so the fss-bb erasure
adapter only needs one shape.
"""
from __future__ import annotations

from typing import Callable, Protocol

import numpy as np
import stim


class BatchDecoder(Protocol):
    """Any decoder exposing decode_batch(shots)->predictions (0/1 ints)."""

    def decode_batch(self, shots: np.ndarray) -> np.ndarray: ...


class BBCodeHarness:
    """Single-implementation decode surface shared across targets."""

    def __init__(self, decoder_factory: Callable[[stim.DetectorErrorModel], BatchDecoder]):
        self._decoder_factory = decoder_factory
        self._compiled: dict[int, BatchDecoder] = {}

    def _get_decoder(self, dem: stim.DetectorErrorModel) -> BatchDecoder:
        key = hash(str(dem))
        if key not in self._compiled:
            self._compiled[key] = self._decoder_factory(dem)
        return self._compiled[key]

    def decode(
        self, dem: stim.DetectorErrorModel, shots: np.ndarray
    ) -> np.ndarray:
        """Batch-decode bool shots (num_shots, num_dets) to int predictions (num_shots, num_obs)."""
        shots = np.ascontiguousarray(np.asarray(shots, dtype=np.bool_))
        decoder = self._get_decoder(dem)
        out = decoder.decode_batch(shots)
        out = np.asarray(out)
        out = out.reshape(len(shots), -1)
        return out.astype(np.int64)

    def decode_timed(
        self, dem: stim.DetectorErrorModel, shots: np.ndarray, component: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        """Per-shot timings: returns (predictions, time_seconds_per_shot)."""
        import time

        decoder = self._get_decoder(dem)
        n = len(shots)
        times = np.empty(n, dtype=np.float64)
        preds = np.empty((n, dem.num_observables), dtype=np.int64)
        for i in range(n):
            t0 = time.perf_counter()
            p = decoder.decode_batch(np.asarray(shots[i : i + 1], dtype=np.bool_))
            times[i] = time.perf_counter() - t0
            preds[i] = np.asarray(p).reshape(-1)
        return preds, times
