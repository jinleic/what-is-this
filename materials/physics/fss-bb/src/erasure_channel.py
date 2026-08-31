"""Quantum erasure channel sampling for CSS codes (code-level, per arXiv:2603.19062v3 Sec 2.2):

    1. each qubit erased independently w.p. p;
    2. erased qubits receive independent Pauli X and Z errors with Pr = 0.5 each
       (Y on an erased qubit w.p. 0.25);
    3. decoding uses erasure-informed priors: erased -> channel prob 0.5,
       non-erased -> 1e-10 (updated per shot);
    4. logical error iff either the X or the Z sector fails.

Why direct sampling rather than stim: the erasure mask differs per shot, so a
stim detector-error model would have to be rebuilt per mask; the paper's own
pipeline updates decoder posteriors "per shot before BP-OSD decoding". Direct
numpy sampling gives bit-identical channel semantics at ~µs/shot, and the
ldpc BpOsdDecoder accepts per-shot priors. Erasure DEMs (stim `ERANGE`/
`ERASE` instruction) are still exported per campaign point for interop with
the shared qldpc-dec harness (metadata/debug use, see export_erasure_dem).
"""

from __future__ import annotations

import numpy as np

NON_ERASED_PRIOR = 1e-10  # paper Sec 2.2 step 3
ERASED_PRIOR = 0.5


def sample_erasure_shots(
    rng: np.random.Generator,
    N: int,
    p: float,
    num_shots: int,
):
    """Sample erasure masks and error components.

    Returns:
        mask : (num_shots, N) bool — erased qubits
        ex   : (num_shots, N) bool — X component of the error
        ez   : (num_shots, N) bool — Z component of the error
    ex/ez are independent Bernoulli(0.5) on erased qubits and identically zero
    elsewhere (pure erasure channel).
    """
    mask = rng.random((num_shots, N)) < p
    erx = rng.random((num_shots, N)) < 0.5
    erz = rng.random((num_shots, N)) < 0.5
    ex = mask & erx
    ez = mask & erz
    return mask, ex, ez


def syndromes(Hx: np.ndarray, Hz: np.ndarray, ex: np.ndarray, ez: np.ndarray):
    """Syndromes for both sectors, vectorized over shots.

    X-checks (rows of Hx) detect Z errors: sx = (Hx @ ez^T)^T.
    Z-checks (rows of Hz) detect X errors: sz = (Hz @ ex^T)^T.
    Returns (sx, sz) uint8 arrays, shape (num_shots, n_xchecks), (num_shots, n_zchecks).
    """
    sx = (ez @ Hx.T % 2).astype(np.uint8)
    sz = (ex @ Hz.T % 2).astype(np.uint8)
    return sx, sz


def sector_failures(
    Lz: np.ndarray,
    Lx: np.ndarray,
    Hx: np.ndarray,
    Hz: np.ndarray,
    ex: np.ndarray,
    ez: np.ndarray,
    cx: np.ndarray,
    cz: np.ndarray,
):
    """Evaluate sector failures given corrections.

    cx/cz: (num_shots, N) bool corrections from the X/Z-sector decoders
    (cx = X-type correction applied against Z errors? Convention here:
      - the X-sector decoder runs on Hx with syndrome sx = Hx ez; its output
        `cz` is a Z-type vector; residual wz = ez + cz. Failure iff
        Hx wz == 0 but Lz wz != 0 (a Z logical was flipped).
      - the Z-sector decoder runs on Hz with syndrome sz = Hz ex; its output
        `cx` is an X-type vector; residual wx = ex + cx. Failure iff
        Hz wx == 0 but Lx wx != 0.
    Residuals with a violated check (Hx wz != 0) count as failures too — the
    decoder failed to return a coset-consistent answer.
    Returns (fail_x, fail_z) bool arrays, shape (num_shots,).
    """
    wz = ez ^ cz
    wx = ex ^ cx
    fail_z = ((wz @ Hx.T % 2).any(axis=1)) | ((wz @ Lz.T % 2).any(axis=1))
    fail_x = ((wx @ Hz.T % 2).any(axis=1)) | ((wx @ Lx.T % 2).any(axis=1))
    return fail_x, fail_z
