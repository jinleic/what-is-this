"""Sector-quotient dimension identity (J.3) and X-centralizer shrink bound (J.1).

Theorem G(ii)'s proof once leaned on the sentence "the X-sector of the quotient
is unchanged" — false as an argument (the pure-X centralizer shrinks, e.g.
phase2_58: 40 -> 12).  The correct support is the identity, true for every
stabilizer code,

    dim(Xcen / S_X) = dim(Zcen / S_Z) = k,

plus the rank identity used in the paper.  These tests lock both, computed
directly on catalogue rows with GF(2) machinery (no stored claims).

Numerics provenance: notes/theorem_j_xsector.md; the sector identity itself was
machine-checked by that session on 368/368 catalogue rows
(results/partial_runs/xsector/catalogue_scan.json) — here we recompute a sample
independently rather than trust the artifact.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
assert _SPEC and _SPEC.loader
E39 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(E39)
E27 = E39.E27

# One row per length, mixed classes, plus all four named rows from the J-note.
SAMPLE_ROWS = ["12_6_0193", "phase2_58", "phase2_60", "phase2_88", "15_6_0219",
               "9_6_0140", "12_6_0217", "12_6_0194"]


def _rows_by_label(catalogue):
    return {E27.catalogue_label(row, i): (i, row) for i, row in enumerate(catalogue)}


def _sector_dims(H: np.ndarray, n: int) -> tuple[int, int, int, int]:
    """(dim Xcen, dim S_X, dim Zcen, dim S_Z) for check matrix H = (HX | HZ)."""

    H = np.asarray(H, dtype=np.uint8)
    Xb, Zb = H[:, :n], H[:, n:]
    # Pure-Z centralizer: right nullspace of the X block (Hx z = 0).
    Zcen = nullspace_np(Xb).shape[0]
    # Pure-X centralizer: right nullspace of the Z block (Hz x = 0).
    Xcen = nullspace_np(Zb).shape[0]
    # S_X(Q): X-parts of row products whose Z-part vanishes.
    lam_z = nullspace_np(Zb.T)  # left nullspace of Zb
    Sx = rank_np((lam_z @ Xb) % 2)
    # S_Z(Q): Z-parts of row products whose X-part vanishes.
    lam_x = nullspace_np(Xb.T)  # left nullspace of Xb
    Sz = rank_np((lam_x @ Zb) % 2)
    return int(Xcen), int(Sx), int(Zcen), int(Sz)


@pytest.fixture(scope="module")
def catalogue():
    return E27.load_catalogue()


def test_sector_quotient_equals_k_all_sampled(catalogue) -> None:
    labels = _rows_by_label(catalogue)
    missing = [lab for lab in SAMPLE_ROWS if lab not in labels]
    assert not missing, f"sample rows missing from catalogue: {missing}"
    for lab in SAMPLE_ROWS:
        _, row = labels[lab]
        _, code = E27.pbb_code(row)
        H = np.asarray(code.H, dtype=np.uint8)
        n = code.H.shape[1] // 2
        Xcen, Sx, Zcen, Sz = _sector_dims(H, n)
        k = n - rank_np(H)
        assert k == Xcen - Sx == Zcen - Sz, (
            f"{lab}: sector identity broken: k={k}, "
            f"Xcen-S_X={Xcen - Sx}, Zcen-S_Z={Zcen - Sz}"
        )


def test_xcentralizer_shrinks_by_at_least_dressing(catalogue) -> None:
    # J.1: Xcen(P) - Xcen(Q) = rho_X >= dim(Delta_bar) = k_P - k_Q.
    labels = _rows_by_label(catalogue)
    for lab in SAMPLE_ROWS:
        _, row = labels[lab]
        _, code = E27.pbb_code(row)
        _, HX, HZ = E27.parent_matrices(row)
        HX = np.asarray(HX, dtype=np.uint8)
        HZ = np.asarray(HZ, dtype=np.uint8)
        n = code.H.shape[1] // 2
        # pure-X centralizers: parent Xcen = right nullspace of its Z block
        # [[0],[HZ]] (first block rows are all-X).
        H_z_parent = np.vstack(
            [np.zeros((HX.shape[0], n), dtype=np.uint8), HZ]
        )
        Xcen_P = nullspace_np(H_z_parent).shape[0]
        Xcen_Q = _sector_dims(np.asarray(code.H, dtype=np.uint8), n)[0]
        k_P = n - rank_np(HX) - rank_np(HZ)
        rho_x = int(Xcen_P) - int(Xcen_Q)
        dim_bar = k_P - code.k
        assert rho_x >= dim_bar >= 0, (
            f"{lab}: rho_X={rho_x} vs dim(Delta_bar)={dim_bar} (k_P={k_P}, k_Q={code.k})"
        )
