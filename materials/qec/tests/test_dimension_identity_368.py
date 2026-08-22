"""Full-catalogue machine check of the PBB dimension identity.

This is the guard named by the attribution block in ``proofs/pbb_structure.md``
(Proposition (ii)): the *proved* identity ``k_P - k_Q = dim(Delta-bar)`` is
re-derived from first principles and machine-verified on **every one of the
368 catalogue rows** -- not a sample.  ``tests/test_pbb_survival.py`` checks
the same identity on 4 rows only; this file covers all 368.

Setup.  A catalogue row carries the cyclic-shift polynomials of a parent BB
code ``H_X = [A B]``, ``H_Z = [B^T A^T]`` and of its perturbed-BB child

    H_Q = [ A  B | C    D   ]
          [ 0  0 | B^T  A^T ].

The dressing space is ``Delta = { lambda [C D] : lambda [A B] = 0 }`` and
``dim(Delta-bar) = dim (Delta + S_Z)/S_Z`` is computed EXACTLY as

    rank(H_Z stacked Delta) - rank(H_Z),

with ``Delta`` built directly from the left nullspace of ``H_X`` -- no
shortcut through ``analyse_pbb`` or the survival module.  Per row:

(a) ``k_P - k_Q == dim(Delta-bar)``, where ``k_Q`` comes from
    ``StabilizerCode.k`` (bitset rank path, independent of ``rank_np``);
(b) ``rank(H_Q) == rank(H_X) + rank(H_Z stacked Delta)``.

Coverage: all 368 rows, including all 50 rows with ``n = 360``.  A full pass
measures ~8 s wall-clock on the reference machine (2026-08-18), far under the
120 s budget, so NO ``n <= 144`` restriction is applied.
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

# Import E27 the same way experiments/exp039_nogo_module.py does: the
# catalogue loader and the parent/PBB builders live there.
_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E27
_SPEC.loader.exec_module(E27)

EXPECTED_ROWS = E27.EXPECTED_CATALOGUE_ROWS  # 368
EXPECTED_DELTA_POSITIVE_ROWS = E27.EXPECTED_DELTA_POSITIVE_ROWS  # 155


def _row_check(row: dict) -> dict:
    """Rebuild parent and PBB from catalogue polynomials; measure both sides."""

    _, HX, HZ = E27.parent_matrices(row)
    _, code = E27.pbb_code(row)
    n = HX.shape[1]
    # Z-part of the dressed first block: H_Q[:dim, n:] == [C D].
    CD = code.H[: HX.shape[0], n:]

    # Delta = { lambda [C D] : lambda [A B] = 0 }, from the left kernel of HX.
    left = nullspace_np(np.ascontiguousarray(HX.T))
    if left.size:
        delta = (left.astype(np.uint8) @ CD.astype(np.uint8)) % 2
    else:
        delta = np.zeros((0, n), dtype=np.uint8)

    rank_hx = rank_np(HX)
    rank_hz = rank_np(HZ)
    rank_hz_delta = rank_np(np.vstack([HZ, delta]))
    dim_delta_bar = rank_hz_delta - rank_hz
    k_parent = n - rank_hx - rank_hz
    return {
        "n": n,
        "k_parent": k_parent,
        "k_pbb": code.k,  # bitset rank path -- independent of rank_np
        "dim_delta_bar": dim_delta_bar,
        "rank_hq": rank_np(code.H),
        "rank_hq_predicted": rank_hx + rank_hz_delta,
    }


@pytest.fixture(scope="module")
def catalogue() -> list[dict]:
    rows = E27.load_catalogue()
    assert len(rows) == EXPECTED_ROWS
    return rows


@pytest.fixture(scope="module")
def checks(catalogue) -> list[dict]:
    """Single pass over the full catalogue; both tests below consume it."""

    return [_row_check(row) for row in catalogue]


@pytest.mark.parametrize("index", range(EXPECTED_ROWS))
def test_dimension_identity_on_row(checks, index):
    """(a) k_parent - k_pbb == dim(Delta-bar), exactly, zero-mismatch."""

    c = checks[index]
    assert (
        c["k_parent"] - c["k_pbb"] == c["dim_delta_bar"]
    ), f"row {index} (n={c['n']}): k_P - k_Q mismatch"


@pytest.mark.parametrize("index", range(EXPECTED_ROWS))
def test_rank_decomposition_on_row(checks, index):
    """(b) rank(H_Q) == rank(H_X) + rank(H_Z stacked Delta), exactly."""

    c = checks[index]
    assert (
        c["rank_hq"] == c["rank_hq_predicted"]
    ), f"row {index} (n={c['n']}): rank(H_Q) decomposition mismatch"


def test_full_368_coverage_and_delta_split(checks):
    """Every catalogue row was exercised; the delta>0 split matches EXP-027."""

    assert len(checks) == EXPECTED_ROWS
    assert sum(c["dim_delta_bar"] > 0 for c in checks) == EXPECTED_DELTA_POSITIVE_ROWS
