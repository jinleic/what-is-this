"""Export the two EXP-066 open-bundle representatives for dist-m4ri preflight."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from qec_research.distance.sat_decide import css_logical_bases  # noqa: E402
SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "39x3.json"
OUT = ROOT / "scratch" / "exp067_m4ri_inputs"


def load_exp066():
    path = ROOT / "experiments" / "exp066_n234_frontier.py"
    spec = importlib.util.spec_from_file_location("exp066_for_m4ri", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matrix_sha256(matrix: np.ndarray) -> str:
    packed = np.packbits(np.asarray(matrix, dtype=np.uint8), axis=None)
    return hashlib.sha256(packed.tobytes()).hexdigest()


def main() -> int:
    exp066 = load_exp066()
    shard = json.loads(SHARD.read_text(encoding="utf-8"))
    residual = [
        record
        for record in shard["records"]
        if record.get("initial_verdict") == "solver_required"
    ]
    bundles = exp066.automorphism_bundles(residual)["records"]
    open_bundles = [
        bundle
        for bundle in bundles
        if all(residual[index]["verdict"] == "undecided" for index in bundle["record_indexes"])
    ]
    if [bundle["bundle_index"] for bundle in open_bundles] != [19, 22]:
        raise RuntimeError("the two EXP-066 open bundles changed")

    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for bundle in open_bundles:
        residual_index = int(bundle["representative_record"])
        record = residual[residual_index]
        hx, hz = exp066._problem_matrices(record)
        if hx.shape != (117, 234) or hz.shape != (117, 234):
            raise RuntimeError(f"unexpected BB matrix shapes: {hx.shape}, {hz.shape}")
        if np.any(hx @ hz.T % 2):
            raise RuntimeError("exported CSS matrices do not commute")

        lx, lz = css_logical_bases(hx, hz)
        if lx.shape != (8, 234) or lz.shape != (8, 234):
            raise RuntimeError(f"unexpected logical shapes: {lx.shape}, {lz.shape}")
        stem = f"bundle{bundle['bundle_index']}"
        swap = np.r_[np.arange(117, 234), np.arange(117)]
        hx_swapped, hz_swapped, lx_swapped = hx[:, swap], hz[:, swap], lx[:, swap]
        hx_path = OUT / f"{stem}_HX.mtx"
        hz_path = OUT / f"{stem}_HZ.mtx"
        lx_path = OUT / f"{stem}_LX.mtx"
        hx_swapped_path = OUT / f"{stem}_block_swapped_HX.mtx"
        hz_swapped_path = OUT / f"{stem}_block_swapped_HZ.mtx"
        lx_swapped_path = OUT / f"{stem}_block_swapped_LX.mtx"
        for path, matrix in (
            (hx_path, hx),
            (hz_path, hz),
            (lx_path, lx),
            (hx_swapped_path, hx_swapped),
            (hz_swapped_path, hz_swapped),
            (lx_swapped_path, lx_swapped),
        ):
            mmwrite(path, sparse.coo_matrix(matrix.astype(np.int8)), field="integer")
        manifest.append(
            {
                "bundle_index": int(bundle["bundle_index"]),
                "bundle_size": int(bundle["size"]),
                "residual_index": residual_index,
                "A": record["A"],
                "B": record["B"],
                "n": int(record["n"]),
                "k": int(record["k_parent"]),
                "known_upper_bound": int(record["witness_bound"]),
                "HX_sha256": matrix_sha256(hx),
                "HZ_sha256": matrix_sha256(hz),
                "HX_path": str(hx_path.relative_to(ROOT)),
                "HZ_path": str(hz_path.relative_to(ROOT)),
                "LX_sha256": matrix_sha256(lx),
                "LX_path": str(lx_path.relative_to(ROOT)),
                "block_swapped_HX_sha256": matrix_sha256(hx_swapped),
                "block_swapped_HZ_sha256": matrix_sha256(hz_swapped),
                "block_swapped_HX_path": str(hx_swapped_path.relative_to(ROOT)),
                "block_swapped_HZ_path": str(hz_swapped_path.relative_to(ROOT)),
                "block_swapped_LX_sha256": matrix_sha256(lx_swapped),
                "block_swapped_LX_path": str(lx_swapped_path.relative_to(ROOT)),
            }
        )

    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
