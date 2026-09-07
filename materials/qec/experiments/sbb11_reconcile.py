"""Exact Gross graph reconciliation; endpoint annotations are not clock times.

This source gate does not construct a memory circuit or compute its distance.
The canonical builder and both existing GF(2) implementations are snapshotted
as inputs to the registered campaign, rather than reimplemented here.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import time
from pathlib import Path

import numpy as np


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def point(u: int, v: int, offset: tuple[int, int]) -> tuple[int, int]:
    return ((offset[0] - 2 * u) % 24, (offset[1] - 2 * v) % 12)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    inputs = args.run / "inputs"
    builder = load_module("canonical_bb", inputs / "bb_codes.py")
    gf2 = load_module("existing_gf2", inputs / "gf2_linalg.py")
    source = json.loads((inputs / "legacy_transcription.json").read_text())
    code = builder.bb_code(12, 6)
    hx, hz = code["Hx"], code["Hz"]
    data_coordinates = [point(u, v, offset) for offset in ((1, 4), (4, 1)) for u in range(12) for v in range(6)]
    qubit_at = {p: q for q, p in enumerate(data_coordinates)}
    assert len(qubit_at) == 144
    assert set(qubit_at) == {(x, y) for y in range(12) for x in range(24) if (x + y) % 2}
    coordinates = {
        "X": [point(u, v, (1, 1)) for u in range(12) for v in range(6)],
        "Z": [point(u, v, (4, 4)) for u in range(12) for v in range(6)],
    }
    figure_checks = {}
    spoke_targets = {}
    colour_classes = {}
    for kind, style in (("X", "filled"), ("Z", "hollow")):
        geometry = source["geometry"][style]
        offsets = {"N": (0, -1), "E": (1, 0), "S": (0, 1), "W": (-1, 0), "NE": geometry["NE"], "SW": geometry["SW"]}
        matrix = np.zeros((72, 144), dtype=np.uint8)
        rows = []
        colours = []
        for i, (x, y) in enumerate(coordinates[kind]):
            assert (x % 2, y % 2) == (geometry["x_parity"], geometry["y_parity"])
            targets = {spoke: qubit_at[((x + dx) % 24, (y + dy) % 12)] for spoke, (dx, dy) in offsets.items()}
            assert len(set(targets.values())) == 6
            matrix[i, list(targets.values())] = 1
            rows.append(targets)
            colours.append((x // 2 + y // 2) % 3)
        figure_checks[kind] = matrix
        spoke_targets[kind] = rows
        colour_classes[kind] = colours
    rank_results = {}
    for kind, matrix in (("X", hx), ("Z", hz)):
        rank_results[kind] = {"canonical_builder": builder.gf2_rank(matrix), "numpy": gf2.rank_np(matrix), "bitset": gf2.rank_bitset(gf2.rows_to_bitsets(matrix), 144)}
        assert set(rank_results[kind].values()) == {66}
        assert matrix.shape == (72, 144)
        assert np.all(matrix.sum(axis=1) == 6)
        assert np.all(matrix.sum(axis=0) == 3)
        assert np.array_equal(matrix, figure_checks[kind])
        for colour in range(3):
            rows = [i for i, c in enumerate(colour_classes[kind]) if c == colour]
            assert len(rows) == 24
            assert np.all(matrix[rows].sum(axis=0) == 1)
    assert not gf2.matmul(hx, hz.T).any()
    assert 144 - rank_results["X"]["bitset"] - rank_results["Z"]["bitset"] == 12
    lx, lz = code["Lx"], code["Lz"]
    assert lx.shape == lz.shape == (12, 144)
    assert not gf2.matmul(hz, lx.T).any()
    assert not gf2.matmul(hx, lz.T).any()
    pairing = gf2.matmul(lx, lz.T)
    assert gf2.rank_np(pairing) == gf2.rank_bitset(gf2.rows_to_bitsets(pairing), 12) == 12
    usage = resource.getrusage(resource.RUSAGE_SELF)
    report = {
        "status": "CANONICAL_GRAPH_EQUIVALENCE_VERIFIED",
        "builder": "inputs/bb_codes.py:bb_code(12,6)",
        "polynomials": {"A": "x^3+y+y^2", "B": "y^3+x+x^2"},
        "shapes": {"Hx": list(hx.shape), "Hz": list(hz.shape)},
        "ranks": rank_results,
        "n": 144,
        "k": 12,
        "nominal_published_code_distance": 12,
        "distance12_recomputed_here": False,
        "CSS_commutation": True,
        "row_weight": 6,
        "column_weight": 3,
        "coordinate_map": {"index": "6u+v, modulo(24,12)", "Xcheck": "(1-2u,1-2v)", "Zcheck": "(4-2u,4-2v)", "Ldata": "(1-2u,4-2v)", "Rdata": "(4-2u,1-2v)"},
        "full_entrywise_figure_equivalence": {"X": True, "Z": True},
        "three_colour_cover": "Each colour has24 checks and touches each data qubit exactly once, separately in X/Z",
        "logical_basis_shape": [12, 144],
        "logical_pairing_rank": 12,
        "monomial_spokes": {"X": {"L_x3": "SW", "L_y": "S", "L_y2": "N", "R_y3": "NE", "R_x": "E", "R_x2": "W"}, "Z": {"L_y_minus3": "SW", "L_x_minus1": "W", "L_x_minus2": "E", "R_x_minus3": "NE", "R_y_minus1": "N", "R_y_minus2": "S"}},
        "endpoint_annotation_semantics": "Not assumed to be absolute times; no clock phase tests performed",
        "circuit_admission": "PENDING exact schedule and direct construction boundary reconciliation",
        "physical_fault_model_derived": False,
        "distance_assay_calls": 0,
        "usage": {"wall_seconds": time.monotonic() - started, "cpu_seconds": usage.ru_utime + usage.ru_stime, "max_rss_bytes_darwin": usage.ru_maxrss},
    }
    np.savez_compressed(args.run / "canonical_gross.npz", Hx=hx, Hz=hz, Lx=lx, Lz=lz, figure_Hx=figure_checks["X"], figure_Hz=figure_checks["Z"], logical_pairing=pairing)
    (args.run / "coordinate_map.json").write_text(json.dumps({"data": data_coordinates, "checks": coordinates, "spoke_targets": spoke_targets, "colour_classes": colour_classes}, indent=2) + "\n")
    (args.run / "canonical_reconciliation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
