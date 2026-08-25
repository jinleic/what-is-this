"""Exact cross-lattice CRT transport checks for EXP-059."""
from __future__ import annotations

import json
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp059", ROOT / "experiments" / "exp059_coprime_transport.py"
)
E59 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E59
assert _SPEC.loader is not None
_SPEC.loader.exec_module(E59)


@pytest.mark.parametrize(
    "source,target",
    [
        ((15, 7), (21, 5)),
        ((15, 7), (35, 3)),
        ((21, 5), (35, 3)),
    ],
)
def test_crt_coordinate_transport_is_bijective(source, target) -> None:
    mapping = E59.coordinate_transport(*source, *target)
    assert mapping.shape == (105,)
    assert np.array_equal(np.sort(mapping), np.arange(105))
    inverse = E59.coordinate_transport(*target, *source)
    assert np.array_equal(inverse[mapping], np.arange(105))


@pytest.mark.parametrize(
    "source,target",
    [
        ((15, 7), (21, 5)),
        ((15, 7), (35, 3)),
        ((21, 5), (35, 3)),
    ],
)
def test_bb_matrices_transport_exactly(source, target) -> None:
    exponent_supports = [
        ([0, 1, 37], [0, 14, 73]),
        ([0, 22, 91], [0, 35, 68]),
    ]
    for A_exponents, B_exponents in exponent_supports:
        A = E59.exponents_to_terms(A_exponents, *source)
        B = E59.exponents_to_terms(B_exponents, *source)
        proof = E59.transport_certificate(*source, *target, A, B)
        assert proof["valid"] is True
        assert proof["HX_exact"] and proof["HZ_exact"]
        assert proof["source_k"] == proof["target_k"]
        assert proof["source_n"] == proof["target_n"] == 210


def test_n210_coprime_family_has_one_exact_search_space() -> None:
    proof = E59.coprime_family_certificate(105)
    assert proof["factorizations"] == [[15, 7], [21, 5], [35, 3], [105, 1]]
    assert proof["screen_lattices"] == [[15, 7], [21, 5], [35, 3]]
    assert proof["all_pairwise_transports_valid"] is True
    assert proof["representative"] == [15, 7]

def test_transported_n210_screen_shards_preserve_every_verdict() -> None:
    source_path = (
        ROOT / "results" / "partial_runs" / "exp055_screen" / "15x7.json"
    )
    source = json.loads(source_path.read_text(encoding="utf-8"))
    for name, lattice in (("21x5.json", [21, 5]), ("35x3.json", [35, 3])):
        target = json.loads(
            (source_path.parent / name).read_text(encoding="utf-8")
        )
        assert [target["ell"], target["m"]] == lattice
        assert target["protocol"] == source["protocol"]
        assert target["candidates_after_symmetry"] == source[
            "candidates_after_symmetry"
        ]
        assert target["orbit_total"] == source["orbit_total"]
        assert target["verdicts"] == source["verdicts"]
        assert target["transport"]["valid"] is True
        assert target["transport"]["all_physical_witnesses_valid"] is True
        assert target["transport"]["records_transported"] == len(source["records"])

def test_transport_missing_proof_support_fails_closed(tmp_path: Path) -> None:
    source_path = (
        ROOT / "results" / "partial_runs" / "exp055_screen" / "15x7.json"
    )
    source = json.loads(source_path.read_text(encoding="utf-8"))

    missing_logical = json.loads(json.dumps(source))
    record = next(
        row for row in missing_logical["records"]
        if row["verdict"] in {
            "dominated_by_witness", "dominated_by_cdcl_witness"
        }
    )
    record.pop("witness_support", None)
    logical_path = tmp_path / "missing-logical.json"
    logical_path.write_text(json.dumps(missing_logical))
    with pytest.raises(RuntimeError, match="witness failed"):
        E59.transport_screen_shard(
            logical_path, 21, 5, tmp_path / "logical-target.json"
        )

    missing_ceiling = json.loads(json.dumps(source))
    record = next(
        row for row in missing_ceiling["records"]
        if row.get("ceiling_witness_support")
    )
    record["verdict"] = "dominated_by_ceiling"
    record.pop("ceiling_witness_support", None)
    ceiling_path = tmp_path / "missing-ceiling.json"
    ceiling_path.write_text(json.dumps(missing_ceiling))
    with pytest.raises(RuntimeError, match="witness failed"):
        E59.transport_screen_shard(
            ceiling_path, 21, 5, tmp_path / "ceiling-target.json"
        )
