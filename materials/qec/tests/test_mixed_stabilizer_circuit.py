"""Validation of the generic mixed-stabilizer circuit compiler.

Gates enforced here:
  G5 (circuit semantics) -- the noiseless circuit produces zero detection
      events and zero observable flips, i.e. it really measures the intended
      stabilizers and preserves the logical subspace.
  G6 (fault propagation) -- Stim's undetectable-logical-error search returns a
      fault set whose size is a *circuit distance* upper bound; on toy codes we
      check it against the known code distance.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import stim

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec, build_memory_circuit, edge_colour_schedule,
    generator_supports, schedule_stats,
)
from qec_research.circuits.bicycle_schedule import pure_z_logical_basis  # noqa: E402
from qec_research.symplectic.core import StabilizerCode  # noqa: E402


def pauli_rows(strings: list[str]) -> np.ndarray:
    n = len(strings[0])
    H = np.zeros((len(strings), 2 * n), dtype=np.uint8)
    for i, s in enumerate(strings):
        for j, ch in enumerate(s):
            if ch in "XY":
                H[i, j] = 1
            if ch in "ZY":
                H[i, n + j] = 1
    return H


FIVE_QUBIT = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]          # [[5,1,3]] non-CSS
STEANE = ["IIIXXXX", "IXXIIXX", "XIXIXIX",
          "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"]                  # [[7,1,3]] CSS


def test_five_qubit_code_is_valid_and_non_css():
    H = pauli_rows(FIVE_QUBIT)
    code = StabilizerCode(H)
    v = code.validate()
    assert v["commutes_numpy"] and v["commutes_bitset"]
    assert v["k"] == 1
    assert not v["is_css"]


def test_edge_colouring_attains_max_degree():
    rng = np.random.default_rng(7)
    for _ in range(25):
        n = int(rng.integers(4, 25))
        r = int(rng.integers(2, 18))
        H = np.zeros((r, 2 * n), dtype=np.uint8)
        for i in range(r):
            k = int(rng.integers(1, min(6, n) + 1))
            cols = rng.choice(n, size=k, replace=False)
            for cq in cols:
                t = rng.integers(0, 3)
                if t == 0:
                    H[i, cq] = 1
                elif t == 1:
                    H[i, n + cq] = 1
                else:
                    H[i, cq] = H[i, n + cq] = 1
        sup = generator_supports(H)
        if not any(s.weight for s in sup):
            continue
        layers = edge_colour_schedule(sup, n)
        deg = np.zeros(n, dtype=int)
        for s in sup:
            for j in s.paulis:
                deg[j] += 1
        delta = max(max((s.weight for s in sup), default=0), int(deg.max()))
        assert len(layers) == delta, (len(layers), delta)
        # every edge scheduled exactly once
        got = sorted((c, q) for L in layers for c, q in L)
        want = sorted((s.index, j) for s in sup for j in s.paulis)
        assert got == want


@pytest.mark.parametrize("strings,name", [(STEANE, "steane"), (FIVE_QUBIT, "five")])
def test_noiseless_circuit_has_no_detection_events(strings, name):
    H = pauli_rows(strings)
    code = StabilizerCode(H)
    logs = pure_z_logical_basis(code)
    spec = CircuitSpec(H=H, observables=logs, rounds=4, p=0.0, basis="Z")
    circ, meta = build_memory_circuit(spec)
    assert meta["n_qubits_total"] == meta["n_data"] + meta["n_ancilla"]
    assert meta["num_observables"] >= 1, "need at least one Z-measurable logical"
    sampler = circ.compile_detector_sampler()
    dets, obs = sampler.sample(2000, separate_observables=True)
    assert not dets.any(), f"{name}: noiseless circuit fired {dets.sum()} detectors"
    assert not obs.any(), f"{name}: noiseless circuit flipped an observable"


@pytest.mark.parametrize("strings,name,dmax", [(STEANE, "steane", 3), (FIVE_QUBIT, "five", 3)])
def test_circuit_distance_upper_bound(strings, name, dmax):
    H = pauli_rows(strings)
    code = StabilizerCode(H)
    spec = CircuitSpec(H=H, observables=pure_z_logical_basis(code), rounds=3, p=1e-3, basis="Z")
    circ, _ = build_memory_circuit(spec)
    err = circ.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=4,
        dont_explore_edges_with_degree_above=4,
        dont_explore_edges_increasing_symptom_degree=False,
    )
    assert 1 <= len(err) <= dmax, f"{name}: circuit distance bound {len(err)}"


def test_detector_error_model_builds_and_is_nontrivial():
    H = pauli_rows(FIVE_QUBIT)
    code = StabilizerCode(H)
    spec = CircuitSpec(H=H, observables=pure_z_logical_basis(code), rounds=3, p=1e-3, basis="Z")
    circ, meta = build_memory_circuit(spec)
    dem = circ.detector_error_model(decompose_errors=False)
    assert dem.num_detectors == meta["num_detectors"]
    assert dem.num_errors > 0


# --------------------------------------------------------------------------
# Regression guards for the two falsified routes (notes/failed_routes.md)
# --------------------------------------------------------------------------
def test_edge_colouring_violates_parity_and_cpsat_does_not():
    """FR-002: a plain edge colouring is conflict-free but measures the wrong
    operators; the parity-enforcing scheduler must fix exactly that."""
    from qec_research.circuits.scheduling import (
        cpsat_schedule, depth_lower_bound, parity_defects, verify_schedule)

    H = pauli_rows(STEANE)
    sup = generator_supports(H)
    n = H.shape[1] // 2

    naive = edge_colour_schedule(sup, n)
    slot = {(c, q): li for li, layer in enumerate(naive) for c, q in layer}
    assert verify_schedule(sup, slot)["ancilla_conflict_free"]
    assert verify_schedule(sup, slot)["qubit_conflict_free"]
    assert parity_defects(sup, slot), "expected the colouring to break the parity rule"

    lb = depth_lower_bound(sup, n)
    good = None
    for T in range(lb, lb + 4):
        r = cpsat_schedule(sup, n, T=T, time_limit_s=60)
        if r.slot is not None and r.verification["valid"]:
            good = r
            break
    assert good is not None
    assert not parity_defects(sup, good.slot)


def test_pure_z_detector_not_same_sector():
    """FR-001: a pure-Z operator's nontriviality must be decided against the
    x-part of a logical basis element, never against a Z-type representative."""
    from qec_research.codes.bicycle import BRAVYI_BB, bb_stabilizer
    from qec_research.distance.sectors import min_pure_z_logical

    code = bb_stabilizer(BRAVYI_BB["[[72,12,6]]"])
    r = min_pure_z_logical(code, time_limit_s=120, workers=4)
    assert r.exact and r.weight == 6, (r.weight, r.exact)
    v = np.zeros(2 * code.n, dtype=np.uint8)
    for j in r.support:
        v[code.n + j] = 1
    assert code.is_logical(v), "returned witness must be a genuine logical"


def test_metadata_reports_the_real_schedule_depth_not_the_colouring_bound():
    """Regression: schedule_stats used to silently report the edge-colouring
    value Delta, which corresponds to an INVALID circuit.  For the Gross code
    Delta = 6 but the minimum valid depth is 7."""
    from qec_research.circuits.bicycle_schedule import (
        bb_supports_and_orbits, pure_z_logical_basis)
    from qec_research.circuits.mixed_stabilizer import schedule_stats
    from qec_research.circuits.scheduling import cpsat_schedule, slots_to_layers
    from qec_research.codes.bicycle import BRAVYI_BB

    code, sup, orb = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
    r = cpsat_schedule(sup, code.n, T=7, time_limit_s=120,
                       symmetry_orbits=orb, random_seed=0)
    assert r.slot is not None and r.verification["valid"]
    layers = slots_to_layers(r.slot, 7)

    # without a schedule, no depth may be invented
    bare = schedule_stats(code.H)
    assert bare["two_qubit_layers"] is None
    assert bare["depth_lower_bound"] == 6

    _, meta = build_memory_circuit(CircuitSpec(
        H=code.H, observables=pure_z_logical_basis(code), rounds=2, p=0.0,
        basis="Z", layers=layers))
    assert meta["two_qubit_layers"] == 7, meta["two_qubit_layers"]
    assert meta["depth_lower_bound"] == 6
    assert meta["total_two_qubit_gates"] == 864
    assert meta["n_qubits_total"] == 288
