"""Circuit-level syndrome extraction for CSS bivariate-bicycle codes.

The seven CNOT layers are the schedule used by the primary artifact for
arXiv:2308.07915.  A layer entry is ``(x_direction, z_direction)``; ``"idle"``
means that check family is being prepared or measured in that layer.  Direction
labels 0, 1, 2 are the three monomials of A (for X checks) / B^T (for Z
checks), and labels 3, 4, 5 are the three monomials of B / A^T.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

import numpy as np
import stim

from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np

__all__ = ["DEPTH7_SCHEDULE", "build_bb_memory_circuit"]

# Source: third_party/BivariateBicycleCodes/decoder_setup.py:45-51.
# Each pair is (sX[t], sZ[t]) for t=0,...,6.
DEPTH7_SCHEDULE: tuple[tuple[int | str, int | str], ...] = (
    ("idle", 3),
    (1, 5),
    (4, 0),
    (3, 1),
    (5, 2),
    (0, 4),
    (2, "idle"),
)


def _as_binary_matrix(matrix: np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(matrix, dtype=np.uint8)
    if result.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix")
    if np.any(result > 1):
        raise ValueError(f"{name} must be binary")
    return result


def _offset_counts(block: np.ndarray, ell: int, m: int) -> Counter[tuple[int, int]]:
    counts: Counter[tuple[int, int]] = Counter()
    rows, cols = np.nonzero(block)
    for row, col in zip(rows.tolist(), cols.tolist(), strict=True):
        ri, rj = divmod(row, m)
        ci, cj = divmod(col, m)
        counts[((ci - ri) % ell, (cj - rj) % m)] += 1
    return counts


def _infer_lattice_and_terms(hx: np.ndarray) -> tuple[int, int, list[tuple[int, int]], list[tuple[int, int]]]:
    """Recover the two torus dimensions and monomial shifts from ``[A B]``."""
    checks, n = hx.shape
    if n != 2 * checks:
        raise ValueError(f"expected HX shape (n/2, n), got {hx.shape}")
    candidates: list[tuple[int, int, list[tuple[int, int]], list[tuple[int, int]]]] = []
    for m in range(1, checks + 1):
        if checks % m:
            continue
        ell = checks // m
        a_counts = _offset_counts(hx[:, :checks], ell, m)
        b_counts = _offset_counts(hx[:, checks:], ell, m)
        if (len(a_counts) == 3 and len(b_counts) == 3
                and set(a_counts.values()) == {checks}
                and set(b_counts.values()) == {checks}):
            candidates.append((ell, m, list(a_counts), list(b_counts)))
    if len(candidates) != 1:
        shapes = [(ell, m) for ell, m, _, _ in candidates]
        raise ValueError(
            "could not uniquely infer the BB torus and its six monomial matchings "
            f"from HX; candidate (ell, m) values: {shapes}"
        )
    return candidates[0]


def _artifact_term_order(
    terms: Sequence[tuple[int, int]], *, polynomial: str
) -> list[tuple[int, int]]:
    """Order recovered shifts according to the A/B convention in the artifact."""
    if polynomial == "A":
        # A = x^a1 + y^a2 + y^a3.
        ordered = sorted(
            terms,
            key=lambda term: (
                0 if term[1] == 0 and term[0] != 0 else
                1 if term[0] == 0 and term[1] != 0 else
                2,
                term[0] + term[1],
                term,
            ),
        )
    elif polynomial == "B":
        # B = y^b1 + x^b2 + x^b3; B=1+x^b2+x^b3 also occurs.
        ordered = sorted(
            terms,
            key=lambda term: (
                0 if term[0] == 0 else 1,
                term[0] + term[1],
                term,
            ),
        )
    else:
        raise ValueError(f"unknown polynomial label {polynomial!r}")
    return ordered


def _shift(index: int, delta: tuple[int, int], ell: int, m: int) -> int:
    i, j = divmod(index, m)
    return ((i + delta[0]) % ell) * m + (j + delta[1]) % m


def _directional_neighbors(hx: np.ndarray, hz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return data-neighbor arrays indexed by direction then check row."""
    h = hx.shape[0]
    ell, m, a_terms, b_terms = _infer_lattice_and_terms(hx)
    a_terms = _artifact_term_order(a_terms, polynomial="A")
    b_terms = _artifact_term_order(b_terms, polynomial="B")

    expected_hz = np.hstack([hx[:, h:].T, hx[:, :h].T])
    if not np.array_equal(hz, expected_hz):
        raise ValueError("HZ is not [B^T A^T] for the supplied HX=[A B]")

    x_neighbors = np.empty((6, h), dtype=np.int64)
    z_neighbors = np.empty((6, h), dtype=np.int64)
    for row in range(h):
        for direction, delta in enumerate(a_terms):
            x_neighbors[direction, row] = _shift(row, delta, ell, m)
            z_neighbors[direction + 3, row] = h + _shift(
                row, (-delta[0], -delta[1]), ell, m
            )
        for offset, delta in enumerate(b_terms):
            direction = offset + 3
            x_neighbors[direction, row] = h + _shift(row, delta, ell, m)
            z_neighbors[offset, row] = _shift(
                row, (-delta[0], -delta[1]), ell, m
            )

    reconstructed_hx = np.zeros_like(hx)
    reconstructed_hz = np.zeros_like(hz)
    rows = np.arange(h)
    for direction in range(6):
        reconstructed_hx[rows, x_neighbors[direction]] ^= 1
        reconstructed_hz[rows, z_neighbors[direction]] ^= 1
    if not np.array_equal(reconstructed_hx, hx) or not np.array_equal(reconstructed_hz, hz):
        raise AssertionError("directional monomial decomposition did not reconstruct HX/HZ")
    return x_neighbors, z_neighbors


def _normalise_schedule(
    schedule: Sequence[Sequence[int | str]] | None,
) -> tuple[tuple[int | str, int | str], ...]:
    layers = DEPTH7_SCHEDULE if schedule is None else tuple(tuple(layer) for layer in schedule)
    if len(layers) != 7 or any(len(layer) != 2 for layer in layers):
        raise ValueError("a depth-7 schedule must contain seven (x_direction, z_direction) pairs")
    for family, position in (("X", 0), ("Z", 1)):
        values = [layer[position] for layer in layers]
        directions = [value for value in values if value != "idle"]
        if sorted(directions) != list(range(6)) or values.count("idle") != 1:
            raise ValueError(
                f"{family} schedule must use directions 0,...,5 exactly once and idle exactly once"
            )
    return layers


def _logical_z_basis(hx: np.ndarray, hz: np.ndarray) -> np.ndarray:
    """Representatives of ker(HX) / row(HZ), i.e. independent pure-Z logicals."""
    stabilizers, _ = rref_np(hz)
    current = stabilizers.copy()
    current_rank = current.shape[0]
    logicals: list[np.ndarray] = []
    for candidate in nullspace_np(hx):
        trial = np.vstack([current, candidate])
        trial_rank = rank_np(trial)
        if trial_rank > current_rank:
            logicals.append(candidate.copy())
            current = trial
            current_rank = trial_rank
    expected = hx.shape[1] - rank_np(hx) - rank_np(hz)
    if len(logicals) != expected:
        raise AssertionError(f"expected {expected} pure-Z logicals, found {len(logicals)}")
    return np.asarray(logicals, dtype=np.uint8)


def _record_target(measurement_index: int, measurement_count: int) -> stim.GateTarget:
    offset = measurement_index - measurement_count
    if offset >= 0:
        raise AssertionError("measurement record target must refer to the past")
    return stim.target_rec(offset)


def _append_noise(circuit: stim.Circuit, name: str, targets: Sequence[int], p: float) -> None:
    if targets:
        circuit.append(name, targets, p)


def build_bb_memory_circuit(
    HX: np.ndarray,
    HZ: np.ndarray,
    rounds: int,
    p: float,
    basis: str = "Z",
    schedule: Sequence[Sequence[int | str]] | None = None,
) -> stim.Circuit:
    """Build a noisy depth-7 BB syndrome-extraction memory circuit.

    Data qubits are ``0..n-1``.  The ``n/2`` X ancillas are ``n..3n/2-1``
    and the ``n/2`` Z ancillas are ``3n/2..2n-1``.  Noise of strength ``p``
    is applied after every two-qubit gate, at every scheduled data idle,
    after reset, and immediately before measurement.

    Only the Z-basis memory experiment is currently defined: data are prepared
    and measured in Z, terminal Z-check detectors close the time boundary, and
    the observables are an independent basis of pure-Z logical operators.
    """
    hx = _as_binary_matrix(HX, "HX")
    hz = _as_binary_matrix(HZ, "HZ")
    if hx.shape != hz.shape:
        raise ValueError(f"HX and HZ shapes differ: {hx.shape} versus {hz.shape}")
    h, n = hx.shape
    if h * 2 != n:
        raise ValueError(f"BB checks must have shape (n/2, n), got {hx.shape}")
    if np.any(hx.sum(axis=1) != 6) or np.any(hz.sum(axis=1) != 6):
        raise ValueError("the depth-7 BB schedule requires weight-6 checks")
    if np.any(hx.sum(axis=0) != 3) or np.any(hz.sum(axis=0) != 3):
        raise ValueError("the depth-7 BB schedule requires data degree three per check family")
    if np.any((hx.astype(np.int64) @ hz.T.astype(np.int64)) & 1):
        raise ValueError("HX and HZ do not commute")
    if not isinstance(rounds, int) or isinstance(rounds, bool) or rounds < 1:
        raise ValueError("rounds must be a positive integer")
    p = float(p)
    if not np.isfinite(p) or p < 0 or p > 1:
        raise ValueError("p must be a finite probability in [0, 1]")
    if basis.upper() != "Z":
        raise NotImplementedError("only the Z-basis BB memory experiment is implemented")

    layers = _normalise_schedule(schedule)
    x_neighbors, z_neighbors = _directional_neighbors(hx, hz)
    x_ancillas = list(range(n, n + h))
    z_ancillas = list(range(n + h, 2 * n))
    data_qubits = list(range(n))

    scheduled_gates: list[list[tuple[int, int]]] = []
    scheduled_idles: list[list[int]] = []
    for layer_index, (x_direction, z_direction) in enumerate(layers):
        gates: list[tuple[int, int]] = []
        if x_direction != "idle":
            gates.extend(
                (x_ancillas[row], int(x_neighbors[int(x_direction), row]))
                for row in range(h)
            )
        if z_direction != "idle":
            gates.extend(
                (int(z_neighbors[int(z_direction), row]), z_ancillas[row])
                for row in range(h)
            )
        used = [qubit for gate in gates for qubit in gate]
        assert len(used) == len(set(used)), (
            f"depth-7 layer {layer_index} is not disjoint; "
            f"{len(used) - len(set(used))} repeated qubit uses"
        )
        used_data = {qubit for qubit in used if qubit < n}
        scheduled_gates.append(gates)
        scheduled_idles.append([qubit for qubit in data_qubits if qubit not in used_data])

    circuit = stim.Circuit()
    # Product-state preparation and the initial Z-ancilla preparation precede
    # the artifact's first CNOT layer.  X ancillas are prepared in layer 0.
    circuit.append("R", data_qubits + z_ancillas)
    _append_noise(circuit, "X_ERROR", data_qubits + z_ancillas, p)
    circuit.append("TICK")

    measurement_count = 0
    previous_x: list[int] | None = None
    previous_z: list[int] | None = None

    for round_index in range(rounds):
        # Layer 0: prepare X ancillas while Z ancillas perform their first CNOT.
        circuit.append("RX", x_ancillas)
        _append_noise(circuit, "Z_ERROR", x_ancillas, p)
        gates = scheduled_gates[0]
        circuit.append("CX", [qubit for gate in gates for qubit in gate])
        _append_noise(circuit, "DEPOLARIZE2", [qubit for gate in gates for qubit in gate], p)
        _append_noise(circuit, "DEPOLARIZE1", scheduled_idles[0], p)
        circuit.append("TICK")

        # Layers 1 through 5 interleave the two check families.
        for layer_index in range(1, 6):
            gates = scheduled_gates[layer_index]
            circuit.append("CX", [qubit for gate in gates for qubit in gate])
            _append_noise(circuit, "DEPOLARIZE2", [qubit for gate in gates for qubit in gate], p)
            _append_noise(circuit, "DEPOLARIZE1", scheduled_idles[layer_index], p)
            circuit.append("TICK")

        # Layer 6: Z ancillas are measured while X ancillas perform their last CNOT.
        _append_noise(circuit, "X_ERROR", z_ancillas, p)
        circuit.append("MZ", z_ancillas)
        current_z = list(range(measurement_count, measurement_count + h))
        measurement_count += h
        gates = scheduled_gates[6]
        circuit.append("CX", [qubit for gate in gates for qubit in gate])
        _append_noise(circuit, "DEPOLARIZE2", [qubit for gate in gates for qubit in gate], p)
        _append_noise(circuit, "DEPOLARIZE1", scheduled_idles[6], p)
        for row, current in enumerate(current_z):
            targets = [_record_target(current, measurement_count)]
            if previous_z is not None:
                targets.append(_record_target(previous_z[row], measurement_count))
            circuit.append("DETECTOR", targets, [float(round_index), float(row), 0.0])
        previous_z = current_z
        circuit.append("TICK")

        # Layer 7 of the artifact timeline: all data idle, X measurement, Z reset.
        _append_noise(circuit, "DEPOLARIZE1", data_qubits, p)
        _append_noise(circuit, "Z_ERROR", x_ancillas, p)
        circuit.append("MX", x_ancillas)
        current_x = list(range(measurement_count, measurement_count + h))
        measurement_count += h
        circuit.append("R", z_ancillas)
        _append_noise(circuit, "X_ERROR", z_ancillas, p)
        if previous_x is not None:
            for row, current in enumerate(current_x):
                circuit.append(
                    "DETECTOR",
                    [
                        _record_target(current, measurement_count),
                        _record_target(previous_x[row], measurement_count),
                    ],
                    [float(round_index), float(row), 1.0],
                )
        previous_x = current_x
        circuit.append("TICK")

    # Destructive Z readout closes the Z-check boundary and reads pure-Z logicals.
    _append_noise(circuit, "X_ERROR", data_qubits, p)
    circuit.append("MZ", data_qubits)
    final_data = list(range(measurement_count, measurement_count + n))
    measurement_count += n
    assert previous_z is not None
    for row, support in enumerate(hz):
        targets = [
            _record_target(final_data[qubit], measurement_count)
            for qubit in np.flatnonzero(support)
        ]
        targets.append(_record_target(previous_z[row], measurement_count))
        circuit.append("DETECTOR", targets, [float(rounds), float(row), 0.0])

    for observable, logical_z in enumerate(_logical_z_basis(hx, hz)):
        circuit.append(
            "OBSERVABLE_INCLUDE",
            [
                _record_target(final_data[qubit], measurement_count)
                for qubit in np.flatnonzero(logical_z)
            ],
            observable,
        )

    assert circuit.num_qubits == 2 * n
    return circuit
