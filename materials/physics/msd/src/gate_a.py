"""Clifford/noise model for the zero-level CCZ Gate-A reproduction.

This module follows the executable pieces described in arXiv:2605.21867:
Fig. 4's [[8,3,2]] encoder, the double-check preparation, three logical-ZZ
teleportation parity measurements, and repeated X/Z syndrome checks.  As in
Sec. IV of that paper, T/T-dagger gates are replaced by identity so that the
circuit is simulable by Stim.  Surface-code patches and their decoder are not
expanded; three verification wires are reused as effective output wires, and
output X-frame bits are formed from the source-block X measurements.

The omitted surface-code boundary and the GHZ-to-parity scheduling scaffold
are [INFERENCE].  They are kept explicit in the metadata and README rather
than presented as an exact reconstruction of the paper's Appendix-A layout.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import stim

try:  # Support imports both as ``src.gate_a`` and with ``src`` on PYTHONPATH.
    from .code_832 import LX, LZ
    from .encoder import apply_encoder, input_basis
    from .layout import ALL_QUBITS, DATA, OUTPUT, SURGERY, SUPERDENSE_CHECKS, SYNDROME, VERIFY
except ImportError:  # pragma: no cover - exercised by direct CLI invocation.
    from code_832 import LX, LZ
    from encoder import apply_encoder, input_basis
    from layout import ALL_QUBITS, DATA, OUTPUT, SURGERY, SUPERDENSE_CHECKS, SYNDROME, VERIFY


@dataclass
class CircuitBuilder:
    """Append Stim operations with the paper's circuit-level Pauli noise.

    ``p`` is used for one- and two-qubit depolarization after unitary gates
    and idle wires.  Reset and measurement errors are represented by a Pauli
    flip immediately before the corresponding Stim operation, which is
    distributionally equivalent to the paper's post-operation flip for the
    measured bit.
    """

    p: float = 0.0

    def __post_init__(self) -> None:
        if not 0.0 <= self.p <= 1.0:
            raise ValueError(f"p must be in [0, 1], got {self.p!r}")
        self.circuit = stim.Circuit()
        self.measurement_count = 0

    def _idle_noise(self, active: Iterable[int]) -> None:
        if not self.p:
            return
        active_set = set(active)
        idle = [q for q in ALL_QUBITS if q not in active_set]
        if idle:
            self.circuit.append("DEPOLARIZE1", idle, self.p)

    def _tick(self, active: Iterable[int]) -> None:
        self._idle_noise(active)
        self.circuit.append("TICK")
    def idle_tick(self) -> None:
        """Advance one all-idle tick, including idle depolarization."""

        self._tick([])


    def reset_many(self, qubits: Sequence[int], *, basis: str) -> None:
        qubits = list(qubits)
        if not qubits:
            return
        if basis not in {"X", "Z"}:
            raise ValueError(f"reset basis must be X or Z, got {basis!r}")
        op = "RX" if basis == "X" else "R"
        error = "Z_ERROR" if basis == "X" else "X_ERROR"
        self.circuit.append(op, qubits)
        if self.p:
            self.circuit.append(error, qubits, self.p)
        self._tick(qubits)

    def h_layer(self, qubits: Sequence[int]) -> None:
        qubits = list(qubits)
        if not qubits:
            return
        self.circuit.append("H", qubits)
        if self.p:
            self.circuit.append("DEPOLARIZE1", qubits, self.p)
        self._tick(qubits)

    def cx_layer(self, pairs: Sequence[tuple[int, int]]) -> None:
        pairs = list(pairs)
        if not pairs:
            return
        flat = [q for pair in pairs for q in pair]
        if len(set(flat)) != len(flat):
            raise ValueError(f"CX layer has overlapping wires: {pairs!r}")
        self.circuit.append("CX", flat)
        if self.p:
            self.circuit.append("DEPOLARIZE2", flat, self.p)
        self._tick(flat)

    def cx_serial(self, pairs: Sequence[tuple[int, int]]) -> None:
        for pair in pairs:
            self.cx_layer([pair])

    def measure_many(self, qubits: Sequence[int], *, basis: str) -> list[int]:
        qubits = list(qubits)
        if not qubits:
            return []
        if basis not in {"X", "Z"}:
            raise ValueError(f"measurement basis must be X or Z, got {basis!r}")
        op = "MX" if basis == "X" else "M"
        error = "Z_ERROR" if basis == "X" else "X_ERROR"
        if self.p:
            self.circuit.append(error, qubits, self.p)
        self.circuit.append(op, qubits)
        records = list(range(self.measurement_count, self.measurement_count + len(qubits)))
        self.measurement_count += len(qubits)
        self._tick(qubits)
        return records

    def _record_target(self, record: int) -> stim.GateTarget:
        # Stim rec offsets are relative to the current end of the record stream.
        return stim.target_rec(record - self.measurement_count)

    def detector(self, record: int) -> None:
        self.circuit.append("DETECTOR", [self._record_target(record)])

    def observable(self, records: Sequence[int], index: int) -> None:
        self.circuit.append(
            "OBSERVABLE_INCLUDE",
            [self._record_target(record) for record in records],
            index,
        )


def _append_ghz_scaffold(builder: CircuitBuilder) -> None:
    """Prepare a four-wire GHZ state before the explicit parity checks.

    The published protocol describes one GHZ ancilla followed by four
    transversal CNOT layers.  Appendix A does not expose a machine-readable
    gate list, so this model records the GHZ preparation and then performs the
    four stabilizer parities with the same four wires after reset.  The reset
    boundary is deliberately visible rather than silently claiming an exact
    simultaneous-GHZ implementation.
    """

    builder.reset_many(SYNDROME, basis="Z")
    builder.h_layer([SYNDROME[0]])
    builder.cx_serial([(SYNDROME[0], q) for q in SYNDROME[1:]])


def _append_stabilizer_check(
    builder: CircuitBuilder, *, ancilla: int, basis: str, support: Iterable[int]
) -> int:
    """Measure one CSS stabilizer with a reset ancilla and return its record."""

    support = sorted(support)
    builder.reset_many([ancilla], basis=basis)
    if basis == "X":
        # |+> control -> data measures the product of data X operators.
        builder.cx_serial([(ancilla, q) for q in support])
    else:
        # data -> |0> target measures the product of data Z operators.
        builder.cx_serial([(q, ancilla) for q in support])
    return builder.measure_many([ancilla], basis=basis)[0]


def _append_syndrome_round(builder: CircuitBuilder) -> None:
    _append_ghz_scaffold(builder)
    for ancilla, (basis, support) in zip(SYNDROME, SUPERDENSE_CHECKS):
        record = _append_stabilizer_check(
            builder, ancilla=ancilla, basis=basis, support=support
        )
        # Every check is +1 for the ideal encoded |+++> state.
        builder.detector(record)


def build_gate_a_circuit(p: float = 0.0) -> stim.Circuit:
    """Build the Clifford-only Gate-A circuit at physical error rate ``p``.

    The returned circuit has three observables, one for each logical output.
    A shot is accepted when all detector bits are zero; an accepted shot is a
    logical failure when any observable bit is one.  The observable relation is
    the ideal teleportation X-frame relation
    ``output_X xor source_logical_X = 0`` for each of the three logicals.
    """

    builder = CircuitBuilder(float(p))

    # Non-fault-tolerant encoding followed by the double-check preparation.
    input_basis(builder, offset=0)
    apply_encoder(builder, offset=0)
    builder.reset_many(VERIFY, basis="X")
    builder.cx_layer(list(zip(DATA, VERIFY)))
    apply_encoder(builder, offset=8, inverse=True)
    verification_qubits = [8 + q for q in (1, 3, 4, 6)]
    verification_records = builder.measure_many(verification_qubits, basis="X")
    for record in verification_records:
        builder.detector(record)
    builder.reset_many(verification_qubits, basis="X")
    apply_encoder(builder, offset=8)
    builder.cx_layer(list(zip(DATA, VERIFY)))

    # T/T-dagger layer is intentionally identity (Sec. IV's Stim model).
    builder.idle_tick()


    # One superdense-style check and two further syndrome rounds.  This is the
    # explicit postselection boundary used by the estimator.
    for _ in range(3):
        _append_syndrome_round(builder)

    # Logical ZZ surgery.  The two reusable ancillas are parity meters; the
    # three verification wires in OUTPUT are the effective d=3 surface-code
    # representatives in this 22-wire model.  These three surgery outcomes are
    # random in the ideal circuit and are intentionally not postselected.
    for k, source_pair in enumerate(LZ):
        source_a, source_b = sorted(source_pair)
        ancilla = SURGERY[k % len(SURGERY)]
        builder.reset_many([ancilla], basis="Z")
        builder.cx_serial(
            [
                (source_a, ancilla),
                (source_b, ancilla),
                (OUTPUT[k], ancilla),
            ]
        )
        builder.measure_many([ancilla], basis="Z")

    # Directly measure the source block in X, as in the paper.  The source
    # logical-X parities provide the classical teleportation frame for the
    # effective output X measurements.
    source_x_records = builder.measure_many(DATA, basis="X")
    output_x_records = builder.measure_many(OUTPUT, basis="X")
    for k, logical_x_support in enumerate(LX):
        records = [output_x_records[k]] + [source_x_records[q] for q in sorted(logical_x_support)]
        builder.observable(records, k)

    return builder.circuit


def circuit_metadata(p: float = 0.0) -> dict[str, object]:
    """Return machine-readable metadata used in campaign artifacts."""

    circuit = build_gate_a_circuit(p)
    return {
        "physical_wires": len(ALL_QUBITS),
        "effective_output_wires": OUTPUT,
        "physical_error_rate": float(p),
        "noise_model": {
            "single_qubit_gate_and_idle": "DEPOLARIZE1(p)",
            "two_qubit_gate": "DEPOLARIZE2(p)",
            "reset": "basis-flip(p)",
            "measurement": "basis-flip(p)",
        },
        "t_layer": "identity (Clifford-only Stim surrogate)",
        "decoder": "detector postselection + effective X-frame parity; not paper decoder",
        "num_qubits": circuit.num_qubits,
        "num_ticks": circuit.num_ticks,
        "num_detectors": circuit.num_detectors,
        "num_observables": circuit.num_observables,
    }
