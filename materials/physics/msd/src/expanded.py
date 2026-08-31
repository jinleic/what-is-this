"""Distance-carrying output patch for the Gate-A surrogate.

The committed Gate-A circuit routes all three observables through distance-1
wires: observable k reads one output MX record plus four source MX records,
so a single measurement or frame flip moves an observable at weight one and
the whole surrogate shows p_L ~ c*p (owner diagnostic: log-log slope 0.992,
coefficient ~167) instead of the paper's reported c ~ 300 * p^2 suppression.

This module keeps ``gate_a.py`` untouched.  It reproduces the
``build_gate_a_circuit`` preamble call-for-call on a builder subclass whose
idle noise also covers the fresh patch wires, then replaces the single-shot
output/source X readouts with a repetition-protected readout:

1. An X-frame repetition codeword per output: the output wire's X frame is
   copied (serial CX, shared control) into two fresh |0> wires, giving a
   3-wire codeword {OUTPUT[k], a, b} whose 3-wire X parity is the frame.
   Per-wire X values are NOT deterministic (the ideal readout is a random
   triplet like (0,1,1)); only the parity is the frame.  The deterministic
   codeword stabilizers are the Z-type pair parities Z_oZ_a and Z_aZ_b.
2. ``output_distance`` rounds of those Z pair checks: round 1 appends
   deterministic (+1) detectors, later rounds consecutive-round comparison
   detectors.  These detect X-type errors; X errors are harmless for an
   X-frame observable, so this layer is scaffolding, not the protection.
3. ``output_distance`` rounds of codeword MX with per-wire consecutive-round
   comparison detectors.  A Z flip (including the Z_ERROR(p) immediately
   before any MX) hits a subset of rounds and fires the comparisons around
   it.  Round 1 has no anchor (the frame is random), so Z flips before
   round 1 remain undetected by design; those live in the narrow
   injection-to-first-MX window.
4. ``output_distance`` rounds of source-block X readout with per-wire
   consecutive-round comparisons; the last round supplies the teleportation
   frame.  Z flips on DATA anywhere between the syndrome rounds and the
   first source round flip every round consistently (no anchor exists
   before the first X measurement of a wire) — this residual weight-1 path
   is inherited from the base model's frame structure, not from the patch.
5. Observable k = XOR(final codeword MX records) XOR XOR(final source MX
   records over the logical-X support) — the same teleportation X-frame
   relation as Gate-A, with the output side spread over 3 records.

Known irreducible weight-1 path ([INFERENCE], stated up front): OUTPUT wires
carry a random, never-anchorable X frame from their re-encoding through the
surgery block to the patch injection point (~100 ticks).  A single Z there
is copied consistently to all three codeword wires and is invisible to every
detector.  Injecting earlier (before surgery) would couple the surgery
ancillas into the codeword (CX(o, anc) maps Z_anc -> Z_o Z_anc), creating
new weight-1 paths, so the patch is deliberately injected post-surgery.
The smoke run therefore measures how much of the slope-1 coefficient is
readout-bound (fixable) versus frame-history-bound (irreducible here).

This is an [INFERENCE] 1D repetition-style patch, not the paper's 2D
surface-code patch: it demonstrates whether detector-protected redundancy on
the readout is sufficient, or a 2D patch boundary is required.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import stim

try:  # Support imports both as ``src.expanded`` and with ``src`` on PYTHONPATH.
    from .code_832 import LX, LZ
    from .encoder import apply_encoder, input_basis
    from .gate_a import CircuitBuilder, _append_stabilizer_check, _append_syndrome_round
    from .layout import ALL_QUBITS, DATA, OUTPUT, SURGERY, VERIFY
except ImportError:  # pragma: no cover - exercised by direct script invocation.
    from code_832 import LX, LZ
    from encoder import apply_encoder, input_basis
    from gate_a import CircuitBuilder, _append_stabilizer_check, _append_syndrome_round
    from layout import ALL_QUBITS, DATA, OUTPUT, SURGERY, VERIFY


N_OUTPUTS = len(OUTPUT)


@dataclass
class ExpandedBuilder(CircuitBuilder):
    """CircuitBuilder whose idle depolarization also covers the patch wires.

    The base class hardcodes idle noise to ``ALL_QUBITS`` (wires 0..21); the
    patch wires appended after the Gate-A preamble must idle-noise too or
    they would silently escape the circuit-level noise model.
    """

    extra_wires: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        super().__post_init__()

    def _idle_noise(self, active: Iterable[int]) -> None:
        if not self.p:
            return
        active_set = set(active)
        pool = list(ALL_QUBITS) + list(self.extra_wires)
        idle = [q for q in pool if q not in active_set]
        if idle:
            self.circuit.append("DEPOLARIZE1", idle, self.p)

    def comparison_detector(self, record_a: int, record_b: int) -> None:
        """Deterministic-parity detector over two measurement records."""

        self.circuit.append(
            "DETECTOR",
            [self._record_target(record_a), self._record_target(record_b)],
        )


def _temporally_compared_readout(
    builder: ExpandedBuilder, qubits: list[int], rounds: int
) -> list[list[int]]:
    """Read ``qubits`` in X for ``rounds`` rounds with temporal comparisons.

    Records are absolute; rounds are outer, wires inner.  Per wire, every
    consecutive round pair gets a comparison detector (a wire's X value is
    time-invariant absent errors).  Round 1 has no temporal anchor because
    the X values are genuinely random.
    """

    all_rounds: list[list[int]] = []
    for round_index in range(rounds):
        records = builder.measure_many(qubits, basis="X")
        if round_index > 0:
            previous = all_rounds[round_index - 1]
            for position in range(len(qubits)):
                builder.comparison_detector(previous[position], records[position])
        all_rounds.append(records)
    return all_rounds


def build_expanded_circuit(
    p: float = 0.0,
    *,
    output_distance: int = 3,
    debug_direct_output_mx: bool = False,
) -> stim.Circuit | tuple[stim.Circuit, dict[str, object]]:
    """Build the patched Gate-A circuit at physical error rate ``p``.

    The preamble (encoding, double-check verification, three syndrome rounds,
    ZZ surgery) uses the identical builder-call sequence as
    ``build_gate_a_circuit``, with the patch injected after the surgery block.
    With ``debug_direct_output_mx`` three extra plain ``MX(OUTPUT)`` records
    are appended after all detector/observable targets and the return value
    becomes ``(circuit, metadata)`` with absolute record indices; otherwise
    only the circuit is returned.
    """

    if output_distance < 3 or output_distance % 2 == 0:
        raise ValueError(
            "output_distance must be an odd integer >= 3 (the X-frame parity "
            "needs an odd number of codeword records), got "
            f"{output_distance!r}"
        )
    p = float(p)

    # Wire plan: per output, one fresh codeword wire pair (a, b) after the
    # representative OUTPUT[k], then two reusable Z-check ancillas.
    extras_per_output = output_distance - 1
    cursor = len(ALL_QUBITS)
    patch_data: dict[int, list[int]] = {}
    for k in range(N_OUTPUTS):
        patch_data[k] = list(range(cursor, cursor + extras_per_output))
        cursor += extras_per_output
    patch_ancillas: dict[int, list[int]] = {}
    for k in range(N_OUTPUTS):
        patch_ancillas[k] = [cursor, cursor + 1]
        cursor += 2
    extra_wires = tuple(range(len(ALL_QUBITS), cursor))
    builder = ExpandedBuilder(p, extra_wires=extra_wires)

    # ---- Gate-A preamble: identical calls to build_gate_a_circuit ----
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

    for _ in range(3):
        _append_syndrome_round(builder)

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

    # ---- distance-carrying output patch ----
    # Fresh |0> codeword wires, then copy the output X frame (serial CX
    # because the control OUTPUT[k] is shared).
    for k in range(N_OUTPUTS):
        extras = patch_data[k]
        builder.reset_many(extras, basis="Z")
        builder.cx_serial([(OUTPUT[k], wire) for wire in extras])

    # output_distance rounds of codeword Z pair checks: Z_oZ_a and Z_aZ_b.
    # Round 1 is deterministic (+1); later rounds compare consecutive rounds.
    check_records: list[dict[int, tuple[int, int]]] = []
    for round_index in range(output_distance):
        round_records: dict[int, tuple[int, int]] = {}
        for k in range(N_OUTPUTS):
            anc_pair, anc_codeword = patch_ancillas[k]
            rep = OUTPUT[k]
            first_extra = patch_data[k][0]
            second_extra = patch_data[k][1]
            rec_pair = _append_stabilizer_check(
                builder, ancilla=anc_pair, basis="Z", support={rep, first_extra}
            )
            rec_codeword = _append_stabilizer_check(
                builder,
                ancilla=anc_codeword,
                basis="Z",
                support={first_extra, second_extra},
            )
            round_records[k] = (rec_pair, rec_codeword)
            if round_index == 0:
                builder.detector(rec_pair)
                builder.detector(rec_codeword)
            else:
                prev_pair, prev_codeword = check_records[round_index - 1][k]
                builder.comparison_detector(prev_pair, rec_pair)
                builder.comparison_detector(prev_codeword, rec_codeword)
        check_records.append(round_records)

    # output_distance rounds of codeword MX: temporal comparisons catch any
    # Z flip that hits a subset of rounds; the last round is the observable.
    patch_rounds: dict[int, list[list[int]]] = {}
    for k in range(N_OUTPUTS):
        codeword = [OUTPUT[k]] + patch_data[k]
        patch_rounds[k] = _temporally_compared_readout(
            builder, codeword, output_distance
        )

    # output_distance rounds of source X readout with per-wire consecutive
    # round comparisons; the last round supplies the teleportation frame.
    source_rounds = _temporally_compared_readout(
        builder, list(DATA), output_distance
    )
    final_source_records = source_rounds[-1]

    for k, logical_x_support in enumerate(LX):
        records = list(patch_rounds[k][-1]) + [
            final_source_records[q] for q in sorted(logical_x_support)
        ]
        builder.observable(records, k)

    direct_output_records: list[int] | None = None
    if debug_direct_output_mx:
        direct_output_records = builder.measure_many(OUTPUT, basis="X")

    metadata: dict[str, object] = {
        "output_distance": output_distance,
        "codeword_wires": {
            k: [OUTPUT[k]] + patch_data[k] for k in range(N_OUTPUTS)
        },
        "patch_check_ancillas": dict(patch_ancillas),
        "final_patch_records": {k: patch_rounds[k][-1] for k in patch_rounds},
        "final_source_records": list(final_source_records),
        "direct_output_records": direct_output_records,
        "num_qubits": builder.circuit.num_qubits,
        "num_detectors": builder.circuit.num_detectors,
        "num_observables": builder.circuit.num_observables,
    }
    if debug_direct_output_mx:
        return builder.circuit, metadata
    return builder.circuit
