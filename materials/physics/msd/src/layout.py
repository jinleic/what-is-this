"""Logical roles for the 22 physical wires in the Gate-A model.

The paper counts 22 physical circuit qubits and three logical surface-code
outputs (arXiv:2605.21867, Sec. III).  The output blocks are not expanded in
the Clifford-only Stim model.  Instead, three wires from the eight-wire
verification block are reused as effective output representatives after the
second verification CNOT.  This keeps the counted circuit wires at 22 while
making the omitted surface-code decoder boundary explicit.

The role partition is therefore:

* 0..7: [[8,3,2]] data block, at cube vertices.
* 8..15: eight |+> verification wires; 8..10 are reused as output
  representatives after verification.
* 16..19: four syndrome wires.
* 20..21: reusable surgery-parity wires.

The role reuse and the unexpanded surface-code patches are [INFERENCE], not
claims that the paper's Appendix-A space-time layout has these exact labels.
"""

DATA = list(range(8))
VERIFY = list(range(8, 16))
OUTPUT = VERIFY[:3]
SYNDROME = list(range(16, 20))
SURGERY = [20, 21]
ALL_QUBITS = list(range(22))
N_PHYSICAL = len(ALL_QUBITS)

# Figure 4(a) has these three depth groups.  The text calls the encoder
# "eight CNOT gates", while the rendered figure contains ten CNOT symbols.
# This list is transcribed from the rendered symbols and is checked by
# preparing all eight +1 stabilizers in smoke.py.
ENCODER_LAYERS = [
    [(4, 0), (6, 7), (1, 5), (3, 2)],
    [(0, 2), (7, 5), (3, 1), (4, 6)],
    [(5, 4), (1, 0)],
]
ENCODER_CX = [pair for layer in ENCODER_LAYERS for pair in layer]

# Four stabilizer checks used by the superdense-syndrome scaffold.  The
# first is X-type and the remaining three are Z-type; supports come directly
# from Eq. (2) of arXiv:2605.21867.
SUPERDENSE_CHECKS = [
    ("X", set(range(8))),
    ("Z", {0, 1, 2, 3}),
    ("Z", {0, 1, 4, 5}),
    ("Z", {0, 2, 4, 6}),
]
