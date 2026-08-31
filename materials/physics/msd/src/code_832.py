"""[[8,3,2]] code constants: stabilizers, logicals, transversal CCZ.

Qubits 0..7 = cube vertices. Code (arXiv:2605.21867 Sec. II):
  S_X1 = X^8, S_Z1 = Z^8,
  S_Z2 = Z0 Z1 Z2 Z3, S_Z3 = Z0 Z1 Z4 Z5, S_Z4 = Z0 Z2 Z4 Z6.
  L_X1 = X0X1X2X3, L_X2 = X0X1X4X5, L_X3 = X0X2X4X6
  L_Z1 = Z0Z4,     L_Z2 = Z0Z2,     L_Z3 = Z0Z1.
Transversal CCZ (Eq. 1): T0 T1^dag T2^dag T3 T4^dag T5 T6 T7^dag.
"""

# (x, y, z) coordinates of cube vertex i: i = 4x + 2y + z (Fig. 1 convention).
CUBE_COORDS = {
    0: (0, 0, 0), 1: (0, 0, 1), 2: (0, 1, 0), 3: (0, 1, 1),
    4: (1, 0, 0), 5: (1, 0, 1), 6: (1, 1, 0), 7: (1, 1, 1),
}

# Z-type stabilizer supports (bits over qubits 0..7).
Z_STABS = [
    set(range(8)),            # S_Z1 = Z^8
    {0, 1, 2, 3},             # S_Z2
    {0, 1, 4, 5},             # S_Z3
    {0, 2, 4, 6},             # S_Z4
]
# X-type stabilizer support.
X_STABS = [set(range(8))]     # S_X1 = X^8

# Logical X operators (weight-4, one per logical qubit, aligned with S_Z j+1).
LX = [{0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6}]
# Logical Z operators (weight 2).
LZ = [{0, 4}, {0, 2}, {0, 1}]

# Transversal CCZ: +1 = T, -1 = T^dagger (Eq. 1 order over qubits 0..7).
CCZ_T_SIGN = {0: +1, 1: -1, 2: -1, 3: +1, 4: -1, 5: +1, 6: +1, 7: -1}

# Phase-polynomial data for the CCZ pattern is kept as constants above.  The
# Clifford-only Gate-A smoke harness verifies the stabilizer state; it does not
# attempt a non-Clifford CCZ simulation.

# The physical order used by the identity-substituted Stim circuit.  The
# actual T/T-dagger gates are intentionally not emitted by gate_a.py.
CCZ_T_ORDER = list(range(8))
