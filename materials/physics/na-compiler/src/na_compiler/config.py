"""Persisted session configuration for na-compiler.

Values recorded 2026-08-29 from first-hand reads:
- MQT QMAP source (github main) + mqt.qmap 3.9.0 wheel behavior.
- arXiv:2505.22715v1 (routing-aware placement) full text.
- arXiv:2604.25478v2 (RSQASM evaluation frame) full text.
"""

# --- MQT QMAP shipped architecture (eval/na/zoned/square_architecture.json) ---
QMAP_SQUARE_ARCH = {
    "name": "na_compiler_scaled_square",
    "operation_duration": {
        "rydberg_gate": 0.36,
        "single_qubit_gate": 52,
        "atom_transfer": 15,
    },
    "operation_fidelity": {
        "rydberg_gate": 0.995,
        "single_qubit_gate": 0.9997,
        "atom_transfer": 0.999,
    },
    "qubit_spec": {"T": 1.5e6},
    "storage_zones": [],
    "entanglement_zones": [],
    "aods": [{"id": 0, "site_separation": 2, "r": 100, "c": 100}],
}

# Physical separations kept identical to QMAP's shipped arch (units um).
STORAGE_SITE_SEPARATION = (4.0, 4.0)
ENTANGLEMENT_SITE_SEPARATION = (12.0, 10.0)

# --- Duration model constants (QMAP official evaluator,
# eval/na/zoned/eval_ids_relaxed_routing.py, read 2026-08-29) ---
TIME_ATOM_TRANSFER_US = 15.0  # per load or store event
# jerk-limited cubic profile for d <= 110 um: t = 2*(4d/j)^(1/3)
T_D_MAX_US = 200.0
D_MAX_UM = 110.0
JERK_UM_PER_US3 = 32.0 * D_MAX_UM / T_D_MAX_US**3
V_MAX_UM_PER_US = D_MAX_UM / T_D_MAX_US * 2.0  # 1.1


def move_batch_duration_us(d_um: float) -> float:
    """QMAP evaluator's move-batch duration for max travel distance d (um)."""
    if d_um <= D_MAX_UM:
        return 2.0 * (4.0 * d_um / JERK_UM_PER_US3) ** (1.0 / 3.0)
    return T_D_MAX_US + (d_um - D_MAX_UM) / V_MAX_UM_PER_US


def transport_duration_us(
    batch_max_distances: "list[float]", n_load_store_events: int
) -> float:
    """Total transport cost: sum of batch durations + 15us per load/store."""
    return (
        sum(move_batch_duration_us(d) for d in batch_max_distances)
        + TIME_ATOM_TRANSFER_US * n_load_store_events
    )
