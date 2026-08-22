"""Generic syndrome-extraction circuit compiler for arbitrary (mixed) Pauli
stabilizer codes, targeting Stim.

Construction
------------
One ancilla per stabilizer generator.  For generator g with Pauli P_j on data
qubit j we prepare the ancilla in |+>, apply the controlled-P_j gates
(CX / CY / CZ with the ancilla as control) and measure the ancilla in the X
basis.  Because controlled-P_j and controlled-P_k act on different targets
they commute, so the ordering inside a generator is free; the outcome is the
eigenvalue of  prod_j P_j.

Scheduling
----------
Gates conflict when they share a qubit, but conflict-freedom is NOT sufficient:
controlled-Paulis from different ancillas onto a shared target fail to commute
when the Paulis anticommute, so an arbitrary edge colouring measures the wrong
operators (see ``circuits/scheduling.py``, Lemma C1, and
``notes/failed_routes.md`` FR-002).  Valid schedules must additionally satisfy
an even-crossing parity rule, which strictly raises the achievable depth: for
weight-6 bivariate-bicycle codes the combinatorial bound is 6 but the minimum
valid depth is 7.

``edge_colour_schedule`` is retained ONLY as a counterexample generator for the
regression tests.  Production schedules come from
``scheduling.cpsat_schedule``.

Noise
-----
Circuit-level uniform depolarising noise, matching the standard convention:
  * DEPOLARIZE2(p) after every two-qubit gate
  * DEPOLARIZE1(p) on every idle data/ancilla location in a two-qubit layer
  * X_ERROR(p) / Z_ERROR(p) after reset and before measurement
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import stim

__all__ = ["PauliSupport", "generator_supports", "edge_colour_schedule",
           "CircuitSpec", "build_memory_circuit", "schedule_stats"]

_PAULI_GATE = {"X": "CX", "Y": "CY", "Z": "CZ"}


@dataclass
class PauliSupport:
    """Support of one stabilizer generator: qubit -> 'X'|'Y'|'Z'."""
    index: int
    paulis: dict[int, str]

    @property
    def weight(self) -> int:
        return len(self.paulis)

    @property
    def is_pure_z(self) -> bool:
        return all(p == "Z" for p in self.paulis.values())

    @property
    def is_pure_x(self) -> bool:
        return all(p == "X" for p in self.paulis.values())

    @property
    def is_mixed(self) -> bool:
        return not (self.is_pure_x or self.is_pure_z)


def generator_supports(H: np.ndarray) -> list[PauliSupport]:
    H = np.asarray(H, dtype=np.uint8) & 1
    n = H.shape[1] // 2
    out = []
    for i, row in enumerate(H):
        d: dict[int, str] = {}
        for j in range(n):
            x, z = int(row[j]), int(row[n + j])
            if x and z:
                d[j] = "Y"
            elif x:
                d[j] = "X"
            elif z:
                d[j] = "Z"
        out.append(PauliSupport(i, d))
    return out


# --------------------------------------------------------------------------
# Konig edge colouring of a bipartite multigraph (checks vs data qubits)
# --------------------------------------------------------------------------
def edge_colour_schedule(supports: list[PauliSupport], n_qubits: int) -> list[list[tuple[int, int]]]:
    """Return layers; each layer is a list of (check_index, qubit_index).

    Uses the standard Konig alternating-path algorithm, which colours a
    bipartite multigraph with exactly Delta = max degree colours.
    """
    edges = [(s.index, j) for s in supports for j in s.paulis]
    if not edges:
        return []
    deg_c: dict[int, int] = {}
    deg_q: dict[int, int] = {}
    for c, q in edges:
        deg_c[c] = deg_c.get(c, 0) + 1
        deg_q[q] = deg_q.get(q, 0) + 1
    delta = max(max(deg_c.values()), max(deg_q.values()))

    at_c: dict[int, dict[int, int]] = {c: {} for c in deg_c}   # check -> colour -> qubit
    at_q: dict[int, dict[int, int]] = {q: {} for q in deg_q}   # qubit -> colour -> check
    edge_colour: dict[tuple[int, int], int] = {}

    def free_colour(tbl: dict[int, int]) -> int:
        for col in range(delta):
            if col not in tbl:
                return col
        raise RuntimeError("degree exceeded Delta")

    cap = 4 * len(edges) + 8
    for (c, q) in edges:
        a = free_colour(at_c[c])
        b = free_colour(at_q[q])
        if a != b:
            # walk the a/b alternating path out of q, collecting its edges
            path: list[tuple[int, int, int]] = []      # (check, qubit, colour)
            side, node, look = "q", q, a
            for _ in range(cap):
                tbl = at_q[node] if side == "q" else at_c[node]
                other = tbl.get(look)
                if other is None:
                    break
                path.append((other, node, look) if side == "q" else (node, other, look))
                node = other
                side = "c" if side == "q" else "q"
                look = b if look == a else a
            else:
                raise RuntimeError("alternating path did not terminate")
            # delete the whole path first, then reinsert with a<->b swapped
            for cc, qq, col in path:
                del at_c[cc][col]
                del at_q[qq][col]
            for cc, qq, col in path:
                new = b if col == a else a
                at_c[cc][new] = qq
                at_q[qq][new] = cc
                edge_colour[(cc, qq)] = new
            if a in at_c[c] or a in at_q[q]:
                raise RuntimeError("alternating-path swap failed to free the colour")
        at_c[c][a] = q
        at_q[q][a] = c
        edge_colour[(c, q)] = a

    layers: list[list[tuple[int, int]]] = [[] for _ in range(delta)]
    for (c, q), col in edge_colour.items():
        layers[col].append((c, q))
    # validate disjointness
    for li, layer in enumerate(layers):
        seen_c, seen_q = set(), set()
        for c, q in layer:
            assert c not in seen_c, f"layer {li} reuses check {c}"
            assert q not in seen_q, f"layer {li} reuses qubit {q}"
            seen_c.add(c)
            seen_q.add(q)
    return layers


def schedule_stats(H: np.ndarray,
                   layers: list[list[tuple[int, int]]] | None = None) -> dict:
    """Static resource counts for a stabilizer code's syndrome circuit.

    ``layers`` must be an *actually valid* schedule when depth is wanted.  With
    ``layers=None`` no depth is reported at all -- only the combinatorial lower
    bound -- because the only depth obtainable without a solver is the
    edge-colouring value, which corresponds to an invalid circuit.
    """
    sup = generator_supports(H)
    n = H.shape[1] // 2
    weights = [s.weight for s in sup]
    deg = np.zeros(n, dtype=int)
    for s in sup:
        for j in s.paulis:
            deg[j] += 1
    lower_bound = max(max(weights, default=0), int(deg.max()) if n else 0)
    out = {
        "num_checks": len(sup),
        "num_data_qubits": n,
        "num_ancillas": len(sup),
        "max_check_weight": int(max(weights)),
        "mean_check_weight": float(np.mean(weights)),
        "total_two_qubit_gates": int(sum(weights)),
        "max_qubit_degree": int(deg.max()),
        "depth_lower_bound": lower_bound,
        "two_qubit_layers": (len(layers) if layers is not None else None),
        "num_mixed_checks": sum(1 for s in sup if s.is_mixed),
        "num_pure_x_checks": sum(1 for s in sup if s.is_pure_x),
        "num_pure_z_checks": sum(1 for s in sup if s.is_pure_z),
        "gate_types": sorted({_PAULI_GATE[p] for s in sup for p in s.paulis.values()}),
    }
    return out


# --------------------------------------------------------------------------
# circuit construction
# --------------------------------------------------------------------------
@dataclass
class CircuitSpec:
    H: np.ndarray
    observables: np.ndarray        # (num_obs, 2n) symplectic logical operators
    rounds: int
    p: float
    basis: str = "Z"               # data initialisation / final measurement basis
    layers: list[list[tuple[int, int]]] | None = None

    def __post_init__(self) -> None:
        self.H = np.asarray(self.H, dtype=np.uint8) & 1
        self.observables = np.asarray(self.observables, dtype=np.uint8) & 1


def build_memory_circuit(spec: CircuitSpec) -> tuple[stim.Circuit, dict]:
    """Build a `rounds`-round memory experiment.

    Data qubits 0..n-1, ancillas n..n+r-1 (one per generator).
    Returns (circuit, metadata).
    """
    H = spec.H
    n = H.shape[1] // 2
    sup = generator_supports(H)
    r = len(sup)
    if spec.layers is not None:
        layers = spec.layers
    else:
        # A plain edge colouring is NOT a valid schedule: controlled-Paulis from
        # different ancillas onto a shared target fail to commute when the
        # Paulis anticommute (see circuits/scheduling.py, Lemma C1).  The
        # default must therefore be the parity-enforcing solver.
        from .scheduling import cpsat_schedule, depth_lower_bound, slots_to_layers

        lb = depth_lower_bound(sup, n)
        layers = None
        for T in range(lb, lb + 12):
            res = cpsat_schedule(sup, n, T=T, time_limit_s=120.0)
            if res.slot is not None and res.verification["valid"]:
                layers = slots_to_layers(res.slot, T)
                break
        if layers is None:
            raise RuntimeError(
                "no valid schedule found; pass `layers=` explicitly or widen the "
                "depth search (see experiments/exp005_schedulability_nogo.py)")
    p = spec.p
    anc = {s.index: n + s.index for s in sup}
    all_q = list(range(n + r))

    if spec.basis == "Z":
        det_ready = [s.is_pure_z for s in sup]
    elif spec.basis == "X":
        det_ready = [s.is_pure_x for s in sup]
    else:
        raise ValueError(spec.basis)

    c = stim.Circuit()
    # --- data initialisation ---
    if spec.basis == "Z":
        c.append("R", list(range(n)))
        if p:
            c.append("X_ERROR", list(range(n)), p)
    else:
        c.append("RX", list(range(n)))
        if p:
            c.append("Z_ERROR", list(range(n)), p)

    meas_count = 0
    meas_index: list[dict[int, int]] = []   # per round: check index -> absolute measurement index

    for rd in range(spec.rounds):
        c.append("RX", [anc[i] for i in range(r)])
        if p:
            c.append("Z_ERROR", [anc[i] for i in range(r)], p)
        for layer in layers:
            by_gate: dict[str, list[int]] = {}
            busy: set[int] = set()
            for ci, q in layer:
                g = _PAULI_GATE[sup[ci].paulis[q]]
                by_gate.setdefault(g, []).extend([anc[ci], q])
                busy.add(anc[ci])
                busy.add(q)
            targets_this_layer: list[int] = []
            for g in sorted(by_gate):
                c.append(g, by_gate[g])
                targets_this_layer.extend(by_gate[g])
            if p:
                if targets_this_layer:
                    c.append("DEPOLARIZE2", targets_this_layer, p)
                idle = [q for q in all_q if q not in busy]
                if idle:
                    c.append("DEPOLARIZE1", idle, p)
            c.append("TICK")
        if p:
            c.append("Z_ERROR", [anc[i] for i in range(r)], p)
        c.append("MX", [anc[i] for i in range(r)])
        this_round = {}
        for i in range(r):
            this_round[i] = meas_count
            meas_count += 1
        meas_index.append(this_round)

        # --- detectors ---
        for i in range(r):
            if rd == 0:
                if det_ready[i]:
                    off = meas_count - meas_index[0][i]
                    c.append("DETECTOR", [stim.target_rec(-off)], [i, rd])
            else:
                o1 = meas_count - meas_index[rd][i]
                o0 = meas_count - meas_index[rd - 1][i]
                c.append("DETECTOR", [stim.target_rec(-o1), stim.target_rec(-o0)], [i, rd])

    # --- final data measurement, in the initialisation basis ---
    if p:
        c.append("X_ERROR" if spec.basis == "Z" else "Z_ERROR", list(range(n)), p)
    c.append("M" if spec.basis == "Z" else "MX", list(range(n)))
    data_meas_base = meas_count
    meas_count += n

    # final detectors: reconstruct each basis-compatible check from the data
    for i, s in enumerate(sup):
        if not det_ready[i]:
            continue
        recs = [stim.target_rec(-(meas_count - (data_meas_base + j))) for j in s.paulis]
        recs.append(stim.target_rec(-(meas_count - meas_index[spec.rounds - 1][i])))
        c.append("DETECTOR", recs, [i, spec.rounds])

    # --- observables ---
    nobs = 0
    for oi, obs in enumerate(spec.observables):
        if spec.basis == "Z":
            if obs[:n].any():
                continue           # not measurable by a Z-basis data readout
            supp = np.flatnonzero(obs[n:])
        else:
            if obs[n:].any():
                continue
            supp = np.flatnonzero(obs[:n])
        recs = [stim.target_rec(-(meas_count - (data_meas_base + int(j)))) for j in supp]
        c.append("OBSERVABLE_INCLUDE", recs, nobs)
        nobs += 1

    meta = {
        **schedule_stats(H, layers),
        "n_data": n, "n_ancilla": r, "n_qubits_total": n + r,
        "rounds": spec.rounds, "p": p, "basis": spec.basis,
        # depth of the schedule actually used, never the edge-colouring bound
        "two_qubit_layers": len(layers),
        "two_qubit_gates_per_round": sum(s.weight for s in sup),
        "num_observables": nobs,
        "num_detectors": c.num_detectors,
    }
    return c, meta
