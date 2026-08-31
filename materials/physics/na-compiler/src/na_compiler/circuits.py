"""Deterministic circuit instances and a minimal QASM interaction loader.

Gate-A QFT inputs are generated with qiskit 2.5.2 and transpiled at
optimization_level=0 to ``{u,cz}``; only the ordered CZ interaction stream is
needed by the transport compiler. Gate-A 3-regular instances are generated
with networkx 3.6.1 and fixed seeds. ``load_qasm`` accepts OpenQASM 2.x and
returns the same ordered interaction graph/layers used by the drivers.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Iterable


@dataclass
class Instance:
    """Bundled circuit instance.

    ``pairs`` is the ordered ASAP-layered CZ interaction stream. ``qasm`` is a
    native-QASM rendering retained in campaign artifacts, making every
    quantitative row replayable without relying on a Python circuit object.
    """

    name: str
    n_qubits: int
    pairs: list[list[tuple[int, int]]]
    family: str
    provenance: dict = field(default_factory=dict)
    qasm: str = ""

    def to_json(self) -> str:
        return json.dumps(
            {
                "name": self.name,
                "n_qubits": self.n_qubits,
                "pairs": self.pairs,
                "family": self.family,
                "provenance": self.provenance,
                "qasm": self.qasm,
            },
            sort_keys=True,
        )

    def interaction_edges(self) -> list[tuple[int, int]]:
        return [pair for layer in self.pairs for pair in layer]


def pairs_to_qasm(n_qubits: int, pairs: Iterable[Iterable[tuple[int, int]]]) -> str:
    """Render a native OpenQASM 2 circuit from an ordered pair stream."""
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n_qubits}];"]
    for layer in pairs:
        for a, b in layer:
            lines.append(f"cz q[{a}],q[{b}];")
    return "\n".join(lines) + "\n"


def load_qasm(text_or_path: str) -> Instance:
    """Load OpenQASM 2.x into an interaction graph + deterministic ASAP layers.

    Supported two-qubit operations are ``cz``, ``cx/cnot``, ``cp``, and
    ``swap`` (all become an interaction edge); one-qubit operations and
    barriers are accepted and ignored for transport. Register references must
    be ``q[index]`` or ``name[index]`` with a single quantum register. The
    loader deliberately rejects unsupported/malformed operations instead of
    silently dropping them.
    """
    import pathlib

    source = pathlib.Path(text_or_path)
    if source.exists():
        text = source.read_text()
        source_name = str(source)
    else:
        text = text_or_path
        source_name = "inline"
    qreg = re.search(r"\bqreg\s+(\w+)\s*\[\s*(\d+)\s*\]\s*;", text)
    if not qreg:
        raise ValueError("QASM must declare one qreg")
    reg_name, n_qubits_s = qreg.groups()
    n_qubits = int(n_qubits_s)
    twoq: list[list[tuple[int, int]]] = []
    op_re = re.compile(r"^\s*([A-Za-z_][\w]*)\s*(?:\([^;]*\))?\s+([^;]+);\s*$")
    qref_re = re.compile(rf"{re.escape(reg_name)}\s*\[\s*(\d+)\s*\]")
    # A qasm line can include // comments; strip only comments outside our
    # tiny grammar (arguments are numeric for the supported operations).
    for raw in text.splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line or line.startswith(("OPENQASM", "include", "qreg", "creg")):
            continue
        m = op_re.match(line)
        if not m:
            raise ValueError(f"unsupported/malformed QASM line: {raw!r}")
        op, args = m.groups()
        refs = [int(x) for x in qref_re.findall(args)]
        if any(x < 0 or x >= n_qubits for x in refs):
            raise ValueError(f"qubit index out of range: {raw!r}")
        op_l = op.lower()
        if op_l in {"barrier", "measure", "reset"}:
            continue
        if op_l in {"cz", "cx", "cnot", "cp", "swap"}:
            if len(refs) != 2 or refs[0] == refs[1]:
                raise ValueError(f"two-qubit operation needs two distinct qrefs: {raw!r}")
            twoq.append([(min(refs), max(refs))])
        elif len(refs) == 1:
            continue
        else:
            raise ValueError(f"unsupported operation {op!r}: {raw!r}")
    layers = asap_merge(twoq)
    return Instance(
        name=source_name,
        n_qubits=n_qubits,
        pairs=layers,
        family="qasm",
        provenance={"source": source_name, "loader": "na_compiler.circuits.load_qasm"},
        qasm=text,
    )


def qft_instance(n: int) -> Instance:
    """Canonical QFT(n), qiskit 2.5.2, no optimization, CZ-native stream."""
    from qiskit import transpile
    from qiskit.circuit.library import QFT

    if n < 2:
        raise ValueError("QFT benchmark requires n >= 2")
    circ = QFT(n, do_swaps=True)
    transpiled = transpile(circ, basis_gates=["u", "cz"], optimization_level=0)
    sequential: list[list[tuple[int, int]]] = []
    for inst in transpiled.data:
        if len(inst.qubits) != 2:
            continue
        qi = [transpiled.find_bit(q).index for q in inst.qubits]
        if inst.operation.name not in {"cz", "cx", "swap"}:
            raise ValueError(f"unexpected QFT two-qubit gate {inst.operation.name}")
        sequential.append([(min(qi), max(qi))])
    layers = asap_merge(sequential)
    return Instance(
        name=f"qft{n}",
        n_qubits=n,
        pairs=layers,
        family="qft",
        provenance={
            "source": "qiskit.circuit.library.QFT",
            "qiskit_version": __import__("qiskit").__version__,
            "transpile": "basis=[u,cz], optimization_level=0, do_swaps=True",
        },
        qasm=pairs_to_qasm(n, layers),
    )


def reg3_instance(n: int, seed: int) -> Instance:
    """Connected 3-regular random graph (networkx) with fixed seed."""
    if n < 4 or n % 2:
        raise ValueError("3-regular benchmark requires even n >= 4")
    import networkx as nx

    # random_regular_graph can produce disconnected graphs; deterministic
    # seed stepping makes the provenance explicit while guaranteeing the
    # requested connected family.
    chosen_seed = seed
    for attempt in range(10_000):
        g = nx.random_regular_graph(3, n, seed=seed + attempt)
        if nx.is_connected(g):
            chosen_seed = seed + attempt
            break
    else:
        raise RuntimeError(f"failed to generate connected 3-regular graph n={n}")
    edges = sorted((min(u, v), max(u, v)) for u, v in g.edges())
    sequential = [[edge] for edge in edges]
    layers = asap_merge(sequential)
    return Instance(
        name=f"reg3_{n}_s{seed}",
        n_qubits=n,
        pairs=layers,
        family="reg3",
        provenance={
            "source": "networkx.random_regular_graph",
            "networkx_version": nx.__version__,
            "n": n,
            "requested_seed": seed,
            "chosen_seed": chosen_seed,
            "edge_order": "sorted lexicographically",
        },
        qasm=pairs_to_qasm(n, layers),
    )


def asap_merge(sequential_pairs: list[list[tuple[int, int]]]) -> list[list[tuple[int, int]]]:
    """Match QMAP ASAPScheduler for CZ-only streams without barriers.

    ``next_layer_for_qubit[q]`` is the first layer after q's last two-qubit
    gate. The entanglement capacity limit is intentionally *not* applied here;
    QMAP's architecture_for() sizes capacity to avoid splitting a benchmark
    layer.
    """
    next_layer: dict[int, int] = {}
    layers: list[list[tuple[int, int]]] = []
    for group in sequential_pairs:
        for a, b in group:
            layer = max(next_layer.get(a, 0), next_layer.get(b, 0))
            while layer >= len(layers):
                layers.append([])
            # The nextLayer values already prevent operand collisions, but the
            # explicit check keeps behavior correct for hand-built inputs.
            while any(a in pair or b in pair for pair in layers[layer]):
                layer += 1
                while layer >= len(layers):
                    layers.append([])
            layers[layer].append((a, b))
            next_layer[a] = layer + 1
            next_layer[b] = layer + 1
    return layers
