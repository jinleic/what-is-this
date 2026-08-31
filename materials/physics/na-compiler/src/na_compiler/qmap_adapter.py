"""mqt.qmap zoned-compiler adapter and .naviz program parser.

``compile_with_qmap`` runs the RoutingAwareCompiler in strict mode on a
scaled ``Architecture`` and keeps the decompiled program in three forms:

* ``NavizProgram`` — ordered events (atom/load/move/store/cz/u/log), the
  machine-truth record;
*紧凑``raw_moves`` — the QMAP router's actual batch sequence, each batch
  carrying atoms and targets;
* ``raw_events_log`` — the raw .naviz string itself for archival.

The parser asserts structural invariants (load-before-move, store-after-
load) rather than silently repairing the record.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class NavizEvent:
    kind: str  # 'atom' | 'load' | 'move' | 'store' | 'cz' | 'u' | 'log'
    atoms: list[str] = field(default_factory=list)
    targets: dict[str, tuple[float, float]] = field(default_factory=dict)
    zone: str = ""


@dataclass
class NavizProgram:
    events: list[NavizEvent]
    atom_locations: dict[str, tuple[float, float]] = field(default_factory=dict)
    raw_moves: list[dict[str, tuple[float, float]]] = field(default_factory=list)
    raw_events_log: str = ""


@dataclass
class QmapRunResult:
    program: NavizProgram
    stats: dict
    qmap_wall_s: float = 0.0


def compile_with_qmap(
    quimap_arch_json: str,
    pairs: list[list[tuple[int, int]]],
    n_qubits: int,
    routing_method: str = "strict",
) -> QmapRunResult:
    """Run MQT QMAP RoutingAwareCompiler (strict routing) on a CZ stream."""
    import mqt.core.ir as mqt_ir
    from mqt.qmap.na.zoned import (
        RoutingAwareCompiler,
        RoutingMethod,
        ZonedNeutralAtomArchitecture,
    )

    arch = ZonedNeutralAtomArchitecture.from_json_string(quimap_arch_json)
    qc = mqt_ir.QuantumComputation(n_qubits)
    if not pairs:
        raise ValueError("benchmark requires at least one two-qubit gate")
    for layer in pairs:
        for a, b in layer:
            if not (0 <= a < n_qubits and 0 <= b < n_qubits):
                raise ValueError(f"gate ({a}, {b}) out of range for {n_qubits} qubits")
            qc.cz(a, b)
    comp = RoutingAwareCompiler(
        arch,
        routing_method=RoutingMethod.strict
        if routing_method == "strict"
        else RoutingMethod.relaxed,
    )
    import time as _time

    t0 = _time.perf_counter()
    code = comp.compile(qc)
    wall = _time.perf_counter() - t0
    stat = comp.stats()
    program = parse_naviz(code)
    return QmapRunResult(
        program=program,
        stats={"raw_stats": {k: str(v) for k, v in stat.items()}},
        qmap_wall_s=wall,
    )


def parse_naviz(code: str) -> NavizProgram:
    """Parse .naviz output into events; structural assertions required.

    The parser never repairs or drops events: any unrecognized line aborts,
    and every move must reference a loaded (in-AOD) atom.
    """
    events: list[NavizEvent] = []
    atom_locations: dict[str, tuple[float, float]] = {}
    lines = code.splitlines()
    i = 0
    in_aod: set[str] = set()
    raw_moves: list[dict[str, tuple[float, float]]] = []
    while i < len(lines):
        line = lines[i]
        m = re.match(r"atom\s+\((-?[\d.]+),\s*(-?[\d.]+)\)\s+(\w+)", line)
        if m:
            x, y, name = float(m.group(1)), float(m.group(2)), m.group(3)
            atom_locations[name] = (x, y)
            events.append(NavizEvent("atom", atoms=[name]))
            i += 1
            continue
        if line.startswith("@+ load") or line.startswith("@+ store"):
            kind = "load" if line.startswith("@+ load") else "store"
            atoms = []
            if "[" in line:
                i += 1
                while "]" not in lines[i]:
                    atoms.append(lines[i].strip())
                    i += 1
                i += 1
            else:
                m2 = re.match(r"@\+ \w+ (\w+)", line)
                assert m2, f"unparseable {kind}: {line!r}"
                atoms = [m2.group(1)]
                i += 1
            if kind == "load":
                for a in atoms:
                    assert a not in in_aod, f"double load {a}"
                    in_aod.add(a)
            else:
                for a in atoms:
                    assert a in in_aod, f"store without load {a}"
                    in_aod.discard(a)
            events.append(NavizEvent(kind, atoms=atoms))
            continue
        if line.startswith("@+ move"):
            targets: dict[str, tuple[float, float]] = {}
            if "[" in line:
                i += 1
                while "]" not in lines[i]:
                    m3 = re.match(r"\((-?[\d.]+),\s*(-?[\d.]+)\)\s+(\w+)", lines[i].strip())
                    assert m3, f"unparseable move line {lines[i]!r}"
                    targets[m3.group(3)] = (float(m3.group(1)), float(m3.group(2)))
                    i += 1
                i += 1
            else:
                m3 = re.match(
                    r"@\+ move\s+\((-?[\d.]+),\s*(-?[\d.]+)\)\s+(\w+)", line
                )
                assert m3, f"unparseable move: {line!r}"
                targets[m3.group(3)] = (float(m3.group(1)), float(m3.group(2)))
                i += 1
            for a in targets:
                assert a in in_aod, f"move of unloaded atom {a}"
            raw_moves.append(targets)
            events.append(NavizEvent("move", targets=targets))
            continue
        if line.startswith("@+ cz"):
            m4 = re.match(r"@\+ cz (\S+)", line)
            assert m4, f"unparseable cz: {line!r}"
            events.append(NavizEvent("cz", zone=m4.group(1)))
            i += 1
            continue
        if line.startswith("@+ u") or line.startswith("@+ rz"):
            events.append(NavizEvent("u"))
            i += 1
            continue
        # blank lines tolerated, anything else is a parse failure
        assert not line.strip(), f"unrecognized naviz line: {line!r}"
        events.append(NavizEvent("log", atoms=[]))
        i += 1
    return NavizProgram(
        events=events,
        atom_locations=atom_locations,
        raw_moves=raw_moves,
        raw_events_log=code,
    )


def extract_placement_trajectory(
    prog: NavizProgram, atom_names: list[str]
) -> list[dict[str, tuple[float, float]]]:
    """Atom positions at each 2Q-gate layer (position during each cz event);

    index 0 = initial storage placement. Coordinates not yet declared for an
    atom are impossible after ``parse_naviz`` (every atom line predates all
    movement), so here every referenced atom must exist — asserted, never
    defaulted.
    """
    trajectory: list[dict[str, tuple[float, float]]] = []
    positions = dict(prog.atom_locations)
    for name in atom_names:
        assert name in positions, f"atom {name} missing from naviz header"
    trajectory.append(dict(positions))
    for ev in prog.events:
        if ev.kind == "move":
            for a, xy in ev.targets.items():
                positions[a] = xy
        elif ev.kind == "cz":
            trajectory.append({name: positions[name] for name in atom_names})
    return trajectory
