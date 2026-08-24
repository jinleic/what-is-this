#!/usr/bin/env python3
"""Exact four-boundary Ising tensors in the high-temperature basis.

[LEMMA] For a finite repository lattice ``G=(V,E)`` and four fixed terminal
spins, this module constructs

    T_b(v) = sum_{sigma on V minus B} product_{uw in E} (1 + v sigma_u sigma_w).

[LEMMA] Its normalized terminal Walsh transform is the exact current signature
``C_S(v)`` that counts edge subsets with odd-degree boundary ``S`` and even
interior degree.  The two constructions are implemented independently.

[UNRESOLVED] No determinant, Pfaffian, or all-size representation is assumed by
this engine; it only supplies exact finite tensors to the bounded tests in
``e210`` and ``e211``.
"""

from __future__ import annotations

import os
import platform
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from ising.lattices import Lattice, hyperrect  # noqa: E402

Poly = list[int]
FOUR_TERMINALS = 4
CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 2_000_000_000
PLUCKER_COEFFICIENTS = (1, -1, 1, -1)


def max_rss_bytes() -> int:
    """Return peak resident bytes under the Darwin/Linux ru_maxrss conventions."""

    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class Budget:
    """Single-process CPU/RSS gate; elapsed wall time is never certifying."""

    def __init__(
        self,
        cpu_seconds: float = CPU_BUDGET_SECONDS,
        rss_bytes: int = RSS_CAP_BYTES,
    ) -> None:
        self.cpu_seconds = float(cpu_seconds)
        self.rss_bytes = int(rss_bytes)
        self.started = time.process_time()
        self.stage = "initialization"

    def used(self) -> float:
        return time.process_time() - self.started

    def check(self, stage: str) -> None:
        self.stage = stage
        used = self.used()
        rss = max_rss_bytes()
        if used >= self.cpu_seconds:
            raise RuntimeError(
                f"process-time budget exceeded at {stage}: {used:.3f}>={self.cpu_seconds}"
            )
        if rss >= self.rss_bytes:
            raise RuntimeError(f"RSS cap exceeded at {stage}: {rss}>={self.rss_bytes}")


def lower_priority(increment: int = 10) -> int | None:
    """Increase niceness when the host permits it and return the resulting value."""

    try:
        os.nice(increment)
        return os.getpriority(os.PRIO_PROCESS, 0)
    except (AttributeError, OSError):
        return None


def trim(poly: Sequence[int]) -> Poly:
    result = [int(value) for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def poly_add(left: Sequence[int], right: Sequence[int]) -> Poly:
    size = max(len(left), len(right))
    result = [0] * size
    for index, value in enumerate(left):
        result[index] += int(value)
    for index, value in enumerate(right):
        result[index] += int(value)
    return trim(result)


def poly_linear_combination(terms: Iterable[tuple[int, Sequence[int]]]) -> Poly:
    materialized = [(int(scale), poly) for scale, poly in terms]
    if not materialized:
        return [0]
    size = max(len(poly) for _, poly in materialized)
    result = [0] * size
    for scale, poly in materialized:
        if not scale:
            continue
        for index, value in enumerate(poly):
            result[index] += scale * int(value)
    return trim(result)


def poly_mul(left: Sequence[int], right: Sequence[int]) -> Poly:
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        if not left_value:
            continue
        for right_degree, right_value in enumerate(right):
            if right_value:
                result[left_degree + right_degree] += int(left_value) * int(right_value)
    return trim(result)


def poly_pow(base: Sequence[int], exponent: int) -> Poly:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    result = [1]
    factor = trim(base)
    power = int(exponent)
    while power:
        if power & 1:
            result = poly_mul(result, factor)
        power >>= 1
        if power:
            factor = poly_mul(factor, factor)
    return result


def monomial(degree: int, coefficient: int = 1) -> Poly:
    if degree < 0:
        raise ValueError("monomial degree must be nonnegative")
    return [0] * degree + [int(coefficient)]


def poly_at_fraction(poly: Sequence[int], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + int(coefficient)
    return result


def _validate_terminals(lat: Lattice, terminals: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(site) for site in terminals)
    if len(result) != FOUR_TERMINALS:
        raise ValueError(f"exactly {FOUR_TERMINALS} terminals are required")
    if len(set(result)) != len(result):
        raise ValueError("terminals must be distinct")
    if any(site < 0 or site >= lat.n_sites for site in result):
        raise ValueError("terminal index is outside the lattice")
    return result


def _alignment_kernels(edge_count: int) -> list[Poly]:
    plus_powers: list[Poly] = [[1]]
    minus_powers: list[Poly] = [[1]]
    for _ in range(edge_count):
        plus_powers.append(poly_mul(plus_powers[-1], [1, 1]))
        minus_powers.append(poly_mul(minus_powers[-1], [1, -1]))
    return [
        poly_mul(plus_powers[aligned], minus_powers[edge_count - aligned])
        for aligned in range(edge_count + 1)
    ]


def boundary_partition_tensor(
    lat: Lattice,
    terminals: Sequence[int],
    budget: Budget | None = None,
) -> list[Poly]:
    """Enumerate the exact fixed-terminal tensor ``T_b(v)`` over integer polynomials.

    Repository convention is retained: a set state bit means spin ``+1``.  The
    terminal tensor index uses the caller's terminal order, with the same bit
    convention.  Only the at-most-4096-state lattices declared in ``e211`` are
    used by this front.
    """

    terminal_sites = _validate_terminals(lat, terminals)
    terminal_position = {site: position for position, site in enumerate(terminal_sites)}
    kernels = _alignment_kernels(lat.n_bonds)
    result = [[0] * (lat.n_bonds + 1) for _ in range(1 << FOUR_TERMINALS)]

    for state in range(1 << lat.n_sites):
        if budget is not None and state % 1024 == 0:
            budget.check(f"spin tensor {lat.shape}: state {state}/{1 << lat.n_sites}")
        boundary_mask = 0
        for site, position in terminal_position.items():
            if (state >> site) & 1:
                boundary_mask |= 1 << position
        aligned = sum(
            ((state >> left) & 1) == ((state >> right) & 1)
            for left, right in lat.bonds
        )
        row = result[boundary_mask]
        for degree, coefficient in enumerate(kernels[aligned]):
            row[degree] += coefficient
    return [trim(poly) for poly in result]


def walsh_currents(tensor: Sequence[Sequence[int]], site_count: int) -> list[Poly]:
    """Return ``2^(-|V|)`` times the terminal-spin Walsh transform of ``tensor``."""

    if len(tensor) != 1 << FOUR_TERMINALS:
        raise ValueError("the four-terminal tensor must have 16 entries")
    divisor = 1 << int(site_count)
    degree_size = max(len(poly) for poly in tensor)
    currents: list[Poly] = []
    for subset in range(1 << FOUR_TERMINALS):
        transformed = [0] * degree_size
        selected = subset.bit_count()
        for boundary_mask, poly in enumerate(tensor):
            selected_plus = (subset & boundary_mask).bit_count()
            sign = -1 if (selected - selected_plus) & 1 else 1
            for degree, coefficient in enumerate(poly):
                transformed[degree] += sign * int(coefficient)
        if any(coefficient % divisor for coefficient in transformed):
            raise ArithmeticError("Walsh coefficients are not divisible by 2^|V|")
        currents.append(trim([coefficient // divisor for coefficient in transformed]))
    return currents


def edge_subset_currents(
    lat: Lattice,
    terminals: Sequence[int],
    budget: Budget | None = None,
) -> list[Poly]:
    """Independently build current polynomials by parity-state edge dynamic programming."""

    terminal_sites = _validate_terminals(lat, terminals)
    terminal_position = {site: position for position, site in enumerate(terminal_sites)}
    terminal_set = set(terminal_sites)
    interior_mask = sum(1 << site for site in range(lat.n_sites) if site not in terminal_set)

    parity_states: dict[int, Poly] = {0: [1]}
    for edge_index, (left, right) in enumerate(lat.bonds):
        if budget is not None:
            budget.check(f"edge-current DP {lat.shape}: edge {edge_index + 1}/{lat.n_bonds}")
        toggle = (1 << left) ^ (1 << right)
        next_states: dict[int, Poly] = {}
        for parity, poly in parity_states.items():
            next_states[parity] = poly_add(next_states.get(parity, [0]), poly)
            included = [0] + poly
            target = parity ^ toggle
            next_states[target] = poly_add(next_states.get(target, [0]), included)
        parity_states = next_states

    currents = [[0] for _ in range(1 << FOUR_TERMINALS)]
    for parity, poly in parity_states.items():
        if parity & interior_mask:
            continue
        boundary_mask = 0
        for site, position in terminal_position.items():
            if (parity >> site) & 1:
                boundary_mask |= 1 << position
        if currents[boundary_mask] != [0]:
            raise AssertionError("two parity states mapped to one boundary current")
        currents[boundary_mask] = trim(poly)
    return currents


def plucker_channels(currents: Sequence[Sequence[int]]) -> list[Poly]:
    """Return the four quadratic channels in the standard four-point order.

    The masks are ``empty|1234, 12|34, 13|24, 14|23``.  Thus the Pfaffian
    Pluecker residual has coefficients ``(+1,-1,+1,-1)``.
    """

    if len(currents) != 1 << FOUR_TERMINALS:
        raise ValueError("the current signature must have 16 entries")
    return [
        poly_mul(currents[0b0000], currents[0b1111]),
        poly_mul(currents[0b0011], currents[0b1100]),
        poly_mul(currents[0b0101], currents[0b1010]),
        poly_mul(currents[0b1001], currents[0b0110]),
    ]


def plucker_residual(
    currents: Sequence[Sequence[int]],
    coefficients: Sequence[int] = PLUCKER_COEFFICIENTS,
) -> Poly:
    if len(coefficients) != 4:
        raise ValueError("the four-channel residual needs four coefficients")
    return poly_linear_combination(zip(coefficients, plucker_channels(currents), strict=True))


def _main() -> int:
    lower_priority()
    budget = Budget()
    planar = hyperrect((2, 3), periodic=False)
    planar_terminals = tuple(
        planar.index(coord) for coord in ((0, 0), (0, 2), (1, 2), (1, 0))
    )
    cube = hyperrect((2, 2, 2), periodic=False)
    cube_terminals = tuple(
        cube.index(coord)
        for coord in ((0, 0, 0), (0, 0, 1), (1, 1, 1), (1, 1, 0))
    )

    planar_tensor = boundary_partition_tensor(planar, planar_terminals, budget)
    planar_walsh = walsh_currents(planar_tensor, planar.n_sites)
    planar_direct = edge_subset_currents(planar, planar_terminals, budget)
    cube_tensor = boundary_partition_tensor(cube, cube_terminals, budget)
    cube_walsh = walsh_currents(cube_tensor, cube.n_sites)
    cube_direct = edge_subset_currents(cube, cube_terminals, budget)

    checks = {
        "[COMPUTATION] planar Walsh/direct agreement": planar_walsh == planar_direct,
        "[COMPUTATION] planar Pluecker residual vanishes": plucker_residual(planar_walsh) == [0],
        "[COMPUTATION] cube Walsh/direct agreement": cube_walsh == cube_direct,
        "[COMPUTATION] skew-cube Pluecker residual is nonzero": plucker_residual(cube_walsh) != [0],
        "[COMPUTATION] resource budget": budget.used() < budget.cpu_seconds
        and max_rss_bytes() < budget.rss_bytes,
    }
    for name, passed in checks.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}", flush=True)
    print(f"[COMPUTATION] cube residual = {plucker_residual(cube_walsh)}", flush=True)
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(_main())
