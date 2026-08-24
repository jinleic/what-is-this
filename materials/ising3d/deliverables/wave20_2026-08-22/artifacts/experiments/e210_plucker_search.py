#!/usr/bin/env python3
"""Bounded exact search for a four-channel Pluecker/Hirota identity.

[COMPUTATION] The declared ansatz contains exactly the normalized coefficient
vectors ``(1,c12,c13,c14)`` with each free integer coefficient in ``[-2,2]``.
Only a 2x3 planar rectangle fixes the coefficients.  A 3x3 planar rectangle is
a fresh control, and no 3D value participates in selection.

[UNRESOLVED] Failure of this 125-element class is not a no-go for arbitrary
determinants, Pfaffians, basis changes, terminal orders, or higher identities.
"""

from __future__ import annotations

import hashlib
import itertools
import sys
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from ising.lattices import Lattice, hyperrect  # noqa: E402

from e209_plucker_tensor import (  # noqa: E402
    Budget,
    PLUCKER_COEFFICIENTS,
    lower_priority,
    max_rss_bytes,
)

MPoly = dict[int, int]
COEFFICIENT_BOUND = 2
MAX_MULTIVARIATE_EDGES = 12


def _validated_terminals(lat: Lattice, terminals: Sequence[int]) -> tuple[int, ...]:
    result = tuple(int(site) for site in terminals)
    if len(result) != 4 or len(set(result)) != 4:
        raise ValueError("four distinct terminals are required")
    if any(site < 0 or site >= lat.n_sites for site in result):
        raise ValueError("terminal index is outside the lattice")
    return result


def multivariate_currents(
    lat: Lattice,
    terminals: Sequence[int],
    budget: Budget | None = None,
) -> list[MPoly]:
    """Return exact currents in independent edge variables using base-three monomial codes.

    A selected edge contributes digit one in its code.  Products of two current
    polynomials are ordinary integer addition of the codes, whose digits then
    lie in ``{0,1,2}``; no carry is possible.
    """

    if lat.n_bonds > MAX_MULTIVARIATE_EDGES:
        raise ValueError(
            f"bounded multivariate search permits at most {MAX_MULTIVARIATE_EDGES} edges"
        )
    terminal_sites = _validated_terminals(lat, terminals)
    positions = {site: position for position, site in enumerate(terminal_sites)}
    terminal_set = set(terminal_sites)
    interior_mask = sum(1 << site for site in range(lat.n_sites) if site not in terminal_set)
    powers_of_three = [3**edge for edge in range(lat.n_bonds)]
    result: list[MPoly] = [{} for _ in range(16)]

    for edge_subset in range(1 << lat.n_bonds):
        if budget is not None and edge_subset % 1024 == 0:
            budget.check(
                f"multivariate currents {lat.shape}: subset {edge_subset}/{1 << lat.n_bonds}"
            )
        parity = 0
        monomial_code = 0
        for edge_index, (left, right) in enumerate(lat.bonds):
            if (edge_subset >> edge_index) & 1:
                parity ^= (1 << left) ^ (1 << right)
                monomial_code += powers_of_three[edge_index]
        if parity & interior_mask:
            continue
        boundary_mask = 0
        for site, position in positions.items():
            if (parity >> site) & 1:
                boundary_mask |= 1 << position
        result[boundary_mask][monomial_code] = (
            result[boundary_mask].get(monomial_code, 0) + 1
        )
    return result


def mpoly_mul(left: MPoly, right: MPoly) -> MPoly:
    result: MPoly = {}
    for left_code, left_coefficient in left.items():
        for right_code, right_coefficient in right.items():
            code = left_code + right_code
            result[code] = (
                result.get(code, 0) + left_coefficient * right_coefficient
            )
    return {code: coefficient for code, coefficient in result.items() if coefficient}


def multivariate_channels(currents: Sequence[MPoly]) -> list[MPoly]:
    if len(currents) != 16:
        raise ValueError("the current signature must have 16 entries")
    return [
        mpoly_mul(currents[0b0000], currents[0b1111]),
        mpoly_mul(currents[0b0011], currents[0b1100]),
        mpoly_mul(currents[0b0101], currents[0b1010]),
        mpoly_mul(currents[0b1001], currents[0b0110]),
    ]


def mpoly_linear_combination(
    channels: Sequence[MPoly], coefficients: Sequence[int]
) -> MPoly:
    if len(channels) != 4 or len(coefficients) != 4:
        raise ValueError("four channels and four coefficients are required")
    result: MPoly = {}
    for scale, channel in zip(coefficients, channels, strict=True):
        for code, coefficient in channel.items():
            result[code] = result.get(code, 0) + int(scale) * coefficient
    return {code: coefficient for code, coefficient in result.items() if coefficient}


def coefficient_candidates(bound: int = COEFFICIENT_BOUND) -> Iterable[tuple[int, int, int, int]]:
    if bound < 0:
        raise ValueError("coefficient bound must be nonnegative")
    for tail in itertools.product(range(-bound, bound + 1), repeat=3):
        yield (1, *tail)


def bounded_search(
    training_channels: Sequence[Sequence[MPoly]],
    bound: int = COEFFICIENT_BOUND,
) -> dict[str, object]:
    """Search only the declared finite normalized integer coefficient box."""

    candidates = list(coefficient_candidates(bound))
    survivors = [
        coefficients
        for coefficients in candidates
        if all(
            not mpoly_linear_combination(channels, coefficients)
            for channels in training_channels
        )
    ]
    return {
        "tag": "[COMPUTATION]",
        "normalization": "coefficient of C_empty*C_1234 is 1",
        "free_integer_coefficient_interval": [-bound, bound],
        "candidate_count": len(candidates),
        "survivors": [list(coefficients) for coefficients in survivors],
    }


def decode_ternary(code: int, edge_count: int) -> list[int]:
    if code < 0:
        raise ValueError("monomial code must be nonnegative")
    digits: list[int] = []
    remaining = int(code)
    for _ in range(edge_count):
        digits.append(remaining % 3)
        remaining //= 3
    if remaining:
        raise ValueError("monomial code exceeds the declared edge count")
    return digits


def uniform_specialization(poly: MPoly, edge_count: int) -> list[int]:
    coefficients = [0] * (2 * edge_count + 1)
    for code, coefficient in poly.items():
        degree = sum(decode_ternary(code, edge_count))
        coefficients[degree] += coefficient
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def mpoly_digest(poly: MPoly) -> str:
    payload = "\n".join(
        f"{code}:{poly[code]}" for code in sorted(poly)
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def first_witness(poly: MPoly, edge_count: int) -> dict[str, object] | None:
    if not poly:
        return None
    degree, code = min(
        (sum(decode_ternary(code, edge_count)), code) for code in poly
    )
    return {
        "monomial_code_base3": code,
        "edge_exponents": decode_ternary(code, edge_count),
        "total_degree": degree,
        "coefficient": poly[code],
    }


def cyclic_corner_terminals(lat: Lattice) -> tuple[int, ...]:
    if len(lat.shape) != 2 or min(lat.shape) < 2:
        raise ValueError("a nondegenerate open rectangle is required")
    rows, cols = lat.shape
    return tuple(
        lat.index(coord)
        for coord in ((0, 0), (0, cols - 1), (rows - 1, cols - 1), (rows - 1, 0))
    )


def _main() -> int:
    lower_priority()
    budget = Budget()

    training = hyperrect((2, 3), periodic=False)
    training_currents = multivariate_currents(
        training, cyclic_corner_terminals(training), budget
    )
    training_channels = multivariate_channels(training_currents)
    search = bounded_search([training_channels])

    holdout = hyperrect((3, 3), periodic=False)
    holdout_currents = multivariate_currents(
        holdout, cyclic_corner_terminals(holdout), budget
    )
    holdout_residual = mpoly_linear_combination(
        multivariate_channels(holdout_currents), PLUCKER_COEFFICIENTS
    )

    cube = hyperrect((2, 2, 2), periodic=False)
    cube_terminals = tuple(
        cube.index(coord)
        for coord in ((0, 0, 0), (0, 0, 1), (1, 1, 1), (1, 1, 0))
    )
    cube_currents = multivariate_currents(cube, cube_terminals, budget)
    cube_residual = mpoly_linear_combination(
        multivariate_channels(cube_currents), PLUCKER_COEFFICIENTS
    )
    joint_survivors = [
        coefficients
        for coefficients in coefficient_candidates(COEFFICIENT_BOUND)
        if not mpoly_linear_combination(training_channels, coefficients)
        and not mpoly_linear_combination(multivariate_channels(cube_currents), coefficients)
    ]

    checks = {
        "[COMPUTATION] planar training selects the Pfaffian signs": search["survivors"]
        == [list(PLUCKER_COEFFICIENTS)],
        "[COMPUTATION] fresh planar multivariate holdout passes": not holdout_residual,
        "[COMPUTATION] skew-cube multivariate residual is nonzero": bool(cube_residual),
        "[COMPUTATION] no bounded candidate passes training and cube": not joint_survivors,
        "[COMPUTATION] resource budget": budget.used() < budget.cpu_seconds
        and max_rss_bytes() < budget.rss_bytes,
    }
    for name, passed in checks.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}", flush=True)
    print(
        f"[COMPUTATION] searched {search['candidate_count']} candidates; "
        f"cube multivariate residual has {len(cube_residual)} terms; "
        f"witness={first_witness(cube_residual, cube.n_bonds)}",
        flush=True,
    )
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(_main())
