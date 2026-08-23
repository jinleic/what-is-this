#!/usr/bin/env python3
"""Exact rational torus Green kernels for the mode-resolved four-point front.

This component computes only finite certificates.  It reuses the audited exact
quadratic-field Fourier engine from e93, but stores and checks the resulting
Green kernels over Q.  It does not write the combined result artifact.
"""

from __future__ import annotations

import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

from e93_mag_floor import Quad, torus_green  # noqa: E402

SCRIPT = "experiments/e215_fourpoint_kernel.py"
SIDES = (4, 6, 8, 12)
CPU_BUDGET_SECONDS = 90.0
RSS_LIMIT_BYTES = 2 * 1024**3


def _as_fraction(value: Fraction | Quad) -> Fraction:
    if isinstance(value, Quad):
        return value.rational()
    return value


def _peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # Darwin reports bytes; Linux and most BSD-derived CI images report KiB.
    return value if sys.platform == "darwin" else value * 1024


def _record(
    checks: list[dict[str, Any]], name: str, condition: bool, detail: str
) -> None:
    passed = bool(condition)
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    checks.append({"name": name, "passed": passed, "detail": detail})


def _defect_equation_holds(
    side: int, green: dict[tuple[int, int, int], Fraction]
) -> bool:
    volume = side**3
    for site, value in green.items():
        lhs = 3 * value
        for axis in range(3):
            for step in (-1, 1):
                neighbour = list(site)
                neighbour[axis] = (neighbour[axis] + step) % side
                lhs -= Fraction(1, 2) * green[tuple(neighbour)]
        rhs = (Fraction(1) if site == (0, 0, 0) else Fraction(0)) - Fraction(
            1, volume
        )
        if lhs != rhs:
            return False
    return True


def build_kernel_data(
    sides: tuple[int, ...] = SIDES,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    started = time.process_time()
    checks: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []

    for side in sides:
        raw_green, orbits, inherited_controls = torus_green(side)
        green = {site: _as_fraction(value) for site, value in raw_green.items()}
        volume = side**3
        c0 = green[(0, 0, 0)]
        minimum = min(green.values())
        maximum = max(green.values())
        diameter = c0 - minimum
        minimizers = sorted(site for site, value in green.items() if value == minimum)

        orbit_rows: list[dict[str, Any]] = []
        for orbit in orbits:
            representative = orbit[0]
            value = green[representative]
            if any(green[site] != value for site in orbit):
                raise AssertionError(f"side {side}: orbit constancy failed")
            orbit_rows.append(
                {
                    "representative": list(representative),
                    "multiplicity": len(orbit),
                    "value": str(value),
                }
            )

        zero_sum = sum(
            (len(orbit) * green[orbit[0]] for orbit in orbits), Fraction(0)
        )
        defect_ok = _defect_equation_holds(side, green)
        coverage = sum(len(orbit) for orbit in orbits)

        _record(
            checks,
            f"L{side}_rational_kernel",
            all(isinstance(value, Fraction) for value in green.values()),
            f"all {volume} values lie in Q",
        )
        _record(
            checks,
            f"L{side}_orbit_partition",
            coverage == volume and len(green) == volume,
            f"{len(orbits)} signed-permutation orbits cover {volume} sites",
        )
        _record(
            checks,
            f"L{side}_poisson_certificate",
            defect_ok and bool(inherited_controls["defect"]),
            "(3-A/2) C = delta_0-1/N at every site",
        )
        _record(
            checks,
            f"L{side}_zero_sum",
            zero_sum == 0 and bool(inherited_controls["zero_sum"]),
            "the zero Fourier mode is removed exactly",
        )
        _record(
            checks,
            f"L{side}_diameter",
            c0 == maximum and minimum < 0 < c0 and diameter > c0,
            f"C(0)={c0}, min C={minimum}, D={diameter}",
        )

        records.append(
            {
                "side": side,
                "volume": volume,
                "field": "Q",
                "orbit_count": len(orbits),
                "orbits": orbit_rows,
                "C0": str(c0),
                "minimum": str(minimum),
                "minimum_sites": [list(site) for site in minimizers],
                "diameter_D": str(diameter),
                "kappa_D_over_2": str(diameter / 2),
                "zero_sum": str(zero_sum),
            }
        )

        cpu = time.process_time() - started
        rss = _peak_rss_bytes()
        if cpu > CPU_BUDGET_SECONDS:
            raise RuntimeError(
                f"kernel CPU budget exceeded after L={side}: {cpu:.3f}s"
            )
        if rss >= RSS_LIMIT_BYTES:
            raise MemoryError(f"kernel RSS limit exceeded after L={side}: {rss}")

    _record(
        checks,
        "kernel_sides_complete",
        tuple(row["side"] for row in records) == tuple(sides),
        f"computed exact kernels for sides {tuple(sides)}",
    )

    data = {
        "definition": {
            "torus": "T_L=(Z/LZ)^3",
            "dispersion": "lambda(k)=3-sum_i cos(k_i)",
            "kernel": "C_L(z)=N^-1 sum_{k!=0} cos(k.z)/lambda(k)",
            "rationality_certificate": (
                "the rational Poisson system (3-A/2)C=delta_0-1/N with "
                "sum_z C(z)=0 has a unique solution"
            ),
        },
        "sides": records,
        "resource_observation": {
            "process_time_seconds": f"{time.process_time() - started:.6f}",
            "peak_rss_bytes": _peak_rss_bytes(),
            "process_time_budget_seconds": CPU_BUDGET_SECONDS,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
        },
    }
    return data, checks


def main() -> int:
    data, checks = build_kernel_data()
    for check in checks:
        print(f"PASS {check['name']} -- {check['detail']}")
    print(
        json.dumps(
            {
                "script": SCRIPT,
                "sides": [row["side"] for row in data["sides"]],
                "check_count": len(checks),
            },
            sort_keys=True,
        )
    )
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
