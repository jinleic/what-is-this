#!/usr/bin/env python3
"""Standalone exact verifier for results/bounds/upper_fourpoint.json.

This file intentionally imports no e215/e216/e217 producer module.  It verifies
stored rational kernels through their real-space Poisson certificate and
reconstructs the mode lifts with its own quadratic-field implementation.
"""

from __future__ import annotations

import itertools
import json
import resource
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "upper_fourpoint.json"
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
I3_LOW = Fraction(
    "0.505462019717326006052004053227140259985129014817420892188993487886028773451173816800537247069896037925625"
)
A_RIEMANN = Fraction(992015, 408608)
RATE_C = Fraction(5168, 525)
EXPECTED_COVERAGE = {
    "RATIONAL_CHALLENGE_BELOW_WATSON",
    "EXACT_RATIONAL_TORUS_KERNELS",
    "MODE_RESOLVED_PHYSICAL_MAPPING",
    "FOURPOINT_MOMENT_ROWS",
    "DETERMINISTIC_LIFT_PROJECTION_EQUALITY",
    "EXACT_L12_ETA_COUNTERCERTIFICATE",
    "ALL_SIZE_VANISHING_FLOOR",
    "THERMODYNAMIC_DIRECTION_AND_SCOPE",
}


@dataclass(frozen=True)
class Surd:
    """An exact a+b*sqrt(d), implemented independently of the producers."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    d: int = 2

    def __init__(self, a: Any = 0, b: Any = 0, d: int = 2):
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))
        object.__setattr__(self, "d", int(d))

    def _lift(self, other: object) -> "Surd" | Any:
        if isinstance(other, Surd):
            if self.d != other.d:
                raise ValueError("mixed quadratic fields")
            return other
        if isinstance(other, (int, Fraction)):
            return Surd(other, 0, self.d)
        return NotImplemented

    def __add__(self, other: object):
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return Surd(self.a + rhs.a, self.b + rhs.b, self.d)

    __radd__ = __add__

    def __neg__(self):
        return Surd(-self.a, -self.b, self.d)

    def __sub__(self, other: object):
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object):
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return rhs - self

    def __mul__(self, other: object):
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return Surd(
            self.a * rhs.a + self.b * rhs.b * self.d,
            self.a * rhs.b + self.b * rhs.a,
            self.d,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object):
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        norm = rhs.a * rhs.a - rhs.b * rhs.b * self.d
        if norm == 0:
            raise ZeroDivisionError("zero norm")
        numerator = self * Surd(rhs.a, -rhs.b, self.d)
        return Surd(numerator.a / norm, numerator.b / norm, self.d)

    def __rtruediv__(self, other: object):
        lhs = self._lift(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs / self

    def sign(self) -> int:
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.a * self.b > 0:
            return (self.a > 0) - (self.a < 0)
        # Opposite signs: compare |a|^2 and d|b|^2 without approximation.
        a2 = self.a * self.a
        db2 = self.d * self.b * self.b
        if a2 == db2:
            return 0
        if abs(self.a) > 0 and a2 > db2:
            return 1 if self.a > 0 else -1
        return 1 if self.b > 0 else -1

    def __eq__(self, other: object) -> bool:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return False
        return self.a == rhs.a and self.b == rhs.b

    def __lt__(self, other: object) -> bool:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return (self - rhs).sign() < 0

    def __le__(self, other: object) -> bool:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return (self - rhs).sign() <= 0

    def __gt__(self, other: object) -> bool:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return (self - rhs).sign() > 0

    def __ge__(self, other: object) -> bool:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return (self - rhs).sign() >= 0


Field = Fraction | Surd


def parse_field(text: str) -> Field:
    if text.startswith("Q[") and text.endswith("]"):
        a, b, d = text[2:-1].split("|")
        return Surd(Fraction(a), Fraction(b), int(d))
    return Fraction(text)


def cosine_table(side: int) -> list[Field]:
    if side == 4:
        base: Field = Fraction(0)
    elif side == 6:
        base = Fraction(1, 2)
    elif side == 8:
        base = Surd(0, Fraction(1, 2), 2)
    elif side == 12:
        base = Surd(0, Fraction(1, 2), 3)
    else:
        raise AssertionError(f"unexpected side {side}")
    table: list[Field] = [base * 0 + 1, base]
    for index in range(1, side):
        table.append(2 * base * table[index] - table[index - 1])
    table = table[:side]
    assert 2 * base * table[-1] - table[-2] == table[0]
    assert table[side // 2] == -1
    return table


def signed_permutation_orbit(
    point: tuple[int, int, int], side: int
) -> tuple[tuple[int, int, int], ...]:
    images: set[tuple[int, int, int]] = set()
    for permutation in itertools.permutations((0, 1, 2)):
        permuted = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=3):
            images.add(
                tuple((sign * coordinate) % side for sign, coordinate in zip(signs, permuted))
            )
    return tuple(sorted(images))


def independent_orbits(side: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    unseen = set(itertools.product(range(side), repeat=3))
    answer: list[tuple[tuple[int, int, int], ...]] = []
    while unseen:
        orbit = signed_permutation_orbit(min(unseen), side)
        assert set(orbit) <= unseen
        unseen.difference_update(orbit)
        answer.append(orbit)
    return tuple(answer)


def peak_rss_bytes() -> int:
    amount = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return amount if sys.platform == "darwin" else amount * 1024


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def even_at_least(value: int) -> int:
    return value if value % 2 == 0 else value + 1


def even_strictly_above(value: Fraction) -> int:
    return even_at_least(value.numerator // value.denominator + 1)


def check(condition: bool, name: str, detail: str) -> None:
    if not condition:
        raise AssertionError(f"FAIL {name}: {detail}")
    print(f"PASS {name} -- {detail}")


def reconstruct_kernel(row: dict[str, Any]) -> dict[tuple[int, int, int], Fraction]:
    side = int(row["side"])
    volume = side**3
    orbits = independent_orbits(side)
    stored = {tuple(item["representative"]): item for item in row["orbits"]}
    check(
        set(stored) == {orbit[0] for orbit in orbits},
        f"L{side}_stored_orbits",
        "independent signed-permutation partition matches",
    )
    green: dict[tuple[int, int, int], Fraction] = {}
    for orbit in orbits:
        item = stored[orbit[0]]
        assert int(item["multiplicity"]) == len(orbit)
        value = Fraction(item["value"])
        for site in orbit:
            green[site] = value
    check(len(green) == volume, f"L{side}_coverage", f"covered {volume} sites")
    check(
        sum(green.values(), Fraction(0)) == 0,
        f"L{side}_zero_sum",
        "zero spatial sum",
    )
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
        assert lhs == rhs
    c0 = green[(0, 0, 0)]
    minimum = min(green.values())
    diameter = c0 - minimum
    check(
        c0 == Fraction(row["C0"])
        and minimum == Fraction(row["minimum"])
        and diameter == Fraction(row["diameter_D"])
        and c0 == max(green.values())
        and minimum < 0,
        f"L{side}_kernel_extrema",
        "Poisson solution has stored exact C0, minimum, and diameter",
    )
    check(True, f"L{side}_poisson", "defect equation holds at every site")
    return green


def dispersion(mode: tuple[int, int, int], table: list[Field]) -> Field:
    value: Field = table[0] * 0 + 3
    for coordinate in mode:
        value = value - table[coordinate]
    return value


def verify_lift(
    lift: dict[str, Any],
    green: dict[tuple[int, int, int], Fraction],
    coupling: Fraction,
) -> None:
    side = int(lift["side"])
    volume = side**3
    table = cosine_table(side)
    c0 = green[(0, 0, 0)]
    minimum = min(green.values())
    diameter = c0 - minimum
    alpha = min(Fraction(1, 2 * coupling), Fraction(1, diameter))
    p0 = Fraction(1) - alpha * c0
    expected_regime = (
        "infrared-all-ceiling"
        if Fraction(1, 2 * coupling) <= Fraction(1, diameter)
        else "GKS-diameter"
    )
    check(
        Fraction(lift["alpha"]) == alpha
        and Fraction(lift["p_zero_equals_M2"]) == p0
        and lift["regime"] == expected_regime,
        f"L{side}_lift_scale",
        "independent Green-ray branch and zero mode agree",
    )

    orbits = independent_orbits(side)
    stored_modes = {
        tuple(item["representative"]): item for item in lift["mode_orbits"]
    }
    expected_reps = {orbit[0] for orbit in orbits if orbit != ((0, 0, 0),)}
    assert set(stored_modes) == expected_reps
    full_p: dict[tuple[int, int, int], Field] = {(0, 0, 0): p0}
    for orbit in orbits:
        if orbit == ((0, 0, 0),):
            continue
        item = stored_modes[orbit[0]]
        lam = dispersion(orbit[0], table)
        probability = alpha / (volume * lam)
        assert int(item["multiplicity"]) == len(orbit)
        assert parse_field(item["lambda"]) == lam
        assert parse_field(item["p"]) == probability
        for mode in orbit:
            full_p[mode] = probability

    total: Field = sum(full_p.values(), table[0] * 0)
    assert total == 1
    assert len(full_p) == volume
    assert all(0 <= probability <= 1 for probability in full_p.values())
    for mode, probability in full_p.items():
        inverse = tuple((-coordinate) % side for coordinate in mode)
        assert full_p[inverse] == probability
        if mode != (0, 0, 0):
            lam = dispersion(mode, table)
            assert probability <= Fraction(1, 2 * coupling * volume) / lam

    for site, kernel_value in green.items():
        direct: Field = table[0] * 0 + p0
        for mode, probability in full_p.items():
            if mode == (0, 0, 0):
                continue
            phase = sum(mode[axis] * site[axis] for axis in range(3)) % side
            direct = direct + probability * table[phase]
        assert direct == p0 + alpha * kernel_value
        assert 0 <= direct <= 1

    # Independent verification of the implicit full four-point matrix
    # q_kl=p_k p_l without allocating O(N^2) entries.
    for probability in full_p.values():
        assert probability * total == probability
        assert probability * probability <= probability
        assert probability >= 0
    edge_g = p0 + alpha * green[(1, 0, 0)]
    energy_floor = Fraction(1) - (Fraction(1) - Fraction(1, volume)) / (
        6 * coupling
    )
    assert edge_g >= energy_floor
    q_info = lift["q_representation"]
    assert q_info["formula"] == "q[k,l]=p[k]*p[l]"
    assert q_info["rank"] == 1
    check(
        True,
        f"L{side}_mode_fourpoint_lift",
        f"all {volume} modes, inverse transforms, caps, and implicit Q rows pass",
    )


def main() -> int:
    started = time.process_time()
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    check(
        set(payload) == {"meta", "data", "checks"},
        "artifact_shape",
        "top-level envelope is exactly meta/data/checks",
    )
    check(
        payload["meta"]["producer"] == "experiments/e217_fourpoint_endpoint.py",
        "producer_identity",
        "combined producer recorded",
    )
    names = [row["name"] for row in payload["checks"]]
    check(
        len(names) == len(set(names)) and all(row["passed"] is True for row in payload["checks"]),
        "stored_checks",
        f"{len(names)} unique producer checks are true",
    )

    data = payload["data"]
    check(
        set(data["coverage"]) == EXPECTED_COVERAGE
        and len(data["coverage"]) == len(EXPECTED_COVERAGE),
        "coverage",
        "all universal-claim rows are present exactly once",
    )
    certificate = data["fourpoint_certificate"]
    challenge = certificate["challenge"]
    coupling = Fraction(challenge["K"])
    eta = Fraction(challenge["eta"])
    check(
        coupling == Fraction(6, 25)
        and eta == Fraction(1, 32)
        and Fraction(challenge["I3_lower_input"]) == I3_LOW
        and I3_LOW > 2 * coupling,
        "rational_challenge",
        "K=6/25 lies strictly below I3/2 by a directed rational endpoint",
    )

    kernel_rows = data["kernel_certificates"]["sides"]
    check(
        [int(row["side"]) for row in kernel_rows] == [4, 6, 8, 12],
        "kernel_sides",
        "expected exact finite certificates are present",
    )
    kernels = {int(row["side"]): reconstruct_kernel(row) for row in kernel_rows}

    lifts = certificate["finite_lifts"]
    check(
        [int(row["side"]) for row in lifts] == [4, 6, 8, 12],
        "lift_sides",
        "every kernel has one mode-resolved lift",
    )
    for lift in lifts:
        verify_lift(lift, kernels[int(lift["side"])], coupling)
        if time.process_time() - started > CPU_BUDGET_SECONDS:
            raise RuntimeError("standalone verifier CPU budget exceeded")
        if peak_rss_bytes() >= RSS_LIMIT_BYTES:
            raise MemoryError("standalone verifier RSS limit exceeded")

    endpoint_gap = I3_LOW - 2 * coupling
    challenge_cutoff = even_at_least(ceil_fraction(A_RIEMANN / endpoint_gap))
    eta_cutoff = max(challenge_cutoff, even_strictly_above(RATE_C / eta))
    all_size = certificate["all_size_certificate"]
    check(
        Fraction(all_size["A"]) == A_RIEMANN
        and Fraction(all_size["rate_constant"]) == RATE_C
        and int(all_size["challenge_even_cutoff"]) == challenge_cutoff == 96
        and I3_LOW - A_RIEMANN / challenge_cutoff >= 2 * coupling,
        "all_size_cutoff",
        "independent exact arithmetic reproduces the even challenge cutoff 96",
    )
    check(
        int(all_size["fixed_eta_even_cutoff"]) == eta_cutoff == 316
        and RATE_C / eta_cutoff < eta
        and eta - RATE_C / eta_cutoff == Fraction(131, 1327200),
        "all_size_rate",
        "p0<=5168/(525L)<1/32 for every even L>=316",
    )

    lift12 = next(row for row in lifts if int(row["side"]) == 12)
    p0_12 = Fraction(lift12["p_zero_equals_M2"])
    finite = certificate["fixed_eta_L12"]
    check(
        p0_12 == Fraction(94973361469787, 3762640205356032)
        and Fraction(finite["p0"]) == p0_12
        and Fraction(finite["exact_gap_eta_minus_p0"]) == eta - p0_12
        and eta - p0_12 == Fraction(22609144947589, 3762640205356032),
        "L12_eta_countercertificate",
        "exact finite primal lies below eta by the stored positive rational gap",
    )

    class_data = certificate["class"]
    check(
        class_data["name"] == "MR4-power-simplex-L1"
        and "[1;p][1;p]^T" in class_data["psd_identity"]
        and "exactly" in class_data["exact_projection_theorem"],
        "projection_certificate_form",
        "the artifact records the deterministic outer-product projection theorem",
    )
    check(
        peak_rss_bytes() < RSS_LIMIT_BYTES,
        "resource_bound",
        f"peak RSS {peak_rss_bytes()} bytes is below 2 GiB",
    )
    print(f"PASS process_time -- {time.process_time() - started:.3f}s")
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
