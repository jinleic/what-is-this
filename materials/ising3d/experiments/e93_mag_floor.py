#!/usr/bin/env python3
"""Exact uniform-in-L collapse of the finite-torus infrared/GKS LP floors.

Producer for results/bounds/mag_floor.json and proofs/mag_floor.md.

Objects (notation of proofs/kc_interval2.md section 2):
  T_L = (Z/LZ)^3, N = L^3, lambda(k) = 3 - sum_i cos k_i,
  C_L(z) = N^-1 sum_{k != 0} cos(k.z)/lambda(k)   (zero-mode-removed torus Green),
  D_L = C_L(0) - min_z C_L(z),  kappa_L = D_L/2,  u_L = (-min_z C_L)/D_L,
  floor(K,L) = 1 - S_*(K,L) = exact optimum of the finite-torus
               infrared/GKS/cap/energy LP maximising S = N^-1 sum_{k!=0} Ghat(k).

Exact results produced here (all arithmetic Fraction or Q(sqrt d), d in {2,3,5}):
  T1  floor(K,L) > 0 strictly for every K > 0 and every L >= 2;
  T2  floor(K,L) = 1 - C_L(0)/(2K) exactly whenever 2K >= D_L (closed form),
      floor(K,L) <= u_L whenever 2K <= D_L (explicit feasible Green profile);
  T3  explicit rational constants (pi < 355/113 etc.) giving
      min_z C_L >= -2 sqrt(J2(L)/N), J2(L) = N^-1 sum_{k!=0} lambda(k)^-2,
      J2(L) <= L*Sigma4/64 with Sigma4 <= 4 pi^2 + pi^4/45,
      C_L(0) >= I3 - (c_Delta+1)/L, and hence
      floor(K,L) <= c*/L,  c* = 5168/525,  for all L >= 4 with
      2K <= I3 - (c_Delta+1)/L;
  T4  falsification of the uniform LP-floor lemma for every K <= I3/2;
  T5  meta-negative: every relaxation of the audited constraint system has an
      even smaller floor, so no momentum-averaged/paired modification can help.

New exact data: torus Green tables at L = 8 (Q(sqrt2)), 10 (Q(sqrt5)),
12 (Q(sqrt3)); exact LP primal/dual certificates at L = 8 in Q(sqrt2);
independent recomputation of the L = 4, 6 LP floors of e85_kc_interval2.

Run:  timeout 3600 .venv/bin/python experiments/e93_mag_floor.py
"""

from __future__ import annotations

import itertools
import json
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Sequence

REPO = Path(__file__).resolve().parents[1]
E85_JSON = REPO / "results" / "bounds" / "kc_interval2.json"
UPPER_JSON = REPO / "results" / "bounds" / "upper_infrared.json"
OUT_JSON = REPO / "results" / "bounds" / "mag_floor.json"

SIDES_GREEN = (4, 6, 8, 10, 12)
SIDES_LP = (4, 6, 8)
K_POINTS = (Fraction(1, 5), Fraction(21, 100), Fraction(11, 50), Fraction(6, 25))

# Certified rational upper bounds for the transcendental constants that enter
# the explicit-rate lemmas.  355/113 > pi and 99/70 > sqrt(2) are classical
# exact rational bounds (Zu Chongzhi; (99/70)^2 = 9801/4900 > 2).
PI_UB = Fraction(355, 113)
SQRT2_UB = Fraction(99, 70)
# Certified lower endpoint of I3, from proofs/upper_infrared.md section 2
# (mpmath.iv at dps=100, outward rounded).  Used only through the inequality
# I3 >= I3_LOW.
I3_LOW = Fraction(
    "0.505462019717326006052004053227140259985129014817420892188993487886028773451173816800537247069896037925625"
)
I3_HALF_UPPER_DECIMAL = "0.2527310098586630030260020266135701299926"


# ---------------------------------------------------------------------------
# Exact quadratic-field arithmetic: value a + b*sqrt(d), a,b Fractions.
# ---------------------------------------------------------------------------


class Quad:
    __slots__ = ("a", "b", "d")

    def __init__(self, a=0, b=0, d=2):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.d = d

    @staticmethod
    def _coerce(other: object, d: int) -> "Quad":
        if isinstance(other, Quad):
            if other.d != d:
                raise ValueError("mixed sqrt fields")
            return other
        if isinstance(other, (int, Fraction)):
            return Quad(other, 0, d)
        return NotImplemented

    def __add__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(self.a + o.a, self.b + o.b, self.d)

    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(self.a - o.a, self.b - o.b, self.d)

    def __rsub__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return o - self

    def __mul__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(
            self.a * o.a + self.b * o.b * o.d,
            self.a * o.b + self.b * o.a,
            self.d,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        norm = o.a * o.a - o.b * o.b * o.d
        if norm == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt d)")
        num = self * Quad(o.a, -o.b, o.d)
        return Quad(num.a / norm, num.b / norm, self.d)

    def __rtruediv__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return o / self

    def __neg__(self):
        return Quad(-self.a, -self.b, self.d)

    def __pos__(self):
        return self

    def __eq__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return self.a == o.a and self.b == o.b

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def sign(self) -> int:
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        ratio = self.a / self.b
        if ratio > 0:
            return 1 if self.b > 0 else -1
        # ratio < 0: a + b*sqrt(d) = b*(ratio + sqrt(d)); ratio + sqrt(d) > 0
        # iff ratio^2 < d (ratio negative; equality impossible for nonsquare d).
        positive = ratio * ratio < self.d
        if self.b > 0:
            return 1 if positive else -1
        return -1 if positive else 1

    def __lt__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() < 0

    def __le__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() <= 0

    def __gt__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() > 0

    def __ge__(self, other):
        o = self._coerce(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() >= 0

    def is_rational(self) -> bool:
        return self.b == 0

    def rational(self) -> Fraction:
        if self.b != 0:
            raise ValueError(f"not rational: {self}")
        return self.a

    def __repr__(self) -> str:
        return f"Quad({self.a},{self.b},sqrt{self.d})"


def field_text(value: Fraction | Quad) -> str:
    if isinstance(value, Quad):
        return f"Q[{value.a}|{value.b}|{value.d}]"
    return str(value)


def parse_field(text: str) -> Fraction | Quad:
    if text.startswith("Q["):
        body = text[2:-1]
        a, b, d = body.split("|")
        return Quad(Fraction(a), Fraction(b), int(d))
    return Fraction(text)


def fraction_text(value: Fraction | Quad) -> str:
    if isinstance(value, Quad):
        return str(value.rational())
    return str(value)


# ---------------------------------------------------------------------------
# Exact cosine tables.  cos(2*pi*j/L) for the quadratic-cosine sides.
# ---------------------------------------------------------------------------

FIELD_DEGREE: dict[int, int | None] = {4: None, 6: None, 8: 2, 10: 5, 12: 3}


def cos_table(side: int) -> tuple[list[Fraction | Quad], Fraction | Quad, int | None]:
    """Exact table cos(2 pi j / L), j = 0..L-1, built from cos(2 pi / L)."""

    d = FIELD_DEGREE[side]
    if side == 4:
        base: Fraction | Quad = Fraction(0)
    elif side == 6:
        base = Fraction(1, 2)
    elif side == 8:
        base = Quad(0, Fraction(1, 2), 2)  # sqrt(2)/2
    elif side == 10:
        base = Quad(Fraction(1, 4), Fraction(1, 4), 5)  # (1+sqrt5)/4
    elif side == 12:
        base = Quad(0, Fraction(1, 2), 3)  # sqrt(3)/2
    else:
        raise ValueError(f"no exact cosine table for side {side}")

    # Algebraic certification of the base value.
    if side == 4:
        assert base == 0
    elif side == 6:
        assert 2 * base == 1
    elif side == 8:
        assert 2 * base * base == 1
    elif side == 10:
        assert 4 * base * base - 2 * base - 1 == 0 and base > 0
    else:
        assert 2 * base * 2 * base == 3 and base > 0

    table: list[Fraction | Quad] = [base - base + 1, base]
    for j in range(1, side):
        table.append(2 * base * table[j] - table[j - 1])
    table = table[: side]  # table[j] for j in 0..side-1
    # Closure controls: periodicity and the exact midpoint value -1.
    wrap = 2 * base * table[side - 1] - table[side - 2]
    assert wrap == table[0]
    assert table[side // 2] == -1
    return table, base, d


# ---------------------------------------------------------------------------
# Signed-permutation orbits (hyperoctahedral B3), as in e85_kc_interval2.
# ---------------------------------------------------------------------------


def signed_permutation_orbit(point: tuple[int, int, int], side: int):
    images = set()
    for permutation in itertools.permutations(range(3)):
        reordered = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=3):
            images.add(
                tuple(
                    (sign * coordinate) % side
                    for sign, coordinate in zip(signs, reordered, strict=True)
                )
            )
    return tuple(sorted(images))

def torus_orbits(side: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    remaining = set(itertools.product(range(side), repeat=3))
    orbits: list[tuple[tuple[int, int, int], ...]] = []

    while remaining:
        point = min(remaining)
        orbit = signed_permutation_orbit(point, side)
        remaining.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


# ---------------------------------------------------------------------------
# Exact torus Green data per side.
# ---------------------------------------------------------------------------


def torus_green(side: int):
    """Exact C_L on every site, via direct mode sums on orbit representatives.

    Returns (green_by_site, orbits, checks).  C_L is invariant under signed
    permutations (lambda and the mode set are), so per-orbit computation
    determines every site.
    """

    table, base, d = cos_table(side)
    volume = side**3
    orbits = torus_orbits(side)
    modes = [m for m in itertools.product(range(side), repeat=3) if m != (0, 0, 0)]
    lam_by_mode: list[Fraction | Quad] = []
    for mode in modes:
        lam = table[0] * 0 + 3
        for axis in range(3):
            lam = lam - table[mode[axis]]
        lam_by_mode.append(lam)
    if any(not (lam > 0) for lam in lam_by_mode):
        raise AssertionError("nonzero mode with nonpositive dispersion")

    rep_values: dict[tuple[int, int, int], Fraction | Quad] = {}
    for orbit in orbits:
        z = orbit[0]
        total = table[0] * 0
        for mode, lam in zip(modes, lam_by_mode, strict=True):
            phase = sum(mode[axis] * z[axis] for axis in range(3)) % side
            total = total + table[phase] / lam
        rep_values[z] = total / volume

    green: dict[tuple[int, int, int], Fraction | Quad] = {}
    for orbit in orbits:
        for site in orbit:
            green[site] = rep_values[orbit[0]]

    # Controls: rational values (Galois invariance), zero spatial sum, defect
    # equation at every orbit representative.
    for value in rep_values.values():
        if isinstance(value, Quad):
            if not value.is_rational():
                raise AssertionError(f"C_L value not rational: {value}")
    zero_sum = table[0] * 0
    for orbit in orbits:
        zero_sum = zero_sum + len(orbit) * rep_values[orbit[0]]
    if zero_sum != 0:
        raise AssertionError("spatial sum of C_L is not zero")
    for orbit in orbits:
        z = orbit[0]
        lap = rep_values[z] * 3
        for axis in range(3):
            for step in (1, -1):
                neighbour = list(z)
                neighbour[axis] = (neighbour[axis] + step) % side
                lap = lap - Fraction(1, 2) * green[tuple(neighbour)]
        target = (table[0] * 0 + 1) if z == (0, 0, 0) else table[0] * 0
        target = target - Fraction(1, volume)
        if lap != target:
            raise AssertionError(f"defect equation failed at {z}")
    return green, orbits, {"zero_sum": True, "defect": True}


def j2_exact(side: int, orbits) -> Fraction:
    """J2(L) = N^-1 sum_{k != 0} lambda(k)^-2, exactly (Galois invariant)."""

    table, _, _ = cos_table(side)
    volume = side**3
    total: Fraction | Quad = table[0] * 0
    for orbit in orbits:
        if orbit == ((0, 0, 0),):
            continue
        k = orbit[0]
        lam = table[0] * 0 + 3
        for axis in range(3):
            lam = lam - table[k[axis]]
        total = total + (table[0] * 0 + len(orbit)) / (lam * lam)
    if isinstance(total, Quad):
        total = total.rational()  # Galois-invariant sum must be rational
    return total / volume


# ---------------------------------------------------------------------------
# Generic exact simplex over Fraction / Quad (Bland-safe Dantzig rule).
# ---------------------------------------------------------------------------


def solve_square(
    matrix: Sequence[Sequence[Fraction | Quad]], rhs: Sequence[Fraction | Quad]
) -> list[Fraction | Quad]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix) or len(rhs) != size:
        raise ValueError("expected a nonempty square system")
    augmented = [list(row) + [rhs[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column] != 0), None)
        if pivot is None:
            raise ValueError("singular matrix")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor != 0:
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(augmented[row], augmented[column], strict=True)
                ]
    return [row[-1] for row in augmented]


def simplex_max(
    a: Sequence[Sequence[Fraction | Quad]],
    b: Sequence[Fraction | Quad],
    c: Sequence[Fraction | Quad],
    iteration_cap: int = 50000,
):
    """Exact primal simplex for max c.x, A x <= b, x >= 0, with dual.

    Entering rule: largest positive reduced cost (Dantzig), falling back to
    Bland's rule after 300 consecutive non-improving iterations, which
    guarantees finite termination.
    """

    row_count = len(a)
    variable_count = len(c)
    if not row_count or any(len(row) != variable_count for row in a):
        raise ValueError("malformed constraint matrix")
    if len(b) != row_count or any(not (value >= 0) for value in b):
        raise ValueError("initial slack basis requires nonnegative RHS")
    zero = b[0] * 0
    one = zero + 1
    full_a = [
        list(row) + [one if column == row_index else zero for column in range(row_count)]
        for row_index, row in enumerate(a)
    ]
    tableau = [list(row) + [b[index]] for index, row in enumerate(full_a)]
    full_c = list(c) + [zero for _ in range(row_count)]
    basis = [variable_count + row for row in range(row_count)]
    stale = 0
    iterations = 0
    objective_prev = zero
    while True:
        iterations += 1
        if iterations > iteration_cap:
            raise RuntimeError("simplex iteration cap reached")
        bland = stale > 300
        c_basis = [full_c[index] for index in basis]
        reduced = [
            full_c[column]
            - sum(c_basis[row] * tableau[row][column] for row in range(row_count))
            for column in range(variable_count + row_count)
        ]
        entering = None
        if bland:
            entering = next((column for column, value in enumerate(reduced) if value > zero), None)
        else:
            best = zero
            for column, value in enumerate(reduced):
                if value > best:
                    best = value
                    entering = column
        if entering is None:
            break
        ratios = [
            (tableau[row][-1] / tableau[row][entering], basis[row], row)
            for row in range(row_count)
            if tableau[row][entering] > zero
        ]
        if not ratios:
            raise ValueError("unbounded LP")
        _, _, leaving = min(ratios, key=lambda item: (item[0], item[1], item[2]))
        divisor = tableau[leaving][entering]
        tableau[leaving] = [value / divisor for value in tableau[leaving]]
        for row in range(row_count):
            if row == leaving:
                continue
            factor = tableau[row][entering]
            if factor != 0:
                tableau[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(tableau[row], tableau[leaving], strict=True)
                ]
        basis[leaving] = entering
        # Cheap progress probe: value of the basic solution objective.
        probe = zero
        for row, index in enumerate(basis):
            probe = probe + full_c[index] * tableau[row][-1]
        if probe == objective_prev:
            stale += 1
        else:
            stale = 0
        objective_prev = probe

    primal_full = [zero for _ in range(variable_count + row_count)]
    for row, column in enumerate(basis):
        primal_full[column] = tableau[row][-1]
    primal = primal_full[:variable_count]
    objective = sum(value * coefficient for value, coefficient in zip(primal, c, strict=True))

    basis_matrix = [[full_a[row][basis[column]] for column in range(row_count)] for row in range(row_count)]
    c_basis = [full_c[column] for column in basis]
    transpose = [[basis_matrix[column][row] for column in range(row_count)] for row in range(row_count)]
    dual = solve_square(transpose, c_basis)
    dual_objective = sum(value * bound for value, bound in zip(dual, b, strict=True))
    primal_feasible = all(value >= zero for value in primal) and all(
        sum(coefficient * value for coefficient, value in zip(row, primal, strict=True)) <= bound
        for row, bound in zip(a, b, strict=True)
    )
    dual_feasible = all(value >= zero for value in dual) and all(
        sum(dual[row] * full_a[row][column] for row in range(row_count)) >= full_c[column]
        for column in range(variable_count + row_count)
    )
    if not (primal_feasible and dual_feasible and objective == dual_objective):
        raise AssertionError("exact primal/dual LP certificate failed")
    return {
        "primal": primal,
        "dual": dual,
        "objective": objective,
        "basis": basis,
        "primal_feasible": primal_feasible,
        "dual_feasible": dual_feasible,
        "iterations": iterations,
    }


# ---------------------------------------------------------------------------
# The finite-torus LP (orbit-symmetrised relaxation, as in e85 route 2).
# ---------------------------------------------------------------------------


def lp_problem(side: int, coupling: Fraction):
    table, _, _ = cos_table(side)
    volume = side**3
    zero = table[0] * 0
    one = zero + 1
    all_orbits = torus_orbits(side)
    mode_orbits = tuple(orbit for orbit in all_orbits if orbit != ((0, 0, 0),))
    site_orbits = tuple(orbit for orbit in all_orbits if orbit != ((0, 0, 0),))
    multiplicities = [len(orbit) for orbit in mode_orbits]
    dispersions: list[Fraction | Quad] = []
    for orbit in mode_orbits:
        lam = zero + 3
        for axis in range(3):
            lam = lam - table[orbit[0][axis]]
        dispersions.append(lam)
    if any(not (value > 0) for value in dispersions):
        raise AssertionError("nonpositive dispersion")

    phase_sums: list[list[Fraction | Quad]] = []
    for site_orbit in site_orbits:
        site = site_orbit[0]
        row = []
        for mode_orbit in mode_orbits:
            total = zero
            for mode in mode_orbit:
                phase = sum(mode[axis] * site[axis] for axis in range(3)) % side
                total = total + table[phase]
            row.append(total)
        phase_sums.append(row)

    rows: list[list[Fraction | Quad]] = []
    rhs: list[Fraction | Quad] = []
    labels: list[str] = []
    for index, dispersion in enumerate(dispersions):
        row = [zero for _ in mode_orbits]
        row[index] = one
        rows.append(row)
        rhs.append(one / (2 * coupling * dispersion))
        labels.append(f"ceiling_mode_{index}")
    for orbit, phase in zip(site_orbits, phase_sums, strict=True):
        site = orbit[0]
        gap = [one * m - coefficient for m, coefficient in zip(multiplicities, phase, strict=True)]
        rows.append(list(gap))
        rhs.append(one * volume)
        labels.append("gks_nonnegative_" + "_".join(map(str, site)))
        rows.append([-value for value in gap])
        rhs.append(zero)
        labels.append("spin_upper_" + "_".join(map(str, site)))
    rows.append([one * value for value in multiplicities])
    rhs.append(one * volume)
    labels.append("zero_mode_nonnegative")
    neighbour_index = next(index for index, orbit in enumerate(site_orbits) if (1, 0, 0) in orbit)
    gap_e = [
        one * m - coefficient
        for m, coefficient in zip(multiplicities, phase_sums[neighbour_index], strict=True)
    ]
    rows.append(list(gap_e))
    rhs.append(one * (volume - 1) / (6 * coupling))
    labels.append("energy_nearest_neighbour")

    return {
        "side": side,
        "coupling": coupling,
        "volume": volume,
        "mode_orbits": mode_orbits,
        "site_orbits": site_orbits,
        "multiplicities": multiplicities,
        "dispersions": dispersions,
        "phase_sums": phase_sums,
        "a": rows,
        "b": rhs,
        "c": [one * value for value in multiplicities],
        "labels": labels,
    }


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    if not passed:
        raise AssertionError(f"FAILED check {name}: {detail}")
    checks.append({"name": name, "passed": True, "detail": detail})


def quad_self_test(checks: list[dict[str, object]]) -> None:
    x = Quad(1, 1, 2)
    y = Quad(3, -2, 2)
    assert (x * x) == Quad(3, 2, 2)
    assert (x + y) == Quad(4, -1, 2)
    assert (x - y) == Quad(-2, 3, 2)
    assert (x / y) * y == x
    assert (Quad(-1, 1, 2)).sign() == 1  # sqrt2 - 1 > 0
    assert (Quad(1, -1, 2)).sign() == -1
    assert (Quad(-3, 1, 2)).sign() == -1  # sqrt2 < 3
    assert (Quad(3, 1, 2)).sign() == 1
    assert (Quad(0, 1, 5)) * Quad(0, 1, 5) == 5
    z = Quad(Fraction(1, 3), Fraction(-2, 7), 3)
    assert (z / Quad(0, 1, 3)) * Quad(0, 1, 3) == z
    _check(checks, "quad_field_self_test", True, "exact Q(sqrt d) arithmetic controls")


def ceil_fraction(value: Fraction) -> int:
    """Exact ceiling for a nonnegative Fraction."""

    return -((-value.numerator) // value.denominator)


def main() -> None:
    started = time.time()
    checks: list[dict[str, object]] = []
    quad_self_test(checks)

    e85 = json.loads(E85_JSON.read_text(encoding="utf-8"))
    upper = json.loads(UPPER_JSON.read_text(encoding="utf-8"))
    e85_floors: dict[tuple[int, str], str] = {}
    for record in e85["data"]["route2_finite_volume_floor"]["records"]:
        for cert in record["finite_floor_certificates"]:
            e85_floors[(record["side"], cert["coupling"])] = cert["magnetization_square_floor"]

    # ------------------------------------------------------------------
    # Explicit constants for the uniform-rate lemmas.
    # ------------------------------------------------------------------
    pi2_ub = PI_UB * PI_UB
    pi4_ub = pi2_ub * pi2_ub
    sigma4_ub = 4 * pi2_ub + pi4_ub / 45  # >= Sigma4 = sum_{j != 0} |j|^-4
    sqrt_sigma4_ub = Fraction(323, 50)
    _check(
        checks,
        "constants_sigma4_sqrt_bound",
        sqrt_sigma4_ub * sqrt_sigma4_ub >= sigma4_ub,
        "Sigma4 <= 4 pi^2 + pi^4/45 <= (323/50)^2 exactly",
    )
    c_star = Fraction(32, 21) * sqrt_sigma4_ub / 1  # (3/2)(64/63) * sqrt(Sigma4)
    # I3 - C_mu(0) is split into {lambda >= mu^2} and {lambda < mu^2}.
    # The elementary bounds lambda <= |k|^2/2 and lambda >= 2|k|^2/pi^2
    # give pi^2*mu/(8 sqrt2) + pi*mu/(4 sqrt2), bounded as below.
    c_delta = SQRT2_UB * (pi2_ub + 2 * PI_UB) / 16
    c_delta_plus_1 = c_delta + 1
    _check(
        checks,
        "constants_c_star_form",
        c_star == Fraction(5168, 525),
        "c* = (3/2)(N/(N-1))_at_L=4 * sqrt(Sigma4ub) = 5168/525 exactly",
    )
    _check(
        checks,
        "constants_endpoint_dominated_by_c_star",
        c_star >= c_delta_plus_1 / I3_LOW,
        "at K=I3/2, the D_L<=I3 case is also bounded by c*/L",
    )

    # ------------------------------------------------------------------
    # Per-side exact Green data, D_L, kappa_L, u_L, J2(L), lemma controls.
    # ------------------------------------------------------------------
    side_records: dict[str, dict[str, object]] = {}
    for side in SIDES_GREEN:
        t0 = time.time()
        green, orbits, green_checks = torus_green(side)
        volume = side**3
        origin = green[(0, 0, 0)]
        minimum = min(green.values())
        minimum_sites = sorted(site for site, value in green.items() if value == minimum)
        d_value = origin - minimum
        kappa = d_value / 2
        u_value = (-minimum) / d_value
        j2 = j2_exact(side, orbits)
        # Lemma M1 control: (min C_L)^2 <= 4 J2 / N  (exact Fractions).
        m1_ok = minimum * minimum <= Fraction(4) * j2 / volume
        # Lemma M2 control: J2(L) <= L * Sigma4ub / 64.
        m2_ok = j2 <= Fraction(side) * sigma4_ub / 64
        # Lemma M3 control: C_L(0) >= I3_LOW - (c_delta+1)/L.
        m3_ok = origin >= I3_LOW - c_delta_plus_1 / side
        # u_L <= c*/L.
        u_rate_ok = u_value <= c_star / side
        # kappa_L < I3/2, decided through the certified rational lower endpoint.
        kappa_ok = d_value < I3_LOW
        _check(
            checks,
            f"green_lemmas_L{side}",
            m1_ok and m2_ok and m3_ok and u_rate_ok and kappa_ok,
            "M1/M2/M3 controls, u_L <= c*/L, and D_L < I3 all hold exactly",
        )
        orbit_table = [
            {
                "representative": list(orbit[0]),
                "multiplicity": len(orbit),
                "C_L": fraction_text(green[orbit[0]]),
            }
            for orbit in orbits
        ]
        side_records[str(side)] = {
            "side": side,
            "sqrt_field": FIELD_DEGREE[side],
            "volume": volume,
            "C_L_at_origin": fraction_text(origin),
            "minimum_C_L": fraction_text(minimum),
            "minimum_sites": [list(site) for site in minimum_sites[:8]],
            "D_L": fraction_text(d_value),
            "kappa_L": fraction_text(kappa),
            "u_L": fraction_text(u_value),
            "u_L_times_L": fraction_text(u_value * side),
            "J2_L": str(j2),
            "M1_control_4J2_over_N": str(Fraction(4) * j2 / volume),
            "green_orbits": orbit_table,
            "green_checks": green_checks,
            "lp": [],
            "wall_seconds": time.time() - t0,
        }
        print(
            f"[side {side}] C0={fraction_text(origin)} min={fraction_text(minimum)} "
            f"kappa={fraction_text(kappa)} u={fraction_text(u_value)} "
            f"J2={j2} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    # ------------------------------------------------------------------
    # LP certificates at sides 4, 6 (Fraction, must match e85) and 8 (Q(sqrt2)).
    # ------------------------------------------------------------------
    for side in SIDES_LP:
        record = side_records[str(side)]
        d_value = Fraction(record["D_L"])
        kappa = Fraction(record["kappa_L"])
        origin = Fraction(record["C_L_at_origin"])
        for coupling in K_POINTS:
            t0 = time.time()
            problem = lp_problem(side, coupling)
            solved = simplex_max(problem["a"], problem["b"], problem["c"])
            objective = solved["objective"]
            s_star = objective / problem["volume"]
            floor = s_star * 0 + 1 - s_star
            wall = time.time() - t0
            regime = "lp"
            closed_form: str | None = None
            if 2 * coupling >= d_value:
                regime = "closed_form"
                expected_floor = Fraction(1) - origin / (2 * coupling)
                closed_form = str(expected_floor)
                if floor != expected_floor:
                    raise AssertionError(
                        f"closed form mismatch L={side} K={coupling}: "
                        f"{expected_floor} vs {floor}"
                    )
            e85_key = (side, str(coupling))
            if e85_key in e85_floors and fraction_text(floor) != e85_floors[e85_key]:
                raise AssertionError(
                    f"e85 floor mismatch L={side} K={coupling}: "
                    f"{e85_floors[e85_key]} vs {floor}"
                )
            # T2 profile certificate: floor <= u_L whenever 2K <= D_L.
            if 2 * coupling <= d_value:
                u_value = Fraction(record["u_L"])
                if not floor <= u_value:
                    raise AssertionError(f"floor exceeds u_L at L={side} K={coupling}")
            # LP floor is strictly positive (T1).
            if floor <= 0:
                raise AssertionError(f"nonpositive floor at L={side} K={coupling}")
            certificate = {
                "coupling": str(coupling),
                "regime": regime,
                "closed_form_floor": closed_form,
                "matches_e85": e85_key in e85_floors,
                "objective_sum_nonzero_modes": field_text(objective),
                "S_star": field_text(s_star),
                "magnetization_square_floor": field_text(floor),
                "iterations": solved["iterations"],
                "mode_orbits": [
                    {
                        "representative": list(orbit[0]),
                        "multiplicity": multiplicity,
                        "lambda": field_text(dispersion),
                    }
                    for orbit, multiplicity, dispersion in zip(
                        problem["mode_orbits"],
                        problem["multiplicities"],
                        problem["dispersions"],
                        strict=True,
                    )
                ],
                "constraint_labels": list(problem["labels"]),
                "primal_x": [field_text(value) for value in solved["primal"]],
                "dual_y": [field_text(value) for value in solved["dual"]],
                "basis": list(solved["basis"]),
                "primal_dual_equal": True,
                "wall_seconds": wall,
            }
            record["lp"].append(certificate)
            print(
                f"[LP L={side} K={coupling}] floor={floor} regime={regime} "
                f"iters={solved['iterations']} ({wall:.1f}s)",
                flush=True,
            )

    # Explicit feasible-profile certificate (T2): x_k = 1/(D_L lambda_k) is
    # feasible for the FULL (non-symmetrised) LP whenever 2K <= D_L; verify at
    # every site for every side, and record G(z) on orbit representatives.
    profile_records = []
    for side in SIDES_GREEN:
        green, orbits, _ = torus_green(side)
        volume = side**3
        d_value = Fraction(fraction_text(green[(0, 0, 0)] - min(green.values())))
        for coupling in K_POINTS:
            if 2 * coupling > d_value:
                continue
            origin = Fraction(fraction_text(green[(0, 0, 0)]))
            ok_sites = True
            for site, value in green.items():
                g_z = Fraction(1) - (origin - Fraction(fraction_text(value))) / d_value
                if not (Fraction(0) <= g_z <= Fraction(1)):
                    ok_sites = False
                    break
            if not ok_sites:
                raise AssertionError(f"Green profile infeasible at L={side} K={coupling}")
            profile_records.append(
                {
                    "side": side,
                    "coupling": str(coupling),
                    "c_coefficient": str(Fraction(1, 1) / d_value),
                    "S_profile": str(origin / d_value),
                    "floor_upper_bound_u_L": str(Fraction(1) - origin / d_value),
                    "all_sites_verified": True,
                }
            )
    _check(
        checks,
        "green_profile_feasibility_all_sides",
        len(profile_records) >= 8,
        "x_k = 1/(D_L lambda_k) verified feasible (all sites, GKS 0<=G<=1, "
        "ceilings, cap) for every stored (side, K) with 2K <= D_L",
    )

    # Falsification table: exact floor sequences along L at each stored K.
    falsification: dict[str, object] = {
        "statement": (
            "[FALSIFIED] For every positive coupling K with 2K <= I3 there is no "
            "m > 0 and no L0 such that the exact finite-torus LP floor satisfies "
            "floor(K,L) >= m for every even L >= L0.  If 2K < I3, then "
            "floor(K,L) <= c*/L for all L >= max(4, ceil((c_delta+1)/(I3-2K))).  "
            "At the endpoint 2K=I3, use the profile c=min(I3^-1,D_L^-1): when "
            "D_L<=I3, M3 gives 1-C_L(0)/I3 <= (c_delta+1)/(I3*L) <= c*/L; when "
            "D_L>=I3, the u_L bound gives floor<=c*/L."
        ),
        "c_star_over_L": {
            str(side): str(c_star / side) for side in SIDES_GREEN
        },
        "floor_sequences": {},
        "u_L_sequence": {str(side): side_records[str(side)]["u_L"] for side in SIDES_GREEN},
        "kappa_sequence": {str(side): side_records[str(side)]["kappa_L"] for side in SIDES_GREEN},
    }
    for coupling in K_POINTS:
        sequence = []
        for side in SIDES_LP:
            cert = next(
                item for item in side_records[str(side)]["lp"] if item["coupling"] == str(coupling)
            )
            sequence.append({"side": side, "floor": cert["magnetization_square_floor"]})
        decreasing = all(
            parse_field(sequence[i]["floor"]) > parse_field(sequence[i + 1]["floor"])
            for i in range(len(sequence) - 1)
        )
        falsification["floor_sequences"][str(coupling)] = {
            "floors": sequence,
            "strictly_decreasing_in_L": decreasing,
        }
        _check(
            checks,
            f"floor_decreasing_K{coupling.numerator}_{coupling.denominator}",
            decreasing,
            "exact LP floor strictly decreases from L=4 to L=6 to L=8 at this K",
        )
    dyadic_f4 = parse_field(
        next(item for item in side_records["4"]["lp"] if item["coupling"] == "1/5")[
            "magnetization_square_floor"
        ]
    )
    dyadic_f8 = parse_field(
        next(item for item in side_records["8"]["lp"] if item["coupling"] == "1/5")[
            "magnetization_square_floor"
        ]
    )
    dyadic_difference = dyadic_f4 - dyadic_f8
    expected_dyadic_difference = Quad(Fraction(28697, 626688), Fraction(1, 816), 2)
    _check(
        checks,
        "dyadic_lp_floor_countercertificate",
        dyadic_difference == expected_dyadic_difference and dyadic_difference > 0,
        "F(1/5,4)-F(1/5,8)=28697/626688+sqrt(2)/816 > 0",
    )
    falsification["dyadic_countercertificate"] = {
        "statement": (
            "[FALSIFIED] Candidate nondecreasing dyadic LP-floor relation "
            "floor(K,2L)>=floor(K,L) fails at K=1/5, L=4."
        ),
        "coupling_K": "1/5",
        "smaller_side_L": 4,
        "larger_side_2L": 8,
        "floor_difference_F_L_minus_F_2L": field_text(dyadic_difference),
        "strictly_positive": True,
    }
    kappa_increasing = all(
        Fraction(side_records[str(a)]["kappa_L"]) < Fraction(side_records[str(b)]["kappa_L"])
        for a, b in zip(SIDES_GREEN, SIDES_GREEN[1:])
    )
    _check(
        checks,
        "kappa_increasing_on_checked_sides",
        kappa_increasing,
        "kappa_L = (C_L(0)-min C_L)/2 strictly increases through L=4,6,8,10,12",
    )
    # Route-(3) conditional, made concrete on the pre-existing rational LP grid.
    # For fixed K,m, H_{K,m}(L) := [floor(K,L) >= m] is decidable from the
    # exact certificate.  If it held for every sufficiently large even L, the
    # finite-torus infrared mechanism would imply K_c <= K.  The three exact
    # checks below pass, but the rate theorem proves H false beyond `first_even`.
    h_coupling = Fraction(11, 50)
    h_floor = Fraction(1, 100)
    h_rows = []
    for side in SIDES_LP:
        cert = next(item for item in side_records[str(side)]["lp"] if item["coupling"] == str(h_coupling))
        value = parse_field(cert["magnetization_square_floor"])
        if not value >= h_floor:
            raise AssertionError(f"H( L={side} ) unexpectedly fails")
        h_rows.append({"side": side, "floor": cert["magnetization_square_floor"], "H_holds": True})
    h_condition_start = ceil_fraction(c_delta_plus_1 / (I3_LOW - 2 * h_coupling))
    h_rate_start = (c_star / h_floor).numerator // (c_star / h_floor).denominator + 1
    h_first_even = max(4, h_condition_start, h_rate_start)
    if h_first_even % 2:
        h_first_even += 1
    _check(
        checks,
        "conditional_H_11_50_1_100_countercertificate",
        2 * h_coupling <= I3_LOW - c_delta_plus_1 / h_first_even
        and c_star / h_first_even < h_floor,
        "for every even L >= first_even, the rate theorem forces H(L) false",
    )
    conditional_h = {
        "tag": "[CONDITIONAL] then [FALSIFIED]",
        "predicate": "H_{K,m}(L): exact LP floor(K,L) >= m",
        "coupling_K": str(h_coupling),
        "floor_m": str(h_floor),
        "implication": (
            "If H_{11/50,1/100}(L) held for every sufficiently large even L, then "
            "the uniform finite-torus magnetization floor would imply K_c <= 11/50."
        ),
        "exact_verified_sides": h_rows,
        "rate_countercertificate": {
            "first_even_L": h_first_even,
            "condition_start_L": h_condition_start,
            "strict_rate_start_L": h_rate_start,
            "conclusion": (
                "For every even L >= first_even_L, 2K <= I3-(c_delta+1)/L and "
                "floor(K,L) <= c*/L < 1/100, so H_{11/50,1/100}(L) is false."
            ),
        },
    }

    data = {
        "model": "nearest-neighbour ferromagnetic Ising model on Z^3, K = beta*J",
        "objects": (
            "floor(K,L) = 1 - S_*(K,L), S_* the maximum of the exact finite-torus "
            "LP over the audited infrared/GKS/cap/energy constraint system of "
            "proofs/kc_interval2.md section 2. Cubic-group averaging makes the "
            "stored orbit form an exact reduction of that LP; the LP itself remains "
            "a relaxation of the physical correlation set. The profile certificates "
            "below are also checked at every site of the unreduced system."
        ),
        "headline": (
            "[THEOREM] The exact finite-torus LP floors are strictly positive at "
            "every finite L but collapse uniformly: floor(K,L) <= u_L = "
            "(-min_z C_L)/D_L <= c*/L with c* = 5168/525 whenever 2K <= D_L, and "
            "2K <= D_L holds for all sufficiently large L whenever 0 < 2K < I3. "
            "The uniform LP-floor lemma needed to move the upper endpoint "
            "K_c <= I3/2 is therefore false for every positive K <= I3/2 within "
            "the audited constraint class."
        ),
        "constants": {
            "pi_upper": str(PI_UB),
            "sqrt2_upper": str(SQRT2_UB),
            "sigma4_upper": str(sigma4_ub),
            "sqrt_sigma4_upper": str(sqrt_sigma4_ub),
            "c_star": str(c_star),
            "c_delta": str(c_delta),
            "c_delta_plus_1": str(c_delta_plus_1),
            "I3_lower_rational": str(I3_LOW),
            "I3_half_upper_decimal": I3_HALF_UPPER_DECIMAL,
            "L1_formula": "ceil((c_delta_plus_1)/(I3_lower - 2K)) for 2K < I3_lower",
        },
        "sides": side_records,
        "green_profile_certificates": profile_records,
        "falsification": falsification,
        "conditional_H": conditional_h,
        "meta_negative": {
            "statement": (
                "[THEOREM] Every modification obtained solely by dropping rows, "
                "orbit/momentum averaging pointwise constraints, or replacing them "
                "by a convex consequence has a feasible set containing the audited "
                "LP feasible set, hence a floor <= floor(K,L) <= u_L -> 0. Such "
                "relaxations cannot yield a useful uniform-in-L floor. A genuinely "
                "stronger paired-momentum constraint would have to be new valid "
                "input not implied by the audited class (for example random-current, "
                "Simon-Lieb, GKS-II four-point, or ABF-type input); it is not "
                "established here."
            )
        },
        "closed_forms": {
            "regime_2K_ge_D": (
                "floor(K,L) = 1 - C_L(0)/(2K) exactly (all-ceiling profile is "
                "feasible and componentwise maximal)"
            ),
            "regime_2K_le_D": (
                "floor(K,L) <= u_L = (-min_z C_L)/D_L (explicit feasible Green "
                "profile 1/(D_L lambda))"
            ),
            "positivity": "floor(K,L) > 0 strictly for every K > 0, every L >= 2",
        },
    }

    elapsed = time.time() - started
    payload = {
        "provenance": {
            "script": "experiments/e93_mag_floor.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": (
                "All decisions are exact: Fraction for rational-cosine sides "
                "(L=4,6) and for all Galois-invariant outputs; Q(sqrt2), "
                "Q(sqrt5), Q(sqrt3) pairs of Fractions for the quadratic-cosine "
                "sides. No floating point enters any decision. The constants "
                "pi < 355/113 and sqrt2 < 99/70 are classical exact rational "
                "bounds; I3 enters only through its certified rational lower "
                "endpoint."
            ),
            "elapsed_seconds": elapsed,
        },
        "data": data,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"\nPASS  ({elapsed:.1f}s)  wrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
