#!/usr/bin/env python3
"""Clean-room verifier for experiments/e93_mag_floor.py.

This test deliberately does not import the producer.  It rebuilds exact
quadratic-field cosine data, zero-mode-removed torus Green functions, every
orbit-symmetrised LP matrix, every stored primal/dual/basis certificate, and
the full-site feasible Green-profile certificates.  It also verifies the
explicit O(1/L) rate arithmetic used in proofs/mag_floor.md.

Run from repo root:
    timeout 3600 .venv/bin/python tests/test_mag_floor.py
"""

from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "mag_floor.json"
E85 = ROOT / "results" / "bounds" / "kc_interval2.json"

SIDES = (4, 6, 8, 10, 12)
LP_SIDES = (4, 6, 8)
K_POINTS = (Fraction(1, 5), Fraction(21, 100), Fraction(11, 50), Fraction(6, 25))
PI_UB = Fraction(355, 113)
SQRT2_UB = Fraction(99, 70)
I3_LOW = Fraction(
    "0.505462019717326006052004053227140259985129014817420892188993487886028773451173816800537247069896037925625"
)


class Quad:
    """Independent Q(sqrt d) implementation, a+b*sqrt(d), exact Fractions."""

    __slots__ = ("a", "b", "d")

    def __init__(self, a=0, b=0, d=2):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.d = d

    @staticmethod
    def _cast(other: object, d: int):
        if isinstance(other, Quad):
            if other.d != d:
                raise AssertionError("mixed quadratic fields")
            return other
        if isinstance(other, (int, Fraction)):
            return Quad(other, 0, d)
        return NotImplemented

    def __add__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(self.a + o.a, self.b + o.b, self.d)

    __radd__ = __add__

    def __sub__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(self.a - o.a, self.b - o.b, self.d)

    def __rsub__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return o - self

    def __mul__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return Quad(self.a * o.a + self.d * self.b * o.b, self.a * o.b + self.b * o.a, self.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        norm = o.a * o.a - self.d * o.b * o.b
        if norm == 0:
            raise ZeroDivisionError("quadratic-field denominator is zero")
        return Quad(
            (self.a * o.a - self.d * self.b * o.b) / norm,
            (self.b * o.a - self.a * o.b) / norm,
            self.d,
        )

    def __rtruediv__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return o / self

    def __neg__(self):
        return Quad(-self.a, -self.b, self.d)

    def __eq__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return self.a == o.a and self.b == o.b

    def __ne__(self, other):
        result = self.__eq__(other)
        return NotImplemented if result is NotImplemented else not result

    def sign(self) -> int:
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        q = self.a / self.b
        # a+b sqrt(d)=b(q+sqrt(d)).  If q>=0 then q+sqrt(d)>0;
        # otherwise its sign is decided by q^2 < d.
        inside_positive = q >= 0 or q * q < self.d
        if self.b > 0:
            return 1 if inside_positive else -1
        return -1 if inside_positive else 1

    def __lt__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() < 0

    def __le__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() <= 0

    def __gt__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() > 0

    def __ge__(self, other):
        o = self._cast(other, self.d)
        if o is NotImplemented:
            return NotImplemented
        return (self - o).sign() >= 0

    def rational(self) -> Fraction:
        assert self.b == 0, f"not rational: {self.a}+({self.b})sqrt({self.d})"
        return self.a


def parse_field(text: str):
    if text.startswith("Q["):
        a, b, d = text[2:-1].split("|")
        return Quad(Fraction(a), Fraction(b), int(d))
    return Fraction(text)


def as_fraction(value) -> Fraction:
    return value.rational() if isinstance(value, Quad) else Fraction(value)


def field_sum(values, zero):
    total = zero
    for value in values:
        total = total + value
    return total


def exact_cosines(side: int):
    """Independent exact cosine construction and algebraic controls."""

    if side == 4:
        base = Fraction(0)
        assert base == 0
    elif side == 6:
        base = Fraction(1, 2)
        assert 2 * base == 1
    elif side == 8:
        base = Quad(0, Fraction(1, 2), 2)
        assert base * base == Fraction(1, 2)
    elif side == 10:
        base = Quad(Fraction(1, 4), Fraction(1, 4), 5)
        assert 4 * base * base - 2 * base - 1 == 0 and base > 0
    elif side == 12:
        base = Quad(0, Fraction(1, 2), 3)
        assert 4 * base * base == 3 and base > 0
    else:
        raise AssertionError(f"unsupported side {side}")
    values = [base * 0 + 1, base]
    for j in range(1, side):
        values.append(2 * base * values[j] - values[j - 1])
    values = values[:side]
    assert values[side // 2] == -1
    assert 2 * base * values[-1] - values[-2] == values[0]
    return values


def signed_orbit(point: tuple[int, int, int], side: int):
    result = set()
    for perm in itertools.permutations(range(3)):
        p = tuple(point[i] for i in perm)
        for signs in itertools.product((-1, 1), repeat=3):
            result.add(tuple((sign * coord) % side for sign, coord in zip(signs, p, strict=True)))
    return tuple(sorted(result))


def all_orbits(side: int):
    pending = set(itertools.product(range(side), repeat=3))
    answer = []
    while pending:
        orbit = signed_orbit(min(pending), side)
        pending.difference_update(orbit)
        answer.append(orbit)
    return tuple(answer)


def exact_green(side: int):
    """Rebuild C_L exactly from raw torus mode sums, one site orbit at a time."""

    cos = exact_cosines(side)
    n = side**3
    orbits = all_orbits(side)
    modes = [k for k in itertools.product(range(side), repeat=3) if k != (0, 0, 0)]
    lambdas = []
    for k in modes:
        lam = cos[0] * 0 + 3
        for axis in range(3):
            lam = lam - cos[k[axis]]
        assert lam > 0
        lambdas.append(lam)

    by_rep = {}
    for orbit in orbits:
        z = orbit[0]
        total = cos[0] * 0
        for k, lam in zip(modes, lambdas, strict=True):
            phase = sum(k[i] * z[i] for i in range(3)) % side
            total = total + cos[phase] / lam
        by_rep[z] = as_fraction(total / n)
    by_site = {}
    for orbit in orbits:
        for z in orbit:
            by_site[z] = by_rep[orbit[0]]

    # All-site exact identities, independent of the producer's implementation.
    assert sum(by_site.values(), Fraction(0)) == 0
    for z, value in by_site.items():
        lap = 3 * value
        for axis in range(3):
            for delta in (-1, 1):
                w = list(z)
                w[axis] = (w[axis] + delta) % side
                lap -= Fraction(1, 2) * by_site[tuple(w)]
        assert lap == (1 if z == (0, 0, 0) else 0) - Fraction(1, n)

    j2_total = cos[0] * 0
    for k, lam in zip(modes, lambdas, strict=True):
        j2_total = j2_total + 1 / (lam * lam)
    j2 = as_fraction(j2_total / n)
    return by_site, orbits, j2


def rebuild_lp(side: int, coupling: Fraction):
    """Independent reconstruction of the exact orbit-reduced LP matrix."""

    cos = exact_cosines(side)
    n = side**3
    zero = cos[0] * 0
    one = zero + 1
    orbits = all_orbits(side)
    mode_orbits = tuple(o for o in orbits if o != ((0, 0, 0),))
    site_orbits = tuple(o for o in orbits if o != ((0, 0, 0),))
    mult = [len(o) for o in mode_orbits]
    dispersions = []
    for orbit in mode_orbits:
        lam = zero + 3
        for i in range(3):
            lam = lam - cos[orbit[0][i]]
        assert lam > zero
        dispersions.append(lam)
    phase_sums = []
    for z_orbit in site_orbits:
        z = z_orbit[0]
        row = []
        for k_orbit in mode_orbits:
            phase = zero
            for k in k_orbit:
                phase = phase + cos[sum(k[i] * z[i] for i in range(3)) % side]
            row.append(phase)
        phase_sums.append(row)

    a, b, labels = [], [], []
    for j, lam in enumerate(dispersions):
        row = [zero for _ in mode_orbits]
        row[j] = one
        a.append(row)
        b.append(one / (2 * coupling * lam))
        labels.append(f"ceiling_mode_{j}")
    for orbit, phase in zip(site_orbits, phase_sums, strict=True):
        z = orbit[0]
        gap = [one * m - q for m, q in zip(mult, phase, strict=True)]
        a.append(list(gap))
        b.append(one * n)
        labels.append("gks_nonnegative_" + "_".join(map(str, z)))
        a.append([-q for q in gap])
        b.append(zero)
        labels.append("spin_upper_" + "_".join(map(str, z)))
    a.append([one * m for m in mult])
    b.append(one * n)
    labels.append("zero_mode_nonnegative")
    e_index = next(j for j, orbit in enumerate(site_orbits) if (1, 0, 0) in orbit)
    gap_e = [one * m - q for m, q in zip(mult, phase_sums[e_index], strict=True)]
    a.append(list(gap_e))
    b.append(one * (n - 1) / (6 * coupling))
    labels.append("energy_nearest_neighbour")
    return {
        "n": n,
        "zero": zero,
        "one": one,
        "mode_orbits": mode_orbits,
        "mult": mult,
        "disp": dispersions,
        "a": a,
        "b": b,
        "c": [one * m for m in mult],
        "labels": labels,
    }


def nonsingular(matrix, zero):
    """Independent exact Gaussian-elimination rank check for a square basis."""

    work = [list(row) for row in matrix]
    size = len(work)
    assert all(len(row) == size for row in work)
    for col in range(size):
        pivot = next((r for r in range(col, size) if work[r][col] != zero), None)
        if pivot is None:
            return False
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
        pivot_value = work[col][col]
        for r in range(col + 1, size):
            factor = work[r][col] / pivot_value
            if factor != zero:
                work[r] = [x - factor * y for x, y in zip(work[r], work[col], strict=True)]
    return True


def verify_lp_certificate(side: int, cert: dict[str, object], stored_side: dict[str, object]):
    coupling = Fraction(cert["coupling"])
    problem = rebuild_lp(side, coupling)
    zero, one = problem["zero"], problem["one"]
    a, b, c = problem["a"], problem["b"], problem["c"]
    x = [parse_field(v) for v in cert["primal_x"]]
    y = [parse_field(v) for v in cert["dual_y"]]
    assert len(x) == len(c)
    assert len(y) == len(a)
    assert cert["constraint_labels"] == problem["labels"]

    # Rebuilt mode-order contract pins the matrix columns and field embedding.
    stored_modes = cert["mode_orbits"]
    assert len(stored_modes) == len(problem["mode_orbits"])
    for saved, orbit, multiplicity, dispersion in zip(
        stored_modes,
        problem["mode_orbits"],
        problem["mult"],
        problem["disp"],
        strict=True,
    ):
        assert saved["representative"] == list(orbit[0])
        assert saved["multiplicity"] == multiplicity
        assert parse_field(saved["lambda"]) == dispersion

    # Exact primal feasibility.
    assert all(value >= zero for value in x)
    row_values = [field_sum((q * value for q, value in zip(row, x, strict=True)), zero) for row in a]
    assert all(value <= bound for value, bound in zip(row_values, b, strict=True))
    primal = field_sum((q * value for q, value in zip(c, x, strict=True)), zero)

    # Exact dual feasibility and equality of objectives.
    assert all(value >= zero for value in y)
    for j in range(len(c)):
        lhs = field_sum((y[i] * a[i][j] for i in range(len(a))), zero)
        assert lhs >= c[j]
    dual = field_sum((value * bound for value, bound in zip(y, b, strict=True)), zero)
    assert primal == dual
    assert primal == parse_field(cert["objective_sum_nonzero_modes"])
    assert primal / problem["n"] == parse_field(cert["S_star"])
    assert one - primal / problem["n"] == parse_field(cert["magnetization_square_floor"])
    assert cert["primal_dual_equal"] is True

    # Stored basis is a genuine nonsingular basis for the stored primal point.
    basis = cert["basis"]
    assert len(basis) == len(a) and len(set(basis)) == len(basis)
    full_a = [list(row) + [one if r == j else zero for j in range(len(a))] for r, row in enumerate(a)]
    slacks = [bound - value for bound, value in zip(b, row_values, strict=True)]
    full_values = x + slacks
    assert all(0 <= index < len(full_values) for index in basis)
    assert all(full_values[index] == zero for index in range(len(full_values)) if index not in basis)
    basis_matrix = [[full_a[r][basis[col]] for col in range(len(basis))] for r in range(len(a))]
    assert nonsingular(basis_matrix, zero)
    for r in range(len(a)):
        assert field_sum((full_a[r][index] * full_values[index] for index in basis), zero) == b[r]

    # Closed-form all-ceiling regime and exact profile upper bound.
    d = Fraction(stored_side["D_L"])
    c0 = Fraction(stored_side["C_L_at_origin"])
    floor = parse_field(cert["magnetization_square_floor"])
    if 2 * coupling >= d:
        assert cert["regime"] == "closed_form"
        assert floor == 1 - c0 / (2 * coupling)
        assert cert["closed_form_floor"] == str(1 - c0 / (2 * coupling))
    else:
        assert cert["regime"] == "lp"
    if 2 * coupling <= d:
        assert floor <= Fraction(stored_side["u_L"])
    assert floor > zero


def verify_side(data: dict[str, object], side: int, constants: dict[str, object]):
    saved = data["sides"][str(side)]
    green, orbits, j2 = exact_green(side)
    n = side**3
    c0 = green[(0, 0, 0)]
    min_c = min(green.values())
    min_sites = sorted(z for z, v in green.items() if v == min_c)
    d = c0 - min_c
    u = -min_c / d

    assert c0 == Fraction(saved["C_L_at_origin"])
    assert min_c == Fraction(saved["minimum_C_L"])
    assert saved["minimum_sites"] == [list(z) for z in min_sites[:8]]
    assert d == Fraction(saved["D_L"])
    assert d / 2 == Fraction(saved["kappa_L"])
    assert u == Fraction(saved["u_L"])
    assert side * u == Fraction(saved["u_L_times_L"])
    assert j2 == Fraction(saved["J2_L"])

    # Independent orbit table, including exact multiplicity and every stored value.
    saved_orbits = saved["green_orbits"]
    assert len(saved_orbits) == len(orbits)
    for row, orbit in zip(saved_orbits, orbits, strict=True):
        assert row["representative"] == list(orbit[0])
        assert row["multiplicity"] == len(orbit)
        assert Fraction(row["C_L"]) == green[orbit[0]]

    # M1, M2, M3, u_L rate, and exact below-I3 check.
    pi2_ub = PI_UB * PI_UB
    pi4_ub = pi2_ub * pi2_ub
    sigma4 = 4 * pi2_ub + pi4_ub / 45
    sqrt_sigma4 = Fraction(323, 50)
    c_star = Fraction(5168, 525)
    c_delta = SQRT2_UB * (pi2_ub + 2 * PI_UB) / 16
    a_delta = c_delta + 1
    assert sqrt_sigma4 * sqrt_sigma4 >= sigma4
    assert min_c * min_c <= 4 * j2 / n
    assert j2 <= side * sigma4 / 64
    assert c0 >= I3_LOW - a_delta / side
    assert u <= c_star / side
    assert d < I3_LOW
    assert c_star >= a_delta / I3_LOW
    assert saved["M1_control_4J2_over_N"] == str(4 * j2 / n)

    # The exact Green profile c=1/D_L is valid whenever 2K<=D_L.
    for profile in data["green_profile_certificates"]:
        if profile["side"] != side:
            continue
        k = Fraction(profile["coupling"])
        assert 2 * k <= d
        coeff = Fraction(profile["c_coefficient"])
        assert coeff == 1 / d
        assert Fraction(profile["S_profile"]) == c0 / d
        assert Fraction(profile["floor_upper_bound_u_L"]) == u
        for z, value in green.items():
            g = 1 - (c0 - value) / d
            assert 0 <= g <= 1
        # Energy row, direct rather than merely invoking the producer's claim.
        ge = 1 - (c0 - green[(1, 0, 0)]) / d
        assert ge >= 1 - Fraction(n - 1, 6 * k * n)

    if side in LP_SIDES:
        certs = saved["lp"]
        assert len(certs) == len(K_POINTS)
        for cert in certs:
            verify_lp_certificate(side, cert, saved)

    return {"side": side, "green": green, "d": d, "u": u}


def main() -> None:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    data = result["data"]
    constants = data["constants"]

    # Exact constant envelope has the intended producer values.
    assert constants["pi_upper"] == "355/113"
    assert constants["sqrt2_upper"] == "99/70"
    assert Fraction(constants["c_star"]) == Fraction(5168, 525)
    pi2_ub = PI_UB * PI_UB
    expected_c_delta = SQRT2_UB * (pi2_ub + 2 * PI_UB) / 16
    assert Fraction(constants["c_delta"]) == expected_c_delta
    assert Fraction(constants["c_delta_plus_1"]) == expected_c_delta + 1
    assert Fraction(constants["I3_lower_rational"]) == I3_LOW

    rebuilt = [verify_side(data, side, constants) for side in SIDES]

    # Cross-artifact independent regression: all previously stored L=4,6 floors match e85.
    e85 = json.loads(E85.read_text(encoding="utf-8"))
    e85_floors = {}
    for record in e85["data"]["route2_finite_volume_floor"]["records"]:
        for cert in record["finite_floor_certificates"]:
            e85_floors[(record["side"], cert["coupling"])] = cert["magnetization_square_floor"]
    for side in (4, 6):
        for cert in data["sides"][str(side)]["lp"]:
            assert cert["matches_e85"] is True
            assert cert["magnetization_square_floor"] == e85_floors[(side, cert["coupling"])]

    # Exact finite-scale trends claimed in the result envelope.
    kappas = [item["d"] / 2 for item in rebuilt]
    assert all(a < b for a, b in zip(kappas, kappas[1:]))
    c_star = Fraction(5168, 525)
    a_delta = expected_c_delta + 1
    for k in K_POINTS:
        floors = []
        for side in LP_SIDES:
            cert = next(c for c in data["sides"][str(side)]["lp"] if c["coupling"] == str(k))
            floors.append(parse_field(cert["magnetization_square_floor"]))
        assert all(a > b for a, b in zip(floors, floors[1:]))
        assert data["falsification"]["floor_sequences"][str(k)]["strictly_decreasing_in_L"] is True

    # Exact countercertificate to the natural nondecreasing dyadic LP-floor
    # proposal F(K, 2L) >= F(K, L), used in proofs/mag_floor.md.
    f4 = parse_field(next(c for c in data["sides"]["4"]["lp"] if c["coupling"] == "1/5")["magnetization_square_floor"])
    f8 = parse_field(next(c for c in data["sides"]["8"]["lp"] if c["coupling"] == "1/5")["magnetization_square_floor"])
    assert f4 - f8 == Fraction(28697, 626688) + Quad(0, Fraction(1, 816), 2)
    assert f4 > f8
    dyadic_saved = data["falsification"]["dyadic_countercertificate"]
    assert dyadic_saved["coupling_K"] == "1/5"
    assert dyadic_saved["smaller_side_L"] == 4 and dyadic_saved["larger_side_2L"] == 8
    assert parse_field(dyadic_saved["floor_difference_F_L_minus_F_2L"]) == f4 - f8
    assert dyadic_saved["strictly_positive"] is True

    # This is the exact algebraic endpoint split used in the proof: at 2K=I3,
    # c=min(I3^-1,D^-1) gives either the M3 bound or the u_L bound, both <=c*/L.
    assert a_delta / I3_LOW <= c_star

    # Route-(3) conditional: exact finite checks pass at L=4,6,8, but the
    # independently recomputed universal rate certificate disproves its all-L premise.
    h = data["conditional_H"]
    h_k = Fraction(h["coupling_K"])
    h_m = Fraction(h["floor_m"])
    assert h_k == Fraction(11, 50) and h_m == Fraction(1, 100)
    assert h["predicate"] == "H_{K,m}(L): exact LP floor(K,L) >= m"
    assert [row["side"] for row in h["exact_verified_sides"]] == [4, 6, 8]
    assert all(parse_field(row["floor"]) >= h_m and row["H_holds"] is True for row in h["exact_verified_sides"])
    h_gap = I3_LOW - 2 * h_k
    h_condition_start = -(-(a_delta / h_gap).numerator // (a_delta / h_gap).denominator)
    h_rate_ratio = c_star / h_m
    h_rate_start = h_rate_ratio.numerator // h_rate_ratio.denominator + 1
    h_first_even = max(4, h_condition_start, h_rate_start)
    if h_first_even % 2:
        h_first_even += 1
    saved_counter = h["rate_countercertificate"]
    assert saved_counter["condition_start_L"] == h_condition_start
    assert saved_counter["strict_rate_start_L"] == h_rate_start
    assert saved_counter["first_even_L"] == h_first_even
    assert 2 * h_k <= I3_LOW - a_delta / h_first_even
    assert c_star / h_first_even < h_m

    print("PASS test_mag_floor: exact Green data, full-site profiles, 12 primal-dual LP certificates, bases, dyadic countercertificate, and rate controls verified")


if __name__ == "__main__":
    main()
