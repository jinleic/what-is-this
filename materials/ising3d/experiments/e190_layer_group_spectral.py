"""Integrate the layer-group spectral obstruction and two exact grid controls.

The primary result is structural and negative: varying conjugators invalidate the
proposed passage from pointwise Gaussian spectra to a fixed quadratic Lie
algebra.  Pair-product calculations are deliberately limited to the established
2x3 and 2x4 controls at t=1/3; they are not used as an all-size argument.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from e188_gaussian_group import run_formalization  # noqa: E402
from e189_varying_conjugator import run_witness  # noqa: E402

ARTIFACT = ROOT / "results" / "spectral" / "layer_group_spectral.json"
PRIME = 1_000_003
CPU_BUDGET_SECONDS = 180.0
RSS_CAP_BYTES = 2_000_000_000
FLOAT64_EXACT = 1 << 53
INT64_MAX = (1 << 63) - 1
EXPECTED_CONTROLS = {
    "2x3": {"sites": 6, "edges": 7, "pair_gcd_degree": 385, "r_lower_bound": 1631},
    "2x4": {"sites": 8, "edges": 10, "pair_gcd_degree": 9329, "r_lower_bound": 23311},
}


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class Budget:
    def __init__(self) -> None:
        self.started = time.process_time()
        self.stage = "initialization"

    def used(self) -> float:
        return time.process_time() - self.started

    def check(self, stage: str | None = None) -> None:
        if stage is not None:
            self.stage = stage
        used = self.used()
        rss = max_rss_bytes()
        if used >= CPU_BUDGET_SECONDS:
            raise RuntimeError(
                f"process-time budget exceeded at {self.stage}: {used:.3f}>={CPU_BUDGET_SECONDS}"
            )
        if rss >= RSS_CAP_BYTES:
            raise RuntimeError(f"RSS cap exceeded at {self.stage}: {rss}>={RSS_CAP_BYTES}")


def check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def grid_bonds(rows: int, cols: int) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for row in range(rows):
        for col in range(cols):
            site = row * cols + col
            if col + 1 < cols:
                bonds.append((site, site + 1))
            if row + 1 < rows:
                bonds.append((site, site + cols))
    return bonds


def build_integer_layer(
    sites: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    q: Fraction,
    budget: Budget,
) -> tuple[int, list[list[int]]]:
    """Build scale*P_t*D_q*P_t over Z by integer Kronecker sweeps."""
    dim = 1 << sites
    scale = t.denominator ** (2 * sites) * q.denominator ** len(bonds)
    diagonal: list[int] = []
    for state in range(dim):
        aligned = 0
        for left, right in bonds:
            left_bit = (state >> (sites - 1 - left)) & 1
            right_bit = (state >> (sites - 1 - right)) & 1
            aligned += left_bit == right_bit
        diagonal.append(q.numerator**aligned * q.denominator ** (len(bonds) - aligned))

    matrix = [[0] * dim for _ in range(dim)]
    for state, value in enumerate(diagonal):
        matrix[state][state] = value

    def apply_left(rows: list[list[int]]) -> None:
        for site in range(sites):
            bit = 1 << (sites - 1 - site)
            for upper in range(dim):
                if upper & bit:
                    continue
                lower = upper | bit
                row0, row1 = rows[upper], rows[lower]
                for col in range(dim):
                    x, y = row0[col], row1[col]
                    row0[col] = t.denominator * x + t.numerator * y
                    row1[col] = t.numerator * x + t.denominator * y
            budget.check(f"integer Kronecker sweep site {site + 1}/{sites}")

    apply_left(matrix)
    matrix = [list(column) for column in zip(*matrix)]
    apply_left(matrix)
    return scale, matrix


def modmul(left: np.ndarray, right: np.ndarray, prime: int) -> np.ndarray:
    length = left.shape[1]
    assert length * (prime - 1) ** 2 < FLOAT64_EXACT
    product = np.asarray(left, dtype=np.float64) @ np.asarray(right, dtype=np.float64)
    return np.mod(product, float(prime)).astype(np.int64)


def moddot(left: np.ndarray, right: np.ndarray, prime: int) -> int:
    assert left.size * (prime - 1) ** 2 < INT64_MAX
    return int(np.dot(left, right)) % prime


def trim(poly: np.ndarray) -> np.ndarray:
    for index in range(min(poly.size, 64)):
        if poly[index]:
            return poly if index == 0 else poly[index:]
    nonzero = np.flatnonzero(poly)
    return np.zeros(1, dtype=np.int64) if nonzero.size == 0 else poly[nonzero[0] :]


def is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def monic(poly: np.ndarray, prime: int) -> np.ndarray:
    value = trim(np.asarray(poly, dtype=np.int64))
    if is_zero(value):
        return value
    lead = int(value[0]) % prime
    return value if lead == 1 else (value * pow(lead, -1, prime)) % prime


def remainder(
    dividend: np.ndarray, divisor: np.ndarray, prime: int, budget: Budget
) -> np.ndarray:
    assert int(divisor[0]) == 1
    work = np.asarray(dividend, dtype=np.int64).copy() % prime
    divisor_degree = divisor.size - 1
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset % 256 == 0:
            budget.check("polynomial Euclidean division")
        lead = int(work[offset])
        if lead:
            window = work[offset : offset + divisor_degree + 1]
            work[offset : offset + divisor_degree + 1] = (window - lead * divisor) % prime
    return trim(work[dividend_degree - divisor_degree + 1 :])


def polynomial_gcd(
    left: np.ndarray, right: np.ndarray, prime: int, budget: Budget
) -> np.ndarray:
    left, right = monic(left, prime), monic(right, prime)
    while not is_zero(right):
        budget.check("polynomial Euclidean gcd")
        left, right = right, monic(remainder(left, right, prime, budget), prime)
    return left


def derivative(poly: np.ndarray, prime: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    return np.ascontiguousarray(
        ((np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % prime)[::-1]
    )


def newton(
    power_sums: np.ndarray,
    degree: int,
    prime: int,
    inverses: list[int],
    budget: Budget,
) -> np.ndarray:
    """Recover a monic descending polynomial by exact chunked Newton identities."""
    chunk = (FLOAT64_EXACT - 1) // ((prime - 1) ** 2)
    assert chunk >= 1
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    sums = np.ascontiguousarray(power_sums[: degree + 1], dtype=np.float64)
    reverse_buffer = np.zeros(degree + 1, dtype=np.float64)
    reverse_buffer[degree] = 1.0
    for k in range(1, degree + 1):
        if k % 1024 == 0:
            budget.check(f"Newton reconstruction degree {degree}: k={k}")
        total = 0
        base = degree - k
        lo = 1
        while lo <= k:
            hi = min(k, lo + chunk - 1)
            total += int(
                np.dot(sums[lo : hi + 1], reverse_buffer[base + lo : base + hi + 1])
            ) % prime
            lo = hi + 1
        value = (-(total % prime) * inverses[k]) % prime
        coefficients[k] = value
        reverse_buffer[degree - k] = float(value)
    return coefficients


def coefficient_digest(descending: np.ndarray) -> str:
    ascending = descending[::-1]
    payload = ",".join(str(int(value)) for value in ascending).encode()
    return hashlib.sha256(payload).hexdigest()


def pair_certificate(
    matrix_rows: list[list[int]], prime: int, budget: Budget
) -> dict[str, object]:
    dim = len(matrix_rows)
    reduced = np.array([[value % prime for value in row] for row in matrix_rows], dtype=np.int64)
    slots = comb(dim, 2)
    needed = 2 * slots

    budget.check(f"{dim}x{dim} power traces")
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for exponent in range(1, dim + 1):
        power = modmul(power, reduced, prime)
        traces.append(int(power.trace()) % prime)
        if exponent % 32 == 0:
            budget.check(f"power traces {exponent}/{dim}")

    inverses = [0] + [pow(k, -1, prime) for k in range(1, needed + 1)]
    charpoly = newton(np.array(traces, dtype=np.int64), dim, prime, inverses, budget)
    tail = charpoly[1:]
    for exponent in range(dim + 1, needed + 1):
        window = np.array(traces[exponent - 1 : exponent - dim - 1 : -1], dtype=np.int64)
        traces.append((-moddot(tail, window, prime)) % prime)
        if exponent % 2048 == 0:
            budget.check(f"Cayley-Hamilton trace recurrence {exponent}/{needed}")

    half = pow(2, -1, prime)
    pair_sums = np.zeros(slots + 1, dtype=np.int64)
    for exponent in range(1, slots + 1):
        pair_sums[exponent] = (
            ((traces[exponent] * traces[exponent] - traces[2 * exponent]) % prime) * half
        ) % prime
    pair_poly = newton(pair_sums, slots, prime, inverses, budget)
    pair_derivative = derivative(pair_poly, prime)
    gcd_poly = polynomial_gcd(pair_poly, pair_derivative, prime, budget)
    divides_pair = is_zero(remainder(pair_poly, gcd_poly, prime, budget))
    divides_derivative = is_zero(remainder(pair_derivative, gcd_poly, prime, budget))
    return {
        "prime": prime,
        "dimension": dim,
        "pair_slots": slots,
        "charpoly_degree": charpoly.size - 1,
        "pair_polynomial_degree": pair_poly.size - 1,
        "pair_gcd_degree": gcd_poly.size - 1,
        "distinct_pair_products_lower_bound": slots - (gcd_poly.size - 1),
        "charpoly_sha256": coefficient_digest(charpoly),
        "pair_polynomial_sha256": coefficient_digest(pair_poly),
        "pair_gcd_sha256": coefficient_digest(gcd_poly),
        "gcd_divides_pair": bool(divides_pair),
        "gcd_divides_derivative": bool(divides_derivative),
    }


def run_control(label: str, cols: int, budget: Budget) -> dict[str, object]:
    started = time.process_time()
    sites = 2 * cols
    bonds = grid_bonds(2, cols)
    t = Fraction(1, 3)
    q = Fraction(5, 3)
    scale, matrix = build_integer_layer(sites, bonds, t, q, budget)
    symmetric = all(matrix[i][j] == matrix[j][i] for i in range(len(matrix)) for j in range(i))
    certificate = pair_certificate(matrix, PRIME, budget)
    return {
        "tag": "[COMPUTATION]",
        "layer": label,
        "sites": sites,
        "bonds": [list(edge) for edge in bonds],
        "edge_count": len(bonds),
        "t": "1/3",
        "q": "5/3",
        "integer_scale": str(scale),
        "integer_matrix_symmetric": symmetric,
        "gaussian_ceiling": 3**sites - 2**sites,
        **certificate,
        "process_time_seconds": round(time.process_time() - started, 3),
    }


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_all() -> tuple[dict[str, object], list[dict[str, object]], Budget]:
    budget = Budget()
    checks: list[dict[str, object]] = []
    formalization = run_formalization()
    witness = run_witness()
    for item in formalization["checks"]:
        check(checks, f"e188: {item['name']}", bool(item["passed"]), str(item["detail"]))
    for item in witness["checks"]:
        check(checks, f"e189: {item['name']}", bool(item["passed"]), str(item["detail"]))

    controls: list[dict[str, object]] = []
    for label, cols in (("2x3", 3), ("2x4", 4)):
        row = run_control(label, cols, budget)
        controls.append(row)
        expected = EXPECTED_CONTROLS[label]
        check(
            checks,
            f"{label}: raw graph dimensions",
            row["sites"] == expected["sites"] and row["edge_count"] == expected["edges"],
            f"n={row['sites']}, edges={row['edge_count']}",
        )
        check(
            checks,
            f"{label}: exact modular pair gcd degree",
            row["pair_gcd_degree"] == expected["pair_gcd_degree"],
            f"degree={row['pair_gcd_degree']} at p={PRIME}",
        )
        check(
            checks,
            f"{label}: one-sided characteristic-zero lower bound",
            row["distinct_pair_products_lower_bound"] == expected["r_lower_bound"]
            and row["distinct_pair_products_lower_bound"] > row["gaussian_ceiling"],
            (
                f"r_Q>={row['distinct_pair_products_lower_bound']}>"
                f"{row['gaussian_ceiling']}"
            ),
        )
        check(
            checks,
            f"{label}: polynomial divisibility audit",
            bool(row["gcd_divides_pair"]) and bool(row["gcd_divides_derivative"]),
            "monic gcd divides both pair polynomial and derivative",
        )
        check(
            checks,
            f"{label}: exact integer construction is symmetric",
            bool(row["integer_matrix_symmetric"]),
            f"scale={row['integer_scale']}",
        )

    budget.check("final resource audit")
    check(
        checks,
        "process-time and RSS budgets",
        budget.used() < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        (
            f"CPU={budget.used():.3f}s<{CPU_BUDGET_SECONDS}; "
            f"RSS={max_rss_bytes()}<{RSS_CAP_BYTES}"
        ),
    )

    data = {
        "formalization": formalization["data"],
        "varying_conjugator_obstruction": witness["data"],
        "physical_curve_non_density": {
            "tag": "[THEOREM]",
            "statement": (
                "Even a hypothetical proof that the ambient two-parameter family is not "
                "contained in the Gaussian conjugacy saturation would not imply that its "
                "physical one-parameter curve is not contained."
            ),
            "curve_equation": "F_phys(x,y)=(y-1)x^2-(y+1)=0",
            "exact_obstruction": (
                "a nonzero ambient obstruction polynomial may be divisible by F_phys and "
                "therefore vanish identically after isotropic specialization"
            ),
            "finite_exception_requirement": (
                "one must exhibit an obstruction whose restriction modulo F_phys is nonzero; "
                "only then does its nonzero univariate numerator have finitely many roots"
            ),
        },
        "fixed_basis_dimension_scope": {
            "tag": "[THEOREM]",
            "statement": (
                "A Lie-dimension or Clifford-grade bound excludes one fixed conjugation of "
                "the whole generated group into quadratics. It does not exclude pointwise "
                "conjugations C(a,b) that vary with the transfer parameters."
            ),
            "vandermonde_scope": (
                "ad_A grade separation remains a valid statement inside Lie(A,B), but supplies "
                "no mechanism that identifies the independently varying spectral conjugators"
            ),
        },
        "control_computations": controls,
        "conclusion": {
            "tag": "[UNRESOLVED]",
            "isotropic_all_size_spectral_nogo_proved": False,
            "exceptional_set": None,
            "exact_result": (
                "the proposed group/algebra differentiation bridge is false, with an exact "
                "Spin(2) varying-conjugator witness and an explicit tangent absorption identity"
            ),
            "scope": (
                "this kills the proposed proof route; it neither produces a Gaussian branching "
                "layer nor disproves a future non-localizing isotropic spectral theorem"
            ),
            "what_remains": (
                "construct a conjugacy-invariant spectral polynomial nonzero modulo the physical "
                "curve for every branching layer, or add hypotheses that force a single "
                "parameter-independent conjugator"
            ),
        },
    }
    return data, checks, budget


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]], budget: Budget
) -> dict[str, object]:
    sources = [
        ROOT / "experiments" / "e188_gaussian_group.py",
        ROOT / "experiments" / "e189_varying_conjugator.py",
        Path(__file__),
    ]
    return {
        "meta": {
            "experiment": "e190_layer_group_spectral",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e190_layer_group_spectral.py",
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy_version": np.__version__,
            "arithmetic": (
                "symbolic Q(r,u,t) identities plus exact F_p power traces, Newton identities, "
                "Cayley-Hamilton recurrence, and monic Euclidean gcd"
            ),
            "modular_direction": (
                "deg gcd_Q <= deg gcd_Fp; controls use only r_Q >= pair_slots-deg gcd_Fp"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "observed_process_time_seconds": round(budget.used(), 3),
                "observed_peak_rss_bytes": max_rss_bytes(),
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in sources
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks, budget = run_all()
    artifact = make_artifact(data, checks, budget)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if all(bool(item["passed"]) for item in checks):
        print("PASS e190 layer-group spectral obstruction", flush=True)
        return 0
    print("FAIL e190 layer-group spectral obstruction", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
