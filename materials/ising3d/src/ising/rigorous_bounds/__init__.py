"""Certified ingredients for elementary bounds on the cubic Ising ``K_c``.

The combinatorial routines use Python integers only.  Transcendental point
values use :mod:`mpmath`; certified enclosures use its directed-rounding
interval context and are returned as decimal endpoint strings so callers do
not silently coerce them to binary floats.
"""

from __future__ import annotations

from collections.abc import Iterable
from itertools import permutations, product
import multiprocessing
from urllib.request import Request, urlopen

import mpmath as mp


OEIS_A001412_URL = "https://oeis.org/A001412/b001412.txt"
OEIS_A001412 = (
    1,
    6,
    30,
    150,
    726,
    3534,
    16926,
    81390,
    387966,
    1853886,
    8809878,
    41934150,
    198842742,
    943974510,
    4468911678,
)

_DIRECTIONS = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
)
_CUBIC_SYMMETRIES = tuple(
    (permutation, signs)
    for permutation in permutations(range(3))
    for signs in product((-1, 1), repeat=3)
)


def fetch_oeis_a001412(
    max_steps: int, url: str = OEIS_A001412_URL, timeout: float = 30.0
) -> tuple[int, ...]:
    """Fetch and parse the OEIS b-file through index ``max_steps``."""
    if max_steps < 0:
        raise ValueError("max_steps must be non-negative")
    request = Request(
        url, headers={"User-Agent": "ising3d-exact/0.1 (research reproducibility)"}
    )
    with urlopen(request, timeout=timeout) as response:  # noqa: S310 - caller URL
        text = response.read().decode("ascii")
    values: dict[int, int] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        index_text, value_text = line.split()[:2]
        index = int(index_text)
        if 0 <= index <= max_steps:
            values[index] = int(value_text)
    missing = [index for index in range(max_steps + 1) if index not in values]
    if missing:
        raise ValueError(f"OEIS b-file is missing indices {missing}")
    return tuple(values[index] for index in range(max_steps + 1))


def _direct_saw_counts(max_steps: int) -> tuple[int, ...]:
    counts = [0] * (max_steps + 1)
    visited = {(0, 0, 0)}

    def visit(point: tuple[int, int, int], steps: int) -> None:
        counts[steps] += 1
        if steps == max_steps:
            return
        x, y, z = point
        for dx, dy, dz in _DIRECTIONS:
            next_point = (x + dx, y + dy, z + dz)
            if next_point in visited:
                continue
            visited.add(next_point)
            visit(next_point, steps + 1)
            visited.remove(next_point)

    visit((0, 0, 0), 0)
    return tuple(counts)


def _transform_path(
    path: tuple[tuple[int, int, int], ...],
    symmetry: tuple[tuple[int, int, int], tuple[int, int, int]],
) -> tuple[tuple[int, int, int], ...]:
    permutation, signs = symmetry
    return tuple(
        (
            signs[0] * point[permutation[0]],
            signs[1] * point[permutation[1]],
            signs[2] * point[permutation[2]],
        )
        for point in path
    )


def _canonical_path(
    path: tuple[tuple[int, int, int], ...],
) -> tuple[tuple[int, int, int], ...]:
    return min(_transform_path(path, symmetry) for symmetry in _CUBIC_SYMMETRIES)


def _prefix_orbits(
    split_steps: int,
) -> tuple[tuple[tuple[tuple[int, int, int], ...], int], ...]:
    """Return one prefix per cubic-symmetry orbit and its exact orbit size."""
    orbit_sizes: dict[tuple[tuple[int, int, int], ...], int] = {}
    visited = {(0, 0, 0)}
    path = [(0, 0, 0)]

    def visit(steps: int) -> None:
        if steps == split_steps:
            canonical = _canonical_path(tuple(path))
            orbit_sizes[canonical] = orbit_sizes.get(canonical, 0) + 1
            return
        x, y, z = path[-1]
        for dx, dy, dz in _DIRECTIONS:
            next_point = (x + dx, y + dy, z + dz)
            if next_point in visited:
                continue
            visited.add(next_point)
            path.append(next_point)
            visit(steps + 1)
            path.pop()
            visited.remove(next_point)

    visit(0)
    return tuple(sorted(orbit_sizes.items()))


def _count_prefix_extensions(
    task: tuple[tuple[tuple[int, int, int], ...], int]
) -> tuple[int, ...]:
    path, max_steps = task
    split_steps = len(path) - 1
    counts = [0] * (max_steps - split_steps + 1)
    visited = set(path)

    def visit(point: tuple[int, int, int], steps: int) -> None:
        counts[steps - split_steps] += 1
        if steps == max_steps:
            return
        x, y, z = point
        for dx, dy, dz in _DIRECTIONS:
            next_point = (x + dx, y + dy, z + dz)
            if next_point in visited:
                continue
            visited.add(next_point)
            visit(next_point, steps + 1)
            visited.remove(next_point)

    visit(path[-1], split_steps)
    return tuple(counts)


def enumerate_saw_counts(
    max_steps: int, *, workers: int = 1, split_steps: int = 6
) -> tuple[int, ...]:
    """Count rooted ``n``-step SAWs on ``Z^3`` exactly for every ``n <= max_steps``.

    The search is ordinary backtracking with a visited-vertex set.  For a long
    search, prefixes are quotiented by all 48 signed coordinate permutations;
    each continuation count is multiplied by the *enumerated* orbit size.
    Distinct prefix orbits may be processed in independent worker processes.
    No tabulated count participates in the computation.
    """
    if max_steps < 0:
        raise ValueError("max_steps must be non-negative")
    if workers < 1:
        raise ValueError("workers must be positive")
    if max_steps <= split_steps or max_steps <= 6:
        return _direct_saw_counts(max_steps)

    split_steps = min(max(1, split_steps), max_steps)
    orbits = _prefix_orbits(split_steps)
    counts = list(_direct_saw_counts(split_steps - 1)) + [0] * (
        max_steps - split_steps + 1
    )
    tasks = tuple((path, max_steps) for path, _ in orbits)

    if workers == 1:
        extensions: Iterable[tuple[int, ...]] = map(_count_prefix_extensions, tasks)
    else:
        from concurrent.futures import ProcessPoolExecutor

        try:
            context = multiprocessing.get_context("fork")
        except ValueError:  # pragma: no cover - Windows has no fork context
            context = multiprocessing.get_context()
        executor = ProcessPoolExecutor(max_workers=workers, mp_context=context)
        chunk_size = max(1, len(tasks) // (workers * 8))
        extensions = executor.map(
            _count_prefix_extensions, tasks, chunksize=chunk_size
        )

    try:
        for (_, orbit_size), extension_counts in zip(orbits, extensions, strict=True):
            for offset, extension_count in enumerate(extension_counts):
                counts[split_steps + offset] += orbit_size * extension_count
    finally:
        if workers != 1:
            executor.shutdown()
    return tuple(counts)


def connective_constant_upper(count: int, steps: int, dps: int = 80) -> mp.mpf:
    """Return ``count**(1/steps)``, an upper bound on the connective constant."""
    _validate_count_and_precision(count, steps, dps)
    with mp.workdps(dps + 10):
        value = mp.power(mp.mpf(count), mp.mpf(1) / steps)
        return +value


def saw_kc_lower(count: int, steps: int, dps: int = 80) -> mp.mpf:
    """Return the SAW-correlation lower bound ``atanh(count**(-1/steps))``."""
    _validate_count_and_precision(count, steps, dps)
    with mp.workdps(dps + 10):
        value = mp.atanh(mp.power(mp.mpf(count), -mp.mpf(1) / steps))
        return +value


def square_lattice_kc(dps: int = 80) -> mp.mpf:
    """Return the exact square-lattice value ``log(1 + sqrt(2))/2``."""
    _validate_precision(dps)
    with mp.workdps(dps + 10):
        return +(mp.log(1 + mp.sqrt(2)) / 2)


def watson_integral_numeric(dimension: int = 3, dps: int = 80) -> mp.mpf:
    r"""Numerically evaluate ``I_d = integral 1/(d-sum(cos(k_i)))``.

    The measure is ``d^d k/(2*pi)^d``.  For ``d >= 3`` the Laplace/Bessel
    representation is integrated after ``t=(s/(1-s))**2``.  The scaled Bessel
    identity ``exp(-t) I_0(t) = 1F1(1/2;1;-2t)`` avoids catastrophic overflow.
    This function is a high-precision numerical cross-check, not the interval
    certificate returned by :func:`watson_i3_interval`.
    """
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("dimension must be a positive integer")
    _validate_precision(dps)
    if dimension <= 2:
        return mp.inf

    with mp.workdps(dps + 20):
        half = mp.mpf("0.5")
        endpoint = (
            2 / mp.power(2 * mp.pi, mp.mpf(3) / 2)
            if dimension == 3
            else mp.mpf(0)
        )

        def transformed_integrand(s: mp.mpf) -> mp.mpf:
            if s == 1:
                return endpoint
            t = mp.power(s / (1 - s), 2)
            scaled_i0 = mp.hyp1f1(half, 1, -2 * t)
            return mp.power(scaled_i0, dimension) * 2 * s / mp.power(1 - s, 3)

        value = mp.quad(
            transformed_integrand,
            [0, mp.mpf("0.5"), mp.mpf("0.8"), mp.mpf("0.95"), mp.mpf("0.995"), 1],
        )
        return +value


def watson_closed_form(dps: int = 80, *, denominator: int = 32) -> mp.mpf:
    r"""Evaluate ``sqrt(6)*prod Gamma(r/24)/(denominator*pi**3)``.

    Glasser--Zucker's simple-cubic Watson integral has ``denominator=32``.
    ``denominator=4`` is accepted solely to quantify the factor-eight error in
    the formula supplied for this research task.
    """
    _validate_precision(dps)
    if not isinstance(denominator, int) or denominator <= 0:
        raise ValueError("denominator must be a positive integer")
    with mp.workdps(dps + 20):
        product_gamma = mp.fprod(
            mp.gamma(mp.mpf(numerator) / 24) for numerator in (1, 5, 7, 11)
        )
        return +(mp.sqrt(6) * product_gamma / (denominator * mp.pi**3))


def watson_i3_interval(dps: int = 80) -> tuple[str, str]:
    r"""Certified interval for ``I_3 = W_sc/3`` via the exact gamma identity."""
    _validate_precision(dps)
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        interval = (
            mp.iv.sqrt(mp.iv.mpf(6))
            * mp.iv.gamma(mp.iv.mpf(1) / 24)
            * mp.iv.gamma(mp.iv.mpf(5) / 24)
            * mp.iv.gamma(mp.iv.mpf(7) / 24)
            * mp.iv.gamma(mp.iv.mpf(11) / 24)
            / (96 * mp.iv.pi**3)
        )
        return _interval_endpoints(interval)
    finally:
        mp.iv.dps = previous_dps


def connective_constant_upper_interval(
    count: int, steps: int, dps: int = 80
) -> tuple[str, str]:
    """Directed-rounding enclosure of ``count**(1/steps)``."""
    _validate_count_and_precision(count, steps, dps)
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        exponent = mp.iv.mpf(1) / steps
        interval = mp.iv.mpf(count) ** exponent
        return _interval_endpoints(interval)
    finally:
        mp.iv.dps = previous_dps


def saw_kc_lower_interval(
    count: int, steps: int, dps: int = 80
) -> tuple[str, str]:
    """Directed-rounding enclosure of the SAW lower-bound constant."""
    _validate_count_and_precision(count, steps, dps)
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        exponent = -mp.iv.mpf(1) / steps
        inverse_mu = mp.iv.mpf(count) ** exponent
        interval = mp.iv.log((1 + inverse_mu) / (1 - inverse_mu)) / 2
        return _interval_endpoints(interval)
    finally:
        mp.iv.dps = previous_dps


def square_lattice_kc_interval(dps: int = 80) -> tuple[str, str]:
    """Directed-rounding enclosure of ``log(1+sqrt(2))/2``."""
    _validate_precision(dps)
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        interval = mp.iv.log(1 + mp.iv.sqrt(mp.iv.mpf(2))) / 2
        return _interval_endpoints(interval)
    finally:
        mp.iv.dps = previous_dps


def infrared_kc_upper(dimension: int = 3, dps: int = 80) -> mp.mpf:
    """Return the reflection-positivity bound ``K_c <= I_d/2``."""
    value = watson_integral_numeric(dimension, dps)
    if mp.isinf(value):
        return value
    with mp.workdps(dps + 10):
        return +(value / 2)


def infrared_kc_upper_interval(dps: int = 80) -> tuple[str, str]:
    """Certified interval for the three-dimensional constant ``I_3/2``."""
    lo, hi = watson_i3_interval(dps)
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        return _interval_endpoints(mp.iv.mpf([lo, hi]) / 2)
    finally:
        mp.iv.dps = previous_dps


def _interval_endpoints(value: object) -> tuple[str, str]:
    text = str(value).strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError(f"unexpected mpmath interval representation: {text}")
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return lower.strip(), upper.strip()


def _validate_count_and_precision(count: int, steps: int, dps: int) -> None:
    if not isinstance(count, int) or count <= 0:
        raise ValueError("count must be a positive integer")
    if not isinstance(steps, int) or steps <= 0:
        raise ValueError("steps must be a positive integer")
    _validate_precision(dps)


def _validate_precision(dps: int) -> None:
    if not isinstance(dps, int) or dps < 30:
        raise ValueError("dps must be an integer of at least 30")


__all__ = [
    "OEIS_A001412",
    "OEIS_A001412_URL",
    "connective_constant_upper",
    "connective_constant_upper_interval",
    "enumerate_saw_counts",
    "fetch_oeis_a001412",
    "infrared_kc_upper",
    "infrared_kc_upper_interval",
    "saw_kc_lower",
    "saw_kc_lower_interval",
    "square_lattice_kc",
    "square_lattice_kc_interval",
    "watson_closed_form",
    "watson_i3_interval",
    "watson_integral_numeric",
]
