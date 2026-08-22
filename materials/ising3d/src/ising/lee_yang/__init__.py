r"""Exact finite-lattice Lee--Yang polynomials and high-precision zeros.

The field convention in this module is

``k = number of down spins``, ``x = exp(-2 K)``, and ``z = exp(-2 H_f)``.

If ``q`` is the number of unsatisfied bonds, the integer array returned by
:func:`field_dos_from_enumeration` or :func:`transfer_field_dos` obeys

.. math::

   Z(K,H_f) = e^{K n_b + H_f N}
       \sum_{k=0}^N \sum_{q=0}^{n_b} c[k,q] x^q z^k.

Thus a row ``c[k, :]`` is the exact integer coefficient polynomial of ``z^k``
in the broken-bond variable.  Spin reversal fixes ``q`` and maps ``k`` to
``N-k``, so the exact symmetry is ``c[k,q] = c[N-k,q]``.
"""

from __future__ import annotations

from math import prod

import mpmath as mp
import numpy as np

from ..exact_enumeration import joint_dos
from ..lattices import Lattice
from ..transfer_matrix import layer_bonds

__all__ = [
    "circle_residual",
    "field_dos_from_enumeration",
    "field_polynomial_at_K",
    "lee_yang_roots",
    "near_edge_density",
    "normalized_root_residual",
    "positive_zero_angles",
    "rational_field_polynomial",
    "rational_field_polynomial_transfer",
    "transfer_field_dos",
]


# Pairwise-distinct primes immediately below 2**30.  They were generated with
# exact integer primality testing.  Keeping every modulus below 2**30 makes a
# product or a sum of two products of residues safe in uint64 arithmetic.
_CRT_PRIMES = (
    1073741789,
    1073741783,
    1073741741,
    1073741723,
    1073741719,
    1073741717,
    1073741689,
    1073741671,
    1073741663,
    1073741651,
    1073741621,
    1073741567,
    1073741561,
    1073741527,
    1073741503,
    1073741477,
    1073741467,
    1073741441,
    1073741419,
    1073741399,
    1073741387,
    1073741381,
    1073741371,
    1073741329,
    1073741311,
    1073741309,
    1073741287,
    1073741237,
    1073741213,
    1073741197,
    1073741189,
    1073741173,
)
_MAX_TRANSFER_ENTRIES = 80_000_000
_MAX_RATIONAL_ENTRIES = 10_000_000


def field_dos_from_enumeration(lat: Lattice) -> np.ndarray:
    """Return exact ``c[k,q]`` using :func:`joint_dos`.

    The last axis of ``joint_dos`` counts up spins and its bond axes count
    satisfied bonds.  This function reverses both conventions to the down-spin
    exponent ``k`` and broken-bond exponent ``q`` used in the module docstring.
    ``joint_dos`` itself enforces the full-enumeration limit ``N <= 30``.
    """

    joint = joint_dos(lat, with_field=True)
    n_sites = lat.n_sites
    n_bonds = lat.n_bonds
    out = np.zeros((n_sites + 1, n_bonds + 1), dtype=np.int64)
    for satisfied_by_direction in np.ndindex(joint.shape[:-1]):
        q = n_bonds - sum(satisfied_by_direction)
        out[:, q] += joint[satisfied_by_direction][::-1]

    assert int(out.sum()) == 1 << n_sites
    assert np.array_equal(out, out[::-1]), "spin-reversal symmetry failed"
    return out


def _periodic_tuple(periodic, dim: int) -> tuple[bool, ...]:
    if isinstance(periodic, (bool, np.bool_)):
        return (bool(periodic),) * dim
    result = tuple(bool(value) for value in periodic)
    if len(result) != dim:
        raise ValueError("periodic must be a bool or have one entry per dimension")
    return result


def _layer_statistics(ns: int, bonds) -> tuple[np.ndarray, np.ndarray]:
    if ns > 31:
        raise ValueError("a dense layer state requires at most 31 sites")
    states = np.arange(1 << ns, dtype=np.uint32)
    broken = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for i, j in bonds:
        broken += (((states >> np.uint32(i)) ^ (states >> np.uint32(j))) & one).astype(
            np.int64
        )
    down = ns - np.bitwise_count(states).astype(np.int64)
    return broken, down


def _shift_field_and_bonds(
    values: np.ndarray, down_shift: np.ndarray, bond_shift: np.ndarray
) -> np.ndarray:
    """Shift the last two axes by state-dependent ``(k,q)`` exponents."""

    n_field = values.shape[-2]
    n_bond = values.shape[-1]
    out = np.zeros_like(values)
    pairs = set(zip(map(int, down_shift), map(int, bond_shift)))
    for dk, dq in pairs:
        selected = (down_shift == dk) & (bond_shift == dq)
        out[selected, ..., dk:, dq:] = values[
            selected, ..., : n_field - dk, : n_bond - dq
        ]
    return out


def _apply_interlayer_bivariate(values: np.ndarray, ns: int) -> np.ndarray:
    """Apply ``prod_i (I + x X_i)`` with the broken-bond axis last."""

    tail = values.shape[1:]
    n_bond = values.shape[-1]
    result = values
    for bit in range(ns):
        hi = 1 << (ns - bit - 1)
        lo = 1 << bit
        reshaped = result.reshape((hi, 2, lo) + tail)
        updated = reshaped.copy()
        updated[..., 1:] += reshaped[:, ::-1][..., : n_bond - 1]
        result = updated.reshape((1 << ns,) + tail)
    return result


def transfer_field_dos(shape, periodic=False) -> np.ndarray:
    """Return exact ``c[k,q]`` by a layer transfer matrix.

    This is the second-polynomial-axis extension of ``ising.transfer_matrix``;
    the shared module is not modified.  Shapes may be two- or three-dimensional,
    with the last component used as the transfer direction.  A fully periodic
    transfer direction is evaluated as an exact trace.  Dense memory is checked
    before allocation, and int64 arithmetic is used only for ``N <= 62``.
    """

    shape = tuple(int(side) for side in shape)
    if len(shape) not in (2, 3) or any(side < 1 for side in shape):
        raise ValueError("shape must be a two- or three-dimensional positive box")
    per = _periodic_tuple(periodic, len(shape))
    cross = shape[:-1]
    length = shape[-1]
    ns = prod(cross)
    n_sites = ns * length
    if n_sites > 62:
        raise ValueError("int64 field DOS requires N <= 62; use rational CRT propagation")

    in_layer = layer_bonds(cross, per[:-1])
    transfer_bonds = 0 if length == 1 else ns * (length if per[-1] else length - 1)
    n_bonds = len(in_layer) * length + transfer_bonds
    n_states = 1 << ns
    trace_factor = n_states if per[-1] and length > 1 else 1
    entries = n_states * trace_factor * (n_sites + 1) * (n_bonds + 1)
    if entries > _MAX_TRANSFER_ENTRIES:
        raise ValueError(
            f"field transfer would allocate {entries:,} int64 entries; "
            "use a smaller cross-section or rational CRT propagation"
        )

    broken, down = _layer_statistics(ns, in_layer)
    if not per[-1] or length == 1:
        values = np.zeros((n_states, n_sites + 1, n_bonds + 1), dtype=np.int64)
        values[:, 0, 0] = 1
        values = _shift_field_and_bonds(values, down, broken)
        for _ in range(1, length):
            values = _apply_interlayer_bivariate(values, ns)
            values = _shift_field_and_bonds(values, down, broken)
        result = values.sum(axis=0)
    else:
        values = np.zeros(
            (n_states, n_states, n_sites + 1, n_bonds + 1), dtype=np.int64
        )
        indices = np.arange(n_states)
        values[indices, indices, 0, 0] = 1
        for _ in range(length):
            values = _shift_field_and_bonds(values, down, broken)
            values = _apply_interlayer_bivariate(values, ns)
        result = values[indices, indices].sum(axis=0)

    assert result.min() >= 0, "int64 overflow"
    assert int(result.sum()) == 1 << n_sites
    assert np.array_equal(result, result[::-1]), "spin-reversal symmetry failed"
    return result


def rational_field_polynomial(
    field_dos: np.ndarray, numerator: int, denominator: int
) -> list[int]:
    """Evaluate an exact field DOS at rational ``x = numerator/denominator``.

    The returned integer coefficients are

    ``A[k] = sum_q c[k,q] numerator**q denominator**(n_b-q)``.

    They equal ``denominator**n_b`` times the coefficient polynomial in the
    module docstring, so the common clearing factor does not change any zero.
    """

    numerator = int(numerator)
    denominator = int(denominator)
    if numerator <= 0 or denominator <= 0:
        raise ValueError("numerator and denominator must be positive")
    values = np.asarray(field_dos)
    if values.ndim != 2:
        raise ValueError("field_dos must have axes (down spins, broken bonds)")
    n_bonds = values.shape[1] - 1
    bond_weights = [
        numerator**q * denominator ** (n_bonds - q) for q in range(n_bonds + 1)
    ]
    coefficients = [
        sum(int(count) * weight for count, weight in zip(row, bond_weights))
        for row in values
    ]
    if coefficients != coefficients[::-1]:
        raise AssertionError("spin-reversal palindrome lost during rational evaluation")
    return coefficients


def field_polynomial_at_K(field_dos: np.ndarray, K, dps: int = 100) -> list[mp.mpf]:
    """Evaluate ``sum_q c[k,q] exp(-2K q)`` with ``dps`` decimal digits."""

    if dps < 50:
        raise ValueError("Lee-Yang calculations require at least 50 decimal digits")
    mp.mp.dps = int(dps)
    coupling = mp.mpf(K)
    x = mp.exp(-2 * coupling)
    values = np.asarray(field_dos)
    if values.ndim != 2:
        raise ValueError("field_dos must have axes (down spins, broken bonds)")
    powers = [mp.mpf(1)]
    for _ in range(1, values.shape[1]):
        powers.append(powers[-1] * x)
    coefficients = [
        mp.fsum(mp.mpf(int(count)) * power for count, power in zip(row, powers))
        for row in values
    ]
    if coefficients != coefficients[::-1]:
        raise AssertionError("spin-reversal palindrome lost during numerical evaluation")
    return coefficients


def _modular_open_polynomial(
    shape: tuple[int, ...],
    per: tuple[bool, ...],
    numerator: int,
    denominator: int,
    modulus: int,
) -> list[int]:
    """One modular residue pass for rational open-boundary propagation."""

    cross = shape[:-1]
    length = shape[-1]
    ns = prod(cross)
    n_states = 1 << ns
    in_layer = layer_bonds(cross, per[:-1])
    broken, down = _layer_statistics(ns, in_layer)
    weights = np.array(
        [
            pow(numerator, int(q), modulus)
            * pow(denominator, len(in_layer) - int(q), modulus)
            % modulus
            for q in broken
        ],
        dtype=np.uint64,
    )
    mod = np.uint64(modulus)
    disagree = np.uint64(numerator % modulus)
    agree = np.uint64(denominator % modulus)

    values = np.zeros((n_states, ns + 1), dtype=np.uint64)
    values[np.arange(n_states), down] = weights
    for _ in range(1, length):
        width = values.shape[1]
        transformed = values
        for bit in range(ns):
            hi = 1 << (ns - bit - 1)
            lo = 1 << bit
            reshaped = transformed.reshape(hi, 2, lo, width)
            updated = np.empty_like(reshaped)
            lower = reshaped[:, 0]
            upper = reshaped[:, 1]
            updated[:, 0] = (agree * lower + disagree * upper) % mod
            updated[:, 1] = (disagree * lower + agree * upper) % mod
            transformed = updated.reshape(n_states, width)

        values = np.zeros((n_states, width + ns), dtype=np.uint64)
        for down_count in range(ns + 1):
            selected = down == down_count
            if np.any(selected):
                values[selected, down_count : down_count + width] = (
                    transformed[selected] * weights[selected, None]
                ) % mod

    # n_states <= 2**16 under the public memory guard, so summing residues
    # below 2**30 cannot overflow uint64 before the final reduction.
    totals = values.sum(axis=0, dtype=np.uint64) % mod
    return [int(value) for value in totals]


def rational_field_polynomial_transfer(
    shape,
    numerator: int,
    denominator: int,
    periodic=False,
) -> list[int]:
    """Exact integer field polynomial via modular transfer and CRT.

    This path avoids the ``N <= 30`` enumeration limit and the bond-polynomial
    axis.  It evaluates a rational ferromagnetic bond fugacity exactly while
    retaining the magnetization polynomial.  The transfer direction must be
    open; this keeps memory proportional to ``2**cross_section * (N+1)``.
    The 4x4x4 open cube therefore needs only about 34 MiB per residue pass.
    """

    shape = tuple(int(side) for side in shape)
    if len(shape) not in (2, 3) or any(side < 1 for side in shape):
        raise ValueError("shape must be a two- or three-dimensional positive box")
    per = _periodic_tuple(periodic, len(shape))
    if per[-1]:
        raise ValueError("rational CRT propagation currently requires an open transfer direction")
    numerator = int(numerator)
    denominator = int(denominator)
    if numerator <= 0 or denominator <= 0:
        raise ValueError("numerator and denominator must be positive")

    cross = shape[:-1]
    length = shape[-1]
    ns = prod(cross)
    n_states = 1 << ns
    n_sites = ns * length
    in_layer = layer_bonds(cross, per[:-1])
    n_bonds = len(in_layer) * length + ns * (length - 1)
    entries = n_states * (n_sites + 1)
    if ns > 16 or entries > _MAX_RATIONAL_ENTRIES:
        raise ValueError(
            f"rational transfer needs {entries:,} residues per final layer; "
            "the supported limit is a 16-site cross-section"
        )

    # Every coefficient is nonnegative and bounded by the total number of
    # configurations times the largest per-bond cleared weight.
    upper_bound = (1 << n_sites) * max(numerator, denominator) ** n_bonds
    coefficients = [0] * (n_sites + 1)
    combined_modulus = 1
    for modulus in _CRT_PRIMES:
        residues = _modular_open_polynomial(
            shape, per, numerator, denominator, modulus
        )
        inverse = pow(combined_modulus % modulus, -1, modulus)
        for k, residue in enumerate(residues):
            correction = ((residue - coefficients[k]) % modulus) * inverse % modulus
            coefficients[k] += combined_modulus * correction
        combined_modulus *= modulus
        if combined_modulus > upper_bound:
            break
    else:
        raise ValueError(
            "the exact coefficient bound exceeds the built-in CRT modulus product; "
            "use a smaller rational numerator and denominator"
        )

    assert all(0 <= value < upper_bound for value in coefficients)
    assert coefficients == coefficients[::-1], "spin-reversal palindrome failed after CRT"
    return coefficients


def _divide_by_z_plus_one(coefficients: list[mp.mpf]) -> list[mp.mpf]:
    # An odd-degree palindrome has z=-1 exactly.  Build only the lower half
    # by synthetic division and reflect it; running the recurrence through
    # every coefficient needlessly accumulates arbitrary-precision rounding.
    quotient = [mp.mpf(0)] * (len(coefficients) - 1)
    middle = (len(quotient) - 1) // 2
    quotient[0] = coefficients[0]
    for index in range(1, middle + 1):
        quotient[index] = coefficients[index] - quotient[index - 1]
    for index in range(middle):
        quotient[-1 - index] = quotient[index]

    reconstructed = (
        [quotient[0]]
        + [
            quotient[index - 1] + quotient[index]
            for index in range(1, len(quotient))
        ]
        + [quotient[-1]]
    )
    scale = max(abs(value) for value in coefficients)
    error = max(
        abs(expected - actual)
        for expected, actual in zip(coefficients, reconstructed)
    )
    if error > 8 * len(coefficients) * mp.eps * scale:
        raise AssertionError("odd palindromic polynomial did not have the root z=-1")
    return quotient


def _chebyshev_initial_roots(coefficients: list[mp.mpf]) -> np.ndarray:
    """Roots in ``t=cos(theta)`` from the palindromic half-polynomial."""

    degree = len(coefficients) - 1
    if degree % 2:
        raise ValueError("Chebyshev reduction requires even degree")
    half = degree // 2
    scale = max(abs(value) for value in coefficients)
    chebyshev = np.array(
        [float(coefficients[half] / scale)]
        + [float(2 * coefficients[half - j] / scale) for j in range(1, half + 1)],
        dtype=np.float64,
    )
    roots = np.polynomial.chebyshev.chebroots(chebyshev)
    roots = roots[np.argsort(roots.real)]
    # This is a detection gate, not a projection of arbitrary roots onto the
    # circle.  A clear complex or out-of-interval reduced root is reported.
    if any(abs(root.imag) > 1e-7 for root in roots):
        raise ArithmeticError("reduced field polynomial has non-real roots at float64 isolation")
    if any(root.real < -1.0000001 or root.real > 1.0000001 for root in roots):
        raise ArithmeticError("reduced field polynomial has a root outside [-1,1]")
    return roots.real


def _reduced_value_derivative(
    coefficients: list[mp.mpf], t: mp.mpf
) -> tuple[mp.mpf, mp.mpf]:
    """Evaluate the Chebyshev half-polynomial and its derivative."""

    half = (len(coefficients) - 1) // 2
    value = coefficients[half]
    derivative = mp.mpf(0)
    if half == 0:
        return value, derivative

    t_previous = mp.mpf(1)
    t_current = t
    u_previous = mp.mpf(1)
    u_current = 2 * t
    value += 2 * coefficients[half - 1] * t_current
    derivative += 2 * coefficients[half - 1]
    for order in range(2, half + 1):
        t_next = 2 * t * t_current - t_previous
        value += 2 * coefficients[half - order] * t_next
        derivative += 2 * coefficients[half - order] * order * u_current
        t_previous, t_current = t_current, t_next
        u_previous, u_current = u_current, 2 * t * u_current - u_previous
    return value, derivative


def _refine_reduced_root(
    coefficients: list[mp.mpf], initial: float, dps: int
) -> mp.mpf:
    scale = max(abs(value) for value in coefficients)
    normalized = [value / scale for value in coefficients]
    root = mp.mpf(repr(float(initial)))
    tolerance = mp.power(10, -dps + 15)
    for _ in range(50):
        value, derivative = _reduced_value_derivative(normalized, root)
        if derivative == 0:
            raise ArithmeticError("multiple reduced root defeated Newton refinement")
        step = value / derivative
        root -= step
        if abs(step) < tolerance:
            break
    else:
        raise ArithmeticError("reduced-root Newton iteration did not converge")
    value, _ = _reduced_value_derivative(normalized, root)
    if abs(value) > mp.power(10, -dps + 25):
        raise ArithmeticError("reduced-root residual is too large")
    if root < -1 or root > 1:
        raise ArithmeticError("high-precision reduced root lies outside [-1,1]")
    return root


def _poly_value_derivative(
    descending: list[mp.mpf], z: mp.mpc
) -> tuple[mp.mpc, mp.mpc]:
    value = mp.mpc(descending[0])
    derivative = mp.mpc(0)
    for coefficient in descending[1:]:
        derivative = derivative * z + value
        value = value * z + coefficient
    return value, derivative


def _complex_newton(
    normalized_descending: list[mp.mpf], initial: mp.mpc, dps: int
) -> mp.mpc:
    # A small radial perturbation makes the final radius an independent output
    # of complex Newton iteration rather than an identity from exp(i theta).
    perturbation = mp.power(10, -min(20, max(8, dps // 4)))
    root = mp.mpc(initial) * (1 + perturbation)
    tolerance = mp.power(10, -dps + 15)
    for _ in range(50):
        value, derivative = _poly_value_derivative(normalized_descending, root)
        if derivative == 0:
            raise ArithmeticError("multiple complex root defeated Newton refinement")
        step = value / derivative
        root -= step
        if abs(step) < tolerance:
            return root
    raise ArithmeticError("complex-root Newton iteration did not converge")


def lee_yang_roots(coefficients, dps: int = 100) -> list[mp.mpc]:
    """Return every field-polynomial root at ``dps >= 50`` digits.

    Palindromy reduces the problem to a real polynomial in
    ``t=(z+z**-1)/2=cos(theta)``.  Float64 Chebyshev roots are used only as
    isolating guesses; each is refined in arbitrary precision, converted to a
    unit-circle guess, radially perturbed, and independently refined against
    the original complex polynomial.  Consequently the reported radial error
    is measured rather than imposed by the parameterization.
    """

    if dps < 50:
        raise ValueError("Lee-Yang root finding requires at least 50 decimal digits")
    mp.mp.dps = int(dps)
    values = [mp.mpf(value) for value in coefficients]
    if len(values) < 2 or values[0] <= 0 or values[-1] <= 0:
        raise ValueError("coefficients must include positive constant and leading terms")
    if any(value < 0 for value in values):
        raise ValueError("ferromagnetic field-polynomial coefficients must be nonnegative")
    if values != values[::-1]:
        raise ValueError("field polynomial must be exactly palindromic")

    original = values
    odd = (len(values) - 1) % 2 == 1
    reduced = _divide_by_z_plus_one(values) if odd else values
    initial_reduced = _chebyshev_initial_roots(reduced)
    refined_t = [
        _refine_reduced_root(reduced, initial, dps) for initial in initial_reduced
    ]

    scale = max(original)
    descending = [value / scale for value in reversed(original)]
    roots: list[mp.mpc] = []
    for t in refined_t:
        theta = mp.acos(t)
        roots.append(_complex_newton(descending, mp.exp(1j * theta), dps))
        roots.append(_complex_newton(descending, mp.exp(-1j * theta), dps))
    if odd:
        roots.append(_complex_newton(descending, mp.mpc(-1), dps))
    roots.sort(key=lambda root: mp.arg(root))
    if len(roots) != len(original) - 1:
        raise AssertionError("root count does not equal polynomial degree")
    return roots


def circle_residual(roots) -> mp.mpf:
    """Return ``max_i ||z_i|-1|``."""

    values = list(roots)
    if not values:
        return mp.mpf(0)
    return max(abs(abs(root) - 1) for root in values)


def normalized_root_residual(coefficients, roots) -> mp.mpf:
    """Return the largest ``|P(z_i)| / sum_k |A_k|``."""

    values = [mp.mpf(value) for value in coefficients]
    scale = mp.fsum(abs(value) for value in values)
    descending = list(reversed(values))
    residuals = []
    for root in roots:
        value, _ = _poly_value_derivative(descending, root)
        residuals.append(abs(value) / scale)
    return max(residuals, default=mp.mpf(0))


def positive_zero_angles(roots) -> list[mp.mpf]:
    """Return positive zero angles in ``(0, pi]``, counting multiplicity."""

    tolerance = mp.power(10, -max(20, mp.mp.dps - 15))
    angles = sorted(mp.arg(root) for root in roots if mp.arg(root) > tolerance)
    return angles


def near_edge_density(roots, n_sites: int, intervals: int = 3) -> list[dict[str, mp.mpf]]:
    """Estimate the positive-arc zero density in the first few intervals.

    For consecutive positive angles, the empirical density is
    ``1 / (N * (theta[j+1]-theta[j]))`` at the interval midpoint.  It is an
    explicitly finite-volume diagnostic, not an infinite-volume exponent.
    """

    if n_sites <= 0:
        raise ValueError("n_sites must be positive")
    angles = positive_zero_angles(roots)
    count = min(max(0, int(intervals)), max(0, len(angles) - 1))
    result = []
    for index in range(count):
        spacing = angles[index + 1] - angles[index]
        result.append(
            {
                "midpoint": (angles[index + 1] + angles[index]) / 2,
                "density": 1 / (n_sites * spacing),
                "spacing": spacing,
            }
        )
    return result
