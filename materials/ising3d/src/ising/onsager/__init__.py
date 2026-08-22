"""Onsager--Kaufman solution and its finite-dimensional Clifford control checks.

The finite-torus formula uses ``m`` rows of ``n`` spins.  ``Kx`` couples
neighbours inside a row and ``Ky`` couples neighbouring rows.  All public
thermodynamic functions return :class:`mpmath.mpf` values and use the caller's
``mp.mp.dps`` setting.

For a periodic row the Jordan--Wigner boundary bond is quadratic only after
fixing the spin-flip (fermion-parity) sector.  ``symmetric_row_transfer_matrix``
therefore distinguishes the physical spin matrix from either of its two
Gaussian sector extensions; see ``notes/onsager_derivation.md``.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
import math
from typing import Sequence

import mpmath as mp
import numpy as np
from scipy.linalg import expm, logm


__all__ = [
    "anisotropic_torus_broken_bond_counts",
    "collapse_bivariate_counts",
    "evaluate_bivariate_partition",
    "free_fermion_diagnostics",
    "kaufman_torus_Z",
    "majorana_operators",
    "onsager_free_energy",
    "reconstruct_torus_Z_from_generators",
    "symmetric_row_transfer_matrix",
]


_I2 = np.eye(2, dtype=np.complex128)
_X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
_Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
_PAULI = {"I": _I2, "X": _X, "Y": _Y, "Z": _Z}


def _positive_sizes(m: int, n: int) -> tuple[int, int]:
    m, n = int(m), int(n)
    if m < 2 or n < 2:
        raise ValueError("the finite-torus formulas require m,n >= 2")
    return m, n


def _nonnegative_couplings(Kx, Ky) -> tuple[mp.mpf, mp.mpf]:
    Kx, Ky = mp.mpf(Kx), mp.mpf(Ky)
    if Kx < 0 or Ky < 0:
        raise ValueError("this implementation uses ferromagnetic couplings Kx,Ky >= 0")
    return Kx, Ky


# ---------------------------------------------------------------------------
# Thermodynamic-limit determinant
# ---------------------------------------------------------------------------


def onsager_free_energy(Kx, Ky) -> mp.mpf:
    r"""Return the infinite-volume reduced free energy per spin.

    This evaluates

    .. math::

       \log 2 + {1\over 8\pi^2}\int_0^{2\pi}\!\int_0^{2\pi}
       \log(C_x C_y-S_x\cos\theta_x-S_y\cos\theta_y)
       \,d\theta_xd\theta_y,

    where ``C_i=cosh(2 K_i)`` and ``S_i=sinh(2 K_i)``.  The inner angular
    integral is performed analytically.  The remaining radicand is evaluated
    in a cancellation-free form, including on the critical line
    ``Sx*Sy == 1``.
    """

    Kx, Ky = _nonnegative_couplings(Kx, Ky)
    sx, sy = mp.sinh(2 * Kx), mp.sinh(2 * Ky)
    cx, cy = mp.cosh(2 * Kx), mp.cosh(2 * Ky)

    # cx*cy - sx - sy is non-negative and vanishes quadratically on the
    # critical line.  The quotient is the same expression without cancellation:
    # (cx*cy-sx-sy)(cx*cy+sx+sy) = (sx*sy-1)^2.
    determinant_gap = (sx * sy - 1) ** 2 / (cx * cy + sx + sy)

    def integrated_angle(theta: mp.mpf) -> mp.mpf:
        # A = cx*cy - sx*cos(theta), and A-sy is formed stably at theta=0.
        a_minus_sy = determinant_gap + 2 * sx * mp.sin(theta / 2) ** 2
        a = a_minus_sy + sy
        radical = mp.sqrt(a_minus_sy * (a_minus_sy + 2 * sy))
        return mp.log((a + radical) / 2)

    integral = mp.quad(integrated_angle, [0, mp.pi / 2, mp.pi])
    return mp.log(2) + integral / (2 * mp.pi)


# ---------------------------------------------------------------------------
# Kaufman's four spin-structure sectors
# ---------------------------------------------------------------------------


def _dual_coupling(K: mp.mpf) -> mp.mpf:
    """Kramers--Wannier dual: tanh(K*) = exp(-2 K), for K > 0."""

    return mp.atanh(mp.exp(-2 * K))


def _kaufman_energies(n: int, Kx: mp.mpf, Ky: mp.mpf) -> tuple[list[mp.mpf], list[mp.mpf]]:
    """Return antiperiodic (odd) and periodic (even) one-particle energies.

    ``gamma_0 = 2*(Ky* - Kx)`` is signed.  This single sign is what changes the
    periodic-sector ``sinh`` product across the critical line.
    """

    ky_star = _dual_coupling(Ky)
    cosh_2ky_star = mp.cosh(2 * ky_star)
    sinh_2ky_star = mp.sinh(2 * ky_star)
    cosh_2kx = mp.cosh(2 * Kx)
    sinh_2kx = mp.sinh(2 * Kx)

    gamma: list[mp.mpf] = []
    for ell in range(2 * n):
        if ell == 0:
            gamma.append(2 * (ky_star - Kx))
            continue
        momentum = mp.pi * ell / n
        argument = (
            cosh_2kx * cosh_2ky_star
            - sinh_2kx * sinh_2ky_star * mp.cos(momentum)
        )
        # Rounding can put an exactly unit argument infinitesimally below one.
        if argument < 1 and 1 - argument <= 32 * mp.eps:
            argument = mp.mpf(1)
        gamma.append(mp.acosh(argument))

    odd = [gamma[2 * r + 1] for r in range(n)]
    even = [gamma[2 * r] for r in range(n)]
    return odd, even


def _cycle_partition(length: int, K: mp.mpf) -> mp.mpf:
    return (2 * mp.cosh(K)) ** length + (2 * mp.sinh(K)) ** length


def kaufman_torus_Z(m: int, n: int, Kx, Ky) -> mp.mpf:
    """Return Kaufman's exact partition function for an ``m x n`` torus.

    The convention fixed against exact integer transfer matrices is

    ``Z = prefactor/2 * (C_odd + S_odd + C_even - S_even)``,

    with momenta ``(2r+1) pi/n`` in the odd sector and ``2r pi/n`` in the even
    sector.  Here ``prefactor=(2 sinh(2 Ky))**(m*n/2)`` and ``gamma_0`` in
    ``S_even`` is signed as ``2*(Ky* - Kx)``.
    """

    m, n = _positive_sizes(m, n)
    Kx, Ky = _nonnegative_couplings(Kx, Ky)

    # The dual-coupling representation is singular at zero coupling although
    # the physical limit is elementary: independent periodic Ising chains.
    if Ky == 0:
        return _cycle_partition(n, Kx) ** m
    if Kx == 0:
        return _cycle_partition(m, Ky) ** n

    odd, even = _kaufman_energies(n, Kx, Ky)
    odd_cosh = mp.fprod(2 * mp.cosh(m * value / 2) for value in odd)
    odd_sinh = mp.fprod(2 * mp.sinh(m * value / 2) for value in odd)
    even_cosh = mp.fprod(2 * mp.cosh(m * value / 2) for value in even)
    even_sinh = mp.fprod(2 * mp.sinh(m * value / 2) for value in even)
    prefactor = (2 * mp.sinh(2 * Ky)) ** (mp.mpf(m * n) / 2)
    return prefactor * (odd_cosh + odd_sinh + even_cosh - even_sinh) / 2


# ---------------------------------------------------------------------------
# Independent exact-integer anisotropic transfer matrix
# ---------------------------------------------------------------------------


def _horizontal_broken_counts(n: int) -> np.ndarray:
    return np.array(
        [
            sum(
                ((state >> site) & 1) != ((state >> ((site + 1) % n)) & 1)
                for site in range(n)
            )
            for state in range(1 << n)
        ],
        dtype=np.int64,
    )


def _shift_horizontal(vec: np.ndarray, shifts: np.ndarray, degree: int) -> np.ndarray:
    out = np.zeros_like(vec)
    for shift_value in np.unique(shifts):
        shift = int(shift_value)
        selected = shifts == shift
        out[selected, :, shift:, :] = vec[selected, :, : degree + 1 - shift, :]
    return out


def _apply_vertical_bonds(vec: np.ndarray, n: int, degree: int) -> np.ndarray:
    """Apply product_j (I + y X_j), with the y-degree on the final axis."""

    dim = 1 << n
    tail = vec.shape[1:]
    current = vec
    for site in range(n):
        high = 1 << (n - site - 1)
        low = 1 << site
        reshaped = current.reshape((high, 2, low) + tail)
        updated = reshaped.copy()
        updated[..., 1:] += reshaped[:, ::-1][..., :degree]
        current = updated.reshape((dim,) + tail)
    return current


@lru_cache(maxsize=None)
def anisotropic_torus_broken_bond_counts(m: int, n: int) -> tuple[tuple[int, ...], ...]:
    """Return exact bivariate broken-bond counts for an ``m x n`` torus.

    Entry ``[qx][qy]`` counts configurations with ``qx`` broken row bonds and
    ``qy`` broken inter-row bonds.  The implementation is an integer transfer
    matrix.  Summing entries on diagonals reproduces
    ``torus_broken_bond_poly((n,m))`` exactly.
    """

    m, n = _positive_sizes(m, n)
    sites = m * n
    if sites > 62:
        raise ValueError("int64 coefficient bound requires m*n <= 62")

    dim = 1 << n
    degree = sites
    horizontal = _horizontal_broken_counts(n)
    trace = np.zeros((degree + 1, degree + 1), dtype=np.int64)
    block = min(8, dim)

    for first in range(0, dim, block):
        stop = min(first + block, dim)
        width = stop - first
        vec = np.zeros((dim, width, degree + 1, degree + 1), dtype=np.int64)
        for initial in range(first, stop):
            vec[initial, initial - first, 0, 0] = 1

        for _ in range(m):
            vec = _shift_horizontal(vec, horizontal, degree)
            vec = _apply_vertical_bonds(vec, n, degree)
            if np.any(vec < 0):
                raise OverflowError("int64 overflow in anisotropic transfer matrix")

        for initial in range(first, stop):
            trace += vec[initial, initial - first]

    total = int(trace.sum())
    if total != 1 << sites:
        raise ArithmeticError(f"configuration count {total} != 2^{sites}")
    return tuple(tuple(int(value) for value in row) for row in trace)


def collapse_bivariate_counts(counts: Sequence[Sequence[int]]) -> list[int]:
    """Collapse ``c[qx][qy]`` to coefficients indexed by ``qx+qy``."""

    size = len(counts)
    if size == 0 or any(len(row) != size for row in counts):
        raise ValueError("counts must be a non-empty square array")
    collapsed = [0] * (2 * size - 1)
    for qx, row in enumerate(counts):
        for qy, coefficient in enumerate(row):
            collapsed[qx + qy] += int(coefficient)
    while len(collapsed) > 1 and collapsed[-1] == 0:
        collapsed.pop()
    return collapsed


def evaluate_bivariate_partition(counts: Sequence[Sequence[int]], Kx, Ky) -> mp.mpf:
    """Evaluate exact bivariate counts as a finite-torus partition function."""

    Kx, Ky = _nonnegative_couplings(Kx, Ky)
    size = len(counts)
    if size == 0 or any(len(row) != size for row in counts):
        raise ValueError("counts must be a non-empty square array")
    sites = size - 1
    x, y = mp.exp(-2 * Kx), mp.exp(-2 * Ky)
    polynomial = mp.fsum(
        mp.mpf(int(coefficient)) * x**qx * y**qy
        for qx, row in enumerate(counts)
        for qy, coefficient in enumerate(row)
        if coefficient
    )
    return mp.exp((Kx + Ky) * sites) * polynomial


# ---------------------------------------------------------------------------
# Explicit spin transfer matrix and Majorana quadratic generator
# ---------------------------------------------------------------------------


def _kron_all(operators: Sequence[np.ndarray]) -> np.ndarray:
    result = np.array([[1]], dtype=np.complex128)
    for operator in operators:
        result = np.kron(result, operator)
    return result


def _site_operator(n: int, site: int, operator: np.ndarray) -> np.ndarray:
    return _kron_all([operator if index == site else _I2 for index in range(n)])


@lru_cache(maxsize=None)
def majorana_operators(n: int) -> tuple[np.ndarray, ...]:
    """Return ``2n`` Jordan--Wigner Majoranas in the X-string convention.

    With zero-based ``j``, ``c[2j] = X_0...X_{j-1} Z_j`` and
    ``c[2j+1] = -X_0...X_{j-1} Y_j``.  Consequently
    ``X_j=-i c[2j]c[2j+1]`` and
    ``Z_j Z_{j+1}=-i c[2j+1]c[2j+2]``.
    """

    n = int(n)
    if n < 1:
        raise ValueError("n must be positive")
    operators: list[np.ndarray] = []
    for site in range(n):
        odd = [_X if index < site else _Z if index == site else _I2 for index in range(n)]
        even = [_X if index < site else _Y if index == site else _I2 for index in range(n)]
        operators.append(_kron_all(odd))
        operators.append(-_kron_all(even))
    return tuple(operators)


def symmetric_row_transfer_matrix(
    n: int,
    Kx,
    Ky,
    *,
    periodic: bool = True,
    parity_sector: int | None = None,
) -> np.ndarray:
    """Build ``V=V2**(1/2) V1 V2**(1/2)`` as a ``2**n`` square matrix.

    ``parity_sector=None`` gives the physical spin transfer matrix.  For a
    periodic row, ``parity_sector=+1`` or ``-1`` instead gives the Gaussian
    extension appropriate to that parity block.  Restricting the extension to
    its named parity reproduces the corresponding block of the physical matrix.
    """

    n = int(n)
    if n < 2:
        raise ValueError("n must be at least two")
    kx, ky = float(Kx), float(Ky)
    if kx < 0 or ky <= 0:
        raise ValueError("the explicit transfer matrix requires Kx >= 0 and Ky > 0")
    if parity_sector not in (None, -1, +1):
        raise ValueError("parity_sector must be None, +1, or -1")
    if parity_sector is not None and not periodic:
        raise ValueError("parity_sector is only meaningful for a periodic row")

    majoranas = majorana_operators(n)
    dim = 1 << n
    horizontal = np.zeros((dim, dim), dtype=np.complex128)

    if periodic and parity_sector is None:
        z_operators = [_site_operator(n, site, _Z) for site in range(n)]
        for site in range(n):
            horizontal += kx * (z_operators[site] @ z_operators[(site + 1) % n])
    else:
        for site in range(n - 1):
            horizontal += kx * (-1j * majoranas[2 * site + 1] @ majoranas[2 * site + 2])
        if periodic:
            horizontal += kx * parity_sector * (1j * majoranas[-1] @ majoranas[0])

    ky_star = float(_dual_coupling(mp.mpf(ky)))
    vertical = np.zeros_like(horizontal)
    for site in range(n):
        vertical += ky_star * (-1j * majoranas[2 * site] @ majoranas[2 * site + 1])

    log_scalar = n * math.log(2 * math.sinh(2 * ky)) / 2
    half_horizontal = expm(horizontal / 2)
    transfer = math.exp(log_scalar) * half_horizontal @ expm(vertical) @ half_horizontal
    # The exact ABA product is Hermitian positive definite.  Remove only the
    # roundoff-level anti-Hermitian part before applying the principal logarithm.
    return (transfer + transfer.conj().T) / 2


def _extract_quadratic_generator(
    transfer: np.ndarray, majoranas: Sequence[np.ndarray]
) -> tuple[complex, np.ndarray, np.ndarray, np.ndarray, float]:
    """Expand log(V)=alpha I + (i/4)c^T A c and return A and its energies."""

    matrix_log = logm(transfer)
    dim = transfer.shape[0]
    alpha = np.trace(matrix_log) / dim
    count = len(majoranas)
    generator = np.zeros((count, count), dtype=np.float64)
    reconstruction = alpha * np.eye(dim, dtype=np.complex128)

    for first in range(count):
        for second in range(first + 1, count):
            bilinear = 1j * majoranas[first] @ majoranas[second]
            coefficient = np.trace(bilinear @ matrix_log) / dim
            real_coefficient = float(coefficient.real)
            generator[first, second] = 2 * real_coefficient
            generator[second, first] = -2 * real_coefficient
            reconstruction += real_coefficient * bilinear

    eigvals = np.linalg.eigvalsh(1j * generator)
    energies = np.maximum(eigvals[count // 2 :].real, 0.0)
    residual = float(np.max(np.abs(matrix_log - reconstruction)))
    return alpha, generator, energies, matrix_log, residual


def _majorana_words(n: int) -> list[tuple[str, ...]]:
    words: list[tuple[str, ...]] = []
    for site in range(n):
        words.append(tuple("X" * site + "Z" + "I" * (n - site - 1)))
        words.append(tuple("X" * site + "Y" + "I" * (n - site - 1)))
    return words


def _allowed_bilinear_words(n: int) -> set[tuple[str, ...]]:
    multiplication = {
        ("I", "I"): "I",
        ("I", "X"): "X",
        ("I", "Y"): "Y",
        ("I", "Z"): "Z",
        ("X", "I"): "X",
        ("Y", "I"): "Y",
        ("Z", "I"): "Z",
        ("X", "X"): "I",
        ("Y", "Y"): "I",
        ("Z", "Z"): "I",
        ("X", "Y"): "Z",
        ("Y", "X"): "Z",
        ("Y", "Z"): "X",
        ("Z", "Y"): "X",
        ("Z", "X"): "Y",
        ("X", "Z"): "Y",
    }
    majorana_words = _majorana_words(n)
    allowed = {tuple("I" * n)}
    for first in range(2 * n):
        for second in range(first + 1, 2 * n):
            allowed.add(
                tuple(
                    multiplication[left, right]
                    for left, right in zip(majorana_words[first], majorana_words[second])
                )
            )
    return allowed


def _max_off_bilinear_coefficient(matrix_log: np.ndarray, n: int) -> tuple[float, str]:
    dim = 1 << n
    allowed = _allowed_bilinear_words(n)
    maximum = 0.0
    maximum_word = ""
    for word in product("IXYZ", repeat=n):
        if word in allowed:
            continue
        pauli_string = _kron_all([_PAULI[letter] for letter in word])
        coefficient = np.vdot(pauli_string, matrix_log) / dim
        magnitude = float(abs(coefficient))
        if magnitude > maximum:
            maximum = magnitude
            maximum_word = "".join(word)
    return maximum, maximum_word


def free_fermion_diagnostics(
    n: int,
    Kx,
    Ky,
    *,
    parity_sector: int | None = +1,
) -> dict[str, object]:
    """Diagnose Majorana-bilinear support and single-particle energies of ``log V``.

    A periodic physical spin matrix is a direct sum of two Gaussian operators,
    not one Gaussian on the full Jordan--Wigner Fock space.  Pass ``+1`` or
    ``-1`` (the default is ``+1``) for the appropriate Gaussian extension.
    Passing ``None`` intentionally diagnoses the unreduced physical matrix.
    """

    n = int(n)
    transfer = symmetric_row_transfer_matrix(
        n, Kx, Ky, periodic=True, parity_sector=parity_sector
    )
    majoranas = majorana_operators(n)
    alpha, generator, energies, matrix_log, reconstruction_residual = (
        _extract_quadratic_generator(transfer, majoranas)
    )
    off_coefficient, off_word = _max_off_bilinear_coefficient(matrix_log, n)

    analytic: list[float] = []
    energy_residual: float | None = None
    if parity_sector is not None:
        kx_mp, ky_mp = mp.mpf(Kx), mp.mpf(Ky)
        odd, even = _kaufman_energies(n, kx_mp, ky_mp)
        selected = odd if parity_sector == +1 else even
        analytic = sorted(float(abs(value)) for value in selected)
        energy_residual = max(abs(found - expected) for found, expected in zip(energies, analytic))

    identity_residual = 0.0
    parity = np.eye(1 << n, dtype=np.complex128)
    for site in range(n):
        spin_x = _site_operator(n, site, _X)
        parity = parity @ spin_x
        majorana_x = -1j * majoranas[2 * site] @ majoranas[2 * site + 1]
        identity_residual = max(identity_residual, float(np.max(np.abs(spin_x - majorana_x))))
        if site + 1 < n:
            spin_bond = _site_operator(n, site, _Z) @ _site_operator(n, site + 1, _Z)
            majorana_bond = -1j * majoranas[2 * site + 1] @ majoranas[2 * site + 2]
            identity_residual = max(
                identity_residual, float(np.max(np.abs(spin_bond - majorana_bond)))
            )
    wrap_bond = _site_operator(n, n - 1, _Z) @ _site_operator(n, 0, _Z)
    majorana_wrap = -1j * majoranas[-1] @ majoranas[0]
    boundary_identity_residual = float(np.max(np.abs(majorana_wrap + parity @ wrap_bond)))

    return {
        "n": n,
        "parity_sector": parity_sector,
        "log_scalar": float(alpha.real),
        "max_log_scalar_imaginary": float(abs(alpha.imag)),
        "single_particle_energies": [float(value) for value in energies],
        "analytic_energies": analytic,
        "max_energy_residual": energy_residual,
        "max_non_bilinear_coefficient": off_coefficient,
        "largest_non_bilinear_word": off_word,
        "max_quadratic_reconstruction_entry": reconstruction_residual,
        "max_generator_identity_residual": identity_residual,
        "boundary_parity_identity_residual": boundary_identity_residual,
        "generator_pfaffian": float(_pfaffian(generator)),
    }


def _pfaffian(matrix: np.ndarray) -> float:
    """Pfaffian by skew-Gaussian elimination (sufficient for the 2n<=12 control)."""

    work = np.array(matrix, dtype=np.float64, copy=True)
    size = work.shape[0]
    if work.shape != (size, size) or size % 2:
        raise ValueError("a Pfaffian requires an even-dimensional square matrix")
    value = 1.0

    for first in range(0, size - 1, 2):
        pivot_index = first + 1 + int(np.argmax(np.abs(work[first, first + 1 :])))
        if abs(work[first, pivot_index]) < 1e-15:
            return 0.0
        if pivot_index != first + 1:
            work[[first + 1, pivot_index], :] = work[[pivot_index, first + 1], :]
            work[:, [first + 1, pivot_index]] = work[:, [pivot_index, first + 1]]
            value = -value

        pivot = work[first, first + 1]
        value *= pivot
        for row in range(first + 2, size):
            for column in range(row + 1, size):
                work[row, column] += (
                    work[first + 1, row] * work[first, column]
                    - work[first, row] * work[first + 1, column]
                ) / pivot
                work[column, row] = -work[row, column]
    return float(value)


def reconstruct_torus_Z_from_generators(m: int, n: int, Kx, Ky) -> dict[str, object]:
    """Reconstruct the torus partition function from the two ``2n`` generators.

    The ordinary Gaussian trace is a product of ``2 cosh(m gamma/2)``.  The
    parity-inserted trace is ``(-1)^n sign(Pf A)`` times the corresponding
    ``sinh`` product.  Projecting with ``(1+p P)/2`` and summing ``p=+1,-1``
    gives the physical spin partition function.
    """

    m, n = _positive_sizes(m, n)
    majoranas = majorana_operators(n)
    total = 0.0
    sectors: list[dict[str, object]] = []

    for parity in (+1, -1):
        transfer = symmetric_row_transfer_matrix(
            n, Kx, Ky, periodic=True, parity_sector=parity
        )
        alpha, generator, energies, _, residual = _extract_quadratic_generator(
            transfer, majoranas
        )
        scalar = math.exp(m * float(alpha.real))
        ordinary_trace = scalar * float(np.prod(2 * np.cosh(m * energies / 2)))
        pfaffian = _pfaffian(generator)
        sinh_product = scalar * float(np.prod(2 * np.sinh(m * energies / 2)))
        if abs(pfaffian) < 1e-14 or sinh_product == 0:
            parity_trace = 0.0
        else:
            parity_trace = ((-1) ** n) * math.copysign(sinh_product, pfaffian)
        contribution = (ordinary_trace + parity * parity_trace) / 2
        total += contribution
        sectors.append(
            {
                "parity": parity,
                "log_scalar": float(alpha.real),
                "energies": [float(value) for value in energies],
                "pfaffian": pfaffian,
                "ordinary_trace": ordinary_trace,
                "parity_inserted_trace": parity_trace,
                "projected_trace": contribution,
                "quadratic_residual": residual,
            }
        )

    return {"Z": total, "sectors": sectors}
