r"""Hilbert transform on R via the Cayley / tan-half-angle map, FFT-based.

Convention (Cordoba-Cordoba-Fontelos, J. Math. Pures Appl. 86 (2006), eq (2)):

    Hf(x) = (1/pi) PV \int_R f(y)/(x-y) dy,   \hat{Hf}(xi) = -i sgn(xi) \hat f(xi).

Map: y = tan(theta/2), theta in (-pi, pi). H commutes with dilation, so no scale enters.

IDENTITY.  f + i Hf is the boundary value of a function analytic in the upper half
plane; the Cayley transform carries this to the disc, so the line transform and the
periodic conjugate function agree up to a real additive constant:

    (Hf)(tan(theta/2)) = Hcirc[g](theta) + C[f],     g(theta) = f(tan(theta/2)),

with Hcirc[e^{ik theta}] = -i sgn(k) e^{ik theta}.  Fixing C by averaging in theta and
using \int (Hf) h = -\int f (Hh) with h = 1/(1+x^2), Hh = x/(1+x^2):

    C[f] = (1/pi) \int_R (Hf)(x) dx/(1+x^2) = -(1/pi) \int_R f(x) x/(1+x^2) dx
         = -(1/2pi) \int_{-pi}^{pi} g(phi) tan(phi/2) dphi,      dx/(1+x^2) = dtheta/2.

THE TRAP.  For f decaying only algebraically, f ~ c|x|^{-p} with 0 < p < 1, that last
integrand behaves like (pi-theta)^{p-1} at the endpoints: integrable, but singular.  A
uniform-grid rule therefore converges at only O(N^{-(1-p)}) -- measured O(N^{-1/2}) for
p ~ 0.4586.  The resulting error is an almost pure CONSTANT offset in Hf, which in the
CCF profile equation is indistinguishable from a shift in the eigenvalue lambda.  Test
data that decays like x^{-2} does not expose this.  Use `hilbert_line_tail` whenever the
data decays algebraically.
"""

import numpy as np
from scipy.special import beta as _beta


def hilbert_circle(g):
    """Periodic conjugate function on a uniform grid over [-pi, pi).

    Multiplier -i sgn(k); the k=0 and (for even N) Nyquist modes are annihilated.
    """
    n = g.shape[-1]
    ghat = np.fft.fft(g)
    k = np.fft.fftfreq(n, d=1.0 / n)
    mult = -1j * np.sign(k)
    if n % 2 == 0:
        mult[n // 2] = 0.0  # Nyquist carries no orientation
    return np.real(np.fft.ifft(ghat * mult))


def make_grid(n, scale=1.0):
    """Half-cell-offset uniform theta grid, the induced y grid, and tan(theta/2).

    theta_j = -pi + 2*pi*(j+1/2)/n keeps every node strictly inside (-pi, pi), so
    tan(theta/2) stays finite (|y| <~ 2n/pi).  A node exactly on -pi gives
    tan(theta/2) ~ 1e16 and the cancellation there wrecks the additive constant
    (observed: H^2+I error 2.6e-6 instead of 1e-16).
    The grid is antisymmetric-paired, theta_j + theta_{n-1-j} = 0, so parity
    constraints hold exactly on the grid.
    """
    theta = -np.pi + 2.0 * np.pi * (np.arange(n) + 0.5) / n
    half = np.tan(theta / 2.0)
    return theta, scale * half, half


def hilbert_line(g, half):
    """H on R for RAPIDLY DECAYING data (faster than |y|^{-1}).

    For algebraically decaying data use `hilbert_line_tail`.
    """
    n = g.shape[-1]
    return hilbert_circle(g) - np.sum(g * half) / n


def pow_tail(y, p):
    """Plain algebraic model m(y) = y (1+y^2)^{-(p+1)/2}, odd, ~ sgn(y)|y|^{-p}.

    Kept because its additive constant is exact (see `const_pow_tail`), but its own
    Hilbert transform is not elementary, so `odd_tail` below is the better model.
    """
    return y * (1.0 + y**2) ** (-(p + 1.0) / 2.0)


def const_pow_tail(p):
    """Exact C for `pow_tail`.

    C[m] = -(1/pi) \\int_R y^2 (1+y^2)^{-(p+3)/2} dy = -(1/pi) B(3/2, p/2),
    from \\int_0^inf x^{a-1}(1+x^2)^{-b} dx = (1/2) B(a/2, b-a/2) with a=3, b=(p+3)/2.
    """
    return -_beta(1.5, p / 2.0) / np.pi


def odd_tail(y, p):
    """Odd model tail with an EXACTLY known Hilbert transform.

    Psi(z) = (1-iz)^{-p} is analytic in the upper half plane and decays, so on R its
    boundary values satisfy Im Psi = H[Re Psi].  With 1-iy = sqrt(1+y^2)e^{-i arctan y},

        M_c(y) = (1+y^2)^{-p/2} cos(p arctan y)   (even)
        M_s(y) = (1+y^2)^{-p/2} sin(p arctan y)   (odd)

    give H[M_c] = M_s, hence H[M_s] = H[H[M_c]] = -M_c.

    M_s ~ sin(p*pi/2) sgn(y)|y|^{-p} as |y| -> inf: the same algebraic decay the CCF
    profile has.  Subtracting it removes the (pi-theta)^p corner that limits the FFT,
    not merely the additive constant.
    """
    r = (1.0 + y**2) ** (-p / 2.0)
    return r * np.sin(p * np.arctan(y))


def odd_tail_hilbert(y, p):
    """H[odd_tail] = -M_c, in closed form."""
    return -((1.0 + y**2) ** (-p / 2.0)) * np.cos(p * np.arctan(y))


def odd_tail_amplitude(p):
    """A in odd_tail(y,p) ~ A sgn(y)|y|^{-p}; A = sin(p*pi/2)."""
    return np.sin(p * np.pi / 2.0)


def hilbert_line_tail(g, half, y, p, c_tail):
    """H on R with the algebraic tail handled in closed form.

    g       : samples of f on the tan grid
    p       : decay exponent
    c_tail  : coefficient of odd_tail(y, p) in f, so f - c_tail*odd_tail decays
              faster than |y|^{-p}

    The exact pair carries the slow part; the FFT and the constant quadrature act only
    on the fast remainder.
    """
    n = g.shape[-1]
    rem = g - c_tail * odd_tail(y, p)
    h_rem = hilbert_circle(rem) - np.sum(rem * half) / n
    return c_tail * odd_tail_hilbert(y, p) + h_rem
