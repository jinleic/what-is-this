"""Self-similar blowup profile for the Cordoba-Cordoba-Fontelos (CCF) equation.

MODEL (Cordoba-Cordoba-Fontelos, J. Math. Pures Appl. 86 (2006) 529-540, eq (1)):
    f_t - (H f) f_x = 0,        H f(x) = (1/pi) PV int f(y)/(x-y) dy.

Vorticity form used by arXiv:2511.22819 eq (2), with omega = f_x, u = H f, u_x = H omega:
    omega_t - u omega_x = omega u_x,      u(0) = 0.

Self-similar ansatz (arXiv:2511.22819 eq (3)):
    omega = (1-t)^{-1} Omega(y),   y = x (1-t)^{-(1+lambda)},   U(y) = int_0^y H Omega.

DERIVATION (done independently here, then checked against the paper's eq (4)):
    omega_t = (1-t)^{-2} [ Omega + (1+lambda) y Omega' ]
    omega_x = (1-t)^{-(2+lambda)} Omega'
    u       = (1-t)^{lambda} U(y),  since H acts in x and commutes with the dilation
    u omega_x = (1-t)^{-2} U Omega'
    omega u_x = (1-t)^{-2} Omega (H Omega)
  giving the profile equation
    (P)   Omega + [ (1+lambda) y - U ] Omega' - Omega (H Omega) = 0
  which matches arXiv:2511.22819 eq (4). Good cross-check on conventions and signs.

EXACT FIRST INTEGRAL (derived here; the form actually solved below).
  Put G = (1+lambda) y - U, so G' = (1+lambda) - H Omega. Then
    (G Omega)' = G' Omega + G Omega'
               = (1+lambda) Omega - Omega H Omega + ( -Omega + Omega H Omega )   [by (P)]
               = lambda Omega.
  Hence, with F(y) = int_0^y Omega and F(0) = 0,
    (Q)   G(y) Omega(y) = lambda F(y),      G = (1+lambda) y - U.
  (Q) is algebraic in Omega: no derivative of the unknown appears. That removes the
  worst conditioning from the collocation system and is what we discretize.

SYMMETRY / NORMALIZATION.
  Parity: (P) is consistent only for Omega odd (then H Omega and F are even, U and G odd).
  The only continuous symmetry left is the dilation Omega(y) -> Omega(b y): under it
  Omega'(0) -> b Omega'(0), while lambda is untouched. (An amplitude rescale
  Omega -> a Omega is NOT a symmetry: it multiplies the nonlinear terms by a.)
  So fixing Omega'(0) = 1 removes the last degree of freedom and makes lambda isolated.

FAR FIELD.  Balancing Omega + (1+lambda) y Omega' ~ 0 gives
    Omega ~ C sgn(y) |y|^{-p},   p = 1/(1+lambda),
so we factor the decay out explicitly and solve for a bounded profile V:
    Omega(y) = y (1 + y^2)^{-(p+1)/2} V(y),   V even, V(0) = Omega'(0) = 1, V(inf) = C.

TARGET.  arXiv:2509.14185 Fig. 2 reports the stable CCF exponent
    lambda = 1.1807776628998
(consistent with Eggers-Fontelos, Nonlinearity 32 (2019), doi 10.1088/1361-6544/ab4e0b,
lambda ~ 1.18078). Recovering that number from an independent discretization is the
validation goal of this file.

STATUS / SCOPE.  This is an inviscid 1D model equation. It is NOT the Clay problem:
adding a constant-viscosity Laplacian and moving to boundary-free R^3 are both
unresolved conceptual gaps. See GAP_MAP.md.
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import hyp2f1

from hilbert import hilbert_line, make_grid


# ---------------------------------------------------------------- far-field models
# Antiderivatives of the two decay models, in closed form, so the slowly-convergent
# tails never go through quadrature. Both are exact.


def int_odd_model(y, p):
    """int_0^y s (1+s^2)^{-(p+1)/2} ds = [ (1+y^2)^{(1-p)/2} - 1 ] / (1-p)."""
    return ((1.0 + y**2) ** ((1.0 - p) / 2.0) - 1.0) / (1.0 - p)


def int_even_model(y, p):
    """int_0^y (1+s^2)^{-p/2} ds = y * 2F1(1/2, p/2; 3/2; -y^2)."""
    return y * hyp2f1(0.5, p / 2.0, 1.5, -(y**2))


def cumulative_theta(vals, y, n):
    """Cumulative int_0^y vals ds on the tan grid, antisymmetric about y=0.

    ds = (1+y^2)/2 dtheta, so we cumulate vals*(1+y^2)/2 in theta with the midpoint
    rule. Nodes are half-cell offset, so y=0 is not a node; we integrate outward from
    the centre of the grid, which sits exactly at y=0 by the antisymmetric pairing.
    """
    dtheta = 2.0 * np.pi / n
    w = vals * (1.0 + y**2) / 2.0 * dtheta
    half = n // 2
    out = np.empty_like(w)
    # positive side: y_{half} is the smallest positive node; midpoint rule from 0.
    out[half:] = np.cumsum(w[half:]) - 0.5 * w[half:]
    # negative side by parity of the integral (odd integrand -> even integral, etc.)
    out[:half] = -(np.cumsum(w[:half][::-1]) - 0.5 * w[:half][::-1])[::-1]
    return out


# ---------------------------------------------------------------- residual


def unpack(z, n):
    """z = [V on the positive half minus its first entry, lambda] -> (V full, lambda)."""
    half = n // 2
    lam = z[-1]
    vpos = np.empty(half)
    vpos[0] = 1.0  # normalization V(0)=1 pins the dilation symmetry
    vpos[1:] = z[:-1]
    v = np.concatenate([vpos[::-1], vpos])  # even extension on the paired grid
    return v, lam


def residual(z, n, y, half_t):
    v, lam = unpack(z, n)
    if not (0.05 < lam < 6.0):
        return np.full(n // 2, 1e6)
    p = 1.0 / (1.0 + lam)

    pref = y * (1.0 + y**2) ** (-(p + 1.0) / 2.0)  # odd, ~ sgn(y)|y|^{-p}
    om = pref * v  # Omega, odd
    hom = hilbert_line(om, half_t)  # H Omega, even

    # F = int_0^y Omega, split as (closed-form model) + (quadrature of remainder).
    c_f = v[-1]  # V(inf) read off at the outermost node
    f_rem = om - c_f * pref
    big_f = c_f * int_odd_model(y, p) + cumulative_theta(f_rem, y, n)

    # U = int_0^y H Omega, same split against the even model (1+y^2)^{-p/2}.
    c_u = hom[-1] * (1.0 + y[-1] ** 2) ** (p / 2.0)
    u_rem = hom - c_u * (1.0 + y**2) ** (-p / 2.0)
    big_u = c_u * int_even_model(y, p) + cumulative_theta(u_rem, y, n)

    # Algebraic profile equation (Q):  [ (1+lambda) y - U ] Omega - lambda F = 0.
    res = ((1.0 + lam) * y - big_u) * om - lam * big_f

    # (Q) is even; evaluate on the positive half only. Scale out the y^{1-p} growth
    # so the far field does not dominate the least-squares weighting.
    w = (1.0 + y**2) ** (-(1.0 - p) / 2.0)
    return (res * w)[n // 2 :]


def solve(n=512, lam0=1.2, verbose=True):
    theta, y, half_t = make_grid(n)
    z0 = np.concatenate([np.ones(n // 2 - 1), [lam0]])
    sol = least_squares(residual, z0, args=(n, y, half_t), xtol=1e-15,
                        ftol=1e-15, gtol=1e-15, method="lm", max_nfev=20000)
    v, lam = unpack(sol.x, n)
    r = residual(sol.x, n, y, half_t)
    if verbose:
        print(f"  N={n:5d}  lambda={lam:.13f}  |res|_inf={np.max(np.abs(r)):.3e}  "
              f"nfev={sol.nfev}")
    return lam, v, y, np.max(np.abs(r))


if __name__ == "__main__":
    print(__doc__.split("STATUS")[0].strip()[:0] or "", end="")
    print("CCF stable self-similar profile: independent spectral-collocation solve")
    print("reference lambda (arXiv:2509.14185 Fig.2) = 1.1807776628998\n")
    lam = 1.2
    for n in (128, 256, 512, 1024, 2048):
        lam, v, y, r = solve(n=n, lam0=lam)
        print(f"        error vs reference = {abs(lam - 1.1807776628998):.3e}")
