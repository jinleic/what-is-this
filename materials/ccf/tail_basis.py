r"""The (1-iz)^{-q} basis: closed under H, d/dy and \int_0^y.

For q > 0, Psi_q(z) = (1-iz)^{-q} is analytic in the upper half plane and decays, so on
R its boundary values satisfy Im Psi = H[Re Psi].  Writing
1 - iy = sqrt(1+y^2) e^{-i arctan y},

    C_q(y) = (1+y^2)^{-q/2} cos(q arctan y)     even
    S_q(y) = (1+y^2)^{-q/2} sin(q arctan y)     odd

ALGEBRA (all exact; from d/dy (1-iy)^{-r} = i r (1-iy)^{-r-1}):

    H[C_q] = S_q                     H[S_q] = -C_q          (H^2 = -I)
    d/dy S_r = r C_{r+1}             d/dy C_r = -r S_{r+1}
    \int_0^y C_q = S_{q-1}/(q-1)     \int_0^y S_q = (1 - C_{q-1})/(q-1)
    C_q(0) = 1                       S_q(0) = 0             S_q'(0) = q

Far field, using arctan y = pi/2 - 1/y + O(y^-3) as y -> +inf:
    S_q ~ y^{-q} [ sin(q pi/2) - cos(q pi/2) q/y + ... ]
    C_q ~ y^{-q} [ cos(q pi/2) + sin(q pi/2) q/y + ... ]

The index may be negative in the antiderivative formulas, where S_{q-1} and C_{q-1}
grow like |y|^{1-q}.  That is the point: for the CCF profile Omega ~ |y|^{-p} forces
U = \int_0^y H(Omega) ~ |y|^{1-p}, and this basis carries that growth in closed form
instead of through a divergent quadrature.  Removing that divergence is why this module
exists.

Dilation: H[f(./s)](y) = (Hf)(y/s) for s > 0, and \int_0^y f(t/s) dt = s \int_0^{y/s} f.
The basis scale s threads through as explicit factors, handled by `Basis`.

SCOPE: this is numerical-method infrastructure for a 1D model equation.  It has no
bearing on the Clay problem.  See ../docs/GAP_MAP.md.
"""

import numpy as np


def S(q, y):
    """S_q(y) = Im (1-iy)^{-q} = (1+y^2)^{-q/2} sin(q arctan y). Odd in y."""
    return (1.0 + y * y) ** (-q / 2.0) * np.sin(q * np.arctan(y))


def C(q, y):
    """C_q(y) = Re (1-iy)^{-q} = (1+y^2)^{-q/2} cos(q arctan y). Even in y."""
    return (1.0 + y * y) ** (-q / 2.0) * np.cos(q * np.arctan(y))


class Basis:
    r"""Odd expansion Omega(y) = sum_k a_k S_{q_k}(y/s), with exact H, \int and d/dy.

    Every method returns an array of shape (len(q), len(y)); contract it with the
    coefficient vector to get the field.
    """

    def __init__(self, q, s=1.0):
        self.q = np.asarray(q, dtype=float)
        self.s = float(s)
        if np.any(np.abs(self.q - 1.0) < 1e-12):
            raise ValueError("q = 1 makes the antiderivative formulas singular")
        if np.any(self.q <= 0.0):
            raise ValueError("need q > 0 for analyticity and decay in the upper half plane")

    def _z(self, y):
        return np.asarray(y, dtype=float) / self.s

    def omega(self, y):
        """S_q(y/s)."""
        z = self._z(y)
        return np.stack([S(qk, z) for qk in self.q])

    def h_omega(self, y):
        """H[S_q(./s)](y) = -C_q(y/s)."""
        z = self._z(y)
        return np.stack([-C(qk, z) for qk in self.q])

    def big_u(self, y):
        r"""\int_0^y H[S_q(./s)] = -s S_{q-1}(y/s)/(q-1). Odd; grows like |y|^{1-q}."""
        z = self._z(y)
        return np.stack([-self.s * S(qk - 1.0, z) / (qk - 1.0) for qk in self.q])

    def big_f(self, y):
        r"""\int_0^y S_q(./s) = s (1 - C_{q-1}(y/s))/(q-1). Even; F(0) = 0."""
        z = self._z(y)
        return np.stack([self.s * (1.0 - C(qk - 1.0, z)) / (qk - 1.0) for qk in self.q])

    def omega_prime(self, y):
        """d/dy S_q(y/s) = (q/s) C_{q+1}(y/s)."""
        z = self._z(y)
        return np.stack([(qk / self.s) * C(qk + 1.0, z) for qk in self.q])

    def omega_prime_at_zero(self):
        """Omega'(0) = sum_k a_k q_k / s: the vector to dot with the coefficients."""
        return self.q / self.s
