#!/usr/bin/env python3
"""Reproducible checker for the single-pair channel split rate (-1.5926).

Configuration: P0 = delta(x*) + delta(x*+d) with shared masses (1/2, 1/2),
P1 = delta(x*) + delta(x*-d) with the SAME shared masses (1/2, 1/2);
paired law (x_1, x_2) = (x*, x*+d) masses (1/2, 1/2) on P0 and
(y_1, y_2) = (x*, x*-d) with the SAME masses (1/2, 1/2) on P1.

Functional (from the c-channel bilinear):
    E_K = sum_{i,j<=2} a_i a_j [K(x_i, x_j) + K(y_i, y_j) - 2K(x_i, y_j)]
with K(s, t) = h(pi(s, t)).

Result (as printed by this file):
    E_K/d^2 = -1.5926157 at d = 1/512, approaching that value from
    -1.6102991 at d = 1/10 monotonically in d; relative spread over the
    whole grid 1.110e-2.  The converged rate is the one quoted.

Rate normalisation, resolved: an earlier note quoted 0.1593 for this
object.  That was the two-slot mass-normalised rate, i.e. 1.5926/10, and
it is superseded.  See Addendum 10 of
uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md, which records the reconciliation:
the rate is -1.5926 and the ratio to the smooth-chart constant is
|rate|/kappa = 12.97, not 1.297.  The mathematics of the rate itself was
never in question; only the normalisation and the stale decimals were.

Consequence: |rate|/kappa = 12.97 >> 1, so the smooth-chart kappa margin
does not dominate this channel rate; a closing argument needs the pinch
band's higher-order behaviour or a per-fiber corner.  Every number in the
CONCLUSION block below is printed from the computed values, so the file
cannot go stale again.

Run: math/.venv/bin/python -I -B math/uc/probe_single_pair_ek.py
"""
import sys
sys.path.insert(0, 'uc')
import mpmath
from fractions import Fraction
from liu9_binding import solve_equation_parameters


def h_mp(u):
    if u == 0 or u == 1:
        return mpmath.mpf(0)
    return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))


def _pi_mp(x, y):
    return x * y * (1 + (1 - x) * (1 - y))


def K(s, t):
    return h_mp(_pi_mp(s, t))


D_GRID = (Fraction(1, 10), Fraction(1, 32), Fraction(1, 128),
          Fraction(1, 256), Fraction(1, 512))


def main():
    params = solve_equation_parameters(200)
    with mpmath.workdps(200):
        beta, m, xstar = params.beta, params.mean, params.x
        print('CONFIGURATION:')
        print(' P0 = delta(x*) + delta(x*+d); masses (1/2, 1/2)')
        print(' P1 = delta(x*) + delta(x*-d); masses (1/2, 1/2)')
        print(' K(s,t) = h(pi(s,t))')
        print()
        rows = []
        for d in D_GRID:
            xm = [xstar, xstar + d]
            ym = [xstar, xstar - d]
            a = [Fraction(1, 2), Fraction(1, 2)]
            E_K = 0
            for i in range(2):
                for j in range(2):
                    E_K += a[i]*a[j]*(K(xm[i], xm[j]) + K(ym[i], ym[j]) - 2*K(xm[i], ym[j]))
            rate2 = E_K / (d*d)
            rate4 = E_K / (d**4)
            rows.append((d, E_K, rate2, rate4))
            print('d = %s:  E_K = %s,  E_K/d^2 = %s,  E_K/d^4 = %s' % (mpmath.nstr(d, 10), mpmath.nstr(E_K, 14), mpmath.nstr(rate2, 14), mpmath.nstr(rate4, 14)))
        rates = [float(r[2]) for r in rows]
        target = rates[-1]
        spread = max(abs(r - target)/abs(target) for r in rates)
        print()
        print('rate spread relative to smallest-d rate: %.3e' % spread)
        kappa = mpmath.mpf(4119063) / mpmath.mpf(33554432)
        print('kappa = %s' % mpmath.nstr(kappa, 16))
        print('|rate| = %s' % mpmath.nstr(abs(target), 14))
        print('|rate|/kappa = %s' % mpmath.nstr(abs(target)/kappa, 12))
        print()
        ratio = abs(target) / kappa
        print('CONCLUSION: converged rate E_K/d^2 = %.10g at d = %s'
              % (target, rows[-1][0]))
        print('channel damage ratio |rate|/kappa = %s: %s than the kappa'
              % (mpmath.nstr(ratio, 6),
                 'LARGER' if ratio > 1 else 'SMALLER'))
        print('F-margin rate, so the smooth-chart margin does NOT dominate')
        print('generally; the induction closing requires the pinch-band')
        print('higher-order behaviour or a per-fiber corner.')
        print('(The 0.1593 of an earlier note was this rate divided by the')
        print(' two-slot mass normalisation; see Addendum 10 of')
        print(' uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md.)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
