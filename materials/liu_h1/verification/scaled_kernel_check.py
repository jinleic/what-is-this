#!/usr/bin/env python3
"""Exact algebra for the scaled H1 corollary; no numerical parameter search.

Claim: for every kappa in [0,1], the residual for
z=(1-s)(1-t)(1+kappa*s*t) is NSD on finite signed measures.
This checker verifies only the new algebra. PSD closure, Taylor's theorem,
and the already-audited reciprocal-kernel lemma remain ordinary mathematics.
Run: nice -n 10 math/.venv/bin/python -I -B
     math/LIU_H1/verification/scaled_kernel_check.py
"""

import json
import resource
import signal

resource.setrlimit(resource.RLIMIT_CPU, (15, 16))
signal.alarm(30)

import sympy as sp


def main():
    a, b, k, r, s, t = sp.symbols("A B kappa r s t", nonnegative=True)
    g = lambda x: x + (1 - x) * sp.log(1 - x)
    phi = g(a * (1 + k * b)) + a * g(-k * b)
    denominator = (1 + k * b) * (1 - a * (1 + k * b))
    obligations = {
        "initial_value": sp.simplify(phi.subs(k, 0) - g(a)),
        "initial_derivative": sp.simplify(sp.diff(phi, k).subs(k, 0) + a * b * sp.log(1 - a)),
        "second_derivative": sp.cancel(sp.diff(phi, k, 2) - a * b**2 / denominator),
    }
    # Independently differentiate the reparameterized Taylor path. This catches
    # the easy-to-miss kappa^2 factor on the integral remainder.
    along_path = phi.subs(k, r * k)
    obligations["scaled_taylor_curvature"] = sp.cancel(
        sp.diff(along_path, r, 2) - k**2 * a * b**2 / denominator.subs(k, r * k)
    )
    av, bv = (1 - s) * (1 - t), s * t
    gap = 1 - av * (1 + r * k * bv)
    obligations["uniform_s_endpoint_bound"] = sp.expand(
        gap - s - (1 - s) * t * (1 - r * k * s * (1 - t))
    )
    obligations["uniform_t_endpoint_bound"] = sp.expand(
        gap - t - (1 - t) * s * (1 - r * k * t * (1 - s))
    )
    obligations["projected_remainder"] = sp.expand(av * (1 + k * bv) - k * av * bv - av)
    assert all(value == 0 for value in obligations.values()), obligations
    # A lost square must fail at a fixed rational point, not merely at a float.
    wrong_curvature = sp.diff(along_path, r, 2) - k * a * b**2 / denominator.subs(k, r * k)
    negative_control = sp.cancel(wrong_curvature).subs(
        {a: sp.Rational(1, 4), b: sp.Rational(1, 4), k: sp.Rational(16, 25), r: sp.Rational(1, 2)}
    )
    assert negative_control != 0
    print(json.dumps({
        "status": "PASS", "exact_identities": list(obligations),
        "lost_square_negative_control": str(negative_control),
        "domain": "kappa,r,s,t in [0,1]; calculus on the nonsingular interior, endpoints by the uniform bound",
        "scope": "Exact symbolic algebra only; continuum PSD conclusion uses the manuscript's stated analytic lemmas.",
        "sympy": sp.__version__,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
