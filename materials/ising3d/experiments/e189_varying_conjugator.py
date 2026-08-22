"""Exact counterexample to fixed-conjugator differentiation rigidity.

For the one-mode scalar spin group, the conjugacy saturation of its diagonal
torus contains every invertible diagonalizable 2x2 matrix.  The family
exp(a H) exp(b X) lies in that saturation on a nonempty complex-open set (and
for every real a,b), although H and X generate sl_2 rather than a conjugate of
the one-dimensional quadratic algebra.
"""

from __future__ import annotations

import sympy as sp


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _zero_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.factor(entry) == 0 for entry in matrix)


def run_witness() -> dict[str, object]:
    checks: list[dict[str, object]] = []
    r, u = sp.symbols("r u", nonzero=True)

    H = sp.Matrix([[1, 0], [0, -1]])
    X = sp.Matrix([[0, 1], [1, 0]])
    C = H * X - X * H
    generated_rank = sp.Matrix.hstack(
        H.reshape(4, 1), X.reshape(4, 1), C.reshape(4, 1)
    ).rank()
    _check(
        checks,
        "H and X generate a three-dimensional Lie algebra",
        generated_rank == 3,
        f"rank(H,X,[H,X])={generated_rank}",
    )
    _check(
        checks,
        "the generated basis closes as sl2",
        _zero_matrix(H * C - C * H - 4 * X)
        and _zero_matrix(X * C - C * X + 4 * H),
        "[H,[H,X]]=4X and [X,[H,X]]=-4H",
    )

    D = sp.diag(r, 1 / r)
    cosh_b = (u + 1 / u) / 2
    sinh_b = (u - 1 / u) / 2
    exp_bX = sp.Matrix([[cosh_b, sinh_b], [sinh_b, cosh_b]])
    M = sp.simplify(D * exp_bX)
    trace = sp.factor(sp.trace(M))
    determinant = sp.factor(M.det())
    discriminant = sp.factor(trace**2 - 4)
    expected_trace = sp.factor((r + 1 / r) * (u + 1 / u) / 2)
    _check(checks, "witness family has determinant one", determinant == 1, f"det={determinant}")
    _check(
        checks,
        "witness trace and discriminant are exact Laurent expressions",
        sp.factor(trace - expected_trace) == 0 and sp.factor(discriminant - (expected_trace**2 - 4)) == 0,
        f"trace={trace}; discriminant={discriminant}",
    )

    sample = sp.simplify(M.subs({r: sp.Rational(2), u: sp.Rational(2)}))
    sample_trace = sp.trace(sample)
    sample_discriminant = sp.factor(sample_trace**2 - 4)
    _check(
        checks,
        "non-diagonal regular-semisimple sample lies in the open set",
        sample[0, 1] != 0 and sample[1, 0] != 0 and sample.det() == 1 and sample_discriminant != 0,
        f"M(2,2)={sample.tolist()}, Delta={sample_discriminant}",
    )

    # At b=0 (u=1), the entire b-velocity is a conjugacy-orbit velocity.
    # If C(b) D C(b)^-1 is differentiated at C(0)=I, its orbit term is [Y,D].
    Y = sp.Matrix([[0, -(r**2) / (r**2 - 1)], [1 / (r**2 - 1), 0]])
    b_velocity = D * X
    absorbed = sp.simplify(Y * D - D * Y - b_velocity)
    _check(
        checks,
        "varying conjugator absorbs the full transverse derivative",
        _zero_matrix(absorbed),
        "[Y(r),diag(r,r^-1)]=diag(r,r^-1)X for r^2!=1",
    )
    Y_at_two = sp.simplify(Y.subs(r, sp.Rational(2)))
    _check(
        checks,
        "conjugator velocity is finite at a regular b=0 point",
        Y_at_two == sp.Matrix([[0, sp.Rational(-4, 3)], [sp.Rational(1, 3), 0]]),
        f"Y(2)={Y_at_two.tolist()}",
    )
    pole_numerator = sp.simplify((r**2 - 1) * Y)
    _check(
        checks,
        "the required conjugator velocity has a pole at the identity",
        pole_numerator.subs(r, 1) != sp.zeros(2),
        "(r^2-1)Y has nonzero limit, so no bounded differentiable choice crosses r=1",
    )

    data = {
        "tag": "[THEOREM]",
        "statement": (
            "The proposed rigidity implication is false: pointwise membership of "
            "exp(aA)exp(bB) in a union of conjugates on an open parameter set does not force "
            "Lie(A,B) into one fixed conjugate of the subgroup Lie algebra."
        ),
        "small_subgroup": {
            "ambient_group": "SL(2,C) on the one-mode Fock space C^2",
            "subgroup": "T=Spin(2,C) image={diag(r,r^-1): r in C^times}",
            "scalar_extension": "C^times T is the group of invertible diagonal 2x2 matrices",
            "conjugacy_saturation": "all invertible diagonalizable 2x2 matrices",
        },
        "generators": {
            "H": [[int(entry) for entry in row] for row in H.tolist()],
            "X": [[int(entry) for entry in row] for row in X.tolist()],
            "commutator": [[int(entry) for entry in row] for row in C.tolist()],
            "generated_lie_algebra": "sl(2,C), dimension 3",
            "quadratic_lie_algebra": "spin(2,C), dimension 1 (dimension 2 after adding scalars)",
        },
        "family": {
            "additive": "M(a,b)=exp(aH)exp(bX)",
            "multiplicative": "r=exp(a), u=exp(b)",
            "matrix": [[str(sp.factor(entry)) for entry in row] for row in M.tolist()],
            "trace": str(trace),
            "determinant": str(determinant),
            "discriminant": str(discriminant),
            "complex_open_set": "Omega={(r,u) in (C^times)^2 : Delta(r,u)!=0}",
            "real_open_statement": (
                "for real a,b, trace=2 cosh(a)cosh(b)>=2; equality only at (0,0), "
                "where M=I, so every real parameter gives a diagonalizable M"
            ),
            "pointwise_conclusion": (
                "on Omega, M has two distinct nonzero eigenvalues and is conjugate to T; "
                "the conjugator varies with (a,b)"
            ),
        },
        "differentiation_at_b_zero": {
            "tag": "[LEMMA]",
            "base": "D=diag(r,r^-1), r^2!=1",
            "derivative": "partial_b M|_{b=0}=D X",
            "conjugator_velocity": [[str(sp.factor(entry)) for entry in row] for row in Y.tolist()],
            "identity": "D X=[Y,D]",
            "interpretation": (
                "the entire B-direction is tangent to the conjugacy orbit and imposes no "
                "condition B in Lie(T); Y has a pole as r approaches 1"
            ),
        },
        "sample": {
            "r": "2",
            "u": "2",
            "matrix": [[str(entry) for entry in row] for row in sample.tolist()],
            "discriminant": str(sample_discriminant),
        },
    }
    return {"data": data, "checks": checks}


def main() -> int:
    result = run_witness()
    for item in result["checks"]:
        print(f"[{'PASS' if item['passed'] else 'FAIL'}] {item['name']}: {item['detail']}")
    if all(bool(item["passed"]) for item in result["checks"]):
        print("PASS e189 varying-conjugator obstruction")
        return 0
    print("FAIL e189 varying-conjugator obstruction")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
