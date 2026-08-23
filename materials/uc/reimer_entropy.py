"""REIMER + ENTROPY: a legitimate new constraint, a tempting new ceiling,
and a sharp impossibility result for coordinatewise use.

STATUS / HEADLINE
=================
PROVED: Reimer's average-set-size theorem gives a valid GLOBAL constraint on
the conditional-probability laws in the Gilmer/Sawin/Cambie entropy argument:

    sum_i E[p_i] = E|A| >= (1/2) log_2|F|
                 = (1/2) sum_i E[h(p_i)],

hence sum_i (2 M_i - L_i) >= 0.  This is a genuine conditioning-stable
constraint, unlike the refuted "support p <= c" idea (p_i is conditional and
can equal 1 even when its mean/frequency is below c).

NUMERICAL / EXACT-EQUATION DISCOVERY: if one incorrectly applies that global
constraint to a SINGLE common law mu, the old c* obstruction is excluded and
the new two-point obstruction is

    mu_R = a_R delta_1 + (1-a_R) delta_{b_R},

where Reimer equality and the entropy functional are both tight:

    (1-a) h(b) = 2[a+(1-a)b],
    (1-alpha)(1-a)^2 h(2b-b^2) + alpha(1-2a) - (1-a)h(b) = 0.

At alpha = 0.0356069 this gives

    a_R = 0.1249081993926099066101029511...
    b_R = 0.2944850724339305712171116825...
    c_R = a_R + (1-a_R)b_R
        = 0.3826096716808159068990890908...,

which exceeds the two-strategy ceiling c* by 2.641383141e-4.  A KKT penalty
lambda_R = 0.00883449108362051295944... makes this obstruction stationary,
and numerical searches over 1--4 pair-orbits return only it.

PROVED REFUTATION: this DOES NOT yield a stronger union-closed theorem by
itself.  Cambie's proof applies the probability-distribution inequality
coordinate by coordinate to the conditional law mu_i of
p_i=Pr(A_i=1|A_{<i}); Reimer constrains only the SUM over coordinates.  The
single-law KKT inequality

    F(mu) - lambda [2M(mu)-L(mu)] >= 0

already fails on every deterministic conditional law

    nu_d = a delta_1 + (1-a) delta_0,   0<a<=c:
    L(nu_d)=Q(nu_d)=C(nu_d)=F(nu_d)=0,
    2M(nu_d)-L(nu_d)=2a>0,

so its left side is -2 lambda a < 0.  This is not a technical defect: such
zero-entropy deterministic coordinates carry exactly the positive Reimer
slack that can compensate entropy-active bad coordinates without contributing
any union-entropy gain.

More strongly, the old c* obstruction mu_* has

    F(mu_*) = 0,
    L(mu_*) - 2M(mu_*) = d = 0.07755434785474821088... > 0.

Choose deterministic coordinates with mean a <= c.  Their Reimer surplus is
2a and their F is zero.  Taking the coordinate-count ratio

    N_det / N_bad = d/(2a)

makes the GLOBAL Reimer inequality an equality while sum_i F(mu_i)=0.  At
a=c_R this ratio is 0.1013491733155...: one deterministic coordinate balances
about 9.8669 bad coordinates.  Both coordinate-law means (c* and c_R) are at
most every target c >= c_R.

CONCLUSION (PROVED at the level of the entropy relaxation): Reimer's scalar
average-size inequality alone cannot break c* in a coordinatewise
Sawin/Cambie proof.  Any genuine gain must use information that couples the
conditional laws across coordinates -- e.g. Reimer's stronger structural
non-interference conditions, realizability of the conditional-law sequence,
or a new multi-coordinate entropy inequality.  The tempting c_R is the
ceiling of an illegitimate single-law relaxation, not a theorem target.

Run: ./.venv/bin/python uc/reimer_entropy.py
"""

from mpmath import mp

mp.dps = 80
ALPHA = mp.mpf("0.0356069")
C_STAR = mp.mpf("0.38234553336670272115")
B_STAR = mp.mpf("0.32945473850303697239")


def h(x):
    x = mp.mpf(x)
    if x <= 0 or x >= 1:
        return mp.mpf(0)
    return -(x * mp.log(x) + (1 - x) * mp.log(1 - x)) / mp.log(2)


def f_two_point(a, b):
    """The pair-orbit functional on a delta_1 + (1-a) delta_b."""
    return ((1 - ALPHA) * (1 - a) ** 2 * h(2 * b - b * b)
            + ALPHA * (1 - 2 * a) - (1 - a) * h(b))


def reimer_slack(a, b):
    """2M-L; Reimer admissibility is nonnegative."""
    mean = a + (1 - a) * b
    entropy = (1 - a) * h(b)
    return 2 * mean - entropy


def solve_new_root():
    equations = (
        lambda a, b: reimer_slack(a, b),
        lambda a, b: f_two_point(a, b),
    )
    return mp.findroot(equations, (mp.mpf("0.125"), mp.mpf("0.2945")),
                       solver="mdnewton", tol=mp.mpf("1e-70"))


def fixed_mean_a(c, b):
    return (c - b) / (1 - b)


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    # Old obstruction.
    a_star = 1 - 1 / (2 - h(B_STAR))
    mean_star = a_star + (1 - a_star) * B_STAR
    entropy_star = (1 - a_star) * h(B_STAR)
    deficit = entropy_star - 2 * mean_star
    assert abs(mean_star - C_STAR) < mp.mpf("2e-19")
    assert abs(f_two_point(a_star, B_STAR)) < mp.mpf("2e-19")
    assert deficit > 0
    print("PROVED/80-digit control: old c* obstruction violates Reimer")
    print("  mean M    =", mp.nstr(mean_star, 40))
    print("  entropy L =", mp.nstr(entropy_star, 40))
    print("  L-2M      =", mp.nstr(deficit, 40), "> 0")
    print()

    # Tempting single-law root.
    a_r, b_r = solve_new_root()
    c_r = a_r + (1 - a_r) * b_r
    assert abs(reimer_slack(a_r, b_r)) < mp.mpf("1e-70")
    assert abs(f_two_point(a_r, b_r)) < mp.mpf("1e-70")
    assert c_r > C_STAR
    print("NUMERICAL (80-digit root of exact equations): tempting single-law root")
    print("  a_R =", mp.nstr(a_r, 60))
    print("  b_R =", mp.nstr(b_r, 60))
    print("  c_R =", mp.nstr(c_r, 60))
    print("  c_R-c* =", mp.nstr(c_r - C_STAR, 35))
    print("  Reimer residual =", mp.nstr(reimer_slack(a_r, b_r), 8))
    print("  F residual       =", mp.nstr(f_two_point(a_r, b_r), 8))
    print()

    # KKT multiplier at fixed mean.
    def f_at_fixed_mean(b):
        return f_two_point(fixed_mean_a(c_r, b), b)

    def s_at_fixed_mean(b):
        return reimer_slack(fixed_mean_a(c_r, b), b)

    lam = mp.diff(f_at_fixed_mean, b_r) / mp.diff(s_at_fixed_mean, b_r)
    penalized = lambda b: f_at_fixed_mean(b) - lam * s_at_fixed_mean(b)
    assert lam > 0
    assert abs(penalized(b_r)) < mp.mpf("1e-70")
    assert abs(mp.diff(penalized, b_r)) < mp.mpf("1e-60")
    assert mp.diff(penalized, b_r, 2) > 0
    print("NUMERICAL (80-digit derivatives): tempting KKT multiplier")
    print("  lambda_R =", mp.nstr(lam, 60))
    print("  H'' at obstruction =", mp.nstr(mp.diff(penalized, b_r, 2), 40), "> 0")
    print()

    # Universal refutation of the local penalty.
    a_det = c_r
    deterministic_penalty = -lam * 2 * a_det
    assert deterministic_penalty < 0
    ratio = deficit / (2 * a_det)
    print("PROVED: deterministic conditional laws refute the local KKT inequality")
    print("  nu_d = a delta_1 + (1-a) delta_0, a=c_R")
    print("  F(nu_d)=0, L(nu_d)=0, 2M(nu_d)=", mp.nstr(2 * a_det, 40))
    print("  F-lambda(2M-L) =", mp.nstr(deterministic_penalty, 40), "< 0")
    print()
    print("PROVED: balanced coordinate-law sequence defeats scalar Reimer use")
    print("  bad-law deficit L-2M =", mp.nstr(deficit, 40))
    print("  deterministic/bad coordinate ratio =", mp.nstr(ratio, 40))
    print("  one deterministic coordinate balances", mp.nstr(1 / ratio, 12),
          "bad coordinates")
    print("  total Reimer slack = 0 and total entropy-functional gain = 0")
    print()
    print("FINAL VERDICT: the Reimer average-size SCALAR alone cannot break c*.")
    print("A stronger cross-coordinate/realizability constraint is required.")
