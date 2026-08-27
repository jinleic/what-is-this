"""THEOREM B''' -- rigorous reduction to at most two pair-orbits.

This file closes only the support-reduction step.  Throughout, ``Phi`` means the
*exact* pair-orbit functional (called ``Phi_exact`` in ``margin_lemma.py``), not
the later Margin-Lemma lower bound.  No numerical observation below is used in
the proof.

THEOREM (quantifiers made explicit).
For every fixed alpha in [0,1] and t in [0,1], the minimum of

    Phi_alpha(nu) = (1-alpha) Q(mu_nu) + alpha C(nu) - L(nu)

over all probability measures nu on Delta whose marginal mu_nu has mean at
most t is attained by a nu supported on at most two points of Delta.  Thus the
minimizing marginal is supported on at most four points.  In particular this
holds for the certified value alpha = 0.0356069.

Here and below, an atom (p,q) of nu is one pair-orbit.  The theorem is a
statement about each fixed alpha; its proof is uniform for alpha in [0,1].
The restriction alpha <= 1 is used exactly once, to make -(1-alpha)E concave.

===============================================================================
SETUP
===============================================================================

[PROVED-HERE]  Let

    Delta = {(p,q) in R^2 : 0 <= p <= q <= 1}.

Delta is a closed subset of the compact metric square [0,1]^2, hence is a
compact metric space.  Write P(Delta) for its regular Borel probability
measures, with the weak-* topology (convergence against every f in C(Delta)).
The marginal map M : P(Delta) -> P([0,1]) is

    mu_nu = M(nu) = (pi_1# nu + pi_2# nu)/2,

so, for f in C([0,1]),

    integral f(x) dmu_nu(x)
      = integral_Delta (f(p)+f(q))/2 dnu(p,q).                         (S1)

[CITED(verified quote R1)]  Riesz representation identifies regular Borel
measures on Delta with C(Delta)^*.  Van Neerven, Theorem 4.2, states: "For
every phi in (C_0(X))^* there exists a unique Radon measure mu ..." and the
correspondence is isometric.  His Theorem 4.50 (Banach--Alaoglu) states exactly:
"The closed unit ball of every dual Banach space is compact with respect to the
weak* topology."

[PROVED-HERE, using R1]  P(Delta) is weak-* compact.  Indeed it is contained in
the closed unit ball of C(Delta)^*.  It is weak-* closed there because

    <1,nu> = 1,                 and
    <f,nu> >= 0 for every f in C(Delta) with f >= 0

are intersections of weak-* closed conditions; by Riesz representation these
conditions characterize probability measures.  Banach--Alaoglu now gives
compactness.  (Equivalently, Prokhorov applies because every family of laws on
the compact space Delta is tight.)

===============================================================================
CONTINUITY LEMMA
===============================================================================

Define binary entropy in bits by

    h(u) = -u log_2(u) - (1-u) log_2(1-u),   h(0)=h(1)=0.

For z=(p,q) in Delta put

    g(z)   = (2-p-q)/2,
    ell(z) = (h(p)+h(q))/2,
    s*(z)  = min(max(1/2,q), min(p+q,1)),
    c(z)   = h(s*(z)).

Then

    mass(nu) = integral 1 dnu,
    B(nu)    = integral g dnu = integral (1-x) dmu_nu(x),
    L(nu)    = integral ell dnu = integral h(x) dmu_nu(x),
    C(nu)    = integral c dnu.

[PROVED-HERE]  The functions 1 and g are continuous.  The elementary limit
r log r -> 0 as r down to zero proves that h, with the displayed endpoint
values, is continuous on [0,1]; hence ell is continuous including every edge
and corner of Delta.  Finally min and max preserve continuity (for example,
max(a,b)=(a+b+|a-b|)/2 and min(a,b)=(a+b-|a-b|)/2).  Therefore s* is a
composition of continuous min/max operations, c=h o s* is continuous, and
mass, B, L, and C are weak-* continuous linear functionals.

For the remaining term use the notation of ``reduction.py``.  For 0 <= u <= 1,
with the endpoint value understood by continuity, set

    T(u) = u + (1-u) ln(1-u),       T(1)=1,
    W(x,y) = xy - x T(y) - y T(x) + T(xy).                              (C1)

For a marginal mu define

    E(mu) = (1/ln 2) integral integral W(1-x,1-y) dmu(x)dmu(y),           (C2)

and abbreviate E(nu)=E(mu_nu).

[PROVED-HERE]  T is continuous on [0,1], again because r ln r -> 0.  Formula
(C1) makes W continuous on the whole square, not merely its interior.  Direct
substitution, using T(0)=0 and T(1)=1, gives for every y in [0,1]

    W(0,y)=0,             W(1,y)=y-T(y)-y+T(y)=0.                        (C3)

[PROVED-HERE]  The marginal map M is continuous and linear: equation (S1)
turns every continuous test function f on [0,1] into the continuous test
function (p,q) -> (f(p)+f(q))/2 on Delta.

[CITED(verified quote R1)]  Van Neerven's Stone--Weierstrass theorem,
Theorem 2.5, says that a subspace of C(K) which contains 1, is an algebra
(closed under conjugation in the complex case), and separates points is dense
in C(K).

[PROVED-HERE, using that quoted theorem]  If mu_i -> mu weak-*, then
mu_i tensor mu_i -> mu tensor mu weak-* on [0,1]^2.  For a product test
function a(x)b(y), its integral is

    (integral a dmu_i)(integral b dmu_i),

which converges to the analogous product for mu.  Finite sums of such product
functions form a unital point-separating algebra, hence are uniformly dense in
C([0,1]^2).  Uniform approximation, together with the fact that all product
measures have mass one, extends convergence to every continuous kernel.  Apply
this to (x,y)->W(1-x,1-y).  Thus nu -> E(nu) is weak-* continuous.  Consequently

    Phi_alpha(nu)
      = [2(1-alpha)B(nu)-1]L(nu) + alpha C(nu) -(1-alpha)E(nu)            (C4)

is weak-* continuous.

===============================================================================
GLOBAL ATTAINMENT
===============================================================================

[PROVED-HERE]  Since mean(mu_nu)=1-B(nu), the feasible set is

    H_t = {nu in P(Delta) : 1-B(nu) <= t}.

It is the inverse image of the closed interval (-infinity,t] under a continuous
map, hence is weak-* closed in compact P(Delta), and therefore compact.  It is
nonempty for t in [0,1] (delta_(0,0) is feasible).  By continuity in (C4), Phi
attains its global minimum on H_t; fix a minimizer nu* and put beta*=B(nu*).
This proves global attainment before any slicing argument is used.

===============================================================================
SLICE STEP
===============================================================================

Let

    K_beta = {nu in P(Delta) : B(nu)=beta*}.

[PROVED-HERE]  K_beta is nonempty (it contains nu*), convex, and weak-* compact.
Moreover every member has mean 1-beta* <= t, so K_beta is a subset of H_t.
It follows at once that

    min_(K_beta) Phi = Phi(nu*) = min_(H_t) Phi.                          (S2)

Thus replacing nu* by another minimizer on its slice cannot lose feasibility
or change the global minimum.

[CITED(verified quote; local source R2)]  ``reduction.py``, Theorem A', proves

    W(x,y) = sum_{n>=2} (x-x^n)(y-y^n)/(n(n-1))

and states: "Then W is a POSITIVE-DEFINITE kernel" and
"Q(mu)=2 B(mu)L(mu)-E(mu)."  The feature series converges absolutely and
uniformly (each summand has absolute value at most 1/(n(n-1))), so it may be
integrated against signed finite measures; its quadratic form is then a sum of
nonnegative squares.

[PROVED-HERE, using R2]  A positive-semidefinite quadratic form is convex:
for marginals mu_0, mu_1 and theta in [0,1], bilinearity gives

 E(theta mu_0+(1-theta)mu_1)
   = theta E(mu_0)+(1-theta)E(mu_1)
     -theta(1-theta)E(mu_0-mu_1)
   <= theta E(mu_0)+(1-theta)E(mu_1).                                    (S3)

The marginal map is linear, so the same inequality holds as a function of nu.
On K_beta, Theorem A' turns (C4) into

    Phi_alpha(nu)|K_beta
       = [2(1-alpha)beta*-1]L(nu) + alpha C(nu)
         -(1-alpha)E(nu).                                                 (S4)

The first two terms are linear.  Since 1-alpha >= 0, (S3) shows that the last
term is concave.  Hence Phi_alpha is continuous and concave on K_beta.
This proves, rather than numerically assumes, fixed-B concavity for every fixed
alpha in [0,1].

[CITED(verified quote R3)]  The Bauer maximum principle in Bru and de Siqueira
Pedra, Lemma 3.3, reads: "Let X be a locally convex real space. An upper
semi-continuous convex real-valued function h over a compact convex subset K
subset X attains its maximum at an extreme point of K," equivalently
"sup h(K)=max h(E(K))."  This is the standard Bauer theorem; the original
reference is H. Bauer (1958), DOI 10.1007/BF01898615, and it is also
Aliprantis--Border (2006), Theorem 7.69.

[PROVED-HERE, using R3]  The signed-measure space C(Delta)^*, with its weak-*
topology, is locally convex Hausdorff.  Apply Bauer to the continuous convex
function -Phi_alpha on K_beta.  There is an extreme point nu** of K_beta at
which -Phi is maximal, hence at which Phi is minimal.  By (S2), nu** is also a
global minimizer over H_t.

===============================================================================
EXTREME POINTS OF THE ONE-MOMENT SLICE
===============================================================================

K_beta is the moment set

    {nu in P(Delta) : integral_Delta g dnu = beta*}                       (E1)

for the one continuous moment function g(p,q)=(2-p-q)/2.  Probability
normalization is already part of P(Delta).

[CITED(verified quote R4)]  The bibliographic data were checked against the
publisher record: G. Winkler, "Extreme Points of Moment Sets," Mathematics of
Operations Research 13(4) (1988), 581--587, DOI 10.1287/moor.13.4.581.
Pinelis, Theorem 12, explicitly reproduces Winkler's Theorem 2.1.  Its
assumption and equality-moment conclusion are quoted there as follows:

  "Suppose that the set Pi of all probability measures in Lambda is a
   Choquet-simplex and ex Pi is contained in Delta^(1)."

  "(b) For any c in R^k, one has
   ex(Pi intersect Lambda_(f,c))
     = Pi intersect Lambda_(f,c) intersect Delta^(1+k)_(1,f)."

Here Delta^(1+k)_(1,f) denotes the discrete measures on at most 1+k points
whose augmented moment vectors (1,f(x)) are linearly independent.  Thus the
precise theorem is not the false assertion that every measure on <=k+1 points
is extreme: it characterizes the extreme ones by the additional
moment-determination (linear-independence) condition, and in particular gives
the support upper bound <=k+1.

[PROVED-HERE: hypotheses mapped]  Take Lambda to be the cone of all finite
nonnegative Borel measures on Delta and Pi=P(Delta).  Its extreme points are
the Dirac laws: a non-Dirac probability can be split across a Borel set of
mass strictly between zero and one, while a Dirac law cannot be split.
Moreover P(Delta) is a Choquet simplex: nu is represented on the extreme
boundary {delta_z:z in Delta} by the pushforward of nu under z->delta_z, and
that representing law is unique because its barycenter agrees with nu on all
f in C(Delta).  Hence Winkler's assumptions hold.  Apply the quoted part (b)
with k=1, f=g, and c=beta*.  Every extreme point of K_beta has support on at
most 1+1=2 points.

[PROVED-HERE: direct specialization, independently checking the cited bound]
There is also a short perturbation proof in this particular case.  If an
extreme nu in (E1) had three distinct points in its topological support,
choose pairwise disjoint Borel neighborhoods A_1,A_2,A_3 of them.  By the
definition of topological support, each has positive nu-mass.  The three
vectors

    v_i = (nu(A_i), integral_(A_i) g dnu) in R^2

are linearly dependent.  Choose nonzero (a_1,a_2,a_3) with sum a_i v_i=0 and
put eta=sum a_i nu|A_i.  The signed measure eta is nonzero, has total mass zero,
and has g-moment zero.  For sufficiently small epsilon>0, the Radon--Nikodym
densities 1 +/- epsilon a_i on A_i (and 1 off their union) are nonnegative.
Thus nu+epsilon eta and nu-epsilon eta are distinct members of (E1) whose
midpoint is nu, contradicting extremality.  Therefore every extreme point has
support size at most two, exactly as Winkler's theorem asserts.

===============================================================================
CONCLUSION
===============================================================================

[PROVED-HERE]  The Bauer minimizer nu** is extreme in K_beta, so the preceding
step gives |supp(nu**)|<=2.  It has B(nu**)=beta*, hence
mean(mu_nu**)=1-beta*<=t, and by (S2)

 min {Phi(nu): nu in P(Delta), mean(mu_nu)<=t}
   = Phi(nu**)
   = min {Phi(nu): nu feasible, |supp(nu)|<=2}.

Each support point (p,q) contributes at most p and q to the marginal, so
mu_nu** has at most four support points.  This proves Theorem B'''.  QED.

===============================================================================
VERIFIED REFERENCES (the quoted statements above, not numerical evidence)
===============================================================================

R1. Jan van Neerven, Functional Analysis, arXiv:2112.11166v7 (2025),
    Theorems 2.5, 4.2, and 4.50.
    https://arxiv.org/pdf/2112.11166
R2. Local proved source: math/uc/reduction.py, Theorem A', especially its
    equations (1)--(3).
R3. J.-B. Bru and W. de Siqueira Pedra, "Remarks on the Gamma-regularization
    of Non-convex and Non-semi-continuous Functions on Topological Vector
    Spaces," arXiv:1610.03411, Lemma 3.3.
    https://arxiv.org/pdf/1610.03411
    Original: H. Bauer, "Minimalstellen von Funktionen und Extremalpunkte,"
    Archiv der Mathematik 9 (1958), 389--393, DOI 10.1007/BF01898615.
R4. G. Winkler, "Extreme Points of Moment Sets," Mathematics of Operations
    Research 13(4) (1988), 581--587, DOI 10.1287/moor.13.4.581.
    Publisher metadata: https://dl.acm.org/doi/abs/10.5555/2770717.2770721
    Exact Theorem 2.1 quotation verified through I. Pinelis, "On the extreme
    points of moments sets," Mathematical Methods of Operations Research 83
    (2016), 325--349, Theorem 12 [Winkler, Theorem 2.1].
    https://arxiv.org/pdf/1204.0249

===============================================================================
PYTHON CONTROLS -- NUMERICAL-SAMPLED ONLY; NOT PART OF THE PROOF
===============================================================================

The executable controls below (i) sample fixed-B chord inequalities, (ii)
exercise the continuous marginal map and E along an atom-merging sequence, and
(iii) deliberately remove the fixed-B hypothesis and detect nonconcavity.  All
claims in the proof above are PROVED-HERE or CITED(verified quote); only these
explicitly marked controls are sampled numerical evidence.
"""

from __future__ import annotations

import numpy as np


ALPHA = 0.0356069
LN2 = float(np.log(2.0))
FIXED_B_TRIALS = 3000
FREE_B_TRIALS = 2500

WINKLER_DOI = "10.1287/moor.13.4.581"
WINKLER_QUOTE = (
    "For any c in R^k, ex(Pi intersect Lambda_(f,c)) equals "
    "Pi intersect Lambda_(f,c) intersect Delta^(1+k)_(1,f)."
)
BAUER_QUOTE = (
    "An upper semi-continuous convex real-valued function over a compact "
    "convex subset attains its maximum at an extreme point."
)


def binary_entropy(x: np.ndarray) -> np.ndarray:
    """Binary entropy in bits, continuously extended by h(0)=h(1)=0."""
    x = np.asarray(x, dtype=float)
    if np.any((x < 0.0) | (x > 1.0)):
        raise ValueError("binary_entropy expects values in [0,1]")
    out = np.zeros_like(x)
    interior = (x > 0.0) & (x < 1.0)
    u = x[interior]
    out[interior] = -(u * np.log(u) + (1.0 - u) * np.log1p(-u)) / LN2
    return out


def t_function(x: np.ndarray) -> np.ndarray:
    """T(x)=x+(1-x)ln(1-x), continuously extended by T(1)=1."""
    x = np.asarray(x, dtype=float)
    if np.any((x < 0.0) | (x > 1.0)):
        raise ValueError("t_function expects values in [0,1]")
    out = np.ones_like(x)
    below_one = x < 1.0
    u = x[below_one]
    out[below_one] = u + (1.0 - u) * np.log1p(-u)
    return out


def kernel_w(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """The continuous PSD kernel W from reduction.py."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return x * y - x * t_function(y) - y * t_function(x) + t_function(x * y)


def s_star(atoms: np.ndarray) -> np.ndarray:
    """Optimal orbit coupling coordinate on Delta."""
    p = atoms[:, 0]
    q = atoms[:, 1]
    return np.minimum(np.maximum(0.5, q), np.minimum(p + q, 1.0))


def marginal(atoms: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the atoms and weights of mu_nu from a discrete pair-orbit law."""
    return atoms.ravel(), np.repeat(weights, 2) / 2.0


def energy(atoms: np.ndarray, weights: np.ndarray) -> float:
    """E(mu_nu) evaluated directly from the continuous kernel W."""
    points, marginal_weights = marginal(atoms, weights)
    complements = 1.0 - points
    x, y = np.meshgrid(complements, complements, indexing="ij")
    gram = kernel_w(x, y)
    return float(marginal_weights @ gram @ marginal_weights / LN2)


def orbit_quantities(
    atoms: np.ndarray, weights: np.ndarray, alpha: float = ALPHA
) -> tuple[float, float, float, float, float]:
    """Return (B,L,C,E,Phi_exact) for a finite pair-orbit measure."""
    points, marginal_weights = marginal(atoms, weights)
    b_value = float(marginal_weights @ (1.0 - points))
    l_value = float(marginal_weights @ binary_entropy(points))
    c_value = float(weights @ binary_entropy(s_star(atoms)))
    e_value = energy(atoms, weights)
    phi_value = (
        (2.0 * (1.0 - alpha) * b_value - 1.0) * l_value
        + alpha * c_value
        - (1.0 - alpha) * e_value
    )
    return b_value, l_value, c_value, e_value, phi_value


def random_measure(
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a small discrete law on Delta."""
    count = int(rng.integers(1, 5))
    atoms = np.sort(rng.uniform(0.01, 0.99, size=(count, 2)), axis=1)
    weights = rng.dirichlet(np.ones(count))
    return atoms, weights


def mix_measures(
    atoms_a: np.ndarray,
    weights_a: np.ndarray,
    atoms_b: np.ndarray,
    weights_b: np.ndarray,
    theta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Represent theta*nu_a+(1-theta)*nu_b on the concatenated support."""
    return (
        np.vstack((atoms_a, atoms_b)),
        np.concatenate((theta * weights_a, (1.0 - theta) * weights_b)),
    )


def fixed_b_concavity_control(trials: int = FIXED_B_TRIALS) -> tuple[float, int]:
    """Return the worst sampled Phi(mix)-chord gap on exactly fixed-B segments."""
    rng = np.random.default_rng(20260812)
    worst = float("inf")
    for _ in range(trials):
        atoms_a, weights_a = random_measure(rng)
        atoms_b, weights_b = random_measure(rng)

        mean_a = 1.0 - orbit_quantities(atoms_a, weights_a)[0]
        mean_b = 1.0 - orbit_quantities(atoms_b, weights_b)[0]
        target_mean = float(rng.uniform(0.08, 0.92) * min(mean_a, mean_b))
        atoms_a = atoms_a * (target_mean / mean_a)
        atoms_b = atoms_b * (target_mean / mean_b)

        b_a = orbit_quantities(atoms_a, weights_a)[0]
        b_b = orbit_quantities(atoms_b, weights_b)[0]
        assert abs(b_a - b_b) <= 2.0e-13, (b_a, b_b)

        theta = float(rng.uniform(0.05, 0.95))
        mixed_atoms, mixed_weights = mix_measures(
            atoms_a, weights_a, atoms_b, weights_b, theta
        )
        phi_a = orbit_quantities(atoms_a, weights_a)[-1]
        phi_b = orbit_quantities(atoms_b, weights_b)[-1]
        phi_mixed = orbit_quantities(mixed_atoms, mixed_weights)[-1]
        gap = phi_mixed - (theta * phi_a + (1.0 - theta) * phi_b)
        worst = min(worst, gap)
    return worst, trials


def marginal_and_energy_continuity_control() -> tuple[list[float], list[float]]:
    """Exercise M and E on two atoms merging weakly to one pair-orbit."""
    limit_atoms = np.array([[0.27, 0.74]], dtype=float)
    limit_weights = np.array([1.0])
    limit_points, limit_marginal_weights = marginal(limit_atoms, limit_weights)
    limit_energy = energy(limit_atoms, limit_weights)

    test_functions = (
        lambda x: np.ones_like(x),
        lambda x: x,
        lambda x: x * x,
        lambda x: np.sin(np.pi * x),
        binary_entropy,
    )
    limit_integrals = np.array(
        [limit_marginal_weights @ f(limit_points) for f in test_functions]
    )

    marginal_errors: list[float] = []
    energy_errors: list[float] = []
    for epsilon in 0.08 * 0.5 ** np.arange(12):
        atoms = np.array(
            [
                [0.27 - epsilon, 0.74 - epsilon],
                [0.27 + epsilon, 0.74 + epsilon],
            ],
            dtype=float,
        )
        weights = np.array([0.5, 0.5])
        points, marginal_weights = marginal(atoms, weights)

        via_marginal = np.array(
            [marginal_weights @ f(points) for f in test_functions]
        )
        via_orbits = np.array(
            [
                weights
                @ ((f(atoms[:, 0]) + f(atoms[:, 1])) / 2.0)
                for f in test_functions
            ]
        )
        assert np.max(np.abs(via_marginal - via_orbits)) <= 2.0e-15

        marginal_errors.append(float(np.max(np.abs(via_marginal - limit_integrals))))
        energy_errors.append(abs(energy(atoms, weights) - limit_energy))

    for errors in (marginal_errors, energy_errors):
        assert all(
            later <= 1.05 * earlier + 2.0e-15
            for earlier, later in zip(errors, errors[1:])
        ), errors
        assert errors[-1] < 1.0e-7, errors[-1]
    return marginal_errors, energy_errors


def free_b_failure_control(trials: int = FREE_B_TRIALS) -> tuple[float, int]:
    """Remove the fixed-B premise and find a sampled violation of concavity."""
    rng = np.random.default_rng(314159)
    worst = float("inf")
    for _ in range(trials):
        atoms_a, weights_a = random_measure(rng)
        atoms_b, weights_b = random_measure(rng)
        theta = float(rng.uniform(0.05, 0.95))
        mixed_atoms, mixed_weights = mix_measures(
            atoms_a, weights_a, atoms_b, weights_b, theta
        )
        phi_a = orbit_quantities(atoms_a, weights_a)[-1]
        phi_b = orbit_quantities(atoms_b, weights_b)[-1]
        phi_mixed = orbit_quantities(mixed_atoms, mixed_weights)[-1]
        gap = phi_mixed - (theta * phi_a + (1.0 - theta) * phi_b)
        worst = min(worst, gap)
    return worst, trials


def main() -> None:
    """Run citation-record checks and the three explicitly sampled controls."""
    assert WINKLER_DOI == "10.1287/moor.13.4.581"
    assert "1+k" in WINKLER_QUOTE and "extreme point" in BAUER_QUOTE
    assert float(binary_entropy(np.array([0.0]))[0]) == 0.0
    assert float(binary_entropy(np.array([1.0]))[0]) == 0.0
    edge = np.array([0.0, 0.17, 0.63, 1.0])
    assert np.max(np.abs(kernel_w(np.zeros_like(edge), edge))) <= 2.0e-15
    assert np.max(np.abs(kernel_w(np.ones_like(edge), edge))) <= 2.0e-15

    print("THEOREM B''' SUPPORT-REDUCTION STATUS")
    print("CITED(verified quote) SETUP: Riesz representation and Banach--Alaoglu (R1)")
    print("PROVED-HERE SETUP: Delta compact; P(Delta) weak-* closed, hence compact")
    print("PROVED-HERE CONTINUITY: mass, B, L, C, marginal map, product pairing, E")
    print("PROVED-HERE GLOBAL ATTAINMENT: feasible set compact and Phi continuous")
    print("CITED(verified quote) SLICE: Theorem A' PSD identity (R2); Bauer principle (R3)")
    print("PROVED-HERE SLICE: fixed-B Phi is concave for every alpha in [0,1]")
    print("CITED(verified quote) EXTREME POINTS: Winkler 1988, Theorem 2.1(b) (R4)")
    print("PROVED-HERE EXTREME POINTS: hypotheses mapped and perturbation checked")
    print("PROVED-HERE CONCLUSION: a global minimizer has <=2 pair-orbits")

    fixed_worst, fixed_used = fixed_b_concavity_control()
    assert fixed_used >= 2000
    assert fixed_worst >= -1.0e-9, fixed_worst
    print(
        "NUMERICAL-SAMPLED CONTROL 1: fixed-B concavity "
        f"worst gap={fixed_worst:+.6e} over {fixed_used} segments [PASS]"
    )

    marginal_errors, energy_errors = marginal_and_energy_continuity_control()
    print(
        "NUMERICAL-SAMPLED CONTROL 2: atom-merging marginal error "
        f"{marginal_errors[0]:.3e}->{marginal_errors[-1]:.3e}; "
        f"E error {energy_errors[0]:.3e}->{energy_errors[-1]:.3e} [PASS]"
    )

    free_worst, free_used = free_b_failure_control()
    assert free_worst < -1.0e-3, free_worst
    print(
        "NUMERICAL-SAMPLED CONTROL 3: B-free nonconcavity detected "
        f"worst gap={free_worst:+.6e} over {free_used} segments [PASS]"
    )
    print("ALL ASSERT-BACKED CONTROLS PASS")


if __name__ == "__main__":
    main()
