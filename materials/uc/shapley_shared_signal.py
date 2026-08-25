"""CERTIFIED nonuniform probability-relaxation failure of direct Shapley A>=0.

This is a boundary result, not a finite-family counterexample.  Let

    c = 1913/5000,
    a = 1913/25000,
    b = 7652/23087,

so c=a+(1-a)b and c<0.38261.  A latent U is Bernoulli(a).  Conditional on U=1,
every X_i equals 1; conditional on U=0, the X_i are iid Bernoulli(b).  Add d
identical coordinate copies of U.  For k=700 and d=1073, all coordinate
marginals obey the cap and even the stronger support-size condition
$2E|X|>=log_2|supp(X)|$ holds, yet the direct random-order functional A is
strictly negative.

The permutation average has an exact O(k) reduction.  Let q be the number of
X-coordinates preceding the first of the d+1 U-copies.  Its negative-
hypergeometric probability is rational.  Before that first U-copy, the local
laws depend only on the number of earlier X's; the first U-copy has one more
exact two-atom law; after it, every remaining X has the same bad law and later
U-copies are deterministic.  Every transport problem therefore has at most two
atoms.  The script checks both rational vertices of each 2x2 transportation
polytope and evaluates all entropy terms with 4096-bit Arb.

SCOPE.  The joint law itself is not uniform and therefore does not refute A>=0
for the original full uniform distribution or a union-closed family.  Because
it is rational, however, `shapley_uniform_lift.py` realizes it as a projection
of an admissible uniform simple family.  It therefore refutes every induction
invariant demanding nonnegative projections.  Only a global argument that
captures the compensating label cost at the full-uniform boundary remains open.

Run: ./.venv/bin/python uc/shapley_shared_signal.py
"""

from fractions import Fraction

from flint import arb, ctx

ctx.prec = 4096
ZERO = arb(0)
ONE = arb(1)
LOG2 = arb(2).log()
ALPHA = arb("0.0356069")

C = Fraction(1913, 5000)
A = Fraction(1913, 25000)
B = Fraction(7652, 23087)
K = 700
D = 1073


def aq(value):
    return arb(value.numerator) / value.denominator


def h(value):
    if value == 0 or value == 1:
        return ZERO
    x = aq(value)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def sstar(p, q):
    return max(p, q, min(p + q, Fraction(1, 2)))


def law(items):
    mass = {}
    for atom, weight in items:
        if weight:
            mass[atom] = mass.get(atom, Fraction()) + weight
    result = tuple(sorted(mass.items()))
    assert sum(weight for _atom, weight in result) == 1
    return result


def min_enclosure(values):
    lower = values[0].lower()
    upper = values[0].upper()
    for value in values[1:]:
        if value.lower() < lower:
            lower = value.lower()
        if value.upper() < upper:
            upper = value.upper()
    return lower.union(upper)


def functional(local_law):
    """Certified F for a one- or two-atom rational law."""
    atoms = tuple(atom for atom, _weight in local_law)
    weights = tuple(weight for _atom, weight in local_law)
    entropy = sum((aq(weight) * h(atom)
                   for atom, weight in local_law), ZERO)
    iid = sum(
        (aq(weights[i] * weights[j])
         * h(atoms[i] + atoms[j] - atoms[i] * atoms[j])
         for i in range(len(atoms)) for j in range(len(atoms))),
        ZERO,
    )

    if len(atoms) == 1:
        coupled = h(sstar(atoms[0], atoms[0]))
    else:
        assert len(atoms) == 2
        w0 = weights[0]
        costs = []
        # Every equal-marginal 2x2 plan is parameterized by pi_00=x.
        for x in (max(Fraction(), 2 * w0 - 1), w0):
            plan = (
                (x, w0 - x),
                (w0 - x, 1 - 2 * w0 + x),
            )
            cost = sum(
                (aq(plan[i][j]) * h(sstar(atoms[i], atoms[j]))
                 for i in range(2) for j in range(2)),
                ZERO,
            )
            costs.append(cost)
        coupled = min_enclosure(costs)
    return (ONE - ALPHA) * iid + ALPHA * coupled - entropy


def certificate(k=K, duplicates=D):
    assert C == A + (1 - A) * B
    bad_law = law(((B, 1 - A), (Fraction(1), A)))
    bad_value = functional(bad_law)

    before_x = []
    first_u = []
    b_power = Fraction(1)
    for q in range(k + 1):
        all_one_mass = A + (1 - A) * b_power
        u_posterior = A / all_one_mass
        first_u.append(functional(law((
            (Fraction(0), 1 - all_one_mass),
            (u_posterior, all_one_mass),
        ))))
        if q < k:
            next_x = (A + (1 - A) * b_power * B) / all_one_mass
            before_x.append(functional(law((
                (B, 1 - all_one_mass),
                (next_x, all_one_mass),
            ))))
        b_power *= B

    n = k + duplicates + 1
    survival = Fraction(1)
    probability_total = Fraction()
    prefix = ZERO
    direct = ZERO
    for q in range(k + 1):
        probability = survival * Fraction(duplicates + 1, n - q)
        probability_total += probability
        direct += aq(probability) * (
            prefix + first_u[q] + (k - q) * bad_value)
        if q < k:
            prefix += before_x[q]
            survival *= Fraction(k - q, n - q)
    assert probability_total == 1

    entropy = h(A) + aq(1 - A) * k * h(B)
    twice_mean = 2 * aq(k * C + (duplicates + 1) * A)
    reimer_slack = twice_mean - entropy
    support_entropy = (arb(2) ** k + ONE).log() / LOG2
    support_slack = twice_mean - support_entropy
    return direct, reimer_slack, support_slack, bad_value, entropy

def compact(value):
    """Short display of an already-certified Arb interval."""
    return "[%.18g, %.18g]" % (float(value.lower()), float(value.upper()))


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    direct, reimer_slack, support_slack, bad_value, entropy = certificate()
    assert C * 100000 < 38261
    assert bad_value.upper() < 0
    assert reimer_slack.lower() > arb("107")
    assert support_slack.lower() > arb("0.0049")
    assert direct.upper() < arb("-0.049")

    print("PROVED [exact rational order weights + 4096-bit Arb]")
    print("  k / duplicate U coordinates / total coordinates =",
          K, "/", D + 1, "/", K + D + 1)
    print("  X marginal c =", C, "= %.18g < 0.38261" % float(C))
    print("  U-copy marginal a =", A, "= %.18g" % float(A))
    print("  bad conditional atom b =", B, "= %.18g" % float(B))
    print("  F(bad law) =", compact(bad_value), "< 0")
    print("  H =", compact(entropy))
    print("  Shannon Reimer slack 2E|X|-H =", compact(reimer_slack), "> 0")
    print("  support Reimer slack 2E|X|-log2|supp| =",
          compact(support_slack), "> 0")
    print("  direct permutation average A =", compact(direct), "< 0")
    print()
    print("FINAL VERDICT")
    print("  REFUTED for general nonuniform binary distributions: direct A>=0.")
    print("  Via uniform_lift: REFUTED as a projection-nonnegativity invariant.")
    print("  NOT REFUTED: A>=0 for the original full uniform or union-closed law.")
