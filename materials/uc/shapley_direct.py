"""DIRECT permutation-average frontier for the Sawin/Cambie functional.

For coordinate i and predecessor set S, let mu_{i,S} be the exact law of
Pr(X_i=1 | X_S) for X uniform on a finite family.  In a uniform random
permutation, S is the predecessor set of i with probability

    |S|! (n-|S|-1)! / n! = 1 / (n * binom(n-1, |S|)).

Consequently the random-order quantity is exactly

    A(X) = sum_i sum_{S subset [n]\\{i}}
           F(mu_{i,S}) / (n * binom(n-1, |S|)).

This is the correct direct average E_pi sum_i F(mu_{i,pi}); it is not the
nonlinear quantity F(E_pi mu_{i,pi}).  The identity is a counting theorem and
requires no convexity of F.  It also reduces n! orders to n*2^(n-1) local laws.

Equivalently, reveal every other coordinate independently with probability
t and call the revealed set R_t.  The beta integral

    integral_0^1 t^|S| (1-t)^(n-1-|S|) dt
      = 1 / (n * binom(n-1, |S|))

gives the exact random-reveal representation

    A(X) = sum_i integral_0^1 E F(mu_{i,R_t}) dt.

The same identity with F replaced by conditional entropy sums to H(X).  This
turns the conjecture into an integral inequality along a posterior-martingale
reveal process; it still does not justify a pointwise inequality in t.

Writing D(t)=sum_i E H(X_i|X_{R_t}), differentiation of the multilinear
extension gives another exact identity:

    D'(t) = -sum_{i != j} E I(X_i;X_j | X_{R_t^{i,j}}) <= 0,

where R_t^{i,j} independently reveals coordinates outside {i,j}.  Thus Reimer
constrains the area integral of a monotone entropy profile whose derivative is
an explicit conditional-mutual-information flow.

The same representation proves exact tensorization.  If X=(Y,Z) is an
independent product, conditioning a Y-coordinate on predecessors from Z changes
nothing, and the relative order induced on the Y-coordinates is uniform (and
likewise for Z).  Hence A(Y,Z)=A(Y)+A(Z), while H(Y,Z)=H(Y)+H(Z).

The surviving conjecture at the current Reimer target scale is

    max_i Pr(X_i=1) <= 0.38261 and 2 E|X| >= H(X)
        ==> A(X) >= 0.

The complete n=4 calculation suggested the stronger normalized bound
A(X)>=gamma_B H(X).  `shapley_n5_complete.py` finds 25 robust n=5 violations
in two exact symmetry orbits, while its complete float64 comparison finds no
negative A among 64,534 nontrivial admissible families (the trivial family has
A=0).  `shapley_shared_signal.py` further proves A<0 if uniformity on the
support is dropped.  Thus uniform finite-family realizability is essential.
The seeded search below is retained as a coverage control: it missed the rare
n=5 normalized counterexamples and cannot replace enumeration.

Evidence labels:
  * PROVED: predecessor-set/random-reveal weighting, independent-product
    tensorization, and the base-family gamma enclosure in
    `shapley_direct_cert.py`.
  * COMPLETE NUMERICAL: every admissible family on [4], with exact Fraction
    conditional laws and float64 HiGHS transport costs.
  * COMPLETE n=5 ENUMERATION / NUMERICAL F: `shapley_n5_complete.py` refutes
    the gamma_B strengthening but finds no negative A among nontrivial records.
  * CERTIFIED BOUNDARY: `shapley_shared_signal.py` refutes A>=0 for general
    nonuniform binary distributions, not for uniform finite families.
  * NUMERICAL / INCOMPLETE: deterministic seeded samples on [5] through [8].

Run: ./.venv/bin/python uc/shapley_direct.py
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import ceil, comb, factorial, floor, log2
import random
import time

from shapley_entropy import functional, law_key

TARGET = 0.38261
BASE4 = (0, 1, 2, 4, 7, 9, 10, 12)


def predecessor_law(mem, coordinate, predecessors):
    """Exact law of Pr(X_coordinate=1 | X_predecessors)."""
    predecessors = tuple(sorted(predecessors))
    groups = {}
    for subset in mem:
        pattern = tuple((subset >> j) & 1 for j in predecessors)
        total, ones = groups.get(pattern, (0, 0))
        groups[pattern] = (total + 1, ones + ((subset >> coordinate) & 1))
    mass = {}
    for total, ones in groups.values():
        atom = Fraction(ones, total)
        mass[atom] = mass.get(atom, Fraction()) + Fraction(total, len(mem))
    assert sum(mass.values(), Fraction()) == 1
    return tuple(sorted(mass.items()))


def predecessor_weight(n, size):
    return Fraction(1, n * comb(n - 1, size))

def reveal_integral_weight(n, size):
    """Exact beta-integral weight for an independently revealed subset."""
    return Fraction(factorial(size) * factorial(n - 1 - size), factorial(n))


def direct_entropy_layers(mem, n):
    """Shapley conditional-entropy contributions by predecessor cardinality."""
    layers = [0.0] * n
    for coordinate in range(n):
        others = tuple(j for j in range(n) if j != coordinate)
        for size in range(n):
            weight = float(predecessor_weight(n, size))
            for predecessors in combinations(others, size):
                entropy = functional(predecessor_law(
                    mem, coordinate, predecessors))[2]
                layers[size] += weight * entropy
    return tuple(layers)


def reveal_profile(layers, t):
    """Evaluate the Bernoulli-reveal Bernstein polynomial at t."""
    n = len(layers)
    return n * sum(
        comb(n - 1, size) * layers[size]
        * t ** size * (1.0 - t) ** (n - 1 - size)
        for size in range(n)
    )

def reveal_entropy_layer_derivative(layers, t):
    """Differentiate the Bernstein representation of D(t)."""
    n = len(layers)
    if n == 1:
        return 0.0
    return n * (n - 1) * sum(
        comb(n - 2, size) * (layers[size + 1] - layers[size])
        * t ** size * (1.0 - t) ** (n - 2 - size)
        for size in range(n - 1)
    )


def reveal_entropy_edge_derivative(mem, n, t):
    """Compute D'(t) as the negative conditional-mutual-information flow."""
    entropy_cache = {}

    def conditional_entropy(coordinate, predecessors):
        key = (coordinate, tuple(sorted(predecessors)))
        if key not in entropy_cache:
            entropy_cache[key] = functional(predecessor_law(
                mem, coordinate, key[1]))[2]
        return entropy_cache[key]

    derivative = 0.0
    for coordinate in range(n):
        for revealed_next in range(n):
            if revealed_next == coordinate:
                continue
            others = tuple(
                j for j in range(n)
                if j != coordinate and j != revealed_next
            )
            for size in range(n - 1):
                probability = t ** size * (1.0 - t) ** (n - 2 - size)
                for predecessors in combinations(others, size):
                    before = conditional_entropy(coordinate, predecessors)
                    after = conditional_entropy(
                        coordinate, predecessors + (revealed_next,))
                    derivative += probability * (after - before)
    return derivative


def direct_average(mem, n):
    """Return A(X) and its n predecessor-cardinality layer contributions."""
    total = 0.0
    layers = [0.0] * n
    for coordinate in range(n):
        others = tuple(j for j in range(n) if j != coordinate)
        for size in range(n):
            weight = float(predecessor_weight(n, size))
            for predecessors in combinations(others, size):
                value = functional(predecessor_law(mem, coordinate,
                                                   predecessors))[0]
                contribution = weight * value
                total += contribution
                layers[size] += contribution
    return total, tuple(layers)


def permutation_average(mem, n):
    total = 0.0
    count = 0
    for order in permutations(range(n)):
        count += 1
        total += sum(functional(law_key(mem, order, position))[0]
                     for position in range(n))
    return total / count


def frequencies(mem, n):
    return tuple(sum((subset >> i) & 1 for subset in mem) / len(mem)
                 for i in range(n))


def reimer_slack(mem):
    average_size = sum(subset.bit_count() for subset in mem) / len(mem)
    return 2.0 * average_size - log2(len(mem))


def admissible(mem, n):
    return (max(frequencies(mem, n)) <= TARGET + 1e-15
            and reimer_slack(mem) >= -1e-12)


def feasible_sizes(n):
    sizes = []
    for family_size in range(2, (1 << n) + 1):
        incidence_cap = n * floor(TARGET * family_size + 1e-15)
        incidence_need = ceil(family_size * log2(family_size) / 2.0 - 1e-12)
        if incidence_cap >= incidence_need:
            sizes.append(family_size)
    return tuple(sizes)


def sample_family(rng, n, family_size, bit_probability):
    members = set()
    attempts = 0
    while len(members) < family_size and attempts < 20 * family_size:
        subset = 0
        for coordinate in range(n):
            if rng.random() < bit_probability:
                subset |= 1 << coordinate
        members.add(subset)
        attempts += 1
    if len(members) != family_size:
        return None
    return tuple(sorted(members))


def falsification_sample(n, wanted, max_draws, seed):
    rng = random.Random(seed)
    sizes = feasible_sizes(n)
    families = {BASE4}  # Embed the known four-coordinate extremizer with dummies.
    draws = 0
    while len(families) < wanted + 1 and draws < max_draws:
        family_size = rng.choice(sizes)
        bit_probability = rng.uniform(0.20, 0.44)
        mem = sample_family(rng, n, family_size, bit_probability)
        draws += 1
        if mem is not None and admissible(mem, n):
            families.add(mem)
    return families, draws


def product_family(left, right, shift):
    return tuple(sorted(a | (b << shift) for a in left for b in right))


def exhaustive_n4():
    records = []
    for family_mask in range(1, 1 << 16):
        mem = tuple(subset for subset in range(16)
                    if (family_mask >> subset) & 1)
        if len(mem) > 1 and admissible(mem, 4):
            value, layers = direct_average(mem, 4)
            records.append((value / log2(len(mem)), value, mem, layers))
    return min(records), len(records)


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    # Exact combinatorial controls.
    for n in range(1, 13):
        for size in range(n):
            assert predecessor_weight(n, size) == reveal_integral_weight(n, size)
        assert sum(comb(n - 1, size) * predecessor_weight(n, size)
                   for size in range(n)) == 1
    base_value, base_layers = direct_average(BASE4, 4)
    explicit_value = permutation_average(BASE4, 4)
    assert abs(base_value - explicit_value) < 2e-12
    gamma = base_value / log2(len(BASE4))
    base_entropy_layers = direct_entropy_layers(BASE4, 4)
    assert abs(sum(base_entropy_layers) - log2(len(BASE4))) < 2e-12
    for t in (0.0, 0.25, 0.5, 0.75, 1.0):
        layer_derivative = reveal_entropy_layer_derivative(
            base_entropy_layers, t)
        edge_derivative = reveal_entropy_edge_derivative(BASE4, 4, t)
        assert abs(layer_derivative - edge_derivative) < 2e-12
        assert edge_derivative <= 2e-12

    product8 = product_family(BASE4, BASE4, 4)
    product_value, _ = direct_average(product8, 8)
    assert admissible(product8, 8)
    assert abs(product_value - 2.0 * base_value) < 5e-12

    print("PROVED [counting + reveal-flow + independent-product identities]")
    print("  beta-integral weights equal predecessor weights for n=1..12")
    print("  D'(t) equals negative conditional-MI flow (direct n=4 control)")
    print("  direct n*2^(n-1) formula agrees with all 4! orders on the base family")
    print("  base A(X) = %.15f" % base_value)
    print("  base gamma=A/H = %.15f" % gamma)
    print("  two-block product A(X) = %.15f = 2*base (float tolerance)" %
          product_value)
    print()
    differences = tuple(
        value - gamma * entropy
        for value, entropy in zip(base_layers, base_entropy_layers)
    )
    print("RANDOM-REVEAL PROFILE OF BASE4 — NUMERICAL F VALUES")
    print("  F predecessor layers:", base_layers)
    print("  H predecessor layers:", base_entropy_layers)
    print("  layer coefficients F-gamma*H:", differences)
    for t in (0.0, 0.5, 1.0):
        value = reveal_profile(base_layers, t)
        entropy = reveal_profile(base_entropy_layers, t)
        print("  t=%.1f: G(t)-gamma*D(t) = %+.15f" %
              (t, value - gamma * entropy))
    print("  D'(1/2) from conditional-MI flow = %.15f" %
          reveal_entropy_edge_derivative(BASE4, 4, 0.5))
    assert reveal_profile(differences, 0.0) < 0
    assert reveal_profile(differences, 0.5) < 0
    assert reveal_profile(differences, 1.0) > 0
    print("  Pointwise-in-t and t<->1-t midpoint bounds are false even for BASE4.")
    print()

    t0 = time.time()
    worst4, count4 = exhaustive_n4()
    print("COMPLETE NUMERICAL n=4")
    print("  nontrivial Reimer-admissible families at max frequency <= %.5f: %d" %
          (TARGET, count4))
    print("  minimum A/H = %.15f" % worst4[0])
    print("  A = %.15f; family = %s" % (worst4[1], worst4[2]))
    print("  predecessor-size layers =", worst4[3])
    print("  matches base gamma within 1e-12:", abs(worst4[0] - gamma) < 1e-12)
    print()

    sample_specs = (
        (5, 1000, 1_000_000, 2026082307),
        (6, 500, 300_000, 2026082308),
        (7, 250, 200_000, 2026082309),
        (8, 100, 100_000, 2026082310),
    )
    sampled_below_gamma = []
    sampled_negative = []
    print("SEEDED COVERAGE CONTROL — NUMERICAL / INCOMPLETE")
    for n, wanted, max_draws, seed in sample_specs:
        families, draws = falsification_sample(n, wanted, max_draws, seed)
        worst = None
        for mem in families:
            value, layers = direct_average(mem, n)
            ratio = value / log2(len(mem))
            record = (ratio, value, mem, layers)
            if worst is None or record[0] < worst[0]:
                worst = record
            if ratio < gamma - 1e-10:
                sampled_below_gamma.append((n, record))
            if value < -1e-10:
                sampled_negative.append((n, record))
        print("  n=%d: %d families after %d draws; min A/H=%.15f; family=%s" %
              (n, len(families), draws, worst[0], worst[2]))
    print("  sampled ratios below the refuted gamma_B benchmark:",
          len(sampled_below_gamma))
    print("  (complete n=5 enumeration finds 25; this sample misses all 25)")
    print("  sampled negative A values:", len(sampled_negative))
    print("  runtime: %.2f s" % (time.time() - t0))
    print()
    print("FINAL FRONTIER")
    print("  REFUTED by complete n=5: A(X) >= gamma_B H(X).")
    print("  REFUTED by a certified shared-signal model for nonuniform laws: A>=0.")
    print("  SURVIVING UNIFORM-FAMILY CONJECTURE: A(X) >= 0 under max marginal")
    print("  <= %.5f and Reimer; complete n=4 and n=5 support it." % TARGET)
