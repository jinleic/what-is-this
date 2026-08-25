"""JOIN-HOMOMORPHIC deterministic dependencies in union-closed families on [4].

`reimer_entropy.py` proves that scalar Reimer slack cannot break the c* ceiling:
zero-entropy deterministic conditional coordinates can carry positive slack at
zero functional cost.  `shapley_gap.py` proves a strict gap when such a
deterministic equality law and the entropy-active c* equality law occur inside
one local law.  The remaining structural question is stronger than role
switching: does random coordinate order force within-law mixing, or otherwise
force a positive direct average of the orderwise local functionals?

This file maps the finite [4] case completely.

If coordinate i is deterministic from a predecessor set S in a union-closed
family, let P be the set of occurring prefix patterns and f:P->{0,1} its truth
function.  P is a join-semilattice, and union closure gives the exact identity

    f(x OR y) = f(x) OR f(y).

Thus f is a join-homomorphism; f^{-1}(0) is a prime join-ideal.  This is the
cross-coordinate information that an arbitrary sink law lacks.
Averaging over orders sometimes places i before an essential predecessor, making
its conditional law entropy-active.  This is necessary order dependence, but
entropy activation alone does not quantify functional cost: order averaging is
E_pi F(mu_{i,pi}), not F(E_pi mu_{i,pi}).

The enumeration:
  * all 65,535 nonempty families on [4]; retain union-closed;
  * identify separating families (distinct incidence columns) with every
    coordinate active (frequency in (0,1)); this is the standard minimal-
    counterexample reduction's relevant class;
  * for every family, coordinate and 4! order, construct the exact Fraction law
    of Pr(A_i=1 | predecessors);
  * for every nontrivial deterministic law, assert the join-homomorphism identity
    on every pair of occurring prefixes;
  * classify the zero-set prime ideal by its truth table, number of minimal
    1-patterns, and predecessor count; measure the fraction of all orders in
    which that coordinate is entropy-active and its Shapley entropy share.

STATUS: the join-homomorphism assertion is a finite exact proof for every
enumerated instance.  Counts and minima are COMPLETE FINITE DATA, not a
universal theorem for arbitrary n.

Run: ./.venv/bin/python uc/join_dependencies.py
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
from math import log2

N = 4
SUBSETS = 1 << N
ORDERS = tuple(permutations(range(N)))


def members(mask):
    return tuple(s for s in range(SUBSETS) if (mask >> s) & 1)


def union_closed(mem):
    p = set(mem)
    return all((x | y) in p for x in mem for y in mem)


def frequencies(mem):
    n = len(mem)
    return tuple(Fraction(sum((s >> i) & 1 for s in mem), n) for i in range(N))


def separating_active(mem):
    columns = [tuple((s >> i) & 1 for s in mem) for i in range(N)]
    freq = frequencies(mem)
    return len(set(columns)) == N and all(0 < x < 1 for x in freq)


def conditional_law(mem, order, position):
    coord = order[position]
    groups = {}
    for subset in mem:
        prefix = tuple((subset >> order[j]) & 1 for j in range(position))
        total, ones = groups.get(prefix, (0, 0))
        groups[prefix] = (total + 1, ones + ((subset >> coord) & 1))
    n = len(mem)
    mass = {}
    for total, ones in groups.values():
        atom = Fraction(ones, total)
        mass[atom] = mass.get(atom, Fraction()) + Fraction(total, n)
    return tuple(sorted(mass.items())), groups


def entropy_fraction_law(key):
    # Float only for the Shapley share; determinism uses exact atoms.
    import math
    total = 0.0
    for atom, weight in key:
        x = float(atom)
        hx = 0.0 if x in (0.0, 1.0) else -x * math.log2(x) - (1-x) * math.log2(1-x)
        total += float(weight) * hx
    return total


def truth_map(groups):
    assert all(ones in (0, total) for total, ones in groups.values())
    return {prefix: int(ones > 0) for prefix, (total, ones) in groups.items()}


def verify_join_homomorphism(truth):
    patterns = set(truth)
    for x in patterns:
        for y in patterns:
            join = tuple(a | b for a, b in zip(x, y))
            assert join in patterns
            assert truth[join] == (truth[x] | truth[y])


def minimal_one_patterns(truth):
    ones = [x for x, value in truth.items() if value]
    return tuple(x for x in ones
                 if not any(y != x and all(a <= b for a, b in zip(y, x))
                            for y in ones))


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    uc = []
    sep = []
    for mask in range(1, 1 << SUBSETS):
        mem = members(mask)
        if union_closed(mem):
            uc.append(mem)
            if separating_active(mem):
                sep.append(mem)

    instance_count = 0
    family_with_dependency = set()
    type_counts = Counter()
    mixing_records = []
    zero_active_order = []

    for family_id, mem in enumerate(sep):
        n = len(mem)
        marginal_h = []
        for f in frequencies(mem):
            x = float(f)
            import math
            marginal_h.append(-x*math.log2(x)-(1-x)*math.log2(1-x))

        # Collect per-coordinate order entropy, and exact deterministic cases.
        by_coordinate = defaultdict(list)
        deterministic_orders = defaultdict(int)
        for order in ORDERS:
            for pos, coord in enumerate(order):
                key, groups = conditional_law(mem, order, pos)
                ent = entropy_fraction_law(key)
                by_coordinate[coord].append(ent)
                if ent < 1e-14:
                    mean = sum(atom * weight for atom, weight in key)
                    if 0 < mean < 1:
                        truth = truth_map(groups)
                        verify_join_homomorphism(truth)
                        mins = minimal_one_patterns(truth)
                        instance_count += 1
                        deterministic_orders[coord] += 1
                        family_with_dependency.add(family_id)
                        # Classification is invariant enough for finite mapping:
                        # predecessor count, zero/one counts, minimal-one count,
                        # and sorted Hamming weights of minimal one-patterns.
                        signature = (
                            pos,
                            sum(v == 0 for v in truth.values()),
                            sum(v == 1 for v in truth.values()),
                            len(mins),
                            tuple(sorted(sum(x) for x in mins)),
                        )
                        type_counts[signature] += 1

        for coord, values in by_coordinate.items():
            if deterministic_orders[coord]:
                active_fraction = sum(v > 1e-14 for v in values) / len(values)
                shapley_entropy = sum(values) / len(values)
                share_ratio = shapley_entropy / marginal_h[coord]
                mixing_records.append((active_fraction, share_ratio,
                                       deterministic_orders[coord], family_id,
                                       coord, mem))
                if active_fraction == 0:
                    zero_active_order.append((family_id, coord, mem))

    print("COMPLETE FINITE COUNTS")
    print("  union-closed families on [4]:", len(uc))
    print("  separating + every coordinate active:", len(sep))
    print("  nontrivial deterministic order-instances:", instance_count)
    print("  separating families with at least one such dependency:",
          len(family_with_dependency))
    print("  distinct dependency signatures:", len(type_counts))
    print()

    print("PROVED [exact finite verification]")
    print("  f(x OR y)=f(x) OR f(y) held in every one of the", instance_count,
          "deterministic instances.")
    print("  zero-set truth classes are prime join-ideals in every instance.")
    print()

    if mixing_records:
        mixing_records.sort(key=lambda row: (row[0], row[1]))
        worst = mixing_records[0]
        print("RANDOM-ORDER ROLE SWITCHING IN THE SEPARATING CLASS")
        print("  minimum fraction of orders making a deterministic-dependent")
        print("  coordinate entropy-active: %.6f" % worst[0])
        print("  corresponding Shapley/marginal entropy share: %.6f" % worst[1])
        print("  deterministic orders: %d / %d" % (worst[2], len(ORDERS)))
        print("  family:", worst[5], "coordinate", worst[4])
        print("  minimum Shapley/marginal entropy ratio over all dependencies: %.6f" %
              min(row[1] for row in mixing_records))
        print("  maximum deterministic-order fraction: %.6f" %
              max(row[2]/len(ORDERS) for row in mixing_records))
        print("  dependent coordinates deterministic in every order:",
              len(zero_active_order))
        assert not zero_active_order
        print()

    print("MOST COMMON JOIN-DEPENDENCY SIGNATURES")
    print("  (predecessors, #zero, #one, #minimal-one, minimal-one weights): count")
    for signature, count in type_counts.most_common(12):
        print(" ", signature, ":", count)
    print()
    print("FINITE VERDICT")
    print("  In every separating active union-closed family on [4], each")
    print("  nontrivial deterministic dependency becomes entropy-active in a")
    print("  positive fraction of random coordinate orders.")
    print("  This proves role switching, not the certified mixture gap:")
    print("  averaging F over orders is not F of the averaged local law.")
    print("  The arbitrary-n frontier is a quantitative within-order mixing")
    print("  theorem or a direct lower bound for the permutation-averaged")
    print("  functional.  This enumeration is evidence, not either proof.")
    print("  No universal conclusion follows from n=4 alone.")
