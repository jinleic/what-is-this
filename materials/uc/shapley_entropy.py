"""RANDOM-ORDER / SHAPLEY ENTROPY experiment on every family over [4].

MOTIVATION
==========
`reimer_entropy.py` proves that Reimer's scalar average-size inequality cannot
break the c* ceiling in a coordinatewise Sawin/Cambie proof: for a FIXED
coordinate order, deterministic conditional laws a*delta_1+(1-a)*delta_0 have
zero entropy-functional cost but positive Reimer slack, and can compensate the
bad c* law globally.

That construction uses order-dependent determinism.  Average the entropy proof
over a uniformly random coordinate permutation.  Conditional entropy shares
then become the Shapley values of the entropy polymatroid: a variable that is
deterministic after its causes can have zero entropy share in that order, but
must inherit shared entropy when it precedes them.  The free slack may vanish.

This file tests that idea exhaustively on every family F subseteq 2^[4]:

  * identify whether F is union-closed;
  * for every one of 4! coordinate orders and every coordinate, compute the
    exact law mu_i of p_i = Pr(A_i=1 | preceding coordinates), A uniform on F;
  * compute Sawin/Cambie's local functional
        F(mu)=(1-alpha) E h(p+q-pq) + alpha min_coupling E h(s*(p,r)) - E h(p),
    where p,q are iid and the p,r term is minimised over ALL equal-marginal
    couplings by an exact finite transportation LP;
  * sum over coordinates; compare the fixed-order range and permutation
    average; record Reimer slack 2 E|A| - log_2|F|.

STATUS.  This is a COMPLETE finite enumeration, but it is NUMERICAL: scipy's
HiGHS solves the finite coupling LPs in float64.  It is a discovery probe, not
a universal proof.  Every local law is represented by exact Fraction weights
and atoms before the LP, and LP results are cached by that exact key.

Run: ./.venv/bin/python uc/shapley_entropy.py
"""

from fractions import Fraction
from itertools import permutations
from math import log2

import numpy as np
from scipy.optimize import linprog

N_GROUND = 4
N_SUBSETS = 1 << N_GROUND
ALPHA = 0.0356069
LOG2 = np.log(2.0)
ORDERS = tuple(permutations(range(N_GROUND)))


def h(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    mask = (x > 0.0) & (x < 1.0)
    y = x[mask]
    out[mask] = -(y * np.log(y) + (1.0 - y) * np.log1p(-y)) / LOG2
    return out


def sstar(p, r):
    return np.maximum(p, np.maximum(r, np.minimum(p + r, 0.5)))


def members(family_mask):
    return tuple(s for s in range(N_SUBSETS) if (family_mask >> s) & 1)


def union_closed(mem):
    present = set(mem)
    return all((a | b) in present for a in mem for b in mem)


def law_key(mem, order, position):
    """Exact law of p_i=Pr(bit i=1 | predecessor pattern)."""
    coordinate = order[position]
    groups = {}
    for subset in mem:
        prefix = tuple((subset >> order[j]) & 1 for j in range(position))
        total, ones = groups.get(prefix, (0, 0))
        groups[prefix] = (total + 1, ones + ((subset >> coordinate) & 1))
    n = len(mem)
    mass = {}
    for total, ones in groups.values():
        atom = Fraction(ones, total)
        mass[atom] = mass.get(atom, Fraction(0)) + Fraction(total, n)
    assert sum(mass.values(), Fraction(0)) == 1
    return tuple(sorted(mass.items()))


_FUNCTIONAL_CACHE = {}


def functional(key):
    """Worst-coupling local functional; cached by exact law key."""
    if key in _FUNCTIONAL_CACHE:
        return _FUNCTIONAL_CACHE[key]
    atoms = np.array([float(atom) for atom, _ in key])
    weights = np.array([float(weight) for _, weight in key])
    m = len(atoms)

    P, Q = np.meshgrid(atoms, atoms, indexing="ij")
    iid = float(weights @ h(P + Q - P * Q) @ weights)
    entropy = float(weights @ h(atoms))

    cost = h(sstar(P, Q)).ravel()
    aeq, beq = [], []
    for i in range(m):
        row = np.zeros((m, m))
        row[i, :] = 1.0
        aeq.append(row.ravel())
        beq.append(weights[i])
    for j in range(m):
        column = np.zeros((m, m))
        column[:, j] = 1.0
        aeq.append(column.ravel())
        beq.append(weights[j])
    result = linprog(cost, A_eq=np.array(aeq), b_eq=np.array(beq),
                     bounds=[(0.0, None)] * (m * m), method="highs")
    assert result.success, result.message
    coupled = float(result.fun)
    value = (1.0 - ALPHA) * iid + ALPHA * coupled - entropy
    mean = float(weights @ atoms)
    ans = (value, mean, entropy, atoms, weights)
    _FUNCTIONAL_CACHE[key] = ans
    return ans


def analyse(mem):
    order_totals = []
    entropy_total = None
    means = []
    deterministic_slacks = []
    for order in ORDERS:
        total = 0.0
        this_entropy = 0.0
        if not means:
            means = [None] * N_GROUND
        for position in range(N_GROUND):
            key = law_key(mem, order, position)
            value, mean, entropy, _, _ = functional(key)
            total += value
            this_entropy += entropy
            means[order[position]] = mean
            if entropy < 1e-14 and mean > 1e-14:
                deterministic_slacks.append(2.0 * mean)
        if entropy_total is None:
            entropy_total = this_entropy
        assert abs(this_entropy - log2(len(mem))) < 2e-12
        order_totals.append(total)
    average_size = sum(bin(s).count("1") for s in mem) / len(mem)
    reimer_slack = 2.0 * average_size - log2(len(mem))
    frequencies = [sum((s >> i) & 1 for s in mem) / len(mem)
                   for i in range(N_GROUND)]
    return {
        "min": min(order_totals),
        "max": max(order_totals),
        "average": sum(order_totals) / len(order_totals),
        "spread": max(order_totals) - min(order_totals),
        "reimer": reimer_slack,
        "max_frequency": max(frequencies),
        "deterministic_slack_count": len(deterministic_slacks),
        "deterministic_slack_total": sum(deterministic_slacks),
    }


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    records = []
    nonempty = 0
    uc_count = 0
    t0 = __import__("time").time()
    for mask in range(1, 1 << N_SUBSETS):
        mem = members(mask)
        nonempty += 1
        uc = union_closed(mem)
        if uc:
            uc_count += 1
        result = analyse(mem)
        result.update(mask=mask, size=len(mem), union_closed=uc, members=mem)
        records.append(result)

    uc_records = [r for r in records if r["union_closed"] and r["members"] != (0,)]
    arbitrary = [r for r in records if r["members"] != (0,)]

    print("NUMERICAL COMPLETE ENUMERATION")
    print("  ground-set size:", N_GROUND)
    print("  nonempty families:", nonempty)
    print("  union-closed families:", uc_count)
    print("  distinct exact conditional-law keys / LPs:", len(_FUNCTIONAL_CACHE))
    print("  runtime: %.2f s" % (__import__("time").time() - t0))
    print()

    def report(name, group):
        worst_avg = min(group, key=lambda r: r["average"])
        worst_min = min(group, key=lambda r: r["min"])
        most_spread = max(group, key=lambda r: r["spread"])
        print(name)
        print("  minimum permutation-AVERAGED sum F = %+.12e" % worst_avg["average"])
        print("    family", worst_avg["members"], "size", worst_avg["size"],
              "maxfreq %.3f" % worst_avg["max_frequency"],
              "Reimer slack %+.6f" % worst_avg["reimer"])
        print("  minimum over ANY FIXED order = %+.12e" % worst_min["min"])
        print("    average for that family = %+.12e" % worst_min["average"],
              "spread %.3e" % worst_min["spread"])
        print("  largest order spread = %.12e" % most_spread["spread"])
        print("    min/avg/max = %+.6e / %+.6e / %+.6e" %
              (most_spread["min"], most_spread["average"], most_spread["max"]))
        print("    family", most_spread["members"])
        print()

    report("ALL FAMILIES (control)", arbitrary)
    report("UNION-CLOSED FAMILIES", uc_records)

    # Compare fixed-order and Shapley averages by frequency bands.
    print("UNION-CLOSED, grouped by maximum-frequency threshold")
    print("  threshold  count    min fixed-order F    min Shapley-avg F")
    for threshold in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        group = [r for r in uc_records if r["max_frequency"] <= threshold + 1e-15]
        if group:
            print("  %-9.2f %-8d %+20.10e %+20.10e" %
                  (threshold, len(group), min(r["min"] for r in group),
                   min(r["average"] for r in group)))
        else:
            print("  %-9.2f %-8d %20s %20s" % (threshold, 0, "--", "--"))
    print()
    print("ALL FAMILIES satisfying the scalar Reimer cut, at target-scale frequencies")
    print("  threshold  count    min fixed-order F    min Shapley-avg F   Shapley gain")
    for threshold in (0.38235, 0.38261, 0.40, 0.45, 0.50):
        group = [r for r in arbitrary
                 if r["max_frequency"] <= threshold + 1e-15
                 and r["reimer"] >= -1e-12]
        if group:
            w = min(group, key=lambda r: r["average"])
            print("  %-9.5f %-8d %+20.10e %+20.10e %+12.4e" %
                  (threshold, len(group), min(r["min"] for r in group),
                   w["average"], w["average"] - w["min"]))
            print("    worst Shapley family", w["members"],
                  "maxfreq %.4f" % w["max_frequency"],
                  "Reimer %+.4e" % w["reimer"])
        else:
            print("  %-9.5f %-8d %20s %20s %12s" %
                  (threshold, 0, "--", "--", "--"))
    print()


    # Does permutation averaging systematically improve the worst order?
    gains = [r["average"] - r["min"] for r in uc_records]
    positive = sum(g > 1e-12 for g in gains)
    print("SHAPLEY EFFECT ON UNION-CLOSED FAMILIES")
    print("  avg(sum F) - min_order(sum F) is nonnegative tautologically;")
    print("  strictly positive (>1e-12) for %d / %d families" %
          (positive, len(uc_records)))
    print("  maximum gain from averaging: %.12e" % max(gains))
    print("  median gain: %.12e" % float(np.median(gains)))
    print()
    print("NUMERICAL ONLY: finite n=4 enumeration + float64 HiGHS couplings.")
    print("A positive Shapley effect identifies a cross-coordinate constraint worth")
    print("formalising; it is not a universal union-closed theorem by itself.")
