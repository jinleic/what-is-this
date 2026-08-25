"""COMPLETE NUMERICAL projection-closure probe for admissible [5] families.

For every one of the 64,534 nontrivial [5] families enumerated by
`shapley_n5_enumerate.cpp`, project the uniform distribution onto every nonempty
proper coordinate subset.  A projection is generally nonuniform: its atom
weights are exact fiber-size multiples of 1/|F|.  This is the first class of
nonuniform laws that an induction on a uniform family must handle.

The finite-data hypothesis tested here was projection nonnegativity. It holds
through [5], but `shapley_uniform_lift.py` refutes it for general admissible
uniform parents by lifting an arbitrary rational negative law. The value of
this script is therefore diagnostic: it maps the small projection lattice,
quantifies negative conditional fibers, and records where fixed-ratio and
uniformity-potential heuristics first looked plausible.

Evidence labels:
  * COMPLETE / EXACT: parent enumeration, projections, rational fiber weights.
  * COMPLETE NUMERICAL: float64 HiGHS transportation costs in A.

Build and run:
  c++ -O3 -std=c++17 uc/shapley_n5_enumerate.cpp -o /tmp/shapley_n5_enum
  ./.venv/bin/python uc/shapley_n5_projections.py /tmp/shapley_n5_enum
"""

from fractions import Fraction
from itertools import combinations
from math import comb, log2
from pathlib import Path
import subprocess
import sys
import time

from shapley_entropy import functional
from shapley_n5_complete import EXPECTED_TOTAL, exact_admissible, members

EXPECTED_UNIQUE = {1: 8, 2: 47, 3: 453, 4: 6035}


def project(mem, kept):
    counts = {}
    for subset in mem:
        image = sum(((subset >> old) & 1) << new
                    for new, old in enumerate(kept))
        counts[image] = counts.get(image, 0) + 1
    states = tuple(sorted(counts))
    probabilities = tuple(Fraction(counts[state], len(mem)) for state in states)
    return states, probabilities

def conditioned_fiber(mem, coordinate, bit):
    kept = tuple(j for j in range(5) if j != coordinate)
    fiber = tuple(subset for subset in mem
                  if ((subset >> coordinate) & 1) == bit)
    if not fiber:
        return None
    states = tuple(sorted(
        sum(((subset >> old) & 1) << new
            for new, old in enumerate(kept))
        for subset in fiber
    ))
    probabilities = tuple(Fraction(1, len(states)) for _state in states)
    return states, probabilities


def conditional_law(states, probabilities, coordinate, predecessors):
    predecessors = tuple(sorted(predecessors))
    groups = {}
    for state, probability in zip(states, probabilities):
        pattern = tuple((state >> j) & 1 for j in predecessors)
        total, ones = groups.get(pattern, (Fraction(), Fraction()))
        groups[pattern] = (
            total + probability,
            ones + probability * ((state >> coordinate) & 1),
        )
    mass = {}
    for total, ones in groups.values():
        atom = ones / total
        mass[atom] = mass.get(atom, Fraction()) + total
    return tuple(sorted(mass.items()))


def direct_value(states, probabilities, dimension):
    total = 0.0
    for coordinate in range(dimension):
        others = tuple(j for j in range(dimension) if j != coordinate)
        for size in range(dimension):
            weight = 1.0 / (dimension * comb(dimension - 1, size))
            for predecessors in combinations(others, size):
                law = conditional_law(
                    states, probabilities, coordinate, predecessors)
                total += weight * functional(law)[0]
    return total

def uniformity_defect(distribution):
    states, probabilities = distribution
    entropy = -sum(float(p) * log2(float(p)) for p in probabilities)
    return log2(len(states)) - entropy


def read_masks(binary):
    result = subprocess.run(
        [str(binary)], capture_output=True, text=True, check=True)
    masks = tuple(int(line, 16) for line in result.stdout.splitlines())
    assert len(masks) == len(set(masks)) == EXPECTED_TOTAL
    return masks


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/shapley_n5_enum")
    if not binary.is_file():
        raise SystemExit("missing enumerator binary; run the documented c++ command")
    masks = read_masks(binary)

    caches = {dimension: {} for dimension in range(1, 5)}
    minima = {}
    projection_count = 0
    negative_count = 0
    t0 = time.time()
    for mask in masks:
        mem = members(mask)
        assert exact_admissible(mem)
        for dimension in range(1, 5):
            for kept in combinations(range(5), dimension):
                projection_count += 1
                states, probabilities = project(mem, kept)
                key = (states, probabilities)
                cache = caches[dimension]
                if key not in cache:
                    cache[key] = direct_value(
                        states, probabilities, dimension)
                value = cache[key]
                if value < -1e-10:
                    negative_count += 1
                if len(states) > 1:
                    record = (value, mask, kept, states, probabilities)
                    if dimension not in minima or record < minima[dimension]:
                        minima[dimension] = record

    assert projection_count == EXPECTED_TOTAL * 30
    assert {dimension: len(cache) for dimension, cache in caches.items()} \
        == EXPECTED_UNIQUE
    assert negative_count == 0

    # Negative control: conditioning is not projection. Individual uniform
    # fibers can have large negative A, so induction needs a cross-fiber term.
    fiber_cache = {}
    fiber_minimum = None
    negative_fiber_count = 0
    for mask in masks:
        mem = members(mask)
        for coordinate in range(5):
            for bit in (0, 1):
                distribution = conditioned_fiber(mem, coordinate, bit)
                if distribution is None:
                    continue
                if distribution not in fiber_cache:
                    fiber_cache[distribution] = direct_value(
                        distribution[0], distribution[1], 4)
                value = fiber_cache[distribution]
                if value < -1e-10:
                    negative_fiber_count += 1
                record = (value, mask, coordinate, bit, distribution)
                if fiber_minimum is None or record < fiber_minimum:
                    fiber_minimum = record
    assert len(fiber_cache) == 8035
    assert negative_fiber_count == 250365
    assert fiber_minimum[0] < -0.7

    # Induction target inside the proper projection lattice: adding a kept
    # coordinate contracts A by at most 1/3. The full-parent level is covered
    # separately by shapley_n5_complete.py and is omitted here for runtime.
    edge_count = 0
    contraction_minimum = None
    ratio_minimum = None
    potential_minimum = None
    global_contraction_minimum = None
    global_potential_minimum = None
    for mask in masks:
        mem = members(mask)
        distributions = {}
        for dimension in range(1, 5):
            for kept in combinations(range(5), dimension):
                distributions[kept] = project(mem, kept)

        for dimension in range(2, 5):
            for kept in combinations(range(5), dimension):
                parent = distributions[kept]
                parent_value = caches[dimension][parent]
                for dropped in kept:
                    child_kept = tuple(j for j in kept if j != dropped)
                    child = distributions[child_kept]
                    child_value = caches[dimension - 1][child]
                    edge_count += 1
                    margin = parent_value - (2.0 / 3.0) * child_value
                    margin_record = (
                        margin, mask, kept, dropped,
                        parent_value, child_value, parent, child)
                    parent_potential = (
                        parent_value + min(uniformity_defect(parent), 0.1))
                    child_potential = (
                        child_value + min(uniformity_defect(child), 0.1))
                    potential_record = (
                        parent_potential - (2.0 / 3.0) * child_potential,
                        mask, kept, dropped, parent, child)
                    if (global_contraction_minimum is None
                            or margin_record < global_contraction_minimum):
                        global_contraction_minimum = margin_record
                    if (global_potential_minimum is None
                            or potential_record < global_potential_minimum):
                        global_potential_minimum = potential_record
                    if child_value <= 1e-12:
                        continue
                    ratio = parent_value / child_value
                    ratio_record = (
                        ratio, mask, kept, dropped,
                        parent_value, child_value, parent, child)
                    if (potential_minimum is None
                            or potential_record < potential_minimum):
                        potential_minimum = potential_record
                    if (contraction_minimum is None
                            or margin_record < contraction_minimum):
                        contraction_minimum = margin_record
                    if ratio_minimum is None or ratio_record < ratio_minimum:
                        ratio_minimum = ratio_record
    assert edge_count == EXPECTED_TOTAL * 70
    assert global_contraction_minimum[0] >= -1e-10
    assert global_potential_minimum[0] >= -1e-10
    assert contraction_minimum[0] > 0
    assert ratio_minimum[0] > 0.7
    assert potential_minimum[0] > 0

    print("COMPLETE NUMERICAL PROJECTION AUDIT")
    print("  parent families:", len(masks))
    print("  coordinate projections:", projection_count)
    print("  unique projected laws by dimension:",
          {dimension: len(cache) for dimension, cache in caches.items()})
    print("  values below -1e-10:", negative_count)
    print("  runtime: %.2f s" % (time.time() - t0))
    print()
    for dimension in range(1, 5):
        value, mask, kept, states, probabilities = minima[dimension]
        entropy = -sum(float(p) * log2(float(p)) for p in probabilities)
        uniformity_defect = log2(len(states)) - entropy
        print("  dimension %d minimum over nonconstant projections" % dimension)
        print("    A = %.15f; parent mask=0x%08x; kept=%s" %
              (value, mask, kept))
        print("    states=%s; probabilities=%s" % (states, probabilities))
        print("    uniformity defect log|supp|-H = %.15f" % uniformity_defect)
    print()
    print("CONDITIONAL-FIBER NEGATIVE CONTROL")
    print("  unique one-coordinate fibers:", len(fiber_cache))
    print("  negative fiber occurrences:", negative_fiber_count)
    print("  minimum A = %.15f; parent mask=0x%08x; conditioned x_%d=%d" %
          (fiber_minimum[0], fiber_minimum[1], fiber_minimum[2],
           fiber_minimum[3]))
    print("  distribution =", fiber_minimum[4])
    print("  Therefore induction cannot prove fibers separately; it needs the")
    print("  cross-fiber mixing/uniformity correction.")
    print()
    print("PROJECTION-LATTICE CONTRACTION PROBE")
    print("  projection edges:", edge_count)
    print("  minimum A(parent)/A(child) for A(child)>0: %.15f" %
          ratio_minimum[0])
    print("    parent / child distributions =", ratio_minimum[6], "/",
          ratio_minimum[7])
    print("  global minimum A(parent)-(2/3)A(child): %.15f" %
          global_contraction_minimum[0])
    print("  minimum for A(child)>0: %.15f" % contraction_minimum[0])
    print("  Truncated potential B=A+min(log|supp|-H,0.1):")
    print("    global minimum corrected margin: %.15f" %
          global_potential_minimum[0])
    print("    minimum for A(child)>0: %.15f" % potential_minimum[0])
    print("  Complete [5] data support both finite-data inequalities.")
    print()
    print("FINITE VERDICT")
    print("  Every proper projection in the complete [5] audit has A>=0.")
    print("  This is numerical small-data evidence only.")
    print("  General projection nonnegativity and the truncated potential are")
    print("  refuted by the balanced-label lift in shapley_uniform_lift.py.")