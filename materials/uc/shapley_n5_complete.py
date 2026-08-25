"""COMPLETE NUMERICAL direct-Shapley enumeration on [5].

The C++ companion enumerates exactly every nontrivial family satisfying

    max_i Pr(X_i=1) <= 0.38261,
    2 E|X| >= log_2 |F|.

This script independently checks every emitted family, brute-force matches the
size-3 and size-4 subcases, and evaluates the direct random-order functional
with exact Fraction conditional laws and float64 HiGHS transport LPs.

The computation refutes the normalized n=4 conjecture A(X)>=gamma_B H(X):
25 robust float64 violations (margin >1e-10) form two exact
coordinate-permutation orbits.  It does not refute the surviving nonnegativity
conjecture A(X)>=0.

Evidence labels:
  * COMPLETE / EXACT: admissible-family enumeration and symmetry classification.
  * COMPLETE NUMERICAL: functional comparison (float64 transport costs).

Build and run:
  c++ -O3 -std=c++17 uc/shapley_n5_enumerate.cpp -o /tmp/shapley_n5_enum
  ./.venv/bin/python uc/shapley_n5_complete.py /tmp/shapley_n5_enum
"""

from collections import Counter
from itertools import combinations, permutations
from math import log2
from pathlib import Path
import subprocess
import sys
import time

from shapley_direct import BASE4, direct_average
from shapley_entropy import union_closed

SPECS = {
    3: (1, 3),
    4: (1, 4),
    6: (2, 8),
    7: (2, 10),
    8: (3, 12),
    9: (3, 15),
    11: (4, 20),
}
EXPECTED_COUNTS = {
    3: 145,
    4: 70,
    6: 2833,
    7: 296,
    8: 44120,
    9: 2910,
    11: 14160,
}
EXPECTED_TOTAL = 64534


def members(family_mask):
    return tuple(subset for subset in range(32)
                 if (family_mask >> subset) & 1)


def exact_admissible(mem):
    cap, incidence_need = SPECS[len(mem)]
    coordinate_counts = [
        sum((subset >> coordinate) & 1 for subset in mem)
        for coordinate in range(5)
    ]
    return max(coordinate_counts) <= cap and sum(coordinate_counts) >= incidence_need


def family_mask(mem):
    return sum(1 << subset for subset in mem)


def columns(mem):
    return tuple(
        tuple((subset >> coordinate) & 1 for subset in mem)
        for coordinate in range(5)
    )


def permute_subset(subset, permutation):
    image = 0
    for old, new in enumerate(permutation):
        if (subset >> old) & 1:
            image |= 1 << new
    return image


def permute_family(mask, permutation):
    image = 0
    for subset in range(32):
        if (mask >> subset) & 1:
            image |= 1 << permute_subset(subset, permutation)
    return image


def orbit(mask):
    return {permute_family(mask, permutation)
            for permutation in permutations(range(5))}


def read_enumerator(binary):
    result = subprocess.run(
        [str(binary)], capture_output=True, text=True, check=True)
    masks = tuple(int(line, 16) for line in result.stdout.splitlines())
    return masks, result.stderr


def brute_small(size):
    found = set()
    for mem in combinations(range(32), size):
        if exact_admissible(mem):
            found.add(family_mask(mem))
    return found


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/shapley_n5_enum")
    if not binary.is_file():
        raise SystemExit("missing enumerator binary; run the documented c++ command")

    masks, enumerator_log = read_enumerator(binary)
    counts = Counter(mask.bit_count() for mask in masks)
    assert len(masks) == len(set(masks)) == EXPECTED_TOTAL
    assert dict(counts) == EXPECTED_COUNTS
    for mask in masks:
        assert exact_admissible(members(mask))
    for size in (3, 4):
        emitted = {mask for mask in masks if mask.bit_count() == size}
        assert emitted == brute_small(size)

    print("COMPLETE EXACT ENUMERATION")
    print(enumerator_log.strip())
    print("  parsed unique families:", len(masks))
    print("  Python constraint check: PASS")
    print("  independent brute checks m=3,4: PASS")
    print()

    base_value, _ = direct_average(BASE4, 4)
    gamma_b = base_value / 3.0
    records = []
    t0 = time.time()
    for index, mask in enumerate(masks, 1):
        mem = members(mask)
        value, layers = direct_average(mem, 5)
        entropy = log2(len(mem))
        incidence = sum(subset.bit_count() for subset in mem)
        reimer = 2.0 * incidence / len(mem) - entropy
        separating = len(set(columns(mem))) == 5
        records.append({
            "mask": mask,
            "members": mem,
            "value": value,
            "ratio": value / entropy,
            "layers": layers,
            "reimer": reimer,
            "separating": separating,
            "union_closed": union_closed(mem),
        })
        if index % 10000 == 0:
            print("  evaluated %d / %d" % (index, len(masks)))

    worst = min(records, key=lambda record: record["ratio"])
    separating_records = [record for record in records if record["separating"]]
    worst_separating = min(separating_records,
                           key=lambda record: record["ratio"])
    equality_records = [record for record in records
                        if abs(record["reimer"]) < 1e-12]
    worst_equality = min(equality_records, key=lambda record: record["ratio"])
    robust_violations = [record for record in records
                         if record["ratio"] < gamma_b - 1e-10]
    boundary_band = [record for record in records
                     if abs(record["ratio"] - gamma_b) <= 1e-10]
    lower_boundary_band = [record for record in boundary_band
                           if record["ratio"] < gamma_b]
    robustly_above = [record for record in records
                      if record["ratio"] > gamma_b + 1e-10]
    nonseparating_violations = [record for record in robust_violations
                                if not record["separating"]]
    separating_violations = [record for record in robust_violations
                             if record["separating"]]

    global_ties = {record["mask"] for record in records
                   if abs(record["ratio"] - worst["ratio"]) < 1e-12}
    separating_ties = {record["mask"] for record in separating_records
                       if abs(record["ratio"]
                              - worst_separating["ratio"]) < 1e-12}
    assert len(robust_violations) == 25
    assert len(nonseparating_violations) == 10
    assert len(separating_violations) == 15
    assert global_ties == orbit(worst["mask"])
    assert separating_ties == orbit(worst_separating["mask"])
    assert {record["mask"] for record in nonseparating_violations} \
        == orbit(worst["mask"])
    assert {record["mask"] for record in separating_violations} \
        == orbit(worst_separating["mask"])
    assert {record["mask"] for record in boundary_band} \
        == orbit(family_mask(BASE4))
    assert len(robust_violations) + len(boundary_band) \
        + len(robustly_above) == len(records)
    assert abs(worst_equality["ratio"] - gamma_b) < 1e-12
    assert min(record["value"] for record in records) > 0
    assert not any(record["union_closed"] for record in records)

    print()
    print("COMPLETE NUMERICAL FUNCTIONAL EVALUATION")
    print("  runtime: %.2f s" % (time.time() - t0))
    print("  float64 LP gamma_B benchmark: %.15f" % gamma_b)
    print("  robust gamma_B violations (margin >1e-10):",
          len(robust_violations))
    print("    nonseparating / separating: %d / %d" %
          (len(nonseparating_violations), len(separating_violations)))
    print("  float64 boundary band |A/H-gamma_B|<=1e-10:",
          len(boundary_band))
    print("    exactly the dummy-coordinate orbit of BASE4: PASS")
    print("    exact A/H=gamma_B follows from deterministic-product tensorization")
    print("    float64 values landing microscopically below benchmark:",
          len(lower_boundary_band), "(classified as exact equality)")
    print("  robust violations are exactly two coordinate-permutation orbits: PASS")
    print()
    print("  global minimum")
    print("    A/H = %.15f; A = %.15f" %
          (worst["ratio"], worst["value"]))
    print("    mask = 0x%08x; family = %s" %
          (worst["mask"], worst["members"]))
    print("    Reimer slack = %.12f; distinct columns = %d" %
          (worst["reimer"], len(set(columns(worst["members"])))))
    print("    predecessor layers =", worst["layers"])
    print()
    print("  separating minimum")
    print("    A/H = %.15f; A = %.15f" %
          (worst_separating["ratio"], worst_separating["value"]))
    print("    mask = 0x%08x; family = %s" %
          (worst_separating["mask"], worst_separating["members"]))
    print("    Reimer slack = %.12f" % worst_separating["reimer"])
    print("    predecessor layers =", worst_separating["layers"])
    print()
    print("  Reimer-equality subcase")
    print("    family count:", len(equality_records))
    print("    minimum A/H = %.15f (matches gamma_B)" %
          worst_equality["ratio"])
    print()
    print("FINAL VERDICT")
    print("  REFUTED: A(X) >= gamma_B H(X) under Reimer + marginal cap.")
    print("  Every enumerated nontrivial family numerically has A(X) > 0.")
    print("  The omitted trivial family {emptyset} is admissible, union-closed,")
    print("  and has A(X)=0; hence the surviving conjecture remains A(X)>=0.")
    print("  No nontrivial admissible [5] family is union-closed; this is a")
    print("  complete entropy-relaxation stress test, not a union-closed theorem.")
