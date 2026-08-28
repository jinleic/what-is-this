#!/usr/bin/env python3
"""Admissibility of every Cartesian power of the cloned-coordinate base.

`n8clone_hi` is the growth base: its certified threshold `-3/80` over dimension
8 gives slope `3/640 = 0.0046875`, above the `n=7` base's `1/250 = 0.004`.  That
claim rests on the powers being admissible, so this module checks it two ways.

Exact integer argument, valid for every `k`
-------------------------------------------
Let `F` have `m` rows on `d` coordinates with every degree equal to `c` and
incidence `I`.  In `F^k` the size is `m^k`, a coordinate lies in `c*m^(k-1)`
rows, and the incidence is `k*I*m^(k-1)`.

*Cap.*  For `m = 70` the cap is `floor(2*70/5) = 28`, and every degree already
equals it.  In the power, `2*m^k/5 = 28*70^(k-1)` exactly, because `70/5 = 14`
is an integer.  So the power's degrees equal its cap exactly, for every `k`, with
no rounding slack consumed.

*Reimer.*  The requirement is `k*I*m^(k-1) >= ceil(m^k * log2(m^k) / 2)`.  Since
the left side is an integer, it suffices that it strictly exceed the real right
side, i.e. that `I/m > log2(m)/2`, i.e. `2I/m > log2 m`, i.e. `2^(2I/m) > m`.
With `I = 224` and `m = 70` that is `2^6.4 > 70`, i.e. the integer inequality
`70^5 = 1,680,700,000 < 2^32 = 4,294,967,296`, which is the base's recorded
`reimer_witness`.  So Reimer holds strictly in every power.

Neither step uses separation, and neither uses the cell structure.  This is the
same computation that `PROOF.md` performs for the 45-row base; the point here is
that a cloned-coordinate base satisfies it too.

Direct instantiation
--------------------
The argument above is then checked against reality by building `F^2` as plain
16-bit masks -- 4,900 rows -- and recomputing size, degrees, cap, incidence,
Reimer threshold, activity, separation and the missing-ordered-join count by
brute force, with no product-aware shortcut.  The join count must match
`1-(1-eps)^2` exactly.  `k=3` is not instantiated: `343,000^2 = 1.18e11` ordered
pairs is out of proportion to what it would add, and the integer argument already
covers every `k`.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/audit_clone_power.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_gate_b_rational import BASES, exact_reimer_threshold  # noqa: E402

DEFAULT_REPORT = HERE / "candidates" / "clone_power_audit.json"


def product_rows(rows: tuple[int, ...], dimension: int, power: int) -> tuple[int, ...]:
    """Rows of F^power, placing each block in its own coordinate range."""
    out = (0,)
    for block in range(power):
        shift = block * dimension
        out = tuple(prefix | (row << shift) for prefix in out for row in rows)
    return tuple(sorted(out))


def missing_ordered_joins(rows: tuple[int, ...], dimension: int) -> int:
    present = bytearray(1 << dimension)
    for row in rows:
        present[row] = 1
    missing = 0
    for left in rows:
        for right in rows:
            if not present[left | right]:
                missing += 1
    return missing


def integer_power_admissibility(
    size: int, incidence: int, degree: int, dimension: int, powers: range
) -> dict[str, object]:
    """Cap and Reimer for F^k by exact integer arithmetic, for each k."""
    rows = []
    for power in powers:
        power_size = size**power
        power_degree = degree * size ** (power - 1)
        power_cap = 2 * power_size // 5
        power_incidence = power * incidence * size ** (power - 1)
        threshold = exact_reimer_threshold(power_size) if power <= 2 else None
        record = {
            "power": power,
            "dimension": dimension * power,
            "size": power_size,
            "degree": power_degree,
            "cap": power_cap,
            "degree_equals_cap": power_degree == power_cap,
            "cap_ok": power_degree <= power_cap,
            "incidence": power_incidence,
        }
        if threshold is not None:
            record["reimer_threshold"] = threshold
            record["reimer_ok"] = power_incidence >= threshold
        else:
            # For larger powers the exact threshold would need size**size, so
            # use the threshold-free form.  Reimer for F^k asks
            # k*I*m^(k-1) >= ceil(k * m^k * log2(m) / 2); the left side is an
            # integer, so it suffices that it dominate the real right side,
            # which cancels k and m^(k-1) down to 2I >= m*log2(m), i.e. to the
            # integer inequality m^m <= 2^(2I).  That is exactly the base's own
            # Reimer condition I >= R_m, since R_m is the least r with
            # 2^(2r) >= m^m.  One base-level check therefore certifies every
            # power at once, for any admissible base.
            record["reimer_ok"] = size**size <= (1 << (2 * incidence))
            record["reimer_threshold_free_witness"] = (
                f"{size}**{size} <= 2**{2 * incidence}"
            )
        rows.append(record)
    return {"powers": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=str, default="n8clone_hi")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-instantiated-power", type=int, default=2)
    parser.add_argument("--max-integer-power", type=int, default=40)
    args = parser.parse_args()

    base = BASES[args.base]
    facts = base.exact_facts()
    rows = base.reconstruct()
    dimension = base.dimension
    size = len(rows)
    degree = facts["coordinate_counts"][0]
    incidence = facts["total_incidence"]
    defect = Fraction(facts["closure_defect"])
    success = 1 - defect

    started = time.monotonic()
    integer_block = integer_power_admissibility(
        size, incidence, degree, dimension, range(1, args.max_integer_power + 1)
    )

    instantiated = []
    for power in range(1, args.max_instantiated_power + 1):
        power_dimension = dimension * power
        built = product_rows(rows, dimension, power)
        counts = tuple(
            sum((row >> coordinate) & 1 for row in built)
            for coordinate in range(power_dimension)
        )
        columns = tuple(
            tuple((row >> coordinate) & 1 for row in built)
            for coordinate in range(power_dimension)
        )
        built_incidence = sum(counts)
        threshold = exact_reimer_threshold(len(built))
        missing = missing_ordered_joins(built, power_dimension)
        observed = Fraction(missing, len(built) ** 2)
        predicted = 1 - success**power
        instantiated.append(
            {
                "power": power,
                "dimension": power_dimension,
                "size": len(built),
                "distinct_rows": len(set(built)) == len(built),
                "coordinate_counts_all_equal": len(set(counts)) == 1,
                "degree": counts[0],
                "cap": 2 * len(built) // 5,
                "degree_equals_cap": counts[0] == 2 * len(built) // 5,
                "total_incidence": built_incidence,
                "reimer_threshold": threshold,
                "reimer_ok": built_incidence >= threshold,
                "active": all(any(column) for column in columns),
                "separating": len(set(columns)) == power_dimension,
                "duplicate_column_pairs": sorted(
                    [i, j]
                    for i in range(power_dimension)
                    for j in range(i + 1, power_dimension)
                    if columns[i] == columns[j]
                ),
                "missing_ordered_join_pairs": missing,
                "closure_defect": str(observed),
                "predicted_closure_defect": str(predicted),
                "defect_identity_holds": observed == predicted,
            }
        )

    report = {
        "base": args.base,
        "base_dimension": dimension,
        "base_size": size,
        "base_degree": degree,
        "base_incidence": incidence,
        "base_closure_defect": str(defect),
        "base_separating": facts["separating"],
        "base_duplicate_column_pairs": facts["duplicate_column_pairs"],
        "reimer_integer_witness": base.reimer_witness,
        "reimer_witness_holds": size**5 < (1 << (4 * dimension)),
        "integer_argument": integer_block,
        "instantiated": instantiated,
        "all_integer_powers_admissible": all(
            record["cap_ok"] and record.get("reimer_ok", True)
            for record in integer_block["powers"]
        ),
        "every_power_degree_equals_cap": all(
            record["degree_equals_cap"] for record in integer_block["powers"]
        ),
        "all_instantiated_admissible": all(
            record["cap_ok"] if "cap_ok" in record else record["degree_equals_cap"]
            for record in instantiated
        )
        and all(record["reimer_ok"] for record in instantiated)
        and all(record["defect_identity_holds"] for record in instantiated),
        # No timing here: the artifact must hash identically on every rerun.
    }
    report["verdict"] = (
        "EVERY_POWER_ADMISSIBLE"
        if report["all_integer_powers_admissible"]
        and report["all_instantiated_admissible"]
        and report["reimer_witness_holds"]
        else "FAILED"
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(json.dumps({k: report[k] for k in (
        "base", "verdict", "reimer_witness_holds",
        "all_integer_powers_admissible", "every_power_degree_equals_cap",
        "all_instantiated_admissible")}, indent=2))
    print(f"wall seconds: {time.monotonic() - started:.2f}", file=sys.stderr)
    for record in instantiated:
        print(
            f"k={record['power']}: n={record['dimension']} m={record['size']} "
            f"deg={record['degree']}=cap {record['degree_equals_cap']} "
            f"inc={record['total_incidence']}>=R={record['reimer_threshold']} "
            f"{record['reimer_ok']} separating={record['separating']} "
            f"defect={record['closure_defect']} identity={record['defect_identity_holds']}"
        )


if __name__ == "__main__":
    main()
