#!/usr/bin/env python3
"""Gate B against the entropy-method ceiling: an exact finite barrier and a hard limit.

The source bridge
-----------------
The frozen Gate B objective

    A_+(F) = (1-alpha) Q(F) + alpha C_+(F) - log2 m,
    alpha = 356069 / 10^7,

comes from the same two-strategy entropy program as Sawin's Question and
Cambie's Question 2 (arXiv:2212.12500v2):

* ``Q`` is exactly their iid-OR term ``E h(p+q-pq)``;
* the lower endpoint of every Gate B Bellman action is exactly Cambie's
  ``s*(p,r) = max(p,r,min(p+r,1/2))``;
* ``C_+`` is **not** just that one greedy term: it maximizes the total entropy
  over every feasible causal continuation.  Hence it is an upper envelope of
  the Sawin-Cambie dependent strategy;
* ``log2 m`` is the chain-rule sum ``E h(p)`` for the uniform law on ``F``.

Thus a valid Sawin-Cambie inequality implies ``A_+ >= 0``, while ``A_+ < 0`` is
the stronger statement that even the causal Bellman upper envelope fails.  The
repository's rational ``alpha=0.0356069`` is the seven-decimal value printed by
Cambie for the optimal mixing weight; the undocumented provenance inside
``gate_b`` is an [INFERENCE], while the formula match is an in-repo derivation.

Status of the published ceiling
-------------------------------
Cambie's Section 2.3 gives a two-atom upper obstruction at

    c* = 0.3823455333667027...

and Section 3 claims the matching lower direction.  The broader ``math/uc``
audit establishes the upper obstruction but labels the matching lower check
**OPEN / COMPUTATIONAL-EVIDENCE**: Cambie explicitly calls its optimization
\"slightly less rigorous\" and uses finite precision plus graphical
confirmation.  This module therefore does **not** promote ``c*`` to a proved
frequency window.  It does three narrower, rigorous things:

1. encloses the defining two-atom obstruction ``(b,a,c*)`` in exact rational
   arithmetic, independently reproducing the reported decimal;
2. verifies that all ten registered Gate B negative bases have maximum
   frequency exactly ``2/5``, beyond the obstruction by at least
   ``0.017654466633``;
3. computes the exact finite failure of the iid entropy chain described next.

The order-extremal finite barrier
--------------------------------
Gilmer's chain lower-bounds ``H(X u Y)`` by the conditional-OR entropy sum for a
*fixed* coordinate order, so it can only prove ``H(X u Y) > H(X)`` through that
sum when ``Q_pi > log2 m`` for some order ``pi``.  This module computes both
order extremes of ``Q_pi`` by dynamic programming over predecessor sets rather
than enumerating ``n!`` orders, and rigorously tests

    max_pi Q_pi < log2 m.

Six of the ten registered bases pass: the one-step iid chain fails at **every**
coordinate order.  The widest exact margin is at least ``0.0393134578``.  This
is a finite uniform-family obstruction at frequency exactly ``2/5``; it does
not rely on Cambie's unverified matching lower check.

Why the requested cap-only defect floor is false
------------------------------------------------
Chase-Lovett (arXiv:2211.11689v1, Example 1.4) prove [REPORTED] that

    F = {x : |x| = psi n + n^(2/3)} u {x : |x| >= (1-psi)n},
    psi = (3-sqrt5)/2,

is ``1-o(1)``-approximately union closed -- exactly ``eps_vee=o(1)`` under the
Gate B ordered-pair convention -- while every frequency is
``psi+o(1) < 2/5``.  It contains the dominant upper layer.  Therefore **no
positive defect floor follows from the cap, even in the presence of a dominant
set.**

This is not a Gate B counterexample because it violates Reimer by a linear
amount: asymptotically its mean set size is ``(psi+o(1))n``, while
``log2|F|/2 = (h(psi)/2+o(1))n = (0.4797...+o(1))n``.  The exact leading-rate
gap enclosed here is at least ``0.0977433528 n``.  Reimer is therefore the only
remaining lever: any floor for the actual Gate B admissible class must use it
essentially.  Chase-Lovett also prove their ``psi`` endpoint optimal for
pair-defect-only arguments, so this is a hard published blocker, not a failed
search.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/barrier_sawin_cambie.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_gate_b_rational import (  # noqa: E402
    ACCUMULATOR_BITS,
    ALPHA,
    BASES,
    decimal_lower,
    decimal_upper,
    entropy_lower,
    entropy_upper,
    log2_lower,
    log2_upper,
    round_down,
    round_up,
)
from bound_local_regime import (  # noqa: E402
    admissible,
    cap,
    defect,
    degrees,
    incidence,
    reimer_threshold,
)

DEFAULT_REPORT = HERE / "candidates" / "sawin_cambie_barrier.json"

# psi = (3 - sqrt 5)/2, the Gilmer/AHS threshold, as a safe rational bracket.
PSI_LOWER = Fraction(3819660112501051, 10**16)
PSI_UPPER = Fraction(3819660112501052, 10**16)

# Exact candidate threshold from the broader math/uc certificate campaign.  Its
# five-dimensional interval certificate is machine-replayed; the universal
# entropy bridge is human-audited, so math/uc labels the resulting theorem
# CANDIDATE rather than promoting Cambie's reported c*.
T_CERT = Fraction(955165028125263, 2500000000000000)


# --- Cambie's ceiling, enclosed without floating point ----------------------


def _entropy_pair(value: Fraction) -> tuple[Fraction, Fraction]:
    return entropy_lower(value), entropy_upper(value)


def _ceiling_equation(point: Fraction) -> tuple[Fraction, Fraction]:
    """Enclosure of f(x) = h(x)(2 - h(x)) - h(2x - x^2), h in bits.

    Cambie's extremal two-point law is supported on {1, b} where f(b) = 0 and b
    is the larger root.  ``t (2 - t)`` is increasing for ``t < 1``, so the product
    term is enclosed by evaluating at the endpoints of the ``h`` enclosure.
    """
    low, high = _entropy_pair(point)
    image = 2 * point - point * point
    image_low, image_high = _entropy_pair(image)
    return low * (2 - low) - image_high, high * (2 - high) - image_low


def cambie_ceiling(steps: int = 80) -> dict[str, object]:
    """Exact rational enclosure of Cambie's ceiling c and its extremal law."""
    low, high = Fraction(1, 4), Fraction(7, 20)
    low_sign = _ceiling_equation(low)
    high_sign = _ceiling_equation(high)
    assert low_sign[1] < 0 < high_sign[0], (low_sign, high_sign)
    for _ in range(steps):
        middle = (low + high) / 2
        enclosure = _ceiling_equation(middle)
        if enclosure[0] > 0:
            high = middle
        elif enclosure[1] < 0:
            low = middle
        else:  # sign unresolvable at this precision; the bracket is the answer
            break
    # h is increasing on [0, 1/2] and b < 1/2, so the enclosure is monotone.
    entropy_low, entropy_high = entropy_lower(low), entropy_upper(high)
    # a = (1 - h(b)) / (2 - h(b)) is decreasing in h(b).
    weight_low = (1 - entropy_high) / (2 - entropy_high)
    weight_high = (1 - entropy_low) / (2 - entropy_low)
    # c = a (1 - b) + b is increasing in both a and b.
    ceiling_low = weight_low * (1 - low) + low
    ceiling_high = weight_high * (1 - high) + high
    return {
        "atom_lower": low,
        "atom_upper": high,
        "weight_lower": weight_low,
        "weight_upper": weight_high,
        "ceiling_lower": ceiling_low,
        "ceiling_upper": ceiling_high,
    }


# --- order-extremal conditional-OR entropy ----------------------------------


def conditional_law(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> tuple[Fraction, list[tuple[Fraction, Fraction]]]:
    """Re-derived conditional law of one coordinate after a predecessor SET.

    Returns ``(zero_mass, impure)`` where ``zero_mass`` is the total mass of
    prefixes forcing the coordinate to 0 and ``impure`` lists ``(p, mass)`` for
    prefixes with ``0 < p < 1``.  Prefixes forcing a 1 are dropped: their OR is
    certain, so they contribute no entropy against any partner.
    """
    groups: dict[int, list[int]] = {}
    mask = 0
    for predecessor in predecessors:
        mask |= 1 << predecessor
    bit = 1 << coordinate
    for row in family:
        record = groups.setdefault(row & mask, [0, 0])
        record[0] += 1
        if row & bit:
            record[1] += 1
    size = len(family)
    zero_mass = Fraction(0)
    impure: list[tuple[Fraction, Fraction]] = []
    for total, ones in groups.values():
        if ones == 0:
            zero_mass += Fraction(total, size)
        elif ones < total:
            impure.append((Fraction(ones, total), Fraction(total, size)))
    impure.sort()
    return zero_mass, impure


def coordinate_term(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> tuple[Fraction, Fraction]:
    """Enclosure of E_{u,v} h(p(u) + p(v) - p(u) p(v)) for one coordinate."""
    zero_mass, impure = conditional_law(family, coordinate, predecessors)
    if not impure:
        return Fraction(0), Fraction(0)
    low = high = Fraction(0)
    if zero_mass:
        twice = 2 * zero_mass
        for probability, mass in impure:
            weight = twice * mass
            low += entropy_lower(probability) * weight
            high += entropy_upper(probability) * weight
    for left, left_mass in impure:
        for right, right_mass in impure:
            union = left + right - left * right
            weight = left_mass * right_mass
            low += entropy_lower(union) * weight
            high += entropy_upper(union) * weight
    return low, high


def order_extremal_q(family: tuple[int, ...], dimension: int) -> dict[str, object]:
    """One-sided bounds on max_pi Q_pi and min_pi Q_pi, by DP over predecessor sets.

    ``Q_pi`` is a sum of per-coordinate terms indexed by (coordinate, predecessor
    SET), so the order extremes are longest/shortest paths through the subset
    lattice: 2^n states rather than n! orders.  The returned ``maximum_upper`` is
    a rigorous upper bound on max_pi Q_pi and ``minimum_lower`` a rigorous lower
    bound on min_pi Q_pi -- each the direction a barrier claim needs.
    """
    full = (1 << dimension) - 1
    upper: list[Fraction | None] = [None] * (1 << dimension)
    lower: list[Fraction | None] = [None] * (1 << dimension)
    choice_up: list[int] = [-1] * (1 << dimension)
    choice_down: list[int] = [-1] * (1 << dimension)
    upper[full] = Fraction(0)
    lower[full] = Fraction(0)
    for state in sorted(range(1 << dimension), key=lambda s: -bin(s).count("1")):
        if state == full:
            continue
        predecessors = tuple(i for i in range(dimension) if (state >> i) & 1)
        best_up: Fraction | None = None
        best_down: Fraction | None = None
        for coordinate in range(dimension):
            if (state >> coordinate) & 1:
                continue
            term_low, term_high = coordinate_term(family, coordinate, predecessors)
            following = state | (1 << coordinate)
            candidate_up = round_up(term_high + upper[following], ACCUMULATOR_BITS)
            candidate_down = round_down(term_low + lower[following], ACCUMULATOR_BITS)
            if best_up is None or candidate_up > best_up:
                best_up, choice_up[state] = candidate_up, coordinate
            if best_down is None or candidate_down < best_down:
                best_down, choice_down[state] = candidate_down, coordinate
        upper[state] = best_up
        lower[state] = best_down

    def walk(choices: list[int]) -> tuple[int, ...]:
        state, order = 0, []
        while state != full:
            coordinate = choices[state]
            order.append(coordinate)
            state |= 1 << coordinate
        return tuple(order)

    return {
        "maximum_upper": upper[0],
        "minimum_lower": lower[0],
        "argmax_order": walk(choice_up),
        "argmin_order": walk(choice_down),
    }


# --- the barrier statement --------------------------------------------------


def barrier_facts(family: tuple[int, ...], dimension: int) -> dict[str, object]:
    """Does the one-step entropy chain fail at every coordinate order?"""
    size = len(family)
    counts = degrees(family, dimension)
    frequency = Fraction(max(counts), size)
    extremes = order_extremal_q(family, dimension)
    size_log_low = log2_lower(Fraction(size))
    size_log_high = log2_upper(Fraction(size))
    maximum = extremes["maximum_upper"]
    minimum = extremes["minimum_lower"]
    margin = size_log_low - maximum
    return {
        "dimension": dimension,
        "family_size": size,
        "maximum_frequency": str(frequency),
        "maximum_frequency_decimal": decimal_upper(frequency, 10),
        "closure_defect": str(defect(family)),
        "admissible": admissible(family, dimension),
        "log2_size_lower_rational": str(size_log_low),
        "log2_size_upper_rational": str(size_log_high),
        "log2_size_lower": decimal_lower(size_log_low, 10),
        "order_maximum_upper_rational": str(maximum),
        "order_maximum_upper": decimal_upper(maximum, 10),
        "order_minimum_lower_rational": str(minimum),
        "order_minimum_lower": decimal_lower(minimum, 10),
        "chain_fails_at_every_order": maximum < size_log_low,
        "barrier_margin_lower_rational": str(margin),
        "barrier_margin_lower": decimal_lower(margin, 10),
        "argmax_order": list(extremes["argmax_order"]),
        "_maximum_upper": maximum,
        "_log2_lower": size_log_low,
        "_log2_upper": size_log_high,
    }


def excluded_sizes(limit: int, threshold: Fraction) -> tuple[int, ...]:
    """Sizes that a valid functional inequality at ``threshold`` would exclude.

    The arithmetic is exact.  The logical status is inherited from the supplied
    threshold: CANDIDATE for ``T_CERT`` and REPORTED/OPEN for Cambie's ``c*``.
    """
    return tuple(
        size
        for size in range(3, limit + 1)
        if Fraction(cap(size), size) <= threshold
    )


def chase_lovett_reimer_gap() -> dict[str, object]:
    """Leading-rate incompatibility of Chase-Lovett Example 1.4 with Reimer."""
    psi_low, psi_high = PSI_LOWER, PSI_UPPER
    entropy_low, entropy_high = entropy_lower(psi_low), entropy_upper(psi_high)
    required_low = entropy_low / 2
    deficit_low = required_low - psi_high
    return {
        "source_status": "REPORTED: Chase-Lovett arXiv:2211.11689v1 Example 1.4",
        "rates_are_asymptotic_leading_terms": True,
        "psi_lower_rational": str(psi_low),
        "psi_upper_rational": str(psi_high),
        "entropy_at_psi_lower_rational": str(entropy_low),
        "entropy_at_psi_upper_rational": str(entropy_high),
        "reimer_required_rate_lower_rational": str(required_low),
        "reimer_deficit_rate_lower_rational": str(deficit_low),
        "psi_upper_decimal": decimal_upper(psi_high, 10),
        "asymptotic_mean_size_rate_upper_decimal": decimal_upper(psi_high, 10),
        "reimer_required_rate_lower_decimal": decimal_lower(required_low, 10),
        "reimer_deficit_rate_lower_decimal": decimal_lower(deficit_low, 10),
        "reimer_satisfied_at_leading_order": psi_high >= required_low,
        "entropy_at_psi_upper_decimal": decimal_upper(entropy_high, 10),
    }

def reimer_endpoint_facts() -> dict[str, object]:
    """Exact logical status of Reimer and the dominant set at zero defect."""
    return {
        "reimer_source": (
            "CITED-DEPENDENCY: D. Reimer, An Average Set Size Theorem, "
            "Combinatorics, Probability and Computing 12(1) (2003), 89-93, "
            "doi:10.1017/S0963548302005230"
        ),
        "zero_defect_means_union_closed": True,
        "reimer_is_automatic_at_zero_defect": True,
        "dominant_set_is_automatic_under_cap": True,
        "dominant_set_derivation": (
            "Union closure contains the union T of all rows. On the active "
            "support, |T|=n and cap gives mean size <=2n/5, hence "
            "|T|>=5*mean/2>=8*mean/5 (strict last step for nonempty support)."
        ),
        "normalization_reduction": (
            "Deleting zero and duplicate columns preserves row identities, "
            "union closure, defect zero, and the frequency cap; Reimer then "
            "holds again by its theorem."
        ),
        "endpoint_consequence": (
            "A positive cap+Reimer+dominant-set defect floor, under either the "
            "displayed or normalized reporting convention, would prove the "
            "2/5 union-closed frequency theorem."
        ),
        "pair_defect_continuous_reimer_bound_refuted": True,
        "refuted_bound": (
            "mean_size >= log2(family_size)/2 - n*g(eps_vee) for every "
            "cap, active, separating family containing the full set, with "
            "g(t)->0 as t->0"
        ),
        "counterexample": (
            "REPORTED: Chase-Lovett arXiv:2211.11689v1 Example 1.4; "
            "the normalized Reimer deficit has positive limiting lower bound"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", type=str, default="")
    parser.add_argument("--size-limit", type=int, default=600)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    arguments = parser.parse_args()

    ceiling = cambie_ceiling()
    ceiling_low = ceiling["ceiling_lower"]
    ceiling_high = ceiling["ceiling_upper"]
    names = (
        [name.strip() for name in arguments.bases.split(",") if name.strip()]
        if arguments.bases
        else list(BASES)
    )

    print("Cambie two-atom upper obstruction, enclosed in exact rationals")
    print("  matching lower direction: OPEN / COMPUTATIONAL-EVIDENCE (math/uc audit)")
    print(f"  extremal atom b in [{decimal_lower(ceiling['atom_lower'], 16)}, "
          f"{decimal_upper(ceiling['atom_upper'], 16)}]")
    print(f"  atom weight   a in [{decimal_lower(ceiling['weight_lower'], 16)}, "
          f"{decimal_upper(ceiling['weight_upper'], 16)}]")
    print(f"  obstruction  c* in [{decimal_lower(ceiling_low, 16)}, "
          f"{decimal_upper(ceiling_high, 16)}]")
    print(f"  2/5 - c* >= {decimal_lower(Fraction(2, 5) - ceiling_high, 12)}")
    print(f"  c* > psi    : {ceiling_low > PSI_UPPER}")
    print(f"  repository alpha = {ALPHA} (Cambie's printed seven-decimal value)")

    entries: dict[str, object] = {}
    all_beyond_obstruction = True
    witnesses: list[str] = []
    witness_margins: dict[str, Fraction] = {}
    for name in names:
        base = BASES[name]
        family = tuple(sorted(base.reconstruct()))
        facts = barrier_facts(family, base.dimension)
        frequency = Fraction(max(degrees(family, base.dimension)), len(family))
        facts["above_cambie_obstruction"] = frequency > ceiling_high
        facts["within_cap_window_beyond_obstruction"] = (
            ceiling_high < frequency <= Fraction(2, 5)
        )
        all_beyond_obstruction &= bool(
            facts["within_cap_window_beyond_obstruction"]
        )
        if facts["chain_fails_at_every_order"]:
            witnesses.append(name)
            witness_margins[name] = facts["_log2_lower"] - facts["_maximum_upper"]
        print(
            f"\n{name}: n={facts['dimension']} m={facts['family_size']} "
            f"maxfreq={facts['maximum_frequency']} defect={facts['closure_defect']}"
        )
        print(
            f"  max_pi Q_pi <= {facts['order_maximum_upper']}   "
            f"log2 m >= {facts['log2_size_lower']}   "
            f"chain fails at every order: {facts['chain_fails_at_every_order']}"
        )
        print(
            f"  min_pi Q_pi >= {facts['order_minimum_lower']}   "
            f"barrier margin >= {facts['barrier_margin_lower']}"
        )
        entries[name] = {
            key: value for key, value in facts.items() if not key.startswith("_")
        }

    candidate_excluded = excluded_sizes(arguments.size_limit, T_CERT)
    reported_excluded = excluded_sizes(arguments.size_limit, ceiling_low)
    exclusions_coincide = candidate_excluded == reported_excluded
    reimer = chase_lovett_reimer_gap()
    endpoint = reimer_endpoint_facts()
    print(
        f"\nsizes excluded conditional on the math/uc candidate theorem "
        f"({len(candidate_excluded)} of them, all <= {max(candidate_excluded)}):"
    )
    print(f"  {list(candidate_excluded)}")
    print(
        "  same list from Cambie's reported threshold: "
        f"{exclusions_coincide}; multiples of five: "
        f"{[size for size in candidate_excluded if size % 5 == 0]}"
    )
    print("\nChase-Lovett Example 1.4 against Reimer (asymptotic leading rates)")
    print(
        f"  mean size rate    <= "
        f"{reimer['asymptotic_mean_size_rate_upper_decimal']} n"
    )
    print(f"  Reimer demands    >= {reimer['reimer_required_rate_lower_decimal']} n")
    print(f"  deficit           >= {reimer['reimer_deficit_rate_lower_decimal']} n")
    print(
        f"  Reimer satisfied  : "
        f"{reimer['reimer_satisfied_at_leading_order']}"
    )

    print("\nZero-defect endpoint")
    print("  Reimer automatic by the cited average-set-size theorem: true")
    print("  dominant set automatic from cap plus union closure: true")
    print("  defect-continuous approximate Reimer inequality: refuted")
    best_margin_fraction = max(witness_margins.values(), default=Fraction(0))
    best_margin = decimal_lower(best_margin_fraction, 10)
    verdict = (
        "CAP_ONLY_DEFECT_FLOOR_FALSE_REIMER_ESSENTIAL_AND_IID_CHAIN_HAS_EXACT_BARRIERS"
        if all_beyond_obstruction
        and witnesses
        and not reimer["reimer_satisfied_at_leading_order"]
        else "INCONCLUSIVE"
    )
    print(
        f"\norder-extremal iid-chain barrier witnesses "
        f"({len(witnesses)} of {len(names)}): {witnesses}"
    )
    print(f"  widest margin among them >= {best_margin}")
    report = {
        "cambie_upper_obstruction": {
            "source_status": (
                "HUMAN-AUDITED UPPER OBSTRUCTION; MATCHING LOWER DIRECTION "
                "OPEN / COMPUTATIONAL-EVIDENCE"
            ),
            "atom_lower_rational": str(ceiling["atom_lower"]),
            "atom_upper_rational": str(ceiling["atom_upper"]),
            "weight_lower_rational": str(ceiling["weight_lower"]),
            "weight_upper_rational": str(ceiling["weight_upper"]),
            "ceiling_lower_rational": str(ceiling_low),
            "ceiling_upper_rational": str(ceiling_high),
            "atom_lower_decimal": decimal_lower(ceiling["atom_lower"], 18),
            "atom_upper_decimal": decimal_upper(ceiling["atom_upper"], 18),
            "weight_lower_decimal": decimal_lower(ceiling["weight_lower"], 18),
            "weight_upper_decimal": decimal_upper(ceiling["weight_upper"], 18),
            "ceiling_lower_decimal": decimal_lower(ceiling_low, 18),
            "ceiling_upper_decimal": decimal_upper(ceiling_high, 18),
            "two_fifths_minus_ceiling_lower_rational": str(
                Fraction(2, 5) - ceiling_high
            ),
            "two_fifths_minus_ceiling_lower_decimal": decimal_lower(
                Fraction(2, 5) - ceiling_high, 12
            ),
            "strictly_above_psi": ceiling_low > PSI_UPPER,
            "repository_alpha": str(ALPHA),
            "alpha_provenance": (
                "INFERENCE: exact seven-decimal value printed by Cambie; "
                "formula match derived in repo"
            ),
        },
        "bases": entries,
        "every_base_beyond_cambie_upper_obstruction": all_beyond_obstruction,
        "barrier_witnesses": witnesses,
        "barrier_witness_count": len(witnesses),
        "widest_barrier_margin_lower_rational": str(best_margin_fraction),
        "widest_barrier_margin_lower_decimal": best_margin,
        "repo_candidate_threshold": {
            "status": (
                "CANDIDATE: machine-replayed interval domain plus "
                "human-audited bridge"
            ),
            "threshold_rational": str(T_CERT),
            "threshold_decimal": decimal_lower(T_CERT, 16),
            "excluded_sizes": list(candidate_excluded),
        },
        "reported_cambie_excluded_sizes": list(reported_excluded),
        "candidate_and_reported_exclusion_lists_coincide": exclusions_coincide,
        "excluded_size_limit": arguments.size_limit,
        "excluded_multiples_of_five": [
            size for size in candidate_excluded if size % 5 == 0
        ],
        "chase_lovett_example": reimer,
        "reimer_zero_defect_endpoint": endpoint,
        "verdict": verdict,
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"\n{verdict}")
    print(f"report: {arguments.report}")


if __name__ == "__main__":
    main()
