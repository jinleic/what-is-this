from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import compare_atlas_classifications as compare_module
from atlas import (
    _bounded_offset_screen,
    carry_trace,
    _carry_positions,
    _input_provenance,
    _near_pattern_key,
    _read_published,
    _record_pattern,
    _serialize_near_patterns,
    _serialize_shift_patterns,
    _shift_pattern_key,
    build_atlas,
)
from compare_atlas_classifications import (
    build_diff,
    carry_positions as independent_carry_positions,
    semantic_sha256,
)
from artifact_io import atomic_write_json
from overlap_density_bound import (
    band_integrand,
    finitary_threshold,
    pair_accounting,
    pair_obstruction,
    band_of,
    band_overlap_lower_bound,
    exponent_tail_bound,
    level_union_bound,
    overlap_family_member,
    overlap_lower_bound,
    CERTIFIED_BANDS,
)
from compensation_structure_analysis import (
    analyze as analyze_compensation_structure,
    blocker_eviction_crt,
    flexible_prefix_residue,
    flexible_prefix_residue_count,
    large_prime_window_obstructions,
    retained_prime_displacement,
    shift_first_bad_term_formula,
    short_cofactor_compensation_formula,
    term_compensation_formula,
)
from compensation_run_search import (
    search as search_compensation_runs,
    segmented_compensation_obstructions,
)
from crt_repair_cover import analyze as analyze_crt_repair_cover
from exact_eviction_certificate import (
    bad_window_positions,
    classify_eviction_candidate,
    surviving_blocker_offsets,
    translation_offset,
    window_large_primes,
)
from good_density_decomposition import (
    classify_term,
    asymptotic_union_bound,
    compensation_count,
    decompose_block,
    union_bound_terms,
)
from forcing_mass_bound import (
    counting_lemma_bound,
    exceeds_target_threshold,
    greedy_size_forcing_modulus,
    modulus_log_lower_bound,
    verify_counting_lemma,
    verify_mertens_constant,
    window_multiplicity_sum,
)
from parity_shift_certificate import (
    classify_shift,
    shift_window_difference,
    small_prime_slacks,
)
from smooth_density_analysis import (
    compensation_levels,
    count_narrow_repair_integers,
    count_short_cofactor_integers,
    measure_compensation_block,
)
from bad_window_near_miss_analysis import analyze as analyze_bad_window_near_misses
from elementary_class_eviction import (
    block_prime,
    covering_prime,
    covering_spans,
    union_measure,
    evict_progression,
    minimal_progression_length,
    progression_prime,
    evict_block,
    evict_initial,
    factor_modulus,
    least_prime_above,
    pinned_split,
    verify_eviction,
)
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    RatioSlackScanner,
    SPFFactorizer,
    TrialFactorizer,
    binomial_valuation,
    carry_count,
    first_witness,
    is_witness,
    is_witness_by_carries,
    is_witness_direct,
    local_floor_contribution,
    natural_shift_obstructions,
    primes_up_to,
    slack,
    witness_certificate,
    slack_floor_contributions,
    zone_contribution,
)
from moving_bad_window_search import (
    search as search_moving_bad_window,
    segmented_large_cofactors,
)
from prime_repair_analysis import analyze as analyze_prime_repairs
from residue_control_analysis import (
    analyze as analyze_residue_controls,
    band_good_residues,
    band_residue_slack,
    carry_free_residue_count,
    crt_coprime,
    fixed_class_single_level_obstruction,
    large_prime_bad_window_formula,
    local_safe_residues,
    prime_power_band,
    natural_shift_spike_formula,
)
from shift_repair_analysis import (
    analyze as analyze_shift_repairs,
    odd_shift_rough_modulus,
    zero_carry_modulus,
)
from tier_density_dp import (
    base_digits,
    carry_bound_terms,
    certified_upper_bound,
    decay_envelope,
    digit_sum_slack,
    envelope_applies,
    failure_count_enumerated,
    failure_count_power,
    failure_count_upto,
    high_digit_count,
)
from zero_carry_corridor_search import search as search_zero_carry_corridor
from zone_analysis import verify_zone_rectangle
from verify_known import input_provenance, read_rows


DATA = Path(__file__).resolve().parents[1] / "data" / "oeis_a375071.csv"
PROJECT = DATA.parent.parent
PRE_FIX_ATLAS = DATA.parent / "atlas_m1_50_k50000_pre_prime_key_fix.json"
ATLAS = DATA.parent / "atlas_m1_50_k50000.json"
SMOKE_PRE_FIX_ATLAS = DATA.parent / "atlas_smoke_m1_8_k500_pre_prime_key_fix.json"
SMOKE_ATLAS = DATA.parent / "atlas_smoke_m1_8_k500.json"
CLASSIFICATION_DIFF = DATA.parent / "atlas_prime_key_classification_diff.json"
CHANGED_CASES = DATA.parent / "atlas_prime_key_changed_cases.jsonl"
SHIFT_REPAIR_ANALYSIS = DATA.parent / "shift_repair_analysis_m1_50_k50000.json"
TIER_DENSITY = DATA.parent / "tier_density_m1_30.json"
FORCING_MASS = DATA.parent / "forcing_mass_bound.json"
OVERLAP_BOUND = DATA.parent / "overlap_density_bound.json"
ELEMENTARY_EVICTION = DATA.parent / "elementary_class_eviction.json"
SURVIVOR_CENSUS = DATA.parent / "survivor_census.json"
GOOD_DECOMPOSITION = DATA.parent / "good_density_decomposition.json"

PRIME_REPAIR_ANALYSIS = DATA.parent / "prime_repair_analysis_r64.json"
ZONE_CERTIFICATE = DATA.parent / "zone_theorem_m1_8_k2000.json"
ZERO_CORRIDOR = DATA.parent / "zero_carry_corridor_m5_18_t2000.json"
ZERO_CORRIDOR_EXTENDED = DATA.parent / "zero_carry_corridor_m8_12_t10000.json"
M27_P7_CORRIDOR = DATA.parent / "partial_zero_carry_m27_p7_t20000.json"
M27_P11_CORRIDOR = DATA.parent / "partial_zero_carry_m27_p11_t20000.json"
RESIDUE_CONTROL_ANALYSIS = (
    DATA.parent / "residue_control_analysis_m1_50_k50000.json"
)
MOVING_BAD_WINDOW_1M = DATA.parent / "moving_bad_window_m27_h1000000.json"
MOVING_BAD_WINDOW_20M = DATA.parent / "moving_bad_window_m27_h20000000.json"
MOVING_BAD_WINDOW_100M = DATA.parent / "moving_bad_window_m27_h100000000.json"
MOVING_BAD_WINDOW_100M_TO_200M = (
    DATA.parent / "moving_bad_window_m27_h100000000_199999999.json"
)
BAD_WINDOW_NEAR_MISS = DATA.parent / "bad_window_near_miss_m27_h1000000.json"
BAD_WINDOW_SURVIVORS = (
    DATA.parent / "bad_window_near_miss_m27_h100000000.json"
)
CRT_REPAIR_COVER = DATA.parent / "crt_repair_cover_m1_50_k50000_r64.json"
COMPENSATION_STRUCTURE = (
    DATA.parent / "compensation_structure_m1_20_k5000_m27_h200000000.json"
)
ADAPTIVE_EVICTION = DATA.parent / "adaptive_eviction_m27_h200000000.json"
COMPENSATION_RUN_FIRST_100M = (
    DATA.parent / "compensation_run_m27_k50001_h100000000.json"
)
COMPENSATION_RUN_NEXT_900M = (
    DATA.parent / "compensation_run_m27_k100050001_h900000000.json"
)
SHIFT_SCALE_RUN_SHARDS = tuple(
    DATA.parent / f"compensation_run_m27_k{start}_h200000000.json"
    for start in (
        5_049_091_644_619,
        5_049_291_644_619,
        5_049_491_644_619,
        5_049_691_644_619,
    )
)

def known_witnesses() -> list[tuple[int, int]]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        return [(int(row["m"]), int(row["k"])) for row in csv.DictReader(stream)]


class ExactArithmeticTests(unittest.TestCase):
    def test_kummer_matches_legendre_on_random_inputs(self) -> None:
        rng = random.Random(389)
        primes = primes_up_to(97)
        for _ in range(1_000):
            a = rng.randrange(0, 2_000)
            b = rng.randrange(0, 2_000)
            p = rng.choice(primes)
            self.assertEqual(
                binomial_valuation(a + b, a, p),
                carry_count(a, b, p),
                (a, b, p),
            )

    def test_all_carry_position_implementations_agree(self) -> None:
        rng = random.Random(389_002)
        primes = primes_up_to(97)
        for _ in range(1_000):
            a = rng.randrange(0, 20_000)
            b = rng.randrange(0, 20_000)
            p = rng.choice(primes)
            expected = tuple(carry_trace(a, b, p)["carry_positions"])
            self.assertEqual(_carry_positions(a, b, p), expected)
            self.assertEqual(independent_carry_positions(a, b, p), expected)

    def test_three_witness_checks_agree_on_small_grid(self) -> None:
        factorizer = TrialFactorizer(140)
        for m in range(13):
            for k in range(1, 121):
                direct = is_witness_direct(m, k, max_n=300)
                self.assertEqual(is_witness(m, k, factorizer), direct, (m, k))
                self.assertEqual(
                    is_witness_by_carries(m, k, factorizer), direct, (m, k)
                )

    def test_ratio_recurrence_matches_prime_valuations(self) -> None:
        factorizer = SPFFactorizer(500)
        for m in range(8):
            scanner = RatioSlackScanner(m, factorizer)
            for k in range(1, 181):
                expected = {
                    p: slack(m, k, p)
                    for p in primes_up_to(m + 2 * k)
                    if slack(m, k, p)
                }
                self.assertEqual(scanner.exponents, expected, (m, k))
                scanner.advance()

    def test_floor_contribution_sign_rule_and_slack_sum(self) -> None:
        rng = random.Random(389_003)
        primes = primes_up_to(97)
        for _ in range(2_000):
            m = rng.randrange(0, 100)
            k = rng.randrange(1, 5_000)
            p = rng.choice(primes)
            x = m + 2 * k
            contributions = slack_floor_contributions(m, k, p)
            self.assertEqual(sum(value for _, value in contributions), slack(m, k, p))
            for q, value in contributions:
                floor_parity = (x // q + m // q) % 2
                expected = (
                    0
                    if floor_parity == 0
                    else (1 if x % q + m % q < q else -1)
                )
                self.assertEqual(value, expected)
                self.assertEqual(value, local_floor_contribution(m, x, q))
                self.assertEqual(value, zone_contribution(m, k, q))
                self.assertIn(value, (-1, 0, 1))

    def test_invalid_arguments_are_rejected(self) -> None:
        factorizer = TrialFactorizer(20)
        with self.assertRaises(ValueError):
            is_witness(-1, 1, factorizer)
        with self.assertRaises(ValueError):
            is_witness(1, 0, factorizer)
        with self.assertRaises(ValueError):
            carry_count(1, 1, 1)
        with self.assertRaises(ValueError):
            TrialFactorizer(0)
        with self.assertRaises(ValueError):
            TrialFactorizer(MAX_TRIAL_FACTORIZER_VALUE + 1)


class ControlLemmaTests(unittest.TestCase):
    def test_prime_power_band_slack_is_an_exact_residue_function(self) -> None:
        for m, p, exponent in ((3, 5, 3), (3, 2, 4), (9, 5, 2)):
            lower, upper = prime_power_band(m, p, exponent)
            good = set(band_good_residues(m, p, exponent))
            modulus = p**exponent
            for k in range(lower, upper + 1):
                expected = band_residue_slack(m, p, exponent, k)
                self.assertEqual(slack(m, k, p), expected)
                self.assertEqual(expected >= 0, k % modulus in good)

    def test_local_safe_residues_extend_to_every_tested_high_part(self) -> None:
        for m in range(1, 13):
            for p in primes_up_to(m):
                modulus, residues = local_safe_residues(m, p)
                self.assertEqual(
                    sum(carry_count(r, m, p) == 0 for r in range(modulus)),
                    carry_free_residue_count(m, p),
                )
                for r in residues:
                    self.assertLess(r + m, modulus)
                    for high in range(8):
                        k = r + high * modulus
                        self.assertGreaterEqual(
                            carry_count(k, k + m, p),
                            carry_count(k, m, p),
                            (m, p, r, high),
                        )

    def test_large_prime_bad_window_formula_matches_exact_slack(self) -> None:
        factorizer = SPFFactorizer(500)
        for m in range(1, 9):
            for k in range(1, 121):
                seen: set[int] = set()
                for position in range(m // 2 + 1, m + 1):
                    for p in factorizer.factor(k + position):
                        if p <= m or p in seen:
                            continue
                        seen.add(p)
                        record = large_prime_bad_window_formula(m, k, p)
                        self.assertEqual(record["predicted_slack"], slack(m, k, p))

    def test_every_negative_prime_power_level_lies_in_the_bad_window(self) -> None:
        for m in range(1, 11):
            first_position = m // 2 + 1
            for k in range(1, 101):
                x = m + 2 * k
                for p in primes_up_to(x):
                    power = p
                    while power <= x:
                        residue_n = (m + k) % power
                        residue_m = m % power
                        contribution = zone_contribution(m, k, power)
                        self.assertEqual(
                            contribution,
                            (2 * residue_n - residue_m) // power,
                        )
                        if contribution == -1:
                            position = m - residue_n
                            self.assertGreaterEqual(position, first_position)
                            self.assertLessEqual(position, m)
                            self.assertEqual((k + position) % power, 0)
                        power *= p
                    if slack(m, k, p) < 0:
                        self.assertTrue(
                            any(
                                (k + position) % p == 0
                                for position in range(first_position, m + 1)
                            )
                        )

    def test_natural_shift_spike_formula_is_exact(self) -> None:
        factorizer = SPFFactorizer(500)
        for m in range(1, 13):
            for k in range(2, 121):
                for p in factorizer.factor(m + 2 * k):
                    if p <= m or p % 2 == 0:
                        continue
                    record = natural_shift_spike_formula(m, k, p)
                    self.assertEqual(
                        record["predicted_source_slack"], slack(m, k, p)
                    )
                    self.assertEqual(
                        record["predicted_shifted_slack"],
                        slack(m + 1, k - 1, p),
                    )

    def test_mirror_lift_binary_exception_and_coprime_crt(self) -> None:
        odd_record = large_prime_bad_window_formula(4, 25 * 24 - 3, 5)
        self.assertEqual(odd_record["predicted_slack"], 0)
        naive_binary = large_prime_bad_window_formula(1, 1, 2)
        repaired_binary = large_prime_bad_window_formula(1, 5, 2)
        self.assertEqual(naive_binary["predicted_slack"], -1)
        self.assertEqual(repaired_binary["predicted_slack"], 0)
        self.assertEqual(crt_coprime(((2, 3), (4, 5))), (14, 15))
        with self.assertRaises(ValueError):
            crt_coprime(((1, 4), (2, 6)))

    def test_every_sampled_fixed_crt_class_contains_an_exact_obstruction(self) -> None:
        cases = (
            ((2, 6, 1, 2, 101), 301),
            ((4, 30, 7, 3, 103), 1_027),
            ((9, 210, 13, 7, 107), 1_063),
            ((5, 12, 5, 3, 101), 401),
        )
        for arguments, expected_k in cases:
            record = fixed_class_single_level_obstruction(*arguments)
            self.assertEqual(record["k"], expected_k)
            self.assertEqual(record["k"] % record["modulus"], record["residue"])
            self.assertEqual(
                slack(record["m"], record["k"], record["prime"]), -1
            )

    def test_segmented_residuals_preserve_every_possible_fatal_prime(self) -> None:
        start = 100
        count = 101
        trial_primes = primes_up_to(14)
        residuals = segmented_large_cofactors(start, count, trial_primes)
        for n, actual in zip(range(start, start + count), residuals, strict=True):
            expected = n
            for p in trial_primes:
                while expected % p == 0:
                    expected //= p
            self.assertEqual(actual, expected)
            self.assertTrue(actual == 1 or actual > 14)

    def test_segmented_compensation_states_match_exact_term_formulas(self) -> None:
        m = 5
        start = 100
        count = 201
        primes, exponents, deficits = segmented_compensation_obstructions(
            m=m,
            start=start,
            count=count,
            trial_primes=primes_up_to(18),
        )
        factorizer = SPFFactorizer(start + count)
        for index, term in enumerate(range(start, start + count)):
            expected = sorted(
                (
                    p,
                    term_compensation_formula(term, p),
                )
                for p in factorizer.factor(term)
                if p > m
                and term_compensation_formula(term, p)["predicted_slack"] < 0
            )
            if expected:
                p, record = expected[0]
                self.assertEqual(primes[index], p)
                self.assertEqual(exponents[index], record["exponent"])
                self.assertEqual(
                    deficits[index],
                    -record["predicted_slack"],
                )
            else:
                self.assertEqual(primes[index], 0)
                self.assertEqual(exponents[index], 0)
                self.assertEqual(deficits[index], 0)

    def test_compensation_run_search_recovers_a_published_minimum(self) -> None:
        result = search_compensation_runs(
            target_m=3,
            start_k=1,
            offset_limit=207,
            chunk_size=67,
        )
        self.assertEqual(
            [record["k"] for record in result["survivors"]["witness_records"]],
            [207],
        )
        self.assertEqual(
            result["counts"]["large_prime_rejected_offsets"]
            + result["counts"]["large_prime_good_offsets"],
            207,
        )
        self.assertTrue(
            all(
                obstruction["prime"] <= 3
                for record in result["survivors"]["reported_records"]
                for obstruction in record["obstructions"]
            )
        )


    def test_large_prime_control_separates_into_bad_window_terms(self) -> None:
        factorizer = SPFFactorizer(500)
        for m in range(1, 9):
            for k in range(1, 121):
                certificate = witness_certificate(m, k, factorizer)
                local = {
                    record["prime"]: record["predicted_slack"]
                    for record in large_prime_window_obstructions(
                        m, k, factorizer
                    )
                }
                exact = {
                    p: value
                    for p, value in certificate.obstructions.items()
                    if p > m
                }
                self.assertEqual(local, exact)
                small_failure = any(
                    p <= m for p in certificate.obstructions
                )
                self.assertEqual(
                    certificate.is_witness,
                    not local and not small_failure,
                )

    def test_flexible_prefix_repair_cone_has_exact_count(self) -> None:
        for p in (3, 5, 7):
            for exponent in (1, 2, 3):
                modulus = p**exponent
                residues = [
                    residue
                    for residue in range(modulus)
                    if flexible_prefix_residue(p, exponent, residue)
                ]
                self.assertEqual(
                    len(residues),
                    flexible_prefix_residue_count(p, exponent),
                )
                for residue in residues:
                    term = p**exponent * residue
                    record = term_compensation_formula(term, p)
                    self.assertEqual(record["exponent"], exponent)
                    self.assertGreaterEqual(record["predicted_slack"], 0)
        self.assertFalse(flexible_prefix_residue(5, 1, 2))
        self.assertTrue(flexible_prefix_residue(5, 1, 3))
        with self.assertRaises(ValueError):
            flexible_prefix_residue_count(2, 1)

    def test_even_shift_spike_is_first_bad_term_compensation(self) -> None:
        factorizer = SPFFactorizer(500)
        for m in range(2, 13, 2):
            for k in range(2, 121):
                for p in factorizer.factor(m + 2 * k):
                    if p <= m + 1:
                        continue
                    record = shift_first_bad_term_formula(m, k, p)
                    self.assertEqual(
                        record["predicted_shifted_slack"],
                        slack(m + 1, k - 1, p),
                    )
                    self.assertEqual(
                        record["predicted_shifted_slack"],
                        record["first_bad_term_slack"],
                    )

    def test_retained_prime_repair_requires_product_displacement(self) -> None:
        self.assertEqual(slack(5, 18, 7), -1)
        self.assertEqual(slack(5, 18, 11), -1)
        record = retained_prime_displacement(18, ((7, 3), (11, 4)))
        self.assertEqual(record["prime_product"], 77)
        self.assertEqual(record["least_forward_delta"], 77)
        self.assertEqual(record["least_forward_candidate"], 95)
        for p, position in ((7, 3), (11, 4)):
            self.assertEqual((record["least_forward_candidate"] + position) % p, 0)

    def test_short_cofactor_spike_formula_is_exact_and_sharp(self) -> None:
        for p in (2, 3, 5, 7, 11):
            for exponent in (1, 2, 3):
                for cofactor in range(1, p):
                    term = p**exponent * cofactor
                    record = short_cofactor_compensation_formula(term, p)
                    indicator = int(2 * cofactor > p)
                    self.assertEqual(record["positive_prefix_levels"], indicator)
                    self.assertEqual(
                        record["predicted_slack"],
                        indicator - exponent,
                    )
                    self.assertEqual(
                        record["predicted_slack"],
                        term_compensation_formula(term, p)["predicted_slack"],
                    )
        self.assertEqual(
            short_cofactor_compensation_formula(15, 5)["predicted_slack"],
            0,
        )
        self.assertEqual(
            short_cofactor_compensation_formula(75, 5)["predicted_slack"],
            -1,
        )
        self.assertEqual(
            short_cofactor_compensation_formula(8, 2)["predicted_slack"],
            -3,
        )
        with self.assertRaises(ValueError):
            short_cofactor_compensation_formula(35, 5)

    def test_odd_prime_minimal_compensating_cofactor_is_sharp(self) -> None:
        for p in (3, 5, 7, 11):
            for exponent in (1, 2, 3):
                power = p**exponent
                minimum = (power + 1) // 2
                boundary = term_compensation_formula(power * minimum, p)
                self.assertEqual(boundary["exponent"], exponent)
                self.assertEqual(
                    boundary["positive_prefix_levels"],
                    exponent,
                )
                self.assertEqual(boundary["predicted_slack"], 0)
                for cofactor in range(1, minimum):
                    if cofactor % p:
                        self.assertLess(
                            term_compensation_formula(power * cofactor, p)[
                                "predicted_slack"
                            ],
                            0,
                        )

        for p in (5, 7):
            for exponent in (1, 2, 3):
                power = p**exponent
                threshold_x = power * (power + 1)
                boundary_k = (threshold_x - 2) // 2
                below_k = (power * (power - 1) - 2) // 2
                self.assertEqual(
                    shift_first_bad_term_formula(
                        2, boundary_k, p
                    )["predicted_shifted_slack"],
                    0,
                )
                self.assertLess(
                    shift_first_bad_term_formula(
                        2, below_k, p
                    )["predicted_shifted_slack"],
                    0,
                )

    def test_zero_carry_crt_evicts_every_prescribed_blocker(self) -> None:
        controls = ((7, 3), (11, 4))
        record = blocker_eviction_crt(5, 18, controls)
        self.assertEqual(record["bad_window_length"], 3)
        self.assertEqual(record["small_prime_zero_carry_modulus"], 1_800)
        self.assertEqual(record["blocker_prime_product"], 77)
        self.assertEqual(record["combined_modulus"], 138_600)
        self.assertEqual(record["combined_eviction_class_count"], 32)

        candidate = record["least_forward_candidate"]
        self.assertGreater(candidate, 18)
        for p in primes_up_to(5):
            self.assertEqual(carry_count(candidate, 5, p), 0)
            self.assertGreaterEqual(slack(5, candidate, p), 0)
        for p, position in controls:
            self.assertEqual((18 + position) % p, 0)
            self.assertEqual(candidate % p, (18 + 3) % p)
            self.assertTrue(
                all((candidate + new_position) % p for new_position in range(3, 6))
            )

        exact_class_count = sum(
            residue % 1_800 == 0
            and all(
                (residue + position) % p
                for p, _ in controls
                for position in range(3, 6)
            )
            for residue in range(record["combined_modulus"])
        )
        self.assertEqual(
            exact_class_count,
            record["combined_eviction_class_count"],
        )

    def test_even_window_eviction_uses_length_not_first_position(self) -> None:
        controls = ((5, 3), (11, 4))
        record = blocker_eviction_crt(4, 7, controls)
        self.assertEqual(record["bad_window_first_position"], 3)
        self.assertEqual(record["bad_window_length"], 2)
        self.assertEqual(record["small_prime_zero_carry_modulus"], 72)
        self.assertEqual(record["blocker_prime_product"], 55)
        self.assertEqual(record["combined_modulus"], 3_960)
        self.assertEqual(record["combined_eviction_class_count"], 27)

        candidate = record["least_forward_candidate"]
        for p, position in controls:
            self.assertEqual((7 + position) % p, 0)
            self.assertEqual(candidate % p, (7 + 2) % p)
            self.assertTrue(
                all((candidate + new_position) % p for new_position in (3, 4))
            )
        exact_class_count = sum(
            residue % 72 == 0
            and all(
                (residue + position) % p
                for p, _ in controls
                for position in (3, 4)
            )
            for residue in range(record["combined_modulus"])
        )
        self.assertEqual(exact_class_count, 27)

    def test_compensation_structure_rejects_mismatched_input_metadata(self) -> None:
        def invoke(
            *,
            moving_path: Path = MOVING_BAD_WINDOW_100M,
            near_path: Path = BAD_WINDOW_SURVIVORS,
            extended_path: Path | None = None,
        ) -> None:
            analyze_compensation_structure(
                max_m=1,
                max_k=1,
                shift_max_m=2,
                shift_max_k=2,
                eviction_max_m=1,
                eviction_max_k=1,
                eviction_candidate_max_x=100,
                moving_path=moving_path,
                near_miss_path=near_path,
                extended_moving_path=extended_path,
            )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            mismatched_near = json.loads(
                BAD_WINDOW_SURVIVORS.read_text(encoding="utf-8")
            )
            mismatched_near["parameters"]["target_m"] = 26
            mismatched_near_path = root / "mismatched_near.json"
            mismatched_near_path.write_text(
                json.dumps(mismatched_near),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "near-miss parameters"):
                invoke(near_path=mismatched_near_path)

            other_moving = json.loads(
                MOVING_BAD_WINDOW_100M.read_text(encoding="utf-8")
            )
            other_near = json.loads(
                BAD_WINDOW_SURVIVORS.read_text(encoding="utf-8")
            )
            other_moving["parameters"]["target_m"] = 26
            other_near["parameters"]["target_m"] = 26
            other_moving_path = root / "other_m_moving.json"
            other_near_path = root / "other_m_near.json"
            other_moving_path.write_text(
                json.dumps(other_moving),
                encoding="utf-8",
            )
            other_near_path.write_text(
                json.dumps(other_near),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "survivor sections require m=27"):
                invoke(
                    moving_path=other_moving_path,
                    near_path=other_near_path,
                )

            for source_path, argument, message in (
                (MOVING_BAD_WINDOW_100M, "moving", "moving schema"),
                (BAD_WINDOW_SURVIVORS, "near", "near-miss schema"),
                (
                    MOVING_BAD_WINDOW_100M_TO_200M,
                    "extended",
                    "extended moving schema",
                ),
            ):
                payload = json.loads(source_path.read_text(encoding="utf-8"))
                payload["schema_version"] = 2
                invalid_path = root / f"{argument}_schema.json"
                invalid_path.write_text(json.dumps(payload), encoding="utf-8")
                with self.subTest(argument=argument):
                    with self.assertRaisesRegex(ValueError, message):
                        if argument == "moving":
                            invoke(moving_path=invalid_path)
                        elif argument == "near":
                            invoke(near_path=invalid_path)
                        else:
                            invoke(extended_path=invalid_path)

    def test_local_formula_boundaries_and_multilevel_compensation(self) -> None:
        self.assertEqual(term_compensation_formula(5, 5)["predicted_slack"], -1)
        self.assertEqual(slack(5, 2, 5), 0)
        factorizer = SPFFactorizer(100)
        self.assertNotIn(
            5,
            {
                record["prime"]
                for record in large_prime_window_obstructions(
                    5, 2, factorizer
                )
            },
        )
        self.assertEqual(
            term_compensation_formula(75, 5)["predicted_slack"],
            -1,
        )
        self.assertEqual(slack(2, 73, 5), -1)
        self.assertEqual(
            term_compensation_formula(600, 5)["predicted_slack"],
            0,
        )
        self.assertEqual(slack(2, 598, 5), 0)

    def test_shift_first_term_identity_excludes_target_small_prime(self) -> None:
        self.assertEqual(slack(5, 2, 5), 0)
        with self.assertRaises(ValueError):
            shift_first_bad_term_formula(4, 3, 5)



class PublishedWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = known_witnesses()
        cls.factorizer = TrialFactorizer(max(m + 2 * k for m, k in cls.rows))

    def test_all_oeis_rows_are_exact_witnesses(self) -> None:
        for m, k in self.rows:
            certificate = witness_certificate(m, k, self.factorizer)
            self.assertTrue(certificate.is_witness, (m, k, certificate.obstructions))
            self.assertEqual(
                is_witness_by_carries(m, k, self.factorizer), True, (m, k)
            )
            for p, left in certificate.left_valuations.items():
                self.assertEqual(left, carry_count(m, k, p), (m, k, p))
                self.assertEqual(
                    certificate.right_valuations[p], carry_count(k, m + k, p)
                )

    def test_each_nontrivial_published_witness_is_tight(self) -> None:
        for m, k in self.rows:
            if m == 0:
                continue
            certificate = witness_certificate(m, k, self.factorizer)
            self.assertGreater(len(certificate.tight_primes), 0, (m, k))
            self.assertEqual(min(certificate.slacks.values()), 0, (m, k))

    def test_bounded_scanner_reproduces_two_published_minima(self) -> None:
        factorizer = SPFFactorizer(5_010)
        self.assertIsNone(first_witness(3, 206, factorizer))
        self.assertEqual(first_witness(3, 207, factorizer), 207)
        self.assertIsNone(first_witness(5, 2_474, factorizer))
        self.assertEqual(first_witness(5, 2_475, factorizer), 2_475)


class NaturalShiftTests(unittest.TestCase):
    def test_shift_slack_identity(self) -> None:
        for m in range(1, 15):
            for k in range(2, 80):
                for p in primes_up_to(m + 2 * k):
                    expected = (
                        slack(m, k, p)
                        + binomial_valuation(m + 1, 1, p)
                        - binomial_valuation(m + 2 * k, 1, p)
                    )
                    self.assertEqual(slack(m + 1, k - 1, p), expected)

    def test_known_odd_shift_counterexample_fails_only_at_three(self) -> None:
        factorizer = TrialFactorizer(1_977)
        self.assertTrue(is_witness(3, 987, factorizer))
        obstructions = natural_shift_obstructions(3, 987, factorizer)
        self.assertEqual([item.prime for item in obstructions], [3])
        self.assertEqual(obstructions[0].shifted_slack, -1)
        self.assertFalse(is_witness(4, 986, factorizer))

    def test_published_odd_m_shifts_have_no_large_prime_obstruction(self) -> None:
        rows = known_witnesses()
        factorizer = TrialFactorizer(max(m + 2 * k for m, k in rows))
        for m, k in rows:
            if m % 2 == 1 and k >= 2:
                obstructions = natural_shift_obstructions(m, k, factorizer)
                self.assertTrue(all(item.prime <= m for item in obstructions), (m, k))


class AtlasContractTests(unittest.TestCase):
    def test_small_atlas_reproduces_first_six_minima_and_boundaries(self) -> None:
        atlas = build_atlas(
            min_m=1,
            max_m=6,
            max_k=2_475,
            repair_radius=4,
            representative_limit=1,
            top_primes=8,
            pattern_limit=16,
            known_path=DATA,
        )
        self.assertEqual(atlas["status"], "EXACT_BOUNDED_COMPUTATION")
        self.assertEqual(
            [row["first_witness"] for row in atlas["m_rows"]],
            [5, 4, 207, 206, 2_475, 984],
        )
        self.assertTrue(
            all(row["published_minimum_reproduced"] for row in atlas["m_rows"])
        )
        self.assertEqual(
            atlas["one_prime_near_witness_patterns"]["total_instances"],
            sum(row["one_prime_near_witness_count"] for row in atlas["m_rows"]),
        )
        self.assertEqual(
            atlas["natural_shift_summary"]["odd_m_large_prime_violations"], 0
        )
        self.assertIn(
            "not CRT-defined repairs",
            atlas["bounded_offset_repair_screen"]["interpretation"],
        )
        expected_digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
        self.assertEqual(atlas["input"]["logical_name"], "data/oeis_a375071.csv")
        self.assertEqual(atlas["input"]["sha256"], expected_digest)

    def test_pattern_classification_never_merges_distinct_prime_bases(self) -> None:
        m, k = 4, 1
        signatures = {
            p: (_carry_positions(k, m, p), _carry_positions(k, m + k, p))
            for p in (2, 3)
        }
        self.assertEqual(signatures[2], signatures[3])

        near_groups: dict[tuple[object, ...], dict[str, object]] = {}
        shift_groups: dict[tuple[object, ...], dict[str, object]] = {}
        for p in (2, 3):
            left_positions, right_positions = signatures[p]
            _record_pattern(
                near_groups,
                _near_pattern_key(p, 1, left_positions, right_positions),
                (m, k, p),
                1,
            )
            _record_pattern(
                shift_groups,
                _shift_pattern_key(
                    p, 0, 0, 1, -1, left_positions, right_positions
                ),
                (m, k, p),
                1,
            )

        near = _serialize_near_patterns(near_groups, 10)
        shifted = _serialize_shift_patterns(shift_groups, 10)
        for summary in (near, shifted):
            self.assertEqual(summary["distinct_patterns"], 2)
            self.assertEqual([row["prime"] for row in summary["patterns"]], [2, 3])
            for row in summary["patterns"]:
                self.assertTrue(
                    all(
                        representative["prime"] == row["prime"]
                        for representative in row["representatives"]
                    )
                )

    def test_dominant_pre_fix_shift_collision_splits_by_prime(self) -> None:
        pre_fix = json.loads(PRE_FIX_ATLAS.read_text(encoding="utf-8"))
        dominant = pre_fix["natural_shift_summary"]["failure_patterns"]["patterns"][0]
        old_pattern = {
            "source_slack": 0,
            "v_p_m_plus_one": 0,
            "v_p_m_plus_two_k": 1,
            "shifted_slack": -1,
            "left_carry_positions": [0],
            "right_carry_positions": [0],
        }
        self.assertEqual(dominant["count"], 6_053)
        self.assertEqual(
            {field: dominant[field] for field in old_pattern},
            old_pattern,
        )
        self.assertEqual(dominant["representatives"][0]["prime"], 5)

        diff = json.loads(CLASSIFICATION_DIFF.read_text(encoding="utf-8"))
        matching_splits = [
            group
            for group in diff["classifications"]["shift"]["split_groups"]
            if group["old_pattern"] == old_pattern
        ]
        self.assertEqual(len(matching_splits), 1)
        split = matching_splits[0]
        self.assertEqual(split["old_count"], 6_053)
        self.assertEqual(split["prime_count"], 1_860)
        self.assertEqual(
            sum(child["count"] for child in split["new_patterns"]),
            split["old_count"],
        )
        child_primes = {
            child["pattern"]["prime"] for child in split["new_patterns"]
        }
        self.assertEqual(len(child_primes), split["prime_count"])
        self.assertEqual(
            len({child["new_pattern_id"] for child in split["new_patterns"]}),
            split["prime_count"],
        )
        for child in split["new_patterns"]:
            self.assertEqual(
                {
                    field: child["pattern"][field]
                    for field in old_pattern
                },
                old_pattern,
            )

    def test_custom_input_provenance_is_neutral_and_content_addressed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            custom = Path(directory) / "candidate.csv"
            custom.write_bytes(DATA.read_bytes())
            snapshot = custom.read_bytes()
            custom.write_text("changed after snapshot", encoding="utf-8")
            for provenance in (
                _input_provenance(custom, snapshot),
                input_provenance(custom, snapshot),
            ):
                self.assertEqual(provenance["logical_name"], "candidate.csv")
                self.assertEqual(
                    provenance["sha256"], hashlib.sha256(snapshot).hexdigest()
                )
                self.assertIn("Caller-supplied", provenance["source"])
            self.assertEqual(len(_read_published(snapshot)), 27)
            self.assertEqual(len(read_rows(snapshot)), 27)

    def test_offset_screen_accounts_for_boundaries_overlap_and_no_repair(self) -> None:
        witness_bits: list[bytearray | None] = [None, bytearray(16), bytearray(16)]
        source = witness_bits[1]
        target = witness_bits[2]
        assert source is not None and target is not None
        for k in (2, 4, 7, 13, 15):
            source[k] = 1
        for k in (4, 7):
            target[k] = 1

        screen = _bounded_offset_screen(witness_bits, 1, 2, 15, 2)
        self.assertEqual(screen["boundary_excluded_failures"], 2)
        self.assertEqual(screen["eligible_natural_shift_failures"], 3)
        self.assertEqual(
            {item["delta"]: item["covered"] for item in screen["coverage_by_offset"]},
            {-2: 1, -1: 0, 0: 0, 1: 2, 2: 0},
        )
        self.assertEqual(
            screen["greedy_offset_menu"], [{"delta": 1, "newly_covered": 2}]
        )
        self.assertEqual(screen["covered_by_greedy_menu"], 2)
        self.assertEqual(screen["uncovered_by_every_tested_offset"], 1)
        self.assertEqual(
            sum(item["newly_covered"] for item in screen["greedy_offset_menu"]),
            screen["covered_by_greedy_menu"],
        )
        self.assertEqual(
            screen["covered_by_greedy_menu"]
            + screen["uncovered_by_every_tested_offset"],
            screen["eligible_natural_shift_failures"],
        )

    def test_small_prime_moduli_match_their_exact_carry_contracts(self) -> None:
        self.assertEqual(zero_carry_modulus(3, primes_up_to(3)), 36)
        self.assertEqual(zero_carry_modulus(5, primes_up_to(5)), 1_800)
        self.assertEqual(odd_shift_rough_modulus(7, primes_up_to(7)), (105, 48))
        for m in range(1, 11):
            primes = primes_up_to(m)
            modulus = zero_carry_modulus(m, primes)
            for p in primes:
                self.assertEqual(carry_count(modulus, m, p), 0)

    def test_small_certificate_builders_execute_behaviorally(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atlas_path = root / "atlas.json"
            small_atlas = build_atlas(
                min_m=1,
                max_m=4,
                max_k=300,
                repair_radius=4,
                representative_limit=1,
                top_primes=8,
                pattern_limit=16,
                known_path=DATA,
            )
            atomic_write_json(atlas_path, small_atlas)
            shift_result = analyze_shift_repairs(
                min_m=1, max_m=4, max_k=300, atlas_path=atlas_path
            )
            prime_result = analyze_prime_repairs(
                min_m=1,
                max_m=4,
                max_k=300,
                radius=4,
                atlas_path=atlas_path,
            )
            zone_result = verify_zone_rectangle(3, 50)
            corridor_result = search_zero_carry_corridor(3, 4, 20)
            control_result = analyze_residue_controls(
                min_m=1,
                max_m=4,
                max_k=300,
                atlas_path=atlas_path,
                formula_max_m=2,
                formula_max_k=20,
            )
            crt_cover_result = analyze_crt_repair_cover(
                min_m=1,
                max_m=4,
                max_k=300,
                radius=4,
                atlas_path=atlas_path,
            )
            moving_result = search_moving_bad_window(
                target_m=4,
                start_k=206,
                offset_limit=20,
                chunk_size=7,
                known_path=DATA,
            )
            near_miss_result = analyze_bad_window_near_misses(
                target_m=4,
                start_k=206,
                offset_limit=20,
                chunk_size=7,
                representative_limit=10,
                known_path=DATA,
            )

            self.assertEqual(shift_result["status"], "EXACT_BOUNDED_ANALYSIS")
            self.assertEqual(
                prime_result["status"], "EXACT_BOUNDED_PRIME_REPAIR_ANALYSIS"
            )
            self.assertEqual(
                zone_result["status"], "PASS_EXACT_FINITE_CERTIFICATE"
            )
            self.assertEqual(
                corridor_result["status"], "EXACT_BOUNDED_CANDIDATE_FAMILY"
            )
            self.assertEqual(
                control_result["status"], "PASS_EXACT_BOUNDED_CONTROL_ANALYSIS"
            )
            self.assertEqual(
                crt_cover_result["status"], "EXACT_BOUNDED_CRT_REPAIR_COVER"
            )
            self.assertEqual(
                moving_result["status"], "EXACT_BOUNDED_MOVING_BAD_WINDOW_SEARCH"
            )
            self.assertEqual(
                near_miss_result["status"],
                "EXACT_BOUNDED_BAD_WINDOW_NEAR_MISS_ANALYSIS",
            )
            self.assertEqual(
                moving_result["counts"]["single_level_rejected_offsets"]
                + moving_result["counts"]["single_level_sieve_survivors"],
                20,
            )
            self.assertEqual(zone_result["counts"]["pairs"], 150)
            self.assertEqual(
                sum(row["candidate_count"] for row in corridor_result["rows"]), 40
            )

    def test_diff_rejects_tampering_and_pins_loaded_snapshots(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pre = root / "pre.json"
            post = root / "post.json"
            cases = root / "cases.jsonl"
            pre.write_bytes(SMOKE_PRE_FIX_ATLAS.read_bytes())
            original_post = SMOKE_ATLAS.read_bytes()
            post.write_bytes(original_post)

            result = build_diff(pre, post, cases, 1, 8, 500)
            self.assertEqual(result["status"], "PASS_EXACT_RECOMPUTATION")
            with self.assertRaises(ValueError):
                build_diff(pre, post, pre, 1, 8, 500)

            tampered = json.loads(original_post)
            tampered["one_prime_near_witness_patterns"]["patterns"][0][
                "prime"
            ] += 2
            atomic_write_json(post, tampered)
            with self.assertRaises(AssertionError):
                build_diff(pre, post, root / "tampered-cases.jsonl", 1, 8, 500)

            post.write_bytes(original_post)
            real_scan = compare_module.scan_events

            def replace_after_snapshot(
                min_m: int, max_m: int, max_k: int
            ) -> dict[str, list[tuple[object, ...]]]:
                post.write_text("{}", encoding="utf-8")
                return real_scan(min_m, max_m, max_k)

            with patch.object(
                compare_module, "scan_events", side_effect=replace_after_snapshot
            ):
                snapshotted = build_diff(
                    pre, post, root / "snapshot-cases.jsonl", 1, 8, 500
                )
            self.assertEqual(
                snapshotted["artifacts"]["repaired"]["sha256"],
                hashlib.sha256(original_post).hexdigest(),
            )

    def test_checked_in_verification_artifact_pins_input_bytes(self) -> None:
        artifact = json.loads(
            (DATA.parent / "known_witness_verification.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact["schema_version"], 1)
        self.assertEqual(
            artifact["input"]["sha256"], hashlib.sha256(DATA.read_bytes()).hexdigest()
        )
        atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
        self.assertEqual(atlas["schema_version"], 2)
        self.assertEqual(
            atlas["input"]["sha256"], hashlib.sha256(DATA.read_bytes()).hexdigest()
        )
        self.assertEqual(
            atlas["implementation_sha256"]["atlas.py"],
            hashlib.sha256((PROJECT / "atlas.py").read_bytes()).hexdigest(),
        )
        for summary in (
            atlas["one_prime_near_witness_patterns"],
            atlas["natural_shift_summary"]["failure_patterns"],
        ):
            self.assertEqual(
                summary["reported_fraction"],
                {
                    "numerator": summary["reported_instances"],
                    "denominator": summary["total_instances"],
                },
            )
        near_patterns = atlas["one_prime_near_witness_patterns"]
        shift_patterns = atlas["natural_shift_summary"]["failure_patterns"]
        self.assertEqual(near_patterns["key_fields"][0], "prime")
        self.assertEqual(shift_patterns["key_fields"][0], "prime")
        for summary in (near_patterns, shift_patterns):
            for pattern in summary["patterns"]:
                self.assertTrue(
                    all(
                        representative["prime"] == pattern["prime"]
                        for representative in pattern["representatives"]
                    )
                )
        diff = json.loads(CLASSIFICATION_DIFF.read_text(encoding="utf-8"))
        self.assertEqual(diff["status"], "PASS_EXACT_RECOMPUTATION")
        self.assertEqual(
            diff["artifacts"]["pre_fix"]["sha256"],
            hashlib.sha256(PRE_FIX_ATLAS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            diff["artifacts"]["repaired"]["sha256"],
            hashlib.sha256(ATLAS.read_bytes()).hexdigest(),
        )
        cases_bytes = CHANGED_CASES.read_bytes()
        self.assertEqual(
            diff["artifacts"]["changed_cases"]["sha256"],
            hashlib.sha256(cases_bytes).hexdigest(),
        )
        self.assertEqual(
            cases_bytes.count(b"\n"),
            diff["artifacts"]["changed_cases"]["case_count"] + 1,
        )
        self.assertEqual(
            diff["classifications"]["near"]["new_distinct_patterns"],
            near_patterns["distinct_patterns"],
        )
        self.assertEqual(
            diff["classifications"]["shift"]["new_distinct_patterns"],
            shift_patterns["distinct_patterns"],
        )
        self.assertIn(
            "regime_event_counts", diff["classifications"]["near"]
        )
        self.assertIn(
            "extremal_patterns", diff["classifications"]["shift"]
        )
        shift_analysis = json.loads(SHIFT_REPAIR_ANALYSIS.read_text(encoding="utf-8"))
        self.assertEqual(shift_analysis["status"], "EXACT_BOUNDED_ANALYSIS")
        self.assertEqual(
            shift_analysis["input"]["atlas_sha256"],
            hashlib.sha256(ATLAS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            shift_analysis["finite_checks"],
            {
                "fatal_prime_power_checks": 4_578,
                "odd_crt_safe_shift_checks": 5_798,
                "zero_carry_prime_checks": 16_826,
            },
        )
        for row in shift_analysis["odd_shift_crt_analysis"]["rows"]:
            self.assertEqual(
                row["crt_safe_source_witnesses"],
                row["crt_safe_successful_shifts"],
            )
        failed_odd_rows = [
            row
            for row in shift_analysis["odd_shift_crt_analysis"]["rows"]
            if row["failed_shifts"]
        ]
        self.assertEqual(
            sum(
                row["failed_sources_repaired_by_safe_witness_within_bound"]
                for row in failed_odd_rows
            ),
            21,
        )
        self.assertEqual(
            max(
                item["periods"]
                for row in failed_odd_rows
                for item in row["safe_witness_repair_period_counts"]
            ),
            559,
        )
        prime_repair = json.loads(PRIME_REPAIR_ANALYSIS.read_text(encoding="utf-8"))
        self.assertEqual(
            prime_repair["status"], "EXACT_BOUNDED_PRIME_REPAIR_ANALYSIS"
        )
        self.assertEqual(
            prime_repair["input"]["atlas_sha256"],
            hashlib.sha256(ATLAS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            prime_repair["counts"],
            {
                "boundary_excluded": 12,
                "exact_original_prime_alignment": 492,
                "failed_natural_shifts": 5_690,
                "new_prime_interference_delay": 3_224,
                "original_primes_repaired_but_no_witness": 1_962,
            },
        )
        for row in shift_analysis["even_shift_spike_analysis"]["rows"]:
            self.assertLessEqual(
                row["failed_shifts_explained_by_fatal_spike"],
                row["failed_shifts"],
            )
        zone = json.loads(ZONE_CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(zone["status"], "PASS_EXACT_FINITE_CERTIFICATE")
        self.assertEqual(zone["counts"]["pairs"], 16_000)
        self.assertEqual(zone["counts"]["level_checks"], 330_023)
        self.assertEqual(zone["counts"]["single_level_obstruction_checks"], 25_536)

        corridor = json.loads(ZERO_CORRIDOR.read_text(encoding="utf-8"))
        extended = json.loads(ZERO_CORRIDOR_EXTENDED.read_text(encoding="utf-8"))
        m27_p7 = json.loads(M27_P7_CORRIDOR.read_text(encoding="utf-8"))
        m27_p11 = json.loads(M27_P11_CORRIDOR.read_text(encoding="utf-8"))
        self.assertEqual(corridor["finite_checks"]["candidate_count"], 24_014)
        self.assertEqual(corridor["finite_checks"]["witness_count"], 5)
        self.assertEqual(extended["finite_checks"]["witness_count"], 0)
        self.assertEqual(m27_p7["finite_checks"]["witness_count"], 0)
        self.assertEqual(m27_p11["finite_checks"]["witness_count"], 0)
        self.assertEqual(
            m27_p7["rows"][0]["failure_candidate_categories"][
                "single_level_large_prime"
            ],
            20_000,
        )
        self.assertEqual(
            m27_p11["rows"][0]["failure_candidate_categories"][
                "single_level_large_prime"
            ],
            20_000,
        )

        control = json.loads(RESIDUE_CONTROL_ANALYSIS.read_text(encoding="utf-8"))
        self.assertEqual(control["status"], "PASS_EXACT_BOUNDED_CONTROL_ANALYSIS")
        self.assertEqual(
            control["input"]["atlas_sha256"],
            hashlib.sha256(ATLAS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            control["implementation_sha256"]["residue_control_analysis.py"],
            hashlib.sha256(
                (PROJECT / "residue_control_analysis.py").read_bytes()
            ).hexdigest(),
        )
        self.assertEqual(
            control["general_large_prime_compensation"]["counts"]["formula_checks"],
            50_360,
        )
        self.assertEqual(
            control["small_prime_local_safe_crt"]["extension_checks"], 523_408
        )
        self.assertEqual(control["prime_power_band_crt"]["periodicity_checks"], 711)
        self.assertEqual(control["fixed_class_obstruction_checks"]["exact_checks"], 4)
        self.assertEqual(
            control["natural_shift_spike_analysis"]["totals"][
                "failed_with_exact_large_spike"
            ],
            5_620,
        )
        self.assertEqual(
            control["natural_shift_spike_analysis"]["totals"][
                "additional_failures_from_multilevel_formula"
            ],
            1_042,
        )
        structure = json.loads(
            COMPENSATION_STRUCTURE.read_text(encoding="utf-8")
        )
        self.assertEqual(
            structure["status"],
            "PASS_EXACT_COMPENSATION_STRUCTURE_ANALYSIS",
        )
        self.assertEqual(structure["schema_version"], 3)
        self.assertEqual(
            structure["input"]["moving_window"]["sha256"],
            hashlib.sha256(MOVING_BAD_WINDOW_100M.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            structure["input"]["near_miss"]["sha256"],
            hashlib.sha256(BAD_WINDOW_SURVIVORS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            structure["input"]["extended_moving_window"]["sha256"],
            hashlib.sha256(
                MOVING_BAD_WINDOW_100M_TO_200M.read_bytes()
            ).hexdigest(),
        )
        for input_name in (
            "moving_window",
            "near_miss",
            "extended_moving_window",
        ):
            self.assertEqual(structure["input"][input_name]["schema_version"], 1)
        self.assertEqual(
            structure["parameters"]["contiguous_moving_offset_count"],
            200_000_000,
        )
        self.assertEqual(
            structure["implementation_sha256"][
                "compensation_structure_analysis.py"
            ],
            hashlib.sha256(
                (PROJECT / "compensation_structure_analysis.py").read_bytes()
            ).hexdigest(),
        )
        for dependency in ("erdos389.py", "residue_control_analysis.py"):
            self.assertEqual(
                structure["implementation_sha256"][dependency],
                hashlib.sha256((PROJECT / dependency).read_bytes()).hexdigest(),
            )
        self.assertEqual(
            structure["small_and_large_prime_separation"]["counts"],
            {
                "large_good_but_small_prime_failed": 190,
                "large_prime_factor_checks": 587_169,
                "large_prime_good_windows": 1_294,
                "large_prime_obstructions": 442_704,
                "pairs": 100_000,
                "short_cofactor_checks": 405_780,
                "short_cofactor_obstructions": 363_792,
                "short_cofactor_prime_power_spikes": 7_173,
                "witnesses": 1_104,
            },
        )
        self.assertEqual(
            structure["even_shift_first_bad_term_identity"]["counts"],
            {
                "prime_checks": 19_654,
                "short_cofactor_checks": 14_510,
                "short_cofactor_prime_power_spikes": 314,
                "short_cofactor_spikes": 12_969,
                "spikes": 15_336,
            },
        )
        self.assertEqual(
            structure["zero_carry_blocker_eviction"]["counts"],
            {
                "canonical_candidate_witnesses": 103,
                "canonical_candidates_beyond_exact_check_bound": 1_109,
                "canonical_candidates_exactly_checked": 2_737,
                "canonical_candidates_with_new_prime_interference": 2_634,
                "new_large_obstruction_prime_checks": 5_570,
                "new_short_cofactor_obstructions": 4_061,
                "new_short_cofactor_prime_power_spikes": 190,
                "old_blocker_prime_checks": 7_782,
                "source_systems": 3_846,
            },
        )
        self.assertEqual(
            structure["zero_carry_blocker_eviction"][
                "maximum_canonical_candidate_decimal_digits"
            ],
            15,
        )
        retained = structure["m27_survivor_retained_prime_bounds"]
        self.assertEqual(
            retained["counts"][
                "all_offsets_with_certified_large_prime_obstruction"
            ],
            100_000_000,
        )
        self.assertEqual(
            [
                row["retained_prime_product_decimal_digits"]
                for row in retained["rows"]
            ],
            [69, 30, 35],
        )
        self.assertEqual(
            [
                row["zero_carry_eviction_candidate_decimal_digits"]
                for row in retained["rows"]
            ],
            [87, 49, 53],
        )
        self.assertTrue(
            all(
                row["zero_carry_eviction_candidate_exceeds_factorizer_bound"]
                for row in retained["rows"]
            )
        )
        extended_retained = structure[
            "m27_extended_survivor_retained_prime_bounds"
        ]
        self.assertEqual(
            extended_retained["counts"],
            {
                "all_offsets_with_certified_large_prime_obstruction": 100_000_000,
                "obstruction_prime_checks": 113,
                "offsets_with_certified_large_prime_obstruction": 100_000_000,
                "single_level_rejected_offsets": 99_999_991,
                "survivors": 9,
                "survivors_with_large_prime_obstructions": 9,
                "survivors_with_small_prime_obstructions": 2,
                "zero_carry_eviction_crt_systems": 9,
            },
        )
        self.assertEqual(
            [
                row["retained_prime_product_decimal_digits"]
                for row in extended_retained["rows"]
            ],
            [36, 74, 79, 45, 79, 78, 42, 48, 48],
        )
        self.assertEqual(
            [
                row["zero_carry_eviction_candidate_decimal_digits"]
                for row in extended_retained["rows"]
            ],
            [55, 94, 98, 64, 98, 97, 61, 67, 67],
        )
        self.assertTrue(
            all(
                row["zero_carry_eviction_candidate_exceeds_factorizer_bound"]
                for row in extended_retained["rows"]
            )
        )
        adaptive = json.loads(ADAPTIVE_EVICTION.read_text(encoding="utf-8"))
        self.assertEqual(
            adaptive["status"],
            "EXACT_BOUNDED_ADAPTIVE_EVICTION_SEARCH",
        )
        self.assertEqual(adaptive["schema_version"], 1)
        self.assertEqual(
            adaptive["input"]["structure_sha256"],
            hashlib.sha256(COMPENSATION_STRUCTURE.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            adaptive["implementation_sha256"]["adaptive_eviction_analysis.py"],
            hashlib.sha256(
                (PROJECT / "adaptive_eviction_analysis.py").read_bytes()
            ).hexdigest(),
        )
        for dependency in ("erdos389.py", "residue_control_analysis.py"):
            self.assertEqual(
                adaptive["implementation_sha256"][dependency],
                hashlib.sha256((PROJECT / dependency).read_bytes()).hexdigest(),
            )
        self.assertEqual(
            adaptive["counts"],
            {
                "candidates_rejected_by_trial_prime": 1,
                "candidates_unclassified_after_trial_scan": 11,
                "old_blocker_prime_checks": 140,
                "old_blocker_window_divisibility_checks": 1_960,
                "safe_candidates_beyond_factorizer_bound": 12,
                "small_prime_slack_checks": 108,
                "source_systems": 12,
                "systems_with_safe_t_in_bound": 12,
                "tested_t_values_through_first_safe": 2_369,
                "trial_prime_divisibility_checks": 12_649_506,
                "trial_prime_divisor_hits": 224,
            },
        )
        self.assertEqual(
            adaptive["uniform_small_prime_system"]["safe_t_count_per_full_period"],
            adaptive["uniform_small_prime_system"]["safe_class_count"],
        )
        self.assertEqual(
            [record["first_safe_t"] for record in adaptive["records"]],
            [706, 63, 10, 70, 119, 423, 153, 22, 24, 409, 268, 90],
        )
        self.assertEqual(
            [
                record["candidate_decimal_digits"]
                for record in adaptive["records"]
            ],
            [72, 32, 36, 37, 77, 81, 47, 80, 80, 45, 51, 50],
        )
        trial_obstructions = [
            record["trial_prime_obstruction"]
            for record in adaptive["records"]
            if record["trial_prime_obstruction"] is not None
        ]
        self.assertEqual(len(trial_obstructions), 1)
        self.assertEqual(
            {
                field: trial_obstructions[0][field]
                for field in ("prime", "position", "predicted_slack")
            },
            {"prime": 137_341, "position": 21, "predicted_slack": -1},
        )
        self.assertTrue(
            all(
                record["candidate_exceeds_factorizer_bound"]
                and all(
                    row["slack"] >= 0 for row in record["small_prime_slacks"]
                )
                for record in adaptive["records"]
            )
        )
        compensation_runs = [
            json.loads(COMPENSATION_RUN_FIRST_100M.read_text(encoding="utf-8")),
            json.loads(COMPENSATION_RUN_NEXT_900M.read_text(encoding="utf-8")),
        ]
        expected_run_counts = (
            (100_000_000, 86_629_164, 13_370_849, 8),
            (900_000_000, 783_699_380, 116_300_633, 9),
        )
        for run, (offsets, bad_terms, good_terms, maximum_run) in zip(
            compensation_runs,
            expected_run_counts,
            strict=True,
        ):
            self.assertEqual(
                run["status"],
                "EXACT_BOUNDED_COMPENSATION_RUN_SEARCH",
            )
            self.assertEqual(
                run["implementation_sha256"]["compensation_run_search.py"],
                hashlib.sha256(
                    (PROJECT / "compensation_run_search.py").read_bytes()
                ).hexdigest(),
            )
            self.assertEqual(
                run["implementation_sha256"][
                    "compensation_structure_analysis.py"
                ],
                hashlib.sha256(
                    (PROJECT / "compensation_structure_analysis.py").read_bytes()
                ).hexdigest(),
            )
            for dependency in ("erdos389.py", "residue_control_analysis.py"):
                self.assertEqual(
                    run["implementation_sha256"][dependency],
                    hashlib.sha256(
                        (PROJECT / dependency).read_bytes()
                    ).hexdigest(),
                )
            self.assertEqual(run["parameters"]["offset_count"], offsets)
            self.assertEqual(
                run["counts"]["large_prime_rejected_offsets"],
                offsets,
            )
            self.assertEqual(run["counts"].get("large_prime_good_offsets", 0), 0)
            self.assertEqual(run["counts"]["compensation_bad_terms"], bad_terms)
            self.assertEqual(run["counts"]["compensation_good_terms"], good_terms)
            self.assertEqual(
                run["counts"]["maximum_consecutive_compensation_good_terms"],
                maximum_run,
            )
            self.assertEqual(
                bad_terms + good_terms,
                offsets + run["parameters"]["bad_window_length"] - 1,
            )
            self.assertEqual(run["survivors"]["reported_records"], [])
            self.assertEqual(run["survivors"]["witness_records"], [])
        m27_atlas_row = next(
            row for row in atlas["m_rows"] if row["m"] == 27
        )
        self.assertEqual(m27_atlas_row["witness_count"], 0)
        self.assertEqual(
            compensation_runs[0]["parameters"]["start_k"],
            atlas["parameters"]["max_k"] + 1,
        )
        self.assertEqual(
            compensation_runs[1]["parameters"]["start_k"],
            compensation_runs[0]["parameters"]["end_k"] + 1,
        )
        self.assertEqual(
            compensation_runs[1]["parameters"]["end_k"],
            1_000_050_000,
        )



        crt_cover = json.loads(CRT_REPAIR_COVER.read_text(encoding="utf-8"))
        self.assertEqual(crt_cover["status"], "EXACT_BOUNDED_CRT_REPAIR_COVER")
        self.assertEqual(
            crt_cover["input"]["atlas_sha256"],
            hashlib.sha256(ATLAS.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            crt_cover["category_counts"],
            {
                "first_small_prime_tier_is_witness": 14,
                "window_blocked_before_later_witness": 3_702,
                "window_blocked_no_witness_within_radius": 1_962,
            },
        )
        self.assertEqual(
            crt_cover["window_blocker_counts"],
            {"compensation_deficit": 1_847, "single_level": 4_759},
        )
        self.assertEqual(crt_cover["counts"]["band_residue_checks"], 905_367)
        self.assertEqual(
            crt_cover["implementation_sha256"]["crt_repair_cover.py"],
            hashlib.sha256((PROJECT / "crt_repair_cover.py").read_bytes()).hexdigest(),
        )

        for moving_path, expected_survivors in (
            (MOVING_BAD_WINDOW_1M, 0),
            (MOVING_BAD_WINDOW_20M, 0),
            (MOVING_BAD_WINDOW_100M, 3),
            (MOVING_BAD_WINDOW_100M_TO_200M, 9),
        ):
            moving = json.loads(moving_path.read_text(encoding="utf-8"))
            self.assertEqual(
                moving["status"], "EXACT_BOUNDED_MOVING_BAD_WINDOW_SEARCH"
            )
            self.assertEqual(
                moving["implementation_sha256"]["moving_bad_window_search.py"],
                hashlib.sha256(
                    (PROJECT / "moving_bad_window_search.py").read_bytes()
                ).hexdigest(),
            )
            survivor_count = moving["counts"].get(
                "single_level_sieve_survivors", 0
            )
            self.assertEqual(
                moving["counts"]["single_level_rejected_offsets"] + survivor_count,
                moving["parameters"]["offset_count"],
            )
            histogram_counts = {
                row["count"]
                for row in moving["single_level_rejections"][
                    "fatal_obstructions_per_window"
                ]
            }
            self.assertEqual(survivor_count, expected_survivors)
            if expected_survivors:
                self.assertIn(0, histogram_counts)
                self.assertEqual(
                    moving["counts"]["survivor_nonwitnesses"],
                    expected_survivors,
                )
            else:
                self.assertNotIn(0, histogram_counts)
            self.assertEqual(
                moving["starting_point"]["published_source"][
                    "natural_shift_obstructions"
                ][0]["prime"],
                5_048_891_644_633,
            )

        near_miss = json.loads(BAD_WINDOW_NEAR_MISS.read_text(encoding="utf-8"))
        self.assertEqual(
            near_miss["status"], "EXACT_BOUNDED_BAD_WINDOW_NEAR_MISS_ANALYSIS"
        )
        self.assertEqual(near_miss["minimum_single_level_obstruction_count"], 1)
        self.assertEqual(near_miss["minimum_window_count"], 5)
        self.assertEqual(
            near_miss["counts"]["mirror_candidates_exactly_checked"], 9
        )
        self.assertEqual(
            near_miss["counts"]["mirror_candidates_with_new_prime_interference"],
            9,
        )
        self.assertEqual(
            near_miss["implementation_sha256"][
                "bad_window_near_miss_analysis.py"
            ],
            hashlib.sha256(
                (PROJECT / "bad_window_near_miss_analysis.py").read_bytes()
            ).hexdigest(),
        )

        survivors = json.loads(BAD_WINDOW_SURVIVORS.read_text(encoding="utf-8"))
        self.assertEqual(
            survivors["status"], "EXACT_BOUNDED_BAD_WINDOW_NEAR_MISS_ANALYSIS"
        )
        self.assertEqual(survivors["minimum_single_level_obstruction_count"], 0)
        self.assertEqual(survivors["minimum_window_count"], 3)
        self.assertEqual(
            survivors["counts"],
            {
                "combined_large_prime_controls_checked": 27,
                "combined_small_prime_controls_checked": 27,
                "mirror_candidates_beyond_factorizer_bound": 0,
                "mirror_candidates_exactly_checked": 27,
                "mirror_candidates_with_new_prime_interference": 27,
                "representative_obstruction_primes": 27,
                "representatives_exactly_checked": 3,
            },
        )
        self.assertTrue(
            all(
                not system["within_factorizer_bound"]
                for system in survivors["combined_crt_systems"]
            )
        )

        corridor_witnesses = [
            (row["m"], witness["k"], witness["natural_shift_is_witness"])
            for row in corridor["rows"]
            for witness in row["witness_representatives"]
        ]
        factorizer = TrialFactorizer(
            max(m + 2 * k for m, k, _ in corridor_witnesses)
        )
        for m, k, shift_is_witness in corridor_witnesses:
            self.assertTrue(is_witness(m, k, factorizer))
            self.assertEqual(
                not natural_shift_obstructions(m, k, factorizer),
                shift_is_witness,
            )
        offset = atlas["bounded_offset_repair_screen"]
        self.assertEqual(
            sum(item["newly_covered"] for item in offset["greedy_offset_menu"]),
            offset["covered_by_greedy_menu"],
        )
        self.assertEqual(
            offset["covered_by_greedy_menu"]
            + offset["uncovered_by_every_tested_offset"],
            offset["eligible_natural_shift_failures"],
        )


class HalfWindowEvictionTests(unittest.TestCase):
    def test_half_window_translation_evicts_every_attached_large_prime(self) -> None:
        max_m, max_k = 14, 400
        factorizer = SPFFactorizer(max_m + max_k + translation_offset(max_m))
        checked = 0
        for m in range(1, max_m + 1):
            offset = translation_offset(m)
            self.assertEqual(offset, -(-m // 2))
            self.assertEqual(len(bad_window_positions(m)), offset)
            for k in range(1, max_k + 1):
                attached = window_large_primes(m, k, factorizer)
                for p, source_position in attached.items():
                    self.assertEqual((k + source_position) % p, 0)
                    for position in bad_window_positions(m):
                        self.assertNotEqual((k + offset + position) % p, 0)
                        checked += 1
        self.assertGreater(checked, 50_000)

    def test_every_smaller_offset_keeps_an_explicit_blocker(self) -> None:
        for m in range(1, 60):
            offset = translation_offset(m)
            first_position = m // 2 + 1
            p = next(q for q in primes_up_to(2 * m + 2) if q > m)
            for delta in range(1, offset):
                source_position = first_position + delta
                self.assertIn(source_position, bad_window_positions(m))
                k = p - source_position
                self.assertGreaterEqual(k, 1)
                self.assertEqual((k + source_position) % p, 0)
                survivor = source_position - delta
                self.assertEqual(survivor, first_position)
                self.assertEqual((k + delta + survivor) % p, 0)
                self.assertIn(delta, surviving_blocker_offsets(m, source_position))

    def test_translated_survivor_candidates_stay_small_and_are_classified(
        self,
    ) -> None:
        structure = json.loads(COMPENSATION_STRUCTURE.read_text(encoding="utf-8"))
        row = structure["m27_survivor_retained_prime_bounds"]["rows"][1]
        k = row["k"]
        candidate = k + translation_offset(27)
        self.assertLessEqual(27 + 2 * candidate, MAX_TRIAL_FACTORIZER_VALUE)
        factorizer = TrialFactorizer(27 + 2 * candidate)
        classification = classify_eviction_candidate(
            27, k, factorizer, tuple(primes_up_to(27))
        )
        self.assertEqual(classification["candidate"], candidate)
        self.assertEqual(classification["candidate_decimal_digits"], 13)
        self.assertTrue(classification["small_prime_tier_ok"])
        self.assertEqual(classification["category"], "new_large_prime_only")
        self.assertTrue(classification["new_large_prime_obstructions"])
        for record in row["obstructions"]:
            for position in bad_window_positions(27):
                self.assertNotEqual(
                    (candidate + position) % record["prime"], 0
                )


class PowersmoothDensityTests(unittest.TestCase):
    def test_compensated_prime_powers_obey_the_minimal_term_bound(self) -> None:
        factorizer = SPFFactorizer(4_000)
        for m in range(1, 13):
            for k in range(1, 300):
                for position in range(m // 2 + 1, m + 1):
                    term = k + position
                    for p, exponent in factorizer.factor(term).items():
                        if p <= m:
                            continue
                        power = p**exponent
                        good = (
                            compensation_levels(term, p, exponent) >= exponent
                        )
                        if good:
                            self.assertLessEqual(power * (power + 1), 2 * term)
                        if power * power > 2 * term:
                            self.assertFalse(good)

    def test_dominant_prime_repairs_only_inside_the_narrow_window(self) -> None:
        factorizer = SPFFactorizer(6_000)
        checked = 0
        for term in range(2, 6_001):
            for p, exponent in factorizer.factor(term).items():
                cofactor = term // p**exponent
                if cofactor >= p:
                    continue
                good = compensation_levels(term, p, exponent) >= exponent
                self.assertEqual(
                    good, exponent == 1 and term < p * p < 2 * term, (term, p)
                )
                checked += 1
        self.assertGreater(checked, 1_000)

    def test_narrow_repair_counts_agree_and_decay(self) -> None:
        small = count_narrow_repair_integers(10_000, 10_000)
        self.assertEqual(small["direct_recount"], small["narrow_repair_integers"])
        larger = count_narrow_repair_integers(100_000, 0)
        self.assertIsNone(larger["direct_recount"])
        self.assertLess(larger["density"], small["density"])

    def test_minimal_compensated_term_is_attained_and_sharp(self) -> None:
        for p in (3, 5, 7, 11):
            for exponent in (1, 2):
                power = p**exponent
                term = power * (power + 1) // 2
                self.assertEqual(
                    compensation_levels(term, p, exponent), exponent
                )
                self.assertEqual(term_compensation_formula(term, p)["exponent"],
                                 exponent)
                smaller = term - power
                if smaller > 0 and smaller % power == 0:
                    self.assertLess(
                        compensation_levels(smaller, p, exponent), exponent
                    )

    def test_short_cofactor_count_agrees_with_largest_prime_factor(self) -> None:
        limit = 20_000
        record = count_short_cofactor_integers(limit)
        self.assertEqual(record["limit"], limit)
        factorizer = SPFFactorizer(limit)
        expected = sum(
            1
            for n in range(2, limit + 1)
            if any(p * p > 2 * n for p in factorizer.factor(n))
        )
        self.assertEqual(record["short_cofactor_integers"], expected)
        self.assertAlmostEqual(record["density"], expected / limit)
        self.assertLess(record["density"], 0.6931471805599453)

    def test_published_witness_window_is_a_thirteen_term_good_run(self) -> None:
        published = dict(known_witnesses())
        source_k = published[25]
        start = source_k + 25 // 2 + 1
        block = measure_compensation_block(
            target_m=27, start=start, count=13, label="unit_control"
        )
        self.assertEqual(block["compensation_good_terms"], 13)
        self.assertEqual(block["longest_run"], 13)
        shift_target_term = published[26] + 26 // 2
        factorizer = TrialFactorizer(shift_target_term)
        self.assertEqual(factorizer.factor(shift_target_term), {shift_target_term: 1})
        self.assertEqual(
            term_compensation_formula(shift_target_term, shift_target_term)[
                "predicted_slack"
            ],
            -1,
        )


class ParityShiftReductionTests(unittest.TestCase):
    def test_window_identity_matches_source_parity(self) -> None:
        for m in range(1, 30):
            for k in (2, 3, 17, 500):
                adjoined = shift_window_difference(m, k)
                if m % 2:
                    self.assertEqual(adjoined, ())
                else:
                    self.assertEqual(adjoined, (k + m // 2,))
                    self.assertEqual(2 * adjoined[0], m + 2 * k)

    def test_odd_source_shifts_never_fail_at_a_large_prime(self) -> None:
        factorizer = SPFFactorizer(2_200)
        checked = 0
        for m in range(1, 12, 2):
            for k in range(2, 1_500):
                if large_prime_window_obstructions(m, k, factorizer):
                    continue
                if any(row["slack"] < 0 for row in small_prime_slacks(m, k)):
                    continue
                row = classify_shift(m, k, factorizer)
                self.assertTrue(row["large_prime_tier_ok"])
                self.assertEqual(row["adjoined_terms"], [])
                self.assertEqual(
                    row["shift_is_witness"], row["small_prime_tier_ok"]
                )
                checked += 1
        self.assertGreater(checked, 100)

    def test_published_even_shifts_fail_at_the_adjoined_term(self) -> None:
        published = dict(known_witnesses())
        m, k = 26, published[26]
        adjoined = shift_window_difference(m, k)
        self.assertEqual(adjoined, (5_048_891_644_633,))
        factorizer = TrialFactorizer(27 + 2 * (k + m))
        row = classify_shift(m, k, factorizer)
        self.assertFalse(row["shift_is_witness"])
        self.assertTrue(row["small_prime_tier_ok"])
        self.assertEqual(
            [record["prime"] for record in row["adjoined_obstructions"]],
            [5_048_891_644_633],
        )
        self.assertEqual(
            row["adjoined_obstructions"][0]["predicted_slack"], -1
        )


class ShiftScaleRunShardTests(unittest.TestCase):
    def test_shift_scale_shards_cover_one_contiguous_billion(self) -> None:
        moving = [
            json.loads(MOVING_BAD_WINDOW_100M.read_text(encoding="utf-8")),
            json.loads(
                MOVING_BAD_WINDOW_100M_TO_200M.read_text(encoding="utf-8")
            ),
        ]
        intervals = [
            (
                payload["parameters"]["start_k"],
                payload["parameters"]["start_k"]
                + payload["parameters"]["offset_count"],
            )
            for payload in moving
        ]
        expected_counts = (
            (24_829_873, 175_170_140, 8),
            (24_661_439, 175_338_574, 9),
            (25_053_346, 174_946_667, 9),
            (24_604_542, 175_395_471, 8),
        )
        for path, (good, bad, longest) in zip(
            SHIFT_SCALE_RUN_SHARDS, expected_counts, strict=True
        ):
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["status"], "EXACT_BOUNDED_COMPENSATION_RUN_SEARCH"
            )
            self.assertEqual(payload["parameters"]["target_m"], 27)
            offsets = payload["parameters"]["offset_count"]
            self.assertEqual(offsets, 200_000_000)
            self.assertEqual(
                payload["counts"]["large_prime_rejected_offsets"], offsets
            )
            self.assertEqual(
                payload["counts"].get("large_prime_good_offsets", 0), 0
            )
            self.assertEqual(payload["counts"]["compensation_good_terms"], good)
            self.assertEqual(payload["counts"]["compensation_bad_terms"], bad)
            self.assertEqual(
                payload["counts"][
                    "maximum_consecutive_compensation_good_terms"
                ],
                longest,
            )
            self.assertLess(longest, 14)
            self.assertEqual(payload["survivors"]["witness_records"], [])
            self.assertEqual(
                payload["implementation_sha256"]["compensation_run_search.py"],
                hashlib.sha256(
                    (PROJECT / "compensation_run_search.py").read_bytes()
                ).hexdigest(),
            )
            start = payload["parameters"]["start_k"]
            intervals.append((start, start + offsets))

        intervals.sort()
        for (_, end), (next_start, _) in zip(intervals, intervals[1:]):
            self.assertEqual(end, next_start)
        self.assertEqual(intervals[0][0], 5_048_891_644_619)
        self.assertEqual(intervals[-1][1], 5_049_891_644_619)
        self.assertEqual(intervals[-1][1] - intervals[0][0], 1_000_000_000)


class TierDensityTests(unittest.TestCase):
    def test_digit_sum_identity_matches_slack(self) -> None:
        rng = random.Random(28)
        for _ in range(400):
            m = rng.randrange(0, 60)
            k = rng.randrange(1, 10**6)
            p = rng.choice(primes_up_to(59))
            self.assertEqual(digit_sum_slack(m, k, p), slack(m, k, p))

    def test_carry_bound_holds_and_is_attained(self) -> None:
        attained = False
        for m in range(1, 14):
            for p in primes_up_to(m):
                for k in range(1, 300):
                    exact = slack(m, k, p)
                    d_count, run, window = carry_bound_terms(m, k, p)
                    margin = exact - (window - d_count - run)
                    self.assertGreaterEqual(margin, 0)
                    attained = attained or margin == 0
        self.assertTrue(attained)

    def test_high_digit_count_matches_definition(self) -> None:
        for p in primes_up_to(50):
            self.assertEqual(
                high_digit_count(p), sum(1 for d in range(p) if 2 * d >= p)
            )

    def test_digit_program_matches_enumeration(self) -> None:
        for m in (0, 1, 5, 12, 27):
            for p in primes_up_to(min(max(m, 2), 11)):
                for digits in (3, 4):
                    if p**digits <= 2 * m or p**digits > 60_000:
                        continue
                    expected = failure_count_enumerated(m, p, p**digits - 1)
                    self.assertEqual(failure_count_power(m, p, digits), expected)
                    self.assertEqual(
                        failure_count_upto(m, p, p**digits - 1), expected
                    )

    def test_prefix_program_matches_enumeration_off_power_boundaries(self) -> None:
        for m in (1, 27):
            for p in (2, 3, 5):
                for limit in (1, 2, 999, 5_000, 12_345):
                    self.assertEqual(
                        failure_count_upto(m, p, limit),
                        failure_count_enumerated(m, p, limit),
                    )

    def test_closed_form_bound_dominates_exact_count(self) -> None:
        for m in (0, 1, 7, 27):
            for p in primes_up_to(min(max(m, 2), 13)):
                for digits in range(len(base_digits(m, p)) + 1, 13):
                    self.assertGreaterEqual(
                        certified_upper_bound(m, p, digits),
                        failure_count_power(m, p, digits),
                    )

    def test_decay_envelope_dominates_bound_where_hypothesis_holds(self) -> None:
        checked = 0
        for m, p, digits in ((1, 2, 32), (1, 2, 40), (2, 3, 33), (27, 2, 96)):
            self.assertTrue(envelope_applies(m, p, digits))
            density = certified_upper_bound(m, p, digits) / p**digits
            self.assertLessEqual(density, decay_envelope(digits))
            checked += 1
        self.assertEqual(checked, 4)

    def test_tier_artifact_pins_inputs_and_records_exact_agreement(self) -> None:
        payload = json.loads(TIER_DENSITY.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "PASS_EXACT_TIER_DENSITY")
        self.assertEqual(
            payload["implementation_sha256"]["tier_density_dp.py"],
            hashlib.sha256((PROJECT / "tier_density_dp.py").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            payload["input"]["density_artifact_sha256"],
            hashlib.sha256(
                (PROJECT / "data" / "smooth_density_m1_20_k3000.json").read_bytes()
            ).hexdigest(),
        )
        self.assertEqual(
            payload["identity_and_inequality"]["minimum_inequality_margin"], 0
        )
        self.assertGreater(payload["envelope_checked_cells"], 0)
        for row in payload["prime_decay"]:
            densities = [cell["density"] for cell in row["cells"]]
            self.assertLess(densities[-1], densities[0])
            for cell in row["cells"]:
                self.assertLessEqual(cell["density"], cell["certified_bound_density"])

    def test_recorded_tier_failure_density_is_small_at_witness_scale(self) -> None:
        payload = json.loads(TIER_DENSITY.read_text(encoding="utf-8"))
        scale = next(
            row for row in payload["scales"] if row["limit"] == 5_048_891_644_621
        )
        self.assertLess(scale["union_failure_density"], 0.02)
        self.assertGreater(scale["surviving_density_lower_bound"], 0.98)
        recomputed = sum(
            failure_count_upto(payload["parameters"]["target_m"], row["p"], scale["limit"])
            for row in scale["per_prime"]
        )
        self.assertEqual(recomputed, scale["union_failure_count"])


class SemanticDigestTests(unittest.TestCase):
    def test_semantic_digest_ignores_exactly_the_environment_blocks(self) -> None:
        payload = json.loads(SMOKE_ATLAS.read_text(encoding="utf-8"))
        baseline = semantic_sha256(payload)
        for field, value in (
            ("generated_at_utc", "1999-12-31T23:59:59Z"),
            ("resources", {"wall_seconds": 12345.0}),
            ("runtime", {"python": "0.0.0"}),
        ):
            mutated = dict(payload)
            mutated[field] = value
            self.assertEqual(semantic_sha256(mutated), baseline)
        substantive = dict(payload)
        substantive["schema_version"] = payload["schema_version"] + 1
        self.assertNotEqual(semantic_sha256(substantive), baseline)

    def test_checked_in_atlas_semantic_digest_is_pinned(self) -> None:
        payload = json.loads(ATLAS.read_text(encoding="utf-8"))
        self.assertEqual(
            semantic_sha256(payload),
            "c1d5f70d134a1544e86f000cf1c9c45df265d5778e04fef6b44cad343589505f",
        )
        diff = json.loads(CLASSIFICATION_DIFF.read_text(encoding="utf-8"))
        self.assertEqual(
            diff["artifacts"]["repaired"]["semantic_sha256"], semantic_sha256(payload)
        )


class ForcingMassTests(unittest.TestCase):
    def test_counting_lemma_holds_for_every_residue(self) -> None:
        report = verify_counting_lemma(8, (2, 3, 5), 3)
        self.assertGreater(report["checked_residues"], 500)
        self.assertGreaterEqual(report["tightest_case"]["slack"], 0)

    def test_counting_lemma_is_nearly_attained(self) -> None:
        self.assertEqual(window_multiplicity_sum(1, 11, 4, 11**4 - 1), 4)
        self.assertAlmostEqual(counting_lemma_bound(1, 11, 4), 4.1)

    def test_mertens_correction_bound_holds(self) -> None:
        report = verify_mertens_constant(50_000)
        self.assertGreater(report["tightest"]["margin"], 0)

    def test_greedy_size_forcing_modulus_meets_hypothesis_and_dominates_bound(self) -> None:
        factorizer = TrialFactorizer(10**14)
        for m, k in ((5, 2475), (13, 7_979_077), (20, 1_019_547_824)):
            greedy = greedy_size_forcing_modulus(m, k, factorizer)
            bound = modulus_log_lower_bound(m, k)
            self.assertGreater(greedy["log_modulus"], bound)
            for row in greedy["terms"]:
                term = row["term"]
                self.assertEqual(term, k + row["position"])
                forced = 10 ** (row["forced_divisor_decimal_digits"] - 1)
                self.assertLess(term / forced, 20 * (2 * term) ** 0.5)

    def test_threshold_exists_only_for_windows_of_five_or_more(self) -> None:
        for m in range(1, 9):
            self.assertIsNone(exceeds_target_threshold(m))
        for m in (9, 13, 27):
            threshold = exceeds_target_threshold(m)
            self.assertIsNotNone(threshold)
            self.assertGreater(modulus_log_lower_bound(m, threshold), math.log(threshold))
            self.assertLessEqual(
                modulus_log_lower_bound(m, threshold - 1), math.log(threshold - 1)
            )

    def test_forcing_artifact_pins_inputs_and_brackets_every_small_witness(self) -> None:
        payload = json.loads(FORCING_MASS.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "PASS_EXACT_FORCING_MASS_BOUND")
        self.assertEqual(
            payload["implementation_sha256"]["forcing_mass_bound.py"],
            hashlib.sha256((PROJECT / "forcing_mass_bound.py").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            payload["input"]["witness_csv_sha256"],
            hashlib.sha256(DATA.read_bytes()).hexdigest(),
        )
        for row in payload["published_witnesses"]:
            self.assertAlmostEqual(
                row["log_modulus_lower_bound"],
                modulus_log_lower_bound(row["m"], row["k"]),
            )
            if "greedy_size_forcing_modulus" in row:
                self.assertGreater(
                    row["greedy_size_forcing_modulus"]["log_modulus"],
                    row["log_modulus_lower_bound"],
                )
        big = [row for row in payload["published_witnesses"] if row["m"] >= 13]
        self.assertTrue(all(row["bound_exceeds_k"] for row in big))


class GoodDensityDecompositionTests(unittest.TestCase):
    def test_compensation_count_matches_independent_implementation(self) -> None:
        factorizer = SPFFactorizer(60_000)
        checked = 0
        for value in range(50_000, 50_400):
            for p, exponent in factorizer.factor(value).items():
                cofactor = value // p**exponent
                self.assertEqual(
                    compensation_count(cofactor, p),
                    compensation_levels(value, p, exponent),
                )
                checked += 1
        self.assertGreater(checked, 400)

    def test_classes_partition_the_block(self) -> None:
        report = decompose_block(100_000, 2_000, 27)
        self.assertEqual(sum(report["counts"].values()), 2_000)
        self.assertAlmostEqual(sum(report["densities"].values()), 1.0)

    def test_independence_prediction_is_refuted_on_the_artifact(self) -> None:
        payload = json.loads(GOOD_DECOMPOSITION.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "PASS_EXACT_GOOD_DENSITY_DECOMPOSITION")
        self.assertEqual(
            payload["implementation_sha256"]["good_density_decomposition.py"],
            hashlib.sha256(
                (PROJECT / "good_density_decomposition.py").read_bytes()
            ).hexdigest(),
        )
        for row in payload["blocks"]:
            self.assertGreater(row["observed_over_predicted_level_only"], 2.0)
            self.assertEqual(sum(row["counts"].values()), row["count"])

    def test_asymptotic_union_bound_exceeds_one(self) -> None:
        report = asymptotic_union_bound()
        self.assertTrue(report["union_bound_exceeds_one"])
        self.assertAlmostEqual(report["level_failure_sum"], 0.32252, places=4)
        self.assertGreater(report["union_bound_total"], 1.0)

    def test_union_bound_terms_refuses_unsieveable_limits(self) -> None:
        with self.assertRaises(ValueError):
            union_bound_terms(10**15, 27)


class ElementaryEvictionTests(unittest.TestCase):
    def test_class_splitting_is_determined_by_the_class(self) -> None:
        rng = random.Random(1534)
        for _ in range(300):
            modulus = rng.randrange(1, 4_000)
            residue = rng.randrange(modulus)
            position = rng.randrange(1, 40)
            factors = factor_modulus(modulus)
            pinned, cofactor_modulus, cofactor_residue = pinned_split(
                residue, position, modulus, factors
            )
            self.assertEqual(pinned * cofactor_modulus, modulus)
            self.assertEqual(math.gcd(cofactor_residue, cofactor_modulus), 1)
            for step in range(1, 6):
                k = residue + step * modulus
                direct = 1
                for prime, exponent in factors.items():
                    level = 0
                    value = k + position
                    while level < exponent and value % prime == 0:
                        value //= prime
                        level += 1
                    direct *= prime**level
                self.assertEqual(direct, pinned)
                self.assertEqual(((k + position) // pinned) % cofactor_modulus,
                                 cofactor_residue)

    def test_factorisation_lemma_kills_the_pair(self) -> None:
        rng = random.Random(35)
        checked = 0
        for _ in range(200):
            m = rng.randrange(1, 40)
            position = rng.choice(list(bad_window_positions(m)))
            pinned = rng.randrange(1, 30)
            cofactor = rng.randrange(1, 30)
            prime = least_prime_above(max(2 * pinned * cofactor, m))
            k = pinned * prime * cofactor - position
            if k < 1:
                continue
            checked += 1
            self.assertLess(slack(m, k, prime), 0)
            if m + 2 * k <= 40_000:
                self.assertFalse(is_witness(m, k, TrialFactorizer(m + 2 * k)))
        self.assertGreater(checked, 150)

    def test_initial_eviction_meets_its_proved_bound(self) -> None:
        for modulus in range(1, 26):
            factors = factor_modulus(modulus)
            for m in range(1, 14):
                for residue in range(modulus):
                    record = evict_initial(m, modulus, residue, factors)
                    verify_eviction(m, modulus, residue, record)
                    self.assertLess(record["k"], 2 * modulus * max(2 * modulus, m))

    def test_block_eviction_holds_inside_its_hypothesis(self) -> None:
        block_start = 10**7
        for modulus in (1, 2, 6, 30, 97, 210, 1000):
            factors = factor_modulus(modulus)
            for m in (1, 5, 13, 27):
                self.assertLessEqual(modulus * block_prime(m, block_start), block_start)
                for residue in range(0, modulus, max(1, modulus // 7)):
                    record = evict_block(m, modulus, residue, block_start, factors)
                    self.assertIsNotNone(record)
                    verify_eviction(m, modulus, residue, record, block_start=block_start)

    def test_block_hypothesis_is_not_vacuous(self) -> None:
        block_start = 10**6
        m = 3
        modulus = 40 * math.isqrt(block_start)
        self.assertGreater(modulus * block_prime(m, block_start), block_start)
        factors = factor_modulus(modulus)
        failures = sum(
            evict_block(m, modulus, residue, block_start, factors) is None
            for residue in range(0, modulus, modulus // 40)
        )
        self.assertGreater(failures, 0)

    def test_progression_prime_respects_the_25_hypothesis(self) -> None:
        # p**2 > 2*top does not imply p > m at small scales; (25) needs both.
        self.assertGreater(progression_prime(5, 10), 5)
        for m in range(1, 60):
            for top in (1, 2, 7, 10, 50, 10**4, 10**9):
                prime = progression_prime(m, top)
                self.assertGreater(prime, m)
                self.assertGreater(prime * prime, 2 * top)

    def test_master_theorem_at_its_minimal_length(self) -> None:
        for start in (1, 97, 10**5):
            for step in range(1, 13):
                for m in range(1, 13):
                    length = minimal_progression_length(m, start, step)
                    record = evict_progression(m, start, step, length)
                    self.assertIsNotNone(record)
                    k = record["k"]
                    self.assertTrue(start <= k <= start + (length - 1) * step)
                    self.assertEqual((k - start) % step, 0)
                    self.assertLess(slack(m, k, record["p"]), 0)
                    self.assertGreater(record["p"], m)
                    self.assertGreater(record["p"] ** 2, 2 * (k + record["i"]))
                    self.assertIsNone(evict_progression(m, start, step, length - 1))

    def test_master_theorem_ignores_the_common_difference(self) -> None:
        rng = random.Random(3663)
        for _ in range(120):
            m = rng.randrange(1, 50)
            start = rng.randrange(10**6, 10**12)
            step = rng.randrange(1, 10**6)
            length = minimal_progression_length(m, start, step)
            record = evict_progression(m, start, step, length)
            self.assertIsNotNone(record)
            self.assertEqual((record["k"] - start) % step, 0)
            self.assertLess(slack(m, record["k"], record["p"]), 0)

    def test_block_prime_carries_the_m_clause_when_m_dominates(self) -> None:
        # If m^2 > 4N + 2m the size condition alone offers a prime <= m, which
        # (25) exempts. Every screened rectangle has N >> m^2/4, where the clause
        # is vacuous, so only a direct test reaches this regime.
        self.assertLessEqual(math.isqrt(4 * 1000 + 2 * 100), 100)
        for m, block_start in ((100, 1000), (60, 500), (31, 200), (12, 30)):
            prime = block_prime(m, block_start)
            self.assertGreater(prime, m)
            self.assertGreater(prime * prime, 2 * (2 * block_start + m))
            modulus = max(1, block_start // prime)
            self.assertLessEqual(modulus * prime, block_start)
            factors = factor_modulus(modulus)
            for residue in range(modulus):
                record = evict_block(m, modulus, residue, block_start, factors)
                self.assertIsNotNone(record)
                verify_eviction(m, modulus, residue, record, block_start=block_start)

    def test_progression_length_bound_is_the_prime_not_the_square_root(self) -> None:
        # Corollary (c): the threshold is P0, which equals 2*sqrt(2Y)(1+o(1))
        # only once 2Y >= m^2.
        for m in (5, 31, 100):
            for top in (m, 2 * m, m * m // 2, m * m, 10 * m * m):
                prime = progression_prime(m, top)
                self.assertGreater(prime, m)
                if 2 * top < m * m:
                    self.assertEqual(prime, least_prime_above(m))

    def test_covering_measure_matches_direct_enumeration(self) -> None:
        for modulus, block_start, m in ((2_000_000, 10**10, 27), (700_000, 10**9, 13),
                                        (999_983, 10**9, 27), (4_096, 10**7, 5)):
            prime = block_prime(m, block_start)
            hit = bytearray(modulus)
            inverse = pow(prime, -1, modulus)
            for position in bad_window_positions(m):
                low = -(-(block_start + position) // prime)
                high = (2 * block_start - 1 + position) // prime
                for step in range(low, high + 1):
                    hit[(step - position * inverse) % modulus] = 1
            self.assertEqual(
                union_measure(covering_spans(m, block_start, modulus, prime), modulus),
                sum(hit),
            )

    def test_covering_certificate_evicts_every_sampled_class(self) -> None:
        rng = random.Random(4620)
        for block_start, m, theta in ((10**9, 27, 0.55), (10**10, 13, 0.52)):
            modulus = int(round(block_start**theta))
            certificate = covering_prime(m, block_start, modulus)
            self.assertIsNotNone(certificate)
            factors = factor_modulus(modulus)
            for _ in range(25):
                residue = rng.randrange(modulus)
                record = evict_block(
                    m, modulus, residue, block_start, factors,
                    first_prime=certificate["p"],
                )
                self.assertIsNotNone(record)
                verify_eviction(m, modulus, residue, record, block_start=block_start)

    def test_covering_cannot_pass_its_measure_ceiling(self) -> None:
        # The union of L intervals of length floor(N/p) cannot exceed L*floor(N/p).
        block_start, m = 10**12, 27
        window = len(list(bad_window_positions(m)))
        ceiling = window * (block_start // block_prime(m, block_start))
        self.assertIsNone(covering_prime(m, block_start, 3 * ceiling, prime_budget=64))
        for row in json.loads(ELEMENTARY_EVICTION.read_text(encoding="utf-8"))[
            "covering"
        ]["rows"]:
            self.assertLessEqual(row["largest_exponent_certified"],
                                 row["measure_ceiling_exponent"])

    def test_smooth_modulus_position_is_permanently_non_fatal(self) -> None:
        # Proposition 16.2: M m-smooth with M^2 >= 2N leaves M | n non-fatal.
        for block_start, modulus in ((10**7, 4620), (10**6, 2310), (10**5, 420)):
            self.assertGreaterEqual(modulus * modulus, block_start)
            largest_prime = max(factor_modulus(modulus))
            for n in range(
                ((block_start + modulus - 1) // modulus) * modulus,
                2 * block_start,
                modulus,
            ):
                remainder, worst = n, 1
                for prime in sorted(factor_modulus(remainder)):
                    power = 1
                    while remainder % prime == 0:
                        remainder //= prime
                        power *= prime
                    if prime > largest_prime:
                        worst = max(worst, power)
                self.assertLessEqual(worst * worst, 2 * n)

    def test_census_artifact_is_exhaustive_and_certifies_large_moduli(self) -> None:
        payload = json.loads(SURVIVOR_CENSUS.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "PASS_EXACT_SURVIVOR_CENSUS")
        self.assertEqual(
            payload["implementation_sha256"]["survivor_census.py"],
            hashlib.sha256((PROJECT / "survivor_census.py").read_bytes()).hexdigest(),
        )
        self.assertAlmostEqual(payload["survivor_density_over_rho2"], 1.0, delta=0.05)
        for m in ("27", "51"):
            self.assertIsNone(payload["multi_position_by_m"][m]["smallest_exponent_failed"])
            self.assertGreater(
                payload["multi_position_by_m"][m]["largest_exponent_certified"], 0.95
            )
        # every single-position failure must be a smooth modulus at or past sqrt(2N)
        for row in payload["rows"]:
            if not row["every_class_evicted"] and row["shape"] != "prime":
                self.assertEqual(row["worst_class"], 0)

    def test_eviction_artifact_pins_its_inputs(self) -> None:
        payload = json.loads(ELEMENTARY_EVICTION.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "PASS_EXACT_ELEMENTARY_CLASS_EVICTION")
        self.assertEqual(
            payload["implementation_sha256"]["elementary_class_eviction.py"],
            hashlib.sha256(
                (PROJECT / "elementary_class_eviction.py").read_bytes()
            ).hexdigest(),
        )
        self.assertTrue(payload["initial_segment"]["exhaustive"])
        self.assertLess(payload["initial_segment"]["worst_ratio_to_bound"], 1.0)
        self.assertTrue(payload["block"]["exhaustive"])
        self.assertTrue(payload["progression"]["exhaustive"])
        self.assertEqual(
            payload["progression"]["one_shorter_still_worked"], 0
        )
        for row in payload["barrier"]["rows"]:
            if row["log_M_over_log_N"] <= 0.75:
                self.assertEqual(row["successes"], row["classes_tested"])



class OverlapDensityTests(unittest.TestCase):
    def test_level_union_bound_tail_dominates_the_discarded_terms(self) -> None:
        short, long = level_union_bound(12), level_union_bound(400)
        self.assertGreaterEqual(short["bound"], long["partial"])
        self.assertLess(short["bound"] - long["bound"], 1e-3)

    def test_exponent_tail_vanishes_with_scale(self) -> None:
        bounds = [exponent_tail_bound(s, 1)["bound"] for s in (1e8, 1e20, 1e40)]
        self.assertTrue(all(a > b for a, b in zip(bounds, bounds[1:])))
        self.assertLess(bounds[-1], 1e-4)

    def test_riemann_sum_is_a_lower_bound_for_the_band_integral(self) -> None:
        # right endpoints of a decreasing integrand under-count; left endpoints
        # over-count. The certificate must sit below the true integral.
        for band in CERTIFIED_BANDS:
            eta, steps = 0.005, 4000
            lo, hi = 1.0 / (band + 2), 1.0 / (band + 1)
            width = (hi - lo) / steps
            left = sum(
                band_integrand(lo + i * width, band, eta) for i in range(steps)
            ) * width * 2.0 ** (-band)
            certified = band_overlap_lower_bound(band, eta, steps)["bound"]
            self.assertLess(certified, left)
            self.assertGreater(certified, left * 0.99)

    def test_shrinking_eta_cannot_shrink_the_overlap_bound(self) -> None:
        bounds = [overlap_lower_bound(eta)["bound"] for eta in (0.02, 0.005, 0.001)]
        self.assertTrue(all(a < b for a, b in zip(bounds, bounds[1:])))

    def test_band_of_matches_its_defining_inequality(self) -> None:
        for value in (10**3, 10**5, 5 * 10**6):
            for prime in (7, 31, 101, 1009, 10007):
                band = band_of(value, prime)
                if band is None:
                    self.assertGreater(prime * prime, 2 * value)
                    continue
                self.assertLessEqual(prime ** (band + 1), 2 * value)
                self.assertGreater(prime ** (band + 2), 2 * value)

    def test_lower_half_residue_count_meets_lemma_17_5(self) -> None:
        for q in (3, 5, 7, 11, 13, 29, 101):
            for r in CERTIFIED_BANDS:
                good = sum(
                    1
                    for v in range(q**r)
                    if all(2 * (v % q**j) <= q**j for j in range(1, r + 1))
                )
                self.assertGreaterEqual(good, q**r * 2.0 ** (-r) * (1 - 1 / q))

    def test_two_band_primes_cannot_share_a_short_cofactor_witness(self) -> None:
        # Lemma 17.4: p > (2w)^(1/2) and two band-r primes with r <= 2 would
        # force (2w)^(1/2 + 2/(r+2)) <= w. Checked as the exponent inequality.
        # The bound is (2w)^(1/2+2/(r+2)) >= 2w > w, so exponent >= 1 suffices.
        for band in CERTIFIED_BANDS:
            self.assertGreaterEqual(0.5 + 2.0 / (band + 2), 1.0)
        self.assertLess(0.5 + 2.0 / (3 + 2), 1.0)  # band 3 is not free

    def test_family_membership_is_empty_when_m_exceeds_every_band_prime(self) -> None:
        # Small-N / large-m corner: no screened block reaches it, so it needs a
        # direct test. Every band prime of w ~ 10^3 is below m = 100, so the
        # family must be empty and no member may be reported.
        factorizer = SPFFactorizer(4000)
        for value in range(1000, 3000):
            self.assertIsNone(
                overlap_family_member(value, 100, factorizer.factor(value))
            )

    def test_family_members_are_inside_both_failure_classes(self) -> None:
        factorizer = SPFFactorizer(60_000)
        seen = 0
        for value in range(50_000, 60_000):
            band = overlap_family_member(value, 1, factorizer.factor(value))
            if band is None:
                continue
            seen += 1
            self.assertIn(band, CERTIFIED_BANDS)
            short, level = classify_term(value, 1, factorizer)
            self.assertTrue(short and level)
        self.assertGreater(seen, 100)


    def test_pair_obstruction_matches_the_measured_accounting(self) -> None:
        # (53) is derived under cross-term independence; the block must confirm
        # it to within the genuine w / w+1 correlation.
        payload = json.loads(OVERLAP_BOUND.read_text(encoding="utf-8"))
        for row, block in zip(payload["pair_accounting"], payload["blocks"]):
            model = pair_obstruction(block["densities"]["both"], row["S1"] / 2 - math.log(2))
            self.assertAlmostEqual(model["depth3_limit"], row["bonferroni_depth3"], delta=3e-3)

    def test_depth_three_route_does_not_close_r2(self) -> None:
        payload = json.loads(OVERLAP_BOUND.read_text(encoding="utf-8"))
        obstruction = payload["pair_obstruction"]
        self.assertFalse(obstruction["closes"])
        self.assertGreater(obstruction["depth3_limit"], 1.0)
        self.assertGreater(obstruction["shortfall_factor"], 5.0)
        # closure would need an overlap below (u-1)/2, which contradicts 17.6
        self.assertLess(
            obstruction["overlap_threshold_for_closure"],
            payload["certificate"]["overlap_bound"],
        )

    def test_four_event_inclusion_exclusion_is_exact_on_the_block(self) -> None:
        row = pair_accounting(200_000, 20_000, 1)
        self.assertAlmostEqual(
            row["S1"] - row["S2"] + row["S3"] - row["S4"], row["exact_union"], places=12
        )
        self.assertGreaterEqual(row["bonferroni_depth3"], row["exact_union"])
        self.assertGreater(row["both_good"], 0.0)

    def test_finitary_threshold_is_vacuous_before_its_own_scale(self) -> None:
        report = finitary_threshold(27, 0.123664**14)
        self.assertEqual(report["constant"], 972)
        # below the crossover the tier bound exceeds X itself, so the criterion
        # cannot be met; above it, it can.
        below = report["log10_scale_nonvacuous"] - 1.0
        above = report["log10_scale_nonvacuous"] + 1.0
        for log_x, vacuous in ((below, True), (above, False)):
            bound_over_x = math.log10(report["constant"]) - report["decay_exponent"] * log_x
            self.assertEqual(bound_over_x > 0.0, vacuous)
        self.assertGreater(
            report["log10_scale_beats_run_density"], report["log10_scale_nonvacuous"]
        )
        with self.assertRaises(ValueError):
            finitary_threshold(27, 0.0)

    def test_certificate_closes_the_first_order_deficit(self) -> None:
        payload = json.loads(OVERLAP_BOUND.read_text(encoding="utf-8"))
        cert = payload["certificate"]
        self.assertGreater(cert["first_order_deficit"], 0.0)
        self.assertGreater(cert["overlap_bound"], cert["first_order_deficit"])
        self.assertGreater(cert["certified_lower_density"], 0.0)
        self.assertAlmostEqual(
            cert["certified_lower_density"],
            1.0
            - cert["first_order_union_total"]
            - cert["limiting_exponent_tail"]
            + cert["overlap_bound"],
            places=12,
        )
        self.assertEqual(
            payload["implementation_sha256"]["overlap_density_bound.py"],
            hashlib.sha256((PROJECT / "overlap_density_bound.py").read_bytes()).hexdigest(),
        )
        for row in payload["blocks"]:
            self.assertTrue(row["identity_holds"])
            self.assertLessEqual(row["family_density"], row["densities"]["both"])


if __name__ == "__main__":
    unittest.main()
