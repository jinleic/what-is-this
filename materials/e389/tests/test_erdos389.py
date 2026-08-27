from __future__ import annotations

import csv
import hashlib
import json
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
)
from artifact_io import atomic_write_json
from compensation_structure_analysis import (
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
from bad_window_near_miss_analysis import analyze as analyze_bad_window_near_misses
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
COMPENSATION_RUN_FIRST_100M = (
    DATA.parent / "compensation_run_m27_k50001_h100000000.json"
)
COMPENSATION_RUN_NEXT_900M = (
    DATA.parent / "compensation_run_m27_k100050001_h900000000.json"
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
        self.assertEqual(
            structure["small_and_large_prime_separation"]["counts"],
            {
                "large_good_but_small_prime_failed": 190,
                "large_prime_factor_checks": 587_169,
                "large_prime_good_windows": 1_294,
                "large_prime_obstructions": 442_704,
                "pairs": 100_000,
                "witnesses": 1_104,
            },
        )
        self.assertEqual(
            structure["even_shift_first_bad_term_identity"]["counts"],
            {"prime_checks": 19_654, "spikes": 15_336},
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
            },
        )
        self.assertEqual(
            [
                row["retained_prime_product_decimal_digits"]
                for row in extended_retained["rows"]
            ],
            [36, 74, 79, 45, 79, 78, 42, 48, 48],
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


if __name__ == "__main__":
    unittest.main()
