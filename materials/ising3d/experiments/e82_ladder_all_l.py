"""Exact wave-9 ladder search: sector diagnosis, boundary obstruction, and an L=8 patch."""

from __future__ import annotations

import gc
import hashlib
import itertools
import json
import signal
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "algebra_growth" / "ladder_all_l.json"
SCRIPT = "experiments/e82_ladder_all_l.py"
TIMEOUT_SECONDS = 3_600
P1 = 2_147_483_647
P2 = 2_147_483_629
MOVES = ("AA", "AB", "BA", "BB")
MOVE_PAIRS = tuple(itertools.combinations(MOVES, 2))
BASE_WORDS = tuple(
    "".join(letters)
    for depth in range(1, 5)
    for letters in itertools.product("AB", repeat=depth)
)
EXPECTED_BASE_HISTOGRAM = {0: 1001, 1: 8269, 2: 12708, 3: 4991, 4: 436}
EXPECTED_STAGE_HISTOGRAMS = {
    3: {3: 8, 4: 88, 5: 520, 6: 920, 7: 672, 8: 408},
    4: {14: 48, 15: 32, 16: 328},
    5: {30: 64, 31: 32, 32: 232},
    6: {61: 40, 63: 48, 64: 144},
    7: {128: 144},
}
EXPECTED_BOUNDARY_HISTOGRAMS = {
    2: {4: 436},
    3: {8: 408},
    4: {16: 328},
    5: {32: 232},
    6: {64: 144},
    7: {128: 144},
    8: {248: 144},
}
EXPECTED_REMOVED_LABELS = (
    "10111011",
    "10111101",
    "10111110",
    "11111010",
    "11111100",
    "11111101",
    "11111110",
    "11111111",
)
PATCH_EXTENSION_LABELS = (
    "000000",
    "100000",
    "010000",
    "110000",
    "001000",
    "000100",
    "100100",
    "010100",
)
EXPECTED_SUPPLEMENT_WORDS = (
    "ABABABABABABBAB",
    "ABABABABABBBBAB",
    "ABABABABBBABBAB",
    "ABABABABBBBBBAB",
    "ABABABBBABABBAB",
    "ABABBBABABABBAB",
    "ABABBBABABBBBAB",
    "ABABBBABBBABBAB",
)
EXPECTED_RELATION_SECTOR_COUNTS = {(0, 0): 3, (0, 1): 0, (1, 0): 3, (1, 1): 2}
EXPECTED_SECTOR_RANKS = {(0, 0): 61, (0, 1): 64, (1, 0): 61, (1, 1): 62}
KNOWN_RELATIONS = json.loads(r"""[{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":187,"pivot_label":"10111011","term_count":19,"terms":[{"coefficient":-7853057352327686058,"label":"00001111"},{"coefficient":25241922304480834876,"label":"00010111"},{"coefficient":-129304695571775854172,"label":"00011011"},{"coefficient":-982063996856384105843,"label":"00011101"},{"coefficient":-1130883532624365454809,"label":"00011110"},{"coefficient":18546535327501495160,"label":"00100111"},{"coefficient":111427725120416217380,"label":"00101011"},{"coefficient":1525779047982068519836,"label":"00101101"},{"coefficient":1854459865558115245710,"label":"00101110"},{"coefficient":-29704159155083841540,"label":"00110011"},{"coefficient":-599664732791886586258,"label":"00110101"},{"coefficient":-892773422107427928438,"label":"00110110"},{"coefficient":64232776834749635564,"label":"00111001"},{"coefficient":181577399782667428260,"label":"00111010"},{"coefficient":-12970915610392028898,"label":"00111100"},{"coefficient":-1370108267480073988,"label":"10011111"},{"coefficient":3907437972129356616,"label":"10101111"},{"coefficient":-7067278625054107992,"label":"10110111"},{"coefficient":5906044059409481224,"label":"10111011"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":189,"pivot_label":"10111101","term_count":19,"terms":[{"coefficient":410018921829350100480,"label":"00001111"},{"coefficient":-597664541827002740440,"label":"00010111"},{"coefficient":-1913599910285479426864,"label":"00011011"},{"coefficient":-12739071720422560410235,"label":"00011101"},{"coefficient":-12834916346092292663508,"label":"00011110"},{"coefficient":-196649649007257787600,"label":"00100111"},{"coefficient":2499565473004361154680,"label":"00101011"},{"coefficient":19808813423773802262940,"label":"00101101"},{"coefficient":21042209757407549874360,"label":"00101110"},{"coefficient":-421253234997969115896,"label":"00110011"},{"coefficient":-7952257389316044495250,"label":"00110101"},{"coefficient":-10149653064136017747000,"label":"00110110"},{"coefficient":866408697840914340860,"label":"00111001"},{"coefficient":2090659304005573094160,"label":"00111010"},{"coefficient":-154511560691192643912,"label":"00111100"},{"coefficient":-20265919654888586812,"label":"10011111"},{"coefficient":5521953758333094000,"label":"10101111"},{"coefficient":42203236470671053064,"label":"10110111"},{"coefficient":11812088118818962448,"label":"10111101"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":190,"pivot_label":"10111110","term_count":6,"terms":[{"coefficient":-1105002572961,"label":"00011110"},{"coefficient":1784495478700,"label":"00101110"},{"coefficient":-823704922418,"label":"00110110"},{"coefficient":151372442200,"label":"00111010"},{"coefficient":-8978829696,"label":"00111100"},{"coefficient":336529200,"label":"10111110"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":250,"pivot_label":"11111010","term_count":20,"terms":[{"coefficient":-310008209268179671040,"label":"01001110"},{"coefficient":2246291314851192125440,"label":"01010110"},{"coefficient":-717021822711198090240,"label":"01011010"},{"coefficient":-186095597210644915200,"label":"01011100"},{"coefficient":-1525708395796177980348,"label":"01011111"},{"coefficient":-979256275959985013760,"label":"01100110"},{"coefficient":-1307951854194832000000,"label":"01101010"},{"coefficient":351773837699595351040,"label":"01101100"},{"coefficient":4198909491492255989735,"label":"01101111"},{"coefficient":1093903526237160826880,"label":"01110010"},{"coefficient":-33289108921224611840,"label":"01110100"},{"coefficient":-3670633689312159538254,"label":"01110111"},{"coefficient":-85172698237052180480,"label":"01111000"},{"coefficient":1059763007412959446740,"label":"01111011"},{"coefficient":-51382448844264700448,"label":"01111101"},{"coefficient":-2953297184266579200,"label":"01111110"},{"coefficient":25271419521566648320,"label":"11011110"},{"coefficient":-37429519266826086400,"label":"11101110"},{"coefficient":-66947928524085790720,"label":"11110110"},{"coefficient":64740219659996108800,"label":"11111010"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":252,"pivot_label":"11111100","term_count":20,"terms":[{"coefficient":-249616476152573742080,"label":"01001110"},{"coefficient":-554295721845334515200,"label":"01010110"},{"coefficient":-591254974616842708480,"label":"01011010"},{"coefficient":-9861852161930437672960,"label":"01011100"},{"coefficient":-21153578227030538231124,"label":"01011111"},{"coefficient":724322930446990533120,"label":"01100110"},{"coefficient":1531397607678326080000,"label":"01101010"},{"coefficient":13941326517406904295680,"label":"01101100"},{"coefficient":58216865451639249159805,"label":"01101111"},{"coefficient":-888643909445963095040,"label":"01110010"},{"coefficient":-5868302381544926794240,"label":"01110100"},{"coefficient":-50892449109922494676602,"label":"01110111"},{"coefficient":623394035275852262400,"label":"01111000"},{"coefficient":14693358010738777990620,"label":"01111011"},{"coefficient":-712405236884304387424,"label":"01111101"},{"coefficient":-40946751808658169600,"label":"01111110"},{"coefficient":-16270393059099182080,"label":"11011110"},{"coefficient":45411620901694361600,"label":"11101110"},{"coefficient":-10391884340236241920,"label":"11110110"},{"coefficient":194220658979988326400,"label":"11111100"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":253,"pivot_label":"11111101","term_count":20,"terms":[{"coefficient":183109289574516060600,"label":"01001111"},{"coefficient":-231825546426065043200,"label":"01010111"},{"coefficient":511275231050049657800,"label":"01011011"},{"coefficient":55279973186880427200,"label":"01011101"},{"coefficient":49137669366357595236,"label":"01011110"},{"coefficient":-196534816297992399500,"label":"01100111"},{"coefficient":-428052330872590774000,"label":"01101011"},{"coefficient":-348273395585543316000,"label":"01101101"},{"coefficient":-158044092144028744545,"label":"01101110"},{"coefficient":190702968708533692600,"label":"01110011"},{"coefficient":318461120889391928000,"label":"01110101"},{"coefficient":173337985719966030978,"label":"01110110"},{"coefficient":-71624707165366338000,"label":"01111001"},{"coefficient":-74016213291987024180,"label":"01111010"},{"coefficient":9221376565097291136,"label":"01111100"},{"coefficient":16782502125013864600,"label":"11011111"},{"coefficient":-50462052833837049000,"label":"11101111"},{"coefficient":59489767439994274600,"label":"11110111"},{"coefficient":-35493082636656585000,"label":"11111011"},{"coefficient":6638092273809956800,"label":"11111101"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":254,"pivot_label":"11111110","term_count":20,"terms":[{"coefficient":1867972045274799266544000,"label":"01001111"},{"coefficient":-2022216625306127047504000,"label":"01010111"},{"coefficient":5123509005506611101312000,"label":"01011011"},{"coefficient":5538321165666545043008000,"label":"01011101"},{"coefficient":1971348384002676448487012,"label":"01011110"},{"coefficient":-3074359569694063667800000,"label":"01100111"},{"coefficient":-4114948090955459765280000,"label":"01101011"},{"coefficient":-10459115007453535362520000,"label":"01101101"},{"coefficient":-6540826779163401765304165,"label":"01101110"},{"coefficient":2563872878195062369264000,"label":"01110011"},{"coefficient":5666231704163240549136000,"label":"01110101"},{"coefficient":6587469774853802724069026,"label":"01110110"},{"coefficient":-911068077592261714920000,"label":"01111001"},{"coefficient":-2291879934442449233882910,"label":"01111010"},{"coefficient":204245069051358921750912,"label":"01111100"},{"coefficient":301598105131555650512000,"label":"11011111"},{"coefficient":-853238891598159201856000,"label":"11101111"},{"coefficient":968075380535011262528000,"label":"11110111"},{"coefficient":-528603147564957019824000,"label":"11111011"},{"coefficient":10977745097813216058000,"label":"11111110"}]},{"coefficient_gcd":1,"exact_residual_support_size":0,"pivot_index":255,"pivot_label":"11111111","term_count":7,"terms":[{"coefficient":-3215763564,"label":"01011111"},{"coefficient":9056323925,"label":"01101111"},{"coefficient":-9030348822,"label":"01110111"},{"coefficient":3800131005,"label":"01111011"},{"coefficient":-668129264,"label":"01111101"},{"coefficient":38382720,"label":"01111110"},{"coefficient":4851000,"label":"11111111"}]}]""")


def generator_terms(length: int) -> dict[str, tuple[int, ...]]:
    if length < 2:
        raise ValueError("the open ladder has at least two rungs")
    n = 2 * length
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return {
        "A": tuple(1 << site for site in range(n)),
        "B": tuple((1 << (n + left)) | (1 << (n + right)) for left, right in bonds),
    }


@lru_cache(maxsize=None)
def cached_terms(length: int) -> dict[str, tuple[int, ...]]:
    return generator_terms(length)


def ladder_bonds(length: int) -> tuple[tuple[int, int], ...]:
    bonds = []
    for row in range(2):
        offset = row * length
        bonds.extend((offset + col, offset + col + 1) for col in range(length - 1))
    bonds.extend((col, length + col) for col in range(length))
    return tuple(bonds)


def sign_half(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    left_right = ((((left >> n) & mask) & (right & mask)).bit_count()) & 1
    right_left = ((((right >> n) & mask) & (left & mask)).bit_count()) & 1
    if left_right == right_left:
        return 0
    return 1 if left_right == 0 else -1


def adjoint(terms: tuple[int, ...], vector: dict[int, int], n: int) -> dict[int, int]:
    result: dict[int, int] = {}
    for term in terms:
        for pauli, coefficient in vector.items():
            sign = sign_half(term, pauli, n)
            if not sign:
                continue
            target = term ^ pauli
            updated = result.get(target, 0) + sign * coefficient
            if updated:
                result[target] = updated
            else:
                result.pop(target, None)
    return result


def vector_for_word(word: str, length: int) -> dict[int, int]:
    terms = cached_terms(length)
    value = {term: 1 for term in terms[word[-1]]}
    for letter in reversed(word[:-1]):
        value = adjoint(terms[letter], value, 2 * length)
    return value


@lru_cache(maxsize=None)
def cached_vector(length: int, word: str) -> dict[int, int]:
    return vector_for_word(word, length)


def rank_mod(rows: list[dict[int, int]], prime: int) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for row in rows:
        value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(value[lead], prime - 2, prime)
                pivots[lead] = {
                    key: coefficient * inverse % prime
                    for key, coefficient in value.items()
                    if coefficient * inverse % prime
                }
                break
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = (value.get(key, 0) - factor * coefficient) % prime
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def exact_rank(rows: list[dict[int, int]]) -> int:
    """Sparse rational echelon rank, used only when a modular full minor is absent."""
    pivots: dict[int, dict[int, Fraction]] = {}
    for row in rows:
        value = {key: Fraction(coefficient) for key, coefficient in row.items() if coefficient}
        while value:
            lead = min(value)
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = 1 / value[lead]
                pivots[lead] = {key: coefficient * inverse for key, coefficient in value.items()}
                break
            factor = value[lead]
            for key, coefficient in pivot.items():
                reduced = value.get(key, Fraction(0)) - factor * coefficient
                if reduced:
                    value[key] = reduced
                else:
                    value.pop(key, None)
    return len(pivots)


def certified_rank(rows: list[dict[int, int]], target: int) -> int:
    modular = rank_mod(rows, P1)
    if modular == target:
        return target
    return exact_rank(rows)


def family_words(base: tuple[str, ...], moves: tuple[str, str], length: int) -> list[str]:
    words = list(base)
    for _ in range(2, length):
        words = [move + word for word in words for move in moves]
    return words


def family_rows(base: tuple[str, ...], moves: tuple[str, str], length: int) -> list[dict[int, int]]:
    words = family_words(base, moves, length)
    assert len(words) == 1 << length
    assert max(map(len, words)) <= 2 * length
    return [cached_vector(length, word) for word in words]


def project_last_rungs(vector: dict[int, int], length: int, width: int) -> dict[int, int]:
    n = 2 * length
    row_mask = ((1 << width) - 1) << (length - width)
    site_mask = row_mask | (row_mask << length)
    pauli_mask = site_mask | (site_mask << n)
    return {code: coefficient for code, coefficient in vector.items() if not code & ~pauli_mask}


def vector_digest(vector: dict[int, int]) -> str:
    digest = hashlib.sha256()
    for code, coefficient in sorted(vector.items()):
        digest.update(f"{code}:{coefficient};".encode("ascii"))
    return digest.hexdigest()


def words_digest(words: list[str]) -> str:
    digest = hashlib.sha256()
    for word in words:
        digest.update(f"{word};".encode("ascii"))
    return digest.hexdigest()


def original_recursive_words(length: int) -> list[tuple[str, str]]:
    family = [("00", "A"), ("01", "BA"), ("10", "ABA"), ("11", "BABA")]
    for _ in range(2, length):
        family = [
            child
            for label, word in family
            for child in ((label + "0", "AB" + word), (label + "1", "BB" + word))
        ]
    return family


def original_word_for_label(label: str) -> str:
    return dict(original_recursive_words(8))[label]


def supplemented_word(extension_label: str) -> str:
    word = "BAB"
    for bit in extension_label:
        word = ("AB" if bit == "0" else "BB") + word
    return word


def run_natural_search() -> tuple[dict[str, object], dict[str, dict[int, int]]]:
    base_histogram: Counter[int] = Counter()
    base_sets: list[tuple[str, ...]] = []
    for base in itertools.combinations(BASE_WORDS, 4):
        rank = certified_rank([cached_vector(2, word) for word in base], 4)
        base_histogram[rank] += 1
        if rank == 4:
            base_sets.append(base)
    assert dict(base_histogram) == EXPECTED_BASE_HISTOGRAM
    assert len(base_sets) == 436

    boundary_histograms: dict[int, Counter[int]] = {2: Counter()}
    for base in base_sets:
        rows = [
            project_last_rungs(cached_vector(2, word), 2, 2)
            for word in base
        ]
        boundary_histograms[2][certified_rank(rows, 4)] += 1
    assert dict(boundary_histograms[2]) == EXPECTED_BOUNDARY_HISTOGRAMS[2]

    stage_records = []
    candidates = [(moves, base) for moves in MOVE_PAIRS for base in base_sets]
    for length, target in ((3, 8), (4, 16), (5, 32), (6, 64), (7, 128)):
        histogram: Counter[int] = Counter()
        survivors: list[tuple[tuple[str, str], tuple[str, ...]]] = []
        for moves, base in candidates:
            rank = certified_rank(family_rows(base, moves, length), target)
            histogram[rank] += 1
            if rank == target:
                survivors.append((moves, base))
        assert dict(histogram) == EXPECTED_STAGE_HISTOGRAMS[length]
        if length == 5:
            assert {moves for moves, _ in survivors} == {("AB", "BB")}
        if length == 7:
            moves, base = survivors[0]
            assert rank_mod(family_rows(base, moves, length), P2) == target

        boundary_histogram: Counter[int] = Counter()
        for moves, base in survivors:
            rows = [
                project_last_rungs(vector, length, min(3, length))
                for vector in family_rows(base, moves, length)
            ]
            boundary_histogram[certified_rank(rows, target)] += 1
        assert dict(boundary_histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[length]
        boundary_histograms[length] = boundary_histogram
        stage_records.append(
            {
                "length": length,
                "rank_histogram": dict(sorted(histogram.items())),
                "survivor_count": len(survivors),
            }
        )
        candidates = survivors

    assert len(candidates) == 144
    seed_words = sorted({word for _, base in candidates for word in base})
    assert len(seed_words) == 15
    boundary_by_word: dict[str, dict[int, int]] = {}
    for seed in seed_words:
        for word in family_words((seed,), ("AB", "BB"), 8):
            boundary_by_word[word] = project_last_rungs(vector_for_word(word, 8), 8, 3)

    boundary_l8_histogram: Counter[int] = Counter()
    survivor_bases = []
    for moves, base in candidates:
        assert moves == ("AB", "BB")
        rows = [boundary_by_word[word] for word in family_words(base, moves, 8)]
        boundary_l8_histogram[exact_rank(rows)] += 1
        survivor_bases.append(base)
    assert dict(boundary_l8_histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[8]
    representative_rows = [
        boundary_by_word[word] for word in family_words(candidates[0][1], candidates[0][0], 8)
    ]
    modular_l8_representative = {
        str(prime): rank_mod(representative_rows, prime) for prime in (P1, P2)
    }
    assert all(rank == 248 for rank in modular_l8_representative.values())
    for length, histogram in boundary_histograms.items():
        assert dict(histogram) == EXPECTED_BOUNDARY_HISTOGRAMS[length]

    record = {
        "base_rank_histogram": dict(sorted(base_histogram.items())),
        "stage_records": stage_records,
        "boundary_rank_histograms": {
            str(length): dict(sorted(histogram.items()))
            for length, histogram in sorted(boundary_histograms.items())
        },
        "boundary_rank_histogram_L8_exact": dict(sorted(boundary_l8_histogram.items())),
        "boundary_rank_L8_modular_representative": modular_l8_representative,
        "survivor_bases": [list(base) for base in survivor_bases],
        "survivor_base_digest": words_digest(["/".join(base) for base in survivor_bases]),
        "seed_words_L8": seed_words,
    }
    return record, boundary_by_word


def permute_pauli(code: int, permutation: tuple[int, ...], n: int) -> int:
    mask = (1 << n) - 1
    x, z = code & mask, (code >> n) & mask
    mapped_x = 0
    mapped_z = 0
    for old, new in enumerate(permutation):
        mapped_x |= ((x >> old) & 1) << new
        mapped_z |= ((z >> old) & 1) << new
    return mapped_x | (mapped_z << n)


def relation_structure(boundary_by_word: dict[str, dict[int, int]]) -> dict[str, object]:
    original = original_recursive_words(8)
    original_rows = {label: boundary_by_word[word] for label, word in original}
    sector_rows: dict[tuple[int, int], list[dict[int, int]]] = {}
    for label, word in original:
        sector = (word.count("A") & 1, word.count("B") & 1)
        sector_rows.setdefault(sector, []).append(original_rows[label])
    assert {sector: len(rows) for sector, rows in sector_rows.items()} == {
        sector: 64 for sector in EXPECTED_RELATION_SECTOR_COUNTS
    }

    relation_sector_counts: Counter[tuple[int, int]] = Counter()
    relation_rows_by_sector: dict[tuple[int, int], list[dict[int, int]]] = {
        sector: [] for sector in EXPECTED_RELATION_SECTOR_COUNTS
    }
    residual_sizes = []
    for relation in KNOWN_RELATIONS:
        terms = relation["terms"]
        sectors = {
            (
                original_word_for_label(term["label"]).count("A") & 1,
                original_word_for_label(term["label"]).count("B") & 1,
            )
            for term in terms
        }
        assert len(sectors) == 1
        sector = sectors.pop()
        relation_sector_counts[sector] += 1
        relation_rows_by_sector[sector].append(
            {int(term["label"], 2): int(term["coefficient"]) for term in terms}
        )
        residual: dict[int, int] = {}
        for term in terms:
            vector = cached_vector(8, original_word_for_label(term["label"]))
            coefficient = int(term["coefficient"])
            for code, value in vector.items():
                total = residual.get(code, 0) + coefficient * value
                if total:
                    residual[code] = total
                else:
                    residual.pop(code, None)
        assert not residual
        residual_sizes.append(len(residual))
    assert {
        sector: relation_sector_counts[sector] for sector in EXPECTED_RELATION_SECTOR_COUNTS
    } == EXPECTED_RELATION_SECTOR_COUNTS

    row_swap = tuple(list(range(8, 16)) + list(range(8)))
    rung_reverse = tuple(list(range(7, -1, -1)) + list(range(15, 7, -1)))
    terms = generator_terms(8)
    for permutation in (row_swap, rung_reverse):
        for generator in terms.values():
            assert {permute_pauli(code, permutation, 16) for code in generator} == set(generator)
    assert all(
        (((left % 8) + (left // 8)) & 1) != (((right % 8) + (right // 8)) & 1)
        for left, right in ladder_bonds(8)
    )

    sector_records = []
    for sector in sorted(EXPECTED_RELATION_SECTOR_COUNTS):
        relation_rank = exact_rank(relation_rows_by_sector[sector])
        relation_count = EXPECTED_RELATION_SECTOR_COUNTS[sector]
        assert relation_rank == relation_count
        modular = {str(prime): rank_mod(sector_rows[sector], prime) for prime in (P1, P2)}
        exact_full_rank = 64 - relation_count
        assert all(rank == exact_full_rank for rank in modular.values())
        assert exact_full_rank == EXPECTED_SECTOR_RANKS[sector]
        sector_records.append(
            {
                "sector": f"{sector[0]}{sector[1]}",
                "family_size": 64,
                "independent_exact_relation_count": relation_count,
                "relation_coefficient_rank_over_Q": relation_rank,
                "boundary_rank_mod": modular,
                "full_family_rank_over_Q": exact_full_rank,
            }
        )
    return {
        "exact_relations": KNOWN_RELATIONS,
        "exact_residual_support_sizes": residual_sizes,
        "relation_sector_counts": {
            f"{sector[0]}{sector[1]}": EXPECTED_RELATION_SECTOR_COUNTS[sector]
            for sector in sorted(EXPECTED_RELATION_SECTOR_COUNTS)
        },
        "sector_records": sector_records,
        "symmetry_certificate": {
            "row_swap_and_rung_reversal": "Each fixes A and B termwise as a set, so every bracket word is fixed.",
            "global_Z_character": "Conjugation by product_v Z_v sends A to -A and B to B.",
            "staggered_X_character": "On the bipartite ladder, conjugation by X on one color class sends A to A and B to -B.",
            "conclusion": "Every recorded relation is homogeneous in (#A mod 2, #B mod 2); reflection and row swap act trivially on every word, so no fixed symmetry sector removes a dependency.",
        },
    }


def echelon_add(pivots: dict[int, dict[int, int]], row: dict[int, int], prime: int) -> int | None:
    value = {key: coefficient % prime for key, coefficient in row.items() if coefficient % prime}
    while value:
        lead = min(value)
        pivot = pivots.get(lead)
        if pivot is None:
            inverse = pow(value[lead], prime - 2, prime)
            pivots[lead] = {
                key: coefficient * inverse % prime
                for key, coefficient in value.items()
                if coefficient * inverse % prime
            }
            return lead
        factor = value[lead]
        for key, coefficient in pivot.items():
            reduced = (value.get(key, 0) - factor * coefficient) % prime
            if reduced:
                value[key] = reduced
            else:
                value.pop(key, None)
    return None


def patched_family(boundary_by_word: dict[str, dict[int, int]]) -> dict[str, object]:
    original = original_recursive_words(8)
    original_rows = [boundary_by_word[word] for _, word in original]
    assert exact_rank(original_rows) == 248
    pivots: dict[int, dict[int, int]] = {}
    observed_removed = []
    kept = []
    for label, word in original:
        if echelon_add(pivots, boundary_by_word[word], P1) is None:
            observed_removed.append(label)
        else:
            kept.append((label, word))
    assert tuple(observed_removed) == EXPECTED_REMOVED_LABELS
    assert len(kept) == 248

    supplements = []
    for extension_label, expected_word in zip(PATCH_EXTENSION_LABELS, EXPECTED_SUPPLEMENT_WORDS):
        word = supplemented_word(extension_label)
        assert word == expected_word
        assert len(word) == 15
        supplements.append(
            {
                "seed": "BAB",
                "extension_label": extension_label,
                "word": word,
                "depth": len(word),
            }
        )
    words = [word for _, word in kept] + [entry["word"] for entry in supplements]
    rows = [boundary_by_word[word] for word in words]
    ranks = {str(prime): rank_mod(rows, prime) for prime in (P1, P2)}
    assert all(rank == 256 for rank in ranks.values())
    return {
        "construction": "Keep the p1-echelon 248 rows of the original AB/BB family and replace its eight dependency pivots by BAB-seed descendants.",
        "removed_original_labels": observed_removed,
        "kept_original_labels": [label for label, _ in kept],
        "supplements": supplements,
        "family_size": len(words),
        "maximum_depth": max(map(len, words)),
        "word_sha256": words_digest(words),
        "boundary_row_digests_sha256": [vector_digest(row) for row in rows],
        "boundary_rank_mod": ranks,
    }


def make_artifact() -> dict[str, object]:
    started = time.monotonic()
    natural_class_search, boundary_by_word = run_natural_search()
    cached_vector.cache_clear()
    gc.collect()
    relations = relation_structure(boundary_by_word)
    cached_vector.cache_clear()
    gc.collect()
    patch = patched_family(boundary_by_word)
    elapsed = time.monotonic() - started
    checks = [
        {
            "name": "exact_depth_compatible_base_enumeration",
            "passed": True,
            "detail": "exact rational fallback proves the complete 30-choose-4 base histogram and 436 independent bases",
        },
        {
            "name": "stationary_two_letter_search_through_L7",
            "passed": True,
            "detail": "the exact staged search leaves 144 AB/BB seed sets, each full through L=7",
        },
        {
            "name": "last_three_rung_boundary_obstruction",
            "passed": True,
            "detail": "all 144 viable stationary families have exact boundary rank 248, not 256, at L=8",
        },
        {
            "name": "primitive_relation_sector_analysis",
            "passed": True,
            "detail": "eight exact full-Pauli residual-zero relations split as 3+0+3+2 under the two character parities",
        },
        {
            "name": "patched_depth_16_family",
            "passed": True,
            "detail": "248 retained words plus eight BAB descendants have boundary rank 256 modulo both primes",
        },
        {
            "name": "honest_scope",
            "passed": True,
            "detail": "the L=8 lower bound is finite; the all-length inequality and non-boundary recursions remain unresolved",
        },
    ]
    return {
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "script": SCRIPT,
            "method": "exact integer Pauli Q_(a|b) commutators; Fraction echelon for all modular-deficient search branches and L=8 boundary ranks; two-prime modular lower bounds for full minors",
            "wall_seconds": round(elapsed, 6),
        },
        "data": {
            "status": "UNRESOLVED",
            "target": "[UNRESOLVED] Prove D_(2L) >= 2^L for every L>=2 by an all-length word recursion.",
            "finite_claim": "[COMPUTATION] The explicit 256-word depth-<=16 patched family has rank 256 over Q at L=8; hence D_16 >= 256.",
            "natural_class_definition": {
                "base_words": "all 30 literal A/B words of depths 1 through 4; use every unordered four-word base set with exact L=2 rank four",
                "extension_moves": "all six unordered distinct pairs from {AA, AB, BA, BB}; W_(L+1) contains mW for every prior W and each move m",
                "depth_control": "every descendant has depth <= 2L",
                "scope": "This is an exact obstruction only to a full-rank last-three-rung boundary-transfer induction in this stationary two-letter class.",
            },
            "natural_class_search": natural_class_search,
            "relation_structure": relations,
            "patched_family": patch,
            "limitations": [
                "Exact boundary rank 248 does not by itself upper-bound the full Pauli rank of every stationary candidate.",
                "The patched L=8 family is a finite certificate, not a recursive all-L construction.",
                "No claim is made about D_(2L) for arbitrary L beyond the certified finite lengths.",
            ],
        },
        "checks": checks,
    }


def write_artifact(artifact: dict[str, object]) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_PATH.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")


def timeout(_signum: int, _frame: object) -> None:
    raise TimeoutError(f"NON-DECISIVE: exceeded {TIMEOUT_SECONDS}s")


def main() -> int:
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        artifact = make_artifact()
        write_artifact(artifact)
        for check in artifact["checks"]:
            print(f"PASS {check['name']}: {check['detail']}")
        print("PASS")
        return 0
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    raise SystemExit(main())
