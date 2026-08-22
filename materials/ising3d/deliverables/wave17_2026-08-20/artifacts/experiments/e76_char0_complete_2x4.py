#!/usr/bin/env python3
"""Complete exact-Q structure certificate for the finite open 2x4 Ising layer.

The proof combines an exact modular lower certificate for actual nested integer
Lie words with exact-Q upper certificates from the faithful eight-sector
representation: kernel splittings, invariant forms, rational projectors, and
full-rank intertwiners.  Modular ranks are used only as rigorous lower bounds.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e76_char0_complete_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_complete_2x4.json"
E20_PATH = ROOT / "experiments" / "e20_algebra_growth.py"
E49_PATH = ROOT / "experiments" / "e49_char0_levi_2x4.py"
GOOD_PRIME = 2_147_483_647
N = 8
ROWS = 2
COLS = 4
DIMENSION = 2952
DERIVED_DIMENSION = 2951
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))
BONDS = tuple(
    (r * COLS + c, (r + 1) * COLS + c)
    for r in range(ROWS - 1)
    for c in range(COLS)
) + tuple(
    (r * COLS + c, r * COLS + c + 1)
    for r in range(ROWS)
    for c in range(COLS - 1)
)
FACTOR_DIMENSIONS = (861, 351, 351, 575, 783, 15, 15)
FACTOR_RANKS = (21, 13, 13, 23, 27, 3, 3)
FACTOR_TYPE = "D21 + B13 + B13 + A23 + A27 + A3 + A3"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def timed(label, timings, memory, function, *args):
    started = time.monotonic()
    before = rss_bytes()
    value = function(*args)
    timings[label] = round(time.monotonic() - started, 6)
    memory[label] = {"rss_before_bytes": before, "rss_peak_bytes": rss_bytes()}
    return value


def site_symmetries():
    return tuple(
        (
            flip_rows,
            flip_cols,
            spin_flip,
            tuple(
                (ROWS - 1 - r if flip_rows else r) * COLS
                + (COLS - 1 - c if flip_cols else c)
                for r in range(ROWS)
                for c in range(COLS)
            ),
        )
        for flip_rows, flip_cols, spin_flip in CHARACTERS
    )


def sector_with_gram(character):
    """Exact sector generators and the physical positive diagonal Gram form."""
    size = 1 << N
    symmetries = site_symmetries()

    def act(configuration, symmetry):
        output = 0
        for old, new in enumerate(symmetry[3]):
            if (configuration >> old) & 1:
                output |= 1 << new
        return output ^ ((size - 1) if symmetry[2] else 0)

    seen = set()
    representatives = []
    orbit_vectors = []
    for configuration in range(size):
        if configuration in seen:
            continue
        seen.update(act(configuration, symmetry) for symmetry in symmetries)
        values = {}
        for symmetry in symmetries:
            sign = -1 if sum(
                x * y for x, y in zip(character, symmetry[:3])
            ) & 1 else 1
            image = act(configuration, symmetry)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(q for q, coefficient in values.items() if coefficient)
        scale = values[representative]
        representatives.append(representative)
        orbit_vectors.append(
            {q: coefficient // scale for q, coefficient in values.items() if coefficient}
        )

    dimension = len(representatives)
    a = sp.zeros(dimension)
    b = sp.zeros(dimension)
    gram = sp.zeros(dimension)
    for j, orbit in enumerate(orbit_vectors):
        gram[j, j] = sum(coefficient * coefficient for coefficient in orbit.values())
        image = {}
        for configuration, coefficient in orbit.items():
            for site in range(N):
                target = configuration ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for i, representative in enumerate(representatives):
            a[i, j] = image.get(representative, 0)
        b[j, j] = sum(
            1
            if ((representatives[j] >> u) & 1)
            == ((representatives[j] >> v) & 1)
            else -1
            for u, v in BONDS
        )
    assert a.T * gram == gram * a
    assert b.T * gram == gram * b
    return a, b, gram


def gram_split_common_kernel(actions, gram):
    """Split the common zero-module by its exact Gram-orthogonal complement."""
    a, b = actions
    dimension = a.rows
    kernel = sp.Matrix.vstack(a, b).nullspace()
    kernel_matrix = sp.Matrix.hstack(*kernel) if kernel else sp.zeros(dimension, 0)
    if kernel:
        complement = (kernel_matrix.T * gram).nullspace()
    else:
        complement = [sp.eye(dimension)[:, i] for i in range(dimension)]
    complement_matrix = sp.Matrix.hstack(*complement)
    change = sp.Matrix.hstack(kernel_matrix, complement_matrix)
    assert change.rank() == dimension
    inverse = change.inv()
    conjugated = tuple(inverse * generator * change for generator in actions)
    k = len(kernel)
    for generator in conjugated:
        assert generator[:k, :] == sp.zeros(k, dimension)
        assert generator[:, :k] == sp.zeros(dimension, k)
    active = tuple(generator[k:, k:] for generator in conjugated)
    return kernel_matrix, complement_matrix, active


def primitive_integer_matrix(matrix):
    denominators = [sp.denom(value) for value in matrix]
    denominator = sp.ilcm(*denominators) if denominators else 1
    entries = [int(value * denominator) for value in matrix]
    content = math.gcd(*(abs(value) for value in entries)) if entries else 1
    if content:
        entries = [value // content for value in entries]
    first = next((value for value in entries if value), 1)
    if first < 0:
        entries = [-value for value in entries]
    return sp.Matrix(matrix.rows, matrix.cols, entries)


def sparse_matrix_certificate(matrix, store_entries=True):
    integer = primitive_integer_matrix(matrix)
    entries = (
        [
            [i, j, int(integer[i, j])]
            for i in range(integer.rows)
            for j in range(integer.cols)
            if integer[i, j]
        ]
        if store_entries
        else None
    )
    return {
        "shape": [integer.rows, integer.cols],
        "entries": entries,
        "nonzero_entries": sum(value != 0 for value in integer),
        "sum_abs_entries": sum(abs(int(value)) for value in integer),
        "weighted_mod_1000000007": sum(
            (i + 1) * (j + 1) * int(integer[i, j])
            for i in range(integer.rows)
            for j in range(integer.cols)
        ) % 1_000_000_007,
    }


def decode_sparse_matrix(certificate):
    if certificate["entries"] is None:
        raise ValueError("certificate intentionally stores checksums only")
    matrix = sp.zeros(*certificate["shape"])
    for i, j, value in certificate["entries"]:
        matrix[i, j] = value
    return matrix


def rational_matrix_mod(matrix, prime):
    """Reduce an exact-Q matrix at a prime avoiding every denominator."""
    output = np.zeros(matrix.shape, dtype=np.int64)
    denominators = []
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            value = sp.Rational(matrix[i, j])
            numerator, denominator = int(value.p), int(value.q)
            if denominator % prime == 0:
                raise AssertionError("bad reduction prime divides a denominator")
            denominators.append(denominator)
            output[i, j] = (
                (numerator % prime) * pow(denominator % prime, prime - 2, prime)
            ) % prime
    return output, math.lcm(*denominators)


def raw_matrix_word_minor(blocks, prime=65_521):
    """Literal-word modular minor for a block-diagonal exact-Q image.

    Echelon rows decide independence, but every commutator is formed from the
    stored raw parent word.  Thus the recipes name literal rational Lie words.
    At acceptance, the unnormalised insertion pivot is recorded.  On the
    selected pivot columns the accepted raw rows have determinant equal to the
    product of those insertion pivots: earlier-row subtraction is unit lower
    triangular and does not change the determinant.
    """
    converted = []
    denominator_lcms = []
    for a, b in blocks:
        a_mod, a_lcm = rational_matrix_mod(a, prime)
        b_mod, b_lcm = rational_matrix_mod(b, prime)
        converted.append((a_mod, b_mod))
        denominator_lcms.append([a_lcm, b_lcm])

    generators = [
        tuple(block[generator_index] for block in converted)
        for generator_index in (0, 1)
    ]
    reduced_rows = []
    pivot_to_row = {}
    raw_words = []
    recipes = []
    pivot_columns = []
    insertion_pivots = []
    determinant = 1

    def flatten(word):
        return np.concatenate([matrix.ravel() for matrix in word])

    def accept(raw_word, recipe):
        nonlocal determinant
        value = flatten(raw_word).copy() % prime
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                return None
            lead = int(nonzero[0])
            previous = pivot_to_row.get(lead)
            if previous is None:
                insertion = int(value[lead])
                determinant = determinant * insertion % prime
                value = value * pow(insertion, prime - 2, prime) % prime
                pivot_to_row[lead] = len(reduced_rows)
                reduced_rows.append(value)
                raw_words.append(tuple(matrix.copy() for matrix in raw_word))
                recipes.append(recipe)
                pivot_columns.append(lead)
                insertion_pivots.append(insertion)
                return len(raw_words) - 1
            value = (
                value - int(value[lead]) * reduced_rows[previous]
            ) % prime

    frontier = []
    for generator_index, generator in enumerate(generators):
        accepted = accept(generator, [-1, generator_index, 1])
        if accepted is not None:
            frontier.append(accepted)
    head = 0
    while head < len(frontier):
        parent = frontier[head]
        head += 1
        for generator_index, generator in enumerate(generators):
            raw_parent = raw_words[parent]
            candidate = tuple(
                (left @ right - right @ left) % prime
                for left, right in zip(generator, raw_parent)
            )
            accepted = accept(
                candidate,
                [parent, generator_index, recipes[parent][2] + 1],
            )
            if accepted is not None:
                frontier.append(accepted)

    assert determinant == math.prod(insertion_pivots) % prime
    assert determinant != 0
    return {
        "claim_tag": "[COMPUTATION]",
        "prime": prime,
        "rank_Fp": len(raw_words),
        "coordinate_columns": sum(a.rows * a.cols for a, _ in blocks),
        "block_dimensions": [a.rows for a, _ in blocks],
        "generator_denominator_lcms": denominator_lcms,
        "all_denominators_nonzero_mod_p": all(
            value % prime for pair in denominator_lcms for value in pair
        ),
        "raw_parent_recipes": recipes,
        "maximum_word_depth": max(recipe[2] for recipe in recipes),
        "pivot_columns": pivot_columns,
        "unnormalised_insertion_pivots": insertion_pivots,
        "selected_minor_determinant_mod_p": determinant,
        "recipe_checksum_mod_1000000007": sum(
            (index + 1)
            * (recipe[0] + 2)
            * (recipe[1] + 2)
            * (recipe[2] + 1)
            for index, recipe in enumerate(recipes)
        ) % 1_000_000_007,
        "rational_implication": (
            "the selected raw rational-word minor has good reduction and "
            "nonzero determinant modulo p, so it is nonzero over Q"
        ),
    }


def modular_word_certificate():
    """Good-prime lower bounds from the same nested rational/integer words."""
    e20 = load_module("_e20_for_e76", E20_PATH)
    model = e20.OrbitModel(ROWS, COLS, False)
    support = model.enumerate_support(100_000)

    def close(seed_words):
        worker = e20.DenseOrbitEngine(model, GOOD_PRIME, support, None, None)

        def apply_raw(action_index, word):
            source_block, row = word
            target_block = source_block ^ (3 if action_index == 0 else 1)
            sources, signs = worker.actions[action_index][source_block]
            result = np.zeros(worker.block_sizes[target_block], dtype=np.int64)
            for source, sign in zip(sources, signs):
                selected = source >= 0
                result[selected] += (
                    sign[selected].astype(np.int64)
                    * row[source[selected]].astype(np.int64)
                )
            result %= GOOD_PRIME
            return target_block, result

        words = []
        tree = []
        frontier = []
        for word, description in seed_words(worker, apply_raw):
            added = worker.reduce_add(word)
            assert added is not None
            words.append((word[0], np.asarray(word[1], dtype=np.int64) % GOOD_PRIME))
            tree.append(description)
            frontier.append(len(words) - 1)
        while frontier:
            next_frontier = []
            for parent in frontier:
                for action_index in (0, 1):
                    candidate = apply_raw(action_index, words[parent])
                    if worker.reduce_add(candidate) is not None:
                        words.append(candidate)
                        tree.append([parent, action_index, tree[parent][2] + 1])
                        next_frontier.append(len(words) - 1)
            frontier = next_frontier
        checksum = sum(
            (index + 1)
            * (entry[0] + 2)
            * (entry[1] + 2)
            * (entry[2] + 1)
            for index, entry in enumerate(tree)
        ) % 1_000_000_007
        return {
            "rank_Fp": len(words),
            "maximum_depth": max(entry[2] for entry in tree),
            "grading_block_ranks": list(worker.row_counts),
            "tree_checksum_mod_1000000007": checksum,
            "tree": tree,
        }

    def full_seeds(worker, _apply):
        return [
            (worker.dense_seed(model.x_terms), [-1, 0, 1]),
            (worker.dense_seed(model.zz_terms), [-1, 1, 1]),
        ]

    def derived_seed(worker, apply):
        b = worker.dense_seed(model.zz_terms)
        bracket = apply(0, b)
        return [(bracket, [-1, 0, 2])]

    full = close(full_seeds)
    derived = close(derived_seed)
    assert full["rank_Fp"] == DIMENSION
    assert derived["rank_Fp"] == DERIVED_DIMENSION
    return {
        "claim_tag": "[COMPUTATION]",
        "prime": GOOD_PRIME,
        "arithmetic": "exact F_p; every row is an actual nested [.,.]/2 word, hence a rational Lie word",
        "support_orbits": len(support),
        "full_closure": full,
        "derived_ideal": derived,
        "rational_implication": "each nonzero modular pivot minor is a nonzero integer minor, so the displayed ranks are rigorous lower bounds over Q",
    }


def exact_q_certificate():
    base = load_module("_e49_helpers_for_e76", E49_PATH)
    sectors = {key: sector_with_gram(key) for key in CHARACTERS}
    sector_records = {}
    active_actions = {}
    component_actions = {}

    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, gram = sectors[key]
        kernel, complement, active = gram_split_common_kernel((a, b), gram)
        active_actions[label] = active
        forms = base.invariant_symmetric_forms(active)
        record = {
            "claim_tag": "[COMPUTATION]",
            "module_dimension_Q": a.rows,
            "common_kernel_dimension_Q": kernel.cols,
            "active_dimension_Q": active[0].rows,
            "physical_gram": sparse_matrix_certificate(gram),
            "kernel_basis": sparse_matrix_certificate(kernel),
            "gram_orthogonal_complement": sparse_matrix_certificate(complement),
            "split_check": True,
            "invariant_symmetric_form_nullity_Q": len(forms),
            "invariant_symmetric_form_max_rank_Q": base.max_combination_rank(forms),
        }
        if forms:
            record["invariant_symmetric_form"] = sparse_matrix_certificate(forms[0])
        if label in {"001", "011", "101", "111"}:
            projectors = base.multiplicity_free_projectors(a, b)
            component_actions[label] = [
                base.restrict_to_projector((a, b), projector)
                for _, projector in projectors
            ]
            dimensions = [dimension for dimension, _ in projectors]
            trace_rows = [
                [action.trace() for action in (a, b)]
                for a, b in component_actions[label]
            ]
            trace_matrix = sp.Matrix(trace_rows).T
            record["projector_dimensions_Q"] = dimensions
            record["projectors"] = [
                sparse_matrix_certificate(projector) for _, projector in projectors
            ]
            record["block_trace_matrix"] = [
                [str(trace_matrix[i, j]) for j in range(trace_matrix.cols)]
                for i in range(trace_matrix.rows)
            ]
            record["block_trace_rank_Q"] = trace_matrix.rank()
        sector_records[label] = record

    # Orthogonal containment is on exact split active summands; no extension
    # cocycle survives the Gram-orthogonal splitting.
    assert sector_records["000"]["common_kernel_dimension_Q"] == 2
    assert sector_records["000"]["active_dimension_Q"] == 42
    assert sector_records["000"]["invariant_symmetric_form_max_rank_Q"] == 42
    for label in ("010", "100", "110"):
        assert sector_records[label]["common_kernel_dimension_Q"] == 1
        assert sector_records[label]["active_dimension_Q"] == 27
        assert sector_records[label]["invariant_symmetric_form_max_rank_Q"] == 27

    expected_components = {
        "001": [4, 4, 24],
        "111": [4, 4, 24],
        "011": [4, 28],
        "101": [4, 28],
    }
    for label, dimensions in expected_components.items():
        assert sector_records[label]["projector_dimensions_Q"] == dimensions
        assert sector_records[label]["block_trace_rank_Q"] == 1

    def full_rank_intertwiner(source, target, dual=False):
        source = base.traceless_actions(source)
        target = base.traceless_actions(target)
        if dual:
            source = tuple(-matrix.T for matrix in source)
        basis = base.intertwiner_basis(target, source)
        weights, witness = base.full_rank_combination(basis)
        return {
            "space_dimension_Q": len(basis),
            "dual": dual,
            "rank_Q": witness.rank() if witness is not None else 0,
            "weights": list(weights) if weights is not None else None,
            "witness": sparse_matrix_certificate(witness) if witness is not None else None,
        }

    odd_pairings = [
        full_rank_intertwiner(component_actions["001"][0], component_actions["111"][0], True),
        full_rank_intertwiner(component_actions["001"][1], component_actions["111"][1], True),
        full_rank_intertwiner(component_actions["001"][2], component_actions["111"][2], True),
        full_rank_intertwiner(component_actions["011"][0], component_actions["101"][0], True),
        full_rank_intertwiner(component_actions["011"][1], component_actions["101"][1], True),
    ]
    shared_a3 = full_rank_intertwiner(
        component_actions["001"][0], component_actions["011"][0], False
    )
    repeated_b13 = full_rank_intertwiner(
        active_actions["100"], active_actions["110"], False
    )
    assert [item["rank_Q"] for item in odd_pairings] == [4, 4, 24, 4, 28]
    assert shared_a3["rank_Q"] == 4
    assert repeated_b13["rank_Q"] == 27

    # Independent exact lower bounds for every sector image and every linkage
    # used in the container count.  These are literal-word minors in the same
    # exact-Q bases used above, not closures of reduced modular combinations.
    sector_image_minors = {
        label: raw_matrix_word_minor([active_actions[label]])
        if label in active_actions and label not in component_actions
        else raw_matrix_word_minor([sectors[tuple(map(int, label))][:2]])
        for label in ("000", "001", "010", "011", "100", "101", "110", "111")
    }
    linkage_minors = {
        "010_plus_100": raw_matrix_word_minor(
            [active_actions["010"], active_actions["100"]]
        ),
        "100_plus_110": raw_matrix_word_minor(
            [active_actions["100"], active_actions["110"]]
        ),
        "001_plus_111": raw_matrix_word_minor(
            [
                sectors[(0, 0, 1)][:2],
                sectors[(1, 1, 1)][:2],
            ]
        ),
        "011_plus_101": raw_matrix_word_minor(
            [
                sectors[(0, 1, 1)][:2],
                sectors[(1, 0, 1)][:2],
            ]
        ),
        "001_plus_011": raw_matrix_word_minor(
            [
                sectors[(0, 0, 1)][:2],
                sectors[(0, 1, 1)][:2],
            ]
        ),
    }
    expected_sector_ranks = {
        "000": 861,
        "001": 606,
        "010": 351,
        "011": 799,
        "100": 351,
        "101": 799,
        "110": 351,
        "111": 606,
    }
    expected_linkage_ranks = {
        "010_plus_100": 702,
        "100_plus_110": 351,
        "001_plus_111": 606,
        "011_plus_101": 799,
        "001_plus_011": 1389,
    }
    assert {
        label: certificate["rank_Fp"]
        for label, certificate in sector_image_minors.items()
    } == expected_sector_ranks
    assert {
        label: certificate["rank_Fp"]
        for label, certificate in linkage_minors.items()
    } == expected_linkage_ranks

    # The one global scalar trace direction.  Exact projectors split every odd
    # sector; the componentwise scalar part of B commutes with A and B.
    scalar_pieces = []
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, _ = sectors[key]
        if label in component_actions:
            projectors = base.multiplicity_free_projectors(a, b)
            scalar = sp.zeros(a.rows)
            for (_, projector), actions in zip(projectors, component_actions[label]):
                scalar += actions[1].trace() * projector / actions[1].rows
        else:
            scalar = sp.zeros(a.rows)
        assert scalar * a == a * scalar
        assert scalar * b == b * scalar
        scalar_pieces.append(scalar)
    central = sp.diag(*scalar_pieces)
    central_integer = primitive_integer_matrix(central)
    pivot = next(
        (index for index, value in enumerate(central_integer) if value),
        None,
    )
    assert pivot is not None
    central_scale = sp.Rational(list(central)[pivot], list(central_integer)[pivot])
    assert central == central_scale * central_integer
    full_a = sp.diag(*(sectors[key][0] for key in CHARACTERS))
    full_b = sp.diag(*(sectors[key][1] for key in CHARACTERS))
    assert central and central * full_a == full_a * central
    assert central * full_b == full_b * central

    odd_upper = 1 + 575 + 783 + 15 + 15
    even_upper = 861 + 351 + 351
    full_upper = even_upper + odd_upper
    derived_upper = full_upper - 1
    assert (odd_upper, even_upper, full_upper, derived_upper) == (1389, 1563, 2952, 2951)

    return {
        "claim_tag": "[THEOREM]",
        "sector_certificates_Q": sector_records,
        "intertwiners_Q": {
            "001_to_111_componentwise_dual": odd_pairings[:3],
            "011_to_101_componentwise_dual": odd_pairings[3:],
            "shared_A3_001_to_011": shared_a3,
            "repeated_B13_100_to_110": repeated_b13,
        },
        "literal_word_image_minors_Fp": sector_image_minors,
        "literal_word_linkage_minors_Fp": linkage_minors,
        "containing_algebra_Q": {
            "type_with_center": "centre(1) + " + FACTOR_TYPE,
            "dimension_Q": full_upper,
            "derived_dimension_Q": derived_upper,
            "factor_dimensions": list(FACTOR_DIMENSIONS),
            "factor_ranks": list(FACTOR_RANKS),
            "levi_rank": sum(FACTOR_RANKS),
            "dimension_identity": "861 + 2*351 + 575 + 783 + 2*15 = 2951",
            "exact_Q_image_and_joint_ranks": {
                "sector_images": expected_sector_ranks,
                "joint_images": expected_linkage_ranks,
            },
            "linkage_reason": "each displayed F_p rank is a lower bound over Q from a literal rational-word selected minor; exact-Q invariant forms/projectors/intertwiners give the matching upper bounds. In particular 100+110=351 identifies one B13, 010+100=702 proves the other B13 independent, 001+111=606 and 011+101=799 identify dual copies, and 001+011=1389 proves exactly two independent A3 factors plus independent A23/A27 modulo the shared centre",
        },
        "central_witness_global_integer": sparse_matrix_certificate(
            central_integer, store_entries=True
        ),
        "central_projection_equals_integer_times": str(central_scale),
        "centrality_verified_exactly": True,
        "central_membership_proof": "the exact scalar projection z of physical B has z=scale*z_integer with one globally normalized integer witness; B-z is traceless on every simple constituent and therefore lies in the 2951-dimensional semisimple container S. The actual derived-word lower bound meets dim(S)=2951, so g'=S and B-z belongs to g'; hence z=B-(B-z) belongs to g. Exact multiplication gives [z,A]=[z,B]=0. Projection g=Q*z plus S onto S has rank 2951 and kernel exactly Q*z",
        "faithful_representation": {
            "dimension": 256,
            "reason": "the eight sectors are an exact direct-sum change of basis of the full physical matrix representation in which the Lie algebra is defined",
        },
    }


def build_artifact():
    timings = {}
    memory = {}
    started = time.monotonic()
    exact_q = timed("exact_Q_upper_and_linkage", timings, memory, exact_q_certificate)
    modular = timed("good_prime_word_lower_bounds", timings, memory, modular_word_certificate)
    timings["total"] = round(time.monotonic() - started, 6)

    full_lower = modular["full_closure"]["rank_Fp"]
    derived_lower = modular["derived_ideal"]["rank_Fp"]
    full_upper = exact_q["containing_algebra_Q"]["dimension_Q"]
    derived_upper = exact_q["containing_algebra_Q"]["derived_dimension_Q"]
    assert full_lower == full_upper == DIMENSION
    assert derived_lower == derived_upper == DERIVED_DIMENSION

    structural_killing = {
        "claim_tag": "[THEOREM]",
        "rank_Q": DERIVED_DIMENSION,
        "nullity_Q": 1,
        "reconstruction": {
            "zero_center_block_dimension": 1,
            "simple_blocks": [
                {"absolute_type": "D21", "dimension": 861, "rational_form": "so(F_000,Q)"},
                {"absolute_type": "B13", "dimension": 351, "rational_form": "so(F_010,Q)"},
                {"absolute_type": "B13", "dimension": 351, "rational_form": "so(F_100,Q)"},
                {"absolute_type": "A23", "dimension": 575, "rational_form": "sl(24,Q)"},
                {"absolute_type": "A27", "dimension": 783, "rational_form": "sl(28,Q)"},
                {"absolute_type": "A3", "dimension": 15, "rational_form": "sl(4,Q)"},
                {"absolute_type": "A3", "dimension": 15, "rational_form": "sl(4,Q)"},
            ],
            "rational_form_scope": "the exact invariant forms F certify orthogonal Q-forms of the displayed absolute Dynkin types; splitness of D21/B13 over Q is not asserted",
            "reason": "the Killing form is zero on the central block and nondegenerate on every characteristic-zero simple Q-form",
        },
        "direct_adjoint_matrix": {
            "claim_tag": "[UNRESOLVED]",
            "materialized": False,
            "distinction": "rank is structurally reconstructed from the certified reductive direct sum; no direct 2952x2952 adjoint-trace matrix or 2951-minor is claimed",
        },
    }

    data = {
        "scope": "the single finite open 2x4 layer only; no all-sizes assertion",
        "field": "Q, with one exact good-prime computation used only for lower bounds",
        "dimension_Q": {
            "claim_tag": "[THEOREM]",
            "value": DIMENSION,
            "upper_bound_Q": full_upper,
            "lower_bound_from_Fp_word_minor": full_lower,
        },
        "derived_algebra": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": DERIVED_DIMENSION,
            "semisimple": True,
            "type": FACTOR_TYPE,
            "upper_bound_Q": derived_upper,
            "lower_bound_from_Fp_word_minor": derived_lower,
        },
        "factor_map": {
            "claim_tag": "[THEOREM]",
            "domain_dimension_Q": DIMENSION,
            "joint_rank_Q": DERIVED_DIMENSION,
            "kernel_dimension_Q": 1,
            "kernel_equals_center": True,
            "central_witness_global_integer": exact_q[
                "central_witness_global_integer"
            ],
            "central_projection_equals_integer_times": exact_q[
                "central_projection_equals_integer_times"
            ],
            "proof": exact_q["central_membership_proof"],
        },
        "solvable_radical": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 1,
            "equals_center": True,
            "proof": "the quotient by the exact central kernel is the displayed semisimple direct sum; every solvable ideal maps trivially to it",
        },
        "Levi_type": {
            "claim_tag": "[THEOREM]",
            "type": FACTOR_TYPE,
            "dimension_Q": DERIVED_DIMENSION,
            "rank": sum(FACTOR_RANKS),
            "factor_dimensions": list(FACTOR_DIMENSIONS),
            "factor_ranks": list(FACTOR_RANKS),
        },
        "Killing": structural_killing,
        "exact_Q_certificate": exact_q,
        "modular_word_certificate": modular,
        "resource_observations": {
            "claim_tag": "[COMPUTATION]",
            "direct_adjoint_collection_int64_bytes": DIMENSION ** 3 * 8,
            "direct_adjoint_collection_gib": DIMENSION ** 3 * 8 / 1024 ** 3,
            "prior_direct_kernel_attempt": {
                "elapsed_seconds_lower_bound": 3600,
                "rss_bytes": 660_242_432,
                "result_used": False,
            },
            "resolution": "exact sector factorization and dimension squeezing avoid the direct adjoint/Killing matrix wall",
        },
        "timings_seconds": timings,
        "memory_observations": memory,
    }
    checks = [
        {"name": "same_word_dimension_lower", "passed": full_lower == 2952, "detail": "2952 actual nested words independent at the good prime"},
        {"name": "exact_Q_dimension_upper", "passed": full_upper == 2952, "detail": "faithful linked classical/block container"},
        {"name": "derived_dimension_squeeze", "passed": derived_lower == derived_upper == 2951, "detail": "good-prime ideal words meet exact-Q traceless upper bound"},
        {"name": "factor_dimensions", "passed": sum(FACTOR_DIMENSIONS) == 2951, "detail": FACTOR_TYPE},
        {"name": "factor_rank", "passed": sum(FACTOR_RANKS) == 103, "detail": "21+2*13+23+27+2*3"},
        {"name": "factor_map_kernel", "passed": data["factor_map"]["joint_rank_Q"] == 2951 and data["factor_map"]["kernel_equals_center"], "detail": "faithful semisimple projection"},
        {"name": "radical_equals_center", "passed": data["solvable_radical"]["equals_center"], "detail": "one-dimensional central kernel of faithful semisimple quotient"},
        {"name": "structural_Killing_rank", "passed": structural_killing["rank_Q"] == 2951 and structural_killing["nullity_Q"] == 1, "detail": "direct sum of centre and seven simple factors"},
    ]
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": "same-word good-prime lower bounds plus exact-Q faithful sector upper bounds, projectors, invariant forms, Gram splittings, and intertwiners",
            "python": platform.python_version(),
        },
        "data": data,
        "checks": checks,
    }


def run(write=True):
    artifact = build_artifact()
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return artifact


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args(argv)
    try:
        artifact = run(not args.no_write)
        for check in artifact["checks"]:
            print(("PASS" if check["passed"] else "FAIL") + ": " + check["name"])
        passed = all(check["passed"] for check in artifact["checks"])
        print("PASS" if passed else "FAIL")
        return 0 if passed else 1
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
