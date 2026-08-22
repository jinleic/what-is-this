#!/usr/bin/env python3
"""Exact finite shadows for the all-size local commutant theorem on Z^3.

The theorem itself is proved in proofs/local_commutant.md.  This producer does
not promote finite calculations into that proof.  It instead checks the two
local coefficient identities used by the induction, exhausts the private-
neighbor geometry in small boxes, and reconstructs exact ordered-Pauli kernels
for representative compact supports with every Z^3 boundary bond retained.
"""
from __future__ import annotations

import hashlib
import json
import platform
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "local_commutant.json"
E52_RESULT = ROOT / "results" / "integrability" / "full_commutant.json"
SCRIPT = "experiments/e115_local_commutant.py"
PRIMES = (2_147_483_647, 2_147_483_629)
COMBINATIONS = ((1, 1), (1, 2), (2, 1), (3, 5), (1, -1), (2, 3))
BOXES = ((2, 2, 2), (3, 2, 2), (2, 2, 3), (3, 3, 1), (3, 3, 2))


# --------------------------------------------------------------------------- #
# Z^3 geometry                                                                  #
# --------------------------------------------------------------------------- #
Site = tuple[int, int, int]
PauliKey = tuple[int, int]  # (X-mask, Z-mask) over one explicit site legend
Image = dict[PauliKey, int | Fraction]


def neighbors(site: Site) -> tuple[Site, ...]:
    x, y, z = site
    return (
        (x + 1, y, z), (x - 1, y, z),
        (x, y + 1, z), (x, y - 1, z),
        (x, y, z + 1), (x, y, z - 1),
    )


def box_sites(shape: tuple[int, int, int]) -> list[Site]:
    a, b, c = shape
    return sorted((x, y, z) for x in range(a) for y in range(b) for z in range(c))


def canonical_peel_site(sites: list[Site]) -> Site:
    """A deterministic maximal-first-coordinate choice; every maximizer works."""
    return max(sorted(sites))


def private_neighbor_census(shape: tuple[int, int, int]) -> dict[str, object]:
    sites = box_sites(shape)
    n = len(sites)
    nonempty = maximizer_checks = maximizer_failures = nonmax_failure_subsets = 0
    witnesses: list[dict[str, object]] = []
    for mask in range(1, 1 << n):
        support = {sites[i] for i in range(n) if (mask >> i) & 1}
        nonempty += 1
        max_first = max(site[0] for site in support)
        for v in support:
            w = (v[0] + 1, v[1], v[2])
            hit = support.intersection(neighbors(w))
            valid = w not in support and hit == {v}
            if v[0] == max_first:
                maximizer_checks += 1
                maximizer_failures += int(not valid)
            elif not valid and len(witnesses) < 3:
                witnesses.append(
                    {
                        "subset": [list(site) for site in sorted(support)],
                        "u": list(v),
                        "u_plus_e1": list(w),
                        "support_neighbors_of_u_plus_e1": [list(site) for site in sorted(hit)],
                    }
                )
        if any(
            u[0] < max_first
            and (
                (u[0] + 1, u[1], u[2]) in support
                or support.intersection(neighbors((u[0] + 1, u[1], u[2]))) != {u}
            )
            for u in support
        ):
            nonmax_failure_subsets += 1
    return {
        "shape": list(shape),
        "n_sites": n,
        "nonempty_subsets": nonempty,
        "maximizer_checks": maximizer_checks,
        "maximizer_failures": maximizer_failures,
        "nonmaximizer_failure_subsets": nonmax_failure_subsets,
        "example_witnesses": witnesses,
    }


# --------------------------------------------------------------------------- #
# Canonical ordered-Pauli maps: exactly the e52/e61 convention                  #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Bond:
    left: Site
    right: Site
    internal: bool


@dataclass
class SupportModel:
    """Finite representation of derivations on Alg(S), with all Z^3 boundary bonds."""

    sites: list[Site]
    ext: list[Site]
    index: dict[Site, int]
    bonds: list[Bond]
    all_bonds: list[Bond]
    columns: list[PauliKey]
    x_images: list[Image]
    z_images: list[Image]

    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def n_ext(self) -> int:
        return len(self.ext)

    @property
    def n_columns(self) -> int:
        return len(self.columns)


def add_term(image: Image, key: PauliKey, value: int | Fraction) -> None:
    if not value:
        return
    next_value = image.get(key, 0) + value
    if next_value:
        image[key] = next_value
    else:
        image.pop(key, None)


def make_model(support_sites: list[Site], internal_only: bool = False) -> SupportModel:
    sites = sorted(support_sites)
    site_set = set(sites)
    # The legend contains S followed by every neighbour reached from S, in a
    # deterministic order.  It is only a finite bookkeeping device for δ_B(Q).
    ext = list(dict.fromkeys(sites + [neighbor for site in sites for neighbor in neighbors(site)]))
    index = {site: i for i, site in enumerate(ext)}
    seen: set[tuple[Site, Site]] = set()
    all_bonds: list[Bond] = []
    for site in sites:
        for neighbor in neighbors(site):
            pair = tuple(sorted((site, neighbor)))
            if pair not in seen:
                seen.add(pair)
                all_bonds.append(Bond(pair[0], pair[1], pair[0] in site_set and pair[1] in site_set))
    bonds = [bond for bond in all_bonds if bond.internal] if internal_only else list(all_bonds)

    columns: list[PauliKey] = []
    for column in range(1 << (2 * len(sites))):
        digits = column
        a_mask = b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            bit = 1 << index[site]
            if label & 1:
                a_mask |= bit
            if label & 2:
                b_mask |= bit
        columns.append((a_mask, b_mask))

    # Repository convention (generator on the LEFT):
    #   1/2 [X_s, Q_(a,b)] = [b_s] Q_(a+e_s,b),
    #   1/2 [Z_u Z_v, Q_(a,b)] = -[a_u xor a_v] Q_(a,b+e_u+e_v).
    x_images: list[Image] = []
    z_images: list[Image] = []
    for a_mask, b_mask in columns:
        x_image: Image = {}
        z_image: Image = {}
        for site in sites:
            bit = 1 << index[site]
            if b_mask & bit:
                add_term(x_image, (a_mask ^ bit, b_mask), 1)
        for bond in bonds:
            left_bit = 1 << index[bond.left]
            right_bit = 1 << index[bond.right]
            if bool(a_mask & left_bit) ^ bool(a_mask & right_bit):
                add_term(z_image, (a_mask, b_mask ^ left_bit ^ right_bit), -1)
        x_images.append(x_image)
        z_images.append(z_image)
    return SupportModel(sites, ext, index, bonds, all_bonds, columns, x_images, z_images)


def linear_combination(x_image: Image, z_image: Image, a: Fraction, b: Fraction) -> Image:
    image: Image = {}
    for key, value in x_image.items():
        add_term(image, key, a * value)
    for key, value in z_image.items():
        add_term(image, key, b * value)
    return image


def pencil_images(model: SupportModel, a: int, b: int) -> list[Image]:
    aq = Fraction(a)
    bq = Fraction(b)
    return [linear_combination(x_image, z_image, aq, bq) for x_image, z_image in zip(model.x_images, model.z_images)]


def simultaneous_images(model: SupportModel) -> list[dict[tuple[int, int, int], int | Fraction]]:
    """Images of Q -> (1/2[A,Q], 1/2[B,Q]) with a disjoint output tag."""
    images: list[dict[tuple[int, int, int], int | Fraction]] = []
    for x_image, z_image in zip(model.x_images, model.z_images):
        image: dict[tuple[int, int, int], int | Fraction] = {}
        for (a_mask, b_mask), value in x_image.items():
            image[(0, a_mask, b_mask)] = value
        for (a_mask, b_mask), value in z_image.items():
            image[(1, a_mask, b_mask)] = value
        images.append(image)
    return images


# --------------------------------------------------------------------------- #
# exact / modular ranks of the image vectors                                   #
# --------------------------------------------------------------------------- #
def rank_q(images: list[dict[object, int | Fraction]]) -> tuple[int, dict[str, int]]:
    """Exact-Q sparse row reduction of image vectors, frequency ordered for low fill."""
    frequency = Counter(key for image in images for key in image)
    order = {
        key: position
        for position, key in enumerate(sorted(frequency, key=lambda key: (frequency[key], key)))
    }
    pivots: dict[object, dict[object, Fraction]] = {}
    max_support = 0
    max_bits = 0
    for source in sorted(images, key=len):
        row = {key: Fraction(value) for key, value in source.items() if value}
        while row:
            lead = min(row, key=order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                coefficient = row[lead]
                if coefficient != 1:
                    row = {key: value / coefficient for key, value in row.items()}
                pivots[lead] = row
                max_support = max(max_support, len(row))
                for value in row.values():
                    max_bits = max(max_bits, value.numerator.bit_length(), value.denominator.bit_length())
                break
            coefficient = row[lead]
            for key, value in old.items():
                reduced = row.get(key, Fraction(0)) - coefficient * value
                if reduced:
                    row[key] = reduced
                else:
                    row.pop(key, None)
    return len(pivots), {
        "rank": len(pivots),
        "maximum_pivot_row_support": max_support,
        "maximum_coefficient_bits": max_bits,
    }


def rank_mod(images: list[dict[object, int | Fraction]], prime: int) -> tuple[int, dict[str, int]]:
    """Independent prime-field cross-check; nullity over Q cannot exceed it."""
    frequency = Counter(key for image in images for key in image)
    order = {
        key: position
        for position, key in enumerate(sorted(frequency, key=lambda key: (frequency[key], key)))
    }
    pivots: dict[object, dict[object, int]] = {}
    max_support = 0
    for source in sorted(images, key=len):
        row: dict[object, int] = {}
        for key, value in source.items():
            if isinstance(value, Fraction):
                residue = value.numerator * pow(value.denominator, prime - 2, prime) % prime
            else:
                residue = value % prime
            if residue:
                row[key] = residue
        while row:
            lead = min(row, key=order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(row[lead], prime - 2, prime)
                pivots[lead] = {key: value * inverse % prime for key, value in row.items()}
                max_support = max(max_support, len(row))
                break
            coefficient = row[lead]
            for key, value in old.items():
                reduced = (row.get(key, 0) - coefficient * value) % prime
                if reduced:
                    row[key] = reduced
                else:
                    row.pop(key, None)
    return len(pivots), {"rank": len(pivots), "maximum_pivot_row_support": max_support}


# --------------------------------------------------------------------------- #
# tensor-slot bookkeeping and exact peeling maps                               #
# --------------------------------------------------------------------------- #
def slot_extract(model: SupportModel, key: PauliKey, site: Site, alpha: int, beta: int) -> PauliKey | None:
    """Project a key to its (X^alpha Z^beta)_site coefficient, clearing that slot."""
    index = model.index[site]
    a_mask, b_mask = key
    if ((a_mask >> index) & 1, (b_mask >> index) & 1) != (alpha, beta):
        return None
    if alpha:
        a_mask ^= 1 << index
    if beta:
        b_mask ^= 1 << index
    return a_mask, b_mask


def slot_part(model: SupportModel, image: Image, site: Site, alpha: int, beta: int) -> Image:
    output: Image = {}
    for key, value in image.items():
        stripped = slot_extract(model, key, site, alpha, beta)
        if stripped is not None:
            add_term(output, stripped, value)
    return output


def single_z_image(model: SupportModel, column: int, site: Site) -> Image:
    """(1/2)[Z_site,Q] in the repository orientation, as a one-term map."""
    a_mask, b_mask = model.columns[column]
    bit = 1 << model.index[site]
    return {(a_mask, b_mask ^ bit): -1} if a_mask & bit else {}


def slot_map_entries(images: list[Image]) -> list[list[int]]:
    entries: list[list[int]] = []
    for source, image in enumerate(images):
        for (a_mask, b_mask), value in sorted(image.items()):
            if isinstance(value, Fraction):
                if value.denominator != 1:
                    raise ArithmeticError(f"nonintegral local map coefficient {value}")
                value = value.numerator
            entries.append([a_mask, b_mask, source, int(value)])
    return entries


def inspect_peeling(model: SupportModel) -> dict[str, object]:
    support = set(model.sites)
    v = canonical_peel_site(model.sites)
    w = (v[0] + 1, v[1], v[2])
    assert w not in support
    assert support.intersection(neighbors(w)) == {v}
    vi = model.index[v]

    wslot_x: list[Image] = []
    wslot_z: list[Image] = []
    vslot_x: list[Image] = []
    vslot_z: list[Image] = []
    w_no_x_or_y = True
    w_identity_exact = True
    v_identity_exact = True
    v_b_zero_on_a_zero = True
    v_a_matches_z_on_a_zero = True
    outside_w_bonds_vanish = True

    for column, (a_mask, b_mask) in enumerate(model.columns):
        x_image = model.x_images[column]
        z_image = model.z_images[column]
        wx = slot_part(model, x_image, w, 0, 1)
        wz = slot_part(model, z_image, w, 0, 1)
        vx = slot_part(model, x_image, v, 1, 1)
        vz = slot_part(model, z_image, v, 1, 1)
        wslot_x.append(wx)
        wslot_z.append(wz)
        vslot_x.append(vx)
        vslot_z.append(vz)
        w_identity_exact &= not wx and wz == single_z_image(model, column, v)
        # The exact Y-slot identity used only after [Q,Z_v]=0, i.e. a_v=0.
        if not (a_mask & (1 << vi)):
            expected = {(a_mask, b_mask ^ (1 << vi)): 1} if b_mask & (1 << vi) else {}
            v_identity_exact &= vx == expected and not vz
            v_b_zero_on_a_zero &= not vz
            v_a_matches_z_on_a_zero &= vx == expected
        for key in list(x_image) + list(z_image):
            a_w = (key[0] >> model.index[w]) & 1
            b_w = (key[1] >> model.index[w]) & 1
            if (a_w, b_w) in ((1, 0), (1, 1)):
                w_no_x_or_y = False

        # These five bonds are present in the formal B but omitted from
        # bonds-meeting-S because both ends lie outside S.  A neighbour u of w
        # need not occur in the finite legend at all; its X-mask is nevertheless
        # zero because every column is supported in S.  Thus a_w xor a_u=0,
        # which is the exact ordered-Pauli criterion for this commutator to
        # vanish.
        for u in neighbors(w):
            if u in support:
                continue
            a_u = bool(a_mask & (1 << model.index[u])) if u in model.index else False
            a_w = bool(a_mask & (1 << model.index[w]))
            if a_u ^ a_w:
                outside_w_bonds_vanish = False

    return {
        "peel_v": list(v),
        "private_w": list(w),
        "support_neighbors_of_w": [list(site) for site in sorted(support.intersection(neighbors(w)))],
        "w_bonds_outside_support": len([u for u in neighbors(w) if u not in support]),
        "w_bonds_outside_support_vanishing": outside_w_bonds_vanish,
        "wslot_identity_exact": w_identity_exact,
        "vslot_identity_exact_after_z_centrality": v_identity_exact,
        "no_Xw_Yw_components": w_no_x_or_y,
        "wslot_x": wslot_x,
        "wslot_z": wslot_z,
        "vslot_x": vslot_x,
        "vslot_z": vslot_z,
        "vslot_b_zero_on_a_v_zero_columns": v_b_zero_on_a_zero,
        "vslot_a_matches_z_slot_on_a_v_zero_columns": v_a_matches_z_on_a_zero,
    }


def local_map_record(sites: list[Site]) -> dict[str, object]:
    model = make_model(sites)
    peeling = inspect_peeling(model)
    return {
        "sites": [list(site) for site in model.sites],
        "sites_legend": [list(site) for site in model.ext],
        "peel_v": peeling["peel_v"],
        "private_w": peeling["private_w"],
        "domain_columns": model.n_columns,
        "wslot_pencil": {
            "a_entries": slot_map_entries(peeling["wslot_x"]),
            "b_entries": slot_map_entries(peeling["wslot_z"]),
        },
        "vslot_pencil": {
            "a_entries": slot_map_entries(peeling["vslot_x"]),
            "b_entries": slot_map_entries(peeling["vslot_z"]),
        },
        "structural_facts": {
            "wslot_a_part_zero": all(not image for image in peeling["wslot_x"]),
            "wslot_b_part_equals_half_ad_Zv": peeling["wslot_identity_exact"],
            "vslot_b_part_zero_on_a_v_zero_columns": peeling["vslot_b_zero_on_a_v_zero_columns"],
            "vslot_a_part_matches_z_slot_on_a_v_zero_columns": peeling["vslot_a_matches_z_slot_on_a_v_zero_columns"],
            "vslot_identity_after_z_centrality": peeling["vslot_identity_exact_after_z_centrality"],
            "vslot_b_part_nonzero_in_general": any(peeling["vslot_z"]),
        },
    }


# --------------------------------------------------------------------------- #
# representative exact kernels                                                 #
# --------------------------------------------------------------------------- #
def analyze_support(
    name: str,
    sites: list[Site],
    internal_exact_nullity: int | None = None,
    internal_graph_label: str | None = None,
) -> dict[str, object]:
    model = make_model(sites)
    peeling = inspect_peeling(model)
    combinations: list[dict[str, object]] = []
    for a, b in COMBINATIONS:
        images = pencil_images(model, a, b)
        rank, diagnostics = rank_q(images)
        nullity = model.n_columns - rank
        assert not images[0], "identity must be in every H_g kernel"
        assert nullity == 1, (name, a, b, rank, model.n_columns)
        combination: dict[str, object] = {
            "a": a,
            "b": b,
            "rank_Q": rank,
            "nullity_Q": nullity,
            "rank_diagnostics": diagnostics,
            "identity_spans_kernel": rank == model.n_columns - 1 and not images[0],
        }
        if (a, b) == (1, 1):
            modular = []
            for prime in PRIMES:
                mod_rank, mod_diagnostics = rank_mod(images, prime)
                modular.append({"prime": prime, "rank_Fp": mod_rank, "nullity_upper_bound_over_Q": model.n_columns - mod_rank, **mod_diagnostics})
                assert mod_rank == rank
            combination["modular_cross_checks"] = modular
        combinations.append(combination)

    paired = simultaneous_images(model)
    simultaneous_rank, simultaneous_diagnostics = rank_q(paired)
    simultaneous_nullity = model.n_columns - simultaneous_rank
    assert simultaneous_nullity == 1, (name, simultaneous_rank, model.n_columns)
    simultaneous_modular = []
    for prime in PRIMES:
        mod_rank, mod_diagnostics = rank_mod(paired, prime)
        simultaneous_modular.append({"prime": prime, "rank_Fp": mod_rank, "nullity_upper_bound_over_Q": model.n_columns - mod_rank, **mod_diagnostics})
        assert mod_rank == simultaneous_rank

    record: dict[str, object] = {
        "name": name,
        "sites": [list(site) for site in model.sites],
        "n_sites": model.n_sites,
        "columns": model.n_columns,
        "bonds_meeting_support": len(model.all_bonds),
        "internal_bonds": sum(bond.internal for bond in model.all_bonds),
        "boundary_bonds": sum(not bond.internal for bond in model.all_bonds),
        "peel_v": peeling["peel_v"],
        "private_w": peeling["private_w"],
        "support_neighbors_of_w": peeling["support_neighbors_of_w"],
        "w_bonds_outside_support": peeling["w_bonds_outside_support"],
        "w_bonds_outside_support_vanishing": peeling["w_bonds_outside_support_vanishing"],
        "combinations": combinations,
        "simultaneous": {
            "rank_Q": simultaneous_rank,
            "nullity_Q": simultaneous_nullity,
            "rank_diagnostics": simultaneous_diagnostics,
            "modular_cross_checks": simultaneous_modular,
            "identity_spans_kernel": simultaneous_rank == model.n_columns - 1 and not paired[0],
        },
        "peeling_wslot_exact": peeling["wslot_identity_exact"],
        "peeling_vslot_exact_after_z_centrality": peeling["vslot_identity_exact_after_z_centrality"],
        "no_Xw_Yw_components": peeling["no_Xw_Yw_components"],
        "identity_column_in_kernel": not pencil_images(model, 1, 1)[0],
    }

    if internal_exact_nullity is not None:
        internal = make_model(sites, internal_only=True)
        paired_internal = simultaneous_images(internal)
        modular = []
        for prime in PRIMES:
            mod_rank, diagnostics = rank_mod(paired_internal, prime)
            modular.append({"prime": prime, "rank_Fp": mod_rank, "nullity_upper_bound_over_Q": internal.n_columns - mod_rank, **diagnostics})
            assert internal.n_columns - mod_rank == internal_exact_nullity
        record["internal_only"] = {
            "graph": internal_graph_label,
            "internal_bonds": len(internal.bonds),
            "simultaneous_nullity_Q": internal_exact_nullity,
            "certificate_source": "exact-Q e52 finite-graph certificate cross-checked here by two-prime ranks",
            "modular_cross_checks": modular,
            "matches_finite_graph_certificate": True,
        }
    return record


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    started = time.monotonic()

    census = [private_neighbor_census(shape) for shape in BOXES]
    total_subsets = sum(int(row["nonempty_subsets"]) for row in census)
    total_maximizers = sum(int(row["maximizer_checks"]) for row in census)
    total_failures = sum(int(row["maximizer_failures"]) for row in census)

    supports = [
        ("single_site", [(0, 0, 0)], None, None),
        ("bond", [(0, 0, 0), (1, 0, 0)], None, None),
        ("path3", [(0, 0, 0), (1, 0, 0), (2, 0, 0)], None, None),
        ("l_tromino", [(0, 0, 0), (1, 0, 0), (0, 1, 0)], None, None),
        ("claw3d", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], None, None),
        ("plaquette", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)], 27, "2x2_open_grid_is_C4"),
        ("star4", [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)], None, None),
        ("corner3d", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1)], None, None),
        ("slab_2x3", [(x, y, 0) for x in range(2) for y in range(3)], 12, "2x3_open_grid"),
    ]
    records = [analyze_support(*support) for support in supports]
    by_name = {record["name"]: record for record in records}

    e52_bytes = E52_RESULT.read_bytes()
    e52 = json.loads(e52_bytes)
    e52_dims = e52["data"]["finite_exact_dimensions"]
    assert e52_dims["2x2"] == 27 and e52_dims["2x3"] == 12
    assert by_name["plaquette"]["internal_only"]["simultaneous_nullity_Q"] == e52_dims["2x2"]
    assert by_name["slab_2x3"]["internal_only"]["simultaneous_nullity_Q"] == e52_dims["2x3"]

    local_maps = {
        "claim_tag": "[COMPUTATION]",
        "definition": (
            "For each support, the entries are exact integer matrices of the coefficient maps "
            "pi_w,01 o (a ad_A + b ad_B) and pi_v,11 o (a ad_A + b ad_B), "
            "with ad_h=(1/2)[h,.] in the repository orientation.  Each tuple is "
            "[target_X_mask,target_Z_mask,source_column,coefficient]; target masks use the "
            "listed finite site legend with the projected slot cleared."
        ),
        "single_site": local_map_record([(0, 0, 0)]),
        "bond": local_map_record([(0, 0, 0), (1, 0, 0)]),
    }

    internal_summary = {
        "claim_tag": "[COMPUTATION]",
        "statement": (
            "The same compact supports have nontrivial finite-graph simultaneous commutants when "
            "boundary Z^3 bonds are removed, but the full Z^3 derivation has nullity one.  The "
            "finite values are exact-Q e52 certificates; this producer independently reproduces "
            "their nullities modulo both listed primes."
        ),
        "source_e52": "results/integrability/full_commutant.json",
        "source_e52_sha256": hashlib.sha256(e52_bytes).hexdigest(),
        "rows": [
            {
                "name": "plaquette",
                "finite_graph": "2x2_open_grid_is_C4",
                "internal_simultaneous_nullity_Q": 27,
                "full_Z3_simultaneous_nullity_Q": by_name["plaquette"]["simultaneous"]["nullity_Q"],
            },
            {
                "name": "slab_2x3",
                "finite_graph": "2x3_open_grid",
                "internal_simultaneous_nullity_Q": 12,
                "full_Z3_simultaneous_nullity_Q": by_name["slab_2x3"]["simultaneous"]["nullity_Q"],
            },
        ],
    }

    checks = [
        check(
            "private_neighbor_exhaustive_boxes",
            total_failures == 0 and total_subsets == 271_099 and total_maximizers > total_subsets,
            "all nonempty subsets of five Z^3 boxes through 3x3x2 satisfy the maximal-first-coordinate private-neighbor lemma",
        ),
        check(
            "nonmaximal_choice_is_not_vacuous",
            all(row["nonmaximizer_failure_subsets"] > 0 and row["example_witnesses"] for row in census),
            "each exhaustive box contains explicit nonmaximal choices for which the private-neighbor conclusion fails",
        ),
        check(
            "Hg_full_boundary_exact_Q_kernels",
            all(
                all(combo["nullity_Q"] == 1 and combo["identity_spans_kernel"] for combo in record["combinations"])
                for record in records
            ),
            "for every representative support and six nonzero rational coefficient pairs, the exact-Q H_g kernel is exactly span{I}",
        ),
        check(
            "simultaneous_corollary_exact_Q_kernels",
            all(record["simultaneous"]["nullity_Q"] == 1 and record["simultaneous"]["identity_spans_kernel"] for record in records),
            "the tagged (ad_A,ad_B) exact-Q kernel is span{I} on every representative support",
        ),
        check(
            "prime_cross_checks",
            all(
                all(row["rank_Fp"] == combo["rank_Q"] for row in combo.get("modular_cross_checks", []))
                and all(row["rank_Fp"] == record["simultaneous"]["rank_Q"] for row in record["simultaneous"]["modular_cross_checks"])
                for record in records
                for combo in record["combinations"]
                if combo["a"] == 1 and combo["b"] == 1
            ),
            "both primes reproduce the stored exact-Q ranks; modular nullity is only an upper bound over Q",
        ),
        check(
            "wslot_isolation_exact",
            all(record["peeling_wslot_exact"] and record["no_Xw_Yw_components"] for record in records),
            "for every ordered-Pauli column, the Z_w slot is exactly b times the single-site ad_Zv coefficient and has no X_w/Y_w contamination",
        ),
        check(
            "outside_w_bonds_explicitly_vanish",
            all(record["w_bonds_outside_support"] == 5 and record["w_bonds_outside_support_vanishing"] for record in records),
            "all five w-incident bonds whose other endpoint lies outside S have zero commutator on every support-restricted Pauli column",
        ),
        check(
            "vslot_peeling_exact_after_Z_centrality",
            all(record["peeling_vslot_exact_after_z_centrality"] for record in records),
            "after the a_v=1 columns are removed by [Q,Z_v]=0, the (1,1)_v/Y_v coefficient is the field coefficient times the stripped Z_v slot",
        ),
        check(
            "local_pencil_maps_exact",
            all(
                local_maps[name]["structural_facts"]["wslot_a_part_zero"]
                and local_maps[name]["structural_facts"]["wslot_b_part_equals_half_ad_Zv"]
                and local_maps[name]["structural_facts"]["vslot_b_part_zero_on_a_v_zero_columns"]
                and local_maps[name]["structural_facts"]["vslot_identity_after_z_centrality"]
                for name in ("single_site", "bond")
            ),
            "stored sparse integer coefficient maps verify both proof projections, including the necessary a_v=0 restriction for the Y_v step",
        ),
        check(
            "finite_graph_boundary_contrast",
            by_name["plaquette"]["internal_only"]["simultaneous_nullity_Q"] == 27
            and by_name["plaquette"]["simultaneous"]["nullity_Q"] == 1
            and by_name["slab_2x3"]["internal_only"]["simultaneous_nullity_Q"] == 12
            and by_name["slab_2x3"]["simultaneous"]["nullity_Q"] == 1,
            "e52 finite-graph commutants are nontrivial only after omitting the boundary bonds that the Z^3 derivation retains",
        ),
    ]

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "method": (
                "exact ordered-Pauli image ranks over Q; exact local tensor-slot maps; exhaustive finite-subset "
                "geometry; two-prime sparse rank cross-checks; finite-graph contrast read only from e52's exact certificate"
            ),
            "certifying_arithmetic": "integers and fractions.Fraction only; elapsed time is noncertifying",
            "total_elapsed_seconds_noncertifying": round(time.monotonic() - started, 6),
        },
        "data": {
            "claim_tags": ["[THEOREM]", "[LEMMA]", "[COMPUTATION]", "[UNRESOLVED]"],
            "theorem": {
                "claim_tag": "[THEOREM]",
                "name": "compact_local_commutant_for_Hg_on_Z3",
                "statement": (
                    "For every a,b in C with a != 0 and b != 0, every finite-support operator Q in the "
                    "quasi-local spin algebra on Z^3 satisfying [Q,H_g]=0 in the finite-commutator-derivation "
                    "sense is a scalar multiple of I."
                ),
                "H_g": "a sum_(x in Z^3) X_x + b sum_(<xy> nearest-neighbour bonds of Z^3) Z_x Z_y",
                "derivation_definition": (
                    "[Q,H_g] := a sum_x [Q,X_x] + b sum_<xy> [Q,Z_x Z_y]; for finite support S both sums "
                    "are finite because only x in S and bonds meeting S can contribute."
                ),
                "peeling_site_rule": "choose v in supp(Q) of maximal first coordinate and w=v+e_1",
                "corollary": {
                    "claim_tag": "[THEOREM]",
                    "statement": "The simultaneous finite-support centralizer {Q:[Q,A]=[Q,B]=0} is C I.",
                    "reason": "[Q,A]=[Q,B]=0 implies [Q,A+B]=0, so apply the theorem with a=b=1.",
                },
                "field_scope": (
                    "Stated over C.  In ordered X^a Z^b coordinates all maps have integer structure constants, "
                    "so the same argument works over any field of characteristic not 2 with a,b nonzero; "
                    "characteristic 2 is excluded because the commutator coefficients 2 vanish."
                ),
            },
            "conventions": {
                "claim_tag": "[LEMMA]",
                "ordered_pauli_basis": "Q_(a,b)=X^a Z^b, with a,b bit masks over a finite ordered support",
                "column_index": "column=sum_i (a_i + 2 b_i) 4^i, with sites in the listed support order",
                "ad_formulas": (
                    "Repository orientation ad_h=(1/2)[h,.]: ad_Xs(Q_(a,b))=[b_s=1]Q_(a+e_s,b); "
                    "ad_ZuZv(Q_(a,b))=-[a_u xor a_v=1]Q_(a,b+e_u+e_v)."
                ),
                "orientation_remark": (
                    "The proof writes [Q,H_g] (operator Q on the left); [Q,h]=-[h,Q].  The exact JSON maps use "
                    "the established e52/e61 generator-left convention, so their signs are recorded explicitly."
                ),
                "slot_projection": (
                    "For a named site u, pi^(u)_(alpha,beta) is the coefficient map in the fixed two-factor order "
                    "M_2(C)_u tensor A_(Z^3 minus {u}): O=sum_(alpha,beta)(X^alpha Z^beta)_u tensor pi^(u)_(alpha,beta)(O)."
                ),
                "private_neighbor_lemma": (
                    "If v maximizes the first coordinate in finite nonempty S and w=v+e_1, then w notin S and "
                    "S intersect N(w)={v}; every other neighbor of w has first coordinate at least v_1+1, "
                    "strictly above the maximum attained in S."
                ),
            },
            "exhaustive_private_neighbor": {
                "claim_tag": "[COMPUTATION]",
                "boxes": census,
                "totals": {"subsets": total_subsets, "maximizer_checks": total_maximizers, "failures": total_failures},
                "scope": "finite exhaustive illustrations only; the all-size lemma is proved geometrically in the proof note",
            },
            "local_coefficient_maps": local_maps,
            "representative_supports": records,
            "internal_contrast_summary": internal_summary,
            "scope_limits": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "This is a compact finite-support theorem for one nondegenerate H_g.  It is stronger than the "
                    "Pauli-string statement but is not a general no-integrability theorem."
                ),
                "excluded": [
                    "Extensive translation sums of local densities have infinite support and are excluded; this is how 1D integrable charges evade the compact-support theorem.",
                    "Quasilocal, weakly local, or any other infinite-support operators are not classified.",
                    "The finite-torus global spin flip P=product_v X_v has full torus support and is not a compact Z^3 counterexample.",
                    "Pure limits a = 0 or b = 0 are excluded; each one-generator centralizer is large.",
                    "Nonlocal transformations, enlarged Hilbert spaces, and alternative representations are not excluded.",
                    "No conclusion about spectra, transfer-operator diagonalization, or thermodynamic integrability follows solely from this compact-local statement.",
                ],
            },
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    failed = [entry["name"] for entry in checks if not entry["passed"]]
    if failed:
        print("FAIL: " + ", ".join(failed))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
