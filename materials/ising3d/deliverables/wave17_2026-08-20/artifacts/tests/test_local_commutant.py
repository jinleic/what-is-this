#!/usr/bin/env python3
"""Independent standalone verifier for the e115 local-commutant certificate.

Re-derives, with code paths independent of experiments/e115_local_commutant.py:
  * the exhaustive private-neighbour census over finite subsets of small boxes,
  * the exact ordered-Pauli commutator maps of the restricted infinite-volume
    derivations (built from operator PRODUCTS, not from the producer's
    column-scan formulas),
  * representative kernel dimensions by dense exact-Q Gaussian elimination
    (small supports) and by reversed-order modular sparse elimination (all),
  * the two slot-projection peeling identities of the proof, their pencil
    decomposition, and their failure under sabotage (boundary bonds removed,
    or a bogus local-basis cancellation).

Fails hard on any mismatch with results/integrability/local_commutant.json.
"""
from __future__ import annotations

import json
import signal
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "local_commutant.json"
E52 = ROOT / "results" / "integrability" / "full_commutant.json"
PRIME = 2_147_483_647
TIMEOUT_SECONDS = 900


# --------------------------------------------------------------------------- #
# geometry of Z^3                                                              #
# --------------------------------------------------------------------------- #
def neighbors(site):
    x, y, z = site
    return [
        (x + 1, y, z), (x - 1, y, z),
        (x, y + 1, z), (x, y - 1, z),
        (x, y, z + 1), (x, y, z - 1),
    ]


def box_sites(shape):
    a, b, c = shape
    return sorted(
        (x, y, z) for x in range(a) for y in range(b) for z in range(c)
    )


# --------------------------------------------------------------------------- #
# independent census of the private-neighbour lemma                            #
# --------------------------------------------------------------------------- #
def census(shape):
    sites = box_sites(shape)
    n = len(sites)
    nonempty = 0
    maximizer_checks = 0
    maximizer_failures = 0
    nonmax_fail = 0
    witnesses = []
    for mask in range(1, 1 << n):
        subset = {sites[i] for i in range(n) if (mask >> i) & 1}
        nonempty += 1
        first_max = max(s[0] for s in subset)
        for u in subset:
            east = (u[0] + 1, u[1], u[2])
            hit = subset.intersection(neighbors(east))
            ok = east not in subset and hit == {u}
            if u[0] == first_max:
                maximizer_checks += 1
                if not ok:
                    maximizer_failures += 1
            elif not ok and len(witnesses) < 3:
                witnesses.append(
                    {
                        "subset": sorted(subset),
                        "u": list(u),
                        "u_plus_e1": list(east),
                        "support_neighbors_of_u_plus_e1": sorted(hit),
                    }
                )
        for u in subset:
            if u[0] < first_max:
                east = (u[0] + 1, u[1], u[2])
                if east in subset or subset.intersection(neighbors(east)) != {u}:
                    nonmax_fail += 1
                    break
    return {
        "shape": list(shape),
        "n_sites": n,
        "nonempty_subsets": nonempty,
        "maximizer_checks": maximizer_checks,
        "maximizer_failures": maximizer_failures,
        "nonmaximizer_failure_subsets": nonmax_fail,
        "example_witnesses": witnesses,
    }


def verify_census(payload):
    stored = payload["data"]["exhaustive_private_neighbor"]
    assert stored["claim_tag"] == "[COMPUTATION]"
    totals = [0, 0, 0]
    for box in stored["boxes"]:
        mine = census(tuple(box["shape"]))
        for key in (
            "n_sites",
            "nonempty_subsets",
            "maximizer_checks",
            "maximizer_failures",
            "nonmaximizer_failure_subsets",
        ):
            assert mine[key] == box[key], (box["shape"], key, mine[key], box[key])
        assert len(box["example_witnesses"]) == len(mine["example_witnesses"])
        for witness, stored_witness in zip(mine["example_witnesses"], box["example_witnesses"]):
            assert witness["u"] == stored_witness["u"]
            assert witness["u_plus_e1"] == stored_witness["u_plus_e1"]
            subset = {tuple(s) for s in stored_witness["subset"]}
            east = tuple(stored_witness["u_plus_e1"])
            u = tuple(stored_witness["u"])
            assert u in subset
            assert east in subset or subset.intersection(neighbors(east)) != {u}
        totals[0] += box["nonempty_subsets"]
        totals[1] += box["maximizer_checks"]
        totals[2] += box["maximizer_failures"]
    assert stored["totals"] == {
        "subsets": totals[0],
        "maximizer_checks": totals[1],
        "failures": totals[2],
    }
    assert totals[2] == 0
    assert totals[0] > 200_000, "exhaustive budget must include the 18-site box"


# --------------------------------------------------------------------------- #
# ordered-Pauli algebra via operator products (independent of the producer)    #
# --------------------------------------------------------------------------- #
def op_mul(left, right):
    """(X^a Z^b)(X^c Z^d) = (-1)^{|b & c|} X^{a^c} Z^{b^d}, exact."""
    out = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            sign = -1 if bin(b & c).count("1") % 2 else 1
            key = (a ^ c, b ^ d)
            out[key] = out.get(key, 0) + sign * x * y
    return {k: v for k, v in out.items() if v}


def op_sub(left, right):
    out = dict(left)
    for k, v in right.items():
        w = out.get(k, 0) - v
        if w:
            out[k] = w
        else:
            out.pop(k, None)
    return out


def op_scale(op, scale):
    return {k: v * scale for k, v in op.items() if v * scale}


class Model:
    """Finite shadow of the infinite-volume derivations on Alg(S)."""

    def __init__(self, sites, internal_only=False):
        self.sites = sorted(sites)
        siteset = set(self.sites)
        ext = list(dict.fromkeys(self.sites + [q for s in self.sites for q in neighbors(s)]))
        self.ext = ext
        self.idx = {p: i for i, p in enumerate(ext)}
        self.nt = len(ext)
        self.columns = []
        for column in range(1 << (2 * len(self.sites))):
            digits = column
            a = b = 0
            for p in self.sites:
                label = digits & 3
                digits >>= 2
                if label & 1:
                    a |= 1 << self.idx[p]
                if label & 2:
                    b |= 1 << self.idx[p]
            self.columns.append((a, b))
        bonds = []
        seen = set()
        for s in self.sites:
            for q in neighbors(s):
                key = tuple(sorted((s, q)))
                if key not in seen:
                    seen.add(key)
                    bonds.append((key[0], key[1], key[0] in siteset and key[1] in siteset))
        self.bonds = [bd for bd in bonds if (bd[2] or not internal_only)] if internal_only else bonds
        self.all_bonds = bonds

    # -- generators as single-term ordered-Pauli operators ------------------- #
    def x_op(self, site):
        return {(1 << self.idx[site], 0): 1}

    def zz_op(self, u, v):
        return {(0, (1 << self.idx[u]) | (1 << self.idx[v])): 1}

    def z_op(self, site):
        return {(0, 1 << self.idx[site]): 1}

    # -- commutators (1/2)[h, P] built from products ------------------------ #
    def half_commutator(self, pauli_term, generator):
        left = op_mul(generator, {pauli_term: 1})
        right = op_mul({pauli_term: 1}, generator)
        half = op_scale(op_sub(left, right), Fraction(1, 2))
        return {k: v for k, v in half.items() if v}

    def image_delta(self, column, a_coef, b_coef, bonds=None):
        """a*(1/2)[sum X, P] + b*(1/2)[sum ZZ, P] with P = columns[column]."""
        term = self.columns[column]
        out = {}
        for site in self.sites:
            piece = self.half_commutator(term, self.x_op(site))
            for k, v in piece.items():
                out[k] = out.get(k, Fraction(0)) + a_coef * v
        for u, v, _internal in (self.bonds if bonds is None else bonds):
            piece = self.half_commutator(term, self.zz_op(u, v))
            for k, w in piece.items():
                out[k] = out.get(k, Fraction(0)) + b_coef * w
        return {k: v for k, v in out.items() if v}

    def images(self, a_coef, b_coef, bonds=None):
        return [self.image_delta(c, a_coef, b_coef, bonds) for c in range(len(self.columns))]

    def tagged_images(self):
        """Images of (ad_A ; ad_B) as disjointly tagged vectors (a=b=1)."""
        out = []
        for c in range(len(self.columns)):
            vec = {}
            for k, v in self.image_delta(c, 1, 0).items():
                vec[(0,) + k] = v
            for k, v in self.image_delta(c, 0, 1).items():
                vec[(1,) + k] = v
            out.append(vec)
        return out


# --------------------------------------------------------------------------- #
# independent ranks                                                            #
# --------------------------------------------------------------------------- #
def dense_rank_q(images):
    keys = sorted({k for img in images for k in img})
    pos = {k: i for i, k in enumerate(keys)}
    rows = [[Fraction(0)] * len(keys) for _ in images]
    for r, img in enumerate(images):
        for k, v in img.items():
            rows[r][pos[k]] = Fraction(v)
    rank = 0
    used = set()
    for col in range(len(keys)):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                f = rows[r][col]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def reverse_rank_mod(images, prime):
    """Sparse modular elimination with reversed row order and greatest pivots."""
    pivots = {}
    rank = 0
    for source in reversed(images):
        vec = {}
        for k, v in source.items():
            r = (v.numerator * pow(v.denominator, prime - 2, prime)) % prime if isinstance(v, Fraction) else v % prime
            if r:
                vec[k] = r
        while vec:
            lead = max(vec)
            old = pivots.get(lead)
            if old is None:
                inv = pow(vec[lead], prime - 2, prime)
                pivots[lead] = {k: v * inv % prime for k, v in vec.items()}
                rank += 1
                break
            f = vec[lead]
            for k, v in old.items():
                r = (vec.get(k, 0) - f * v) % prime
                if r:
                    vec[k] = r
                else:
                    vec.pop(k, None)
    return rank


# --------------------------------------------------------------------------- #
# slot projections                                                             #
# --------------------------------------------------------------------------- #
def slot_key(model, key, site, alpha, beta):
    """Return stripped key if the (X^alpha Z^beta)_site slot matches, else None."""
    i = model.idx[site]
    a_bit = (key[0] >> i) & 1
    b_bit = (key[1] >> i) & 1
    if (a_bit, b_bit) != (alpha, beta):
        return None
    if beta:
        key = (key[0], key[1] ^ (1 << i))
    if alpha:
        key = (key[0] ^ (1 << i), key[1])
    return key


def slot_part(model, image, site, alpha, beta):
    out = {}
    for k, v in image.items():
        stripped = slot_key(model, k, site, alpha, beta)
        if stripped is not None:
            out[stripped] = out.get(stripped, Fraction(0)) + v
    return {k: v for k, v in out.items() if v}


# --------------------------------------------------------------------------- #
# verification of the stored representative supports                           #
# --------------------------------------------------------------------------- #
def canonical_peel(sites):
    return max(sorted(sites))


def verify_support(record):
    sites = [tuple(s) for s in record["sites"]]
    model = Model(sites)
    n_columns = 1 << (2 * len(sites))
    assert record["columns"] == n_columns
    assert record["bonds_meeting_support"] == len(model.all_bonds)
    assert record["internal_bonds"] == sum(1 for bd in model.all_bonds if bd[2])
    assert record["boundary_bonds"] == len(model.all_bonds) - record["internal_bonds"]

    v = canonical_peel(sites)
    w = (v[0] + 1, v[1], v[2])
    assert tuple(record["peel_v"]) == v and tuple(record["private_w"]) == w
    assert w not in set(sites)
    sset = set(sites)
    assert sset.intersection(neighbors(w)) == {v}, "private-neighbour property"

    # kernel dimensions: independent modular rank for every stored combination
    for combo in record["combinations"]:
        a_c, b_c = Fraction(combo["a"]), Fraction(combo["b"])
        images = model.images(a_c, b_c)
        nullity = n_columns - reverse_rank_mod(images, PRIME)
        assert nullity == combo["nullity_Q"] == 1, (record["name"], combo)
    tagged = model.tagged_images()
    assert n_columns - reverse_rank_mod(tagged, PRIME) == record["simultaneous"]["nullity_Q"] == 1

    # dense exact-Q rank for small supports
    if len(sites) <= 4:
        images = model.images(Fraction(1), Fraction(1))
        rank = dense_rank_q(images)
        combo11 = next(c for c in record["combinations"] if c["a"] == 1 and c["b"] == 1)
        assert rank == combo11["rank_Q"] == n_columns - 1, (record["name"], rank)

    # identity column spans the kernel: delta(I) = 0 and rank = n-1
    assert not model.image_delta(0, Fraction(1), Fraction(1))

    # peeling identity I1: pi^{(w)}_{01} of delta_g = b * (1/2)[Z_v, .]
    for combo in record["combinations"]:
        a_c, b_c = Fraction(combo["a"]), Fraction(combo["b"])
        for c in range(n_columns):
            image = model.image_delta(c, a_c, b_c)
            lhs = slot_part(model, image, w, 0, 1)
            rhs = op_scale(model.half_commutator(model.columns[c], model.z_op(v)), b_c)
            assert lhs == rhs, ("I1", record["name"], combo, c)
            # no X_w or Y_w components ever appear
            assert not slot_part(model, image, w, 1, 0)
            assert not slot_part(model, image, w, 1, 1)

    # peeling identity I2b: after [Q,Z_v]=0, only a_v=0 columns remain;
    # pi^{(v)}_{11} of (1/2)[A,Q] is exactly the stripped Z_v coefficient.
    a_c = Fraction(1)
    vi = model.idx[v]
    for c in range(n_columns):
        a_mask, b_mask = model.columns[c]
        image = model.image_delta(c, a_c, Fraction(0))
        lhs = slot_part(model, image, v, 1, 1)
        if (a_mask >> vi) & 1:
            continue
        if (b_mask >> vi) & 1:
            assert lhs == {(a_mask, b_mask ^ (1 << vi)): Fraction(1)}, ("I2b", record["name"], c)
            # No local-basis cancellation is possible on Z_v-slotted columns.
            full = slot_part(model, model.image_delta(c, Fraction(1), Fraction(1)), v, 1, 1)
            assert full, ("cancellation-sensitivity", record["name"], c)
        else:
            assert not lhs, ("I2b-empty", record["name"], c)

    # ZZ contributions to the (1,1)@v slot exist only from a_v = 1 columns
    for c in range(n_columns):
        a_mask, _ = model.columns[c]
        zz = slot_part(model, model.image_delta(c, Fraction(0), Fraction(1)), v, 1, 1)
        if (a_mask >> vi) & 1:
            continue
        assert not zz, ("psi-support", record["name"], c)

    # Every w-bond whose *other* endpoint is outside S has zero commutator:
    # Q has zero X-mask at w and at that endpoint.  The latter may lie outside
    # the finite bookkeeping legend, so use the ordered-Pauli criterion rather
    # than materializing an unnecessary tensor factor.
    for other in neighbors(w):
        if other in sset:
            continue
        other_index = model.idx.get(other)
        for a_mask, _ in model.columns:
            a_w = bool(a_mask & (1 << model.idx[w]))
            a_other = bool(a_mask & (1 << other_index)) if other_index is not None else False
            assert not (a_w ^ a_other)

    # sabotage: dropping boundary bonds breaks identity I1
    internal_bonds = [bd for bd in model.all_bonds if bd[2]]
    broken = False
    for c in range(n_columns):
        lhs = slot_part(model, model.image_delta(c, Fraction(1), Fraction(1), bonds=internal_bonds), w, 0, 1)
        rhs = model.half_commutator(model.columns[c], model.z_op(v))
        if lhs != rhs:
            broken = True
            break
    assert broken, ("sabotage-boundary", record["name"])

    # internal-only contrast record
    if "internal_only" in record:
        inside = Model(sites, internal_only=True)
        tagged_int = inside.tagged_images()
        stored = record["internal_only"]
        if "simultaneous_nullity_Q" in stored:
            assert n_columns - reverse_rank_mod(tagged_int, PRIME) == stored["simultaneous_nullity_Q"]


def verify_local_maps(payload):
    maps = payload["data"]["local_coefficient_maps"]
    assert maps["claim_tag"] == "[COMPUTATION]"
    for name in ("single_site", "bond"):
        record = maps[name]
        sites = [tuple(s) for s in record["sites"]]
        model = Model(sites)
        v = canonical_peel(sites)
        w = (v[0] + 1, v[1], v[2])
        assert tuple(record["peel_v"]) == v and tuple(record["private_w"]) == w
        legend = record["sites_legend"]
        assert [tuple(s) for s in legend] == model.ext
        n_columns = 1 << (2 * len(sites))
        assert record["domain_columns"] == n_columns

        def entries_to_dict(entries):
            out = {}
            for target_a, target_b, source, coeff in entries:
                out[((target_a, target_b), source)] = Fraction(coeff)
            return out

        for slot_name, site, alpha, beta, part in (
            ("wslot_pencil", w, 0, 1, "b"),
            ("wslot_pencil", w, 0, 1, "a"),
            ("vslot_pencil", v, 1, 1, "a"),
            ("vslot_pencil", v, 1, 1, "b"),
        ):
            stored = entries_to_dict(record[slot_name][f"{part}_entries"])
            # recompute: the 'a'/'b' component of the slot map as a map of the pencil
            if part == "a":
                images = model.images(Fraction(1), Fraction(0))
            else:
                images = model.images(Fraction(0), Fraction(1))
            recomputed = {}
            for c in range(n_columns):
                proj = slot_part(model, images[c], site, alpha, beta)
                for k, val in proj.items():
                    recomputed[((k[0], k[1]), c)] = val
            assert stored == recomputed, (name, slot_name, part)

        facts = record["structural_facts"]
        assert facts["wslot_a_part_zero"]
        assert facts["vslot_b_part_zero_on_a_v_zero_columns"]
        assert facts["vslot_a_part_matches_z_slot_on_a_v_zero_columns"]
        # sensitivity: the b-part of the v-slot map is NOT identically zero
        b_entries = entries_to_dict(record["vslot_pencil"]["b_entries"])
        assert b_entries, (name, "vslot b-part must be nonzero in general")


def verify_meta(payload):
    prov = payload["provenance"]
    assert prov["script"] == "experiments/e115_local_commutant.py"
    assert all(check["passed"] for check in payload["checks"]), [
        c["name"] for c in payload["checks"] if not c["passed"]
    ]
    theorem = payload["data"]["theorem"]
    assert theorem["claim_tag"] == "[THEOREM]"
    for phrase in ("finite-support", "a", "scalar"):
        assert phrase in theorem["statement"] or phrase == "a"
    assert "corollary" in theorem and theorem["corollary"]["claim_tag"] == "[THEOREM]"
    scope = payload["data"]["scope_limits"]
    excluded = " ".join(scope["excluded"]).lower()
    for phrase in ("translation", "quasilocal", "torus", "a = 0"):
        assert phrase in excluded, phrase
    # The recorded private-neighbour lemma must match the geometry: w+e_1 is a
    # neighbour of w with first coordinate v_1+2, so a claim of equality would be
    # false.  Verify the true "at least" form both textually and geometrically.
    lemma = payload["data"]["conventions"]["private_neighbor_lemma"]
    assert "at least v_1+1" in lemma, lemma
    v = (0, 0, 0)
    w = (v[0] + 1, v[1], v[2])
    neighbours = [
        tuple(w[i] + s * (1 if i == axis else 0) for i in range(3))
        for axis in range(3)
        for s in (1, -1)
    ]
    others = [p for p in neighbours if p != v]
    assert len(others) == 5, others
    assert all(p[0] >= v[0] + 1 for p in others), others
    assert any(p[0] == v[0] + 2 for p in others), others


def verify_e52_cross(payload):
    stored = json.loads(E52.read_text(encoding="utf-8"))
    dims = stored["data"]["finite_exact_dimensions"]
    summary = payload["data"]["internal_contrast_summary"]
    rows = {row["name"]: row for row in summary["rows"]}
    assert rows["plaquette"]["internal_simultaneous_nullity_Q"] == dims["2x2"] == 27
    assert rows["slab_2x3"]["internal_simultaneous_nullity_Q"] == dims["2x3"] == 12


def main() -> None:
    signal.alarm(TIMEOUT_SECONDS)
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    verify_meta(payload)
    verify_census(payload)
    for record in payload["data"]["representative_supports"]:
        verify_support(record)
    verify_local_maps(payload)
    verify_e52_cross(payload)
    print("PASS")


if __name__ == "__main__":
    main()
