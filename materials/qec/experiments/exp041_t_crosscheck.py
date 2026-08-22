"""EXP-041: independent, SAT-free cross-check of EXP-039's T(P) and d_Z(P).

Theorem H's strength depends on T(P) being exact.  EXP-039 obtained T by
iterated SAT queries (``minimum_weight_module`` in
``src/qec_research/codes/pbb_nogo.py``).  This experiment re-derives T with a
completely independent algorithm that never touches a SAT solver:

1. **Exhaustive meet-in-the-middle enumeration.**  A Z-sector vector is
   ``z = (z1 | z2)`` with each block ``L = ell*m`` bits; the pure-Z centralizer
   condition is ``A z1 = B z2`` over GF(2).  For every weight split
   ``(w1, w2)`` with ``w1 + w2 <= d`` we enumerate all ``C(L, w1)`` masks,
   compute each syndrome ``A z1`` by XORing prepacked column bitmasks, and hash
   it against the ``B z2`` syndromes of all ``C(L, w2)`` masks.  Equality of
   uint64 syndrome bitmasks is exact, so the matched set is *exactly*
   ``{ z in ker[A B] : wt(z) <= d }`` -- no approximation, no solver.

2. **Nontriviality.**  Each kept vector is tested against ``S_Z =
   rowspace(H_Z)`` by echelon (rank) reduction: ``z`` survives iff
   ``rank(S_Z + z) > rank(S_Z)``.  The smallest weight with a survivor is the
   brute-force ``d_Z_bf``; it must equal the certificate's ``d_z_parent``, and
   every certified SAT witness must appear in the survivor set at exactly that
   weight.

3. **T from translation orbits.**  For each minimum-weight survivor, all
   ``ell*m`` cyclic translates ``x^a y^b`` (own implementation) are inserted
   into an echelon basis that starts as ``S_Z``; the number of inserted rows is
   ``T_bf = rank(M + S_Z) - rank(S_Z)``.  Since ``M subset ker[A B]`` and
   ``dim ker[A B] - rank S_Z = k_parent``, accumulation stops early (soundly)
   once ``T_bf = k_parent``.

Scope rule (per assignment): enumeration cost ``sum_{w<=d_Z} C(L, w)``;
parents whose cost exceeds 2e7 are recorded as skipped (the Gross code's
parent among them).  Self-tests before the main loop: own orbit function vs
``pbb_survival.translation_orbit`` (bit-exact, random vectors); own echelon
rank vs ``rank_np`` (random matrices); the MITM enumerator vs full
``nullspace_np`` span enumeration on every parent with ker dimension <= 20
(all three n=36 parents: bit-exact set equality).

No SAT solver is used anywhere; ``pysat``/``sat_decide`` are never imported.
Pure GF(2) linear algebra and integer bitmask arithmetic only.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCHEMA = "exp041-t-crosscheck-v1"
COST_BUDGET = 2e7
CERT_DIR = ROOT / "results" / "partial_runs" / "exp039"
OUT = ROOT / "results" / "processed" / "exp041_t_crosscheck.json"

# Catalogue/parent-matrix construction must be identical to the prime path, so
# the loader is shared (read-only).  Everything *computed* from the matrices
# below is independent: no sat_decide, no pbb_nogo, no minimum_weight_module.
_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E27
_SPEC.loader.exec_module(E27)

from qec_research.codes.pbb_survival import (  # noqa: E402  # orbit VALIDATION only
    translation_orbit as reference_translation_orbit,
)
from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ---------------------------------------------------------------------------
# Own GF(2) machinery (independent of the SAT path)
# ---------------------------------------------------------------------------
def rows_to_ints(M: np.ndarray) -> list[int]:
    """Row -> Python int bitmask, bit j == column j."""
    return [
        int.from_bytes(np.packbits(row, bitorder="little").tobytes(), "little")
        for row in M
    ]


def int_to_row(z: int, ncols: int) -> np.ndarray:
    nbytes = (ncols + 7) // 8
    return np.unpackbits(
        np.frombuffer(z.to_bytes(nbytes, "little"), dtype=np.uint8),
        bitorder="little",
    )[:ncols]


def column_masks(M: np.ndarray) -> np.ndarray:
    """uint64 array; entry j is the row-bitmask of M[:, j] (bit r == row r)."""
    nrows, ncols = M.shape
    if nrows > 64 or ncols > 64:
        raise ValueError("column_masks assumes dimension <= 64 here")
    out = np.zeros(ncols, dtype=np.uint64)
    for j in range(ncols):
        bits = np.flatnonzero(M[:, j])
        v = np.uint64(0)
        for r in bits:
            v |= np.uint64(1) << np.uint64(r)
        out[j] = v
    return out


class Echelon:
    """Row-echelon GF(2) basis over Python ints; pivot = highest set bit.

    Invariant: distinct pivots.  ``reduce`` XORs away every bit that is an
    existing pivot (descending), so ``reduce(t) == 0`` iff ``t`` lies in the
    span.  ``insert`` returns True exactly when the span grew.
    """

    def __init__(self, rows: list[int] = ()) -> None:
        self.pivots: dict[int, int] = {}
        for row in rows:
            self.insert(row)

    def copy(self) -> "Echelon":
        e = Echelon()
        e.pivots = dict(self.pivots)
        return e

    def reduce(self, t: int) -> int:
        while t:
            p = t.bit_length() - 1
            b = self.pivots.get(p)
            if b is None:
                break
            t ^= b
        return t

    def insert(self, t: int) -> bool:
        r = self.reduce(t)
        if not r:
            return False
        self.pivots[r.bit_length() - 1] = r
        return True

    def contains(self, t: int) -> bool:
        return self.reduce(t) == 0

    @property
    def rank(self) -> int:
        return len(self.pivots)


def orbit_masks(z: int, ell: int, m: int) -> list[int]:
    """Own implementation of all ``ell*m`` translates of a 2-block vector.

    Bit layout: block 0 occupies bits [0, L), block 1 bits [L, 2L),
    position p = i*m + j inside a block.  The (a, b) translate moves
    (i, j) -> (i + a mod ell, j + b mod m), matching np.roll on each axis.
    """
    L = ell * m
    blocks = np.empty((2, L), dtype=np.uint8)
    half = int_to_row(z & ((1 << L) - 1), L)
    blocks[0] = half
    blocks[1] = int_to_row(z >> L, L)
    arr = blocks.reshape(2, ell, m)
    out: list[int] = []
    for a in range(ell):
        rolled_a = np.roll(arr, a, axis=1)
        for b in range(m):
            flat = np.roll(rolled_a, b, axis=2).reshape(2 * L)
            out.append(
                int.from_bytes(
                    np.packbits(flat, bitorder="little").tobytes(), "little"
                )
            )
    return out


# ---------------------------------------------------------------------------
# Meet-in-the-middle exact enumerator
# ---------------------------------------------------------------------------
def side_enumeration(
    L: int, w: int, cols: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """All weight-w masks of L bits and their syndrome bitmasks.

    ``cols[j]`` is the syndrome bitmask contributed by setting position j;
    the syndrome of a mask is the XOR of its positions' entries (exact GF(2)).
    """
    total = _comb(L, w)
    masks = np.empty(total, dtype=np.uint64)
    synd = np.empty(total, dtype=np.uint64)
    if w == 0:
        masks[0] = 0
        synd[0] = 0
        return masks, synd
    pos = np.fromiter(
        itertools.chain.from_iterable(itertools.combinations(range(L), w)),
        dtype=np.int64,
        count=total * w,
    ).reshape(total, w)
    CHUNK = 262144
    for lo in range(0, total, CHUNK):
        hi = min(total, lo + CHUNK)
        p = pos[lo:hi].astype(np.uint64)
        masks[lo:hi] = np.bitwise_or.reduce(  # positions are distinct => OR==SUM
            np.uint64(1) << p, axis=1
        )
        synd[lo:hi] = np.bitwise_xor.reduce(cols[p], axis=1)
    return masks, synd


def _comb(n: int, k: int) -> int:
    from math import comb

    return comb(n, k)


def enumeration_cost(L: int, d: int) -> int:
    return sum(_comb(L, w) for w in range(d + 1))


class SideCache:
    """Cached (masks, syndromes, sorted unique syndromes) per block side/weight."""

    def __init__(self, L: int, d: int, cols_a: np.ndarray, cols_b: np.ndarray):
        self.L = L
        self.data: dict[tuple[str, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for w in range(d + 1):
            for side, cols in (("A", cols_a), ("B", cols_b)):
                masks, synd = side_enumeration(L, w, cols)
                self.data[(side, w)] = (masks, synd, np.unique(synd))

    def split_matches(self, w1: int, w2: int) -> list[int]:
        """All z=(z1|z2) with wt(z1)=w1, wt(z2)=w2 and A z1 = B z2."""
        masks_a, synd_a, uniq_a = self.data[("A", w1)]
        masks_b, synd_b, uniq_b = self.data[("B", w2)]
        keys = np.intersect1d(uniq_a, uniq_b)
        if keys.size == 0:
            return []
        left_idx = np.nonzero(np.isin(synd_a, keys))[0]
        right_idx = np.nonzero(np.isin(synd_b, keys))[0]
        right_map: dict[int, list[int]] = {}
        for j in right_idx:
            right_map.setdefault(int(synd_b[j]), []).append(int(masks_b[j]))
        out: list[int] = []
        L = self.L
        for i in left_idx:
            lm = int(masks_a[i])
            for rm in right_map[int(synd_a[i])]:
                out.append(lm | (rm << L))
        return out

    def all_kernel_vectors(self, d: int) -> dict[int, list[int]]:
        """weight -> every z in ker[A B] of that weight, for weight <= d."""
        by_weight: dict[int, list[int]] = {w: [] for w in range(d + 1)}
        for w1 in range(d + 1):
            for w2 in range(d + 1 - w1):
                matches = self.split_matches(w1, w2)
                if matches:
                    by_weight[w1 + w2].extend(matches)
        return by_weight


# ---------------------------------------------------------------------------
# Brute-force T for one parent
# ---------------------------------------------------------------------------
def brute_force_parent(
    HX: np.ndarray,
    HZ: np.ndarray,
    ell: int,
    m: int,
    d_cert: int,
    k_parent: int,
    cert_witnesses: list[list[int]],
) -> tuple[dict[str, Any], dict[int, list[int]]]:
    L = ell * m
    n = 2 * L
    t0 = time.perf_counter()

    A = np.ascontiguousarray(HX[:, :L], dtype=np.uint8) & 1
    B = np.ascontiguousarray(HX[:, L:], dtype=np.uint8) & 1
    cols_a = column_masks(A)
    cols_b = column_masks(B)

    rank_hx = int(rank_np(np.ascontiguousarray(HX)))
    rank_hz = int(rank_np(np.ascontiguousarray(HZ)))
    k_from_matrices = n - rank_hx - rank_hz

    cache = SideCache(L, d_cert, cols_a, cols_b)
    by_weight = cache.all_kernel_vectors(d_cert)

    total_matches = sum(len(v) for v in by_weight.values())
    distinct = set()
    for v in by_weight.values():
        distinct.update(v)
    matches_are_distinct = len(distinct) == total_matches

    # S_Z membership is decided by echelon reduction == rank test.
    sz_basis = Echelon(rows_to_ints(np.ascontiguousarray(HZ)))
    if sz_basis.rank != rank_hz:
        raise AssertionError("echelon rank of H_Z disagrees with rank_np")

    # Nontrivial vectors per weight; d_Z_bf = first weight with a survivor.
    nontriv: dict[int, list[int]] = {}
    for w in range(d_cert + 1):
        keep = [z for z in by_weight[w] if not sz_basis.contains(z)]
        if keep:
            nontriv[w] = sorted(keep)
    d_z_bf = min(nontriv) if nontriv else None

    # Certificate witnesses must lie in the survivor set at exactly d_cert.
    witness_masks = {
        int.from_bytes(
            np.packbits(np.asarray(wv, dtype=np.uint8), bitorder="little").tobytes(),
            "little",
        )
        for wv in cert_witnesses
    }
    witness_weights_ok = all(
        bin(w).count("1") == d_cert for w in witness_masks
    )
    at_d = set(nontriv.get(d_cert, []))
    witnesses_contained = witness_masks <= at_d

    # T_bf: span of orbits of ALL minimum-weight logicals, modulo S_Z.
    t_bf = 0
    orbit_rows_seen = 0
    translates_in_kernel_checked = True
    if d_z_bf is not None:
        logicals = nontriv[d_z_bf]
        basis_t = sz_basis.copy()
        ker_check_budget = 25  # direct per-translate ker checks on a prefix
        for li, z in enumerate(logicals):
            orb = orbit_masks(z, ell, m)
            orbit_rows_seen += len(orb)
            if li < ker_check_budget:
                for t in orb:
                    t1, t2 = t & ((1 << L) - 1), t >> L
                    s = 0
                    x = t1
                    while x:
                        p = (x & -x).bit_length() - 1
                        s ^= int(cols_a[p])
                        x &= x - 1
                    y = t2
                    while y:
                        p = (y & -y).bit_length() - 1
                        s ^= int(cols_b[p])
                        y &= y - 1
                    if s != 0:
                        translates_in_kernel_checked = False
            for t in orb:
                if basis_t.insert(t):
                    t_bf += 1
            if t_bf >= k_parent:
                break
        if t_bf > k_parent:  # impossible if everything above is correct
            raise AssertionError("T exceeded k_parent: quotient bound violated")

    return {
        "d_z_bf": d_z_bf,
        "t_bf": t_bf,
        "k_from_matrices": int(k_from_matrices),
        "rank_hx": rank_hx,
        "rank_hz": rank_hz,
        "num_enumerated_kernel_vectors": total_matches,
        "matches_are_distinct": matches_are_distinct,
        "num_nontrivial_by_weight": {str(w): len(v) for w, v in sorted(nontriv.items())},
        "num_min_weight_logicals": len(nontriv[d_z_bf]) if d_z_bf is not None else 0,
        "num_cert_witnesses": len(witness_masks),
        "witness_weights_ok": witness_weights_ok,
        "witnesses_contained_in_bf_set": witnesses_contained,
        "orbit_rows_seen": orbit_rows_seen,
        "translates_in_kernel_checked": translates_in_kernel_checked,
        "wall_time_s": time.perf_counter() - t0,
    }, by_weight


# ---------------------------------------------------------------------------
# Self-tests (independence protocol)
# ---------------------------------------------------------------------------
def self_tests(shapes: set[tuple[int, int]]) -> dict[str, Any]:
    rng = np.random.default_rng(20260817)
    out: dict[str, Any] = {}

    # 1. Own orbit implementation, bit-exact against the shared reference.
    orbit_ok = True
    for ell, m in sorted(shapes):
        L = ell * m
        for _ in range(20):
            v = (rng.random(2 * L) < 0.2).astype(np.uint8)
            z = int.from_bytes(
                np.packbits(v, bitorder="little").tobytes(), "little"
            )
            ref = reference_translation_orbit(v, ell, m, blocks=2)
            ref_masks = sorted(rows_to_ints(ref))
            mine = sorted(orbit_masks(z, ell, m))
            if ref_masks != mine:
                orbit_ok = False
    if not orbit_ok:
        raise AssertionError("self-test: orbit_masks disagrees with reference")
    out["orbit_matches_reference"] = True

    # 2. Own echelon rank vs rank_np on random matrices.
    rng2 = np.random.default_rng(7)
    for rows, cols in ((30, 36), (40, 72), (20, 72), (48, 108)):
        for _ in range(25):
            M = (rng2.random((rows, cols)) < 0.35).astype(np.uint8)
            e = Echelon(rows_to_ints(M))
            if e.rank != int(rank_np(M)):
                raise AssertionError("self-test: echelon rank disagrees with rank_np")
    out["echelon_rank_matches_rank_np"] = True

    # 3. Column-XOR syndrome equals dense matmul.
    rng3 = np.random.default_rng(11)
    M = (rng3.random((36, 36)) < 0.3).astype(np.uint8)
    cols = column_masks(M)
    for _ in range(200):
        z = int(rng3.integers(0, 1 << 36))
        s = 0
        x = z
        while x:
            p = (x & -x).bit_length() - 1
            s ^= int(cols[p])
            x &= x - 1
        dense = (M.astype(np.int64) @ int_to_row(z, 36).astype(np.int64)) % 2
        want = int.from_bytes(np.packbits(dense, bitorder="little").tobytes(), "little")
        if s != want:
            raise AssertionError("self-test: column-XOR syndrome wrong")
    out["column_syndrome_matches_dense"] = True
    return out


def nullspace_cross_check(
    HX: np.ndarray, by_weight: dict[int, list[int]], d_cert: int
) -> dict[str, Any] | None:
    """Enumerate the FULL span of nullspace_np(HX) via Gray code and compare
    the weight<=d set bit-exactly with the MITM enumerator's output."""
    N = nullspace_np(np.ascontiguousarray(HX))
    kdim = int(N.shape[0])
    if kdim > 20:
        return None
    basis = rows_to_ints(N)
    found: set[int] = set()
    z = 0
    if z.bit_count() <= d_cert:
        found.add(z)
    for i in range(1, 1 << kdim):
        g = i ^ (i >> 1)
        prev = (i - 1) ^ ((i - 1) >> 1)
        z ^= basis[(g ^ prev).bit_length() - 1]
        if z.bit_count() <= d_cert:
            found.add(z)
    mine: set[int] = set()
    for v in by_weight.values():
        mine.update(v)
    return {
        "nullspace_dim": kdim,
        "span_enumerated": 1 << kdim,
        "set_equal": found == mine,
        "count_span_side": len(found),
        "count_mitm_side": len(mine),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def load_certificates_local() -> dict[str, dict[str, Any]]:
    """Same source as exp039_nogo_module.load_certificates(): the atomic
    per-parent JSON certificates.  Reimplemented to keep this file free of
    any import from the SAT path."""
    out: dict[str, dict[str, Any]] = {}
    for path in sorted(CERT_DIR.glob("parent_*.json")):
        rec = json.loads(path.read_text(encoding="utf-8"))
        rec["_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        out[rec["fingerprint"]] = rec
    return out


def catalogue_members() -> dict[str, dict[str, Any]]:
    """fingerprint -> {spec, catalogue labels, matrices} for every catalogue parent."""
    rows = E27.load_catalogue()
    out: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(rows):
        spec, HX, HZ = E27.parent_matrices(row)
        fp = E27.matrix_fingerprint(HX, HZ)
        entry = out.setdefault(
            fp,
            {
                "fingerprint": fp,
                "ell": int(spec.ell),
                "m": int(spec.m),
                "n": int(HX.shape[1]),
                "HX": HX,
                "HZ": HZ,
                "labels": [],
            },
        )
        entry["labels"].append(E27.catalogue_label(row, index))
    return out


def run() -> dict[str, Any]:
    started = time.perf_counter()
    certs = load_certificates_local()
    parents = catalogue_members()

    processed_records: list[dict[str, Any]] = []
    skipped_records: list[dict[str, Any]] = []
    uncertified: list[dict[str, Any]] = []
    shapes = {
        (c["ell"], c["m"])
        for c in certs.values()
        if enumeration_cost(c["ell"] * c["m"], c["d_z_parent"]) <= COST_BUDGET
    }
    checks = self_tests(shapes)

    for fp, cert in sorted(certs.items(), key=lambda kv: (kv[1]["ell"] * kv[1]["m"], kv[0])):
        ell, m, d_cert = int(cert["ell"]), int(cert["m"]), int(cert["d_z_parent"])
        L = ell * m
        cost = enumeration_cost(L, d_cert)
        parent = parents.get(fp)
        if parent is None:
            skipped_records.append(
                {"fingerprint": fp, "reason": "no catalogue parent matches fingerprint"}
            )
            continue
        if cost > COST_BUDGET:
            skipped_records.append(
                {
                    "fingerprint": fp,
                    "ell": ell,
                    "m": m,
                    "n": int(cert["n"]),
                    "d_z_cert": d_cert,
                    "enumeration_cost": cost,
                    "reason": f"enumeration cost {cost} exceeds budget {int(COST_BUDGET)}",
                }
            )
            continue

        print(
            f"[exp041] brute force ell={ell} m={m} n={2 * L} d={d_cert} "
            f"cost={cost} fp={fp[:12]} ...",
            flush=True,
        )
        bf, by_weight = brute_force_parent(
            parent["HX"],
            parent["HZ"],
            ell,
            m,
            d_cert,
            int(cert["k_parent"]),
            cert.get("witness_vectors") or [],
        )
        # Full-span cross-check of the MITM enumerator when affordable.
        ns_check = nullspace_cross_check(parent["HX"], by_weight, d_cert)
        if ns_check is not None and not ns_check["set_equal"]:
            bf["nullspace_cross_check"] = ns_check
            raise AssertionError(
                f"MITM enumerator disagrees with nullspace span on {fp[:16]}"
            )
        bf["nullspace_cross_check"] = ns_check

        d_z_bf = bf["d_z_bf"]
        t_bf = bf["t_bf"]
        d_z_match = d_z_bf == d_cert
        t_cert = int(cert["T"])
        t_match = t_bf == t_cert
        k_match = bf["k_from_matrices"] == int(cert["k_parent"])
        integrity = (
            bf["matches_are_distinct"]
            and bf["witness_weights_ok"]
            and bf["witnesses_contained_in_bf_set"]
            and bf["translates_in_kernel_checked"]
            and (ns_check is None or ns_check["set_equal"])
            and k_match
        )
        status = "OK" if (d_z_match and t_match and integrity) else "MISMATCH"
        processed_records.append(
            {
                "fingerprint": fp,
                "cert_sha256": cert["_sha256"],
                "labels": parent["labels"],
                "ell": ell,
                "m": m,
                "n": int(cert["n"]),
                "enumeration_cost": cost,
                "d_Z_bf": d_z_bf,
                "d_Z_cert": d_cert,
                "d_z_match": d_z_match,
                "T_bf": t_bf,
                "T_cert": t_cert,
                "T_match": t_match,
                "T_is_exact_cert": bool(cert["T_is_exact"]),
                "family_closed_cert": bool(cert["family_closed"]),
                "k_parent_cert": int(cert["k_parent"]),
                "k_parent_match": k_match,
                "brute_force": bf,
                "status": status,
            }
        )
        print(
            f"[exp041]   -> d_Z_bf={d_z_bf} (cert {d_cert}) "
            f"T_bf={t_bf} (cert {t_cert}) status={status} "
            f"({bf['wall_time_s']:.1f}s)",
            flush=True,
        )

    # Parents of the catalogue the SAT path never certified (informational).
    for fp, parent in sorted(parents.items()):
        if fp not in certs:
            uncertified.append(
                {
                    "fingerprint": fp,
                    "ell": parent["ell"],
                    "m": parent["m"],
                    "n": parent["n"],
                    "labels": parent["labels"],
                    "reason": "no EXP-039 certificate exists; out of scope for comparison",
                }
            )

    n_mismatch = sum(1 for r in processed_records if r["status"] != "OK")
    if n_mismatch:
        bad = [r["fingerprint"][:16] for r in processed_records if r["status"] != "OK"]
        verdict = "MISMATCH"
        verdict_detail = (
            f"{n_mismatch} of {len(processed_records)} processed parents disagree "
            f"with the EXP-039 certificates: {bad}"
        )
    elif skipped_records:
        verdict = "PARTIAL"
        verdict_detail = (
            f"{len(skipped_records)} rows skipped (enumeration cost > 2e7 budget), "
            f"all {len(processed_records)} processed parents match"
        )
    else:
        verdict = "MATCH_ALL"
        verdict_detail = f"all {len(processed_records)} certified parents match"

    return {
        "schema": SCHEMA,
        "utc": utc_now(),
        "budget": int(COST_BUDGET),
        "cert_dir": str(CERT_DIR),
        "num_certificates": len(certs),
        "num_catalogue_parents": len(parents),
        "self_tests": checks,
        "processed": processed_records,
        "skipped": skipped_records,
        "uncertified_parents": uncertified,
        "num_processed": len(processed_records),
        "num_skipped": len(skipped_records),
        "num_mismatch": n_mismatch,
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "sat_solver_used": False,
        "pysat_imported": "pysat" in sys.modules,
        "sat_decide_imported": "qec_research.distance.sat_decide" in sys.modules,
        "wall_time_s": time.perf_counter() - started,
    }


def main() -> int:
    payload = run()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(f".json.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(OUT)
    print(f"[exp041] wrote {OUT}")
    print(f"[exp041] verdict: {payload['verdict']} -- {payload['verdict_detail']}")
    return 0 if payload["verdict"] != "MISMATCH" else 1


if __name__ == "__main__":
    raise SystemExit(main())
