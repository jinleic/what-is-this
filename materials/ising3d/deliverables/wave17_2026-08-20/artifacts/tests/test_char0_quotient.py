"""Standalone verifier for experiments/e143_char0_quotient.py.

Recomputes, from the raw integers stored in results/algebra/char0_quotient.json and
independently implemented code (no import of the producer):
  * the Pauli-string basis closure (own BFS) and its size 1056,
  * the DIAGONAL Killing form entries kappa_vv (own popcount/dict implementation),
    the Killing rank / nullity / radical dimension, and the signature (544, 512),
  * perfection (derived support = whole basis),
  * the even-Z-weight parity fact behind the V+/V- split,
  * the invariant alternating forms B_+- : antisymmetry, exact invariance on all 13
    generator blocks (rebuilt here), nonvanishing integer determinant (own Bareiss),
  * the rank certificate: mod-p determinants of the stored 528x528 pivot minors
    (own mod-p LU) at both stored primes => rank_Q >= 528, plus the ambient bound
    dim sp(32) = 528 => rank_Q = 528 = saturation,
  * the quotient-dimension exclusion logic for sp(6)/gl(6)/sp(14),
  * sha256 checksums of the stored raw arrays.

Run:  nice -n 19 env PYTHONPATH=src .venv/bin/python tests/test_char0_quotient.py
Prints OK on success; raises AssertionError on any mismatch.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "results" / "algebra" / "char0_quotient.json"
P_PRIMES = [2147483647, 2147483629]


# ---------------------------------------------------------------- independent core

def sympl(v: int, w: int, n: int) -> int:
    mask = (1 << n) - 1
    return (bin((v >> n) & (w & mask)).count("1")
            + bin((w >> n) & (v & mask)).count("1")) & 1


def bfs_closure(gens, n):
    """Worklist closure under v,w -> v^w for anticommuting pairs (independent of producer)."""
    seen = dict.fromkeys(gens)
    todo = list(gens)
    while todo:
        v = todo.pop()
        for w in list(seen):
            if sympl(v, w, n):
                s = v ^ w
                if s not in seen:
                    seen[s] = None
                    todo.append(s)
    return sorted(seen)


def sign_xz(va, wa, n):
    """(-1)^{b_v . a_w}: bit j of the string id holds a_j, bit n+j holds b_j."""
    bv = va >> n
    return 1 if bin(bv & wa).count("1") % 2 == 0 else -1


def killing_diag_independent(S, n):
    """kappa_vv = sum_j c(v,j) c(v, j^v), c(v,j) = s1 - s2 in {0, +-2}; pure python."""
    kap = []
    for v in S:
        total = 0
        for j in S:
            c1 = sign_xz(v, j, n) - sign_xz(j, v, n)
            if c1 == 0:
                continue
            t = v ^ j
            c2 = sign_xz(v, t, n) - sign_xz(t, v, n)
            total += c1 * c2
        kap.append(total)
    return kap


def pauli_64(v, n):
    """Q_v = X^a Z^b as a 2^n x 2^n integer matrix (row-major numpy, own kron order)."""
    X = np.array([[0, 1], [1, 0]], dtype=np.int64)
    Z = np.array([[1, 0], [0, -1]], dtype=np.int64)
    J = np.array([[0, -1], [1, 0]], dtype=np.int64)
    I2 = np.array([[1, 0], [0, 1]], dtype=np.int64)
    letters = {0: I2, 1: X, 2: Z, 3: J}
    M = np.array([[1]], dtype=np.int64)
    for site in range(n):
        a = (v >> site) & 1
        b = (v >> (n + site)) & 1
        M = np.kron(M, letters[a | (b << 1)])
    return M


def parity_pair_basis(n):
    """Representatives (s, s^comp) of the complement involution on {0,1}^n."""
    comp = (1 << n) - 1
    reps = [s for s in range(1 << n) if s < (s ^ comp)]
    return reps


def block_of(M, n, which):
    """Parity block in the SAME convention as the artifact: with pair vectors
    u_i = e_{s_i} + e_{s_i^comp} (plus) or e_{s_i} - e_{s_i^comp} (minus),
    block[i,j] = (u_i^T M u_j) / 2."""
    reps = parity_pair_basis(n)
    comp = (1 << n) - 1
    k = len(reps)
    U = np.zeros((1 << n, k), dtype=np.int64)
    for j, s in enumerate(reps):
        U[s, j] = 1
        U[s ^ comp, j] = 1 if which == "plus" else -1
    return (U.T @ M @ U) // 2


def det_bareiss(M):
    """Exact integer determinant, fraction-free Bareiss (independent of sympy)."""
    A = [[int(x) for x in row] for row in M]
    nn = len(A)
    sign = 1
    prev = 1
    for k in range(nn - 1):
        if A[k][k] == 0:
            swap = next((r for r in range(k + 1, nn) if A[r][k] != 0), None)
            if swap is None:
                return 0
            A[k], A[swap] = A[swap], A[k]
            sign = -sign
        for i in range(k + 1, nn):
            for j in range(k + 1, nn):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[nn - 1][nn - 1]


def det_mod_p(M, p):
    """Determinant mod p by fraction-free-free LU with modular inverses (own code)."""
    A = (np.array(M, dtype=np.int64) % p).copy()
    nn = A.shape[0]
    det = 1
    for i in range(nn):
        nz = np.nonzero(A[i:, i] % p)[0]
        if len(nz) == 0:
            return 0
        r = i + int(nz[0])
        if r != i:
            A[[i, r]] = A[[r, i]]
            det = (-det) % p
        piv = int(A[i, i])
        det = (det * piv) % p
        inv = pow(piv, p - 2, p)
        fac = (A[i + 1:, i] * inv) % p
        A[i + 1:] = (A[i + 1:] - np.outer(fac, A[i])) % p
    return int(det % p)


def rank_mod_p_own(rows, p):
    A = (np.array(rows, dtype=np.int64) % p).copy()
    r = 0
    m, ncols = A.shape
    for c in range(ncols):
        nz = np.nonzero(A[r:, c] % p)[0]
        if len(nz) == 0:
            continue
        piv = r + int(nz[0])
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy()
        col[r] = 0
        mask = col % p != 0
        A[mask] = (A[mask] - np.outer(col[mask], A[r])) % p
        r += 1
        if r == m:
            break
    return r


# ---------------------------------------------------------------- main verification

def main():
    art = json.loads(ART.read_text())
    main_case = art["cases"]["main_grid_2x3"]
    n = main_case["n_sites"]
    assert n == 6
    S = main_case["string_basis"]
    d = main_case["dimension_Q"]

    # 1. closure re-derivation
    gens = main_case["generator_strings"]
    S2 = bfs_closure(gens, n)
    assert S2 == S, "string basis mismatch under independent BFS"
    assert len(S) == len(set(S)) == 1056 == d, "dim != 1056"

    # closure property re-check (pairwise, spot full sweep over generators x basis)
    Sset = set(S)
    for v in gens:
        for w in S:
            if sympl(v, w, n) and (v ^ w) not in Sset:
                raise AssertionError("closure violated")

    # 2. Killing form re-derivation: diagonal entries, rank, nullity, radical, signature
    kap = killing_diag_independent(S, n)
    assert kap == main_case["killing"]["diagonal"], "killing diagonal mismatch"
    nonzero = sum(1 for x in kap if x != 0)
    assert nonzero == d, "killing form degenerate: rank < dim"
    killing_rank, killing_nullity = nonzero, d - nonzero
    assert killing_nullity == 0 and main_case["radical_dimension_Q"] == 0
    sig = (sum(1 for x in kap if x > 0), sum(1 for x in kap if x < 0))
    assert sig == tuple(main_case["killing"]["signature"]) == (544, 512), sig
    # split sp_32 + sp_32 signature prediction: 2*((528+16)/2, (528-16)/2)
    assert sig == (2 * (528 + 16) // 2, 2 * (528 - 16) // 2)

    # 3. perfection re-derivation
    Tset = set()
    for v in S:
        for w in S:
            if sympl(v, w, n):
                Tset.add(v ^ w)
    assert Tset == Sset, "[g,g] support != g support (not perfect)"

    # 4. parity fact: even Z-weight for every basis string
    for v in S:
        z = bin(v >> n).count("1")
        assert z % 2 == 0, "odd Z-weight string in basis"

    # 5. invariant forms: rebuild generator blocks, exact invariance, Bareiss det
    for half_tag in ("plus", "minus"):
        f = main_case["invariant_forms"][half_tag]
        B = np.array(f["matrix"], dtype=np.int64)
        assert f["form_type"] == "alternating"
        assert np.array_equal(B, -B.T), "B not alternating"
        assert np.abs(B).max() <= 2, "B entries unexpectedly large"
        detB = det_bareiss(B.tolist())
        assert detB == f["det"] and detB != 0, f"B_{half_tag} determinant"
        which = "plus" if half_tag == "plus" else "minus"
        for gv in gens:
            X = block_of(pauli_64(gv, n), n, which)
            assert np.array_equal(X.T @ B + B @ X, np.zeros_like(B)), \
                f"invariance fails for generator {gv} on {half_tag}"

    # 6. rank certificates: rebuild the pivot-row blocks, verify minors mod both primes
    reps_needed = main_case["image_ranks"]
    for half_tag in ("plus", "minus"):
        cert = reps_needed[half_tag]
        which = "plus" if half_tag == "plus" else "minus"
        prow, pcol = cert["pivot_rows"], cert["pivot_cols"]
        assert len(prow) == len(pcol) == 528
        assert cert["rank_mod_p"][str(P_PRIMES[0])] == 528
        assert cert["rank_mod_p"][str(P_PRIMES[1])] == 528
        rows = []
        for ridx in prow:
            rows.append(block_of(pauli_64(S[ridx], n), n, which).reshape(-1))
        stack = np.array(rows, dtype=np.int64)
        minor = stack[:, pcol]
        for p in P_PRIMES:
            dm = det_mod_p(minor, p)
            assert dm != 0 and dm == cert["minor_det_mod_p"][str(p)], \
                f"pivot minor determinant mismatch {half_tag} p={p}"
        # independent full-stack rank at one prime
        r = rank_mod_p_own(stack, P_PRIMES[0])
        assert r == 528, f"full stack rank {r} != 528 on {half_tag}"

    # 7. ambient saturation arithmetic: dim sp(32) = m(2m+1) with m = 16
    m = 16
    assert m * (2 * m + 1) == 528
    amb = main_case["identification"]["plus"]["ambient_dim"]
    assert amb == 528 and main_case["identification"]["plus"]["saturates"] is True

    # 8. quotient-dimension exclusion logic
    quot = main_case["quotient_certificate"]
    assert quot["all_quotient_dimensions"] == [0, 528, 1056]
    assert quot["perfect"] is True
    for bad in (21, 36, 105):
        assert bad not in quot["all_quotient_dimensions"]

    # 9. checksums
    def sha(x):
        return hashlib.sha256(json.dumps(x).encode()).hexdigest()
    ck = art["checksums"]
    assert ck["main_string_basis_sha256"] == sha(S)
    assert ck["main_killing_diagonal_sha256"] == sha(main_case["killing"]["diagonal"])
    assert ck["main_B_plus_sha256"] == sha(main_case["invariant_forms"]["plus"]["matrix"])
    assert ck["main_B_minus_sha256"] == sha(main_case["invariant_forms"]["minus"]["matrix"])

    # 10. theorem field present and checks all pass
    assert art["theorem"]["claim_tag"].startswith("THEOREM")
    assert all(c["passed"] for c in art["checks"]), "stored checks not all passed"

    print("OK: char0_quotient artifact verified")
    print(f"  dim_Q g = {d}; Killing rank {killing_rank}, nullity {killing_nullity}, "
          f"signature {sig}; radical 0; perfect")
    print("  g = I_1 (+) I_2, dims 528/528, each ~= sp_32(Q) "
          "(invariant alternating forms dets "
          f"{det_bareiss(main_case['invariant_forms']['plus']['matrix'])}/"
          f"{det_bareiss(main_case['invariant_forms']['minus']['matrix'])}; "
          "528x528 pivot minors nonzero mod both primes)")
    print("  quotient dims {0,528,1056}: sp(6)/gl(6)/sp(14) excluded")


if __name__ == "__main__":
    sys.exit(main())
