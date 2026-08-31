#!/usr/bin/env python3
"""Gate C orientation sweep.

Swept set (per pre_statement.md of campaign 2026-08-30T113838Z...):
  - 5 decompositions (paper55, perminov58, sun56, mws59, stapleton60), all
    Z-ternary rank-23, each 729/729 Brent-anchored (verify_anchors.py).
  - Orientations: (X, Y, Z) over the 3x3 signed permutation matrices (48 each,
    110592 triples) x sigma-orbit {Id, sigma, sigma^2} (order 3) = 331776 per
    decomposition. (U,V)-swap orientations are NOT valid A*B decompositions
    (they compute B*A; gate C session 2 verified Brent-failure) and are excluded.
  - Objective per orientation: d(U'), d(V'), d(W'-factor) = counts of distinct
    non-input sign-classes of the 23 target vectors per side (rigorous LB on
    the per-side addition counts, per gate B's certified lemma), plus the
    favorable-order transposition bound output >= d(Wfac) + 23 - 9 - ... 
    No solver in the sweep: d-counts only, exact ±1 arithmetic.

  Total objective reported = dL + dR + dO where
    dL = d(X U Y^-1), dR = d(Y V Z^-1), dO = d(Z W X^-1) + (23 - 9)
  (the +23-9 = +14 term is the transposition gap for the output stage; every
  swept orientation's OUTPUT circuit cost >= d(W'-factor) + 14, matching the
  certified model used for the anchors: paper55 output 28 = d 14 + 14 ... per
  side bookkeeping below).

Anchor cross-check model (must reproduce 55 for paper55 under Id):
  - left:  C >= d(U) = 12, paper witness 13
  - right: C >= d(V) = 13, paper witness 14
  - out:   C >= d(Wfac) + 14 = 13 + 14 = 27, paper witness 28
  So the sweep total is dL + dR + dO + 14 and compares against 55/57/57.

Checkpointing: appends one JSON line per (decomposition, sigma-power) triple
index block to sweep_checkpoint.jsonl after every 4096 triples; crash-safe.
"""
import sys
import json
import time
from pathlib import Path
from itertools import product, permutations

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")

from gatec_decomps import LOADERS, META, R

CAMPAIGN = Path("/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/"
                "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56")
CKPT = CAMPAIGN / "sweep_checkpoint.jsonl"
LANDSCAPE = CAMPAIGN / "landscape.json"

N = 9


def canon(v):
    t = tuple(v)
    n = tuple(-x for x in t)
    return t if t <= n else n


def d_count(targets):
    """Distinct non-input sign-classes of 23 ternary target vectors in Z^9."""
    s = set()
    for t in targets:
        if any(t):
            s.add(canon(t))
    in_dirs = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    return len(s - in_dirs)


def build_sp_matrices():
    """48 3x3 signed permutation matrices as tuples of 3 tuples."""
    out = []
    for p in permutations(range(3)):
        for signs in product((1, -1), repeat=3):
            M = [[0] * 3 for _ in range(3)]
            for row in range(3):
                M[row][p[row]] = signs[row]
            out.append(tuple(tuple(r) for r in M))
    return tuple(out)


SP = build_sp_matrices()
assert len(SP) == 48


def sandwich(U, V, W, X, Y, Z):
    """Verified isotropy action on the factor triples. The tensor
    M333 = Σ e_(i,k) ⊗ e_(k,j) ⊗ e_(i,j) is fixed by the substitution
    A -> X A Y⁻¹, B -> Y B Z⁻¹, C -> ... acting on the factor linear forms as
      U'[r][(i2,k2)] = Σ_{i,k} (X⁻¹)[i][i2] · Y[k2][k] · U[r][(i,k)]
      V'[r][(k2,j2)] = Σ_{k,j} (Y⁻¹)[k][k2] · Z[j2][j] · V[r][(k,j)]
      W'[r][(i2,j2)] = Σ_{i,j} (X⁻¹)[i][i2] · Z[j2][j] · W[r][(i,j)]
    Verified: brent 729/729 over Z on ALL five anchored decompositions for
    sampled triples (campaign log; also the invariant check per sweep run)."""
    Xi, Yi = inv_sp(X), inv_sp(Y)
    U2, V2, W2 = [], [], []
    for r in range(R):
        u33 = [list(U[r][0:3]), list(U[r][3:6]), list(U[r][6:9])]
        u2 = [[sum(Xi[i][i2] * Y[k2][k] * u33[i][k]
                   for i in range(3) for k in range(3))
               for k2 in range(3)] for i2 in range(3)]
        U2.append(tuple(u2[0] + u2[1] + u2[2]))
        v33 = [list(V[r][0:3]), list(V[r][3:6]), list(V[r][6:9])]
        v2 = [[sum(Yi[k][k2] * Z[j2][j] * v33[k][j]
                   for k in range(3) for j in range(3))
               for j2 in range(3)] for k2 in range(3)]
        V2.append(tuple(v2[0] + v2[1] + v2[2]))
        w33 = [list(W[r][0:3]), list(W[r][3:6]), list(W[r][6:9])]
        w2 = [[sum(Xi[i][i2] * Z[j2][j] * w33[i][j]
                   for i in range(3) for j in range(3))
               for j2 in range(3)] for i2 in range(3)]
        W2.append(tuple(w2[0] + w2[1] + w2[2]))
    return U2, V2, W2



def inv_sp(M):
    """Inverse of a 3x3 signed permutation matrix."""
    out = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            if M[i][j]:
                out[j][i] = M[i][j]   # transpose of monomial gives inverse
    return out


def sigma_orbit(U, V, W):
    """{Id, sigma, sigma^2} with sigma: (U,V,W) -> (V, W^T, U^T). W-factor
    transpose: Wfac 23x9 -> 9x23 -> rows. We keep W in Wfac rows (23 x 9);
    the transposition acts by reshaping each product's 3x3 and transposing."""
    def T(blocks):
        out = []
        for r in range(R):
            b = [blocks[r][0:3], blocks[r][3:6], blocks[r][6:9]]
            bt = [list(x) for x in zip(*b)]
            out.append(tuple(bt[0] + bt[1] + bt[2]))
        return out
    o0 = (U, V, W)
    o1 = (V, T(W), T(U))          # sigma:  (U,V,W) -> (V, W^T, U^T)
    o2 = (T(W), U, T(V))          # sigma^2: apply sigma to o1 -> (W^T, U, V^T)
    return [o0, o1, o2]


def col_order_invariance_test(U, V, W, tries=8, seed=17):
    """Random column permutation of the 23 products must leave all d() fixed."""
    import random
    rng = random.Random(seed)
    base = (d_count(U), d_count(V), d_count(W))
    for _ in range(tries):
        pm = list(range(R))
        rng.shuffle(pm)
        U2 = [U[i] for i in pm]
        V2 = [V[i] for i in pm]
        W2 = [W[i] for i in pm]
        got = (d_count(U2), d_count(V2), d_count(W2))
        if got != base:
            return False, base, got
    return True, base, None


def sandwich_validity_sanity(U, V, W, X, Y, Z):
    """The sandwich action preserves the tensor mathematically; we verify once
    per (decomposition, sigma-power, one representative triple) with a full
    729/729 Brent check. Cheap enough and catches transcription errors."""
    U2, V2, W2 = sandwich(U, V, W, X, Y, Z)
    f = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ip in range(3):
                    for jp in range(3):
                        for kp in range(3):
                            s = 0
                            for r in range(R):
                                ur, vr, wr = U2[r][3*i+k], V2[r][3*kp+j], W2[r][3*ip+jp]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            expect = 1 if (i == ip and j == jp and k == kp) else 0
                            if s != expect:
                                f += 1
    return f


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else sorted(LOADERS)
    results = {name: {"hits_le_54": [], "min": None, "dist": {}} for name in only}
    t0 = time.time()
    ckpt_f = open(CKPT, "a")

    for name in only:
        U0, V0, W0 = LOADERS[name]()
        ok, base, got = col_order_invariance_test(U0, V0, W0)
        if not ok:
            print(f"COUNTER BUG: {name} column-order variance {base} vs {got}")
            sys.exit(1)
        print(f"{name}: column-order invariance OK, base d(L,R,O) = {base}")
        # one full Brent check per sigma power on a representative sandwich
        X0, Y0, Z0 = SP[0], SP[1], SP[2]
        for sp_idx, (U, V, W) in enumerate(sigma_orbit(U0, V0, W0)):
            f = sandwich_validity_sanity(U, V, W, X0, Y0, Z0)
            if f:
                print(f"VALIDITY FAIL {name} sigma^{sp_idx}: {f} brent failures on sandwich")
                sys.exit(1)
        print(f"{name}: sandwich validity OK (one triple per sigma power, 729/729)")

        dmin = None
        dist = {}
        hits = []
        cnt = 0
        for sp_idx, (U, V, W) in enumerate(sigma_orbit(U0, V0, W0)):
            for xi, X in enumerate(SP):
                for yi, Y in enumerate(SP):
                    for zi, Z in enumerate(SP):
                        U2, V2, W2 = sandwich(U, V, W, X, Y, Z)
                        dl, dr, do = d_count(U2), d_count(V2), d_count(W2)
                        total = dl + dr + do + 14   # + transposition gap
                        dist[total] = dist.get(total, 0) + 1
                        if dmin is None or total < dmin:
                            dmin = total
                            best = (sp_idx, xi, yi, zi, dl, dr, do)
                        if total <= 54:
                            hits.append((sp_idx, xi, yi, zi, dl, dr, do, total))
                        cnt += 1
                        if cnt % 4096 == 0:
                            ckpt_f.write(json.dumps(
                                {"decomp": name, "done": cnt, "minsofar": dmin,
                                 "t": round(time.time() - t0, 1)}) + "\n")
                            ckpt_f.flush()
        results[name] = {"min": dmin, "best": best, "dist": dist, "hits_le_54": hits}
        print(f"{name}: SWEEP DONE min_total={dmin} best={best} hits<=54={len(hits)} "
              f"({time.time()-t0:.0f}s)")
        ckpt_f.write(json.dumps({"decomp": name, "done": "FINAL", "min": dmin,
                                 "dist": dist}) + "\n")
        ckpt_f.flush()
        LANDSCAPE.write_text(json.dumps(results, indent=1))

    ckpt_f.close()
    print("ALL SWEEPS COMPLETE", json.dumps({k: v["min"] for k, v in results.items()}))


if __name__ == "__main__":
    main()
