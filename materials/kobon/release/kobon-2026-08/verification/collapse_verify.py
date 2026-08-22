#!/usr/bin/env python3
"""Machine-check the corner-descent lemma and the convention collapse.

Lemma (corner descent).  Every triangle T bounded by three lines of an
arrangement contains a triangular face of the arrangement.

Theorem (collapse).  face(A) = broad(A) for every arrangement A, hence
Kgen(n) = K(n).

This checks, on exact-rational arrangements with planted degeneracies:

  (L1)  every valid triangle contains a triangular FACE (verified by
        actually running the descent: repeatedly cut by a crossing line and
        keep a triangular piece, and confirm the endpoint is an uncrossed
        triangle contained in the original);
  (L2)  the descent's containment is exact (piece interior inside T);
  (T1)  face(A) equals broad(A), where face(A) is the number of triangular
        faces and broad(A) is the exact maximum interior-disjoint family
        with crossings allowed.

Any failure prints a counterexample and exits nonzero.

Usage:  collapse_verify.py N TRIALS [SEED]
"""
from __future__ import annotations

import json
import math
import random
import sys
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "math" / "kobon"))
import engine  # noqa: E402

sys.path.insert(0, str(ROOT / "scratch" / "kobon"))
from convention_sweep import MODES, hom, mis, sample  # noqa: E402


def side(ms, bs, r, p):
    """Sign of line r at point p: >0 above-line value, i.e. b+mx-y."""
    x, y = p
    return bs[r] + ms[r] * x - y


def crossing_lines(n, ms, bs, t, verts):
    out = []
    for r in range(n):
        if r in t:
            continue
        s = [side(ms, bs, r, v) for v in verts]
        if any(v > 0 for v in s) and any(v < 0 for v in s):
            out.append(r)
    return out


def descend(n, ms, bs, Xp, t):
    """Run the corner descent from triangle t.  Returns (face_triple, path)."""
    path = [t]
    cur = t
    guard = 0
    while True:
        guard += 1
        if guard > 3 * n + 5:
            raise AssertionError(f"descent did not terminate from {t}")
        verts = engine.tri_verts(Xp, cur)
        cross = crossing_lines(n, ms, bs, cur, verts)
        if not cross:
            return cur, path
        c = cross[0]
        # candidate pieces: replace one line of cur by c
        cands = []
        for drop in cur:
            piece = tuple(sorted(set(cur) - {drop} | {c}))
            if len(piece) != 3 or not engine.tri_ok(Xp, piece):
                continue
            pv = engine.tri_verts(Xp, piece)
            # piece must sit inside cur: every vertex weakly inside, and the
            # piece must be nondegenerate
            inside = True
            for r in cur:
                if r in piece:
                    continue
            # weak containment test: each piece vertex on the correct side of
            # every side-line of cur
            cv = verts
            for r in cur:
                ref = [side(ms, bs, r, v) for v in cv if side(ms, bs, r, v) != 0]
                if not ref:
                    continue
                sgn = 1 if ref[0] > 0 else -1
                for v in pv:
                    s = side(ms, bs, r, v)
                    if s != 0 and (1 if s > 0 else -1) != sgn:
                        inside = False
                        break
                if not inside:
                    break
            if not inside:
                continue
            nxt = crossing_lines(n, ms, bs, piece, pv)
            if c in nxt:
                continue
            if not set(nxt) <= set(cross) - {c}:
                continue
            cands.append(piece)
        if not cands:
            raise AssertionError(
                f"no triangular piece for {cur} cut by {c} (from {t})")
        cur = cands[0]
        path.append(cur)


def main() -> int:
    n = int(sys.argv[1])
    trials = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = random.Random(seed)
    checked_tri = checked_arr = 0
    fails = []
    depth_hist = {}
    for i in range(trials):
        s = sample(n, rng, MODES[i % len(MODES)])
        if s is None:
            continue
        ms, bs = s
        Xp = engine.crossings(n, ms, bs)
        tris = [t for t in combinations(range(n), 3) if engine.tri_ok(Xp, t)]
        if not tris:
            continue
        checked_arr += 1
        H = {t: tuple(hom(v) for v in engine.tri_verts(Xp, t)) for t in tris}
        faces = []
        for t in tris:
            verts = engine.tri_verts(Xp, t)
            if not crossing_lines(n, ms, bs, t, verts):
                faces.append(t)
            # (L1)+(L2): descend and check the endpoint
            try:
                f, path = descend(n, ms, bs, Xp, t)
            except AssertionError as e:
                fails.append(("descent", str(e)))
                continue
            checked_tri += 1
            depth_hist[len(path) - 1] = depth_hist.get(len(path) - 1, 0) + 1
            if crossing_lines(n, ms, bs, f, engine.tri_verts(Xp, f)):
                fails.append(("endpoint-crossed", str((t, f))))
            if f != t and engine.interiors_disjoint_h(H[t], H[f]):
                fails.append(("not-contained", str((t, f))))
        conf = {t: set() for t in tris}
        for a, b in combinations(tris, 2):
            if not engine.interiors_disjoint_h(H[a], H[b]):
                conf[a].add(b)
                conf[b].add(a)
        broad, _ = mis(tris, conf)
        # faces are pairwise disjoint, so face(A) = number of triangular faces
        fcount = len(faces)
        fmax, _ = mis(faces, conf) if faces else (0, ())
        if fmax != fcount:
            fails.append(("faces-conflict", str((fcount, fmax))))
        if broad != fcount:
            fails.append(("GAP", json.dumps(
                {"face": fcount, "broad": broad,
                 "lines": [[str(m), str(b)] for m, b in zip(ms, bs)]})))
    print(json.dumps({"n": n, "arrangements": checked_arr,
                      "triangles_descended": checked_tri,
                      "descent_depths": dict(sorted(depth_hist.items())),
                      "failures": len(fails)}))
    for kind, detail in fails[:5]:
        print("FAIL", kind, detail)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
