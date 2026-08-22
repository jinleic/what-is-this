#!/usr/bin/env python3
"""Search for ARRANGEMENT-LEVEL convention separation.

For an arrangement A of n distinct affine lines define

  face(A)  = max size of a pairwise interior-disjoint family of triangles
             whose open interiors meet no line   (= triangular faces)
  broad(A) = max size of a pairwise interior-disjoint family of triangles
             (crossings allowed)

Always face(A) <= broad(A).  A strict gap would show the broad convention is
strictly stronger as a functional on arrangements, and is the mechanism by
which Kgen(n) > K(n) could ever happen.  This sweeps exact-rational
arrangements -- generic, with parallel classes, and with planted multiple
points -- and reports every gap it finds.

Usage:  convention_sweep.py N TRIALS [SEED]
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


def hom(p):
    x, y = p
    d = 1
    for v in (x, y):
        d = d * v.denominator // math.gcd(d, v.denominator)
    return (int(x * d), int(y * d), d)


def mis(nodes, conf):
    """Exact maximum independent set with a simple bound."""
    nodes = sorted(nodes, key=lambda t: -len(conf[t]))
    best = [0, None]

    def rec(i, cur, banned):
        if len(cur) + len(nodes) - i <= best[0]:
            return
        if i == len(nodes):
            if len(cur) > best[0]:
                best[0], best[1] = len(cur), tuple(cur)
            return
        t = nodes[i]
        if t not in banned:
            rec(i + 1, cur + [t], banned | conf[t])
        rec(i + 1, cur, banned)

    rec(0, [], set())
    return best[0], best[1]


def analyse(n, ms, bs):
    Xp = engine.crossings(n, ms, bs)
    tris = [t for t in combinations(range(n), 3) if engine.tri_ok(Xp, t)]
    if not tris:
        return None
    H = {t: tuple(hom(v) for v in engine.tri_verts(Xp, t)) for t in tris}
    crossed = {}
    for t in tris:
        vs = engine.tri_verts(Xp, t)
        c = []
        for r in range(n):
            if r in t:
                continue
            s = [bs[r] + ms[r] * x - y for (x, y) in vs]
            if any(v > 0 for v in s) and any(v < 0 for v in s):
                c.append(r)
        crossed[t] = c
    conf = {t: set() for t in tris}
    for a, b in combinations(tris, 2):
        if not engine.interiors_disjoint_h(H[a], H[b]):
            conf[a].add(b)
            conf[b].add(a)
    broad, bfam = mis(tris, conf)
    faces = [t for t in tris if not crossed[t]]
    face, ffam = mis(faces, conf) if faces else (0, ())
    return {"broad": broad, "face": face, "broad_family": bfam,
            "face_family": ffam, "n_tri": len(tris), "n_face": len(faces)}


def sample(n, rng, mode):
    ms = [Fr(rng.randint(-30, 30), rng.randint(1, 7)) for _ in range(n)]
    if mode == "parallel" and n >= 4:            # one parallel class of size 2
        ms[1] = ms[0]
    if mode == "parallel3" and n >= 6:           # a class of size 3
        ms[1] = ms[0]
        ms[2] = ms[0]
    if len(set(ms)) < n - (2 if mode == "parallel3" else
                           1 if mode == "parallel" else 0):
        return None
    ms = sorted(ms)
    bs = [Fr(rng.randint(-30, 30), rng.randint(1, 7)) for _ in range(n)]
    if len(set(zip(ms, bs))) < n:
        return None
    if mode in ("triple", "quad", "two_triples"):
        px = Fr(rng.randint(-6, 6), rng.randint(1, 3))
        py = Fr(rng.randint(-6, 6), rng.randint(1, 3))
        grp = {"triple": 3, "quad": 4, "two_triples": 3}[mode]
        idx = rng.sample(range(n), grp)
        for r in idx:
            bs[r] = py - ms[r] * px
        if mode == "two_triples":
            qx = Fr(rng.randint(-6, 6), rng.randint(1, 3))
            qy = Fr(rng.randint(-6, 6), rng.randint(1, 3))
            rest = [r for r in range(n) if r not in idx]
            if len(rest) >= 3:
                for r in rng.sample(rest, 3):
                    bs[r] = qy - ms[r] * qx
    if len(set(zip(ms, bs))) < n:
        return None
    return ms, bs


MODES = ("generic", "parallel", "parallel3", "triple", "quad", "two_triples")


def main() -> int:
    n = int(sys.argv[1])
    trials = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = random.Random(seed)
    gaps, tally, done = [], {}, 0
    for i in range(trials):
        mode = MODES[i % len(MODES)]
        s = sample(n, rng, mode)
        if s is None:
            continue
        res = analyse(n, *s)
        if res is None:
            continue
        done += 1
        key = (mode, res["face"], res["broad"])
        tally[key] = tally.get(key, 0) + 1
        if res["broad"] > res["face"]:
            ms, bs = s
            gaps.append({"mode": mode, "face": res["face"],
                         "broad": res["broad"],
                         "lines_frac": [[str(m), str(b)] for m, b in zip(ms, bs)],
                         "broad_family": [list(t) for t in res["broad_family"]],
                         "face_family": [list(t) for t in res["face_family"]]})
            print(json.dumps({"stage": "GAP", "n": n, **gaps[-1]}), flush=True)
            if len(gaps) >= 5:
                break
    out = {"stage": "SUMMARY", "n": n, "arrangements": done,
           "gaps": len(gaps),
           "by_mode_face_broad": {f"{m}:{f}->{b}": c
                                  for (m, f, b), c in sorted(tally.items())}}
    print(json.dumps(out), flush=True)
    if gaps:
        p = ROOT / "scratch" / "kobon" / "discoveries" / f"convention_gap_n{n}.json"
        p.write_text(json.dumps(gaps, indent=1))
        print(json.dumps({"stage": "written", "path": str(p)}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
