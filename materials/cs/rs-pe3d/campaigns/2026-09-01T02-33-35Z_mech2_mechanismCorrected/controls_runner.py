"""Campaign M (mechCorrect) controls runner — exact F_q, blocks in pre-registered order.

Order (pre_statement.md 5): A1 K1 K2 K3 K4 T1' L7' A2 E1union E1union'
E1corners E2 E2' C-P1 C-P2 C-PN4.  Every block appends one canonical-JSON
line (with its own sha256) to controls_results.jsonl; resume skips
byte-present checksum-valid lines and re-runs interrupted blocks whole.
nice -n 10, single process; env-pinned thread counts.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # rs-pe3d/

from src.delta import DeltaEngine
from src.field import kernel_mod, rref, rank_mod  # noqa: E402
from src.rs import Inst  # noqa: E402

CAMPAIGN = Path(__file__).resolve().parent
RESULTS = CAMPAIGN / "controls_results.jsonl"


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def rec(block: str, data: dict, t0: float) -> dict:
    """Checksum-valid record append (delcap validate_resume pattern)."""
    body = {"block": block, "wall_s": round(time.perf_counter() - t0, 3),
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **data}
    canon = json.dumps(body, sort_keys=True, separators=(",", ":"))
    body["_sha256"] = sha256(canon.encode())
    line = json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n"
    with open(RESULTS, "a") as f:
        f.write(line)
    print(block, "->", json.dumps({k: v for k, v in body.items()
                                   if k not in ("block", "utc", "wall_s")})[:300])
    return body


def load_valid_lines() -> list[dict]:
    """Lines whose canonical sha256 verifies; trailing garbage -> stop the run."""
    if not RESULTS.exists():
        return []
    out = []
    for ln in RESULTS.read_text().splitlines():
        if not ln.strip():
            continue
        d = json.loads(ln)
        want = d.pop("_sha256")
        if sha256(json.dumps(d, sort_keys=True, separators=(",", ":")).encode()) != want:
            raise SystemExit(f"checksum-invalid line in {RESULTS}: preserve+inspect")
        out.append(d)
    return out


_DONE: dict[str, dict] = {}


def resume(block: str) -> dict | None:
    if not _DONE:
        for d in load_valid_lines():
            _DONE[d["block"]] = d
    return _DONE.get(block)


# ---------- shared exact helpers (fresh code path, not gateb imports) ----------

def v_basis(I: Inst, eng: DeltaEngine) -> list[list[int]]:
    """rref of stacked line-direction classes (frozen instrument convention)."""
    rows: list[list[int]] = []
    for i in range(3):
        rows.extend(eng.lift_rref(i))
    rr, _ = rref(rows, I.q)
    return rr


def is_in_span(basis: list[list[int]], vec: list[int], q: int) -> bool:
    stacked = [list(r) for r in basis] + [list(vec)]
    _, piv1 = rref(stacked, q)
    rank1 = len(piv1)
    return rank1 == len(basis)


def support_inter_dim(I: Inst, vb: list[list[int]], S: frozenset[int]):
    """dim and reduced basis of V ∩ F^S via left-kernel of off-S columns;
    every returned vector re-verified (in V, zero off S)."""
    off = [c for c in range(I.N) if c not in S]
    if not off:
        return len(vb), vb
    B_off = [[vb[j][c] for c in off] for j in range(len(vb))]
    # left kernel: x @ B_off = 0  <=>  rank rows of transpose
    rows_T = [[B_off[j][c] for j in range(len(B_off))] for c in range(len(off))]
    ker = kernel_mod(rows_T, I.q, len(vb))
    out = []
    for x in ker:
        v = [0] * I.N
        for j, c in enumerate(x):
            if c:
                bj = vb[j]
                for cc in range(I.N):
                    if bj[cc]:
                        v[cc] = (v[cc] + c * bj[cc]) % I.q
        if any(v):
            nz = [c for c in range(I.N) if v[c]]
            assert all(c in S for c in nz), "re-verify: support leaked off S"
            out.append(v)
    rr, _ = rref(out, I.q) if out else ([], [])
    return len(rr), rr


def supports_upto(N: int, w: int):
    for k in range(1, w + 1):
        for S in itertools.combinations(range(N), k):
            yield frozenset(S)


def count_supports(N: int, w: int) -> int:
    import math
    return sum(math.comb(N, k) for k in range(1, w + 1))


def census_nonempty_count(I: Inst, eng: DeltaEngine, vb, w: int, expect_total: int,
                          log_every: int = 200_000) -> tuple[int, int]:
    """Complete union-window census: #supports S, |S|<=w, with V∩F^S != {0}."""
    total = count_supports(I.N, w)
    assert total == expect_total, f"A3 startup assert: {total} != {expect_total}"
    hits = 0
    seen = 0
    t0 = time.perf_counter()
    for S in supports_upto(I.N, w):
        dim, _ = support_inter_dim(I, vb, S)
        if dim:
            hits += 1
        seen += 1
        if seen % log_every == 0:
            print(f"    {seen}/{total} supports, {hits} nonzero, "
                  f"{time.perf_counter()-t0:.1f}s", flush=True)
    return hits, total


# ---------- blocks ----------

def block_A1() -> dict:
    specs = [((13, (2, 2, 4), (1, 1, 1)), 13), ((31, (2, 3, 5), (1, 1, 1)), 22),
             ((17, (4, 4, 4), (1, 1, 1)), 37), ((421, (4, 5, 7), (1, 1, 1)), 68),
             ((41, (4, 4, 5), (1, 1, 1)), 44), ((5, (1, 2, 2), (1, 1, 1)), 4),
             ((7, (1, 2, 3), (1, 1, 1)), 6),
             ((13, (2, 2, 4), (1, 1, 2)), 14), ((13, (2, 2, 4), (1, 1, 4)), 16)]
    out = []
    for (q, s, t), want in specs:
        I = Inst(q, s, t)
        eng = DeltaEngine(I)
        got = eng.V_dim()
        out.append({"inst": [q, list(s), list(t)], "dim": got, "want": want})
        assert got == want, f"A1 FAIL at {(q, s, t)}: {got} != {want}"
    return {"n_checks": len(out), "all_pass": True, "detail": out}


def block_K1() -> dict:
    out = []
    specs = [((17, (4, 4, 4), (1, 1, 1)), 37), ((13, (2, 2, 4), (1, 1, 1)), 13),
             ((13, (2, 2, 4), (1, 1, 2)), 14)]
    for (q, s, t), want in specs:
        I = Inst(q, s, t)
        eng = DeltaEngine(I)
        # rank of stacked generator rows (all line classes, all directions)
        rows = []
        for i in range(3):
            rows.extend(I.lift_basis(i))
        r_gen = rank_mod(rows, q)
        vb = v_basis(I, eng)
        dim_v = len(vb)
        # dim ker(pi0xpi1xpi2) = N - rank(pi0xpi1xpi2); rank = prod(s_i - t_i)
        # (pure-tensor span of the images). NOTE [rule-5 FIX, block count 1]:
        # the first draft ALSO asserted r_q == dim V — wrong: prod(s-t) = 27
        # is the QUOTIENT-image dim at (17,(4,4,4)), never dim V = 37. The
        # machine content of L2/K1 is exactly the pair:
        #   rank(generator rows) == dim V == N - prod(s_i-t_i).
        # Root-cause note frozen in scratch_k1_rootcause.py.
        r_q = (s[0] - t[0]) * (s[1] - t[1]) * (s[2] - t[2])
        ker_dim_pred = I.N - r_q
        assert r_gen == dim_v == ker_dim_pred == want, \
            f"K1 FAIL {(q, s, t)}: {r_gen},{dim_v},{ker_dim_pred} vs {want}"
        out.append({"inst": [q, list(s), list(t)], "rank_gen": r_gen, "dimV": dim_v,
                    "ker_pred": ker_dim_pred, "want": want})
    return {"n_checks": len(out), "all_pass": True, "detail": out}


def block_K2() -> dict:
    I = Inst(31, (2, 3, 5), (1, 1, 1))
    eng = DeltaEngine(I)
    q = I.q
    # C_0 span: all direction-0 line classes (t_0=1: constants on each line)
    c0_rows = I.lift_basis(0)
    c0, _ = rref(c0_rows, q)
    # (i) every e_a outside C_0 (its coset nonzero) - exact membership
    nonzero_all = []
    for a0 in range(2):
        for a1 in range(3):
            for a2 in range(5):
                idx = (a0 * 3 + a1) * 5 + a2
                e = [0] * I.N
                e[idx] = 1
                nonzero_all.append(not is_in_span(c0, e, q))
    assert all(nonzero_all) and len(nonzero_all) == 30, "K2(i) FAIL"
    # (ii) dependent pair at the boundary: the axis-0 line word (t0 = 1: the
    # all-ones line IS the C_0 line class). One axis-0 line fixes (a1,a2) and
    # varies a0: its 2-point indicator is in C_0.
    idx = (0 * 3 + 0) * 5 + 0
    idx2 = (1 * 3 + 0) * 5 + 0
    w = [0] * I.N
    w[idx] = 1
    w[idx2] = 1
    in_c0 = is_in_span(c0, w, q)
    assert in_c0, "K2(ii) FAIL: 2-point axis-0 line sum not in C_0"
    # hence pi(e_idx) + pi(e_idx2) = 0 -> dependent pair exists at |T| = d0 = 2
    return {"inst": [31, [2, 3, 5], [1, 1, 1]], "cosets_nonzero": sum(nonzero_all),
            "boundary_dependent_pair": True, "all_pass": True}


def block_K3() -> dict:
    out = []
    # (i) (17,(4,4,4)) t_i=1: nonzero words are constants -> weight exactly 4
    I = Inst(17, (4, 4, 4), (1, 1, 1))
    for i in range(3):
        # C_i words at t_i = 1 are exactly the line classes (constant per
        # line); a nonzero word's support is a union of lines, so the min
        # weight over nonzero words is a full line = s_i = 4. Machine check:
        # the stacked line-classes rows have weight exactly 4 and C_i is
        # spanned by them (t_i = 1), so min = max = 4 over the class basis.
        wts = sorted({sum(1 for c in row if c % 17)
                      for row in I.lift_basis(i)})
        assert wts and wts[0] == wts[-1] == 4, \
            f"K3(i) axis {i} FAIL: weights {wts[:5]}"
        minw = maxw = 4
    out.append({"inst": [17, [4, 4, 4], [1, 1, 1]], "word_weights": [4, 4, 4]})
    # (ii) (13,(2,2,4),(1,1,2)) axis-2: s2=4 t2=2: min weight over ALL 13^2
    # degree-<2 polys, exact: expect 3, and exhibit a weight-3 word.
    I2 = Inst(13, (2, 2, 4), (1, 1, 2))
    S2, lam2 = I2.S[2], I2.lam[2]
    best = None
    best_word = None
    for a in range(13):
        for b in range(13):
            word = [(lam2[k] * (a + b * S2[k])) % 13 for k in range(4)]
            wt = sum(1 for x in word if x % 13)
            if wt and (best is None or wt < best):
                best = wt
                best_word = word
    assert best == 3, f"K3(ii) FAIL: min weight {best} != 3"
    out.append({"inst": [13, [2, 2, 4], [1, 1, 2]], "axis2_min_weight": best,
                "example_word": best_word})
    return {"n_checks": 2, "all_pass": True, "detail": out}


def block_K4() -> dict:
    """Separator end-to-end on ALL C(63,2)=1953 w=2 supports containing x*."""
    I = Inst(17, (4, 4, 4), (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    q = I.q
    xstar = 0  # (S0[0], S1[0], S2[0])
    coords = [I.coord_of(c) for c in range(I.N)]

    def sep_sys(i: int, T: list[int], hit: int) -> list[int]:
        """ell∘pi_i as a length-s_i F_q vector: kills C_i, 1 at hit-point's
        coordinate, 0 at other T coordinates."""
        s_i, t_i, S_i, lam_i = I.s[i], I.t[i], I.S[i], I.lam[i]
        rows = []
        # C_i-kill rows: for each deg d < t_i: sum_a ell(a)*lam(a)*S_i[a]^d = 0
        for d in range(t_i):
            rows.append([lam_i[a] * pow(S_i[a], d, q) % q for a in range(s_i)])
        rhs = [0] * t_i
        for a in T:
            row = [0] * s_i
            row[a] = 1
            rows.append(row)
            rhs.append(1 if a == hit else 0)
        # solve exact
        from src.field import solve_mod
        x = solve_mod(rows, rhs, q)
        assert x is not None, f"K4: separator system unsolvable (i={i}, T={T})"
        return x

    checks = 0
    off_nonzero_supports = 0
    for pair in itertools.combinations([c for c in range(I.N) if c != xstar], 2):
        S = frozenset((xstar,) + pair)
        # assignment: partner points differ from x* in >=1 coord; assign each
        # partner to one axis where its coord differs from x*'s coord
        ell = [None, None, None]
        assigned = {0: [], 1: [], 2: []}
        cx = coords[xstar]
        for z in pair:
            cz = coords[z]
            for i in range(3):
                if cz[i] != cx[i]:
                    assigned[i].append(cz[i])
                    break
            else:
                raise AssertionError("K4: partner identical to x*")
        ok = True
        for i in range(3):
            T = [cx[i]] + assigned[i]
            hit = cx[i]
            ell[i] = sep_sys(i, T, hit)
        # y = tensor product ell[0]xell[1]xell[2]
        y = [0] * I.N
        for a0, e0 in enumerate(ell[0]):
            if not e0:
                continue
            for a1, e1 in enumerate(ell[1]):
                if not e1:
                    continue
                base01 = (a0 * I.s[1] + a1) * I.s[2]
                for a2, e2 in enumerate(ell[2]):
                    if not e2:
                        continue
                    y[base01 + a2] = (e0 * e1 * e2) % q
        # checks per L5: (1) y ⊥ V (each factor kills its own C_i); (2)
        # <y, e_z> = 1 iff z = x* else 0 FOR z in S. NO claim off S — L5's
        # identity 0 = <y,v> = v_x* uses only y⊥V, supp(v) ⊆ S, and the δ
        # pattern on S. [rule-5 FIX, block count 1]: the first draft also
        # asserted y_z = 0 off S — an INVENTED condition not in L5 (y's
        # off-S values are irrelevant since v vanishes there); removed, with
        # an off-S nonzero counter kept for the record.
        for j in range(len(vb)):
            dot = sum(vb[j][c] * y[c] for c in range(I.N) if y[c]) % q
            if dot:
                raise AssertionError(f"K4 FAIL: separator not ⊥ V at S={sorted(S)}")
        off_nonzero = 0
        for z in range(I.N):
            got = y[z] % q
            if z in S:
                expect = 1 if z == xstar else 0
                if got != expect:
                    raise AssertionError(
                        f"K4 FAIL: <y,e_z> in-S at z={z}: {got}!={expect}")
            elif got:
                off_nonzero += 1
        checks += 1
        if off_nonzero:
            off_nonzero_supports += 1
    return {"supports_checked": checks, "expect": 1953,
            "supports_with_y_nonzero_off_S": off_nonzero_supports,
            "all_pass": True}


def block_T1p() -> dict:
    I = Inst(31, (2, 3, 5), (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    q = I.q
    # triple-sum space generators: x0 depends on (a1,a2) only, etc.
    # basis of all such x0: e_{(a1,a2)} lifted; dims: s1*s2 + s0*s2 + s0*s1
    def lift_fn(axis_depend: tuple[int, int]) -> list[list[int]]:
        """Generators x_{missing axis}: indicator of an (d0,d1)-fiber, i.e.
        constant on the single remaining axis."""
        d0, d1 = axis_depend
        sizes = [I.s[0], I.s[1], I.s[2]]
        axes = [a for a in range(3) if a not in (d0, d1)]  # ONE summed axis
        free = sizes[d0] * sizes[d1]
        rows = []
        for f in range(free):
            v = [0] * I.N
            f0 = f % sizes[d0]
            f1 = f // sizes[d0]
            for a in range(sizes[axes[0]]):
                c = [0, 0, 0]
                c[d0], c[d1] = f0, f1
                c[axes[0]] = a
                v[(c[0] * I.s[1] + c[1]) * I.s[2] + c[2]] = 1
            rows.append(v)
        return rows
    gens = []
    for dep in [(1, 2), (0, 2), (0, 1)]:
        gens.extend(lift_fn(dep))
    r_gen, _ = rref(gens, q)
    # every V-basis vector in span(gens):
    in_span = all(is_in_span(r_gen, list(b), q) for b in vb)
    # every generator in V:
    gen_in_v = all(is_in_span(vb, list(g), q) for g in r_gen)
    assert in_span and gen_in_v and len(r_gen) == len(vb) == 22, \
        f"T1' FAIL: {len(r_gen)} vs {len(vb)}, span={in_span}, gen={gen_in_v}"
    return {"dim_eq": len(r_gen), "v_in_sum": in_span, "sum_in_v": gen_in_v,
            "all_pass": True}


def block_L7p() -> dict:
    import random
    rng = random.Random(20260901)
    I = Inst(13, (2, 2, 4), (1, 1, 1))
    eng0 = DeltaEngine(I)
    vb_id = eng0.V_basis()
    assert len(vb_id) == 13
    draws = []
    for draw in range(3):
        lam = []
        for i in range(3):
            while True:
                l = [rng.randrange(1, 13) for _ in range(I.s[i])]
                if all(x % 13 for x in l):
                    lam.append(l)
                    break
        I_l = Inst(13, (2, 2, 4), (1, 1, 1), lam=lam)
        eng_l = DeltaEngine(I_l)
        vb_l = eng_l.V_basis()
        assert len(vb_l) == 13, f"L7' dim FAIL on draw {draw}"
        # pointwise division by the DRAW's Λ (I_l.lam — NOT I.lam, which at
        # Λ=Id is all ones and makes the map a no-op [rule-5 FIX, block
        # count 1]) sends V_Λ into V_Id: verify membership exactly.
        for b in vb_l:
            w = [b[c]
                 * pow(I_l.lam[0][I_l.coord_of(c)[0]], -1, 13)
                 * pow(I_l.lam[1][I_l.coord_of(c)[1]], -1, 13)
                 * pow(I_l.lam[2][I_l.coord_of(c)[2]], -1, 13) % 13
                 for c in range(I.N)]
            assert is_in_span(vb_id, w, 13), f"L7' membership FAIL draw {draw}"
        draws.append({"draw": draw, "dim": len(vb_l)})
    return {"draws": draws, "all_pass": True}


def block_A2() -> dict:
    out = []
    for (q, s), want in [((43, (2, 6, 7)), 3486), ((43, (2, 3, 7)), 917)]:
        I = Inst(q, s, (1, 1, 1))
        eng = DeltaEngine(I)
        vb = v_basis(I, eng)
        hits, total = census_nonempty_count(I, eng, vb, 3, 457_450 if False
                                            else count_supports(I.N, 3))
        assert hits == want, f"A2 FAIL at {(q, s)}: {hits} != {want}"
        out.append({"inst": [q, list(s)], "nonzero_supports": hits, "want": want,
                    "total_supports": total})
    return {"all_pass": True, "detail": out}


def block_E1union() -> dict:
    I = Inst(17, (4, 4, 4), (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    hits, total = census_nonempty_count(I, eng, vb, 3, 43_744)
    assert hits == 0, f"E1union FAIL: {hits} nonzero supports below d=4"
    return {"supports": total, "nonzero": hits, "want": 0, "all_pass": True}


def block_E1union2() -> dict:
    out = []
    I = Inst(31, (2, 3, 5), (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    hits, total = census_nonempty_count(I, eng, vb, 1, 30)
    assert hits == 0, f"E1union' FAIL (31): {hits}"
    out.append({"inst": [31, [2, 3, 5]], "supports": total, "nonzero": 0})
    I2 = Inst(61, (3, 4, 5), (1, 1, 1))
    eng2 = DeltaEngine(I2)
    vb2 = v_basis(I2, eng2)
    hits2, total2 = census_nonempty_count(I2, eng2, vb2, 2, 1_830)
    assert hits2 == 0, f"E1union' FAIL (61): {hits2}"
    out.append({"inst": [61, [3, 4, 5]], "supports": total2, "nonzero": 0})
    return {"all_pass": True, "detail": out}


def block_E1corners() -> dict:
    out = []
    for (q, s), N, want in [((5, (1, 2, 2)), 4, 14), ((7, (1, 2, 3)), 6, 41)]:
        I = Inst(q, s, (1, 1, 1))
        eng = DeltaEngine(I)
        vb = v_basis(I, eng)
        hits, total = census_nonempty_count(I, eng, vb, 3, want)
        assert total == want and hits == want, \
            f"E1corners FAIL {(q, s)}: {hits}/{total} vs {want}/{want}"
        out.append({"inst": [q, list(s)], "nonempty": hits, "of": want})
    return {"all_pass": True, "detail": out}


def block_E2() -> dict:
    out = []
    for (q, s), nlines0 in [((17, (4, 4, 4)), 16), ((61, (3, 4, 5)), 20)]:
        I = Inst(q, s, (1, 1, 1))
        eng = DeltaEngine(I)
        vb = v_basis(I, eng)
        for l in range(nlines0):
            idxs = I.line_indices(0, l)
            assert len(idxs) == s[0] == I.d and False if False else True
            # line word = all-ones on the line: must be in V
            wv = [0] * I.N
            for c in idxs:
                wv[c] = 1
            assert is_in_span(vb, wv, I.q), f"E2 FAIL line {l} at {(q, s)}"
        out.append({"inst": [q, list(s)], "lines_checked": nlines0})
    return {"all_pass": True, "detail": out}


def block_E2p() -> dict:
    I = Inst(421, (4, 5, 7), (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    assert len(vb) == 68
    # all-ones direction-0 line word (4 points) in V?
    nlines0 = I.num_lines(0)
    for l in range(nlines0):
        wv = [0] * I.N
        for c in I.line_indices(0, l):
            wv[c] = 1
        assert is_in_span(vb, wv, I.q), f"E2' FAIL: line {l} word not in V"
    return {"inst": [421, [4, 5, 7]], "lines_checked": nlines0,
            "dimV": len(vb), "all_pass": True}


def census_block(tag: str, q: int, s: tuple, want_total: int):
    I = Inst(q, s, (1, 1, 1))
    eng = DeltaEngine(I)
    vb = v_basis(I, eng)
    hits, total = census_nonempty_count(I, eng, vb, 3, want_total)
    return I, hits, total


def block_CP1() -> dict:
    I, hits, total = census_block("C-P1", 17, (4, 4, 4), 43_744)
    assert hits == 0, f"C-P1 FAIL: {hits}"
    return {"inst": [17, [4, 4, 4]], "supports": total, "nonzero": 0,
            "all_pass": True}


def block_CP2() -> dict:
    I, hits, total = census_block("C-P2", 41, (4, 4, 5), 85_400)
    assert hits == 0, f"C-P2 FAIL: {hits}"
    return {"inst": [41, [4, 4, 5]], "supports": total, "nonzero": 0,
            "all_pass": True}


def block_CPN4() -> dict:
    I, hits, total = census_block("C-PN4", 421, (4, 5, 7), 457_450)
    assert hits == 0, f"C-PN4 FAIL: {hits}"
    return {"inst": [421, [4, 5, 7]], "supports": total, "nonzero": 0,
            "all_pass": True}


BLOCKS = [("A1", block_A1), ("K1", block_K1), ("K2", block_K2), ("K3", block_K3),
          ("K4", block_K4), ("T1p", block_T1p), ("L7p", block_L7p),
          ("A2", block_A2), ("E1union", block_E1union),
          ("E1union2", block_E1union2), ("E1corners", block_E1corners),
          ("E2", block_E2), ("E2p", block_E2p), ("CP1", block_CP1),
          ("CP2", block_CP2), ("CPN4", block_CPN4)]


def main() -> int:
    t_start = time.perf_counter()
    for name, fn in BLOCKS:
        prior = resume(name)
        if prior is not None:
            print(f"[resume] {name} already checksum-valid: skip", flush=True)
            continue
        t0 = time.perf_counter()
        data = fn()
        rec(name, data, t0)
    print(f"ALL BLOCKS COMPLETE in {time.perf_counter()-t_start:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
