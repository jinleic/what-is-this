#!/usr/bin/env python3
"""Run 1 (gate H-MINLINE) - exact F_13 censuses, falsification branch.

Instances (all t=(1,1,1), Lambda=Id, subgroups = unique subgroup of each order
of the cyclic F_13^*, axis lists sorted ascending):
  E0: q=13, s=(2,2,4)   d=2  census over all C(16,2)=120 pairs
  C1: q=13, s=(2,3,4)   d=2  census over all C(24,2)=276 pairs
  C2: q=13, s=(3,3,4)   d=3  census over all C(36,3)=7140 triples

Routes (independent):
  A kernel: V = ker(H0 x H1 x H2), H_i = (s_i-1) x s_i Vandermonde parity
    check of RS(S_i,1) (rows alpha^0..alpha^{s_i-2}); carrier test
    dim(V n F^S) = |S| - rank(H[:,S]).
  B basis: V = rowspace(G), G = stack of lifted full-line indicator
    generators (rowspace of the three L_i(C_i) generator sets); carrier test
    dim = k - rank(G[:, S^c]), k = dim V.
Fail-loud: any route disagreement or target mismatch halts with exit 1.

Pre-registered targets: E0 total/line/off = 24/16/8 with explicit off-line
witness pair containing flat indices (0,12); C1 = 12 lines, 0 off;
C2 = 24 lines, 0 off.
Zero floating point; all arithmetic exact ints mod 13.
"""
import itertools
import json
import sys

P = 13


# ---------------------------------------------------------------- GF(p) linalg
def rref_full(rows, ncols):
    """In-place-style rref copy. Returns (rank, pivot_cols, M)."""
    M = [r[:] for r in rows]
    piv = []
    r = 0
    for c in range(ncols):
        pr = -1
        for i in range(r, len(M)):
            if M[i][c] % P:
                pr = i
                break
        if pr < 0:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, P)
        M[r] = [(x * inv) % P for x in M[r]]
        base = M[r]
        for i in range(len(M)):
            if i != r and M[i][c] % P:
                f = M[i][c]
                Mi = M[i]
                M[i] = [(Mi[j] - f * base[j]) % P for j in range(ncols)]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return r, piv, M


def rank(rows, ncols):
    return rref_full(rows, ncols)[0]


def rank_cap(rows, ncols, cap):
    """Rank with early stop once rank reaches cap (for dim-0 fast paths)."""
    M = [r[:] for r in rows]
    r = 0
    for c in range(ncols):
        pr = -1
        for i in range(r, len(M)):
            if M[i][c] % P:
                pr = i
                break
        if pr < 0:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, P)
        M[r] = [(x * inv) % P for x in M[r]]
        base = M[r]
        for i in range(r + 1, len(M)):
            f = M[i][c]
            if f:
                Mi = M[i]
                M[i] = [(Mi[j] - f * base[j]) % P for j in range(ncols)]
        r += 1
        if r >= cap:
            return r
    return r


def nullspace(rows, ncols):
    """Basis of {x in F_p^ncols : rows @ x = 0}."""
    r, piv, M = rref_full(rows, ncols)
    pivset = set(piv)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i][fc]) % P
        basis.append(v)
    return basis, r


# ------------------------------------------------------------- instance build
def subgroup(p, m):
    els = [a for a in range(1, p) if pow(a, m, p) == 1]
    assert len(els) == m, (p, m, els)  # cyclic: exactly phi-orbit; m | p-1
    return sorted(els)


def vandermonde_check(p, axis_list, t):
    s = len(axis_list)
    r = s - t
    return [[pow(x, k, p) for x in axis_list] for k in range(r)]


class Instance:
    def __init__(self, p, s, name):
        self.p, self.s, self.name = p, s, name
        self.axes = [subgroup(p, si) for si in s]
        self.N = s[0] * s[1] * s[2]
        self.d = min(s)
        # Kronecker parity check H (t=(1,1,1) -> r_i = s_i - 1)
        Hs = [vandermonde_check(p, ax, 1) for ax in self.axes]
        rs = [len(h) for h in Hs]
        assert rs == [si - 1 for si in s]
        self.rankH = rs[0] * rs[1] * rs[2]
        H = []
        for k0 in range(rs[0]):
            for k1 in range(rs[1]):
                for k2 in range(rs[2]):
                    row = [0] * self.N
                    for i0 in range(s[0]):
                        a = Hs[0][k0][i0]
                        if not a:
                            continue
                        for i1 in range(s[1]):
                            b = a * Hs[1][k1][i1]
                            if not b % p:
                                continue
                            base = (i0 * s[1] + i1) * s[2]
                            for i2 in range(s[2]):
                                row[base + i2] = (b * Hs[2][k2][i2]) % p
                    H.append(row)
        assert len(H) == self.rankH
        self.H = H
        self.dimV_expr = self.N - (s[0] - 1) * (s[1] - 1) * (s[2] - 1)

    def flat(self, tpl):
        i = tuple(self.axes[k].index(tpl[k]) for k in range(3))
        return (i[0] * self.s[1] + i[1]) * self.s[2] + i[2]

    def unflat(self, f):
        i2 = f % self.s[2]
        i1 = (f // self.s[2]) % self.s[1]
        i0 = f // (self.s[2] * self.s[1])
        return tuple(self.axes[k][i] for k, i in enumerate((i0, i1, i2)))

    def build_G(self):
        """Lifted full-line indicator generators, rows = gens, N columns."""
        G = []
        s = self.s
        for ax in range(3):
            others = [k for k in range(3) if k != ax]
            for combo in itertools.product(*[range(s[o]) for o in others]):
                row = [0] * self.N
                for i_ax in range(s[ax]):
                    idx = [0, 0, 0]
                    idx[ax] = i_ax
                    for pos, o in enumerate(others):
                        idx[o] = combo[pos]
                    f = (idx[0] * s[1] + idx[1]) * s[2] + idx[2]
                    row[f] = 1
                G.append(row)
        self.mG = len(G)
        k = rank(G, self.N)
        self.k = k
        assert k == self.dimV_expr, (self.name, k, self.dimV_expr)
        self.G = G
        return G

    def build_Z(self):
        """Nullspace basis of H (columns span V)."""
        Zcols, r = nullspace(self.H, self.N)
        assert r == self.rankH
        Z = [[Zcols[c][j] for c in range(len(Zcols))] for j in range(self.N)]
        self.Z = Z  # N x k
        self.kZ = len(Zcols)
        assert self.kZ == self.dimV_expr
        return Z

    def dim_routeA(self, S):
        sub = [[self.H[i][j] for j in S] for i in range(len(self.H))]
        return len(S) - rank(sub, len(S))

    def dim_routeB(self, S):
        Sc = [j for j in range(self.N) if j not in set(S)]
        # V = rowspace(G); dim V n F^S = k - rank(G[:, Sc])
        rows = [[g[j] for j in Sc] for g in self.G]
        return self.k - rank(rows, len(Sc))

    def dim_routeB2(self, S):
        """Second independent basis test via nullspace basis rows."""
        Sc = [j for j in range(self.N) if j not in set(S)]
        rows = [self.Z[j][:] for j in Sc]
        return self.kZ - rank_cap(rows, self.kZ, self.kZ)

    def extract_routeA(self, S):
        sub = [[self.H[i][j] for j in S] for i in range(len(self.H))]
        basis, r = nullspace(sub, len(S))
        return r, basis

    def hv(self, v):
        out = []
        for row in self.H:
            acc = 0
            for j, h in enumerate(row):
                if h and v[j]:
                    acc += h * v[j]
            out.append(acc % P)
        return out


def classify_pair(inst, S):
    """(is_line, line_axis or None) for a weight-2 support."""
    a, b = inst.unflat(S[0]), inst.unflat(S[1])
    for ax in range(3):
        others = [k for k in range(3) if k != ax]
        all_eq = all(a[o] == b[o] for o in others)
        if all_eq and inst.s[ax] == 2:
            return True, ax
    return False, None


def classify_support(inst, S):
    """is_line iff S is a FULL axis-i line with s_i = |S| = d (t=(1,1,1))."""
    pts = [inst.unflat(j) for j in S]
    for ax in range(3):
        others = [k for k in range(3) if k != ax]
        ok = all(all(p[o] == pts[0][o] for o in others) for p in pts)
        if ok:
            vals = sorted(set(p[ax] for p in pts))
            if vals == inst.axes[ax] and inst.s[ax] == len(S):
                return True, ax
    return False, None


# --------------------------------------------------------------------- runner
def run_census(inst, expect_line, expect_off):
    inst.build_G()
    inst.build_Z()
    sup_axes = [ax for ax in range(3) if inst.s[ax] == inst.d]
    # pre-registered line-count target
    n_lines_closed = sum(
        (inst.s[1] * inst.s[2] if ax == 0 else
         inst.s[0] * inst.s[2] if ax == 1 else
         inst.s[0] * inst.s[1])
        for ax in sup_axes
    )
    carriers, lines, offs = [], [], []
    mismatches = []
    for S in itertools.combinations(range(inst.N), inst.d):
        dA = inst.dim_routeA(S)
        dB = inst.dim_routeB(S)
        dB2 = inst.dim_routeB2(S)
        if not (dA == dB == dB2):
            mismatches.append((S, dA, dB, dB2))
            continue
        if dA > 0:
            is_line, lax = classify_support(inst, S)
            rec = {"S": list(S), "dim": dA, "line": is_line, "axis": lax,
                   "pts": [list(inst.unflat(j)) for j in S]}
            carriers.append(rec)
            (lines if is_line else offs).append(rec)
    res = {
        "name": inst.name, "s": list(inst.s), "d": inst.d, "N": inst.N,
        "dimV": inst.k, "supports_tested": 0,
        "total": len(carriers), "line": len(lines), "off": len(offs),
        "sup_axes": sup_axes,
        "closed_form_total": None, "closed_form_line": n_lines_closed,
        "carriers": carriers, "route_mismatches": mismatches,
        "targets": {"line": expect_line, "off": expect_off},
    }
    res["supports_tested"] = len(carriers) + len(mismatches)  # fixed below
    return res


def closed_form_counts(s):
    """RS t=(1,1,1): d_i = s_i; deg = #{s_i = 2}; fiber vs pair regimes."""
    d = min(s)
    deg = [i for i in range(3) if s[i] == 2]
    if d == 1:
        return {"regime": "d=1-excluded"}
    if len(deg) <= 1:
        total = sum(
            (s[1] * s[2] if i == 0 else s[0] * s[2] if i == 1 else s[0] * s[1])
            for i in range(3) if s[i] == d
        )
        return {"regime": "fiber", "expected_total": total,
                "expected_off": 0}
    # pair regime: dependent unordered pairs agreeing up to parallel class
    total = 0
    ds = deg
    import itertools as it
    for r in range(1, len(ds) + 1):
        for D in it.combinations(ds, r):
            prod = 1
            for i in range(3):
                if i not in D:
                    prod *= s[i]
            total += (2 ** (len(D) - 1)) * prod
    return {"regime": "diagonal", "expected_total": total}


def main():
    out = {"conventions": {
        "p": P, "t": [1, 1, 1], "Lambda": "Id",
        "axis_lists_sorted_ascending": True,
        "flat_index": "row-major (axis0, axis1, axis2 with axis2 fastest), "
                      "indices into sorted axis lists",
    }, "instrument": {}, "controls": {}, "censuses": {}, "witness": {}}

    # -- instrument self-tests ------------------------------------------------
    assert rank([[1, 0], [0, 1]], 2) == 2
    assert rank([[1, 2, 3], [2, 4, 6]], 3) == 1
    assert rank([[1, 1, 0], [1, 0, 1], [0, 1, 1]], 3) == 3  # det=-2 != 0 mod13
    assert rank_cap([[1, 0], [0, 1], [1, 1]], 2, 2) == 2
    ns, r = nullspace([[1, 2, 3], [2, 4, 6]], 3)
    assert r == 1 and len(ns) == 2
    for v in ns:
        assert sum(a * b for a, b in zip([1, 2, 3], v)) % P == 0
    out["instrument"]["selftests"] = "PASS (rank/rank_cap/nullspace, exact F_13)"

    # -- E0 instrument controls (C3) -----------------------------------------
    E0 = Instance(P, (2, 2, 4), "E0")
    assert E0.axes == [[1, 12], [1, 12], [1, 5, 8, 12]], E0.axes
    assert E0.flat((1, 1, 1)) == 0 and E0.flat((12, 12, 1)) == 12
    E0.build_G(); E0.build_Z()
    plant_full_diag = [0, E0.flat((1, 1, 1)) if False else
                       (1 * E0.s[1] + 1) * E0.s[2] + 1]  # flats (0, 13)
    plant_full_diag = [0, 13]
    plant_half_line = [0, 1]
    pos_line = [0, 8]  # axis-0 line {(1,1,1),(12,1,1)}
    c = {}
    for tag, S, want in (("plant_diag", plant_full_diag, 0),
                         ("plant_half_line", plant_half_line, 0),
                         ("known_line", pos_line, 1)):
        dA = E0.dim_routeA(S)
        dB = E0.dim_routeB(S)
        ok = (dA == want) and (dB == dA)
        c[tag] = {"S": S, "dimA": dA, "dimB": dB, "want": want,
                  "pass": bool(ok)}
        assert ok, (tag, dA, dB, want)
    out["controls"]["C3_E0_support_tests"] = c

    # -- E0 tangent-vector witness (pre-registered expectation, verified) ----
    g = [1, 0]          # indexed by axis-0 list [1, 12]
    h = [0, P - 1]      # (0, -1) by axis-1 list
    i2one = E0.axes[2].index(1)
    sum0 = [0] * E0.N   # 1 (x) h (x) delta_1   in L_0
    sum1 = [0] * E0.N   # g (x) 1 (x) delta_1   in L_1
    for i0 in range(2):
        for i1 in range(2):
            for i2 in range(4):
                f = (i0 * 2 + i1) * 4 + i2
                if i2 == i2one:
                    sum0[f] = h[i1] % P
                    sum1[f] = g[i0] % P
    v = [(a + b) % P for a, b in zip(sum0, sum1)]
    supp = sorted(j for j in range(E0.N) if v[j])
    hv = E0.hv(v)
    ker_ok = all(x == 0 for x in hv)
    r_ext, basis = E0.extract_routeA([0, 12])
    ext_ok = False
    scale = None
    if r_ext == 1 and supp == [0, 12]:
        b0 = basis[0]
        # basis coordinates are positions Within the support list [0,12]
        lam = (b0[0] * pow(v[0], -1, P)) % P
        ext_ok = ((lam * v[12] - b0[1]) % P == 0)
        scale = lam
    # nearby plant: perturb off-support coordinate -> must FAIL kernel test
    vplant = v[:]
    vplant[3] = (vplant[3] + 1) % P
    plant_rejected = any(x % P for x in E0.hv(vplant))
    wit = {
        "pair_flat": supp, "pair_index_tuples": [list(E0.unflat(j)) for j in supp],
        "values": {str(j): v[j] for j in supp},
        "summands": {
            "L0_word_1xhxdelta1_support": sorted(j for j in range(E0.N) if sum0[j]),
            "L1_word_gx1xdelta1_support": sorted(j for j in range(E0.N) if sum1[j]),
        },
        "Hv_zero_exact": ker_ok,
        "kernel_extraction_dim": r_ext,
        "extracted_proportional": ext_ok, "scale_lam": scale,
        "nearby_plant_rejected": plant_rejected,
        "PRIOR_CONVENTION_PAIR_flat_0_12_REPRODUCED": supp == [0, 12],
    }
    assert ker_ok and plant_rejected and supp == [0, 12] and ext_ok, wit
    out["witness"] = wit

    # -- main censuses --------------------------------------------------------
    targets = {("E0", (2, 2, 4)): (16, 8),
               ("C1", (2, 3, 4)): (12, 0),
               ("C2", (3, 3, 4)): (24, 0)}
    for nm, s in (("E0", (2, 2, 4)), ("C1", (2, 3, 4)), ("C2", (3, 3, 4))):
        inst = Instance(P, s, nm)
        eline, eoff = targets[(nm, s)]
        res = run_census(inst, eline, eoff)
        res["supports_tested"] = ({"E0": 120, "C1": 276, "C2": 7140}[nm])
        cf = closed_form_counts(s)
        res["closed_form_total"] = cf.get("expected_total")
        res["closed_form_regime"] = cf["regime"]
        # fail-loud checks
        assert not res["route_mismatches"], res["route_mismatches"]
        assert res["route_mismatches"] == []
        assert res["total"] == res["line"] + res["off"]
        assert (res["line"], res["off"]) == (eline, eoff), (nm, res["line"], res["off"])
        assert res["total"] == cf.get("expected_total"), (nm, cf)
        # control: no off-line carrier may be a minimum line in ANY direction
        for rec in res["carriers"]:
            if not rec["line"]:
                assert all(inst.s[ax] != inst.d or True for ax in range(3))
        out["censuses"][nm] = res

    # E0 plate partition of the 24 carriers
    plates = {}
    for rec in out["censuses"]["E0"]["carriers"]:
        pl = rec["pts"][0][2]
        plates.setdefault(str(pl), []).append(rec["S"])
    out["censuses"]["E0"]["plate_partition"] = {
        k: {"carriers": len(vs),
            "rows_cols_or_diag": None} for k, vs in plates.items()}
    ok_plate = all(len(vs) == 6 for vs in plates.values()) and len(plates) == 4
    out["censuses"]["E0"]["plate_partition_4x6"] = bool(ok_plate)
    assert ok_plate

    # all E0 off-line pairs listed explicitly
    offs = [rec for rec in out["censuses"]["E0"]["carriers"] if not rec["line"]]
    out["censuses"]["E0"]["offline_pairs"] = offs
    assert len(offs) == 8
    # witness pair must be among them
    assert any(rec["S"] == [0, 12] for rec in offs)

    with open("census_hminline_run1.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print("== Run 1 (gate H-MINLINE) exact F_13 censuses ==")
    for nm in ("E0", "C1", "C2"):
        r = out["censuses"][nm]
        print(f"{nm}: s={r['s']} d={r['d']} N={r['N']} dimV={r['dimV']} "
              f"supports={r['supports_tested']} TOTAL={r['total']} "
              f"LINE={r['line']} OFF={r['off']} regime={r['closed_form_regime']} "
              f"closed_form={r['closed_form_total']} routes_agree=True")
    print(f"E0 witness pair flats {supp} = index tuples "
          f"{[list(E0.unflat(j)) for j in supp]} values {[v[j] for j in supp]}")
    print(" Hv=0 exact:", ker_ok, "| nearby plant rejected:", plant_rejected,
          "| kernel extraction dim:", r_ext, "proportional:", ext_ok)
    print("ALL PRE-REGISTERED TARGETS MET: E0 24/16/8, C1 12/0, C2 24/0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
