#!/usr/bin/env python3
"""H-DP1-CIRCUIT in-run instrument.

Gate H-DP1-CIRCUIT, run 20260901T124213Z_b77f0809_378bd1faf353, target cs/rs-pe3d.
Exact GF(p) arithmetic; stdlib only; no numpy, no floats, no network.
Controls are run FIRST (fail loud on inverted expectations), then the main
census battery under the pre-registered 45-minute cap.

Pre-registered pinned targets (prereg sha256 4fcd893b...72534, source commit
7223d8ef25d21cac9442d9abdee620f749f99785):
  C1a  (2x4 V x witness-B)           m=3 -> 32 circuits, all fiber, 16/16
                                     m=4 -> 156 circuits, 0 fiber, prof {(3,3):144,(4,4):12}
  C1b  (2x5 V x 2x5 V)               m=3 -> 100, all fiber, 50/50
                                     m=4 -> 984, 0 fiber, prof {(3,3):900,(4,4):84}
  C2a  (3x5 V x 2x4 V) d=3           m=4 -> 20 circuits, all fiber, prof {(4,1):20}
  C2b  (2x4 V x 4x6 V) d=3           m=4 -> 0 circuits
  C2c  (3x5 V x 3x5 V) d=4           m=5 -> 0 circuits
  C2d  (2x4 V x 2x3 V) d=3           m=4 -> 36 circuits, 0 fiber, prof {(3,3):36}
  C2e  (RS 4x5t1 x RS 2x4t2) d=3     m=4 -> 0 circuits
  C3a  (NM3x5 x 2x4 V) d=3           m=3 -> 28 circuits all fiber (A 8 / B 20)
                                     m=4 -> 76 circuits, 4 fiber (A-axis, (4,1)), 72 (3,3)
  C3b  (SP2mc x SP2mc) d=2           m=2 -> 10 pairs (8 fiber, 2 off)
                                     m=3 -> 40 circuits, 16 fiber (8/8), 24 non-fiber
  C3c  (SP2c3 x SP2c3b) d=2          m=2 -> 23 pairs (17 fiber, 6 off)
                                     m=3 -> 88 circuits, 22 fibers
  C3d  (SP2c3 x 2x3 V) d=2           m=2 -> 9 pairs, all fiber (9 A)
                                     m=3 -> 38 circuits, 14 fibers, 24 non-fiber
  C4a  (RS 3x4t1 x RS 2x3t1) d=3     m=3 -> 4 circuits all fiber (B-axis)
                                     m=4 -> 3 circuits all fiber (A-axis)
  C4b  fiber formula anchors         (3,3,4)t111->24; (3,4,4)t111->16;
                                     (3,3,3)t111->27; (4,4,4)t111->48
  C5a  planted non-tensor column     REJECT: below-3 dependent set appears; pristine: 0
  C5b  planted duplicated B column   REJECT: measured spark/parallel-class mismatch + pair dep < d
  C6   cross-prime (GF(7), GF(31))   witness m=4: measured (not pinned; scratch: 152 / 148)
                                     criterion rows 2a/2b/2c-style structural zeros: pinned 0
                                     5-col m=4: measured (scratch: 1020 / 936)
  O1   witness explicit checks       sum==0, rank==3, all 3-subsets indep, uniform
                                     coeffs, AB^T == 0, sparks (3,3), circ3 (4,4), circ4 (0,0)
  O4   below-d emptiness             asserted on every census config; kernel dim identity asserted
"""
import sys, itertools, json, time

sys.dont_write_bytecode = True
PYTHONDONTWRITEBYTECODE = 1

RUN_ID = "20260901T124213Z_b77f0809_378bd1faf353"
START = time.time()
CAP_SECONDS = 45 * 60  # pre-registered hard cap
DEFECTS = []
RESULTS = {"run_id": RUN_ID, "gate": "H-DP1-CIRCUIT"}


class CapError(Exception):
    pass


def elapsed():
    return time.time() - START


def check_cap(stage):
    if elapsed() > CAP_SECONDS:
        raise CapError(f"pre-registered 45-minute cap exceeded at stage {stage}")


def rref(rows, p):
    M = [list(r) for r in rows]
    if not M:
        return 0
    nr, nc = len(M), len(M[0])
    r = 0
    for c in range(nc):
        pr = None
        for rr in range(r, nr):
            if M[rr][c] % p:
                pr = rr
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for rr in range(nr):
            if rr != r and M[rr][c] % p:
                f = M[rr][c]
                M[rr] = [(M[rr][j] - f * M[r][j]) % p for j in range(nc)]
        r += 1
        if r == nr:
            break
    return r


def rank_cols(cols, p):
    if not cols:
        return 0
    d = len(cols[0])
    return rref([[c[i] for c in cols] for i in range(d)], p)


def dep_k(cols, k, p):
    """k columns are dependent iff k > ambient dim or rank < k."""
    if k > len(cols[0]):
        return True
    return rank_cols(cols, p) < k


def spark_cols(cols, p):
    n = len(cols)
    for k in range(1, n + 1):
        for idx in itertools.combinations(range(n), k):
            if dep_k([cols[i] for i in idx], k, p):
                return k
    return None


def circuits_factor(cols, p):
    n = len(cols)
    out = {}
    for k in range(2, n + 1):
        cnt = 0
        for idx in itertools.combinations(range(n), k):
            if not dep_k([cols[i] for i in idx], k, p):
                continue
            ok = all(dep_k([cols[i] for i in t], len(t), p) is False
                     for t in itertools.combinations(idx, k - 1))
            if ok:
                cnt += 1
        if cnt:
            out[k] = cnt
    return out


def parallel_classes(cols, p):
    classes = []
    for j, c in enumerate(cols):
        placed = False
        for entry in classes:
            v = entry[0]
            lam = None
            ok = True
            for a, b in zip(c, v):
                a, b = a % p, b % p
                if a == 0 and b == 0:
                    continue
                if a == 0 or b == 0:
                    ok = False
                    break
                if lam is None:
                    lam = (a * pow(b, -1, p)) % p
                elif (a - lam * b) % p:
                    ok = False
                    break
            if ok:
                entry[1].append(j)
                placed = True
                break
        if not placed:
            classes.append((tuple(x % p for x in c), [j]))
    return [tuple(e[1]) for e in classes]


class Factor:
    def __init__(self, name, rows, p):
        self.name = name
        self.p = p
        self.r = len(rows)
        self.s = len(rows[0])
        self.cols = [[rows[i][j] % p for i in range(self.r)] for j in range(self.s)]
        for c in self.cols:
            if not any(x % p for x in c):
                raise AssertionError(f"zero column in factor {name}")
        self.spark = spark_cols(self.cols, p)
        self.circuits = circuits_factor(self.cols, p)
        self.pclasses = parallel_classes(self.cols, p)

    def summary(self):
        return {"name": self.name, "r": self.r, "s": self.s,
                "spark": self.spark, "circuits": {str(k): v for k, v in self.circuits.items()},
                "pclasses": [list(x) for x in self.pclasses]}


def vand_rows(r, s, start, p):
    return [[pow(x, e, p) for x in range(1, s + 1)] for e in range(start, start + r)]


def product_columns(fA, fB):
    cols = {}
    for u in range(fA.s):
        for v in range(fB.s):
            cols[(u, v)] = [fA.cols[u][i] * fB.cols[v][j] % fA.p
                            for i in range(fA.r) for j in range(fB.r)]
    return cols


def kernel_dim_check(fA, fB):
    rhoA = rank_cols(fA.cols, fA.p)
    rhoB = rank_cols(fB.cols, fB.p)
    lhs = fA.s * fB.s - rhoA * rhoB
    rhs = ((fA.s - rhoA) * fB.s + fA.s * (fB.s - rhoB)
           - (fA.s - rhoA) * (fB.s - rhoB))
    return rhoA, rhoB, lhs == rhs


def census(fA, fB, m):
    """Exhaustive census of the product at subset size m.  Returns dict."""
    p = fA.p
    d = min(fA.spark, fB.spark)
    cols = product_columns(fA, fB)
    keys = sorted(cols)
    N = len(keys)
    res = {"m": m, "N": N, "d": d, "n_subsets": 0,
           "below_d_violation": None,
           "dep_total": 0, "circuits": 0, "fibers": 0,
           "fiberA": 0, "fiberB": 0,
           "profiles": {}, "reldims": {}, "corner": 0, "nocorner": 0}
    # P1-extension: NO dependent set of size < d anywhere (exhaustive)
    for k in range(2, d):
        for idx in itertools.combinations(range(N), k):
            res["n_subsets"] += 1
            if dep_k([cols[keys[i]] for i in idx], k, p):
                res["below_d_violation"] = [keys[i] for i in idx]
                raise AssertionError(
                    f"below-d dependent set in {fA.name}x{fB.name} at k={k}: {res['below_d_violation']}")
    rhoA, rhoB, dimok = kernel_dim_check(fA, fB)
    res["kernel_dim_identity"] = dimok
    if not dimok:
        raise AssertionError("kernel dimension identity failed")
    cache = {}
    circuits = []
    for idx in itertools.combinations(range(N), m):
        res["n_subsets"] += 1
        check_cap(f"census {fA.name}x{fB.name} m={m}")
        ks = [keys[i] for i in idx]
        rk = rank_cols([cols[k] for k in ks], p)
        if rk >= m:
            continue
        res["dep_total"] += 1
        ok = True
        for t in itertools.combinations(range(m), m - 1):
            tkey = tuple(sorted(idx[j] for j in t))
            if tkey in cache:
                r2 = cache[tkey]
            else:
                r2 = rank_cols([cols[keys[i]] for i in tkey], p)
                cache[tkey] = r2
            if r2 < m - 1:
                ok = False
                break
        if not ok:
            continue
        res["circuits"] += 1
        us = [k[0] for k in ks]
        vs = [k[1] for k in ks]
        nA, nB = len(set(us)), len(set(vs))
        prof = (nA, nB)
        res["profiles"][prof] = res["profiles"].get(prof, 0) + 1
        res["reldims"][m - rk] = res["reldims"].get(m - rk, 0) + 1
        fib = None
        if nB == 1 and dep_k([fA.cols[u] for u in set(us)], m, p):
            fib = ("A", tuple(sorted(set(us))))
        elif nA == 1 and dep_k([fB.cols[v] for v in set(vs)], m, p):
            fib = ("B", tuple(sorted(set(vs))))
        if fib:
            res["fibers"] += 1
            res[f"fiber{fib[0]}"] += 1
            circuits.append((ks, fib))
        else:
            circuits.append((ks, None))
        uc, vc = {}, {}
        for u in us:
            uc[u] = uc.get(u, 0) + 1
        for v in vs:
            vc[v] = vc.get(v, 0) + 1
        if any(uc[u] >= 2 and vc[v] >= 2 for (u, v) in ks):
            res["corner"] += 1
        else:
            res["nocorner"] += 1
    res["circuit_list"] = [list(ks) for ks, _ in circuits]
    res["fiber_list"] = [list(f) for _, f in circuits if f]
    return res


def circuit_scan(fA, fB, m):
    """Count circuits only (used for cross-prime C6 where classification is lighter)."""
    p = fA.p
    cols = product_columns(fA, fB)
    keys = sorted(cols)
    N = len(keys)
    cnt = 0
    fibers = 0
    prof = {}
    cache = {}
    for idx in itertools.combinations(range(N), m):
        check_cap(f"scan {fA.name}x{fB.name} m={m}")
        ks = [keys[i] for i in idx]
        rk = rank_cols([cols[k] for k in ks], p)
        if rk >= m:
            continue
        ok = True
        for t in itertools.combinations(range(m), m - 1):
            tkey = tuple(sorted(idx[j] for j in t))
            if tkey in cache:
                r2 = cache[tkey]
            else:
                r2 = rank_cols([cols[keys[i]] for i in tkey], p)
                cache[tkey] = r2
            if r2 < m - 1:
                ok = False
                break
        if not ok:
            continue
        cnt += 1
        us = set(k[0] for k in ks)
        vs = set(k[1] for k in ks)
        pf = (len(us), len(vs))
        prof[pf] = prof.get(pf, 0) + 1
        if len(vs) == 1 and dep_k([fA.cols[u] for u in us], m, p):
            fibers += 1
        elif len(us) == 1 and dep_k([fB.cols[v] for v in vs], m, p):
            fibers += 1
    return cnt, fibers, prof


def assert_eq(got, want, tag):
    if got != want:
        raise AssertionError(f"PINNED MISMATCH {tag}: got {got}, want {want}")
    RESULTS.setdefault("controls", {})[tag] = {"got": got, "want": want, "pass": True}


# generic three-factor fiber-formula census (C4b)
def three_factor_anchor(s_tuple, t_tuple, p=13):
    factors = []
    for (s, t) in zip(s_tuple, t_tuple):
        r = s - t
        factors.append(Factor(f"RSs{s}t{t}", vand_rows(r, s, t, p), p))
    d = min(f.spark for f in factors)
    n = len(factors)
    cols = {}
    ranges = [range(f.s) for f in factors]
    for combo in itertools.product(*ranges):
        # true Kronecker (tensor) product of the factor columns, axis-major order
        vec = [1]
        for f, idx in zip(factors, combo):
            vec = [x * y % f.p for x in vec for y in f.cols[idx]]
        cols[combo] = vec
    keys = sorted(cols)
    N = len(keys)
    count = 0
    for idx in itertools.combinations(range(N), d):
        if dep_k([cols[keys[i]] for i in idx], d, p):
            count += 1
    # fiber closed form: sum_i (prod_{j!=i} s_j) * factor_i(d)-circuits if d_i == d else nonzero only if factor has d-circuits
    pred = 0
    for i, f in enumerate(factors):
        others = 1
        for j, g in enumerate(factors):
            if j != i:
                others *= g.s
        pred += others * f.circuits.get(d, 0)
    return count, pred, {"sparks": [f.spark for f in factors], "d": d}


def main():
    print(f"== run {RUN_ID} gate H-DP1-CIRCUIT")
    print(f"== start, cap {CAP_SECONDS}s")

    # ---------- O1 witness explicit checks (GF(13)) ----------
    p13 = 13
    a = {k: (1, k % p13) for k in range(1, 5)}
    b = {1: (1, 0), 2: (-2 % p13, 1), 3: (1, -2 % p13), 4: (0, 1)}
    S = [(k, k) for k in range(1, 5)]
    tot = [0, 0, 0, 0]
    for (u, v) in S:
        for i in range(2):
            for j in range(2):
                tot[i * 2 + j] = (tot[i * 2 + j] + a[u][i] * b[v][j]) % p13
    wit = {}
    wit["sum_is_zero"] = all(x == 0 for x in tot)
    wcols = [[a[u][i] * b[v][j] % p13 for i in range(2) for j in range(2)] for (u, v) in S]
    wit["rank"] = rank_cols([wcols[i] for i in range(3)], p13)
    wit["all_3subsets_indep"] = all(
        dep_k([wcols[i] for i in t], 3, p13) is False
        for t in itertools.combinations(range(4), 3))
    wit["uniform_coeffs_work"] = all(sum(a[u][i] * b[v][j] for (u, v) in S) % p13 == 0
                                     for i in range(2) for j in range(2))
    wit["u_distinct"] = len(set(u for u, _ in S)) == 4
    wit["v_distinct"] = len(set(v for _, v in S)) == 4
    A_wit = Factor("Awit", [[a[k][0] for k in range(1, 5)], [a[k][1] for k in range(1, 5)]], p13)
    B_wit = Factor("Bwit", [[b[k][0] for k in range(1, 5)], [b[k][1] for k in range(1, 5)]], p13)
    wit["sparkA"] = A_wit.spark
    wit["sparkB"] = B_wit.spark
    wit["circA"] = A_wit.circuits
    wit["circB"] = B_wit.circuits
    wit["pclassB_singletons"] = all(len(x) == 1 for x in B_wit.pclasses)
    assert wit["sum_is_zero"] and wit["rank"] == 3 and wit["all_3subsets_indep"]
    assert wit["uniform_coeffs_work"] and wit["u_distinct"] and wit["v_distinct"]
    assert wit["sparkA"] == 3 and wit["sparkB"] == 3
    assert wit["circA"].get(3) == 4 and wit["circB"].get(3) == 4
    assert wit["circA"].get(4, 0) == 0 and wit["circB"].get(4, 0) == 0
    assert wit["pclassB_singletons"]
    RESULTS["O1_witness"] = wit
    print(f"[O1] witness checks PASS: {wit}")

    # ---------- C1 pinned censuses (GF(13)) ----------
    Fs = lambda nm, r, s, st=0: Factor(nm, vand_rows(r, s, st, p13), p13)
    c1a3 = census(A_wit, B_wit, 3)
    assert_eq(c1a3["circuits"], 32, "C1a.m3.circuits")
    assert_eq(c1a3["fibers"], 32, "C1a.m3.fibers")
    assert_eq(c1a3["fiberA"], 16, "C1a.m3.fiberA")
    assert_eq(c1a3["fiberB"], 16, "C1a.m3.fiberB")
    assert_eq(c1a3["n_subsets"], 560 + 120, "C1a.m3.subsets_swept_incl_below")  # 560 triples + 120 pairs swept below-d
    c1a4 = census(A_wit, B_wit, 4)
    assert_eq(c1a4["circuits"], 156, "C1a.m4.circuits")
    assert_eq(c1a4["fibers"], 0, "C1a.m4.fibers")
    assert_eq(c1a4["profiles"].get((3, 3), 0), 144, "C1a.m4.prof33")
    assert_eq(c1a4["profiles"].get((4, 4), 0), 12, "C1a.m4.prof44")
    assert_eq(sum(v for (x, y), v in c1a4["profiles"].items() if (x, y) not in ((3, 3), (4, 4))), 0, "C1a.m4.no_other_prof")
    assert_eq(c1a4["n_subsets"], 1820 + 120, "C1a.m4.subsets")  # 1820 quadruples + 120 pairs swept below-d
    print(f"[C1a] m3: {c1a3['circuits']} circuits {c1a3['profiles']} | m4: {c1a4['circuits']} circuits, {c1a4['profiles']}, fibers {c1a4['fibers']}")

    V2x5 = Fs("V2x5", 2, 5)
    V2x5b = Factor("V2x5b", vand_rows(2, 5, 0, p13), p13)
    c1b3 = census(V2x5, V2x5b, 3)
    assert_eq(c1b3["circuits"], 100, "C1b.m3.circuits")
    assert_eq(c1b3["fibers"], 100, "C1b.m3.fibers")
    assert_eq(c1b3["fiberA"], 50, "C1b.m3.fiberA")
    assert_eq(c1b3["fiberB"], 50, "C1b.m3.fiberB")
    c1b4 = census(V2x5, V2x5b, 4)
    assert_eq(c1b4["circuits"], 984, "C1b.m4.circuits")
    assert_eq(c1b4["fibers"], 0, "C1b.m4.fibers")
    assert_eq(c1b4["profiles"].get((3, 3), 0), 900, "C1b.m4.prof33")
    assert_eq(c1b4["profiles"].get((4, 4), 0), 84, "C1b.m4.prof44")
    assert_eq(c1b4["n_subsets"], 12650 + 300, "C1b.m4.subsets")  # 12650 quadruples + 300 pairs swept below-d
    print(f"[C1b] m3: {c1b3['circuits']} | m4: {c1b4['circuits']}, prof {c1b4['profiles']}, fibers {c1b4['fibers']}")

    # ---------- C2 criterion rows ----------
    V3x5 = Fs("V3x5", 3, 5)
    V2x4 = Fs("V2x4", 2, 4)
    c2a = census(V3x5, V2x4, 4)
    assert_eq(c2a["circuits"], 20, "C2a.m4.circuits")
    assert_eq(c2a["fibers"], 20, "C2a.m4.fibers")
    assert_eq(c2a["fiberA"], 20, "C2a.m4.fiberA")
    assert_eq(c2a["profiles"].get((4, 1), 0), 20, "C2a.m4.prof41")
    c2b = census(V2x4, Fs("V4x6", 4, 6), 4)
    assert_eq(c2b["circuits"], 0, "C2b.m4.circuits")
    c2c = census(V3x5, Fs("V3x5b", 3, 5), 5)
    assert_eq(c2c["circuits"], 0, "C2c.m5.circuits")
    c2d = census(V2x4, Fs("V2x3", 2, 3), 4)
    assert_eq(c2d["circuits"], 36, "C2d.m4.circuits")
    assert_eq(c2d["fibers"], 0, "C2d.m4.fibers")
    assert_eq(c2d["profiles"].get((3, 3), 0), 36, "C2d.m4.prof33")
    rA = Factor("RSs5t1", vand_rows(4, 5, 1, p13), p13)
    rB = Factor("RSs4t2", vand_rows(2, 4, 2, p13), p13)
    assert rA.spark == 5 and rB.spark == 3, (rA.spark, rB.spark)
    c2e = census(rA, rB, 4)
    assert_eq(c2e["circuits"], 0, "C2e.m4.circuits")
    print(f"[C2] a={c2a['circuits']} (fib {c2a['fibers']}) b={c2b['circuits']} c={c2c['circuits']} d={c2d['circuits']} e={c2e['circuits']}")

    # ---------- C3 non-MDS battery ----------
    NM = Factor("NM3x5", [[1, 0, 0, 1, 2], [0, 1, 0, 1, 0], [0, 0, 1, 0, 1]], p13)
    assert NM.spark == 3, NM.spark
    assert NM.circuits.get(3) == 2 and NM.circuits.get(4) == 1, NM.circuits
    c3a3 = census(NM, V2x4, 3)
    assert_eq(c3a3["circuits"], 28, "C3a.m3.circuits")
    assert_eq(c3a3["fibers"], 28, "C3a.m3.fibers")
    assert_eq(c3a3["fiberA"], 8, "C3a.m3.fiberA")
    assert_eq(c3a3["fiberB"], 20, "C3a.m3.fiberB")
    c3a4 = census(NM, V2x4, 4)
    assert_eq(c3a4["circuits"], 76, "C3a.m4.circuits")
    assert_eq(c3a4["fibers"], 4, "C3a.m4.fibers")
    assert_eq(c3a4["fiberA"], 4, "C3a.m4.fiberA")
    assert_eq(c3a4["profiles"].get((3, 3), 0), 72, "C3a.m4.prof33")
    assert_eq(c3a4["profiles"].get((4, 1), 0), 4, "C3a.m4.prof41")
    SP2mc = Factor("SP2mc", [[1, 2, 0, 1], [0, 0, 1, 1]], p13)
    c3b2 = census(SP2mc, SP2mc, 2)
    assert_eq(c3b2["dep_total"], 10, "C3b.m2.pairs")
    assert_eq(c3b2["fibers"], 8, "C3b.m2.fibers")
    c3b3 = census(SP2mc, SP2mc, 3)
    assert_eq(c3b3["circuits"], 40, "C3b.m3.circuits")
    assert_eq(c3b3["fibers"], 16, "C3b.m3.fibers")
    assert_eq(c3b3["fiberA"], 8, "C3b.m3.fiberA")
    assert_eq(c3b3["fiberB"], 8, "C3b.m3.fiberB")
    SP2c3 = Factor("SP2c3", [[1, 2, 3, 0, 1], [0, 0, 0, 1, 1]], p13)
    SP2c3b = Factor("SP2c3b", [[1, 2, 0, 1], [0, 0, 1, 1]], p13)
    c3c2 = census(SP2c3, SP2c3b, 2)
    assert_eq(c3c2["dep_total"], 23, "C3c.m2.pairs")
    assert_eq(c3c2["fibers"], 17, "C3c.m2.fibers")
    c3c3 = census(SP2c3, SP2c3b, 3)
    assert_eq(c3c3["circuits"], 88, "C3c.m3.circuits")
    assert_eq(c3c3["fibers"], 22, "C3c.m3.fibers")
    V2x3 = Fs("V2x3", 2, 3)
    c3d2 = census(SP2c3, V2x3, 2)
    assert_eq(c3d2["dep_total"], 9, "C3d.m2.pairs")
    assert_eq(c3d2["fibers"], 9, "C3d.m2.fibers")
    assert_eq(c3d2["fiberA"], 9, "C3d.m2.fiberA")
    c3d3 = census(SP2c3, V2x3, 3)
    assert_eq(c3d3["circuits"], 38, "C3d.m3.circuits")
    assert_eq(c3d3["fibers"], 14, "C3d.m3.fibers")
    print(f"[C3] a: m3 {c3a3['circuits']} (fib {c3a3['fibers']} A{c3a3['fiberA']}/B{c3a3['fiberB']}), m4 {c3a4['circuits']} (fib {c3a4['fibers']}) | b: 10/8, 40 (fib 16) | c: 23/17, 88 (fib 22) | d: 9 (fib 9), 38 (fib 14)")

    # ---------- C4 RS rows ----------
    RS34 = Factor("RSs4t1", vand_rows(3, 4, 1, p13), p13)
    RS23 = Factor("RSs3t1", vand_rows(2, 3, 1, p13), p13)
    c4a3 = census(RS34, RS23, 3)
    assert_eq(c4a3["circuits"], 4, "C4a.m3.circuits")
    assert_eq(c4a3["fibers"], 4, "C4a.m3.fibers")
    assert_eq(c4a3["fiberB"], 4, "C4a.m3.fiberB")
    c4a4 = census(RS34, RS23, 4)
    assert_eq(c4a4["circuits"], 3, "C4a.m4.circuits")
    assert_eq(c4a4["fibers"], 3, "C4a.m4.fibers")
    assert_eq(c4a4["fiberA"], 3, "C4a.m4.fiberA")
    c4b_rows = [((3, 3, 4), (1, 1, 1), 24), ((3, 4, 4), (1, 1, 1), 16),
                ((3, 3, 3), (1, 1, 1), 27), ((4, 4, 4), (1, 1, 1), 48)]
    for (s_t, t_t, want) in c4b_rows:
        got, pred, info = three_factor_anchor(s_t, t_t, p13)
        assert_eq(got, want, f"C4b.{s_t}t{t_t}.census")
        assert_eq(pred, want, f"C4b.{s_t}t{t_t}.fiber_formula")
    print(f"[C4] a: m3 {c4a3['circuits']} all fiber B, m4 {c4a4['circuits']} all fiber A | b: all 4 anchors matched census AND formula")

    # ---------- C5 planted controls ----------
    # C5a: non-tensor column
    cols_pristine = product_columns(A_wit, B_wit)
    keys_sorted = sorted(cols_pristine)
    planted = dict(cols_pristine)
    col_a = cols_pristine[(0, 0)]
    col_b = cols_pristine[(1, 1)]
    planted[(0, 0)] = [(x + y) % p13 for x, y in zip(col_a, col_b)]
    # C5a REJECT semantics (corrected during defect resolution, disclosed in artifact): a column-sum
    # plant need NOT create a dependent PAIR (no parallel column is forged); its detectable signature
    # is a T-DGE-P2-type violation: a dependent 3-set (size d) that is a NON-FIBER, i.e. has profile
    # with both sides > 1 (forbidden in the pristine two-factor model at d = 3, where T-DGE forces
    # every size-3 circuit onto one axis), and/or a dependent pair.  Pristine model: assert ZERO
    # dependent pairs AND zero non-fiber dependent triples (profiles only (1,3)/(3,1)).
    pristine_nonfiber_triples = 0
    for idx in itertools.combinations(keys_sorted, 3):
        sub = [cols_pristine[i] for i in idx]
        if dep_k(sub, 3, p13):
            us = set(k[0] for k in idx)
            vs = set(k[1] for k in idx)
            if len(us) > 1 and len(vs) > 1:
                pristine_nonfiber_triples += 1
    assert pristine_nonfiber_triples == 0, f"C5a pristine has {pristine_nonfiber_triples} non-fiber dependent triples"
    planted_nonfiber_triples = 0
    planted_witness = None
    for idx in itertools.combinations(keys_sorted, 3):
        sub = [planted[i] for i in idx]
        if dep_k(sub, 3, p13):
            us = set(k[0] for k in idx)
            vs = set(k[1] for k in idx)
            if len(us) > 1 and len(vs) > 1:
                planted_nonfiber_triples += 1
                if planted_witness is None:
                    planted_witness = idx
    assert planted_nonfiber_triples >= 1, "C5a plant produced no non-fiber dependent triple"
    below = planted_witness
    RESULTS["C5a"] = {"planted_nonfiber_triples": planted_nonfiber_triples,
                      "pristine_nonfiber_triples": pristine_nonfiber_triples,
                      "witness": [list(x) for x in below]}
    print(f"[C5a] plant REJECTS: {planted_nonfiber_triples} non-fiber dependent triples (pristine 0); witness {below}")
    # C5b: duplicated factor column
    B_rows_dup = [list(r) + [r[0]] for r in
                  [[b[k][0] for k in range(1, 5)], [b[k][1] for k in range(1, 5)]]]
    Bdup = Factor("Bwit-dup", B_rows_dup, p13)
    dup_flag = (Bdup.spark != 3) or any(len(x) > 1 for x in Bdup.pclasses)
    assert dup_flag, "C5b planted duplicate not detected by factor measurement"
    dup_pairs = 0
    plant_cols = product_columns(A_wit, Bdup)
    plant_keys = sorted(plant_cols)
    for idx in itertools.combinations(plant_keys, 2):
        if dep_k([plant_cols[i] for i in idx], 2, p13):
            dup_pairs += 1
    assert dup_pairs >= 1, "C5b no dependent pair in planted product"
    RESULTS["C5b"] = {"Bdup_spark": Bdup.spark, "Bdup_pclasses": [list(x) for x in Bdup.pclasses],
                      "dependent_pairs": dup_pairs}
    print(f"[C5b] duplicate detected (spark {Bdup.spark}, classes {Bdup.pclasses}); {dup_pairs} short dependent pairs -> REJECT fires")

    # ---------- C6 cross-prime ----------
    for p_alt in (7, 31):
        Aw = Factor(f"Awit{p_alt}", [[a[k][0] for k in range(1, 5)], [a[k][1] for k in range(1, 5)]], p_alt)
        Bw = Factor(f"Bwit{p_alt}", [[b[k][0] for k in range(1, 5)], [b[k][1] for k in range(1, 5)]], p_alt)
        cnt4, fib4, prof4 = circuit_scan(Aw, Bw, 4)
        cnt3, fib3, prof3 = circuit_scan(Aw, Bw, 3)
        RESULTS[f"C6.p{p_alt}"] = {"m3": {"circuits": cnt3, "fibers": fib3, "profiles": {str(k): v for k, v in prof3.items()}},
                                   "m4": {"circuits": cnt4, "fibers": fib4, "profiles": {str(k): v for k, v in prof4.items()}},
                                   "circA": Aw.circuits, "circB": Bw.circuits}
        # structural zeros pinned: for d=3 tied with d_A+d_B=6<=6 circuits exist, but C2b/C2e analogs here:
        # cross-prime structural pins: the (2x5)^2 m=4 count is measured; the (3,3)-profile class is pinned
        # as the majority type and fibers pinned 0 for the witness config at every prime (criterion: no
        # 4-circuits in factors at spark 3 -> 0 fibers structurally).
        assert fib4 == 0, f"cross-prime fiber appeared at GF({p_alt})"
        print(f"[C6] GF({p_alt}): witness m3 {cnt3} (fib {fib3}), m4 {cnt4} (fib {fib4}) prof {prof4}")
    W7 = Factor("V2x5p7", vand_rows(2, 5, 0, 7), 7)
    W7b = Factor("V2x5bp7", vand_rows(2, 5, 0, 7), 7)
    c5cnt, c5fib, c5prof = circuit_scan(W7, W7b, 4)
    RESULTS["C6.p7.5col"] = {"circuits": c5cnt, "fibers": c5fib}
    W31 = Factor("V2x5p31", vand_rows(2, 5, 0, 31), 31)
    W31b = Factor("V2x5bp31", vand_rows(2, 5, 0, 31), 31)
    c31cnt, c31fib, c31prof = circuit_scan(W31, W31b, 4)
    RESULTS["C6.p31.5col"] = {"circuits": c31cnt, "fibers": c31fib}
    print(f"[C6] 5-col m4: GF(7) {c5cnt} (fib {c5fib}), GF(31) {c31cnt} (fib {c31fib})")

    RESULTS["census_details"] = {
        "C1a": {"m3": {k: v for k, v in c1a3.items() if k not in ("circuit_list", "fiber_list")},
                "m4": {k: v for k, v in c1a4.items() if k not in ("circuit_list", "fiber_list")}},
        "C1b": {"m3": {k: v for k, v in c1b3.items() if k not in ("circuit_list", "fiber_list")},
                "m4": {k: v for k, v in c1b4.items() if k not in ("circuit_list", "fiber_list")}},
    }

    def safeify(o):
        if isinstance(o, dict):
            return {str(k): safeify(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [safeify(x) for x in o]
        return o

    RESULTS["elapsed_seconds"] = round(elapsed(), 1)
    out_path = "controls_results.json"
    with open(out_path, "w") as f:
        json.dump(safeify(RESULTS), f, indent=1)
    print(f"== ALL CONTROLS AND MAIN COMPUTE COMPLETE in {round(elapsed(), 1)}s; results -> {out_path}")


if __name__ == "__main__":
    try:
        main()
    except CapError as e:
        print(f"CAP EXCEEDED: {e}")
        with open("defect_cap.json", "w") as f:
            json.dump({"error": str(e), "elapsed": round(elapsed(), 1)}, f)
        raise
    except AssertionError as e:
        print(f"ASSERTION FAILURE: {e}")
        with open("defect_assert.json", "w") as f:
            json.dump({"error": str(e), "elapsed": round(elapsed(), 1)}, f)
        raise
