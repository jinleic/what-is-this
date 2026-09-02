#!/usr/bin/env python3
"""H-MIX-CROSSING in-run instrument.

Gate: H-MIX-CROSSING (run 5, target cs/rs-pe3d, agent RsPe3dCrossing).
Prereg: prereg/H_CROSSING_PREREG_2026-09-01.md, commit cd2939676258c8be164f227cf7a391e1e1a4a9b4,
sha256 8c24ac710b7559f02d02acc1f57bdf022951c25ce14e58794f09b1eddd993ef5.

Verifies Theorem H-MIX-CROSSING (crossing/mixed-profile channel of A (x) B at
profile (d_A,d_B), size d_A+d_B-2): construction, circuitness, bijection with
closed form d_A*d_B*C_A(d_A)*C_B(d_B), field independence, reach iff
d_A+d_B = d+3, converse, plus controls T1-T9 per prereg section 3.

Exact GF(p) arithmetic; stdlib only; no numpy/floats. Every claim is checked
in-run AFTER campaign init. Defect policy: broken outputs preserved; any fix
is disclosed in defect_*.json and in theorem_crossing.md.
"""
import sys, itertools, json, time

sys.dont_write_bytecode = True
PYTHONDONTWRITEBYTECODE = 1

RUN_ID = "20260901T132343Z_4f169cbb_2f81dd8f2d8d"
START = time.time()
CAP_SECONDS = 60 * 60  # pre-registered hard cap (prereg section 4)
DEFECTS = []
RESULTS = {"run_id": RUN_ID, "gate": "H-MIX-CROSSING"}


class CapError(Exception):
    pass


def elapsed():
    return time.time() - START


def check_cap(stage):
    if elapsed() > CAP_SECONDS:
        raise CapError(f"pre-registered 60-minute cap exceeded at stage {stage}")


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


def circuit_relations_factor(f, k):
    n = f.s
    cols = f.cols
    out = []
    for idx in itertools.combinations(range(n), k):
        if dep_k([cols[i] for i in idx], k, f.p) is False:
            continue
        if any(dep_k([cols[i] for i in t], len(t), f.p) for t in itertools.combinations(idx, k - 1)):
            continue
        # circuit: relation space on idx is 1-dimensional; solve
        # sum_{i in idx} c_i cols[i] = 0, c_{idx[0]} = 1
        i0 = idx[0]
        M = [[cols[i][r_] % f.p for i in idx[1:]] for r_ in range(f.r)]
        # augment with -col_{i0}
        rhs = [(-cols[i0][r_]) % f.p for r_ in range(f.r)]
        # solve M c' = rhs over GF(p) via Gauss-Jordan on [M | rhs]
        A = [list(M[r_]) + [rhs[r_]] for r_ in range(f.r)]
        nr, nc = f.r, len(idx)  # nc = (k-1) coeffs + 1 rhs
        r = 0
        pivcols = []
        for c in range(len(idx) - 1):
            pr = None
            for rr in range(r, nr):
                if A[rr][c] % f.p:
                    pr = rr
                    break
            if pr is None:
                continue
            A[r], A[pr] = A[pr], A[r]
            inv = pow(A[r][c], -1, f.p)
            A[r] = [(x * inv) % f.p for x in A[r]]
            for rr in range(nr):
                if rr != r and A[rr][c] % f.p:
                    fct = A[rr][c]
                    A[rr] = [(A[rr][j] - fct * A[r][j]) % f.p for j in range(nc)]
            pivcols.append(c)
            r += 1
            if r == nr:
                break
        # consistency: no row 0=nonzero
        for rr in range(r, nr):
            if all(A[rr][c] % f.p == 0 for c in range(len(idx) - 1)) and A[rr][len(idx) - 1] % f.p:
                raise AssertionError("circuit relation solve inconsistent")
        sol = {}
        for c, i in enumerate(idx[1:]):
            v = 0
            for ri, pc in enumerate(pivcols):
                if pc == c:
                    v = A[ri][len(idx) - 1]
                    break
            sol[i] = v % f.p
        full = [0] * n
        full[i0] = 1
        for i in idx[1:]:
            full[i] = sol[i]
        assert all(full[i] % f.p for i in idx), "circuit relation not full support"
        out.append((idx, full))
    return out


def crossing_relations(fA, fB):
    """Theorem X construction: for each d_A-circuit R of A, i0 in R, d_B-circuit J of B,
    j0 in J: C = c e_{j0}^T + e_{i0} delta^T with delta scaled so c_{i0} + delta_{j0} = 0.
    Returns dict param -> (support cells, nonzero cell coefficients)."""
    p = fA.p
    dA, dB = fA.spark, fB.spark
    cirelsA = circuit_relations_factor(fA, dA)  # list of (idx, full vector)
    cirelsB = circuit_relations_factor(fB, dB)
    out = {}
    for (R, cvec) in cirelsA:
        for (J, dvec) in cirelsB:
            for i0 in R:
                for j0 in J:
                    # scale delta (full-support circuit relation of B on J) so that c_{i0} + s*d_{j0} = 0
                    s = (-cvec[i0] * pow(dvec[j0], -1, p)) % p
                    key = (R, i0, J, j0)
                    cells = {}
                    for i in R:
                        if i == i0:
                            continue
                        cells[(i, j0)] = cvec[i] % p
                    for j in J:
                        if j == j0:
                            continue
                        cells[(i0, j)] = (s * dvec[j]) % p
                    cells[(i0, j0)] = (cvec[i0] + s * dvec[j0]) % p  # == 0
                    assert cells[(i0, j0)] == 0
                    del cells[(i0, j0)]
                    out[key] = cells
    return out


def set_independence_sanity(fA, fB, ks, coeffs):
    """Check sum_i coeff_i * h(k_i) == 0 (exact)."""
    p = fA.p
    cols = product_columns(fA, fB)
    acc = [0] * (fA.r * fB.r)
    for k, c in zip(ks, coeffs):
        if c == 0:
            continue
        col = cols[k]
        for t in range(len(acc)):
            acc[t] = (acc[t] + c * col[t]) % p
    return all(x == 0 for x in acc)


def profile_of(ks):
    us = set(k[0] for k in ks)
    vs = set(k[1] for k in ks)
    return (len(us), len(vs))


def census(fA, fB, m):
    """Exhaustive census of the product at subset size m. Returns dict."""
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
    circuit_sets = []
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
        prof = profile_of(ks)
        res["profiles"][prof] = res["profiles"].get(prof, 0) + 1
        res["reldims"][m - rk] = res["reldims"].get(m - rk, 0) + 1
        fib = None
        if prof[1] == 1 and dep_k([fA.cols[u] for u in set(k[0] for k in ks)], m, p):
            fib = ("A", tuple(sorted(set(k[0] for k in ks))))
        elif prof[0] == 1 and dep_k([fB.cols[v] for v in set(k[1] for k in ks)], m, p):
            fib = ("B", tuple(sorted(set(k[1] for k in ks))))
        if fib:
            res["fibers"] += 1
            res[f"fiber{fib[0]}"] += 1
            circuit_sets.append((ks, fib))
        else:
            circuit_sets.append((ks, None))
    res["circuit_sets"] = circuit_sets
    return res


def assert_eq(got, want, tag):
    if got != want:
        raise AssertionError(f"PINNED MISMATCH {tag}: got {got}, want {want}")
    RESULTS.setdefault("controls", {})[tag] = {"got": got, "want": want, "pass": True}


# ---------------- Theorem X machinery ----------------

def mixed_prediction(fA, fB):
    dA, dB = fA.spark, fB.spark
    return dA * dB * fA.circuits.get(dA, 0) * fB.circuits.get(dB, 0)


def check_set_equality(fA, fB, cen, cls, tag):
    """The census's profile-(dA,dB) circuit sets must EXACTLY equal the constructed
    supports S(R,i0,J,j0). Also checks relation-space dim 1 and full-support relation."""
    p = fA.p
    dA, dB = fA.spark, fB.spark
    cols = product_columns(fA, fB)
    crel = crossing_relations(fA, fB)
    constructed = set()
    relinfo = {}
    for key, cells in crel.items():
        sup = tuple(sorted(cells.keys()))
        constructed.add(sup)
        relinfo[sup] = cells
    measured = set()
    m = cen["m"]
    for ks, fib in cen["circuit_sets"]:
        if profile_of(ks) == (dA, dB):
            measured.add(tuple(sorted(ks)))
            assert fib is None, f"{tag}: profile-(dA,dB) circuit flagged as fiber {fib}"
    miss = measured - constructed
    extra = constructed - measured
    assert_eq(len(measured), len(constructed), f"{tag}.cardinality")
    assert not miss, f"{tag}: measured-but-not-constructed {sorted(miss)[:5]}"
    assert not extra, f"{tag}: constructed-but-not-measured {sorted(extra)[:5]}"
    # relation-space dim 1 and explicit relation check on EVERY constructed support
    for sup, cells in list(relinfo.items()):
        ks = list(sup)
        coeffs = [cells[k] for k in ks]
        assert set_independence_sanity(fA, fB, ks, coeffs), f"{tag}: relation not in kernel for {sup}"
        rk = rank_cols([cols[k] for k in ks], p)
        assert rk == m - 1, f"{tag}: rank {rk} != m-1 for {sup}"
    RESULTS.setdefault("set_equality", {})[tag] = {
        "measured": len(measured), "constructed": len(constructed),
        "reldim_checked": len(relinfo)}
    return len(measured)


def main():
    print(f"== run {RUN_ID} gate H-MIX-CROSSING")
    print(f"== start, cap {CAP_SECONDS}s; prereg commit cd2939676258c8be164f227cf7a391e1e1a4a9b4")
    p13 = 13

    # ---------- T1 (ACCEPT block): witness + O1 explicit instance ----------
    a_rows = [[1, 1, 1, 1], [1, 2, 3, 4]]          # a_k = (1,k), k=1..4
    b_rows = [[1, -2 % 13, 1, 0], [0, 1, -2 % 13, 1]]  # cols (1,0),(-2,1),(1,-2),(0,1)
    Aw = Factor("Awit", a_rows, p13)
    Bw = Factor("Bwit", b_rows, p13)
    assert Aw.spark == 3 and Bw.spark == 3, (Aw.spark, Bw.spark)
    assert Aw.circuits == {3: 4} and Bw.circuits == {3: 4}, (Aw.circuits, Bw.circuits)
    assert Aw.circuits.get(4, 0) == 0 and Bw.circuits.get(4, 0) == 0
    assert all(len(x) == 1 for x in Bw.pclasses), Bw.pclasses
    _, _, idok = kernel_dim_check(Aw, Bw)
    assert idok
    RESULTS["T1_witness"] = {"sparkA": Aw.spark, "sparkB": Bw.spark,
                             "circA": Aw.circuits, "circB": Bw.circuits,
                             "pclassB_singletons": True, "kernel_dim_identity": True}
    print(f"[T1] witness sparks (3,3), circuits {{3:4}}/{{3:4}}, pclasses singleton, kernel-dim identity OK")
    # O1 explicit crossing instance: R={0,1,2}, i0=0; J={0,1,2}, j0=1
    cirelsA = circuit_relations_factor(Aw, 3)
    cirelsB = circuit_relations_factor(Bw, 3)
    R0 = cirelsA[0][0]  # {0,1,2}
    J0 = cirelsB[0][0]
    i0, j0 = 0, 1
    cvec = dict(zip(R0, [cirelsA[0][1][i] for i in R0]))
    dvec = dict(zip(J0, [cirelsB[0][1][j] for j in J0]))
    s_scale = (-cvec[i0] * pow(dvec[j0], -1, p13)) % p13
    cells = {}
    for i in R0:
        if i != i0:
            cells[(i, j0)] = cvec[i]
    for j in J0:
        if j != j0:
            cells[(i0, j)] = (s_scale * dvec[j]) % p13
    # center cancels: c_{i0} + s*d_{j0} == 0
    assert (cvec[i0] + s_scale * dvec[j0]) % p13 == 0
    o1 = {}
    o1["center_cancels"] = True
    o1["support_size"] = len(cells)
    o1["profile"] = profile_of(sorted(cells.keys()))
    o1["in_kernel"] = set_independence_sanity(Aw, Bw, sorted(cells.keys()),
                                              [cells[k] for k in sorted(cells.keys())])
    o1["support_all_cells_nonzero"] = all(v % p13 for v in cells.values())
    colsP = product_columns(Aw, Bw)
    o1["rank"] = rank_cols([colsP[k] for k in sorted(cells.keys())], p13)
    sup = sorted(cells.keys())
    o1["all_3subsets_indep"] = all(
        dep_k([colsP[k] for k in t], 3, p13) is False
        for t in itertools.combinations(sup, 3))
    assert o1["support_size"] == 4 and o1["profile"] == (3, 3)
    assert o1["in_kernel"] and o1["rank"] == 3 and o1["all_3subsets_indep"]
    assert o1["support_all_cells_nonzero"]
    RESULTS["O1_instance"] = o1
    print(f"[O1] explicit instance: support {sup}, rank {o1['rank']}, kernel-y, 3-subsets indep, profile (3,3)")

    # ---------- T2: witness census GF(13) m in {3,4} with set equality ----------
    c_m3 = census(Aw, Bw, 3)
    assert_eq(c_m3["circuits"], 32, "T2.m3.circuits")
    assert_eq(c_m3["fibers"], 32, "T2.m3.fibers")
    assert_eq(c_m3["fiberA"], 16, "T2.m3.fiberA")
    assert_eq(c_m3["fiberB"], 16, "T2.m3.fiberB")
    assert_eq(c_m3["n_subsets"], 560 + 120, "T2.m3.subsets_swept")
    print(f"[T2] m3: {c_m3['circuits']} circuits all fibers 16/16")
    c_m4 = census(Aw, Bw, 4)
    assert_eq(c_m4["circuits"], 156, "T2.m4.circuits")
    assert_eq(c_m4["fibers"], 0, "T2.m4.fibers")
    assert_eq(c_m4["profiles"].get((3, 3), 0), 144, "T2.m4.prof33")
    assert_eq(c_m4["profiles"].get((4, 4), 0), 12, "T2.m4.prof44")
    assert_eq(sum(v for (x, y), v in c_m4["profiles"].items() if (x, y) not in ((3, 3), (4, 4))),
              0, "T2.m4.no_other_prof")
    assert_eq(c_m4["n_subsets"], 1820 + 120, "T2.m4.subsets_swept")
    assert mixed_prediction(Aw, Bw) == 144
    print(f"[T2] m4: {c_m4['circuits']} circuits, fibers 0, prof {c_m4['profiles']}")
    # R,J index conventions: census keys are (u,v) 0-based; circuit_relations uses 0-based too
    n_mix = check_set_equality(Aw, Bw, c_m4, 144, "T2.witness")
    print(f"[T2] set equality witness m4: {n_mix} mixed <-> constructed, reldim 1 all")

    # ---------- T3a: cancellation sensitivity (perturb any one coefficient) ----------
    # For EVERY constructed crossing relation: pristine rank m*-1, reldim 1; perturbing
    # any single coefficient destroys the dependence (rank -> m*).
    def perturb_check(fA, fB, tag):
        """T3a: for EVERY constructed crossing relation: pristine rank m*-1;
        perturbing ANY single coefficient by ANY nonzero delta destroys the
        dependence (perturbed coefficients not a kernel vector)."""
        p = fA.p
        cr = crossing_relations(fA, fB)
        cols = product_columns(fA, fB)
        checked = 0
        pert_count = 0
        for key, cells in cr.items():
            sup = sorted(cells.keys())
            checked += 1
            ks = sup
            base = [cells[k] for k in ks]
            rk0 = rank_cols([cols[k] for k in ks], p)
            assert rk0 == len(sup) - 1, f"{tag}: pristine rank {rk0} != m*-1 at {key}"
            for pos in range(len(sup)):
                for delta in range(1, p):
                    pert = list(base)
                    pert[pos] = (pert[pos] + delta) % p
                    if pert[pos] == 0:
                        continue  # relation must keep full support
                    pert_count += 1
                    if set_independence_sanity(fA, fB, ks, pert):
                        raise AssertionError(
                            f"{tag}: perturbed relation stayed in-kernel at {key} pos {pos} delta {delta}")
        assert checked > 0, f"{tag}: no relations constructed"
        RESULTS.setdefault("T3a", {})[tag] = {"relations_checked": checked,
                                              "perturbations_swept": pert_count,
                                              "all_destroy_dependence": True}
        print(f"[T3a] {tag}: {checked} relations x {pert_count} perturbations, all destroy dependence")

    perturb_check(Aw, Bw, "witness.full")
    # ---------- T3b: spectrum-corruption quantitative row: 2x3 V x witness B ----------
    V2x3 = Factor("V2x3", vand_rows(2, 3, 0, p13), p13)
    assert V2x3.spark == 3 and V2x3.circuits == {3: 1}, (V2x3.spark, V2x3.circuits)
    c_tb = census(V2x3, Bw, 4)
    # prediction: mixed = d_A*d_B*C_A(3)*C_B(3) = 3*3*1*4 = 36; no all-distinct (s_A''=3 < 4);
    # no fibers (no 4-circuits on either factor: sparks (3,3), 4-circuits 0)
    assert V2x3.circuits.get(4, 0) == 0 and Bw.circuits.get(4, 0) == 0
    assert_eq(c_tb["circuits"], 36, "T3b.m4.total")
    assert_eq(c_tb["fibers"], 0, "T3b.m4.fibers")
    assert_eq(c_tb["profiles"].get((3, 3), 0), 36, "T3b.m4.prof33")
    assert_eq(sum(v for (x, y), v in c_tb["profiles"].items() if (x, y) != (3, 3)), 0,
              "T3b.m4.no_other_prof")
    check_set_equality(V2x3, Bw, c_tb, 36, "T3b.2x3V")

    # ---------- T3c: duplicated factor column REJECT ----------
    b_rows_dup = [list(r) + [r[0]] for r in b_rows]
    Bdup = Factor("Bwit-dup", b_rows_dup, p13)
    dup_flag = (Bdup.spark != 3) or any(len(x) > 1 for x in Bdup.pclasses)
    assert dup_flag, "T3c planted duplicate not detected by factor measurement"
    dup_pairs = 0
    plant_cols = product_columns(Aw, Bdup)
    plant_keys = sorted(plant_cols)
    for idx in itertools.combinations(plant_keys, 2):
        if dep_k([plant_cols[i] for i in idx], 2, p13):
            dup_pairs += 1
    assert dup_pairs >= 1, "T3c no dependent pair in planted product"
    RESULTS["T3c"] = {"Bdup_spark": Bdup.spark, "Bdup_pclasses_len": [len(x) for x in Bdup.pclasses],
                      "dependent_pairs": dup_pairs, "reject_fired": True}
    print(f"[T3c] duplicate REJECT: spark {Bdup.spark}, {dup_pairs} short dependent pairs")

    # ---------- T3d: non-tensor column REJECT ----------
    cols_pristine = product_columns(Aw, Bw)
    keys_sorted = sorted(cols_pristine)
    planted = dict(cols_pristine)
    col_a = cols_pristine[(0, 0)]
    col_b = cols_pristine[(1, 1)]
    planted[(0, 0)] = [(x + y) % p13 for x, y in zip(col_a, col_b)]
    def nonfiber_triples(cols):
        cnt = 0
        w = None
        for idx in itertools.combinations(keys_sorted, 3):
            sub = [cols[i] for i in idx]
            if dep_k(sub, 3, p13):
                us = set(k[0] for k in idx)
                vs = set(k[1] for k in idx)
                if len(us) > 1 and len(vs) > 1:
                    cnt += 1
                    if w is None:
                        w = idx
        return cnt, w
    p_cnt, _ = nonfiber_triples(cols_pristine)
    pl_cnt, pl_wit = nonfiber_triples(planted)
    assert p_cnt == 0, f"T3d pristine has {p_cnt} non-fiber dependent triples"
    assert pl_cnt >= 1, "T3d plant produced no non-fiber dependent triple"
    RESULTS["T3d"] = {"pristine_nonfiber_triples": p_cnt, "planted_nonfiber_triples": pl_cnt,
                      "witness": [list(x) for x in pl_wit]}
    print(f"[T3d] non-tensor REJECT: pristine 0, planted {pl_cnt}")

    # ---------- T4: 2x5 pair GF(13) ----------
    V2x5 = Factor("V2x5", vand_rows(2, 5, 0, p13), p13)
    V2x5b = Factor("V2x5b", vand_rows(2, 5, 0, p13), p13)
    assert V2x5.spark == 3 and V2x5b.spark == 3
    assert V2x5.circuits == {3: 10} and V2x5b.circuits == {3: 10}, (V2x5.circuits, V2x5b.circuits)
    c_2x5_4 = census(V2x5, V2x5b, 4)
    assert_eq(c_2x5_4["circuits"], 984, "T4.m4.circuits")
    assert_eq(c_2x5_4["fibers"], 0, "T4.m4.fibers")
    assert_eq(c_2x5_4["profiles"].get((3, 3), 0), 900, "T4.m4.prof33")
    assert_eq(c_2x5_4["profiles"].get((4, 4), 0), 84, "T4.m4.prof44_anchor")
    assert_eq(sum(v for (x, y), v in c_2x5_4["profiles"].items() if (x, y) not in ((3, 3), (4, 4))),
              0, "T4.m4.no_other_prof")
    assert_eq(c_2x5_4["n_subsets"], 12650 + 300, "T4.m4.subsets_swept")
    assert mixed_prediction(V2x5, V2x5b) == 900
    check_set_equality(V2x5, V2x5b, c_2x5_4, 900, "T4.2x5")
    c_2x5_3 = census(V2x5, V2x5b, 3)
    assert_eq(c_2x5_3["circuits"], 100, "T4.m3.circuits")
    assert_eq(c_2x5_3["fibers"], 100, "T4.m3.fibers")
    assert_eq(c_2x5_3["fiberA"], 50, "T4.m3.fiberA")
    assert_eq(c_2x5_3["fiberB"], 50, "T4.m3.fiberB")
    print(f"[T4] 2x5: m3 {c_2x5_3['circuits']} (fib 50/50), m4 {c_2x5_4['circuits']}, prof {c_2x5_4['profiles']}")

    # ---------- T5: 2x6 pair GF(13), full census ----------
    V2x6 = Factor("V2x6", vand_rows(2, 6, 0, p13), p13)
    V2x6b = Factor("V2x6b", vand_rows(2, 6, 0, p13), p13)
    assert V2x6.spark == 3 and V2x6.circuits == {3: 20}, (V2x6.spark, V2x6.circuits)
    c_2x6 = census(V2x6, V2x6b, 4)
    assert_eq(c_2x6["profiles"].get((3, 3), 0), 3600, "T5.m4.prof33")
    assert_eq(c_2x6["fibers"], 0, "T5.m4.fibers")
    assert_eq(sum(v for (x, y), v in c_2x6["profiles"].items() if (x, y) not in ((3, 3), (4, 4))),
              0, "T5.m4.no_other_prof")
    assert_eq(c_2x6["n_subsets"], 58905 + 630, "T5.m4.subsets_swept")
    assert mixed_prediction(V2x6, V2x6b) == 3600
    n26 = check_set_equality(V2x6, V2x6b, c_2x6, 3600, "T5.2x6")
    RESULTS["T5_total_and_prof44"] = {"total": c_2x6["circuits"],
                                      "prof44": c_2x6["profiles"].get((4, 4), 0)}
    print(f"[T5] 2x6: total {c_2x6['circuits']}, prof {c_2x6['profiles']}, set equality {n26}")

    # ---------- T6: 3x4 MDS pair GF(13): m=5 -> 0; m=6 -> 16 all (4,4) ----------
    M34 = Factor("M34", vand_rows(3, 4, 1, p13), p13)
    M34b = Factor("M34b", vand_rows(3, 4, 1, p13), p13)
    assert M34.spark == 4 and M34b.spark == 4, (M34.spark, M34b.spark)
    assert M34.circuits == {4: 1} and M34b.circuits == {4: 1}, (M34.circuits, M34b.circuits)
    assert M34.circuits.get(5, 0) == 0 and M34b.circuits.get(5, 0) == 0
    c_m5 = census(M34, M34b, 5)
    assert_eq(c_m5["circuits"], 0, "T6.m5.circuits")
    assert_eq(c_m5["dep_total"], 96, "T6.m5.dep_total (8 fiber-4-circuits x 12 other cols)")
    assert_eq(c_m5["n_subsets"], 4368 + 680, "T6.m5.subsets (C(16,5) + below-d 680)")
    c_m6 = census(M34, M34b, 6)
    assert_eq(c_m6["circuits"], 16, "T6.m6.circuits")
    assert_eq(c_m6["fibers"], 0, "T6.m6.fibers")
    assert_eq(c_m6["profiles"].get((4, 4), 0), 16, "T6.m6.prof44")
    assert_eq(sum(v for (x, y), v in c_m6["profiles"].items() if (x, y) != (4, 4)), 0,
              "T6.m6.no_other_prof")
    assert_eq(c_m6["n_subsets"], 8008 + 120 + 560, "T6.m6.subsets_swept")
    assert mixed_prediction(M34, M34b) == 16
    n36 = check_set_equality(M34, M34b, c_m6, 16, "T6.3x4")
    print(f"[T6] 3x4: m5 0/{c_m5['n_subsets']}, m6 16 prof (4,4) set equality {n36}")

    # ---------- T7: cross-prime field independence ----------
    a_rows_i = [[1, 1, 1, 1], [1, 2, 3, 4]]
    b_rows_i = [[1, -2, 1, 0], [0, 1, -2, 1]]
    totals_wanted = {7: 152, 11: 148, 13: 156, 17: 148, 31: 148}
    prof44_wanted = {7: 8, 11: 4, 13: 12, 17: 4, 31: 4}
    for p_alt in (7, 11, 13, 17, 31):
        Ap = Factor(f"Awit{p_alt}", a_rows_i, p_alt)
        Bp = Factor(f"Bwit{p_alt}", b_rows_i, p_alt)
        assert Ap.spark == 3 and Bp.spark == 3, (p_alt, Ap.spark, Bp.spark)
        assert Ap.circuits == {3: 4} and Bp.circuits == {3: 4}, (p_alt, Ap.circuits, Bp.circuits)
        assert Ap.circuits.get(4, 0) == 0 and Bp.circuits.get(4, 0) == 0
        cnt4, fib4, prof4 = circuit_scan(Ap, Bp, 4)
        assert prof4.get((3, 3), 0) == 144, (p_alt, prof4)
        assert fib4 == 0, (p_alt, fib4)
        assert_eq(cnt4, totals_wanted[p_alt], f"T7.p{p_alt}.total")
        assert_eq(prof4.get((4, 4), 0), prof44_wanted[p_alt], f"T7.p{p_alt}.prof44")
        RESULTS[f"T7.p{p_alt}"] = {"m4_total": cnt4, "prof33": prof4.get((3, 3)),
                                   "prof44": prof4.get((4, 4), 0),
                                   "fibers": fib4, "circA": {str(k): v for k, v in Ap.circuits.items()},
                                   "circB": {str(k): v for k, v in Bp.circuits.items()}}
        print(f"[T7] GF({p_alt}): spectra {{3:4}} both, (3,3)=144, total {cnt4} (4,4)={prof4.get((4,4),0)}")

    # ---------- T8: concat guard ----------
    def three_factor_true(s_tuple, t_tuple, p=13):
        factors = []
        for (s, t) in zip(s_tuple, t_tuple):
            r = s - t
            factors.append(Factor(f"RSs{s}t{t}", vand_rows(r, s, t, p), p))
        d = min(f.spark for f in factors)
        cols = {}
        ranges = [range(f.s) for f in factors]
        for combo in itertools.product(*ranges):
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
        dim = len(next(iter(cols.values())))
        return count, factors, d, dim

    def three_factor_concat_buggy(s_tuple, t_tuple, p=13):
        """The run-3 defect path: block concatenation instead of Kronecker product."""
        factors = []
        for (s, t) in zip(s_tuple, t_tuple):
            r = s - t
            factors.append(Factor(f"RSs{s}t{t}", vand_rows(r, s, t, p), p))
        d = min(f.spark for f in factors)
        cols = {}
        for combo in itertools.product(*[range(f.s) for f in factors]):
            vec = [0] * sum(f.r for f in factors)  # WRONG: concat dims
            off = 0
            for f, idx in zip(factors, combo):
                for i in range(f.r):
                    vec[off + i] = f.cols[idx][i]
                off += f.r
            cols[combo] = vec
        keys = sorted(cols)
        N = len(keys)
        count = 0
        for idx in itertools.combinations(range(N), d):
            if dep_k([cols[keys[i]] for i in idx], d, p):
                count += 1
        dim = len(next(iter(cols.values())))
        return count, d, dim

    anchors = [((3, 3, 4), (1, 1, 1), 24), ((3, 4, 4), (1, 1, 1), 16),
               ((3, 3, 3), (1, 1, 1), 27), ((4, 4, 4), (1, 1, 1), 48)]
    concat_rows = {}
    for (s_t, t_t, want) in anchors:
        got, factors, d, dim_true = three_factor_true(s_t, t_t, p13)
        assert_eq(got, want, f"T8.{s_t}t{t_t}.census")
        # fiber formula from factor data
        pred = 0
        for i, f in enumerate(factors):
            others = 1
            for j, g in enumerate(factors):
                if j != i:
                    others *= g.s
            pred += others * f.circuits.get(d, 0)
        assert_eq(pred, want, f"T8.{s_t}t{t_t}.fiber_formula")
        got_b, d_b, dim_buggy = three_factor_concat_buggy(s_t, t_t, p13)
        assert dim_buggy != dim_true, f"T8 concat guard: dims equal at {s_t}"
        concat_rows[str(s_t)] = {"true_dim": dim_true, "buggy_dim": dim_buggy,
                                 "true_count": got, "buggy_count": got_b}
        print(f"[T8] {s_t}t{t_t}: census {got} = formula {want}; concat-buggy dim {dim_buggy} != true {dim_true} (count {got_b})")
    RESULTS["T8_concat_guard"] = concat_rows

    # ---------- T9: regrouping demo ----------
    # A'' = 2x3 V (spark 3, C(3)=1); B_comp = H2 (x) H3 = 2x3 V (x) 3x4 MDS
    # composite spark must be min(3,4)=3 by P1-IH
    H2 = Factor("H2", vand_rows(2, 3, 0, p13), p13)
    H3 = Factor("H3", vand_rows(3, 4, 1, p13), p13)
    # build composite factor B_comp explicitly: columns are h2 (x) h3
    comp_cols = {}
    for u in range(H2.s):
        for v in range(H3.s):
            comp_cols[(u, v)] = [H2.cols[u][i] * H3.cols[v][j] % p13
                                 for i in range(H2.r) for j in range(H3.r)]
    # composite as rows for Factor: need r = H2.r*H3.r = 2*3 = 6, s = H2.s*H3.s = 12
    comp_rows = [[comp_cols[(u, v)][i] for (u, v) in sorted(comp_cols)]
                 for i in range(H2.r * H3.r)]
    Bcomp = Factor("Bcomp", comp_rows, p13)
    assert Bcomp.spark == 3, Bcomp.spark
    CB = Bcomp.circuits.get(3, 0)
    CB4 = Bcomp.circuits.get(4, 0)
    A2 = Factor("A2", vand_rows(2, 3, 0, p13), p13)
    c_t9 = census(A2, Bcomp, 4)
    pred = mixed_prediction(A2, Bcomp)
    assert_eq(c_t9["profiles"].get((3, 3), 0), pred, "T9.m4.prof33_predicted")
    assert_eq(c_t9["fibers"], 3 * CB4, "T9.m4.fibers_predicted")
    assert_eq(sum(v for (x, y), v in c_t9["profiles"].items() if (x, y) != (3, 3)), 3 * CB4,
              "T9.m4.no_other_prof_besides_fibers")
    assert_eq(c_t9["circuits"], pred + 3 * CB4, "T9.m4.total_predicted")
    n_t9 = check_set_equality(A2, Bcomp, c_t9, pred, "T9.composite")
    assert_eq(c_t9["n_subsets"], 58905 + 630, "T9.m4.subsets_swept")
    RESULTS["T9"] = {"Bcomp_spark": Bcomp.spark, "Bcomp_circuits": {str(k): v for k, v in Bcomp.circuits.items()},
                     "pred_mixed": pred, "measured": n_t9}
    print(f"[T9] composite-B spark {Bcomp.spark} (=min=P1-IH), C_B(3)={CB}, mixed {n_t9} = predicted {pred}, fibers {c_t9['fibers']} = 3*C_B(4)")

    # ---------- finalize ----------
    RESULTS["census_details"] = {
        "T2.m3": {k: v for k, v in c_m3.items() if k != "circuit_sets"},
        "T2.m4": {k: v for k, v in c_m4.items() if k != "circuit_sets"},
        "T3b.m4": {k: v for k, v in c_tb.items() if k != "circuit_sets"},
        "T4.m4": {k: v for k, v in c_2x5_4.items() if k != "circuit_sets"},
        "T5.m4": {k: v for k, v in c_2x6.items() if k != "circuit_sets"},
        "T6.m5": {k: v for k, v in c_m5.items() if k != "circuit_sets"},
        "T6.m6": {k: v for k, v in c_m6.items() if k != "circuit_sets"},
        "T9.m4": {k: v for k, v in c_t9.items() if k != "circuit_sets"},
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


def circuit_scan(fA, fB, m):
    """Count/profile circuits only (cross-prime, no below-d pass)."""
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
