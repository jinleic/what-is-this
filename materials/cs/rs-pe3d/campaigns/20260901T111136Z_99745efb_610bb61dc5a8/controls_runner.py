#!/usr/bin/env python3
"""Run 3 (gate H-MINLINE-DGE3-PROOF) exhaustive matrix-model controls, exact GF(p).

Self-contained stdlib instrument. All arithmetic exact (Python ints mod p,
pow(x,-1,p)); NO floats, NO numpy. Runs INSIDE the minted run dir after
init; the campaign verdict rests only on THIS in-run execution (prereg
section 0). Emits controls_results.json.

Families (prereg section 3):
  M1 fiber-regime structured GF(13)    (6 configs, expected carrier counts)
  M2 diagonal-regime structured GF(13) (7 configs, P4 pair formula + split)
  M3 randomized non-MDS GF(7)          (21 configs, measured sparks)
  M4 planted controls                  (non-tensor column, duplicated factor column)
  M5 RS census reconciliation          (15 rows through the parity-check model)
  M6 boundary sharpness                (>=2 spark-2 factors => >=1 diagonal;
                                        <=1 spark-2 factor => 0 diagonals)
"""
import itertools
import json
import math
import random
import sys
import time

P13 = 13
P7 = 7
RESULT_FILE = "controls_results.json"


# ---------------- exact GF(p) linear algebra (self-contained) -------------
def rref(rows, ncols, p):
    M = [list(r) for r in rows]
    r = 0
    for c in range(ncols):
        pr = -1
        for i in range(r, len(M)):
            if M[i][c] % p:
                pr = i
                break
        if pr < 0:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % p for j in range(ncols)]
        r += 1
        if r == len(M):
            break
    return r, M


def rank(rows, ncols, p):
    return rref(rows, ncols, p)[0] if rows else 0


def columns_dependent(cols, p):
    """cols: list of coordinate vectors (same length). Dependent iff rank < #cols."""
    if not cols:
        return False
    n = len(cols)
    rr = rank([[c[i] for c in cols] for i in range(len(cols[0]))], n, p)
    return rr < n


def spark_exhaustive(M, ncols, p):
    """Least size of a dependent column set; None if all independent."""
    rows_n = len(M)
    for k in range(2, ncols + 1):
        for S in itertools.combinations(range(ncols), k):
            sub = [[M[i][j] for j in S] for i in range(rows_n)]
            if rank(sub, k, p) < k:
                return k
    return None


def dependent_subsets_below(M, ncols, p, kmax):
    """All dependent column subsets of size <= kmax (for plant detection)."""
    rows_n = len(M)
    out = []
    for k in range(2, kmax + 1):
        for S in itertools.combinations(range(ncols), k):
            sub = [[M[i][j] for j in S] for i in range(rows_n)]
            if rank(sub, k, p) < k:
                out.append(S)
    return out


# ---------------- product model ------------------------------------------
def kron(A, B, p):
    ra, ca = len(A), len(A[0])
    rb, cb = len(B), len(B[0])
    return [[A[i][a] * B[j][b] % p for a in range(ca) for b in range(cb)]
            for i in range(ra) for j in range(rb)]


def kron_many(mats, p):
    out = mats[0]
    for B in mats[1:]:
        out = kron(out, B, p)
    return out


def flat_to_multi(f, dims):
    out = []
    for dd in reversed(dims):
        out.append(f % dd)
        f //= dd
    return tuple(reversed(out))


def prod_col(mats, mu, p):
    v = [1]
    for m, j in zip(mats, mu):
        col = [m[i][j] % p for i in range(len(m))]
        v = [a * b % p for a in v for b in col]
    return v


def measure_sparks(mats, p):
    return [spark_exhaustive(m, len(m[0]), p) for m in mats]


def parallel_classes(M, p):
    """Partition of column indices into parallel classes (nonzero columns)."""
    ncols = len(M[0])
    cols = [[M[i][j] % p for i in range(len(M))] for j in range(ncols)]
    assert all(any(c) for c in cols), "zero column in factor (hypothesis!)"
    classes = []
    assigned = [False] * ncols
    for j in range(ncols):
        if assigned[j]:
            continue
        cls = [j]
        assigned[j] = True
        for k in range(j + 1, ncols):
            if not assigned[k] and columns_dependent([cols[j], cols[k]], p):
                cls.append(k)
                assigned[k] = True
        classes.append(cls)
    return classes


def is_fiber(mats, dims, sparks, p, mus, d):
    """mus: list of d tuples. Fiber iff agrees off some axis ax with a
    size-d dependent (circuit) column set of factor ax (distinct axvals)."""
    n = len(mats)
    for ax in range(n):
        others = [q for q in range(n) if q != ax]
        agree = all(all(mu[q] == mus[0][q] for q in others) for mu in mus)
        if not agree or sparks[ax] != d:
            continue
        axvals = tuple(mu[ax] for mu in mus)
        if len(set(axvals)) != len(axvals):
            continue
        colvecs = [[mats[ax][i][j] for i in range(len(mats[ax]))]
                   for j in axvals]  # list of COLUMN vectors
        if columns_dependent(colvecs, p):
            return True
    return False


def census_planted(columns, p, name, d):
    """Exhaustive weight-d census over an EXPLICIT column list (arbitrary
    vectors, not assumed tensor-structured); classifies fibers via the
    factor model to expose off-fiber carriers. Used by the M4 plant."""
    N = len(columns)
    carriers, fibers, offs = [], [], []
    ncols = len(columns[0])
    for S in itertools.combinations(range(N), d):
        cols = [columns[j] for j in S]
        if not columns_dependent(cols, p):
            continue
        carriers.append(S)
        # fiber test (per factor model of the PRISTINE grid the columns came
        # from): use is_fiber on the two pristine 5x5 factors by mapping
        # flat -> tuple; the planted column at flat 0 breaks tensorhood, but
        # the fiber test only reads POSITIONS, so it still applies.
        mus = [flat_to_multi(j, [5, 5]) for j in S]
        sparks = [3, 3]
        (fibers if is_fiber(FACTOR_CTX["mats"], [5, 5], sparks, p, mus, d)
         else offs).append(S)
    return {"name": name, "carriers": len(carriers), "fiber": len(fibers),
            "off": len(offs),
            "off_sets": [tuple(s) for s in offs]}


FACTOR_CTX = {"mats": None}


def census_and_check(mats, p, name, d=None):
    """Exhaustive weight-d carrier census + regime checks for a product.

    Returns dict with sparks, d, carriers, fibers, off, and classification
    per the theorem. Asserts spark(product)=d via below-d emptiness (full)
    plus at least one dependent d-set.
    """
    prod = kron_many(mats, p)
    dims = [len(m[0]) for m in mats]
    N = 1
    for dd in dims:
        N *= dd
    sparks = measure_sparks(mats, p)
    assert all(s is not None for s in sparks), (name, "factor sparks None")
    if d is None:
        d = min(sparks)
    assert min(sparks) == d, (name, sparks, d)
    # (A): no dependent subset of size < d at all — exhaustive
    for k in range(2, d):
        for S in itertools.combinations(range(N), k):
            cols = [prod_col(mats, flat_to_multi(j, dims), p) for j in S]
            if columns_dependent(cols, p):
                return {"name": name, "FAIL": f"below-d dependent set {S}",
                        "sparks": sparks, "d": d}
    # at least one dependent d-set => spark(product) = d exactly
    found = None

    for S in itertools.combinations(range(N), d):
        cols = [prod_col(mats, flat_to_multi(j, dims), p) for j in S]
        if columns_dependent(cols, p):
            found = S
            break
    assert found is not None, (name, "no dependent d-set; spark> d?")
    carriers, fibers, offs = [], [], []
    for S in itertools.combinations(range(N), d):
        mus = [flat_to_multi(j, dims) for j in S]
        cols = [prod_col(mats, mu, p) for mu in mus]
        if not columns_dependent(cols, p):
            continue
        carriers.append(S)
        (fibers if is_fiber(mats, dims, sparks, p, mus, d) else offs).append(S)
    return {"name": name, "dims": dims, "sparks": sparks, "d": d,
            "carriers": len(carriers), "fiber": len(fibers), "off": len(offs),
            "carrier_sets": [list(s) for s in carriers],
            "fiber_sets": [list(s) for s in fibers],
            "off_sets": [list(s) for s in offs]}


def pair_census(mats, p, name):
    """d=2 exhaustive census + P4 formula check with MEASURED multiplicities."""
    res = census_and_check(mats, p, name, d=2)
    dims = res["dims"]
    sparks = res["sparks"]
    cls = [parallel_classes(m, p) for m in mats]
    mult = [[len(c) for c in cc] for cc in cls]
    A = [sum(mm * mm for mm in mm_list) for mm_list in mult]
    N = 1
    for dd in dims:
        N *= dd
    expected_pairs = (math.prod(A) - N) // 2
    expected_fibers = 0
    for j in range(len(dims)):
        rest = math.prod(dims[i] for i in range(len(dims)) if i != j)
        expected_fibers += rest * sum(m * (m - 1) // 2 for m in mult[j])
    assert res["carriers"] == expected_pairs, (
        name, res["carriers"], expected_pairs, "P4 pair formula MISMATCH")
    # fiber split must match the theorem's split exactly
    assert res["fiber"] == expected_fibers, (
        name, res["fiber"], expected_fibers, "fiber subcount MISMATCH")
    res["A"] = A
    res["class_multiplicities"] = mult
    res["expected_pairs"] = expected_pairs
    res["expected_fibers"] = expected_fibers
    res["expected_off"] = expected_pairs - expected_fibers
    return res


def fiber_regime_census(mats, p, name, d):
    """Fiber regime: expected carriers = sum over argmin axes of
    (# size-d circuits) x prod_{j!=i} s_j; ALL must be fibers (off=0)."""
    res = census_and_check(mats, p, name, d=d)
    dims = res["dims"]
    expected = 0
    for i in range(len(dims)):
        if res["sparks"][i] != d:
            continue
        n_circ = 0
        # circuit test = dependent AND every proper subset independent
        for S in itertools.combinations(range(dims[i]), d):
            # list of COLUMN vectors of the chosen factor columns
            colvecs = [[mats[i][r][j] for r in range(len(mats[i]))] for j in S]
            if not columns_dependent(colvecs, p):
                continue
            minimal = True
            for k in range(1, d):
                for T in itertools.combinations(S, k):
                    sub = [[mats[i][r][j] for r in range(len(mats[i]))]
                           for j in T]
                    if columns_dependent(sub, p):
                        minimal = False
                        break
                if not minimal:
                    break
            if minimal:
                n_circ += 1
        rest = math.prod(dims[j] for j in range(len(dims)) if j != i)
        expected += n_circ * rest
    assert res["carriers"] == expected, (
        name, res["carriers"], expected, "P4 fiber-regime formula MISMATCH")
    assert res["off"] == 0, (name, res["off"], "off-fiber carrier in fiber regime")
    res["expected_carriers"] = expected
    return res


# ---------------- factor builders ----------------------------------------
def vandermonde(p, ncols, nrows=None):
    """nrows x ncols: col j = (1, a_j, ..., a_j^{nrows-1}); a_j distinct, nonzero."""
    if nrows is None:
        nrows = ncols - 1
    alphas = [(j % (p - 1)) + 1 for j in range(ncols)]
    assert len(set(alphas)) == ncols, (p, ncols)
    return [[pow(a, k, p) for a in alphas] for k in range(nrows)]


def row1_distinct(p, ncols):
    """1 x ncols all-distinct nonzero => spark 2, ONE class of size ncols."""
    row = [(j % (p - 1)) + 1 for j in range(ncols)]
    assert len(set(row)) == ncols
    return [row]


def row1_multiclass(p):
    """Pinned layout [c, 2c, e]: classes {0,1} (size 2) and {2} (size 1)."""
    c, e = 1, 3
    assert (2 * c) % p != c % p and e % p not in (c % p, (2 * c) % p)
    return [[c, (2 * c) % p, e % p]]


def rs_parity_check(p, s, t):
    """(s-t) x s Vandermonde block, row exponents t..s-1, points 1..s."""
    alphas = list(range(1, s + 1))
    return [[pow(a, k, p) for a in alphas] for k in range(t, s)]


# ---------------- control families ---------------------------------------
def run():
    rng = random.Random(20260901)
    out = {"meta": {
        "protocol": "T-DGE exhaustive matrix-model controls, in-run",
        "gate": "H-MINLINE-DGE3-PROOF",
        "run_id": "20260901T111136Z_99745efb_610bb61dc5a8",
        "p_structured": P13, "p_random": P7, "seed": 20260901,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "note": "Assignment-supplied targets pinned in the committed prereg; "
                "in-run results are the only evidential record (prereg s0)."
    }, "M1": [], "M2": [], "M3": [], "M4": [], "M5": [], "M6": [],
        "verdict_fields": {}}
    all_ok = True

    # ---- M1: fiber regime structured, GF(13) -----------------------------
    m1 = [([(3, 5), (3, 5)], 100), ([(3, 5), (4, 5)], 50),
          ([(3, 4), (3, 4), (3, 4)], 192),
          ([(3, 4), (4, 4), (4, 4)], 64),
          ([(4, 6), (4, 6)], 180), ([(4, 6), (5, 6)], 90)]
    for spec, expect in m1:
        mats = [vandermonde(P13, c, sp - 1) for sp, c in spec]
        d = min(sp for sp, _ in spec)
        res = fiber_regime_census(mats, P13, f"M1 spec={spec}", d)
        assert res["carriers"] == expect, (spec, res["carriers"], expect)
        assert res["off"] == 0 and not res["off_sets"]
        out["M1"].append({k: v for k, v in res.items()
                          if k not in ("carrier_sets", "fiber_sets", "off_sets")})
        print(f"M1 {spec}: d={d} carriers={res['carriers']} (expect {expect}) "
              f"fiber={res['fiber']} off=0 OK")
        all_ok = all_ok and res["carriers"] == expect

    # ---- M2: diagonal regime structured, GF(13), P4 pair formula ---------
    m2 = [
        ("a", [row1_distinct(P13, 5), row1_distinct(P13, 5)],
         "one class of 5 per factor"),
        ("b", [row1_distinct(P13, 4), row1_distinct(P13, 4)],
         "Main anchor [2,2](4,4): 120 = 48 fibers + 72 off"),
        ("c", [row1_distinct(P13, 5), row1_distinct(P13, 5),
               vandermonde(P13, 4, 3)],
         "two degenerate + one clean"),
        ("d", [row1_distinct(P13, 5), vandermonde(P13, 5, 3)],
         "single degenerate: fibers only"),
        ("e", [row1_distinct(P13, 4), vandermonde(P13, 4, 3),
               vandermonde(P13, 4, 3)],
         "single degenerate among three"),
        ("f", [row1_multiclass(P13), vandermonde(P13, 4, 3)],
         "multi-class spark-2 factor, single degenerate"),
        ("g", [row1_multiclass(P13), row1_distinct(P13, 4)],
         "multi-class x one-class, both degenerate: 34 pairs, off > 0"),
    ]
    for tag, mats, note in m2:
        res = pair_census(mats, P13, f"M2{tag} {note}")
        ndeg = sum(1 for m in mats if spark_exhaustive(m, len(m[0]), P13) == 2)
        res["num_spark2_factors"] = ndeg
        if ndeg <= 1:
            assert res["off"] == 0, (tag, res["off"])
        out["M2"].append({k: v for k, v in res.items()
                          if k not in ("carrier_sets", "fiber_sets", "off_sets")})
        print(f"M2{tag}: N={math.prod(res['dims'])} pairs={res['carriers']} "
              f"(formula {res['expected_pairs']}) fiber={res['fiber']} "
              f"off={res['off']} ndeg={ndeg} OK")
        all_ok = all_ok and res["carriers"] == res["expected_pairs"]

    # ---- M3: randomized non-MDS battery, GF(7) ---------------------------
    violations = 0
    n_fiber = n_diag = 0
    draws = 0
    trial = 0
    while trial < 21 and draws < 200:
        draws += 1
        nf = rng.choice([2, 3, 4])
        dims = [rng.randint(2, 4) for _ in range(nf)]
        mats = []
        for dd in dims:
            nr = rng.randint(1, dd)
            m = [[rng.randrange(P7) for _ in range(dd)] for _ in range(nr)]
            for j in range(dd):
                while not any(m[i][j] % P7 for i in range(nr)):
                    m[rng.randrange(nr)][j] = rng.randrange(1, P7)
            if nr < dd:
                for _ in range(4000):
                    if spark_exhaustive(m, dd, P7) == nr + 1:
                        break
                    m = [[rng.randrange(P7) for _ in range(dd)]
                         for _ in range(nr)]
                    for j in range(dd):
                        while not any(m[i][j] % P7 for i in range(nr)):
                            m[rng.randrange(nr)][j] = rng.randrange(1, P7)
                else:
                    raise AssertionError("could not sample spark")
            mats.append(m)
        sparks = measure_sparks(mats, P7)
        # A factor with nr == dd (square, full rank) can be ALL-independent
        # (spark None = infinity) in a random draw: the theorem's hypothesis
        # wants sparks that are reached; exclude None-flagged matrices by
        # resampling until every factor spark is finite (recorded honestly).
        if any(s is None for s in sparks):
            continue  # invalid draw (infinite-spark factor): resample
        trial += 1  # a valid (finite-spark) configuration: recorded
        dmin = min(sparks)
        assert dmin is not None and dmin >= 2, (trial, sparks)
        if dmin == 2:
            res = pair_census(mats, P7, f"M3 trial{trial}")
            ndeg = sum(1 for s in sparks if s == 2)
            res["num_spark2_factors"] = ndeg
            if ndeg >= 2:
                n_diag += 1
                if res["off"] < 1:
                    violations += 1  # M6: must have >= 1 diagonal
            else:
                n_fiber += 1
                if res["off"] != 0:
                    violations += 1
        else:
            res = census_and_check(mats, P7, f"M3 trial{trial}", d=dmin)
            assert res["off"] == 0, (trial, res["off"])
            n_fiber += 1
        out["M3"].append({k: v for k, v in res.items()
                          if k not in ("carrier_sets", "fiber_sets", "off_sets")})
        print(f"M3 trial{trial}: nf={nf} dims={dims} sparks={sparks} "
              f"d={res['d']} carriers={res['carriers']} "
              f"fiber={res['fiber']} off={res['off']}")
    assert trial == 21, trial
    assert violations == 0, violations

    # ---- M4: planted controls --------------------------------------------
    m4 = {}
    # (a) non-tensor column plant: replace the product column at flat 0 by
    # the SUM of two OTHER product columns (flats 13 and 21). The planted
    # matrix is no longer a tensor product; the triple {0,13,21} is a
    # dependent d-set that is NOT an axis fiber, so the machinery (which
    # asserts off == 0 in the fiber regime) must REJECT it, while the
    # pristine matrix passes the identical path.
    base = [vandermonde(P13, 5, 2), vandermonde(P13, 5, 2)]
    prod = kron_many(base, P13)
    dims = [5, 5]
    clean_below = dependent_subsets_below(prod, 25, P13, 2)
    assert not clean_below, "pristine matrix must have no dependent pairs"
    pristine = fiber_regime_census(base, P13, "M4 pristine", 3)
    assert pristine["off"] == 0 and pristine["carriers"] == 100, (
        pristine["carriers"], pristine["off"])  # 10 circuits x 5 x 2 axes
    q = prod_col(base, (2, 3), P13)
    r = prod_col(base, (4, 1), P13)
    planted = [list(row) for row in prod]
    for i in range(len(planted)):
        planted[i][0] = (q[i] + r[i]) % P13  # overwrite flat-0 column
    planted_below = dependent_subsets_below(planted, 25, P13, 2)
    assert not planted_below, "plant must not corrupt below-d emptiness"
    # classification of the planted matrix through the SAME machinery: the
    # weight-d census must now see at least one off-fiber dependent triple
    # ({0,13,21}: planted(0) - planted(13) - planted(21) = 0), violating
    # the fiber regime's off == 0 => REJECT.
    FACTOR_CTX["mats"] = base
    plant_cols = [[row[j] for row in planted] for j in range(25)]  # TRUE columns
    post = census_planted(plant_cols, P13, "M4 planted", 3)
    assert post["off"] >= 1, (post["off"], "planted non-tensor NOT rejected")
    assert [0, 13, 21] in [list(s) for s in post["off_sets"]], \
        post["off_sets"][:4]
    m4["nontensor_plant"] = {
        "pristine_below_d_sets": len(clean_below),
        "pristine_carriers": pristine["carriers"], "pristine_off": 0,
        "planted_below_d_sets": len(planted_below),
        "planted_carriers": post["carriers"], "planted_off": post["off"],
        "off_witness": [0, 13, 21],
        "verdict": "REJECT planted as expected (off-fiber carrier appears)"}
    # (b) duplicated factor column plant
    dup = [list(r) for r in vandermonde(P13, 5, 2)]
    dup[0].append(dup[0][0]); dup[1].append(dup[1][0])  # duplicate col 0
    dupspark = spark_exhaustive(dup, 6, P13)
    assert dupspark == 2, dupspark
    m4["duplicated_factor_column"] = {
        "registered_spark": 3, "measured_spark_after_plant": dupspark,
        "verdict": "REJECT: measured spark 2 != registered 3, "
                   "config mismatch fires"}
    # fail-loud inversion check: the pristine anchor passes pair_census
    pristine_res = pair_census([row1_distinct(P13, 4), row1_distinct(P13, 4)],
                               P13, "M4 pristine anchor")
    m4["pristine_control"] = {"carriers": pristine_res["carriers"],
                              "expected": pristine_res["expected_pairs"],
                              "verdict": "ACCEPT as expected"}
    out["M4"] = m4
    print(f"M4 plants: non-tensor REJECT (off={post['off']} >= 1, witness "
          f"[0,13,21]), duplicated-column REJECT (spark 2 != 3), "
          f"pristine ACCEPT")

    # ---- M5: RS census reconciliation, GF(13) ----------------------------
    rows5 = [
        ((2, 3, 4), (1, 1, 1), 12, "fiber"),
        ((2, 4, 4), (1, 1, 1), 16, "fiber"),
        ((3, 3, 3), (1, 1, 1), 27, "fiber"),
        ((3, 3, 4), (1, 1, 1), 24, "fiber"),
        ((3, 4, 4), (1, 1, 1), 16, "fiber"),
        ((4, 4, 4), (1, 1, 1), 48, "fiber"),
        ((3, 3, 3), (2, 1, 1), 27, "fiber_circuit_form"),
        ((4, 4, 4), (2, 1, 1), 64, "fiber_circuit_form"),
        ((2, 2, 4), (2, 1, 1), 16, "trivial_zero_map_count_N"),
        ((2, 3, 4), (2, 2, 2), 24, "trivial_zero_map_count_N"),
        ((2, 2, 3), (1, 1, 1), 18, "diagonal"),
        ((2, 2, 2), (1, 1, 1), 28, "diagonal"),
        ((3, 3, 3), (2, 2, 1), 108, "diagonal_circuit_form"),
    ]
    for s, t, expect, kind in rows5:
        mats = [rs_parity_check(P13, si, ti) for si, ti in zip(s, t)]
        # d_i = s_i - t_i + 1; an EMPTY parity check (s_i = t_i) is the
        # zero-map clause: spark 1, count N (prereg sections 1 and 3).
        # Detected BEFORE measuring sparks (a 0-row matrix has no columns).
        if any(si == ti for si, ti in zip(s, t)):
            N = 1
            for si in s:
                N *= si
            carriers, fiber, off, dd = N, N, 0, 1
            assert N == expect, (s, t, N, expect)
            sparks = [-1 if si == ti else si - ti + 1
                      for si, ti in zip(s, t)]  # -1 marks the zero-map axis
        else:
            sparks = measure_sparks(mats, P13)
            dmin = min(sparks)
            assert dmin >= 2, (s, t, sparks)
            if dmin == 2:
                res = pair_census(mats, P13, f"M5 {s}t{t}")
                carriers, fiber, off, dd = (res["carriers"], res["fiber"],
                                            res["off"], 2)
                assert carriers == expect, (s, t, carriers, expect)
            else:
                res = fiber_regime_census(mats, P13, f"M5 {s}t{t}", dmin)
                carriers, fiber, off, dd = (res["carriers"], res["fiber"],
                                            res["off"], dmin)
                assert carriers == expect, (s, t, carriers, expect)
        out["M5"].append({"s": list(s), "t": list(t),
                          "sparks": sparks, "d": dd, "carriers": carriers,
                          "fiber": fiber, "off": off, "expected": expect,
                          "kind": kind, "pass": carriers == expect})
        print(f"M5 s={s} t={t}: sparks={sparks} carriers={carriers} "
              f"(expect {expect}) off={off} OK")

    out["M5"].append({
        "s": [2, 2, 4], "t": [1, 1, 1], "row": "Run1_anchor",
        "check": "count 24 = 16 lines + 8 diagonals; regime diagonal"})
    mats = [rs_parity_check(P13, si, ti)
            for si, ti in zip((2, 2, 4), (1, 1, 1))]
    res = pair_census(mats, P13, "M5 Run1 anchor (2,2,4)t111")
    assert res["carriers"] == 24, res["carriers"]
    assert res["fiber"] == 16 and res["off"] == 8, (res["fiber"], res["off"])
    # exhibit the 8 off-fiber pairs, map to Run-1 flat-(0,12) pair
    dims = res["dims"]
    off_pairs = []
    for S in itertools.combinations(range(16), 2):
        mus = [flat_to_multi(j, dims) for j in S]
        cols = [prod_col(mats, mu, P13) for mu in mus]
        if columns_dependent(cols, P13) and not is_fiber(
                mats, dims, res["sparks"], P13, mus, 2):
            off_pairs.append(S)
    assert len(off_pairs) == 8, off_pairs
    assert [0, 12] in [list(s) for s in off_pairs], off_pairs[:4]
    out["M5"][-1].update({"carriers": 24, "fiber": 16, "off": 8,
                          "off_pairs_flat": [list(s) for s in off_pairs],
                          "pass": True})
    print(f"M5 Run1 anchor (2,2,4)t111: 24 = 16 lines + 8 diagonals, "
          f"flat-(0,12) pair reproduced, OK")

    # ---- M6: boundary sharpness ------------------------------------------
    m6rows = []
    for tag, mats, note, ndeg_expect in [
        ("a", [row1_distinct(P13, 5), row1_distinct(P13, 5)], "2x degenerate", 2),
        ("b", [row1_distinct(P13, 4), row1_distinct(P13, 4)], "Main anchor", 2),
        ("c", [row1_distinct(P13, 5), row1_distinct(P13, 5),
               vandermonde(P13, 4, 3)], "2x + clean", 2),
        ("g", [row1_multiclass(P13), row1_distinct(P13, 4)], "multi x one", 2),
        ("d", [row1_distinct(P13, 5), vandermonde(P13, 5, 3)],
         "single degenerate", 1),
        ("f", [row1_multiclass(P13), vandermonde(P13, 4, 3)],
         "multi-class single degenerate", 1),
    ]:
        res = pair_census(mats, P13, f"M6{tag}")
        ndeg = sum(1 for m in mats if spark_exhaustive(m, len(m[0]), P13) == 2)
        assert ndeg == ndeg_expect, (tag, ndeg, ndeg_expect)
        if ndeg >= 2:
            assert res["off"] >= 1, (tag, "no diagonal in 2-degenerate config")
            witness = res["off_sets"][0] if res["off_sets"] else None
        else:
            assert res["off"] == 0, (tag, res["off"], "diagonal in single-degenerate")
            witness = None
        m6rows.append({"tag": tag, "note": note, "ndeg": ndeg,
                       "off": res["off"], "witness": witness,
                       "sharp": True, "pass": True})
        print(f"M6{tag}: ndeg={ndeg} off={res['off']} "
              f"{'witness ' + str(witness) if witness else 'off=0'} OK")
    out["M6"] = m6rows

    out["meta"]["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    out["verdict_fields"] = {
        "M1_pass": all(r.get("carriers") == r.get("expected_carriers")
                       and r.get("off") == 0 for r in out["M1"]),
        "M2_pass": all(r.get("carriers") == r.get("expected_pairs")
                       and r.get("fiber") == r.get("expected_fibers")
                       for r in out["M2"]),
        "M3_pass": violations == 0 and len(out["M3"]) == 21,
        "M4_pass": True,
        "M5_pass": all(r.get("pass", True) for r in out["M5"]),
        "M6_pass": all(r.get("pass") for r in out["M6"]),
        "all_ok": all_ok,
    }
    out["verdict_fields"]["all_families_pass"] = all(
        v for k, v in out["verdict_fields"].items() if k != "all_ok")
    assert out["verdict_fields"]["all_families_pass"], out["verdict_fields"]

    with open(RESULT_FILE, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print("ALL SIX CONTROL FAMILIES PASS —", RESULT_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(run())
