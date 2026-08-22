"""Standalone verification of proofs/alternation_law.md (wave-16 lead takeover).

Re-derives, with no producer imports: the closure dims at L=2,3; the fiber law
(lambda_z = J xor z, affine hyperplanes) at L=2,3 in one encoding and at L=4 from the
STORED certificate basis in results/algebra/char0_2x4.json (second encoding); the c(z)
checkerboard fits; the z=J exclusion/inclusion rule; and the dimension law at all
certified points. Run: PYTHONPATH=src python tests/test_alternation_law.py
"""
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def generators(L):
    n = 2 * L
    gens = [1 << (n + i) for i in range(n)]
    for c in range(L - 1):
        for row in (0, 1):
            gens.append((1 << (c * 2 + row)) | (1 << ((c + 1) * 2 + row)))
    for c in range(L):
        gens.append((1 << (c * 2)) | (1 << (c * 2 + 1)))
    return gens


def closure(gens, n):
    S = set(gens)
    frontier = sorted(S)
    while frontier:
        Slist = sorted(S)
        new = set()
        xw = [((v >> n) & ((1 << n) - 1)) for v in Slist]
        zw = [(v & ((1 << n) - 1)) for v in Slist]
        for v in frontier:
            xv, zv = (v >> n) & ((1 << n) - 1), v & ((1 << n) - 1)
            for xo, zo, w in zip(xw, zw, Slist):
                if (bin(xv & zo).count("1") + bin(zv & xo).count("1")) % 2 == 1:
                    s = v ^ w
                    if s not in S:
                        new.add(s)
        S |= new
        frontier = sorted(new)
    return S


def runs():
    ok = []

    def check(name, cond):
        ok.append((name, bool(cond)))
        print(("  ok: " if cond else "FAIL: ") + name)

    # --- L=2,3 exact closures and dim law -------------------------------------
    sets = {}
    for L, dim_expect in ((2, 56), (3, 1056)):
        n = 2 * L
        S = closure(generators(L), n)
        sets[L] = S
        check(f"L={L} closure dim == {dim_expect}", len(S) == dim_expect)
        check(f"L={L} dim law 2^(n-1)(2^(n-1)-(-1)^L)",
              len(S) == (1 << (n - 1)) * ((1 << (n - 1)) - ((-1) ** L)))

    # --- fiber law at L=2,3 (this test's encoding: x high, z low) --------------
    for L in (2, 3):
        n = 2 * L
        Jm = (1 << n) - 1
        byz = defaultdict(list)
        for v in sets[L]:
            byz[v & Jm].append(v >> n)
        # z=J rule
        check(f"L={L} z=J present iff L odd", (Jm in byz) == (L % 2 == 1))
        if L % 2 == 1:
            check("L=3 z=J fiber is FULL (all 64 x)", len(byz.get(Jm, [])) == (1 << n))
        # all other fibers: affine hyperplane with lambda_z = J xor z, size 2^(n-1)
        good = 0
        for z, xs in byz.items():
            if z == Jm:
                continue
            lam = Jm ^ z
            c = bin(lam & xs[0]).count("1") % 2
            if len(xs) == (1 << (n - 1)) and all(
                    bin(lam & x).count("1") % 2 == c for x in xs):
                good += 1
        nfib = (1 << (n - 1)) - 1  # even-z patterns minus z=J
        check(f"L={L} fiber law on {good}/{nfib} fibers", good == nfib)
        # c(z) = 1 + z.d with checkerboard d (d and J xor d coincide on even z)
        cand = set()
        for dvec in range(1 << n):
            if all((1 ^ (bin(z & dvec).count("1") % 2)) ==
                   bin((Jm ^ z) & xs0).count("1") % 2
                   for z in byz if z != Jm
                   for xs0 in [byz[z][0]]):
                cand.add(dvec)
        # bipartition color classes of the ladder in this test's labeling:
        # site i = col i//2, row i%2; colour = (col + row) mod 2
        cls_a = sum(1 << i for i in range(n) if ((i // 2) + (i % 2)) % 2 == 0)
        check(f"L={L} some affine d fits and a checkerboard class is among them",
              bool(cand) and bool({cls_a, Jm ^ cls_a} & cand))

    # --- fiber law at L=4 from the STORED certificate basis (second encoding) ---
    art = os.path.join(ROOT, "results", "algebra", "char0_2x4.json")
    data = json.load(open(art))
    basis = data["cases"]["main_grid_2x4"]["string_basis"]
    check("stored 2x4 basis has 16256 strings", len(basis) == 16256)
    byz = defaultdict(list)
    for v in basis:  # stored encoding: z high 8, x low 8
        byz[(v >> 8) & 0xFF].append(v & 0xFF)
    check("stored 2x4: 127 distinct even-z, z=255 absent",
          len(byz) == 127 and 255 not in byz
          and all(bin(z).count("1") % 2 == 0 for z in byz))
    good = 0
    for z, xs in byz.items():
        lam = 255 ^ z
        c = bin(lam & xs[0]).count("1") % 2
        if len(xs) == 128 and all(bin(lam & x).count("1") % 2 == c for x in xs):
            good += 1
    check("stored 2x4 fibers all affine lambda_z = J xor z, |fiber|=128", good == 127)
    check("stored 2x4 dim law 128*127", 128 * 127 == 16256)

    # --- c(z) = 1 + z.d checkerboard consistency at L=4 (existence over GF2) ----
    # affine-linear system on (1, z0..z7); require at least one consistent solution
    # and that the checkerboard functional is one of them
    zs = sorted(byz)
    # checkerboard classes for the stored (e148) labeling, derived from ITS bond list
    bonds148 = [[0, 4], [0, 1], [1, 5], [1, 2], [2, 6], [2, 3], [3, 7],
                [4, 5], [5, 6], [6, 7]]
    color = [None] * 8
    color[0] = 0
    stack = [0]
    while stack:
        u = stack.pop()
        for a, b in bonds148:
            if a == u or b == u:
                v = b if a == u else a
                if color[v] is None:
                    color[v] = 1 - color[u]
                    stack.append(v)
    classes = [sum(1 << i for i in range(8) if color[i] == t) for t in (0, 1)]
    fits = any(
        all((1 ^ (bin(z & dvec).count("1") % 2)) ==
            (bin((255 ^ z) & byz[z][0]).count("1") % 2) for z in zs)
        for dvec in classes)
    check("stored 2x4: c(z) = 1 + z.d for a checkerboard d", fits)

    # --- certified type data cross-check against artifacts -----------------------
    cc = data["cases"]
    check("2x4 type so128+so128 per certificate",
          "so(128)" in cc["main_grid_2x4"]["decomposition_statement"])
    check("2x3 type sp32+sp32 per certificate",
          "sp(32)" in cc["control_grid_2x3"]["decomposition_statement"])
    check("2x2 type so8+so8 per certificate",
          "so(8)" in cc["control_grid_2x2"]["decomposition_statement"])

    # --- prediction table recorded in the artifact -------------------------------
    law = json.load(open(os.path.join(ROOT, "results", "algebra",
                                      "alternation_law.json")))
    pts = {p["L"]: p for p in law["points"]}
    check("L=5 point = 262656 (modular probe) recorded as PROBE",
          pts[5]["dim"] == 262656 and "PROBE" in pts[5]["status"])
    check("L=6 prediction 4192256 = 2*dim so_2048 recorded",
          pts[6]["dim"] == 4192256 and 2048 * 2047 == 4192256)
    check("L=7 prediction 67117056 = 2*dim sp_8192 recorded",
          pts[7]["dim"] == 67117056 and 2 * 4096 * 8193 == 67117056)

    n_ok = sum(1 for _, c in ok if c)
    print(f"\n{n_ok}/{len(ok)} checks passed")
    assert n_ok == len(ok), [name for name, c in ok if not c]


if __name__ == "__main__":
    runs()
    print("PASS")
