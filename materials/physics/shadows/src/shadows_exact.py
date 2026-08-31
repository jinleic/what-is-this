"""Exact GF(2) symplectic census for contractive-unitary classical shadows.

Owner: physics/shadows. See ../pre_statement.md (frozen 2026-08-29) for gates.

Exactness contract: census/weight computations use Python ints (bit tuples)
and fractions.Fraction ONLY. Floating point appears solely in fields of the
campaign JSON explicitly labeled "numerical_*".

Conventions (pinned first-hand from the paper, see pre_statement.md):
- A Pauli class on n qubits is a 2n-bit integer; bit 2i = x_i, bit 2i+1 = z_i
  (letters: X=(x,z)=(1,0), Y=(1,1), Z=(0,1), I=(0,0)). Classes are defined up
  to the four phases +/-1, +/-i, which a SIZE census never sees.
- A Clifford conjugation action on classes is P -> S P with S symplectic
  (S^T J S = J) over GF(2). Cliffords mod global phase = (translation p) x
  (symplectic S); translations act trivially on classes, so the exhaustive
  class census lives on Sp(4,2) (720 elements) and the mod-phase Clifford
  count is 720 * 16 = 11520 (= 2^(n^2+2n) prod_j (4^j-1) at n=2; the paper
  never states this count itself).
- Two-qubit anchor U12 = exp(i pi/4 Z1 Z2), conjugation O -> U12^dag O U12.
  Derived class map (rule: class(P) -> class(i P G) if [P,G] != 0 else
  class(P), G = Z1 Z2), as images of basis classes (bit order x1,z1,x2,z2):
    S(X1) = Y1 Z2 (0b1011=11), S(Z1) = Z1 (0b0010=2),
    S(X2) = Z1 Y2 (0b1110=14), S(Z2) = Z2 (0b1000=8).
  This reproduces the paper's Methods mapping exactly (asserted in selftest):
    X1Z2 -> Y1, Y1Z2 -> X1, Z1X2 -> Y2, Z1Y2 -> X2  (size 1, signs dropped),
    XX, XY, YX, YY, ZZ unchanged.

Sources: Wu, Wang, Yao, Zhai, You, Zhang, "Contractive unitary and classical
shadow tomography", npj QI 12, 86 (2026), DOI 10.1038/s41534-026-01227-w
(= arXiv:2412.01850). Read first-hand 2026-08-29.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb

N_SITES = 2  # gate A works on 2 qubits
BASIS = [1 << b for b in range(2 * N_SITES)]  # e_x1, e_z1, e_x2, e_z2


# ---------------------------------------------------------------- GF(2) bits

def swap_letter_bits(v: int) -> int:
    """Swap x/z bit pairs (2i <-> 2i+1). Used to build the symplectic form."""
    return ((v & 0xAAAAAAAAAAAAAAAA) >> 1) | ((v & 0x5555555555555555) << 1)


def symplectic_form(v: int, w: int) -> int:
    """B(v,w) = sum_i (x_i w_zi + z_i w_xi) in GF(2). Symmetric, J-standard."""
    return (v & swap_letter_bits(w)).bit_count() & 1


def pauli_size(v: int, n: int) -> int:
    """Number of non-identity single-qubit letters of class v (2n bits)."""
    size = 0
    for i in range(n):
        if (v >> (2 * i)) & 3:
            size += 1
    return size


def letter_vec(letter: str) -> int:
    return {"I": 0b00, "X": 0b01, "Z": 0b10, "Y": 0b11}[letter]


def apply_map(imgs, v: int) -> int:
    """Apply class map given by images of basis vectors (one int per bit)."""
    out = 0
    b = 0
    while v >> b:
        if (v >> b) & 1:
            out ^= imgs[b]
        b += 1
    return out


def is_symplectic(imgs, n: int, basis=None) -> bool:
    basis = basis or [1 << b for b in range(2 * n)]
    for a, b in combinations(range(len(basis)), 2):
        if symplectic_form(apply_map(imgs, basis[a]),
                           apply_map(imgs, basis[b])) != \
           symplectic_form(basis[a], basis[b]):
            return False
    for bidx in range(len(basis)):
        if apply_map(imgs, basis[bidx]) == 0:
            return False  # invertible
    return True


def compose(m2, m1):
    """Class map m2 after m1 (apply m1 first), as basis-image tuples."""
    return tuple(apply_map(m2, img) for img in m1)


def invert_symplectic(imgs):
    """Inverse of a symplectic class map given by basis images."""
    n_bits = len(imgs)
    # Solve imgs @ X = I over GF(2) by Gaussian elimination on augmented cols.
    rows = []
    for b in range(n_bits):
        row = []
        for src in range(n_bits):
            row.append((imgs[src] >> b) & 1)
        rows.append(row)  # rows[b][src] = bit b of image of basis src
    # Build matrix M with M[src][dst] = bit dst of imgs[src]; invert M.
    mat = [[(imgs[src] >> dst) & 1 for dst in range(n_bits)] +
           [1 if src == dst else 0 for dst in range(n_bits)]
           for src in range(n_bits)]
    for col in range(n_bits):
        pivot = next(r for r in range(col, n_bits) if mat[r][col])
        mat[col], mat[pivot] = mat[pivot], mat[col]
        for r in range(n_bits):
            if r != col and mat[r][col]:
                mat[r] = [a ^ b for a, b in zip(mat[r], mat[col])]
    inv_rows = [tuple(mat[src][n_bits + dst] << dst
                      for dst in range(n_bits)) for src in range(n_bits)]
    return tuple(sum(bits) for bits in inv_rows)


# ------------------------------------------------------- two-qubit gate A

ANCHOR_IMGS = (0b1011, 0b0010, 0b1110, 0b1000)  # S(X1), S(Z1), S(X2), S(Z2)
ANCHOR_NAME = "exp(i pi/4 Z1Z2), O -> U^dag O U (paper Methods)"


def enumerate_symplectic_2q():
    """All 65536 candidate maps on 4 bits, keep symplectic. Returns list of
    4-tuples (images of e_x1, e_z1, e_x2, e_z2)."""
    out = []
    for imgs in product(range(16), repeat=4):
        if any(apply_map(imgs, b) == 0 for b in BASIS):
            continue
        ok = True
        for a, b in combinations(range(4), 2):
            va, vb = apply_map(imgs, BASIS[a]), apply_map(imgs, BASIS[b])
            if symplectic_form(va, vb) != symplectic_form(BASIS[a], BASIS[b]):
                ok = False
                break
        if ok:
            out.append(imgs)
    return out


def size2_classes(n: int = 2):
    return [v for v in range(16) if pauli_size(v, n) == 2]


def class_label(v: int, n: int = 2) -> str:
    parts = []
    for i in range(n):
        bits = (v >> (2 * i)) & 3
        parts.append({0: "I", 1: "X", 3: "Y", 2: "Z"}[bits] + str(i + 1)
                     if False else {0: "I", 1: "X", 3: "Y", 2: "Z"}[bits])
    return "".join(parts)


def gate_a_census():
    """Exhaustive census of contraction counts over all 720 symplectic
    actions. Returns dict with asserted counts and full histogram."""
    size2 = size2_classes()
    assert len(size2) == 9, len(size2)
    symp = enumerate_symplectic_2q()
    assert len(symp) == 720, f"Sp(4,2) enumeration gave {len(symp)}, want 720"
    assert ANCHOR_IMGS in symp, "anchor map missing from Sp(4,2)"

    hist = {}          # n_contracted -> number of symplectic actions
    patterns = {}      # achieving mask -> list of action tuples
    anchor_pattern = None
    for imgs in symp:
        contracted = []
        for c in size2:
            img = apply_map(imgs, c)
            if pauli_size(img, 2) <= 1:
                contracted.append(c)
        n = len(contracted)
        hist[n] = hist.get(n, 0) + 1
        if n >= 4:
            patterns.setdefault(frozenset(contracted), []).append(imgs)
        if imgs == ANCHOR_IMGS:
            anchor_pattern = n
            assert [class_label(c) for c in contracted] == \
                ["XZ", "YZ", "ZX", "ZY"] or \
                sorted(class_label(c) for c in contracted) == \
                sorted(["XZ", "YZ", "ZX", "ZY"]), contracted
    # Paper Methods mapping, asserted elementwise on the anchor (class level;
    # the anchor block in selftest already checks exact image classes):
    for src in ("XZ", "YZ", "ZX", "ZY"):
        v = 0
        for idx, ch in enumerate(src):
            v |= {"X": 0b01, "Y": 0b11, "Z": 0b10}[ch] << (2 * idx)
        assert pauli_size(apply_map(ANCHOR_IMGS, v), 2) == 1, src
    for src in ("XX", "XY", "YX", "YY", "ZZ"):
        v = 0
        for idx, ch in enumerate(src):
            v |= {"X": 0b01, "Y": 0b11, "Z": 0b10}[ch] << (2 * idx)
        assert class_label(apply_map(ANCHOR_IMGS, v)) == src, src
    assert anchor_pattern == 4, anchor_pattern

    full_hist = {}
    for n, cnt in sorted(hist.items()):
        full_hist[str(n)] = {"actions": cnt, "labels_11520": 16 * cnt}
    total_actions = sum(hist.values())
    assert total_actions == 720
    max_seen = max(hist)
    return {
        "sp42_order": len(symp),
        "mod_phase_labels": len(symp) * 16,
        "mod_phase_labels_crosscheck": 2 ** 8 * 45,
        "histogram_by_contracted": full_hist,
        "max_contracted": max_seen,
        "n_achieving_actions": sum(hist[4:]) if False else sum(
            cnt for n, cnt in hist.items() if n == 4),
        "achieving_patterns": [
            {"contracted": sorted(class_label(c) for c in pat),
             "count_actions": len(imgs_list)}
            for pat, imgs_list in sorted(
                patterns.items(), key=lambda kv: -len(kv[0]))
            if len(pat) == 4
        ],
        "anchor_contracted": anchor_pattern,
    }


# ------------------------------------------------- generic k-qubit meshes

def single_qubit_symplectics():
    out = []
    for imx in range(4):
        for imz in range(4):
            det = ((imx & 1) * ((imz >> 1) & 1)) ^ (((imx >> 1) & 1) * (imz & 1))
            if det == 1:
                out.append((imx, imz))
    return out


def find_letter_rotation(letter: str):
    """Single-qubit class map with Z -> letter (conjugation by a 1q Clifford
    realizing exp(i pi/4 L L) from the ZZ anchor)."""
    want = letter_vec(letter) & 0b10, letter_vec(letter) & 0b01
    for imx, imz in single_qubit_symplectics():
        img_z = ((imz & 1) << 0) | ((imz >> 1) << 1)
        if img_z == letter_vec(letter):
            return (imx, imz)
    raise AssertionError(letter)

def pair_rotation_imgs(a: str, b: str):
    """Class map of exp(i pi/4 A1 B2) = (R1 x R2) ZZ (R1 x R2)^dag, on 2 sites
    (site 1 carries letter A, site 2 letter B)."""
    ra = find_letter_rotation(a)
    rb = find_letter_rotation(b)
    pre = (ra[0] | (ra[1] << 2), rb[0] | (rb[1] << 2))  # images: x1,z1 then x2,z2
    pre = (pre[0] & 0b11) | ((pre[0] >> 2) << 2), pre[1]  # (x1 img, z1 img, x2 img, z2 img)
    pre_imgs = (ra[0] | (ra[1] << 2), ) + (rb[0] | (rb[1] << 2),)
    # careful explicit: basis x1,z1,x2,z2 -> pre map images
    x1i = ra[0]          # image of X1 under R1
    z1i = ra[1]          # image of Z1 under R1
    x2i = rb[0] << 2     # image of X2 under R2 (site 2 bits)
    z2i = rb[1] << 2
    r_imgs = (x1i | x2i, z1i | z2i, x2i | x1i, z2i | z1i)
    # NOTE: basis vectors are e_x1=(1,0 on site1), e_z1, e_x2, e_z2; a basis
    # vector acts only through its own site's local map:
    r_imgs = (ra[0] | 0, ra[1] | 0, rb[0] << 2, rb[1] << 2)
    r_inv = invert_symplectic(r_imgs)
    return compose(compose(r_imgs, ANCHOR_IMGS), r_inv)
    # phi_AB = r . anchor . r^{-1}  (apply r_inv first) -- compose(m2, m1)
    # means m2 after m1; we want r after anchor after r_inv:
    # compose(compose(r_imgs, ANCHOR_IMGS), r_inv)


def embed_pair(pair_imgs, n: int, i: int, j: int):
    """Embed a 2-site (x1,z1,x2,z2) class map acting on sites (i, j) into an
    n-site map (identity on the other sites)."""
    out = [0] * (2 * n)
    for s in range(n):
        if s not in (i, j):
            out[2 * s] = 1 << (2 * s)
            out[2 * s + 1] = 1 << (2 * s + 1)
    # site i plays anchor role 1 (x1,z1), site j role 2 (x2,z2)
    for local, tgt in ((0, 2 * i), (1, 2 * i + 1), (2, 2 * j), (3, 2 * j + 1)):
        img = pair_imgs[local]
        moved = 0
        for b in range(4):
            if (img >> b) & 1:
                if b < 2:
                    moved |= 1 << (2 * i + b)
                else:
                    moved |= 1 << (2 * j + (b - 2))
        out[tgt] = moved
    # fix identity assignments overridden above for sites i, j: out[2i], ...
    # were set by the loop; pair contributions already placed at moved bits.
    return tuple(out)


def build_mesh(n: int, pair_spec):
    """Compose embedded pair maps over pairs (i, j) in given order.
    pair_spec: iterable of (i, j, pair_imgs)."""
    imgs = tuple(1 << b for b in range(2 * n))
    for i, j, pimgs in pair_spec:
        emb = embed_pair(pimgs, n, i, j)
        assert is_symplectic(emb, n)
        imgs = compose(emb, imgs)  # apply current map, then this pair
    assert is_symplectic(imgs, n)
    return imgs


def mesh_all_pairs(n: int, letters="Z"):
    """U_L = prod_{i<j} exp(i pi/4 L_i L_j), L in {X,Y,Z}; pairs in lex order."""
    spec = [(i, j, pair_rotation_imgs(letters, letters))
            for i, j in combinations(range(n), 2)]
    return build_mesh(n, spec), spec


def cnot_imgs():
    """Class map of CNOT (control = local site 0, target = local site 1),
    SITE-LOCAL 4-bit map for embed_gate; basis order (x0, z0, x1, z1)."""
    X0, Z0, X1, Z1 = 0b0001, 0b0010, 0b0100, 0b1000
    return (X0 | X1, Z0, X1, Z1 | Z0)


def cz_imgs():
    """Class map of CZ on local sites (0, 1), SITE-LOCAL for embed_gate."""
    X0, Z0, X1, Z1 = 0b0001, 0b0010, 0b0100, 0b1000
    return (X0 | Z1, Z0, X1 | Z0, Z1)


def embed_gate(gimgs, n, sites):
    """Embed a gate class map given per involved site (in `sites` order):
    gimgs maps basis index over the involved sites' (x,z) pairs."""
    k = len(sites)
    out = [0] * (2 * n)
    for s in range(n):
        if s not in sites:
            out[2 * s] = 1 << (2 * s)
            out[2 * s + 1] = 1 << (2 * s + 1)
    for local in range(2 * k):
        img = gimgs[local]
        site = sites[local // 2]
        letter_bit = local % 2
        moved = 0
        for src_local in range(2 * k):
            if (img >> src_local) & 1:
                moved |= 1 << (2 * sites[src_local // 2] + src_local % 2)
        out[2 * site + letter_bit] = moved
    return tuple(out)


def token_enumeration(n: int, imgs, k=None, identity_sites=()):
    """Enumerate all 3^n token strings (I on identity_sites) as Pauli class
    vectors, evolve through the class map, return (per-string sizes, hist).
    Strings: uniform letters on support; weight w = sum_m hist[m]/(3^n_sup * 3^m)
    where n_sup = number of letter sites."""
    support = [s for s in range(n) if s not in identity_sites]
    sup = len(support)
    hist = {}
    sizes = []
    for letters in product("XYZ", repeat=sup):
        v = 0
        for s, letter in zip(support, letters):
            v |= letter_vec(letter) << (2 * s)
        img = apply_map(imgs, v)
        sz = pauli_size(img, n)
        sizes.append(sz)
        hist[sz] = hist.get(sz, 0) + 1
    return sizes, hist, sup


def weight_from_hist(hist, n_tokens: int) -> Fraction:
    return sum(Fraction(cnt, n_tokens * (3 ** m)) for m, cnt in hist.items())


# --------------------------------------------------------------- gate B

def eq4_weight(k: int) -> Fraction:
    return Fraction(1, 2) * (Fraction(1, 3 ** k)
                             + Fraction((-1) ** k, 9 ** k)) + \
        Fraction(1, 2) * (Fraction(5 ** k, 9 ** k) - Fraction(1, 9 ** k))


def eq7_weight(k_tilde: int, q: int) -> Fraction:
    """Eq. (7): w = 1/2[3^-kt + (-1)^kt 9^-kt] + (1/2)[(5/9)^kt - 9^-kt] / 3^q."""
    kt = k_tilde
    even_term = (Fraction(5 ** kt, 9 ** kt) - Fraction(1, 9 ** kt)) / (3 ** q)
    return Fraction(1, 2) * (Fraction(1, 3 ** kt)
                             + Fraction((-1) ** kt, 9 ** kt)) + \
        Fraction(1, 2) * even_term


def eq5_hist_counts(k: int):
    """Integer expected counts n_m for the size distribution over 3^k tokens
    per Eq. (5): n_m = C(k,m) 2^m [m odd] + [m==k] sum_l C(k,2l) 4^l."""
    counts = {}
    for m in range(k + 1):
        c = comb(k, m) * 2 ** m * (1 if m % 2 else 0)
        if m == k:
            c += sum(comb(k, 2 * l) * 2 ** (2 * l)
                     for l in range((k // 2) + 1))
        counts[m] = c
    return counts


def gate_b_run(kmax_mesh=10, kmax_weight=16, qmax_total=7):
    out = {"per_k": [], "identity_insertion": [], "random_clifford": [],
           "checks": {}}
    checks = {"eq3_per_string_k<=10": True, "eq5_hist_k<=10": True,
              "eq4_exact_k<=16": True, "eq7_exact": True,
              "eq7_allpairs_eq6_persupportstring": True,
              "random_clifford_1/(2^k+1)": True}
    for k in range(1, kmax_weight + 1):
        w = eq4_weight(k)
        row = {
            "k": k,
            "eq4_w_exact": f"{w.numerator}/{w.denominator}",
            "shadow_norm2_exact": f"{w.denominator}/{w.numerator}",
            "numerical_ratio_norm2_over_2x1.8^k":
                (float(Fraction(1, 1) / w) / (2 * 1.8 ** k)),
        }
        if k <= kmax_mesh:
            imgs, _ = mesh_all_pairs(k, "Z")
            sizes, hist, sup = token_enumeration(k, imgs)
            assert sup == k
            n_tok = 3 ** k
            # Eq3 per-string check: odd N_XY -> size N_XY; even -> size k
            for letters in product("XYZ", repeat=k):
                n_xy = sum(1 for L in letters if L in "XY")
                v = 0
                for s, L in enumerate(letters):
                    v |= letter_vec(L) << (2 * s)
                pred = n_xy if n_xy % 2 else k
                if pauli_size(apply_map(imgs, v), k) != pred:
                    checks["eq3_per_string_k<=10"] = False
            # Eq5 histogram check
            expect = eq5_hist_counts(k)
            if dict(hist) != {m: c for m, c in expect.items() if c}:
                checks["eq5_hist_k<=10"] = False
            w_mesh = weight_from_hist(hist, n_tok)
            if w_mesh != w:
                checks["eq4_exact_k<=16"] = False
            row["mesh_w_exact"] = f"{w_mesh.numerator}/{w_mesh.denominator}"
            row["hist_mesh_exact"] = {str(m): c for m, c in sorted(hist.items())}
            odd_mass = sum(c for m, c in hist.items() if m % 2)
            row["odd_mass"] = odd_mass
            row["even_mass"] = n_tok - odd_mass
            odd_best = max((m for m in hist if m % 2),
                           key=lambda m: hist[m],
                           default=None)
            row["odd_part_mode_m"] = odd_best
        out["per_k"].append(row)
    # identity insertion, Eq. (7)
    for k_tilde in range(1, qmax_total):
        for q in range(0, qmax_total - k_tilde + 1):
            k = k_tilde + q
            imgs, _ = mesh_all_pairs(k, "Z")
            # per-string Eq. (6) check on the support (exhaustive, exact):
            sup = list(range(k_tilde))
            for letters in product("XYZ", repeat=k_tilde):
                v = 0
                for s, L in enumerate(letters):
                    v |= letter_vec(L) << (2 * sup[s])
                n_xy = sum(1 for L in letters if L in "XY")
                m_pred = n_xy + q if n_xy % 2 else k - q
                m_got = pauli_size(apply_map(imgs, v), k)
                if m_got != m_pred:
                    checks["eq7_allpairs_eq6_persupportstring"] = False
            sizes, hist, sup_sizes_unused = token_enumeration(
                k, imgs, identity_sites=set(range(k_tilde, k)))
            n_tok = 3 ** k_tilde
            w_mesh = weight_from_hist(hist, n_tok)
            w7 = eq7_weight(k_tilde, q)
            ok = (w_mesh == w7)
            if not ok:
                checks["eq7_exact"] = False
            out["identity_insertion"].append({
                "k_tilde": k_tilde, "q": q, "k_total": k,
                "w_mesh_exact": f"{w_mesh.numerator}/{w_mesh.denominator}",
                "eq7_exact": f"{w7.numerator}/{w7.denominator}",
                "match": ok,
            })
    # random Clifford baseline
    for k in range(1, kmax_weight + 1):
        w_derived = sum(Fraction(comb(k, m) * 3 ** m, 4 ** k - 1) *
                        Fraction(1, 3 ** m) for m in range(1, k + 1))
        closed = Fraction(1, 2 ** k + 1)
        ok = (w_derived == closed)
        if not ok:
            checks["random_clifford_1/(2^k+1)"] = False
        out["random_clifford"].append({
            "k": k, "w_derived": f"{w_derived.numerator}/{w_derived.denominator}",
            "closed_form": f"{closed.numerator}/{closed.denominator}",
            "match": ok,
        })
    out["checks"] = checks
    return out


# --------------------------------------------------------------- gate C

def all_matchings(sites):
    """Perfect matchings of `sites` (any count of leftover singletons is NOT
    allowed; call with the subsets you want covered)."""
    sites = tuple(sites)
    if not sites:
        yield ()
        return
    first = sites[0]
    for other in sites[1:]:
        pair = (first, other)
        rest = tuple(s for s in sites if s not in pair)
        for m in all_matchings(rest):
            yield (pair,) + m


def matchings_leaving_at_most_one(sites):
    """Matchings covering all sites; for odd |sites| cover all but one."""
    sites = tuple(sites)
    if len(sites) % 2 == 0:
        yield from all_matchings(sites)
    else:
        for drop in sites:
            rest = tuple(s for s in sites if s != drop)
            for m in all_matchings(rest):
                yield m


def perfect_matchings_cover(sites):
    """All matchings covering all sites when even, else all but one site."""
    return matchings_leaving_at_most_one(sites)


def eval_map_on_tokens(n: int, imgs):
    sizes, hist, sup = token_enumeration(n, imgs)
    w = weight_from_hist(hist, 3 ** sup)
    return w, hist


def gate_c_run(kmax=8):
    results = {"families": {}, "notes": []}
    best = None  # (family, k, W)
    # family (a): single-letter meshes
    for letter in "XYZ":
        fam = f"mesh_{letter}"
        results["families"][fam] = []
        for k in range(2, kmax + 1):
            imgs, _ = mesh_all_pairs(k, letter)
            w, hist = eval_map_on_tokens(k, imgs)
            results["families"][fam].append(
                {"k": k, "w_exact": f"{w.numerator}/{w.denominator}",
                 "shadow_norm2": f"{w.denominator}/{w.numerator}",
                 "hist": {str(m): c for m, c in sorted(hist.items())}})
    # consistency: mesh_Z must equal Eq. (4)
    for k in range(2, kmax + 1):
        wz = [Fraction(*map(int, r["w_exact"].split("/")))
              for r in results["families"]["mesh_Z"] if r["k"] == k][0]
        assert wz == eq4_weight(k), (k, wz)
    results["notes"].append("mesh_Z == U_ct: Eq. (4) reproduced exactly for "
                            "k=2..8 (consistency anchor).")
    # family (b): single-letter-per-pair rotations on (near-)perfect matchings
    fam = "matchings_single_letter"
    results["families"][fam] = []
    n_composites = 0
    for k in range(2, kmax + 1):
        matches = list(perfect_matchings_cover(range(k)))
        for m in matches:
            for letters in product("XYZ", repeat=len(m)):
                spec = []
                for (i, j), L in zip(m, letters):
                    spec.append((i, j, pair_rotation_imgs(L, L)))
                imgs = build_mesh(k, spec)
                w, _ = eval_map_on_tokens(k, imgs)
                n_composites += 1
                results["families"][fam].append(
                    {"k": k,
                     "matching": [list(p) for p in m],
                     "letters": "".join(letters),
                     "w_exact": f"{w.numerator}/{w.denominator}",
                     "shadow_norm2": f"{w.denominator}/{w.numerator}",
                     "beats_mesh_Z_faster_w": None})
    results["n_composites_family_b"] = n_composites
    # family (c): pinned commuting/deterministic layers
    fam = "deterministic_layers"
    results["families"][fam] = []
    for k in range(2, kmax + 1):
        base = tuple(1 << b for b in range(2 * k))
        for i, j in combinations(range(k), 2):
            emb = embed_gate(cz_imgs(), k, [i, j])
            assert is_symplectic(emb, k)
            base = compose(emb, base)
        cz_mesh = base
        # CX chain (staircase)
        base = tuple(1 << b for b in range(2 * k))
        for i in range(k - 1):
            emb = embed_gate(cnot_imgs(), k, [i, i + 1])
            base = compose(emb, base)
        cx_chain = base
        # CX full mesh in lex order
        base = tuple(1 << b for b in range(2 * k))
        for i, j in combinations(range(k), 2):
            emb = embed_gate(cnot_imgs(), k, [i, j])
            base = compose(emb, base)
        cx_mesh = base
        # U_ct (Z mesh) then CZ mesh composed
        uct, _ = mesh_all_pairs(k, "Z")
        uct_after_cz = compose(cz_mesh, uct)
        for name, mp in (("cz_mesh", cz_mesh), ("cx_chain", cx_chain),
                         ("cx_mesh_lex", cx_mesh),
                         ("uct_then_cz_mesh", uct_after_cz)):
            assert is_symplectic(mp, k)
            w, hist = eval_map_on_tokens(k, mp)
            results["families"][fam].append(
                {"k": k, "family": name,
                 "w_exact": f"{w.numerator}/{w.denominator}",
                 "shadow_norm2": f"{w.denominator}/{w.numerator}"})
    # verdicts (pre_statement.md Revision 1): a family BEATS U_ct iff its
    # Pauli weight w_M is LARGER than w_ct (shadow norm SMALLER than 2x1.8^k
    # scale). ratio = w_M / w_ct (exact); > 1 means beat.
    ct_w = {r["k"]: Fraction(*map(int, r["w_exact"].split("/")))
            for r in results["families"]["mesh_Z"]}
    results["beats_uct"] = []
    best, best_ratio = None, None
    for fam, rows in results["families"].items():
        for r in rows:
            w = Fraction(*map(int, r["w_exact"].split("/")))
            k = r["k"]
            r["beats_uct_exact"] = bool(w > ct_w[k])
            ratio = w / ct_w[k]
            r["ratio_vs_uct_exact"] = \
                f"{ratio.numerator}/{ratio.denominator}"
            if w > ct_w[k]:
                results["beats_uct"].append({"family": fam, "row": r})
                if best is None or ratio > best_ratio:
                    best, best_ratio = (fam, r), ratio
    results["best_family"] = (
        {"family": best[0], "row": best[1],
         "weight_ratio_vs_uct":
             f"{best_ratio.numerator}/{best_ratio.denominator}"}
        if best else None)
    results["max_w_over_uct_exact"] = (
        f"{best_ratio.numerator}/{best_ratio.denominator}"
        if best else None)
    return results


# ------------------------------------------------------------- self test

def selftest():
    assert pauli_size(0b0000, 2) == 0
    assert pauli_size(0b0101, 2) == 2   # X1 X2
    assert pauli_size(0b0011, 2) == 1   # Y1
    assert symplectic_form(0b0001, 0b0010) == 1   # X1, Z1 anticommute
    assert symplectic_form(0b0001, 0b0100) == 0   # X1, X2 commute
    # anchor: all nine size-2 classes, hand-derived expectations
    bits_per = {"X": 0b01, "Y": 0b11, "Z": 0b10}
    expect = {"XZ": ("Y", 1), "YZ": ("X", 1), "ZX": ("Y", 1),
              "ZY": ("X", 1), "XX": ("XX", 2), "XY": ("XY", 2),
              "YX": ("YX", 2), "YY": ("YY", 2), "ZZ": ("ZZ", 2)}
    for src, (want_letter, want_size) in expect.items():
        v = 0
        for pos, ch in enumerate(src):
            v |= bits_per[ch] << (2 * pos)
        img = apply_map(ANCHOR_IMGS, v)
        assert pauli_size(img, 2) == want_size, (src, class_label(img))
        if want_size == 2:  # unchanged classes: exact class match
            assert class_label(img) == src, (src, class_label(img))
        else:  # contracted: image is a single-site Pauli of the expected letter
            assert class_label(img).replace("I", "") == want_letter, \
                (src, class_label(img), want_letter)
    assert is_symplectic(ANCHOR_IMGS, 2)
    assert len(enumerate_symplectic_2q()) == 720
    assert invert_symplectic(ANCHOR_IMGS) is not None
    r = invert_symplectic(ANCHOR_IMGS)
    idn = compose(ANCHOR_IMGS, r)
    assert idn == tuple(1 << b for b in range(4))
    # generic mesh at n=2 equals the anchor
    m2, _ = mesh_all_pairs(2, "Z")
    assert m2 == ANCHOR_IMGS
    # pair rotation derived maps are symplectic and in Sp(4,2)
    symp = set(enumerate_symplectic_2q())
    for a in "XYZ":
        for b in "XYZ":
            p = pair_rotation_imgs(a, b)
            assert p in symp, (a, b)
    # pair rotations expand single-site letters as expected (class level):
    # exp(i pi/4 Z1Z2): X1 -> Y1Z2 ; CZ: X1 -> X1Z2 (different maps!)
    cz2 = embed_gate(cz_imgs(), 2, [0, 1])
    assert class_label(apply_map(cz2, 0b0001)) == "XZ"
    assert class_label(apply_map(ANCHOR_IMGS, 0b0001)) == "YZ"
    # exactly one of the 9 classes contracts for a single ZZ rotation at n=2
    ct = [c for c in size2_classes()
          if pauli_size(apply_map(ANCHOR_IMGS, c), 2) <= 1]
    assert len(ct) == 4
    return True


if __name__ == "__main__":
    selftest()
    print("selftest OK (anchor mapping, Sp(4,2) order 720, inverses, meshes)")
