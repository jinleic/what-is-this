"""Gate A driver — see module docstring of gate_a.py for the check list."""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

from flint import fmpz

from tensor_data import (
    C_ALIASES,
    EXPANDED_PRODUCTS_TEXT,
    LEFT_SLP,
    OUTPUT_SLP,
    PRODUCT_PAIRS,
    RIGHT_SLP,
    U_BLOCK_PRINTED,
    V_BLOCK_PRINTED,
    W_BLOCK_PRINTED,
)

HERE = Path(__file__).resolve().parent
R = 23


# --------------------------------------------------------------------------
# SLP expansion
# --------------------------------------------------------------------------
def expand_slp_terms(slp, allowed_atoms):
    """Return (env, gate_count). env[name] = {atom: coef} over allowed atoms.
    Each gate must have exactly 2 terms (cost-1 addition/subtraction)."""
    env = {}
    for name, terms in slp:
        assert len(terms) == 2, f"gate {name}: expected a 2-term addition gate"
        acc: dict[str, int] = {}
        for sign, atom in terms:
            src = env.get(atom)
            if src is None and atom.startswith(allowed_atoms):
                acc[atom] = acc.get(atom, 0) + sign
            elif src is not None:
                for k, v in src.items():
                    acc[k] = acc.get(k, 0) + sign * v
            else:
                raise AssertionError(f"gate {name}: atom {atom} not available")
        env[name] = acc
    return env, len(slp)


def slp_factors():
    """Expand the paper's three networks into factor tensors (9x23 blocks)."""
    left_env, n_left = expand_slp_terms(LEFT_SLP, "A")
    right_env, n_right = expand_slp_terms(RIGHT_SLP, "B")

    for lu, rv in PRODUCT_PAIRS:
        assert lu in left_env or lu.startswith("A"), f"left operand {lu} undefined"
        assert rv in right_env or rv.startswith("B"), f"right operand {rv} undefined"

    U = [[0] * R for _ in range(9)]
    V = [[0] * R for _ in range(9)]
    for r, (lu, rv) in enumerate(PRODUCT_PAIRS):
        lterms = {lu: 1} if lu.startswith("A") else left_env[lu]
        rterms = {rv: 1} if rv.startswith("B") else right_env[rv]
        for a, c in lterms.items():
            U[int(a[1:])][r] += c
        for b, c in rterms.items():
            V[int(b[1:])][r] += c
        # ternary check
        for x in list(lterms.values()) + list(rterms.values()):
            assert x in (-1, 0, 1)

    # Output network over symbols M_r
    out_env: dict[str, dict[int, int]] = {}
    for name, terms in OUTPUT_SLP:
        assert len(terms) == 2, f"output gate {name}: expected 2-term gate"
        acc: dict[int, int] = {}
        for sign, atom in terms:
            if atom.startswith("M"):
                idx = int(atom[1:])
                acc[idx] = acc.get(idx, 0) + sign
            else:
                for k, v in out_env[atom].items():
                    acc[k] = acc.get(k, 0) + sign * v
        out_env[name] = acc

    W = [[0] * R for _ in range(9)]
    for k, wname in enumerate(C_ALIASES):
        for midx, coef in out_env[wname].items():
            W[k][midx] += coef

    return (U, V, W), (n_left, n_right, len(OUTPUT_SLP))


def to_fmpz_block(block):
    return [[fmpz(x) for x in row] for row in block]


def brent_check(U, V, W):
    """All 729 Brent identities over Z (exact fmpz). Returns (failures, kinds)."""
    failures = []
    kinds = {"nonzero_pass": 0, "zero_pass": 0, "nonzero_fail": 0, "zero_fail": 0}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # pull the row once
                Ur = U[3 * i + k]
                for ip in range(3):
                    for jp in range(3):
                        Wr = W[3 * ip + jp]
                        for kp in range(3):
                            Vr = V[3 * kp + j]
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = Ur[r], Vr[r], Wr[r]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            expect = 1 if (i == ip and j == jp and k == kp) else 0
                            sv = int(s)
                            if sv == expect:
                                kinds["nonzero_pass" if expect else "zero_pass"] += 1
                            else:
                                failures.append((i, j, k, ip, jp, kp, sv, expect))
                                kinds["nonzero_fail" if expect else "zero_fail"] += 1
    return failures, kinds


def check_printed_blocks(U, V, W):
    u_ok = all(U[i][r] == U_BLOCK_PRINTED[i][r] for i in range(9) for r in range(R))
    v_ok = all(V[i][r] == V_BLOCK_PRINTED[i][r] for i in range(9) for r in range(R))
    w_ok = all(W[i][r] == W_BLOCK_PRINTED[i][r] for i in range(9) for r in range(R))
    return u_ok, v_ok, w_ok


def check_termary(U, V, W):
    return all(
        x in (-1, 0, 1)
        for block in (U, V, W)
        for row in block
        for x in row
    )


def check_expanded_products(U, V, W):
    """Section 5.1: p_r+1 = M_r must match factor columns; c_ij must match W."""
    import re
    amap = {f"a{i}{j}": 3 * (i - 1) + (j - 1) for i in (1, 2, 3) for j in (1, 2, 3)}
    bmap = {f"b{i}{j}": 3 * (i - 1) + (j - 1) for i in (1, 2, 3) for j in (1, 2, 3)}

    def parse_side(expr, vmap):
        terms = re.findall(r"([+-]?)\s*(\d*)\s*\*?\s*([ab]\d\d)", expr)
        coef = [0] * 9
        for sign, num, var in terms:
            c = int(num) if num else 1
            if sign == "-":
                c = -c
            coef[vmap[var]] += c
        return coef

    lines = [ln.strip() for ln in EXPANDED_PRODUCTS_TEXT.strip().splitlines() if ln.strip()]
    p_lines = [ln for ln in lines if ln.startswith("p")]
    c_lines = [ln for ln in lines if ln.startswith("c")]
    assert len(p_lines) == R and len(c_lines) == 9
    mismatches = []
    for r, ln in enumerate(p_lines):
        lhs, rhs = ln.split("=", 1)
        left_expr, right_expr = rhs.split("*")
        lcoef = parse_side(left_expr, amap)
        rcoef = parse_side(right_expr, bmap)
        for i in range(9):
            if lcoef[i] != U[i][r]:
                mismatches.append(f"p{r+1:02d} left A{i}: {lcoef[i]} vs {U[i][r]}")
            if rcoef[i] != V[i][r]:
                mismatches.append(f"p{r+1:02d} right B{i}: {rcoef[i]} vs {V[i][r]}")
    cmap = {}
    for ln in c_lines:
        lhs, rhs = ln.split("=", 1)
        terms = re.findall(r"([+-]?)\s*p(\d\d)", rhs)
        acc = {}
        for sign, num in terms:
            idx = int(num) - 1
            s = -1 if sign == "-" else 1
            acc[idx] = acc.get(idx, 0) + s
        i, j = int(lhs[1]), int(lhs[2])
        cmap[3 * (i - 1) + (j - 1)] = acc
    for k in range(9):
        acc = cmap[k]
        for r in range(R):
            expect = acc.get(r, 0)
            if W[k][r] != expect:
                mismatches.append(f"c-row {k} product {r}: W={W[k][r]} vs expanded {expect}")
    return mismatches


def check_perminov(json_path, U, V, W):
    """Decode Perminov's cr58_cn122 JSON using HIS conventions, independently
    re-implemented from his public loader source (Scheme.from_reduced,
    commit 98ba522):
      - THREE SEPARATE 0-based namespaces: u-space (A-entries 0..8 + fresh
        u-vars from id 9), v-space (B-entries 0..8 + fresh v-vars from id 9),
        w-space (product slots 0..22 + fresh w-vars from id 23) — directly
        from __parse_reduced_vars(fresh_vars, real_variables);
        are COLUMN-major (3x3); the paper's involution CPERM maps Perminov's
        C indexing to the paper's row-major indexing.
      - data["u"][r], data["v"][r]: expansion of product r into (A|B)-entries.
    Returns (mismatch_list, complexity_dict)."""
    data = json.loads(Path(json_path).read_text())
    assert data["n"] == [3, 3, 3] and data["m"] == 23

    CPERM = [0, 3, 6, 1, 4, 7, 2, 5, 8]  # paper Section 3.2, an involution

    nu, nv, nw = (len(data[k]) for k in ("u_fresh", "v_fresh", "w_fresh"))
    # Perminov's __parse_reduced_vars: fresh variable ids are 0-based
    # (real_variables + i) in THREE SEPARATE namespaces:
    #   u-space: real vars 0..8 = A-entries,   fresh ids 9..9+nu-1
    #   v-space: real vars 0..8 = B-entries,   fresh ids 9..9+nv-1
    #   w-space: real vars 0..22 = products,   fresh ids 23..23+nw-1
    known_u = {9 + i: data["u_fresh"][i] for i in range(nu)}
    known_v = {9 + i: data["v_fresh"][i] for i in range(nv)}
    known_w = {23 + i: data["w_fresh"][i] for i in range(nw)}

    def expand(expr, known, nreal, side_of):
        """Resolve ids: < nreal -> (side_of, direct id); in `known` -> recurse."""
        out: list[tuple[int, int, int]] = []
        for e in expr:
            ix, val = e["index"], e["value"]
            if ix in known:
                for s2, i2, v2 in expand(known[ix], known, nreal, side_of):
                    out.append((s2, i2, v2 * val))
            elif 0 <= ix < nreal:
                out.append((side_of, ix, val))
            else:
                raise AssertionError(f"unknown variable index {ix} in {side_of}-space")
        return out

    Un = [[0] * R for _ in range(9)]
    Vn = [[0] * R for _ in range(9)]
    Wn = [[0] * R for _ in range(9)]
    for r, expr in enumerate(data["u"]):
        for s, i, c in expand(expr, known_u, 9, 0):
            Un[i][r] += c
    for r, expr in enumerate(data["v"]):
        for s, i, c in expand(expr, known_v, 9, 1):
            Vn[i][r] += c
    for kcol, expr in enumerate(data["w"]):
        k = CPERM[kcol]
        for s, r, c in expand(expr, known_w, 23, 2):
            Wn[k][r] += c

    tern = all(x in (-1, 0, 1) for blk in (Un, Vn, Wn) for row in blk for x in row)
    assert tern, "decoded tensor not ternary - decode error"

    mism = []
    for i in range(9):
        for r in range(R):
            if Un[i][r] != U[i][r]:
                mism.append(("U", i, r, Un[i][r], U[i][r]))
            if Vn[i][r] != V[i][r]:
                mism.append(("V", i, r, Vn[i][r], V[i][r]))
            if Wn[i][r] != W[i][r]:
                mism.append(("W", i, r, Wn[i][r], W[i][r]))
    return mism, data["complexity"]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--perminov-json", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = {}

    # A1: SLP expansion + counts
    (U, V, W), counts = slp_factors()
    n_left, n_right, n_out = counts
    results["slp_gate_counts"] = {"left": n_left, "right": n_right, "output": n_out,
                                  "total": n_left + n_right + n_out}
    results["counts_pass"] = (n_left, n_right, n_out) == (13, 14, 28)

    Uf, Vf, Wf = to_fmpz_block(U), to_fmpz_block(V), to_fmpz_block(W)

    # A4: Brent (exact)
    failures, kinds = brent_check(Uf, Vf, Wf)
    results["brent_kinds"] = kinds
    results["brent_failures"] = failures[:20]
    results["brent_all_pass"] = (len(failures) == 0 and kinds["nonzero_pass"] == 27
                                 and kinds["zero_pass"] == 702)

    # A2: printed blocks
    u_ok, v_ok, w_ok = check_printed_blocks(U, V, W)
    results["printed_blocks_match"] = {"U": u_ok, "V": v_ok, "W": w_ok}

    # ternary alphabet
    results["ternary"] = check_termary(U, V, W)

    # A3: expanded products / outputs (Section 5.1)
    mism_exp = check_expanded_products(U, V, W)
    results["expanded_products_mismatches"] = mism_exp[:20]
    results["expanded_products_pass"] = len(mism_exp) == 0
    if args.perminov_json:
        try:
            mism, complexity = check_perminov(args.perminov_json, U, V, W)
            results["perminov_mismatches"] = mism[:20]
            results["perminov_pass"] = len(mism) == 0
        except Exception as e:  # noqa: BLE001
            results["perminov_error"] = repr(e)

    verdict = (
        results["counts_pass"]
        and results["brent_all_pass"]
        and all(results["printed_blocks_match"].values())
        and results["ternary"]
        and results["expanded_products_pass"]
        and results.get("perminov_pass", False)
    )
    results["VERDICT_GATE_A"] = "PASS" if verdict else "FAIL"

    out = json.dumps(results, indent=2, default=int)
    print(out)
    if args.out:
        Path(args.out).write_text(out)
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
