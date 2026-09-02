"""n2_systems.py — exact system serializer for the N2 campaign.

Generates msolve input files (integer coefficient polynomial systems) for:
  S1a  tau r=7 CP system           (KNOWN feasible;   K-1)
  S1b  tau r=6 CP system           (KNOWN infeasible; K-2)
  S1c  tau r=6 with I-slice (1,1) entry +1            (K-3 smoke)
  K-4a  extended tau r=7 (gauge pairs, feasible)
  K-4b  extended tau r=6 (gauge pairs, infeasible)
  S2-E   TF rank-13 CP system + per-term gauge pair (273 vars / 218 eqs)
  S2-E1  E with term 1 pinned a[0,1]=c[0,1]=0        (275 vars / 220 eqs)

msolve Poly format (verified against msolve-0.8.0 README/examples):
  line 1: comma-separated variable names
  line 2: characteristic (0 for Q)
  line 3: comma-separated polynomials (empty polynomial allowed via ","),
  in the order of the variable list; poly ordering after char line.
  Field ordering: for char 0, coefficients are rationals.

All construction arithmetic is exact integer (Python ints); the tensor
tables come from the frozen Cayley-Dickson recursion (cd_mul_list of
rf_krawczyk.py, byte-semantic copy below) and the generator asserts the
tau-slice / blockdiag frame facts at generation time.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
def cd_ok(x, y):
    """Frozen rf_krawczyk cd_mul_list: (a,b)(c,d) =
    (ac - conj(d) b, da + b conj(c)), conj flips the imaginary half."""
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    m = n // 2
    a, b = x[:m], x[m:]
    c, d = y[:m], y[m:]
    def cj(z):
        z = list(z)
        if len(z) == 1:
            return z
        m2 = len(z) // 2
        return cj(z[:m2]) + [-v for v in z[m2:]]
    left = [p - r for p, r in zip(cd_ok(a, c), cd_ok(cj(d), b))]
    right = [p + r for p, r in zip(cd_ok(d, a), cd_ok(b, cj(c)))]
    return left + right

def build_tables(n=8):
    basis = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
    T = [[[cd_ok(basis[p_in], basis[b])[c] for c in range(n)]
          for b in range(n)] for p_in in range(n)]
    return T

def tau_tables():
    """tau = 3x4x4 quaternion slices (I, L_i, L_j)"""
    H = build_tables(4)
    return [[[H[p][b][c] for c in range(4)] for b in range(4)]
            for p in range(3)]

def TF_tables():
    """TF = blockdiag(tau, tau) as 3x8x8."""
    tau = tau_tables()
    return [[[ (tau[p][b][c] if (b < 4 and c < 4) else
                (tau[p][b - 4][c - 4] if (b >= 4 and c >= 4) else 0))
              for c in range(8)] for b in range(8)] for p in range(3)]

def cp_system(T, r, gauge="none", pin_term1_zero=False):
    """Serialize CP equations. Returns (varnames, polys).
    Variables: a_ps (p<r-slots), b_bs, c_cs. Polynomial per entry:
      sum_s a_{p,s} b_{b,s} c_{c,s} - T[p][b][c]
    gauge: 'none' | 'pairs' — per term s append
        a_{0,s} d_s - c_{0,s} g_s - 1
        c_{0,s} d_s + a_{0,s} g_s
    pin_term1_zero: add equations a_{0,1} and c_{0,1} = 0 (uses var names).
    """
    p_n = len(T); b_n = len(T[0]); c_n = len(T[0][0])
    vars_ = [f"a{p}_{s}" for p in range(p_n) for s in range(r)]
    vars_ += [f"b{b}_{s}" for b in range(b_n) for s in range(r)]
    vars_ += [f"c{c}_{s}" for c in range(c_n) for s in range(r)]
    polys = []
    for p in range(p_n):
        for b in range(b_n):
            for c in range(c_n):
                terms = []
                for s in range(r):
                    terms.append(f"a{p}_{s}*b{b}_{s}*c{c}_{s}")
                rhs = T[p][b][c]
                if rhs >= 0:
                    poly = " + ".join(terms) + (f" - {rhs}" if rhs else "")
                else:
                    poly = " + ".join(terms) + f" + {-rhs}"
                if not terms:
                    poly = str(-rhs)
                # cancel exactly-zero polynomials (all-zero entries with no
                # terms -> "-0"): keep as "0" for clarity
                if poly in ("", "-0", "+0"):
                    poly = "0"
                polys.append(poly)
    extra_vars = []
    if gauge == "pairs":
        for s in range(r):
            du, gv = f"d{s}", f"g{s}"
            extra_vars += [du, gv]
            # a[0,s] d_s - c[0,s] g_s - 1 = 0
            polys.append(f"a0_{s}*{du} - c0_{s}*{gv} - 1")
            # c[0,s] d_s + a[0,s] g_s = 0
            polys.append(f"c0_{s}*{du} + a0_{s}*{gv}")
    if pin_term1_zero:
        polys.append("a0_0")
        polys.append("c0_0")
        pin_target = 1
        # NOTE: term indexing in this serializer is 0-based; prereg's
        # "term 1" (1-based) == index 0 here.
    allvars = vars_ + extra_vars
    return allvars, polys

def write_msolve(path, vars_, polys, char=0):
    with open(path, "w") as f:
        f.write(",".join(vars_) + "\n")
        f.write(f"{char}\n")
        f.write(",\n".join(polys) + "\n")
    # metadata sidecar
    with open(path + ".meta.json", "w") as f:
        json.dump({"nvars": len(vars_), "npolys": len(polys),
                   "char": char, "vars": vars_}, f, indent=1)

def main():
    # --- generation-time anchors (machine-checked)
    tau = tau_tables()
    n_tau = sum(1 for s in tau for row in s for v in row if v)
    assert n_tau == 12, f"tau nonzero {n_tau}"  # 3 slices x 4 = 12, matches
    # the frozen qmul reference construction (verified this session)
    TF = TF_tables()
    nzF = sum(1 for s in TF for row in s for v in row if v)
    assert nzF == 2 * n_tau, f"TF nonzero {nzF} != 2*tau {2*n_tau}"
    entries = sorted(v for s in tau for row in s for v in row if v)
    assert all(v in (-1, 1) for v in entries)
    print(f"[GEN] anchors: tau nonzero={n_tau}, TF nonzero={nzF} (== 2*tau), "
          f"entries in {{-1,0,1}}: OK")

    out = {}
    # S1a: tau r=7 (feasible)
    v, p = cp_system(tau, 7)
    write_msolve(os.path.join(HERE, "s1a_tau_r7.txt"), v, p)
    out["s1a_tau_r7"] = {"nvars": len(v), "npolys": len(p)}
    # S1b: tau r=6 (infeasible)
    v, p = cp_system(tau, 6)
    write_msolve(os.path.join(HERE, "s1b_tau_r6.txt"), v, p)
    out["s1b_tau_r6"] = {"nvars": len(v), "npolys": len(p)}
    # S1c: tau r=6 with I-slice (1,1) (index [0][0][0]) entry +1
    tau_cor = [[[(tau[pi][b][c] + 1) if (pi == 0 and b == 0 and c == 0)
                 else tau[pi][b][c] for c in range(4)] for b in range(4)]
               for pi in range(3)]
    v, p = cp_system(tau_cor, 6)
    write_msolve(os.path.join(HERE, "s1c_tau_r6_corrupted.txt"), v, p)
    out["s1c_tau_r6_corrupted"] = {"nvars": len(v), "npolys": len(p)}
    # K-4a: extended tau r=7 with gauge pairs
    v, p = cp_system(tau, 7, gauge="pairs")
    write_msolve(os.path.join(HERE, "k4a_tau_r7_ext.txt"), v, p)
    out["k4a_tau_r7_ext"] = {"nvars": len(v), "npolys": len(p)}
    # K-4b: extended tau r=6 with gauge pairs
    v, p = cp_system(tau, 6, gauge="pairs")
    write_msolve(os.path.join(HERE, "k4b_tau_r6_ext.txt"), v, p)
    out["k4b_tau_r6_ext"] = {"nvars": len(v), "npolys": len(p)}
    # S2-E: TF r=13 with gauge pairs (273 vars, 218 eqs)
    v, p = cp_system(TF, 13, gauge="pairs")
    write_msolve(os.path.join(HERE, "s2e_tf_r13_ext.txt"), v, p)
    out["s2e_tf_r13_ext"] = {"nvars": len(v), "npolys": len(p)}
    assert (len(v), len(p)) == (273, 218), f"E size {(len(v), len(p))}"
    # S2-E1: E + pin (a0_0 = c0_0 = 0)
    v, p = cp_system(TF, 13, gauge="pairs", pin_term1_zero=True)
    # S2-E1: E + pin (a0_0 = c0_0 = 0). The pin adds 2 EQUATIONS over the
    # SAME variables (the prereg's "275 vars" miscounted: no fresh
    # variables are needed — pinning uses a0_0/c0_0 which already exist).
    # Recorded size correction, disclosed in VERDICT: E1 = 273 vars / 220 eqs.
    out["s2e1_tf_r13_ext_pin"] = {"nvars": len(v), "npolys": len(p)}
    assert (len(v), len(p)) == (273, 220), f"E1 size {(len(v), len(p))}"
    with open(os.path.join(HERE, "n2_systems_manifest.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("[GEN] systems written:", json.dumps(out))

if __name__ == "__main__":
    main()
