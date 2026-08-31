"""Gate A census: waterfall measurement on binary Goppa instances.

Per instance (m, n, t, seed) everything below is EXACT finite-field
arithmetic over E = F_{2^m} (MUL/INV table engine, cross-validated vs
python-flint):

  QUANTITIES (per pre_statement.md):
  N_eval(c)   = dim{A in E[Z]_{<=2t+2} : A(a)=0 for a in S_c}
                = 2t+3 - rank(Vandermonde_c)   (family evaluation system)
  N_flag(c)   = nullity of the FULL derivative-flag Step-3 system
                Sol(S_c) = {H in E[Z]^k_{<=D} : (d^[j]H)(a) in V_a^[j],
                                                 a in S_c, 0 <= j < r_a}
                measured by incremental projected echelon over the ambient
                coefficient space of dimension k(D+1).  This is Apon's
                Table-1 protocol ("perfect oracle") verbatim.
  N_fam(c)    = nullity of the flag system RESTRICTED to the (2t+3)-
                dimensional family T(A)  (Apon's Lemma-7 claim: equals
                N_eval(c); a contradiction here is the refutation event)
  rank(a,j)   = rank of the (a, j) flag block restricted to the family
                (predicted <= 1, = 1 iff j even and the (j+1)-st F-jet at
                a extends V_a^[j]; j odd rows vanish by the (j+1) = 0
                parity in (j+1) d^[j+1]F)
  jetprof(a)  = dims (r_0 <= ... <= r_R) of V_a^[j] (the F-jet flag at a)
  degenerate(a) = [ exists j < R with r_{j+1} = r_j ]  (Remark-8 failure)

  FLAG BLOCK STRUCTURE: order-j conditions at point a on ambient H are
  H_(a,j) (x) w_j(a)^T   (Kronecker; H_(a,j) = complement basis of
  V_a^[j]; w_j(d) = binom(d,j) a^{d-j}, Lucas parity exact).

  CERTIFICATES (frozen per instance):
   - G, support, seeds, guard outcomes (incl. Delta_{p,q} pair)
   - jet profiles per held point
   - N_eval/N_flag/N_fam for c = 1..2t+4
   - echelon snapshot (first rows) of the full system at c = 2t+3
   - lower-bound certificate: A_S = prod_{a in S_{2t+2}} (Z-a): machine
     check that T(A_S) jets lie in the flags at every held point (depth R)
"""
from __future__ import annotations

import hashlib
import json
import random

import numpy as np

from fastfield import EField
from gfield import GF, _factor
from instance import Instance


# --------------------------------------------------------------------------- #
# instance build (vectorized; reuses guard code from instance.py at small n)
# --------------------------------------------------------------------------- #

def build_instance(m, n, t, seed):
    inst = Instance(m, n, t, seed)
    return inst


# --------------------------------------------------------------------------- #
# jet machinery
# --------------------------------------------------------------------------- #

def lucas_w(ef, D, j, a):
    """w[d] = binom(d, j) * a^(d-j) for d = 0..D  (0 when j not <= d or
    binom even)."""
    d = np.arange(D + 1, dtype=np.int64)
    par = ((d & j) == j) & (d >= j)          # binom(d,j) odd  (Lucas)
    w = np.zeros(D + 1, dtype=np.uint16)
    idx = np.nonzero(par)[0]
    if j == 0:
        # a^(d-0) = a^d ; zero point handled
        ap = ef.powers(a, D)
        w[idx] = ap[idx]
        return w
    # a != 0: a^(d-j) with d-j in [0, D-j]; requires D >= j
    if a != 0 and D >= j:
        ap = ef.powers(a, D - j)
        w[idx] = ap[idx - j]
    else:
        # a == 0 (or D < j): only d = j term (a^0 = 1) survives
        sel = idx[idx == j] if a == 0 else idx[:0]
        w[sel] = 1
    return w


def poly_jets(ef, coeffs, a, order):
    """jets 0..order of poly (coeff array) at a -> list of scalars."""
    cf = np.asarray(coeffs, dtype=np.uint16)
    out = []
    for j in range(order + 1):
        w = lucas_w(ef, len(cf) - 1, j, a)
        terms = ef.MUL[cf, w]
        out.append(int(np.bitwise_xor.reduce(terms)) if terms.size else 0)
    return out


# --------------------------------------------------------------------------- #
# Apon Def-4 map pieces (per instance, vectorized where it matters)
# --------------------------------------------------------------------------- #

def compute_R_G(inst, ef):
    """R_G = (Pi' * Pi^{-1}) mod G^2, deg < 2t."""
    gf = inst.gf
    G2 = gf.pmul(inst.G, inst.G)              # deg 2t, monic
    # Pi mod G^2
    Pi_mod = gf.pmod(inst.Pi, G2)
    # inverse of Pi_mod mod G^2 (Pi invertible mod G2 since gcd(Pi, G) = 1)
    Inv = gf.pinv_mod(Pi_mod, G2)
    RG = gf.pmod(gf.pmul(gf.ptrim(inst.PiD), Inv), G2)
    return gf.ptrim(RG)


def B_of(ef, gf, RG, G2, A_coeffs):
    """B_A = (A * R_G) mod G^2."""
    return gf.pmod(gf.pmul(np.asarray(A_coeffs, dtype=np.uint16).tolist(),
                           RG), G2)


# --------------------------------------------------------------------------- #
# per-instance census
# --------------------------------------------------------------------------- #

def jets_F_at_points(inst, ef, points, order):
    """F-jets at given support points: returns dict a -> (k, order+1) uint16.

    F's coordinate f satisfies f(a_i) = Y[j,i]/lam_i.  Jets are linear in
    the value vector v: (d^[e] f)(a) = sum_i v_i (d^[e] L_i)(a), where
    L_i = Pi/(Z-a_i)/Pi'(a_i).

    (d^[e] L_i)(a):
      a == a_i : d^[e+1] Pi(a_i) / Pi'(a_i)           [Pi has simple root]
      a != a_i : jet of Pi/(Z-a_i) at a divided by Pi'(a_i):
                 Pi(a+u) = (a-a_i+u) * Pi_i(a+u)  =>  series division.
    """
    gf = inst.gf
    k, n, D = inst.k, inst.n, inst.D
    sup = inst.support
    sup_of = {a: i for i, a in enumerate(sup)}
    out = {}
    # Pi jets at each held point once
    Pi_cf = np.array(inst.Pi, dtype=np.uint16)
    PiDa = {a: int(ef.MUL[Pi_cf, lucas_w(ef, len(Pi_cf)-1, 1, a)].sum() % ef.q) for a in []}
    for a in points:
        PiJ = poly_jets(ef, Pi_cf, a, order + 1)   # need order+1 for L i=a case
        Lj = np.zeros((n, order + 1), dtype=np.uint16)
        for i, ai in enumerate(sup):
            if ai == a:
                # d^[e] L_i(a) = PiJ[e+1] / Pi'(a_i)
                pii = int(0)
                w1 = lucas_w(ef, len(Pi_cf) - 1, 1, ai)
                pii = int(np.bitwise_xor.reduce(ef.MUL[Pi_cf, w1]))
                for e in range(order + 1):
                    Lj[i, e] = ef.MUL[PiJ[e + 1], ef.INV[pii]]
            else:
                d = a ^ ai   # a - a_i (char 2)
                # series division: Pi_i(a+u) = Pi(a+u) / (d + u), both known
                # jets to `order`; divisor unit series [d, 1, 0, ...]
                num = np.array(PiJ[: order + 1], dtype=np.uint16)
                den = np.zeros(order + 1, dtype=np.uint16)
                den[0] = d
                den[1] = 1
                # solve (den * q = num): q_0 = num_0/d; q_e = num_e + den_1*q_{e-1}
                q = np.zeros(order + 1, dtype=np.uint16)
                acc = num[0]
                q[0] = ef.MUL[acc, ef.INV[d]]
                for e in range(1, order + 1):
                    # num_e = d*q_e + q_{e-1}  => q_e = (num_e + q_{e-1})/d
                    lhs = num[e] ^ q[e - 1]
                    q[e] = ef.MUL[lhs, ef.INV[d]]
                pii = int(np.bitwise_xor.reduce(ef.MUL[Pi_cf, lucas_w(ef, len(Pi_cf) - 1, 1, ai)]))
                for e in range(order + 1):
                    Lj[i, e] = ef.MUL[q[e], ef.INV[pii]]
        # F jets: (k, order+1) = V(binary) @ (Lj * inv(lam_i))
        Vm = inst.Y.astype(np.uint16)                     # k x n binary
        Linv = np.zeros((n, order + 1), dtype=np.uint16)
        for i in range(n):
            Linv[i] = ef.MUL[ef.INV[inst.lam[i]], Lj[i]]
        # (k,n) x (n,R+1): FJ[k,e] = XOR_i MUL[V[k,i], Linv[i,e]]
        FJ = np.zeros((k, order + 1), dtype=np.uint16)
        for e in range(order + 1):
            col = ef.MUL[Vm, Linv[:, e][None, :]]        # (k, n)
            FJ[:, e] = np.bitwise_xor.reduce(col, axis=1)
        out[a] = FJ
    return out


def flag_data(inst, ef, points, R):
    """per held point: F-jet flag bases, complements H_(a,j), rank profile."""
    FJ = jets_F_at_points(inst, ef, points, R)
    data = {}
    for a in points:
        M = FJ[a]                       # (k, R+1)
        prof = []
        bases = []                      # flag row bases V^[j] (as row stacks)
        cur = np.zeros((0, inst.k), dtype=np.uint16)
        for j in range(R + 1):
            cur = np.vstack([cur, M[:, j][None, :]])
            red, piv = ef.rref(cur)
            prof.append(len(piv))
            bases.append(red[: len(piv)].copy())
        comps = []
        for j in range(R + 1):
            comps.append(complement(ef, bases[j], inst.k))
        data[a] = {"FJ": M, "prof": prof, "flagbase": bases, "H": comps}
    return data


def complement(ef, W, k):
    """rows spanning {u : u·w = 0 for all rows w of W}, over E."""
    if W.shape[0] == 0:
        return np.eye(k, dtype=np.uint16)
    K = ef.row_nullspace(W)
    if not K:
        return np.zeros((0, k), dtype=np.uint16)
    return np.array(K, dtype=np.uint16)


# --- family restricted system ----------------------------------------------

def family_restricted_system(inst, ef, fd, R):
    """rows: (a, j<R, b) -> h_b . jets of T(Z^m)(a).  Built EXPLICITLY
    (naive series route). Returns dict (a, j) -> matrix rows x (2t+3)."""
    gf = inst.gf
    t2 = 2 * inst.t + 2
    cols = t2 + 1
    RG = compute_R_G(inst, ef)
    G2 = gf.pmul(inst.G, inst.G)
    G2c = np.array(gf.ptrim(G2), dtype=np.uint16)
    RGc = np.array(RG, dtype=np.uint16)
    # family monomials Z^m, m = 0..2t+2
    blocks = {}
    points = list(fd.keys())
    for a in points:
        d = fd[a]
        # G series at a to order R+1 ; invG2 series to order R
        GJ = poly_jets(ef, np.array(inst.G, dtype=np.uint16), a, R + 1)
        Ga = GJ[0]
        assert Ga != 0
        G2J = conv(ef, GJ, GJ, R + 1)
        invG2 = series_inv_jets(ef, G2J, R)
        # B_m and jet rows for each m
        # jets of B_m at a: (2t+3, R+1) scalars; jets of Z^m at a: scalars
        ZJ = np.zeros((cols, R + 2), dtype=np.uint16)   # jets of Z^m
        BJ = np.zeros((cols, R + 1), dtype=np.uint16)   # jets of B_m
        rho = np.zeros(cols, dtype=np.uint16)           # rho_A
        for mm in range(cols):
            ZJ[mm] = jets_monomial(ef, mm, a, R + 1)
            if mm <= 2 * inst.t + 2:
                Am = np.zeros(mm + 1, dtype=np.uint16)
                Am[mm] = 1
                Bm = gf.pmod(gf.pmul(Am.tolist(), RG), G2)
                BJ[mm] = poly_jets_vec(ef, np.array(gf.ptrim(Bm), dtype=np.uint16), a, R)
            if mm == t2 and (inst.D % 2 == 0):
                rho[mm] = ef.MUL[1, ef.INV[inst.G[inst.t]]] if False else rho_of(inst, ef)
        # F' jets: (k, R+1) = (q+1 mod 2) * FJ_{q+1}
        FJp = np.zeros((inst.k, R + 1), dtype=np.uint16)
        for q in range(R):
            FJp[:, q] = ef.MUL[(q + 1) % 2, d["FJ"][:, q + 1]]
        # per order j: Psi (cols x k) then restrict with H
        for j in range(R):
            Psi = np.zeros((cols, inst.k), dtype=np.uint16)
            for e in range(j + 1):
                for p in range(e + 1):
                    q = e - p
                    Psi ^= ef.MUL[ZJ[:, p][:, None], FJp[:, q][None, :]]
                    Psi ^= ef.MUL[BJ[:, p][:, None], d["FJ"][:, q][None, :]]
            # multiply by invG2 series: jet_j(Q) = sum_{e+q'=j} invG2_e * numjet_q'
            # (we accumulated numjets per e already into Psi? NO - restructure)
            Psi = np.zeros((cols, inst.k), dtype=np.uint16)
            numjets = []
            for e in range(j + 1):
                N = np.zeros((cols, inst.k), dtype=np.uint16)
                for p in range(e + 1):
                    q = e - p
                    N ^= ef.MUL[ZJ[:, p][:, None], FJp[:, q][None, :]]
                    N ^= ef.MUL[BJ[:, p][:, None], d["FJ"][:, q][None, :]]
                numjets.append(N)
            for e in range(j + 1):
                qp = j - e
                Psi ^= ef.MUL[invG2[e], numjets[qp]]
            # rho_m * Z F term: jets: d^[j](Z F) = a * FJ_j + FJ_{j-1}
            if rho.any():
                ZFj = ef.MUL[a, d["FJ"][:, j]]
                if j >= 1:
                    ZFj = ZFj ^ d["FJ"][:, j - 1]
                Psi ^= ef.MUL[rho, ZFj[None, :]]
            # restrict: rows h . Psi
            H = d["H"][j]
            if H.shape[0] == 0:
                blocks[(a, j)] = np.zeros((0, cols), dtype=np.uint16)
                continue
            rows = np.zeros((H.shape[0], cols), dtype=np.uint16)
            for b in range(H.shape[0]):
                # rows[b, m] = sum_cc MUL[H[b, cc], Psi[m, cc]]
                prod = ef.MUL[H[b][None, :], Psi]        # (cols, k)
                rows[b] = np.bitwise_xor.reduce(prod, axis=1)
            blocks[(a, j)] = rows
    return blocks


def jets_monomial(ef, mm, a, order):
    out = np.zeros(order + 1, dtype=np.uint16)
    d = np.arange(0, mm + 1)
    for e in range(order + 1):
        # jet_e of Z^mm at a = binom(mm, e) a^{mm-e}
        if e > mm:
            out[e] = 0
        elif ((mm & e) == e) and e <= mm:
            if mm - e == 0:
                out[e] = 1
            elif a == 0:
                out[e] = 0
            else:
                ap = ef.powers(a, mm - e)
                out[e] = ap[mm - e]
    return out


def poly_jets_vec(ef, cf, a, order):
    if len(cf) == 0:
        return np.zeros(order + 1, dtype=np.uint16)
    out = np.zeros(order + 1, dtype=np.uint16)
    for j in range(order + 1):
        w = lucas_w(ef, len(cf) - 1, j, a)
        out[j] = np.bitwise_xor.reduce(ef.MUL[cf, w]) if cf.size else 0
    return out


def conv(ef, u, v, order):
    """jet-coefficient convolution: out[e] = sum_{i+j=e} u_i v_j (capped)."""
    u = np.asarray(u, dtype=np.uint16)
    v = np.asarray(v, dtype=np.uint16)
    out = np.zeros(order + 1, dtype=np.uint16)
    for e in range(order + 1):
        i0 = max(0, e - len(v) + 1)
        i1 = min(e, len(u) - 1)
        if i1 < i0:
            continue
        uu = np.asarray(u[i0:i1 + 1], dtype=np.uint16)
        vv = np.asarray(v[e - i1:e - i0 + 1][::-1], dtype=np.uint16)
        out[e] = np.bitwise_xor.reduce(ef.MUL[uu, vv])
    return out


def series_inv_jets(ef, s, order):
    out = np.zeros(order + 1, dtype=np.uint16)
    out[0] = ef.INV[s[0]]
    for e in range(1, order + 1):
        i_idx = np.arange(1, e + 1)
        i_hi = min(e, len(s) - 1)
        i_idx = np.arange(1, i_hi + 1)
        terms = ef.MUL[s[i_idx], out[e - i_idx]]
        acc = np.bitwise_xor.reduce(terms) if terms.size else 0
        out[e] = ef.MUL[acc, ef.INV[s[0]]]
    return out


def rho_of(inst, ef):
    """rho_{Z^{2t+2}} = 1 * (D mod 2) * gamma^{-2}; (D mod 2) handled by
    caller."""
    gamma = inst.G[inst.t]
    return ef.MUL[1, ef.INV[ef.MUL[gamma, gamma]]]


# --- full ambient system (small t) ------------------------------------------

def full_flag_blocks(inst, ef, fd, points, R):
    """Kronecker blocks H_(a,j) (x) w_j(a)^T over the ambient coefficient
    space; returns list in point-major order with metadata."""
    blocks = []
    k, D = inst.k, inst.D
    for a in points:
        for j in range(R):
            H = fd[a]["H"][j]
            if H.shape[0] == 0:
                continue
            w = lucas_w(ef, D, j, a)
            # block[(b, cc), (cc', d)] = H[b, cc'] * w[d] * [cc == cc']
            # implemented as: for each cc: rows H[:, cc] outer w  placed at
            # coordinate block cc
            nb = H.shape[0]
            B = np.zeros((nb * k, k * (D + 1)), dtype=np.uint16)
            for cc in range(k):
                seg = H[:, cc][None, :] & np.ones((1,), dtype=np.uint16)
                B[cc * nb:(cc + 1) * nb, cc * (D + 1):(cc + 1) * (D + 1)] = \
                    ef.MUL[H[:, cc][:, None], w[None, :]]
            blocks.append((a, j, B))
    return blocks


def incremental_sol_dims(inst, ef, blocks, c_values):
    """kernel dimension of the stacked system after the blocks of the first
    c points; blocks are (a, j, B) point-major. Returns dict c -> dim."""
    kD = inst.k * (inst.D + 1)
    dims = {}
    # maintain kernel basis K (rows) of the stack so far
    K = None   # rows of kernel of ALL blocks so far, in coefficient space
    nxt = 0
    pts = sorted(set(a for a, j, B in blocks), key=inst.support.index)
    per_point = {}
    for a, j, B in blocks:
        per_point.setdefault(a, []).append(B)
    K = np.eye(kD, dtype=np.uint16)  # trivial kernel of empty system
    ci = 0
    for c in c_values:
        # add points up to c
        while ci < min(c, len(pts)):
            a = pts[ci]
            Bs = per_point.get(a, [])
            if Bs:
                Bst = np.vstack(Bs)
                K = kernel_with_fixed_kernel(ef, K, Bst)
            ci += 1
        dims[c] = K.shape[0]
    return dims


def kernel_with_fixed_kernel(ef, K, Bnew):
    """kernel of the stacked system {x : B x = 0 (all blocks so far)} with
    K the CURRENT kernel basis rows.  Bnew's rows are expressed in
    K-coordinates via the bilinear pairing, then row-eliminated.
    Returns the new kernel basis rows (ambient coordinates).
    """
    d, N = K.shape
    rb = Bnew.shape[0]
    # C[b, s] = <Bnew[b], K[s]> = XOR_i MUL[Bnew[b,i], K[s,i]]
    Bn = np.asarray(Bnew, dtype=np.uint16)
    # (rb, d, N) is too big; do it row by row
    C = np.zeros((rb, d), dtype=np.uint16)
    for b in range(rb):
        prod = ef.MUL[Bn[b][None, :], K]          # (d, N)
        C[b] = np.bitwise_xor.reduce(prod, axis=1)
    # constraints in K-coordinates: for x = sum_s u_s K[s],
    # (Bnew[b])(x) = sum_s u_s C[b, s]  = (C u)_b.  Kernel condition: C u = 0.
    Kc_list = ef.row_nullspace(C)
    Kc = np.array(Kc_list, dtype=np.uint16) if Kc_list else np.zeros((0, d), dtype=np.uint16)
    out = np.zeros((Kc.shape[0], N), dtype=np.uint16)
    for s in range(Kc.shape[0]):
        acc = np.zeros(N, dtype=np.uint16)
        nz = np.nonzero(Kc[s])[0]
        for i in nz:
            acc ^= ef.MUL[int(Kc[s, i]), K[i]]
        out[s] = acc
    return out


# --------------------------------------------------------------------------- #
# main census run
# --------------------------------------------------------------------------- #

def census(m, n, t, seed, R=4, full_sol=True, record_dir=None):
    inst = build_instance(m, n, t, seed)
    ef = EField(m)
    t2p3 = 2 * t + 3
    cmax = min(2 * t + 4, n)
    S = inst.support[:cmax]

    fd = flag_data(inst, ef, S, R)

    guards = inst.guards
    res = {
        "m": m, "n": n, "t": t, "k": inst.k, "D": inst.D, "seed": seed,
        "prim": hex(ef.prim), "G": inst.G, "g_tries": inst.g_tries,
        "support": inst.support,
        "lam_first10": inst.lam[:10],
        "R": R,
        "guards": {kk: vv for kk, vv in guards.items() if kk not in ("degs",)},
        "n_degenerate_positions_instances": None,
    }

    # ---- jet profiles & degeneracy
    jetprof = {a: fd[a]["prof"] for a in S}
    deg = {}
    for a in S:
        pr = fd[a]["prof"]
        deg[a] = any(pr[j + 1] == pr[j] for j in range(R))
    res["jetprof"] = {str(a): v for a, v in jetprof.items()}
    res["degenerate_held"] = {str(a): bool(v) for a, v in deg.items()}
    res["n_degenerate_held"] = int(sum(deg.values()))

    # ---- family evaluation nullities N_eval(c)
    Nf = {}
    for c in range(1, cmax + 1):
        rows = []
        for a in S[:c]:
            w0 = lucas_w(ef, t2p3 - 1, 0, a)
            rows.append(w0)
        rk = ef.rank(np.array(rows, dtype=np.uint16)) if rows else 0
        Nf[c] = t2p3 - rk
    res["N_eval"] = Nf

    # ---- family restricted system (direct build): ranks and nullities
    blocks_fam = family_restricted_system(inst, ef, fd, R)
    per_point_rank = {}
    rows_all = []
    for a in S:
        rk = 0
        for j in range(R):
            B = blocks_fam[(a, j)]
            r = ef.rank(B) if B.shape[0] else 0
            per_point_rank[(a, j)] = r
            rk = max(rk, r)
            rows_all.append(B)
    res["rank_point_order"] = {f"{a},{j}": int(v) for (a, j), v in per_point_rank.items()}
    # escalation check: any rank >= 2
    bad = {kk: vv for kk, vv in per_point_rank.items() if vv >= 2}
    res["REFUTATION_EVENT"] = bool(bad)
    res["refuting_blocks"] = {f"{a},{j}": int(v) for (a, j), v in bad.items()}
    # N_fam(c): nullity of stacked family-restricted rows for first c points
    Nfam = {}
    for c in range(1, cmax + 1):
        stack = []
        for a in S[:c]:
            for j in range(R):
                stack.append(blocks_fam[(a, j)])
        if stack:
            M = np.vstack(stack)
            Nfam[c] = t2p3 - (ef.rank(M) if M.shape[0] else 0)
        else:
            Nfam[c] = t2p3
    res["N_fam"] = Nfam

    # ---- derived prediction for cross-check:
    # contribution of a = 1 if any even j < R has prof[j+1] > prof[j]
    contrib = {}
    for a in S:
        pr = fd[a]["prof"]
        contrib[a] = any((j % 2 == 0) and (j + 1 <= R) and pr[j + 1] > pr[j]
                         for j in range(R))
    res["derived_contribution"] = {str(a): bool(v) for a, v in contrib.items()}
    # N_fam derived = 2t+3 - #{contributing a in S_c}
    Nfam_pred = {}
    for c in range(1, cmax + 1):
        Nfam_pred[c] = t2p3 - sum(contrib[a] for a in S[:c])
    res["N_fam_derived"] = Nfam_pred
    res["N_fam_matches_derived"] = all(Nfam[c] == Nfam_pred[c] for c in Nfam)

    # ---- full ambient Sol dims (small t only)
    if full_sol and inst.k * (inst.D + 1) <= 60000:
        fblocks = full_flag_blocks(inst, ef, fd, S, R)
        dims = incremental_sol_dims(inst, ef, fblocks, list(range(1, cmax + 1)))
        res["N_flag_full"] = dims
        # waterfall increments
        inc = {}
        prev = inst.k * (inst.D + 1)
        for c in sorted(dims):
            inc[c] = prev - dims[c]
            prev = dims[c]
        res["N_flag_increments"] = inc
        # c_need: least c with N_flag == 1
        cn = None
        for c in sorted(dims):
            if dims[c] == 1:
                cn = c
                break
        res["c_need_full"] = cn
        res["c_need_eq_2t3"] = (cn == t2p3) if cn is not None else None

    # ---- lower-bound certificate at c = 2t+2: T(A_S) in flags
    cert = lower_bound_certificate(inst, ef, fd, S, R)
    res["lb_cert"] = cert

    return res, inst, ef, fd


def lower_bound_certificate(inst, ef, fd, S, R):
    """A_S = prod_{a in S_{2t+2}} (Z - a). Verify T(A_S) jets in flags at
    every held point (depth R). Returns structured result."""
    gf = inst.gf
    t2p2 = 2 * inst.t + 2
    if len(S) < t2p2:
        return {"skipped": "need 2t+2 held points"}
    Sc = S[:t2p2]
    AS = [1]
    for a in Sc:
        AS = gf.ptrim(gf.pmul(AS, [a, 1]))
    okdeg = gf.pdeg(AS) == t2p2
    # T(A_S): B_A, Q_A, rho
    RG = compute_R_G(inst, ef)
    G2 = gf.pmul(inst.G, inst.G)
    BA = gf.pmod(gf.pmul(AS, RG), G2)
    rho = 0
    if inst.D % 2 == 0:
        # [Z^{2t+2}] A_S = lc = 1 ; rho = gamma^{-2}
        gamma = inst.G[inst.t]
        rho = ef.INV[ef.MUL[gamma, gamma]]
    # jets of T(A_S) at each held point and membership in flags
    checks = {}
    pts = Sc if len(Sc) <= 40 else S[:t2p2]  # check at all c = 2t+2 pts
    for j in range(R):
        pass
    # F jets at needed points (S subset) already in fd; need at Sc points
    need = [a for a in Sc] 
    FJ = jets_F_at_points(inst, ef, need, R)
    # T(A_S) jets at each a: num = A_S F' + BA F ; jets; times invG2; + rho ZF
    ASof = np.array(AS, dtype=np.uint16)
    BAof = np.array(gf.ptrim(BA), dtype=np.uint16)
    all_ok = True
    detail = {}
    for a in need:
        d = fd[a]
        FJm = FJ[a]
        GJ = poly_jets(ef, np.array(inst.G, dtype=np.uint16), a, R + 1)
        G2J = conv(ef, GJ, GJ, R + 1)
        invG2 = series_inv_jets(ef, G2J, R)
        ASJ = poly_jets_vec(ef, ASof, a, R)
        BAJ = poly_jets_vec(ef, BAof, a, R)
        ok = True
        for j in range(R):
            # num jets at j: sum_{p+q=j} ASJ_p FJp_q + BAJ_p FJ_q
            num = np.zeros(inst.k, dtype=np.uint16)
            for p in range(j + 1):
                q = j - p
                # F' jet q = (q+1 mod 2) FJ_{q+1}
                if q <= R - 1:
                    s = (q + 1) % 2
                    fpq = ef.MUL[s, FJm[:, q + 1]] if s else np.zeros(inst.k, dtype=np.uint16)
                    num ^= ef.MUL[ASJ[p], fpq]
                num ^= ef.MUL[BAJ[p], FJm[:, q]]
            # times invG2
            qjet = np.zeros(inst.k, dtype=np.uint16)
            for e in range(j + 1):
                qjet ^= ef.MUL[invG2[e], num] if e == j else qjet  # CAREFUL below
            # fix: conv: jet_j(Q) = sum_{e+q=j} invG2_e * num_q
            numjets = []
            for ee in range(j + 1):
                Nv = np.zeros(inst.k, dtype=np.uint16)
                for p in range(ee + 1):
                    q = ee - p
                    if q <= R - 1:
                        s = (q + 1) % 2
                        fpq = ef.MUL[s, FJm[:, q + 1]] if s else np.zeros(inst.k, dtype=np.uint16)
                        Nv ^= ef.MUL[ASJ[p], fpq]
                    Nv ^= ef.MUL[BAJ[p], FJm[:, q]]
                numjets.append(Nv)
            qjet = np.zeros(inst.k, dtype=np.uint16)
            for e in range(j + 1):
                qjet ^= ef.MUL[invG2[e], numjets[j - e]]
            # rho Z F: a FJ_j + FJ_{j-1}
            tjet = ef.MUL[rho, (ef.MUL[a, FJm[:, j]] ^ (FJm[:, j - 1] if j >= 1 else 0))]
            tj = qjet ^ tjet
            # membership in V^[j]: rank of flagbase + tj unchanged
            W = d["flagbase"][j]
            aug = np.vstack([W, tj[None, :]]).astype(np.uint16)
            r2 = ef.rank(aug)
            r1 = ef.rank(W) if W.shape[0] else 0
            if r2 != r1:
                ok = False
                detail[f"{a},{j}"] = "OUTSIDE"
                break
        checks[str(a)] = ok
        all_ok = all_ok and ok
    return {"A_S": AS, "deg_ok": okdeg, "T_AS_in_flags_all_points": all_ok,
            "per_point": checks}
