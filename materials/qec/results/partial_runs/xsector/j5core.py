"""Shared GF(2) machinery for the J.5 attempt (Theorem J.5 of xsector notes).

Conventions: index k = i*m + j on the ell x m torus.  Row 0 of a polynomial
matrix P is the coefficient vector of its polynomial p (row-lambda = lambda * p,
convolution).  Reversal J flips the lattice coordinates: poly(P^T) = p[J-idx].

Key objects:
  M(C,D) = A C^T + B D^T          (validity: M symmetric, i.e. poly m = m[J])
  W      = rowspace(HX) cap ker[CD] = {lam HX : M lam^T = 0}
  S_X(Q) = {lam HX : lam CD in S_Z + Delta},  Delta = L CD, L = left-ker H_X
  demoted x = lam HX in W not in S_X(Q);  partner z = lam CD in ker H_X
"""
from __future__ import annotations

import itertools

import numpy as np

from qec_research.codes.bicycle import monomial_matrix, poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np, matmul as gf2matmul  # noqa: E402


class Lattice:
    def __init__(self, ell, m):
        self.ell, self.m, self.dim = ell, m, ell * m
        X = np.zeros((self.dim, self.dim, self.dim), dtype=np.uint8)
        for k, (a, b) in enumerate(itertools.product(range(ell), range(m))):
            X[k] = monomial_matrix(ell, m, a, b)
        self.X = X
        # reversal permutation J: coordinate k=(i,j) -> (-i,-j)
        J = np.zeros((self.dim, self.dim), dtype=np.uint8)
        for k, (a, b) in enumerate(itertools.product(range(ell), range(m))):
            J[((-a) % ell) * m + ((-b) % m), k] = 1
        self.J = J

    def mat(self, termset):
        return poly_matrix(self.ell, self.m, [tuple(t) for t in termset])

    def sparse(self, coeff_vec):
        out = np.zeros((self.dim, self.dim), np.uint8)
        for k in np.flatnonzero(coeff_vec):
            out ^= self.X[k]
        return out

    def terms_of(self, M):
        c = M[0]
        return sorted(((int(k) // self.m) % self.ell, int(k) % self.m)
                      for k in np.flatnonzero(c))


class ParentData:
    """A BB parent plus all the derived spaces used by the J.5 analysis."""

    def __init__(self, lat: Lattice, A_terms, B_terms):
        self.lat = lat
        self.A_terms = tuple(map(tuple, A_terms))
        self.B_terms = tuple(map(tuple, B_terms))
        dim, N = lat.dim, 2 * lat.dim
        self.dim, self.N = dim, N
        A = lat.mat(self.A_terms)
        B = lat.mat(self.B_terms)
        self.A, self.B = A, B
        HX = np.hstack([A, B])
        HZ = np.hstack([B.T, A.T])
        self.HX, self.HZ = HX, HZ
        self.rHX, self.rHZ = rank_np(HX), rank_np(HZ)
        self.kP = N - self.rHX - self.rHZ
        self.L = nullspace_np(HX.T)          # left kernel of H_X, dim = dim - rHX
        self.row0_maps = self._row0_maps()   # (c,d) -> poly-row of M
        self.validity = self._validity()     # V basis (rows are (c,d) coeff vectors)
        self.syzygy = nullspace_np(self.row0_maps)  # {M = 0} perturbations
        self.sigma = self.syzygy.shape[0] - self.rHZ  # nontrivial excess

    def _row0_maps(self):
        """(c,d) -> row 0 of M = A C^T + B D^T, as an dim x 2*dim matrix.

        row0(A C^T) = J A c  (machine-pinned identity; J = coordinate reversal).
        """
        lat, dim = self.lat, self.dim
        J = lat.J
        PA = (J @ self.A) % 2
        PB = (J @ self.B) % 2
        return np.hstack([PA, PB]).astype(np.uint8)

    def _validity(self):
        """Basis of V = {(c,d) : M symmetric} = ker(c,d -> m + J m) with
        m = row0(M).  (I+J)(JA c + JB d) = (I+J)A c + (I+J)B d."""
        lat, dim = self.lat, self.dim
        J = lat.J
        IJ = (np.eye(dim, dtype=np.uint8) + J) % 2
        C1 = (IJ @ self.A) % 2
        C2 = (IJ @ self.B) % 2
        N = np.hstack([C1, C2]).astype(np.uint8)
        return nullspace_np(N)

    def M_of(self, C, D):
        return (gf2matmul(self.A, C.T) ^ gf2matmul(self.B, D.T)).astype(np.uint8)

    def CD_of_coeff(self, cd):
        dim = self.dim
        return np.hstack([self.lat.sparse(cd[:dim]), self.lat.sparse(cd[dim:])])


def light_stabilizers(parent: ParentData, d_lb, max_lam=3, structured=True):
    """Enumerate x = lam H_X of weight < d_lb with |lam| <= max_lam (+ structured).

    Returns list of dicts {x(bytes), w, lam(list of row indices)} sorted by w.
    """
    HX, dim = parent.HX, parent.dim
    R = HX.astype(np.uint8)
    out = {}
    wR = R.sum(axis=1).astype(np.int64)

    def record(lam_idx, xvec):
        w = int(xvec.sum())
        if 0 < w < d_lb:
            key = xvec.tobytes()
            if key not in out:
                out[key] = {"x": xvec.copy(), "w": w, "lam": tuple(int(i) for i in lam_idx)}

    for i in range(dim):
        record((i,), R[i])
    # pairs
    if max_lam >= 2:
        intR = R.astype(np.int32)
        ov = (intR @ intR.T)  # overlaps
        for i in range(dim):
            for j in range(i + 1, dim):
                w = wR[i] + wR[j] - 2 * ov[i, j]
                if 0 < w < d_lb:
                    record((i, j), R[i] ^ R[j])
    # triples
    if max_lam >= 3:
        for i in range(dim):
            for j in range(i + 1, dim):
                xij = R[i] ^ R[j]
                wij = wR[i] + wR[j] - 2 * int((R[i].astype(np.int32) @ R[j]))
                for k in range(j + 1, dim):
                    wk = wij + wR[k] - 2 * int(xij.astype(np.int32) @ R[k])
                    if 0 < wk < d_lb:
                        record((i, j, k), xij ^ R[k])
    if structured:
        ell, m, lat = parent.lat.ell, parent.lat.m, parent.lat
        gen = []
        gen.append(np.ones(dim, np.uint8))  # full torus sum (rank drop check)
        for i in range(ell):  # x-line at y=i
            v = np.zeros(dim, np.uint8)
            for a in range(ell):
                v[a * m + i % m] = 1
            gen.append(v)
        for j in range(ell):  # y-line at x=j
            v = np.zeros(dim, np.uint8)
            for b in range(m):
                v[j * m + b] = 1
            gen.append(v)
        if ell % 2 == 0:
            h = np.zeros(dim, np.uint8)  # (1 + x^{ell/2})
            for b in range(m):
                h[0 * m + b] = 1
                h[(ell // 2) * m + b] = 1
            gen.append(h)
        if m % 2 == 0:
            h = np.zeros(dim, np.uint8)  # (1 + y^{m/2})
            h[0] = 1
            h[m // 2] = 1
            gen.append(h)
        for v in gen:
            w = int(v.sum())
            if w == 0:
                continue
            # all translates
            Mv = lat.sparse(v)
            for k in range(dim):
                lam = Mv[k]
                x = (lam @ HX) % 2
                record(tuple(np.flatnonzero(lam).tolist()), x)
    res = sorted(out.values(), key=lambda d: d["w"])
    return res


def demotion_spaces(parent: ParentData, C, D):
    """Return dict with kerM, Delta-basis, SZplusDelta nullspace probe, S_X(Q) probe."""
    dim = parent.dim
    CD = np.hstack([C, D])
    M = parent.M_of(C, D)
    kerM = nullspace_np(M.T)
    L = parent.L
    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, parent.N), np.uint8)
    SZD = np.vstack([parent.HZ, Delta]) if Delta.shape[0] else parent.HZ.copy()
    probe_SZD = nullspace_np(SZD)  # z in S_Z + Delta  <=>  z . probe^T == 0
    return {"kerM": kerM, "M": M, "CD": CD, "SZD": SZD, "probe_SZD": probe_SZD,
            "delta_bar": int(rank_np(SZD) - parent.rHZ)}


def z_in_SZplusDelta(z, probe_SZD):
    return not bool(((z @ probe_SZD.T) % 2).any())


def sharp_probe(parent: ParentData, lights, d_lb, rng, max_lams=64, samples_per=8):
    """For each light stabilizer x0 = lam0 HX with w < d_lb, build
    V'(lam0) = V cap {M lam0^T = 0} and test demotion of x0 on samples.
    Returns (probes, witness_candidates)."""
    dim = parent.dim
    V = parent.validity
    P0 = parent.row0_maps
    lam_rows = parent.HX  # dim x N
    probes = []
    cands = []
    for rec in lights[:max_lams]:
        lam0 = np.zeros(dim, np.uint8)
        lam0[list(rec["lam"])] = 1
        x0 = rec["x"]
        w0 = rec["w"]
        # constraint: lambda0 * m = 0 where m = row0(M); row0 map P0 (dim x 2dim)
        Lam0 = parent.lat.sparse(lam0)
        K0 = (Lam0 @ P0) % 2            # (c,d) -> lambda0 * m
        K0V = (K0 @ V.T) % 2            # restricted to V coordinates
        selker = nullspace_np(K0V)
        dimVp = selker.shape[0]
        probes.append({"w": int(w0), "lam": list(rec["lam"]), "dim_Vp": int(dimVp)})
        if dimVp == 0:
            continue
        # sample: basis vectors + random combos
        todo = []
        for r in range(dimVp):
            todo.append(selker[r])
        for _ in range(samples_per):
            sel = (rng.integers(0, 2, dimVp) % 2).astype(np.uint8)
            if sel.any():
                todo.append((sel @ selker) % 2)
        for sel in todo:
            cd = (sel @ V) % 2
            z1 = (lam0 @ CD_of(parent, cd)) % 2
            # demoted iff z not in S_Z + Delta for THIS (C,D)
            sp = demotion_spaces(parent, *split_CD(parent, cd))
            if not z_in_SZplusDelta(z1, sp["probe_SZD"]):
                cands.append({"w": int(w0), "lam": list(rec["lam"]),
                              "cd": cd.tolist(), "w_z": int(z1.sum())})
                break  # one candidate per lambda0 is enough for the census
    return probes, cands


def CD_of(parent, cd):
    dim = parent.dim
    return np.hstack([parent.lat.sparse(cd[:dim]), parent.lat.sparse(cd[dim:])])


def split_CD(parent, cd):
    dim = parent.dim
    return parent.lat.sparse(cd[:dim]), parent.lat.sparse(cd[dim:])


def syzygy_demotion_scan(parent: ParentData, lights, d_lb, rng, samples=256):
    """Scan the M=0 perturbation space (syzygy kernel) for demoted light
    stabilizers.  Returns (tested, candidates)."""
    S = parent.syzygy
    dimS = S.shape[0]
    cands = []
    tested = 0
    if dimS == 0 or not lights:
        return tested, cands
    todo = [S[r] for r in range(dimS)]
    for _ in range(samples):
        sel = (rng.integers(0, 2, dimS)).astype(np.uint8)
        if sel.any():
            todo.append((sel @ S) % 2)
    for cd in todo:
        C, D = split_CD(parent, cd)
        sp = demotion_spaces(parent, C, D)
        M = sp["M"]
        if M.any():
            continue  # paranoia: syzygy must give M=0 exactly
        tested += 1
        for rec in lights:
            if rec["w"] >= d_lb:
                continue
            lam0 = np.zeros(parent.dim, np.uint8)
            lam0[list(rec["lam"])] = 1
            z1 = (lam0 @ sp["CD"]) % 2
            if not z_in_SZplusDelta(z1, sp["probe_SZD"]):
                cands.append({"w": int(rec["w"]), "lam": list(rec["lam"]),
                              "cd": cd.tolist(), "w_z": int(z1.sum())})
                break
    return tested, cands


def batch_demotion_test(parent: ParentData, C, D, lights):
    """Given a valid (C,D), test demotion of every light stabilizer at once.

    A light x = lam HX demotes iff (i) lam M = 0  (x in ker[CD], i.e. x in W)
    and (ii) partner z = lam CD not in S_Z + Delta.
    Returns (list of dicts {w, lam, w_z}, spaces-dict)."""
    dim = parent.dim
    sp = demotion_spaces(parent, C, D)
    CD = sp["CD"]
    probe = sp["probe_SZD"]
    M = sp["M"]
    if not lights:
        return [], sp
    lam_rows = np.zeros((len(lights), dim), np.uint8)
    for i, rec in enumerate(lights):
        lam_rows[i, list(rec["lam"])] = 1
    in_W = ~((lam_rows @ M) % 2).any(axis=1)      # x in ker[CD]
    if not in_W.any():
        return [], sp
    Z = (lam_rows @ CD) % 2          # partner vectors
    dem = in_W & ((Z @ probe.T) % 2).any(axis=1)
    out = []
    for i in np.nonzero(dem)[0]:
        x = (lam_rows[i] @ parent.HX) % 2
        w = int(x.sum())
        if w:
            out.append({"w": w, "lam": list(lights[i]["lam"]), "w_z": int(Z[i].sum())})
    return out, sp


def verify_full(parent: ParentData, cd, lam_sup, d_lb):
    """Dual-path verification of a single demotion event. Returns dict."""
    dim = parent.dim
    cd = np.asarray(cd, dtype=np.uint8)
    C, D = split_CD(parent, cd)
    CD = np.hstack([C, D])
    lam0 = np.zeros(dim, np.uint8)
    lam0[list(lam_sup)] = 1
    x0 = (lam0 @ parent.HX) % 2
    z0 = (lam0 @ CD) % 2
    from qec_research.gf2.linalg import matmul as gf2matmul
    M = (gf2matmul(parent.A, C.T) ^ gf2matmul(parent.B, D.T)).astype(np.uint8)
    sp = demotion_spaces(parent, C, D)
    # path A: z not in S_Z + Delta
    inS_A = z_in_SZplusDelta(z0, sp["probe_SZD"])
    # path B (exp-044 formulation): x0 in S_X(Q) via Wperp
    Wperp = nullspace_np(sp["SZD"])
    Lam = nullspace_np(((CD @ Wperp.T) % 2).T)
    SXQ = (Lam @ parent.HX) % 2
    ns = nullspace_np(SXQ)
    x0_in_SXQ = not bool(((x0 @ ns.T) % 2).sum())
    rhoX = int(rank_np(np.vstack([parent.HZ, CD])) - parent.rHZ)
    res = {
        "M_symmetric": not bool((M ^ M.T).any()),
        "M_is_zero": not bool(M.any()),
        "w_x": int(x0.sum()), "w_z": int(z0.sum()),
        "x_in_XcenQ": (not bool(((x0 @ parent.HZ.T) % 2).any())
                       and not bool(((x0 @ CD.T) % 2).any())),
        "x_in_rowspace_HX": not bool(((x0 @ nullspace_np(parent.HX).T) % 2).sum()),
        "z_in_kerHX": not bool(((z0 @ parent.HX.T) % 2).any()),
        "demoted_A": not inS_A,
        "demoted_B": not x0_in_SXQ,
        "paths_agree": (not inS_A) == (not x0_in_SXQ),
        "rho_X": rhoX,
        "XcenQ_dim": int(2 * dim - parent.rHZ - rhoX),
        "delta_bar": sp["delta_bar"],
        "k_Q": int(parent.kP - sp["delta_bar"]),
        "below_bound": int(x0.sum()) < (d_lb if d_lb is not None else 1 << 60),
        "A_terms": [list(t) for t in parent.A_terms],
        "B_terms": [list(t) for t in parent.B_terms],
        "C_terms": [list(t) for t in parent.lat.terms_of(C)],
        "D_terms": [list(t) for t in parent.lat.terms_of(D)],
        "lam_support": [int(s) for s in lam_sup],
        "x0_support": sorted(int(i) for i in np.flatnonzero(x0)),
        "z0_support": sorted(int(i) for i in np.flatnonzero(z0)),
    }
    res["verified"] = bool(
        res["M_symmetric"] and res["x_in_XcenQ"] and res["x_in_rowspace_HX"]
        and res["z_in_kerHX"] and res["demoted_A"] and res["demoted_B"]
        and res["paths_agree"] and res["w_x"] > 0
        and (d_lb is None or res["w_x"] < d_lb))
    return res


def verify_candidate(parent: ParentData, cand, d_lb):
    """Full GF(2) verification of a violation witness candidate."""
    cd = np.array(cand["cd"], dtype=np.uint8)
    C, D = split_CD(parent, cd)
    lam0 = np.zeros(parent.dim, np.uint8)
    lam0[list(cand["lam"])] = 1
    x0 = (lam0 @ parent.HX) % 2
    z1 = (lam0 @ np.hstack([C, D])) % 2
    sp = demotion_spaces(parent, C, D)
    M = sp["M"]
    checks = {
        "M_symmetric": not bool((M ^ M.T).any()),
        "x_in_XcenQ": (not bool(((x0 @ parent.HZ.T) % 2).any())
                       and not bool(((x0 @ np.hstack([C, D]).T) % 2).any())),
        "x_in_rowspace_HX": not bool(((x0 @ nullspace_np(parent.HX).T) % 2).sum()),
        "z_not_in_SZ_plus_Delta": not z_in_SZplusDelta(z1, sp["probe_SZD"]),
        "z_in_kerHX": not bool(((z1 @ parent.HX.T) % 2).any()),
        "w_x": int(x0.sum()),
        "w_z": int(z1.sum()),
        "below_bound": int(x0.sum()) < d_lb,
        "x_nonzero": bool(x0.any()),
    }
    checks["verified"] = all([
        checks["M_symmetric"], checks["x_in_XcenQ"], checks["x_in_rowspace_HX"],
        checks["z_not_in_SZ_plus_Delta"], checks["below_bound"], checks["x_nonzero"],
    ])
    checks["x_support"] = sorted(int(i) for i in np.flatnonzero(x0))
    checks["C_terms"] = parent.lat.terms_of(C)
    checks["D_terms"] = parent.lat.terms_of(D)
    return checks
