"""Gate A instance construction for binary Goppa codes (mceliece target).

Builds, for (m, n, t):
  E   = F_{2^m}          (log/exp tables, byte-exact)
  G   = seeded random monic irreducible degree-t poly over E (Rabin-tested;
          irreducible over a field => squarefree)
  L   = support: n distinct elements of E taken as a prefix of a seeded
          shuffle of a fixed enumeration of E (full enumeration must have
          >= n elements; n <= 2^m always holds)
  Pi  = prod (Z - a_i);  D = n - 2t - 1
  lam_i = G(a_i)^2 / Pi'(a_i)          (Apon eq. (1) with kappa = 1;
          GIJS Fact 2.3 lambda')
  code = binary Goppa code Gamma(L, G) via parity rows h_i = a_i^j / G(a_i)
  Y    = k x n binary generator matrix (k = n - mt) from the nullspace
  F    = (f_1..f_k), f_j the unique poly deg <= D with
          lam_i * f_j(a_i) = Y[j][i] in F2  (equivalently f_j(a_i) = Y[j][i]/lam_i)

Guards (all exact; from pre_statement.md):
  (alpha) Pi F' + Pi' F = G^2 F^{(2)}   coefficientwise (Apon Lemma 3)
  (beta)  Delta_{p,q} = f_p f_q' + f_q f_p' != 0 for an exhibited pair
  (gamma) k = n - mt
  (delta) lam_i F(a_i) in F2^k exactly, for every i
  (eps)   gcd(f_1..f_k) is a nonzero constant; max_j deg f_j = D exactly

If any guard fails the instance is DEGENERATE under the seed stream rule
of pre_statement.md section 6: log + advance seed + regenerate.
"""
from __future__ import annotations

import hashlib
import json
import random

import numpy as np

from gfield import GF, _factor


# --------------------------------------------------------------------------
# F2 linear algebra (numpy, exact)
# --------------------------------------------------------------------------

def gf2_nullspace(A: np.ndarray):
    """columns basis of nullspace of binary matrix A (rows x cols)."""
    A = (np.asarray(A, dtype=np.uint8) % 2).copy()
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        pr = None
        for rr in range(r, rows):
            if A[rr, c]:
                pr = rr
                break
        if pr is None:
            continue
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        A[nzr] ^= A[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    free = np.array([c for c in range(cols) if c not in set(piv)], dtype=int)
    basis = []
    for fc in free:
        v = np.zeros(cols, dtype=np.uint8)
        v[fc] = 1
        for ri, pc in enumerate(piv):
            if A[ri, fc]:
                v[pc] = 1
        basis.append(v)
    return basis, piv


# --------------------------------------------------------------------------
# instance
# --------------------------------------------------------------------------

class Instance:
    def __init__(self, m: int, n: int, t: int, seed: int):
        self.m, self.n, self.t, self.seed = m, n, t, seed
        self.gf = GF(m)
        self.k = n - m * t
        self.D = n - 2 * t - 1
        assert self.k > 0
        self.build()

    # ---- construction -----------------------------------------------------
    def build(self):
        gf = self.gf
        m, n, t = self.m, self.n, self.t
        self.log_events = []

        # 1. irreducible monic degree-t G (seeded)
        rng = random.Random(self.seed)
        tries = 0
        while True:
            tries += 1
            G = [rng.randrange(gf.q) for _ in range(t)] + [1]
            if gf.pis_irreducible(G):
                break
        self.G, self.g_tries = G, tries

        # 2. support: prefix of a seeded shuffle of [0..2^m)
        idx = list(range(gf.q))
        rng.shuffle(idx)
        support = idx[:n]
        self.support = support

        # 3. Pi = prod (Z - a_i) = prod (Z + a_i) in char 2
        Pi = [1]
        for a in support:
            Pi = gf.ptrim(gf.pmul(Pi, [a, 1]))
        self.Pi = Pi
        self.PiD = gf.pderiv(Pi)
        assert gf.pdeg(self.PiD) >= 0  # Pi' nonzero since support distinct

        PiDa = [gf.peval(self.PiD, a) for a in support]
        assert all(v != 0 for v in PiDa), "support must be distinct"
        lam = []
        for i, a in enumerate(support):
            ga = gf.peval(G, a)
            assert ga != 0, "G must not vanish on support"
            lam.append(gf.mul(gf.mul(ga, ga), gf.inv(PiDa[i])))
        self.lam = lam
        self.Ga = [gf.peval(G, a) for a in support]

        # 4. parity checks h_{j,i} = a_i^j / G(a_i), j < t; lift to F2^m rows
        Hrows = []
        for j in range(t):
            h = [gf.mul(gf.pow(a, j), gf.inv(self.Ga[i])) for i, a in enumerate(support)]
            for b in range(m):
                Hrows.append([(h[i] >> b) & 1 for i in range(n)])
        ns, piv = gf2_nullspace(np.array(Hrows, dtype=np.uint8))
        assert len(ns) >= self.k, (len(ns), self.k)
        # use the first k basis vectors (rows are in RREF; pivots chosen
        # lexicographically) — deterministic
        self.Y = np.array([list(v) for v in ns[: self.k]], dtype=np.uint8)  # k x n

        # 5. interpolate F: f_j(a_i) = Y[j][i] / lam_i  (deg <= D; unique since D < n)
        Ls = self._lagrange()
        F = []
        for j in range(self.k):
            row = self.Y[j]
            f = []
            for i in range(n):
                if row[i]:
                    f = gf.padd(f, gf.pscale(Ls[i], gf.inv(lam[i])))
            F.append(gf.ptrim(f))
        self.F = F

        # 6. guards
        self.guards = self._guards()

    def _lagrange(self):
        """L_i = Pi/(Z - a_i)/Pi'(a_i); L_i(a_j) = delta_ij (exact)."""
        gf = self.gf
        Ls = []
        for i, a in enumerate(self.support):
            q, r = gf.pdivmod(self.Pi, [self.support[i], 1])
            assert gf.pdeg(r) < 0
            denom = self.gf.peval(self.PiD, a)
            Ls.append(gf.pscale(q, gf.inv(denom)))
        return Ls

    # ---- guards ------------------------------------------------------------
    def _guards(self):
        gf = self.gf
        g = {}

        # (delta) lam_i * F(a_i) equals Y[:, i] in F_2^k exactly
        exact = True
        evals = []
        for i, a in enumerate(self.support):
            col = [gf.peval(f, a) for f in self.F]
            evals.append(col)
            li = self.lam[i]
            for j in range(self.k):
                w = gf.mul(li, col[j])
                if w >= (1 << self.m):
                    exact = False
                    break
                if int(self.Y[j][i]) != (w & 1):
                    exact = False
                    break
            if not exact:
                break
        g["delta_yi_binary"] = exact
        self.F_eval = evals  # F(a_i) per support point (k-vectors)

        # (alpha) Pi F' + Pi' F = G^2 F^(2) coordinatewise
        Gsq = gf.pmul(self.G, self.G)
        alph = True
        for f in self.F:
            fp = gf.pderiv(f)
            lhs = gf.ptrim(gf.padd(gf.pmul(self.Pi, fp), gf.pmul(self.PiD, f)))
            rhs = gf.ptrim(gf.pmul(Gsq, gf.pmul(f, f)))
            if lhs != rhs:
                alph = False
                break
        g["alpha_identity"] = alph

        # (beta) Delta_{p,q} != 0 for some pair; exhibit first found
        found = None
        deltapq = None
        for p in range(self.k):
            fp = gf.pderiv(self.F[p])
            if fp == []:
                continue
            for q in range(p + 1, self.k):
                fq = self.F[q]
                fqp = gf.pderiv(fq)
                d = gf.ptrim(gf.padd(gf.pmul(self.F[p], fqp), gf.pmul(fq, fp)))
                if d != []:
                    found, deltapq = (p, q), gf.pdeg(d)
                    break
            if found:
                break
        g["beta_delta_nonzero"] = found is not None
        g["beta_pair"] = found
        g["beta_delta_deg"] = deltapq

        # (gamma) k = n - mt
        g["gamma_k"] = self.k == self.n - self.m * self.t

        # (eps) gcd of coordinates is a constant; max degree exactly D
        gcur = self.F[0][:]
        for f in self.F[1:]:
            gcur = gf.ptrim(gf.pgcd(gcur, f))
            if gf.pdeg(gcur) == 0:
                gcur = [1]
                break
        g["eps_gcd_const"] = gf.pdeg(gcur) <= 0
        degs = [gf.pdeg(f) for f in self.F]
        g["eps_maxdeg_is_D"] = max(degs) == self.D
        g["eps_degs_min_max"] = (min(degs), max(degs))

        # extra: record degree histogram summary
        g["degs"] = degs
        g["d_uniform_input"] = None
        return g

    # ---- serialization ------------------------------------------------------
    def summary(self) -> dict:
        return {
            "m": self.m, "n": self.n, "t": self.t, "k": self.k, "D": self.D,
            "seed": self.seed, "g_tries": self.g_tries,
            "G": self.G, "G_hex": [hex(c) for c in self.G],
            "support_first20": self.support[:20],
            "guards": {k: v for k, v in self.guards.items() if k != "degs"},
            "prim": hex(self.gf.prim),
        }

    def full_record(self) -> dict:
        """records G, support, lam (enough to rebuild F bit-exactly via the
        deterministic construction in src/instance.py at the recorded seed)."""
        return {
            "m": self.m, "n": self.n, "t": self.t, "seed": self.seed,
            "prim": hex(self.gf.prim),
            "G": self.G,
            "support": self.support,
            "lam": self.lam,
            "k": self.k, "D": self.D,
        }
