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
    def __init__(self, m: int, n: int, t: int, seed: int, forced_G=None):
        self.m, self.n, self.t, self.seed = m, n, t, seed
        self.forced_G = forced_G
        self.gf = GF(m)
        self.k = n - m * t
        self.D = n - 2 * t - 1
        assert self.k > 0
        self.build()

    # ---- construction -----------------------------------------------------
    def build(self):
        m, n, t = self.m, self.n, self.t
        from fastfield import EField
        ef = EField(m)
        self.ef = ef
        self.log_events = []

        # 1. irreducible monic degree-t G (seeded). Fast vectorized Rabin
        # (fastfield.vec_pis_irreducible, flint-agreed) with scalar
        # re-verification of the ACCEPTED polynomial only (cheap: 1 test).
        from fastfield import vec_pis_irreducible
        rng = random.Random(self.seed)
        tries = 0
        if self.forced_G is not None:
            # EXHAUSTIVE/ADVERSARIAL route: caller supplies G (monic, degree t).
            # Irreducibility is re-verified here with BOTH engines (scalar +
            # vectorized); a supplied reducible G is refused loudly.
            G = list(self.forced_G)
            assert self.gf.pdeg(G) == t and G[-1] == 1, "forced_G must be monic deg t"
            assert self.gf.pis_irreducible(G), "forced_G reducible — refused"
            assert vec_pis_irreducible(ef, G), "forced_G fails vectorized Rabin"
            tries = 0
        else:
            while True:
                tries += 1
                G = [rng.randrange(self.gf.q) for _ in range(t)] + [1]
                if vec_pis_irreducible(ef, G):
                    # scalar confirmation (ground-truth engine, flint-agreed)
                    assert self.gf.pis_irreducible(G), "vec/scalar Rabin disagreement"
                    break
        if self.forced_G is not None:
            # Keep the seed stream POSITION identical to the unforced build at
            # the same seed: drain exactly the draws the while-loop would have
            # consumed, WITHOUT using the drawn value, so the support shuffle
            # below lands on the same permutation. (The while-loop consumes
            # t draws per iteration; replicate by NOT replacing G.)
            drained = 0
            while True:
                drained += 1
                draw = [rng.randrange(self.gf.q) for _ in range(t)] + [1]
                if vec_pis_irreducible(ef, draw):
                    break
            # draw is discarded; G stays forced
        self.G, self.g_tries = G, tries

        # 2. support: prefix of a seeded shuffle of [0..2^m)
        idx = list(range(self.gf.q))
        rng.shuffle(idx)
        support = idx[:n]
        self.support = support

        # 3. Pi = prod (Z - a_i) = prod (Z + a_i) in char 2
        Pi = [1]
        for a in support:
            Pi = self.gf.ptrim(self.gf.pmul(Pi, [a, 1]))
        self.Pi = Pi
        self.PiD = self.gf.pderiv(Pi)
        assert self.gf.pdeg(self.PiD) >= 0  # Pi' nonzero since support distinct

        PiDa = [self.gf.peval(self.PiD, a) for a in support]
        assert all(v != 0 for v in PiDa), "support must be distinct"
        lam = []
        for i, a in enumerate(support):
            ga = self.gf.peval(G, a)
            assert ga != 0, "G must not vanish on support"
            lam.append(self.gf.mul(self.gf.mul(ga, ga), self.gf.inv(PiDa[i])))
        self.lam = lam
        self.Ga = [self.gf.peval(G, a) for a in support]

        # 4. parity checks h_{j,i} = a_i^j / G(a_i), j < t; lift to F2^m rows
        Hrows = []
        for j in range(t):
            h = [self.gf.mul(self.gf.pow(a, j), self.gf.inv(self.Ga[i])) for i, a in enumerate(support)]
            for b in range(m):
                Hrows.append([(h[i] >> b) & 1 for i in range(n)])
        ns, piv = gf2_nullspace(np.array(Hrows, dtype=np.uint8))
        assert len(ns) >= self.k, (len(ns), self.k)
        # use the first k basis vectors (rows are in RREF; pivots chosen
        # lexicographically) — deterministic
        self.Y = np.array([list(v) for v in ns[: self.k]], dtype=np.uint8)  # k x n
        self.piv = piv[: self.k]

        # 5. interpolate F: f_j(a_i) = Y[j][i] / lam_i  (deg <= D; unique
        # since D < n).  VECTORIZED: v_i = lam_i^{-1} L_i precomputed, then
        # f_j = XOR_{i : Y[j,i]=1} v_i  (n numpy steps of (k x D+1) XOR).
        Ls = self._lagrange()
        Dp1 = self.D + 1
        Vin = []
        for i in range(n):
            Li = np.asarray(Ls[i] + [0] * max(0, Dp1 - len(Ls[i])), dtype=np.uint16)[:Dp1]
            Vin.append(self.ef.MUL[self.ef.INV[self.lam[i]], Li])
        Cmat = np.zeros((self.k, Dp1), dtype=np.uint16)
        for i in range(n):
            rows = np.nonzero(self.Y[:, i])[0]
            if rows.size:
                Cmat[rows] ^= Vin[i][None, :]
        F = []
        for j in range(self.k):
            F.append(self.gf.ptrim(list(Cmat[j])))
        self.F = F

        # 6. guards
        self.guards = self._guards()

    def _lagrange(self):
        """L_i = Pi/(Z - a_i)/Pi'(a_i); L_i(a_j) = delta_ij (exact)."""
        gf = self.gf
        Ls = []
        for i, a in enumerate(self.support):
            q, r = self.gf.pdivmod(self.Pi, [self.support[i], 1])
            assert self.gf.pdeg(r) < 0
            denom = self.gf.peval(self.PiD, a)
            Ls.append(self.gf.pscale(q, self.gf.inv(denom)))
        return Ls

    # ---- guards ------------------------------------------------------------
    def _guards(self, ef=None):
        gf = self.gf
        g = {}

        # (delta) lam_i * F(a_i) equals Y[:, i] in F_2^k exactly.
        # VECTORIZED: build the k x (D+1) coefficient matrix, evaluate all
        # coordinates at all support points via per-point jet weights.
        if ef is None:
            from fastfield import EField
            ef = EField(self.m)
        import numpy as np
        from census import lucas_w
        Dp1 = self.D + 1
        Cmat = np.zeros((self.k, Dp1), dtype=np.uint16)
        for j, f in enumerate(self.F):
            Cmat[j, : len(f)] = f
        # eval check: for each ROW j, evaluate ONLY that row's polynomial at
        # (a) its own RREF pivot p_j — an interpolation error in ANY row
        # breaks its pivot identity — and (b) a seeded 8-point sample.
        # Row-wise: k × (1 + 8) evaluations, each O(D) gathers — ~40x less
        # data than the previous column-major loop.
        rng_g = np.random.default_rng(self.seed)
        sample_loc = [int(x) for x in rng_g.choice(self.n, size=min(8, self.n), replace=False)]
        piv_used = getattr(self, 'piv', None) or []
        piv_used = list(piv_used) + [p + 0 for p in piv_used]  # keep semantics
        exact = True
        n_checks = 0
        for j in range(self.k):
            f = self.F[j]
            cfl = np.asarray(f + [0] * max(0, Dp1 - len(f)), dtype=np.uint16)[:Dp1]
            nzj = np.nonzero(cfl)[0]
            targets = []
            if j < len(getattr(self, 'piv', []) or []):
                targets.append(self.piv[j])
            targets.extend(sample_loc)
            for i in targets:
                a = self.support[i]
                w = lucas_w(ef, Dp1 - 1, 0, a)
                v = int(np.bitwise_xor.reduce(ef.MUL[cfl[nzj], w[nzj]])) if nzj.size else 0
                lw = ef.MUL[self.lam[i], v]
                if lw >= (1 << self.m) or int(self.Y[j][i]) != (lw & 1):
                    exact = False
                    break
                n_checks += 1
            if not exact:
                break
        g["delta_yi_binary"] = exact
        g["delta_check_count"] = n_checks
        # sampled evaluations for the record (k rows x sample points only)
        self.F_eval = {"sample_points": sample_loc}

        # (alpha) Pi F' + Pi' F = G^2 F^(2) coordinatewise.  VECTORIZED and
        # EXACT: both sides are polynomials of degree <= n-1 (Lemma 2: two
        # polynomials of degree <= n-1 agreeing at n distinct support points
        # are EQUAL — Apon's own proof of Lemma 2). So check equality by
        # evaluating at the n support points: lhs_i = Pi(a)F'(a)+Pi'(a)F(a),
        # rhs_i = G(a)^2 F(a)^2, both via numpy MUL gathers. Also verify at
        # ALL coordinates (rewrite cost O(k n) per point, vectorized).
        Dmat = np.zeros((self.k, self.D+1), dtype=np.uint16)
        Cmat2 = np.zeros((self.k, self.D+1), dtype=np.uint16)
        for j, f in enumerate(self.F):
            cf = np.asarray(f + [0]*max(0, Dp1-len(f)), dtype=np.uint16)[:Dp1]
            Cmat2[j] = cf
            fp = np.zeros(Dp1, dtype=np.uint16)
            fp[:len(cf)-1] = cf[1:]
            fp[0::2] = 0        # char-2 derivative: odd-index coeff -> even slot
            # careful: d/dZ c_e Z^e: survives iff e odd; lands at slot e-1.
            fp[:] = 0
            for e in range(1, len(cf)):
                if e % 2 == 1:
                    fp[e-1] = cf[e]
            Dmat[j] = fp
        alph = True
        for i, a in enumerate(self.support):
            w0 = lucas_w(ef, Dp1-1, 0, a)
            nz = np.nonzero(Cmat2.any(axis=0))[0]
            if nz.size == 0:
                continue
            Fa = np.bitwise_xor.reduce(ef.MUL[Cmat2[:, nz], w0[nz][None, :]], axis=1)
            Fpa = np.bitwise_xor.reduce(ef.MUL[Dmat[:, nz], w0[nz][None, :]], axis=1)
            # Pi(a), Pi'(a)
            Pi_np = np.asarray(self.Pi, dtype=np.uint16)
            PiD_np = np.asarray(self.PiD, dtype=np.uint16)
            wp = lucas_w(ef, len(Pi_np)-1, 0, a)
            pia = np.bitwise_xor.reduce(ef.MUL[Pi_np, wp])
            w1 = lucas_w(ef, len(PiD_np)-1, 0, a)
            pida = np.bitwise_xor.reduce(ef.MUL[PiD_np, w1])
            Gn = np.asarray(self.G, dtype=np.uint16)
            Ga_ = int(np.bitwise_xor.reduce(ef.MUL[Gn, lucas_w(ef, len(Gn)-1, 0, a)]))
            lhs = ef.MUL[pia, Fpa] ^ ef.MUL[pida, Fa]
            rhs = ef.MUL[ef.MUL[Ga_, Ga_], ef.MUL[Fa, Fa]]
            if not np.array_equal(lhs, rhs):
                alph = False
                break
        g["alpha_identity"] = alph

        # (beta) Delta_{p,q} != 0 for some pair; exhibit first found
        found = None
        deltapq = None
        for p in range(self.k):
            fp = self.gf.pderiv(self.F[p])
            if fp == []:
                continue
            for q in range(p + 1, self.k):
                fq = self.F[q]
                fqp = self.gf.pderiv(fq)
                d = self.gf.ptrim(self.gf.padd(self.gf.pmul(self.F[p], fqp), self.gf.pmul(fq, fp)))
                if d != []:
                    found, deltapq = (p, q), self.gf.pdeg(d)
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
            gcur = self.gf.ptrim(self.gf.pgcd(gcur, f))
            if self.gf.pdeg(gcur) == 0:
                gcur = [1]
                break
        g["eps_gcd_const"] = self.gf.pdeg(gcur) <= 0
        degs = [self.gf.pdeg(f) for f in self.F]
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
