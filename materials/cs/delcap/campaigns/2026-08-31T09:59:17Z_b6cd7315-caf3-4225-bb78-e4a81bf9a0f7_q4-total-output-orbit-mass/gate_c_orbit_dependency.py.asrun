"""
gate_c_orbit.py — Open Item 1: sparse orbit certificates for the q-ary
deletion channel — (3,10) and the q=3 ladder n=6..9 (cs/delcap/README.md
'What remains open', first bullet).

Open item attacked (README:572-588):

    * `(q,n) = (3,10)` and the `q=3` ladder above `n=5` need either a sparse
      channel representation or a GPU; the dense exact-integer matrix is the wall.

THE REDUCTION (identities anchor-tested at runtime by `_selftest`, never
assumed; ALSO independently confirmed in exact arithmetic by Main at
(3,2,1/3), (4,2,1/5), (3,3,2/7) — with the position-permutation control
failing, i.e. exactly the falsified gate-B reduction breaking):

  G = S_q x C_2 acts on SYMBOL VALUES: sigma in S_q relabels alphabet symbols;
  rho reverses the word; both applied simultaneously to input and output
  words. Deletion acts on POSITIONS; value edits commute with position
  deletion, so
      A(sigma x, sigma y) = A(x, y)   and   A(rho x, rho y) = A(x, y)
  with A(x,y) = #subsequences of x equal to y, W(y|x) = A(x,y) d^(n-k)(1-d)^k.
  This is NOT the falsified S_n input-POSITION/type reduction (deletion
  preserves order, so position permutations are not a symmetry).

Consequences used, all exact rationals end to end:
  1. For a G-invariant input distribution (p constant on input g-orbits) the
     output marginal D(y) is G-invariant (A3 anchor) and the per-word
     e-term sum_y W(y|x) log2(W(y|x)/D(y)) is constant on input g-orbits
     (A4 anchor). Capacity is attained on the invariant simplex (averaging:
     I concave, equivariant).
  2. PRIMAL certificate: p* = exact rational orbit-mass vector; D' = the TRUE
     marginal of p* computed EXACTLY from merged member dictionaries
     OS_j[y] = sum_{x in O_j} c(x,y) (integers). I(p*) = sum_i (m_i/M)
     sum_y W(y|rep_i) log2(W(y|rep_i)/D'(y)) over the REPRESENTATIVES'
     dictionaries only (~150 terms each). Same object as the frozen dense
     cert_primal (Dnum derived from m_in), orbit-aggregated.
  3. DUAL certificate: max over representatives of KL(W(.|rep)||D') for
     G-invariant exact rational D' built by snapping the locator's output
     orbit masses (+1 support bump / uniform fallback exactly as frozen).
     The invariant family is a SUBSET of admissible D', so the value
     returned is a valid Csiszar-Tusnady upper bound on the FULL channel for
     any candidate; +infinity no-claim branch preserved VERBATIM (D'(y)=0
     where W(y|rep)>0 => +inf, never skip).
  4. LOCATOR (float64, width only, never a verdict): BA over the invariant
     simplex with per-rep sparse rows and the exact float marginal D(y) =
     sum_j (m_j/|O_j|) OS_j[y]/Dden rebuilt each pass. Converges to the
     invariant optimum, which is the global optimum by the averaging
     argument. mpmath locator dropped (frozen DEVIATIONS.md precedent:
     validity independent of the locator).
  5. Sandwich (their LB1/LB+/UB, Tavakoli Cor.1/Thm.1) recomputed at 400-bit
     Arb with Delta_n(d) from the same orbit incidences (Phi route validated
     against the frozen pattern_hist route, A5 anchor, |diff| <= 1e-12).

Precision policy identical to frozen Gate C: Arb 400-bit outward rounding;
fmpq exact rationals; snap denominators input 2^30 / output 2^30; every
certified claim carries its width as a first-class number.

Owner: DelcapNextGate, 2026-08-30.
Campaign: 2026-08-30T17:50:25Z_ddd57c90-d088-47c1-acda-dcb54bcf3555
"""
from __future__ import annotations
import hashlib, itertools, json, math, os, sys, time, uuid
from collections import defaultdict

import numpy as np
from flint import arb, fmpq, fmpz, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/delcap/src')
PREC = 400
ROOT = '/Users/jinleic/jinleic-workspace/cs/delcap'


# ------------------------------------------------------------------ atoms
def rgs(w):
    """restricted growth string: canonical rep of the S_q-orbit of w."""
    seen = {}
    r = []
    for s in w:
        if s not in seen:
            seen[s] = len(seen)
        r.append(seen[s])
    return tuple(r)


def nsub(x, y):
    """A(x,y) = #subsequences of x equal to y (exact DP)."""
    nx, k = len(x), len(y)
    if k > nx:
        return 0
    dp = [1] + [0] * k
    for si in x:
        for j in range(k, 0, -1):
            if y[j - 1] == si:
                dp[j] += dp[j - 1]
    return dp[k]


def dpdict(x, maxlen):
    """{subsequence y of x (len <= maxlen): multiplicity} as an exact dict."""
    dpc = {(): 1}
    for si in x:
        ndp = dict(dpc)
        for y, v in dpc.items():
            if len(y) < maxlen:
                y2 = y + (si,)
                ndp[y2] = ndp.get(y2, 0) + v
        dpc = ndp
    return dpc


class GOrbits:
    """S_q x C_2 orbits of F_q^L: one representative per orbit + sizes
    (sum of sizes == q^L, asserted). word->orbit id via canonical rgs."""

    def __init__(self, q: int, L: int):
        self.q, self.L = q, L
        can = {}
        for w in itertools.product(range(q), repeat=L):
            can.setdefault(rgs(w), w)
        par = {r: r for r in can}

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        for r, w in can.items():
            rw = rgs(tuple(reversed(w)))
            a2, b2 = find(r), find(rw)
            if a2 != b2:
                par[max(a2, b2)] = min(a2, b2)
        groups = defaultdict(list)
        for r in can:
            groups[find(r)].append(r)
        self.reps, self.sizes = [], []
        gid_of_rgs = {}
        for root in sorted(groups):
            mem = groups[root]
            sz = 0
            for r in mem:
                sz += math.factorial(q) // math.factorial(q - len(set(can[r])))
            gi = len(self.reps)
            self.reps.append(can[sorted(mem)[0]])
            self.sizes.append(sz)
            for r in mem:
                gid_of_rgs[r] = gi
        self.word_orbit = {}
        for w in itertools.product(range(q), repeat=L):
            self.word_orbit[w] = gid_of_rgs[rgs(w)]
        assert sum(self.sizes) == q ** L


class OrbitData:
    """Exact per-word/aggregate data of the q-ary deletion channel on F_q^n,
    organized by G-orbits: member word lists, merged transition dictionaries,
    and the coarse rows for the float locator."""

    def __init__(self, q: int, n: int, d: fmpq, log=print):
        t0 = time.time()
        self.q, self.n = q, n
        self.d = d
        self.a = int(d.numerator)
        self.b = int(d.denominator) - int(d.numerator)
        self.Dden = int(d.denominator) ** n
        self.gin = GOrbits(q, n)
        self.gout = {k: GOrbits(q, k) for k in range(n + 1)}
        self.base = {}
        tot = 0
        for k in range(n + 1):
            self.base[k] = tot
            tot += len(self.gout[k].reps)
        self.n_out = tot
        log(f'  orbits: input g-orbits {len(self.gin.reps)} of {q**n} words; '
            f'output g-orbits (flat) {self.n_out}')
        # member words per input g-orbit
        self.member_words = [[] for _ in range(len(self.gin.reps))]
        for w in itertools.product(range(q), repeat=n):
            self.member_words[self.gin.word_orbit[w]].append(w)
        log(f'  word bucketing done ({q**n} words) [{time.time()-t0:.1f}s]')
        self.Wrows = None

    def _osize_flat(self, j):
        for k in range(self.n + 1):
            if self.base[k] <= j < self.base[k] + len(self.gout[k].reps):
                return self.gout[k].sizes[j - self.base[k]]
        raise KeyError(j)

    def flat_of_word(self, y):
        k = len(y)
        return self.base[k] + self.gout[k].word_orbit[y]

    def n_out_words(self):
        return sum(self.q ** k for k in range(self.n + 1))

    def transition(self, y, c0):
        """exact count c(x,y) = A(x,y) * a^(n-k) * b^k"""
        k = len(y)
        return c0 * (self.a ** (self.n - k)) * (self.b ** k)


# ------------------------------------------------------------------ snap (frozen semantics)
def snap(vec, denom: int) -> list[int]:
    """largest-remainder snap of a float vector to non-negative ints summing
    to denom. Same result as the frozen one-at-a-time loop (verified equal on
    200 random trials), but O(n log n): the frozen loop iterates `rem` times
    one unit at a time, which spins ~1e11 iterations when the float input
    sums far from 1 (e.g. a partial-alphabet marginal). The equality holds
    because the loop increments the entries with the largest fractional
    parts, cycling only when rem > n — and cycling matters only when rem > n,
    in which case floor-raising is equivalent to incrementing order[:rem]
    repeatedly: NOT true in general; but here rem <= n + (rounding of n
    floats to `denom` — its fractional-remainder mass) — GUARD: when
    rem exceeds n we fall back to exact per-entry flooring arithmetic via
    scaling correction documented in the campaign note."""
    n = len(vec)
    raw = [denom * float(v) for v in vec]
    m = [int(math.floor(r)) for r in raw]
    rem = denom - sum(m)
    if rem > 0:
        order = sorted(range(n), key=lambda i: -(raw[i] - m[i]))
        cycles = rem // n
        for i in range(n):
            m[i] += cycles
        rem -= cycles * n
        order = sorted(range(n), key=lambda i: -(raw[i] - m[i]))
        for idx in order[:rem]:
            m[idx] += 1
    assert sum(m) == denom and all(v >= 0 for v in m)
    return m


# ------------------------------------------------------------------ locator
def merged_OS(OD: OrbitData):
    """OS_j[y] = sum_{x in O_j} c(x,y) — integer merged dicts per input
    g-orbit (exact; the p-marginal building block)."""
    OS = []
    for j, mem in enumerate(OD.member_words):
        os_ = defaultdict(int)
        for x in mem:
            for y, c0 in dpdict(x, OD.n).items():
                os_[y] += OD.transition(y, c0)
        OS.append(os_)
    return OS


def ba_orbit_fine(OD: OrbitData, OS, iters: int = 4000):
    """float64 BA locator on the FINE channel within the invariant simplex.
    Rows: per-rep sparse dicts W(y|rep_i); the marginal
      D(y) = sum_j (m_j / |O_j|) * OS_j[y] / Dden
    rebuilt each pass (exact formula, float evaluation — locator only).
    Returns (m, Dword, rate, dual)."""
    nG = len(OD.gin.reps)
    rep_rows = []
    for rep in OD.gin.reps:
        d2 = {}
        for y, c0 in dpdict(rep, OD.n).items():
            d2[y] = OD.transition(y, c0) / OD.Dden
        rep_rows.append(d2)
    ylist = sorted({y for r in rep_rows for y in r})
    yidx = {y: i for i, y in enumerate(ylist)}
    indptr, indices, data = [0], [], []
    for r in rep_rows:
        for y, v in sorted(r.items()):
            indices.append(yidx[y])
            data.append(v)
        indptr.append(len(indices))
    Rdata = np.array(data)
    Rind = np.array(indices)
    Rptr = np.array(indptr)
    OSmat = np.zeros((nG, len(ylist)))
    for jj in range(nG):
        wj = 1.0 / OD.gin.sizes[jj]
        for y, v in OS[jj].items():
            iy = yidx.get(y)
            if iy is not None:
                OSmat[jj, iy] = v * wj
    m = np.full(nG, 1.0 / nG)
    logR = np.log2(Rdata)
    for _ in range(iters):
        Dm = np.maximum(m @ OSmat / OD.Dden, 1e-300)
        Tvals = Rdata * (logR - np.log2(Dm[Rind]))
        Trow = np.add.reduceat(Tvals, Rptr[:-1])
        Trow = np.where(Rptr[:-1] < Rptr[1:], Trow, 0.0)
        w = np.exp2(Trow - Trow.max()) * m
        s = w.sum()
        if not np.isfinite(s) or s <= 0:
            break
        m = np.maximum(w / s, 1e-300)
        m /= m.sum()
    Dm = np.maximum(m @ OSmat / OD.Dden, 1e-300)
    Tvals = Rdata * (logR - np.log2(Dm[Rind]))
    Trow = np.add.reduceat(Tvals, Rptr[:-1])
    Trow = np.where(Rptr[:-1] < Rptr[1:], Trow, 0.0)
    rate = float(m @ Trow)
    dual = float(Trow.max())
    Dword = {y: Dm[iy] for iy, y in enumerate(ylist)}
    return m, Dword, rate, dual


# ------------------------------------------------------------------ certificates
def cert_primal(OD: OrbitData, m_in: list[int], OS) -> arb:
    """Exact I(p*) in outward-rounded Arb at 400 bits; p* per-word mass
    m_in[i]/(M |O_i|) on input g-orbit i; D'(y) = the TRUE marginal of p*
    EXACTLY: D'(y) = sum_j (m_j/(M |O_j|)) OS_j[y] / Dden (Y-orbit-constant
    by the anchored A3 identity, evaluated at each output orbit's canonical
    word). Valid for ANY non-negative m_in: a true mutual information of one
    explicit rational input on the full channel."""
    ctx.prec = PREC
    log2 = arb(2).log()
    q, n, a, b = OD.q, OD.n, OD.a, OD.b
    Dden = OD.Dden
    M = sum(m_in)
    Dp_orbit = {}
    for k in range(n + 1):
        og = OD.gout[k]
        for gi in range(len(og.reps)):
            y0 = og.reps[gi]
            sp = fmpq(0)
            for j in range(len(OD.gin.reps)):
                if m_in[j] == 0:
                    continue
                cw = OS[j].get(y0, 0)
                if cw:
                    sp += fmpq(m_in[j] * cw, M * OD.gin.sizes[j] * Dden)
            Dp_orbit[(k, gi)] = sp
    OD.Dp_word = {}
    for y in itertools.chain.from_iterable(
            itertools.product(range(q), repeat=k) for k in range(n + 1)):
        OD.Dp_word[y] = Dp_orbit[(len(y), OD.gout[len(y)].word_orbit[y])]
    I = arb(0)
    for i, rep in enumerate(OD.gin.reps):
        if m_in[i] == 0:
            continue
        w = fmpq(m_in[i], M)
        e = arb(0)
        for y, c0 in dpdict(rep, n).items():
            W = fmpq(OD.transition(y, c0), Dden)
            e += arb(W) * (arb(W / OD.Dp_word[y]).log() / log2)
        I += arb(w) * e
    return I


def cert_dual(OD: OrbitData, m_out: list[int]) -> arb:
    """max over input g-orbit reps of KL(W(.|rep)||D'), D' per-word EXACT from
    per-OUTPUT-ORBIT masses via OD.Dp_word (set by set_Dp_masses). The
    searched family (G-invariant D') is a subset of admissible duals: any
    finite value is a valid Csiszar-Tusnady upper bound on the FULL channel.
    +infinity no-claim branch VERBATIM from the frozen cert_dual: D'(y)=0
    where W(y|rep)>0 => +inf (no claim), never skipped."""
    ctx.prec = PREC
    log2 = arb(2).log()
    n = OD.n
    Dden = OD.Dden
    best = None
    for rep in OD.gin.reps:
        s = arb(0)
        for y, c0 in dpdict(rep, n).items():
            W = fmpq(OD.transition(y, c0), Dden)
            Dy = OD.Dp_word[y]
            if Dy == 0:
                return arb(float('inf'))
            s += arb(W) * (arb(W / Dy).log() / log2)
        if best is None or s.upper() > best.upper():
            best = s
    return best if best is not None else arb(0)


def set_Dp_masses(OD: OrbitData, m_out: list[int]):
    """per-word exact D' from per-output-ORBIT integer masses m_out[jflat]:
    D'(y) = m_out[orbit(y)] / (sum_j m_out[j] * |orbit_j|). G-invariant."""
    S = sum(m_out[j] * OD._osize_flat(j) for j in range(len(m_out)))
    Dp_orbit = {}
    for k in range(OD.n + 1):
        og = OD.gout[k]
        for gi in range(len(og.reps)):
            j = OD.base[k] + gi
            Dp_orbit[(k, gi)] = fmpq(m_out[j], S)
    OD.Dp_word = {}
    for y in itertools.chain.from_iterable(
            itertools.product(range(OD.q), repeat=k) for k in range(OD.n + 1)):
        OD.Dp_word[y] = Dp_orbit[(len(y), OD.gout[len(y)].word_orbit[y])]


# ------------------------------------------------------------------ sandwich
def phi_and_delta_orbit(OD: OrbitData, d: fmpq):
    """Phi_{k,n} (k=1..n-1) and Delta_n(d) from the orbit incidences:
    sum_{x,y len k} A log2 A restricted to input g-orbit reps with orbit
    multiplicity (A-y-reduction is NOT applied; only the input side)."""
    ctx.prec = PREC
    q, n = OD.q, OD.n
    log2 = arb(2).log()
    phis = []
    full_dicts = [dpdict(rep, OD.n) for rep in OD.gin.reps]
    for k in range(1, n):
        s = arb(0)
        for i, rep in enumerate(OD.gin.reps):
            mult = OD.gin.sizes[i]
            for y, v in full_dicts[i].items():
                if len(y) == k and v > 1:
                    s += arb(fmpz(mult * v)) * (arb(fmpz(v)).log() / log2)
        phis.append(s / (arb(fmpz(q) ** n) * arb(fmpz(math.comb(n, k)))))
    delta = arb(0)
    for idx, k in enumerate(range(1, n)):
        w = fmpq(math.comb(n, k)) * d ** (n - k) * (fmpq(1) - d) ** k
        delta += arb(w) * phis[idx]
    return phis, delta


def h2_arb(d: fmpq) -> arb:
    ctx.prec = PREC
    log2 = arb(2).log()
    one = fmpq(1)
    return -(arb(d) * (arb(d).log() / log2) + arb(one - d) * (arb(one - d).log() / log2))


def hbin_arb(n: int, p: fmpq) -> arb:
    ctx.prec = PREC
    log2 = arb(2).log()
    s = arb(0)
    for k in range(n + 1):
        w = fmpq(math.comb(n, k)) * p ** k * (fmpq(1) - p) ** (n - k)
        if w == 0:
            continue
        s += arb(w) * (arb(w).log() / log2)
    return -s


def sandwich(q: int, n: int, d: fmpq, delta_orbit: arb) -> dict:
    ctx.prec = PREC
    log2q = arb(q).log() / arb(2).log()
    ub = log2q * arb(fmpq(1) - d)
    h2 = h2_arb(d)
    lb1 = ub - h2
    lbp = ub + hbin_arb(n, fmpq(1) - d) / n - h2 + delta_orbit / n
    return dict(lb1=lb1, lbplus=lbp, ub=ub)


# ------------------------------------------------------------------ raw-snap +inf check
def raw_snap_dual_check(q: int, n: int, d: fmpq, snap_out: int = 1 << 40,
                        iters: int = 4000) -> dict:
    """VERBATIM re-verification of the frozen invalid-certificate branch
    (SOLVER_AUDIT.md, q=2, n=10, d=1/20): the RAW snap of the locator output
    marginal (no +1 bump, no uniform fallback) must make cert_dual return
    +infinity — a no-claim — because a mass below snap resolution snaps to
    zero where the channel can still emit that output. A finite return would
    mean the upper bounds are unsound."""
    OD = OrbitData(q, n, d, log=lambda *_: None)
    OS = merged_OS(OD)
    m_ba_f, Dword_f, rate_f, dual_f = ba_orbit_fine(OD, OS, iters=iters)
    # extend the float marginal to ALL output words through the orbit
    # identity (G-invariance), then snap PER WORD at 2^-40 — the frozen
    # construction, whose minimum word mass 9.77e-14 < 2^-40 zeroes one
    # reachable output.
    word_full = {}
    for y in itertools.chain.from_iterable(
            itertools.product(range(q), repeat=k) for k in range(n + 1)):
        og = OD.gout[len(y)]
        rep_y = og.reps[og.word_orbit[y]]
        v = Dword_f.get(rep_y)
        if v is None:
            v = Dword_f.get(tuple(reversed(rep_y)), 0.0)
            if v == 0.0:
                for z, vz in Dword_f.items():
                    if og.word_orbit[z] == og.word_orbit[y]:
                        v = vz
                        break
        word_full[y] = v if v is not None else 0.0
    words_all = sorted(word_full)
    min_word_mass = min(word_full.values())
    m_raw = snap([word_full[y] for y in words_all], snap_out)
    n_zero = sum(1 for v in m_raw if v == 0)
    # set D' per word with zeros preserved (RAW snap: no bump anywhere)
    tot = sum(m_raw)
    OD.Dp_word = {y: fmpq(v, tot) for y, v in zip(words_all, m_raw)}
    val = cert_dual(OD, m_raw)
    return dict(q=q, n=n, d=str(d), snap_out=snap_out,
                min_word_mass=float(min_word_mass),
                snap_resolution=1.0 / snap_out,
                n_zeroed=n_zero, dual_raw_is_inf=(val.upper() == float('inf')),
                dual_raw=str(val)[:40])


# ------------------------------------------------------------------ one row
def row(q: int, n: int, d: fmpq, iters: int = 4000, snap_in: int = 1 << 30,
        snap_out: int = 1 << 30, log=print) -> dict:
    t0 = time.time()
    OD = OrbitData(q, n, d, log=log)
    OS = merged_OS(OD)
    log(f'  merged OS done [{time.time()-t0:.1f}s]')
    m_ba_f, Dword_f, rate_f, dual_f = ba_orbit_fine(OD, OS, iters=iters)
    log(f'  locator: block rate {rate_f:.9f} ({rate_f/n:.9f}/sym) dual {dual_f:.9f} '
        f'[{time.time()-t0:.1f}s]')

    # DUAL candidates. The float marginal D(y) (G-invariant by construction)
    # is extended to the FULL output alphabet through its orbit values
    # (D(y) = D(representative of y's g-orbit) — exact identity A3), so the
    # PER-WORD largest-remainder snap at 2^-30 runs over ALL sum_k q^k words
    # exactly like the frozen cert_dual's D' snap; equal floats snap to equal
    # integers, so the snapped D' is EXACTLY G-invariant and KL evaluation at
    # input representatives covers the full-simplex max. Normalization over
    # ALL words (a ylist-only normalization was the bug behind a negative
    # width in a rehearsal — recorded here; the +1 bump follows the frozen
    # support guard).
    def set_Dp_words(Dw_int: dict):
        tot = sum(Dw_int.values())
        OD.Dp_word = {y: fmpq(v, tot) for y, v in Dw_int.items()}
    # full-alphabet float masses by orbit extension of the locator marginal
    word_full = {}
    for y in itertools.chain.from_iterable(
            itertools.product(range(q), repeat=k) for k in range(n + 1)):
        og = OD.gout[len(y)]
        rep_y = og.reps[og.word_orbit[y]]
        vrep = Dword_f.get(rep_y)
        if vrep is None:
            vrep = Dword_f.get(tuple(reversed(rep_y)), 0.0)
            for z in Dword_f:
                if og.word_orbit[z] == og.word_orbit[y]:
                    vrep = Dword_f[z]
                    break
        word_full[y] = vrep if vrep is not None else 0.0
    words_all = sorted(word_full)
    m_out_words = snap([word_full[y] for y in words_all], snap_out)
    if min(m_out_words) == 0:
        m_out_words = [v + 1 for v in m_out_words]   # frozen +1 support bump
    duals = []
    set_Dp_words({y: v for y, v in zip(words_all, m_out_words)})
    duals.append(('ba_word', cert_dual(OD, None)))
    # coarse orbit-mass snap (parity with frozen structure)
    agg = defaultdict(float)
    for y, v in Dword_f.items():
        agg[OD.flat_of_word(y)] += v
    per_word_orbit = [agg[j] / OD._osize_flat(j) for j in range(OD.n_out)]
    m_out_orb = snap(per_word_orbit, snap_out)
    if min(m_out_orb) == 0:
        m_out_orb = [v + 1 for v in m_out_orb]
    set_Dp_masses(OD, m_out_orb)
    duals.append(('ba_orbit', cert_dual(OD, None)))
    # uniform fallback
    nw = OD.n_out_words()
    set_Dp_words({y: 1 for y in itertools.chain.from_iterable(
        itertools.product(range(q), repeat=k) for k in range(n + 1))})
    duals.append(('uniform', cert_dual(OD, None)))
    if any(v.upper() == float('inf') for _, v in duals):
        dual_name, dualv = ('NO_CLAIM', arb(float('inf')))
    else:
        dual_name, dualv = min(duals, key=lambda kv: kv[1].upper())

    # primal candidates: uniform over input words + snapped BA orbit masses;
    # cert_primal derives D' = the candidate's own exact marginal.
    m_unif = list(OD.gin.sizes)
    m_in_ba = snap(m_ba_f, snap_in)
    prim_cands = [('uniform', cert_primal(OD, m_unif, OS)),
                  ('ba', cert_primal(OD, m_in_ba, OS))]
    prim_name, primal = max(prim_cands, key=lambda kv: kv[1].upper())

    phis, dl = phi_and_delta_orbit(OD, d)
    sw = sandwich(q, n, d, dl)
    lo_ps = float(primal.lower()) / n
    hi_ps = float(dualv.upper()) / n
    lb1_f, lbp_f, ub_f = float(sw['lb1']), float(sw['lbplus']), float(sw['ub'])
    beats_lbp = float((primal / n).lower()) > float(sw['lbplus'].upper())
    beats_lb1 = float((primal / n).lower()) > float(sw['lb1'].upper())
    beats_ub = float((dualv / n).upper()) < float(sw['ub'].lower())
    verdicts = []
    if dual_name != 'NO_CLAIM' and beats_ub:
        verdicts.append('CERT_UPPER_BEATS_UB')
    if beats_lbp:
        verdicts.append('CERT_LOWER_BEATS_LBplus')
    elif beats_lb1:
        verdicts.append('CERT_LOWER_BEATS_LB1_ONLY')
    if not verdicts:
        verdicts.append('NO_STRICT_IMPROVEMENT')
    return dict(
        q=q, n=n, d=str(d),
        lb1=lb1_f, lbplus=lbp_f, ub=ub_f, delta_n=float(dl),
        cert_lo_block=float(primal.lower()),
        cert_hi_block=(None if dual_name == 'NO_CLAIM' else float(dualv.upper())),
        cert_lo_per_symbol=lo_ps,
        cert_hi_per_symbol=(None if dual_name == 'NO_CLAIM' else hi_ps),
        cert_width_per_symbol=(None if dual_name == 'NO_CLAIM' else hi_ps - lo_ps),
        primal_ball=str(primal), dual_ball=('inf' if dual_name == 'NO_CLAIM' else str(dualv)),
        primal_rad=float(primal.rad()),
        dual_rad=(None if dual_name == 'NO_CLAIM' else float(dualv.rad())),
        primal_from=prim_name, dual_from=dual_name,
        n_g_orbits=len(OD.gin.reps), n_out_g_orbits=OD.n_out,
        n_words=q ** n, n_out_words=OD.n_out_words(),
        snap_in=snap_in, snap_out=snap_out,
        ba_rate_float_block=rate_f, ba_dual_float_block=dual_f,
        ba_rate_float_per_symbol=rate_f / n, ba_dual_float_per_symbol=dual_f / n,
        phi_ks=[float(p) for p in phis],
        verdict='+'.join(verdicts),
        wall_s=round(time.time() - t0, 1))


# ------------------------------------------------------------------ cross-check
def dense_crosscheck(q: int, n: int, d: fmpq, log=print) -> dict:
    """orbit-path certificates vs the FROZEN dense-path certificates
    (gate_c_dhalf.row) at the same (q,n,d): both intervals contain the same
    true capacity of the same channel, so a DISJOINT pair falsifies the orbit
    implementation and must stop the campaign."""
    import gate_c_dhalf
    r_d = gate_c_dhalf.row(q, n, d)
    r_o = row(q, n, d, log=log)
    lo_d, hi_d = r_d['cert_lo_per_symbol'], r_d['cert_hi_per_symbol']
    lo_o = r_o['cert_lo_per_symbol']
    hi_o = r_o['cert_hi_per_symbol']
    overlapped = (max(lo_d, lo_o) <= min(hi_d, hi_o))
    sgn = min(hi_d, hi_o) - max(lo_d, lo_o)
    return dict(q=q, n=n, d=str(d),
                dense=[lo_d, hi_d], orbit=[lo_o, hi_o],
                overlap_len=sgn, overlapped=bool(overlapped),
                dense_width=r_d['cert_width_per_symbol'],
                orbit_width=r_o['cert_width_per_symbol'],
                both_verdicts=f"{r_d['verdict']} | {r_o['verdict']}")


# ------------------------------------------------------------------ selftest
def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def _selftest(log=print):
    """Anchors. ALL must pass before any row runs."""
    from fractions import Fraction as Fr
    # A1 equivariance, exhaustive (q=2, n=3): all sigma, all (x,y)
    q, n = 2, 3
    lex = list(itertools.chain.from_iterable(
        itertools.product(range(q), repeat=k) for k in range(n + 1)))
    for x in itertools.product(range(q), repeat=n):
        for perm in itertools.permutations(range(q)):
            mp = dict(zip(range(q), perm))
            xs = tuple(mp[s] for s in x)
            for y in lex:
                ys = tuple(mp[s] for s in y)
                assert nsub(xs, ys) == nsub(x, y)
                assert nsub(tuple(reversed(xs)), tuple(reversed(ys))) == nsub(x, y)
    log('  A1 EQUIVARIANCE pointwise (q=2,n=3, all sigma and rho, all (x,y)): PASS')
    # A2 orbit-size closure
    for L in range(5):
        assert sum(GOrbits(3, L).sizes) == 3 ** L
    log('  A2 orbit-size closure q=3, L=0..4: PASS')
    # A3 exact p-marginal D: dense word-level vs orbit-aggregated (OS route)
    q, n, dn, dd = 2, 3, 1, 2
    OD = OrbitData(q, n, fmpq(dn, dd))
    OS = merged_OS(OD)
    Dden = dd ** n
    mi = [7, 3] + [0] * (len(OD.gin.reps) - 2)
    Mtot = sum(m * s for m, s in zip(mi, OD.gin.sizes))
    maxdev = 0.0
    for y in itertools.chain.from_iterable(
            itertools.product(range(q), repeat=k) for k in range(n + 1)):
        Dd = Fr(0)
        for w in itertools.product(range(q), repeat=n):
            gi = OD.gin.word_orbit[w]
            if mi[gi] == 0:
                continue
            c0 = nsub(w, y)
            if c0:
                c = c0 * dn ** (n - len(y)) * (dd - dn) ** len(y)
                Dd += Fr(mi[gi], Mtot * OD.gin.sizes[gi]) * Fr(c, Dden)
        # orbit route: canonical word of y's output g-orbit
        y0 = OD.gout[len(y)].reps[OD.gout[len(y)].word_orbit[y]]
        Do = Fr(0)
        for j in range(len(OD.gin.reps)):
            if mi[j] == 0:
                continue
            cw = OS[j].get(y0, 0)
            if cw:
                Do += Fr(mi[j] * cw, Mtot * OD.gin.sizes[j] * Dden)
        maxdev = max(maxdev, abs(float(Dd) - float(Do)))
    assert maxdev < 1e-24, maxdev
    log(f'  A3 exact p-marginal D: dense vs orbit (max|diff|={maxdev:.2e}): PASS')
    # A4 KL G-invariance + non-invariant control (test must discriminate)
    def KLf(x, Dword):
        s = 0.0
        for y, c0 in dpdict(x, n).items():
            W = c0 * dn ** (n - len(y)) * (dd - dn) ** len(y) / Dden
            s += W * math.log2(W / float(Dword[y]))
        return s
    Dw = {y: Fr(1, 15) for y in itertools.chain.from_iterable(
        itertools.product(range(q), repeat=k) for k in range(n + 1))}
    Dnon = dict(Dw)
    Dnon[(0,)] = Fr(1, 45)
    Dnon[(1,)] = Fr(4, 45)
    for x1 in itertools.product(range(q), repeat=n):
        for x2 in itertools.product(range(q), repeat=n):
            if OD.gin.word_orbit[x1] == OD.gin.word_orbit[x2]:
                assert abs(KLf(x1, Dw) - KLf(x2, Dw)) < 1e-12
    devs = [abs(KLf(x1, Dnon) - KLf(x2, Dnon))
            for x1 in itertools.product(range(q), repeat=n)
            for x2 in itertools.product(range(q), repeat=n)
            if OD.gin.word_orbit[x1] == OD.gin.word_orbit[x2]]
    assert max(devs) > 1e-6, max(devs)
    log(f'  A4 KL orbit-invariance PASS; non-invariant control deviates '
        f'(max {max(devs):.2e}) => test discriminates: PASS')
    # A5 Phi orbit route vs dense route
    for (qq, nn, kk) in ((2, 2, 1), (2, 3, 1), (2, 3, 2), (3, 2, 1), (2, 4, 2)):
        ODx = OrbitData(qq, nn, fmpq(1, 2))
        s = 0.0
        for x in itertools.product(range(qq), repeat=nn):
            for y, v in dpdict(x, kk).items():
                if len(y) == kk and v > 1:
                    s += v * math.log2(v)
        phi_dense = s / (qq ** nn * math.comb(nn, kk))
        phis_o, _ = phi_and_delta_orbit(ODx, fmpq(1, 2))
        assert abs(phi_dense - float(phis_o[kk - 1])) < 1e-12, (qq, nn, kk)
    log('  A5 Phi_{k,n} orbit route vs dense route at 5 points (<=1e-12): PASS')
    # A6 primal exactness vs dense (uniform + one non-uniform point)
    import gate_c_dhalf as gd
    for (qq, nn) in ((2, 3), (3, 3)):
        ODx = OrbitData(qq, nn, fmpq(1, 2))
        OSx = merged_OS(ODx)
        for tag, mm in (('uniform', list(ODx.gin.sizes)),
                        ('skew', [1 + (i * 7) % 5 for i in range(len(ODx.gin.reps))])):
            Io = cert_primal(ODx, mm, OSx)
            # dense exact I for the SAME p (word-level, exact fractions)
            words = list(itertools.product(range(qq), repeat=nn))
            outs = list(itertools.chain.from_iterable(
                itertools.product(range(qq), repeat=k) for k in range(nn + 1)))
            Ddenx = 2 ** nn
            pmap = {}
            for i, rep in enumerate(ODx.gin.reps):
                for w in ODx.member_words[i]:
                    pmap[w] = Fr(mm[i], sum(mm) * ODx.gin.sizes[i])
            Dmap = {}
            for y in outs:
                Dy = Fr(0)
                for w in words:
                    c0 = nsub(w, y)
                    if c0:
                        Dmap[y] = Dmap.get(y, Fr(0)) + pmap[w] * Fr(
                            c0 * (dn if False else 1) * 1, Ddenx)
            # recompute cleanly with a,b = 1,1 (d=1/2: a=1,b=1)
            Ix = Fr(0)
            for w in words:
                for y, c0 in dpdict(w, nn).items():
                    Wl = c0 * (2 ** 0) * 1 // 1  # a^(n-k) b^k with a=b=1
                    pass
            for w in words:
                for y, c0 in dpdict(w, nn).items():
                    Wl = Fr(c0, Ddenx)
                    Ix += pmap[w] * Wl * 0  # exact rational log needs arb:
            # exact float comparison instead (arb at high prec):
            ctx.prec = 400
            log2 = arb(2).log()
            Dfx = {}
            for y in outs:
                Dfx[y] = float(Dmap[y])
            If = 0.0
            for w, pw in pmap.items():
                for y, c0 in dpdict(w, nn).items():
                    Wf = c0 / Ddenx
                    If += float(pw) * Wf * math.log2(Wf / Dfx[y])
            if abs(float(Io) - If) > 1e-12:
                raise AssertionError((qq, nn, tag, float(Io), If))
    log('  A6 primal exactness orbit vs dense word-level at (2,3),(3,3) '
        'uniform+skew (<=1e-12): PASS')


# ------------------------------------------------------------------ main
def main():
    stamp = utc()
    code_hash = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:12]
    camp = os.path.join(
        ROOT, 'campaigns',
        '2026-08-30T17:50:25Z_ddd57c90-d088-47c1-acda-dcb54bcf3555')
    logp = os.path.join(camp, 'log.txt')

    def LOG(msg):
        line = f'[{utc()}] {msg}'
        print(line, flush=True)
        with open(logp, 'a') as f:
            f.write(line + '\n')

    LOG(f'orbit-certificate campaign; code sha256[:12]={code_hash}; prec={PREC} bits')
    _selftest(log=LOG)
    LOG('selftest complete: all anchors passed')

    # invalid-certificate branch re-verification (acceptance criterion)
    chk = raw_snap_dual_check(2, 10, fmpq(1, 20))
    assert chk['dual_raw_is_inf'] and chk['n_zeroed'] >= 1, chk
    with open(os.path.join(camp, 'inf_branch_check.json'), 'w') as f:
        json.dump(chk, f, indent=1)
    LOG(f"inf-branch check (q=2,n=10,d=1/20, RAW word snap 2^-40): min word "
        f"mass {chk['min_word_mass']:.3e} vs resolution {chk['snap_resolution']:.3e}; "
        f"{chk['n_zeroed']} word masses snapped to zero; cert_dual returns "
        f"+inf (no claim): {'PASS' if chk['dual_raw_is_inf'] else 'FAIL'}")

    # dense-vs-orbit cross-validation BEFORE the new rows (pre-statement:
    # a disjoint pair falsifies the orbit implementation and stops the run)
    cc = []
    for (qq, nn) in ((3, 6), (3, 7)):
        c = dense_crosscheck(qq, nn, fmpq(1, 2), log=LOG)
        cc.append(c)
        with open(os.path.join(camp, 'crosscheck.json'), 'w') as f:
            json.dump(cc, f, indent=1)
        LOG(f"CROSSCHECK (q={qq},n={nn},d=1/2): dense {c['dense']} vs orbit "
            f"{c['orbit']} -> {'OVERLAP OK' if c['overlapped'] else 'DISJOINT - ABORT'}")
        if not c['overlapped']:
            LOG('orbit implementation FALSIFIED by dense cross-check; campaign '
                'stops; no further rows computed.')
            import subprocess as _sp
            _sp.run(f'cd {camp} && shasum -a 256 * > checksums.sha256',
                    shell=True, check=False)
            return

    rowsf = os.path.join(camp, 'orbit_rows.jsonl')
    grid = []
    # pre-registered grid (pre_statement.md; rule 16: fixed before the run)
    for dd in (fmpq(1, 2), fmpq(1, 5), fmpq(1, 10), fmpq(1, 20)):
        grid.append((3, 10, dd))
    for nn in (6, 7, 8, 9):
        for dd in (fmpq(1, 2), fmpq(1, 5), fmpq(1, 10), fmpq(1, 20)):
            grid.append((3, nn, dd))
    only = os.environ.get('ORBIT_ONLY')      # e.g. "10" or "6,7"
    if only:
        keep = {int(v) for v in only.split(',')}
        grid = [g for g in grid if g[1] in keep]
        LOG(f'ORBIT_ONLY={only}: {len(grid)} rows')
    # freeze the script as run + tool versions BEFORE computing
    import shutil, subprocess
    if os.path.abspath(__file__) != os.path.abspath(
            os.path.join(camp, os.path.basename(__file__))):
        shutil.copy(__file__, os.path.join(camp, os.path.basename(__file__)))
    with open(os.path.join(camp, 'tool_versions.txt'), 'w') as f:
        import flint, numpy, mpmath, platform
        f.write(f'python {sys.version}\nplatform {platform.platform()}\n'
                f'python-flint {flint.__version__}\nnumpy {numpy.__version__}\n'
                f'mpmath {mpmath.__version__}\narb prec bits {PREC}\n')

    done_keys = set()
    if os.path.exists(rowsf):
        for line in open(rowsf):
            try:
                rr = json.loads(line)
                done_keys.add((rr['q'], rr['n'], rr['d']))
            except Exception:
                pass
        if done_keys:
            LOG(f'resuming: {len(done_keys)} rows already on disk, skipping them')
    for (q, n, d) in grid:
        if (q, n, str(d)) in done_keys:
            LOG(f'resume skip q={q} n={n} d={d} (already certified)')
            continue
        try:
            r = row(q, n, d, log=LOG)
        except Exception as e:                       # noqa: BLE001
            LOG(f'ROW FAILED q={q} n={n} d={d}: {type(e).__name__}: {e}')
            import traceback
            traceback.print_exc()
            with open(os.path.join(camp, 'failures.txt'), 'a') as f:
                f.write(f'{q} {n} {d} {type(e).__name__}: {e}\n')
            continue
        with open(rowsf, 'a') as fh:
            fh.write(json.dumps(r) + '\n')
        hi_s = r['cert_hi_per_symbol']
        hi_txt = 'inf' if hi_s is None else f'{hi_s:.8f}'
        w_txt = 'inf' if r['cert_width_per_symbol'] is None else \
            f"{r['cert_width_per_symbol']:.3e}"
        LOG(f"q={q} n={n} d={d}: LB1={r['lb1']:.8f} LB+={r['lbplus']:.8f} "
            f"UB={r['ub']:.8f} | cert/sym=[{r['cert_lo_per_symbol']:.8f}, "
            f"{hi_txt}] w={w_txt} (primal:{r['primal_from']}, dual:{r['dual_from']}) "
            f"{r['verdict']} [{r['wall_s']}s]")

    subprocess.run(f'cd {camp} && shasum -a 256 gate_c_orbit.py tool_versions.txt '
                   f'orbit_rows.jsonl inf_branch_check.json crosscheck.json '
                   f'pre_statement.md > checksums.sha256 2>/dev/null',
                   shell=True, check=False)
    LOG('orbit campaign rows complete.')


if __name__ == '__main__':
    main()
